#!/usr/bin/env python3
"""Freeze GSE143157 cohorts and stream-copy the paired primary Clean-77 HDF5."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from step2b_core import (
    as_bool,
    build_manifest,
    cohort_counts,
    subset_array_metadata,
    subset_donor_metadata,
)


VERSION = "1.0.0"
PRIMARY_H5_NAME = "GSE143157_primary_clean77.h5"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="GSE143157 Step 2B: freeze Full-78/Clean-77 cohorts and build primary HDF5",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--step1-dir", required=True, type=Path)
    parser.add_argument("--step2-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--chunk-rows", type=int, default=5000)
    parser.add_argument("--gzip-level", type=int, default=4)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate CSV/JSON evidence and print cohort counts without writing files",
    )
    return parser.parse_args()


def require_files(root: Path, names: list[str], label: str) -> dict[str, Path]:
    paths = {name: root / name for name in names}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise ValueError(f"{label} files are missing: {missing}")
    return paths


def decode_strings(values) -> list[str]:
    return [item.decode("utf-8") if isinstance(item, bytes) else str(item) for item in values]


def hash_small_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_ready(value):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def prepare_metadata(
    step1: dict[str, Path], step2: dict[str, Path]
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame], dict[str, object]]:
    donors = pd.read_csv(step1["metadata_donors_aligned.csv"])
    arrays = pd.read_csv(step1["metadata_arrays_aligned.csv"])
    outliers = pd.read_csv(step2["outlier_review.csv"])
    step1_summary = json.loads(step1["qc_summary.json"].read_text(encoding="utf-8"))
    step2_summary = json.loads(step2["step2_summary.json"].read_text(encoding="utf-8"))
    if step2_summary.get("donors") != 78:
        raise ValueError("Step 2 summary does not report 78 donors.")
    if step2_summary.get("review_candidate_count") != int(
        outliers["overall_review_candidate"].map(as_bool).sum()
    ):
        raise ValueError("Step 2 summary and outlier_review.csv candidate counts disagree.")

    manifest = build_manifest(donors, outliers)
    full78 = subset_donor_metadata(
        donors, manifest, "include_sensitivity_full78", "sensitivity_full78_column"
    )
    clean77 = subset_donor_metadata(
        donors, manifest, "include_primary_clean77", "primary_clean77_column"
    )
    no_two76 = subset_donor_metadata(
        donors, manifest, "include_sensitivity_no2_76", "sensitivity_no2_76_column"
    )
    tables = {
        "metadata_donors_sensitivity_full78.csv": full78,
        "metadata_donors_primary_clean77.csv": clean77,
        "metadata_donors_sensitivity_no_MSA30_Control11_76.csv": no_two76,
        "metadata_arrays_sensitivity_full78.csv": subset_array_metadata(arrays, full78),
        "metadata_arrays_primary_clean77.csv": subset_array_metadata(arrays, clean77),
        "metadata_arrays_sensitivity_no_MSA30_Control11_76.csv": subset_array_metadata(
            arrays, no_two76
        ),
    }
    context = {
        "step1_summary": step1_summary,
        "step2_summary": step2_summary,
        "outliers": outliers,
        "cohort_counts": {
            "sensitivity_full78": cohort_counts(full78),
            "primary_clean77": cohort_counts(clean77),
            "sensitivity_no_MSA30_Control11_76": cohort_counts(no_two76),
        },
    }
    return manifest, tables, context


def copy_primary_h5(
    source_path: Path,
    destination_path: Path,
    donor_metadata: pd.DataFrame,
    chunk_rows: int,
    gzip_level: int,
) -> dict[str, object]:
    try:
        import h5py
    except ImportError as exc:
        raise ValueError("h5py is not installed. Run the Windows runner so dependencies are installed.") from exc

    keep_indices = donor_metadata["original_hdf5_column"].astype(int).to_numpy()
    expected_ids = donor_metadata["donor_id"].astype(str).tolist()
    partial = destination_path.with_suffix(destination_path.suffix + ".partial")
    if partial.exists():
        raise ValueError(f"Incomplete temporary HDF5 already exists: {partial}")

    matrix_stats: dict[str, object] = {}
    with h5py.File(source_path, "r") as source:
        required = ["donor_ids", "5mC/values", "5mC/cpg_ids", "5hmC/values", "5hmC/cpg_ids"]
        missing = [name for name in required if name not in source]
        if missing:
            raise ValueError(f"Source HDF5 is missing datasets: {missing}")
        source_ids = decode_strings(source["donor_ids"][:])
        if len(source_ids) != 78:
            raise ValueError(f"Source HDF5 has {len(source_ids)} donors instead of 78.")
        if [source_ids[index] for index in keep_indices] != expected_ids:
            raise ValueError("Primary metadata order does not match source HDF5 donor order.")

        with h5py.File(partial, "w", libver="latest") as dest:
            dest.attrs["pipeline"] = "GSE143157 Step 2B cohort freeze"
            dest.attrs["pipeline_version"] = VERSION
            dest.attrs["cohort"] = "primary_clean77"
            dest.attrs["source_hdf5"] = str(source_path.resolve())
            dest.attrs["excluded_primary_donor"] = "MSA 30"
            dest.attrs["construction"] = "Column subset only; CpG rows and values otherwise unchanged"
            dest.create_dataset(
                "donor_ids",
                data=np.asarray(expected_ids, dtype=object),
                dtype=h5py.string_dtype(encoding="utf-8"),
            )
            dest.create_dataset("source_donor_indices", data=keep_indices.astype(np.int64))

            for modality in ["5mC", "5hmC"]:
                src_group = source[modality]
                src_values = src_group["values"]
                if src_values.shape[1] != 78:
                    raise ValueError(f"{modality} source matrix does not have 78 donor columns.")
                group = dest.create_group(modality)
                for key, value in src_group.attrs.items():
                    group.attrs[key] = value
                source.copy(src_group["cpg_ids"], group, name="cpg_ids")

                n_rows = int(src_values.shape[0])
                h5_chunk_rows = min(max(1, chunk_rows), n_rows)
                out_values = group.create_dataset(
                    "values",
                    shape=(n_rows, len(keep_indices)),
                    dtype=src_values.dtype,
                    chunks=(h5_chunk_rows, len(keep_indices)),
                    compression="gzip",
                    compression_opts=gzip_level,
                    shuffle=True,
                )
                for key, value in src_values.attrs.items():
                    out_values.attrs[key] = value

                missing_cells = 0
                finite_min = np.inf
                finite_max = -np.inf
                for start in range(0, n_rows, chunk_rows):
                    stop = min(start + chunk_rows, n_rows)
                    block = np.asarray(src_values[start:stop, keep_indices])
                    out_values[start:stop, :] = block
                    finite = np.isfinite(block)
                    missing_cells += int((~finite).sum())
                    if finite.any():
                        finite_min = min(finite_min, float(block[finite].min()))
                        finite_max = max(finite_max, float(block[finite].max()))
                    if stop % 100000 < chunk_rows or stop == n_rows:
                        print(f"[{modality}] copied {stop:,}/{n_rows:,} CpGs", flush=True)

                total_cells = n_rows * len(keep_indices)
                matrix_stats[modality] = {
                    "shape": [n_rows, len(keep_indices)],
                    "dtype": str(src_values.dtype),
                    "missing_cells": missing_cells,
                    "missing_fraction": missing_cells / total_cells,
                    "finite_min": None if not np.isfinite(finite_min) else finite_min,
                    "finite_max": None if not np.isfinite(finite_max) else finite_max,
                }
            dest.flush()

    os.replace(partial, destination_path)
    return matrix_stats


def write_html_report(path: Path, decision: dict[str, object]) -> None:
    counts = decision["cohorts"]
    matrices = decision["primary_hdf5"]["matrix_stats"]
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>GSE143157 Step 2B cohort report</title>
<style>
body{{font-family:Arial,sans-serif;max-width:1000px;margin:28px auto;padding:0 18px;color:#222;line-height:1.55}}
h1,h2{{color:#17365d}} .good{{background:#eaf6ec;border-left:5px solid #2e7d32;padding:12px 16px}}
.warn{{background:#fff4ce;border-left:5px solid #d6a800;padding:12px 16px}}
table{{border-collapse:collapse;width:100%}}th,td{{padding:8px;border-bottom:1px solid #ddd;text-align:left}}
code{{background:#f4f4f4;padding:2px 4px}}
</style></head><body>
<h1>GSE143157 Step 2B: frozen analysis cohorts</h1>
<p class="good"><strong>Primary paired cohort:</strong> 77 donors (37 Control, 40 MSA). MSA 30 is removed from both 5mC and 5hmC so the modality comparison uses the same people.</p>
<p class="warn"><strong>Control 11 remains in the primary cohort.</strong> PCA evidence alone is not treated as proof of technical failure. A separate 76-donor stress-test definition is supplied.</p>
<h2>Cohorts</h2><table><tr><th>Cohort</th><th>Donors</th><th>Control</th><th>MSA</th><th>Purpose</th></tr>
<tr><td>Primary Clean-77</td><td>{counts['primary_clean77']['donors']}</td><td>37</td><td>40</td><td>Main paired 5mC/5hmC analysis</td></tr>
<tr><td>Full-78</td><td>{counts['sensitivity_full78']['donors']}</td><td>37</td><td>41</td><td>Full-cohort sensitivity analysis</td></tr>
<tr><td>No-two 76</td><td>{counts['sensitivity_no_MSA30_Control11_76']['donors']}</td><td>36</td><td>40</td><td>Optional stress-test only</td></tr></table>
<h2>Why MSA 30 is excluded from the paired primary cohort</h2>
<ul><li>BS-array Detection-P failure exceeded the documented threshold.</li><li>5hmC missingness and negative BS−oxBS differences were robust QC outliers.</li><li>5hmC PCA was an extreme outlier, while 5mC was not.</li></ul>
<p>This rule is based on converging technical evidence, not diagnosis or PCA distance alone.</p>
<h2>Primary HDF5</h2>
<p><code>{html.escape(decision['primary_hdf5']['file'])}</code></p>
<ul><li>5mC shape: {matrices['5mC']['shape'][0]:,} × {matrices['5mC']['shape'][1]}</li>
<li>5hmC shape: {matrices['5hmC']['shape'][0]:,} × {matrices['5hmC']['shape'][1]}</li>
<li>No CpG row was removed or recalculated in Step 2B; only the MSA 30 donor column was omitted.</li></ul>
<h2>Analysis rule</h2><p>Use Clean-77 for the main paired comparison. Repeat final models on Full-78. The 76-donor definition is a stress-test and must not silently replace the primary cohort.</p>
</body></html>"""
    path.write_text(document, encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    started = time.time()
    if args.chunk_rows < 1:
        raise ValueError("--chunk-rows must be positive.")
    if not 0 <= args.gzip_level <= 9:
        raise ValueError("--gzip-level must be between 0 and 9.")

    step1_dir = args.step1_dir.resolve()
    step2_dir = args.step2_dir.resolve()
    output = args.output.resolve()
    step1 = require_files(
        step1_dir,
        [
            "SUCCESS.txt",
            "GSE143157_step1_matrices.h5",
            "metadata_donors_aligned.csv",
            "metadata_arrays_aligned.csv",
            "qc_summary.json",
        ],
        "Step 1",
    )
    step2 = require_files(
        step2_dir,
        ["SUCCESS.txt", "outlier_review.csv", "step2_summary.json"],
        "Step 2",
    )
    manifest, tables, context = prepare_metadata(step1, step2)
    print("Validated Step 1 and Step 2 evidence.")
    for name, details in context["cohort_counts"].items():
        print(f"{name}: {details['donors']} donors, {details['condition_counts']}")
    if args.dry_run:
        print("DRY RUN SUCCESS: no output files were written.")
        return
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Output directory is not empty; nothing will be overwritten: {output}")
    output.mkdir(parents=True, exist_ok=True)

    manifest.to_csv(output / "cohort_manifest.csv", index=False)
    for name, table in tables.items():
        table.to_csv(output / name, index=False)

    clean77 = tables["metadata_donors_primary_clean77.csv"]
    primary_h5 = output / PRIMARY_H5_NAME
    matrix_stats = copy_primary_h5(
        step1["GSE143157_step1_matrices.h5"],
        primary_h5,
        clean77,
        args.chunk_rows,
        args.gzip_level,
    )

    indexed_outliers = context["outliers"].set_index("donor_id")
    decision = {
        "pipeline": {
            "name": "GSE143157 Step 2B cohort freeze",
            "version": VERSION,
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": time.time() - started,
        },
        "inputs": {
            "step1_dir": str(step1_dir),
            "step2_dir": str(step2_dir),
            "step1_summary_sha256": hash_small_file(step1["qc_summary.json"]),
            "step2_summary_sha256": hash_small_file(step2["step2_summary.json"]),
            "outlier_review_sha256": hash_small_file(step2["outlier_review.csv"]),
        },
        "frozen_decision": {
            "primary_cohort": "primary_clean77",
            "primary_excluded_donors": ["MSA 30"],
            "primary_exclusion_basis": "Converging objective BS/5hmC technical QC evidence",
            "primary_same_donors_for_both_modalities": True,
            "control_11_primary_status": "retained",
            "control_11_reason": "Cross-modality PCA evidence without objective Detection-P/missingness technical failure",
            "automatic_pca_only_exclusion": False,
        },
        "MSA_30_step2_evidence": {
            "evidence_count": int(indexed_outliers.loc["MSA 30", "evidence_count"]),
            "evidence": str(indexed_outliers.loc["MSA 30", "outlier_evidence"]),
            "5mC_pca_robust_distance_z": float(indexed_outliers.loc["MSA 30", "5mC_pca_robust_distance_z"]),
            "5hmC_pca_robust_distance_z": float(indexed_outliers.loc["MSA 30", "5hmC_pca_robust_distance_z"]),
            "5hmC_negative_fraction_before_filter": float(indexed_outliers.loc["MSA 30", "5hmC_negative_difference_fraction_before_filter"]),
            "5hmC_missing_fraction_after_filter": float(indexed_outliers.loc["MSA 30", "5hmC_missing_fraction_after_probe_filter"]),
            "worst_detection_failure_fraction": float(indexed_outliers.loc["MSA 30", "worst_detection_failure_fraction"]),
        },
        "cohorts": context["cohort_counts"],
        "primary_hdf5": {
            "file": PRIMARY_H5_NAME,
            "construction": "Memory-safe row-chunked donor-column subset; CpG values are not recalculated",
            "source_donor_indices_dataset": "source_donor_indices",
            "matrix_stats": matrix_stats,
        },
        "downstream_policy": [
            "Use primary_clean77 for the main paired 5mC-only, 5hmC-only, and combined comparison.",
            "Repeat final findings on sensitivity_full78.",
            "Use sensitivity_no_MSA30_Control11_76 only as a stress-test.",
            "Do not reuse Step 2 top-variable CpGs as supervised ML features.",
            "Any ML feature selection must occur inside each training fold.",
        ],
    }
    with (output / "cohort_decision.json").open("w", encoding="utf-8") as handle:
        json.dump(decision, handle, indent=2, ensure_ascii=False, allow_nan=False, default=json_ready)
    write_html_report(output / "STEP2B_REPORT.html", decision)
    (output / "PRIMARY_COHORT.txt").write_text(
        "PRIMARY COHORT: Clean-77\n"
        "Included: 37 Control + 40 MSA = 77 donors\n"
        "Primary exclusion: MSA 30\n"
        "Reason: objective BS/5hmC technical QC failure\n"
        "Control 11: retained in primary; optional 76-donor stress-test supplied\n"
        "Full-78: required sensitivity analysis\n",
        encoding="utf-8",
    )
    (output / "SUCCESS.txt").write_text(
        "SUCCESS\n"
        f"Completed UTC: {decision['pipeline']['completed_utc']}\n"
        "Primary cohort: 77 donors (37 Control, 40 MSA)\n"
        "Primary excluded donor: MSA 30\n"
        "Full-78 sensitivity cohort preserved\n"
        "Control 11 retained in primary cohort\n",
        encoding="utf-8",
    )
    print("\nSTEP 2B SUCCESS")
    print(f"Primary HDF5: {primary_h5}")
    print(f"Open report: {output / 'STEP2B_REPORT.html'}")


def main() -> int:
    try:
        run(parse_args())
    except (ValueError, OSError, KeyError, pd.errors.ParserError) as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
