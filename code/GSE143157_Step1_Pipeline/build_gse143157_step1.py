#!/usr/bin/env python3
"""Build donor-aligned 5mC and 5hmC matrices from GSE143157.

The large interleaved input is read in row chunks.  It is never loaded in full.
The resulting matrices are stored as compressed, chunked HDF5 datasets so Step 2
can read them incrementally.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pandas as pd


VERSION = "1.1.0"

SAMPLE_ALIASES = (
    "matrix_array_column",
    "matrix_column",
    "sample_column",
    "beta_column",
    "array_column",
    "sample_number",
    "sample",
    "array",
)
DONOR_ALIASES = (
    "biological_subject_id",
    "biological_subject",
    "subject_id",
    "donor_id",
    "individual_id",
    "subject",
    "donor",
)
CONVERSION_ALIASES = (
    "conversion",
    "conversion_type",
    "assay",
    "bisulfite_conversion",
    "treatment",
)
CONDITION_ALIASES = (
    "condition",
    "diagnosis",
    "disease_status",
    "phenotype",
    "group",
    "disease",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Stream GSE143157_BMIQ_normalized_beta.txt and construct paired "
            "donor-level 5mC (oxBS) and 5hmC (BS - oxBS) matrices."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", required=True, type=Path, help="Large .txt or .txt.gz matrix")
    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="GSE143157_family.soft.gz directly, or an existing 156-row donor mapping CSV/TSV",
    )
    parser.add_argument("--output", required=True, type=Path, help="New or empty output directory")
    parser.add_argument("--chunk-rows", type=int, default=5000, help="CpG rows read per chunk")
    parser.add_argument("--missing-threshold", type=float, default=0.20, help="Remove 5hmC CpGs above this missing fraction")
    parser.add_argument("--detection-p-threshold", type=float, default=0.01, help="Detection-P failure cutoff for QC")
    parser.add_argument(
        "--mask-detection-failures",
        action="store_true",
        help="Also set beta values to missing when Detection-P is non-finite or above the cutoff",
    )
    parser.add_argument("--sample-field", help="Metadata field holding SAMPLE 1 ... SAMPLE 156")
    parser.add_argument("--donor-field", help="Metadata donor/biological-subject field")
    parser.add_argument("--conversion-field", help="Metadata BS/oxBS conversion field")
    parser.add_argument("--condition-field", help="Metadata MSA/Control field")
    parser.add_argument("--expected-cpgs", type=int, default=731661)
    parser.add_argument("--expected-arrays", type=int, default=156)
    parser.add_argument("--expected-donors", type=int, default=78)
    parser.add_argument("--expected-msa", type=int, default=41)
    parser.add_argument("--expected-control", type=int, default=37)
    parser.add_argument(
        "--allow-unexpected-counts",
        action="store_true",
        help="Permit non-GSE143157 counts (intended only for testing or deliberate subsets)",
    )
    parser.add_argument(
        "--allow-out-of-range-beta",
        action="store_true",
        help="Convert finite beta values outside [0,1] to missing instead of stopping",
    )
    parser.add_argument(
        "--skip-duplicate-cpg-check",
        action="store_true",
        help="Use less RAM by not keeping an exact set of probe IDs",
    )
    parser.add_argument("--progress-every", type=int, default=50000, help="Print progress after roughly this many rows")
    return parser.parse_args()


def fail(message: str) -> None:
    raise ValueError(message)


def slug(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


def clean_text(value: object) -> str:
    if value is None:
        return ""
    return re.sub(r"\s+", " ", str(value).strip())


def sample_key(value: object) -> str:
    text = clean_text(value)
    numeric = re.fullmatch(r"(?:sample\s*)?(\d+)(?:\.0+)?", text, flags=re.IGNORECASE)
    if numeric:
        return f"sample {int(numeric.group(1))}"
    return text.casefold()


def normalize_conversion(value: object) -> str:
    compact = re.sub(r"[^a-z0-9]+", "", clean_text(value).lower())
    if compact.startswith("oxbs") or "oxidativebisulfite" in compact:
        return "oxBS"
    if compact in {"bs", "bisulfite", "standardbisulfite", "bisulphite"}:
        return "BS"
    return clean_text(value)


def normalize_condition(value: object) -> str:
    text = clean_text(value)
    compact = re.sub(r"[^a-z0-9]+", "", text.lower())
    if "multiplesystematrophy" in compact or compact.startswith("msa"):
        return "MSA"
    if compact in {"control", "ctrl", "hc", "healthy", "healthycontrol", "normal"}:
        return "Control"
    return text


def open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("rt", encoding="utf-8-sig", newline="")


def read_matrix_header(path: Path) -> list[str]:
    with open_text(path) as handle:
        line = handle.readline()
    if not line:
        fail(f"Input matrix is empty: {path}")
    return next(csv.reader([line.rstrip("\r\n")], delimiter="\t"))


def read_family_soft_metadata(path: Path) -> pd.DataFrame:
    """Extract the 156 array records and donor pairing from GEO family SOFT."""
    records: list[dict[str, object]] = []
    current: dict[str, object] | None = None

    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    records.append(current)
                current = {
                    "gsm": clean_text(line.split("=", 1)[1]),
                    "characteristics": [],
                }
                continue
            if current is None or not line.startswith("!Sample_"):
                continue
            key, separator, value = line.partition(" = ")
            if not separator:
                continue
            key = key.removeprefix("!Sample_")
            value = clean_text(value)
            if key == "characteristics_ch1":
                characteristics = current["characteristics"]
                assert isinstance(characteristics, list)
                characteristics.append(value)
            elif key in {"title", "geo_accession", "source_name_ch1"}:
                current[key] = value
    if current is not None:
        records.append(current)

    rows: list[dict[str, str]] = []
    for index, record in enumerate(records, start=1):
        title = clean_text(record.get("title", ""))
        match = re.fullmatch(
            r"(SAMPLE\s+\d+)\s*:\s*(.+?)\s*,\s*(oxBS|BS)",
            title,
            flags=re.IGNORECASE,
        )
        if not match:
            fail(f"Could not parse GEO sample title #{index}: {title!r}")
        sample_column = re.sub(r"\s+", " ", match.group(1).upper())
        donor_id = clean_text(match.group(2))
        title_conversion = normalize_conversion(match.group(3))

        characteristics: dict[str, str] = {}
        raw_characteristics = record.get("characteristics", [])
        assert isinstance(raw_characteristics, list)
        for item in raw_characteristics:
            key, separator, value = str(item).partition(":")
            if separator:
                characteristics[slug(key)] = clean_text(value)

        characteristic_conversion = normalize_conversion(characteristics.get("conversion", ""))
        if characteristic_conversion and characteristic_conversion != title_conversion:
            fail(
                f"Conversion mismatch in {sample_column}: title={title_conversion}, "
                f"characteristics={characteristic_conversion}"
            )
        condition = normalize_condition(characteristics.get("condition", ""))
        rows.append(
            {
                "sample_column": sample_column,
                "sample_number": sample_column.split()[-1],
                "gsm": clean_text(record.get("geo_accession", record.get("gsm", ""))),
                "biological_subject_id": donor_id,
                "condition": condition,
                "conversion": title_conversion,
                "sex": characteristics.get("sex", ""),
                "age": characteristics.get("age", ""),
                "pmi_hours": characteristics.get("pmi_hours", ""),
                "brain_bank": characteristics.get("brain_bank", ""),
                "identifier": characteristics.get("identifier", ""),
                "sample_type": characteristics.get("sample_type", ""),
                "tissue": characteristics.get("tissue", ""),
                "source": clean_text(record.get("source_name_ch1", "")),
                "geo_title": title,
            }
        )

    if not rows:
        fail(f"No ^SAMPLE records were found in GEO family SOFT: {path}")
    return pd.DataFrame(rows)


def read_metadata(path: Path) -> pd.DataFrame:
    lower_name = path.name.lower()
    if lower_name.endswith(".soft") or lower_name.endswith(".soft.gz"):
        return read_family_soft_metadata(path)
    if path.suffix.lower() in {".xlsx", ".xls"}:
        fail("Metadata must be saved as CSV or TSV, not Excel. Use Save As -> CSV UTF-8 in Excel.")
    try:
        frame = pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)
    except Exception as exc:
        fail(f"Could not read metadata as CSV/TSV: {exc}")
    frame.columns = [clean_text(c).lstrip("\ufeff") for c in frame.columns]
    if frame.empty:
        fail("Metadata contains no data rows.")
    if len(set(frame.columns)) != len(frame.columns):
        fail("Metadata contains duplicate column names.")
    return frame


def resolve_named_field(frame: pd.DataFrame, override: str | None, aliases: Iterable[str], label: str) -> str | None:
    by_slug = {slug(c): c for c in frame.columns}
    if override:
        if override in frame.columns:
            return override
        candidate = by_slug.get(slug(override))
        if candidate:
            return candidate
        fail(f"--{label}-field '{override}' is not a metadata column. Available: {list(frame.columns)}")
    for alias in aliases:
        if alias in by_slug:
            return by_slug[alias]
    return None


def infer_sample_field(frame: pd.DataFrame, beta_headers: list[str], override: str | None) -> str:
    if override:
        field = resolve_named_field(frame, override, SAMPLE_ALIASES, "sample")
        assert field is not None
        return field

    beta_keys = set(map(sample_key, beta_headers))
    scores: list[tuple[int, int, str]] = []
    aliases = {name: index for index, name in enumerate(SAMPLE_ALIASES)}
    for column in frame.columns:
        keys = [sample_key(v) for v in frame[column]]
        score = sum(k in beta_keys for k in keys)
        priority = aliases.get(slug(column), len(aliases) + 1)
        scores.append((score, -priority, column))
    scores.sort(reverse=True)
    best_score, _, best_field = scores[0]
    if best_score != len(beta_headers):
        top = [(column, score) for score, _, column in scores[:5]]
        fail(
            "Could not align all matrix SAMPLE columns to metadata. "
            f"Best candidate '{best_field}' matched {best_score}/{len(beta_headers)}. "
            f"Top scores: {top}. Supply --sample-field explicitly if needed."
        )
    return best_field


def infer_donor_field(frame: pd.DataFrame, override: str | None) -> str:
    named = resolve_named_field(frame, override, DONOR_ALIASES, "donor")
    if named:
        return named
    candidates = []
    for column in frame.columns:
        values = frame[column].map(clean_text)
        counts = values.value_counts(dropna=False)
        if len(counts) * 2 == len(frame) and (counts == 2).all() and "" not in counts.index:
            candidates.append(column)
    if len(candidates) == 1:
        return candidates[0]
    fail(
        "Could not uniquely identify the donor field. "
        f"Pair-shaped candidates: {candidates}. Supply --donor-field explicitly."
    )


def infer_conversion_field(frame: pd.DataFrame, override: str | None) -> str:
    named = resolve_named_field(frame, override, CONVERSION_ALIASES, "conversion")
    if named:
        return named
    candidates = []
    for column in frame.columns:
        values = frame[column].map(normalize_conversion)
        if set(values) == {"BS", "oxBS"}:
            candidates.append(column)
    if len(candidates) == 1:
        return candidates[0]
    fail(f"Could not uniquely identify BS/oxBS field. Candidates: {candidates}. Use --conversion-field.")


def infer_condition_field(frame: pd.DataFrame, override: str | None) -> str:
    named = resolve_named_field(frame, override, CONDITION_ALIASES, "condition")
    if named:
        return named
    candidates = []
    for column in frame.columns:
        values = set(frame[column].map(normalize_condition))
        if values == {"MSA", "Control"}:
            candidates.append(column)
    if len(candidates) == 1:
        return candidates[0]
    fail(f"Could not uniquely identify MSA/Control field. Candidates: {candidates}. Use --condition-field.")


def check_expected(observed: int, expected: int, label: str, allow: bool) -> None:
    if observed != expected and not allow:
        fail(f"Expected {expected} {label}, observed {observed}. Check the files or use --allow-unexpected-counts deliberately.")


def align_metadata(
    frame: pd.DataFrame,
    beta_headers: list[str],
    detection_headers: list[str],
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray, dict[str, str]]:
    sample_field = infer_sample_field(frame, beta_headers, args.sample_field)
    donor_field = infer_donor_field(frame, args.donor_field)
    conversion_field = infer_conversion_field(frame, args.conversion_field)
    condition_field = infer_condition_field(frame, args.condition_field)

    fields = {
        "sample": sample_field,
        "donor": donor_field,
        "conversion": conversion_field,
        "condition": condition_field,
    }

    work = frame.copy()
    work["_array_key"] = work[sample_field].map(sample_key)
    work["_donor_id"] = work[donor_field].map(clean_text)
    work["_conversion"] = work[conversion_field].map(normalize_conversion)
    work["_condition"] = work[condition_field].map(normalize_condition)

    if (work["_donor_id"] == "").any():
        fail("Metadata has an empty donor ID.")
    if work["_array_key"].duplicated().any():
        dupes = work.loc[work["_array_key"].duplicated(keep=False), sample_field].tolist()[:10]
        fail(f"Metadata maps an array more than once. Examples: {dupes}")

    beta_keys = [sample_key(v) for v in beta_headers]
    missing = [beta_headers[i] for i, key in enumerate(beta_keys) if key not in set(work["_array_key"])]
    extra = work.loc[~work["_array_key"].isin(beta_keys), sample_field].tolist()
    if missing or extra:
        fail(f"Metadata/matrix mismatch. Missing metadata for: {missing[:10]}; extra metadata rows: {extra[:10]}")

    ordered = work.set_index("_array_key", drop=False).loc[beta_keys].reset_index(drop=True)
    # Assignment (rather than DataFrame.insert) deliberately supports mapping
    # files that already use one of these canonical output column names.
    ordered["matrix_beta_index"] = np.arange(len(ordered), dtype=int)
    ordered["matrix_array_column"] = beta_headers
    ordered["matrix_detection_p_column"] = detection_headers

    conversions = Counter(ordered["_conversion"])
    if set(conversions) != {"BS", "oxBS"}:
        fail(f"Conversion values must normalize to only BS and oxBS; observed {dict(conversions)}")

    pairs = []
    for donor_id, group in ordered.groupby("_donor_id", sort=False):
        counts = Counter(group["_conversion"])
        if counts != {"BS": 1, "oxBS": 1}:
            fail(f"Donor '{donor_id}' must have exactly one BS and one oxBS row; observed {dict(counts)}")
        conditions = set(group["_condition"])
        if len(conditions) != 1:
            fail(f"Donor '{donor_id}' has inconsistent condition labels: {conditions}")
        bs_row = group.loc[group["_conversion"] == "BS"].iloc[0]
        ox_row = group.loc[group["_conversion"] == "oxBS"].iloc[0]
        pairs.append((int(ox_row["matrix_beta_index"]), donor_id, bs_row, ox_row))
    pairs.sort(key=lambda item: item[0])

    donor_rows = []
    original_columns = list(frame.columns)
    excluded = {sample_field, donor_field, conversion_field, condition_field}
    for _, donor_id, bs_row, ox_row in pairs:
        row: dict[str, object] = {
            "donor_index": len(donor_rows),
            "donor_id": donor_id,
            "condition": ox_row["_condition"],
            "oxBS_array_column": ox_row["matrix_array_column"],
            "BS_array_column": bs_row["matrix_array_column"],
            "oxBS_beta_index": int(ox_row["matrix_beta_index"]),
            "BS_beta_index": int(bs_row["matrix_beta_index"]),
        }
        for column in original_columns:
            if column in excluded:
                continue
            ox_value = clean_text(ox_row[column])
            bs_value = clean_text(bs_row[column])
            if ox_value == bs_value:
                row[column] = ox_value
            else:
                row[f"oxBS__{column}"] = ox_value
                row[f"BS__{column}"] = bs_value
        donor_rows.append(row)

    donor_frame = pd.DataFrame(donor_rows)
    bs_indices = donor_frame["BS_beta_index"].to_numpy(dtype=int)
    ox_indices = donor_frame["oxBS_beta_index"].to_numpy(dtype=int)
    return ordered, donor_frame, bs_indices, ox_indices, fields


def append_rows(dataset: h5py.Dataset, values: np.ndarray) -> None:
    if values.shape[0] == 0:
        return
    start = dataset.shape[0]
    dataset.resize(start + values.shape[0], axis=0)
    dataset[start : start + values.shape[0]] = values


def append_strings(dataset: h5py.Dataset, values: list[str]) -> None:
    if not values:
        return
    start = dataset.shape[0]
    dataset.resize(start + len(values), axis=0)
    dataset[start : start + len(values)] = np.asarray(values, dtype=object)


def finite_min(values: np.ndarray, axis: int = 0) -> np.ndarray:
    masked = np.where(np.isfinite(values), values, np.nan)
    with np.errstate(all="ignore"):
        result = np.nanmin(masked, axis=axis)
    return result


def finite_max(values: np.ndarray, axis: int = 0) -> np.ndarray:
    masked = np.where(np.isfinite(values), values, np.nan)
    with np.errstate(all="ignore"):
        result = np.nanmax(masked, axis=axis)
    return result


def json_number(value: float | np.floating) -> float | None:
    numeric = float(value)
    return numeric if math.isfinite(numeric) else None


def run(args: argparse.Namespace) -> None:
    started = time.time()
    input_path = args.input.resolve()
    metadata_path = args.metadata.resolve()
    output_dir = args.output.resolve()

    if not input_path.is_file():
        fail(f"Input matrix not found or not a file: {input_path}")
    if not metadata_path.is_file():
        fail(f"Metadata not found or not a file: {metadata_path}")
    if args.chunk_rows < 1:
        fail("--chunk-rows must be positive.")
    if not 0 <= args.missing_threshold < 1:
        fail("--missing-threshold must be in [0, 1).")
    if not 0 < args.detection_p_threshold <= 1:
        fail("--detection-p-threshold must be in (0, 1].")
    if output_dir.exists() and any(output_dir.iterdir()):
        fail(f"Output directory is not empty: {output_dir}. Choose a new folder; existing results are never overwritten.")
    output_dir.mkdir(parents=True, exist_ok=True)

    header = read_matrix_header(input_path)
    if len(header) < 3 or (len(header) - 1) % 2 != 0:
        fail(f"Expected ID_REF plus interleaved beta/Detection-P pairs; observed {len(header)} total columns.")
    beta_headers = header[1::2]
    detection_headers = header[2::2]
    n_arrays = len(beta_headers)
    check_expected(n_arrays, args.expected_arrays, "beta arrays", args.allow_unexpected_counts)
    if len(detection_headers) != n_arrays:
        fail("Beta and Detection-P column counts differ.")
    detection_name_hits = sum("detect" in name.lower() or "pval" in name.lower() for name in detection_headers)
    if detection_name_hits < max(1, int(0.90 * n_arrays)):
        fail(
            "The positional columns do not look interleaved: fewer than 90% of the expected "
            "Detection-P headers contain 'detect' or 'pval'."
        )

    metadata = read_metadata(metadata_path)
    check_expected(len(metadata), n_arrays, "metadata rows matching arrays", args.allow_unexpected_counts)
    array_meta, donor_meta, bs_indices, ox_indices, resolved_fields = align_metadata(
        metadata, beta_headers, detection_headers, args
    )
    n_donors = len(donor_meta)
    check_expected(n_donors, args.expected_donors, "paired donors", args.allow_unexpected_counts)
    condition_counts = Counter(donor_meta["condition"])
    check_expected(condition_counts.get("MSA", 0), args.expected_msa, "MSA donors", args.allow_unexpected_counts)
    check_expected(condition_counts.get("Control", 0), args.expected_control, "Control donors", args.allow_unexpected_counts)

    print(f"Validated layout: {n_arrays} beta + {n_arrays} Detection-P columns")
    print(f"Validated metadata: {n_donors} paired donors; conditions={dict(condition_counts)}")
    print(f"Chunk size: {args.chunk_rows:,} CpGs; detection masking: {args.mask_detection_failures}")

    numeric_names = [f"value_{i:03d}" for i in range(1, 2 * n_arrays + 1)]
    parser_names = ["probe_id", *numeric_names]
    dtypes = {"probe_id": str, **{name: np.float32 for name in numeric_names}}

    array_beta_valid = np.zeros(n_arrays, dtype=np.int64)
    array_beta_nonfinite = np.zeros(n_arrays, dtype=np.int64)
    array_beta_out_of_range = np.zeros(n_arrays, dtype=np.int64)
    array_beta_sum = np.zeros(n_arrays, dtype=np.float64)
    array_beta_min = np.full(n_arrays, np.inf, dtype=np.float64)
    array_beta_max = np.full(n_arrays, -np.inf, dtype=np.float64)
    array_det_valid = np.zeros(n_arrays, dtype=np.int64)
    array_det_nonfinite = np.zeros(n_arrays, dtype=np.int64)
    array_det_out_of_range = np.zeros(n_arrays, dtype=np.int64)
    array_det_gt_threshold = np.zeros(n_arrays, dtype=np.int64)
    array_det_gt_005 = np.zeros(n_arrays, dtype=np.int64)

    donor_5mc_missing = np.zeros(n_donors, dtype=np.int64)
    donor_5mc_sum = np.zeros(n_donors, dtype=np.float64)
    donor_5mc_valid = np.zeros(n_donors, dtype=np.int64)
    donor_hmc_negative_all = np.zeros(n_donors, dtype=np.int64)
    donor_hmc_missing_retained = np.zeros(n_donors, dtype=np.int64)
    donor_hmc_valid_retained = np.zeros(n_donors, dtype=np.int64)
    donor_hmc_sum_retained = np.zeros(n_donors, dtype=np.float64)
    donor_detection_masked_5mc = np.zeros(n_donors, dtype=np.int64)
    donor_detection_masked_hmc_pair = np.zeros(n_donors, dtype=np.int64)

    total_rows = 0
    retained_hmc_rows = 0
    removed_hmc_rows = 0
    total_negative_hmc_cells = 0
    total_hmc_input_missing_cells = 0
    total_detection_masked_beta_cells = 0
    duplicate_count = 0
    duplicate_examples: list[str] = []
    seen_probe_ids: set[str] | None = None if args.skip_duplicate_cpg_check else set()

    h5_tmp = output_dir / "GSE143157_step1_matrices.incomplete.h5"
    h5_final = output_dir / "GSE143157_step1_matrices.h5"
    probe_qc_tmp = output_dir / "5hmC_probe_qc.incomplete.tsv.gz"
    probe_qc_final = output_dir / "5hmC_probe_qc.tsv.gz"

    string_dtype = h5py.string_dtype(encoding="utf-8")
    h5_row_chunk = min(max(256, args.chunk_rows), 4096)
    next_progress = args.progress_every

    with h5py.File(h5_tmp, "w") as h5, gzip.open(
        probe_qc_tmp, "wt", encoding="utf-8", newline=""
    ) as probe_qc:
        h5.attrs["pipeline_version"] = VERSION
        h5.attrs["source_accession"] = "GSE143157"
        h5.attrs["construction_5mC"] = "oxBS beta"
        h5.attrs["construction_5hmC"] = "BS beta - oxBS beta"
        h5.attrs["negative_5hmC_handling"] = "set to NaN"
        h5.attrs["missing_fraction_filter"] = float(args.missing_threshold)
        h5.attrs["detection_p_threshold"] = float(args.detection_p_threshold)
        h5.attrs["detection_p_masking_enabled"] = bool(args.mask_detection_failures)
        h5.attrs["matrix_orientation"] = "rows=CpGs, columns=donors"
        h5.create_dataset("donor_ids", data=np.asarray(donor_meta["donor_id"], dtype=object), dtype=string_dtype)
        h5.create_dataset("conditions", data=np.asarray(donor_meta["condition"], dtype=object), dtype=string_dtype)

        group_5mc = h5.create_group("5mC")
        values_5mc = group_5mc.create_dataset(
            "values",
            shape=(0, n_donors),
            maxshape=(None, n_donors),
            dtype=np.float32,
            chunks=(h5_row_chunk, n_donors),
            compression="gzip",
            compression_opts=4,
            shuffle=True,
            fillvalue=np.nan,
        )
        ids_5mc = group_5mc.create_dataset(
            "cpg_ids", shape=(0,), maxshape=(None,), dtype=string_dtype, chunks=(h5_row_chunk,), compression="gzip"
        )

        group_hmc = h5.create_group("5hmC")
        values_hmc = group_hmc.create_dataset(
            "values",
            shape=(0, n_donors),
            maxshape=(None, n_donors),
            dtype=np.float32,
            chunks=(h5_row_chunk, n_donors),
            compression="gzip",
            compression_opts=4,
            shuffle=True,
            fillvalue=np.nan,
        )
        ids_hmc = group_hmc.create_dataset(
            "cpg_ids", shape=(0,), maxshape=(None,), dtype=string_dtype, chunks=(h5_row_chunk,), compression="gzip"
        )

        probe_qc.write(
            "cpg_id\tinput_missing_pairs\tnegative_differences\tmissing_after_negative\tmissing_fraction\tretained\n"
        )

        reader = pd.read_csv(
            input_path,
            sep="\t",
            header=None,
            names=parser_names,
            skiprows=1,
            dtype=dtypes,
            chunksize=args.chunk_rows,
            compression="infer",
            encoding="utf-8-sig",
            keep_default_na=True,
            na_values=["NA", "NaN", "nan", "N/A", "null", "NULL", ""],
            on_bad_lines="error",
            engine="c",
            memory_map=not input_path.name.lower().endswith(".gz"),
        )

        for chunk_number, chunk in enumerate(reader, start=1):
            probe_ids = [clean_text(v) for v in chunk.iloc[:, 0].tolist()]
            if any(not probe for probe in probe_ids):
                fail(f"Found an empty probe ID in chunk {chunk_number}.")
            if seen_probe_ids is not None:
                for probe in probe_ids:
                    if probe in seen_probe_ids:
                        duplicate_count += 1
                        if len(duplicate_examples) < 20:
                            duplicate_examples.append(probe)
                    else:
                        seen_probe_ids.add(probe)

            beta = chunk.iloc[:, 1::2].to_numpy(dtype=np.float32, copy=False)
            detection = chunk.iloc[:, 2::2].to_numpy(dtype=np.float32, copy=False)
            if beta.shape[1] != n_arrays or detection.shape[1] != n_arrays:
                fail(f"Internal positional parse error in chunk {chunk_number}: beta={beta.shape}, detection={detection.shape}")

            beta_finite = np.isfinite(beta)
            beta_out = beta_finite & ((beta < 0.0) | (beta > 1.0))
            if beta_out.any() and not args.allow_out_of_range_beta:
                row_i, col_i = np.argwhere(beta_out)[0]
                fail(
                    f"Finite beta outside [0,1] at probe {probe_ids[int(row_i)]}, "
                    f"array {beta_headers[int(col_i)]}: {float(beta[row_i, col_i])}. "
                    "Inspect the file; use --allow-out-of-range-beta only if deliberate."
                )

            array_beta_valid += beta_finite.sum(axis=0)
            array_beta_nonfinite += (~beta_finite).sum(axis=0)
            array_beta_out_of_range += beta_out.sum(axis=0)
            array_beta_sum += np.where(beta_finite, beta, 0.0).sum(axis=0, dtype=np.float64)
            chunk_min = finite_min(beta)
            chunk_max = finite_max(beta)
            array_beta_min = np.minimum(array_beta_min, np.where(np.isfinite(chunk_min), chunk_min, np.inf))
            array_beta_max = np.maximum(array_beta_max, np.where(np.isfinite(chunk_max), chunk_max, -np.inf))

            det_finite = np.isfinite(detection)
            det_out = det_finite & ((detection < 0.0) | (detection > 1.0))
            array_det_valid += det_finite.sum(axis=0)
            array_det_nonfinite += (~det_finite).sum(axis=0)
            array_det_out_of_range += det_out.sum(axis=0)
            array_det_gt_threshold += (det_finite & (detection > args.detection_p_threshold)).sum(axis=0)
            array_det_gt_005 += (det_finite & (detection > 0.05)).sum(axis=0)

            ox = beta[:, ox_indices].astype(np.float32, copy=True)
            bs = beta[:, bs_indices].astype(np.float32, copy=True)
            ox[~np.isfinite(ox)] = np.nan
            bs[~np.isfinite(bs)] = np.nan
            if args.allow_out_of_range_beta:
                ox[(ox < 0.0) | (ox > 1.0)] = np.nan
                bs[(bs < 0.0) | (bs > 1.0)] = np.nan

            if args.mask_detection_failures:
                det_ox = detection[:, ox_indices]
                det_bs = detection[:, bs_indices]
                bad_ox = (~np.isfinite(det_ox)) | (det_ox > args.detection_p_threshold)
                bad_bs = (~np.isfinite(det_bs)) | (det_bs > args.detection_p_threshold)
                donor_detection_masked_5mc += bad_ox.sum(axis=0)
                donor_detection_masked_hmc_pair += (bad_ox | bad_bs).sum(axis=0)
                total_detection_masked_beta_cells += int(bad_ox.sum() + bad_bs.sum())
                ox[bad_ox] = np.nan
                bs[bad_bs] = np.nan

            five_mc = ox
            input_missing = ~np.isfinite(bs) | ~np.isfinite(ox)
            difference = bs - ox
            negative = np.isfinite(difference) & (difference < 0.0)
            five_hmc = difference.astype(np.float32, copy=True)
            five_hmc[negative | input_missing] = np.nan

            input_missing_count = input_missing.sum(axis=1).astype(np.int16)
            negative_count = negative.sum(axis=1).astype(np.int16)
            missing_count = (~np.isfinite(five_hmc)).sum(axis=1).astype(np.int16)
            missing_fraction = missing_count.astype(np.float64) / n_donors
            retained = missing_fraction <= args.missing_threshold

            append_rows(values_5mc, five_mc)
            append_strings(ids_5mc, probe_ids)
            retained_values = five_hmc[retained]
            retained_ids = [probe_ids[i] for i in np.flatnonzero(retained)]
            append_rows(values_hmc, retained_values)
            append_strings(ids_hmc, retained_ids)

            qc_lines = [
                f"{probe_ids[i]}\t{int(input_missing_count[i])}\t{int(negative_count[i])}\t"
                f"{int(missing_count[i])}\t{missing_fraction[i]:.8f}\t{int(retained[i])}"
                for i in range(len(probe_ids))
            ]
            probe_qc.write("\n".join(qc_lines) + "\n")

            donor_5mc_missing += (~np.isfinite(five_mc)).sum(axis=0)
            donor_5mc_valid += np.isfinite(five_mc).sum(axis=0)
            donor_5mc_sum += np.nansum(five_mc, axis=0, dtype=np.float64)
            donor_hmc_negative_all += negative.sum(axis=0)
            donor_hmc_missing_retained += (~np.isfinite(retained_values)).sum(axis=0)
            donor_hmc_valid_retained += np.isfinite(retained_values).sum(axis=0)
            donor_hmc_sum_retained += np.nansum(retained_values, axis=0, dtype=np.float64)

            rows_here = len(probe_ids)
            total_rows += rows_here
            kept_here = int(retained.sum())
            retained_hmc_rows += kept_here
            removed_hmc_rows += rows_here - kept_here
            total_negative_hmc_cells += int(negative.sum())
            total_hmc_input_missing_cells += int(input_missing.sum())

            if total_rows >= next_progress:
                elapsed = time.time() - started
                print(
                    f"Processed {total_rows:,} CpGs | retained 5hmC {retained_hmc_rows:,} | "
                    f"elapsed {elapsed / 60:.1f} min",
                    flush=True,
                )
                while next_progress <= total_rows:
                    next_progress += args.progress_every
            if chunk_number % 20 == 0:
                h5.flush()

        h5.attrs["total_5mC_cpgs"] = int(total_rows)
        h5.attrs["retained_5hmC_cpgs"] = int(retained_hmc_rows)
        h5.attrs["donor_count"] = int(n_donors)
        h5.flush()

    check_expected(total_rows, args.expected_cpgs, "CpG rows", args.allow_unexpected_counts)
    if duplicate_count:
        fail(f"Found {duplicate_count} duplicate probe IDs. Examples: {duplicate_examples}")

    h5_tmp.replace(h5_final)
    probe_qc_tmp.replace(probe_qc_final)

    array_export = array_meta.copy()
    array_export["donor_id"] = array_export["_donor_id"]
    array_export["condition_normalized"] = array_export["_condition"]
    array_export["conversion_normalized"] = array_export["_conversion"]
    array_export = array_export.drop(columns=["_array_key", "_donor_id", "_condition", "_conversion"])
    canonical_front = [
        "matrix_beta_index",
        "matrix_array_column",
        "matrix_detection_p_column",
        "donor_id",
        "condition_normalized",
        "conversion_normalized",
    ]
    array_export = array_export[canonical_front + [c for c in array_export.columns if c not in canonical_front]]
    array_export.to_csv(output_dir / "metadata_arrays_aligned.csv", index=False)
    donor_meta.to_csv(output_dir / "metadata_donors_aligned.csv", index=False)

    per_array_rows = []
    for i, beta_name in enumerate(beta_headers):
        meta_row = array_meta.iloc[i]
        beta_mean = array_beta_sum[i] / array_beta_valid[i] if array_beta_valid[i] else math.nan
        per_array_rows.append(
            {
                "matrix_array_column": beta_name,
                "matrix_detection_p_column": detection_headers[i],
                "donor_id": meta_row["_donor_id"],
                "condition": meta_row["_condition"],
                "conversion": meta_row["_conversion"],
                "total_cpgs": total_rows,
                "beta_valid_count": int(array_beta_valid[i]),
                "beta_nonfinite_count": int(array_beta_nonfinite[i]),
                "beta_out_of_range_count": int(array_beta_out_of_range[i]),
                "beta_mean": json_number(beta_mean),
                "beta_min": json_number(array_beta_min[i]),
                "beta_max": json_number(array_beta_max[i]),
                "detection_p_valid_count": int(array_det_valid[i]),
                "detection_p_nonfinite_count": int(array_det_nonfinite[i]),
                "detection_p_out_of_range_count": int(array_det_out_of_range[i]),
                "detection_p_gt_0_01_count": int(array_det_gt_threshold[i]),
                "detection_p_gt_0_01_fraction_of_all_cpgs": float(array_det_gt_threshold[i] / total_rows),
                "detection_p_gt_0_05_count": int(array_det_gt_005[i]),
                "detection_p_gt_0_05_fraction_of_all_cpgs": float(array_det_gt_005[i] / total_rows),
            }
        )
    pd.DataFrame(per_array_rows).to_csv(output_dir / "qc_per_array.csv", index=False)

    per_donor = donor_meta[["donor_index", "donor_id", "condition", "oxBS_array_column", "BS_array_column"]].copy()
    per_donor["5mC_missing_count"] = donor_5mc_missing
    per_donor["5mC_missing_fraction"] = donor_5mc_missing / total_rows
    per_donor["5mC_mean"] = np.divide(
        donor_5mc_sum,
        donor_5mc_valid,
        out=np.full(n_donors, np.nan, dtype=float),
        where=donor_5mc_valid > 0,
    )
    per_donor["5hmC_negative_difference_count_before_filter"] = donor_hmc_negative_all
    per_donor["5hmC_negative_difference_fraction_before_filter"] = donor_hmc_negative_all / total_rows
    per_donor["5hmC_missing_count_after_probe_filter"] = donor_hmc_missing_retained
    per_donor["5hmC_missing_fraction_after_probe_filter"] = donor_hmc_missing_retained / max(retained_hmc_rows, 1)
    per_donor["5hmC_mean_after_probe_filter"] = np.divide(
        donor_hmc_sum_retained,
        donor_hmc_valid_retained,
        out=np.full(n_donors, np.nan, dtype=float),
        where=donor_hmc_valid_retained > 0,
    )
    per_donor["detection_masked_5mC_count"] = donor_detection_masked_5mc
    per_donor["detection_masked_5hmC_pair_count"] = donor_detection_masked_hmc_pair
    per_donor.to_csv(output_dir / "qc_per_donor.csv", index=False)

    elapsed = time.time() - started
    beta_global_min = min((v for v in array_beta_min if math.isfinite(v)), default=math.nan)
    beta_global_max = max((v for v in array_beta_max if math.isfinite(v)), default=math.nan)
    summary = {
        "pipeline": {
            "name": "GSE143157 Step 1 streaming 5mC/5hmC builder",
            "version": VERSION,
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": elapsed,
        },
        "inputs": {
            "matrix_path": str(input_path),
            "matrix_size_bytes": input_path.stat().st_size,
            "metadata_path": str(metadata_path),
            "resolved_metadata_fields": resolved_fields,
        },
        "parameters": {
            "chunk_rows": args.chunk_rows,
            "missing_threshold": args.missing_threshold,
            "detection_p_threshold": args.detection_p_threshold,
            "detection_p_masking_enabled": args.mask_detection_failures,
            "allow_out_of_range_beta": args.allow_out_of_range_beta,
            "duplicate_cpg_check_performed": not args.skip_duplicate_cpg_check,
        },
        "observed_structure": {
            "total_input_columns": len(header),
            "cpg_rows": total_rows,
            "beta_arrays": n_arrays,
            "detection_p_arrays": len(detection_headers),
            "paired_donors": n_donors,
            "condition_counts": dict(condition_counts),
            "duplicate_cpg_ids": duplicate_count if not args.skip_duplicate_cpg_check else None,
            "beta_global_min": json_number(beta_global_min),
            "beta_global_max": json_number(beta_global_max),
            "beta_nonfinite_cells": int(array_beta_nonfinite.sum()),
            "beta_out_of_range_cells": int(array_beta_out_of_range.sum()),
            "detection_p_nonfinite_cells": int(array_det_nonfinite.sum()),
            "detection_p_out_of_range_cells": int(array_det_out_of_range.sum()),
            "detection_p_failures_gt_threshold": int(array_det_gt_threshold.sum()),
        },
        "construction": {
            "5mC": "paired donor oxBS beta value",
            "5hmC": "paired donor BS beta minus oxBS beta",
            "negative_5hmC": "set to missing (NaN)",
            "5hmC_probe_filter": "retained when missing_fraction <= threshold; removed only when > threshold",
            "imputation": "not performed in Step 1",
            "detection_p_use": (
                "masked beta when Detection-P was non-finite or above threshold"
                if args.mask_detection_failures
                else "reported for QC only; did not alter beta values"
            ),
        },
        "outputs_summary": {
            "5mC_cpgs": total_rows,
            "5mC_donors": n_donors,
            "5hmC_cpgs_retained": retained_hmc_rows,
            "5hmC_cpgs_removed_gt_missing_threshold": removed_hmc_rows,
            "5hmC_donors": n_donors,
            "5hmC_negative_cells_set_missing": total_negative_hmc_cells,
            "5hmC_input_missing_pair_cells": total_hmc_input_missing_cells,
            "detection_masked_beta_cells": total_detection_masked_beta_cells,
        },
        "files": {
            "matrices_hdf5": h5_final.name,
            "array_metadata": "metadata_arrays_aligned.csv",
            "donor_metadata": "metadata_donors_aligned.csv",
            "per_array_qc": "qc_per_array.csv",
            "per_donor_qc": "qc_per_donor.csv",
            "per_probe_5hmC_qc": probe_qc_final.name,
        },
    }
    with (output_dir / "qc_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False, allow_nan=False)

    success_text = (
        "SUCCESS\n"
        f"Completed UTC: {summary['pipeline']['completed_utc']}\n"
        f"5mC shape: {total_rows} CpGs x {n_donors} donors\n"
        f"5hmC shape after >{args.missing_threshold:.0%} missing filter: {retained_hmc_rows} CpGs x {n_donors} donors\n"
        "Read qc_summary.json before Step 2.\n"
    )
    (output_dir / "SUCCESS.txt").write_text(success_text, encoding="utf-8")

    print("\nSUCCESS")
    print(f"5mC:  {total_rows:,} CpGs x {n_donors} donors")
    print(f"5hmC: {retained_hmc_rows:,} CpGs x {n_donors} donors after missingness filter")
    print(f"QC summary: {output_dir / 'qc_summary.json'}")
    print(f"Matrices:   {h5_final}")


def main() -> int:
    args = parse_args()
    try:
        run(args)
    except (ValueError, OSError, pd.errors.ParserError) as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
