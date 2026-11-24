#!/usr/bin/env python3
"""Integrity checks for GSE143157 Step 2B outputs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def decode(values) -> list[str]:
    return [item.decode("utf-8") if isinstance(item, bytes) else str(item) for item in values]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def check(output: Path) -> None:
    required = [
        "SUCCESS.txt",
        "PRIMARY_COHORT.txt",
        "STEP2B_REPORT.html",
        "cohort_decision.json",
        "cohort_manifest.csv",
        "metadata_donors_sensitivity_full78.csv",
        "metadata_donors_primary_clean77.csv",
        "metadata_donors_sensitivity_no_MSA30_Control11_76.csv",
        "metadata_arrays_sensitivity_full78.csv",
        "metadata_arrays_primary_clean77.csv",
        "metadata_arrays_sensitivity_no_MSA30_Control11_76.csv",
        "GSE143157_primary_clean77.h5",
    ]
    missing = [name for name in required if not (output / name).is_file()]
    require(not missing, f"Missing Step 2B files: {missing}")

    decision = json.loads((output / "cohort_decision.json").read_text(encoding="utf-8"))
    require(decision["cohorts"]["primary_clean77"]["donors"] == 77, "Primary count is not 77.")
    require(decision["cohorts"]["sensitivity_full78"]["donors"] == 78, "Full count is not 78.")
    require(
        decision["cohorts"]["sensitivity_no_MSA30_Control11_76"]["donors"] == 76,
        "No-two sensitivity count is not 76.",
    )
    primary_meta = pd.read_csv(output / "metadata_donors_primary_clean77.csv")
    manifest = pd.read_csv(output / "cohort_manifest.csv")
    require(len(primary_meta) == 77, "Primary metadata does not have 77 rows.")
    require("MSA 30" not in set(primary_meta["donor_id"]), "MSA 30 is still in primary metadata.")
    require("Control 11" in set(primary_meta["donor_id"]), "Control 11 is missing from primary metadata.")
    require(primary_meta["condition"].value_counts().to_dict() == {"MSA": 40, "Control": 37}, "Primary labels are wrong.")
    require(len(manifest) == 78, "Manifest does not have 78 rows.")

    source_path = Path(decision["inputs"]["step1_dir"]) / "GSE143157_step1_matrices.h5"
    destination_path = output / "GSE143157_primary_clean77.h5"
    keep = primary_meta["original_hdf5_column"].astype(int).to_numpy()
    expected_ids = primary_meta["donor_id"].astype(str).tolist()
    with h5py.File(source_path, "r") as source, h5py.File(destination_path, "r") as dest:
        require(decode(dest["donor_ids"][:]) == expected_ids, "Primary HDF5 donor IDs/order are wrong.")
        require(np.array_equal(dest["source_donor_indices"][:], keep), "Source donor indices are wrong.")
        for modality in ["5mC", "5hmC"]:
            src_values = source[f"{modality}/values"]
            dst_values = dest[f"{modality}/values"]
            require(dst_values.shape == (src_values.shape[0], 77), f"{modality} output shape is wrong.")
            rows = sorted({0, src_values.shape[0] // 2, src_values.shape[0] - 1})
            for row in rows:
                expected = np.asarray(src_values[row : row + 1, keep])
                observed = np.asarray(dst_values[row : row + 1, :])
                require(np.array_equal(expected, observed, equal_nan=True), f"{modality} value mismatch at row {row}.")
            src_ids = source[f"{modality}/cpg_ids"]
            dst_ids = dest[f"{modality}/cpg_ids"]
            require(src_ids.shape == dst_ids.shape, f"{modality} CpG ID shape differs.")
            for row in rows:
                require(decode(src_ids[row : row + 1]) == decode(dst_ids[row : row + 1]), f"{modality} CpG ID mismatch at row {row}.")
    print("STEP 2B INTEGRITY CHECK PASSED")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    try:
        check(args.output.resolve())
    except (ValueError, OSError, KeyError, json.JSONDecodeError, pd.errors.ParserError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
