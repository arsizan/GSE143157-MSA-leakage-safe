#!/usr/bin/env python3
"""Quick integrity check for the GSE143157 Step 1 output folder."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate GSE143157 Step 1 outputs")
    parser.add_argument("output", type=Path, help="Step 1 result folder")
    args = parser.parse_args()
    folder = args.output.resolve()

    required = [
        "SUCCESS.txt",
        "qc_summary.json",
        "GSE143157_step1_matrices.h5",
        "metadata_arrays_aligned.csv",
        "metadata_donors_aligned.csv",
        "qc_per_array.csv",
        "qc_per_donor.csv",
        "5hmC_probe_qc.tsv.gz",
    ]
    missing = [name for name in required if not (folder / name).is_file()]
    if missing:
        print(f"FAILED: missing files: {missing}", file=sys.stderr)
        return 2

    with (folder / "qc_summary.json").open(encoding="utf-8") as handle:
        summary = json.load(handle)
    donors = pd.read_csv(folder / "metadata_donors_aligned.csv")
    arrays = pd.read_csv(folder / "metadata_arrays_aligned.csv")

    with h5py.File(folder / "GSE143157_step1_matrices.h5", "r") as h5:
        shape_5mc = tuple(h5["5mC/values"].shape)
        shape_hmc = tuple(h5["5hmC/values"].shape)
        ids_5mc = h5["5mC/cpg_ids"].shape[0]
        ids_hmc = h5["5hmC/cpg_ids"].shape[0]
        donor_ids = h5["donor_ids"].shape[0]

    expected_5mc = (
        summary["outputs_summary"]["5mC_cpgs"],
        summary["outputs_summary"]["5mC_donors"],
    )
    expected_hmc = (
        summary["outputs_summary"]["5hmC_cpgs_retained"],
        summary["outputs_summary"]["5hmC_donors"],
    )
    errors = []
    if shape_5mc != expected_5mc or ids_5mc != shape_5mc[0]:
        errors.append(f"5mC mismatch: HDF5={shape_5mc}, IDs={ids_5mc}, summary={expected_5mc}")
    if shape_hmc != expected_hmc or ids_hmc != shape_hmc[0]:
        errors.append(f"5hmC mismatch: HDF5={shape_hmc}, IDs={ids_hmc}, summary={expected_hmc}")
    if len(donors) != donor_ids or len(donors) != shape_5mc[1]:
        errors.append(f"Donor mismatch: CSV={len(donors)}, HDF5 IDs={donor_ids}, matrix={shape_5mc[1]}")
    if len(arrays) != 2 * len(donors):
        errors.append(f"Array pairing mismatch: arrays={len(arrays)}, donors={len(donors)}")

    if errors:
        print("FAILED")
        for error in errors:
            print(f"- {error}")
        return 2

    print("PASS: Step 1 outputs are internally consistent.")
    print(f"5mC shape:  {shape_5mc[0]:,} CpGs x {shape_5mc[1]} donors")
    print(f"5hmC shape: {shape_hmc[0]:,} CpGs x {shape_hmc[1]} donors")
    print(f"Condition counts: {donors['condition'].value_counts().to_dict()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
