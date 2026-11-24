#!/usr/bin/env python3
"""Integrity checks for GSE143157 Step 3 outputs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


COHORTS = {
    "primary_clean77": (77, {"MSA": 40, "Control": 37}),
    "sensitivity_full78": (78, {"MSA": 41, "Control": 37}),
    "sensitivity_no_MSA30_Control11_76": (76, {"MSA": 40, "Control": 36}),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def partition_signature(frame: pd.DataFrame) -> tuple[str, ...]:
    parts = []
    for fold in sorted(frame["outer_test_fold"].unique()):
        parts.append(",".join(sorted(frame.loc[frame["outer_test_fold"] == fold, "donor_id"])))
    return tuple(sorted(parts))


def check_cohort(output: Path, cohort: str, donors: int, conditions: dict[str, int], summary: dict) -> None:
    audit = pd.read_csv(output / f"covariate_audit_{cohort}.csv")
    outer = pd.read_csv(output / f"outer_folds_{cohort}.csv")
    inner = pd.read_csv(output / f"inner_folds_{cohort}.csv")
    balance = pd.read_csv(output / f"outer_fold_balance_{cohort}.csv")
    require(set(audit["covariate"]) == {"age", "sex", "brain_bank", "pmi_hours"}, f"{cohort}: audit covariates are wrong.")
    require(audit["fdr_bh"].between(0, 1).all(), f"{cohort}: invalid FDR values.")
    repeats = int(summary["parameters"]["outer_repeats"])
    outer_folds = int(summary["parameters"]["outer_folds"])
    inner_folds = int(summary["parameters"]["inner_folds"])
    require(len(outer) == donors * repeats, f"{cohort}: wrong outer assignment row count.")
    require(len(balance) == repeats * outer_folds, f"{cohort}: wrong outer balance row count.")
    signatures = set()
    all_donors = set(outer["donor_id"])
    require(len(all_donors) == donors, f"{cohort}: wrong unique donor count.")
    require(outer.drop_duplicates("donor_id")["condition"].value_counts().to_dict() == conditions, f"{cohort}: wrong labels.")
    for repeat in range(1, repeats + 1):
        current = outer.loc[outer["repeat"] == repeat]
        require(len(current) == donors, f"{cohort}: repeat {repeat} row count is wrong.")
        require(current["donor_id"].nunique() == donors, f"{cohort}: repeat {repeat} duplicates/misses donors.")
        require(set(current["outer_test_fold"]) == set(range(1, outer_folds + 1)), f"{cohort}: repeat {repeat} fold labels are wrong.")
        signatures.add(partition_signature(current))
        for fold in range(1, outer_folds + 1):
            test = current.loc[current["outer_test_fold"] == fold]
            require(set(test["condition"]) == {"Control", "MSA"}, f"{cohort}: an outer test fold lacks a class.")
            inner_current = inner.loc[(inner["outer_repeat"] == repeat) & (inner["outer_fold"] == fold)]
            expected_train = all_donors - set(test["donor_id"])
            require(set(inner_current["donor_id"]) == expected_train, f"{cohort}: inner donors do not equal outer training donors.")
            require(inner_current["donor_id"].nunique() == len(inner_current), f"{cohort}: duplicate inner donor assignment.")
            require(set(inner_current["inner_validation_fold"]) == set(range(1, inner_folds + 1)), f"{cohort}: inner fold labels are wrong.")
            for inner_fold in range(1, inner_folds + 1):
                validation = inner_current.loc[inner_current["inner_validation_fold"] == inner_fold]
                require(set(validation["condition"]) == {"Control", "MSA"}, f"{cohort}: an inner fold lacks a class.")
    require(len(signatures) == repeats, f"{cohort}: repeated outer partitions are not unique.")


def check(output: Path) -> None:
    required = [
        "SUCCESS.txt",
        "STEP3_REPORT.html",
        "step3_summary.json",
        "model_input_policy.json",
        "covariate_balance_primary.png",
        "outer_fold_balance_primary.png",
    ]
    for cohort in COHORTS:
        required.extend(
            [
                f"covariate_audit_{cohort}.csv",
                f"outer_folds_{cohort}.csv",
                f"inner_folds_{cohort}.csv",
                f"outer_fold_balance_{cohort}.csv",
            ]
        )
    missing = [name for name in required if not (output / name).is_file()]
    require(not missing, f"Missing Step 3 files: {missing}")
    summary = json.loads((output / "step3_summary.json").read_text(encoding="utf-8"))
    policy = json.loads((output / "model_input_policy.json").read_text(encoding="utf-8"))
    require(policy["primary_cohort"] == "primary_clean77", "Primary policy cohort is wrong.")
    require(set(policy["same_outer_and_inner_donor_folds_for_modalities"]) == {"5mC", "5hmC", "combined"}, "Modality fold policy is wrong.")
    for cohort, (donors, conditions) in COHORTS.items():
        check_cohort(output, cohort, donors, conditions, summary)
    print("STEP 3 INTEGRITY CHECK PASSED")


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
