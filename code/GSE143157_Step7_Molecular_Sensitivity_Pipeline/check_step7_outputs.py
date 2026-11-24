#!/usr/bin/env python3
"""Strict integrity checks for GSE143157 Step 7 outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


FULL78 = "sensitivity_full78"
CLEAN76 = "sensitivity_no_MSA30_Control11_76"
COHORTS = {FULL78: 78, CLEAN76: 76}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    return parser.parse_args()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read_csv(path: Path) -> pd.DataFrame:
    require(path.is_file(), f"Missing CSV: {path}")
    frame = pd.read_csv(path)
    require(len(frame) > 0, f"Empty CSV: {path}")
    return frame


def check_cohort(root: Path, cohort: str, donors: int, mode: str) -> None:
    folder = root / cohort
    require((folder / "SUCCESS.txt").is_file(), f"Missing cohort SUCCESS: {cohort}")
    summary_path = folder / "cohort_summary.json"
    require(summary_path.is_file(), f"Missing cohort summary: {cohort}")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    require(summary["parameters"]["analysis_mode"] == mode, f"Analysis mode mismatch: {cohort}")
    require(summary["cohort"]["name"] == cohort, f"Cohort-name mismatch: {cohort}")
    require(int(summary["cohort"]["donors"]) == donors, f"Donor-count mismatch: {cohort}")
    require(set(summary["models"]) == {
        "5mC__molecular_plus_covariates",
        "5hmC__molecular_plus_covariates",
        "combined__molecular_plus_covariates",
    }, f"Unexpected models: {cohort}")

    predictions = read_csv(folder / f"predictions_{cohort}.csv")
    donor_means = read_csv(folder / f"donor_mean_predictions_{cohort}.csv")
    repeats = read_csv(folder / f"repeat_metrics_{cohort}.csv")
    aggregate = read_csv(folder / f"aggregate_metrics_{cohort}.csv")
    choices = read_csv(folder / f"hyperparameter_choices_{cohort}.csv")
    tuning = read_csv(folder / f"tuning_summary_{cohort}.csv")
    stability = read_csv(folder / f"feature_stability_{cohort}.csv")
    candidate_pool = read_csv(folder / f"candidate_pool_summary_{cohort}.csv")
    read_csv(folder / f"selected_features_{cohort}.csv")
    read_csv(folder / f"outer_coefficients_{cohort}.csv")
    read_csv(folder / f"auc_comparisons_{cohort}.csv")

    expected_repeats = 10 if mode == "FULL_RESEARCH_RUN" else 1
    require(len(predictions) == donors * expected_repeats * 3, f"Prediction row count is wrong: {cohort}")
    require(not predictions.duplicated(["model_id", "repeat", "donor_id"]).any(), f"Duplicate prediction: {cohort}")
    require(predictions.groupby(["model_id", "repeat"])["donor_id"].nunique().eq(donors).all(), f"Incomplete donor predictions: {cohort}")
    require(set(predictions["outer_fold"]) == {1, 2, 3, 4, 5}, f"Outer folds incomplete: {cohort}")
    require(len(donor_means) == donors * 3, f"Donor-mean row count is wrong: {cohort}")
    require(len(repeats) == expected_repeats * 3, f"Repeat metrics count is wrong: {cohort}")
    require(len(aggregate) == 3 and aggregate["model_id"].nunique() == 3, f"Aggregate models incomplete: {cohort}")
    require(len(choices) == expected_repeats * 5 * 3, f"Hyperparameter choice count is wrong: {cohort}")
    require(choices.groupby(["repeat", "outer_fold", "model_id"]).size().eq(1).all(), f"Duplicate hyperparameter choice: {cohort}")
    require(tuning.groupby(["repeat", "outer_fold", "model_id"])["selected"].apply(lambda values: values.astype(str).str.lower().eq("true").sum()).eq(1).all(), f"Tuning selections are invalid: {cohort}")
    require(len(candidate_pool) == expected_repeats * 5 * 2, f"Candidate-pool count is wrong: {cohort}")
    require(stability["selection_frequency"].between(0, 1).all(), f"Invalid stability frequency: {cohort}")
    require(not stability.duplicated(["model_id", "feature_name"]).any(), f"Duplicate stability feature: {cohort}")
    for column in ("roc_auc_donor_mean_oof", "pr_auc_donor_mean_oof", "balanced_accuracy_donor_mean_oof"):
        values = pd.to_numeric(aggregate[column], errors="coerce")
        require(np.isfinite(values).all() and values.between(0, 1).all(), f"Invalid {column}: {cohort}")
    for name in (f"performance_{cohort}.png", f"roc_{cohort}.png", f"incremental_auc_{cohort}.png", f"feature_stability_{cohort}.png", "COHORT_REPORT.html"):
        path = folder / name
        require(path.is_file() and path.stat().st_size > 1000, f"Missing or tiny cohort artifact: {path}")


def run(output: Path) -> None:
    root = output.resolve()
    required_root = [
        "SUCCESS.txt", "RUN_CONFIGURATION.json", "step7_summary.json", "STEP7_REPORT.html",
        "aggregate_metrics_across_cohorts.csv", "primary_candidate_sensitivity.csv", "stable_set_overlap.csv",
        "performance_across_cohorts.png", "primary_candidate_selection_heatmap.png",
        "candidate_effect_concordance.png", "stable_set_overlap.png",
    ]
    for name in required_root:
        path = root / name
        require(path.is_file() and path.stat().st_size > 0, f"Missing root output: {path}")
    summary = json.loads((root / "step7_summary.json").read_text(encoding="utf-8"))
    mode = summary["analysis_mode"]
    require(mode in {"FULL_RESEARCH_RUN", "SMOKE_TEST_NOT_FOR_RESEARCH"}, "Invalid Step 7 analysis mode.")
    for cohort, donors in COHORTS.items():
        check_cohort(root, cohort, donors, mode)

    performance = read_csv(root / "aggregate_metrics_across_cohorts.csv")
    tracking = read_csv(root / "primary_candidate_sensitivity.csv")
    overlap = read_csv(root / "stable_set_overlap.csv")
    require(len(performance) == 12, "Across-cohort performance must contain 3 cohorts x 4 models.")
    require(performance.groupby("cohort")["model_id"].nunique().eq(4).all(), "Across-cohort models are incomplete.")
    require(len(tracking) == 20 and tracking["feature_name"].nunique() == 20, "Primary candidate tracking must contain 20 unique CpGs.")
    for column in ("primary_selection_frequency", "full78_selection_frequency", "clean76_selection_frequency", "minimum_sensitivity_selection_frequency"):
        require(pd.to_numeric(tracking[column], errors="coerce").between(0, 1).all(), f"Invalid tracking frequency: {column}")
    calculated = tracking[["full78_selection_frequency", "clean76_selection_frequency"]].min(axis=1)
    require(np.allclose(calculated, tracking["minimum_sensitivity_selection_frequency"]), "Minimum sensitivity frequency is inconsistent.")
    expected_support = (
        (tracking["minimum_sensitivity_selection_frequency"] >= 0.60)
        & tracking["effect_direction_matches_both"].astype(str).str.lower().eq("true")
    )
    observed_support = tracking["sensitivity_supported_both"].astype(str).str.lower().eq("true")
    require((expected_support == observed_support).all(), "Sensitivity-support classification is inconsistent.")
    require(len(overlap) == 9 and overlap["jaccard"].between(0, 1).all(), "Stable-set overlap table is invalid.")
    require(int(summary["counts"]["primary_candidates"]) == 20, "Summary candidate count is invalid.")
    require(int(summary["counts"]["primary_candidates_supported_in_both"]) == int(observed_support.sum()), "Summary support count is inconsistent.")
    print("STEP 7 INTEGRITY CHECK PASSED")
    if mode != "FULL_RESEARCH_RUN":
        print("SMOKE TEST ONLY - DO NOT USE THESE NUMBERS IN THE PAPER")


def main() -> int:
    try:
        run(parse_args().output)
    except (ValueError, OSError, KeyError, json.JSONDecodeError, pd.errors.ParserError) as exc:
        print(f"STEP 7 INTEGRITY CHECK FAILED: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
