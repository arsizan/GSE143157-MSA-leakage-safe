#!/usr/bin/env python3
"""Integrity checks for GSE143157 Step 4 outputs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


COHORTS = {
    "primary_clean77": 77,
    "sensitivity_full78": 78,
    "sensitivity_no_MSA30_Control11_76": 76,
}
MODELS = {"age_only", "no_age_covariates", "all_covariates"}
RESULT_TYPES = [
    "predictions",
    "hyperparameter_choices",
    "tuning_scores",
    "outer_coefficients",
    "coefficient_stability",
    "repeat_metrics",
    "donor_mean_predictions",
    "aggregate_metrics",
    "pairwise_auc_comparisons",
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def check_cohort(output: Path, step3: Path, cohort: str, donors: int) -> None:
    frames = {
        result: pd.read_csv(output / f"{result}_{cohort}.csv") for result in RESULT_TYPES
    }
    predictions = frames["predictions"]
    expected_prediction_rows = donors * 10 * len(MODELS)
    require(len(predictions) == expected_prediction_rows, f"{cohort}: wrong prediction row count.")
    require(set(predictions["model"]) == MODELS, f"{cohort}: wrong models in predictions.")
    require(predictions["probability_msa"].between(0, 1).all(), f"{cohort}: invalid probabilities.")
    require(predictions[["model", "repeat", "donor_id"]].duplicated().sum() == 0, f"{cohort}: duplicate donor prediction within a repeat/model.")
    counts = predictions.groupby(["model", "repeat"])["donor_id"].nunique()
    require((counts == donors).all(), f"{cohort}: not every donor was predicted once per repeat/model.")
    require((predictions["predicted_label"] == np.where(predictions["probability_msa"] >= 0.5, "MSA", "Control")).all(), f"{cohort}: prediction threshold mismatch.")

    outer = pd.read_csv(step3 / f"outer_folds_{cohort}.csv")
    expected_fold = outer.set_index(["repeat", "donor_id"])["outer_test_fold"]
    observed = predictions.copy()
    observed["expected_fold"] = [
        expected_fold.loc[(repeat, donor)]
        for repeat, donor in zip(observed["repeat"], observed["donor_id"])
    ]
    require((observed["outer_fold"] == observed["expected_fold"]).all(), f"{cohort}: predictions do not use frozen outer folds.")

    choices = frames["hyperparameter_choices"]
    require(len(choices) == 3 * 10 * 5, f"{cohort}: wrong hyperparameter-choice row count.")
    require(set(choices["selected_C"]) <= {0.01, 0.1, 1.0, 10.0, 100.0}, f"{cohort}: unexpected C selected.")
    require(set(choices["selected_class_weight"]) <= {"none", "balanced"}, f"{cohort}: unexpected class weight.")
    tuning = frames["tuning_scores"]
    require(len(tuning) == 3 * 10 * 5 * 10, f"{cohort}: wrong tuning-grid row count.")
    require(tuning.groupby(["model", "repeat", "outer_fold"])["selected"].sum().eq(1).all(), f"{cohort}: tuning selections are not one per outer fit.")
    require(len(frames["repeat_metrics"]) == 3 * 10, f"{cohort}: wrong repeat-metric row count.")
    require(len(frames["donor_mean_predictions"]) == donors * 3, f"{cohort}: wrong donor-mean row count.")
    require((frames["donor_mean_predictions"]["repeats"] == 10).all(), f"{cohort}: donor means do not use 10 repeats.")
    require(len(frames["aggregate_metrics"]) == 3, f"{cohort}: wrong aggregate-metric row count.")
    require(len(frames["pairwise_auc_comparisons"]) == 3, f"{cohort}: wrong pairwise comparison count.")
    require(not frames["outer_coefficients"].empty, f"{cohort}: coefficients are empty.")
    require(not frames["coefficient_stability"].empty, f"{cohort}: coefficient stability is empty.")


def check(output: Path, step3: Path) -> None:
    required = [
        "SUCCESS.txt",
        "STEP4_REPORT.html",
        "step4_summary.json",
        "roc_primary.png",
        "pr_primary.png",
        "metrics_primary.png",
        "coefficients_primary.png",
    ]
    for cohort in COHORTS:
        required.extend(f"{result}_{cohort}.csv" for result in RESULT_TYPES)
    missing = [name for name in required if not (output / name).is_file()]
    require(not missing, f"Missing Step 4 files: {missing}")
    summary = json.loads((output / "step4_summary.json").read_text(encoding="utf-8"))
    require(summary["parameters"]["outer_folds"] == 5, "Outer-fold count is wrong.")
    require(summary["parameters"]["outer_repeats"] == 10, "Outer-repeat count is wrong.")
    require(summary["parameters"]["inner_folds"] == 4, "Inner-fold count is wrong.")
    require(summary["parameters"]["classification_threshold"] == 0.5, "Threshold is wrong.")
    for cohort, donors in COHORTS.items():
        check_cohort(output, step3, cohort, donors)
    print("STEP 4 INTEGRITY CHECK PASSED")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--step3-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        check(args.output.resolve(), args.step3_dir.resolve())
    except (ValueError, OSError, KeyError, json.JSONDecodeError, pd.errors.ParserError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
