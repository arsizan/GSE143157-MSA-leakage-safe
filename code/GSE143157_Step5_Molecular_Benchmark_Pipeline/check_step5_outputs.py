#!/usr/bin/env python3
"""Integrity checks for GSE143157 Step 5 outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


MODALITIES = {"5mC", "5hmC", "combined"}
VIEWS = {"molecular_only", "molecular_plus_covariates", "residualized_molecular"}
MODEL_IDS = {f"{modality}__{view}" for modality in MODALITIES for view in VIEWS}
RESULT_TABLES = [
    "predictions",
    "hyperparameter_choices",
    "tuning_summary",
    "outer_coefficients",
    "selected_features",
    "candidate_pool_summary",
    "repeat_metrics",
    "donor_mean_predictions",
    "aggregate_metrics",
    "feature_stability",
    "auc_comparisons",
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def check(output: Path, step3: Path, allow_smoke: bool) -> None:
    required = [
        "RUN_CONFIGURATION.json",
        "step5_summary.json",
        "SUCCESS.txt",
        "STEP5_REPORT.html",
        "performance_primary.png",
        "roc_primary.png",
        "incremental_auc_primary.png",
        "feature_stability_primary.png",
    ] + [f"{name}_primary_clean77.csv" for name in RESULT_TABLES]
    missing = [name for name in required if not (output / name).is_file()]
    require(not missing, f"Missing Step 5 files: {missing}")

    summary = json.loads((output / "step5_summary.json").read_text(encoding="utf-8"))
    parameters = summary["parameters"]
    repeats = int(parameters["outer_repeats"])
    mode = parameters["analysis_mode"]
    if not allow_smoke:
        require(mode == "FULL_RESEARCH_RUN", "This is not a full research run.")
        require(repeats == 10, "A full run must contain 10 outer repeats.")
    else:
        require(mode in {"FULL_RESEARCH_RUN", "SMOKE_TEST_NOT_FOR_RESEARCH"}, "Unknown analysis mode.")
    require(parameters["outer_folds"] == 5, "Outer-fold count is wrong.")
    require(parameters["inner_folds"] == 4, "Inner-fold count is wrong.")
    require(parameters["classification_threshold"] == 0.5, "Classification threshold is wrong.")

    frames = {name: pd.read_csv(output / f"{name}_primary_clean77.csv") for name in RESULT_TABLES}
    predictions = frames["predictions"]
    require(len(predictions) == 77 * repeats * 9, "Prediction row count is wrong.")
    require(set(predictions["model_id"]) == MODEL_IDS, "The nine molecular configurations are incomplete.")
    require(predictions["probability_msa"].between(0, 1).all(), "Invalid probabilities detected.")
    require(
        predictions[["model_id", "repeat", "donor_id"]].duplicated().sum() == 0,
        "A donor has duplicate predictions within a model/repeat.",
    )
    donor_counts = predictions.groupby(["model_id", "repeat"])["donor_id"].nunique()
    require((donor_counts == 77).all(), "Not every donor was predicted once per configuration/repeat.")
    expected_labels = np.where(predictions["probability_msa"] >= 0.5, "MSA", "Control")
    require((predictions["predicted_label"] == expected_labels).all(), "Prediction threshold mismatch.")

    outer = pd.read_csv(step3 / "outer_folds_primary_clean77.csv")
    outer = outer.loc[outer["repeat"].between(1, repeats)]
    expected_fold = outer.set_index(["repeat", "donor_id"])["outer_test_fold"]
    observed_fold = [
        expected_fold.loc[(repeat, donor)]
        for repeat, donor in zip(predictions["repeat"], predictions["donor_id"])
    ]
    require(
        np.array_equal(predictions["outer_fold"].to_numpy(), np.asarray(observed_fold)),
        "Predictions do not follow the frozen Step 3 outer folds.",
    )

    choices = frames["hyperparameter_choices"]
    require(len(choices) == 9 * repeats * 5, "Hyperparameter-choice row count is wrong.")
    require(set(choices["model_id"]) == MODEL_IDS, "Hyperparameter choices omit a configuration.")
    tuning = frames["tuning_summary"]
    grid_size = len(parameters["feature_counts"]) * len(parameters["c_grid"]) * len(parameters["class_weights"])
    require(len(tuning) == 9 * repeats * 5 * grid_size, "Tuning-grid row count is wrong.")
    require(
        tuning.groupby(["model_id", "repeat", "outer_fold"])["selected"].sum().eq(1).all(),
        "Tuning does not select exactly one setting per outer fit.",
    )
    require(len(frames["candidate_pool_summary"]) == 2 * repeats * 5, "Candidate-pool summary count is wrong.")
    require(set(frames["candidate_pool_summary"]["modality"]) == {"5mC", "5hmC"}, "Candidate modalities are wrong.")
    require(
        (frames["candidate_pool_summary"]["candidate_pool_size"] == parameters["candidate_pool"]).all(),
        "A candidate pool has the wrong size.",
    )
    expected_selected = int(choices["selected_total_molecular_features"].sum())
    require(len(frames["selected_features"]) == expected_selected, "Selected-feature row count is wrong.")
    require(not frames["outer_coefficients"].empty, "Outer coefficients are empty.")
    require(len(frames["repeat_metrics"]) == 9 * repeats, "Repeat-metric row count is wrong.")
    require(len(frames["donor_mean_predictions"]) == 77 * 9, "Donor-mean prediction count is wrong.")
    require(
        (frames["donor_mean_predictions"]["repeats"] == repeats).all(),
        "Donor means do not use every outer repeat.",
    )
    require(len(frames["aggregate_metrics"]) == 9, "Aggregate-metric row count is wrong.")
    require(len(frames["auc_comparisons"]) == 24, "Planned AUC-comparison count is wrong.")
    require(not frames["feature_stability"].empty, "Feature-stability table is empty.")
    print("STEP 5 INTEGRITY CHECK PASSED")
    if mode == "SMOKE_TEST_NOT_FOR_RESEARCH":
        print("SMOKE TEST ONLY - DO NOT USE THESE NUMBERS IN THE PAPER")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--step3-dir", required=True, type=Path)
    parser.add_argument("--allow-smoke", action="store_true")
    args = parser.parse_args()
    try:
        check(args.output.resolve(), args.step3_dir.resolve(), args.allow_smoke)
    except (ValueError, OSError, KeyError, json.JSONDecodeError, pd.errors.ParserError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
