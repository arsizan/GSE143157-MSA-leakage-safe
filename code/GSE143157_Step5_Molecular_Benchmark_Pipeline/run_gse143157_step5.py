#!/usr/bin/env python3
"""Run leakage-safe primary-cohort molecular benchmarks for GSE143157."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import html
import json
import os
from pathlib import Path
import pickle
import sys
import time

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from sklearn.metrics import precision_recall_curve, roc_curve, auc

from step5_core import (
    MODEL_SPECS,
    MODALITIES,
    VIEWS,
    benjamini_hochberg,
    binary_metrics,
    build_design,
    decode_strings,
    fit_covariate_transform,
    fit_logistic,
    paired_auc_bootstrap,
    prepare_candidates,
    stratified_bootstrap_metrics,
    stream_variance_candidates,
)


VERSION = "1.0.0"
PRIMARY_DONORS = 77
OUTER_FOLDS = 5
INNER_FOLDS = 4
MODALITY_COLORS = {"5mC": "#2387bd", "5hmC": "#1aa77b", "combined": "#dc6d11"}
VIEW_LABELS = {
    "molecular_only": "Molecular only",
    "molecular_plus_covariates": "Molecular + covariates",
    "residualized_molecular": "Covariate-residualized molecular",
}


def parse_csv_numbers(text: str, cast) -> list:
    values = [cast(item.strip()) for item in text.split(",") if item.strip()]
    if not values:
        raise ValueError("A tuning grid cannot be empty.")
    return sorted(set(values))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step2b-dir", required=True, type=Path)
    parser.add_argument("--step3-dir", required=True, type=Path)
    parser.add_argument("--step4-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--candidate-pool", type=int, default=5000)
    parser.add_argument("--feature-counts", default="25,50,100")
    parser.add_argument("--c-grid", default="0.01,0.1,1,10")
    parser.add_argument("--class-weights", default="none,balanced")
    parser.add_argument("--chunk-rows", type=int, default=5000)
    parser.add_argument("--missing-threshold", type=float, default=0.20)
    parser.add_argument("--minimum-variance", type=float, default=1e-8)
    parser.add_argument("--bootstrap-iterations", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=143157)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Small one-repeat software test only; never use smoke output as a research result.",
    )
    return parser.parse_args()


def file_sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def json_clean(value):
    if isinstance(value, dict):
        return {str(key): json_clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_clean(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def model_id(modality: str, view: str) -> str:
    return f"{modality}__{view}"


def validate_inputs(step2b: Path, step3: Path, step4: Path) -> dict[str, object]:
    paths = {
        "h5": step2b / "GSE143157_primary_clean77.h5",
        "metadata": step2b / "metadata_donors_primary_clean77.csv",
        "cohort_decision": step2b / "cohort_decision.json",
        "step2b_success": step2b / "SUCCESS.txt",
        "outer": step3 / "outer_folds_primary_clean77.csv",
        "inner": step3 / "inner_folds_primary_clean77.csv",
        "policy": step3 / "model_input_policy.json",
        "step3_success": step3 / "SUCCESS.txt",
        "step4_predictions": step4 / "donor_mean_predictions_primary_clean77.csv",
        "step4_metrics": step4 / "aggregate_metrics_primary_clean77.csv",
        "step4_success": step4 / "SUCCESS.txt",
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise ValueError(f"Required input files are missing: {missing}")

    metadata = pd.read_csv(paths["metadata"])
    outer = pd.read_csv(paths["outer"])
    inner = pd.read_csv(paths["inner"])
    step4_predictions = pd.read_csv(paths["step4_predictions"])
    policy = json.loads(paths["policy"].read_text(encoding="utf-8"))

    if len(metadata) != PRIMARY_DONORS or metadata["donor_id"].nunique() != PRIMARY_DONORS:
        raise ValueError("Primary metadata does not contain 77 unique donors.")
    if "MSA 30" in set(metadata["donor_id"]):
        raise ValueError("MSA 30 is unexpectedly present in the primary Clean-77 metadata.")
    if set(metadata["condition"]) != {"Control", "MSA"}:
        raise ValueError("Primary condition labels are invalid.")
    if sorted(outer["repeat"].unique()) != list(range(1, 11)):
        raise ValueError("Step 3 outer folds do not contain repeats 1-10.")
    if sorted(outer["outer_test_fold"].unique()) != list(range(1, OUTER_FOLDS + 1)):
        raise ValueError("Step 3 outer fold labels are invalid.")
    if outer.groupby(["repeat", "donor_id"]).size().ne(1).any():
        raise ValueError("A donor is duplicated or absent within a Step 3 outer repeat.")
    expected_views = {
        "covariate_only_baseline",
        "molecular_only",
        "molecular_plus_covariates",
        "training_fold_residualized_molecular_sensitivity",
    }
    if not expected_views.issubset(set(policy.get("required_model_views", []))):
        raise ValueError("Step 3 model policy does not contain the required Step 5 views.")
    needed_step4 = step4_predictions[step4_predictions["model"].isin(["age_only", "all_covariates"])]
    if len(needed_step4) != PRIMARY_DONORS * 2:
        raise ValueError("Step 4 donor-mean benchmark predictions are incomplete.")

    with h5py.File(paths["h5"], "r") as h5:
        required_h5 = ["donor_ids", "5mC/values", "5mC/cpg_ids", "5hmC/values", "5hmC/cpg_ids"]
        missing_h5 = [name for name in required_h5 if name not in h5]
        if missing_h5:
            raise ValueError(f"Primary HDF5 is missing datasets: {missing_h5}")
        h5_donors = decode_strings(h5["donor_ids"][:])
        if h5_donors != metadata["donor_id"].astype(str).tolist():
            raise ValueError("HDF5 donor order does not match primary metadata.")
        if h5["5mC/values"].shape != (731661, PRIMARY_DONORS):
            raise ValueError(f"Unexpected 5mC shape: {h5['5mC/values'].shape}")
        if h5["5hmC/values"].shape != (405408, PRIMARY_DONORS):
            raise ValueError(f"Unexpected 5hmC shape: {h5['5hmC/values'].shape}")

    return {
        "paths": paths,
        "metadata": metadata.reset_index(drop=True),
        "outer": outer,
        "inner": inner,
        "step4_predictions": step4_predictions,
    }


def prepare_split(
    pools,
    metadata: pd.DataFrame,
    labels: np.ndarray,
    train_indices: np.ndarray,
    evaluation_indices: np.ndarray,
    maximum_features: int,
    missing_threshold: float,
    minimum_variance: float,
) -> dict[str, object]:
    cov_train, cov_evaluate, cov_names = fit_covariate_transform(
        metadata, train_indices, evaluation_indices
    )
    raw = {}
    residualized = {}
    for modality in ("5mC", "5hmC"):
        raw[modality] = prepare_candidates(
            pools[modality],
            train_indices,
            evaluation_indices,
            labels[train_indices],
            metadata,
            maximum_features,
            False,
            missing_threshold,
            minimum_variance,
        )
        residualized[modality] = prepare_candidates(
            pools[modality],
            train_indices,
            evaluation_indices,
            labels[train_indices],
            metadata,
            maximum_features,
            True,
            missing_threshold,
            minimum_variance,
        )
    return {
        "train_indices": train_indices,
        "evaluation_indices": evaluation_indices,
        "labels_train": labels[train_indices],
        "labels_evaluate": labels[evaluation_indices],
        "raw": raw,
        "residualized": residualized,
        "covariate_train": cov_train,
        "covariate_evaluate": cov_evaluate,
        "covariate_names": cov_names,
    }


def design_from_split(split: dict[str, object], modality: str, view: str, feature_count: int):
    return build_design(
        modality,
        view,
        feature_count,
        split["raw"],
        split["residualized"],
        split["covariate_train"],
        split["covariate_evaluate"],
        split["covariate_names"],
    )


def select_hyperparameters(
    inner_splits: list[dict[str, object]],
    feature_counts: list[int],
    c_grid: list[float],
    class_weights: list[str],
    repeat: int,
    outer_fold: int,
    seed: int,
) -> tuple[list[dict[str, object]], dict[tuple[str, str], dict[str, object]]]:
    tuning_rows: list[dict[str, object]] = []
    selected: dict[tuple[str, str], dict[str, object]] = {}
    for spec_index, (modality, view) in enumerate(MODEL_SPECS):
        start = len(tuning_rows)
        for feature_count in feature_counts:
            for c_value in c_grid:
                for class_weight in class_weights:
                    fold_scores = []
                    for inner_index, split in enumerate(inner_splits, start=1):
                        train, evaluate, _ = design_from_split(split, modality, view, feature_count)
                        model = fit_logistic(
                            train,
                            split["labels_train"],
                            c_value,
                            class_weight,
                            seed + spec_index * 10000 + inner_index,
                        )
                        probabilities = model.predict_proba(evaluate)[:, 1]
                        fold_scores.append(
                            binary_metrics(split["labels_evaluate"], probabilities)["roc_auc"]
                        )
                    tuning_rows.append(
                        {
                            "repeat": repeat,
                            "outer_fold": outer_fold,
                            "modality": modality,
                            "view": view,
                            "model_id": model_id(modality, view),
                            "feature_count_per_modality": feature_count,
                            "total_molecular_features": feature_count * (2 if modality == "combined" else 1),
                            "C": c_value,
                            "class_weight": class_weight,
                            "mean_inner_roc_auc": float(np.mean(fold_scores)),
                            "sd_inner_roc_auc": float(np.std(fold_scores, ddof=1)),
                            "selected": False,
                        }
                    )
        candidates = list(range(start, len(tuning_rows)))
        best_index = min(
            candidates,
            key=lambda index: (
                -tuning_rows[index]["mean_inner_roc_auc"],
                tuning_rows[index]["sd_inner_roc_auc"],
                tuning_rows[index]["feature_count_per_modality"],
                tuning_rows[index]["C"],
                0 if tuning_rows[index]["class_weight"] == "none" else 1,
            ),
        )
        tuning_rows[best_index]["selected"] = True
        selected[(modality, view)] = tuning_rows[best_index]
    return tuning_rows, selected


def run_outer_fit(
    h5: h5py.File,
    metadata: pd.DataFrame,
    labels: np.ndarray,
    outer: pd.DataFrame,
    inner: pd.DataFrame,
    repeat: int,
    outer_fold: int,
    parameters: dict[str, object],
) -> dict[str, list[dict[str, object]]]:
    started = time.time()
    donor_to_index = {donor: index for index, donor in enumerate(metadata["donor_id"].astype(str))}
    test_donors = outer.loc[
        (outer["repeat"] == repeat) & (outer["outer_test_fold"] == outer_fold), "donor_id"
    ].astype(str).tolist()
    test_indices = np.asarray([donor_to_index[donor] for donor in test_donors], dtype=np.int64)
    train_indices = np.asarray(
        [index for index in range(len(metadata)) if index not in set(test_indices)], dtype=np.int64
    )
    if set(metadata.iloc[train_indices]["donor_id"]) & set(test_donors):
        raise ValueError("Outer train/test donor overlap detected.")

    print(f"  selecting outer-training variance pools", flush=True)
    pools = {}
    candidate_rows = []
    for modality in ("5mC", "5hmC"):
        pools[modality] = stream_variance_candidates(
            h5,
            modality,
            train_indices,
            parameters["candidate_pool"],
            parameters["chunk_rows"],
            parameters["missing_threshold"],
            parameters["minimum_variance"],
        )
        pool = pools[modality]
        candidate_rows.append(
            {
                "repeat": repeat,
                "outer_fold": outer_fold,
                "modality": modality,
                "total_cpgs": pool.total_features,
                "eligible_training_cpgs": pool.eligible_features,
                "candidate_pool_size": len(pool.row_indices),
                "largest_training_variance": float(pool.outer_train_variances.max()),
                "smallest_selected_training_variance": float(pool.outer_train_variances.min()),
            }
        )

    inner_table = inner.loc[
        (inner["outer_repeat"] == repeat) & (inner["outer_fold"] == outer_fold)
    ].copy()
    expected_train_donors = set(metadata.iloc[train_indices]["donor_id"].astype(str))
    if set(inner_table["donor_id"].astype(str)) != expected_train_donors:
        raise ValueError(f"Inner assignments do not equal outer training donors for repeat {repeat}, fold {outer_fold}.")
    inner_splits = []
    maximum_features = max(parameters["feature_counts"])
    for inner_fold in range(1, INNER_FOLDS + 1):
        validation_donors = inner_table.loc[
            inner_table["inner_validation_fold"] == inner_fold, "donor_id"
        ].astype(str).tolist()
        validation_indices = np.asarray([donor_to_index[donor] for donor in validation_donors], dtype=np.int64)
        inner_train_indices = np.asarray(
            [index for index in train_indices if index not in set(validation_indices)], dtype=np.int64
        )
        inner_splits.append(
            prepare_split(
                pools,
                metadata,
                labels,
                inner_train_indices,
                validation_indices,
                maximum_features,
                parameters["missing_threshold"],
                parameters["minimum_variance"],
            )
        )
    print("  tuning 9 molecular configurations inside frozen inner folds", flush=True)
    tuning_rows, selected = select_hyperparameters(
        inner_splits,
        parameters["feature_counts"],
        parameters["c_grid"],
        parameters["class_weights"],
        repeat,
        outer_fold,
        parameters["seed"] + repeat * 100000 + outer_fold * 1000,
    )

    outer_split = prepare_split(
        pools,
        metadata,
        labels,
        train_indices,
        test_indices,
        maximum_features,
        parameters["missing_threshold"],
        parameters["minimum_variance"],
    )
    prediction_rows: list[dict[str, object]] = []
    choice_rows: list[dict[str, object]] = []
    coefficient_rows: list[dict[str, object]] = []
    selected_feature_rows: list[dict[str, object]] = []
    for config_index, (modality, view) in enumerate(MODEL_SPECS):
        choice = selected[(modality, view)]
        feature_count = int(choice["feature_count_per_modality"])
        train, test, feature_records = design_from_split(outer_split, modality, view, feature_count)
        model = fit_logistic(
            train,
            labels[train_indices],
            choice["C"],
            choice["class_weight"],
            parameters["seed"] + repeat * 100000 + outer_fold * 1000 + config_index,
        )
        probabilities = model.predict_proba(test)[:, 1]
        for donor_index, probability in zip(test_indices, probabilities):
            prediction_rows.append(
                {
                    "repeat": repeat,
                    "outer_fold": outer_fold,
                    "modality": modality,
                    "view": view,
                    "model_id": model_id(modality, view),
                    "donor_index": int(donor_index),
                    "donor_id": str(metadata.iloc[donor_index]["donor_id"]),
                    "condition": str(metadata.iloc[donor_index]["condition"]),
                    "y_true": int(labels[donor_index]),
                    "probability_msa": float(probability),
                    "predicted_label": "MSA" if probability >= 0.5 else "Control",
                }
            )
        choice_rows.append(
            {
                "repeat": repeat,
                "outer_fold": outer_fold,
                "modality": modality,
                "view": view,
                "model_id": model_id(modality, view),
                "selected_feature_count_per_modality": feature_count,
                "selected_total_molecular_features": feature_count * (2 if modality == "combined" else 1),
                "selected_C": float(choice["C"]),
                "selected_class_weight": choice["class_weight"],
                "selected_mean_inner_roc_auc": float(choice["mean_inner_roc_auc"]),
                "selected_sd_inner_roc_auc": float(choice["sd_inner_roc_auc"]),
            }
        )
        coefficients = model.coef_[0]
        if len(coefficients) != len(feature_records):
            raise ValueError("Coefficient count does not match final design columns.")
        for record, coefficient in zip(feature_records, coefficients):
            row = {
                "repeat": repeat,
                "outer_fold": outer_fold,
                "modality": modality,
                "view": view,
                "model_id": model_id(modality, view),
                **record,
                "coefficient": float(coefficient),
            }
            coefficient_rows.append(row)
            if record["feature_type"] == "molecular":
                selected_feature_rows.append(row.copy())

    runtime = time.time() - started
    for row in candidate_rows:
        row["outer_fit_elapsed_seconds"] = runtime
    return {
        "predictions": prediction_rows,
        "hyperparameter_choices": choice_rows,
        "tuning_summary": tuning_rows,
        "outer_coefficients": coefficient_rows,
        "selected_features": selected_feature_rows,
        "candidate_pool_summary": candidate_rows,
    }


def make_repeat_metrics(predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (modality, view, repeat), frame in predictions.groupby(["modality", "view", "repeat"]):
        result = binary_metrics(frame["y_true"].to_numpy(), frame["probability_msa"].to_numpy())
        rows.append({"modality": modality, "view": view, "model_id": model_id(modality, view), "repeat": repeat, **result})
    return pd.DataFrame(rows).sort_values(["view", "modality", "repeat"]).reset_index(drop=True)


def make_donor_means(predictions: pd.DataFrame) -> pd.DataFrame:
    frame = (
        predictions.groupby(["modality", "view", "model_id", "donor_id", "condition", "y_true"], as_index=False)
        .agg(probability_msa_mean=("probability_msa", "mean"), probability_msa_sd=("probability_msa", "std"), repeats=("repeat", "nunique"))
    )
    frame["predicted_label"] = np.where(frame["probability_msa_mean"] >= 0.5, "MSA", "Control")
    return frame.sort_values(["view", "modality", "donor_id"]).reset_index(drop=True)


def make_aggregate_metrics(
    donor_means: pd.DataFrame,
    repeat_metrics: pd.DataFrame,
    bootstrap_iterations: int,
    seed: int,
) -> pd.DataFrame:
    rows = []
    metric_names = list(binary_metrics(np.array([0, 1]), np.array([0.25, 0.75])))
    for config_index, ((modality, view), frame) in enumerate(donor_means.groupby(["modality", "view"])):
        labels = frame["y_true"].to_numpy()
        probabilities = frame["probability_msa_mean"].to_numpy()
        point = binary_metrics(labels, probabilities)
        intervals = stratified_bootstrap_metrics(
            labels, probabilities, bootstrap_iterations, seed + config_index * 100003
        )
        repeat_frame = repeat_metrics.loc[
            (repeat_metrics["modality"] == modality) & (repeat_metrics["view"] == view)
        ]
        row = {"modality": modality, "view": view, "model_id": model_id(modality, view), "donors": len(frame), "outer_repeats": int(frame["repeats"].iloc[0])}
        for name in metric_names:
            row[f"{name}_repeat_mean"] = float(repeat_frame[name].mean())
            row[f"{name}_repeat_sd"] = float(repeat_frame[name].std(ddof=1)) if len(repeat_frame) > 1 else 0.0
            row[f"{name}_donor_mean_oof"] = point[name]
            row[f"{name}_bootstrap_low"] = intervals[name][0]
            row[f"{name}_bootstrap_high"] = intervals[name][1]
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["view", "modality"]).reset_index(drop=True)


def make_feature_stability(selected_features: pd.DataFrame, total_outer_fits: int) -> pd.DataFrame:
    rows = []
    grouped = selected_features.groupby(["modality", "view", "model_id", "feature_name", "source_modality"])
    for keys, frame in grouped:
        coefficients = frame["coefficient"].to_numpy(dtype=float)
        rows.append(
            {
                "modality": keys[0],
                "view": keys[1],
                "model_id": keys[2],
                "feature_name": keys[3],
                "source_modality": keys[4],
                "selected_outer_fits": len(frame),
                "total_outer_fits": total_outer_fits,
                "selection_frequency": len(frame) / total_outer_fits,
                "coefficient_mean_when_selected": float(coefficients.mean()),
                "coefficient_sd_when_selected": float(coefficients.std(ddof=1)) if len(coefficients) > 1 else 0.0,
                "positive_coefficient_fraction": float((coefficients > 0).mean()),
                "mean_selection_score": float(frame["selection_score"].mean()),
                "mean_effect_msa_minus_control": float(frame["effect_msa_minus_control"].mean()),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["view", "modality", "selection_frequency", "feature_name"],
        ascending=[True, True, False, True],
    ).reset_index(drop=True)


def make_auc_comparisons(
    donor_means: pd.DataFrame,
    step4_predictions: pd.DataFrame,
    bootstrap_iterations: int,
    seed: int,
) -> pd.DataFrame:
    records = []
    comparison_index = 0
    for (modality, view), molecular in donor_means.groupby(["modality", "view"]):
        for baseline in ("age_only", "all_covariates"):
            covariate = step4_predictions.loc[step4_predictions["model"] == baseline]
            merged = molecular.merge(
                covariate[["donor_id", "y_true", "probability_msa_mean"]],
                on=["donor_id", "y_true"],
                suffixes=("_molecular", "_baseline"),
                validate="one_to_one",
            )
            result = paired_auc_bootstrap(
                merged["y_true"].to_numpy(),
                merged["probability_msa_mean_molecular"].to_numpy(),
                merged["probability_msa_mean_baseline"].to_numpy(),
                bootstrap_iterations,
                seed + comparison_index * 1009,
            )
            records.append(
                {
                    "comparison_family": "molecular_vs_covariate_baseline",
                    "model_a": model_id(modality, view),
                    "model_b": baseline,
                    **result,
                    "bootstrap_iterations": bootstrap_iterations,
                }
            )
            comparison_index += 1
    for view in VIEWS:
        combined = donor_means.loc[(donor_means["modality"] == "combined") & (donor_means["view"] == view)]
        for single in ("5mC", "5hmC"):
            other = donor_means.loc[(donor_means["modality"] == single) & (donor_means["view"] == view)]
            merged = combined.merge(
                other[["donor_id", "y_true", "probability_msa_mean"]],
                on=["donor_id", "y_true"],
                suffixes=("_combined", "_single"),
                validate="one_to_one",
            )
            result = paired_auc_bootstrap(
                merged["y_true"].to_numpy(),
                merged["probability_msa_mean_combined"].to_numpy(),
                merged["probability_msa_mean_single"].to_numpy(),
                bootstrap_iterations,
                seed + comparison_index * 1009,
            )
            records.append(
                {
                    "comparison_family": "combined_vs_single_modality",
                    "model_a": model_id("combined", view),
                    "model_b": model_id(single, view),
                    **result,
                    "bootstrap_iterations": bootstrap_iterations,
                }
            )
            comparison_index += 1
    frame = pd.DataFrame(records)
    frame["fdr_bh_within_family"] = np.nan
    for family, indices in frame.groupby("comparison_family").groups.items():
        frame.loc[indices, "fdr_bh_within_family"] = benjamini_hochberg(
            frame.loc[indices, "bootstrap_two_sided_p"]
        )
    return frame


def plot_performance(aggregate: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), sharey=True)
    for axis, view in zip(axes, VIEWS):
        subset = aggregate.loc[aggregate["view"] == view].set_index("modality").loc[list(MODALITIES)]
        x = np.arange(len(MODALITIES))
        values = subset["roc_auc_donor_mean_oof"].to_numpy()
        lower = values - subset["roc_auc_bootstrap_low"].to_numpy()
        upper = subset["roc_auc_bootstrap_high"].to_numpy() - values
        axis.bar(x, values, color=[MODALITY_COLORS[item] for item in MODALITIES], alpha=0.9)
        axis.errorbar(x, values, yerr=np.vstack([lower, upper]), fmt="none", ecolor="black", capsize=4)
        axis.axhline(0.5, color="#777777", linestyle="--", linewidth=1)
        axis.set_xticks(x, MODALITIES)
        axis.set_ylim(0, 1.03)
        axis.set_title(VIEW_LABELS[view])
    axes[0].set_ylabel("Donor-mean out-of-fold ROC-AUC (95% CI)")
    fig.suptitle("Primary Clean-77 molecular benchmarks")
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_roc_curves(donor_means: pd.DataFrame, step4: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), sharex=True, sharey=True)
    baseline = step4.loc[step4["model"] == "all_covariates"]
    for axis, view in zip(axes, VIEWS):
        for modality in MODALITIES:
            frame = donor_means.loc[(donor_means["view"] == view) & (donor_means["modality"] == modality)]
            fpr, tpr, _ = roc_curve(frame["y_true"], frame["probability_msa_mean"])
            axis.plot(fpr, tpr, color=MODALITY_COLORS[modality], linewidth=2, label=f"{modality} (AUC={auc(fpr, tpr):.3f})")
        fpr, tpr, _ = roc_curve(baseline["y_true"], baseline["probability_msa_mean"])
        axis.plot(fpr, tpr, color="#555555", linestyle=":", linewidth=2, label=f"Step 4 all covariates ({auc(fpr, tpr):.3f})")
        axis.plot([0, 1], [0, 1], color="#888888", linestyle="--", linewidth=1)
        axis.set_title(VIEW_LABELS[view])
        axis.set_xlabel("False-positive rate")
        axis.legend(frameon=False, fontsize=8)
    axes[0].set_ylabel("True-positive rate")
    fig.suptitle("Primary Clean-77 molecular ROC curves")
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_incremental(comparisons: pd.DataFrame, output: Path) -> None:
    frame = comparisons.loc[
        (comparisons["comparison_family"] == "molecular_vs_covariate_baseline")
        & (comparisons["model_b"] == "all_covariates")
    ].copy()
    frame = frame.sort_values("roc_auc_difference_a_minus_b")
    values = frame["roc_auc_difference_a_minus_b"].to_numpy()
    lower = values - frame["bootstrap_low"].to_numpy()
    upper = frame["bootstrap_high"].to_numpy() - values
    fig, axis = plt.subplots(figsize=(9.5, 6.2))
    y = np.arange(len(frame))
    axis.errorbar(values, y, xerr=np.vstack([lower, upper]), fmt="o", color="#245a91", ecolor="#444444", capsize=3)
    axis.axvline(0, color="#777777", linestyle="--", linewidth=1)
    axis.set_yticks(y, frame["model_a"])
    axis.set_xlabel("ROC-AUC difference versus Step 4 all-covariate baseline (95% CI)")
    axis.set_title("Incremental molecular discrimination")
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_feature_stability(stability: pd.DataFrame, output: Path) -> None:
    frame = stability.loc[
        (stability["modality"] == "combined") & (stability["view"] == "molecular_plus_covariates")
    ].nlargest(20, ["selection_frequency", "mean_selection_score"])
    frame = frame.sort_values(["selection_frequency", "mean_selection_score"])
    fig, axis = plt.subplots(figsize=(10, 7.2))
    colors = [MODALITY_COLORS[source] for source in frame["source_modality"]]
    axis.barh(np.arange(len(frame)), frame["selection_frequency"], color=colors)
    axis.set_yticks(np.arange(len(frame)), frame["feature_name"])
    axis.set_xlim(0, 1)
    axis.set_xlabel("Selection frequency across outer fits")
    axis.set_title("Top repeated CpGs: combined molecular + covariate model")
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_report(output: Path, aggregate: pd.DataFrame, comparisons: pd.DataFrame, stability: pd.DataFrame, parameters: dict[str, object]) -> None:
    metrics = aggregate[["modality", "view", "roc_auc_donor_mean_oof", "roc_auc_bootstrap_low", "roc_auc_bootstrap_high", "pr_auc_donor_mean_oof", "balanced_accuracy_donor_mean_oof", "mcc_donor_mean_oof"]]
    benchmark = comparisons.loc[(comparisons["comparison_family"] == "molecular_vs_covariate_baseline") & (comparisons["model_b"] == "all_covariates")]
    top = stability.sort_values("selection_frequency", ascending=False).head(30)
    document = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>GSE143157 Step 5 report</title>
<style>body{{font-family:Segoe UI,Arial,sans-serif;max-width:1200px;margin:32px auto;padding:0 20px;color:#1d2630}}h1,h2{{color:#173b5f}}.good{{background:#eaf7ee;border-left:5px solid #238b45;padding:12px}}.warn{{background:#fff6df;border-left:5px solid #d79000;padding:12px}}table{{border-collapse:collapse;width:100%;font-size:13px}}th,td{{border:1px solid #d7dee5;padding:7px;text-align:left}}th{{background:#eef3f7}}img{{max-width:100%;margin:12px 0}}</style></head><body>
<h1>GSE143157 Step 5 — Primary Clean-77 molecular benchmark</h1>
<p class="good"><strong>Leakage guardrail:</strong> outer-test donors did not contribute to variance pools, missingness rules, imputation, scaling, supervised CpG selection, residualization, or tuning.</p>
<p class="warn"><strong>Interpretation:</strong> this report estimates prediction, not causality. CpGs require stability and biological annotation before being described as signatures.</p>
<h2>Design</h2><p>Three modalities × three model views; frozen 5×{parameters['outer_repeats']} outer CV and four inner folds. Combined models select the stated feature count separately from each modality.</p>
<h2>Performance</h2>{metrics.to_html(index=False,border=0,float_format=lambda x:f'{x:.4f}')}
<img src="performance_primary.png" alt="Performance"><img src="roc_primary.png" alt="ROC curves">
<h2>Incremental value versus Step 4 all-covariate baseline</h2>{benchmark.to_html(index=False,border=0,float_format=lambda x:f'{x:.4f}')}
<img src="incremental_auc_primary.png" alt="Incremental AUC">
<h2>Feature stability</h2>{top[["modality","view","feature_name","selected_outer_fits","selection_frequency","coefficient_mean_when_selected","positive_coefficient_fraction"]].to_html(index=False,border=0,float_format=lambda x:f'{x:.4f}')}
<img src="feature_stability_primary.png" alt="Feature stability">
<h2>Next decision</h2><p>Only after reviewing paired uncertainty and stability should one workflow be locked for Full-78 sensitivity analysis. Do not select a single lucky repeat or outer fold.</p>
</body></html>"""
    output.write_text(document, encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    started = time.time()
    step2b = args.step2b_dir.resolve()
    step3 = args.step3_dir.resolve()
    step4 = args.step4_dir.resolve()
    output = args.output.resolve()
    inputs = validate_inputs(step2b, step3, step4)

    if args.smoke:
        outer_repeats = 1
        candidate_pool = 300
        feature_counts = [10, 25]
        c_grid = [0.1, 1.0]
        class_weights = ["none"]
        bootstrap_iterations = 200
        analysis_mode = "SMOKE_TEST_NOT_FOR_RESEARCH"
    else:
        outer_repeats = 10
        candidate_pool = int(args.candidate_pool)
        feature_counts = parse_csv_numbers(args.feature_counts, int)
        c_grid = parse_csv_numbers(args.c_grid, float)
        class_weights = [item.strip() for item in args.class_weights.split(",") if item.strip()]
        bootstrap_iterations = int(args.bootstrap_iterations)
        analysis_mode = "FULL_RESEARCH_RUN"
    if not set(class_weights).issubset({"none", "balanced"}) or not class_weights:
        raise ValueError("Class weights must be none and/or balanced.")
    if candidate_pool < max(feature_counts):
        raise ValueError("Candidate pool must be at least the largest feature count.")
    if args.chunk_rows < 100 or not 0 <= args.missing_threshold < 1:
        raise ValueError("Invalid chunk size or missingness threshold.")
    if not args.smoke and bootstrap_iterations < 1000:
        raise ValueError("Full research runs require at least 1000 bootstrap iterations.")

    parameters = {
        "analysis_mode": analysis_mode,
        "outer_folds": OUTER_FOLDS,
        "outer_repeats": outer_repeats,
        "inner_folds": INNER_FOLDS,
        "candidate_pool": candidate_pool,
        "candidate_pool_method": "outer-training-only variance ranking; no diagnosis labels",
        "feature_counts": feature_counts,
        "combined_feature_count_rule": "feature_count separately from 5mC and 5hmC",
        "supervised_selection": "absolute Welch statistic inside each inner-training fold; refit on outer training only",
        "c_grid": c_grid,
        "class_weights": class_weights,
        "inner_selection_metric": "ROC-AUC",
        "missing_threshold": float(args.missing_threshold),
        "minimum_variance": float(args.minimum_variance),
        "chunk_rows": int(args.chunk_rows),
        "classification_threshold": 0.5,
        "bootstrap_iterations": bootstrap_iterations,
        "seed": int(args.seed),
    }
    paths = inputs["paths"]
    h5_stat = paths["h5"].stat()
    configuration = {
        "pipeline_version": VERSION,
        "parameters": parameters,
        "inputs": {
            "h5": str(paths["h5"].resolve()),
            "h5_size": h5_stat.st_size,
            "h5_mtime_ns": h5_stat.st_mtime_ns,
            "metadata_sha256": file_sha256(paths["metadata"]),
            "outer_folds_sha256": file_sha256(paths["outer"]),
            "inner_folds_sha256": file_sha256(paths["inner"]),
            "model_policy_sha256": file_sha256(paths["policy"]),
            "step4_predictions_sha256": file_sha256(paths["step4_predictions"]),
        },
    }
    output.mkdir(parents=True, exist_ok=True)
    if (output / "SUCCESS.txt").exists():
        raise ValueError(f"A successful Step 5 result already exists: {output}")
    config_path = output / "RUN_CONFIGURATION.json"
    if config_path.exists():
        previous = json.loads(config_path.read_text(encoding="utf-8"))
        if previous != configuration:
            raise ValueError("Incomplete output uses a different configuration. Rename it before rerunning.")
    else:
        config_path.write_text(json.dumps(configuration, indent=2, ensure_ascii=False), encoding="utf-8")
    checkpoint_dir = output / "_checkpoints"
    checkpoint_dir.mkdir(exist_ok=True)

    metadata = inputs["metadata"]
    labels = (metadata["condition"].astype(str) == "MSA").astype(int).to_numpy()
    completed_results = []
    with h5py.File(paths["h5"], "r") as h5:
        for repeat in range(1, outer_repeats + 1):
            for outer_fold in range(1, OUTER_FOLDS + 1):
                checkpoint = checkpoint_dir / f"repeat_{repeat:02d}_fold_{outer_fold}.pkl"
                if checkpoint.exists():
                    print(f"[repeat {repeat}/{outer_repeats}, fold {outer_fold}/5] loading checkpoint", flush=True)
                    with checkpoint.open("rb") as handle:
                        result = pickle.load(handle)
                else:
                    print(f"[repeat {repeat}/{outer_repeats}, fold {outer_fold}/5] starting", flush=True)
                    result = run_outer_fit(
                        h5,
                        metadata,
                        labels,
                        inputs["outer"],
                        inputs["inner"],
                        repeat,
                        outer_fold,
                        parameters,
                    )
                    partial = checkpoint.with_suffix(".pkl.partial")
                    with partial.open("wb") as handle:
                        pickle.dump(result, handle, protocol=pickle.HIGHEST_PROTOCOL)
                    os.replace(partial, checkpoint)
                    print(f"[repeat {repeat}/{outer_repeats}, fold {outer_fold}/5] checkpoint saved", flush=True)
                completed_results.append(result)

    result_names = ["predictions", "hyperparameter_choices", "tuning_summary", "outer_coefficients", "selected_features", "candidate_pool_summary"]
    frames = {
        name: pd.DataFrame([row for result in completed_results for row in result[name]])
        for name in result_names
    }
    predictions = frames["predictions"].sort_values(["view", "modality", "repeat", "donor_index"])
    repeat_metrics = make_repeat_metrics(predictions)
    donor_means = make_donor_means(predictions)
    aggregate = make_aggregate_metrics(donor_means, repeat_metrics, bootstrap_iterations, args.seed + 7000001)
    total_outer_fits = outer_repeats * OUTER_FOLDS
    stability = make_feature_stability(frames["selected_features"], total_outer_fits)
    comparisons = make_auc_comparisons(donor_means, inputs["step4_predictions"], bootstrap_iterations, args.seed + 9000001)

    final_frames = {
        **frames,
        "repeat_metrics": repeat_metrics,
        "donor_mean_predictions": donor_means,
        "aggregate_metrics": aggregate,
        "feature_stability": stability,
        "auc_comparisons": comparisons,
    }
    for name, frame in final_frames.items():
        frame.to_csv(output / f"{name}_primary_clean77.csv", index=False)

    plot_performance(aggregate, output / "performance_primary.png")
    plot_roc_curves(donor_means, inputs["step4_predictions"], output / "roc_primary.png")
    plot_incremental(comparisons, output / "incremental_auc_primary.png")
    plot_feature_stability(stability, output / "feature_stability_primary.png")
    write_report(output / "STEP5_REPORT.html", aggregate, comparisons, stability, parameters)

    summary = {
        "pipeline": {
            "name": "GSE143157 Step 5 primary molecular nested-CV benchmark",
            "version": VERSION,
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": time.time() - started,
            "scikit_learn_version": sklearn.__version__,
            "h5py_version": h5py.__version__,
        },
        "parameters": parameters,
        "inputs": configuration["inputs"],
        "cohort": {"name": "primary_clean77", "donors": PRIMARY_DONORS, "control": 37, "msa": 40, "excluded": ["MSA 30"]},
        "models": [model_id(modality, view) for modality, view in MODEL_SPECS],
        "aggregate_metrics": aggregate.to_dict(orient="records"),
        "auc_comparisons": comparisons.to_dict(orient="records"),
        "guardrails": [
            "Outer-test donors never influenced variance pools, missingness filters, imputation, scaling, supervised selection, residualization, or tuning.",
            "Outer-training-only unsupervised variance pools were fixed before inner tuning and used no diagnosis labels.",
            "All diagnosis-based CpG selection was fitted separately inside each inner-training fold.",
            "The final outer model refitted every transformation using outer-training donors only.",
            "The same frozen Step 3 donor folds were used for 5mC, 5hmC, and combined models.",
            "Smoke-test output is never a research result.",
        ],
    }
    (output / "step5_summary.json").write_text(
        json.dumps(json_clean(summary), indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8"
    )
    success_text = (
        "SUCCESS\n"
        f"Completed UTC: {summary['pipeline']['completed_utc']}\n"
        f"Analysis mode: {analysis_mode}\n"
        f"Outer predictions: {len(predictions)}\n"
        f"Outer repeats: {outer_repeats}\n"
        "Models: 3 modalities x 3 views\n"
        "All outer transformations and supervised selections were training-only\n"
    )
    (output / "SUCCESS.txt").write_text(success_text, encoding="utf-8")
    print("\nSTEP 5 SUCCESS")
    if args.smoke:
        print("SMOKE TEST ONLY - DO NOT USE THESE NUMBERS IN THE PAPER")
    print(f"Open report: {output / 'STEP5_REPORT.html'}")


def main() -> int:
    try:
        run(parse_args())
    except (ValueError, OSError, KeyError, json.JSONDecodeError, pd.errors.ParserError) as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
