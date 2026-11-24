#!/usr/bin/env python3
"""Core utilities for the GSE143157 Step 7 molecular sensitivity analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable
import warnings

import h5py
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    log_loss,
    matthews_corrcoef,
    roc_auc_score,
)


MODALITIES = ("5mC", "5hmC", "combined")
VIEWS = ("molecular_only", "molecular_plus_covariates", "residualized_molecular")
MODEL_SPECS = tuple((modality, view) for view in VIEWS for modality in MODALITIES)
NUMERIC_COVARIATES = ("age", "pmi_hours")
CATEGORICAL_COVARIATES = ("sex", "brain_bank")


@dataclass
class CandidatePool:
    modality: str
    row_indices: np.ndarray
    cpg_ids: np.ndarray
    outer_train_variances: np.ndarray
    values: np.ndarray
    eligible_features: int
    total_features: int


@dataclass
class PreparedFeatures:
    modality: str
    feature_names: np.ndarray
    scores: np.ndarray
    effects: np.ndarray
    missing_fractions: np.ndarray
    train: np.ndarray
    evaluate: np.ndarray


def decode_strings(values: Iterable[object]) -> list[str]:
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def _best_positions(scores: np.ndarray, row_indices: np.ndarray, count: int) -> np.ndarray:
    if count <= 0:
        return np.empty(0, dtype=np.int64)
    count = min(count, len(scores))
    if count < len(scores):
        provisional = np.argpartition(scores, len(scores) - count)[-count:]
    else:
        provisional = np.arange(len(scores), dtype=np.int64)
    order = np.lexsort((row_indices[provisional], -scores[provisional]))
    return provisional[order]


def _read_cpg_ids(dataset: h5py.Dataset, indices: np.ndarray) -> np.ndarray:
    order = np.argsort(indices)
    sorted_indices = indices[order]
    sorted_values = np.asarray(decode_strings(dataset[sorted_indices]), dtype=object)
    restore = np.argsort(order)
    return sorted_values[restore]


def stream_variance_candidates(
    h5: h5py.File,
    modality: str,
    training_columns: np.ndarray,
    pool_size: int,
    chunk_rows: int,
    missing_threshold: float,
    minimum_variance: float,
    cohort_columns: np.ndarray | None = None,
) -> CandidatePool:
    """Select an outer-training-only variance pool while retaining values for later application."""
    values_ds = h5[f"{modality}/values"]
    ids_ds = h5[f"{modality}/cpg_ids"]
    total_features, source_donors = values_ds.shape
    if cohort_columns is None:
        cohort_columns = np.arange(source_donors, dtype=np.int64)
    else:
        cohort_columns = np.asarray(cohort_columns, dtype=np.int64)
    if cohort_columns.ndim != 1 or not len(cohort_columns):
        raise ValueError("Cohort HDF5 columns must be a non-empty one-dimensional array.")
    if len(np.unique(cohort_columns)) != len(cohort_columns):
        raise ValueError("Cohort HDF5 columns contain duplicates.")
    if cohort_columns.min() < 0 or cohort_columns.max() >= source_donors:
        raise ValueError("A cohort HDF5 column is outside the source matrix.")
    if np.any(np.diff(cohort_columns) <= 0):
        raise ValueError("Cohort HDF5 columns must be strictly increasing for deterministic reads.")
    cohort_donors = len(cohort_columns)
    training_columns = np.asarray(sorted(map(int, training_columns)), dtype=np.int64)
    if training_columns.min() < 0 or training_columns.max() >= cohort_donors:
        raise ValueError("A training column is outside the cohort matrix.")
    minimum_finite = int(np.ceil((1.0 - missing_threshold) * len(training_columns)))

    best_scores = np.empty(0, dtype=np.float64)
    best_indices = np.empty(0, dtype=np.int64)
    best_values = np.empty((0, cohort_donors), dtype=np.float32)
    eligible_total = 0

    for start in range(0, total_features, chunk_rows):
        stop = min(start + chunk_rows, total_features)
        block_all = np.asarray(values_ds[start:stop, cohort_columns], dtype=np.float32)
        block_train = block_all[:, training_columns].astype(np.float64, copy=False)
        finite = np.isfinite(block_train)
        finite_counts = finite.sum(axis=1)
        eligible = finite_counts >= minimum_finite
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            medians = np.nanmedian(block_train, axis=1)
        imputed = np.where(finite, block_train, medians[:, None])
        variances = np.var(imputed, axis=1, ddof=1)
        eligible &= np.isfinite(variances) & (variances > minimum_variance)
        eligible_total += int(eligible.sum())
        local_positions = np.flatnonzero(eligible)
        if local_positions.size:
            local_rows = start + local_positions
            local_scores = variances[local_positions]
            keep_local = _best_positions(local_scores, local_rows, pool_size)
            local_positions = local_positions[keep_local]
            local_rows = local_rows[keep_local]
            local_scores = local_scores[keep_local]
            merged_scores = np.concatenate([best_scores, local_scores])
            merged_indices = np.concatenate([best_indices, local_rows])
            merged_values = np.vstack([best_values, block_all[local_positions]])
            keep_merged = _best_positions(merged_scores, merged_indices, pool_size)
            best_scores = merged_scores[keep_merged]
            best_indices = merged_indices[keep_merged]
            best_values = merged_values[keep_merged]
        if stop % 100000 < chunk_rows or stop == total_features:
            print(
                f"    [{modality}] streamed {stop:,}/{total_features:,} CpGs; "
                f"eligible={eligible_total:,}",
                flush=True,
            )

    if len(best_indices) < pool_size:
        raise ValueError(
            f"{modality}: only {len(best_indices)} eligible CpGs were available for a pool of {pool_size}."
        )
    final_order = np.lexsort((best_indices, -best_scores))
    best_indices = best_indices[final_order]
    best_scores = best_scores[final_order]
    best_values = best_values[final_order]
    cpg_ids = _read_cpg_ids(ids_ds, best_indices)
    return CandidatePool(
        modality=modality,
        row_indices=best_indices,
        cpg_ids=cpg_ids,
        outer_train_variances=best_scores,
        values=best_values,
        eligible_features=eligible_total,
        total_features=total_features,
    )


def fit_covariate_transform(
    metadata: pd.DataFrame,
    training_indices: np.ndarray,
    evaluation_indices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    train_frame = metadata.iloc[training_indices]
    eval_frame = metadata.iloc[evaluation_indices]
    train_columns: list[np.ndarray] = []
    eval_columns: list[np.ndarray] = []
    names: list[str] = []

    for column in NUMERIC_COVARIATES:
        train_values = pd.to_numeric(train_frame[column], errors="coerce").to_numpy(dtype=float)
        eval_values = pd.to_numeric(eval_frame[column], errors="coerce").to_numpy(dtype=float)
        median = float(np.nanmedian(train_values))
        train_values = np.where(np.isfinite(train_values), train_values, median)
        eval_values = np.where(np.isfinite(eval_values), eval_values, median)
        mean = float(train_values.mean())
        scale = float(train_values.std(ddof=0))
        if not np.isfinite(scale) or scale < 1e-12:
            scale = 1.0
        train_columns.append((train_values - mean) / scale)
        eval_columns.append((eval_values - mean) / scale)
        names.append(column)

    for column in CATEGORICAL_COVARIATES:
        train_values = train_frame[column].fillna("Missing").astype(str)
        eval_values = eval_frame[column].fillna("Missing").astype(str)
        mode = train_values.mode(dropna=False)
        replacement = str(mode.iloc[0]) if len(mode) else "Missing"
        train_values = train_values.replace("nan", replacement)
        eval_values = eval_values.replace("nan", replacement)
        levels = sorted(train_values.unique().tolist())
        for level in levels[1:]:
            train_columns.append((train_values.to_numpy() == level).astype(float))
            eval_columns.append((eval_values.to_numpy() == level).astype(float))
            names.append(f"{column}_{level}")

    return (
        np.column_stack(train_columns).astype(np.float64),
        np.column_stack(eval_columns).astype(np.float64),
        names,
    )


def _welch_scores(matrix: np.ndarray, labels: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    group0 = matrix[labels == 0]
    group1 = matrix[labels == 1]
    mean0 = group0.mean(axis=0)
    mean1 = group1.mean(axis=0)
    var0 = group0.var(axis=0, ddof=1)
    var1 = group1.var(axis=0, ddof=1)
    denominator = np.sqrt(var0 / len(group0) + var1 / len(group1) + 1e-12)
    effects = mean1 - mean0
    scores = np.abs(effects) / denominator
    scores[~np.isfinite(scores)] = 0.0
    return scores, effects


def prepare_candidates(
    pool: CandidatePool,
    training_indices: np.ndarray,
    evaluation_indices: np.ndarray,
    labels_training: np.ndarray,
    metadata: pd.DataFrame,
    maximum_features: int,
    residualize: bool,
    missing_threshold: float,
    minimum_variance: float,
) -> PreparedFeatures:
    train = pool.values[:, training_indices].T.astype(np.float64, copy=False)
    evaluate = pool.values[:, evaluation_indices].T.astype(np.float64, copy=False)
    finite_train = np.isfinite(train)
    missing_fraction = 1.0 - finite_train.mean(axis=0)
    valid = missing_fraction <= missing_threshold
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        medians = np.nanmedian(train, axis=0)
    valid &= np.isfinite(medians)
    train = np.where(finite_train, train, medians[None, :])
    evaluate = np.where(np.isfinite(evaluate), evaluate, medians[None, :])
    variances = train.var(axis=0, ddof=1)
    valid &= np.isfinite(variances) & (variances > minimum_variance)
    valid_positions = np.flatnonzero(valid)
    if len(valid_positions) < maximum_features:
        raise ValueError(
            f"{pool.modality}: only {len(valid_positions)} valid candidates remain; "
            f"need at least {maximum_features}."
        )
    train = train[:, valid_positions]
    evaluate = evaluate[:, valid_positions]
    feature_names = pool.cpg_ids[valid_positions]
    selected_missing = missing_fraction[valid_positions]

    if residualize:
        cov_train, cov_evaluate, _ = fit_covariate_transform(
            metadata, training_indices, evaluation_indices
        )
        design_train = np.column_stack([np.ones(len(training_indices)), cov_train])
        design_evaluate = np.column_stack([np.ones(len(evaluation_indices)), cov_evaluate])
        coefficients, *_ = np.linalg.lstsq(design_train, train, rcond=None)
        train = train - design_train @ coefficients
        evaluate = evaluate - design_evaluate @ coefficients

    scores, effects = _welch_scores(train, labels_training)
    tie_rows = pool.row_indices[valid_positions]
    order = np.lexsort((tie_rows, -scores))[:maximum_features]
    train = train[:, order]
    evaluate = evaluate[:, order]
    feature_names = feature_names[order]
    scores = scores[order]
    effects = effects[order]
    selected_missing = selected_missing[order]

    means = train.mean(axis=0)
    scales = train.std(axis=0, ddof=0)
    scales[~np.isfinite(scales) | (scales < 1e-12)] = 1.0
    train = (train - means) / scales
    evaluate = (evaluate - means) / scales
    return PreparedFeatures(
        modality=pool.modality,
        feature_names=feature_names,
        scores=scores,
        effects=effects,
        missing_fractions=selected_missing,
        train=train,
        evaluate=evaluate,
    )


def build_design(
    modality: str,
    view: str,
    feature_count: int,
    raw: dict[str, PreparedFeatures],
    residualized: dict[str, PreparedFeatures],
    covariate_train: np.ndarray,
    covariate_evaluate: np.ndarray,
    covariate_names: list[str],
) -> tuple[np.ndarray, np.ndarray, list[dict[str, object]]]:
    source = residualized if view == "residualized_molecular" else raw
    source_modalities = ("5mC", "5hmC") if modality == "combined" else (modality,)
    train_parts: list[np.ndarray] = []
    eval_parts: list[np.ndarray] = []
    records: list[dict[str, object]] = []
    for source_modality in source_modalities:
        prepared = source[source_modality]
        train_parts.append(prepared.train[:, :feature_count])
        eval_parts.append(prepared.evaluate[:, :feature_count])
        for rank in range(feature_count):
            records.append(
                {
                    "feature_name": f"{source_modality}::{prepared.feature_names[rank]}",
                    "feature_type": "molecular",
                    "source_modality": source_modality,
                    "within_modality_rank": rank + 1,
                    "selection_score": float(prepared.scores[rank]),
                    "effect_msa_minus_control": float(prepared.effects[rank]),
                    "training_missing_fraction": float(prepared.missing_fractions[rank]),
                }
            )

    if view == "molecular_plus_covariates":
        train_parts.append(covariate_train)
        eval_parts.append(covariate_evaluate)
        for name in covariate_names:
            records.append(
                {
                    "feature_name": f"covariate::{name}",
                    "feature_type": "covariate",
                    "source_modality": "covariate",
                    "within_modality_rank": np.nan,
                    "selection_score": np.nan,
                    "effect_msa_minus_control": np.nan,
                    "training_missing_fraction": np.nan,
                }
            )
    return np.column_stack(train_parts), np.column_stack(eval_parts), records


def fit_logistic(
    train: np.ndarray,
    labels: np.ndarray,
    c_value: float,
    class_weight_name: str,
    random_seed: int,
) -> LogisticRegression:
    class_weight = None if class_weight_name == "none" else "balanced"
    model = LogisticRegression(
        C=float(c_value),
        class_weight=class_weight,
        solver="liblinear",
        l1_ratio=0.0,
        max_iter=5000,
        random_state=int(random_seed),
    )
    model.fit(train, labels)
    return model


def binary_metrics(labels: np.ndarray, probabilities: np.ndarray) -> dict[str, float]:
    labels = np.asarray(labels, dtype=int)
    probabilities = np.clip(np.asarray(probabilities, dtype=float), 1e-8, 1 - 1e-8)
    predictions = (probabilities >= 0.5).astype(int)
    true_positive = int(((predictions == 1) & (labels == 1)).sum())
    false_negative = int(((predictions == 0) & (labels == 1)).sum())
    true_negative = int(((predictions == 0) & (labels == 0)).sum())
    false_positive = int(((predictions == 1) & (labels == 0)).sum())
    sensitivity = true_positive / max(1, true_positive + false_negative)
    specificity = true_negative / max(1, true_negative + false_positive)
    return {
        "roc_auc": float(roc_auc_score(labels, probabilities)),
        "pr_auc": float(average_precision_score(labels, probabilities)),
        "balanced_accuracy": float(balanced_accuracy_score(labels, predictions)),
        "sensitivity": float(sensitivity),
        "specificity": float(specificity),
        "mcc": float(matthews_corrcoef(labels, predictions)),
        "accuracy": float(accuracy_score(labels, predictions)),
        "brier_score": float(brier_score_loss(labels, probabilities)),
        "log_loss": float(log_loss(labels, probabilities, labels=[0, 1])),
    }


def stratified_bootstrap_metrics(
    labels: np.ndarray,
    probabilities: np.ndarray,
    iterations: int,
    seed: int,
) -> dict[str, tuple[float, float]]:
    labels = np.asarray(labels, dtype=int)
    probabilities = np.asarray(probabilities, dtype=float)
    group0 = np.flatnonzero(labels == 0)
    group1 = np.flatnonzero(labels == 1)
    rng = np.random.default_rng(seed)
    samples: dict[str, list[float]] = {name: [] for name in binary_metrics(labels, probabilities)}
    for _ in range(iterations):
        indices = np.concatenate(
            [rng.choice(group0, len(group0), replace=True), rng.choice(group1, len(group1), replace=True)]
        )
        result = binary_metrics(labels[indices], probabilities[indices])
        for name, value in result.items():
            samples[name].append(value)
    return {
        name: (float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975)))
        for name, values in samples.items()
    }


def paired_auc_bootstrap(
    labels: np.ndarray,
    probabilities_a: np.ndarray,
    probabilities_b: np.ndarray,
    iterations: int,
    seed: int,
) -> dict[str, float]:
    labels = np.asarray(labels, dtype=int)
    probabilities_a = np.asarray(probabilities_a, dtype=float)
    probabilities_b = np.asarray(probabilities_b, dtype=float)
    observed = float(roc_auc_score(labels, probabilities_a) - roc_auc_score(labels, probabilities_b))
    group0 = np.flatnonzero(labels == 0)
    group1 = np.flatnonzero(labels == 1)
    rng = np.random.default_rng(seed)
    differences = np.empty(iterations, dtype=float)
    for index in range(iterations):
        selected = np.concatenate(
            [rng.choice(group0, len(group0), replace=True), rng.choice(group1, len(group1), replace=True)]
        )
        differences[index] = roc_auc_score(labels[selected], probabilities_a[selected]) - roc_auc_score(
            labels[selected], probabilities_b[selected]
        )
    lower = float(np.quantile(differences, 0.025))
    upper = float(np.quantile(differences, 0.975))
    left = (float((differences <= 0).sum()) + 1.0) / (iterations + 1.0)
    right = (float((differences >= 0).sum()) + 1.0) / (iterations + 1.0)
    return {
        "roc_auc_difference_a_minus_b": observed,
        "bootstrap_low": lower,
        "bootstrap_high": upper,
        "bootstrap_two_sided_p": min(1.0, 2.0 * min(left, right)),
    }


def benjamini_hochberg(values: Iterable[float]) -> np.ndarray:
    p_values = np.asarray(list(values), dtype=float)
    order = np.argsort(p_values)
    adjusted = np.empty(len(p_values), dtype=float)
    running = 1.0
    for rank_index in range(len(p_values) - 1, -1, -1):
        original_position = order[rank_index]
        rank = rank_index + 1
        running = min(running, p_values[original_position] * len(p_values) / rank)
        adjusted[original_position] = min(1.0, running)
    return adjusted
