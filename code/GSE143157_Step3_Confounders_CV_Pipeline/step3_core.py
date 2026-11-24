"""Statistics and deterministic nested-CV helpers for GSE143157 Step 3."""

from __future__ import annotations

import json
import math
import re
from collections.abc import Iterable

import numpy as np
import pandas as pd


CONTINUOUS_COVARIATES = ["age", "pmi_hours"]
CATEGORICAL_COVARIATES = ["sex", "brain_bank"]
REQUIRED_COLUMNS = [
    "donor_index",
    "original_hdf5_column",
    "donor_id",
    "condition",
    "sex",
    "age",
    "pmi_hours",
    "brain_bank",
]


def require_columns(frame: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def validate_metadata(
    frame: pd.DataFrame,
    expected_donors: int,
    expected_conditions: dict[str, int],
    label: str,
) -> pd.DataFrame:
    require_columns(frame, REQUIRED_COLUMNS, label)
    result = frame.copy()
    if len(result) != expected_donors:
        raise ValueError(f"{label}: expected {expected_donors} donors, found {len(result)}.")
    if result["donor_id"].duplicated().any():
        raise ValueError(f"{label}: duplicate donor IDs found.")
    if result["condition"].value_counts().to_dict() != expected_conditions:
        raise ValueError(
            f"{label}: expected condition counts {expected_conditions}, "
            f"found {result['condition'].value_counts().to_dict()}."
        )
    indices = pd.to_numeric(result["donor_index"], errors="raise").astype(int).to_numpy()
    if not np.array_equal(indices, np.arange(expected_donors)):
        raise ValueError(f"{label}: donor_index is not sequential 0..{expected_donors - 1}.")
    for column in CONTINUOUS_COVARIATES:
        result[column] = pd.to_numeric(result[column], errors="coerce")
        if result[column].notna().sum() < expected_donors // 2:
            raise ValueError(f"{label}: too many missing/non-numeric values in {column}.")
    for column in CATEGORICAL_COVARIATES:
        result[column] = result[column].fillna("Unknown").astype(str)
    return result


def hedges_g(msa: np.ndarray, control: np.ndarray) -> float:
    msa = msa[np.isfinite(msa)]
    control = control[np.isfinite(control)]
    n1, n0 = len(msa), len(control)
    if n1 < 2 or n0 < 2:
        return float("nan")
    pooled_var = ((n1 - 1) * np.var(msa, ddof=1) + (n0 - 1) * np.var(control, ddof=1)) / (
        n1 + n0 - 2
    )
    if pooled_var <= 0:
        return 0.0
    d = (float(np.mean(msa)) - float(np.mean(control))) / math.sqrt(pooled_var)
    correction = 1.0 - 3.0 / (4.0 * (n1 + n0) - 9.0)
    return correction * d


def continuous_permutation_p(
    values: np.ndarray,
    labels: np.ndarray,
    permutations: int,
    seed: int,
) -> tuple[float, float]:
    finite = np.isfinite(values)
    values = values[finite]
    labels = labels[finite]
    msa = labels == "MSA"
    observed = float(np.mean(values[msa]) - np.mean(values[~msa]))
    rng = np.random.default_rng(seed)
    extreme = 0
    for _ in range(permutations):
        permuted = rng.permutation(labels)
        group = permuted == "MSA"
        difference = float(np.mean(values[group]) - np.mean(values[~group]))
        if abs(difference) >= abs(observed) - 1e-15:
            extreme += 1
    return observed, (extreme + 1.0) / (permutations + 1.0)


def contingency_chi_square(labels: np.ndarray, values: np.ndarray) -> tuple[float, np.ndarray, list[str]]:
    categories = sorted({str(value) for value in values})
    table = np.zeros((2, len(categories)), dtype=float)
    for row, condition in enumerate(["Control", "MSA"]):
        for column, category in enumerate(categories):
            table[row, column] = np.sum((labels == condition) & (values == category))
    row_totals = table.sum(axis=1, keepdims=True)
    col_totals = table.sum(axis=0, keepdims=True)
    expected = row_totals @ col_totals / table.sum()
    valid = expected > 0
    chi_square = float(np.sum(((table - expected) ** 2)[valid] / expected[valid]))
    return chi_square, table, categories


def categorical_permutation_p(
    values: np.ndarray,
    labels: np.ndarray,
    permutations: int,
    seed: int,
) -> tuple[float, float, np.ndarray, list[str]]:
    observed, table, categories = contingency_chi_square(labels, values)
    rng = np.random.default_rng(seed)
    extreme = 0
    for _ in range(permutations):
        permuted = rng.permutation(labels)
        statistic, _, _ = contingency_chi_square(permuted, values)
        if statistic >= observed - 1e-15:
            extreme += 1
    return observed, (extreme + 1.0) / (permutations + 1.0), table, categories


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    output = np.empty_like(adjusted)
    output[order] = adjusted
    return output


def effect_magnitude(effect_name: str, effect: float) -> str:
    value = abs(effect)
    if effect_name == "Hedges g":
        return "large" if value >= 0.8 else "moderate" if value >= 0.5 else "small" if value >= 0.2 else "negligible"
    return "large" if value >= 0.5 else "moderate" if value >= 0.3 else "small" if value >= 0.1 else "negligible"


def covariate_audit(metadata: pd.DataFrame, permutations: int, seed: int) -> pd.DataFrame:
    labels = metadata["condition"].astype(str).to_numpy()
    rows: list[dict[str, object]] = []
    for offset, covariate in enumerate(CONTINUOUS_COVARIATES):
        values = pd.to_numeric(metadata[covariate], errors="coerce").to_numpy(dtype=float)
        msa = values[labels == "MSA"]
        control = values[labels == "Control"]
        observed, p_value = continuous_permutation_p(
            values, labels, permutations, seed + 1009 * (offset + 1)
        )
        effect = hedges_g(msa, control)
        rows.append(
            {
                "covariate": covariate,
                "covariate_type": "continuous",
                "n_control_nonmissing": int(np.isfinite(control).sum()),
                "n_msa_nonmissing": int(np.isfinite(msa).sum()),
                "control_summary": f"mean={np.nanmean(control):.3f}; sd={np.nanstd(control, ddof=1):.3f}; median={np.nanmedian(control):.3f}",
                "msa_summary": f"mean={np.nanmean(msa):.3f}; sd={np.nanstd(msa, ddof=1):.3f}; median={np.nanmedian(msa):.3f}",
                "difference_msa_minus_control": observed,
                "effect_name": "Hedges g",
                "effect_size": effect,
                "permutation_p_two_sided": p_value,
                "category_counts_control": "",
                "category_counts_msa": "",
            }
        )
    for offset, covariate in enumerate(CATEGORICAL_COVARIATES):
        values = metadata[covariate].fillna("Unknown").astype(str).to_numpy()
        statistic, p_value, table, categories = categorical_permutation_p(
            values, labels, permutations, seed + 5003 * (offset + 1)
        )
        denominator = len(values) * min(table.shape[0] - 1, table.shape[1] - 1)
        effect = 0.0 if denominator <= 0 else math.sqrt(statistic / denominator)
        control_counts = {categories[i]: int(table[0, i]) for i in range(len(categories))}
        msa_counts = {categories[i]: int(table[1, i]) for i in range(len(categories))}
        control_props = table[0] / max(1.0, table[0].sum())
        msa_props = table[1] / max(1.0, table[1].sum())
        rows.append(
            {
                "covariate": covariate,
                "covariate_type": "categorical",
                "n_control_nonmissing": int(table[0].sum()),
                "n_msa_nonmissing": int(table[1].sum()),
                "control_summary": json.dumps(control_counts, ensure_ascii=False, sort_keys=True),
                "msa_summary": json.dumps(msa_counts, ensure_ascii=False, sort_keys=True),
                "difference_msa_minus_control": float(np.max(np.abs(msa_props - control_props))),
                "effect_name": "Cramer's V",
                "effect_size": effect,
                "permutation_p_two_sided": p_value,
                "category_counts_control": json.dumps(control_counts, ensure_ascii=False, sort_keys=True),
                "category_counts_msa": json.dumps(msa_counts, ensure_ascii=False, sort_keys=True),
            }
        )
    result = pd.DataFrame(rows)
    result["fdr_bh"] = benjamini_hochberg(result["permutation_p_two_sided"].to_numpy())
    result["effect_magnitude"] = [
        effect_magnitude(name, value)
        for name, value in zip(result["effect_name"], result["effect_size"])
    ]
    result["model_policy"] = "Prespecified adjustment plus covariate-only baseline"
    return result


def _stratified_assignment(labels: np.ndarray, folds: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    assignment = np.full(len(labels), -1, dtype=int)
    for label in sorted(set(labels)):
        indices = np.flatnonzero(labels == label)
        rng.shuffle(indices)
        fold_order = rng.permutation(folds)
        for position, index in enumerate(indices):
            assignment[index] = int(fold_order[position % folds])
    if np.any(assignment < 0):
        raise ValueError("Fold assignment failed.")
    return assignment


def _category_proportions(values: pd.Series, categories: list[str]) -> np.ndarray:
    strings = values.fillna("Unknown").astype(str)
    return np.asarray([(strings == category).mean() for category in categories], dtype=float)


def fold_balance_score(metadata: pd.DataFrame, assignment: np.ndarray, folds: int) -> float:
    deviations: list[float] = []
    for covariate in CONTINUOUS_COVARIATES:
        values = pd.to_numeric(metadata[covariate], errors="coerce")
        global_mean = float(values.mean())
        global_sd = float(values.std(ddof=1))
        scale = global_sd if np.isfinite(global_sd) and global_sd > 0 else 1.0
        global_missing = float(values.isna().mean())
        for fold in range(folds):
            subset = values.iloc[np.flatnonzero(assignment == fold)]
            deviations.append(abs(float(subset.mean()) - global_mean) / scale)
            deviations.append(abs(float(subset.isna().mean()) - global_missing))
    for covariate in CATEGORICAL_COVARIATES:
        categories = sorted(metadata[covariate].fillna("Unknown").astype(str).unique())
        global_props = _category_proportions(metadata[covariate], categories)
        for fold in range(folds):
            subset = metadata.iloc[np.flatnonzero(assignment == fold)][covariate]
            deviations.extend(np.abs(_category_proportions(subset, categories) - global_props))
    values = np.asarray(deviations, dtype=float)
    return float(values.mean() + 0.75 * values.max())


def canonical_signature(assignment: np.ndarray, donor_ids: np.ndarray, folds: int) -> tuple[str, ...]:
    parts = [
        ",".join(sorted(donor_ids[assignment == fold].astype(str).tolist()))
        for fold in range(folds)
    ]
    return tuple(sorted(parts))


def optimized_stratified_assignment(
    metadata: pd.DataFrame,
    folds: int,
    seed: int,
    candidates: int,
    forbidden_signatures: set[tuple[str, ...]] | None = None,
) -> tuple[np.ndarray, float, int, tuple[str, ...]]:
    labels = metadata["condition"].astype(str).to_numpy()
    donor_ids = metadata["donor_id"].astype(str).to_numpy()
    forbidden_signatures = forbidden_signatures or set()
    best: tuple[np.ndarray, float, int, tuple[str, ...]] | None = None
    for candidate in range(candidates):
        candidate_seed = seed + candidate * 104729
        assignment = _stratified_assignment(labels, folds, candidate_seed)
        signature = canonical_signature(assignment, donor_ids, folds)
        if signature in forbidden_signatures:
            continue
        score = fold_balance_score(metadata, assignment, folds)
        if best is None or score < best[1]:
            best = (assignment, score, candidate_seed, signature)
    if best is None:
        raise ValueError("Could not find a distinct balanced fold assignment.")
    return best


def _safe_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def summarize_outer_fold(
    metadata: pd.DataFrame,
    assignment: np.ndarray,
    repeat: int,
    fold: int,
    partition_score: float,
) -> dict[str, object]:
    test = metadata.iloc[np.flatnonzero(assignment == fold)]
    train = metadata.iloc[np.flatnonzero(assignment != fold)]
    row: dict[str, object] = {
        "repeat": repeat,
        "outer_fold": fold + 1,
        "partition_balance_score": partition_score,
        "test_donors": len(test),
        "train_donors": len(train),
        "test_control": int((test["condition"] == "Control").sum()),
        "test_msa": int((test["condition"] == "MSA").sum()),
    }
    for covariate in CONTINUOUS_COVARIATES:
        train_values = pd.to_numeric(train[covariate], errors="coerce")
        test_values = pd.to_numeric(test[covariate], errors="coerce")
        row[f"{covariate}_train_mean"] = float(train_values.mean())
        row[f"{covariate}_test_mean"] = float(test_values.mean())
        row[f"{covariate}_test_minus_train"] = float(test_values.mean() - train_values.mean())
        row[f"{covariate}_test_missing_fraction"] = float(test_values.isna().mean())
    for covariate in CATEGORICAL_COVARIATES:
        categories = sorted(metadata[covariate].fillna("Unknown").astype(str).unique())
        for category in categories:
            name = _safe_name(category)
            row[f"{covariate}_{name}_train_fraction"] = float(
                (train[covariate].fillna("Unknown").astype(str) == category).mean()
            )
            row[f"{covariate}_{name}_test_fraction"] = float(
                (test[covariate].fillna("Unknown").astype(str) == category).mean()
            )
    return row


def build_nested_cv(
    metadata: pd.DataFrame,
    repeats: int,
    outer_folds: int,
    inner_folds: int,
    candidates: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    outer_rows: list[pd.DataFrame] = []
    inner_rows: list[pd.DataFrame] = []
    balance_rows: list[dict[str, object]] = []
    used_signatures: set[tuple[str, ...]] = set()
    chosen_seeds: list[int] = []
    for repeat_index in range(repeats):
        print(f"  building outer repeat {repeat_index + 1}/{repeats}...", flush=True)
        assignment, score, chosen_seed, signature = optimized_stratified_assignment(
            metadata,
            outer_folds,
            seed + repeat_index * 1_000_003,
            candidates,
            used_signatures,
        )
        used_signatures.add(signature)
        chosen_seeds.append(chosen_seed)
        outer = metadata.copy()
        outer.insert(0, "repeat", repeat_index + 1)
        outer.insert(1, "outer_test_fold", assignment + 1)
        outer_rows.append(outer)

        for fold in range(outer_folds):
            balance_rows.append(
                summarize_outer_fold(metadata, assignment, repeat_index + 1, fold, score)
            )
            train_positions = np.flatnonzero(assignment != fold)
            train = metadata.iloc[train_positions].reset_index(drop=True)
            inner_assignment, inner_score, inner_seed, _ = optimized_stratified_assignment(
                train,
                inner_folds,
                seed + repeat_index * 1_000_003 + (fold + 1) * 10_007,
                max(50, candidates // 2),
            )
            inner = train.copy()
            inner.insert(0, "outer_repeat", repeat_index + 1)
            inner.insert(1, "outer_fold", fold + 1)
            inner.insert(2, "inner_validation_fold", inner_assignment + 1)
            inner.insert(3, "inner_partition_balance_score", inner_score)
            inner.insert(4, "inner_seed", inner_seed)
            inner_rows.append(inner)

    outer_frame = pd.concat(outer_rows, ignore_index=True)
    inner_frame = pd.concat(inner_rows, ignore_index=True)
    balance_frame = pd.DataFrame(balance_rows)
    details = {
        "outer_partition_seeds": chosen_seeds,
        "unique_outer_partitions": len(used_signatures),
        "mean_outer_partition_balance_score": float(balance_frame["partition_balance_score"].mean()),
        "max_abs_outer_age_test_minus_train": float(balance_frame["age_test_minus_train"].abs().max()),
        "max_abs_outer_pmi_test_minus_train": float(balance_frame["pmi_hours_test_minus_train"].abs().max()),
    }
    return outer_frame, inner_frame, balance_frame, details
