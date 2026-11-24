"""Leakage-safe covariate-only nested-CV helpers for GSE143157 Step 4."""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    log_loss,
    matthews_corrcoef,
    roc_auc_score,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


@dataclass(frozen=True)
class ModelSpec:
    name: str
    continuous: tuple[str, ...]
    categorical: tuple[str, ...]
    description: str


MODEL_SPECS = {
    "age_only": ModelSpec(
        "age_only",
        ("age",),
        (),
        "Diagnostic model using age alone",
    ),
    "no_age_covariates": ModelSpec(
        "no_age_covariates",
        ("pmi_hours",),
        ("sex", "brain_bank"),
        "Sex, brain bank, and PMI without age",
    ),
    "all_covariates": ModelSpec(
        "all_covariates",
        ("age", "pmi_hours"),
        ("sex", "brain_bank"),
        "Prespecified age, sex, brain bank, and PMI benchmark",
    ),
}

METRIC_NAMES = [
    "roc_auc",
    "pr_auc",
    "balanced_accuracy",
    "sensitivity",
    "specificity",
    "mcc",
    "accuracy",
    "brier_score",
    "log_loss",
]


def encode_condition(values: pd.Series) -> np.ndarray:
    mapping = {"Control": 0, "MSA": 1}
    if not set(values.astype(str)).issubset(mapping):
        raise ValueError("Condition values must be Control or MSA.")
    return values.astype(str).map(mapping).to_numpy(dtype=int)


def build_pipeline(spec: ModelSpec, c_value: float, class_weight: str | None, seed: int) -> Pipeline:
    transformers = []
    if spec.continuous:
        numeric = Pipeline(
            [
                ("impute", SimpleImputer(strategy="median")),
                ("scale", StandardScaler()),
            ]
        )
        transformers.append(("continuous", numeric, list(spec.continuous)))
    if spec.categorical:
        categorical = Pipeline(
            [
                ("impute", SimpleImputer(strategy="most_frequent")),
                (
                    "onehot",
                    OneHotEncoder(
                        drop="if_binary",
                        handle_unknown="ignore",
                        sparse_output=False,
                    ),
                ),
            ]
        )
        transformers.append(("categorical", categorical, list(spec.categorical)))
    preprocess = ColumnTransformer(
        transformers,
        remainder="drop",
        sparse_threshold=0.0,
        verbose_feature_names_out=False,
    )
    classifier = LogisticRegression(
        C=float(c_value),
        class_weight=class_weight,
        solver="liblinear",
        max_iter=5000,
        random_state=int(seed),
    )
    return Pipeline([("preprocess", preprocess), ("classifier", classifier)])


def positive_probability(model: Pipeline, frame: pd.DataFrame) -> np.ndarray:
    classifier = model.named_steps["classifier"]
    class_position = int(np.flatnonzero(classifier.classes_ == 1)[0])
    return model.predict_proba(frame)[:, class_position]


def metric_values(y_true: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    y_true = np.asarray(y_true, dtype=int)
    probability = np.asarray(probability, dtype=float)
    prediction = (probability >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, prediction, labels=[0, 1]).ravel()
    sensitivity = tp / (tp + fn) if tp + fn else float("nan")
    specificity = tn / (tn + fp) if tn + fp else float("nan")
    clipped = np.clip(probability, 1e-7, 1 - 1e-7)
    return {
        "roc_auc": float(roc_auc_score(y_true, probability)),
        "pr_auc": float(average_precision_score(y_true, probability)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, prediction)),
        "sensitivity": float(sensitivity),
        "specificity": float(specificity),
        "mcc": float(matthews_corrcoef(y_true, prediction)),
        "accuracy": float(accuracy_score(y_true, prediction)),
        "brier_score": float(brier_score_loss(y_true, probability)),
        "log_loss": float(log_loss(y_true, clipped, labels=[0, 1])),
    }


def select_hyperparameters(
    outer_train: pd.DataFrame,
    inner_assignments: pd.DataFrame,
    spec: ModelSpec,
    c_grid: list[float],
    class_weights: list[str | None],
    seed: int,
) -> tuple[dict[str, object], pd.DataFrame]:
    if set(outer_train["donor_id"]) != set(inner_assignments["donor_id"]):
        raise ValueError("Inner-fold donors do not exactly match outer-training donors.")
    inner_map = inner_assignments.set_index("donor_id")["inner_validation_fold"]
    frame = outer_train.copy()
    frame["inner_validation_fold"] = frame["donor_id"].map(inner_map)
    if frame["inner_validation_fold"].isna().any():
        raise ValueError("Missing inner-fold assignment.")
    y = encode_condition(frame["condition"])
    rows: list[dict[str, object]] = []
    for combination_index, (c_value, class_weight) in enumerate(
        itertools.product(c_grid, class_weights)
    ):
        fold_aucs = []
        for inner_fold in sorted(frame["inner_validation_fold"].unique()):
            validation_mask = frame["inner_validation_fold"].eq(inner_fold).to_numpy()
            train_mask = ~validation_mask
            model = build_pipeline(
                spec,
                c_value,
                class_weight,
                seed + combination_index * 101 + int(inner_fold),
            )
            model.fit(frame.loc[train_mask], y[train_mask])
            probability = positive_probability(model, frame.loc[validation_mask])
            fold_aucs.append(float(roc_auc_score(y[validation_mask], probability)))
        rows.append(
            {
                "model": spec.name,
                "C": float(c_value),
                "class_weight": "none" if class_weight is None else str(class_weight),
                "mean_inner_roc_auc": float(np.mean(fold_aucs)),
                "sd_inner_roc_auc": float(np.std(fold_aucs, ddof=1)),
                "min_inner_roc_auc": float(np.min(fold_aucs)),
                "inner_fold_roc_aucs": ";".join(f"{value:.8f}" for value in fold_aucs),
            }
        )
    tuning = pd.DataFrame(rows)
    tuning["class_weight_preference"] = tuning["class_weight"].ne("none").astype(int)
    tuning = tuning.sort_values(
        [
            "mean_inner_roc_auc",
            "sd_inner_roc_auc",
            "C",
            "class_weight_preference",
        ],
        ascending=[False, True, True, True],
        kind="stable",
    ).reset_index(drop=True)
    tuning["selected"] = False
    tuning.loc[0, "selected"] = True
    best = tuning.iloc[0]
    return {
        "C": float(best["C"]),
        "class_weight": None if best["class_weight"] == "none" else str(best["class_weight"]),
        "mean_inner_roc_auc": float(best["mean_inner_roc_auc"]),
        "sd_inner_roc_auc": float(best["sd_inner_roc_auc"]),
    }, tuning.drop(columns=["class_weight_preference"])


def extract_coefficients(
    fitted_model: Pipeline,
    model_name: str,
    repeat: int,
    outer_fold: int,
    c_value: float,
    class_weight: str | None,
) -> pd.DataFrame:
    names = fitted_model.named_steps["preprocess"].get_feature_names_out()
    coefficients = fitted_model.named_steps["classifier"].coef_[0]
    return pd.DataFrame(
        {
            "model": model_name,
            "repeat": repeat,
            "outer_fold": outer_fold,
            "feature": names,
            "coefficient": coefficients,
            "odds_ratio_per_transformed_unit": np.exp(coefficients),
            "selected_C": c_value,
            "selected_class_weight": "none" if class_weight is None else class_weight,
        }
    )


def repeat_metrics(predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (model, repeat), group in predictions.groupby(["model", "repeat"], sort=False):
        metrics = metric_values(group["y_true"].to_numpy(), group["probability_msa"].to_numpy())
        rows.append({"model": model, "repeat": int(repeat), "donors": len(group), **metrics})
    return pd.DataFrame(rows)


def donor_mean_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        predictions.groupby(["model", "donor_id", "condition", "y_true"], as_index=False)
        .agg(
            probability_msa_mean=("probability_msa", "mean"),
            probability_msa_sd=("probability_msa", "std"),
            repeats=("repeat", "nunique"),
        )
    )
    grouped["predicted_label"] = np.where(
        grouped["probability_msa_mean"] >= 0.5, "MSA", "Control"
    )
    return grouped


def stratified_bootstrap_metrics(
    y_true: np.ndarray,
    probability: np.ndarray,
    iterations: int,
    seed: int,
) -> dict[str, tuple[float, float]]:
    y_true = np.asarray(y_true, dtype=int)
    probability = np.asarray(probability, dtype=float)
    control = np.flatnonzero(y_true == 0)
    msa = np.flatnonzero(y_true == 1)
    rng = np.random.default_rng(seed)
    values = {metric: [] for metric in METRIC_NAMES}
    for _ in range(iterations):
        sampled = np.concatenate(
            [
                rng.choice(control, size=len(control), replace=True),
                rng.choice(msa, size=len(msa), replace=True),
            ]
        )
        metrics = metric_values(y_true[sampled], probability[sampled])
        for metric in METRIC_NAMES:
            values[metric].append(metrics[metric])
    return {
        metric: (
            float(np.nanpercentile(metric_values_, 2.5)),
            float(np.nanpercentile(metric_values_, 97.5)),
        )
        for metric, metric_values_ in values.items()
    }


def aggregate_metrics(
    repeats: pd.DataFrame,
    donor_means: pd.DataFrame,
    bootstrap_iterations: int,
    seed: int,
) -> pd.DataFrame:
    rows = []
    for model_index, (model, donor_group) in enumerate(donor_means.groupby("model", sort=False)):
        repeat_group = repeats.loc[repeats["model"] == model]
        aggregate = metric_values(
            donor_group["y_true"].to_numpy(), donor_group["probability_msa_mean"].to_numpy()
        )
        intervals = stratified_bootstrap_metrics(
            donor_group["y_true"].to_numpy(),
            donor_group["probability_msa_mean"].to_numpy(),
            bootstrap_iterations,
            seed + model_index * 1009,
        )
        row: dict[str, object] = {
            "model": model,
            "donors": len(donor_group),
            "outer_repeats": int(repeat_group["repeat"].nunique()),
        }
        for metric in METRIC_NAMES:
            row[f"{metric}_repeat_mean"] = float(repeat_group[metric].mean())
            row[f"{metric}_repeat_sd"] = float(repeat_group[metric].std(ddof=1))
            row[f"{metric}_donor_mean_oof"] = aggregate[metric]
            row[f"{metric}_bootstrap_low"] = intervals[metric][0]
            row[f"{metric}_bootstrap_high"] = intervals[metric][1]
        rows.append(row)
    return pd.DataFrame(rows)


def pairwise_auc_comparisons(
    donor_means: pd.DataFrame,
    iterations: int,
    seed: int,
) -> pd.DataFrame:
    wide = donor_means.pivot(index=["donor_id", "y_true"], columns="model", values="probability_msa_mean").reset_index()
    y = wide["y_true"].to_numpy(dtype=int)
    control = np.flatnonzero(y == 0)
    msa = np.flatnonzero(y == 1)
    rng = np.random.default_rng(seed)
    rows = []
    pairs = [
        ("all_covariates", "age_only"),
        ("all_covariates", "no_age_covariates"),
        ("age_only", "no_age_covariates"),
    ]
    for first, second in pairs:
        first_probability = wide[first].to_numpy(dtype=float)
        second_probability = wide[second].to_numpy(dtype=float)
        observed = float(
            roc_auc_score(y, first_probability) - roc_auc_score(y, second_probability)
        )
        differences = []
        for _ in range(iterations):
            sampled = np.concatenate(
                [
                    rng.choice(control, size=len(control), replace=True),
                    rng.choice(msa, size=len(msa), replace=True),
                ]
            )
            differences.append(
                roc_auc_score(y[sampled], first_probability[sampled])
                - roc_auc_score(y[sampled], second_probability[sampled])
            )
        differences = np.asarray(differences)
        two_sided = min(
            1.0,
            2.0
            * min(
                (np.sum(differences <= 0) + 1) / (iterations + 1),
                (np.sum(differences >= 0) + 1) / (iterations + 1),
            ),
        )
        rows.append(
            {
                "model_a": first,
                "model_b": second,
                "roc_auc_difference_a_minus_b": observed,
                "bootstrap_low": float(np.percentile(differences, 2.5)),
                "bootstrap_high": float(np.percentile(differences, 97.5)),
                "bootstrap_two_sided_p": float(two_sided),
                "bootstrap_iterations": iterations,
            }
        )
    return pd.DataFrame(rows)


def coefficient_stability(coefficients: pd.DataFrame) -> pd.DataFrame:
    return (
        coefficients.groupby(["model", "feature"], as_index=False)
        .agg(
            outer_fits=("coefficient", "size"),
            coefficient_mean=("coefficient", "mean"),
            coefficient_sd=("coefficient", "std"),
            coefficient_median=("coefficient", "median"),
            positive_fraction=("coefficient", lambda values: float((values > 0).mean())),
            odds_ratio_median=("odds_ratio_per_transformed_unit", "median"),
        )
    )
