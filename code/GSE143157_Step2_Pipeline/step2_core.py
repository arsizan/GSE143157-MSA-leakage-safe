"""Numerical core for GSE143157 Step 2 (no file or plotting dependencies)."""

from __future__ import annotations

import math

import numpy as np


def robust_z(values: np.ndarray) -> np.ndarray:
    """Median/MAD z-score; deterministic standard-deviation fallback."""
    array = np.asarray(values, dtype=np.float64)
    median = np.nanmedian(array)
    mad = np.nanmedian(np.abs(array - median))
    if np.isfinite(mad) and mad > 0:
        return 0.6744897501960817 * (array - median) / mad
    std = np.nanstd(array, ddof=1)
    if np.isfinite(std) and std > 0:
        return (array - np.nanmean(array)) / std
    return np.zeros_like(array)


def median_impute_rows(matrix: np.ndarray) -> tuple[np.ndarray, int]:
    """Return a float32 copy with non-finite cells replaced by row medians."""
    values = np.asarray(matrix, dtype=np.float32).copy()
    missing = ~np.isfinite(values)
    missing_count = int(missing.sum())
    if not missing_count:
        return values, 0
    with np.errstate(all="ignore"):
        medians = np.nanmedian(values, axis=1)
    if (~np.isfinite(medians)).any():
        bad_rows = np.flatnonzero(~np.isfinite(medians))[:10].tolist()
        raise ValueError(f"Cannot median-impute all-missing feature rows; examples: {bad_rows}")
    row_index, _ = np.where(missing)
    values[missing] = medians[row_index]
    return values, missing_count


def pca_from_features_by_samples(
    matrix: np.ndarray,
    n_components: int = 10,
    outlier_pcs: int = 5,
) -> dict[str, np.ndarray | int]:
    """PCA via a small sample Gram matrix.

    Input rows are CpGs/features and columns are donors/samples. Each feature is
    centered across donors, but not variance-scaled.
    """
    values = np.asarray(matrix, dtype=np.float64)
    if values.ndim != 2 or values.shape[0] < 2 or values.shape[1] < 3:
        raise ValueError(f"PCA requires a 2-D feature x sample matrix; observed {values.shape}")
    if not np.isfinite(values).all():
        raise ValueError("PCA matrix still contains non-finite values after imputation.")

    centered = values - values.mean(axis=1, keepdims=True)
    gram = (centered.T @ centered) / (values.shape[0] - 1)
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = np.clip(eigenvalues[order], 0.0, None)
    eigenvectors = eigenvectors[:, order]

    total = eigenvalues.sum()
    if not np.isfinite(total) or total <= 0:
        raise ValueError("PCA total variance is not positive.")
    explained = eigenvalues / total
    component_count = min(int(n_components), values.shape[1] - 1)
    scores = eigenvectors[:, :component_count] * np.sqrt(
        eigenvalues[:component_count] * (values.shape[0] - 1)
    )

    # Eigenvector signs are arbitrary. Orient each component deterministically.
    for component in range(scores.shape[1]):
        anchor = int(np.argmax(np.abs(scores[:, component])))
        if scores[anchor, component] < 0:
            scores[:, component] *= -1

    used_pcs = min(int(outlier_pcs), scores.shape[1])
    score_robust_z = np.column_stack([robust_z(scores[:, i]) for i in range(used_pcs)])
    robust_distance = np.sqrt(np.sum(score_robust_z**2, axis=1))
    robust_distance_z = robust_z(robust_distance)
    max_abs_pc_z = np.max(np.abs(score_robust_z), axis=1)

    return {
        "scores": scores,
        "eigenvalues": eigenvalues[:component_count],
        "explained_ratio": explained[:component_count],
        "cumulative_ratio": np.cumsum(explained[:component_count]),
        "robust_distance": robust_distance,
        "robust_distance_z": robust_distance_z,
        "max_abs_pc_z": max_abs_pc_z,
        "outlier_pcs_used": used_pcs,
    }


def sample_correlations(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Pearson sample correlation, off-diagonal means, and robust z-scores."""
    values = np.asarray(matrix, dtype=np.float64)
    correlations = np.corrcoef(values.T)
    if not np.isfinite(correlations).all():
        raise ValueError("Sample-correlation matrix contains non-finite values.")
    n_samples = correlations.shape[0]
    mean_to_others = (correlations.sum(axis=1) - 1.0) / (n_samples - 1)
    return correlations, mean_to_others, robust_z(mean_to_others)


def threshold_flags(
    distance_z: np.ndarray,
    max_abs_pc_z: np.ndarray,
    correlation_z: np.ndarray,
    threshold: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return PCA, correlation, and combined modality candidate flags."""
    pca_flag = (np.asarray(distance_z) > threshold) | (np.asarray(max_abs_pc_z) > threshold)
    correlation_flag = np.asarray(correlation_z) < -threshold
    return pca_flag, correlation_flag, pca_flag | correlation_flag


def safe_json_number(value: float) -> float | None:
    number = float(value)
    return number if math.isfinite(number) else None
