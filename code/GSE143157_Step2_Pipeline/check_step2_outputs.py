#!/usr/bin/env python3
"""Validate the generated GSE143157 Step 2 result folder."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser(description="Check GSE143157 Step 2 outputs")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    root = args.output.resolve()
    required = [
        "SUCCESS.txt",
        "step2_summary.json",
        "STEP2_REPORT.html",
        "outlier_review.csv",
        "pca_scores_5mC.csv",
        "pca_scores_5hmC.csv",
        "pca_variance_5mC.csv",
        "pca_variance_5hmC.csv",
        "sample_correlations_5mC.csv",
        "sample_correlations_5hmC.csv",
        "top_variable_cpgs_5mC.csv",
        "top_variable_cpgs_5hmC.csv",
        "pca_5mC.png",
        "pca_5hmC.png",
        "correlation_heatmap_5mC.png",
        "correlation_heatmap_5hmC.png",
        "outlier_diagnostics_5mC.png",
        "outlier_diagnostics_5hmC.png",
        "pca_covariates_5mC.png",
        "pca_covariates_5hmC.png",
        "scree_5mC.png",
        "scree_5hmC.png",
        "qc_5hmC_donor_missingness.png",
    ]
    missing = [name for name in required if not (root / name).is_file()]
    if missing:
        print(f"FAILED: missing files: {missing}", file=sys.stderr)
        return 2

    summary = json.loads((root / "step2_summary.json").read_text(encoding="utf-8"))
    outliers = pd.read_csv(root / "outlier_review.csv")
    errors = []
    if len(outliers) != 78 or outliers["donor_id"].nunique() != 78:
        errors.append(f"outlier_review.csv must contain 78 unique donors; observed {len(outliers)} rows")
    if outliers["automatic_exclusion"].astype(bool).any():
        errors.append("Step 2 must not automatically exclude any donor")
    if int(outliers["overall_review_candidate"].astype(bool).sum()) != summary["review_candidate_count"]:
        errors.append("Review-candidate count differs between CSV and JSON")

    for label in ["5mC", "5hmC"]:
        scores = pd.read_csv(root / f"pca_scores_{label}.csv")
        variance = pd.read_csv(root / f"pca_variance_{label}.csv")
        correlations = pd.read_csv(root / f"sample_correlations_{label}.csv")
        top = pd.read_csv(root / f"top_variable_cpgs_{label}.csv")
        if len(scores) != 78 or scores["donor_id"].nunique() != 78:
            errors.append(f"{label} PCA scores do not contain 78 unique donors")
        if scores[["PC1", "PC2"]].isna().any().any():
            errors.append(f"{label} PC1/PC2 contains missing values")
        if len(correlations) != 78 or correlations.shape[1] != 79:
            errors.append(f"{label} correlation table is not 78 x 78 plus donor_id")
        if len(top) != summary["feature_selection"][label]["selected_cpgs"]:
            errors.append(f"{label} top-variable CpG count differs from JSON")
        explained = variance["explained_variance_ratio"]
        if (explained < 0).any() or explained.sum() > 1.000001:
            errors.append(f"{label} explained-variance ratios are invalid")

    small_figures = [name for name in required if name.endswith(".png") and (root / name).stat().st_size < 5000]
    if small_figures:
        errors.append(f"Figures appear empty or too small: {small_figures}")

    if errors:
        print("FAILED")
        for error in errors:
            print(f"- {error}")
        return 2
    print("PASS: Step 2 outputs are internally consistent.")
    print(f"Review candidates: {summary['review_candidate_count']}")
    print(f"High-priority review: {summary['high_priority_review_count']}")
    print("No donor was automatically excluded.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
