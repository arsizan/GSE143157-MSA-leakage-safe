#!/usr/bin/env python3
"""Step 2: donor-level PCA, correlations, and outlier evidence for GSE143157."""

from __future__ import annotations

import argparse
import gc
import html
import json
import math
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import h5py
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from step2_core import (
    median_impute_rows,
    pca_from_features_by_samples,
    robust_z,
    safe_json_number,
    sample_correlations,
    threshold_flags,
)


VERSION = "1.0.0"
CONTROL_COLOR = "#0072B2"
MSA_COLOR = "#D55E00"
NEUTRAL_COLOR = "#666666"
FLAG_COLOR = "#CC0000"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="GSE143157 Step 2: memory-safe PCA, sample correlation, and outlier evidence",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--step1-dir", required=True, type=Path, help="Completed Step1_results_full78 folder")
    parser.add_argument("--output", required=True, type=Path, help="New or empty Step 2 output folder")
    parser.add_argument("--top-variable-cpgs", type=int, default=50000)
    parser.add_argument("--chunk-rows", type=int, default=5000)
    parser.add_argument("--pca-components", type=int, default=10)
    parser.add_argument("--outlier-pcs", type=int, default=5)
    parser.add_argument("--outlier-z-threshold", type=float, default=3.5)
    parser.add_argument("--detection-failure-threshold", type=float, default=0.01)
    return parser.parse_args()


def fail(message: str) -> None:
    raise ValueError(message)


def decode_strings(values) -> list[str]:
    return [item.decode("utf-8") if isinstance(item, bytes) else str(item) for item in values]


def validate_inputs(step1_dir: Path) -> dict[str, Path]:
    files = {
        "h5": step1_dir / "GSE143157_step1_matrices.h5",
        "donors": step1_dir / "metadata_donors_aligned.csv",
        "donor_qc": step1_dir / "qc_per_donor.csv",
        "array_qc": step1_dir / "qc_per_array.csv",
        "summary": step1_dir / "qc_summary.json",
        "success": step1_dir / "SUCCESS.txt",
    }
    missing = [str(path) for path in files.values() if not path.is_file()]
    if missing:
        fail(f"Step 1 files are missing: {missing}")
    return files


def select_top_variable_features(
    values_ds: h5py.Dataset,
    ids_ds: h5py.Dataset,
    top_k: int,
    chunk_rows: int,
    label: str,
) -> tuple[np.ndarray, pd.DataFrame, dict[str, int | float]]:
    n_features, n_donors = values_ds.shape
    variances = np.empty(n_features, dtype=np.float32)
    missing_fractions = np.empty(n_features, dtype=np.float32)
    print(f"[{label}] Pass 1/2: calculating variance for {n_features:,} CpGs...", flush=True)
    for start in range(0, n_features, chunk_rows):
        stop = min(start + chunk_rows, n_features)
        raw = values_ds[start:stop, :]
        missing_fractions[start:stop] = (~np.isfinite(raw)).mean(axis=1)
        imputed, _ = median_impute_rows(raw)
        variances[start:stop] = np.var(imputed, axis=1, ddof=1)
        if stop % 100000 < chunk_rows or stop == n_features:
            print(f"[{label}] variance progress: {stop:,}/{n_features:,}", flush=True)

    eligible = np.flatnonzero(np.isfinite(variances) & (variances > 0))
    if eligible.size < 2:
        fail(f"{label}: fewer than two variable CpGs are available.")
    selected_count = min(int(top_k), int(eligible.size))
    eligible_variances = variances[eligible]
    if selected_count == eligible.size:
        selected = eligible
    else:
        local = np.argpartition(eligible_variances, -selected_count)[-selected_count:]
        selected = eligible[local]
    selected = np.sort(selected)
    selected_mask = np.zeros(n_features, dtype=bool)
    selected_mask[selected] = True

    matrix = np.empty((selected_count, n_donors), dtype=np.float32)
    selected_ids: list[str] = []
    selected_indices: list[int] = []
    write_at = 0
    imputed_cells = 0
    print(f"[{label}] Pass 2/2: loading top {selected_count:,} variable CpGs...", flush=True)
    for start in range(0, n_features, chunk_rows):
        stop = min(start + chunk_rows, n_features)
        local_positions = np.flatnonzero(selected_mask[start:stop])
        if local_positions.size == 0:
            continue
        raw = values_ds[start:stop, :][local_positions]
        imputed, missing_here = median_impute_rows(raw)
        count = imputed.shape[0]
        matrix[write_at : write_at + count] = imputed
        global_positions = (local_positions + start).astype(int)
        selected_indices.extend(global_positions.tolist())
        selected_ids.extend(decode_strings(ids_ds[global_positions]))
        imputed_cells += missing_here
        write_at += count
    if write_at != selected_count:
        fail(f"{label}: selected feature load mismatch ({write_at} vs {selected_count}).")

    table = pd.DataFrame(
        {
            "hdf5_row_index": selected_indices,
            "cpg_id": selected_ids,
            "variance_after_pca_imputation": variances[selected_indices],
            "missing_fraction_before_pca_imputation": missing_fractions[selected_indices],
        }
    ).sort_values("variance_after_pca_imputation", ascending=False, kind="stable")
    table.insert(0, "variance_rank", np.arange(1, len(table) + 1))
    details = {
        "input_cpgs": int(n_features),
        "positive_variance_cpgs": int(eligible.size),
        "selected_cpgs": int(selected_count),
        "selected_missing_cells_median_imputed_for_qc_only": int(imputed_cells),
    }
    return matrix, table, details


def build_modality_analysis(
    label: str,
    matrix: np.ndarray,
    metadata: pd.DataFrame,
    n_components: int,
    outlier_pcs: int,
    threshold: float,
) -> dict[str, object]:
    pca = pca_from_features_by_samples(matrix, n_components=n_components, outlier_pcs=outlier_pcs)
    correlations, mean_correlation, correlation_z = sample_correlations(matrix)
    pca_flag, corr_flag, modality_flag = threshold_flags(
        pca["robust_distance_z"], pca["max_abs_pc_z"], correlation_z, threshold
    )

    scores = metadata.copy()
    score_values = np.asarray(pca["scores"])
    for i in range(score_values.shape[1]):
        scores[f"PC{i + 1}"] = score_values[:, i]
    scores["pca_robust_distance"] = pca["robust_distance"]
    scores["pca_robust_distance_z"] = pca["robust_distance_z"]
    scores["pca_max_abs_component_z"] = pca["max_abs_pc_z"]
    scores["mean_sample_correlation"] = mean_correlation
    scores["mean_sample_correlation_z"] = correlation_z
    scores["flag_pca"] = pca_flag
    scores["flag_low_correlation"] = corr_flag
    scores["flag_modality_review"] = modality_flag

    variance = pd.DataFrame(
        {
            "component": [f"PC{i + 1}" for i in range(len(pca["eigenvalues"]))],
            "eigenvalue": pca["eigenvalues"],
            "explained_variance_ratio": pca["explained_ratio"],
            "cumulative_explained_variance_ratio": pca["cumulative_ratio"],
        }
    )
    return {
        "label": label,
        "scores": scores,
        "variance": variance,
        "correlations": correlations,
        "outlier_pcs_used": int(pca["outlier_pcs_used"]),
    }


def colors_for_condition(conditions: pd.Series) -> list[str]:
    return [MSA_COLOR if value == "MSA" else CONTROL_COLOR for value in conditions]


def annotate_candidates(ax, frame: pd.DataFrame, x: str, y: str, flag: str) -> None:
    labels = frame.loc[frame[flag] | frame["donor_id"].eq("MSA 30")]
    for _, row in labels.iterrows():
        ax.annotate(
            row["donor_id"],
            (row[x], row[y]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
            color=FLAG_COLOR if bool(row[flag]) else NEUTRAL_COLOR,
        )


def plot_pca(scores: pd.DataFrame, variance: pd.DataFrame, label: str, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 6.3))
    for condition, color in [("Control", CONTROL_COLOR), ("MSA", MSA_COLOR)]:
        subset = scores.loc[scores["condition"] == condition]
        ax.scatter(subset["PC1"], subset["PC2"], s=48, c=color, label=condition, alpha=0.82, edgecolors="white", linewidths=0.6)
    flagged = scores["flag_modality_review"]
    ax.scatter(scores.loc[flagged, "PC1"], scores.loc[flagged, "PC2"], s=105, facecolors="none", edgecolors=FLAG_COLOR, linewidths=1.5, label="Review candidate")
    annotate_candidates(ax, scores, "PC1", "PC2", "flag_modality_review")
    pc1 = 100 * variance.loc[0, "explained_variance_ratio"]
    pc2 = 100 * variance.loc[1, "explained_variance_ratio"]
    ax.set_xlabel(f"PC1 ({pc1:.2f}% variance)")
    ax.set_ylabel(f"PC2 ({pc2:.2f}% variance)")
    ax.set_title(f"{label} donor PCA — top variable CpGs")
    ax.axhline(0, color="#DDDDDD", linewidth=0.8, zorder=0)
    ax.axvline(0, color="#DDDDDD", linewidth=0.8, zorder=0)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_scree(variance: pd.DataFrame, label: str, output: Path) -> None:
    shown = variance.head(10)
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    ax.bar(shown["component"], 100 * shown["explained_variance_ratio"], color=CONTROL_COLOR, alpha=0.85)
    ax.plot(shown["component"], 100 * shown["cumulative_explained_variance_ratio"], color=MSA_COLOR, marker="o", label="Cumulative")
    ax.set_ylabel("Explained variance (%)")
    ax.set_xlabel("Principal component")
    ax.set_title(f"{label} PCA scree plot")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_correlation_heatmap(correlations: np.ndarray, donors: list[str], label: str, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.2, 8.0))
    image = ax.imshow(correlations, cmap="viridis", vmin=float(np.nanpercentile(correlations, 2)), vmax=1.0, aspect="equal")
    tick_positions = list(range(0, len(donors), 5))
    ax.set_xticks(tick_positions, [donors[i] for i in tick_positions], rotation=90, fontsize=7)
    ax.set_yticks(tick_positions, [donors[i] for i in tick_positions], fontsize=7)
    ax.set_title(f"{label} sample-correlation heatmap")
    ax.set_xlabel("Donor")
    ax.set_ylabel("Donor")
    if "MSA 30" in donors:
        pos = donors.index("MSA 30")
        ax.axhline(pos, color=FLAG_COLOR, linewidth=0.8)
        ax.axvline(pos, color=FLAG_COLOR, linewidth=0.8)
    fig.colorbar(image, ax=ax, label="Pearson correlation", fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_outlier_diagnostics(scores: pd.DataFrame, label: str, threshold: float, output: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(10.0, 7.3), sharex=True)
    x = np.arange(len(scores))
    colors = colors_for_condition(scores["condition"])
    axes[0].scatter(x, scores["pca_robust_distance_z"], c=colors, s=34)
    axes[0].axhline(threshold, color=FLAG_COLOR, linestyle="--", linewidth=1.0, label=f"Flag threshold {threshold:g}")
    axes[0].set_ylabel("Robust PCA-distance z")
    axes[0].legend(frameon=False)
    axes[1].scatter(x, scores["mean_sample_correlation_z"], c=colors, s=34)
    axes[1].axhline(-threshold, color=FLAG_COLOR, linestyle="--", linewidth=1.0)
    axes[1].set_ylabel("Mean-correlation robust z")
    axes[1].set_xlabel("Donors in aligned matrix order")
    flagged = scores["flag_modality_review"] | scores["donor_id"].eq("MSA 30")
    for idx in np.flatnonzero(flagged):
        axes[0].annotate(scores.iloc[idx]["donor_id"], (idx, scores.iloc[idx]["pca_robust_distance_z"]), xytext=(3, 4), textcoords="offset points", fontsize=7)
        axes[1].annotate(scores.iloc[idx]["donor_id"], (idx, scores.iloc[idx]["mean_sample_correlation_z"]), xytext=(3, 4), textcoords="offset points", fontsize=7)
    fig.suptitle(f"{label} outlier evidence (flags are not exclusions)")
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_covariates(scores: pd.DataFrame, variance: pd.DataFrame, label: str, output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 9.0))
    x, y = scores["PC1"], scores["PC2"]
    pc1 = 100 * variance.loc[0, "explained_variance_ratio"]
    pc2 = 100 * variance.loc[1, "explained_variance_ratio"]

    for condition, color in [("Control", CONTROL_COLOR), ("MSA", MSA_COLOR)]:
        mask = scores["condition"].eq(condition)
        axes[0, 0].scatter(x[mask], y[mask], c=color, s=36, label=condition, alpha=0.82)
    axes[0, 0].set_title("Condition")
    axes[0, 0].legend(frameon=False)

    sex_colors = {"Male": "#009E73", "Female": "#CC79A7"}
    for value in scores["sex"].dropna().unique():
        mask = scores["sex"].eq(value)
        axes[0, 1].scatter(x[mask], y[mask], c=sex_colors.get(value, NEUTRAL_COLOR), s=36, label=value, alpha=0.82)
    axes[0, 1].set_title("Sex")
    axes[0, 1].legend(frameon=False)

    bank_values = list(scores["brain_bank"].dropna().unique())
    bank_colors = ["#56B4E9", "#E69F00", "#999999"]
    for i, value in enumerate(bank_values):
        mask = scores["brain_bank"].eq(value)
        label_value = "Bispebjerg" if "Bispebjerg" in value else "MRC London" if "MRC" in value else value
        axes[1, 0].scatter(x[mask], y[mask], c=bank_colors[i % len(bank_colors)], s=36, label=label_value, alpha=0.82)
    axes[1, 0].set_title("Brain bank")
    axes[1, 0].legend(frameon=False)

    ages = pd.to_numeric(scores["age"], errors="coerce")
    scatter = axes[1, 1].scatter(x, y, c=ages, cmap="plasma", s=38, alpha=0.85)
    axes[1, 1].set_title("Age")
    fig.colorbar(scatter, ax=axes[1, 1], label="Age (years)")

    for ax in axes.flat:
        ax.set_xlabel(f"PC1 ({pc1:.2f}%)")
        ax.set_ylabel(f"PC2 ({pc2:.2f}%)")
    fig.suptitle(f"{label} PCA by biological and technical covariates")
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_qc_vs_hmc(outliers: pd.DataFrame, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.3, 6.1))
    for condition, color in [("Control", CONTROL_COLOR), ("MSA", MSA_COLOR)]:
        subset = outliers.loc[outliers["condition"] == condition]
        ax.scatter(
            100 * subset["5hmC_negative_difference_fraction_before_filter"],
            100 * subset["5hmC_missing_fraction_after_probe_filter"],
            c=color,
            s=46,
            label=condition,
            alpha=0.82,
        )
    labels = outliers.loc[outliers["overall_review_candidate"] | outliers["donor_id"].eq("MSA 30")]
    for _, row in labels.iterrows():
        ax.annotate(
            row["donor_id"],
            (100 * row["5hmC_negative_difference_fraction_before_filter"], 100 * row["5hmC_missing_fraction_after_probe_filter"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )
    ax.set_xlabel("Negative BS−oxBS differences before filtering (%)")
    ax.set_ylabel("5hmC missingness after probe filtering (%)")
    ax.set_title("Donor-level 5hmC QC")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def relative_img(path: Path) -> str:
    return html.escape(path.name)


def write_html_report(
    output: Path,
    summary: dict,
    outliers: pd.DataFrame,
    figures: dict[str, Path],
) -> None:
    candidate_table = outliers.loc[outliers["overall_review_candidate"], [
        "donor_id",
        "condition",
        "evidence_count",
        "outlier_evidence",
        "high_priority_review",
    ]].copy()
    if candidate_table.empty:
        table_html = "<p>No automatic review candidates at the documented thresholds.</p>"
    else:
        table_html = candidate_table.to_html(index=False, border=0, classes="data-table", escape=True)
    figure_html = "\n".join(
        f'<section><h2>{html.escape(title)}</h2><img src="{relative_img(path)}" alt="{html.escape(title)}"></section>'
        for title, path in figures.items()
    )
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>GSE143157 Step 2 QC report</title>
<style>
body{{font-family:Arial,sans-serif;max-width:1080px;margin:28px auto;padding:0 18px;color:#222;line-height:1.5}}
h1,h2{{color:#17365d}} .notice{{background:#fff4ce;padding:12px 16px;border-left:5px solid #d6a800}}
.metrics{{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px}}
.metric{{background:#f3f6fa;padding:12px}} .metric strong{{display:block;font-size:1.35rem}}
img{{max-width:100%;height:auto}} section{{margin:34px 0}} .data-table{{border-collapse:collapse;width:100%;font-size:.9rem}}
.data-table th,.data-table td{{padding:7px;border-bottom:1px solid #ddd;text-align:left}}
code{{background:#f4f4f4;padding:2px 4px}}
</style></head><body>
<h1>GSE143157 Step 2: PCA and outlier review</h1>
<p class="notice"><strong>Important:</strong> A review candidate is not an excluded donor. No sample was deleted in Step 2.</p>
<div class="metrics">
<div class="metric"><span>Donors</span><strong>{summary['donors']}</strong></div>
<div class="metric"><span>Variable CpGs per modality</span><strong>{summary['parameters']['top_variable_cpgs']:,}</strong></div>
<div class="metric"><span>Review candidates</span><strong>{summary['review_candidate_count']}</strong></div>
<div class="metric"><span>High-priority review</span><strong>{summary['high_priority_review_count']}</strong></div>
</div>
<h2>Review-candidate table</h2>{table_html}
<h2>Method summary</h2>
<ul><li>Separate 5mC and 5hmC PCA using the top variable CpGs.</li>
<li>Remaining 5hmC missing values were median-imputed per CpG for unsupervised QC only.</li>
<li>CpGs were centered across donors and were not variance-scaled.</li>
<li>Outlier evidence combines robust PCA distance, low mean sample correlation, Detection-P QC, and 5hmC missingness.</li>
<li>These full-cohort QC selections must not be reused as supervised ML feature selection.</li></ul>
{figure_html}
</body></html>"""
    output.write_text(document, encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    started = time.time()
    step1_dir = args.step1_dir.resolve()
    output_dir = args.output.resolve()
    if args.top_variable_cpgs < 100:
        fail("--top-variable-cpgs must be at least 100.")
    if args.chunk_rows < 1:
        fail("--chunk-rows must be positive.")
    if args.pca_components < 5:
        fail("--pca-components must be at least 5 for the documented summary.")
    if args.outlier_pcs < 2 or args.outlier_pcs > args.pca_components:
        fail("--outlier-pcs must be at least 2 and no greater than --pca-components.")
    if args.outlier_z_threshold <= 0:
        fail("--outlier-z-threshold must be positive.")
    files = validate_inputs(step1_dir)
    if output_dir.exists() and any(output_dir.iterdir()):
        fail(f"Output directory is not empty: {output_dir}. Existing analyses are never overwritten.")
    output_dir.mkdir(parents=True, exist_ok=True)

    donor_meta = pd.read_csv(files["donors"])
    donor_qc = pd.read_csv(files["donor_qc"])
    array_qc = pd.read_csv(files["array_qc"])
    step1_summary = json.loads(files["summary"].read_text(encoding="utf-8"))
    if len(donor_meta) != 78 or set(donor_meta["condition"]) != {"Control", "MSA"}:
        fail("Expected 78 donors with Control/MSA labels in Step 1 metadata.")
    if donor_meta["donor_id"].duplicated().any():
        fail("Step 1 donor metadata has duplicate donor IDs.")

    metadata_columns = [
        column for column in ["donor_index", "donor_id", "condition", "sex", "age", "pmi_hours", "brain_bank"]
        if column in donor_meta.columns
    ]
    analysis_meta = donor_meta[metadata_columns].copy()
    modality_results: dict[str, dict[str, object]] = {}
    selection_details: dict[str, dict[str, int | float]] = {}

    with h5py.File(files["h5"], "r") as h5:
        h5_donors = decode_strings(h5["donor_ids"][:])
        if h5_donors != donor_meta["donor_id"].astype(str).tolist():
            fail("HDF5 donor order does not match metadata_donors_aligned.csv.")
        for label, group_name in [("5mC", "5mC"), ("5hmC", "5hmC")]:
            matrix, top_features, details = select_top_variable_features(
                h5[f"{group_name}/values"],
                h5[f"{group_name}/cpg_ids"],
                args.top_variable_cpgs,
                args.chunk_rows,
                label,
            )
            top_features.to_csv(output_dir / f"top_variable_cpgs_{label}.csv", index=False)
            result = build_modality_analysis(
                label,
                matrix,
                analysis_meta,
                args.pca_components,
                args.outlier_pcs,
                args.outlier_z_threshold,
            )
            scores = result["scores"]
            variance = result["variance"]
            correlations = result["correlations"]
            assert isinstance(scores, pd.DataFrame)
            assert isinstance(variance, pd.DataFrame)
            assert isinstance(correlations, np.ndarray)
            scores.to_csv(output_dir / f"pca_scores_{label}.csv", index=False)
            variance.to_csv(output_dir / f"pca_variance_{label}.csv", index=False)
            correlation_frame = pd.DataFrame(correlations, columns=h5_donors)
            correlation_frame.insert(0, "donor_id", h5_donors)
            correlation_frame.to_csv(output_dir / f"sample_correlations_{label}.csv", index=False)

            figures = {
                f"pca_{label}": output_dir / f"pca_{label}.png",
                f"scree_{label}": output_dir / f"scree_{label}.png",
                f"correlation_{label}": output_dir / f"correlation_heatmap_{label}.png",
                f"outlier_{label}": output_dir / f"outlier_diagnostics_{label}.png",
                f"covariates_{label}": output_dir / f"pca_covariates_{label}.png",
            }
            plot_pca(scores, variance, label, figures[f"pca_{label}"])
            plot_scree(variance, label, figures[f"scree_{label}"])
            plot_correlation_heatmap(correlations, h5_donors, label, figures[f"correlation_{label}"])
            plot_outlier_diagnostics(scores, label, args.outlier_z_threshold, figures[f"outlier_{label}"])
            plot_covariates(scores, variance, label, figures[f"covariates_{label}"])

            modality_results[label] = result
            selection_details[label] = details
            del matrix, correlations
            gc.collect()

    outliers = donor_meta.copy()
    for label in ["5mC", "5hmC"]:
        scores = modality_results[label]["scores"]
        assert isinstance(scores, pd.DataFrame)
        columns = [
            "donor_id",
            "pca_robust_distance",
            "pca_robust_distance_z",
            "pca_max_abs_component_z",
            "mean_sample_correlation",
            "mean_sample_correlation_z",
            "flag_pca",
            "flag_low_correlation",
            "flag_modality_review",
        ]
        rename = {column: f"{label}_{column}" for column in columns if column != "donor_id"}
        outliers = outliers.merge(scores[columns].rename(columns=rename), on="donor_id", how="left", validate="one_to_one")

    outliers = outliers.merge(
        donor_qc[
            [
                "donor_id",
                "5hmC_negative_difference_fraction_before_filter",
                "5hmC_missing_fraction_after_probe_filter",
                "5hmC_mean_after_probe_filter",
            ]
        ],
        on="donor_id",
        how="left",
        validate="one_to_one",
    )
    detection_per_donor = (
        array_qc.sort_values("detection_p_gt_0_01_fraction_of_all_cpgs", ascending=False)
        .groupby("donor_id", as_index=False)
        .first()[["donor_id", "matrix_array_column", "conversion", "detection_p_gt_0_01_fraction_of_all_cpgs"]]
        .rename(
            columns={
                "matrix_array_column": "worst_detection_array",
                "conversion": "worst_detection_conversion",
                "detection_p_gt_0_01_fraction_of_all_cpgs": "worst_detection_failure_fraction",
            }
        )
    )
    outliers = outliers.merge(detection_per_donor, on="donor_id", how="left", validate="one_to_one")

    outliers["5hmC_missing_robust_z"] = robust_z(outliers["5hmC_missing_fraction_after_probe_filter"].to_numpy())
    outliers["5hmC_negative_robust_z"] = robust_z(outliers["5hmC_negative_difference_fraction_before_filter"].to_numpy())
    outliers["flag_5hmC_missing_qc"] = outliers["5hmC_missing_robust_z"] > args.outlier_z_threshold
    outliers["flag_5hmC_negative_qc"] = outliers["5hmC_negative_robust_z"] > args.outlier_z_threshold
    outliers["flag_detection_qc"] = outliers["worst_detection_failure_fraction"] > args.detection_failure_threshold

    evidence_columns = [
        "5mC_flag_pca",
        "5mC_flag_low_correlation",
        "5hmC_flag_pca",
        "5hmC_flag_low_correlation",
        "flag_5hmC_missing_qc",
        "flag_5hmC_negative_qc",
        "flag_detection_qc",
    ]
    outliers["evidence_count"] = outliers[evidence_columns].astype(bool).sum(axis=1)
    evidence_labels = {
        "5mC_flag_pca": "5mC PCA",
        "5mC_flag_low_correlation": "5mC low correlation",
        "5hmC_flag_pca": "5hmC PCA",
        "5hmC_flag_low_correlation": "5hmC low correlation",
        "flag_5hmC_missing_qc": "5hmC missingness",
        "flag_5hmC_negative_qc": "negative 5hmC",
        "flag_detection_qc": "Detection-P",
    }
    outliers["outlier_evidence"] = outliers.apply(
        lambda row: "; ".join(label for column, label in evidence_labels.items() if bool(row[column])),
        axis=1,
    )
    outliers["overall_review_candidate"] = outliers["evidence_count"] > 0
    outliers["high_priority_review"] = outliers["evidence_count"] >= 2
    outliers["automatic_exclusion"] = False
    outliers = outliers.sort_values(
        ["high_priority_review", "evidence_count", "5hmC_missing_fraction_after_probe_filter"],
        ascending=[False, False, False],
        kind="stable",
    )
    outliers.to_csv(output_dir / "outlier_review.csv", index=False)

    qc_figure = output_dir / "qc_5hmC_donor_missingness.png"
    plot_qc_vs_hmc(outliers, qc_figure)
    figure_paths = {
        "5mC PCA": output_dir / "pca_5mC.png",
        "5hmC PCA": output_dir / "pca_5hmC.png",
        "5mC outlier diagnostics": output_dir / "outlier_diagnostics_5mC.png",
        "5hmC outlier diagnostics": output_dir / "outlier_diagnostics_5hmC.png",
        "5hmC donor QC": qc_figure,
        "5mC correlations": output_dir / "correlation_heatmap_5mC.png",
        "5hmC correlations": output_dir / "correlation_heatmap_5hmC.png",
        "5mC covariates": output_dir / "pca_covariates_5mC.png",
        "5hmC covariates": output_dir / "pca_covariates_5hmC.png",
        "5mC scree plot": output_dir / "scree_5mC.png",
        "5hmC scree plot": output_dir / "scree_5hmC.png",
    }

    candidate_rows = outliers.loc[outliers["overall_review_candidate"]]
    msa30 = outliers.loc[outliers["donor_id"].eq("MSA 30")].iloc[0]
    summary = {
        "pipeline": {
            "name": "GSE143157 Step 2 PCA and outlier review",
            "version": VERSION,
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": time.time() - started,
        },
        "inputs": {
            "step1_dir": str(step1_dir),
            "step1_pipeline_version": step1_summary["pipeline"]["version"],
            "5mC_shape": [step1_summary["outputs_summary"]["5mC_cpgs"], 78],
            "5hmC_shape": [step1_summary["outputs_summary"]["5hmC_cpgs_retained"], 78],
        },
        "parameters": {
            "top_variable_cpgs": args.top_variable_cpgs,
            "chunk_rows": args.chunk_rows,
            "pca_components_saved": args.pca_components,
            "outlier_pcs_used": args.outlier_pcs,
            "outlier_robust_z_threshold": args.outlier_z_threshold,
            "detection_failure_fraction_threshold": args.detection_failure_threshold,
            "feature_centering": "per CpG across donors",
            "feature_variance_scaling": False,
            "5hmC_pca_imputation": "per-CpG median, unsupervised QC only",
            "automatic_donor_exclusion": False,
        },
        "donors": len(donor_meta),
        "condition_counts": donor_meta["condition"].value_counts().to_dict(),
        "feature_selection": selection_details,
        "explained_variance": {
            label: {
                "PC1": safe_json_number(modality_results[label]["variance"].loc[0, "explained_variance_ratio"]),
                "PC2": safe_json_number(modality_results[label]["variance"].loc[1, "explained_variance_ratio"]),
                "first_5_cumulative": safe_json_number(modality_results[label]["variance"].loc[4, "cumulative_explained_variance_ratio"]),
            }
            for label in ["5mC", "5hmC"]
        },
        "review_candidate_count": int(outliers["overall_review_candidate"].sum()),
        "high_priority_review_count": int(outliers["high_priority_review"].sum()),
        "review_candidates": candidate_rows[
            ["donor_id", "condition", "evidence_count", "outlier_evidence", "high_priority_review"]
        ].to_dict(orient="records"),
        "MSA_30": {
            "evidence_count": int(msa30["evidence_count"]),
            "evidence": msa30["outlier_evidence"],
            "overall_review_candidate": bool(msa30["overall_review_candidate"]),
            "high_priority_review": bool(msa30["high_priority_review"]),
            "automatic_exclusion": False,
        },
        "interpretation_guardrails": [
            "Outlier flags request review; they do not prove technical failure or justify deletion by themselves.",
            "Full-cohort variable CpG selection is for unsupervised QC only and must not be reused as supervised ML feature selection.",
            "Any donor removal must be reported and accompanied by a full-cohort sensitivity analysis.",
        ],
    }
    with (output_dir / "step2_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False, allow_nan=False)
    write_html_report(output_dir / "STEP2_REPORT.html", summary, outliers, figure_paths)
    (output_dir / "SUCCESS.txt").write_text(
        "SUCCESS\n"
        f"Completed UTC: {summary['pipeline']['completed_utc']}\n"
        f"Donors analyzed: {len(donor_meta)}\n"
        f"Review candidates: {summary['review_candidate_count']}\n"
        "No donor was automatically excluded.\n",
        encoding="utf-8",
    )

    print("\nSTEP 2 SUCCESS")
    print(f"Review candidates: {summary['review_candidate_count']}")
    print(f"High-priority review: {summary['high_priority_review_count']}")
    print(f"MSA 30 evidence: {summary['MSA_30']['evidence'] or 'none at automatic thresholds'}")
    print(f"Open report: {output_dir / 'STEP2_REPORT.html'}")


def main() -> int:
    args = parse_args()
    try:
        run(args)
    except (ValueError, OSError, KeyError, pd.errors.ParserError) as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
