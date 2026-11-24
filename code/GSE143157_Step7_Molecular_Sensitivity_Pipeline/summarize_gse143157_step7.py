#!/usr/bin/env python3
"""Summarize the two prespecified GSE143157 Step 7 sensitivity cohorts."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import html
import json
from pathlib import Path
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


VERSION = "1.0.0"
PRIMARY_COHORT = "primary_clean77"
FULL78 = "sensitivity_full78"
CLEAN76 = "sensitivity_no_MSA30_Control11_76"
COHORT_ORDER = (PRIMARY_COHORT, FULL78, CLEAN76)
COHORT_LABELS = {
    PRIMARY_COHORT: "Primary Clean-77",
    FULL78: "Sensitivity Full-78",
    CLEAN76: "Sensitivity Clean-76",
}
MODEL_ORDER = (
    "all_covariates",
    "5mC__molecular_plus_covariates",
    "5hmC__molecular_plus_covariates",
    "combined__molecular_plus_covariates",
)
MODEL_LABELS = {
    "all_covariates": "All covariates",
    "5mC__molecular_plus_covariates": "5mC + covariates",
    "5hmC__molecular_plus_covariates": "5hmC + covariates",
    "combined__molecular_plus_covariates": "Combined + covariates",
}
MODEL_COLORS = {
    "all_covariates": "#6c737a",
    "5mC__molecular_plus_covariates": "#2387bd",
    "5hmC__molecular_plus_covariates": "#1aa77b",
    "combined__molecular_plus_covariates": "#dc6d11",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step4-dir", required=True, type=Path)
    parser.add_argument("--step5-dir", required=True, type=Path)
    parser.add_argument("--step6-dir", required=True, type=Path)
    parser.add_argument("--full78-dir", required=True, type=Path)
    parser.add_argument("--clean76-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
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


def require_files(paths: dict[str, Path]) -> None:
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise ValueError(f"Required Step 7 inputs are missing: {missing}")


def load_inputs(args: argparse.Namespace) -> dict[str, object]:
    step4 = args.step4_dir.resolve()
    step5 = args.step5_dir.resolve()
    step6 = args.step6_dir.resolve()
    full78 = args.full78_dir.resolve()
    clean76 = args.clean76_dir.resolve()
    paths = {
        "step4_success": step4 / "SUCCESS.txt",
        "step5_success": step5 / "SUCCESS.txt",
        "step5_config": step5 / "RUN_CONFIGURATION.json",
        "primary_aggregate": step5 / "aggregate_metrics_primary_clean77.csv",
        "primary_stability": step5 / "feature_stability_primary_clean77.csv",
        "step6_success": step6 / "SUCCESS.txt",
        "primary_candidates": step6 / "primary_stable_cpgs_primary_clean77.csv",
        "full78_success": full78 / "SUCCESS.txt",
        "full78_summary": full78 / "cohort_summary.json",
        "full78_aggregate": full78 / f"aggregate_metrics_{FULL78}.csv",
        "full78_stability": full78 / f"feature_stability_{FULL78}.csv",
        "clean76_success": clean76 / "SUCCESS.txt",
        "clean76_summary": clean76 / "cohort_summary.json",
        "clean76_aggregate": clean76 / f"aggregate_metrics_{CLEAN76}.csv",
        "clean76_stability": clean76 / f"feature_stability_{CLEAN76}.csv",
    }
    for cohort in COHORT_ORDER:
        paths[f"step4_{cohort}"] = step4 / f"aggregate_metrics_{cohort}.csv"
    require_files(paths)

    primary_config = json.loads(paths["step5_config"].read_text(encoding="utf-8"))
    sensitivity_summaries = {
        FULL78: json.loads(paths["full78_summary"].read_text(encoding="utf-8")),
        CLEAN76: json.loads(paths["clean76_summary"].read_text(encoding="utf-8")),
    }
    modes = {summary["parameters"]["analysis_mode"] for summary in sensitivity_summaries.values()}
    if len(modes) != 1:
        raise ValueError("The two sensitivity cohorts were not run in the same analysis mode.")
    mode = next(iter(modes))
    if mode == "FULL_RESEARCH_RUN":
        reference = primary_config["parameters"]
        keys = ("outer_folds", "outer_repeats", "inner_folds", "candidate_pool", "feature_counts", "c_grid", "class_weights", "missing_threshold", "minimum_variance")
        for cohort, summary in sensitivity_summaries.items():
            mismatched = [key for key in keys if summary["parameters"].get(key) != reference.get(key)]
            if mismatched:
                raise ValueError(f"{cohort} does not match the locked Step 5 policy: {mismatched}")
    elif mode != "SMOKE_TEST_NOT_FOR_RESEARCH":
        raise ValueError(f"Unknown sensitivity analysis mode: {mode}")

    candidates = pd.read_csv(paths["primary_candidates"])
    if len(candidates) != 20 or candidates["feature_name"].nunique() != 20:
        raise ValueError("Step 6 primary candidate table is not the expected 20 unique CpGs.")
    return {
        "paths": paths,
        "mode": mode,
        "candidates": candidates,
        "primary_config": primary_config,
        "sensitivity_summaries": sensitivity_summaries,
    }


def load_performance(paths: dict[str, Path]) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    molecular_paths = {
        PRIMARY_COHORT: paths["primary_aggregate"],
        FULL78: paths["full78_aggregate"],
        CLEAN76: paths["clean76_aggregate"],
    }
    for cohort, path in molecular_paths.items():
        frame = pd.read_csv(path)
        frame = frame.loc[frame["view"] == "molecular_plus_covariates"].copy()
        if set(frame["model_id"]) != set(MODEL_ORDER[1:]):
            raise ValueError(f"Molecular metrics are incomplete for {cohort}.")
        frame["cohort"] = cohort
        frame["cohort_label"] = COHORT_LABELS[cohort]
        rows.append(frame)

        baseline = pd.read_csv(paths[f"step4_{cohort}"])
        baseline = baseline.loc[baseline["model"] == "all_covariates"].copy()
        if len(baseline) != 1:
            raise ValueError(f"Step 4 all-covariate benchmark is incomplete for {cohort}.")
        baseline["modality"] = "covariates"
        baseline["view"] = "covariate_only"
        baseline["model_id"] = "all_covariates"
        baseline["cohort"] = cohort
        baseline["cohort_label"] = COHORT_LABELS[cohort]
        rows.append(baseline)

    combined = pd.concat(rows, ignore_index=True, sort=False)
    combined["model_label"] = combined["model_id"].map(MODEL_LABELS)
    primary_auc = combined.loc[combined["cohort"] == PRIMARY_COHORT].set_index("model_id")["roc_auc_donor_mean_oof"]
    combined["delta_roc_auc_vs_primary_clean77"] = combined.apply(
        lambda row: float(row["roc_auc_donor_mean_oof"] - primary_auc.loc[row["model_id"]]), axis=1
    )
    return combined.sort_values(["cohort", "model_id"]).reset_index(drop=True)


def coefficient_consistency(frame: pd.DataFrame) -> pd.Series:
    fraction = pd.to_numeric(frame["positive_coefficient_fraction"], errors="coerce")
    return np.maximum(fraction, 1.0 - fraction)


def direction(value: float) -> str:
    if not np.isfinite(value) or value == 0:
        return "not_observed"
    return "MSA_higher" if value > 0 else "MSA_lower"


def load_stability(paths: dict[str, Path]) -> dict[str, pd.DataFrame]:
    tables = {
        PRIMARY_COHORT: pd.read_csv(paths["primary_stability"]),
        FULL78: pd.read_csv(paths["full78_stability"]),
        CLEAN76: pd.read_csv(paths["clean76_stability"]),
    }
    for cohort, frame in tables.items():
        frame = frame.loc[frame["model_id"] == "combined__molecular_plus_covariates"].copy()
        frame["coefficient_direction_consistency"] = coefficient_consistency(frame)
        frame["effect_direction"] = frame["mean_effect_msa_minus_control"].map(direction)
        if frame["feature_name"].duplicated().any():
            raise ValueError(f"Duplicate combined-model stability feature in {cohort}.")
        tables[cohort] = frame
    return tables


def make_candidate_tracking(candidates: pd.DataFrame, stability: dict[str, pd.DataFrame]) -> pd.DataFrame:
    keep = [
        "feature_name", "source_modality", "cpg_id", "gene_symbols", "stability_tier",
        "manifest_qc_status", "interpretation_eligible", "selection_frequency",
        "coefficient_direction_consistency", "mean_effect_msa_minus_control", "effect_direction",
    ]
    result = candidates[keep].copy().rename(columns={
        "selection_frequency": "primary_selection_frequency",
        "coefficient_direction_consistency": "primary_direction_consistency",
        "mean_effect_msa_minus_control": "primary_mean_effect_msa_minus_control",
        "effect_direction": "primary_effect_direction",
    })
    for cohort in (FULL78, CLEAN76):
        lookup = stability[cohort].set_index("feature_name")
        short = "full78" if cohort == FULL78 else "clean76"
        result[f"{short}_selected_outer_fits"] = result["feature_name"].map(lookup["selected_outer_fits"]).fillna(0).astype(int)
        result[f"{short}_selection_frequency"] = result["feature_name"].map(lookup["selection_frequency"]).fillna(0.0)
        result[f"{short}_direction_consistency"] = result["feature_name"].map(lookup["coefficient_direction_consistency"])
        result[f"{short}_mean_effect_msa_minus_control"] = result["feature_name"].map(lookup["mean_effect_msa_minus_control"])
        result[f"{short}_effect_direction"] = result[f"{short}_mean_effect_msa_minus_control"].map(direction)
        result[f"{short}_effect_direction_matches_primary"] = (
            result[f"{short}_effect_direction"] == result["primary_effect_direction"]
        ) & (result[f"{short}_effect_direction"] != "not_observed")

    result["minimum_sensitivity_selection_frequency"] = result[[
        "full78_selection_frequency", "clean76_selection_frequency"
    ]].min(axis=1)
    result["effect_direction_matches_both"] = (
        result["full78_effect_direction_matches_primary"]
        & result["clean76_effect_direction_matches_primary"]
    )
    result["sensitivity_supported_both"] = (
        (result["minimum_sensitivity_selection_frequency"] >= 0.60)
        & result["effect_direction_matches_both"]
    )
    one_supported = (
        (result[["full78_selection_frequency", "clean76_selection_frequency"]].max(axis=1) >= 0.60)
        & (
            result["full78_effect_direction_matches_primary"]
            | result["clean76_effect_direction_matches_primary"]
        )
    )
    result["sensitivity_support_class"] = np.where(
        result["sensitivity_supported_both"],
        "supported_in_both_at_prespecified_threshold",
        np.where(one_supported, "supported_in_one_cohort_only", "not_supported_at_0.60_in_either"),
    )
    return result.sort_values(
        ["sensitivity_supported_both", "minimum_sensitivity_selection_frequency", "primary_selection_frequency"],
        ascending=[False, False, False],
    ).reset_index(drop=True)


def make_overlap(stability: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, dict[str, set[str]]]:
    sets: dict[str, set[str]] = {}
    for cohort, frame in stability.items():
        selected = frame.loc[
            (frame["selection_frequency"] >= 0.60)
            & (frame["coefficient_direction_consistency"] >= 0.80),
            "feature_name",
        ]
        sets[cohort] = set(selected.astype(str))
    rows = []
    for cohort_a in COHORT_ORDER:
        for cohort_b in COHORT_ORDER:
            union = sets[cohort_a] | sets[cohort_b]
            intersection = sets[cohort_a] & sets[cohort_b]
            rows.append({
                "cohort_a": cohort_a,
                "cohort_b": cohort_b,
                "stable_set_a": len(sets[cohort_a]),
                "stable_set_b": len(sets[cohort_b]),
                "intersection": len(intersection),
                "union": len(union),
                "jaccard": len(intersection) / len(union) if union else 1.0,
            })
    return pd.DataFrame(rows), sets


def plot_performance(performance: pd.DataFrame, output: Path) -> None:
    fig, axis = plt.subplots(figsize=(12.5, 6.5))
    x = np.arange(len(COHORT_ORDER), dtype=float)
    width = 0.18
    offsets = np.linspace(-1.5 * width, 1.5 * width, len(MODEL_ORDER))
    for offset, model in zip(offsets, MODEL_ORDER):
        frame = performance.loc[performance["model_id"] == model].set_index("cohort").loc[list(COHORT_ORDER)]
        values = frame["roc_auc_donor_mean_oof"].to_numpy(dtype=float)
        lower = values - frame["roc_auc_bootstrap_low"].to_numpy(dtype=float)
        upper = frame["roc_auc_bootstrap_high"].to_numpy(dtype=float) - values
        axis.bar(x + offset, values, width=width, color=MODEL_COLORS[model], label=MODEL_LABELS[model])
        axis.errorbar(x + offset, values, yerr=np.vstack([lower, upper]), fmt="none", ecolor="#222222", capsize=3, linewidth=1)
    axis.axhline(0.5, color="#777777", linestyle="--", linewidth=1)
    axis.set_xticks(x, [COHORT_LABELS[item] for item in COHORT_ORDER])
    axis.set_ylim(0, 1.03)
    axis.set_ylabel("Donor-mean out-of-fold ROC-AUC (95% bootstrap CI)")
    axis.set_title("Prespecified molecular sensitivity analysis across overlapping cohorts")
    axis.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_candidate_heatmap(tracking: pd.DataFrame, output: Path) -> None:
    columns = ["primary_selection_frequency", "full78_selection_frequency", "clean76_selection_frequency"]
    matrix = tracking.set_index("feature_name")[columns].to_numpy(dtype=float)
    labels = []
    for _, row in tracking.iterrows():
        gene = str(row["gene_symbols"]) if pd.notna(row["gene_symbols"]) and str(row["gene_symbols"]).strip() else "IGR"
        labels.append(f"{row['source_modality']} {row['cpg_id']} {gene}")
    fig, axis = plt.subplots(figsize=(8.8, 10.5))
    image = axis.imshow(matrix, aspect="auto", cmap="YlGnBu", vmin=0, vmax=1)
    axis.set_xticks(np.arange(3), ["Primary Clean-77", "Full-78", "Clean-76"])
    axis.set_yticks(np.arange(len(labels)), labels)
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            color = "white" if matrix[row, column] >= 0.62 else "black"
            axis.text(column, row, f"{matrix[row, column]:.2f}", ha="center", va="center", color=color, fontsize=8)
    axis.set_title("Selection frequency of the 20 primary CpG candidates")
    fig.colorbar(image, ax=axis, label="Selection frequency across outer fits", shrink=0.8)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_effects(tracking: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.2), sharex=True, sharey=True)
    limits = tracking[[
        "primary_mean_effect_msa_minus_control",
        "full78_mean_effect_msa_minus_control",
        "clean76_mean_effect_msa_minus_control",
    ]].to_numpy(dtype=float)
    finite = limits[np.isfinite(limits)]
    bound = max(abs(finite.min()), abs(finite.max())) * 1.12 if len(finite) else 0.2
    for axis, short, title in zip(axes, ("full78", "clean76"), ("Full-78", "Clean-76")):
        x = tracking["primary_mean_effect_msa_minus_control"].to_numpy(dtype=float) * 100
        y = tracking[f"{short}_mean_effect_msa_minus_control"].to_numpy(dtype=float) * 100
        colors = np.where(tracking["source_modality"] == "5mC", "#2387bd", "#1aa77b")
        axis.scatter(x, y, c=colors, s=48, alpha=0.9)
        axis.plot([-bound * 100, bound * 100], [-bound * 100, bound * 100], linestyle="--", color="#777777")
        axis.axhline(0, color="#aaaaaa", linewidth=0.8)
        axis.axvline(0, color="#aaaaaa", linewidth=0.8)
        axis.set_xlim(-bound * 100, bound * 100)
        axis.set_ylim(-bound * 100, bound * 100)
        axis.set_title(title)
        axis.set_xlabel("Primary MSA − Control effect (beta percentage points)")
    axes[0].set_ylabel("Sensitivity-cohort MSA − Control effect (beta percentage points)")
    fig.suptitle("Direction and magnitude of primary candidate effects")
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_overlap(overlap: pd.DataFrame, output: Path) -> None:
    matrix = overlap.pivot(index="cohort_a", columns="cohort_b", values="jaccard").loc[list(COHORT_ORDER), list(COHORT_ORDER)].to_numpy(dtype=float)
    fig, axis = plt.subplots(figsize=(7.4, 6.2))
    image = axis.imshow(matrix, cmap="Blues", vmin=0, vmax=1)
    labels = [COHORT_LABELS[item] for item in COHORT_ORDER]
    axis.set_xticks(np.arange(3), labels, rotation=20, ha="right")
    axis.set_yticks(np.arange(3), labels)
    for row in range(3):
        for column in range(3):
            axis.text(column, row, f"{matrix[row, column]:.2f}", ha="center", va="center", color="white" if matrix[row, column] > 0.55 else "black")
    axis.set_title("Jaccard overlap of stable combined-model CpG sets")
    fig.colorbar(image, ax=axis, label="Jaccard index", shrink=0.82)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_report(output: Path, performance: pd.DataFrame, tracking: pd.DataFrame, overlap: pd.DataFrame, mode: str) -> None:
    primary_model = performance.loc[performance["model_id"] == "combined__molecular_plus_covariates", [
        "cohort_label", "donors", "roc_auc_donor_mean_oof", "roc_auc_bootstrap_low", "roc_auc_bootstrap_high",
        "pr_auc_donor_mean_oof", "balanced_accuracy_donor_mean_oof", "delta_roc_auc_vs_primary_clean77",
    ]]
    supported = tracking.loc[tracking["sensitivity_supported_both"], [
        "source_modality", "cpg_id", "gene_symbols", "stability_tier",
        "primary_selection_frequency", "full78_selection_frequency", "clean76_selection_frequency",
        "primary_effect_direction",
    ]]
    pairwise = overlap.loc[overlap["cohort_a"] < overlap["cohort_b"]]
    warning = "SMOKE TEST ONLY — numbers in this report are not research results." if mode != "FULL_RESEARCH_RUN" else "The two sensitivity cohorts overlap the primary cohort and are robustness analyses, not external validation."
    document = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>GSE143157 Step 7 sensitivity report</title>
<style>body{{font-family:Segoe UI,Arial,sans-serif;max-width:1250px;margin:32px auto;padding:0 20px;color:#1d2630}}h1,h2{{color:#173b5f}}.good{{background:#eaf7ee;border-left:5px solid #238b45;padding:12px}}.warn{{background:#fff2df;border-left:5px solid #d46b08;padding:12px}}table{{border-collapse:collapse;width:100%;font-size:13px;margin:10px 0 18px}}th,td{{border:1px solid #d7dee5;padding:7px;text-align:left}}th{{background:#eef3f7}}img{{max-width:100%;margin:10px 0 24px}}</style></head><body>
<h1>GSE143157 Step 7 — Prespecified sensitivity validation</h1>
<p class="good"><strong>Design locked before review:</strong> only 5mC, 5hmC and combined molecular-plus-covariate models were re-estimated using the frozen Step 3 nested-CV folds. All molecular preprocessing and selection remained outer-training-only.</p>
<p class="warn"><strong>Interpretation:</strong> {html.escape(warning)} Candidate sensitivity support requires selection frequency ≥0.60 and matching mean-effect direction in both alternative cohorts. This threshold is descriptive and is not a p-value.</p>
<h2>Performance robustness</h2>{primary_model.to_html(index=False,border=0,float_format=lambda value:f'{value:.4f}')}
<img src="performance_across_cohorts.png" alt="Performance across cohorts">
<h2>Primary candidate tracking</h2><p>{len(supported)} of 20 primary candidates met the prespecified sensitivity-support rule in both alternative cohorts.</p>
{supported.to_html(index=False,border=0,float_format=lambda value:f'{value:.3f}')}
<img src="primary_candidate_selection_heatmap.png" alt="Primary candidate stability heatmap">
<img src="candidate_effect_concordance.png" alt="Candidate effect concordance">
<h2>Stable-set overlap</h2>{pairwise.to_html(index=False,border=0,float_format=lambda value:f'{value:.3f}')}
<img src="stable_set_overlap.png" alt="Stable-set overlap">
<h2>Reporting guardrail</h2><p>These results assess sensitivity to the inclusion or exclusion of MSA 30 and Control 11. They cannot establish generalization to a new brain bank or an independent cohort.</p>
</body></html>"""
    output.write_text(document, encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    started = time.time()
    output = args.output.resolve()
    inputs = load_inputs(args)
    if (output / "SUCCESS.txt").exists():
        raise ValueError(f"A successful Step 7 result already exists: {output}")
    output.mkdir(parents=True, exist_ok=True)
    paths = inputs["paths"]
    performance = load_performance(paths)
    stability = load_stability(paths)
    tracking = make_candidate_tracking(inputs["candidates"], stability)
    overlap, stable_sets = make_overlap(stability)

    performance.to_csv(output / "aggregate_metrics_across_cohorts.csv", index=False)
    tracking.to_csv(output / "primary_candidate_sensitivity.csv", index=False)
    overlap.to_csv(output / "stable_set_overlap.csv", index=False)

    plot_performance(performance, output / "performance_across_cohorts.png")
    plot_candidate_heatmap(tracking, output / "primary_candidate_selection_heatmap.png")
    plot_effects(tracking, output / "candidate_effect_concordance.png")
    plot_overlap(overlap, output / "stable_set_overlap.png")
    write_report(output / "STEP7_REPORT.html", performance, tracking, overlap, inputs["mode"])

    configuration = {
        "pipeline_version": VERSION,
        "analysis_mode": inputs["mode"],
        "locked_model_view": "molecular_plus_covariates",
        "sensitivity_support_rule": {
            "minimum_selection_frequency_in_each_alternative_cohort": 0.60,
            "mean_effect_direction_must_match_primary_in_each_alternative_cohort": True,
            "descriptive_not_p_value": True,
        },
        "input_sha256": {name: sha256(path) for name, path in paths.items() if path.suffix in {".csv", ".json", ".txt"}},
    }
    (output / "RUN_CONFIGURATION.json").write_text(json.dumps(configuration, indent=2, ensure_ascii=False), encoding="utf-8")
    summary = {
        "pipeline": {
            "name": "GSE143157 Step 7 prespecified overlapping-cohort sensitivity validation",
            "version": VERSION,
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": time.time() - started,
        },
        "analysis_mode": inputs["mode"],
        "cohorts": {
            PRIMARY_COHORT: {"donors": 77, "control": 37, "msa": 40, "excluded": ["MSA 30"]},
            FULL78: {"donors": 78, "control": 37, "msa": 41, "excluded": []},
            CLEAN76: {"donors": 76, "control": 36, "msa": 40, "excluded": ["MSA 30", "Control 11"]},
        },
        "counts": {
            "primary_candidates": int(len(tracking)),
            "primary_candidates_supported_in_both": int(tracking["sensitivity_supported_both"].sum()),
            "stable_combined_model_features_primary": len(stable_sets[PRIMARY_COHORT]),
            "stable_combined_model_features_full78": len(stable_sets[FULL78]),
            "stable_combined_model_features_clean76": len(stable_sets[CLEAN76]),
        },
        "combined_model_performance": performance.loc[
            performance["model_id"] == "combined__molecular_plus_covariates",
            ["cohort", "donors", "roc_auc_donor_mean_oof", "roc_auc_bootstrap_low", "roc_auc_bootstrap_high", "pr_auc_donor_mean_oof", "balanced_accuracy_donor_mean_oof", "delta_roc_auc_vs_primary_clean77"],
        ].to_dict(orient="records"),
        "guardrails": [
            "The sensitivity cohorts were prespecified in Step 2B and their donor folds were frozen in Step 3.",
            "Only the Step 6-locked molecular-plus-covariate model view was re-estimated.",
            "All variance filtering, imputation, scaling, supervised CpG selection and tuning were nested inside outer training data.",
            "The cohorts overlap and are not independent external validation datasets.",
            "Sensitivity-support thresholds are descriptive stability rules, not statistical significance tests.",
        ],
    }
    (output / "step7_summary.json").write_text(json.dumps(json_clean(summary), indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    success = (
        "SUCCESS\n"
        f"Completed UTC: {summary['pipeline']['completed_utc']}\n"
        f"Analysis mode: {inputs['mode']}\n"
        "Cohorts: Primary Clean-77, Sensitivity Full-78, Sensitivity Clean-76\n"
        "Locked view: molecular_plus_covariates\n"
        f"Primary candidates supported in both sensitivity cohorts: {summary['counts']['primary_candidates_supported_in_both']}/20\n"
        "Interpretation: overlapping-cohort robustness analysis, not independent external validation\n"
    )
    (output / "SUCCESS.txt").write_text(success, encoding="utf-8")
    print("\nSTEP 7 SUCCESS")
    if inputs["mode"] != "FULL_RESEARCH_RUN":
        print("SMOKE TEST ONLY - DO NOT USE THESE NUMBERS IN THE PAPER")
    print(f"Open report: {output / 'STEP7_REPORT.html'}")


def main() -> int:
    try:
        run(parse_args())
    except (ValueError, OSError, KeyError, json.JSONDecodeError, pd.errors.ParserError) as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
