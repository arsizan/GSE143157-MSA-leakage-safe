from __future__ import annotations

import argparse
import hashlib
import html
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


PIPELINE_VERSION = "1.0.0"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Freeze the verified GSE143157 Steps 1-8 into paper-ready tables, figures, and guardrails."
    )
    for step in range(1, 9):
        parser.add_argument(f"--step{step}-dir", required=True, type=Path)
    parser.add_argument("--step2b-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def require_success(directory: Path, label: str) -> None:
    if not directory.is_dir():
        raise FileNotFoundError(f"{label} directory not found: {directory}")
    success = directory / "SUCCESS.txt"
    if not success.is_file() or "SUCCESS" not in success.read_text(encoding="utf-8"):
        raise RuntimeError(f"{label} is not complete: {success} is missing or invalid")


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def format_ci(point: float, low: float, high: float) -> str:
    return f"{point:.3f} ({low:.3f}-{high:.3f})"


def copy_figure(source: Path, destination: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(f"Required figure not found: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def build_workflow_figure(output_path: Path, counts: dict[str, Any]) -> None:
    stages = [
        ("Public paired EPIC data", "78 donors; 156 arrays\n37 Control / 41 MSA"),
        ("5mC and 5hmC construction", f"5mC: {counts['cpg_5mc']:,} CpGs\n5hmC: {counts['cpg_5hmc']:,} retained CpGs"),
        ("QC and cohort freeze", "Primary Clean-77\nMSA 30 excluded for convergent technical QC"),
        ("Leakage-safe nested CV", "5 outer folds x 10 repeats\n4 inner folds; training-only feature operations"),
        ("Molecular benchmarking", "5mC, 5hmC, combined\nwith prespecified covariates"),
        ("Annotation and sensitivity", f"20 primary CpGs\n{counts['supported_candidates']}/20 supported in both sensitivity cohorts"),
        ("Bias-aware pathway analysis", "GOmeth; EPIC probe-number and\nmulti-gene bias correction"),
    ]

    fig, ax = plt.subplots(figsize=(10, 13))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, len(stages) * 2)
    ax.axis("off")
    colors = ["#E8F1FA", "#EAF6F1", "#FFF4E5", "#F1ECF8", "#E8F1FA", "#EAF6F1", "#FFF4E5"]
    for idx, ((title, body), color) in enumerate(zip(stages, colors)):
        y = len(stages) * 2 - 1.2 - idx * 2
        box = plt.Rectangle((1.2, y - 0.65), 7.6, 1.3, facecolor=color, edgecolor="#375A7F", linewidth=1.5)
        ax.add_patch(box)
        ax.text(1.55, y + 0.22, title, fontsize=13, fontweight="bold", va="center")
        ax.text(8.45, y - 0.12, body, fontsize=10.5, ha="right", va="center")
        if idx < len(stages) - 1:
            ax.annotate("", xy=(5, y - 1.35), xytext=(5, y - 0.7), arrowprops=dict(arrowstyle="-|>", color="#375A7F", lw=1.6))
    ax.set_title("GSE143157 analysis workflow", fontsize=18, pad=18)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def add_audit(rows: list[dict[str, Any]], item: str, expected: str, observed: Any, passed: bool, evidence: str, notes: str = "") -> None:
    rows.append(
        {
            "audit_item": item,
            "expected": expected,
            "observed": str(observed),
            "status": "PASS" if passed else "FAIL",
            "evidence_file": evidence,
            "notes": notes,
        }
    )


def main() -> None:
    args = parse_args()
    step_dirs = {f"step{i}": getattr(args, f"step{i}_dir").resolve() for i in range(1, 9)}
    step_dirs["step2b"] = args.step2b_dir.resolve()
    output = args.output.resolve()

    for label, directory in step_dirs.items():
        require_success(directory, label.upper())
    if (output / "SUCCESS.txt").is_file():
        raise RuntimeError(f"A completed Step 9 result already exists. Nothing was overwritten: {output}")
    output.mkdir(parents=True, exist_ok=True)
    figures_dir = output / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc)

    step1 = load_json(step_dirs["step1"] / "qc_summary.json")
    step2 = load_json(step_dirs["step2"] / "step2_summary.json")
    step2b = load_json(step_dirs["step2b"] / "cohort_decision.json")
    step3 = load_json(step_dirs["step3"] / "step3_summary.json")
    step4 = load_json(step_dirs["step4"] / "step4_summary.json")
    step5 = load_json(step_dirs["step5"] / "step5_summary.json")
    step6 = load_json(step_dirs["step6"] / "step6_summary.json")
    step7 = load_json(step_dirs["step7"] / "step7_summary.json")
    step8 = load_json(step_dirs["step8"] / "step8_summary.json")
    step8_config = load_json(step_dirs["step8"] / "RUN_CONFIGURATION.json")

    if step7.get("analysis_mode") != "FULL_RESEARCH_RUN":
        raise RuntimeError("Step 7 is not a FULL_RESEARCH_RUN")
    if step8_config.get("prior_prob") is not True or step8_config.get("fractional_multi_gene_counts") is not True:
        raise RuntimeError("Step 8 bias-correction locks are not present")

    # Table 1: cohort definition.
    cohort_order = ["primary_clean77", "sensitivity_full78", "sensitivity_no_MSA30_Control11_76"]
    cohort_labels = {
        "primary_clean77": "Primary Clean-77",
        "sensitivity_full78": "Sensitivity Full-78",
        "sensitivity_no_MSA30_Control11_76": "Sensitivity Clean-76",
    }
    cohort_roles = {
        "primary_clean77": "Primary analysis",
        "sensitivity_full78": "Outlier-inclusion sensitivity",
        "sensitivity_no_MSA30_Control11_76": "Strict outlier stress-test",
    }
    cohort_rows: list[dict[str, Any]] = []
    for cohort in cohort_order:
        frozen = step2b["cohorts"][cohort]
        s7 = step7["cohorts"][cohort]
        cohort_rows.append(
            {
                "cohort": cohort_labels[cohort],
                "role": cohort_roles[cohort],
                "donors_n": int(s7["donors"]),
                "control_n": int(s7["control"]),
                "msa_n": int(s7["msa"]),
                "female_n": int(frozen["sex_counts"]["Female"]),
                "male_n": int(frozen["sex_counts"]["Male"]),
                "age_mean_years": float(frozen["age_mean"]),
                "age_range_years": f"{int(frozen['age_min'])}-{int(frozen['age_max'])}",
                "excluded_donors": "; ".join(s7.get("excluded", [])) or "None",
            }
        )
    cohort_table = pd.DataFrame(cohort_rows)
    cohort_table.to_csv(output / "TABLE1_cohorts.csv", index=False)

    # Table 2: locked model-view performance across all prespecified cohorts.
    performance = pd.read_csv(step_dirs["step7"] / "aggregate_metrics_across_cohorts.csv")
    required_perf = {
        "cohort", "cohort_label", "model_label", "model_id", "donors",
        "roc_auc_donor_mean_oof", "roc_auc_bootstrap_low", "roc_auc_bootstrap_high",
        "pr_auc_donor_mean_oof", "pr_auc_bootstrap_low", "pr_auc_bootstrap_high",
        "balanced_accuracy_donor_mean_oof", "balanced_accuracy_bootstrap_low", "balanced_accuracy_bootstrap_high",
    }
    missing_perf = required_perf.difference(performance.columns)
    if missing_perf:
        raise RuntimeError(f"Step 7 performance table is missing columns: {sorted(missing_perf)}")

    baseline_auc = (
        performance.loc[performance["model_id"] == "all_covariates", ["cohort", "roc_auc_donor_mean_oof"]]
        .set_index("cohort")["roc_auc_donor_mean_oof"]
        .to_dict()
    )
    perf_out = pd.DataFrame(
        {
            "cohort": performance["cohort_label"],
            "donors_n": performance["donors"].astype(int),
            "model": performance["model_label"],
            "model_id": performance["model_id"],
            "roc_auc": performance["roc_auc_donor_mean_oof"],
            "roc_auc_ci_low": performance["roc_auc_bootstrap_low"],
            "roc_auc_ci_high": performance["roc_auc_bootstrap_high"],
            "roc_auc_95ci": [
                format_ci(p, lo, hi)
                for p, lo, hi in zip(
                    performance["roc_auc_donor_mean_oof"],
                    performance["roc_auc_bootstrap_low"],
                    performance["roc_auc_bootstrap_high"],
                )
            ],
            "pr_auc": performance["pr_auc_donor_mean_oof"],
            "pr_auc_ci_low": performance["pr_auc_bootstrap_low"],
            "pr_auc_ci_high": performance["pr_auc_bootstrap_high"],
            "balanced_accuracy": performance["balanced_accuracy_donor_mean_oof"],
            "balanced_accuracy_ci_low": performance["balanced_accuracy_bootstrap_low"],
            "balanced_accuracy_ci_high": performance["balanced_accuracy_bootstrap_high"],
            "delta_roc_auc_vs_all_covariates": [
                float(auc) - float(baseline_auc[cohort])
                for cohort, auc in zip(performance["cohort"], performance["roc_auc_donor_mean_oof"])
            ],
            "validation_scope": "Overlapping-cohort internal sensitivity; not external validation",
        }
    )
    perf_out.to_csv(output / "TABLE2_model_performance.csv", index=False)

    # Table 3: annotated primary candidates joined to the prespecified sensitivity evidence.
    stable = pd.read_csv(step_dirs["step6"] / "primary_stable_cpgs_primary_clean77.csv")
    sensitivity = pd.read_csv(step_dirs["step7"] / "primary_candidate_sensitivity.csv")
    candidates = stable.merge(
        sensitivity,
        on=["feature_name", "source_modality", "cpg_id"],
        how="inner",
        suffixes=("_step6", "_step7"),
        validate="one_to_one",
    )
    if len(candidates) != int(step6["counts"]["primary_stable_cpgs"]):
        raise RuntimeError("Step 6/7 candidate join did not preserve all primary candidates")
    candidates["sensitivity_supported_both"] = as_bool(candidates["sensitivity_supported_both"])
    candidates["effect_direction_matches_both"] = as_bool(candidates["effect_direction_matches_both"])

    def candidate_priority(row: pd.Series) -> str:
        manifest_pass = str(row["manifest_qc_status_step6"]) == "pass_manifest_audit"
        supported = bool(row["sensitivity_supported_both"])
        tier = "Tier A" if str(row["stability_tier_step6"]).startswith("A_") else "Tier B"
        if supported and manifest_pass:
            return f"Core robust, manifest-audited ({tier})"
        if supported:
            return f"Robust but manifest-SNP review required ({tier})"
        if manifest_pass:
            return f"Secondary, cohort-sensitive ({tier})"
        return f"Secondary, cohort-sensitive and manifest review required ({tier})"

    candidates["paper_priority"] = candidates.apply(candidate_priority, axis=1)
    candidates["mean_effect_percentage_points"] = candidates["mean_effect_msa_minus_control"] * 100.0
    candidates["gene_symbols_final"] = candidates["gene_symbols_step6"].fillna("")
    candidates["paper_interpretation"] = "Exploratory CpG association candidate; not an independently validated biomarker"
    candidate_columns = [
        "source_modality", "cpg_id", "gene_symbols_final", "stability_tier_step6",
        "selection_frequency", "coefficient_direction_consistency", "effect_direction",
        "mean_effect_percentage_points", "chromosome", "mapinfo", "ucsc_refgene_group",
        "relation_to_cpg_island", "manifest_qc_status_step6", "primary_selection_frequency",
        "full78_selection_frequency", "clean76_selection_frequency",
        "minimum_sensitivity_selection_frequency", "effect_direction_matches_both",
        "sensitivity_supported_both", "sensitivity_support_class", "paper_priority",
        "paper_interpretation",
    ]
    candidate_table = candidates[candidate_columns].rename(
        columns={
            "gene_symbols_final": "gene_symbols",
            "stability_tier_step6": "stability_tier",
            "effect_direction": "primary_effect_direction",
            "manifest_qc_status_step6": "manifest_qc_status",
            "selection_frequency": "primary_model_selection_frequency",
            "coefficient_direction_consistency": "primary_coefficient_direction_consistency",
        }
    )
    candidate_table = candidate_table.sort_values(
        ["sensitivity_supported_both", "minimum_sensitivity_selection_frequency", "primary_model_selection_frequency"],
        ascending=[False, False, False],
    )
    candidate_table.to_csv(output / "TABLE3_robust_cpg_candidates.csv", index=False)

    # Table 4: bias-aware pathway results. FDR status is explicit to prevent nominal overclaiming.
    pathways = pd.read_csv(step_dirs["step8"] / "top_terms_all.csv")
    pathways["FDR"] = pd.to_numeric(pathways["FDR"], errors="coerce")
    pathways["P.DE"] = pd.to_numeric(pathways["P.DE"], errors="coerce")
    pathways["reporting_status"] = pathways["FDR"].apply(
        lambda value: "FDR-significant" if pd.notna(value) and value < 0.05 else "Exploratory ranking only"
    )
    pathways = (
        pathways.sort_values(["modality", "collection", "P.DE"])
        .groupby(["modality", "collection"], as_index=False, group_keys=False)
        .head(10)
    )
    pathway_columns = [
        "modality", "collection", "term_id", "display_term", "Ont", "N", "DE",
        "P.DE", "FDR", "SigGenesInSet", "reporting_status",
    ]
    pathways[pathway_columns].to_csv(output / "TABLE4_bias_aware_pathways.csv", index=False)
    shutil.copy2(step_dirs["step8"] / "top_terms_all.csv", output / "SUPPLEMENT_all_top_pathways.csv")

    # Paired primary comparison used for the central performance claim.
    comparisons = pd.read_csv(step_dirs["step5"] / "auc_comparisons_primary_clean77.csv")
    target_comparison = comparisons.loc[
        (comparisons["model_a"] == "combined__molecular_plus_covariates")
        & (comparisons["model_b"] == "all_covariates")
    ]
    if len(target_comparison) != 1:
        raise RuntimeError("Could not uniquely locate combined + covariates versus all-covariates comparison")
    target_comparison = target_comparison.iloc[0]

    # Final methods audit.
    audit_rows: list[dict[str, Any]] = []
    add_audit(audit_rows, "Paired donors", "78", step1["observed_structure"]["paired_donors"], step1["observed_structure"]["paired_donors"] == 78, "Step1/qc_summary.json")
    add_audit(audit_rows, "Disease/control counts", "41 MSA / 37 Control", f"{step1['observed_structure']['condition_counts']['MSA']} MSA / {step1['observed_structure']['condition_counts']['Control']} Control", step1["observed_structure"]["condition_counts"] == {"Control": 37, "MSA": 41}, "Step1/qc_summary.json")
    add_audit(audit_rows, "5hmC construction", "BS - oxBS; negatives to missing", f"{step1['construction']['5hmC']}; {step1['construction']['negative_5hmC']}", step1["construction"]["5hmC"] == "paired donor BS beta minus oxBS beta" and "missing" in step1["construction"]["negative_5hmC"].lower(), "Step1/qc_summary.json")
    add_audit(audit_rows, "5hmC missingness filter", "remove only >20% missing", step1["parameters"]["missing_threshold"], float(step1["parameters"]["missing_threshold"]) == 0.2, "Step1/qc_summary.json")
    add_audit(audit_rows, "Outlier exclusion policy", "No PCA-only automatic exclusion", step2["parameters"]["automatic_donor_exclusion"], step2["parameters"]["automatic_donor_exclusion"] is False, "Step2/step2_summary.json")
    add_audit(audit_rows, "Primary cohort freeze", "77 donors; exclude MSA 30", f"{step2b['cohorts']['primary_clean77']['donors']} donors; {step2b['frozen_decision']['primary_excluded_donors']}", step2b["cohorts"]["primary_clean77"]["donors"] == 77 and step2b["frozen_decision"]["primary_excluded_donors"] == ["MSA 30"], "Step2B/cohort_decision.json")
    add_audit(audit_rows, "Nested CV", "5 outer x 10 repeats; 4 inner", f"{step3['parameters']['outer_folds']} x {step3['parameters']['outer_repeats']}; {step3['parameters']['inner_folds']} inner", step3["parameters"]["outer_folds"] == 5 and step3["parameters"]["outer_repeats"] == 10 and step3["parameters"]["inner_folds"] == 4, "Step3/step3_summary.json")
    add_audit(audit_rows, "Fold construction leakage", "No CpG values used", step3["guardrails"][0], "No CpG values" in step3["guardrails"][0], "Step3/step3_summary.json")
    step5_mode = step5.get(
        "analysis_mode",
        step5.get("parameters", {}).get("analysis_mode", step5.get("pipeline", {}).get("analysis_mode", "")),
    )
    add_audit(audit_rows, "Step 5 research run", "FULL_RESEARCH_RUN", step5_mode, step5_mode == "FULL_RESEARCH_RUN", "Step5/step5_summary.json")
    add_audit(audit_rows, "Stable candidates", "20", step6["counts"]["primary_stable_cpgs"], step6["counts"]["primary_stable_cpgs"] == 20, "Step6/step6_summary.json")
    add_audit(audit_rows, "Sensitivity run", "FULL_RESEARCH_RUN", step7["analysis_mode"], step7["analysis_mode"] == "FULL_RESEARCH_RUN", "Step7/step7_summary.json")
    add_audit(audit_rows, "Candidates supported in both", "15", step7["counts"]["primary_candidates_supported_in_both"], step7["counts"]["primary_candidates_supported_in_both"] == 15, "Step7/step7_summary.json")
    add_audit(audit_rows, "GOmeth probe-number correction", "TRUE", step8_config["prior_prob"], step8_config["prior_prob"] is True, "Step8/RUN_CONFIGURATION.json")
    add_audit(audit_rows, "GOmeth multi-gene correction", "TRUE", step8_config["fractional_multi_gene_counts"], step8_config["fractional_multi_gene_counts"] is True, "Step8/RUN_CONFIGURATION.json")
    audit = pd.DataFrame(audit_rows)
    audit.to_csv(output / "FINAL_METHODS_AUDIT.csv", index=False)
    if not (audit["status"] == "PASS").all():
        failed = audit.loc[audit["status"] != "PASS", "audit_item"].tolist()
        raise RuntimeError("Final methods audit failed: " + ", ".join(failed))

    pathway_sig_count = int(step8.get("total_fdr_significant_terms", 0))
    claims = pd.DataFrame(
        [
            ["Association language", "ALLOWED", "Use 'MSA-associated methylation/hydroxymethylation signatures'."],
            ["Classification performance", "ALLOWED WITH UNCERTAINTY", "Report donor-mean out-of-fold estimates with 95% donor-bootstrap intervals."],
            ["Molecular improvement over covariates", "NOT ESTABLISHED", "Primary paired comparison CI crosses zero and multiplicity-adjusted evidence is not significant."],
            ["Independent validation", "NOT ALLOWED", "Full-78 and Clean-76 overlap the primary cohort; they are sensitivity analyses."],
            ["Diagnostic biomarker", "NOT ALLOWED", "No independent clinical validation or prospective diagnostic evaluation was performed."],
            ["Stable CpGs as statistically significant", "NOT ALLOWED", "Selection frequency is a stability measure, not a p-value."],
            ["Causal mechanism", "NOT ALLOWED", "Observational post-mortem data cannot establish causality."],
            ["Pathway significance", "ALLOWED" if pathway_sig_count > 0 else "NO FDR-SIGNIFICANT PATHWAY", f"GOmeth FDR<0.05 count = {pathway_sig_count}; nominal rankings remain exploratory."],
            ["Generalization", "NOT ESTABLISHED", "Single public dataset; external replication is required."],
        ],
        columns=["claim_topic", "status", "paper_wording_guardrail"],
    )
    claims.to_csv(output / "PAPER_CLAIMS_GUARDRAILS.csv", index=False)

    decision_log = pd.DataFrame(
        [
            [1, "Construct 5mC and 5hmC", "5mC=oxBS; 5hmC=BS-oxBS; negative differences set missing; >20% missing CpGs removed"],
            [2, "Unsupervised QC", "PCA/correlation/missingness flags requested review; flags were not automatic exclusions"],
            ["2B", "Freeze cohorts", "Primary excludes MSA 30; Full-78 and Clean-76 prespecified as sensitivity cohorts"],
            [3, "Freeze confounders and CV", "Age, sex, brain bank, PMI; donor-level repeated nested CV"],
            [4, "Measure non-molecular predictability", "Covariate-only baselines established before molecular comparison"],
            [5, "Molecular benchmark", "All preprocessing, selection and tuning restricted to training folds"],
            [6, "Interpret stable candidates", "Selection stability separated from statistical significance; EPIC manifest audited"],
            [7, "Outlier-decision sensitivity", "Locked molecular-plus-covariate view repeated on overlapping cohorts"],
            [8, "Bias-aware pathways", "GOmeth with separate modality backgrounds; FDR<0.05 reporting rule"],
            [9, "Final freeze", "Verified tables, figures, hashes and paper-claim guardrails generated"],
        ],
        columns=["step", "decision", "locked_rule"],
    )
    decision_log.to_csv(output / "FINAL_ANALYSIS_DECISION_LOG.csv", index=False)

    counts = {
        "cpg_5mc": int(step1["outputs_summary"]["5mC_cpgs"]),
        "cpg_5hmc": int(step1["outputs_summary"]["5hmC_cpgs_retained"]),
        "supported_candidates": int(step7["counts"]["primary_candidates_supported_in_both"]),
    }
    build_workflow_figure(figures_dir / "Figure0_analysis_workflow.png", counts)
    copy_figure(step_dirs["step5"] / "performance_primary.png", figures_dir / "Figure1_primary_molecular_benchmarks.png")
    copy_figure(step_dirs["step7"] / "performance_across_cohorts.png", figures_dir / "Figure2_sensitivity_performance.png")
    copy_figure(step_dirs["step7"] / "primary_candidate_selection_heatmap.png", figures_dir / "Figure3_candidate_robustness.png")
    copy_figure(step_dirs["step7"] / "candidate_effect_concordance.png", figures_dir / "FigureS1_candidate_effect_concordance.png")
    copy_figure(step_dirs["step7"] / "stable_set_overlap.png", figures_dir / "FigureS2_stable_set_overlap.png")
    copy_figure(step_dirs["step6"] / "stable_cpg_annotation_primary.png", figures_dir / "FigureS3_candidate_annotation.png")
    copy_figure(step_dirs["step8"] / "top_enrichment_5mC.png", figures_dir / "FigureS4_5mC_pathway_ranking.png")
    copy_figure(step_dirs["step8"] / "top_enrichment_5hmC.png", figures_dir / "FigureS5_5hmC_pathway_ranking.png")
    copy_figure(step_dirs["step2"] / "qc_5hmC_donor_missingness.png", figures_dir / "FigureS6_5hmC_QC.png")

    primary_perf = perf_out.loc[perf_out["cohort"] == "Primary Clean-77"].set_index("model")
    combined_sensitivity = perf_out.loc[perf_out["model"] == "Combined + covariates"].set_index("cohort")
    paper_numbers = {
        "dataset": {
            "accession": "GSE143157",
            "paired_donors_full": 78,
            "control_full": 37,
            "msa_full": 41,
            "primary_donors": 77,
            "primary_control": 37,
            "primary_msa": 40,
            "primary_excluded": ["MSA 30"],
            "5mC_cpgs": counts["cpg_5mc"],
            "5hmC_cpgs_after_missingness_filter": counts["cpg_5hmc"],
        },
        "primary_roc_auc": {
            label: {
                "estimate": float(primary_perf.loc[label, "roc_auc"]),
                "ci_low": float(primary_perf.loc[label, "roc_auc_ci_low"]),
                "ci_high": float(primary_perf.loc[label, "roc_auc_ci_high"]),
            }
            for label in ["All covariates", "5mC + covariates", "5hmC + covariates", "Combined + covariates"]
        },
        "combined_sensitivity_roc_auc": {
            cohort: {
                "estimate": float(combined_sensitivity.loc[cohort, "roc_auc"]),
                "ci_low": float(combined_sensitivity.loc[cohort, "roc_auc_ci_low"]),
                "ci_high": float(combined_sensitivity.loc[cohort, "roc_auc_ci_high"]),
            }
            for cohort in ["Primary Clean-77", "Sensitivity Full-78", "Sensitivity Clean-76"]
        },
        "combined_vs_all_covariates_primary": {
            "roc_auc_difference": float(target_comparison["roc_auc_difference_a_minus_b"]),
            "bootstrap_ci_low": float(target_comparison["bootstrap_low"]),
            "bootstrap_ci_high": float(target_comparison["bootstrap_high"]),
            "bootstrap_two_sided_p": float(target_comparison["bootstrap_two_sided_p"]),
            "fdr_bh_within_family": float(target_comparison["fdr_bh_within_family"]),
        },
        "candidate_stability": {
            "primary_candidates": int(step7["counts"]["primary_candidates"]),
            "supported_in_both_sensitivity_cohorts": counts["supported_candidates"],
            "all_effect_directions_match_both": bool(candidate_table["effect_direction_matches_both"].all()),
            "jaccard_primary_vs_full78": float(pd.read_csv(step_dirs["step7"] / "stable_set_overlap.csv").query("cohort_a == 'primary_clean77' and cohort_b == 'sensitivity_full78'")["jaccard"].iloc[0]),
            "jaccard_primary_vs_clean76": float(pd.read_csv(step_dirs["step7"] / "stable_set_overlap.csv").query("cohort_a == 'primary_clean77' and cohort_b == 'sensitivity_no_MSA30_Control11_76'")["jaccard"].iloc[0]),
        },
        "pathway_enrichment": {
            "method": "missMethyl::gometh",
            "selected_5mC_cpgs": int(step6["counts"]["gometh_selected_5mC"]),
            "selected_5hmC_cpgs": int(step6["counts"]["gometh_selected_5hmC"]),
            "fdr_significant_terms": pathway_sig_count,
            "interpretation": "exploratory; selected CpGs were ML-stable, not EWAS-significant",
        },
    }
    with (output / "PAPER_NUMBERS.json").open("w", encoding="utf-8") as handle:
        json.dump(paper_numbers, handle, indent=2, ensure_ascii=False)

    methods_results = f"""# GSE143157 locked methods and results numbers

This file is a data-backed drafting aid. It does not replace author review.

## Cohort and preprocessing

The analysis used 78 paired donors (41 MSA and 37 controls) with matched BS and oxBS EPIC profiles. We represented 5mC by the oxBS beta value and estimated 5hmC as BS minus oxBS. Negative 5hmC differences were treated as missing, and CpGs with more than 20% missing 5hmC values were removed. This retained {counts['cpg_5hmc']:,} 5hmC CpGs; the 5mC matrix contained {counts['cpg_5mc']:,} CpGs.

The primary cohort contained 77 donors (40 MSA and 37 controls) after excluding MSA 30 based on convergent technical QC evidence. Full-78 and Clean-76 were prespecified overlapping-cohort sensitivity analyses, not independent external validation cohorts.

## Validation design

Models used five outer folds repeated ten times and four inner folds. Missing-value handling, variance filtering, imputation, scaling, supervised CpG selection, residualization, and hyperparameter tuning were fitted using training donors only. Age, sex, brain bank, and post-mortem interval were prespecified covariates.

## Primary discrimination

In Primary Clean-77, donor-mean out-of-fold ROC-AUC was {primary_perf.loc['All covariates', 'roc_auc']:.3f} (95% donor-bootstrap CI {primary_perf.loc['All covariates', 'roc_auc_ci_low']:.3f}-{primary_perf.loc['All covariates', 'roc_auc_ci_high']:.3f}) for all covariates, {primary_perf.loc['5mC + covariates', 'roc_auc']:.3f} for 5mC plus covariates, {primary_perf.loc['5hmC + covariates', 'roc_auc']:.3f} for 5hmC plus covariates, and {primary_perf.loc['Combined + covariates', 'roc_auc']:.3f} for combined molecular plus covariates. The combined-versus-all-covariates ROC-AUC difference was {float(target_comparison['roc_auc_difference_a_minus_b']):.3f} (95% paired donor-bootstrap CI {float(target_comparison['bootstrap_low']):.3f} to {float(target_comparison['bootstrap_high']):.3f}; within-family BH-adjusted value {float(target_comparison['fdr_bh_within_family']):.3f}). Thus, incremental molecular discrimination beyond covariates was not established.

## Sensitivity and candidate stability

Combined plus covariate ROC-AUC was {combined_sensitivity.loc['Sensitivity Full-78', 'roc_auc']:.3f} in Full-78 and {combined_sensitivity.loc['Sensitivity Clean-76', 'roc_auc']:.3f} in Clean-76. Of 20 primary stable CpG candidates, {counts['supported_candidates']} met the prespecified selection-frequency and effect-direction rule in both sensitivity cohorts. All 20 retained the primary effect direction in both overlapping sensitivity cohorts. These analyses support robustness to the cohort decision but do not constitute external replication.

## Pathway analysis

Bias-aware exploratory enrichment used missMethyl::gometh with EPIC annotation, separate tested-probe backgrounds for 5mC and 5hmC, probe-number prior probabilities, equivalent CpG counts, and fractional counting for multi-gene probes. The input comprised seven manifest-audited stable CpGs per modality. The number of FDR-significant terms was {pathway_sig_count}. Nominally ranked pathways must be described as exploratory rather than confirmed mechanisms.

## Mandatory limitations

- Single public post-mortem dataset.
- Small donor sample and wide uncertainty intervals.
- Overlapping sensitivity cohorts are not external validation.
- Selection frequency is not a p-value.
- Stable CpGs are exploratory associations, not validated diagnostic biomarkers.
- Observational data do not establish causality.
"""
    (output / "METHODS_RESULTS_NUMBERS.md").write_text(methods_results, encoding="utf-8")

    # Record hashes for the critical upstream evidence used in this freeze.
    key_files = {
        "Step1/qc_summary.json": step_dirs["step1"] / "qc_summary.json",
        "Step2/step2_summary.json": step_dirs["step2"] / "step2_summary.json",
        "Step2B/cohort_decision.json": step_dirs["step2b"] / "cohort_decision.json",
        "Step3/step3_summary.json": step_dirs["step3"] / "step3_summary.json",
        "Step3/model_input_policy.json": step_dirs["step3"] / "model_input_policy.json",
        "Step4/step4_summary.json": step_dirs["step4"] / "step4_summary.json",
        "Step4/aggregate_metrics_primary_clean77.csv": step_dirs["step4"] / "aggregate_metrics_primary_clean77.csv",
        "Step5/step5_summary.json": step_dirs["step5"] / "step5_summary.json",
        "Step5/aggregate_metrics_primary_clean77.csv": step_dirs["step5"] / "aggregate_metrics_primary_clean77.csv",
        "Step5/auc_comparisons_primary_clean77.csv": step_dirs["step5"] / "auc_comparisons_primary_clean77.csv",
        "Step5/feature_stability_primary_clean77.csv": step_dirs["step5"] / "feature_stability_primary_clean77.csv",
        "Step6/step6_summary.json": step_dirs["step6"] / "step6_summary.json",
        "Step6/primary_stable_cpgs_primary_clean77.csv": step_dirs["step6"] / "primary_stable_cpgs_primary_clean77.csv",
        "Step7/step7_summary.json": step_dirs["step7"] / "step7_summary.json",
        "Step7/aggregate_metrics_across_cohorts.csv": step_dirs["step7"] / "aggregate_metrics_across_cohorts.csv",
        "Step7/primary_candidate_sensitivity.csv": step_dirs["step7"] / "primary_candidate_sensitivity.csv",
        "Step8/step8_summary.json": step_dirs["step8"] / "step8_summary.json",
        "Step8/enrichment_summary.csv": step_dirs["step8"] / "enrichment_summary.csv",
        "Step8/top_terms_all.csv": step_dirs["step8"] / "top_terms_all.csv",
    }
    manifest_rows = []
    for role, path in key_files.items():
        if not path.is_file():
            raise FileNotFoundError(f"Critical upstream evidence is missing: {path}")
        manifest_rows.append(
            {
                "evidence_role": role,
                "absolute_path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    pd.DataFrame(manifest_rows).to_csv(output / "FINAL_INPUT_FILE_MANIFEST.csv", index=False)

    completed = datetime.now(timezone.utc)
    run_config = {
        "pipeline": "GSE143157 Step 9 final audit and paper package",
        "version": PIPELINE_VERSION,
        "completed_utc": completed.isoformat(),
        "step_directories": {key: str(value) for key, value in step_dirs.items()},
        "output": str(output),
        "primary_cohort": "Primary Clean-77",
        "primary_model_view": "molecular_plus_covariates",
        "candidate_support_rule": "selection frequency >=0.60 and primary effect direction matched in both Full-78 and Clean-76",
        "pathway_significance_rule": "GOmeth Benjamini-Hochberg FDR <0.05",
    }
    with (output / "RUN_CONFIGURATION.json").open("w", encoding="utf-8") as handle:
        json.dump(run_config, handle, indent=2, ensure_ascii=False)

    report = f"""<!doctype html><html><head><meta charset='utf-8'><title>GSE143157 Step 9</title>
<style>body{{font-family:Segoe UI,Arial,sans-serif;max-width:1200px;margin:32px auto;line-height:1.45;color:#202124}}h1,h2{{color:#17365d}}.ok{{background:#e8f5e9;border-left:5px solid #2e7d32;padding:12px}}.warn{{background:#fff8e1;border-left:5px solid #f9a825;padding:12px}}table{{border-collapse:collapse;width:100%;font-size:13px;margin:14px 0}}th,td{{border:1px solid #ccc;padding:6px;text-align:left}}th{{background:#eef3f8}}img{{max-width:100%;height:auto;border:1px solid #ddd}}</style></head><body>
<h1>GSE143157 Step 9: Final verified analysis freeze</h1>
<div class='ok'><b>All locked method checks passed.</b> The tables and figures below are the paper-number source of truth.</div>
<div class='warn'><b>Central interpretation:</b> Molecular signatures showed moderate discrimination and candidate-direction stability, but incremental performance beyond covariates and external generalization were not established.</div>
<h2>Table 1. Cohorts</h2>{cohort_table.to_html(index=False, border=0, escape=True)}
<h2>Table 2. Locked performance</h2>{perf_out[['cohort','model','roc_auc_95ci','pr_auc','balanced_accuracy','delta_roc_auc_vs_all_covariates']].to_html(index=False, border=0, escape=True, float_format=lambda x: f'{x:.3f}')}
<h2>Table 3. Primary CpG candidates</h2>{candidate_table[['source_modality','cpg_id','gene_symbols','stability_tier','primary_model_selection_frequency','minimum_sensitivity_selection_frequency','sensitivity_supported_both','manifest_qc_status','paper_priority']].to_html(index=False, border=0, escape=True, float_format=lambda x: f'{x:.3f}')}
<h2>Pathway status</h2><p>GOmeth FDR-significant terms: <b>{pathway_sig_count}</b>. Selected CpGs were ML-stable candidates, not EWAS-significant loci.</p>
<h2>Workflow</h2><img src='figures/Figure0_analysis_workflow.png' alt='Analysis workflow'>
<h2>Sensitivity performance</h2><img src='figures/Figure2_sensitivity_performance.png' alt='Sensitivity performance'>
<h2>Candidate robustness</h2><img src='figures/Figure3_candidate_robustness.png' alt='Candidate robustness'>
<h2>Final files</h2><ul><li>TABLE1-4: paper-ready source tables</li><li>PAPER_NUMBERS.json: exact machine-readable numbers</li><li>METHODS_RESULTS_NUMBERS.md: data-backed drafting aid</li><li>PAPER_CLAIMS_GUARDRAILS.csv: claims that are and are not supported</li><li>FINAL_FREEZE.json: cryptographic hashes of this package</li></ul>
</body></html>"""
    (output / "STEP9_REPORT.html").write_text(report, encoding="utf-8")

    # Freeze every generated output except the freeze record and success marker themselves.
    freeze_files = []
    for path in sorted(output.rglob("*")):
        if not path.is_file() or path.name in {"FINAL_FREEZE.json", "SUCCESS.txt"}:
            continue
        freeze_files.append(
            {
                "relative_path": path.relative_to(output).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    freeze_seed = "\n".join(f"{item['relative_path']}:{item['sha256']}" for item in freeze_files).encode("utf-8")
    freeze = {
        "freeze_id": hashlib.sha256(freeze_seed).hexdigest(),
        "created_utc": completed.isoformat(),
        "pipeline_version": PIPELINE_VERSION,
        "files": freeze_files,
        "instructions": "Do not edit frozen files in place. Create a versioned revision and regenerate FINAL_FREEZE.json.",
    }
    with (output / "FINAL_FREEZE.json").open("w", encoding="utf-8") as handle:
        json.dump(freeze, handle, indent=2, ensure_ascii=False)

    (output / "SUCCESS.txt").write_text(
        "\n".join(
            [
                "SUCCESS",
                f"Completed UTC: {completed.isoformat()}",
                "Analysis mode: FINAL_VERIFIED_FREEZE",
                f"Freeze ID: {freeze['freeze_id']}",
                f"Primary candidates supported in both sensitivity cohorts: {counts['supported_candidates']}/20",
                f"GOmeth FDR-significant pathways: {pathway_sig_count}",
                "Interpretation: internal donor-level analysis with overlapping-cohort sensitivity; not external validation",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print("STEP 9 BUILD COMPLETED")
    print(f"Freeze ID: {freeze['freeze_id']}")


if __name__ == "__main__":
    main()
