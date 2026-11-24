#!/usr/bin/env python3
"""Step 4: covariate-only leakage-safe nested-CV benchmarks."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import sklearn

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, roc_curve

from step4_core import (
    METRIC_NAMES,
    MODEL_SPECS,
    aggregate_metrics,
    build_pipeline,
    coefficient_stability,
    donor_mean_predictions,
    encode_condition,
    extract_coefficients,
    pairwise_auc_comparisons,
    positive_probability,
    repeat_metrics,
    select_hyperparameters,
)


VERSION = "1.0.0"
COHORTS = {
    "primary_clean77": 77,
    "sensitivity_full78": 78,
    "sensitivity_no_MSA30_Control11_76": 76,
}
MODEL_COLORS = {
    "age_only": "#0072B2",
    "no_age_covariates": "#009E73",
    "all_covariates": "#D55E00",
}
MODEL_LABELS = {
    "age_only": "Age only",
    "no_age_covariates": "Sex + bank + PMI",
    "all_covariates": "All covariates",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="GSE143157 Step 4 covariate-only nested-CV benchmark",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--step2b-dir", required=True, type=Path)
    parser.add_argument("--step3-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--c-grid", default="0.01,0.1,1,10,100")
    parser.add_argument("--bootstrap-iterations", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=143157)
    return parser.parse_args()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_c_grid(text: str) -> list[float]:
    try:
        values = sorted({float(item.strip()) for item in text.split(",")})
    except ValueError as exc:
        raise ValueError("--c-grid must be comma-separated positive numbers.") from exc
    if not values or any(value <= 0 for value in values):
        raise ValueError("Every C value must be positive.")
    return values


def validate_inputs(
    step2b: Path, step3: Path
) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame], dict[str, object], dict[str, Path]]:
    required = {
        "step2b_success": step2b / "SUCCESS.txt",
        "cohort_decision": step2b / "cohort_decision.json",
        "step3_success": step3 / "SUCCESS.txt",
        "step3_summary": step3 / "step3_summary.json",
        "model_policy": step3 / "model_input_policy.json",
    }
    for cohort in COHORTS:
        required[f"outer_{cohort}"] = step3 / f"outer_folds_{cohort}.csv"
        required[f"inner_{cohort}"] = step3 / f"inner_folds_{cohort}.csv"
    missing = [str(path) for path in required.values() if not path.is_file()]
    if missing:
        raise ValueError(f"Required Step 2B/Step 3 files are missing: {missing}")
    summary = json.loads(required["step3_summary"].read_text(encoding="utf-8"))
    policy = json.loads(required["model_policy"].read_text(encoding="utf-8"))
    if policy["primary_cohort"] != "primary_clean77":
        raise ValueError("Step 3 policy primary cohort is not Clean-77.")
    if set(policy["same_outer_and_inner_donor_folds_for_modalities"]) != {
        "5mC",
        "5hmC",
        "combined",
    }:
        raise ValueError("Step 3 modality-fold policy has changed.")
    if summary["parameters"]["outer_folds"] != 5 or summary["parameters"]["outer_repeats"] != 10:
        raise ValueError("Expected the frozen 5-fold x 10-repeat outer CV.")
    if summary["parameters"]["inner_folds"] != 4:
        raise ValueError("Expected the frozen 4-fold inner CV.")

    outer_tables: dict[str, pd.DataFrame] = {}
    inner_tables: dict[str, pd.DataFrame] = {}
    for cohort, donors in COHORTS.items():
        outer = pd.read_csv(required[f"outer_{cohort}"])
        inner = pd.read_csv(required[f"inner_{cohort}"])
        if len(outer) != donors * 10 or outer["donor_id"].nunique() != donors:
            raise ValueError(f"{cohort}: outer-fold table has the wrong size.")
        outer_tables[cohort] = outer
        inner_tables[cohort] = inner
    return outer_tables, inner_tables, {"summary": summary, "policy": policy}, required


def run_cohort(
    cohort: str,
    outer: pd.DataFrame,
    inner: pd.DataFrame,
    c_grid: list[float],
    bootstrap_iterations: int,
    seed: int,
) -> dict[str, pd.DataFrame]:
    metadata = outer.loc[outer["repeat"] == 1].copy()
    metadata = metadata.sort_values("donor_index", kind="stable").drop(columns=["repeat", "outer_test_fold"])
    metadata = metadata.set_index("donor_id", drop=False)
    prediction_rows: list[dict[str, object]] = []
    choice_rows: list[dict[str, object]] = []
    tuning_frames: list[pd.DataFrame] = []
    coefficient_frames: list[pd.DataFrame] = []
    class_weights: list[str | None] = [None, "balanced"]

    for model_index, (model_name, spec) in enumerate(MODEL_SPECS.items()):
        print(f"[{cohort}] model: {model_name}", flush=True)
        for repeat in range(1, 11):
            print(f"  repeat {repeat}/10", flush=True)
            repeat_outer = outer.loc[outer["repeat"] == repeat]
            for outer_fold in range(1, 6):
                test_ids = repeat_outer.loc[
                    repeat_outer["outer_test_fold"] == outer_fold, "donor_id"
                ].astype(str).tolist()
                train_ids = repeat_outer.loc[
                    repeat_outer["outer_test_fold"] != outer_fold, "donor_id"
                ].astype(str).tolist()
                outer_train = metadata.loc[train_ids].reset_index(drop=True)
                outer_test = metadata.loc[test_ids].reset_index(drop=True)
                inner_current = inner.loc[
                    (inner["outer_repeat"] == repeat) & (inner["outer_fold"] == outer_fold)
                ].copy()
                fold_seed = seed + model_index * 10_000_019 + repeat * 100_003 + outer_fold
                best, tuning = select_hyperparameters(
                    outer_train,
                    inner_current,
                    spec,
                    c_grid,
                    class_weights,
                    fold_seed,
                )
                tuning.insert(0, "cohort", cohort)
                tuning.insert(1, "repeat", repeat)
                tuning.insert(2, "outer_fold", outer_fold)
                tuning_frames.append(tuning)
                choice_rows.append(
                    {
                        "cohort": cohort,
                        "model": model_name,
                        "repeat": repeat,
                        "outer_fold": outer_fold,
                        "outer_train_donors": len(outer_train),
                        "outer_test_donors": len(outer_test),
                        "selected_C": best["C"],
                        "selected_class_weight": "none" if best["class_weight"] is None else best["class_weight"],
                        "selected_mean_inner_roc_auc": best["mean_inner_roc_auc"],
                        "selected_sd_inner_roc_auc": best["sd_inner_roc_auc"],
                    }
                )
                y_train = encode_condition(outer_train["condition"])
                y_test = encode_condition(outer_test["condition"])
                fitted = build_pipeline(
                    spec,
                    float(best["C"]),
                    best["class_weight"],
                    fold_seed + 99991,
                )
                fitted.fit(outer_train, y_train)
                probability = positive_probability(fitted, outer_test)
                for row_index, donor in outer_test.reset_index(drop=True).iterrows():
                    prediction_rows.append(
                        {
                            "cohort": cohort,
                            "model": model_name,
                            "repeat": repeat,
                            "outer_fold": outer_fold,
                            "donor_id": donor["donor_id"],
                            "condition": donor["condition"],
                            "y_true": int(y_test[row_index]),
                            "probability_msa": float(probability[row_index]),
                            "predicted_label": "MSA" if probability[row_index] >= 0.5 else "Control",
                            "selected_C": best["C"],
                            "selected_class_weight": "none" if best["class_weight"] is None else best["class_weight"],
                        }
                    )
                coefficient_frames.append(
                    extract_coefficients(
                        fitted,
                        model_name,
                        repeat,
                        outer_fold,
                        float(best["C"]),
                        best["class_weight"],
                    )
                )

    predictions = pd.DataFrame(prediction_rows)
    choices = pd.DataFrame(choice_rows)
    tuning_scores = pd.concat(tuning_frames, ignore_index=True)
    coefficients = pd.concat(coefficient_frames, ignore_index=True)
    repeats = repeat_metrics(predictions)
    donor_means = donor_mean_predictions(predictions)
    aggregate = aggregate_metrics(
        repeats, donor_means, bootstrap_iterations, seed + 30_000_091
    )
    comparisons = pairwise_auc_comparisons(
        donor_means, bootstrap_iterations, seed + 40_000_093
    )
    stability = coefficient_stability(coefficients)
    return {
        "predictions": predictions,
        "hyperparameter_choices": choices,
        "tuning_scores": tuning_scores,
        "outer_coefficients": coefficients,
        "coefficient_stability": stability,
        "repeat_metrics": repeats,
        "donor_mean_predictions": donor_means,
        "aggregate_metrics": aggregate,
        "pairwise_auc_comparisons": comparisons,
    }


def plot_roc(donor_means: pd.DataFrame, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for model, group in donor_means.groupby("model", sort=False):
        fpr, tpr, _ = roc_curve(group["y_true"], group["probability_msa_mean"])
        from sklearn.metrics import roc_auc_score

        auc_value = roc_auc_score(group["y_true"], group["probability_msa_mean"])
        ax.plot(fpr, tpr, linewidth=2, color=MODEL_COLORS[model], label=f"{MODEL_LABELS[model]} (AUC={auc_value:.3f})")
    ax.plot([0, 1], [0, 1], "--", color="#777777", label="Chance")
    ax.set_xlabel("False-positive rate")
    ax.set_ylabel("True-positive rate")
    ax.set_title("Primary Clean-77 covariate-only ROC")
    ax.legend(frameon=False)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_pr(donor_means: pd.DataFrame, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    prevalence = float(donor_means.drop_duplicates("donor_id")["y_true"].mean())
    for model, group in donor_means.groupby("model", sort=False):
        precision, recall, _ = precision_recall_curve(
            group["y_true"], group["probability_msa_mean"]
        )
        from sklearn.metrics import average_precision_score

        ap = average_precision_score(group["y_true"], group["probability_msa_mean"])
        ax.plot(recall, precision, linewidth=2, color=MODEL_COLORS[model], label=f"{MODEL_LABELS[model]} (AP={ap:.3f})")
    ax.axhline(prevalence, linestyle="--", color="#777777", label=f"Prevalence={prevalence:.3f}")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Primary Clean-77 covariate-only precision–recall")
    ax.legend(frameon=False)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_metrics(aggregate: pd.DataFrame, output: Path) -> None:
    metrics = ["roc_auc", "pr_auc", "balanced_accuracy"]
    x = np.arange(len(metrics))
    width = 0.24
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    for model_index, model in enumerate(MODEL_SPECS):
        row = aggregate.loc[aggregate["model"] == model].iloc[0]
        values = [row[f"{metric}_donor_mean_oof"] for metric in metrics]
        lower = [row[f"{metric}_bootstrap_low"] for metric in metrics]
        upper = [row[f"{metric}_bootstrap_high"] for metric in metrics]
        errors = np.asarray([np.asarray(values) - np.asarray(lower), np.asarray(upper) - np.asarray(values)])
        ax.bar(
            x + (model_index - 1) * width,
            values,
            width,
            yerr=errors,
            capsize=3,
            color=MODEL_COLORS[model],
            label=MODEL_LABELS[model],
            alpha=0.9,
        )
    ax.axhline(0.5, linestyle="--", color="#777777", linewidth=1)
    ax.set_xticks(x, ["ROC-AUC", "PR-AUC", "Balanced accuracy"])
    ax.set_ylim(0, 1.03)
    ax.set_ylabel("Out-of-fold performance")
    ax.set_title("Primary Clean-77 covariate benchmarks (95% donor bootstrap CI)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_coefficients(coefficients: pd.DataFrame, output: Path) -> None:
    subset = coefficients.loc[coefficients["model"] == "all_covariates"].copy()
    features = sorted(subset["feature"].unique())
    values = [subset.loc[subset["feature"] == feature, "coefficient"].to_numpy() for feature in features]
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.boxplot(values, showfliers=False)
    ax.set_xticks(range(1, len(features) + 1))
    ax.set_xticklabels(features, rotation=30, ha="right")
    ax.axhline(0, color="#777777", linestyle="--", linewidth=1)
    ax.set_ylabel("Regularized logistic coefficient")
    ax.set_title("All-covariate coefficient stability across 50 outer fits")
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_report(output: Path, aggregate: pd.DataFrame, comparisons: pd.DataFrame, stability: pd.DataFrame, parameters: dict[str, object]) -> None:
    display_metrics = aggregate[
        [
            "model",
            "roc_auc_donor_mean_oof",
            "roc_auc_bootstrap_low",
            "roc_auc_bootstrap_high",
            "pr_auc_donor_mean_oof",
            "balanced_accuracy_donor_mean_oof",
            "mcc_donor_mean_oof",
            "brier_score_donor_mean_oof",
        ]
    ].copy()
    metrics_html = display_metrics.to_html(index=False, border=0, classes="data-table", float_format=lambda value: f"{value:.4f}")
    comparisons_html = comparisons.to_html(index=False, border=0, classes="data-table", float_format=lambda value: f"{value:.4f}")
    stability_html = stability.loc[stability["model"] == "all_covariates", ["feature", "coefficient_mean", "coefficient_sd", "positive_fraction", "odds_ratio_median"]].to_html(index=False, border=0, classes="data-table", float_format=lambda value: f"{value:.4f}")
    document = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>GSE143157 Step 4 report</title><style>
body{{font-family:Arial,sans-serif;max-width:1080px;margin:28px auto;padding:0 18px;line-height:1.5;color:#222}}h1,h2{{color:#17365d}}
.notice{{background:#fff4ce;border-left:5px solid #d6a800;padding:12px 16px}}.good{{background:#eaf6ec;border-left:5px solid #2e7d32;padding:12px 16px}}
.data-table{{border-collapse:collapse;width:100%;font-size:.87rem}}.data-table th,.data-table td{{padding:7px;border-bottom:1px solid #ddd;text-align:left}}
img{{max-width:100%;height:auto}}code{{background:#f4f4f4;padding:2px 4px}}</style></head><body>
<h1>GSE143157 Step 4: covariate-only nested-CV benchmark</h1>
<p class="good"><strong>Leakage check:</strong> Every preprocessing fit and hyperparameter decision used inner/outer training donors only. Each donor received one outer-test prediction per repeat.</p>
<p class="notice"><strong>Interpretation:</strong> These are confounding benchmarks, not molecular biomarkers. A later methylation model must improve over the all-covariate benchmark to claim added molecular predictive value.</p>
<h2>Primary Clean-77 aggregate out-of-fold metrics</h2>{metrics_html}
<h2>Paired donor-bootstrap ROC-AUC comparisons</h2>{comparisons_html}
<h2>ROC curve</h2><img src="roc_primary.png" alt="Primary ROC curves">
<h2>Precision–recall curve</h2><img src="pr_primary.png" alt="Primary precision recall curves">
<h2>Metric comparison</h2><img src="metrics_primary.png" alt="Primary performance metrics">
<h2>All-covariate coefficient stability</h2>{stability_html}<img src="coefficients_primary.png" alt="Coefficient stability">
<h2>Frozen design</h2><ul><li>Outer CV: 5 folds × 10 repeats.</li><li>Inner CV: 4 folds.</li><li>Inner selection metric: ROC-AUC.</li><li>C grid: {html.escape(str(parameters['c_grid']))}; class-weight options: none and balanced.</li><li>Prediction threshold: 0.5, not optimized on outer test data.</li></ul>
</body></html>"""
    output.write_text(document, encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    started = time.time()
    c_grid = parse_c_grid(args.c_grid)
    if args.bootstrap_iterations < 500:
        raise ValueError("--bootstrap-iterations must be at least 500.")
    step2b = args.step2b_dir.resolve()
    step3 = args.step3_dir.resolve()
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Output directory is not empty; nothing will be overwritten: {output}")
    outer_tables, inner_tables, context, required = validate_inputs(step2b, step3)
    output.mkdir(parents=True, exist_ok=True)

    all_results: dict[str, dict[str, pd.DataFrame]] = {}
    for cohort_index, cohort in enumerate(COHORTS):
        all_results[cohort] = run_cohort(
            cohort,
            outer_tables[cohort],
            inner_tables[cohort],
            c_grid,
            args.bootstrap_iterations,
            args.seed + cohort_index * 100_000_007,
        )
        for result_name, frame in all_results[cohort].items():
            frame.to_csv(output / f"{result_name}_{cohort}.csv", index=False)

    primary = all_results["primary_clean77"]
    plot_roc(primary["donor_mean_predictions"], output / "roc_primary.png")
    plot_pr(primary["donor_mean_predictions"], output / "pr_primary.png")
    plot_metrics(primary["aggregate_metrics"], output / "metrics_primary.png")
    plot_coefficients(primary["outer_coefficients"], output / "coefficients_primary.png")

    parameters = {
        "outer_folds": 5,
        "outer_repeats": 10,
        "inner_folds": 4,
        "inner_selection_metric": "ROC-AUC",
        "c_grid": c_grid,
        "class_weight_grid": ["none", "balanced"],
        "classification_threshold": 0.5,
        "bootstrap_iterations": args.bootstrap_iterations,
        "bootstrap_unit": "donor; stratified by diagnosis",
        "random_seed": args.seed,
    }
    write_report(
        output / "STEP4_REPORT.html",
        primary["aggregate_metrics"],
        primary["pairwise_auc_comparisons"],
        primary["coefficient_stability"],
        parameters,
    )
    summary = {
        "pipeline": {
            "name": "GSE143157 Step 4 covariate-only nested-CV benchmark",
            "version": VERSION,
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": time.time() - started,
            "scikit_learn_version": sklearn.__version__,
        },
        "inputs": {
            "step2b_dir": str(step2b),
            "step3_dir": str(step3),
            "step3_summary_sha256": file_sha256(required["step3_summary"]),
            "model_input_policy_sha256": file_sha256(required["model_policy"]),
        },
        "parameters": parameters,
        "models": {
            name: {
                "continuous": list(spec.continuous),
                "categorical": list(spec.categorical),
                "description": spec.description,
            }
            for name, spec in MODEL_SPECS.items()
        },
        "cohorts": {
            cohort: {
                "donors": COHORTS[cohort],
                "prediction_rows": len(results["predictions"]),
                "primary_metrics": results["aggregate_metrics"].to_dict(orient="records"),
                "pairwise_auc_comparisons": results["pairwise_auc_comparisons"].to_dict(orient="records"),
            }
            for cohort, results in all_results.items()
        },
        "interpretation_guardrails": [
            "This step measures prediction from baseline covariates only; it is not a molecular biomarker result.",
            "Molecular models must use the identical outer and inner donor folds.",
            "Added molecular value requires comparison against the all_covariates benchmark.",
            "No outer-test prediction was used for preprocessing, hyperparameter selection, or threshold tuning.",
            "Report all repeats and uncertainty; do not select the best repeat.",
        ],
    }
    with (output / "step4_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False, allow_nan=False)
    (output / "SUCCESS.txt").write_text(
        "SUCCESS\n"
        f"Completed UTC: {summary['pipeline']['completed_utc']}\n"
        "Benchmark: covariate-only nested CV\n"
        "Models: age_only, no_age_covariates, all_covariates\n"
        "Cohorts: primary Clean-77, Full-78 sensitivity, optional 76 stress-test\n"
        "Every donor received one outer-test prediction per repeat and model\n",
        encoding="utf-8",
    )
    print("\nSTEP 4 SUCCESS")
    print(f"Open report: {output / 'STEP4_REPORT.html'}")


def main() -> int:
    try:
        run(parse_args())
    except (ValueError, OSError, KeyError, json.JSONDecodeError, pd.errors.ParserError) as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
