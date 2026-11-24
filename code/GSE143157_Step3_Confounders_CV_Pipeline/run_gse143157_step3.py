#!/usr/bin/env python3
"""Step 3: confounder audit and frozen leakage-safe nested-CV folds."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from step3_core import build_nested_cv, covariate_audit, validate_metadata


VERSION = "1.0.0"
COHORTS = {
    "primary_clean77": {
        "file": "metadata_donors_primary_clean77.csv",
        "donors": 77,
        "conditions": {"MSA": 40, "Control": 37},
    },
    "sensitivity_full78": {
        "file": "metadata_donors_sensitivity_full78.csv",
        "donors": 78,
        "conditions": {"MSA": 41, "Control": 37},
    },
    "sensitivity_no_MSA30_Control11_76": {
        "file": "metadata_donors_sensitivity_no_MSA30_Control11_76.csv",
        "donors": 76,
        "conditions": {"MSA": 40, "Control": 36},
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="GSE143157 Step 3: covariate audit and deterministic nested CV",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--step2b-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--outer-folds", type=int, default=5)
    parser.add_argument("--outer-repeats", type=int, default=10)
    parser.add_argument("--inner-folds", type=int, default=4)
    parser.add_argument("--fold-candidates", type=int, default=200)
    parser.add_argument("--permutations", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=143157)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_args(args: argparse.Namespace) -> None:
    if args.outer_folds < 3 or args.outer_folds > 10:
        raise ValueError("--outer-folds must be between 3 and 10.")
    if args.outer_repeats < 2:
        raise ValueError("--outer-repeats must be at least 2.")
    if args.inner_folds < 3 or args.inner_folds >= args.outer_folds + 1:
        raise ValueError("--inner-folds must be at least 3 and suitable for outer training sets.")
    if args.fold_candidates < 50:
        raise ValueError("--fold-candidates must be at least 50.")
    if args.permutations < 1000:
        raise ValueError("--permutations must be at least 1000.")


def load_inputs(root: Path) -> tuple[dict[str, pd.DataFrame], dict[str, object], dict[str, Path]]:
    required_names = [
        "SUCCESS.txt",
        "cohort_decision.json",
        "cohort_manifest.csv",
        "GSE143157_primary_clean77.h5",
    ] + [spec["file"] for spec in COHORTS.values()]
    paths = {name: root / name for name in required_names}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise ValueError(f"Step 2B files are missing: {missing}")
    with paths["GSE143157_primary_clean77.h5"].open("rb") as handle:
        signature = handle.read(8)
    if signature != b"\x89HDF\r\n\x1a\n":
        raise ValueError("Primary Clean-77 file does not have a valid HDF5 signature.")
    decision = json.loads(paths["cohort_decision.json"].read_text(encoding="utf-8"))
    if decision["frozen_decision"]["primary_cohort"] != "primary_clean77":
        raise ValueError("Step 2B primary cohort is not the frozen Clean-77 cohort.")
    if decision["primary_hdf5"]["matrix_stats"]["5mC"]["shape"][1] != 77:
        raise ValueError("Step 2B primary 5mC matrix does not report 77 donors.")
    if decision["primary_hdf5"]["matrix_stats"]["5hmC"]["shape"][1] != 77:
        raise ValueError("Step 2B primary 5hmC matrix does not report 77 donors.")

    metadata: dict[str, pd.DataFrame] = {}
    for cohort, spec in COHORTS.items():
        frame = pd.read_csv(paths[spec["file"]])
        metadata[cohort] = validate_metadata(
            frame, spec["donors"], spec["conditions"], cohort
        )
    if set(metadata["sensitivity_full78"]["donor_id"]) - set(metadata["primary_clean77"]["donor_id"]) != {"MSA 30"}:
        raise ValueError("Primary-vs-Full donor difference is not exactly MSA 30.")
    if set(metadata["sensitivity_full78"]["donor_id"]) - set(metadata["sensitivity_no_MSA30_Control11_76"]["donor_id"]) != {"MSA 30", "Control 11"}:
        raise ValueError("The optional 76-donor sensitivity definition is wrong.")
    return metadata, decision, paths


def preprocessing_policy(args: argparse.Namespace) -> dict[str, object]:
    return {
        "policy_version": "1.0.0",
        "primary_cohort": "primary_clean77",
        "same_outer_and_inner_donor_folds_for_modalities": ["5mC", "5hmC", "combined"],
        "covariates": ["age", "sex", "brain_bank", "pmi_hours"],
        "required_model_views": [
            "covariate_only_baseline",
            "molecular_only",
            "molecular_plus_covariates",
            "training_fold_residualized_molecular_sensitivity",
        ],
        "nested_cv": {
            "outer_folds": args.outer_folds,
            "outer_repeats": args.outer_repeats,
            "inner_folds": args.inner_folds,
            "purpose_outer": "unbiased performance estimation only",
            "purpose_inner": "hyperparameter and feature-count selection only",
            "split_optimization_inputs": "condition labels and baseline covariates only; no CpG values",
        },
        "inside_each_outer_training_fold_only": [
            "5hmC missingness filtering using training donors only",
            "per-CpG median imputation fitted on training donors only",
            "zero/near-zero variance filtering",
            "unsupervised variance prefiltering",
            "all supervised feature selection",
            "feature scaling",
            "covariate encoder fitting",
            "molecular residualization sensitivity model fitting",
            "hyperparameter tuning through the predefined inner folds",
        ],
        "application_to_outer_test": [
            "reuse training-fold CpG list in the same order",
            "reuse training-fold medians",
            "reuse training-fold scaling parameters",
            "reuse training-fold covariate encoding",
            "reuse training-fold residualization coefficients when that sensitivity model is evaluated",
        ],
        "prohibited": [
            "global imputation before cross-validation",
            "global scaling before cross-validation",
            "using Step 2 top-variable CpGs as supervised-model features",
            "selecting CpGs using all diagnosis labels before cross-validation",
            "tuning on outer-test performance",
            "reporting only the best repeat or fold",
            "using different donor folds for 5mC, 5hmC, and combined views",
        ],
        "molecular_scale": {
            "5mC": "Use beta values for the primary predictive pipeline; any M-value analysis is a separately documented sensitivity analysis.",
            "5hmC": "Use non-negative BS-oxBS difference values; do not apply a beta-value logit transform.",
        },
        "reporting": [
            "aggregate predictions from every outer test fold and report all repeats",
            "report balanced accuracy, ROC-AUC, PR-AUC, sensitivity, specificity, and MCC with uncertainty",
            "compare against the covariate-only baseline",
            "repeat final locked workflow on Full-78",
            "use the 76-donor cohort only as an optional stress-test",
        ],
    }


def save_plot_covariates(metadata: pd.DataFrame, output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"Control": "#0072B2", "MSA": "#D55E00"}
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.2))
    rng = np.random.default_rng(143157)
    for ax, covariate, title in [
        (axes[0, 0], "age", "Age"),
        (axes[0, 1], "pmi_hours", "Post-mortem interval"),
    ]:
        groups = []
        for position, condition in enumerate(["Control", "MSA"], start=1):
            values = pd.to_numeric(
                metadata.loc[metadata["condition"] == condition, covariate], errors="coerce"
            ).dropna()
            groups.append(values.to_numpy())
            jitter = rng.normal(0, 0.045, len(values))
            ax.scatter(
                np.full(len(values), position) + jitter,
                values,
                color=colors[condition],
                alpha=0.65,
                s=28,
            )
        ax.boxplot(groups, positions=[1, 2], widths=0.45, showfliers=False)
        ax.set_xticks([1, 2], ["Control", "MSA"])
        ax.set_title(title)
    for ax, covariate, title in [
        (axes[1, 0], "sex", "Sex"),
        (axes[1, 1], "brain_bank", "Brain bank"),
    ]:
        table = pd.crosstab(metadata["condition"], metadata[covariate], normalize="index")
        table = table.reindex(["Control", "MSA"])
        bottom = np.zeros(2)
        palette = ["#56B4E9", "#E69F00", "#009E73", "#CC79A7"]
        for index, category in enumerate(table.columns):
            values = table[category].to_numpy()
            label = "MRC London" if "MRC" in str(category) else "Bispebjerg" if "Bispebjerg" in str(category) else str(category)
            ax.bar([0, 1], values, bottom=bottom, label=label, color=palette[index % len(palette)])
            bottom += values
        ax.set_xticks([0, 1], ["Control", "MSA"])
        ax.set_ylim(0, 1)
        ax.set_ylabel("Proportion")
        ax.set_title(title)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Primary Clean-77 covariate balance")
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def save_plot_fold_balance(balance: pd.DataFrame, output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    measures = [
        ("age_test_minus_train", "Age mean difference"),
        ("pmi_hours_test_minus_train", "PMI mean difference"),
    ]
    female_test = next(column for column in balance.columns if column.startswith("sex_female_test_fraction"))
    female_train = female_test.replace("_test_", "_train_")
    balance = balance.copy()
    balance["female_fraction_difference"] = balance[female_test] - balance[female_train]
    mrc_test = next(column for column in balance.columns if column.startswith("brain_bank_mrc_") and column.endswith("_test_fraction"))
    mrc_train = mrc_test.replace("_test_", "_train_")
    balance["mrc_fraction_difference"] = balance[mrc_test] - balance[mrc_train]
    measures += [
        ("female_fraction_difference", "Female fraction difference"),
        ("mrc_fraction_difference", "MRC-bank fraction difference"),
    ]
    fig, axes = plt.subplots(4, 1, figsize=(11.0, 7.8), sharex=True)
    for ax, (column, title) in zip(axes, measures):
        matrix = balance.pivot(index="outer_fold", columns="repeat", values=column).to_numpy()
        limit = max(abs(float(np.nanmin(matrix))), abs(float(np.nanmax(matrix))), 1e-6)
        image = ax.imshow(matrix, cmap="coolwarm", vmin=-limit, vmax=limit, aspect="auto")
        ax.set_yticks(range(matrix.shape[0]), range(1, matrix.shape[0] + 1))
        ax.set_ylabel("Fold")
        ax.set_title(title, fontsize=10)
        fig.colorbar(image, ax=ax, fraction=0.02, pad=0.02)
    axes[-1].set_xticks(range(10), range(1, 11))
    axes[-1].set_xlabel("Outer repeat")
    fig.suptitle("Primary Clean-77 outer-test versus training covariate balance")
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_report(path: Path, summary: dict[str, object], audit: pd.DataFrame) -> None:
    table = audit[
        [
            "covariate",
            "control_summary",
            "msa_summary",
            "effect_name",
            "effect_size",
            "effect_magnitude",
            "permutation_p_two_sided",
            "fdr_bh",
        ]
    ].to_html(index=False, border=0, classes="data-table", float_format=lambda value: f"{value:.4g}")
    document = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>GSE143157 Step 3 report</title><style>
body{{font-family:Arial,sans-serif;max-width:1080px;margin:28px auto;padding:0 18px;line-height:1.5;color:#222}}h1,h2{{color:#17365d}}
.notice{{background:#fff4ce;border-left:5px solid #d6a800;padding:12px 16px}}.good{{background:#eaf6ec;border-left:5px solid #2e7d32;padding:12px 16px}}
.data-table{{border-collapse:collapse;width:100%;font-size:.88rem}}.data-table th,.data-table td{{padding:7px;border-bottom:1px solid #ddd;text-align:left}}
img{{max-width:100%;height:auto}}code{{background:#f4f4f4;padding:2px 4px}}</style></head><body>
<h1>GSE143157 Step 3: confounders and leakage-safe nested CV</h1>
<p class="good"><strong>Frozen primary cohort:</strong> Clean-77. The exact same donor folds must be used for 5mC, 5hmC, and combined models.</p>
<p class="notice"><strong>Important:</strong> Covariate balance does not make global preprocessing safe. Imputation, variance filtering, feature selection, scaling, residualization, and tuning must still be fitted inside the outer-training data only.</p>
<h2>Primary covariate audit</h2><p>Two-sided permutation tests used {summary['parameters']['permutations']:,} permutations; four p-values were adjusted by Benjamini–Hochberg FDR. Effect sizes, not p-values alone, determine interpretation.</p>{table}
<h2>Frozen CV design</h2><ul><li>{summary['parameters']['outer_folds']}-fold outer CV repeated {summary['parameters']['outer_repeats']} times.</li><li>{summary['parameters']['inner_folds']}-fold inner CV inside every outer-training set.</li><li>Partitions were diagnosis-stratified and selected for age, sex, brain-bank, and PMI balance without reading CpG values.</li><li>Covariate-only, molecular-only, and molecular-plus-covariate models are required.</li></ul>
<h2>Primary covariate distribution</h2><img src="covariate_balance_primary.png" alt="Primary covariate balance">
<h2>Outer-fold balance</h2><img src="outer_fold_balance_primary.png" alt="Outer fold balance">
<h2>Leakage rule</h2><p>Do not use the Step 2 top-variable CpG list for supervised modeling. Every molecular preprocessing decision must be learned again inside each training fold.</p>
</body></html>"""
    path.write_text(document, encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    validate_args(args)
    started = time.time()
    root = args.step2b_dir.resolve()
    output = args.output.resolve()
    metadata, step2b_decision, paths = load_inputs(root)
    print("Validated Step 2B cohorts and primary HDF5 signature.")

    audits: dict[str, pd.DataFrame] = {}
    outer_tables: dict[str, pd.DataFrame] = {}
    inner_tables: dict[str, pd.DataFrame] = {}
    balance_tables: dict[str, pd.DataFrame] = {}
    cv_details: dict[str, object] = {}
    for index, (cohort, frame) in enumerate(metadata.items()):
        print(f"[{cohort}] covariate audit...", flush=True)
        audits[cohort] = covariate_audit(
            frame, args.permutations, args.seed + index * 100_003
        )
        print(f"[{cohort}] optimized repeated nested-CV folds...", flush=True)
        outer, inner, balance, details = build_nested_cv(
            frame,
            args.outer_repeats,
            args.outer_folds,
            args.inner_folds,
            args.fold_candidates,
            args.seed + index * 10_000_019,
        )
        outer_tables[cohort] = outer
        inner_tables[cohort] = inner
        balance_tables[cohort] = balance
        cv_details[cohort] = details

    print("Primary covariate audit:")
    print(
        audits["primary_clean77"][
            ["covariate", "effect_name", "effect_size", "effect_magnitude", "permutation_p_two_sided", "fdr_bh"]
        ].to_string(index=False)
    )
    if args.dry_run:
        print("DRY RUN SUCCESS: no output files were written.")
        return
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Output directory is not empty; nothing will be overwritten: {output}")
    output.mkdir(parents=True, exist_ok=True)

    for cohort in COHORTS:
        audits[cohort].to_csv(output / f"covariate_audit_{cohort}.csv", index=False)
        outer_tables[cohort].to_csv(output / f"outer_folds_{cohort}.csv", index=False)
        inner_tables[cohort].to_csv(output / f"inner_folds_{cohort}.csv", index=False)
        balance_tables[cohort].to_csv(output / f"outer_fold_balance_{cohort}.csv", index=False)

    policy = preprocessing_policy(args)
    with (output / "model_input_policy.json").open("w", encoding="utf-8") as handle:
        json.dump(policy, handle, indent=2, ensure_ascii=False, allow_nan=False)

    summary = {
        "pipeline": {
            "name": "GSE143157 Step 3 confounder audit and nested-CV freeze",
            "version": VERSION,
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": time.time() - started,
        },
        "input": {
            "step2b_dir": str(root),
            "cohort_decision_sha256": sha256(paths["cohort_decision.json"]),
            "primary_hdf5_file": "GSE143157_primary_clean77.h5",
            "primary_hdf5_shapes": step2b_decision["primary_hdf5"]["matrix_stats"],
        },
        "parameters": {
            "outer_folds": args.outer_folds,
            "outer_repeats": args.outer_repeats,
            "inner_folds": args.inner_folds,
            "fold_candidates_per_outer_repeat": args.fold_candidates,
            "permutations": args.permutations,
            "random_seed": args.seed,
            "p_value_adjustment": "Benjamini-Hochberg across four covariates within each cohort",
        },
        "cohorts": {
            cohort: {
                "donors": len(frame),
                "condition_counts": {
                    str(key): int(value) for key, value in frame["condition"].value_counts().items()
                },
                "outer_assignment_rows": len(outer_tables[cohort]),
                "inner_assignment_rows": len(inner_tables[cohort]),
                "outer_balance_rows": len(balance_tables[cohort]),
                "cv_details": cv_details[cohort],
            }
            for cohort, frame in metadata.items()
        },
        "primary_covariate_audit": audits["primary_clean77"].to_dict(orient="records"),
        "frozen_policy_file": "model_input_policy.json",
        "guardrails": [
            "No CpG values were read to create or optimize folds.",
            "The same frozen donor folds must be used for 5mC, 5hmC, and combined models.",
            "All imputation, feature filtering/selection, scaling, residualization, and tuning occur inside training data only.",
            "Outer-test data are used once for final prediction and never for model selection.",
            "Full-78 is a required sensitivity analysis; the 76-donor cohort is an optional stress-test.",
        ],
    }
    with (output / "step3_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False, allow_nan=False)

    save_plot_covariates(metadata["primary_clean77"], output / "covariate_balance_primary.png")
    save_plot_fold_balance(
        balance_tables["primary_clean77"], output / "outer_fold_balance_primary.png"
    )
    write_report(output / "STEP3_REPORT.html", summary, audits["primary_clean77"])
    (output / "SUCCESS.txt").write_text(
        "SUCCESS\n"
        f"Completed UTC: {summary['pipeline']['completed_utc']}\n"
        "Primary cohort: Clean-77\n"
        f"Nested CV: {args.outer_folds} outer folds x {args.outer_repeats} repeats; {args.inner_folds} inner folds\n"
        "Same donor folds frozen for 5mC, 5hmC, and combined models\n"
        "No CpG values were used to create folds\n",
        encoding="utf-8",
    )
    print("\nSTEP 3 SUCCESS")
    print(f"Open report: {output / 'STEP3_REPORT.html'}")


def main() -> int:
    try:
        run(parse_args())
    except (ValueError, OSError, KeyError, json.JSONDecodeError, pd.errors.ParserError) as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
