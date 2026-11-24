#!/usr/bin/env python3
"""Strict integrity checks for GSE143157 Step 6 outputs."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


PRIMARY_MODEL = "combined__molecular_plus_covariates"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read_gzip_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def check(output: Path, step5: Path) -> None:
    required = [
        "RUN_CONFIGURATION.json",
        "step6_summary.json",
        "SUCCESS.txt",
        "PATHWAY_ENRICHMENT_STATUS.txt",
        "STEP6_REPORT.html",
        "annotated_feature_stability_primary_clean77.csv",
        "primary_stable_cpgs_primary_clean77.csv",
        "primary_stable_cpg_gene_map_primary_clean77.csv",
        "candidate_gene_evidence_primary_clean77.csv",
        "published_signal_overlap_primary_clean77.csv",
        "annotation_qc_summary_primary_clean77.csv",
        "gometh_selected_5mC.txt.gz",
        "gometh_selected_5hmC.txt.gz",
        "gometh_background_5mC.txt.gz",
        "gometh_background_5hmC.txt.gz",
        "stable_cpg_annotation_primary.png",
        "effect_stability_primary.png",
        "genomic_context_primary.png",
        "gene_evidence_primary.png",
    ]
    missing = [name for name in required if not (output / name).is_file()]
    require(not missing, f"Missing Step 6 files: {missing}")
    require(all((output / name).stat().st_size > 0 for name in required), "An output file is empty.")

    summary = json.loads((output / "step6_summary.json").read_text(encoding="utf-8"))
    policy = summary["analysis_policy"]
    counts = summary["counts"]
    require(policy["primary_model"] == PRIMARY_MODEL, "Wrong primary model.")
    require(policy["stability_is_not_a_p_value"] is True, "Stability interpretation flag is wrong.")
    require(policy["automatic_pathway_enrichment_run"] is False, "Unexpected pathway claims found.")

    step5_input = pd.read_csv(step5 / "feature_stability_primary_clean77.csv")
    annotated = pd.read_csv(output / "annotated_feature_stability_primary_clean77.csv")
    primary = pd.read_csv(output / "primary_stable_cpgs_primary_clean77.csv")
    gene_map = pd.read_csv(output / "primary_stable_cpg_gene_map_primary_clean77.csv")
    genes = pd.read_csv(output / "candidate_gene_evidence_primary_clean77.csv")
    published = pd.read_csv(output / "published_signal_overlap_primary_clean77.csv")

    require(len(annotated) == len(step5_input), "Annotated feature table row count changed.")
    require(
        annotated[["model_id", "feature_name"]].duplicated().sum() == 0,
        "Annotated model-feature rows are duplicated.",
    )
    require(len(primary) == counts["primary_stable_cpgs"], "Primary candidate count mismatch.")
    require(primary["feature_name"].is_unique, "Primary candidate CpGs are duplicated.")
    require(set(primary["model_id"]) == {PRIMARY_MODEL}, "Primary candidates came from another model.")
    require(
        set(primary["stability_tier"]).issubset({"A_high_stability", "B_moderate_stability"}),
        "An undeclared stability tier entered the primary list.",
    )
    expected_mask = (
        (annotated["model_id"] == PRIMARY_MODEL)
        & (
            (
                (annotated["selection_frequency"] >= 0.80)
                & (annotated["coefficient_direction_consistency"] >= 0.90)
            )
            | (
                (annotated["selection_frequency"] >= 0.60)
                & (annotated["coefficient_direction_consistency"] >= 0.80)
            )
        )
    )
    expected = set(annotated.loc[expected_mask, "feature_name"])
    require(set(primary["feature_name"]) == expected, "Primary stability filter was not reproduced.")
    require(primary["manifest_found"].all(), "A primary candidate lacks manifest annotation.")
    require(primary["selection_frequency"].between(0, 1).all(), "Invalid selection frequency.")
    require(
        primary["coefficient_direction_consistency"].between(0.5, 1).all(),
        "Invalid direction consistency.",
    )
    require(
        int(primary["stability_tier"].eq("A_high_stability").sum()) == counts["tier_a"],
        "Tier A count mismatch.",
    )
    require(
        int(primary["stability_tier"].eq("B_moderate_stability").sum()) == counts["tier_b"],
        "Tier B count mismatch.",
    )
    require(gene_map["gene_symbol"].fillna("").ne("").all(), "Blank gene symbols in gene map.")
    require(gene_map[["gene_symbol", "feature_name"]].duplicated().sum() == 0, "Duplicate gene-CpG map.")
    require(genes["gene_symbol"].is_unique, "Gene evidence contains duplicate genes.")
    require(len(published) == counts["published_reference_rows"], "Published-reference row mismatch.")

    for modality in ("5mC", "5hmC"):
        selected = read_gzip_lines(output / f"gometh_selected_{modality}.txt.gz")
        expected_selected = sorted(
            primary.loc[
                (primary["source_modality"] == modality) & primary["interpretation_eligible"],
                "cpg_id",
            ].astype(str).unique()
        )
        require(selected == expected_selected, f"GOmeth selected list is wrong for {modality}.")
        require(
            len(selected) == counts[f"gometh_selected_{modality}"],
            f"GOmeth selected count mismatch for {modality}.",
        )
        with gzip.open(output / f"gometh_background_{modality}.txt.gz", "rt", encoding="utf-8") as handle:
            background_count = sum(1 for line in handle if line.strip())
        require(
            background_count == counts[f"gometh_background_{modality}"],
            f"GOmeth background count mismatch for {modality}.",
        )
        require(set(selected).issubset(set(primary["cpg_id"])), "Selected GOmeth CpG not in primary table.")

    require(np.isfinite(primary["mean_effect_msa_minus_control"]).all(), "Primary effects are not finite.")
    print("STEP 6 INTEGRITY CHECK PASSED")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--step5-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        check(args.output.resolve(), args.step5_dir.resolve())
    except (ValueError, OSError, KeyError, json.JSONDecodeError, pd.errors.ParserError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
