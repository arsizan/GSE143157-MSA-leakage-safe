#!/usr/bin/env python3
"""GSE143157 Step 6: annotate and audit stable CpG features.

This step is intentionally descriptive. Selection frequency is not treated as a
p-value, and pathway enrichment is not performed with a naive gene-level test.
Instead, the script creates inputs for a separate bias-aware missMethyl/GOmeth
analysis.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import html
import json
from pathlib import Path
import sys
import time

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PIPELINE_VERSION = "1.0.0"
PRIMARY_MODEL = "combined__molecular_plus_covariates"
TIER_A_FREQ = 0.80
TIER_A_DIRECTION = 0.90
TIER_B_FREQ = 0.60
TIER_B_DIRECTION = 0.80
COMMON_SNP_MAF = 0.05

PAPER_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC7063795/"
MISSMETHYL_URL = (
    "https://bioconductor.org/packages/release/bioc/vignettes/"
    "missMethyl/inst/doc/missMethyl.html"
)

MANIFEST_FIELDS = {
    "Genome_Build": "genome_build",
    "CHR": "chromosome",
    "MAPINFO": "mapinfo",
    "UCSC_RefGene_Name": "ucsc_refgene_name",
    "UCSC_RefGene_Accession": "ucsc_refgene_accession",
    "UCSC_RefGene_Group": "ucsc_refgene_group",
    "UCSC_CpG_Islands_Name": "ucsc_cpg_island_name",
    "Relation_to_UCSC_CpG_Island": "relation_to_cpg_island",
    "Phantom5_Enhancers": "phantom5_enhancers",
    "DMR": "manifest_dmr",
    "450k_Enhancer": "enhancer_450k",
    "HMM_Island": "hmm_island",
    "Regulatory_Feature_Name": "regulatory_feature_name",
    "Regulatory_Feature_Group": "regulatory_feature_group",
    "GencodeBasicV12_NAME": "gencode_basic_name",
    "GencodeBasicV12_Group": "gencode_basic_group",
    "DNase_Hypersensitivity_NAME": "dnase_hypersensitivity_name",
    "DNase_Hypersensitivity_Evidence_Count": "dnase_evidence_count",
    "OpenChromatin_NAME": "open_chromatin_name",
    "OpenChromatin_Evidence_Count": "open_chromatin_evidence_count",
    "TFBS_NAME": "tfbs_name",
    "TFBS_Evidence_Count": "tfbs_evidence_count",
    "SNP_ID": "snp_id",
    "SNP_DISTANCE": "snp_distance",
    "SNP_MinorAlleleFrequency": "snp_minor_allele_frequency",
    "Random_Loci": "random_loci",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--step5-dir", required=True, type=Path)
    parser.add_argument("--step2b-dir", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--published-signals", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--chunk-size", type=int, default=50000)
    return parser.parse_args()


def file_sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def clean_text(value) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return str(value).strip()


def unique_semicolon_values(value) -> list[str]:
    seen: set[str] = set()
    result: list[str] = []
    for item in clean_text(value).split(";"):
        item = item.strip()
        if item and item not in seen:
            seen.add(item)
            result.append(item)
    return result


def max_manifest_maf(value) -> float:
    values: list[float] = []
    for item in clean_text(value).split(";"):
        try:
            values.append(float(item))
        except ValueError:
            continue
    return max(values) if values else np.nan


def feature_to_cpg(value: str) -> str:
    text = clean_text(value)
    return text.split("::", 1)[-1]


def load_manifest_subset(path: Path, target_cpgs: set[str]) -> tuple[pd.DataFrame, int]:
    records: list[dict[str, object]] = []
    assay_rows = 0
    found: set[str] = set()
    print(f"Scanning EPIC manifest for {len(target_cpgs):,} unique CpGs...")
    with gzip.open(path, "rt", encoding="utf-8-sig", errors="replace", newline="") as handle:
        for line in handle:
            if line.startswith("[Assay]"):
                break
        else:
            raise ValueError("The manifest does not contain an [Assay] section.")
        header = next(csv.reader([next(handle)]))
        reader = csv.DictReader(handle, fieldnames=header)
        for row in reader:
            section_marker = clean_text(row.get("IlmnID"))
            if section_marker.startswith("["):
                break
            assay_rows += 1
            cpg_id = clean_text(row.get("Name")) or clean_text(row.get("IlmnID"))
            if cpg_id in target_cpgs and cpg_id not in found:
                record: dict[str, object] = {"cpg_id": cpg_id}
                for source, destination in MANIFEST_FIELDS.items():
                    record[destination] = clean_text(row.get(source))
                records.append(record)
                found.add(cpg_id)
            if assay_rows % 100000 == 0:
                print(
                    f"  manifest rows {assay_rows:,}; targets found "
                    f"{len(found):,}/{len(target_cpgs):,}"
                )
    print(f"Manifest scan complete: {len(found):,}/{len(target_cpgs):,} targets found.")
    return pd.DataFrame(records), assay_rows


def add_manifest_qc(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    result["manifest_found"] = result["genome_build"].fillna("").astype(str).ne("")
    result["chromosome"] = result["chromosome"].fillna("").astype(str).str.replace("chr", "", regex=False)
    result["sex_chromosome_flag"] = result["chromosome"].str.upper().isin({"X", "Y"})
    result["random_locus_flag"] = (
        result["random_loci"].fillna("").astype(str).str.upper().isin({"TRUE", "1", "YES"})
    )
    result["max_manifest_snp_maf"] = result["snp_minor_allele_frequency"].map(max_manifest_maf)
    result["common_manifest_snp_flag"] = result["max_manifest_snp_maf"].fillna(-1).ge(COMMON_SNP_MAF)
    result["gene_symbols"] = result.apply(
        lambda row: ";".join(
            unique_semicolon_values(row.get("ucsc_refgene_name"))
            or unique_semicolon_values(row.get("gencode_basic_name"))
        ),
        axis=1,
    )
    result["gene_annotation_source"] = np.where(
        result["ucsc_refgene_name"].fillna("").astype(str).ne(""),
        "UCSC_RefGene_Name",
        np.where(
            result["gencode_basic_name"].fillna("").astype(str).ne(""),
            "GencodeBasicV12_NAME",
            "intergenic_or_unmapped",
        ),
    )
    result["gene_count"] = result["gene_symbols"].map(lambda value: len(unique_semicolon_values(value)))
    result["intergenic_flag"] = result["gene_count"].eq(0)
    result["multi_gene_flag"] = result["gene_count"].gt(1)

    def qc_status(row: pd.Series) -> str:
        if not bool(row["manifest_found"]):
            return "fail_missing_manifest_annotation"
        if bool(row["sex_chromosome_flag"]):
            return "review_sex_chromosome"
        if bool(row["random_locus_flag"]):
            return "review_random_locus"
        if bool(row["common_manifest_snp_flag"]):
            return "review_common_manifest_snp"
        return "pass_manifest_audit"

    result["manifest_qc_status"] = result.apply(qc_status, axis=1)
    result["cross_reactive_probe_status"] = (
        "upstream_processed_matrix_reported_ChAMP_standard_filtering;_not_independently_retested_here"
    )
    result["interpretation_eligible"] = result["manifest_qc_status"].eq("pass_manifest_audit")
    return result


def assign_stability_tier(row: pd.Series) -> str:
    frequency = float(row["selection_frequency"])
    direction = float(row["coefficient_direction_consistency"])
    if frequency >= TIER_A_FREQ and direction >= TIER_A_DIRECTION:
        return "A_high_stability"
    if frequency >= TIER_B_FREQ and direction >= TIER_B_DIRECTION:
        return "B_moderate_stability"
    if frequency >= 0.50 and direction >= 0.70:
        return "C_review_only"
    return "not_primary_stable"


def enrich_stability(stability: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    frame = stability.copy()
    frame["cpg_id"] = frame["feature_name"].map(feature_to_cpg)
    frame["positive_coefficient_fraction"] = pd.to_numeric(
        frame["positive_coefficient_fraction"], errors="coerce"
    )
    frame["selection_frequency"] = pd.to_numeric(frame["selection_frequency"], errors="coerce")
    frame["coefficient_mean_when_selected"] = pd.to_numeric(
        frame["coefficient_mean_when_selected"], errors="coerce"
    )
    frame["mean_effect_msa_minus_control"] = pd.to_numeric(
        frame["mean_effect_msa_minus_control"], errors="coerce"
    )
    frame["coefficient_direction_consistency"] = np.maximum(
        frame["positive_coefficient_fraction"], 1.0 - frame["positive_coefficient_fraction"]
    )
    frame["effect_direction"] = np.select(
        [frame["mean_effect_msa_minus_control"] > 0, frame["mean_effect_msa_minus_control"] < 0],
        ["MSA_higher", "MSA_lower"],
        default="no_mean_difference",
    )
    frame["coefficient_direction"] = np.select(
        [frame["coefficient_mean_when_selected"] > 0, frame["coefficient_mean_when_selected"] < 0],
        ["toward_MSA", "toward_Control"],
        default="zero",
    )
    frame["coefficient_effect_direction_agree"] = (
        np.sign(frame["coefficient_mean_when_selected"])
        == np.sign(frame["mean_effect_msa_minus_control"])
    )
    frame["stability_tier"] = frame.apply(assign_stability_tier, axis=1)
    result = frame.merge(manifest, on="cpg_id", how="left", validate="many_to_one")
    result = add_manifest_qc(result)
    return result


def make_gene_map(primary: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, row in primary.iterrows():
        for gene in unique_semicolon_values(row["gene_symbols"]):
            rows.append(
                {
                    "gene_symbol": gene,
                    "feature_name": row["feature_name"],
                    "cpg_id": row["cpg_id"],
                    "source_modality": row["source_modality"],
                    "stability_tier": row["stability_tier"],
                    "selection_frequency": row["selection_frequency"],
                    "coefficient_direction_consistency": row["coefficient_direction_consistency"],
                    "effect_direction": row["effect_direction"],
                    "mean_effect_msa_minus_control": row["mean_effect_msa_minus_control"],
                    "chromosome": row["chromosome"],
                    "mapinfo": row["mapinfo"],
                    "ucsc_refgene_group": row["ucsc_refgene_group"],
                    "relation_to_cpg_island": row["relation_to_cpg_island"],
                    "manifest_qc_status": row["manifest_qc_status"],
                    "interpretation_eligible": row["interpretation_eligible"],
                }
            )
    columns = [
        "gene_symbol",
        "feature_name",
        "cpg_id",
        "source_modality",
        "stability_tier",
        "selection_frequency",
        "coefficient_direction_consistency",
        "effect_direction",
        "mean_effect_msa_minus_control",
        "chromosome",
        "mapinfo",
        "ucsc_refgene_group",
        "relation_to_cpg_island",
        "manifest_qc_status",
        "interpretation_eligible",
    ]
    return pd.DataFrame(rows, columns=columns)


def make_gene_evidence(gene_map: pd.DataFrame) -> pd.DataFrame:
    if gene_map.empty:
        return pd.DataFrame(
            columns=[
                "rank",
                "gene_symbol",
                "stable_cpg_count",
                "cpg_ids",
                "modalities",
                "best_stability_tier",
                "max_selection_frequency",
                "mean_selection_frequency",
                "minimum_direction_consistency",
                "effect_directions",
                "all_manifest_qc_pass",
            ]
        )
    rows: list[dict[str, object]] = []
    tier_priority = {"A_high_stability": 2, "B_moderate_stability": 1}
    for gene, group in gene_map.groupby("gene_symbol", sort=False):
        tiers = list(group["stability_tier"].astype(str))
        best_tier = max(tiers, key=lambda value: tier_priority.get(value, 0))
        rows.append(
            {
                "gene_symbol": gene,
                "stable_cpg_count": int(group["cpg_id"].nunique()),
                "cpg_ids": ";".join(sorted(group["cpg_id"].astype(str).unique())),
                "modalities": ";".join(sorted(group["source_modality"].astype(str).unique())),
                "best_stability_tier": best_tier,
                "max_selection_frequency": float(group["selection_frequency"].max()),
                "mean_selection_frequency": float(group["selection_frequency"].mean()),
                "minimum_direction_consistency": float(
                    group["coefficient_direction_consistency"].min()
                ),
                "effect_directions": ";".join(sorted(group["effect_direction"].astype(str).unique())),
                "all_manifest_qc_pass": bool(group["interpretation_eligible"].all()),
                "tier_priority": tier_priority.get(best_tier, 0),
            }
        )
    result = pd.DataFrame(rows).sort_values(
        ["tier_priority", "max_selection_frequency", "stable_cpg_count", "gene_symbol"],
        ascending=[False, False, False, True],
    )
    result.insert(0, "rank", range(1, len(result) + 1))
    return result.drop(columns="tier_priority")


def make_published_overlap(
    published: pd.DataFrame,
    annotated: pd.DataFrame,
    primary: pd.DataFrame,
    manifest: pd.DataFrame,
    h5_presence: dict[str, set[str]],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    primary_keys = set(zip(primary["source_modality"], primary["cpg_id"]))
    for _, ref in published.iterrows():
        modality = clean_text(ref["modality"])
        cpg_id = clean_text(ref["cpg_id"])
        matching = annotated.loc[
            (annotated["source_modality"] == modality) & (annotated["cpg_id"] == cpg_id)
        ]
        if matching.empty:
            selected_outer_fits = 0
            max_frequency = 0.0
            observed_effect = np.nan
            observed_direction = "not_selected"
        else:
            best = matching.sort_values("selection_frequency", ascending=False).iloc[0]
            selected_outer_fits = int(pd.to_numeric(matching["selected_outer_fits"]).max())
            max_frequency = float(matching["selection_frequency"].max())
            observed_effect = float(best["mean_effect_msa_minus_control"])
            observed_direction = clean_text(best["effect_direction"])
        manifest_row = manifest.loc[manifest["cpg_id"] == cpg_id]
        manifest_gene = ""
        manifest_chr = ""
        manifest_pos = ""
        if not manifest_row.empty:
            manifest_gene = clean_text(manifest_row.iloc[0].get("ucsc_refgene_name"))
            manifest_chr = clean_text(manifest_row.iloc[0].get("chromosome"))
            manifest_pos = clean_text(manifest_row.iloc[0].get("mapinfo"))
        reported_effect = pd.to_numeric(pd.Series([ref["reported_effect_msa_minus_control"]]), errors="coerce").iloc[0]
        direction_match = False
        if pd.notna(reported_effect) and pd.notna(observed_effect) and observed_effect != 0:
            direction_match = bool(np.sign(float(reported_effect)) == np.sign(float(observed_effect)))
        rows.append(
            {
                **ref.to_dict(),
                "present_in_step2b_modality_background": cpg_id in h5_presence.get(modality, set()),
                "manifest_gene": manifest_gene,
                "manifest_chromosome": manifest_chr,
                "manifest_position": manifest_pos,
                "selected_in_any_step5_outer_fit": selected_outer_fits > 0,
                "maximum_selected_outer_fits_any_model": selected_outer_fits,
                "maximum_selection_frequency_any_model": max_frequency,
                "primary_stable_exact_overlap": (modality, cpg_id) in primary_keys,
                "observed_step5_mean_effect": observed_effect,
                "observed_step5_effect_direction": observed_direction,
                "reported_observed_direction_match": direction_match,
            }
        )
    return pd.DataFrame(rows)


def decode_ids(values) -> list[str]:
    decoded: list[str] = []
    for value in values:
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8"))
        else:
            decoded.append(str(value))
    return decoded


def export_backgrounds(
    h5_path: Path,
    output: Path,
    published_ids: set[str],
    chunk_size: int,
) -> tuple[dict[str, int], dict[str, set[str]]]:
    counts: dict[str, int] = {}
    presence: dict[str, set[str]] = {"5mC": set(), "5hmC": set()}
    with h5py.File(h5_path, "r") as handle:
        for modality in ("5mC", "5hmC"):
            dataset = handle[f"{modality}/cpg_ids"]
            destination = output / f"gometh_background_{modality}.txt.gz"
            count = 0
            with gzip.open(destination, "wt", encoding="utf-8", newline="\n") as writer:
                for start in range(0, len(dataset), chunk_size):
                    ids = decode_ids(dataset[start : start + chunk_size])
                    if ids:
                        writer.write("\n".join(ids) + "\n")
                    count += len(ids)
                    presence[modality].update(set(ids).intersection(published_ids))
            counts[modality] = count
            print(f"Exported GOmeth background for {modality}: {count:,} CpGs")
    return counts, presence


def write_selected_gometh_inputs(primary: pd.DataFrame, output: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    for modality in ("5mC", "5hmC"):
        selected = sorted(
            primary.loc[
                (primary["source_modality"] == modality) & primary["interpretation_eligible"],
                "cpg_id",
            ].astype(str).unique()
        )
        with gzip.open(
            output / f"gometh_selected_{modality}.txt.gz", "wt", encoding="utf-8", newline="\n"
        ) as writer:
            if selected:
                writer.write("\n".join(selected) + "\n")
        counts[modality] = len(selected)
    return counts


def plot_stable_candidates(primary: pd.DataFrame, output: Path) -> None:
    frame = primary.iloc[::-1].copy()
    frame["display_gene"] = frame["gene_symbols"].replace("", "IGR/unmapped").str.split(";").str[0]
    frame["label"] = frame.apply(
        lambda row: f"{row['source_modality']}  {row['cpg_id']}  {row['display_gene']}", axis=1
    )
    colors = frame["source_modality"].map({"5mC": "#2b8cbe", "5hmC": "#1aa77e"})
    fig, ax = plt.subplots(figsize=(12, 9))
    ax.barh(frame["label"], frame["selection_frequency"], color=colors)
    ax.axvline(TIER_B_FREQ, color="#777777", linestyle="--", linewidth=1)
    for index, (_, row) in enumerate(frame.iterrows()):
        tier = "A" if row["stability_tier"].startswith("A_") else "B"
        ax.text(min(float(row["selection_frequency"]) + 0.012, 0.98), index, tier, va="center", fontsize=9)
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("Selection frequency across 50 outer fits")
    ax.set_title("Primary stable CpGs with EPIC gene annotation")
    ax.grid(axis="x", alpha=0.2)
    ax.legend(
        handles=[
            plt.Line2D([0], [0], color="#2b8cbe", lw=8, label="5mC"),
            plt.Line2D([0], [0], color="#1aa77e", lw=8, label="5hmC"),
        ],
        loc="lower right",
    )
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_effect_stability(primary: pd.DataFrame, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 7))
    colors = {"5mC": "#2b8cbe", "5hmC": "#1aa77e"}
    for modality, group in primary.groupby("source_modality"):
        ax.scatter(
            group["mean_effect_msa_minus_control"] * 100,
            group["selection_frequency"],
            s=65,
            color=colors.get(modality, "#777777"),
            label=modality,
            alpha=0.9,
        )
        ordered = group.sort_values(
            ["selection_frequency", "mean_effect_msa_minus_control", "cpg_id"],
            ascending=[False, True, True],
        )
        offsets = [(5, 5), (5, 15), (5, -13), (5, 25), (5, -23)]
        for label_index, (_, row) in enumerate(ordered.iterrows()):
            gene = unique_semicolon_values(row["gene_symbols"])
            label = gene[0] if gene else row["cpg_id"]
            if gene and (group["gene_symbols"] == row["gene_symbols"]).sum() > 1:
                label = f"{label} ({row['cpg_id']})"
            ax.annotate(
                label,
                (row["mean_effect_msa_minus_control"] * 100, row["selection_frequency"]),
                xytext=offsets[label_index % len(offsets)],
                textcoords="offset points",
                fontsize=8,
                arrowprops={"arrowstyle": "-", "color": "#999999", "lw": 0.5},
            )
    ax.axvline(0, color="#777777", linewidth=1)
    ax.axhline(TIER_B_FREQ, color="#777777", linestyle="--", linewidth=1)
    ax.set_xlabel("Mean outer-training MSA − Control effect (beta percentage points)")
    ax.set_ylabel("Selection frequency")
    ax.set_ylim(0.55, 0.98)
    ax.set_title("Direction and stability of primary CpG candidates")
    ax.legend()
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def first_annotation(value: str, default: str) -> str:
    items = unique_semicolon_values(value)
    return items[0] if items else default


def plot_genomic_context(primary: pd.DataFrame, output: Path) -> None:
    frame = primary.copy()
    frame["gene_region"] = frame["ucsc_refgene_group"].map(
        lambda value: first_annotation(value, "Intergenic/unmapped")
    )
    frame["gene_region"] = frame["gene_region"].replace(
        {"5URT": "5'UTR", "5UTR": "5'UTR", "3URT": "3'UTR", "3UTR": "3'UTR"}
    )
    frame["island_relation"] = frame["relation_to_cpg_island"].replace("", "Open_Sea").fillna("Open_Sea")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, column, title in [
        (axes[0], "gene_region", "Gene-relative annotation"),
        (axes[1], "island_relation", "CpG-island relation"),
    ]:
        table = pd.crosstab(frame[column], frame["source_modality"])
        table = table.reindex(columns=["5mC", "5hmC"], fill_value=0)
        table.plot.barh(ax=ax, color=["#2b8cbe", "#1aa77e"])
        ax.set_title(title)
        ax.set_xlabel("Number of primary stable CpGs")
        ax.set_ylabel("")
        ax.grid(axis="x", alpha=0.2)
    fig.suptitle("Genomic context of primary stable CpGs")
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_gene_evidence(gene_evidence: pd.DataFrame, output: Path) -> None:
    if gene_evidence.empty:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No gene-mapped primary candidates", ha="center", va="center")
        ax.axis("off")
    else:
        frame = gene_evidence.head(15).iloc[::-1]
        colors = frame["modalities"].map(
            lambda value: "#7b6fd0" if ";" in value else ("#2b8cbe" if value == "5mC" else "#1aa77e")
        )
        fig, ax = plt.subplots(figsize=(9, 7))
        ax.barh(frame["gene_symbol"], frame["max_selection_frequency"], color=colors)
        ax.axvline(TIER_B_FREQ, color="#777777", linestyle="--", linewidth=1)
        ax.set_xlim(0, 1.0)
        ax.set_xlabel("Maximum CpG selection frequency")
        ax.set_title("Top gene-mapped evidence (descriptive, not pathway significance)")
        ax.grid(axis="x", alpha=0.2)
        ax.legend(
            handles=[
                plt.Line2D([0], [0], color="#2b8cbe", lw=8, label="5mC"),
                plt.Line2D([0], [0], color="#1aa77e", lw=8, label="5hmC"),
                plt.Line2D([0], [0], color="#7b6fd0", lw=8, label="both"),
            ],
            loc="lower right",
        )
        fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def dataframe_html(frame: pd.DataFrame, columns: list[str], limit: int = 30) -> str:
    table = frame.loc[:, columns].head(limit).copy()
    for column in table.select_dtypes(include=["float", "float64"]).columns:
        table[column] = table[column].map(lambda value: "" if pd.isna(value) else f"{value:.3f}")
    return table.to_html(index=False, escape=True, border=0, classes="results")


def write_report(
    output: Path,
    primary: pd.DataFrame,
    genes: pd.DataFrame,
    published_overlap: pd.DataFrame,
    summary: dict[str, object],
) -> None:
    counts = summary["counts"]
    candidate_table = dataframe_html(
        primary.sort_values("selection_frequency", ascending=False),
        [
            "source_modality",
            "cpg_id",
            "gene_symbols",
            "stability_tier",
            "selection_frequency",
            "coefficient_direction_consistency",
            "effect_direction",
            "mean_effect_msa_minus_control",
            "ucsc_refgene_group",
            "relation_to_cpg_island",
            "manifest_qc_status",
        ],
        limit=30,
    )
    published_table = dataframe_html(
        published_overlap,
        [
            "modality",
            "cpg_id",
            "gene",
            "publication_status",
            "primary_stable_exact_overlap",
            "maximum_selection_frequency_any_model",
        ],
        limit=20,
    )
    gene_table = dataframe_html(
        genes,
        [
            "rank",
            "gene_symbol",
            "stable_cpg_count",
            "cpg_ids",
            "modalities",
            "best_stability_tier",
            "max_selection_frequency",
            "effect_directions",
            "all_manifest_qc_pass",
        ],
        limit=30,
    )
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>GSE143157 Step 6 report</title>
<style>
body{{font-family:Arial,sans-serif;max-width:1200px;margin:28px auto;padding:0 18px;color:#222;line-height:1.45}}
h1,h2{{color:#17365d}} .card{{display:inline-block;border:1px solid #ddd;border-radius:8px;padding:12px 18px;margin:5px;background:#f7fafc}}
.warning{{border-left:5px solid #d97706;background:#fff7ed;padding:12px;margin:16px 0}}
.ok{{border-left:5px solid #16866b;background:#effaf6;padding:12px;margin:16px 0}}
table.results{{border-collapse:collapse;width:100%;font-size:13px;margin:12px 0 25px}}
table.results th,table.results td{{border:1px solid #ddd;padding:6px;text-align:left;vertical-align:top}}
table.results th{{background:#eaf1f8}} img{{max-width:100%;border:1px solid #ddd;margin:10px 0}}
code{{background:#f1f3f5;padding:2px 4px}} a{{color:#0b5cab}}
</style></head><body>
<h1>GSE143157 Step 6 — CpG annotation and interpretation audit</h1>
<div class="ok"><b>Completed.</b> The primary source is the combined molecular + covariate model from the leakage-safe Step 5 analysis.</div>
<div class="card"><b>Primary stable CpGs</b><br>{counts['primary_stable_cpgs']}</div>
<div class="card"><b>Tier A / Tier B</b><br>{counts['tier_a']} / {counts['tier_b']}</div>
<div class="card"><b>Gene mapped</b><br>{counts['gene_mapped_primary']}</div>
<div class="card"><b>Unique candidate genes</b><br>{counts['unique_candidate_genes']}</div>
<div class="card"><b>Published exact overlaps</b><br>{counts['published_primary_exact_overlap']}</div>
<h2>What the stability tiers mean</h2>
<p>Tier A requires selection in at least 80% of outer fits and at least 90% coefficient-direction consistency. Tier B requires at least 60% selection and 80% direction consistency. These are descriptive reproducibility rules, not hypothesis-test p-values.</p>
<div class="warning"><b>Important:</b> A stable feature is not automatically a causal locus, a clinically validated biomarker, or a statistically significant EWAS probe. Gene and pathway interpretation must remain exploratory.</div>
<h2>Primary stable CpGs</h2>{candidate_table}
<img src="stable_cpg_annotation_primary.png" alt="Stable CpG annotation">
<img src="effect_stability_primary.png" alt="Effect and stability">
<h2>Genomic context</h2><img src="genomic_context_primary.png" alt="Genomic context">
<h2>Gene-mapped evidence</h2>{gene_table}
<img src="gene_evidence_primary.png" alt="Gene evidence">
<h2>Exact comparison with the original GSE143157 paper</h2>
<p>The original report identified five genome-wide significant 5mC probes and highlighted a 5mC-to-5hmC shift at AREL1/cg08753407. This table checks exact CpG-and-modality overlap without treating absence as a contradiction, because the published EWAS and the present predictive stability analysis answer different questions.</p>
{published_table}
<h2>Pathway-analysis policy</h2>
<p>No naive gene-level enrichment was run. Methylation arrays have probe-number and multi-gene annotation bias. Files named <code>gometh_selected_*.txt.gz</code> and <code>gometh_background_*.txt.gz</code> were prepared for a separate bias-aware missMethyl/GOmeth analysis. Any future enrichment must be labelled exploratory because the input CpGs are ML-stable features rather than EWAS-significant probes.</p>
<h2>Sources</h2>
<ul><li><a href="{html.escape(PAPER_URL)}">Original GSE143157 publication</a></li>
<li><a href="{html.escape(MISSMETHYL_URL)}">missMethyl bias-aware gene-set guidance</a></li>
<li>Illumina MethylationEPIC v1.0 manifest supplied locally with GPL21145; Genome Build 37 coordinates.</li></ul>
</body></html>"""
    output.write_text(document, encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    started = time.time()
    step5 = args.step5_dir.resolve()
    step2b = args.step2b_dir.resolve()
    manifest_path = args.manifest.resolve()
    published_path = args.published_signals.resolve()
    output = args.output.resolve()
    required = [
        step5 / "SUCCESS.txt",
        step5 / "step5_summary.json",
        step5 / "feature_stability_primary_clean77.csv",
        step2b / "GSE143157_primary_clean77.h5",
        manifest_path,
        published_path,
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise ValueError(f"Required inputs are missing: {missing}")
    if output.joinpath("SUCCESS.txt").exists():
        raise ValueError(f"A successful Step 6 result already exists: {output}")
    output.mkdir(parents=True, exist_ok=True)

    step5_summary = json.loads((step5 / "step5_summary.json").read_text(encoding="utf-8"))
    if step5_summary["parameters"]["analysis_mode"] != "FULL_RESEARCH_RUN":
        raise ValueError("Step 5 is not a full research run.")
    if step5_summary["cohort"]["name"] != "primary_clean77":
        raise ValueError("Step 5 does not use the expected primary_clean77 cohort.")

    stability = pd.read_csv(step5 / "feature_stability_primary_clean77.csv")
    published = pd.read_csv(published_path)
    target_cpgs = set(stability["feature_name"].map(feature_to_cpg)) | set(
        published["cpg_id"].astype(str)
    )
    manifest, manifest_assay_rows = load_manifest_subset(manifest_path, target_cpgs)
    annotated = enrich_stability(stability, manifest)

    primary_model_rows = annotated.loc[annotated["model_id"] == PRIMARY_MODEL].copy()
    primary = primary_model_rows.loc[
        primary_model_rows["stability_tier"].isin({"A_high_stability", "B_moderate_stability"})
    ].copy()
    primary = primary.sort_values(
        ["selection_frequency", "coefficient_direction_consistency", "feature_name"],
        ascending=[False, False, True],
    )
    if primary.empty:
        raise ValueError("No primary stable candidates passed the declared thresholds.")

    h5_path = step2b / "GSE143157_primary_clean77.h5"
    published_ids = set(published["cpg_id"].astype(str))
    background_counts, h5_presence = export_backgrounds(
        h5_path, output, published_ids, args.chunk_size
    )
    selected_gometh_counts = write_selected_gometh_inputs(primary, output)

    gene_map = make_gene_map(primary)
    gene_evidence = make_gene_evidence(gene_map)
    published_overlap = make_published_overlap(
        published, annotated, primary, manifest, h5_presence
    )

    annotated.to_csv(output / "annotated_feature_stability_primary_clean77.csv", index=False)
    primary.to_csv(output / "primary_stable_cpgs_primary_clean77.csv", index=False)
    gene_map.to_csv(output / "primary_stable_cpg_gene_map_primary_clean77.csv", index=False)
    gene_evidence.to_csv(output / "candidate_gene_evidence_primary_clean77.csv", index=False)
    published_overlap.to_csv(
        output / "published_signal_overlap_primary_clean77.csv", index=False
    )

    tier_a = int(primary["stability_tier"].eq("A_high_stability").sum())
    tier_b = int(primary["stability_tier"].eq("B_moderate_stability").sum())
    counts = {
        "step5_feature_stability_rows": int(len(stability)),
        "unique_step5_cpgs": int(stability["feature_name"].map(feature_to_cpg).nunique()),
        "manifest_assay_rows_scanned": int(manifest_assay_rows),
        "manifest_targets_requested": int(len(target_cpgs)),
        "manifest_targets_found": int(len(manifest)),
        "primary_model_feature_rows": int(len(primary_model_rows)),
        "primary_stable_cpgs": int(len(primary)),
        "tier_a": tier_a,
        "tier_b": tier_b,
        "gene_mapped_primary": int((~primary["intergenic_flag"]).sum()),
        "intergenic_primary": int(primary["intergenic_flag"].sum()),
        "manifest_qc_pass_primary": int(primary["interpretation_eligible"].sum()),
        "unique_candidate_genes": int(gene_map["gene_symbol"].nunique()) if not gene_map.empty else 0,
        "published_reference_rows": int(len(published_overlap)),
        "published_primary_exact_overlap": int(
            published_overlap["primary_stable_exact_overlap"].sum()
        ),
        "gometh_background_5mC": background_counts["5mC"],
        "gometh_background_5hmC": background_counts["5hmC"],
        "gometh_selected_5mC": selected_gometh_counts["5mC"],
        "gometh_selected_5hmC": selected_gometh_counts["5hmC"],
    }
    qc_summary = pd.DataFrame(
        [{"metric": key, "value": value} for key, value in counts.items()]
    )
    qc_summary.to_csv(output / "annotation_qc_summary_primary_clean77.csv", index=False)

    plot_stable_candidates(primary, output / "stable_cpg_annotation_primary.png")
    plot_effect_stability(primary, output / "effect_stability_primary.png")
    plot_genomic_context(primary, output / "genomic_context_primary.png")
    plot_gene_evidence(gene_evidence, output / "gene_evidence_primary.png")

    completed = datetime.now(timezone.utc).isoformat()
    summary = {
        "pipeline": {
            "name": "GSE143157 Step 6 CpG annotation and interpretation audit",
            "version": PIPELINE_VERSION,
            "completed_utc": completed,
            "elapsed_seconds": time.time() - started,
        },
        "analysis_policy": {
            "primary_model": PRIMARY_MODEL,
            "tier_A": {
                "minimum_selection_frequency": TIER_A_FREQ,
                "minimum_coefficient_direction_consistency": TIER_A_DIRECTION,
            },
            "tier_B": {
                "minimum_selection_frequency": TIER_B_FREQ,
                "minimum_coefficient_direction_consistency": TIER_B_DIRECTION,
            },
            "common_manifest_snp_maf_flag_threshold": COMMON_SNP_MAF,
            "stability_is_not_a_p_value": True,
            "automatic_pathway_enrichment_run": False,
            "pathway_reason": (
                "Naive enrichment was prohibited because EPIC gene-set tests require correction "
                "for probe-number and multi-gene bias. GOmeth-ready files were produced."
            ),
        },
        "inputs": {
            "step5_dir": str(step5),
            "step2b_dir": str(step2b),
            "manifest": str(manifest_path),
            "feature_stability_sha256": file_sha256(
                step5 / "feature_stability_primary_clean77.csv"
            ),
            "manifest_sha256": file_sha256(manifest_path),
            "h5_sha256": file_sha256(h5_path),
            "published_signals_sha256": file_sha256(published_path),
        },
        "counts": counts,
        "sources": {
            "original_GSE143157_publication": PAPER_URL,
            "bias_aware_enrichment_guidance": MISSMETHYL_URL,
            "manifest_coordinate_build": "Genome Build 37 / hg19",
        },
        "top_primary_candidates": primary.head(20)[
            [
                "source_modality",
                "cpg_id",
                "gene_symbols",
                "stability_tier",
                "selection_frequency",
                "coefficient_direction_consistency",
                "effect_direction",
                "mean_effect_msa_minus_control",
                "manifest_qc_status",
            ]
        ].to_dict(orient="records"),
    }
    (output / "step6_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (output / "RUN_CONFIGURATION.json").write_text(
        json.dumps(
            {
                "pipeline_version": PIPELINE_VERSION,
                "primary_model": PRIMARY_MODEL,
                "tier_A_frequency": TIER_A_FREQ,
                "tier_A_direction": TIER_A_DIRECTION,
                "tier_B_frequency": TIER_B_FREQ,
                "tier_B_direction": TIER_B_DIRECTION,
                "common_manifest_snp_maf": COMMON_SNP_MAF,
                "chunk_size": args.chunk_size,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    (output / "PATHWAY_ENRICHMENT_STATUS.txt").write_text(
        "NOT RUN IN THE PRIMARY PYTHON STEP\n"
        "Reason: naive gene enrichment is biased for methylation arrays.\n"
        "Bias-aware GOmeth input files were created for separate exploratory analysis.\n"
        f"Guidance: {MISSMETHYL_URL}\n",
        encoding="utf-8",
    )
    write_report(output / "STEP6_REPORT.html", primary, gene_evidence, published_overlap, summary)
    (output / "SUCCESS.txt").write_text(
        "SUCCESS\n"
        f"Completed UTC: {completed}\n"
        "Analysis: CpG annotation and descriptive stability audit\n"
        f"Primary model: {PRIMARY_MODEL}\n"
        f"Primary stable CpGs: {len(primary)} (Tier A={tier_a}, Tier B={tier_b})\n"
        f"Manifest targets found: {len(manifest)}/{len(target_cpgs)}\n"
        "Naive pathway enrichment: NOT RUN\n"
        "GOmeth-ready selected and background CpG files: CREATED\n",
        encoding="utf-8",
    )
    print("STEP 6 SUCCESS")
    print(f"Primary stable CpGs: {len(primary)} (Tier A={tier_a}, Tier B={tier_b})")
    print(f"Open report: {output / 'STEP6_REPORT.html'}")


def main() -> int:
    args = parse_args()
    try:
        run(args)
    except (ValueError, OSError, KeyError, csv.Error, json.JSONDecodeError, pd.errors.ParserError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
