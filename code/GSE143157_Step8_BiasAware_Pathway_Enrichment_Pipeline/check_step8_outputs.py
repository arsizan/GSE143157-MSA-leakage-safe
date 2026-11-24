from __future__ import annotations

import csv
import gzip
import json
import math
import sys
from pathlib import Path


REQUIRED = [
    "SUCCESS.txt",
    "step8_summary.json",
    "RUN_CONFIGURATION.json",
    "gometh_input_qc.csv",
    "enrichment_summary.csv",
    "top_terms_all.csv",
    "significant_terms_FDR05.csv",
    "gometh_GO_5mC.csv",
    "gometh_KEGG_5mC.csv",
    "gometh_GO_5hmC.csv",
    "gometh_KEGG_5hmC.csv",
    "top_enrichment_5mC.png",
    "top_enrichment_5hmC.png",
    "package_versions.csv",
    "R_SESSION_INFO.txt",
    "step8_console_log.txt",
    "STEP8_REPORT.html",
]


def fail(message: str) -> None:
    raise SystemExit(f"STEP 8 INTEGRITY CHECK FAILED: {message}")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(value: str, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        fail(f"{label} is not numeric: {value!r}")
        raise exc
    if not math.isfinite(result):
        fail(f"{label} is not finite")
    return result


def count_gzip_lines(path: Path) -> int:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip())


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: python check_step8_outputs.py <Step6 directory> <Step8 directory>")

    step6 = Path(sys.argv[1]).resolve()
    output = Path(sys.argv[2]).resolve()
    if not step6.is_dir() or not output.is_dir():
        fail("Step 6 or Step 8 directory does not exist")

    missing = [name for name in REQUIRED if not (output / name).is_file()]
    if missing:
        fail("missing output files: " + ", ".join(missing))

    success = (output / "SUCCESS.txt").read_text(encoding="utf-8")
    if "BIAS_AWARE_EXPLORATORY_ENRICHMENT" not in success:
        fail("SUCCESS.txt does not identify the locked analysis mode")

    summary = json.loads((output / "step8_summary.json").read_text(encoding="utf-8"))
    config = json.loads((output / "RUN_CONFIGURATION.json").read_text(encoding="utf-8"))
    if config.get("array_type") != "EPIC":
        fail("array_type is not EPIC")
    for key in ("prior_prob", "equivalent_cpg_counts", "fractional_multi_gene_counts"):
        if config.get(key) is not True:
            fail(f"{key} was not locked to TRUE")
    if config.get("selected_cpgs_are_ewas_significant") is not False:
        fail("EWAS-significance guardrail is missing")

    input_qc = read_csv(output / "gometh_input_qc.csv")
    if {row["modality"] for row in input_qc} != {"5mC", "5hmC"}:
        fail("input QC does not contain exactly 5mC and 5hmC")
    for row in input_qc:
        modality = row["modality"]
        selected_path = step6 / f"gometh_selected_{modality}.txt.gz"
        background_path = step6 / f"gometh_background_{modality}.txt.gz"
        if not selected_path.is_file() or not background_path.is_file():
            fail(f"Step 6 GOmeth inputs are missing for {modality}")
        if int(row["selected_cpgs"]) != count_gzip_lines(selected_path):
            fail(f"selected CpG count mismatch for {modality}")
        if int(row["background_cpgs"]) != count_gzip_lines(background_path):
            fail(f"background CpG count mismatch for {modality}")
        if row["selected_subset_of_background"].lower() != "true":
            fail(f"selected CpGs are not a subset of the {modality} background")

    enrichment = read_csv(output / "enrichment_summary.csv")
    expected_pairs = {(m, c) for m in ("5mC", "5hmC") for c in ("GO", "KEGG")}
    observed_pairs = {(row["modality"], row["collection"]) for row in enrichment}
    if observed_pairs != expected_pairs:
        fail("enrichment_summary.csv does not contain all four modality/collection pairs")

    result_lookup: dict[tuple[str, str], list[dict[str, str]]] = {}
    for modality, collection in expected_pairs:
        path = output / f"gometh_{collection}_{modality}.csv"
        rows = read_csv(path)
        if not rows:
            fail(f"{path.name} contains no tested terms")
        for row in rows:
            p = as_float(row["P.DE"], f"{path.name} P.DE")
            fdr = as_float(row["FDR"], f"{path.name} FDR")
            if not (0.0 <= p <= 1.0 and 0.0 <= fdr <= 1.0):
                fail(f"{path.name} contains P/FDR outside [0,1]")
        result_lookup[(modality, collection)] = rows

    for row in enrichment:
        key = (row["modality"], row["collection"])
        rows = result_lookup[key]
        if int(row["tested_terms"]) != len(rows):
            fail(f"tested-term count mismatch for {key}")
        expected_sig = sum(float(item["FDR"]) < 0.05 for item in rows)
        if int(row["fdr_lt_0_05"]) != expected_sig:
            fail(f"FDR-significant count mismatch for {key}")

    significant = read_csv(output / "significant_terms_FDR05.csv")
    if any(float(row["FDR"]) >= 0.05 for row in significant):
        fail("significant_terms_FDR05.csv contains FDR >= 0.05")
    expected_total = sum(
        float(row["FDR"]) < 0.05
        for rows in result_lookup.values()
        for row in rows
    )
    if len(significant) != expected_total:
        fail("combined significant-term count is inconsistent")
    if int(summary.get("total_fdr_significant_terms", -1)) != expected_total:
        fail("step8_summary.json significant-term count is inconsistent")

    for image in ("top_enrichment_5mC.png", "top_enrichment_5hmC.png"):
        if (output / image).stat().st_size < 10_000:
            fail(f"{image} appears incomplete")

    print("STEP 8 INTEGRITY CHECK PASSED")


if __name__ == "__main__":
    main()

