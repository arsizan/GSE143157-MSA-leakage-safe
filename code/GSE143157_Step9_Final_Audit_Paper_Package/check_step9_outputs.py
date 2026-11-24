from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path


REQUIRED = [
    "SUCCESS.txt",
    "FINAL_FREEZE.json",
    "RUN_CONFIGURATION.json",
    "PAPER_NUMBERS.json",
    "TABLE1_cohorts.csv",
    "TABLE2_model_performance.csv",
    "TABLE3_robust_cpg_candidates.csv",
    "TABLE4_bias_aware_pathways.csv",
    "SUPPLEMENT_all_top_pathways.csv",
    "FINAL_METHODS_AUDIT.csv",
    "FINAL_INPUT_FILE_MANIFEST.csv",
    "FINAL_ANALYSIS_DECISION_LOG.csv",
    "PAPER_CLAIMS_GUARDRAILS.csv",
    "METHODS_RESULTS_NUMBERS.md",
    "STEP9_REPORT.html",
    "figures/Figure0_analysis_workflow.png",
    "figures/Figure1_primary_molecular_benchmarks.png",
    "figures/Figure2_sensitivity_performance.png",
    "figures/Figure3_candidate_robustness.png",
]


def fail(message: str) -> None:
    raise SystemExit(f"STEP 9 INTEGRITY CHECK FAILED: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python check_step9_outputs.py <Step9 output directory>")
    output = Path(sys.argv[1]).resolve()
    if not output.is_dir():
        fail("output directory does not exist")
    missing = [name for name in REQUIRED if not (output / name).is_file()]
    if missing:
        fail("missing files: " + ", ".join(missing))

    success = (output / "SUCCESS.txt").read_text(encoding="utf-8")
    if "FINAL_VERIFIED_FREEZE" not in success:
        fail("SUCCESS.txt is not a final verified freeze")

    cohorts = read_csv(output / "TABLE1_cohorts.csv")
    if len(cohorts) != 3 or {int(row["donors_n"]) for row in cohorts} != {76, 77, 78}:
        fail("cohort table does not contain the locked 76/77/78 donor cohorts")

    performance = read_csv(output / "TABLE2_model_performance.csv")
    if len(performance) != 12:
        fail("performance table must contain four models across three cohorts")
    if {row["model"] for row in performance} != {
        "All covariates", "5mC + covariates", "5hmC + covariates", "Combined + covariates"
    }:
        fail("performance model set is not locked")
    for row in performance:
        for column in ("roc_auc", "roc_auc_ci_low", "roc_auc_ci_high", "pr_auc", "balanced_accuracy"):
            value = float(row[column])
            if not 0 <= value <= 1:
                fail(f"{column} outside [0,1]")

    candidates = read_csv(output / "TABLE3_robust_cpg_candidates.csv")
    if len(candidates) != 20:
        fail("candidate table must contain 20 primary stable CpGs")
    supported = sum(row["sensitivity_supported_both"].lower() == "true" for row in candidates)
    if supported != 15:
        fail(f"expected 15 sensitivity-supported candidates, observed {supported}")
    if not all(row["effect_direction_matches_both"].lower() == "true" for row in candidates):
        fail("at least one candidate effect direction does not match both sensitivity cohorts")

    audit = read_csv(output / "FINAL_METHODS_AUDIT.csv")
    if not audit or any(row["status"] != "PASS" for row in audit):
        fail("final methods audit contains a failure")

    numbers = json.loads((output / "PAPER_NUMBERS.json").read_text(encoding="utf-8"))
    if numbers["dataset"]["primary_donors"] != 77:
        fail("PAPER_NUMBERS primary donor count is incorrect")
    if numbers["candidate_stability"]["supported_in_both_sensitivity_cohorts"] != 15:
        fail("PAPER_NUMBERS candidate support count is incorrect")
    if numbers["combined_vs_all_covariates_primary"]["bootstrap_ci_low"] > 0 or numbers["combined_vs_all_covariates_primary"]["bootstrap_ci_high"] < 0:
        fail("central comparison CI unexpectedly excludes zero; inspect source comparison")

    freeze = json.loads((output / "FINAL_FREEZE.json").read_text(encoding="utf-8"))
    records = freeze.get("files", [])
    if not records:
        fail("FINAL_FREEZE.json contains no file records")
    for record in records:
        path = output / record["relative_path"]
        if not path.is_file():
            fail(f"frozen file is missing: {record['relative_path']}")
        if path.stat().st_size != int(record["size_bytes"]):
            fail(f"frozen file size changed: {record['relative_path']}")
        if sha256(path) != record["sha256"]:
            fail(f"frozen file hash changed: {record['relative_path']}")
    seed = "\n".join(f"{item['relative_path']}:{item['sha256']}" for item in records).encode("utf-8")
    if hashlib.sha256(seed).hexdigest() != freeze.get("freeze_id"):
        fail("freeze ID does not match the frozen file list")

    for image in (output / "figures").glob("*.png"):
        if image.stat().st_size < 10_000:
            fail(f"figure appears incomplete: {image.name}")

    print("STEP 9 INTEGRITY CHECK PASSED")


if __name__ == "__main__":
    main()

