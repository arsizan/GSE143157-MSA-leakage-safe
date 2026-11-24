# Leakage-safe nested validation of paired 5mC and 5hmC profiles in multiple system atrophy (GSE143157)

Analysis code and frozen outputs for:

> M. A. R. Sizan, N. Aman, and M. K. Samy, "Stable Multiple System Atrophy-Associated CpG Prioritization from 5mC and 5hmC Profiles Using Leakage-Safe Nested Validation," in *Proc. 2026 29th International Conference on Computer and Information Technology (ICCIT)*, Cox's Bazar, Bangladesh, Dec. 2026.

<!-- Keep the title above identical to the camera-ready title, and add the IEEE Xplore DOI when it is assigned. -->

## Contents

| Path | What it holds |
|---|---|
| `code/GSE143157_StepN_*/` | The script for each analysis step, its integrity check (`check_stepN_outputs.py`), its Python requirements, the Windows runner used (`RUN_STEPN_WINDOWS.cmd`) and a step guide in Bangla (`README_STEPN_BN.md`). |
| `code/GSE143157_Step8_Label_Repair/` | The display-label repair applied to the Step 8 output after the GOmeth run (labels only; no statistic changed). |
| `results/StepN_*/` | Every output of the nine analysis steps exactly as frozen for the paper: CSV tables, JSON summaries, HTML reports and PNG figures. |
| `results/*/RUN_CONFIGURATION.json`, `*_summary.json` | The parameters, seeds and input fingerprints of each step. |
| `results/Step9_Final_Audit_Paper_Package/` | Final audit: paper tables, `PAPER_NUMBERS.json`, `METHODS_RESULTS_NUMBERS.md`, claim guardrails, figures, and `FINAL_INPUT_FILE_MANIFEST.csv` (SHA-256 of the 19 files behind the paper's numbers). |

Each step's check script passes when run against the matching `results/` folder for Steps 2–9 (for example `python code/GSE143157_Step5_Molecular_Benchmark_Pipeline/check_step5_outputs.py results/Step5_Molecular_Benchmark_Primary77 --step3-dir results/Step3_Confounders_CV`). The checks for Steps 1 and 2B also need the HDF5 matrices listed under *Not included*.

## Input data (not redistributed)

GEO series [GSE143157](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143157) (Rydbirk et al., *Acta Neuropathol. Commun.*, 2020): post-mortem prefrontal cortex from 41 MSA cases and 37 controls, each with one bisulfite (BS) and one oxidative-bisulfite (oxBS) Illumina MethylationEPIC array.

- `GSE143157_BMIQ_normalized_beta.txt`: BMIQ-normalized beta matrix (about 2.35 GB; 731,661 CpGs, 156 beta and 156 detection-P columns)
- `GSE143157_family.soft.gz`: sample metadata

## Pipeline

| Step | Folder in `results/` | What it did |
|---|---|---|
| 1 | `Step1_results_full78` | Donor-aligned 5mC (oxBS beta) and 5hmC (BS − oxBS; negative values set to missing; CpGs with > 20% missing removed) views; detection P < 0.01 summarized without masking beta values. |
| 2 | `Step2_PCA_Outliers_full78` | PCA robust distance, mean sample correlation, 5hmC missingness and negative-difference QC for all 78 donors; flags triggered review, not automatic exclusion. |
| 2B | `Step2B_Cohorts` | Frozen cohorts: primary Clean-77 (MSA 30 excluded), and Full-78 and Clean-76 (MSA 30 and Control 11 excluded) as overlapping sensitivity cohorts. |
| 3 | `Step3_Confounders_CV` | Covariate audit (age, sex, brain bank, PMI) and frozen donor folds: 5 outer folds × 10 repeats, 4 inner folds, seed 143157; no CpG values were used to build the folds. |
| 4 | `Step4_Covariate_Benchmark` | Covariate-only benchmarks: age only; sex, brain bank and PMI; all four covariates. |
| 5 | `Step5_Molecular_Benchmark_Primary77` | L2-regularized logistic regression (scikit-learn 1.8.0, liblinear) for 5mC, 5hmC and combined views: molecular-only, molecular + covariates, covariate-residualized. Inside training folds only: missingness and variance filtering, median imputation, scaling, a label-blind 5,000-CpG variance screen, Welch-statistic selection of 25, 50 or 100 CpGs per modality, and tuning of C ∈ {0.01, 0.1, 1, 10} and class weights by mean inner ROC-AUC. Donor-mean out-of-fold metrics with 2,000 stratified donor-bootstrap replicates; Benjamini–Hochberg FDR. |
| 6 | `Step6_Annotation_Interpretation_Primary77` | Stability tiers (A: selection frequency ≥ 0.80 and direction consistency ≥ 0.90; B: ≥ 0.60 and ≥ 0.80), EPIC hg19 manifest annotation, overlap with the source study's reported CpGs. |
| 7 | `Step7_Molecular_Sensitivity` | The locked molecular + covariates workflow on Full-78 and Clean-76; a candidate is supported when its selection frequency is ≥ 0.60 with the primary effect direction in both cohorts. |
| 8 | `Step8_BiasAware_Pathway_Enrichment` | GOmeth (R 4.6.1, missMethyl 1.46.0; probe-number and multi-gene corrections) GO and KEGG enrichment of the seven manifest-passing stable CpGs per modality against modality-specific backgrounds. A later presentation-label repair changed display labels only (`label_repair_audit.csv`; the pre-repair files are kept in `_backup_before_label_repair_20260821_024317`). |
| 9 | `Step9_Final_Audit_Paper_Package` | Final freeze, audit checks, tables, figures and paper numbers. |

## Running the code

Each step was run on Windows with its `RUN_STEPN_WINDOWS.cmd`. The runner creates a private Python environment (`.venv`), installs `requirements_stepN.txt`, runs the step and then its check script. The runners point to the original folders under `D:\Paper\MSA_Research\02_Epigenetics\GSE143157\`, so change those paths at the top of each runner before use. Steps run in order (1, 2, 2B, 3, 4, 5, 6, 7, 8, 9), and each step reads the outputs of the steps before it. Step 8 is an R script (`run_gse143157_step8.R`) that needs R 4.6.1 with missMethyl 1.46.0. Steps 5 and 7, the nested cross-validation runs, take the longest.

## Where the paper's numbers come from

| Paper item | File |
|---|---|
| Molecular models (Table II, Fig. 2) | `results/Step5_Molecular_Benchmark_Primary77/aggregate_metrics_primary_clean77.csv` |
| Age-only and covariate benchmarks | `results/Step4_Covariate_Benchmark/aggregate_metrics_primary_clean77.csv` |
| Paired ROC-AUC differences, including against age alone | `results/Step5_Molecular_Benchmark_Primary77/auc_comparisons_primary_clean77.csv` |
| Full-78 and Clean-76 comparisons | `results/Step7_Molecular_Sensitivity/sensitivity_*/auc_comparisons_*.csv` |
| 20 stable CpGs, tiers and cohort support | `results/Step9_Final_Audit_Paper_Package/TABLE3_robust_cpg_candidates.csv` |
| Stability after covariate residualization | `results/Step5_Molecular_Benchmark_Primary77/feature_stability_primary_clean77.csv` (model `combined__residualized_molecular`) |
| Overlap with the source study's FDR < 0.05 CpGs | `results/Step6_Annotation_Interpretation_Primary77/published_signal_overlap_primary_clean77.csv` |
| GO and KEGG results | `results/Step8_BiasAware_Pathway_Enrichment/enrichment_summary.csv` and `gometh_*.csv` |
| Headline numbers | `results/Step9_Final_Audit_Paper_Package/PAPER_NUMBERS.json` |

## Not included

- The processed HDF5 matrices `Step1_results_full78/GSE143157_step1_matrices.h5` and `Step2B_Cohorts/GSE143157_primary_clean77.h5` (about 300 MB each), which exceed GitHub's 100 MB file limit.
- Step 5 and Step 7 resume checkpoints (`_checkpoints/*.pkl`).
- The private Python environments (`.venv`) that the runners create.

## Citation

See `CITATION.cff`. If you use this code or these outputs, please also cite the original dataset: R. Rydbirk et al., "Epigenetic modulation of AREL1 and increased HLA expression in brains of multiple system atrophy patients," *Acta Neuropathol. Commun.*, vol. 8, Art. no. 29, 2020, doi: 10.1186/s40478-020-00908-7.
