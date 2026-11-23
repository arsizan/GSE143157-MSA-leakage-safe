# GSE143157 locked methods and results numbers

This file is a data-backed drafting aid. It does not replace author review.

## Cohort and preprocessing

The analysis used 78 paired donors (41 MSA and 37 controls) with matched BS and oxBS EPIC profiles. We represented 5mC by the oxBS beta value and estimated 5hmC as BS minus oxBS. Negative 5hmC differences were treated as missing, and CpGs with more than 20% missing 5hmC values were removed. This retained 405,408 5hmC CpGs; the 5mC matrix contained 731,661 CpGs.

The primary cohort contained 77 donors (40 MSA and 37 controls) after excluding MSA 30 based on convergent technical QC evidence. Full-78 and Clean-76 were prespecified overlapping-cohort sensitivity analyses, not independent external validation cohorts.

## Validation design

Models used five outer folds repeated ten times and four inner folds. Missing-value handling, variance filtering, imputation, scaling, supervised CpG selection, residualization, and hyperparameter tuning were fitted using training donors only. Age, sex, brain bank, and post-mortem interval were prespecified covariates.

## Primary discrimination

In Primary Clean-77, donor-mean out-of-fold ROC-AUC was 0.693 (95% donor-bootstrap CI 0.560-0.809) for all covariates, 0.632 for 5mC plus covariates, 0.663 for 5hmC plus covariates, and 0.653 for combined molecular plus covariates. The combined-versus-all-covariates ROC-AUC difference was -0.041 (95% paired donor-bootstrap CI -0.205 to 0.124; within-family BH-adjusted value 0.703). Thus, incremental molecular discrimination beyond covariates was not established.

## Sensitivity and candidate stability

Combined plus covariate ROC-AUC was 0.705 in Full-78 and 0.685 in Clean-76. Of 20 primary stable CpG candidates, 15 met the prespecified selection-frequency and effect-direction rule in both sensitivity cohorts. All 20 retained the primary effect direction in both overlapping sensitivity cohorts. These analyses support robustness to the cohort decision but do not constitute external replication.

## Pathway analysis

Bias-aware exploratory enrichment used missMethyl::gometh with EPIC annotation, separate tested-probe backgrounds for 5mC and 5hmC, probe-number prior probabilities, equivalent CpG counts, and fractional counting for multi-gene probes. The input comprised seven manifest-audited stable CpGs per modality. The number of FDR-significant terms was 0. Nominally ranked pathways must be described as exploratory rather than confirmed mechanisms.

## Mandatory limitations

- Single public post-mortem dataset.
- Small donor sample and wide uncertainty intervals.
- Overlapping sensitivity cohorts are not external validation.
- Selection frequency is not a p-value.
- Stable CpGs are exploratory associations, not validated diagnostic biomarkers.
- Observational data do not establish causality.
