# GSE143157 — Step 4: covariate-only benchmark

## কেন এই ধাপ প্রয়োজন

Primary cohort-এ Control-এর mean age 72.97 বছর এবং MSA-এর mean age 66.00 বছর। তাই সরাসরি methylation model ভালো score করলে সেটি disease signal নাকি age signal—তা আলাদা করে বোঝা দরকার।

Step 4 তিনটি model চালাবে:

1. `age_only` — শুধু age
2. `no_age_covariates` — sex + brain bank + PMI
3. `all_covariates` — age + sex + brain bank + PMI

এগুলো molecular biomarker model নয়। এগুলো হলো confounding benchmark। পরবর্তী 5mC/5hmC model-কে `all_covariates` model-এর চেয়ে অতিরিক্ত predictive value দেখাতে হবে।

## Leakage-safe নিয়ম

- Step 3-এর exact outer ও inner donor folds ব্যবহার হবে।
- Continuous imputation ও scaling প্রতিটি inner/outer training set-এ আলাদাভাবে fit হবে।
- C এবং class-weight শুধু inner folds দিয়ে নির্বাচন হবে।
- Outer-test donor hyperparameter selection বা threshold tuning-এ ব্যবহার হবে না।
- প্রতিটি donor প্রতি model ও repeat-এ ঠিক একবার outer-test prediction পাবে।
- Prediction threshold স্থির `0.5`।

## কোথায় extract করবেন

`GSE143157_Step4_Covariate_Benchmark_Pipeline.zip` extract করুন:

```text
D:\Paper\MSA_Research
```

তারপর folder হবে:

```text
D:\Paper\MSA_Research\GSE143157_Step4_Covariate_Benchmark_Pipeline
```

## চালানোর নিয়ম

Folder খুলে double-click করুন:

```text
RUN_STEP4_WINDOWS.cmd
```

Black window বন্ধ করবেন না। অনেকগুলো ছোট inner-CV model fit হবে, তাই কয়েক মিনিট লাগতে পারে। শেষ হলে browser-এ report খুলবে এবং terminal-এ দেখাবে:

```text
STEP 4 FINISHED SUCCESSFULLY
STEP 4 INTEGRITY CHECK PASSED
```

## Command Prompt দিয়ে চালাতে চাইলে

```bat
cd /d D:\Paper\MSA_Research\GSE143157_Step4_Covariate_Benchmark_Pipeline
RUN_STEP4_WINDOWS.cmd
```

## Input

```text
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step2B_Cohorts
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step3_Confounders_CV
```

## Output

```text
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step4_Covariate_Benchmark
```

প্রধান output:

- `STEP4_REPORT.html`
- `step4_summary.json`
- `aggregate_metrics_primary_clean77.csv`
- `repeat_metrics_primary_clean77.csv`
- `donor_mean_predictions_primary_clean77.csv`
- `pairwise_auc_comparisons_primary_clean77.csv`
- `hyperparameter_choices_primary_clean77.csv`
- `coefficient_stability_primary_clean77.csv`
- Full-78 ও optional-76 cohort-এর একই ধরনের files
- `roc_primary.png`, `pr_primary.png`
- `metrics_primary.png`, `coefficients_primary.png`
- `SUCCESS.txt`

## Result পাওয়ার পরে কী দেখব

- Age-only ROC-AUC কত
- Sex/bank/PMI age ছাড়া কতটা predict করে
- All-covariate model age-only-এর চেয়ে সত্যিই ভালো কি না
- Confidence interval chance level অতিক্রম করছে কি না
- Age coefficient 50 outer fits-এ একই direction-এ আছে কি না

## Run শেষে আমাকে যা পাঠাবেন

```text
SUCCESS.txt
step4_summary.json
aggregate_metrics_primary_clean77.csv
repeat_metrics_primary_clean77.csv
pairwise_auc_comparisons_primary_clean77.csv
coefficient_stability_primary_clean77.csv
```

আর figures:

```text
roc_primary.png
pr_primary.png
metrics_primary.png
coefficients_primary.png
```

বড় prediction/tuning CSV upload করার প্রয়োজন নেই, যদি না পরে চাই। Error হলে black window-এর সম্পূর্ণ screenshot পাঠাবেন।
