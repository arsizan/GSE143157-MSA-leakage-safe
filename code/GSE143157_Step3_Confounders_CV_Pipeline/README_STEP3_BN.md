# GSE143157 — Step 3: confounder audit এবং leakage-safe CV

## এই ধাপে কী হবে

Step 3 এখনো disease-classification model train করবে না। এটি model training-এর নিয়ম ও donor folds আগে থেকে স্থির করবে। এতে ফল দেখার পরে সুবিধামতো sample split পরিবর্তন করা যাবে না।

Pipeline তিনটি কাজ করবে:

1. Age, sex, brain bank এবং post-mortem interval (PMI)-এর MSA–Control imbalance effect size ও permutation test দিয়ে মাপবে।
2. Primary Clean-77, Full-78 sensitivity এবং optional 76-donor stress-test-এর জন্য fixed nested cross-validation folds তৈরি করবে।
3. কোন preprocessing training fold-এর ভেতরে করতে হবে তার immutable policy লিখবে।

## Frozen cross-validation design

```text
Outer CV: 5 folds × 10 repeats
Inner CV: প্রতিটি outer-training set-এর ভেতরে 4 folds
Random seed: 143157
```

- Outer test fold শুধু unbiased performance measurement-এর জন্য।
- Inner folds feature count ও hyperparameter নির্বাচন করবে।
- একই donor folds 5mC-only, 5hmC-only এবং combined model-এ ব্যবহার করতে হবে।
- Fold তৈরি করতে কোনো CpG value পড়া হবে না। Diagnosis label দিয়ে stratification এবং age/sex/brain-bank/PMI দিয়ে balance করা হবে।

## কোথায় extract করবেন

`GSE143157_Step3_Confounders_CV_Pipeline.zip` extract করুন:

```text
D:\Paper\MSA_Research
```

Folder হবে:

```text
D:\Paper\MSA_Research\GSE143157_Step3_Confounders_CV_Pipeline
```

## সবচেয়ে সহজে চালানোর নিয়ম

Folder খুলে double-click করুন:

```text
RUN_STEP3_WINDOWS.cmd
```

প্রথমবার package install-এর জন্য internet connection লাগতে পারে। শেষ হলে report browser-এ খুলবে এবং black window-তে দেখাবে:

```text
STEP 3 FINISHED SUCCESSFULLY
STEP 3 INTEGRITY CHECK PASSED
```

## Command Prompt দিয়ে চালাতে চাইলে

```bat
cd /d D:\Paper\MSA_Research\GSE143157_Step3_Confounders_CV_Pipeline
RUN_STEP3_WINDOWS.cmd
```

## Input

```text
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step2B_Cohorts
```

## Output

```text
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step3_Confounders_CV
```

প্রধান output:

- `STEP3_REPORT.html` — সহজ summary ও figures
- `step3_summary.json` — authoritative parameters/results
- `model_input_policy.json` — leakage prevention rules
- `covariate_audit_primary_clean77.csv` — primary confounder results
- `outer_folds_primary_clean77.csv` — প্রতিটি donor কোন repeat-এর কোন outer-test fold-এ যাবে
- `inner_folds_primary_clean77.csv` — প্রতিটি outer-training set-এর inner validation folds
- `outer_fold_balance_primary_clean77.csv` — fold-wise covariate balance
- Full-78 ও 76 cohort-এর জন্য একই ধরনের files
- `covariate_balance_primary.png`
- `outer_fold_balance_primary.png`
- `SUCCESS.txt`

## Leakage prevention-এর সহজ অর্থ

প্রতিটি outer fold-এ:

```text
Outer training donors
        ↓
Missingness filter / imputation / variance filter
        ↓
Feature selection / scaling / hyperparameter tuning
        ↓
Locked fitted pipeline
        ↓
Outer test donors-এ একবার prediction
```

Outer-test donor-এর information দিয়ে median, feature list, scaler বা model parameter বানানো যাবে না। Step 2-এর top-variable CpG list supervised ML-এ সরাসরি ব্যবহার করা যাবে না।

## Run শেষে আমাকে যা দেবেন

```text
SUCCESS.txt
step3_summary.json
model_input_policy.json
covariate_audit_primary_clean77.csv
outer_fold_balance_primary_clean77.csv
```

আর এই দুইটি figure:

```text
covariate_balance_primary.png
outer_fold_balance_primary.png
```

বড় outer/inner fold CSV upload করার দরকার নেই, যদি না পরে চাই। Error হলে black window-এর সম্পূর্ণ screenshot পাঠাবেন।
