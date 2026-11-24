# GSE143157 — Step 5: Primary molecular benchmark

## এই ধাপে কী হবে

Primary Clean-77 cohort ব্যবহার করে নয়টি model configuration পরীক্ষা হবে:

```text
3 modality: 5mC, 5hmC, combined
3 view: molecular only, molecular + covariates, covariate-residualized molecular
```

এখানে `combined` model প্রতি modality থেকে সমান সংখ্যক CpG নেয়। যেমন feature count 50 হলে 50টি 5mC + 50টি 5hmC CpG ব্যবহৃত হবে।

Step 4-এর age-only এবং all-covariate predictions-এর সঙ্গে paired comparison হবে। প্রধান প্রশ্ন:

> Molecular data কি age/sex/brain-bank/PMI baseline-এর বাইরে অতিরিক্ত prediction দিতে পারে?

## Leakage কীভাবে আটকানো হয়েছে

- Step 3-এর একই frozen donor folds ব্যবহার হয়।
- Outer-test donor variance ranking, missingness rule, imputation, scaling বা CpG selection-এ অংশ নেয় না।
- Diagnosis দেখে CpG নির্বাচন প্রতিটি inner-training fold-এর ভেতরে নতুন করে হয়।
- Hyperparameter শুধু inner validation দিয়ে নির্বাচন হয়।
- Outer-test donor একবার শুধু final prediction-এর জন্য ব্যবহৃত হয়।
- 5mC beta values অপরিবর্তিত থাকে; 5hmC non-negative BS−oxBS values ব্যবহৃত হয়।

## Memory-safe ব্যবস্থা

পুরো HDF5 matrix RAM-এ load হয় না। প্রতি outer fold-এ 5,000-row chunk করে matrix পড়া হয়। Outer-training donor দিয়ে সবচেয়ে variable 5,000 CpG-এর candidate pool তৈরি হয়। এরপর inner-training donor দিয়ে supervised selection হয়।

## কোথায় extract করবেন

ZIP extract করুন:

```text
D:\Paper\MSA_Research
```

তারপর folder হবে:

```text
D:\Paper\MSA_Research\GSE143157_Step5_Molecular_Benchmark_Pipeline
```

## কীভাবে চালাবেন

Folder খুলে double-click করুন:

```text
RUN_STEP5_WINDOWS.cmd
```

Black window বন্ধ করবেন না। Computer charger-এ রাখুন এবং sleep বন্ধ রাখুন। Computer অনুযায়ী সাধারণত প্রায় 15–60 মিনিট লাগতে পারে; ধীর computer-এ আরও বেশি সময় লাগতে পারে।

প্রতি outer fold শেষ হলে checkpoint save হবে। মাঝখানে বন্ধ হয়ে গেলে একই `RUN_STEP5_WINDOWS.cmd` আবার double-click করলে completed folds থেকে resume করবে।

শেষে দেখাবে:

```text
STEP 5 SUCCESS
STEP 5 INTEGRITY CHECK PASSED
STEP 5 FINISHED SUCCESSFULLY
```

## Command Prompt দিয়ে চালাতে চাইলে

```bat
cd /d D:\Paper\MSA_Research\GSE143157_Step5_Molecular_Benchmark_Pipeline
RUN_STEP5_WINDOWS.cmd
```

## Input folders

```text
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step2B_Cohorts
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step3_Confounders_CV
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step4_Covariate_Benchmark
```

## Result folder

```text
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step5_Molecular_Benchmark_Primary77
```

প্রধান output:

- `STEP5_REPORT.html`
- `step5_summary.json`
- `aggregate_metrics_primary_clean77.csv`
- `auc_comparisons_primary_clean77.csv`
- `repeat_metrics_primary_clean77.csv`
- `donor_mean_predictions_primary_clean77.csv`
- `hyperparameter_choices_primary_clean77.csv`
- `feature_stability_primary_clean77.csv`
- `selected_features_primary_clean77.csv`
- `performance_primary.png`
- `roc_primary.png`
- `incremental_auc_primary.png`
- `feature_stability_primary.png`
- `SUCCESS.txt`

`_checkpoints` folder delete করবেন না যতক্ষণ full run সফল না হয়।

## Run শেষে আমাকে যা পাঠাবেন

```text
SUCCESS.txt
step5_summary.json
aggregate_metrics_primary_clean77.csv
auc_comparisons_primary_clean77.csv
repeat_metrics_primary_clean77.csv
hyperparameter_choices_primary_clean77.csv
feature_stability_primary_clean77.csv
```

এবং চারটি PNG figure পাঠাবেন। বড় prediction, tuning ও selected-feature files পরে দরকার হলে চাইব।

## গুরুত্বপূর্ণ সতর্কতা

Code-এর `--smoke` option শুধু software test-এর জন্য। Smoke output-এর numbers paper-এ ব্যবহার করা যাবে না। Windows runner নিজে full research run চালায়।
