# GSE143157 — Step 2B: analysis cohort তৈরি

## এই ধাপে কী হবে

এই pipeline তিনটি donor list তৈরি করবে:

1. **Primary Clean-77:** 37 Control + 40 MSA। MSA 30 বাদ থাকবে।
2. **Sensitivity Full-78:** সব 78 donor থাকবে। Final result এই cohort-এ আবার পরীক্ষা করতে হবে।
3. **Optional stress-test 76:** MSA 30 এবং Control 11 বাদ থাকবে। এটি primary dataset নয়।

Primary analysis-এর জন্য নতুন `GSE143157_primary_clean77.h5` তৈরি হবে। বড় Step 1 HDF5 পুরোটা RAM-এ load হবে না; ছোট row-chunk ধরে copy হবে। Original Step 1 ও Step 2 files বদলানো বা overwrite করা হবে না।

## কেন MSA 30 primary cohort-এ থাকবে না

MSA 30-এর 5mC মোটামুটি স্বাভাবিক হলেও BS-dependent 5hmC-তে চারটি evidence একসঙ্গে পাওয়া গেছে:

- Detection-P failure
- অত্যধিক negative `BS − oxBS`
- filtering-এর পরে অত্যধিক missing 5hmC
- extreme 5hmC PCA outlier

Paired 5mC-only, 5hmC-only এবং combined model-এর ন্যায্য তুলনার জন্য তিনটিতেই একই 77 donor ব্যবহার করা হবে। তাই primary HDF5 থেকে MSA 30-এর donor column দুই modality থেকেই বাদ যাবে।

Control 11 PCA-তে outlier হলেও Detection-P/missingness technical failure নেই। তাই সে primary Clean-77-এ থাকবে। শুধু optional 76-donor stress-test-এ তাকে বাদ দিয়ে result কতটা বদলায় দেখা যাবে।

## কোথায় extract করবেন

`GSE143157_Step2B_Cohorts_Pipeline.zip` extract করুন:

```text
D:\Paper\MSA_Research
```

Extract হওয়ার পরে folder হওয়া উচিত:

```text
D:\Paper\MSA_Research\GSE143157_Step2B_Cohorts_Pipeline
```

## সবচেয়ে সহজে চালানোর নিয়ম

উপরের folder খুলে এই file-এ double-click করুন:

```text
RUN_STEP2B_WINDOWS.cmd
```

Black window চলার সময় বন্ধ করবেন না এবং computer sleep-এ যেতে দেবেন না। HDF5 copy ও compression-এর জন্য কয়েক মিনিট লাগতে পারে। প্রথমবার Python packages install করার জন্য internet connection লাগতে পারে।

শেষে browser-এ report নিজে থেকে খুলবে এবং black window-তে দেখাবে:

```text
STEP 2B FINISHED SUCCESSFULLY
STEP 2B INTEGRITY CHECK PASSED
```

## Command Prompt দিয়ে চালাতে চাইলে

```bat
cd /d D:\Paper\MSA_Research\GSE143157_Step2B_Cohorts_Pipeline
RUN_STEP2B_WINDOWS.cmd
```

PowerShell বা Python editor-এ আলাদা আলাদা code paste করার প্রয়োজন নেই।

## Input folders

```text
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step1_results_full78
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step2_PCA_Outliers_full78
```

## Output folder

```text
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step2B_Cohorts
```

গুরুত্বপূর্ণ output:

- `STEP2B_REPORT.html` — beginner-friendly cohort report
- `PRIMARY_COHORT.txt` — primary decision-এর ছোট summary
- `cohort_decision.json` — paper/audit-এর authoritative decision record
- `cohort_manifest.csv` — 78 donor-এর কোন cohort-এ inclusion আছে
- `metadata_donors_primary_clean77.csv` — Clean-77 donor metadata
- `metadata_donors_sensitivity_full78.csv` — Full-78 sensitivity metadata
- `metadata_donors_sensitivity_no_MSA30_Control11_76.csv` — optional stress-test metadata
- `metadata_arrays_*.csv` — প্রতিটি cohort-এর paired array metadata
- `GSE143157_primary_clean77.h5` — Step 3-এর primary paired matrices
- `SUCCESS.txt` — completed-run marker

## Run শেষ হলে আমাকে যা পাঠাবেন

```text
SUCCESS.txt
PRIMARY_COHORT.txt
cohort_decision.json
cohort_manifest.csv
metadata_donors_primary_clean77.csv
```

`GSE143157_primary_clean77.h5` upload করার দরকার নেই। Error হলে black window-এর সম্পূর্ণ screenshot পাঠাবেন।

## গুরুত্বপূর্ণ নিয়ম

- Step 2-এর top-variable CpG list পরবর্তী ML feature selection হিসেবে ব্যবহার করা যাবে না।
- Supervised feature selection প্রতিটি training fold-এর ভেতরে করতে হবে।
- Main result Clean-77-এ হবে এবং final finding Full-78-এ sensitivity analysis হিসেবে পুনরায় পরীক্ষা করতে হবে।
- Original study-এর 63-donor cohort exact excluded-donor list না পাওয়া পর্যন্ত অনুমান করে তৈরি করা যাবে না।
