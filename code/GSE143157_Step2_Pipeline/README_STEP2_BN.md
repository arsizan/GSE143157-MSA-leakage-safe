# GSE143157 — Step 2: PCA, correlation এবং outlier review

## Step 2-এর উদ্দেশ্য

Step 2 কোনো disease-classification model বানাবে না এবং কোনো donor delete করবে না। এটি যাচাই করবে:

- 5mC data-তে কোনো donor অন্যদের থেকে অস্বাভাবিক দূরে কি না
- 5hmC data-তে একই সমস্যা আছে কি না
- sample-to-sample correlation কম কি না
- age, sex বা brain bank PCA pattern explain করছে কি না
- MSA 30-এর সমস্যা শুধু Detection-P-তে, নাকি পুরো methylation profile-এও দেখা যায়

## MSA 30-এর সহজ ব্যাখ্যা

MSA 30-এর oxBS measurement `SAMPLE 67` clean, কিন্তু BS measurement `SAMPLE 145`-এ Detection-P failure অনেক বেশি। ফলে `BS − oxBS` difference অনেক জায়গায় negative হয়েছে এবং missing করতে হয়েছে। এটি technical measurement problem হতে পারে; donor-এর biology ভুল—এটি প্রমাণিত নয়। তাই Step 2 evidence ছাড়া donor বাদ দেওয়া হবে না।

## Analysis method

5mC ও 5hmC আলাদাভাবে:

1. সবচেয়ে variable 50,000 CpG নির্বাচন
2. 5hmC-এর remaining missing values শুধু PCA/QC-এর জন্য per-CpG median দিয়ে পূরণ
3. প্রতি CpG donor-দের মধ্যে mean-center করা; variance scaling করা হবে না
4. donor-level PCA
5. mean sample-correlation
6. robust PCA-distance ও robust correlation flags
7. condition, sex, brain bank এবং age দিয়ে PCA plots
8. Detection-P এবং 5hmC missingness evidence একত্র করে review list

Top-variable CpGs পুরো cohort দেখে নির্বাচন করা হচ্ছে বলে এগুলো শুধু unsupervised QC-এর জন্য। এগুলো পরবর্তী supervised ML feature selection হিসেবে ব্যবহার করা যাবে না। ML feature selection training fold-এর ভেতরে করতে হবে।

## আপনার exact directory অনুযায়ী চালানোর নিয়ম

`GSE143157_Step2_Pipeline.zip` extract করুন:

```text
D:\Paper\MSA_Research
```

Extract হওয়ার পরে folder:

```text
D:\Paper\MSA_Research\GSE143157_Step2_Pipeline
```

Folder খুলে double-click করুন:

```text
RUN_STEP2_WINDOWS.cmd
```

প্রথমবার packages install করার জন্য internet connection লাগবে। Black window বন্ধ করবেন না এবং computer sleep-এ যেতে দেবেন না।

শেষে দেখাবে:

```text
STEP 2 FINISHED SUCCESSFULLY
```

## Input

Runner নিজে এই completed Step 1 folder ব্যবহার করবে:

```text
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step1_results_full78
```

## Output

সব result যাবে:

```text
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step2_PCA_Outliers_full78
```

সবচেয়ে আগে browser-এ খুলবেন:

```text
STEP2_REPORT.html
```

গুরুত্বপূর্ণ output:

- `step2_summary.json` — authoritative summary ও review candidates
- `outlier_review.csv` — donor-wise সব evidence ও flags
- `pca_scores_5mC.csv`, `pca_scores_5hmC.csv` — PCA coordinates
- `pca_variance_5mC.csv`, `pca_variance_5hmC.csv` — explained variance
- `sample_correlations_*.csv` — 78×78 correlation matrices
- `top_variable_cpgs_*.csv` — PCA/QC-only CpG lists
- `pca_*.png` — main PCA figures
- `pca_covariates_*.png` — age/sex/brain-bank checks
- `correlation_heatmap_*.png` — sample correlations
- `outlier_diagnostics_*.png` — robust flag evidence
- `qc_5hmC_donor_missingness.png` — MSA 30-এর missingness context
- `SUCCESS.txt` — completed-run marker

## Run শেষ হলে আমাকে যা দেবেন

ছোট files:

```text
SUCCESS.txt
step2_summary.json
outlier_review.csv
pca_scores_5mC.csv
pca_scores_5hmC.csv
pca_variance_5mC.csv
pca_variance_5hmC.csv
```

আর এই PNG figures:

```text
pca_5mC.png
pca_5hmC.png
correlation_heatmap_5mC.png
correlation_heatmap_5hmC.png
outlier_diagnostics_5mC.png
outlier_diagnostics_5hmC.png
qc_5hmC_donor_missingness.png
pca_covariates_5mC.png
pca_covariates_5hmC.png
```

`sample_correlations_*.csv` এবং `top_variable_cpgs_*.csv` upload করার দরকার নেই, যদি না পরে চাই। Error হলে black window-এর সম্পূর্ণ screenshot পাঠাবেন।
