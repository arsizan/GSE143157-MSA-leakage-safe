# GSE143157 Step 7 — Sensitivity Validation

এই ধাপে Primary Clean-77 result দুইটি আগে থেকে নির্ধারিত alternative cohort-এ পরীক্ষা করা হবে:

- **Full-78:** সব 78 donor; MSA 30-সহ।
- **Strict Clean-76:** MSA 30 এবং Control 11—দুজনকেই বাদ দেওয়া।

এগুলো independent external validation cohort নয়। একই dataset-এর প্রায় একই donor ব্যবহার করে result outlier decision-এর প্রতি sensitive কি না পরীক্ষা করা হচ্ছে।

## কোথায় extract করবে

ZIP file-টি এই folder-এ extract করো:

```text
D:\Paper\MSA_Research
```

তারপর folder হবে:

```text
D:\Paper\MSA_Research\GSE143157_Step7_Molecular_Sensitivity_Pipeline
```

## কীভাবে চালাবে

এই file-এ double-click করো:

```text
RUN_STEP7_WINDOWS.cmd
```

প্রথমবার প্রয়োজনীয় Python packages install হবে। Internet connection লাগতে পারে।

## কত সময় লাগতে পারে

সাধারণত প্রায় 30–120 মিনিট লাগতে পারে। ধীর computer-এ আরও বেশি সময় লাগতে পারে।

- Laptop charger লাগিয়ে রাখবে।
- Sleep mode বন্ধ রাখবে।
- কালো window বন্ধ করবে না।
- অন্য ভারী program বন্ধ রাখলে ভালো।

প্রতিটি outer fold-এর পরে checkpoint save হয়। মাঝখানে computer বন্ধ হয়ে গেলে একই CMD file আবার চালালে completed checkpoint থেকে resume করবে।

## কোন model চলবে

Step 6-এর পরে model view lock করা হয়েছে। Step 7 নতুন করে best view খুঁজবে না। শুধু:

```text
5mC + covariates
5hmC + covariates
combined 5mC/5hmC + covariates
```

প্রতিটি cohort-এ একই Step 3 nested-CV folds এবং Step 5 tuning grid ব্যবহার করা হবে। Outer-test donor কোনো feature filtering, imputation, scaling, CpG selection বা tuning-এ ব্যবহার হবে না।

## Output কোথায় হবে

```text
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step7_Molecular_Sensitivity
```

সফল হলে কালো window-তে দেখা যাবে:

```text
STEP 7 FINISHED SUCCESSFULLY
```

এবং report নিজে থেকে browser-এ খুলবে।

## প্রধান output

- `SUCCESS.txt` — পুরো Step 7 সফল হয়েছে কি না।
- `step7_summary.json` — exact counts ও প্রধান ফল।
- `aggregate_metrics_across_cohorts.csv` — তিন cohort-এর performance।
- `primary_candidate_sensitivity.csv` — Step 6-এর 20 CpG দুই sensitivity cohort-এ কতবার selected হয়েছে এবং direction মিলেছে কি না।
- `stable_set_overlap.csv` — stable feature-set overlap।
- `performance_across_cohorts.png` — ROC-AUC comparison।
- `primary_candidate_selection_heatmap.png` — candidate stability heatmap।
- `candidate_effect_concordance.png` — effect direction comparison।
- `stable_set_overlap.png` — stable-set Jaccard overlap।
- `STEP7_REPORT.html` — সহজে দেখার report।

দুইটি sensitivity subfolder-এ predictions, hyperparameter choices, coefficients, selected features এবং checkpoints থাকবে। এগুলো audit/reproducibility-এর জন্য রাখা হবে।

## Sensitivity-supported candidate-এর নিয়ম

একটি Primary CpG-কে `supported_in_both` বলা হবে যদি:

1. Full-78-এ selection frequency কমপক্ষে 0.60 হয়;
2. Clean-76-এ selection frequency কমপক্ষে 0.60 হয়; এবং
3. দুই cohort-এই mean MSA–Control effect direction Primary Clean-77-এর সঙ্গে মিলে।

এটি descriptive stability rule—p-value বা independent replication নয়।

## শেষ হলে কী পাঠাবে

এই files পাঠাবে:

```text
SUCCESS.txt
step7_summary.json
aggregate_metrics_across_cohorts.csv
primary_candidate_sensitivity.csv
stable_set_overlap.csv
performance_across_cohorts.png
primary_candidate_selection_heatmap.png
candidate_effect_concordance.png
stable_set_overlap.png
```

`SUCCESS.txt` তৈরি হওয়ার আগে কোনো সংখ্যা paper-এ ব্যবহার করবে না।
