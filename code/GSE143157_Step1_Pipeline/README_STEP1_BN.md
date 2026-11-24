# GSE143157 — Step 1: 5mC/5hmC matrix construction

এই package Windows-এ 2.3 GB `GSE143157_BMIQ_normalized_beta.txt` file-টি অল্প অল্প CpG row করে পড়ে। পুরো file RAM-এ একসাথে ওঠে না। এটি 156 array measurement-কে metadata দিয়ে 78 donor-এ pair করে এবং Step 2 PCA/outlier analysis-এর জন্য compressed HDF5 matrix ও QC files তৈরি করে।

## আপনার বর্তমান directory-এর জন্য সবচেয়ে সহজ পদ্ধতি

আপনার screenshot ও directory সরাসরি যাচাই করে exact paths পাওয়া গেছে:

```text
Matrix:
D:\Paper\MSA_Research\GSE143157_BMIQ_normalized_beta.txt\GSE143157_BMIQ_normalized_beta.txt

Metadata:
D:\Paper\MSA_Research\GSE143157_BMIQ_normalized_beta.txt\GSE143157_family.soft.gz

Step 1 result যাবে:
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step1_results_full78
```

আপনার জন্য recommended procedure:

1. `GSE143157_Step1_Pipeline.zip` extract করে পুরো folder-টি এখানে রাখুন:

   ```text
   D:\Paper\MSA_Research\GSE143157_Step1_Pipeline
   ```

2. ওই folder খুলে `RUN_STEP1_WINDOWS.cmd`-এ double-click করুন।
3. Windows warning দেখালে **More info → Run anyway** শুধু তখনই চাপবেন, যখন file-টি এই package থেকেই এসেছে।
4. কালো window বন্ধ করবেন না এবং computer sleep-এ যেতে দেবেন না। শেষে `STEP 1 FINISHED SUCCESSFULLY` লেখা দেখাবে।

এই one-click runner environment তৈরি, packages install, main analysis এবং output integrity check—সব নিজে করবে। `GSE143157_RAW.tar` বা `GSE143157_RAW`-এর IDAT files Step 1-এর জন্য ব্যবহার/আরও extract করতে হবে না। `.txt.gz`-ও আবার extract করতে হবে না, কারণ 2.3 GB text file ইতিমধ্যে পাওয়া গেছে।

## এই Step-এ বৈজ্ঞানিকভাবে কী হচ্ছে

প্রতি donor-এর দুটি measurement আছে:

- `oxBS beta` = 5mC estimate
- `BS beta` = 5mC + 5hmC signal
- `5hmC = BS beta - oxBS beta`

তারপর:

1. Negative 5hmC difference biological level হিসেবে গ্রহণ করা হয় না; `NaN`/missing করা হয়।
2. কোনো CpG-তে 78 donor-এর 20%-এর বেশি missing থাকলে সেই CpG 5hmC matrix থেকে বাদ যায়। 78 donor হলে সর্বোচ্চ 15টি missing রাখা যায়; 16 বা বেশি হলে বাদ।
3. Step 1-এ imputation করা হয় না। Step 2-এর unsupervised PCA-এর জন্য আলাদা documented imputation করা হবে; পরে ML-এ imputation শুধু training fold-এ fit করতে হবে, যাতে leakage না হয়।
4. Default run-এ Detection-P values QC report-এ যায়, কিন্তু beta values বদলায় না। কারণ inputটি authors-এর BMIQ-normalized, probe-QC-passed matrix। Detection-P দিয়ে আবার masking করতে চাইলে explicit `--mask-detection-failures` option আছে; এই alternate run-টি primary run-এর সঙ্গে mix করা যাবে না।

এই construction GSE143157-এর primary study-র Methods অনুসরণ করে: BS–oxBS delta, negative values missing, এবং `NA fraction > 0.2` probe removal। Primary source: Rydbirk et al., 2020, DOI: https://doi.org/10.1186/s40478-020-00908-7

## আগে folder প্রস্তুত করুন

ধরা হলো আপনার research folder:

```text
D:\Paper\MSA_Research
```

এই package folder-টি সেখানে copy/extract করে এমন রাখুন:

```text
D:\Paper\MSA_Research\GSE143157_Step1_Pipeline\
```

Raw matrix ও metadata যেখানে আছে সেখানেই থাকতে পারে; 2.3 GB file copy করার প্রয়োজন নেই। নিচের command-এ তাদের full path দেবেন।

Script আপনার `GSE143157_family.soft.gz` সরাসরি parse করে 156-row donor mapping তৈরি করতে পারে। চাইলে আগে থেকে তৈরি CSV/TSV mapping-ও ব্যবহার করা যায়। CSV/TSV-এর সাধারণ column names script নিজে চিনতে পারে, যেমন:

```text
sample_number / sample_column
biological_subject / biological_subject_id / donor_id
condition
conversion
```

প্রতি donor-এর ঠিক একটি `BS` ও একটি `oxBS` row থাকতে হবে। যদি column name আলাদা হয়, নিচে override command দেখানো হয়েছে।

## একদম exact Windows commands

### 1. Command Prompt খুলে package folder-এ যান

```bat
cd /d "D:\Paper\MSA_Research\GSE143157_Step1_Pipeline"
```

### 2. আলাদা Python environment তৈরি করুন

```bat
python -m venv .venv
```

যদি `python` command না চলে কিন্তু `py` চলে:

```bat
py -3 -m venv .venv
```

### 3. প্রয়োজনীয় packages install করুন

```bat
.venv\Scripts\python.exe -m pip install --upgrade pip
.venv\Scripts\python.exe -m pip install -r requirements_step1.txt
```

### 4. Input paths আছে কি না যাচাই করুন

আপনার ক্ষেত্রে metadata source হলো `GSE143157_family.soft.gz`:

```bat
dir "D:\Paper\MSA_Research\GSE143157_BMIQ_normalized_beta.txt\GSE143157_BMIQ_normalized_beta.txt"
dir "D:\Paper\MSA_Research\GSE143157_BMIQ_normalized_beta.txt\GSE143157_family.soft.gz"
```

দুই command-এই file দেখালে পরের ধাপে যান।

### 5. Primary Step 1 run করুন

```bat
.venv\Scripts\python.exe build_gse143157_step1.py --input "D:\Paper\MSA_Research\GSE143157_BMIQ_normalized_beta.txt\GSE143157_BMIQ_normalized_beta.txt" --metadata "D:\Paper\MSA_Research\GSE143157_BMIQ_normalized_beta.txt\GSE143157_family.soft.gz" --output "D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step1_results_full78"
```

Compressed `.txt.gz` file সরাসরি ব্যবহার করলে:

```bat
.venv\Scripts\python.exe build_gse143157_step1.py --input "D:\Paper\MSA_Research\GSE143157_BMIQ_normalized_beta.txt.gz" --metadata "D:\Paper\MSA_Research\GSE143157_BMIQ_normalized_beta.txt\GSE143157_family.soft.gz" --output "D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step1_results_full78"
```

একই output folder non-empty থাকলে script ইচ্ছাকৃতভাবে থামবে। পুরোনো result overwrite না করে নতুন নাম দিন, যেমন:

```text
GSE143157_STEP1_RESULTS_v2
```

### Metadata column auto-detection ব্যর্থ হলে

আপনার mapping CSV-তে column names যদি হয় `Sample`, `Biological subject`, `Conversion`, `Condition`, exact override command:

```bat
.venv\Scripts\python.exe build_gse143157_step1.py --input "D:\Paper\MSA_Research\GSE143157_BMIQ_normalized_beta.txt\GSE143157_BMIQ_normalized_beta.txt" --metadata "D:\Paper\MSA_Research\GSE143157_donor_mapping.csv" --output "D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step1_results_full78" --sample-field "Sample" --donor-field "Biological subject" --conversion-field "Conversion" --condition-field "Condition"
```

### Optional sensitivity run: Detection-P masking

এটি primary run নয়; আলাদা output folder ব্যবহার করুন:

```bat
.venv\Scripts\python.exe build_gse143157_step1.py --input "D:\Paper\MSA_Research\GSE143157_BMIQ_normalized_beta.txt\GSE143157_BMIQ_normalized_beta.txt" --metadata "D:\Paper\MSA_Research\GSE143157_BMIQ_normalized_beta.txt\GSE143157_family.soft.gz" --output "D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step1_results_detp_masked" --mask-detection-failures --detection-p-threshold 0.01
```

## Run চলার সময় কী দেখবেন

প্রথমে script layout ও metadata validate করবে। তারপর প্রায় প্রতি 50,000 CpG-এ progress দেখাবে:

```text
Validated layout: 156 beta + 156 Detection-P columns
Validated metadata: 78 paired donors; conditions={'Control': 37, 'MSA': 41}
Processed 50,000 CpGs | retained 5hmC ... | elapsed ... min
...
SUCCESS
```

`SUCCESS` না দেখা পর্যন্ত result complete ধরে নেবেন না। Laptop sleep বন্ধ রাখা ভালো। Script-এর default chunk 5,000 rows; RAM কম হলে `--chunk-rows 2000` যোগ করতে পারেন। এতে RAM কম লাগবে, সময় একটু বাড়বে।

## Output file-এর সহজ অর্থ

Result folder-এ থাকবে:

### `SUCCESS.txt`

Run সম্পূর্ণ হয়েছে—এই marker। এটি না থাকলে Step 2 শুরু করবেন না।

### `qc_summary.json`

সবচেয়ে গুরুত্বপূর্ণ machine-readable summary। এতে input dimensions, 78 donor pairing, MSA/control count, beta/Detection-P QC, কত negative 5hmC missing হয়েছে, এবং কত CpG >20% missing-এর জন্য বাদ গেছে থাকবে। Paper-এ কোনো সংখ্যা লেখার আগে এই file থেকে নিতে হবে।

### `GSE143157_step1_matrices.h5`

Compressed HDF5 data file। ভিতরে:

```text
/5mC/values     rows = CpGs, columns = 78 donors
/5mC/cpg_ids
/5hmC/values    rows = retained CpGs, columns = 78 donors
/5hmC/cpg_ids
/donor_ids
/conditions
```

এটি Excel-এ খোলার file নয়। Step 2 Python script chunk করে পড়বে। `float32` storage যথেষ্ট precision রাখে এবং file size/RAM কমায়।

### `metadata_arrays_aligned.csv`

156 array-এর matrix order, donor ID, condition, এবং BS/oxBS conversion। Pairing audit করার file।

### `metadata_donors_aligned.csv`

78 donor-এর matrix column order। Step 2 plots ও পরবর্তী donor-aware ML split এই order ব্যবহার করবে। একই donor-এর BS/oxBS আর আলাদা sample নয়।

### `qc_per_array.csv`

প্রতি 156 array-এর mean/range beta এবং Detection-P failure counts/fractions। আগে flagged `SAMPLE 145` এখানে independently দেখা যাবে; script কোনো outlier নিজে delete করে না।

### `qc_per_donor.csv`

প্রতি 78 donor-এর 5mC/5hmC missingness, mean values এবং negative 5hmC difference fraction। Step 2 sample-level QC-তে কাজে লাগবে।

### `5hmC_probe_qc.tsv.gz`

প্রতি input CpG-এর missing-pair count, negative-difference count, final missing fraction এবং retained/removed flag। Filtering সম্পূর্ণ auditable করে। Compressed text; Excel-এ পুরোটা খোলার দরকার নেই।

## Output integrity check

Primary run সফল হওয়ার পরে:

```bat
.venv\Scripts\python.exe check_step1_outputs.py "D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step1_results_full78"
```

Expected শেষ message:

```text
PASS: Step 1 outputs are internally consistent.
```

এই check scientific validity প্রমাণ করে না; files/shapes/metadata একে অপরের সঙ্গে internally consistent কি না যাচাই করে।

## কোন result আগে থেকে ধরে নেওয়া হয়নি

- Full 78-donor run-এ retained 5hmC CpG count script চালানোর আগে জানা নয়।
- Published paper-এ 15 outlier সরিয়ে 63 sample-এর পরে 405,408 5hmC probe ছিল। আমাদের Step 1 full 78 donor দিয়ে PCA/outlier discovery প্রস্তুত করছে, তাই retained count একই হতেই হবে—এমন expectation রাখা যাবে না।
- Step 1 কোনো donor/outlier বাদ দেয় না, কোনো missing value impute করে না, এবং কোনো ML classification করে না।

## Run শেষ হলে পরের কাজ

`SUCCESS.txt`, `qc_summary.json`, `qc_per_array.csv`, `qc_per_donor.csv`, এবং `metadata_donors_aligned.csv` আগে review করতে হবে। তারপর Step 2-তে donor-level PCA, robust sample-distance, correlation heatmap, এবং outlier sensitivity plan করা যাবে।
