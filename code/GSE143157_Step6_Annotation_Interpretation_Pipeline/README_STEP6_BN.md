# GSE143157 Step 6 — সহজ নির্দেশনা

এই ধাপে Step 5-এর stable CpG-গুলোকে Illumina EPIC manifest ব্যবহার করে gene, chromosome, genomic region এবং CpG-island context-এর সঙ্গে যুক্ত করা হবে।

## এখন কী করবেন

1. ZIP file-টি `D:\Paper\MSA_Research`-এ Extract করুন।
2. Extract করা `GSE143157_Step6_Annotation_Interpretation_Pipeline` folder খুলুন।
3. `RUN_STEP6_WINDOWS.cmd` file-এ double-click করুন।
4. কালো window বন্ধ করবেন না। সাধারণত 5–20 মিনিট লাগতে পারে।
5. কাজ শেষ হলে HTML report নিজে থেকে খুলবে।

Final output folder:

```text
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step6_Annotation_Interpretation_Primary77
```

## Stable candidate কীভাবে নির্ধারিত হয়েছে?

Primary model:

```text
combined 5mC + 5hmC + covariates
```

- Tier A: কমপক্ষে 80% outer fit-এ selected এবং coefficient direction কমপক্ষে 90% consistent।
- Tier B: কমপক্ষে 60% outer fit-এ selected এবং coefficient direction কমপক্ষে 80% consistent।

এগুলো stability rule—p-value বা genome-wide significance নয়।

## Output file-এর অর্থ

- `primary_stable_cpgs_primary_clean77.csv`: আমাদের primary CpG candidates।
- `annotated_feature_stability_primary_clean77.csv`: Step 5-এর সব selected-feature evidence ও manifest annotation।
- `primary_stable_cpg_gene_map_primary_clean77.csv`: এক CpG একাধিক gene-এ map করলে আলাদা row।
- `candidate_gene_evidence_primary_clean77.csv`: gene-wise descriptive evidence summary।
- `published_signal_overlap_primary_clean77.csv`: original GSE143157 paper-এর reported CpG-এর exact overlap।
- `annotation_qc_summary_primary_clean77.csv`: annotation ও QC count।
- চারটি PNG: candidate stability, direction, genomic context এবং gene evidence।
- `STEP6_REPORT.html`: সহজ report।
- `gometh_selected_*.txt.gz` এবং `gometh_background_*.txt.gz`: bias-aware pathway analysis-এর জন্য প্রস্তুত input।
- `SUCCESS.txt`: Step 6 ঠিকভাবে শেষ হয়েছে কি না।

## Pathway analysis এখনই automatic কেন নয়?

EPIC array-তে প্রতিটি gene-এর probe সংখ্যা সমান নয়। সাধারণ gene enrichment করলে বেশি probe থাকা gene ভুলভাবে বেশি গুরুত্বপূর্ণ দেখাতে পারে। তাই main Python step কোনো naïve pathway result তৈরি করে না। Bias-aware `missMethyl/GOmeth`-এর input তৈরি করে।

`RUN_OPTIONAL_MISSMETHYL_ENRICHMENT.cmd` এখন চালানোর দরকার নেই। Main Step 6 result আগে review করতে দিন। R আগে থেকেই install থাকলে পরে এটি ব্যবহার করা যাবে।

## শেষ হলে যে files দেবেন

```text
SUCCESS.txt
step6_summary.json
annotation_qc_summary_primary_clean77.csv
primary_stable_cpgs_primary_clean77.csv
candidate_gene_evidence_primary_clean77.csv
published_signal_overlap_primary_clean77.csv
stable_cpg_annotation_primary.png
effect_stability_primary.png
genomic_context_primary.png
gene_evidence_primary.png
```

## গুরুত্বপূর্ণ বৈজ্ঞানিক সতর্কতা

- Stable CpG মানেই causal gene নয়।
- Stable CpG মানেই clinically validated biomarker নয়।
- Gene mapping regulatory relationship প্রমাণ করে না।
- Original EWAS এবং আমাদের predictive stability analysis আলাদা প্রশ্ন করে; exact overlap না থাকলে সেটি software error নয়।
- Coordinates Illumina EPIC v1.0 manifest-এর Genome Build 37/hg19 অনুযায়ী।
