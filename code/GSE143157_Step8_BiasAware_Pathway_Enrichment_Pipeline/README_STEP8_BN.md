# GSE143157 Step 8 — Bias-aware Pathway Enrichment

এই ধাপে Step 6-এর stable, manifest-audited CpG candidate ব্যবহার করে `missMethyl::gometh` pathway analysis হবে। সাধারণ gene-list enrichment করা হবে না, কারণ methylation array-তে একেক gene-এর probe সংখ্যা আলাদা এবং একটি CpG একাধিক gene-এ map করতে পারে।

## আগে একবার R install করো

যদি computer-এ R আগে থেকে না থাকে:

1. এই official page খোলো: <https://cran.r-project.org/bin/windows/base/>
2. **Download R for Windows** installer download করো।
3. Installer-এ default options ও default folder রেখে install করো।
4. Computer restart দরকার নেই; এই CMD আবার চালাও।

RStudio দরকার নেই। শুধু R লাগবে।

## ZIP কোথায় extract করবে

ZIP file-টি এখানে extract করো:

```text
D:\Paper\MSA_Research
```

তাহলে folder হবে:

```text
D:\Paper\MSA_Research\GSE143157_Step8_BiasAware_Pathway_Enrichment_Pipeline
```

## কীভাবে চালাবে

এই file-এ double-click করো:

```text
RUN_STEP8_WINDOWS.cmd
```

প্রথম run-এ R/Bioconductor packages download হবে। প্রায় 20–60 মিনিট লাগতে পারে এবং internet connection লাগবে। Package install হয়ে গেলে enrichment নিজে থেকে চলবে।

Script version 1.0.4 থেকে packageগুলো `C:\Program Files\R`-এ নয়, ব্যবহারকারীর writable folder-এ install হয়:

```text
C:\Users\<your-name>\AppData\Local\GSE143157_R_library
```

তাই Administrator হিসেবে CMD চালানোর প্রয়োজন নেই। আগের failed run-এর partial Step 8 folder মুছতেও হবে না।
Source annotation packages একটির পর একটি install হবে (`Ncpus=1`), তাই আলাদা Rtools বা `make` command দরকার নেই। আগের failed installation-এর `00LOCK` folder থাকলে script শুধু project-specific R library থেকে সেটি নিজে পরিষ্কার করবে।
`KEGGREST` dependency-টি আগে আলাদাভাবে install ও load-test হবে; তারপর EPIC annotation এবং `org.Hs.eg.db` install হবে। আগের run-এ অর্ধেক install হওয়া target package প্রয়োজন হলে স্বয়ংক্রিয়ভাবে force-reinstall হবে।
GO ও KEGG-এর আলাদা raw column names একটি fixed common result schema-তে map হয়। Error হলেও log sink সঠিকভাবে বন্ধ হবে, ফলে R process আর অনন্তকাল চলবে না।

## Method lock

```text
Function: missMethyl::gometh
Array: Illumina MethylationEPIC v1
Annotation: ilm10b4.hg19
Probe-number correction: ON
Multi-gene fractional counting: ON
5mC background: 731,661 tested CpGs
5hmC background: 405,408 retained/tested CpGs
Selected CpGs: 7 per modality
Significance: Benjamini-Hochberg FDR < 0.05
```

Selected CpG মাত্র 7টি করে হওয়ায় power সীমিত। `FDR < 0.05` pathway না পাওয়া সম্পূর্ণ বৈধ result। শুধু nominal P-value দেখে pathway-কে significant বলা যাবে না।

## Output কোথায় হবে

```text
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step8_BiasAware_Pathway_Enrichment
```

## শেষ হলে কী পাঠাবে

```text
SUCCESS.txt
step8_summary.json
enrichment_summary.csv
top_terms_all.csv
significant_terms_FDR05.csv
top_enrichment_5mC.png
top_enrichment_5hmC.png
STEP8_REPORT.html
```

`SUCCESS.txt` এবং `STEP 8 INTEGRITY CHECK PASSED` হওয়ার আগে pathway result paper-এ ব্যবহার করবে না।
