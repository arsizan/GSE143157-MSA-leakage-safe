# GSE143157 Step 9 — Final Audit and Paper Package

Step 8 সফল হওয়ার পরে এই ধাপ চালাবে। এটি নতুন model train করবে না। Steps 1–8-এর verified result একসাথে মিলিয়ে paper-এর জন্য final source-of-truth package তৈরি করবে।

## ZIP কোথায় extract করবে

```text
D:\Paper\MSA_Research
```

তাহলে folder হবে:

```text
D:\Paper\MSA_Research\GSE143157_Step9_Final_Audit_Paper_Package
```

## কীভাবে চালাবে

এই file-এ double-click করো:

```text
RUN_STEP9_WINDOWS.cmd
```

সাধারণত 2–10 মিনিট লাগবে। এটি Steps 1–8-এর `SUCCESS.txt` না পেলে নিজে থেকেই থেমে যাবে।

## Output কোথায় হবে

```text
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step9_Final_Audit_Paper_Package
```

## প্রধান output

- `TABLE1_cohorts.csv` — primary ও sensitivity cohort-এর exact counts।
- `TABLE2_model_performance.csv` — চারটি locked model, তিনটি cohort এবং confidence interval।
- `TABLE3_robust_cpg_candidates.csv` — annotation, effect, manifest QC ও sensitivity supportসহ 20 CpG।
- `TABLE4_bias_aware_pathways.csv` — GOmeth top terms এবং explicit FDR status।
- `PAPER_NUMBERS.json` — paper-এ ব্যবহারযোগ্য exact machine-readable numbers।
- `METHODS_RESULTS_NUMBERS.md` — data-backed English drafting aid।
- `PAPER_CLAIMS_GUARDRAILS.csv` — কোন দাবি করা যাবে এবং কোনটি যাবে না।
- `FINAL_METHODS_AUDIT.csv` — প্রতিটি locked rule pass করেছে কি না।
- `FINAL_INPUT_FILE_MANIFEST.csv` — upstream evidence file-এর SHA-256।
- `FINAL_FREEZE.json` — final output-এর cryptographic freeze record।
- `figures` — paper-ready নাম দেওয়া main ও supplementary figures।
- `STEP9_REPORT.html` — সহজে দেখার report।

## খুব গুরুত্বপূর্ণ

`FINAL_FREEZE.json` তৈরি হওয়ার পরে output files সরাসরি edit করবে না। Paper-এর figure crop বা table formatting দরকার হলে মূল freeze অক্ষত রেখে আলাদা manuscript folder-এ copy করবে।

## শেষ হলে কী পাঠাবে

```text
SUCCESS.txt
PAPER_NUMBERS.json
TABLE1_cohorts.csv
TABLE2_model_performance.csv
TABLE3_robust_cpg_candidates.csv
TABLE4_bias_aware_pathways.csv
FINAL_METHODS_AUDIT.csv
PAPER_CLAIMS_GUARDRAILS.csv
FINAL_FREEZE.json
STEP9_REPORT.html
```

সাথে `figures` folder-এর screenshot বা প্রধান PNG files পাঠাবে। `STEP 9 INTEGRITY CHECK PASSED` না দেখা পর্যন্ত paper লেখা শুরু করবে না।

