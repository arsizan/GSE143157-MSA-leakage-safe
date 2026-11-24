# Step 8 label-only repair

এই patch শুধু repeated GO/KEGG নাম ঠিক করবে। `gometh` আবার চালাবে না এবং term ID, P-value, FDR বা gene result বদলাবে না।

## কী করবে

1. ZIP-টি extract করো।
2. Extracted folder-এর ভেতর `APPLY_STEP8_LABEL_REPAIR.cmd`-এ double-click করো।
3. শেষে সবুজ `STEP 8 LABEL REPAIR FINISHED SUCCESSFULLY` লেখা দেখবে এবং corrected HTML report খুলবে।
4. তারপর Step 9 যদি আগে শেষ হয়ে থাকে, Step 9 আবার চালাবে—কারণ তার Table 4 পুরোনো repeated labels পড়ে থাকতে পারে।

Patch নিজে থেকে বর্তমান original files-এর SHA-256 hash মিলিয়ে দেখবে। সব না মিললে কিছু overwrite করবে না। Apply হওয়ার আগে original files এই folder-এ backup হবে:

```text
D:\Paper\MSA_Research\02_Epigenetics\GSE143157\02_Preprocessing\Step8_BiasAware_Pathway_Enrichment\_backup_before_label_repair_YYYYMMDD_HHMMSS
```

Repair-এর পর `LABEL_REPAIR_APPLIED.txt`, `LABEL_REPAIR_LOCAL_APPLY.txt` এবং `label_repair_audit.csv` তৈরি হবে।
