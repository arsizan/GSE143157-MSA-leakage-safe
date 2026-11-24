args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript run_bias_aware_enrichment.R <Step6 output directory>")
}
out_dir <- normalizePath(args[[1]], mustWork = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
needed <- c(
  "missMethyl",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "org.Hs.eg.db"
)
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages(library(missMethyl))

for (modality in c("5mC", "5hmC")) {
  selected_path <- file.path(out_dir, paste0("gometh_selected_", modality, ".txt.gz"))
  background_path <- file.path(out_dir, paste0("gometh_background_", modality, ".txt.gz"))
  selected <- readLines(gzfile(selected_path))
  background <- readLines(gzfile(background_path))
  if (length(selected) < 5) {
    warning(paste("Too few selected CpGs for", modality, "- skipped"))
    next
  }
  go <- gometh(
    sig.cpg = selected,
    all.cpg = background,
    collection = "GO",
    array.type = "EPIC",
    prior.prob = TRUE,
    fract.counts = TRUE,
    plot.bias = FALSE,
    sig.genes = TRUE
  )
  kegg <- gometh(
    sig.cpg = selected,
    all.cpg = background,
    collection = "KEGG",
    array.type = "EPIC",
    prior.prob = TRUE,
    fract.counts = TRUE,
    plot.bias = FALSE,
    sig.genes = TRUE
  )
  write.csv(go, file.path(out_dir, paste0("gometh_GO_", modality, ".csv")), row.names = TRUE)
  write.csv(kegg, file.path(out_dir, paste0("gometh_KEGG_", modality, ".csv")), row.names = TRUE)
}

writeLines(
  c(
    "SUCCESS",
    paste("Completed UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    "Method: missMethyl::gometh with EPIC probe-number and multi-gene bias correction",
    "Interpretation: exploratory; input CpGs were ML-stable candidates, not EWAS-significant probes"
  ),
  file.path(out_dir, "MISSMETHYL_SUCCESS.txt")
)
