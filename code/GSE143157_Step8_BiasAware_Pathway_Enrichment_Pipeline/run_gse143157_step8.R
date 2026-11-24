args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "Rscript run_gse143157_step8.R <Step6 output directory> <Step8 output directory>"
)
if (length(args) != 2) {
  stop(usage)
}

step6_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)

if (file.exists(file.path(output_dir, "SUCCESS.txt"))) {
  stop("A completed Step 8 result already exists. Nothing was overwritten.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

options(timeout = 3600)
options(repos = c(CRAN = "https://cloud.r-project.org"))

configure_user_library <- function() {
  local_app_data <- Sys.getenv("LOCALAPPDATA", unset = "")
  if (!nzchar(local_app_data)) {
    local_app_data <- path.expand("~")
  }
  r_minor <- strsplit(R.version$minor, ".", fixed = TRUE)[[1]][[1]]
  library_dir <- file.path(
    local_app_data,
    "GSE143157_R_library",
    paste0("R-", R.version$major, ".", r_minor)
  )
  dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(library_dir) || file.access(library_dir, 2) != 0) {
    stop(paste("Cannot create a writable user R library:", library_dir))
  }
  .libPaths(unique(c(library_dir, .libPaths())))
  Sys.setenv(R_LIBS_USER = library_dir)
  normalizePath(library_dir, winslash = "/", mustWork = TRUE)
}

user_library <- configure_user_library()
cat("Using writable user R package library:\n", user_library, "\n\n")

clear_stale_install_locks <- function() {
  entries <- list.files(user_library, full.names = TRUE, all.files = TRUE)
  locks <- entries[grepl("^00LOCK", basename(entries))]
  if (length(locks) > 0L) {
    cat("Removing stale package-install lock(s) left by the previous failed run:\n")
    cat(paste0("  ", locks, collapse = "\n"), "\n")
    unlink(locks, recursive = TRUE, force = TRUE)
  }
}

clear_stale_install_locks()

install_required_packages <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages(
      "BiocManager",
      repos = "https://cloud.r-project.org",
      lib = user_library
    )
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    install.packages(
      "jsonlite",
      repos = "https://cloud.r-project.org",
      lib = user_library
    )
  }

  # KEGGREST must load before AnnotationDbi, org.Hs.eg.db, bumphunter,
  # and the EPIC annotation package can complete their installation tests.
  # A prior interrupted parallel install may leave the higher-level packages
  # present but unloadable, so bootstrap this dependency first and verify it.
  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    cat("Bootstrapping missing dependency KEGGREST before annotation packages...\n")
    BiocManager::install(
      "KEGGREST",
      ask = FALSE,
      update = FALSE,
      force = TRUE,
      lib = user_library,
      Ncpus = 1L
    )
  }
  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    stop("KEGGREST could not be installed or loaded; dependent annotation packages were not attempted.")
  }

  bioc_needed <- c(
    "missMethyl",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "org.Hs.eg.db"
  )
  missing <- bioc_needed[
    !vapply(bioc_needed, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing) > 0) {
    BiocManager::install(
      missing,
      ask = FALSE,
      update = FALSE,
      force = TRUE,
      lib = user_library,
      Ncpus = 1L
    )
  }

  still_missing <- c("missMethyl", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "org.Hs.eg.db", "jsonlite")[
    !vapply(
      c("missMethyl", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "org.Hs.eg.db", "jsonlite"),
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]
  if (length(still_missing) > 0) {
    stop(paste("Required R packages are still missing:", paste(still_missing, collapse = ", ")))
  }
}

install_required_packages()
suppressPackageStartupMessages(library(missMethyl))

run_analysis <- function() {
log_path <- file.path(output_dir, "step8_console_log.txt")
log_con <- file(log_path, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)
on.exit({
  if (sink.number(type = "message") != 2L) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  if (isOpen(log_con)) close(log_con)
}, add = TRUE)

started <- Sys.time()
cat("GSE143157 Step 8 started:", format(started, tz = "UTC", usetz = TRUE), "\n")
cat("Step 6 input:", step6_dir, "\n")
cat("Step 8 output:", output_dir, "\n")

required_step6 <- c(
  "SUCCESS.txt",
  "step6_summary.json",
  "gometh_selected_5mC.txt.gz",
  "gometh_background_5mC.txt.gz",
  "gometh_selected_5hmC.txt.gz",
  "gometh_background_5hmC.txt.gz"
)
missing_step6 <- required_step6[!file.exists(file.path(step6_dir, required_step6))]
if (length(missing_step6) > 0) {
  stop(paste("Missing Step 6 files:", paste(missing_step6, collapse = ", ")))
}

read_cpgs <- function(path) {
  values <- trimws(readLines(gzfile(path), warn = FALSE))
  values <- unique(values[nzchar(values)])
  invalid <- values[!grepl("^cg[0-9]+$", values)]
  if (length(invalid) > 0) {
    stop(paste("Invalid CpG identifiers in", basename(path)))
  }
  values
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

standardize_result <- function(x, modality, collection) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x$term_id <- rownames(x)
  rownames(x) <- NULL

  if (collection == "GO") {
    if ("TERM" %in% colnames(x)) x$Term <- as.character(x$TERM)
    if ("ONTOLOGY" %in% colnames(x)) x$Ont <- as.character(x$ONTOLOGY)
  }
  if (collection == "KEGG" && "Description" %in% colnames(x)) {
    x$Pathway <- as.character(x$Description)
  }

  for (column in c("Term", "Ont", "Pathway", "SigGenesInSet")) {
    if (!column %in% colnames(x)) x[[column]] <- ""
    x[[column]][is.na(x[[column]])] <- ""
  }
  for (column in c("N", "DE", "P.DE", "FDR")) {
    if (!column %in% colnames(x)) x[[column]] <- NA_real_
    x[[column]] <- safe_numeric(x[[column]])
  }

  x$modality <- modality
  x$collection <- collection
  x$display_term <- ifelse(collection == "GO", x$Term, x$Pathway)
  blank <- !nzchar(trimws(x$display_term))
  x$display_term[blank] <- x$term_id[blank]
  x <- x[order(x$P.DE, x$FDR, na.last = TRUE), , drop = FALSE]
  x$rank_by_p <- seq_len(nrow(x))
  x$fdr_significant_0_05 <- !is.na(x$FDR) & x$FDR < 0.05
  x$nominal_p_lt_0_05 <- !is.na(x$P.DE) & x$P.DE < 0.05

  preferred <- c(
    "modality", "collection", "term_id", "display_term", "Term", "Ont",
    "Pathway", "N", "DE", "P.DE", "FDR", "SigGenesInSet", "rank_by_p",
    "fdr_significant_0_05", "nominal_p_lt_0_05"
  )
  x[, preferred, drop = FALSE]
}

plot_top_terms <- function(results, modality, selected_count, output_path) {
  top <- do.call(rbind, lapply(results, function(x) head(x, 10)))
  top <- top[order(top$P.DE, na.last = TRUE), , drop = FALSE]
  top <- head(top, 15)

  png(output_path, width = 2200, height = 1500, res = 190)
  on.exit(dev.off(), add = TRUE)

  if (nrow(top) == 0 || all(is.na(top$P.DE))) {
    plot.new()
    text(0.5, 0.5, paste("No", modality, "enrichment results were returned"), cex = 1.4)
    return(invisible(NULL))
  }

  labels <- paste0(top$collection, ": ", top$display_term)
  labels <- ifelse(nchar(labels) > 70, paste0(substr(labels, 1, 67), "..."), labels)
  scores <- -log10(pmax(top$P.DE, .Machine$double.xmin))
  colours <- ifelse(top$collection == "GO", "#2C7FB8", "#F28E2B")

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(5.5, 24, 4.5, 2), xpd = FALSE)
  bp <- barplot(
    rev(scores),
    names.arg = rev(labels),
    horiz = TRUE,
    las = 1,
    col = rev(colours),
    border = NA,
    xlab = expression(-log[10]("bias-aware GOmeth P-value")),
    main = paste0(modality, " exploratory bias-aware pathway ranking (", selected_count, " CpGs)")
  )
  abline(v = -log10(0.05), lty = 2, col = "#666666")
  legend(
    "bottomright",
    legend = c("GO", "KEGG", "nominal P=0.05"),
    fill = c("#2C7FB8", "#F28E2B", NA),
    border = c(NA, NA, NA),
    lty = c(NA, NA, 2),
    col = c(NA, NA, "#666666"),
    bty = "n"
  )
  invisible(bp)
}

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

modalities <- c("5mC", "5hmC")
all_results <- list()
summary_rows <- list()
input_rows <- list()

for (modality in modalities) {
  selected_path <- file.path(step6_dir, paste0("gometh_selected_", modality, ".txt.gz"))
  background_path <- file.path(step6_dir, paste0("gometh_background_", modality, ".txt.gz"))
  selected <- read_cpgs(selected_path)
  background <- read_cpgs(background_path)

  missing_from_background <- setdiff(selected, background)
  if (length(missing_from_background) > 0) {
    stop(paste(modality, "selected CpGs were absent from the tested background"))
  }
  if (length(selected) < 5) {
    stop(paste(modality, "has fewer than 5 selected CpGs; enrichment is not run"))
  }

  input_rows[[length(input_rows) + 1]] <- data.frame(
    modality = modality,
    selected_cpgs = length(selected),
    background_cpgs = length(background),
    selected_unique = !anyDuplicated(selected),
    background_unique = !anyDuplicated(background),
    selected_subset_of_background = length(missing_from_background) == 0,
    selected_file_md5 = unname(tools::md5sum(selected_path)),
    background_file_md5 = unname(tools::md5sum(background_path)),
    stringsAsFactors = FALSE
  )

  modality_results <- list()
  for (collection in c("GO", "KEGG")) {
    cat("Running", modality, collection, "GOmeth...\n")
    raw <- missMethyl::gometh(
      sig.cpg = selected,
      all.cpg = background,
      collection = collection,
      array.type = "EPIC",
      plot.bias = FALSE,
      prior.prob = TRUE,
      equiv.cpg = TRUE,
      fract.counts = TRUE,
      genomic.features = "ALL",
      sig.genes = TRUE
    )
    result <- standardize_result(raw, modality, collection)
    modality_results[[collection]] <- result
    all_results[[paste(modality, collection, sep = "_")]] <- result

    write.csv(
      result,
      file.path(output_dir, paste0("gometh_", collection, "_", modality, ".csv")),
      row.names = FALSE,
      na = ""
    )

    valid_p <- result$P.DE[is.finite(result$P.DE)]
    valid_fdr <- result$FDR[is.finite(result$FDR)]
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      modality = modality,
      collection = collection,
      selected_cpgs = length(selected),
      background_cpgs = length(background),
      tested_terms = nrow(result),
      nominal_p_lt_0_05 = sum(result$nominal_p_lt_0_05, na.rm = TRUE),
      fdr_lt_0_05 = sum(result$fdr_significant_0_05, na.rm = TRUE),
      minimum_p_value = if (length(valid_p) > 0) min(valid_p) else NA_real_,
      minimum_fdr = if (length(valid_fdr) > 0) min(valid_fdr) else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  plot_top_terms(
    modality_results,
    modality,
    length(selected),
    file.path(output_dir, paste0("top_enrichment_", modality, ".png"))
  )
}

input_qc <- do.call(rbind, input_rows)
enrichment_summary <- do.call(rbind, summary_rows)
all_terms <- do.call(rbind, all_results)
rownames(all_terms) <- NULL

top_terms <- do.call(
  rbind,
  lapply(split(all_terms, interaction(all_terms$modality, all_terms$collection)), function(x) {
    head(x[order(x$P.DE, x$FDR, na.last = TRUE), , drop = FALSE], 20)
  })
)
rownames(top_terms) <- NULL

significant_terms <- all_terms[
  !is.na(all_terms$FDR) & all_terms$FDR < 0.05,
  ,
  drop = FALSE
]

write.csv(input_qc, file.path(output_dir, "gometh_input_qc.csv"), row.names = FALSE, na = "")
write.csv(enrichment_summary, file.path(output_dir, "enrichment_summary.csv"), row.names = FALSE, na = "")
write.csv(top_terms, file.path(output_dir, "top_terms_all.csv"), row.names = FALSE, na = "")
write.csv(
  significant_terms,
  file.path(output_dir, "significant_terms_FDR05.csv"),
  row.names = FALSE,
  na = ""
)

package_names <- c(
  "R",
  "BiocManager",
  "missMethyl",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "org.Hs.eg.db",
  "jsonlite"
)
package_versions <- data.frame(
  package = package_names,
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("BiocManager")),
    as.character(packageVersion("missMethyl")),
    as.character(packageVersion("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")),
    as.character(packageVersion("org.Hs.eg.db")),
    as.character(packageVersion("jsonlite"))
  ),
  stringsAsFactors = FALSE
)
write.csv(package_versions, file.path(output_dir, "package_versions.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_SESSION_INFO.txt"))

completed <- Sys.time()
configuration <- list(
  pipeline_version = "1.0.4",
  step6_dir = step6_dir,
  output_dir = output_dir,
  user_r_library = user_library,
  array_type = "EPIC",
  annotation = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  genome_build = "hg19",
  collections = c("GO", "KEGG"),
  package_install_ncpus = 1,
  dependency_bootstrap = "KEGGREST_then_targets",
  prior_prob = TRUE,
  equivalent_cpg_counts = TRUE,
  fractional_multi_gene_counts = TRUE,
  genomic_features = "ALL",
  fdr_method = "Benjamini-Hochberg as returned by missMethyl::gometh",
  fdr_reporting_threshold = 0.05,
  interpretation = "exploratory_hypothesis_generating",
  selected_cpg_definition = "Step 6 stable candidates that passed the manifest audit",
  selected_cpgs_are_ewas_significant = FALSE
)
jsonlite::write_json(
  configuration,
  file.path(output_dir, "RUN_CONFIGURATION.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  digits = 16
)

summary_json <- list(
  pipeline = list(
    name = "GSE143157 Step 8 bias-aware exploratory pathway enrichment",
    version = "1.0.4",
    completed_utc = format(completed, tz = "UTC", usetz = TRUE),
    elapsed_seconds = as.numeric(difftime(completed, started, units = "secs"))
  ),
  method = list(
    function_name = "missMethyl::gometh",
    array_type = "EPIC",
    annotation = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    genome_build = "hg19",
    prior_prob = TRUE,
    equiv_cpg = TRUE,
    fract_counts = TRUE,
    collections = c("GO", "KEGG"),
    fdr_threshold = 0.05
  ),
  inputs = lapply(seq_len(nrow(input_qc)), function(i) as.list(input_qc[i, , drop = FALSE])),
  results = lapply(seq_len(nrow(enrichment_summary)), function(i) as.list(enrichment_summary[i, , drop = FALSE])),
  total_fdr_significant_terms = nrow(significant_terms),
  interpretation_guardrails = c(
    "The input CpGs were ML-stable candidates, not EWAS-significant CpGs.",
    "GOmeth corrects for differing probes per gene and multi-gene annotated probes.",
    "5mC and 5hmC used separate tested-probe backgrounds.",
    "Only FDR below 0.05 is called statistically significant.",
    "Nominally ranked pathways are exploratory and must not be presented as confirmed mechanisms.",
    "The selected sets contain only seven CpGs per modality, so power is limited."
  )
)
jsonlite::write_json(
  summary_json,
  file.path(output_dir, "step8_summary.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  digits = 16,
  na = "null"
)

summary_html_rows <- apply(enrichment_summary, 1, function(row) {
  paste0(
    "<tr><td>", html_escape(row[["modality"]]), "</td><td>",
    html_escape(row[["collection"]]), "</td><td>", row[["selected_cpgs"]],
    "</td><td>", row[["background_cpgs"]], "</td><td>", row[["tested_terms"]],
    "</td><td>", row[["nominal_p_lt_0_05"]], "</td><td>", row[["fdr_lt_0_05"]],
    "</td><td>", format(as.numeric(row[["minimum_fdr"]]), digits = 4), "</td></tr>"
  )
})

status_text <- if (nrow(significant_terms) > 0) {
  paste(nrow(significant_terms), "term(s) passed FDR < 0.05.")
} else {
  "No term passed FDR < 0.05. This is a valid, reportable result."
}

html <- paste0(
  "<!doctype html><html><head><meta charset='utf-8'><title>GSE143157 Step 8</title>",
  "<style>body{font-family:Segoe UI,Arial,sans-serif;max-width:1100px;margin:32px auto;line-height:1.5;color:#202124}",
  "h1,h2{color:#17365d}.ok{background:#e8f5e9;border-left:5px solid #2e7d32;padding:12px}",
  ".warn{background:#fff8e1;border-left:5px solid #f9a825;padding:12px}",
  "table{border-collapse:collapse;width:100%;margin:14px 0}th,td{border:1px solid #ccc;padding:7px;text-align:left}",
  "th{background:#eef3f8}img{max-width:100%;height:auto;border:1px solid #ddd}</style></head><body>",
  "<h1>GSE143157 Step 8: Bias-aware pathway enrichment</h1>",
  "<div class='ok'><b>Run completed.</b> Method: missMethyl::gometh with EPIC probe-number and multi-gene bias correction.</div>",
  "<div class='warn'><b>Interpretation:</b> ", html_escape(status_text),
  " The tested CpGs are stable ML candidates, not EWAS-significant loci; all pathway results are exploratory.</div>",
  "<h2>QC and result counts</h2><table><tr><th>Modality</th><th>Collection</th><th>Selected CpGs</th>",
  "<th>Background CpGs</th><th>Tested terms</th><th>Nominal P&lt;0.05</th><th>FDR&lt;0.05</th><th>Minimum FDR</th></tr>",
  paste(summary_html_rows, collapse = ""), "</table>",
  "<h2>Top-ranked 5mC terms</h2><img src='top_enrichment_5mC.png' alt='5mC top enrichment'>",
  "<h2>Top-ranked 5hmC terms</h2><img src='top_enrichment_5hmC.png' alt='5hmC top enrichment'>",
  "<h2>Files to use</h2><ul><li><code>enrichment_summary.csv</code>: exact counts and minimum P/FDR.</li>",
  "<li><code>top_terms_all.csv</code>: top 20 terms per modality and collection.</li>",
  "<li><code>significant_terms_FDR05.csv</code>: only FDR-significant terms; it may correctly contain zero rows.</li>",
  "<li><code>gometh_GO_*.csv</code> and <code>gometh_KEGG_*.csv</code>: complete result tables.</li></ul>",
  "</body></html>"
)
writeLines(html, file.path(output_dir, "STEP8_REPORT.html"), useBytes = TRUE)

writeLines(
  c(
    "SUCCESS",
    paste("Completed UTC:", format(completed, tz = "UTC", usetz = TRUE)),
    "Analysis mode: BIAS_AWARE_EXPLORATORY_ENRICHMENT",
    "Method: missMethyl::gometh; EPIC; prior.prob=TRUE; equiv.cpg=TRUE; fract.counts=TRUE",
    paste("FDR-significant terms:", nrow(significant_terms)),
    "Interpretation: exploratory; ML-stable candidates are not EWAS-significant probes"
  ),
  file.path(output_dir, "SUCCESS.txt")
)

cat("GSE143157 Step 8 completed successfully.\n")
cat("FDR-significant terms:", nrow(significant_terms), "\n")
}

run_analysis()
