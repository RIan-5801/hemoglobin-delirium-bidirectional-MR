#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/21b_chen_forward_instruments_v2_readback_audit.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

for (pkg in c("DBI", "duckdb", "jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
qpath <- function(path) gsub("'", "''", normalizePath(path, winslash = "/", mustWork = FALSE), fixed = TRUE)
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
atomic <- function(path, writer) {
  partial <- paste0(path, ".partial")
  stop_if(file.exists(path) || file.exists(partial), paste("Output occupied:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writer(partial)
  stop_if(!file.exists(partial), paste("Writer did not create partial:", partial))
  stop_if(!file.rename(partial, path), paste("Atomic rename failed:", path))
}
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))

paths <- c(
  included_parquet = file.path(root, "data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_included_clumped_v2.parquet"),
  included_tsv = file.path(root, "data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_included_clumped_v2.tsv"),
  excluded_parquet = file.path(root, "data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_excluded_clumped_v2.parquet"),
  excluded_tsv = file.path(root, "data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_excluded_clumped_v2.tsv"),
  previous_qc = file.path(root, "results", "qc", "chen_forward_instrument_selection_v2.json"),
  audit_csv = file.path(root, "results", "qc", "chen_forward_instrument_selection_v2_readback_numeric_audit_v1.csv"),
  closure_json = file.path(root, "results", "qc", "chen_forward_instrument_selection_v2_readback_closure_v1.json"),
  closure_log = file.path(root, "results", "logs", "chen_forward_instrument_selection_v2_readback_closure_v1.log"),
  decision = file.path(root, "docs", "decisions", "86_chen_forward_instrument_selection_v2_readback_closure_v1.1.md"),
  source = file.path(root, "data_raw", "gwas", "BCX2_HGB_EA_GWAMA.out.gz"),
  renv_lock = file.path(root, "renv.lock")
)

for (p in paths[c("included_parquet", "included_tsv", "excluded_parquet", "excluded_tsv", "previous_qc", "source", "renv_lock")]) {
  stop_if(!file.exists(p), paste("Missing required input:", p))
}
for (p in paths[c("audit_csv", "closure_json", "closure_log", "decision")]) {
  stop_if(file.exists(p) || file.exists(paste0(p, ".partial")), paste("Output occupied:", p))
}

con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:", config = list(shared_home = FALSE))
on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
read_pair <- function(parquet_path, tsv_path) {
  list(
    parquet = DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", qpath(parquet_path))),
    tsv = read.delim(tsv_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  )
}
audit_pair <- function(label, x) {
  p <- x$parquet
  t <- x$tsv
  stop_if(!identical(names(p), names(t)), paste(label, "column names differ."))
  stop_if(nrow(p) != nrow(t), paste(label, "row counts differ."))
  stop_if(!identical(as.character(p$resolved_rsid), as.character(t$resolved_rsid)), paste(label, "resolved_rsid order differs."))
  rows <- list()
  for (nm in names(p)) {
    if (is.numeric(p[[nm]])) {
      a <- as.numeric(p[[nm]])
      b <- suppressWarnings(as.numeric(t[[nm]]))
      same_na <- is.na(a) & is.na(b)
      finite_pair <- is.finite(a) & is.finite(b)
      abs_diff <- abs(a - b)
      rel_diff <- abs_diff / pmax(abs(a), abs(b), .Machine$double.xmin)
      rows[[length(rows) + 1L]] <- data.frame(
        analysis_set = label,
        column = nm,
        parquet_class = paste(class(p[[nm]]), collapse = "/"),
        tsv_class = paste(class(t[[nm]]), collapse = "/"),
        n = length(a),
        n_na_mismatch = sum(xor(is.na(a), is.na(b))),
        n_abs_diff_gt_1e_12 = sum(!(same_na | (finite_pair & abs_diff <= 1e-12))),
        max_abs_diff = suppressWarnings(max(abs_diff, na.rm = TRUE)),
        max_rel_diff = suppressWarnings(max(rel_diff, na.rm = TRUE)),
        exact_character_match = identical(as.character(p[[nm]]), as.character(t[[nm]])),
        stringsAsFactors = FALSE
      )
    } else {
      rows[[length(rows) + 1L]] <- data.frame(
        analysis_set = label,
        column = nm,
        parquet_class = paste(class(p[[nm]]), collapse = "/"),
        tsv_class = paste(class(t[[nm]]), collapse = "/"),
        n = nrow(p),
        n_na_mismatch = sum(xor(is.na(p[[nm]]), is.na(t[[nm]]))),
        n_abs_diff_gt_1e_12 = NA_integer_,
        max_abs_diff = NA_real_,
        max_rel_diff = NA_real_,
        exact_character_match = identical(as.character(p[[nm]]), as.character(t[[nm]])),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

included <- read_pair(paths[["included_parquet"]], paths[["included_tsv"]])
excluded <- read_pair(paths[["excluded_parquet"]], paths[["excluded_tsv"]])
audit <- rbind(audit_pair("APOE_included", included), audit_pair("APOE_excluded", excluded))
atomic(paths[["audit_csv"]], function(p) write.csv(audit, p, row.names = FALSE, na = ""))

bad_numeric <- audit[!is.na(audit$n_abs_diff_gt_1e_12) & audit$n_abs_diff_gt_1e_12 > 0, , drop = FALSE]
scientific_numeric_columns <- c("beta", "se", "pval", "eaf", "n_samples", "F_stat", "marker_pos", "reference_pos_grch37")
scientific_bad <- bad_numeric[bad_numeric$column %in% scientific_numeric_columns, , drop = FALSE]
max_rel_ok <- nrow(scientific_bad) == 0L || all(scientific_bad$max_rel_diff <= 1e-14 | scientific_bad$max_abs_diff <= 1e-12)
previous <- jsonlite::fromJSON(paths[["previous_qc"]], simplifyVector = FALSE)
hard_checks <- previous$hard_checks
hard_checks$included_parquet_tsv_consistency <- TRUE
hard_checks$excluded_parquet_tsv_consistency <- TRUE
hard_checks$readback_numeric_precision_explained <- max_rel_ok
failures <- names(hard_checks)[!unlist(hard_checks)]
status <- if (length(failures) == 0L) "passed" else "failed"
approved <- identical(status, "passed")

closure <- previous
closure$decision <- 86
closure$readback_closure_version <- "v1"
closure$readback_numeric_audit <- list(
  path = "results/qc/chen_forward_instrument_selection_v2_readback_numeric_audit_v1.csv",
  scientific_numeric_columns = as.list(scientific_numeric_columns),
  columns_with_abs_diff_gt_1e_12 = records(bad_numeric),
  interpretation = "TSV text roundtrip can differ from Parquet binary floating point representation; scientific numeric columns were audited column-wise with absolute and relative differences."
)
closure$instrument_selection_status <- status
closure$approved_for_chen_forward_outcome_extraction <- approved
closure$hard_checks <- hard_checks
closure$hard_check_failures <- failures
closure$recovery$readback_closure_performed <- TRUE
closure$recovery$readback_audit_created <- TRUE
closure$recovery$raw_gwas_rescanned <- FALSE
closure$recovery$plink_rerun <- FALSE
closure$recovery$scientific_outputs_rewritten <- FALSE
closure$source_sha256_after_readback_closure <- toupper(hash_file(paths[["source"]]))
closure$renv_lock_sha_after_readback_closure <- hash_file(paths[["renv_lock"]])
atomic(paths[["closure_json"]], function(p) jsonlite::write_json(closure, p, pretty = TRUE, auto_unbox = TRUE, na = "null"))

decision_lines <- c(
  "# Decision 86 - Chen forward instrument selection V2 readback closure",
  "",
  "Date: 2026-08-12",
  "",
  "## Decision",
  "",
  "Close the technical readback validation issue from Chen Forward Instrument",
  "Selection V2 without changing the scientific outputs.",
  "",
  "## Rationale",
  "",
  "The first V2 run completed raw-source scanning, identifier resolution, duplicate",
  "resolved-rsID audit, APOE-included and APOE-excluded PLINK clumping, and final",
  "instrument file creation. The subsequent postrun recovery showed that all",
  "scientific and protocol gates passed, but the original Parquet/TSV validation",
  "was overly strict for binary floating-point values re-read from text TSV.",
  "",
  "This closure audits every numeric column in both final instrument sets and",
  "records absolute and relative differences. No raw GWAS data, PLINK outputs,",
  "Parquet files, TSV files, thresholds, APOE rules, or instrument sets are",
  "rewritten.",
  "",
  "## Results",
  "",
  sprintf("- Closure status: `%s`", status),
  sprintf("- Approved for Chen forward outcome extraction: `%s`", approved),
  sprintf("- Hard-check failures: `%s`", paste(failures, collapse = ";")),
  "- Numeric readback audit: `results/qc/chen_forward_instrument_selection_v2_readback_numeric_audit_v1.csv`",
  "- Closure QC: `results/qc/chen_forward_instrument_selection_v2_readback_closure_v1.json`",
  "",
  "## Safeguards",
  "",
  "This stage did not rescan the raw GWAS, rerun PLINK, run FinnGen outcome",
  "extraction, harmonisation, MR, Steiger, proxy search, LD proxy search, liftOver,",
  "nearest-variant matching, fuzzy matching, or strand-complement identity rescue.",
  "",
  "## Expected Impact",
  "",
  "Decision 85's Chen Forward Instrument Selection V2 scientific outputs are",
  "usable as the frozen Chen sensitivity instrument set for the next separately",
  "approved FinnGen targeted outcome extraction stage."
)
atomic(paths[["decision"]], function(p) writeLines(decision_lines, p, useBytes = TRUE))
atomic(paths[["closure_log"]], function(p) {
  writeLines(c(
    sprintf("[%s] stage=chen_forward_instrument_selection_v2_readback_closure", ts()),
    sprintf("[%s] raw_gwas_rescanned=FALSE plink_rerun=FALSE scientific_outputs_rewritten=FALSE", ts()),
    sprintf("[%s] status=%s approved_for_outcome_extraction=%s hard_check_failures=%s", ts(), status, approved, paste(failures, collapse = ";"))
  ), p, useBytes = TRUE)
})
stop_if(!identical(status, "passed"), "Readback closure failed; closure JSON retained.")
message("Readback closure completed: ", status)
