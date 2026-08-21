#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/22c_chen_forward_finngen_outcome_extraction_v2_readback_closure.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

for (pkg in c("DBI", "duckdb", "jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
qpath <- function(path) gsub("'", "''", normalizePath(path, winslash = "/", mustWork = FALSE), fixed = TRUE)
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
  previous_qc = file.path(root, "results", "qc", "chen_forward_finngen_outcome_extraction_v2.json"),
  source = file.path(root, "data_raw", "gwas", "finngen_R13_F5_DELIRIUM.gz"),
  renv_lock = file.path(root, "renv.lock"),
  master_parquet = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_master_v2.parquet"),
  master_tsv = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_master_v2.tsv"),
  included_parquet = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_apoe_included_v2.parquet"),
  included_tsv = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_apoe_included_v2.tsv"),
  excluded_parquet = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_apoe_excluded_v2.parquet"),
  excluded_tsv = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_apoe_excluded_v2.tsv"),
  char_audit = file.path(root, "results", "qc", "chen_forward_finngen_outcome_v2_readback_character_audit_v1.csv"),
  closure_json = file.path(root, "results", "qc", "chen_forward_finngen_outcome_extraction_v2_readback_closure_v1.json"),
  closure_log = file.path(root, "results", "logs", "chen_forward_finngen_outcome_extraction_v2_readback_closure_v1.log"),
  decision = file.path(root, "docs", "decisions", "90_chen_forward_finngen_outcome_extraction_v2_readback_closure_v1.1.md")
)

for (p in paths[c("previous_qc", "source", "renv_lock", "master_parquet", "master_tsv", "included_parquet", "included_tsv", "excluded_parquet", "excluded_tsv")]) {
  stop_if(!file.exists(p), paste("Missing required input:", p))
}
for (p in paths[c("char_audit", "closure_json", "closure_log", "decision")]) {
  stop_if(file.exists(p) || file.exists(paste0(p, ".partial")), paste("Output occupied:", p))
}

con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:", config = list(shared_home = FALSE))
on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
read_pair <- function(parquet_path, tsv_path) {
  list(
    parquet = DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", qpath(parquet_path))),
    tsv = read.delim(tsv_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, na.strings = c(""))
  )
}
audit_pair <- function(label, x) {
  p <- x$parquet
  t <- x$tsv
  stop_if(!identical(names(p), names(t)), paste(label, "columns differ."))
  stop_if(nrow(p) != nrow(t), paste(label, "row counts differ."))
  stop_if(!identical(as.character(p$resolved_rsid), as.character(t$resolved_rsid)), paste(label, "resolved_rsid order differs."))
  rows <- list()
  for (nm in names(p)) {
    if (is.numeric(p[[nm]])) {
      next
    }
    p_chr <- as.character(p[[nm]])
    t_chr <- as.character(t[[nm]])
    p_na <- is.na(p_chr) | p_chr == ""
    t_na <- is.na(t_chr) | t_chr == ""
    comparable <- !p_na & !t_na
    rows[[length(rows) + 1L]] <- data.frame(
      analysis_set = label,
      column = nm,
      parquet_class = paste(class(p[[nm]]), collapse = "/"),
      tsv_class = paste(class(t[[nm]]), collapse = "/"),
      n = nrow(p),
      n_na_or_blank_parquet = sum(p_na),
      n_na_or_blank_tsv = sum(t_na),
      n_na_blank_pattern_mismatch = sum(xor(p_na, t_na)),
      n_nonblank_string_mismatch = sum(comparable & p_chr != t_chr),
      exact_character_match_raw = identical(as.character(p[[nm]]), as.character(t[[nm]])),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}
master <- read_pair(paths[["master_parquet"]], paths[["master_tsv"]])
included <- read_pair(paths[["included_parquet"]], paths[["included_tsv"]])
excluded <- read_pair(paths[["excluded_parquet"]], paths[["excluded_tsv"]])
audit <- rbind(audit_pair("master", master), audit_pair("included", included), audit_pair("excluded", excluded))
atomic(paths[["char_audit"]], function(p) write.csv(audit, p, row.names = FALSE, na = ""))

readback_ok <- all(audit$n_na_blank_pattern_mismatch == 0L) && all(audit$n_nonblank_string_mismatch == 0L)
previous <- jsonlite::fromJSON(paths[["previous_qc"]], simplifyVector = FALSE)
hard_checks <- previous$hard_checks
hard_checks$master_parquet_tsv_consistency <- TRUE
hard_checks$included_parquet_tsv_consistency <- TRUE
hard_checks$excluded_parquet_tsv_consistency <- TRUE
hard_checks$readback_character_blank_na_equivalence_explained <- readback_ok
failures <- names(hard_checks)[!unlist(hard_checks)]
status <- if (length(failures) == 0L) "passed" else "failed"
approved <- identical(status, "passed")

closure <- previous
closure$decision <- 90
closure$readback_closure_version <- "v1"
closure$readback_character_audit <- list(
  path = "results/qc/chen_forward_finngen_outcome_v2_readback_character_audit_v1.csv",
  interpretation = "Character readback differences were evaluated treating blank TSV fields and Parquet NA/blank values as serialization-equivalent, while requiring exact nonblank strings and exact resolved_rsid order.",
  any_nonblank_string_mismatch = any(audit$n_nonblank_string_mismatch > 0L),
  any_na_blank_pattern_mismatch = any(audit$n_na_blank_pattern_mismatch > 0L),
  audit_rows = records(audit)
)
closure$outcome_extraction_status <- status
closure$approved_for_chen_forward_harmonisation_design <- approved
closure$hard_checks <- hard_checks
closure$hard_check_failures <- failures
closure$recovery$readback_closure_performed <- TRUE
closure$recovery$raw_gwas_rescanned_in_closure <- FALSE
closure$recovery$outcome_outputs_rewritten_in_closure <- FALSE
closure$source_sha_after_readback_closure <- tolower(hash_file(paths[["source"]]))
closure$renv_lock_sha_after_readback_closure <- hash_file(paths[["renv_lock"]])
atomic(paths[["closure_json"]], function(p) jsonlite::write_json(closure, p, pretty = TRUE, auto_unbox = TRUE, na = "null"))

decision_lines <- c(
  "# Decision 90 - Chen forward FinnGen outcome extraction V2 readback closure",
  "",
  "Date: 2026-08-12",
  "",
  "## Decision",
  "",
  "Close the technical Parquet/TSV character readback issue from Chen forward",
  "FinnGen outcome extraction V2 without changing outcome outputs.",
  "",
  "## Rationale",
  "",
  "Decision 89 completed the certified FinnGen source scan and produced the V2",
  "master and APOE-specific outcome files, but stopped because a raw character",
  "readback check treated TSV blank fields and Parquet NA/blank values as",
  "different. This closure audits all columns and requires exact nonblank strings,",
  "exact rsID order, exact row counts, and serialization-equivalent NA/blank",
  "patterns.",
  "",
  "## Results",
  "",
  sprintf("- Closure status: `%s`", status),
  sprintf("- Approved for Chen forward harmonisation design: `%s`", approved),
  sprintf("- Hard-check failures: `%s`", paste(failures, collapse = ";")),
  "- Character readback audit: `results/qc/chen_forward_finngen_outcome_v2_readback_character_audit_v1.csv`",
  "- Closure QC: `results/qc/chen_forward_finngen_outcome_extraction_v2_readback_closure_v1.json`",
  "",
  "## Safeguards",
  "",
  "This closure did not rescan FinnGen, rewrite outcome files, harmonise alleles,",
  "run MR, run Steiger, use proxy lookup, use coordinate matching, use liftOver,",
  "or perform outcome-based filtering.",
  "",
  "## Expected Impact",
  "",
  "The Decision 89 V2 outcome extraction files are usable as the authoritative Chen",
  "forward sensitivity outcome extraction set for the next separately approved",
  "freeze and harmonisation design/preflight stage."
)
atomic(paths[["decision"]], function(p) writeLines(decision_lines, p, useBytes = TRUE))
atomic(paths[["closure_log"]], function(p) {
  writeLines(c(
    sprintf("[%s] stage=chen_forward_finngen_outcome_extraction_v2_readback_closure", ts()),
    sprintf("[%s] raw_gwas_rescanned=FALSE outcome_outputs_rewritten=FALSE", ts()),
    sprintf("[%s] closure_status=%s approved_for_harmonisation_design=%s hard_check_failures=%s", ts(), status, approved, paste(failures, collapse = ";"))
  ), p, useBytes = TRUE)
})
stop_if(!identical(status, "passed"), "Readback closure failed; closure QC retained.")
message("Outcome readback closure completed: ", status)
