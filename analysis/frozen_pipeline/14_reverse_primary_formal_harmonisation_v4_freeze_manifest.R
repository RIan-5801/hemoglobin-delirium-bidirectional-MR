#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/14_reverse_primary_formal_harmonisation_v4_freeze_manifest.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

for (pkg in c("digest", "jsonlite")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
hash_file <- function(path) digest::digest(file = path, algo = "sha256")
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

out <- c(
  manifest = file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_formal_harmonisation_v4_freeze_manifest.csv"),
  json = file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_formal_harmonisation_v4_freeze.json"),
  log = file.path(root, "results", "logs", "vuckovic_hb_reverse_primary_formal_harmonisation_v4_freeze.log")
)
stop_if(any(file.exists(c(out, paste0(out, ".partial")))), "A freeze final or partial target exists; refusing to overwrite.")
dir.create(dirname(out[["manifest"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["log"]]), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = out[["log"]], append = TRUE)

main <- function() {
  log_line("stage=reverse_primary_formal_harmonisation_v4_freeze")
  rel <- c(
    "docs/decisions/55_vuckovic_hb_reverse_primary_formal_harmonisation_v4_tsv_schema_fix_v1.1.md",
    "docs/decisions/56_vuckovic_hb_reverse_primary_formal_harmonisation_v4_freeze_v1.1.md",
    "R/14_reverse_primary_formal_harmonisation_v4.R",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_primary_harmonised_apoe_included_v4.parquet",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_primary_harmonised_apoe_included_v4.tsv",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_primary_harmonised_apoe_excluded_v4.parquet",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_primary_harmonised_apoe_excluded_v4.tsv",
    "results/qc/vuckovic_hb_reverse_primary_formal_harmonisation_counts_v4.csv",
    "results/qc/vuckovic_hb_reverse_primary_formal_harmonisation_v4.json",
    "results/logs/vuckovic_hb_reverse_primary_formal_harmonisation_v4.log",
    "renv.lock"
  )
  paths <- file.path(root, rel)
  stop_if(any(!file.exists(paths)), paste("Missing freeze input(s):", paste(rel[!file.exists(paths)], collapse = "; ")))

  qc <- jsonlite::fromJSON(file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_formal_harmonisation_v4.json"), simplifyVector = FALSE)
  failures <- character()
  add_fail <- function(x) failures <<- unique(c(failures, x))
  if (!identical(qc$harmonisation_version, "v4")) add_fail("harmonisation_version_not_v4")
  if (!identical(qc$harmonisation_status, "passed")) add_fail("harmonisation_status_not_passed")
  if (length(qc$hard_check_failures) != 0L) add_fail("hard_check_failures_present")
  if (!identical(qc$analysis_direction, "delirium_to_Hb")) add_fail("analysis_direction_mismatch")
  if (as.integer(qc$included_final_valid_instrument_count) != 1L) add_fail("included_final_valid_count_mismatch")
  if (as.integer(qc$excluded_final_valid_instrument_count) != 1L) add_fail("excluded_final_valid_count_mismatch")
  if (!isTRUE(qc$overall_relaxed_threshold_trigger)) add_fail("relaxed_threshold_trigger_not_true")
  if (!isTRUE(qc$approved_for_reverse_relaxed_threshold_branch)) add_fail("relaxed_branch_not_approved")
  if (!identical(qc$relaxed_threshold_status, "triggered_not_started")) add_fail("relaxed_threshold_status_mismatch")
  if (abs(as.numeric(qc$relaxed_threshold) - 5e-6) > 0) add_fail("relaxed_threshold_mismatch")
  if (!isTRUE(qc$renv_lock_unchanged)) add_fail("renv_lock_changed_in_source_qc")
  if (!isTRUE(qc$included_parquet_tsv_consistency)) add_fail("included_parquet_tsv_consistency_not_true")
  if (!isTRUE(qc$excluded_parquet_tsv_consistency)) add_fail("excluded_parquet_tsv_consistency_not_true")

  manifest <- data.frame(
    file_role = c("decision_source", "freeze_decision", "script_source", "harmonised_included_parquet",
                  "harmonised_included_tsv", "harmonised_excluded_parquet", "harmonised_excluded_tsv",
                  "counts_source", "qc_source", "log_source", "renv_lock"),
    relative_path = rel,
    file_size_bytes = as.numeric(file.info(paths)$size),
    sha256 = vapply(paths, hash_file, character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(manifest, paste0(out[["manifest"]], ".partial"), row.names = FALSE)
  manifest_sha <- hash_file(paste0(out[["manifest"]], ".partial"))
  freeze_status <- if (length(failures) == 0L) "passed" else "failed"
  result <- list(
    freeze_version = "v1",
    decision = 56,
    authoritative_formal_harmonisation_version = "v4",
    source_decision = 55,
    analysis_direction = "delirium_to_Hb",
    analysis_role = "secondary_reverse_primary",
    harmonisation_status = qc$harmonisation_status,
    included_final_valid_instrument_count = qc$included_final_valid_instrument_count,
    excluded_final_valid_instrument_count = qc$excluded_final_valid_instrument_count,
    primary_relaxed_threshold_trigger_threshold = qc$primary_relaxed_threshold_trigger_threshold,
    overall_relaxed_threshold_trigger = qc$overall_relaxed_threshold_trigger,
    relaxed_threshold = qc$relaxed_threshold,
    relaxed_threshold_status = qc$relaxed_threshold_status,
    relaxed_threshold_trigger_reason = qc$relaxed_threshold_trigger_reason,
    approved_for_reverse_relaxed_threshold_branch = qc$approved_for_reverse_relaxed_threshold_branch,
    no_relaxed_threshold_analysis_started_by_freeze = TRUE,
    no_reverse_mr = TRUE,
    no_proxy = TRUE,
    no_liftover = TRUE,
    manifest_sha256 = manifest_sha,
    freeze_status = freeze_status,
    approved_for_reverse_relaxed_threshold_instrument_selection = identical(freeze_status, "passed"),
    hard_check_failures = failures
  )
  jsonlite::write_json(result, paste0(out[["json"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  if (!identical(freeze_status, "passed")) stop("Freeze checks failed; partial outputs retained.", call. = FALSE)
  for (path in out[c("manifest", "json")]) {
    stop_if(file.exists(path), paste("Output appeared during run:", path))
    stop_if(!file.rename(paste0(path, ".partial"), path), paste("Atomic rename failed:", path))
  }
  log_line("freeze_status=passed manifest_sha256=", manifest_sha)
}

tryCatch(main(), error = function(e) {
  log_line("status=failed error=", conditionMessage(e))
  quit(status = 1L)
})
