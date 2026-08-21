#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
pkgs <- c("digest", "DBI", "duckdb", "jsonlite")
if (any(!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))) stop("Required installed package missing; no automatic installation.", call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") stop("Usage: Rscript.exe R/13_reverse_primary_harmonisation_preflight_v1_freeze_manifest.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

stop_if <- function(x, m) if (isTRUE(x)) stop(m, call. = FALSE)
hash_file <- function(p) digest::digest(file = p, algo = "sha256", serialize = FALSE)
sql_string <- function(p, must_work = TRUE) paste0("'", gsub("'", "''", normalizePath(p, winslash = "/", mustWork = must_work), fixed = TRUE), "'")
num_equal <- function(a, b, tol = 1e-12) {
  a <- as.numeric(a); b <- as.numeric(b)
  (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tol)
}

paths <- list(
  script = file.path(root, "R", "13_reverse_primary_harmonisation_preflight_v1.R"),
  pq = file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_primary_harmonisation_preflight_v1.parquet"),
  tsv = file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_primary_harmonisation_preflight_v1.tsv"),
  json = file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_harmonisation_preflight_v1.json"),
  log = file.path(root, "results", "logs", "vuckovic_hb_reverse_primary_harmonisation_preflight_v1.log"),
  decision50 = file.path(root, "docs", "decisions", "50_vuckovic_hb_reverse_primary_harmonisation_preflight_v1_v1.1.md"),
  extraction_freeze = file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_outcome_extraction_v1_freeze.json")
)
out <- c(
  manifest = file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_harmonisation_preflight_v1_freeze_manifest.csv"),
  json = file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_harmonisation_preflight_v1_freeze.json"),
  log = file.path(root, "results", "logs", "vuckovic_hb_reverse_primary_harmonisation_preflight_v1_freeze.log")
)
stop_if(any(!file.exists(unlist(paths))), "Required preflight or upstream authority file is missing.")
stop_if(any(file.exists(c(out, paste0(out, ".partial")))), "A preflight freeze final or partial target exists; refusing to overwrite.")
dir.create(dirname(out[["manifest"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["log"]]), recursive = TRUE, showWarnings = FALSE)
log_line <- function(x) cat(sprintf("%s %s\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), x), file = out[["log"]], append = TRUE)

main <- function() {
  log_line("stage=reverse_primary_harmonisation_preflight_v1_freeze")
  qc <- jsonlite::fromJSON(paths$json, simplifyVector = FALSE)
  extraction_freeze <- jsonlite::fromJSON(paths$extraction_freeze, simplifyVector = FALSE)
  hard_check_failures <- character()
  add_fail <- function(x) hard_check_failures <<- unique(c(hard_check_failures, x))

  if (!identical(qc$preflight_status, "passed")) add_fail("preflight_status_not_passed")
  if (!isTRUE(qc$approved_for_reverse_primary_formal_harmonisation)) add_fail("not_approved_for_formal_harmonisation")
  if (length(qc$hard_check_failures) != 0L) add_fail("preflight_hard_check_failures_not_empty")
  if (!identical(qc$harmonisation_performed, FALSE)) add_fail("harmonisation_already_performed")
  if (as.integer(qc$beta_flip_performed_count) != 0L) add_fail("beta_flip_performed_count_nonzero")
  if (as.integer(qc$eaf_flip_performed_count) != 0L) add_fail("eaf_flip_performed_count_nonzero")
  if (!identical(qc$relaxed_threshold_status, "not_triggered")) add_fail("relaxed_threshold_not_triggered_failed")
  if (!identical(extraction_freeze$freeze_status, "passed")) add_fail("outcome_extraction_freeze_not_passed")

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  pq <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(paths$pq)))
  tsv <- read.delim(paths$tsv, check.names = FALSE)
  if (!identical(names(pq), names(tsv))) add_fail("parquet_tsv_schema_mismatch")
  if (nrow(pq) != nrow(tsv)) add_fail("parquet_tsv_row_count_mismatch")
  for (col in intersect(names(pq), names(tsv))) {
    if (is.numeric(pq[[col]]) || is.numeric(tsv[[col]])) {
      if (!all(num_equal(pq[[col]], tsv[[col]]))) add_fail(paste0("numeric_mismatch_", col))
    } else if (!identical(as.character(pq[[col]]), as.character(tsv[[col]]))) {
      add_fail(paste0("character_mismatch_", col))
    }
  }

  required_cols <- c("analysis_set", "target_rsid", "exposure_effect_allele", "exposure_other_allele",
                     "exposure_beta", "exposure_se", "exposure_pval", "exposure_eaf",
                     "outcome_effect_allele_raw", "outcome_other_allele_raw", "outcome_beta_raw",
                     "outcome_se_raw", "outcome_pval_raw", "outcome_eaf_raw", "palindromic_snp",
                     "raw_orientation_class", "beta_flip_required", "eaf_flip_required",
                     "strand_flip_required", "eligible_for_formal_harmonisation")
  if (!all(required_cols %in% names(pq))) add_fail("required_preflight_columns_missing")

  counts <- list(
    included_preflight_count = sum(pq$analysis_set == "APOE included"),
    excluded_preflight_count = sum(pq$analysis_set == "APOE excluded"),
    included_eligible_for_formal_harmonisation_count = sum(pq$analysis_set == "APOE included" & pq$eligible_for_formal_harmonisation),
    excluded_eligible_for_formal_harmonisation_count = sum(pq$analysis_set == "APOE excluded" & pq$eligible_for_formal_harmonisation),
    palindromic_count = sum(pq$palindromic_snp),
    exact_match_count = sum(pq$raw_orientation_class == "exact_match"),
    swapped_match_count = sum(pq$raw_orientation_class == "swapped_match"),
    strand_exact_match_count = sum(pq$raw_orientation_class == "strand_exact_match"),
    strand_swapped_match_count = sum(pq$raw_orientation_class == "strand_swapped_match"),
    incompatible_count = sum(pq$raw_orientation_class == "incompatible"),
    invalid_count = sum(pq$raw_orientation_class == "invalid"),
    beta_flip_required_count = sum(pq$beta_flip_required),
    eaf_flip_required_count = sum(pq$eaf_flip_required)
  )
  hard_checks <- list(
    preflight_status_passed = identical(qc$preflight_status, "passed"),
    approved_for_reverse_primary_formal_harmonisation = isTRUE(qc$approved_for_reverse_primary_formal_harmonisation),
    hard_check_failures_empty = length(qc$hard_check_failures) == 0L,
    no_harmonisation_performed = identical(qc$harmonisation_performed, FALSE),
    no_beta_flip_performed = as.integer(qc$beta_flip_performed_count) == 0L,
    no_eaf_flip_performed = as.integer(qc$eaf_flip_performed_count) == 0L,
    relaxed_threshold_not_started = identical(qc$relaxed_threshold_status, "not_triggered"),
    parquet_tsv_consistency = length(grep("mismatch", hard_check_failures)) == 0L,
    outcome_extraction_freeze_passed = identical(extraction_freeze$freeze_status, "passed"),
    no_source_scan = TRUE,
    no_reverse_mr = TRUE
  )
  freeze_status <- if (length(hard_check_failures) == 0L && all(unlist(hard_checks))) "passed" else "failed"
  approved <- identical(freeze_status, "passed")

  rel <- c(
    "R/13_reverse_primary_harmonisation_preflight_v1.R",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_primary_harmonisation_preflight_v1.parquet",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_primary_harmonisation_preflight_v1.tsv",
    "results/qc/vuckovic_hb_reverse_primary_harmonisation_preflight_v1.json",
    "results/logs/vuckovic_hb_reverse_primary_harmonisation_preflight_v1.log",
    "docs/decisions/50_vuckovic_hb_reverse_primary_harmonisation_preflight_v1_v1.1.md",
    "results/qc/vuckovic_hb_reverse_primary_outcome_extraction_v1_freeze.json"
  )
  manifest <- data.frame(file_role = "authority", relative_path = rel,
                         file_size_bytes = file.info(file.path(root, rel))$size,
                         sha256 = vapply(file.path(root, rel), hash_file, character(1)),
                         stringsAsFactors = FALSE)
  write.csv(manifest, paste0(out[["manifest"]], ".partial"), row.names = FALSE)
  manifest_sha256 <- hash_file(paste0(out[["manifest"]], ".partial"))
  js <- c(list(
    freeze_version = "v1",
    authoritative_reverse_primary_harmonisation_preflight_version = "v1",
    analysis_direction = "delirium_to_Hb",
    included_target_rsid = pq$target_rsid[pq$analysis_set == "APOE included"],
    excluded_target_rsid = pq$target_rsid[pq$analysis_set == "APOE excluded"],
    outcome_extraction_freeze_manifest_sha256 = extraction_freeze$manifest_sha256,
    manifest_sha256 = manifest_sha256,
    freeze_status = freeze_status,
    approved_for_reverse_primary_formal_harmonisation = approved,
    hard_checks = hard_checks,
    hard_check_failures = hard_check_failures
  ), counts)
  jsonlite::write_json(js, paste0(out[["json"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  for (p in out[c("manifest", "json")]) {
    stop_if(file.exists(p), paste("Output appeared during run:", p))
    stop_if(!file.rename(paste0(p, ".partial"), p), paste("Atomic rename failed:", p))
  }
  log_line(sprintf("SUCCESS freeze_status=%s manifest_sha256=%s hard_check_failures=%s", freeze_status, manifest_sha256, paste(hard_check_failures, collapse = ";")))
  if (!approved) stop("Preflight freeze did not pass; formal harmonisation must not proceed.", call. = FALSE)
}

tryCatch(main(), error = function(e) {
  log_line(paste0("TERMINATED_PRIMARY: ", conditionMessage(e)))
  quit(status = 1L)
})

