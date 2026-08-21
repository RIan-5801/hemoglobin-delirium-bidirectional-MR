#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/16_reverse_relaxed_outcome_extraction_v1_freeze_manifest.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

for (pkg in c("DBI", "duckdb", "jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
hash_file <- function(path) digest::digest(file = path, algo = "sha256")
sql_string <- function(path, must_work = TRUE) paste0("'", gsub("'", "''", normalizePath(path, winslash = "/", mustWork = must_work), fixed = TRUE), "'")
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
num_equal <- function(a, b, tol = 1e-12) {
  a <- as.numeric(a); b <- as.numeric(b)
  (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tol)
}

out <- c(
  manifest = file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_outcome_extraction_v1_freeze_manifest.csv"),
  json = file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_outcome_extraction_v1_freeze.json"),
  log = file.path(root, "results", "logs", "vuckovic_hb_reverse_relaxed_outcome_extraction_v1_freeze.log")
)
stop_if(any(file.exists(c(out, paste0(out, ".partial")))), "A relaxed extraction freeze final or partial target exists; refusing to overwrite.")
dir.create(dirname(out[["manifest"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["log"]]), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = out[["log"]], append = TRUE)

main <- function() {
  log_line("stage=reverse_relaxed_outcome_extraction_v1_freeze")
  rel <- c(
    "docs/decisions/59_finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze_v1.1.md",
    "docs/decisions/60_vuckovic_hb_reverse_relaxed_outcome_extraction_v1_v1.1.md",
    "docs/decisions/61_vuckovic_hb_reverse_relaxed_outcome_extraction_v1_freeze_v1.1.md",
    "R/16_extract_vuckovic_hb_reverse_relaxed_outcomes_v1.R",
    "data_derived/reverse_outcome_extraction/finngen_r13_delirium_reverse_relaxed_targets_v1.tsv",
    "data_derived/reverse_outcome_extraction/vuckovic_hb_reverse_relaxed_all_matches_v1.parquet",
    "data_derived/reverse_outcome_extraction/vuckovic_hb_reverse_relaxed_all_matches_v1.tsv",
    "data_derived/reverse_outcome_extraction/vuckovic_hb_reverse_relaxed_unique_matches_v1.parquet",
    "data_derived/reverse_outcome_extraction/vuckovic_hb_reverse_relaxed_unique_matches_v1.tsv",
    "results/qc/vuckovic_hb_reverse_relaxed_target_match_status_v1.csv",
    "results/qc/vuckovic_hb_reverse_relaxed_outcome_extraction_v1.json",
    "results/logs/vuckovic_hb_reverse_relaxed_outcome_extraction_v1.log",
    "results/qc/finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze.json",
    "results/qc/finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze_manifest.csv",
    "renv.lock"
  )
  paths <- file.path(root, rel)
  stop_if(any(!file.exists(paths)), paste("Missing freeze input(s):", paste(rel[!file.exists(paths)], collapse = "; ")))

  qc <- jsonlite::fromJSON(file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_outcome_extraction_v1.json"), simplifyVector = FALSE)
  instrument_freeze <- jsonlite::fromJSON(file.path(root, "results", "qc", "finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze.json"), simplifyVector = FALSE)
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  read_pq <- function(path) DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(path)))
  targets <- read.delim(file.path(root, "data_derived", "reverse_outcome_extraction", "finngen_r13_delirium_reverse_relaxed_targets_v1.tsv"), check.names = FALSE)
  status <- read.csv(file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_target_match_status_v1.csv"), check.names = FALSE)
  all_pq <- read_pq(file.path(root, "data_derived", "reverse_outcome_extraction", "vuckovic_hb_reverse_relaxed_all_matches_v1.parquet"))
  unique_pq <- read_pq(file.path(root, "data_derived", "reverse_outcome_extraction", "vuckovic_hb_reverse_relaxed_unique_matches_v1.parquet"))
  all_tsv <- read.delim(file.path(root, "data_derived", "reverse_outcome_extraction", "vuckovic_hb_reverse_relaxed_all_matches_v1.tsv"), check.names = FALSE)
  unique_tsv <- read.delim(file.path(root, "data_derived", "reverse_outcome_extraction", "vuckovic_hb_reverse_relaxed_unique_matches_v1.tsv"), check.names = FALSE)

  failures <- character()
  add_fail <- function(x) failures <<- unique(c(failures, x))
  if (!identical(qc$extraction_status, "passed")) add_fail("extraction_status_not_passed")
  if (!isTRUE(qc$approved_for_reverse_relaxed_harmonisation_preflight)) add_fail("not_approved_for_preflight")
  if (length(qc$hard_check_failures) != 0L) add_fail("extraction_hard_check_failures_not_empty")
  if (!identical(qc$analysis_direction, "delirium_to_Hb") || !identical(qc$analysis_role, "secondary_reverse_exploratory_relaxed")) add_fail("direction_or_role_mismatch")
  if (!identical(qc$proxy_used, FALSE) || !identical(qc$liftover_used, FALSE) || !identical(qc$harmonisation_performed, FALSE)) add_fail("prohibited_action_flag_not_false")
  if (!identical(instrument_freeze$freeze_status, "passed")) add_fail("upstream_instrument_freeze_not_passed")
  if (!setequal(targets$target_rsid, status$target_rsid) || !setequal(targets$target_rsid, unique_pq$target_rsid)) add_fail("target_set_mismatch")
  if (nrow(targets) != as.integer(qc$union_target_count)) add_fail("union_target_count_mismatch")
  if (sum(targets$included_member) != as.integer(qc$included_target_count)) add_fail("included_target_count_mismatch")
  if (sum(targets$excluded_member) != as.integer(qc$excluded_target_count)) add_fail("excluded_target_count_mismatch")
  if (sum(targets$included_member & targets$excluded_member) != as.integer(qc$shared_target_count)) add_fail("shared_target_count_mismatch")
  if (sum(status$match_class == "unique_match") != as.integer(qc$unique_match_count)) add_fail("unique_match_count_mismatch")
  if (sum(status$match_class == "missing_match") != as.integer(qc$missing_match_count)) add_fail("missing_match_count_mismatch")
  if (sum(status$match_class == "multiple_match") != as.integer(qc$multiple_match_count)) add_fail("multiple_match_count_mismatch")
  if (sum(status$included_member & status$match_class == "unique_match") != as.integer(qc$included_unique_match_count)) add_fail("included_unique_match_count_mismatch")
  if (sum(status$excluded_member & status$match_class == "unique_match") != as.integer(qc$excluded_unique_match_count)) add_fail("excluded_unique_match_count_mismatch")
  if (sum(status$included_member & status$excluded_member & status$match_class == "unique_match") != as.integer(qc$shared_unique_match_count)) add_fail("shared_unique_match_count_mismatch")
  if (!identical(names(all_pq), names(all_tsv)) || nrow(all_pq) != nrow(all_tsv)) add_fail("all_matches_parquet_tsv_consistency_failed")
  if (!identical(names(unique_pq), names(unique_tsv)) || nrow(unique_pq) != nrow(unique_tsv)) add_fail("unique_matches_parquet_tsv_consistency_failed")
  for (col in c("target_rsid", "outcome_beta_raw", "outcome_se_raw", "outcome_lp_raw", "outcome_pval_raw", "outcome_eaf_raw")) {
    if (col %in% names(unique_pq)) {
      if (is.numeric(unique_pq[[col]]) || is.numeric(unique_tsv[[col]])) {
        if (!all(num_equal(unique_pq[[col]], unique_tsv[[col]]))) add_fail(paste0("unique_numeric_mismatch_", col))
      } else if (!identical(as.character(unique_pq[[col]]), as.character(unique_tsv[[col]]))) {
        add_fail(paste0("unique_character_mismatch_", col))
      }
    }
  }

  manifest <- data.frame(
    file_role = c("instrument_freeze_decision", "extraction_decision", "extraction_freeze_decision", "extraction_script",
                  "target_registry", "all_matches_parquet", "all_matches_tsv", "unique_matches_parquet",
                  "unique_matches_tsv", "target_status_csv", "extraction_qc_json", "extraction_log",
                  "instrument_freeze_json", "instrument_freeze_manifest", "renv_lock"),
    relative_path = rel,
    file_size_bytes = as.numeric(file.info(paths)$size),
    sha256 = vapply(paths, hash_file, character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(manifest, paste0(out[["manifest"]], ".partial"), row.names = FALSE)
  manifest_sha <- hash_file(paste0(out[["manifest"]], ".partial"))
  freeze_status <- if (length(failures) == 0L) "passed" else "failed"
  freeze <- list(
    freeze_version = "v1",
    decision = 61,
    authoritative_reverse_relaxed_outcome_extraction_version = "v1",
    source_relaxed_instrument_version = "v2",
    analysis_direction = "delirium_to_Hb",
    analysis_role = "secondary_reverse_exploratory_relaxed",
    branch_type = "protocol_prespecified_fallback",
    p_threshold = 5e-6,
    included_target_count = sum(targets$included_member),
    excluded_target_count = sum(targets$excluded_member),
    shared_target_count = sum(targets$included_member & targets$excluded_member),
    union_target_count = nrow(targets),
    unique_match_count = sum(status$match_class == "unique_match"),
    missing_match_count = sum(status$match_class == "missing_match"),
    multiple_match_count = sum(status$match_class == "multiple_match"),
    included_unique_match_count = sum(status$included_member & status$match_class == "unique_match"),
    excluded_unique_match_count = sum(status$excluded_member & status$match_class == "unique_match"),
    shared_unique_match_count = sum(status$included_member & status$excluded_member & status$match_class == "unique_match"),
    palindromic_count = sum(status$palindromic_snp),
    outcome_source_path = qc$outcome_source_path,
    outcome_source_sha256 = qc$source_sha$certified,
    outcome_genome_build = qc$outcome_genome_build,
    outcome_effect_scale = qc$outcome_effect_scale,
    outcome_n_study = qc$outcome_n_study,
    manifest_sha256 = manifest_sha,
    freeze_status = freeze_status,
    approved_for_reverse_relaxed_harmonisation_preflight = identical(freeze_status, "passed"),
    hard_checks = list(
      extraction_status_passed = identical(qc$extraction_status, "passed"),
      extraction_approved_for_preflight = isTRUE(qc$approved_for_reverse_relaxed_harmonisation_preflight),
      extraction_hard_check_failures_empty = length(qc$hard_check_failures) == 0L,
      no_proxy = identical(qc$proxy_used, FALSE),
      no_liftover = identical(qc$liftover_used, FALSE),
      no_harmonisation = identical(qc$harmonisation_performed, FALSE),
      target_registry_status_unique_consistent = !("target_set_mismatch" %in% failures),
      all_matches_consistency = !("all_matches_parquet_tsv_consistency_failed" %in% failures),
      unique_matches_consistency = !("unique_matches_parquet_tsv_consistency_failed" %in% failures),
      no_vuckovic_rescan = TRUE,
      no_reverse_mr = TRUE
    ),
    hard_check_failures = failures
  )
  jsonlite::write_json(freeze, paste0(out[["json"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  if (!identical(freeze_status, "passed")) stop("Relaxed outcome extraction freeze failed; partial outputs retained.", call. = FALSE)
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
