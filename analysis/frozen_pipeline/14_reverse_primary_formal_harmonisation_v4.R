#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
pkgs <- c("digest", "DBI", "duckdb", "jsonlite")
if (any(!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))) stop("Required installed package missing; no automatic installation.", call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") stop("Usage: Rscript.exe R/14_reverse_primary_formal_harmonisation_v4.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

stop_if <- function(x, m) if (isTRUE(x)) stop(m, call. = FALSE)
hash_file <- function(p) digest::digest(file = p, algo = "sha256", serialize = FALSE)
sql_string <- function(p, must_work = TRUE) paste0("'", gsub("'", "''", normalizePath(p, winslash = "/", mustWork = must_work), fixed = TRUE), "'")
sql_ident <- function(con, x) as.character(DBI::dbQuoteIdentifier(con, x))
num_equal <- function(a, b, tol = 1e-12) {
  a <- as.numeric(a); b <- as.numeric(b)
  (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tol)
}
canon_char <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}
records <- function(x) lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))

preflight_freeze_json <- file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_harmonisation_preflight_v1_freeze.json")
preflight_pq <- file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_primary_harmonisation_preflight_v1.parquet")
renv_lock <- file.path(root, "renv.lock")
out <- c(
  inc_pq = file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_primary_harmonised_apoe_included_v4.parquet"),
  inc_tsv = file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_primary_harmonised_apoe_included_v4.tsv"),
  exc_pq = file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_primary_harmonised_apoe_excluded_v4.parquet"),
  exc_tsv = file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_primary_harmonised_apoe_excluded_v4.tsv"),
  counts = file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_formal_harmonisation_counts_v4.csv"),
  json = file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_formal_harmonisation_v4.json"),
  log = file.path(root, "results", "logs", "vuckovic_hb_reverse_primary_formal_harmonisation_v4.log")
)
stop_if(any(!file.exists(c(preflight_freeze_json, preflight_pq, renv_lock))), "Required preflight freeze or renv.lock input missing.")
stop_if(any(file.exists(c(out, paste0(out, ".partial")))), "A formal harmonisation V4 final or partial target exists; refusing to overwrite.")
dir.create(dirname(out[["inc_pq"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["json"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["log"]]), recursive = TRUE, showWarnings = FALSE)
log_line <- function(x) cat(sprintf("%s %s\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), x), file = out[["log"]], append = TRUE)

write_pq <- function(con, x, path) {
  nm <- paste0("tmp_", digest::digest(path, algo = "xxhash32", serialize = FALSE))
  DBI::dbWriteTable(con, nm, x, temporary = TRUE, overwrite = TRUE)
  DBI::dbExecute(con, sprintf("COPY %s TO %s (FORMAT PARQUET, COMPRESSION ZSTD)", sql_ident(con, nm), sql_string(path, must_work = FALSE)))
  DBI::dbRemoveTable(con, nm)
}
check_pq_tsv <- function(con, pq_path, tsv_path) {
  pq <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(pq_path)))
  tsv <- read.delim(tsv_path, check.names = FALSE, colClasses = "character", na.strings = character())
  ok <- identical(names(pq), names(tsv)) && nrow(pq) == nrow(tsv) && !anyDuplicated(names(pq))
  if (!ok) return(FALSE)
  for (col in names(pq)) {
    if (is.numeric(pq[[col]])) {
      ok <- ok && all(num_equal(pq[[col]], tsv[[col]]))
    } else if (is.logical(pq[[col]])) {
      ok <- ok && identical(toupper(canon_char(pq[[col]])), toupper(canon_char(tsv[[col]])))
    } else {
      ok <- ok && identical(canon_char(pq[[col]]), canon_char(tsv[[col]]))
    }
  }
  ok
}

main <- function() {
  log_line("stage=reverse_primary_formal_harmonisation_v4")
  renv_before <- hash_file(renv_lock)
  freeze <- jsonlite::fromJSON(preflight_freeze_json, simplifyVector = FALSE)
  stop_if(!identical(freeze$freeze_status, "passed"), "Preflight freeze gate is not passed.")
  stop_if(!isTRUE(freeze$approved_for_reverse_primary_formal_harmonisation), "Preflight freeze did not approve formal harmonisation.")

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  pre <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(preflight_pq)))
  stop_if(nrow(pre) != 2L, "Expected exactly two preflight rows.")
  stop_if(sum(pre$analysis_set == "APOE included") != 1L || sum(pre$analysis_set == "APOE excluded") != 1L, "Included/excluded preflight rows incomplete.")

  h <- pre
  h$rsid <- h$target_rsid
  h$orientation_class <- h$raw_orientation_class
  h$outcome_effect_allele_harmonised <- h$outcome_effect_allele_raw
  h$outcome_other_allele_harmonised <- h$outcome_other_allele_raw
  h$outcome_beta_harmonised <- h$outcome_beta_raw
  h$outcome_se_harmonised <- h$outcome_se_raw
  h$outcome_pval_harmonised <- h$outcome_pval_raw
  h$outcome_eaf_harmonised <- h$outcome_eaf_raw
  h$beta_flipped <- FALSE
  h$eaf_flipped <- FALSE
  h$strand_flipped <- FALSE
  h$record_excluded <- FALSE
  h$exclusion_reason <- ""

  swapped <- h$orientation_class %in% c("swapped_match", "strand_swapped_match")
  strand <- h$orientation_class %in% c("strand_exact_match", "strand_swapped_match")
  incompatible <- h$orientation_class %in% c("incompatible", "invalid")
  if (any(swapped)) {
    h$outcome_effect_allele_harmonised[swapped] <- h$outcome_other_allele_raw[swapped]
    h$outcome_other_allele_harmonised[swapped] <- h$outcome_effect_allele_raw[swapped]
    h$outcome_beta_harmonised[swapped] <- -h$outcome_beta_raw[swapped]
    h$outcome_eaf_harmonised[swapped] <- ifelse(is.na(h$outcome_eaf_raw[swapped]), NA_real_, 1 - h$outcome_eaf_raw[swapped])
    h$beta_flipped[swapped] <- TRUE
    h$eaf_flipped[swapped] <- TRUE
  }
  if (any(strand)) h$strand_flipped[strand] <- TRUE
  if (any(incompatible)) {
    h$record_excluded[incompatible] <- TRUE
    h$exclusion_reason[incompatible] <- h$orientation_class[incompatible]
  }
  h$final_valid_instrument <- h$source_match_unique & h$valid_exposure_alleles & h$valid_outcome_alleles &
    h$orientation_class %in% c("exact_match", "swapped_match", "strand_exact_match", "strand_swapped_match") &
    !h$record_excluded

  keep <- c("analysis_set", "rsid", "exposure_effect_allele", "exposure_other_allele", "exposure_beta", "exposure_se",
            "exposure_pval", "exposure_eaf", "outcome_effect_allele_raw", "outcome_other_allele_raw",
            "outcome_beta_raw", "outcome_se_raw", "outcome_pval_raw", "outcome_eaf_raw",
            "outcome_effect_allele_harmonised", "outcome_other_allele_harmonised",
            "outcome_beta_harmonised", "outcome_se_harmonised", "outcome_pval_harmonised",
            "outcome_eaf_harmonised", "orientation_class", "palindromic_snp", "beta_flipped",
            "eaf_flipped", "strand_flipped", "record_excluded", "exclusion_reason", "final_valid_instrument")
  h <- h[, keep]
  inc <- h[h$analysis_set == "APOE included", , drop = FALSE]
  exc <- h[h$analysis_set == "APOE excluded", , drop = FALSE]

  hard_check_failures <- character()
  add_fail <- function(x) hard_check_failures <<- unique(c(hard_check_failures, x))
  if (!all(pre$eligible_for_formal_harmonisation == h$final_valid_instrument)) add_fail("all_preflight_eligible_records_accounted_for_failed")
  if (!identical(pre$raw_orientation_class, h$orientation_class)) add_fail("orientation_matches_frozen_preflight_failed")
  exact <- h$orientation_class == "exact_match"
  if (any(exact)) {
    if (!all(num_equal(h$outcome_beta_harmonised[exact], h$outcome_beta_raw[exact]))) add_fail("exact_match_beta_invariance_failed")
    if (!all(num_equal(h$outcome_se_harmonised[exact], h$outcome_se_raw[exact]))) add_fail("exact_match_se_invariance_failed")
    if (!all(num_equal(h$outcome_pval_harmonised[exact], h$outcome_pval_raw[exact]))) add_fail("exact_match_pval_invariance_failed")
    if (!all(num_equal(h$outcome_eaf_harmonised[exact], h$outcome_eaf_raw[exact]))) add_fail("exact_match_eaf_invariance_failed")
    if (!all(h$outcome_effect_allele_harmonised[exact] == h$exposure_effect_allele[exact] &
             h$outcome_other_allele_harmonised[exact] == h$exposure_other_allele[exact])) add_fail("exact_match_allele_invariance_failed")
  }
  if (!identical(h$exposure_beta, pre$exposure_beta) || !identical(h$exposure_se, pre$exposure_se) ||
      !identical(h$exposure_pval, pre$exposure_pval) || !identical(h$exposure_eaf, pre$exposure_eaf)) add_fail("exposure_fields_unchanged_failed")

  write_pq(con, inc, paste0(out[["inc_pq"]], ".partial"))
  write.table(inc, paste0(out[["inc_tsv"]], ".partial"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  write_pq(con, exc, paste0(out[["exc_pq"]], ".partial"))
  write.table(exc, paste0(out[["exc_tsv"]], ".partial"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  inc_consistency <- check_pq_tsv(con, paste0(out[["inc_pq"]], ".partial"), paste0(out[["inc_tsv"]], ".partial"))
  exc_consistency <- check_pq_tsv(con, paste0(out[["exc_pq"]], ".partial"), paste0(out[["exc_tsv"]], ".partial"))
  if (!inc_consistency) add_fail("included_parquet_tsv_consistency_failed")
  if (!exc_consistency) add_fail("excluded_parquet_tsv_consistency_failed")

  included_final <- sum(inc$final_valid_instrument)
  excluded_final <- sum(exc$final_valid_instrument)
  included_lt3 <- included_final < 3L
  excluded_lt3 <- excluded_final < 3L
  overall_trigger <- included_lt3
  relaxed_status <- if (overall_trigger) "triggered_not_started" else "not_triggered"
  trigger_reason <- if (overall_trigger) "final_valid_independent_primary_instruments_after_harmonisation_lt_3" else ""
  approved_relaxed <- overall_trigger
  counts <- data.frame(analysis_set = c("APOE included", "APOE excluded"),
                       input_count = c(nrow(inc), nrow(exc)),
                       final_valid_instrument_count = c(included_final, excluded_final),
                       lt3 = c(included_lt3, excluded_lt3),
                       stringsAsFactors = FALSE)
  write.csv(counts, paste0(out[["counts"]], ".partial"), row.names = FALSE)

  hard_checks <- list(
    preflight_freeze_gate = identical(freeze$freeze_status, "passed"),
    included_input_present = nrow(inc) == 1L,
    excluded_input_present = nrow(exc) == 1L,
    all_preflight_eligible_records_accounted_for = all(pre$eligible_for_formal_harmonisation == h$final_valid_instrument),
    orientation_matches_frozen_preflight = identical(pre$raw_orientation_class, h$orientation_class),
    included_harmonisation_completed = nrow(inc) == 1L,
    excluded_harmonisation_completed = nrow(exc) == 1L,
    exact_match_invariance = length(grep("exact_match", hard_check_failures)) == 0L,
    outcome_se_unchanged = all(num_equal(h$outcome_se_harmonised, h$outcome_se_raw)),
    outcome_pval_unchanged = all(num_equal(h$outcome_pval_harmonised, h$outcome_pval_raw)),
    exposure_fields_unchanged = !("exposure_fields_unchanged_failed" %in% hard_check_failures),
    final_valid_counts_computed = !is.na(included_final) && !is.na(excluded_final),
    included_parquet_tsv_consistency = inc_consistency,
    excluded_parquet_tsv_consistency = exc_consistency,
    protocol_trigger_evaluated_after_formal_harmonisation = TRUE,
    no_relaxed_threshold_analysis_started = TRUE,
    no_reverse_mr = TRUE,
    no_proxy = TRUE,
    no_liftover = TRUE,
    no_source_scan = TRUE,
    no_clumping = TRUE
  )
  harmonisation_status <- if (length(hard_check_failures) == 0L && all(unlist(hard_checks))) "passed" else "failed"
  renv_after <- hash_file(renv_lock)
  renv_lock_unchanged <- identical(renv_before, renv_after)
  if (!renv_lock_unchanged) {
    hard_check_failures <- unique(c(hard_check_failures, "renv_lock_changed"))
    harmonisation_status <- "failed"
  }
  qc <- list(
    harmonisation_version = "v4",
    supersedes_failed_harmonisation_versions = c("v1", "v2", "v3"),
    analysis_direction = "delirium_to_Hb",
    analysis_role = "secondary_reverse_primary",
    source_instrument_version = "v4",
    source_outcome_extraction_version = "v1",
    source_preflight_version = "v1",
    included_input_count = nrow(inc),
    excluded_input_count = nrow(exc),
    included_final_valid_instrument_count = included_final,
    excluded_final_valid_instrument_count = excluded_final,
    exact_match_count = sum(h$orientation_class == "exact_match"),
    swapped_match_count = sum(h$orientation_class == "swapped_match"),
    strand_exact_match_count = sum(h$orientation_class == "strand_exact_match"),
    strand_swapped_match_count = sum(h$orientation_class == "strand_swapped_match"),
    palindromic_count = sum(h$palindromic_snp),
    incompatible_count = sum(h$orientation_class == "incompatible"),
    invalid_count = sum(h$orientation_class == "invalid"),
    beta_flip_count = sum(h$beta_flipped),
    eaf_flip_count = sum(h$eaf_flipped),
    strand_flip_count = sum(h$strand_flipped),
    excluded_record_count = sum(h$record_excluded),
    harmonisation_rows = records(h),
    included_parquet_tsv_consistency = inc_consistency,
    excluded_parquet_tsv_consistency = exc_consistency,
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after,
    renv_lock_unchanged = renv_lock_unchanged,
    informational_findings = c("renv out-of-sync warning may be emitted by project activation; renv.lock was not modified"),
    primary_relaxed_threshold_trigger_threshold = 3L,
    included_lt3 = included_lt3,
    excluded_lt3 = excluded_lt3,
    overall_relaxed_threshold_trigger = overall_trigger,
    relaxed_threshold = 5e-6,
    relaxed_threshold_status = relaxed_status,
    relaxed_threshold_trigger_reason = trigger_reason,
    approved_for_reverse_relaxed_threshold_branch = approved_relaxed,
    harmonisation_status = harmonisation_status,
    hard_checks = hard_checks,
    hard_check_failures = hard_check_failures
  )
  jsonlite::write_json(qc, paste0(out[["json"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  if (!identical(harmonisation_status, "passed")) stop("Formal harmonisation V4 did not pass; partial outputs retained.", call. = FALSE)
  for (p in out[c("inc_pq", "inc_tsv", "exc_pq", "exc_tsv", "counts", "json")]) {
    stop_if(file.exists(p), paste("Output appeared during run:", p))
    stop_if(!file.rename(paste0(p, ".partial"), p), paste("Atomic rename failed:", p))
  }
  log_line(sprintf("SUCCESS harmonisation_status=%s included_final=%s excluded_final=%s relaxed_threshold_status=%s hard_check_failures=%s", harmonisation_status, included_final, excluded_final, relaxed_status, paste(hard_check_failures, collapse = ";")))
}

tryCatch(main(), error = function(e) {
  log_line(paste0("TERMINATED_PRIMARY: ", conditionMessage(e)))
  quit(status = 1L)
})

