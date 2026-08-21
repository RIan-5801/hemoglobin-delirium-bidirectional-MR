#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/29s_chen_reverse_outcome_extraction_v1_readback_recovery.R [--project-root <path>]", call. = FALSE)
}
setwd(root)

for (pkg in c("DBI", "duckdb", "jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

rel <- function(...) file.path(root, ...)
norm <- function(x) gsub("\\\\", "/", x)
relpath <- function(x) norm(sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", norm(normalizePath(x, winslash = "/", mustWork = FALSE))))
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
atomic_write <- function(path, writer) {
  partial <- paste0(path, ".partial")
  stop_if(file.exists(path) || file.exists(partial), paste("Output occupied:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writer(partial)
  stop_if(!file.exists(partial), paste("Writer did not create partial:", partial))
  stop_if(!file.rename(partial, path), paste("Atomic rename failed:", path))
}
write_json <- function(x, path) atomic_write(path, function(p) {
  jsonlite::write_json(x, p, pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)
})
write_text <- function(lines, path) atomic_write(path, function(p) writeLines(lines, p, useBytes = TRUE))
write_csv <- function(x, path) atomic_write(path, function(p) {
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.csv(x, p, row.names = FALSE, na = "")
})
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))

paths <- list(
  script = rel("R", "29s_chen_reverse_outcome_extraction_v1_readback_recovery.R"),
  failed_script = rel("R", "29_chen_reverse_outcome_extraction_v1.R"),
  extraction_script = rel("R", "29r_chen_reverse_outcome_extraction_v1_technical_recovery.R"),
  extraction_qc = rel("results", "qc", "chen_reverse_outcome_extraction_v1.json"),
  extraction_decision = rel("docs", "decisions", "110_chen_reverse_outcome_extraction_v1_v1.1.md"),
  initial_failed_log = rel("results", "logs", "chen_reverse_outcome_extraction_v1.log"),
  extraction_recovery_log = rel("results", "logs", "chen_reverse_outcome_extraction_v1_technical_recovery.log"),
  union_targets = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_union_targets_v1.tsv"),
  master_parquet = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_master_v1.parquet"),
  master_tsv = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_master_v1.tsv"),
  strict_inc_parquet = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_strict_apoe_included_v1.parquet"),
  strict_inc_tsv = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_strict_apoe_included_v1.tsv"),
  strict_exc_parquet = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_strict_apoe_excluded_v1.parquet"),
  strict_exc_tsv = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_strict_apoe_excluded_v1.tsv"),
  relaxed_inc_parquet = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_relaxed_apoe_included_v1.parquet"),
  relaxed_inc_tsv = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_relaxed_apoe_included_v1.tsv"),
  relaxed_exc_parquet = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_relaxed_apoe_excluded_v1.parquet"),
  relaxed_exc_tsv = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_relaxed_apoe_excluded_v1.tsv"),
  match_audit = rel("results", "qc", "chen_reverse_outcome_match_audit_v1.csv"),
  missing = rel("results", "qc", "chen_reverse_outcome_missing_v1.tsv"),
  qc_json = rel("results", "qc", "chen_reverse_outcome_extraction_v1_readback_recovery.json"),
  mismatch_csv = rel("results", "qc", "chen_reverse_outcome_extraction_v1_readback_recovery_mismatch_audit.csv"),
  log = rel("results", "logs", "chen_reverse_outcome_extraction_v1_readback_recovery.log"),
  decision = rel("docs", "decisions", "111_chen_reverse_outcome_extraction_v1_readback_recovery_v1.1.md"),
  renv_lock = rel("renv.lock")
)

outputs <- unlist(paths[c("qc_json", "mismatch_csv", "log", "decision")])
occupied <- outputs[file.exists(outputs) | file.exists(paste0(outputs, ".partial"))]
stop_if(length(occupied) > 0L, paste("Output or partial exists:", paste(occupied, collapse = "; ")))

required_inputs <- unlist(paths[c(
  "failed_script", "extraction_script", "extraction_qc", "extraction_decision", "initial_failed_log",
  "extraction_recovery_log", "union_targets", "master_parquet", "master_tsv", "strict_inc_parquet",
  "strict_inc_tsv", "strict_exc_parquet", "strict_exc_tsv", "relaxed_inc_parquet", "relaxed_inc_tsv",
  "relaxed_exc_parquet", "relaxed_exc_tsv", "match_audit", "missing", "renv_lock"
)])
missing_inputs <- required_inputs[!file.exists(required_inputs)]
stop_if(length(missing_inputs) > 0L, paste("Required input missing:", paste(missing_inputs, collapse = "; ")))

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 111L), paste0("Expected next decision 111, found ", next_decision, "; no outputs written."))

dir.create(dirname(paths$log), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = paths$log, append = TRUE)
log_line("stage=chen_reverse_outcome_extraction_v1_readback_recovery_start")

extraction_qc <- read_json(paths$extraction_qc)
stop_if(!identical(extraction_qc$outcome_extraction_status, "failed"), "Decision 111 recovery expects Decision 110 extraction QC to be failed.")
stop_if(!setequal(extraction_qc$hard_check_failures, c("master_parquet_tsv_consistency", "branch_parquet_tsv_consistency")), "Decision 110 failure pattern is not the approved readback-only pattern.")
stop_if(!identical(extraction_qc$union_unique_exact_match_count, 12L), "Decision 110 exact match count is not 12.")
stop_if(!identical(extraction_qc$union_missing_count, 0L), "Decision 110 missing count is not zero.")
stop_if(!identical(extraction_qc$union_multiple_exact_match_count, 0L), "Decision 110 multiple count is not zero.")
stop_if(!identical(extraction_qc$union_marker_effect_allele_incompatible_count, 0L), "Decision 110 incompatible count is not zero.")
stop_if(!identical(extraction_qc$source_sha_before, extraction_qc$source_sha_after), "Decision 110 source SHA changed.")

con_db <- DBI::dbConnect(duckdb::duckdb(), read_only = FALSE)
on.exit(try(DBI::dbDisconnect(con_db, shutdown = TRUE), silent = TRUE), add = TRUE)
read_parquet <- function(path) DBI::dbGetQuery(con_db, sprintf("SELECT * FROM read_parquet(%s)", DBI::dbQuoteString(con_db, norm(path))))
read_tsv_chr <- function(path) utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
is_missing_token <- function(x) is.na(x) | identical(x, "") | (!is.na(x) & x == "")
validate_pair <- function(label, parquet_path, tsv_path, key = "target_rsid") {
  p <- read_parquet(parquet_path)
  t <- read_tsv_chr(tsv_path)
  same_cols <- identical(names(p), names(t))
  same_n <- identical(nrow(p), nrow(t))
  same_order <- same_n && identical(as.character(p[[key]]), as.character(t[[key]]))
  char_cols <- names(p)[!vapply(p, is.numeric, logical(1))]
  num_cols <- names(p)[vapply(p, is.numeric, logical(1))]
  mismatches <- data.frame(pair = character(), column = character(), row = integer(), parquet = character(), tsv = character(), mismatch_type = character(), stringsAsFactors = FALSE)
  for (k in char_cols) {
    a <- as.character(p[[k]])
    b <- as.character(t[[k]])
    for (i in seq_along(a)) {
      missing_equiv <- isTRUE(is_missing_token(a[[i]]) && is_missing_token(b[[i]]))
      exact <- !is.na(a[[i]]) && !is.na(b[[i]]) && identical(a[[i]], b[[i]])
      if (!missing_equiv && !exact) {
        mismatches <- rbind(mismatches, data.frame(
          pair = label, column = k, row = i,
          parquet = ifelse(is.na(a[[i]]), "<NA>", a[[i]]),
          tsv = ifelse(is.na(b[[i]]), "<NA>", b[[i]]),
          mismatch_type = "character_value_mismatch",
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  numeric_stats <- lapply(num_cols, function(k) {
    a <- as.numeric(p[[k]])
    b <- suppressWarnings(as.numeric(t[[k]]))
    both_missing <- (is.na(a) & is_missing_token(as.character(t[[k]])))
    both_num <- is.finite(a) & is.finite(b)
    abs_diff <- rep(NA_real_, length(a))
    rel_diff <- rep(NA_real_, length(a))
    abs_diff[both_num] <- abs(a[both_num] - b[both_num])
    rel_diff[both_num] <- abs_diff[both_num] / pmax(abs(a[both_num]), abs(b[both_num]), .Machine$double.xmin)
    data.frame(
      column = k,
      all_cells_ok = all(both_missing | both_num),
      max_abs = if (any(both_num)) max(abs_diff[both_num], na.rm = TRUE) else 0,
      max_rel = if (any(both_num)) max(rel_diff[both_num], na.rm = TRUE) else 0,
      stringsAsFactors = FALSE
    )
  })
  numeric_stats <- if (length(numeric_stats)) do.call(rbind, numeric_stats) else data.frame(column = character(), all_cells_ok = logical(), max_abs = numeric(), max_rel = numeric())
  numeric_ok <- all(numeric_stats$all_cells_ok) && all(numeric_stats$max_abs <= 1e-12 | numeric_stats$max_rel <= 1e-12)
  list(
    summary = data.frame(
      pair = label,
      parquet = relpath(parquet_path),
      tsv = relpath(tsv_path),
      row_count = nrow(p),
      same_cols = same_cols,
      same_n = same_n,
      same_order = same_order,
      character_ok_with_missing_equivalence = nrow(mismatches) == 0L,
      numeric_ok = numeric_ok,
      pair_ok = same_cols && same_n && same_order && nrow(mismatches) == 0L && numeric_ok,
      stringsAsFactors = FALSE
    ),
    mismatches = mismatches,
    numeric_stats = numeric_stats
  )
}

pairs <- list(
  master = validate_pair("master", paths$master_parquet, paths$master_tsv),
  strict_included = validate_pair("strict_included", paths$strict_inc_parquet, paths$strict_inc_tsv),
  strict_excluded = validate_pair("strict_excluded", paths$strict_exc_parquet, paths$strict_exc_tsv),
  relaxed_included = validate_pair("relaxed_included", paths$relaxed_inc_parquet, paths$relaxed_inc_tsv),
  relaxed_excluded = validate_pair("relaxed_excluded", paths$relaxed_exc_parquet, paths$relaxed_exc_tsv)
)
pair_summary <- do.call(rbind, lapply(pairs, `[[`, "summary"))
mismatch_audit <- do.call(rbind, lapply(pairs, `[[`, "mismatches"))
if (is.null(mismatch_audit) || nrow(mismatch_audit) == 0L) {
  mismatch_audit <- data.frame(pair = character(), column = character(), row = integer(), parquet = character(), tsv = character(), mismatch_type = character(), stringsAsFactors = FALSE)
}
write_csv(mismatch_audit, paths$mismatch_csv)

master <- read_tsv_chr(paths$master_tsv)
match_audit <- utils::read.csv(paths$match_audit, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
empty_reason_all_exact <- all(master$match_status == "unique_exact_match") && all(master$reason == "" | is.na(master$reason))
match_audit_counts_ok <- all(as.integer(match_audit$raw_match_count) == 1L) && all(match_audit$match_status == "unique_exact_match")

hard_checks <- list(
  decision_110_failure_preserved = file.exists(paths$extraction_decision) && identical(extraction_qc$outcome_extraction_status, "failed"),
  failure_pattern_limited_to_readback_consistency = setequal(extraction_qc$hard_check_failures, c("master_parquet_tsv_consistency", "branch_parquet_tsv_consistency")),
  source_scan_completed_in_decision_110 = isTRUE(extraction_qc$source_scan_completed),
  source_sha_unchanged_in_decision_110 = identical(extraction_qc$source_sha_before, extraction_qc$source_sha_after),
  all_targets_exact_in_decision_110 = identical(extraction_qc$union_unique_exact_match_count, 12L) && identical(extraction_qc$union_missing_count, 0L) && identical(extraction_qc$union_multiple_exact_match_count, 0L) && identical(extraction_qc$union_marker_effect_allele_incompatible_count, 0L),
  exact_matching_method_preserved = identical(extraction_qc$matching_method, "exact_BCX2_marker_id_string_identity"),
  no_proxy = isFALSE(extraction_qc$proxy_used),
  no_liftover = isFALSE(extraction_qc$liftover_used),
  no_nearest_variant = isFALSE(extraction_qc$nearest_variant_used),
  no_strand_complement_identity_rescue = isFALSE(extraction_qc$strand_complement_identity_rescue_used),
  no_harmonisation = isFALSE(extraction_qc$harmonisation_performed),
  no_mr = isFALSE(extraction_qc$mr_run),
  no_steiger = isFALSE(extraction_qc$steiger_run),
  parquet_tsv_readback_ok_with_missing_equivalence = all(pair_summary$pair_ok),
  no_true_value_mismatches_after_missing_equivalence = nrow(mismatch_audit) == 0L,
  empty_reason_column_explains_prior_char_failure = empty_reason_all_exact,
  match_audit_counts_ok = match_audit_counts_ok,
  branch_counts_preserved = identical(vapply(extraction_qc$branch_counts, function(x) x$exact, integer(1)), c(1L, 1L, 11L, 10L)),
  renv_lock_unchanged = identical(hash_file(paths$renv_lock), "253471c02e1e47a40d0f68b296d4ae2b1df471d757ac86328b3e974018d039f3")
)
hard_check_failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
recovery_status <- if (length(hard_check_failures) == 0L) "passed" else "failed"

qc <- list(
  readback_recovery_version = "v1",
  decision = 111,
  date = format(Sys.Date()),
  recovery_status = recovery_status,
  extraction_decision = 110,
  extraction_qc = relpath(paths$extraction_qc),
  approved_for_freeze = identical(recovery_status, "passed"),
  recovery_scope = "Parquet/TSV readback validation only; TSV empty strings are treated as missing-value encodings equivalent to Parquet NA.",
  extraction_scientific_outputs_changed = FALSE,
  raw_source_scan_performed = FALSE,
  harmonisation_performed = FALSE,
  mr_run = FALSE,
  steiger_run = FALSE,
  pair_summary = records(pair_summary),
  mismatch_audit = relpath(paths$mismatch_csv),
  hard_checks = hard_checks,
  hard_check_failures = hard_check_failures,
  informational_findings = list(
    prior_failed_checks = extraction_qc$hard_check_failures,
    prior_failure_explanation = "Parquet readback represented reason as NA while TSV represented reason as empty string because write.table used na=''. All exact-match rows legitimately have no failure reason.",
    source_rows_scanned_in_decision_110 = extraction_qc$source_rows_scanned,
    union_unique_exact_match_count = extraction_qc$union_unique_exact_match_count,
    union_missing_count = extraction_qc$union_missing_count,
    union_multiple_exact_match_count = extraction_qc$union_multiple_exact_match_count,
    union_marker_effect_allele_incompatible_count = extraction_qc$union_marker_effect_allele_incompatible_count
  )
)
write_json(qc, paths$qc_json)

decision_lines <- c(
  "# Decision 111: Chen Reverse Outcome Extraction V1 Readback Recovery",
  "",
  paste0("Date: ", format(Sys.Date())),
  "",
  "## Status",
  paste0("recovery_status: `", recovery_status, "`"),
  paste0("approved_for_freeze: `", identical(recovery_status, "passed"), "`"),
  "",
  "## Decision",
  "Classify the Decision 110 Parquet/TSV consistency failure as a technical readback encoding issue limited to missing values in the `reason` column.",
  "",
  "This recovery does not alter target variants, matching rules, source data, effect alleles, outcome estimates, branch membership, or any scientific parameter.",
  "",
  "## Evidence",
  paste0("- Decision 110 status preserved as: `", extraction_qc$outcome_extraction_status, "`."),
  paste0("- Decision 110 failed checks preserved as: `", paste(extraction_qc$hard_check_failures, collapse = ";"), "`."),
  paste0("- Source rows scanned in Decision 110: `", extraction_qc$source_rows_scanned, "`."),
  paste0("- Union exact / missing / multiple / incompatible: `", extraction_qc$union_unique_exact_match_count, " / ", extraction_qc$union_missing_count, " / ", extraction_qc$union_multiple_exact_match_count, " / ", extraction_qc$union_marker_effect_allele_incompatible_count, "`."),
  paste0("- Parquet/TSV pair checks passed with missing-value equivalence: `", all(pair_summary$pair_ok), "`."),
  paste0("- True mismatches after missing-value equivalence: `", nrow(mismatch_audit), "`."),
  "",
  "## Outputs",
  paste0("- `", relpath(paths$qc_json), "`"),
  paste0("- `", relpath(paths$mismatch_csv), "`"),
  paste0("- `", relpath(paths$log), "`")
)
write_text(decision_lines, paths$decision)

log_line("recovery_status=", recovery_status, "; hard_check_failures=", paste(hard_check_failures, collapse = ","))
cat("recovery_status=", recovery_status, "\n", sep = "")
cat("approved_for_freeze=", identical(recovery_status, "passed"), "\n", sep = "")
cat("hard_check_failures=", paste(hard_check_failures, collapse = ","), "\n", sep = "")
cat("true_mismatches_after_missing_equivalence=", nrow(mismatch_audit), "\n", sep = "")
