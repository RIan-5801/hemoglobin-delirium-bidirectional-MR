#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/31b_chen_reverse_formal_harmonisation_v1_readback_recovery.R [--project-root <path>]", call. = FALSE)
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
records <- function(x) {
  if (!is.data.frame(x) || nrow(x) == 0L) return(list())
  lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
}
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
write_csv <- function(x, path) atomic_write(path, function(p) {
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.csv(x, p, row.names = FALSE, na = "")
})
write_text <- function(lines, path) atomic_write(path, function(p) writeLines(lines, p, useBytes = TRUE))
is_missing_token <- function(x) is.na(x) | (!is.na(x) & as.character(x) == "")

paths <- list(
  script = rel("R", "31b_chen_reverse_formal_harmonisation_v1_readback_recovery.R"),
  formal_script = rel("R", "31_chen_reverse_formal_harmonisation_v1.R"),
  formal_qc = rel("results", "qc", "chen_reverse_formal_harmonisation_v1.json"),
  formal_decision = rel("docs", "decisions", "114_chen_reverse_formal_harmonisation_v1_v1.1.md"),
  formal_log = rel("results", "logs", "chen_reverse_formal_harmonisation_v1.log"),
  counts = rel("results", "qc", "chen_reverse_formal_harmonisation_counts_v1.csv"),
  transform_audit = rel("results", "qc", "chen_reverse_formal_harmonisation_transform_audit_v1.csv"),
  excluded = rel("results", "qc", "chen_reverse_formal_harmonisation_excluded_snps_v1.tsv"),
  master_parquet = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_master_v1.parquet"),
  master_tsv = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_master_v1.tsv"),
  strict_inc_parquet = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_strict_apoe_included_v1.parquet"),
  strict_inc_tsv = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_strict_apoe_included_v1.tsv"),
  strict_exc_parquet = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_strict_apoe_excluded_v1.parquet"),
  strict_exc_tsv = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_strict_apoe_excluded_v1.tsv"),
  relaxed_inc_parquet = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_relaxed_apoe_included_v1.parquet"),
  relaxed_inc_tsv = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_relaxed_apoe_included_v1.tsv"),
  relaxed_exc_parquet = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_relaxed_apoe_excluded_v1.parquet"),
  relaxed_exc_tsv = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_relaxed_apoe_excluded_v1.tsv"),
  recovery_qc = rel("results", "qc", "chen_reverse_formal_harmonisation_v1_readback_recovery.json"),
  mismatch_csv = rel("results", "qc", "chen_reverse_formal_harmonisation_v1_readback_recovery_mismatch_audit.csv"),
  log = rel("results", "logs", "chen_reverse_formal_harmonisation_v1_readback_recovery.log"),
  decision = rel("docs", "decisions", "115_chen_reverse_formal_harmonisation_v1_readback_recovery_v1.1.md"),
  renv_lock = rel("renv.lock")
)

outputs <- unlist(paths[c("recovery_qc", "mismatch_csv", "log", "decision")])
occupied <- outputs[file.exists(outputs) | file.exists(paste0(outputs, ".partial"))]
stop_if(length(occupied) > 0L, paste("Output or partial exists:", paste(occupied, collapse = "; ")))
inputs <- unlist(paths[c("formal_script", "formal_qc", "formal_decision", "formal_log", "counts", "transform_audit", "excluded", "master_parquet", "master_tsv", "strict_inc_parquet", "strict_inc_tsv", "strict_exc_parquet", "strict_exc_tsv", "relaxed_inc_parquet", "relaxed_inc_tsv", "relaxed_exc_parquet", "relaxed_exc_tsv", "renv_lock")])
missing_inputs <- inputs[!file.exists(inputs)]
stop_if(length(missing_inputs) > 0L, paste("Required input missing:", paste(missing_inputs, collapse = "; ")))

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 115L), paste0("Expected next decision 115, found ", next_decision, "; no outputs written."))

dir.create(dirname(paths$log), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = paths$log, append = TRUE)
log_line("stage=chen_reverse_formal_harmonisation_v1_readback_recovery_start")

formal_qc <- read_json(paths$formal_qc)
stop_if(!identical(formal_qc$harmonisation_status, "failed"), "Decision 115 recovery expects Decision 114 QC to be failed.")
stop_if(!setequal(formal_qc$hard_check_failures, c("master_parquet_tsv_consistency", "all_four_branch_parquet_tsv_consistency")), "Decision 114 failure pattern is not readback-only.")

nullable_character_fields <- c(
  "strict_included_projected_exclusion_reason",
  "strict_excluded_projected_exclusion_reason",
  "relaxed_included_projected_exclusion_reason",
  "relaxed_excluded_projected_exclusion_reason",
  "outcome_effect_allele_harmonised",
  "outcome_other_allele_harmonised",
  "exclusion_reason"
)

con_db <- DBI::dbConnect(duckdb::duckdb(), read_only = FALSE)
on.exit(try(DBI::dbDisconnect(con_db, shutdown = TRUE), silent = TRUE), add = TRUE)
read_parquet <- function(path) DBI::dbGetQuery(con_db, sprintf("SELECT * FROM read_parquet(%s)", DBI::dbQuoteString(con_db, norm(path))))
validate_pair <- function(label, parquet, tsv, key = "rsid") {
  p <- read_parquet(parquet)
  t <- utils::read.delim(tsv, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
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
      ok <- if (k %in% nullable_character_fields) {
        (is_missing_token(a[[i]]) && is_missing_token(b[[i]])) || (!is.na(a[[i]]) && !is.na(b[[i]]) && a[[i]] == b[[i]])
      } else {
        !is.na(a[[i]]) && !is.na(b[[i]]) && a[[i]] == b[[i]]
      }
      if (!ok) {
        mismatches <- rbind(mismatches, data.frame(pair = label, column = k, row = i, parquet = ifelse(is.na(a[[i]]), "<NA>", a[[i]]), tsv = ifelse(is.na(b[[i]]), "<NA>", b[[i]]), mismatch_type = "character_value_mismatch", stringsAsFactors = FALSE))
      }
    }
  }
  numeric_stats <- lapply(num_cols, function(k) {
    a <- as.numeric(p[[k]])
    b <- suppressWarnings(as.numeric(t[[k]]))
    both_missing <- is.na(a) & is_missing_token(as.character(t[[k]]))
    both_num <- is.finite(a) & is.finite(b)
    abs_diff <- rep(NA_real_, length(a))
    rel_diff <- rep(NA_real_, length(a))
    abs_diff[both_num] <- abs(a[both_num] - b[both_num])
    rel_diff[both_num] <- abs_diff[both_num] / pmax(abs(a[both_num]), abs(b[both_num]), .Machine$double.xmin)
    data.frame(column = k, all_cells_ok = all(both_missing | both_num), max_abs = if (any(both_num)) max(abs_diff[both_num], na.rm = TRUE) else 0, max_rel = if (any(both_num)) max(rel_diff[both_num], na.rm = TRUE) else 0, stringsAsFactors = FALSE)
  })
  numeric_stats <- if (length(numeric_stats)) do.call(rbind, numeric_stats) else data.frame(column = character(), all_cells_ok = logical(), max_abs = numeric(), max_rel = numeric())
  numeric_ok <- all(numeric_stats$all_cells_ok) && all(numeric_stats$max_abs <= 1e-12 | numeric_stats$max_rel <= 1e-12)
  summary <- data.frame(pair = label, row_count = nrow(p), same_cols = same_cols, same_n = same_n, same_order = same_order, character_ok_with_field_specific_missing_equivalence = nrow(mismatches) == 0L, numeric_ok = numeric_ok, pair_ok = same_cols && same_n && same_order && nrow(mismatches) == 0L && numeric_ok, stringsAsFactors = FALSE)
  list(summary = summary, mismatches = mismatches, numeric_stats = numeric_stats)
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

master <- read_parquet(paths$master_parquet)
final_rows <- master[master$final_valid, , drop = FALSE]
hard_checks <- list(
  decision_114_failure_preserved = identical(formal_qc$harmonisation_status, "failed"),
  failure_pattern_limited_to_readback_consistency = setequal(formal_qc$hard_check_failures, c("master_parquet_tsv_consistency", "all_four_branch_parquet_tsv_consistency")),
  all_non_readback_hard_checks_passed_in_decision_114 = all(unlist(formal_qc$hard_checks[setdiff(names(formal_qc$hard_checks), formal_qc$hard_check_failures)])),
  parquet_tsv_readback_ok_with_field_specific_missing_equivalence = all(pair_summary$pair_ok),
  no_true_value_mismatches_after_missing_equivalence = nrow(mismatch_audit) == 0L,
  final_valid_total_reverified = nrow(final_rows) == 11L,
  no_palindromic_final_valid = sum(master$final_valid & master$palindromic_recomputed) == 0L,
  final_alleles_aligned = all(final_rows$outcome_effect_allele_harmonised == final_rows$exposure_effect_allele) && all(final_rows$outcome_other_allele_harmonised == final_rows$exposure_other_allele),
  no_mr = isFALSE(formal_qc$mr_run),
  no_steiger = isFALSE(formal_qc$steiger_run),
  renv_lock_unchanged = identical(hash_file(paths$renv_lock), "253471c02e1e47a40d0f68b296d4ae2b1df471d757ac86328b3e974018d039f3")
)
hard_check_failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
recovery_status <- if (length(hard_check_failures) == 0L) "passed" else "failed"

qc <- list(
  readback_recovery_version = "v1",
  decision = 115,
  date = format(Sys.Date()),
  recovery_status = recovery_status,
  formal_harmonisation_decision = 114,
  formal_qc = relpath(paths$formal_qc),
  approved_for_chen_reverse_mr_input_freeze = identical(recovery_status, "passed"),
  recovery_scope = "Parquet/TSV readback validation only; TSV empty strings are treated as missing-value encodings only for explicitly listed nullable character fields.",
  nullable_character_fields_for_tsv_blank_normalization = nullable_character_fields,
  formal_scientific_outputs_changed = FALSE,
  harmonisation_rerun = FALSE,
  mr_run = FALSE,
  steiger_run = FALSE,
  pair_summary = records(pair_summary),
  mismatch_audit = relpath(paths$mismatch_csv),
  hard_checks = hard_checks,
  hard_check_failures = hard_check_failures,
  informational_findings = list(
    decision_114_status_preserved = formal_qc$harmonisation_status,
    decision_114_failures_preserved = formal_qc$hard_check_failures,
    prior_failure_explanation = "Decision 114 omitted four branch-level projected exclusion reason fields from nullable character readback normalization; no true value mismatch was observed.",
    final_valid_total = nrow(final_rows)
  )
)
write_json(qc, paths$recovery_qc)

decision_lines <- c(
  "# Decision 115: Chen Reverse Formal Harmonisation V1 Readback Recovery",
  "",
  paste0("Date: ", format(Sys.Date())),
  "",
  "## Status",
  paste0("recovery_status: `", recovery_status, "`"),
  paste0("approved_for_chen_reverse_mr_input_freeze: `", identical(recovery_status, "passed"), "`"),
  "",
  "## Decision",
  "Classify the Decision 114 Parquet/TSV consistency failure as a technical readback encoding issue limited to explicitly nullable character fields.",
  "",
  "This recovery does not rerun harmonisation, alter alleles, change final-valid membership, perform MR, or run Steiger.",
  "",
  "## Evidence",
  paste0("- Decision 114 status preserved as: `", formal_qc$harmonisation_status, "`."),
  paste0("- Decision 114 failed checks preserved as: `", paste(formal_qc$hard_check_failures, collapse = ";"), "`."),
  paste0("- Parquet/TSV pair checks passed with field-specific missing-value equivalence: `", all(pair_summary$pair_ok), "`."),
  paste0("- True mismatches after missing-value equivalence: `", nrow(mismatch_audit), "`."),
  paste0("- Final-valid total reverified: `", nrow(final_rows), "`."),
  paste0("- Hard-check failures: `", if (length(hard_check_failures) == 0L) "none" else paste(hard_check_failures, collapse = ";"), "`."),
  "",
  "## Outputs",
  paste0("- `", relpath(paths$recovery_qc), "`"),
  paste0("- `", relpath(paths$mismatch_csv), "`"),
  paste0("- `", relpath(paths$log), "`")
)
write_text(decision_lines, paths$decision)

log_line("recovery_status=", recovery_status, "; hard_check_failures=", paste(hard_check_failures, collapse = ","))
cat("recovery_status=", recovery_status, "\n", sep = "")
cat("approved_for_chen_reverse_mr_input_freeze=", identical(recovery_status, "passed"), "\n", sep = "")
cat("hard_check_failures=", paste(hard_check_failures, collapse = ","), "\n", sep = "")
cat("true_mismatches_after_missing_equivalence=", nrow(mismatch_audit), "\n", sep = "")
