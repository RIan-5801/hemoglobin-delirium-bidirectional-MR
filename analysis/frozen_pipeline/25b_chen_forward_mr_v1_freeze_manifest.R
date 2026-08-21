#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/25b_chen_forward_mr_v1_freeze_manifest.R [--project-root <path>]", call. = FALSE)
}
setwd(root)

for (pkg in c("jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

rel <- function(...) file.path(root, ...)
norm <- function(x) gsub("\\\\", "/", x)
relpath <- function(x) norm(sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", norm(normalizePath(x, winslash = "/", mustWork = FALSE))))
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
write_csv_precise <- function(x, path) {
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}
write_json <- function(x, path) {
  partial <- paste0(path, ".partial")
  jsonlite::write_json(x, partial, pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}
write_text <- function(lines, path) {
  partial <- paste0(path, ".partial")
  writeLines(lines, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))

paths <- list(
  decision94 = rel("docs", "decisions", "94_chen_forward_harmonised_mr_inputs_v1_freeze_v1.1.md"),
  decision98 = rel("docs", "decisions", "98_chen_forward_mr_method_alignment_amendment_v1_readback_audit_v1.1.md"),
  decision99 = rel("docs", "decisions", "99_chen_forward_mr_analysis_contract_v2_v1.1.md"),
  decision100 = rel("docs", "decisions", "100_chen_forward_mr_v1_v1.1.md"),
  decision101 = rel("docs", "decisions", "101_chen_forward_mr_presso_technical_recovery_v1_v1.1.md"),
  decision102 = rel("docs", "decisions", "102_chen_forward_mr_v1_presso_completion_closure_v1.1.md"),
  script25 = rel("R", "25_chen_forward_mr_v1.R"),
  script25a = rel("R", "25a_chen_forward_mr_presso_recovery_v1.R"),
  script25a2 = rel("R", "25a2_chen_forward_mr_presso_recovery_v1_finalize_bounded_attempt.R"),
  script25a3 = rel("R", "25a3_chen_forward_mr_v1_presso_completion_closure.R"),
  script25b = rel("R", "25b_chen_forward_mr_v1_freeze_manifest.R"),
  estimates = rel("results", "tables", "chen_forward_mr_estimates_v1.csv"),
  heterogeneity = rel("results", "tables", "chen_forward_heterogeneity_v1.csv"),
  egger = rel("results", "tables", "chen_forward_egger_intercept_v1.csv"),
  original_presso = rel("results", "tables", "chen_forward_mr_presso_v1.csv"),
  recovery_presso = rel("results", "tables", "chen_forward_mr_presso_recovery_v1.csv"),
  loo = rel("results", "tables", "chen_forward_leave_one_out_v1.csv"),
  comparison = rel("results", "tables", "chen_forward_vuckovic_comparison_v1.csv"),
  mr_qc = rel("results", "qc", "chen_forward_mr_v1.json"),
  recovery_qc = rel("results", "qc", "chen_forward_mr_presso_recovery_v1.json"),
  closure_qc = rel("results", "qc", "chen_forward_mr_v1_presso_completion_closure.json"),
  mr_log = rel("results", "logs", "chen_forward_mr_v1.log"),
  recovery_log = rel("results", "logs", "chen_forward_mr_presso_recovery_v1.log"),
  closure_log = rel("results", "logs", "chen_forward_mr_v1_presso_completion_closure.log"),
  renv_lock = rel("renv.lock"),
  manifest = rel("results", "qc", "chen_forward_mr_v1_freeze_manifest.csv"),
  freeze_json = rel("results", "qc", "chen_forward_mr_v1_freeze.json"),
  freeze_log = rel("results", "logs", "chen_forward_mr_v1_freeze.log"),
  decision = rel("docs", "decisions", "103_chen_forward_mr_v1_freeze_v1.1.md")
)

required <- unlist(paths[setdiff(names(paths), c("manifest", "freeze_json", "freeze_log", "decision"))])
missing <- required[!file.exists(required)]
if (length(missing) > 0L) stop("Missing required input(s): ", paste(missing, collapse = "; "), call. = FALSE)

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 103L), paste0("Expected next decision 103, found ", next_decision, "; no outputs written."))

targets <- unlist(paths[c("manifest", "freeze_json", "freeze_log", "decision")])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(occupied, collapse = "; ")))

renv_before <- hash_file(paths$renv_lock)
mr_qc <- read_json(paths$mr_qc)
recovery_qc <- read_json(paths$recovery_qc)
closure_qc <- read_json(paths$closure_qc)
est <- read.csv(paths$estimates, stringsAsFactors = FALSE, check.names = FALSE)
het <- read.csv(paths$heterogeneity, stringsAsFactors = FALSE, check.names = FALSE)
egger <- read.csv(paths$egger, stringsAsFactors = FALSE, check.names = FALSE)
original_presso <- read.csv(paths$original_presso, stringsAsFactors = FALSE, check.names = FALSE)
recovery_presso <- read.csv(paths$recovery_presso, stringsAsFactors = FALSE, check.names = FALSE)
loo <- read.csv(paths$loo, stringsAsFactors = FALSE, check.names = FALSE)
comparison <- read.csv(paths$comparison, stringsAsFactors = FALSE, check.names = FALSE)

loo_summary <- do.call(rbind, lapply(split(loo, loo$analysis_set), function(x) {
  i <- which.max(as.numeric(x$absolute_shift))
  data.frame(
    analysis_set = x$analysis_set[[1L]],
    total_removed_snp_rows = nrow(x),
    max_absolute_shift = as.numeric(x$absolute_shift[[i]]),
    max_shift_rsid = x$removed_rsid[[i]],
    any_sign_change = any(x$sign_change == "True" | x$sign_change == TRUE),
    any_nominal_significance_change = any(x$nominal_significance_change == "True" | x$nominal_significance_change == TRUE),
    stringsAsFactors = FALSE
  )
}))

manifest_records <- data.frame(
  relative_path = relpath(unlist(paths[setdiff(names(paths), c("manifest", "freeze_json", "freeze_log", "decision"))])),
  file_role = setdiff(names(paths), c("manifest", "freeze_json", "freeze_log", "decision")),
  scientific_authority = TRUE,
  file_size_bytes = as.numeric(file.info(unlist(paths[setdiff(names(paths), c("manifest", "freeze_json", "freeze_log", "decision"))]))$size),
  sha256 = vapply(unlist(paths[setdiff(names(paths), c("manifest", "freeze_json", "freeze_log", "decision"))]), hash_file, character(1)),
  stringsAsFactors = FALSE
)

partial_manifest <- paste0(paths$manifest, ".partial")
on.exit(unlink(partial_manifest, force = TRUE), add = TRUE)
write_csv_precise(manifest_records, partial_manifest)
manifest_sha <- hash_file(partial_manifest)
renv_after <- hash_file(paths$renv_lock)

hard_checks <- list(
  decision_100_core_mr_gate = identical(mr_qc$mr_status, "passed") &&
    length(mr_qc$hard_check_failures) == 0L &&
    isTRUE(mr_qc$approved_for_chen_forward_results_interpretation),
  mr_presso_completion_closure_gate = identical(closure_qc$closure_status, "passed") &&
    isTRUE(closure_qc$approved_for_chen_forward_mr_results_freeze) &&
    length(closure_qc$hard_check_failures) == 0L,
  main_results_reverified = nrow(est) == 10L &&
    all(c("beta", "se", "ci_lower", "ci_upper", "pval", "OR", "OR_ci_lower", "OR_ci_upper", "method_id", "nsnp") %in% names(est)) &&
    all(abs(est$OR - exp(est$beta)) <= 1e-12),
  heterogeneity_reverified = nrow(het) == 4L && all(c("Q", "df", "pval") %in% names(het)),
  egger_intercept_reverified = nrow(egger) == 2L && all(c("intercept", "se", "pval") %in% names(egger)),
  loo_reverified = nrow(loo) == as.integer(mr_qc$included_nsnp) + as.integer(mr_qc$excluded_nsnp) &&
    all(c("removed_rsid", "absolute_shift", "sign_change", "nominal_significance_change") %in% names(loo)),
  mr_presso_status_truthful = identical(closure_qc$mr_presso_final_status, "technically_unavailable_under_frozen_configuration") &&
    all(original_presso$value == "failed_or_timeout") &&
    all(grepl("not_completed|not_started", recovery_presso$value)),
  no_single_snp = isFALSE(mr_qc$single_snp_run) && !file.exists(rel("results", "tables", "chen_forward_single_snp_v1.csv")),
  no_steiger = isFALSE(mr_qc$steiger_run),
  no_posthoc_filtering = isTRUE(mr_qc$hard_checks$no_posthoc_filtering),
  independent_replication_false = isFALSE(mr_qc$independent_replication),
  no_mr_rerun = TRUE,
  renv_lock_unchanged = identical(renv_before, renv_after)
)
failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
freeze_status <- if (length(failures) == 0L) "passed" else "failed"

included_results <- est[est$analysis_set == "APOE_included", , drop = FALSE]
excluded_results <- est[est$analysis_set == "APOE_excluded", , drop = FALSE]

freeze <- list(
  freeze_version = "v1",
  date = "2026-08-13",
  authoritative_chen_forward_mr_version = "v1",
  analysis_direction = "Hb_to_delirium",
  analysis_role = "forward_alternative_hb_gwas_sensitivity",
  independent_replication = FALSE,
  source_mr_input_freeze_decision = 94,
  mr_contract_decision = 99,
  mr_execution_decision = 100,
  mr_presso_recovery_decision = 101,
  mr_presso_completion_closure_decision = 102,
  mr_presso_completion_status = closure_qc$mr_presso_final_status,
  included_nsnp = mr_qc$included_nsnp,
  excluded_nsnp = mr_qc$excluded_nsnp,
  included_results = records(included_results),
  excluded_results = records(excluded_results),
  heterogeneity_results = records(het),
  egger_intercept_results = records(egger),
  mr_presso_results_or_status = list(
    original_decision100 = records(original_presso),
    recovery_decision101 = records(recovery_presso),
    final_status = closure_qc$mr_presso_final_status
  ),
  leave_one_out_summary = records(loo_summary),
  single_snp_run = FALSE,
  steiger_run = FALSE,
  vuckovic_robustness_comparison = records(comparison),
  manifest_path = relpath(paths$manifest),
  manifest_sha256 = manifest_sha,
  freeze_status = freeze_status,
  approved_for_chen_reverse_sensitivity_design = identical(freeze_status, "passed"),
  hard_checks = hard_checks,
  hard_check_failures = failures,
  informational_findings = list(
    chen_analysis_is_sensitivity_not_independent_replication = TRUE,
    mr_presso_not_negative_and_no_outlier_statement_not_supported = TRUE,
    reverse_chen_sensitivity_not_started = TRUE,
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after
  )
)

decision_lines <- c(
  "# Decision 103: Chen Forward MR V1 Freeze",
  "",
  "Date: 2026-08-13",
  "",
  "## Status",
  paste0("freeze_status: `", freeze_status, "`"),
  paste0("approved_for_chen_reverse_sensitivity_design: `", identical(freeze_status, "passed"), "`"),
  "",
  "## Decision",
  "Chen Forward MR V1 results are frozen after Decision 102 MR-PRESSO completion closure.",
  "",
  "MR-PRESSO final status is `technically_unavailable_under_frozen_configuration`; no Global Test P value or no-outlier conclusion is inferred.",
  "",
  "Chen Forward MR remains an alternative-Hb-GWAS robustness sensitivity analysis and is not an independent replication.",
  "",
  "## Hard Check Failures",
  if (length(failures) == 0L) "- none" else paste0("- `", failures, "`"),
  "",
  "## Outputs",
  paste0("- `", relpath(paths$manifest), "`"),
  paste0("- `", relpath(paths$freeze_json), "`"),
  paste0("- `", relpath(paths$freeze_log), "`"),
  paste0("- `", relpath(paths$decision), "`")
)

log_lines <- c(
  "[2026-08-13] Chen Forward MR V1 freeze",
  paste0("freeze_status=", freeze_status),
  paste0("manifest_sha256=", manifest_sha),
  paste0("approved_for_chen_reverse_sensitivity_design=", identical(freeze_status, "passed")),
  paste0("hard_check_failures=", if (length(failures) == 0L) "[]" else paste(failures, collapse = ";")),
  "No MR was rerun."
)

stop_if(length(failures) > 0L, paste("Freeze hard checks failed:", paste(failures, collapse = "; ")))
if (!file.rename(partial_manifest, paths$manifest)) stop("Atomic rename failed: ", paths$manifest, call. = FALSE)
write_json(freeze, paths$freeze_json)
write_text(log_lines, paths$freeze_log)
write_text(decision_lines, paths$decision)

cat("freeze_status=", freeze_status, "\n", sep = "")
cat("manifest_sha256=", manifest_sha, "\n", sep = "")
cat("approved_for_chen_reverse_sensitivity_design=", identical(freeze_status, "passed"), "\n", sep = "")
cat("hard_check_failures=[]\n")
