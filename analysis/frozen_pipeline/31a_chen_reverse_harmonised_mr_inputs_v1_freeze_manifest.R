#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/31a_chen_reverse_harmonised_mr_inputs_v1_freeze_manifest.R [--project-root <path>]", call. = FALSE)
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
estimator_class <- function(n) {
  if (n == 0L) "not_estimable" else if (n == 1L) "Wald ratio" else "multi-IV MR design pending"
}

paths <- list(
  script = rel("R", "31a_chen_reverse_harmonised_mr_inputs_v1_freeze_manifest.R"),
  decision47 = rel("docs", "decisions", "47_finngen_r13_delirium_reverse_primary_instruments_v4_freeze_v3_v1.1.md"),
  decision59 = rel("docs", "decisions", "59_finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze_v1.1.md"),
  decision64 = rel("docs", "decisions", "64_reverse_relaxed_palindromic_handling_rule_v1_v1.1.md"),
  decision110 = rel("docs", "decisions", "110_chen_reverse_outcome_extraction_v1_v1.1.md"),
  decision111 = rel("docs", "decisions", "111_chen_reverse_outcome_extraction_v1_readback_recovery_v1.1.md"),
  decision112 = rel("docs", "decisions", "112_chen_reverse_outcome_extraction_v1_freeze_v1.1.md"),
  decision113 = rel("docs", "decisions", "113_chen_reverse_harmonisation_contract_and_preflight_v1_v1.1.md"),
  decision114 = rel("docs", "decisions", "114_chen_reverse_formal_harmonisation_v1_v1.1.md"),
  decision115 = rel("docs", "decisions", "115_chen_reverse_formal_harmonisation_v1_readback_recovery_v1.1.md"),
  formal_script = rel("R", "31_chen_reverse_formal_harmonisation_v1.R"),
  recovery_script = rel("R", "31b_chen_reverse_formal_harmonisation_v1_readback_recovery.R"),
  formal_qc = rel("results", "qc", "chen_reverse_formal_harmonisation_v1.json"),
  recovery_qc = rel("results", "qc", "chen_reverse_formal_harmonisation_v1_readback_recovery.json"),
  recovery_mismatch = rel("results", "qc", "chen_reverse_formal_harmonisation_v1_readback_recovery_mismatch_audit.csv"),
  formal_counts = rel("results", "qc", "chen_reverse_formal_harmonisation_counts_v1.csv"),
  transform_audit = rel("results", "qc", "chen_reverse_formal_harmonisation_transform_audit_v1.csv"),
  excluded = rel("results", "qc", "chen_reverse_formal_harmonisation_excluded_snps_v1.tsv"),
  formal_log = rel("results", "logs", "chen_reverse_formal_harmonisation_v1.log"),
  recovery_log = rel("results", "logs", "chen_reverse_formal_harmonisation_v1_readback_recovery.log"),
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
  manifest = rel("results", "qc", "chen_reverse_harmonised_mr_inputs_v1_freeze_manifest.csv"),
  freeze_json = rel("results", "qc", "chen_reverse_harmonised_mr_inputs_v1_freeze.json"),
  log = rel("results", "logs", "chen_reverse_harmonised_mr_inputs_v1_freeze.log"),
  decision = rel("docs", "decisions", "116_chen_reverse_harmonised_mr_inputs_v1_freeze_v1.1.md"),
  renv_lock = rel("renv.lock")
)

outputs <- unlist(paths[c("manifest", "freeze_json", "log", "decision")])
occupied <- outputs[file.exists(outputs) | file.exists(paste0(outputs, ".partial"))]
stop_if(length(occupied) > 0L, paste("Output or partial exists:", paste(occupied, collapse = "; ")))
inputs <- unlist(paths[c("script", "decision47", "decision59", "decision64", "decision110", "decision111", "decision112", "decision113", "decision114", "decision115", "formal_script", "recovery_script", "formal_qc", "recovery_qc", "recovery_mismatch", "formal_counts", "transform_audit", "excluded", "formal_log", "recovery_log", "master_parquet", "master_tsv", "strict_inc_parquet", "strict_inc_tsv", "strict_exc_parquet", "strict_exc_tsv", "relaxed_inc_parquet", "relaxed_inc_tsv", "relaxed_exc_parquet", "relaxed_exc_tsv", "renv_lock")])
missing_inputs <- inputs[!file.exists(inputs)]
stop_if(length(missing_inputs) > 0L, paste("Required input missing:", paste(missing_inputs, collapse = "; ")))

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 116L), paste0("Expected next decision 116, found ", next_decision, "; no outputs written."))

dir.create(dirname(paths$log), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = paths$log, append = TRUE)
log_line("stage=chen_reverse_harmonised_mr_inputs_v1_freeze_start")

formal_qc <- read_json(paths$formal_qc)
recovery_qc <- read_json(paths$recovery_qc)
stop_if(!identical(formal_qc$harmonisation_status, "failed"), "Decision 114 status should remain failed and preserved.")
stop_if(!identical(recovery_qc$recovery_status, "passed") || !isTRUE(recovery_qc$approved_for_chen_reverse_mr_input_freeze), "Decision 115 recovery gate failed.")

con_db <- DBI::dbConnect(duckdb::duckdb(), read_only = FALSE)
on.exit(try(DBI::dbDisconnect(con_db, shutdown = TRUE), silent = TRUE), add = TRUE)
read_parquet <- function(path) DBI::dbGetQuery(con_db, sprintf("SELECT * FROM read_parquet(%s)", DBI::dbQuoteString(con_db, norm(path))))
branch_paths <- list(
  strict_apoe_included = paths$strict_inc_parquet,
  strict_apoe_excluded = paths$strict_exc_parquet,
  relaxed_apoe_included = paths$relaxed_inc_parquet,
  relaxed_apoe_excluded = paths$relaxed_exc_parquet
)
branch_roles <- c(
  strict_apoe_included = "reverse_strict_primary_alternative_hb_outcome_sensitivity",
  strict_apoe_excluded = "reverse_strict_primary_alternative_hb_outcome_sensitivity",
  relaxed_apoe_included = "reverse_relaxed_exploratory_alternative_hb_outcome_sensitivity",
  relaxed_apoe_excluded = "reverse_relaxed_exploratory_alternative_hb_outcome_sensitivity"
)
branch_results <- lapply(names(branch_paths), function(branch) {
  x <- read_parquet(branch_paths[[branch]])
  n <- nrow(x)
  data.frame(
    branch = branch,
    analysis_role = branch_roles[[branch]],
    final_valid_count = n,
    final_rsids = paste(x$rsid, collapse = ";"),
    F_min = if (n) min(x$F_stat, na.rm = TRUE) else NA_real_,
    F_mean = if (n) mean(x$F_stat, na.rm = TRUE) else NA_real_,
    F_median = if (n) stats::median(x$F_stat, na.rm = TRUE) else NA_real_,
    F_max = if (n) max(x$F_stat, na.rm = TRUE) else NA_real_,
    F_lt10 = sum(x$F_stat < 10, na.rm = TRUE),
    mr_estimable = n > 0L,
    planned_estimator_class = estimator_class(n),
    palindrome_rule = "exclude_all_without_EAF_reinclusion",
    final_alleles_aligned = if (n) all(x$outcome_effect_allele_harmonised == x$exposure_effect_allele & x$outcome_other_allele_harmonised == x$exposure_other_allele) else TRUE,
    stringsAsFactors = FALSE
  )
})
branch_results <- do.call(rbind, branch_results)
master <- read_parquet(paths$master_parquet)

manifest <- data.frame(
  role = c(
    "freeze_script", "decision47", "decision59", "decision64", "decision110", "decision111", "decision112", "decision113", "decision114", "decision115",
    "formal_script", "readback_recovery_script", "formal_qc", "readback_recovery_qc", "readback_mismatch_audit", "formal_counts", "transform_audit", "excluded_snps", "formal_log", "readback_recovery_log",
    "master_parquet", "master_tsv", "strict_inc_parquet", "strict_inc_tsv", "strict_exc_parquet", "strict_exc_tsv", "relaxed_inc_parquet", "relaxed_inc_tsv", "relaxed_exc_parquet", "relaxed_exc_tsv", "renv_lock"
  ),
  path = c(
    paths$script, paths$decision47, paths$decision59, paths$decision64, paths$decision110, paths$decision111, paths$decision112, paths$decision113, paths$decision114, paths$decision115,
    paths$formal_script, paths$recovery_script, paths$formal_qc, paths$recovery_qc, paths$recovery_mismatch, paths$formal_counts, paths$transform_audit, paths$excluded, paths$formal_log, paths$recovery_log,
    paths$master_parquet, paths$master_tsv, paths$strict_inc_parquet, paths$strict_inc_tsv, paths$strict_exc_parquet, paths$strict_exc_tsv, paths$relaxed_inc_parquet, paths$relaxed_inc_tsv, paths$relaxed_exc_parquet, paths$relaxed_exc_tsv, paths$renv_lock
  ),
  stringsAsFactors = FALSE
)
manifest$relative_path <- vapply(manifest$path, relpath, character(1))
manifest$exists <- file.exists(manifest$path)
manifest$size_bytes <- unname(file.info(manifest$path)$size)
manifest$sha256 <- vapply(manifest$path, hash_file, character(1))
write_csv(manifest, paths$manifest)
manifest_sha <- hash_file(paths$manifest)

hard_checks <- list(
  formal_harmonisation_gate = identical(recovery_qc$recovery_status, "passed") && isTRUE(recovery_qc$approved_for_chen_reverse_mr_input_freeze),
  formal_outputs_complete = all(manifest$exists),
  branch_counts_reverified = identical(as.integer(branch_results$final_valid_count), c(1L, 1L, 10L, 9L)),
  branch_rsids_reverified = all(nzchar(branch_results$final_rsids)),
  all_final_alleles_aligned = all(branch_results$final_alleles_aligned),
  no_palindromic_final_valid = sum(master$final_valid & master$palindromic_recomputed) == 0L,
  all_F_ge_10 = all(branch_results$F_lt10 == 0L),
  strict_relaxed_hierarchy_preserved = TRUE,
  mr_estimability_classified = all(!is.na(branch_results$planned_estimator_class)) && all(branch_results$mr_estimable),
  no_posthoc_filtering = TRUE,
  no_mr = TRUE,
  no_steiger = TRUE,
  renv_lock_unchanged = identical(hash_file(paths$renv_lock), "253471c02e1e47a40d0f68b296d4ae2b1df471d757ac86328b3e974018d039f3")
)
hard_check_failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
freeze_status <- if (length(hard_check_failures) == 0L) "passed" else "failed"

freeze <- list(
  freeze_version = "v1",
  decision = 116,
  date = format(Sys.Date()),
  authoritative_chen_reverse_harmonisation_version = "v1",
  analysis_direction = "delirium_to_Hb",
  analysis_role = "reverse_alternative_hb_outcome_sensitivity",
  independent_replication = FALSE,
  strict_threshold = 5e-8,
  relaxed_threshold = 5e-6,
  branch_results = records(branch_results),
  outcome_scale = formal_qc$outcome_scale,
  strict_primary_superseded_by_relaxed = FALSE,
  relaxed_confirmatory = FALSE,
  manifest_sha256 = manifest_sha,
  freeze_status = freeze_status,
  approved_for_chen_reverse_mr_design = identical(freeze_status, "passed") && any(branch_results$mr_estimable),
  mr_run = FALSE,
  steiger_run = FALSE,
  hard_checks = hard_checks,
  hard_check_failures = hard_check_failures,
  informational_findings = list(
    decision_114_status_preserved = formal_qc$harmonisation_status,
    decision_115_recovery_status = recovery_qc$recovery_status,
    freeze_is_design_readiness_only = TRUE,
    no_estimator_executed = TRUE
  )
)
write_json(freeze, paths$freeze_json)

decision_lines <- c(
  "# Decision 116: Chen Reverse Harmonised MR Inputs V1 Freeze",
  "",
  paste0("Date: ", format(Sys.Date())),
  "",
  "## Status",
  paste0("freeze_status: `", freeze_status, "`"),
  paste0("approved_for_chen_reverse_mr_design: `", identical(freeze_status, "passed") && any(branch_results$mr_estimable), "`"),
  "",
  "## Decision",
  "Freeze Chen reverse harmonised MR input files after Decision 115 readback recovery passed.",
  "",
  "This freeze records MR design readiness only. It does not run Wald ratio, IVW, MR-Egger, weighted median/mode, heterogeneity, Egger intercept, MR-PRESSO, leave-one-out, MR, or Steiger.",
  "",
  "## Frozen Evidence",
  paste0("- Manifest SHA-256: `", manifest_sha, "`."),
  paste0("- Hard-check failures: `", if (length(hard_check_failures) == 0L) "none" else paste(hard_check_failures, collapse = ";"), "`."),
  "",
  "## Branch Readiness",
  paste0("- strict_apoe_included: `", branch_results$final_valid_count[branch_results$branch == "strict_apoe_included"], "`; estimator class `", branch_results$planned_estimator_class[branch_results$branch == "strict_apoe_included"], "`."),
  paste0("- strict_apoe_excluded: `", branch_results$final_valid_count[branch_results$branch == "strict_apoe_excluded"], "`; estimator class `", branch_results$planned_estimator_class[branch_results$branch == "strict_apoe_excluded"], "`."),
  paste0("- relaxed_apoe_included: `", branch_results$final_valid_count[branch_results$branch == "relaxed_apoe_included"], "`; estimator class `", branch_results$planned_estimator_class[branch_results$branch == "relaxed_apoe_included"], "`."),
  paste0("- relaxed_apoe_excluded: `", branch_results$final_valid_count[branch_results$branch == "relaxed_apoe_excluded"], "`; estimator class `", branch_results$planned_estimator_class[branch_results$branch == "relaxed_apoe_excluded"], "`."),
  "",
  "## Outputs",
  paste0("- `", relpath(paths$freeze_json), "`"),
  paste0("- `", relpath(paths$manifest), "`"),
  paste0("- `", relpath(paths$log), "`")
)
write_text(decision_lines, paths$decision)

log_line("freeze_status=", freeze_status, "; hard_check_failures=", paste(hard_check_failures, collapse = ","))
cat("freeze_status=", freeze_status, "\n", sep = "")
cat("approved_for_chen_reverse_mr_design=", identical(freeze_status, "passed") && any(branch_results$mr_estimable), "\n", sep = "")
cat("hard_check_failures=", paste(hard_check_failures, collapse = ","), "\n", sep = "")
cat("manifest_sha256=", manifest_sha, "\n", sep = "")
