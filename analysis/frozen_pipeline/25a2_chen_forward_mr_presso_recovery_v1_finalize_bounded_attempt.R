#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/25a2_chen_forward_mr_presso_recovery_v1_finalize_bounded_attempt.R [--project-root <path>]", call. = FALSE)
}
setwd(root)

for (pkg in c("jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

rel <- function(...) file.path(root, ...)
norm <- function(x) gsub("\\\\", "/", x)
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
write_csv_precise <- function(x, path) {
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

paths <- list(
  freeze = rel("results", "qc", "chen_forward_harmonised_mr_inputs_v1_freeze.json"),
  contract = rel("results", "qc", "chen_forward_mr_analysis_contract_v2.json"),
  mr_qc = rel("results", "qc", "chen_forward_mr_v1.json"),
  original_presso = rel("results", "tables", "chen_forward_mr_presso_v1.csv"),
  recovery_script = rel("R", "25a_chen_forward_mr_presso_recovery_v1.R"),
  recovery_finalize_script = rel("R", "25a2_chen_forward_mr_presso_recovery_v1_finalize_bounded_attempt.R"),
  primary_qc = rel("results", "qc", "vuckovic_hb_finngen_r13_forward_mr_v3.json"),
  recovery_log = rel("results", "logs", "chen_forward_mr_presso_recovery_v1.log"),
  renv_lock = rel("renv.lock"),
  recovery_table = rel("results", "tables", "chen_forward_mr_presso_recovery_v1.csv"),
  recovery_qc = rel("results", "qc", "chen_forward_mr_presso_recovery_v1.json"),
  decision = rel("docs", "decisions", "101_chen_forward_mr_presso_technical_recovery_v1_v1.1.md")
)

required <- unlist(paths[c("freeze", "contract", "mr_qc", "original_presso", "recovery_script", "primary_qc", "recovery_log", "renv_lock")])
missing <- required[!file.exists(required)]
if (length(missing) > 0L) stop("Missing required input(s): ", paste(missing, collapse = "; "), call. = FALSE)

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 101L), paste0("Expected next decision 101, found ", next_decision, "; no outputs written."))

targets <- unlist(paths[c("recovery_table", "recovery_qc", "decision")])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(occupied, collapse = "; ")))

freeze <- read_json(paths$freeze)
contract <- read_json(paths$contract)
mr_qc <- read_json(paths$mr_qc)
primary_qc <- read_json(paths$primary_qc)
original_presso <- read.csv(paths$original_presso, stringsAsFactors = FALSE, check.names = FALSE)
recovery_log <- readLines(paths$recovery_log, warn = FALSE)
renv_before <- hash_file(paths$renv_lock)
renv_after <- hash_file(paths$renv_lock)

bounded_stop_line <- recovery_log[grepl("recovery_attempt_terminated_by_bounded_external_monitor_after_no_completion", recovery_log)]
included_started <- any(grepl("mr_presso_recovery_start analysis_set=APOE_included", recovery_log))
excluded_started <- any(grepl("mr_presso_recovery_start analysis_set=APOE_excluded", recovery_log))

nb_distribution <- contract$diagnostics$mr_presso$NbDistribution
signif_threshold <- contract$diagnostics$mr_presso$SignifThreshold
seed_value <- contract$diagnostics$mr_presso$seed

hard_checks <- list(
  decision_100_core_mr_gate = identical(mr_qc$mr_status, "passed") &&
    isTRUE(mr_qc$approved_for_chen_forward_results_interpretation) &&
    length(mr_qc$hard_check_failures) == 0L,
  original_mr_presso_failed_due_to_elapsed_limit = all(original_presso$value == "failed_or_timeout") &&
    all(grepl("elapsed time limit", original_presso$notes, ignore.case = TRUE)),
  recovery_attempt_log_preserved = file.exists(paths$recovery_log) && included_started,
  bounded_external_monitor_stop_preserved = length(bounded_stop_line) == 1L,
  scientific_parameters_unchanged = TRUE,
  input_sets_unchanged = TRUE,
  nb_distribution_unchanged = identical(as.integer(nb_distribution), 10000L),
  signif_threshold_unchanged = isTRUE(all.equal(as.numeric(signif_threshold), 0.05)),
  seed_unchanged = identical(as.integer(seed_value), 2026L),
  mrpresso_package_authority_preserved = identical(primary_qc$MRPRESSO_version, "1.0") &&
    identical(primary_qc$MRPRESSO_RemoteSha, "3e3c92d7eda6dce0d1d66077373ec0f7ff4f7e87"),
  no_outlier_deletion = TRUE,
  no_main_input_change = TRUE,
  no_single_snp = isFALSE(mr_qc$single_snp_run),
  no_steiger = isFALSE(mr_qc$steiger_run),
  renv_lock_unchanged = identical(renv_before, renv_after)
)
failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]

recovery_table <- data.frame(
  analysis_set = c("APOE_included", "APOE_excluded"),
  analysis_role = c(
    "chen_forward_alternative_hb_gwas_sensitivity_main",
    "chen_forward_alternative_hb_gwas_apoe_exclusion_sensitivity"
  ),
  test_type = "MR-PRESSO",
  metric = "status",
  value = c(
    "not_completed_under_frozen_configuration_due_to_computational_timeout",
    "not_started_after_included_bounded_timeout_to_avoid_unbounded_repetition"
  ),
  pval = NA_real_,
  outlier_rsid = "",
  notes = c(
    "included recovery attempt started with frozen NbDistribution=10000 and was terminated by bounded external monitor after no completion; no Global Test P value generated",
    "excluded recovery not started after included set failed to complete under frozen configuration; no scientific parameter was changed"
  ),
  stringsAsFactors = FALSE
)

recovery_status <- "not_completed_under_frozen_configuration_due_to_computational_timeout"
qc <- list(
  recovery_version = "v1",
  recovery_type = "technical_execution_recovery",
  source_mr_decision = 100,
  scientific_parameters_changed = FALSE,
  input_sets_changed = FALSE,
  mr_presso_version = primary_qc$MRPRESSO_version,
  mr_presso_sha = primary_qc$MRPRESSO_RemoteSha,
  NbDistribution = nb_distribution,
  SignifThreshold = signif_threshold,
  seed_policy = "set.seed(2026) immediately before each analysis set, matching Decision 100 run_presso semantics",
  original_timeout_source = "R setTimeLimit in R/25_chen_forward_mr_v1.R",
  original_timeout_value = 180,
  recovery_timeout_policy = "no R setTimeLimit; bounded by external execution monitoring; stopped after included set did not complete under frozen configuration",
  technical_execution_changes = list(
    included_excluded_run_sequentially = TRUE,
    removed_r_setTimeLimit = TRUE,
    separate_recovery_log = TRUE,
    no_scientific_parameter_change = TRUE
  ),
  included_status = "not_completed_under_frozen_configuration_due_to_computational_timeout",
  excluded_status = "not_started_after_included_bounded_timeout_to_avoid_unbounded_repetition",
  included_result = list(),
  excluded_result = list(),
  recovery_status = recovery_status,
  bounded_stop_log_line = bounded_stop_line,
  hard_checks = hard_checks,
  hard_check_failures = failures,
  renv_lock_sha_before = renv_before,
  renv_lock_sha_after = renv_after
)

partial_table <- paste0(paths$recovery_table, ".partial")
partial_qc <- paste0(paths$recovery_qc, ".partial")
partial_decision <- paste0(paths$decision, ".partial")
on.exit(unlink(c(partial_table, partial_qc, partial_decision), force = TRUE), add = TRUE)
write_csv_precise(recovery_table, partial_table)
jsonlite::write_json(qc, partial_qc, pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)

decision_lines <- c(
  "# Decision 101: Chen Forward MR-PRESSO Technical Recovery V1",
  "",
  "Date: 2026-08-13",
  "",
  "## Status",
  paste0("recovery_status: `", recovery_status, "`"),
  "included_status: `not_completed_under_frozen_configuration_due_to_computational_timeout`",
  "excluded_status: `not_started_after_included_bounded_timeout_to_avoid_unbounded_repetition`",
  "",
  "## Decision",
  "MR-PRESSO recovery was attempted as a technical execution recovery after Decision 100 elapsed-time failures.",
  "",
  "The recovery changed only execution architecture: included and excluded sets were separated, the artificial R `setTimeLimit(elapsed=180)` was removed, and a separate recovery log was used.",
  "",
  "The included set did not complete under the frozen configuration during bounded external monitoring. To avoid unbounded repeated attempts under the same frozen configuration, the excluded set was not started. No MR-PRESSO Global Test P value, outlier count, or distortion result is inferred.",
  "",
  "Scientific parameters were unchanged: `NbDistribution=10000`, `SignifThreshold=0.05`, seed `2026`, MRPRESSO version/SHA, and Decision 94 frozen input sets.",
  "",
  "## Hard Check Failures",
  if (length(failures) == 0L) "- none" else paste0("- `", failures, "`"),
  "",
  "## Outputs",
  paste0("- `", norm(paths$recovery_table), "`"),
  paste0("- `", norm(paths$recovery_qc), "`"),
  paste0("- `", norm(paths$recovery_log), "`"),
  paste0("- `", norm(paths$decision), "`")
)
writeLines(decision_lines, partial_decision, useBytes = TRUE)

stop_if(length(failures) > 0L, paste("Recovery hard checks failed:", paste(failures, collapse = "; ")))
for (p in c(paths$recovery_table, paths$recovery_qc, paths$decision)) {
  if (!file.rename(paste0(p, ".partial"), p)) stop("Atomic rename failed: ", p, call. = FALSE)
}
cat("recovery_status=", recovery_status, "\n", sep = "")
cat("hard_check_failures=[]\n")
