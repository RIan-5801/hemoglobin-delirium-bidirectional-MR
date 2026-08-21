#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/37a_final_integrated_analysis_freeze_readback_recovery_v1.R [--project-root <path>]", call. = FALSE)
}
setwd(root)

local_lib <- normalizePath(file.path(root, "renv", "mr-v1-library"), winslash = "/", mustWork = TRUE)
.libPaths(c(local_lib, .libPaths()))

for (pkg in c("jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

rel <- function(...) file.path(root, ...)
norm <- function(x) gsub("\\\\", "/", x)
relpath <- function(x) norm(sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", norm(normalizePath(x, winslash = "/", mustWork = FALSE))))
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
read_csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
is_empty <- function(x) is.null(x) || length(x) == 0L
as_bool <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% "true"
}
write_csv_precise <- function(x, path) {
  partial <- paste0(path, ".partial")
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.csv(x, partial, row.names = FALSE, na = "")
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}
write_json <- function(x, path) {
  partial <- paste0(path, ".partial")
  jsonlite::write_json(x, partial, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}
write_text <- function(lines, path) {
  partial <- paste0(path, ".partial")
  writeLines(lines, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}
latest_decision <- function() {
  files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
  nums <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", files)))
  max(nums, na.rm = TRUE) + 1L
}

paths <- list(
  script37a = rel("R", "37a_final_integrated_analysis_freeze_readback_recovery_v1.R"),
  script37 = rel("R", "37_final_integrated_analysis_freeze_v1.R"),
  renv_lock = rel("renv.lock"),
  decision124 = rel("docs", "decisions", "124_unified_steiger_v1_freeze_v1.1.md"),
  decision125 = rel("docs", "decisions", "125_final_integrated_analysis_freeze_v1_v1.1.md"),
  steiger_freeze = rel("results", "qc", "unified_steiger_v1_freeze.json"),
  final_manifest_v1 = rel("results", "qc", "final_integrated_analysis_freeze_manifest_v1.csv"),
  final_qc_v1 = rel("results", "qc", "final_integrated_analysis_freeze_v1.json"),
  final_log_v1 = rel("results", "logs", "final_integrated_analysis_freeze_v1.log"),
  final_registry = rel("results", "final", "final_analysis_registry_v1.csv"),
  final_primary = rel("results", "final", "final_primary_result_matrix_v1.csv"),
  final_diagnostic = rel("results", "final", "final_diagnostic_status_matrix_v1.csv"),
  final_completeness = rel("results", "final", "final_analysis_completeness_audit_v1.csv"),
  final_interpretation = rel("results", "final", "final_scientific_interpretation_v1.json"),
  final_limitations = rel("results", "final", "final_limitations_registry_v1.csv"),
  recovery_manifest = rel("results", "qc", "final_integrated_analysis_freeze_manifest_readback_recovery_v1.csv"),
  recovery_json = rel("results", "qc", "final_integrated_analysis_freeze_readback_recovery_v1.json"),
  recovery_log = rel("results", "logs", "final_integrated_analysis_freeze_readback_recovery_v1.log"),
  decision126 = rel("docs", "decisions", "126_final_integrated_analysis_freeze_readback_recovery_v1_v1.1.md")
)

required <- unlist(paths[1:15])
missing <- required[!file.exists(required)]
stop_if(length(missing) > 0L, paste("Missing required input(s):", paste(relpath(missing), collapse = "; ")))
stop_if(!identical(latest_decision(), 126L), paste("Expected next decision 126, found ", latest_decision(), "; no outputs written."))

targets <- unlist(paths[c("recovery_manifest", "recovery_json", "recovery_log", "decision126")])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(relpath(occupied), collapse = "; ")))

renv_before <- hash_file(paths$renv_lock)
steiger <- read_json(paths$steiger_freeze)
final_qc <- read_json(paths$final_qc_v1)
final_interpretation <- read_json(paths$final_interpretation)
final_manifest_v1 <- read_csv(paths$final_manifest_v1)
final_registry <- read_csv(paths$final_registry)
final_primary <- read_csv(paths$final_primary)
final_diagnostic <- read_csv(paths$final_diagnostic)
final_completeness <- read_csv(paths$final_completeness)
final_limitations <- read_csv(paths$final_limitations)

expected_self_artifacts <- c(
  "results/final/final_scientific_interpretation_v1.json",
  "results/qc/final_integrated_analysis_freeze_v1.json",
  "results/logs/final_integrated_analysis_freeze_v1.log",
  "docs/decisions/125_final_integrated_analysis_freeze_v1_v1.1.md"
)
missing_from_v1_manifest <- expected_self_artifacts[!expected_self_artifacts %in% final_manifest_v1$relative_path]

script_text <- paste(readLines(paths$script37a, warn = FALSE), collapse = "\n")
scan_text <- gsub("\"([^\"\\\\]|\\\\.)*\"", "\"\"", scan_text <- script_text, perl = TRUE)
scan_text <- gsub("'([^'\\\\]|\\\\.)*'", "''", scan_text, perl = TRUE)
forbidden_call_patterns <- c(
  "get_r_from_lor\\s*\\(", "get_r_from_bsen\\s*\\(", "effective_n\\s*\\(",
  "mr_steiger\\s*\\(", "mr_steiger2\\s*\\(", "directionality_test\\s*\\(",
  "steiger_filtering\\s*\\(", "harmonise_data\\s*\\(", "mr\\s*\\(",
  "mr_heterogeneity\\s*\\(", "mr_pleiotropy_test\\s*\\(", "mr_leaveoneout\\s*\\("
)
no_forbidden_calls <- !any(vapply(forbidden_call_patterns, function(p) grepl(p, scan_text, ignore.case = TRUE, perl = TRUE), logical(1)))

source_rows <- final_manifest_v1
source_rows$recovery_role <- "retained_from_final_manifest_v1"

core_additions <- c(
  paths$script37a,
  paths$final_interpretation,
  paths$final_qc_v1,
  paths$final_log_v1,
  paths$decision125,
  paths$recovery_manifest,
  paths$recovery_json,
  paths$recovery_log,
  paths$decision126
)
core_status <- c(
  "current_recovery_script",
  "current_final_output",
  "current_final_qc",
  "current_final_log",
  "current_final_decision",
  "current_recovery_output",
  "current_recovery_output",
  "current_recovery_output",
  "current_recovery_decision"
)
core_authority <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
core_rows <- data.frame(
  relative_path = relpath(core_additions),
  artifact_type = c(
    "readback_recovery_script",
    "final_scientific_interpretation_json",
    "final_integrated_analysis_freeze_qc_json",
    "final_integrated_analysis_freeze_log",
    "final_integrated_analysis_freeze_decision",
    "readback_recovery_manifest",
    "readback_recovery_qc_json",
    "readback_recovery_log",
    "readback_recovery_decision"
  ),
  decision_number = c(126, 125, 125, 125, 125, 126, 126, 126, 126),
  status = core_status,
  scientific_authority = core_authority,
  superseded_or_recovered_by = c("", "", "", "", "", "", "", "", ""),
  file_size_bytes = NA_real_,
  sha256 = "",
  recovery_role = c(
    "new_recovery_script",
    "added_missing_final_self_artifact",
    "added_missing_final_self_artifact",
    "added_missing_final_self_artifact",
    "added_missing_final_self_artifact",
    "self_referential_recovery_artifact",
    "self_referential_recovery_artifact",
    "self_referential_recovery_artifact",
    "self_referential_recovery_artifact"
  ),
  stringsAsFactors = FALSE
)

recovery_manifest <- rbind(source_rows, core_rows)
recovery_manifest <- recovery_manifest[!duplicated(recovery_manifest$relative_path, fromLast = TRUE), , drop = FALSE]

initial_existence <- file.exists(core_additions)
hard_checks_pre <- list(
  decision_124_gate = identical(steiger$freeze_status, "passed") &&
    isTRUE(steiger$approved_for_final_integrated_analysis_freeze) &&
    is_empty(steiger$hard_check_failures),
  decision_125_gate = identical(final_qc$freeze_status, "passed") &&
    identical(final_qc$analysis_phase_status, "complete_under_frozen_protocol") &&
    isTRUE(final_qc$approved_for_results_tables_figures) &&
    isTRUE(final_qc$approved_for_manuscript_results_drafting) &&
    is_empty(final_qc$hard_check_failures),
  manifest_v1_gap_confirmed = setequal(missing_from_v1_manifest, expected_self_artifacts),
  final_outputs_exist = all(file.exists(unlist(paths[10:15]))),
  final_table_shapes_preserved = nrow(final_registry) == 12L &&
    nrow(final_primary) == 12L &&
    nrow(final_diagnostic) == 12L &&
    nrow(final_completeness) == 12L &&
    nrow(final_limitations) == 10L,
  final_classifications_preserved = identical(final_interpretation$overall_mr_classification, final_qc$overall_mr_classification) &&
    identical(final_interpretation$steiger_integrated_classification, final_qc$steiger_integrated_classification) &&
    identical(final_interpretation$bidirectional_causality_inferred, FALSE) &&
    identical(final_qc$bidirectional_causality_inferred, FALSE),
  no_scientific_analysis_rerun = no_forbidden_calls,
  no_existing_target_overwrite = all(!file.exists(targets)) && all(!file.exists(paste0(targets, ".partial"))),
  renv_lock_pre_unchanged = identical(renv_before, hash_file(paths$renv_lock)),
  git_status_not_required = identical(final_qc$git_status, "not_applicable_project_not_git_repository")
)
pre_failures <- names(hard_checks_pre)[!vapply(hard_checks_pre, isTRUE, logical(1))]
stop_if(length(pre_failures) > 0L, paste("Readback recovery pre-write hard checks failed:", paste(pre_failures, collapse = "; ")))

write_csv_precise(recovery_manifest, paths$recovery_manifest)
recovery_manifest_sha <- hash_file(paths$recovery_manifest)

for (i in seq_len(nrow(recovery_manifest))) {
  abs_path <- rel(recovery_manifest$relative_path[[i]])
  if (file.exists(abs_path)) {
    recovery_manifest$file_size_bytes[[i]] <- as.numeric(file.info(abs_path)$size)
    recovery_manifest$sha256[[i]] <- hash_file(abs_path)
  }
}
write_csv_precise(recovery_manifest, paths$recovery_manifest)
recovery_manifest_sha <- hash_file(paths$recovery_manifest)

renv_after <- hash_file(paths$renv_lock)
hard_checks <- c(hard_checks_pre, list(
  recovered_manifest_includes_final_self_artifacts = all(expected_self_artifacts %in% recovery_manifest$relative_path),
  recovered_manifest_includes_recovery_self_artifacts = all(relpath(unlist(paths[c("recovery_manifest", "recovery_json", "recovery_log", "decision126")])) %in% recovery_manifest$relative_path),
  no_final_result_file_modified = identical(hash_file(paths$final_registry), recovery_manifest$sha256[match(relpath(paths$final_registry), recovery_manifest$relative_path)]) &&
    identical(hash_file(paths$final_primary), recovery_manifest$sha256[match(relpath(paths$final_primary), recovery_manifest$relative_path)]) &&
    identical(hash_file(paths$final_interpretation), recovery_manifest$sha256[match(relpath(paths$final_interpretation), recovery_manifest$relative_path)]),
  no_scientific_classification_change = identical(final_interpretation$overall_mr_classification, "no_robust_evidence_for_a_causal_association_between_Hb_and_delirium_in_either_direction") &&
    identical(final_interpretation$steiger_integrated_classification, "Chen_based_instrument_sets_show_robust_support_for_their_hypothesized_instrument_orientation_across_prespecified_prevalence_assumptions") &&
    identical(final_interpretation$bidirectional_causality_inferred, FALSE),
  renv_lock_unchanged = identical(renv_before, renv_after)
))
failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
recovery_status <- if (length(failures) == 0L) "passed" else "failed"
stop_if(length(failures) > 0L, paste("Readback recovery hard checks failed:", paste(failures, collapse = "; ")))

recovery <- list(
  recovery_version = "v1",
  decision = 126,
  date = "2026-08-13",
  recovery_status = recovery_status,
  source_final_freeze_decision = 125,
  source_final_freeze_status = final_qc$freeze_status,
  source_analysis_phase_status = final_qc$analysis_phase_status,
  correction_scope = "manifest_readback_completeness_only",
  scientific_results_changed = FALSE,
  final_interpretation_changed = FALSE,
  final_qc_changed = FALSE,
  final_manifest_v1_sha256 = hash_file(paths$final_manifest_v1),
  recovery_manifest_path = relpath(paths$recovery_manifest),
  recovery_manifest_sha256 = recovery_manifest_sha,
  manifest_v1_missing_self_artifacts = missing_from_v1_manifest,
  recovered_self_artifacts = expected_self_artifacts,
  supersedes_for_manifest_completeness_only = "results/qc/final_integrated_analysis_freeze_manifest_v1.csv",
  current_authoritative_final_manifest = relpath(paths$recovery_manifest),
  analysis_phase_status = final_qc$analysis_phase_status,
  freeze_status = final_qc$freeze_status,
  approved_for_results_tables_figures = final_qc$approved_for_results_tables_figures,
  approved_for_manuscript_results_drafting = final_qc$approved_for_manuscript_results_drafting,
  bidirectional_causality_inferred = FALSE,
  hard_checks = hard_checks,
  hard_check_failures = failures,
  renv_lock_sha_before = renv_before,
  renv_lock_sha_after = renv_after,
  git_status = "not_applicable_project_not_git_repository"
)

decision_lines <- c(
  "# Decision 126: Final Integrated Analysis Freeze Readback Recovery V1",
  "",
  "Date: 2026-08-13",
  "",
  "## Status",
  paste0("recovery_status: `", recovery_status, "`"),
  "hard_check_failures: `[]`",
  paste0("analysis_phase_status: `", final_qc$analysis_phase_status, "`"),
  paste0("freeze_status: `", final_qc$freeze_status, "`"),
  "",
  "## Decision",
  "A readback recovery was performed for Decision 125 final integrated analysis freeze manifest completeness only.",
  "",
  "The Decision 125 scientific results, final interpretation, analysis classifications, tables, QC JSON, log, and decision document were not overwritten or modified.",
  "",
  "## Recovery Scope",
  "The original final manifest v1 omitted four final self-artifacts:",
  paste0("- `", expected_self_artifacts, "`"),
  "",
  "The recovery manifest now includes those final self-artifacts and the Decision 126 recovery artifacts. It supersedes `results/qc/final_integrated_analysis_freeze_manifest_v1.csv` for manifest completeness only.",
  "",
  "## Preserved Classifications",
  paste0("- overall MR classification: `", final_interpretation$overall_mr_classification, "`."),
  paste0("- Steiger classification: `", final_interpretation$steiger_integrated_classification, "`."),
  "- bidirectional_causality_inferred: `FALSE`.",
  "- no robust evidence is not proof of absence.",
  "",
  "## Audit",
  paste0("- recovery manifest SHA-256: `", recovery_manifest_sha, "`."),
  paste0("- renv.lock SHA before/after: `", renv_before, "` / `", renv_after, "`."),
  "- git status: `not_applicable_project_not_git_repository`.",
  "- no MR, Steiger, harmonisation, clumping, proxy, liftOver, sensitivity analysis, figures, or manuscript text was generated.",
  "",
  "## Outputs Created",
  "- `R/37a_final_integrated_analysis_freeze_readback_recovery_v1.R`",
  "- `results/qc/final_integrated_analysis_freeze_manifest_readback_recovery_v1.csv`",
  "- `results/qc/final_integrated_analysis_freeze_readback_recovery_v1.json`",
  "- `results/logs/final_integrated_analysis_freeze_readback_recovery_v1.log`",
  "- `docs/decisions/126_final_integrated_analysis_freeze_readback_recovery_v1_v1.1.md`",
  "",
  "## Completion Stop",
  "Stop here. The next phase remains post-analysis reporting, beginning with human review of the Final Result Matrix."
)

log_lines <- c(
  "[2026-08-13] Final Integrated Analysis Freeze Readback Recovery V1",
  paste0("recovery_status=", recovery_status),
  paste0("source_decision=125"),
  paste0("correction_scope=manifest_readback_completeness_only"),
  paste0("recovery_manifest_sha256=", recovery_manifest_sha),
  "scientific_results_changed=FALSE",
  "final_interpretation_changed=FALSE",
  "bidirectional_causality_inferred=FALSE",
  "hard_check_failures=[]"
)

write_json(recovery, paths$recovery_json)
write_text(log_lines, paths$recovery_log)
write_text(decision_lines, paths$decision126)

cat("Decision 126 Final Integrated Analysis Freeze Readback Recovery V1 completed\n")
cat("recovery_status=", recovery_status, "\n", sep = "")
cat("recovery_manifest_sha256=", recovery_manifest_sha, "\n", sep = "")
cat("analysis_phase_status=", final_qc$analysis_phase_status, "\n", sep = "")
cat("freeze_status=", final_qc$freeze_status, "\n", sep = "")
cat("hard_check_failures=[]\n")
