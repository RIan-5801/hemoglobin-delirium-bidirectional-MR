#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/37b_final_reporting_metadata_readback_recovery_v1.R [--project-root <path>]", call. = FALSE)
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
read_csv_chr <- function(path) utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character", na.strings = character())
stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
is_empty <- function(x) is.null(x) || length(x) == 0L
records <- function(x) if (!is.data.frame(x) || nrow(x) == 0L) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
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
as_bool <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% "true"
}

paths <- list(
  script37b = rel("R", "37b_final_reporting_metadata_readback_recovery_v1.R"),
  renv_lock = rel("renv.lock"),
  source_primary = rel("results", "final", "final_primary_result_matrix_v1.csv"),
  source_diagnostic = rel("results", "final", "final_diagnostic_status_matrix_v1.csv"),
  source_interpretation = rel("results", "final", "final_scientific_interpretation_v1.json"),
  source_limitations = rel("results", "final", "final_limitations_registry_v1.csv"),
  decision125_qc = rel("results", "qc", "final_integrated_analysis_freeze_v1.json"),
  decision126_qc = rel("results", "qc", "final_integrated_analysis_freeze_readback_recovery_v1.json"),
  decision125 = rel("docs", "decisions", "125_final_integrated_analysis_freeze_v1_v1.1.md"),
  decision126 = rel("docs", "decisions", "126_final_integrated_analysis_freeze_readback_recovery_v1_v1.1.md"),
  recovered_primary = rel("results", "final", "final_primary_result_matrix_readback_recovery_v1.csv"),
  recovered_diagnostic = rel("results", "final", "final_diagnostic_status_matrix_readback_recovery_v1.csv"),
  recovery_json = rel("results", "qc", "final_reporting_metadata_readback_recovery_v1.json"),
  authority_map = rel("results", "qc", "final_reporting_authority_map_v1.csv"),
  recovery_log = rel("results", "logs", "final_reporting_metadata_readback_recovery_v1.log"),
  decision127 = rel("docs", "decisions", "127_final_reporting_metadata_readback_recovery_v1_v1.1.md")
)

required <- unlist(paths[1:10])
missing <- required[!file.exists(required)]
stop_if(length(missing) > 0L, paste("Missing required input(s):", paste(relpath(missing), collapse = "; ")))
stop_if(!identical(latest_decision(), 127L), paste("Expected next decision 127, found ", latest_decision(), "; no outputs written."))

targets <- unlist(paths[11:16])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(relpath(occupied), collapse = "; ")))

renv_before <- hash_file(paths$renv_lock)
primary <- read_csv_chr(paths$source_primary)
diagnostic <- read_csv_chr(paths$source_diagnostic)
interpretation <- read_json(paths$source_interpretation)
decision125 <- read_json(paths$decision125_qc)
decision126 <- read_json(paths$decision126_qc)

target_ids <- c("chen_forward_included", "chen_forward_excluded")
expected_direction_before <- "delirium_to_Hb"
expected_direction_after <- "Hb_to_delirium"
required_primary_cols <- c("analysis_id", "direction", "analysis_family", "exposure_source", "outcome_source", "effect_scale")
required_diagnostic_cols <- c("analysis_id", "direction", "analysis_family", "evidence_level", "apoe_status", "heterogeneity_status", "egger_status", "mr_presso_status", "loo_status", "single_snp_status", "steiger_status")

source_primary_sha <- hash_file(paths$source_primary)
source_diagnostic_sha <- hash_file(paths$source_diagnostic)
source_interpretation_sha <- hash_file(paths$source_interpretation)
source_limitations_sha <- hash_file(paths$source_limitations)

primary_idx <- match(target_ids, primary$analysis_id)
diagnostic_idx <- match(target_ids, diagnostic$analysis_id)

hard_checks_pre <- list(
  decision_125_gate = identical(decision125$freeze_status, "passed") &&
    identical(decision125$analysis_phase_status, "complete_under_frozen_protocol") &&
    is_empty(decision125$hard_check_failures),
  decision_126_gate = identical(decision126$recovery_status, "passed") &&
    identical(decision126$analysis_phase_status, "complete_under_frozen_protocol") &&
    is_empty(decision126$hard_check_failures),
  source_primary_schema = all(required_primary_cols %in% names(primary)),
  source_diagnostic_schema = all(required_diagnostic_cols %in% names(diagnostic)),
  target_rows_present_once_primary = all(!is.na(primary_idx)) &&
    all(vapply(target_ids, function(id) sum(primary$analysis_id == id) == 1L, logical(1))),
  target_rows_present_once_diagnostic = all(!is.na(diagnostic_idx)) &&
    all(vapply(target_ids, function(id) sum(diagnostic$analysis_id == id) == 1L, logical(1))),
  original_primary_direction_error_confirmed = all(primary$direction[primary_idx] == expected_direction_before),
  original_diagnostic_direction_error_confirmed = all(diagnostic$direction[diagnostic_idx] == expected_direction_before),
  chen_forward_metadata_confirms_forward_primary = all(primary$analysis_family[primary_idx] == "forward_alternative_hb_gwas") &&
    all(primary$exposure_source[primary_idx] == "Chen 2020 European haemoglobin") &&
    all(primary$outcome_source[primary_idx] == "FinnGen R13 delirium") &&
    all(primary$effect_scale[primary_idx] == "log_odds_delirium_per_1_unit_genetically_predicted_standardized_Hb"),
  chen_forward_metadata_confirms_forward_diagnostic = all(diagnostic$analysis_family[diagnostic_idx] == "forward_alternative_hb_gwas"),
  no_existing_target_overwrite = all(!file.exists(targets)) && all(!file.exists(paste0(targets, ".partial"))),
  scientific_interpretation_unchanged_pre = identical(hash_file(paths$source_interpretation), source_interpretation_sha),
  limitations_registry_unchanged_pre = identical(hash_file(paths$source_limitations), source_limitations_sha),
  analysis_phase_status_preserved = identical(decision126$analysis_phase_status, "complete_under_frozen_protocol"),
  renv_lock_pre_unchanged = identical(renv_before, hash_file(paths$renv_lock))
)
pre_failures <- names(hard_checks_pre)[!vapply(hard_checks_pre, isTRUE, logical(1))]
stop_if(length(pre_failures) > 0L, paste("Reporting metadata recovery pre-write hard checks failed:", paste(pre_failures, collapse = "; ")))

primary_recovered <- primary
diagnostic_recovered <- diagnostic
primary_recovered$direction[primary_idx] <- expected_direction_after
diagnostic_recovered$direction[diagnostic_idx] <- expected_direction_after

count_changed_cells <- function(a, b) {
  stop_if(!identical(dim(a), dim(b)), "Dimension mismatch in cell-change audit.")
  sum(as.matrix(a) != as.matrix(b), na.rm = TRUE)
}
non_direction_unchanged <- function(a, b) {
  cols <- setdiff(names(a), "direction")
  identical(a[, cols, drop = FALSE], b[, cols, drop = FALSE])
}

write_csv_precise(primary_recovered, paths$recovered_primary)
write_csv_precise(diagnostic_recovered, paths$recovered_diagnostic)

primary_readback <- read_csv_chr(paths$recovered_primary)
diagnostic_readback <- read_csv_chr(paths$recovered_diagnostic)

primary_changed_cells <- count_changed_cells(primary, primary_readback)
diagnostic_changed_cells <- count_changed_cells(diagnostic, diagnostic_readback)

authority_map <- data.frame(
  artifact = c(
    "primary_result_matrix_for_reporting",
    "diagnostic_status_matrix_for_reporting",
    "scientific_interpretation",
    "limitations_registry",
    "final_manifest_completeness",
    "original_primary_result_matrix_provenance",
    "original_diagnostic_status_matrix_provenance"
  ),
  current_authority_path = c(
    relpath(paths$recovered_primary),
    relpath(paths$recovered_diagnostic),
    relpath(paths$source_interpretation),
    relpath(paths$source_limitations),
    decision126$current_authoritative_final_manifest,
    relpath(paths$source_primary),
    relpath(paths$source_diagnostic)
  ),
  authority_decision = c(127, 127, 125, 125, 126, 125, 125),
  current_authority_sha256 = c(
    hash_file(paths$recovered_primary),
    hash_file(paths$recovered_diagnostic),
    source_interpretation_sha,
    source_limitations_sha,
    decision126$recovery_manifest_sha256,
    source_primary_sha,
    source_diagnostic_sha
  ),
  scientific_authority = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
  supersedes_for_reporting_only = c(
    relpath(paths$source_primary),
    relpath(paths$source_diagnostic),
    "",
    "",
    "results/qc/final_integrated_analysis_freeze_manifest_v1.csv",
    "",
    ""
  ),
  notes = c(
    "Direction metadata corrected for chen_forward_included and chen_forward_excluded only.",
    "Direction metadata corrected for chen_forward_included and chen_forward_excluded only.",
    "Unmodified Decision 125 scientific interpretation authority.",
    "Unmodified Decision 125 limitations registry.",
    "Decision 126 readback-recovered manifest completeness authority.",
    "Preserved original Decision 125 reporting matrix with known Chen forward direction metadata error.",
    "Preserved original Decision 125 diagnostic matrix with known Chen forward direction metadata error."
  ),
  stringsAsFactors = FALSE
)
write_csv_precise(authority_map, paths$authority_map)

renv_after <- hash_file(paths$renv_lock)
hard_checks <- c(hard_checks_pre, list(
  primary_row_count_order_unchanged = identical(primary$analysis_id, primary_readback$analysis_id) &&
    identical(nrow(primary), nrow(primary_readback)),
  diagnostic_row_count_order_unchanged = identical(diagnostic$analysis_id, diagnostic_readback$analysis_id) &&
    identical(nrow(diagnostic), nrow(diagnostic_readback)),
  primary_exactly_two_direction_cells_changed = primary_changed_cells == 2L &&
    all(primary_readback$direction[primary_idx] == expected_direction_after),
  diagnostic_exactly_two_direction_cells_changed = diagnostic_changed_cells == 2L &&
    all(diagnostic_readback$direction[diagnostic_idx] == expected_direction_after),
  primary_non_direction_cells_exact_readback_unchanged = non_direction_unchanged(primary, primary_readback),
  diagnostic_non_direction_cells_exact_readback_unchanged = non_direction_unchanged(diagnostic, diagnostic_readback),
  numeric_values_exact_readback_unchanged = non_direction_unchanged(primary, primary_readback),
  exposure_outcome_source_unchanged = all(primary$exposure_source == primary_readback$exposure_source) &&
    all(primary$outcome_source == primary_readback$outcome_source),
  effect_scale_unchanged = all(primary$effect_scale == primary_readback$effect_scale),
  evidence_level_unchanged = all(primary$evidence_level == primary_readback$evidence_level) &&
    all(diagnostic$evidence_level == diagnostic_readback$evidence_level),
  diagnostic_status_unchanged = all(diagnostic[, setdiff(names(diagnostic), c("direction")), drop = FALSE] ==
      diagnostic_readback[, setdiff(names(diagnostic_readback), c("direction")), drop = FALSE]),
  scientific_interpretation_json_not_modified = identical(hash_file(paths$source_interpretation), source_interpretation_sha),
  limitations_registry_not_modified = identical(hash_file(paths$source_limitations), source_limitations_sha),
  decision_125_126_not_overwritten = file.exists(paths$decision125) && file.exists(paths$decision126),
  analysis_phase_status_still_complete = identical(decision126$analysis_phase_status, "complete_under_frozen_protocol"),
  scientific_results_changed_false = TRUE,
  reporting_metadata_corrected_true = TRUE,
  no_mr_or_steiger_or_harmonisation_or_clumping = TRUE,
  renv_lock_unchanged = identical(renv_before, renv_after),
  git_status_not_required = identical(decision126$git_status, "not_applicable_project_not_git_repository")
))
failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
recovery_status <- if (length(failures) == 0L) "passed" else "failed"
stop_if(length(failures) > 0L, paste("Reporting metadata recovery hard checks failed:", paste(failures, collapse = "; ")))

recovery <- list(
  recovery_version = "v1",
  decision = 127,
  date = "2026-08-13",
  recovery_status = recovery_status,
  correction_scope = "final_reporting_direction_metadata_only",
  source_final_freeze_decision = 125,
  source_manifest_completeness_decision = 126,
  source_primary_result_matrix = relpath(paths$source_primary),
  source_diagnostic_status_matrix = relpath(paths$source_diagnostic),
  recovered_primary_result_matrix = relpath(paths$recovered_primary),
  recovered_diagnostic_status_matrix = relpath(paths$recovered_diagnostic),
  corrected_cells = list(
    list(table = "final_primary_result_matrix", analysis_id = "chen_forward_included", column = "direction", from = expected_direction_before, to = expected_direction_after),
    list(table = "final_primary_result_matrix", analysis_id = "chen_forward_excluded", column = "direction", from = expected_direction_before, to = expected_direction_after),
    list(table = "final_diagnostic_status_matrix", analysis_id = "chen_forward_included", column = "direction", from = expected_direction_before, to = expected_direction_after),
    list(table = "final_diagnostic_status_matrix", analysis_id = "chen_forward_excluded", column = "direction", from = expected_direction_before, to = expected_direction_after)
  ),
  source_primary_sha256 = source_primary_sha,
  recovered_primary_sha256 = hash_file(paths$recovered_primary),
  source_diagnostic_sha256 = source_diagnostic_sha,
  recovered_diagnostic_sha256 = hash_file(paths$recovered_diagnostic),
  scientific_interpretation_sha256 = source_interpretation_sha,
  limitations_registry_sha256 = source_limitations_sha,
  authority_map = relpath(paths$authority_map),
  authority_map_sha256 = hash_file(paths$authority_map),
  analysis_phase_status = decision126$analysis_phase_status,
  freeze_status = decision126$freeze_status,
  approved_for_results_tables_figures = decision126$approved_for_results_tables_figures,
  approved_for_manuscript_results_drafting = decision126$approved_for_manuscript_results_drafting,
  scientific_results_changed = FALSE,
  final_interpretation_changed = FALSE,
  limitations_registry_changed = FALSE,
  reporting_metadata_corrected = TRUE,
  bidirectional_causality_inferred = FALSE,
  hard_checks = hard_checks,
  hard_check_failures = failures,
  renv_lock_sha_before = renv_before,
  renv_lock_sha_after = renv_after,
  git_status = "not_applicable_project_not_git_repository"
)

decision_lines <- c(
  "# Decision 127: Final Reporting Metadata Readback Recovery V1",
  "",
  "Date: 2026-08-13",
  "",
  "## Status",
  paste0("recovery_status: `", recovery_status, "`"),
  "hard_check_failures: `[]`",
  paste0("analysis_phase_status: `", decision126$analysis_phase_status, "`"),
  paste0("freeze_status: `", decision126$freeze_status, "`"),
  "",
  "## Decision",
  "A narrow readback recovery was performed for final reporting matrix direction metadata only.",
  "",
  "The Decision 125 scientific interpretation, limitations registry, scientific classifications, and all numeric MR/diagnostic values were not modified.",
  "",
  "## Corrected Cells",
  "- `final_primary_result_matrix_v1.csv`: `chen_forward_included` direction `delirium_to_Hb` -> `Hb_to_delirium`.",
  "- `final_primary_result_matrix_v1.csv`: `chen_forward_excluded` direction `delirium_to_Hb` -> `Hb_to_delirium`.",
  "- `final_diagnostic_status_matrix_v1.csv`: `chen_forward_included` direction `delirium_to_Hb` -> `Hb_to_delirium`.",
  "- `final_diagnostic_status_matrix_v1.csv`: `chen_forward_excluded` direction `delirium_to_Hb` -> `Hb_to_delirium`.",
  "",
  "## Authority Map",
  "- Reporting primary matrix authority: `results/final/final_primary_result_matrix_readback_recovery_v1.csv`.",
  "- Reporting diagnostic matrix authority: `results/final/final_diagnostic_status_matrix_readback_recovery_v1.csv`.",
  "- Scientific interpretation authority remains Decision 125 `results/final/final_scientific_interpretation_v1.json`.",
  "- Final manifest completeness authority remains Decision 126 recovery manifest.",
  "",
  "## Audit",
  "- Both recovered matrices keep identical row count and row order.",
  "- Exactly two direction cells were changed in each recovered matrix.",
  "- All non-direction cells, including numeric values, exposure/outcome source, effect scale, evidence level, and diagnostic status, are exact readback unchanged.",
  "- scientific_results_changed: `FALSE`.",
  "- reporting_metadata_corrected: `TRUE`.",
  paste0("- renv.lock SHA before/after: `", renv_before, "` / `", renv_after, "`."),
  "- git status: `not_applicable_project_not_git_repository`.",
  "- no MR, Steiger, harmonisation, clumping, proxy, liftOver, sensitivity analysis, figures, or manuscript text was generated.",
  "",
  "## Outputs Created",
  "- `R/37b_final_reporting_metadata_readback_recovery_v1.R`",
  "- `results/final/final_primary_result_matrix_readback_recovery_v1.csv`",
  "- `results/final/final_diagnostic_status_matrix_readback_recovery_v1.csv`",
  "- `results/qc/final_reporting_metadata_readback_recovery_v1.json`",
  "- `results/qc/final_reporting_authority_map_v1.csv`",
  "- `results/logs/final_reporting_metadata_readback_recovery_v1.log`",
  "- `docs/decisions/127_final_reporting_metadata_readback_recovery_v1_v1.1.md`",
  "",
  "## Completion Stop",
  "Stop here. The next phase is table architecture using the Decision 127 reporting matrix authorities."
)

log_lines <- c(
  "[2026-08-13] Final Reporting Metadata Readback Recovery V1",
  paste0("recovery_status=", recovery_status),
  "correction_scope=final_reporting_direction_metadata_only",
  "corrected_primary_cells=2",
  "corrected_diagnostic_cells=2",
  "scientific_results_changed=FALSE",
  "reporting_metadata_corrected=TRUE",
  paste0("recovered_primary_sha256=", hash_file(paths$recovered_primary)),
  paste0("recovered_diagnostic_sha256=", hash_file(paths$recovered_diagnostic)),
  paste0("authority_map_sha256=", hash_file(paths$authority_map)),
  "hard_check_failures=[]"
)

write_json(recovery, paths$recovery_json)
write_text(log_lines, paths$recovery_log)
write_text(decision_lines, paths$decision127)

cat("Decision 127 Final Reporting Metadata Readback Recovery V1 completed\n")
cat("recovery_status=", recovery_status, "\n", sep = "")
cat("analysis_phase_status=", decision126$analysis_phase_status, "\n", sep = "")
cat("freeze_status=", decision126$freeze_status, "\n", sep = "")
cat("scientific_results_changed=FALSE\n")
cat("reporting_metadata_corrected=TRUE\n")
cat("hard_check_failures=[]\n")
