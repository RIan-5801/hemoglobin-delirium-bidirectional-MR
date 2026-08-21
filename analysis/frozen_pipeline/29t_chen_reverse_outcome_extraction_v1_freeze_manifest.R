#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/29t_chen_reverse_outcome_extraction_v1_freeze_manifest.R [--project-root <path>]", call. = FALSE)
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
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE))

paths <- list(
  script = rel("R", "29t_chen_reverse_outcome_extraction_v1_freeze_manifest.R"),
  decision109_qc = rel("results", "qc", "chen_reverse_outcome_extraction_execution_contract_v1_readback_closure.json"),
  target_authority = rel("results", "qc", "chen_reverse_outcome_target_authority_audit_v1_readback_recovery.csv"),
  source_cert = rel("results", "qc", "chen_2020_hb_source_certification_v1.json"),
  chen_raw = rel("data_raw", "gwas", "BCX2_HGB_EA_GWAMA.out.gz"),
  failed_script = rel("R", "29_chen_reverse_outcome_extraction_v1.R"),
  extraction_script = rel("R", "29r_chen_reverse_outcome_extraction_v1_technical_recovery.R"),
  readback_script = rel("R", "29s_chen_reverse_outcome_extraction_v1_readback_recovery.R"),
  extraction_qc = rel("results", "qc", "chen_reverse_outcome_extraction_v1.json"),
  readback_qc = rel("results", "qc", "chen_reverse_outcome_extraction_v1_readback_recovery.json"),
  readback_mismatch = rel("results", "qc", "chen_reverse_outcome_extraction_v1_readback_recovery_mismatch_audit.csv"),
  initial_failed_log = rel("results", "logs", "chen_reverse_outcome_extraction_v1.log"),
  extraction_recovery_log = rel("results", "logs", "chen_reverse_outcome_extraction_v1_technical_recovery.log"),
  readback_log = rel("results", "logs", "chen_reverse_outcome_extraction_v1_readback_recovery.log"),
  decision110 = rel("docs", "decisions", "110_chen_reverse_outcome_extraction_v1_v1.1.md"),
  decision111 = rel("docs", "decisions", "111_chen_reverse_outcome_extraction_v1_readback_recovery_v1.1.md"),
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
  manifest_csv = rel("results", "qc", "chen_reverse_outcome_extraction_v1_freeze_manifest.csv"),
  freeze_json = rel("results", "qc", "chen_reverse_outcome_extraction_v1_freeze.json"),
  log = rel("results", "logs", "chen_reverse_outcome_extraction_v1_freeze.log"),
  decision = rel("docs", "decisions", "112_chen_reverse_outcome_extraction_v1_freeze_v1.1.md"),
  renv_lock = rel("renv.lock")
)

outputs <- unlist(paths[c("manifest_csv", "freeze_json", "log", "decision")])
occupied <- outputs[file.exists(outputs) | file.exists(paste0(outputs, ".partial"))]
stop_if(length(occupied) > 0L, paste("Output or partial exists:", paste(occupied, collapse = "; ")))

inputs <- unlist(paths[c(
  "decision109_qc", "target_authority", "source_cert", "chen_raw", "failed_script", "extraction_script",
  "readback_script", "extraction_qc", "readback_qc", "readback_mismatch", "initial_failed_log",
  "extraction_recovery_log", "readback_log", "decision110", "decision111", "union_targets",
  "master_parquet", "master_tsv", "strict_inc_parquet", "strict_inc_tsv", "strict_exc_parquet",
  "strict_exc_tsv", "relaxed_inc_parquet", "relaxed_inc_tsv", "relaxed_exc_parquet",
  "relaxed_exc_tsv", "match_audit", "missing", "renv_lock"
)])
missing_inputs <- inputs[!file.exists(inputs)]
stop_if(length(missing_inputs) > 0L, paste("Required input missing:", paste(missing_inputs, collapse = "; ")))

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 112L), paste0("Expected next decision 112, found ", next_decision, "; no outputs written."))

dir.create(dirname(paths$log), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = paths$log, append = TRUE)
log_line("stage=chen_reverse_outcome_extraction_v1_freeze_start")

contract <- read_json(paths$decision109_qc)
source_cert <- read_json(paths$source_cert)
extraction_qc <- read_json(paths$extraction_qc)
readback_qc <- read_json(paths$readback_qc)
target_authority <- utils::read.csv(paths$target_authority, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
master <- utils::read.delim(paths$master_tsv, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
match_audit <- utils::read.csv(paths$match_audit, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
missing_tsv <- utils::read.delim(paths$missing, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
mismatch_audit <- utils::read.csv(paths$readback_mismatch, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")

branch_counts <- data.frame(
  branch = c("strict_apoe_included", "strict_apoe_excluded", "relaxed_apoe_included", "relaxed_apoe_excluded"),
  expected_targets = c(1L, 1L, 11L, 10L),
  exact_targets = c(
    sum(master$strict_included_member == "TRUE" & master$match_status == "unique_exact_match"),
    sum(master$strict_excluded_member == "TRUE" & master$match_status == "unique_exact_match"),
    sum(master$relaxed_included_member == "TRUE" & master$match_status == "unique_exact_match"),
    sum(master$relaxed_excluded_member == "TRUE" & master$match_status == "unique_exact_match")
  ),
  ready_for_harmonisation_preflight = TRUE,
  stringsAsFactors = FALSE
)
branch_counts$ready_for_harmonisation_preflight <- branch_counts$exact_targets > 0L

manifest_items <- data.frame(
  category = c(
    "script", "script", "script", "script",
    "authority", "authority", "authority", "raw_source",
    "qc", "qc", "qc", "qc", "qc",
    "log", "log", "log",
    "decision", "decision", "decision",
    "data", "data", "data", "data", "data", "data", "data", "data", "data", "data", "data", "data",
    "environment"
  ),
  role = c(
    "initial_failed_extraction_script", "technical_recovery_extraction_script", "readback_recovery_script", "freeze_manifest_script",
    "decision_109_contract", "corrected_target_authority", "chen_source_certification", "chen_raw_source_read_only",
    "decision_110_extraction_qc_failed_preserved", "decision_111_readback_recovery_qc", "decision_111_mismatch_audit", "match_audit", "missing_audit",
    "initial_failed_scan_log", "technical_recovery_extraction_log", "readback_recovery_log",
    "decision_110_extraction", "decision_111_readback_recovery", "decision_112_freeze",
    "union_targets", "master_parquet", "master_tsv", "strict_included_parquet", "strict_included_tsv", "strict_excluded_parquet", "strict_excluded_tsv",
    "relaxed_included_parquet", "relaxed_included_tsv", "relaxed_excluded_parquet", "relaxed_excluded_tsv", "freeze_manifest_csv",
    "renv_lock"
  ),
  path = c(
    paths$failed_script, paths$extraction_script, paths$readback_script, paths$script,
    paths$decision109_qc, paths$target_authority, paths$source_cert, paths$chen_raw,
    paths$extraction_qc, paths$readback_qc, paths$readback_mismatch, paths$match_audit, paths$missing,
    paths$initial_failed_log, paths$extraction_recovery_log, paths$readback_log,
    paths$decision110, paths$decision111, paths$decision,
    paths$union_targets, paths$master_parquet, paths$master_tsv, paths$strict_inc_parquet, paths$strict_inc_tsv, paths$strict_exc_parquet, paths$strict_exc_tsv,
    paths$relaxed_inc_parquet, paths$relaxed_inc_tsv, paths$relaxed_exc_parquet, paths$relaxed_exc_tsv, paths$manifest_csv,
    paths$renv_lock
  ),
  stringsAsFactors = FALSE
)
manifest_items$relative_path <- vapply(manifest_items$path, relpath, character(1))
manifest_items$exists <- file.exists(manifest_items$path)
manifest_items$size_bytes <- unname(file.info(manifest_items$path)$size)
manifest_items$sha256 <- vapply(seq_len(nrow(manifest_items)), function(i) {
  if (manifest_items$exists[[i]] && !identical(manifest_items$relative_path[[i]], relpath(paths$manifest_csv)) && !identical(manifest_items$relative_path[[i]], relpath(paths$decision))) {
    hash_file(manifest_items$path[[i]])
  } else {
    NA_character_
  }
}, character(1))

write_csv(manifest_items, paths$manifest_csv)
manifest_sha256 <- hash_file(paths$manifest_csv)

hard_checks <- list(
  decision_109_contract_frozen = identical(contract$authoritative_contract_status, "frozen") && isTRUE(contract$approved_for_chen_reverse_outcome_extraction),
  strict_exposure_authority_preserved = identical(contract$strict_exposure_authority$freeze_status, "passed") && identical(contract$strict_exposure_authority$included_nsnp, 1L) && identical(contract$strict_exposure_authority$excluded_nsnp, 1L),
  relaxed_exposure_authority_preserved = identical(contract$relaxed_exposure_authority$freeze_status, "passed") && identical(contract$relaxed_exposure_authority$included_nsnp, 11L) && identical(contract$relaxed_exposure_authority$excluded_nsnp, 10L),
  corrected_target_audit_used = identical(relpath(paths$target_authority), contract$corrected_target_audit_authority),
  target_counts_preserved = nrow(target_authority) == 12L && nrow(master) == 12L,
  extraction_decision_110_preserved_as_failed = identical(extraction_qc$outcome_extraction_status, "failed"),
  extraction_failure_is_readback_only = setequal(extraction_qc$hard_check_failures, c("master_parquet_tsv_consistency", "branch_parquet_tsv_consistency")),
  readback_recovery_decision_111_passed = identical(readback_qc$recovery_status, "passed") && isTRUE(readback_qc$approved_for_freeze),
  no_true_readback_mismatches = nrow(mismatch_audit) == 0L,
  source_rows_scanned_once_in_successful_extraction = identical(extraction_qc$source_rows_scanned, 50638016L),
  chen_source_sha_matches_certification = identical(tolower(hash_file(paths$chen_raw)), tolower(source_cert$source_sha256)) && identical(extraction_qc$source_sha_before, extraction_qc$source_sha_after),
  all_targets_unique_exact = all(master$match_status == "unique_exact_match") && all(as.integer(match_audit$raw_match_count) == 1L),
  no_missing_targets = nrow(missing_tsv) == 0L,
  no_multiple_match_file_created = !file.exists(rel("results", "qc", "chen_reverse_outcome_multiple_matches_v1.tsv")),
  no_incompatible_targets = identical(extraction_qc$union_marker_effect_allele_incompatible_count, 0L),
  branch_counts_preserved = all(branch_counts$expected_targets == branch_counts$exact_targets),
  outcome_scale_preserved = identical(extraction_qc$outcome_scale, "standardized_quantitative_Hb_effect"),
  exact_marker_id_matching_only = identical(extraction_qc$matching_method, "exact_BCX2_marker_id_string_identity"),
  no_proxy = isFALSE(extraction_qc$proxy_used),
  no_liftover = isFALSE(extraction_qc$liftover_used),
  no_nearest_variant = isFALSE(extraction_qc$nearest_variant_used),
  no_strand_complement_identity_rescue = isFALSE(extraction_qc$strand_complement_identity_rescue_used),
  no_exposure_reselection = isTRUE(extraction_qc$hard_checks$no_exposure_reselection),
  no_reclumping = isTRUE(extraction_qc$hard_checks$no_reclumping),
  no_harmonisation = isFALSE(extraction_qc$harmonisation_performed) && isFALSE(readback_qc$harmonisation_performed),
  no_mr = isFALSE(extraction_qc$mr_run) && isFALSE(readback_qc$mr_run),
  no_steiger = isFALSE(extraction_qc$steiger_run) && isFALSE(readback_qc$steiger_run),
  all_manifest_items_exist_except_self_hash_placeholders = all(manifest_items$exists),
  manifest_sha256_created = nzchar(manifest_sha256),
  renv_lock_unchanged = identical(hash_file(paths$renv_lock), "253471c02e1e47a40d0f68b296d4ae2b1df471d757ac86328b3e974018d039f3")
)
hard_check_failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
freeze_status <- if (length(hard_check_failures) == 0L) "passed" else "failed"

freeze <- list(
  freeze_version = "v1",
  decision = 112,
  date = format(Sys.Date()),
  freeze_status = freeze_status,
  extraction_decision = 110,
  readback_recovery_decision = 111,
  approved_for_chen_reverse_harmonisation_preflight = identical(freeze_status, "passed"),
  analysis_direction = "delirium_to_Hb",
  analysis_role = "reverse_alternative_hb_outcome_sensitivity",
  independent_replication = FALSE,
  outcome_source = "Chen_2020_Hb_BCX2",
  outcome_scale = "standardized_quantitative_Hb_effect",
  source_rows_scanned = extraction_qc$source_rows_scanned,
  source_sha_before = extraction_qc$source_sha_before,
  source_sha_after = extraction_qc$source_sha_after,
  union_target_count = nrow(master),
  union_unique_exact_match_count = extraction_qc$union_unique_exact_match_count,
  union_missing_count = extraction_qc$union_missing_count,
  union_multiple_exact_match_count = extraction_qc$union_multiple_exact_match_count,
  union_marker_effect_allele_incompatible_count = extraction_qc$union_marker_effect_allele_incompatible_count,
  branch_counts = records(branch_counts),
  matching_method = "exact_BCX2_marker_id_string_identity",
  proxy_used = FALSE,
  liftover_used = FALSE,
  nearest_variant_used = FALSE,
  strand_complement_identity_rescue_used = FALSE,
  harmonisation_performed = FALSE,
  mr_run = FALSE,
  steiger_run = FALSE,
  manifest_csv = relpath(paths$manifest_csv),
  manifest_sha256 = manifest_sha256,
  hard_checks = hard_checks,
  hard_check_failures = hard_check_failures,
  informational_findings = list(
    decision_110_status_preserved = extraction_qc$outcome_extraction_status,
    decision_110_failures_preserved = extraction_qc$hard_check_failures,
    decision_111_recovery_status = readback_qc$recovery_status,
    decision_111_recovery_scope = readback_qc$recovery_scope,
    initial_failed_log_preserved = relpath(paths$initial_failed_log),
    extraction_recovery_log = relpath(paths$extraction_recovery_log),
    readback_recovery_log = relpath(paths$readback_log)
  )
)
write_json(freeze, paths$freeze_json)

decision_lines <- c(
  "# Decision 112: Chen Reverse Outcome Extraction V1 Freeze",
  "",
  paste0("Date: ", format(Sys.Date())),
  "",
  "## Status",
  paste0("freeze_status: `", freeze_status, "`"),
  paste0("approved_for_chen_reverse_harmonisation_preflight: `", identical(freeze_status, "passed"), "`"),
  "",
  "## Decision",
  "Freeze the Chen reverse outcome extraction V1 artifacts after Decision 111 readback recovery classified the Decision 110 Parquet/TSV failure as a missing-value encoding artifact.",
  "",
  "This freeze locks targeted Chen outcome data for reverse sensitivity harmonisation preflight only. It does not perform harmonisation, MR, Steiger, proxy search, liftOver, nearest-variant matching, reclumping, or exposure reselection.",
  "",
  "## Frozen Evidence",
  paste0("- Manifest: `", relpath(paths$manifest_csv), "`."),
  paste0("- Manifest SHA-256: `", manifest_sha256, "`."),
  paste0("- Chen source SHA before/after extraction: `", extraction_qc$source_sha_before, "` / `", extraction_qc$source_sha_after, "`."),
  paste0("- Source rows scanned: `", extraction_qc$source_rows_scanned, "`."),
  paste0("- Union exact / missing / multiple / incompatible: `", extraction_qc$union_unique_exact_match_count, " / ", extraction_qc$union_missing_count, " / ", extraction_qc$union_multiple_exact_match_count, " / ", extraction_qc$union_marker_effect_allele_incompatible_count, "`."),
  paste0("- Hard-check failures: `", if (length(hard_check_failures) == 0L) "none" else paste(hard_check_failures, collapse = ";"), "`."),
  "",
  "## Branch Readiness",
  paste0("- strict_apoe_included: `", branch_counts$exact_targets[branch_counts$branch == "strict_apoe_included"], "/", branch_counts$expected_targets[branch_counts$branch == "strict_apoe_included"], "` exact targets."),
  paste0("- strict_apoe_excluded: `", branch_counts$exact_targets[branch_counts$branch == "strict_apoe_excluded"], "/", branch_counts$expected_targets[branch_counts$branch == "strict_apoe_excluded"], "` exact targets."),
  paste0("- relaxed_apoe_included: `", branch_counts$exact_targets[branch_counts$branch == "relaxed_apoe_included"], "/", branch_counts$expected_targets[branch_counts$branch == "relaxed_apoe_included"], "` exact targets."),
  paste0("- relaxed_apoe_excluded: `", branch_counts$exact_targets[branch_counts$branch == "relaxed_apoe_excluded"], "/", branch_counts$expected_targets[branch_counts$branch == "relaxed_apoe_excluded"], "` exact targets."),
  "",
  "## Outputs",
  paste0("- `", relpath(paths$freeze_json), "`"),
  paste0("- `", relpath(paths$manifest_csv), "`"),
  paste0("- `", relpath(paths$log), "`")
)
write_text(decision_lines, paths$decision)

log_line("freeze_status=", freeze_status, "; hard_check_failures=", paste(hard_check_failures, collapse = ","))
cat("freeze_status=", freeze_status, "\n", sep = "")
cat("approved_for_chen_reverse_harmonisation_preflight=", identical(freeze_status, "passed"), "\n", sep = "")
cat("hard_check_failures=", paste(hard_check_failures, collapse = ","), "\n", sep = "")
cat("manifest_sha256=", manifest_sha256, "\n", sep = "")
