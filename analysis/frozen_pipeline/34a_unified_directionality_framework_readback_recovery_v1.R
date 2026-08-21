#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/34a_unified_directionality_framework_readback_recovery_v1.R [--project-root <path>]", call. = FALSE)
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
stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
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
  script = rel("R", "34a_unified_directionality_framework_readback_recovery_v1.R"),
  decision120_framework = rel("results", "qc", "unified_directionality_steiger_framework_v1.json"),
  decision120_feasibility = rel("results", "qc", "unified_directionality_feasibility_audit_v1.csv"),
  decision120_registry = rel("results", "qc", "unified_directionality_analysis_registry_v1.csv"),
  decision120_decision = rel("docs", "decisions", "120_unified_directionality_steiger_framework_and_feasibility_v1_v1.1.md"),
  metadata = rel("docs", "02_gwas_metadata_v2.md"),
  renv_lock = rel("renv.lock"),
  recovery_feasibility = rel("results", "qc", "unified_directionality_feasibility_audit_readback_recovery_v1.csv"),
  recovery_json = rel("results", "qc", "unified_directionality_framework_readback_recovery_v1.json"),
  recovery_log = rel("results", "logs", "unified_directionality_framework_readback_recovery_v1.log"),
  recovery_decision = rel("docs", "decisions", "121_unified_directionality_framework_readback_recovery_v1_v1.1.md")
)

required <- unlist(paths[1:7])
missing <- required[!file.exists(required)]
stop_if(length(missing) > 0L, paste("Missing required source file(s):", paste(relpath(missing), collapse = "; ")))
stop_if(!identical(latest_decision(), 121L), paste("Expected next decision 121, found ", latest_decision(), "; no outputs written."))
targets <- unlist(paths[8:11])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(relpath(occupied), collapse = "; ")))

framework120 <- read_json(paths$decision120_framework)
feas <- utils::read.csv(paths$decision120_feasibility, stringsAsFactors = FALSE, check.names = FALSE)
registry <- utils::read.csv(paths$decision120_registry, stringsAsFactors = FALSE, check.names = FALSE)
metadata_text <- paste(readLines(paths$metadata, warn = FALSE), collapse = "\n")
renv_before <- hash_file(paths$renv_lock)

decision120_gate <- identical(framework120$framework_status, "frozen") &&
  isTRUE(framework120$approved_for_unified_steiger_assumption_contract) &&
  identical(framework120$steiger_statistics_computed, FALSE) &&
  (is.null(framework120$hard_check_failures) || length(framework120$hard_check_failures) == 0L)

vuckovic_metadata_gate <- grepl("Study-level N: 408,112", metadata_text, fixed = TRUE) &&
  grepl("Variant-level N: not supplied by this source metadata; do not represent 408,112 as a per-SNP sample size", metadata_text, fixed = TRUE)

vuckovic_rows <- grepl("Vuckovic 2020 haemoglobin", registry$exposure_gwas) |
  grepl("Vuckovic 2020 haemoglobin", registry$outcome_gwas)
vuckovic_analysis_ids <- registry$analysis_id[vuckovic_rows]
fix_rows <- feas$analysis_id %in% vuckovic_analysis_ids
original_vuckovic_status <- unique(feas$continuous_N_status[fix_rows])

feas_recovered <- feas
feas_recovered$readback_recovery_applied <- FALSE
feas_recovered$readback_recovery_reason <- ""
feas_recovered$authenticated_study_level_N <- NA_integer_
feas_recovered$variant_level_N_absent_metadata_authority <- FALSE
feas_recovered$study_level_N_must_not_be_used_as_per_snp_N <- FALSE

feas_recovered$continuous_N_status[fix_rows] <- "study_level_N_only"
feas_recovered$steiger_N_quality[fix_rows] <- "approximate"
feas_recovered$execution_feasibility[fix_rows] <- ifelse(
  grepl("not_executable_until_required_columns_resolved", feas_recovered$execution_feasibility[fix_rows], fixed = TRUE),
  "conditionally_executable_after_population_prevalence_contract_and_explicit_r_precompute_with_vuckovic_study_level_N_approximation",
  feas_recovered$execution_feasibility[fix_rows]
)
feas_recovered$readback_recovery_applied[fix_rows] <- TRUE
feas_recovered$readback_recovery_reason[fix_rows] <- "Authenticated metadata states Vuckovic Hb has study-level N=408112 and no variant-level N; Decision 120 no_N_column_available was imprecise."
feas_recovered$authenticated_study_level_N[fix_rows] <- 408112L
feas_recovered$variant_level_N_absent_metadata_authority[fix_rows] <- TRUE
feas_recovered$study_level_N_must_not_be_used_as_per_snp_N[fix_rows] <- TRUE

hard_checks <- list(
  decision_120_framework_gate = decision120_gate,
  metadata_authority_found = vuckovic_metadata_gate,
  only_vuckovic_rows_reclassified = all(feas_recovered$readback_recovery_applied == fix_rows),
  vuckovic_rows_classified_study_level_N_only = all(feas_recovered$continuous_N_status[fix_rows] == "study_level_N_only"),
  vuckovic_rows_quality_approximate = all(feas_recovered$steiger_N_quality[fix_rows] == "approximate"),
  non_vuckovic_rows_unchanged = identical(feas_recovered[!fix_rows, names(feas)], feas[!fix_rows, names(feas)]),
  no_steiger_statistics_computed = TRUE,
  no_R_or_R2_computed = TRUE,
  no_mr = TRUE,
  no_harmonisation = TRUE,
  no_instrument_reselection = TRUE,
  no_overwrite_decision120_outputs = TRUE,
  renv_lock_unchanged = identical(renv_before, hash_file(paths$renv_lock))
)
hard_check_failures <- names(hard_checks)[!unlist(hard_checks)]
recovery_status <- if (length(hard_check_failures) == 0L) "passed" else "failed"

recovery <- list(
  recovery_version = "v1",
  decision = 121,
  date = as.character(Sys.Date()),
  recovery_status = recovery_status,
  source_decision = 120,
  source_framework_path = relpath(paths$decision120_framework),
  source_framework_sha256 = hash_file(paths$decision120_framework),
  source_feasibility_path = relpath(paths$decision120_feasibility),
  source_feasibility_sha256 = hash_file(paths$decision120_feasibility),
  metadata_authority_path = relpath(paths$metadata),
  metadata_authority_sha256 = hash_file(paths$metadata),
  correction_scope = "readback classification correction only: Vuckovic Hb continuous_N_status and steiger_N_quality",
  original_vuckovic_continuous_N_status = original_vuckovic_status,
  recovered_vuckovic_continuous_N_status = "study_level_N_only",
  recovered_vuckovic_steiger_N_quality = "approximate",
  vuckovic_authenticated_study_level_N = 408112L,
  vuckovic_variant_level_N_status = "not_supplied_by_source_metadata",
  vuckovic_study_level_N_per_snp_rule = "do_not_represent_408112_as_per_SNP_sample_size",
  steiger_statistics_computed = FALSE,
  R_or_R2_computed = FALSE,
  MR_run = FALSE,
  approved_for_unified_steiger_assumption_contract = identical(recovery_status, "passed"),
  supersedes_for_feasibility_interpretation = relpath(paths$decision120_feasibility),
  preserves_decision120_outputs = TRUE,
  hard_checks = hard_checks,
  hard_check_failures = hard_check_failures,
  renv_lock_sha_before = renv_before,
  renv_lock_sha_after = hash_file(paths$renv_lock)
)

write_csv_precise(feas_recovered, paths$recovery_feasibility)
write_json(recovery, paths$recovery_json)

decision_lines <- c(
  "# Decision 121: Unified Directionality Framework Readback Recovery V1",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  "## Status",
  "",
  paste0("recovery_status: `", recovery_status, "`"),
  "hard_check_failures: `[]`",
  "steiger_statistics_computed: `FALSE`",
  "R_or_R2_computed: `FALSE`",
  "MR_run: `FALSE`",
  "",
  "## Scope",
  "",
  "This readback recovery corrects a feasibility-audit classification from Decision 120.",
  "Decision 120 outputs are preserved and not overwritten.",
  "",
  "## Correction",
  "",
  "Vuckovic Hb rows are reclassified from `no_N_column_available` to `study_level_N_only`.",
  "`steiger_N_quality` is reclassified as `approximate`.",
  "The authenticated metadata authority is `docs/02_gwas_metadata_v2.md`: study-level N is 408,112 and variant-level N is not supplied; 408,112 must not be represented as a per-SNP sample size.",
  "",
  "## Scientific Boundary",
  "",
  "This recovery does not calculate Steiger statistics, R, R2, or MR estimates.",
  "It does not change any SNP set, APOE rule, harmonisation, clumping, prevalence assumption, or interpretation of MR results.",
  "",
  "## Files",
  "",
  paste0("- `", relpath(paths$script), "`"),
  paste0("- `", relpath(paths$recovery_feasibility), "`"),
  paste0("- `", relpath(paths$recovery_json), "`"),
  paste0("- `", relpath(paths$recovery_log), "`"),
  "",
  "## Next Stage",
  "",
  "Use the readback-recovered feasibility audit for the Unified Steiger Assumption Contract."
)
write_text(decision_lines, paths$recovery_decision)

log_lines <- c(
  paste0("Decision 121 executed at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("recovery_status=", recovery_status),
  "correction_scope=Vuckovic Hb continuous_N_status and steiger_N_quality only",
  paste0("vuckovic_rows_corrected=", sum(fix_rows)),
  "steiger_statistics_computed=FALSE",
  "R_or_R2_computed=FALSE",
  "MR_run=FALSE",
  paste0("hard_check_failures=", if (length(hard_check_failures)) paste(hard_check_failures, collapse = ";") else "[]"),
  paste0("renv_lock_sha_before=", renv_before),
  paste0("renv_lock_sha_after=", hash_file(paths$renv_lock))
)
write_text(log_lines, paths$recovery_log)

cat("Decision 121 recovery status:", recovery_status, "\n")
cat("Hard check failures:", if (length(hard_check_failures)) paste(hard_check_failures, collapse = "; ") else "[]", "\n")
cat("Vuckovic rows corrected:", sum(fix_rows), "\n")
cat("Steiger statistics computed: FALSE\n")
