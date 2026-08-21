#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/27a_chen_reverse_outcome_extraction_preflight_plan_v1_readback_recovery.R [--project-root <path>]", call. = FALSE)
}
setwd(root)

for (pkg in c("jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

rel <- function(...) file.path(root, ...)
norm <- function(x) gsub("\\\\", "/", x)
relpath <- function(x) {
  norm(sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", norm(normalizePath(x, winslash = "/", mustWork = FALSE))))
}
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
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

paths <- list(
  script = rel("R", "27a_chen_reverse_outcome_extraction_preflight_plan_v1_readback_recovery.R"),
  failed_script = rel("R", "27_chen_reverse_outcome_extraction_preflight_plan_v1.R"),
  failed_qc = rel("results", "qc", "chen_reverse_outcome_extraction_preflight_plan_v1.json"),
  failed_log = rel("results", "logs", "chen_reverse_outcome_extraction_preflight_plan_v1.log"),
  failed_decision = rel("docs", "decisions", "105_chen_reverse_outcome_extraction_preflight_plan_v1_v1.1.md"),
  decision104 = rel("docs", "decisions", "104_chen_reverse_sensitivity_design_v1_v1.1.md"),
  design104_qc = rel("results", "qc", "chen_reverse_sensitivity_design_v1.json"),
  decision78 = rel("docs", "decisions", "78_chen_2020_hb_source_certification_v1_v1.1.md"),
  decision83 = rel("docs", "decisions", "83_chen_2020_hb_official_source_dictionary_audit_v1_v1.1.md"),
  decision84 = rel("docs", "decisions", "84_chen_identifier_resolution_amendment_v1_v1.1.md"),
  chen_source_cert = rel("results", "qc", "chen_2020_hb_source_certification_v1.json"),
  chen_readme = rel("docs", "source_metadata", "readme_BCX2_meta_analyses.txt"),
  reverse_primary_instrument_freeze = rel("results", "qc", "finngen_r13_delirium_reverse_primary_instruments_v4_freeze_v3.json"),
  reverse_strict_mr_freeze = rel("results", "qc", "reverse_strict_primary_mr_v1_freeze.json"),
  reverse_relaxed_mr_freeze = rel("results", "qc", "reverse_relaxed_mr_v1_freeze.json"),
  relaxed_harmonisation_freeze = rel("results", "qc", "vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_freeze.json"),
  reverse_relaxed_palindrome_rule = rel("results", "qc", "reverse_relaxed_palindromic_handling_rule_v1.json"),
  renv_lock = rel("renv.lock"),
  recovery_qc = rel("results", "qc", "chen_reverse_outcome_extraction_preflight_plan_v1_readback_recovery.json"),
  recovery_log = rel("results", "logs", "chen_reverse_outcome_extraction_preflight_plan_v1_readback_recovery.log"),
  decision106 = rel("docs", "decisions", "106_chen_reverse_outcome_extraction_preflight_plan_v1_readback_recovery_v1.1.md")
)

required_inputs <- unlist(paths[c(
  "script", "failed_script", "failed_qc", "failed_log", "failed_decision",
  "decision104", "design104_qc", "decision78", "decision83", "decision84",
  "chen_source_cert", "chen_readme", "reverse_primary_instrument_freeze",
  "reverse_strict_mr_freeze", "reverse_relaxed_mr_freeze",
  "relaxed_harmonisation_freeze", "reverse_relaxed_palindrome_rule", "renv_lock"
)])
missing <- required_inputs[!file.exists(required_inputs)]
if (length(missing) > 0L) stop("Missing required input(s): ", paste(missing, collapse = "; "), call. = FALSE)

targets <- unlist(paths[c("recovery_qc", "recovery_log", "decision106")])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(occupied, collapse = "; ")))

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 106L), paste0("Expected next decision 106, found ", next_decision, "; no outputs written."))

text_file <- function(path) paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
failed_qc <- read_json(paths$failed_qc)
design104 <- read_json(paths$design104_qc)
chen_cert <- read_json(paths$chen_source_cert)
reverse_primary <- read_json(paths$reverse_primary_instrument_freeze)
reverse_strict <- read_json(paths$reverse_strict_mr_freeze)
reverse_relaxed <- read_json(paths$reverse_relaxed_mr_freeze)
relaxed_harmonisation_freeze <- read_json(paths$relaxed_harmonisation_freeze)
readme_text <- text_file(paths$chen_readme)
decision84_text <- text_file(paths$decision84)

known_false_positive_failures <- c(
  "chen_source_certification_passed",
  "documented_vs_inferred_scale_separated",
  "outcome_extraction_performed"
)

failure_set_expected <- setequal(unlist(failed_qc$hard_check_failures), known_false_positive_failures)

forbidden_units <- unlist(chen_cert$phenotype$forbidden_interpretations, use.names = FALSE)
contains_all <- function(text, terms) all(vapply(terms, grepl, logical(1), x = text, fixed = TRUE))

corrected_hard_checks <- list(
  failed_decision_105_preserved = identical(failed_qc$plan_status, "failed") && failure_set_expected,
  failure_modes_are_technical_false_positives = failure_set_expected,
  decision_104_design_gate = identical(design104$design_status, "passed") &&
    isTRUE(design104$approved_for_chen_reverse_sensitivity_outcome_extraction_plan),
  chen_source_certification_passed = identical(chen_cert$certification_status, "passed") &&
    identical(as.integer(chen_cert$decision), 78L),
  chen_official_dictionary_present = contains_all(readme_text, c(
    "reference_allele - Effect allele",
    "other_allele - Non effect allele",
    "eaf - effect allele frequency",
    "beta - Overall effect size/beta value for meta-analysis",
    "n_samples - Number of samples with marker present"
  )),
  chen_other_allele_reverified = identical(chen_cert$allele_convention$other_allele_field, "other_allele") &&
    grepl("other_allele - Non effect allele", readme_text, fixed = TRUE),
  chen_beta_scale_reverified = identical(chen_cert$phenotype$effect_scale, "standardized_quantitative_Hb_effect"),
  documented_vs_inferred_scale_separated = isFALSE(chen_cert$phenotype$physical_unit_claim_allowed) &&
    setequal(forbidden_units, c("g/dL", "g/L", "clinical haemoglobin concentration units")),
  strict_reverse_instruments_frozen_and_reused = identical(reverse_primary$freeze_status, "passed") &&
    isTRUE(reverse_primary$approved_for_reverse_primary_outcome_extraction),
  strict_reverse_mr_frozen = identical(reverse_strict$freeze_status, "passed") &&
    identical(reverse_strict$instrument_threshold, 5e-08),
  relaxed_reverse_instruments_frozen_and_reused = identical(relaxed_harmonisation_freeze$freeze_status, "passed"),
  relaxed_reverse_exploratory_frozen = identical(reverse_relaxed$freeze_status, "passed") &&
    identical(reverse_relaxed$p_threshold, 5e-06) &&
    identical(reverse_relaxed$branch_type, "protocol_triggered_exploratory_fallback"),
  no_exposure_reselection = TRUE,
  cross_build_identity_plan_defined = !is.null(failed_qc$identity_bridge_plan),
  no_direct_cross_build_coordinate_matching = TRUE,
  decision_84_scope_checked = grepl("It does not apply to Vuckovic, FinnGen, reverse analyses", decision84_text, fixed = TRUE),
  decision_84_reverse_adoption_is_new_decision = TRUE,
  no_proxy = TRUE,
  no_liftover = TRUE,
  no_nearest_variant = TRUE,
  no_fuzzy_matching = TRUE,
  no_strand_complement_rescue = TRUE,
  missing_handling_defined = identical(failed_qc$missing_multiple_incompatible_handling$missing$category, "missing"),
  multiple_handling_defined = identical(failed_qc$missing_multiple_incompatible_handling$multiple$category, "multiple_exact_match"),
  incompatible_handling_defined = identical(failed_qc$missing_multiple_incompatible_handling$incompatible$category, "marker_effect_allele_incompatible"),
  no_chen_raw_scan = isFALSE(failed_qc$prohibited_now$chen_raw_scan_performed),
  outcome_extraction_not_performed = isFALSE(failed_qc$prohibited_now$outcome_extraction_performed),
  no_harmonisation = isFALSE(failed_qc$prohibited_now$harmonisation_performed),
  no_mr = isFALSE(failed_qc$prohibited_now$mr_run),
  no_steiger = isFALSE(failed_qc$prohibited_now$steiger_run)
)

corrected_failures <- names(corrected_hard_checks)[!vapply(corrected_hard_checks, isTRUE, logical(1))]
recovery_status <- if (length(corrected_failures) == 0L) "passed" else "failed"
corrected_plan_status <- if (identical(recovery_status, "passed")) "frozen" else "failed"

manifest_inputs <- unlist(paths[c(
  "script", "failed_script", "failed_qc", "failed_log", "failed_decision",
  "decision104", "design104_qc", "decision78", "decision83", "decision84",
  "chen_source_cert", "chen_readme", "reverse_primary_instrument_freeze",
  "reverse_strict_mr_freeze", "reverse_relaxed_mr_freeze",
  "relaxed_harmonisation_freeze", "reverse_relaxed_palindrome_rule", "renv_lock"
)])
input_manifest <- lapply(names(manifest_inputs), function(nm) {
  list(
    file_role = nm,
    relative_path = relpath(manifest_inputs[[nm]]),
    sha256 = hash_file(manifest_inputs[[nm]]),
    file_size_bytes = unname(file.info(manifest_inputs[[nm]])$size)
  )
})

recovery <- list(
  recovery_version = "v1",
  decision = 106,
  date = format(Sys.Date()),
  recovery_status = recovery_status,
  source_plan_decision = 105,
  source_plan_status = failed_qc$plan_status,
  source_hard_check_failures = failed_qc$hard_check_failures,
  false_positive_classification = list(
    chen_source_certification_passed = "integer/numeric strict-comparison artifact; source certification status is passed and decision is 78",
    documented_vs_inferred_scale_separated = "JSON list/vector strict-comparison artifact; physical unit claims remain forbidden",
    outcome_extraction_performed = "boolean polarity artifact; source field confirms outcome extraction was not performed"
  ),
  corrected_plan_status = corrected_plan_status,
  analysis_direction = failed_qc$analysis_direction,
  analysis_role = failed_qc$analysis_role,
  independent_replication = failed_qc$independent_replication,
  chen_outcome_semantics_reverification = failed_qc$chen_outcome_semantics_reverification,
  exposure_instrument_authority_plan = failed_qc$exposure_instrument_authority_plan,
  identity_bridge_plan = failed_qc$identity_bridge_plan,
  decision_84_scope_reuse_status = failed_qc$decision_84_scope_reuse_status,
  missing_multiple_incompatible_handling = failed_qc$missing_multiple_incompatible_handling,
  prohibited_now = failed_qc$prohibited_now,
  corrected_hard_checks = corrected_hard_checks,
  corrected_hard_check_failures = corrected_failures,
  input_manifest = input_manifest,
  raw_file_hash_recomputation = "not_recomputed_in_readback_recovery; no Chen raw GWAS scan was performed",
  approved_for_chen_reverse_outcome_extraction_execution = identical(corrected_plan_status, "frozen")
)

write_json(recovery, paths$recovery_qc)

decision_lines <- c(
  "# Decision 106: Chen Reverse Outcome Extraction Preflight Plan V1 Readback Recovery",
  "",
  paste0("Date: ", format(Sys.Date())),
  "",
  "## Status",
  paste0("recovery_status: `", recovery_status, "`"),
  paste0("corrected_plan_status: `", corrected_plan_status, "`"),
  paste0("approved_for_chen_reverse_outcome_extraction_execution: `", identical(corrected_plan_status, "frozen"), "`"),
  "",
  "## Decision",
  "Classify the Decision 105 hard-check failures as technical false positives after readback against frozen source authorities.",
  "",
  "Decision 105 is preserved as failed historical evidence. Decision 106 records the corrected readback status for the same plan without rerunning or overwriting Decision 105 outputs.",
  "",
  "## False-Positive Findings",
  "- `chen_source_certification_passed`: integer/numeric strict-comparison artifact; Chen source certification is passed and decision is 78.",
  "- `documented_vs_inferred_scale_separated`: JSON list/vector comparison artifact; g/dL, g/L, and clinical-unit claims remain forbidden.",
  "- `outcome_extraction_performed`: boolean polarity artifact; outcome extraction was not performed.",
  "",
  "## Preserved Scientific Plan",
  "- Chen reverse remains `FinnGen R13 delirium -> Chen 2020 haemoglobin`.",
  "- Role remains reverse alternative-Hb outcome sensitivity, not reverse primary and not independent replication.",
  "- Strict and relaxed reverse exposure authorities must be reused; no exposure reselection or reclumping is allowed.",
  "- Identity matching remains exact GRCh37 coordinate plus unordered allele-set after canonical-rsID bridge.",
  "- Proxy, liftOver, nearest-variant, fuzzy matching, and strand-complement rescue remain forbidden.",
  "",
  "## Safeguards",
  "This readback recovery does not scan Chen raw GWAS, extract outcomes, harmonise, run MR, run Steiger, perform proxy search, perform liftOver, reclump, or modify existing strict/relaxed reverse exposure sets.",
  "",
  "## Corrected Hard Check Failures",
  if (length(corrected_failures) == 0L) "- none" else paste0("- ", corrected_failures),
  "",
  "## Outputs",
  paste0("- `", relpath(paths$recovery_qc), "`"),
  paste0("- `", relpath(paths$recovery_log), "`"),
  paste0("- `", relpath(paths$decision106), "`")
)
write_text(decision_lines, paths$decision106)

log_lines <- c(
  paste0("[", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "] stage=chen_reverse_outcome_extraction_preflight_plan_v1_readback_recovery_start"),
  "decision=106",
  paste0("source_plan_decision=105"),
  paste0("source_plan_status=", failed_qc$plan_status),
  paste0("source_hard_check_failures=", paste(unlist(failed_qc$hard_check_failures), collapse = ",")),
  paste0("recovery_status=", recovery_status),
  paste0("corrected_plan_status=", corrected_plan_status),
  paste0("corrected_hard_check_failures=", paste(corrected_failures, collapse = ",")),
  paste0("approved_for_chen_reverse_outcome_extraction_execution=", identical(corrected_plan_status, "frozen")),
  "chen_raw_scan_performed=FALSE",
  "outcome_extraction_performed=FALSE",
  "harmonisation_performed=FALSE",
  "mr_run=FALSE",
  "steiger_run=FALSE"
)
write_text(log_lines, paths$recovery_log)

cat("recovery_status=", recovery_status, "\n", sep = "")
cat("corrected_plan_status=", corrected_plan_status, "\n", sep = "")
cat("corrected_hard_check_failures=", paste(corrected_failures, collapse = ","), "\n", sep = "")
cat("approved_for_chen_reverse_outcome_extraction_execution=", identical(corrected_plan_status, "frozen"), "\n", sep = "")
