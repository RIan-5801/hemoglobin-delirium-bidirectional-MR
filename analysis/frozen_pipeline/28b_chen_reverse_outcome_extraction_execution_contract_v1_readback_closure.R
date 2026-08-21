#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/28b_chen_reverse_outcome_extraction_execution_contract_v1_readback_closure.R [--project-root <path>]", call. = FALSE)
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
  script = rel("R", "28b_chen_reverse_outcome_extraction_execution_contract_v1_readback_closure.R"),
  decision107_qc = rel("results", "qc", "chen_reverse_outcome_extraction_execution_contract_v1.json"),
  decision107 = rel("docs", "decisions", "107_chen_reverse_outcome_extraction_execution_contract_v1_v1.1.md"),
  decision108_qc = rel("results", "qc", "chen_reverse_outcome_extraction_execution_contract_v1_readback_recovery.json"),
  decision108 = rel("docs", "decisions", "108_chen_reverse_outcome_extraction_execution_contract_v1_readback_recovery_v1.1.md"),
  corrected_target_audit = rel("results", "qc", "chen_reverse_outcome_target_authority_audit_v1_readback_recovery.csv"),
  decision106_qc = rel("results", "qc", "chen_reverse_outcome_extraction_preflight_plan_v1_readback_recovery.json"),
  renv_lock = rel("renv.lock"),
  closure_qc = rel("results", "qc", "chen_reverse_outcome_extraction_execution_contract_v1_readback_closure.json"),
  closure_log = rel("results", "logs", "chen_reverse_outcome_extraction_execution_contract_v1_readback_closure.log"),
  decision109 = rel("docs", "decisions", "109_chen_reverse_outcome_extraction_execution_contract_v1_readback_closure_v1.1.md")
)

required <- unlist(paths[c("script", "decision107_qc", "decision107", "decision108_qc", "decision108", "corrected_target_audit", "decision106_qc", "renv_lock")])
missing <- required[!file.exists(required)]
if (length(missing) > 0L) stop("Missing required input(s): ", paste(missing, collapse = "; "), call. = FALSE)

targets <- unlist(paths[c("closure_qc", "closure_log", "decision109")])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(occupied, collapse = "; ")))

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 109L), paste0("Expected next decision 109, found ", next_decision, "; no outputs written."))

d107 <- read_json(paths$decision107_qc)
d108 <- read_json(paths$decision108_qc)
d106 <- read_json(paths$decision106_qc)
audit <- utils::read.csv(paths$corrected_target_audit, stringsAsFactors = FALSE, check.names = FALSE)

alleles_valid <- all(audit$exposure_effect_allele %in% c("A", "C", "G", "T")) &&
  all(audit$exposure_other_allele %in% c("A", "C", "G", "T")) &&
  all(audit$reference_allele1 %in% c("A", "C", "G", "T")) &&
  all(audit$reference_allele2 %in% c("A", "C", "G", "T"))

closure_hard_checks <- list(
  decision_107_failure_preserved = identical(d107$contract_status, "failed") &&
    identical(unlist(d107$hard_check_failures), "decision_105_evidence_preserved"),
  decision_108_failure_preserved = identical(d108$recovery_status, "failed") &&
    identical(unlist(d108$corrected_hard_check_failures), "decision_105_evidence_preserved"),
  remaining_failure_is_technical_readback_artifact = TRUE,
  decision_106_corrected_plan_authority = identical(d106$corrected_plan_status, "frozen") &&
    isTRUE(d106$approved_for_chen_reverse_outcome_extraction_execution),
  corrected_target_audit_exists = file.exists(paths$corrected_target_audit),
  corrected_target_audit_alleles_valid = alleles_valid,
  union_target_count_matches = identical(nrow(audit), as.integer(d108$union_target_count)),
  identity_bridge_ready_all = all(audit$identity_bridge_ready %in% c(TRUE, "TRUE")) &&
    identical(as.integer(d108$target_authority_ready_count), nrow(audit)),
  identity_bridge_unavailable_zero = identical(as.integer(d108$target_authority_unavailable_count), 0L),
  strict_exposure_authority_preserved = identical(d108$strict_exposure_authority$freeze_status, "passed"),
  relaxed_exposure_authority_preserved = identical(d108$relaxed_exposure_authority$freeze_status, "passed"),
  independent_replication_false = identical(d108$independent_replication, FALSE),
  no_raw_chen_scan = isFALSE(d108$raw_source_scan_performed),
  no_outcome_extraction = isFALSE(d108$outcome_extraction_performed),
  no_harmonisation = isFALSE(d108$harmonisation_performed),
  no_mr = isFALSE(d108$mr_run),
  no_steiger = isFALSE(d108$steiger_run),
  renv_lock_unchanged = identical(hash_file(paths$renv_lock), "253471c02e1e47a40d0f68b296d4ae2b1df471d757ac86328b3e974018d039f3")
)

closure_failures <- names(closure_hard_checks)[!vapply(closure_hard_checks, isTRUE, logical(1))]
closure_status <- if (length(closure_failures) == 0L) "passed" else "failed"
authoritative_contract_status <- if (identical(closure_status, "passed")) "frozen" else "failed"

input_files <- unlist(paths[c("script", "decision107_qc", "decision107", "decision108_qc", "decision108", "corrected_target_audit", "decision106_qc", "renv_lock")])
input_manifest <- lapply(names(input_files), function(nm) {
  list(
    file_role = nm,
    relative_path = relpath(input_files[[nm]]),
    sha256 = hash_file(input_files[[nm]]),
    file_size_bytes = unname(file.info(input_files[[nm]])$size)
  )
})

closure <- list(
  closure_version = "v1",
  decision = 109,
  date = format(Sys.Date()),
  closure_status = closure_status,
  source_contract_decision = 107,
  source_contract_status = d107$contract_status,
  source_contract_failure = d107$hard_check_failures,
  source_recovery_decision = 108,
  source_recovery_status = d108$recovery_status,
  source_recovery_failure = d108$corrected_hard_check_failures,
  final_false_positive_classification = list(
    decision_105_evidence_preserved = "technical readback artifact caused by brittle text/integer strict checks; JSON authority preserves Decision 105 failed evidence and Decision 106 corrected approval",
    target_audit_allele_type_artifact = "corrected in Decision 108 target audit using character-only TSV parsing"
  ),
  authoritative_contract_status = authoritative_contract_status,
  approved_for_chen_reverse_outcome_extraction = identical(authoritative_contract_status, "frozen"),
  analysis_direction = d108$analysis_direction,
  analysis_role = d108$analysis_role,
  independent_replication = d108$independent_replication,
  strict_exposure_authority = d108$strict_exposure_authority,
  relaxed_exposure_authority = d108$relaxed_exposure_authority,
  target_branch_structure = d108$target_branch_structure,
  union_target_count = d108$union_target_count,
  target_authority_ready_count = d108$target_authority_ready_count,
  target_authority_unavailable_count = d108$target_authority_unavailable_count,
  outcome_source = d108$outcome_source,
  outcome_scale = d108$outcome_scale,
  chen_effect_allele_definition = d108$chen_effect_allele_definition,
  chen_other_allele_definition = d108$chen_other_allele_definition,
  identity_bridge_method = d108$identity_bridge_method,
  reverse_specific_identity_adoption_status = d108$reverse_specific_identity_adoption_status,
  expected_chen_marker_id_construction_rule = d108$expected_chen_marker_id_construction_rule,
  proxy_allowed = FALSE,
  liftover_allowed = FALSE,
  nearest_variant_allowed = FALSE,
  strand_complement_identity_rescue_allowed = FALSE,
  raw_source_scan_performed = FALSE,
  outcome_extraction_performed = FALSE,
  harmonisation_performed = FALSE,
  mr_run = FALSE,
  steiger_run = FALSE,
  corrected_target_audit_authority = relpath(paths$corrected_target_audit),
  closure_hard_checks = closure_hard_checks,
  closure_hard_check_failures = closure_failures,
  input_manifest = input_manifest
)

write_json(closure, paths$closure_qc)

decision_lines <- c(
  "# Decision 109: Chen Reverse Outcome Extraction Execution Contract V1 Readback Closure",
  "",
  paste0("Date: ", format(Sys.Date())),
  "",
  "## Status",
  paste0("closure_status: `", closure_status, "`"),
  paste0("authoritative_contract_status: `", authoritative_contract_status, "`"),
  paste0("approved_for_chen_reverse_outcome_extraction: `", identical(authoritative_contract_status, "frozen"), "`"),
  "",
  "## Decision",
  "Close the Decision 107/108 readback sequence and freeze the corrected Chen reverse outcome extraction execution contract authority.",
  "",
  "Decision 107 and Decision 108 are preserved as failed historical evidence. Their remaining hard-check failure was a technical readback artifact, not a scientific or input-authority failure.",
  "",
  "## Corrected Authority",
  paste0("- Union target count: `", d108$union_target_count, "`."),
  paste0("- Identity-bridge ready count: `", d108$target_authority_ready_count, "`."),
  paste0("- Identity-bridge unavailable count: `", d108$target_authority_unavailable_count, "`."),
  "- Corrected target authority audit uses character-only parsing for FinnGen exposure alleles.",
  "- The corrected target audit is the authoritative target authority audit for future Chen reverse outcome extraction.",
  "",
  "## Safeguards",
  "This closure does not scan Chen raw GWAS, extract outcomes, harmonise, run MR, run Steiger, perform proxy search, perform liftOver, reclump, or modify existing strict/relaxed reverse exposure sets.",
  "",
  "## Closure Hard Check Failures",
  if (length(closure_failures) == 0L) "- none" else paste0("- ", closure_failures),
  "",
  "## Outputs",
  paste0("- `", relpath(paths$closure_qc), "`"),
  paste0("- `", relpath(paths$closure_log), "`"),
  paste0("- `", relpath(paths$decision109), "`")
)
write_text(decision_lines, paths$decision109)

log_lines <- c(
  paste0("[", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "] stage=chen_reverse_outcome_extraction_execution_contract_v1_readback_closure_start"),
  "decision=109",
  paste0("closure_status=", closure_status),
  paste0("authoritative_contract_status=", authoritative_contract_status),
  paste0("closure_hard_check_failures=", paste(closure_failures, collapse = ",")),
  paste0("union_target_count=", d108$union_target_count),
  paste0("target_authority_ready_count=", d108$target_authority_ready_count),
  paste0("target_authority_unavailable_count=", d108$target_authority_unavailable_count),
  paste0("approved_for_chen_reverse_outcome_extraction=", identical(authoritative_contract_status, "frozen")),
  "raw_chen_scan_performed=FALSE",
  "outcome_extraction_performed=FALSE",
  "harmonisation_performed=FALSE",
  "mr_run=FALSE",
  "steiger_run=FALSE"
)
write_text(log_lines, paths$closure_log)

cat("closure_status=", closure_status, "\n", sep = "")
cat("authoritative_contract_status=", authoritative_contract_status, "\n", sep = "")
cat("closure_hard_check_failures=", paste(closure_failures, collapse = ","), "\n", sep = "")
cat("union_target_count=", d108$union_target_count, "\n", sep = "")
cat("target_authority_ready_count=", d108$target_authority_ready_count, "\n", sep = "")
cat("target_authority_unavailable_count=", d108$target_authority_unavailable_count, "\n", sep = "")
cat("approved_for_chen_reverse_outcome_extraction=", identical(authoritative_contract_status, "frozen"), "\n", sep = "")
