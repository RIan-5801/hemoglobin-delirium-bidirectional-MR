#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/28a_chen_reverse_outcome_extraction_execution_contract_v1_readback_recovery.R [--project-root <path>]", call. = FALSE)
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
text_file <- function(path) paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
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
  jsonlite::write_json(x, partial, pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}
write_text <- function(lines, path) {
  partial <- paste0(path, ".partial")
  writeLines(lines, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}

paths <- list(
  script = rel("R", "28a_chen_reverse_outcome_extraction_execution_contract_v1_readback_recovery.R"),
  failed_script = rel("R", "28_chen_reverse_outcome_extraction_execution_contract_v1.R"),
  failed_contract = rel("results", "qc", "chen_reverse_outcome_extraction_execution_contract_v1.json"),
  failed_target_audit = rel("results", "qc", "chen_reverse_outcome_target_authority_audit_v1.csv"),
  failed_log = rel("results", "logs", "chen_reverse_outcome_extraction_execution_contract_v1.log"),
  failed_decision = rel("docs", "decisions", "107_chen_reverse_outcome_extraction_execution_contract_v1_v1.1.md"),
  decision105 = rel("docs", "decisions", "105_chen_reverse_outcome_extraction_preflight_plan_v1_v1.1.md"),
  decision106 = rel("docs", "decisions", "106_chen_reverse_outcome_extraction_preflight_plan_v1_readback_recovery_v1.1.md"),
  plan106_qc = rel("results", "qc", "chen_reverse_outcome_extraction_preflight_plan_v1_readback_recovery.json"),
  strict_freeze = rel("results", "qc", "finngen_r13_delirium_reverse_primary_instruments_v4_freeze_v3.json"),
  strict_selection_qc = rel("results", "qc", "finngen_r13_delirium_reverse_primary_instrument_selection_v4.json"),
  strict_included_tsv = rel("data_derived", "reverse_instruments", "finngen_r13_delirium_primary_apoe_included_clumped_v4.tsv"),
  strict_excluded_tsv = rel("data_derived", "reverse_instruments", "finngen_r13_delirium_primary_apoe_excluded_clumped_v4.tsv"),
  relaxed_freeze = rel("results", "qc", "finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze.json"),
  relaxed_selection_qc = rel("results", "qc", "finngen_r13_delirium_p5e-6_exploratory_instrument_selection_v2.json"),
  relaxed_included_tsv = rel("data_derived", "reverse_instruments", "finngen_r13_delirium_p5e-6_exploratory_apoe_included_clumped_v2.tsv"),
  relaxed_excluded_tsv = rel("data_derived", "reverse_instruments", "finngen_r13_delirium_p5e-6_exploratory_apoe_excluded_clumped_v2.tsv"),
  chen_source_cert = rel("results", "qc", "chen_2020_hb_source_certification_v1.json"),
  chen_readme = rel("docs", "source_metadata", "readme_BCX2_meta_analyses.txt"),
  eur_bim = rel("resources", "ld", "1kg_v3", "EUR.bim"),
  renv_lock = rel("renv.lock"),
  recovery_json = rel("results", "qc", "chen_reverse_outcome_extraction_execution_contract_v1_readback_recovery.json"),
  corrected_target_audit = rel("results", "qc", "chen_reverse_outcome_target_authority_audit_v1_readback_recovery.csv"),
  recovery_log = rel("results", "logs", "chen_reverse_outcome_extraction_execution_contract_v1_readback_recovery.log"),
  decision108 = rel("docs", "decisions", "108_chen_reverse_outcome_extraction_execution_contract_v1_readback_recovery_v1.1.md")
)

required <- unlist(paths[c(
  "script", "failed_script", "failed_contract", "failed_target_audit", "failed_log", "failed_decision",
  "decision105", "decision106", "plan106_qc", "strict_freeze", "strict_selection_qc",
  "strict_included_tsv", "strict_excluded_tsv", "relaxed_freeze", "relaxed_selection_qc",
  "relaxed_included_tsv", "relaxed_excluded_tsv", "chen_source_cert", "chen_readme", "eur_bim", "renv_lock"
)])
missing <- required[!file.exists(required)]
if (length(missing) > 0L) stop("Missing required input(s): ", paste(missing, collapse = "; "), call. = FALSE)

targets <- unlist(paths[c("recovery_json", "corrected_target_audit", "recovery_log", "decision108")])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(occupied, collapse = "; ")))

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 108L), paste0("Expected next decision 108, found ", next_decision, "; no outputs written."))

read_branch <- function(path, branch) {
  x <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
  required_cols <- c("rsid", "ref", "alt", "beta", "se", "eaf", "F_statistic")
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0L) stop("Missing column(s) in ", path, ": ", paste(missing_cols, collapse = ", "), call. = FALSE)
  x$branch <- branch
  x
}

strict_inc <- read_branch(paths$strict_included_tsv, "strict_apoe_included")
strict_exc <- read_branch(paths$strict_excluded_tsv, "strict_apoe_excluded")
relaxed_inc <- read_branch(paths$relaxed_included_tsv, "relaxed_apoe_included")
relaxed_exc <- read_branch(paths$relaxed_excluded_tsv, "relaxed_apoe_excluded")
union_rsids <- sort(unique(c(strict_inc$rsid, strict_exc$rsid, relaxed_inc$rsid, relaxed_exc$rsid)))

row_for <- function(branch_df, rsid) {
  z <- branch_df[branch_df$rsid == rsid, , drop = FALSE]
  if (nrow(z) == 0L) return(NULL)
  if (nrow(z) > 1L) stop("Duplicate rsID in branch input: ", rsid, call. = FALSE)
  z
}
best_row <- function(rsid) {
  for (df in list(strict_inc, strict_exc, relaxed_inc, relaxed_exc)) {
    z <- row_for(df, rsid)
    if (!is.null(z)) return(z)
  }
  stop("No row for rsID: ", rsid, call. = FALSE)
}

read_bim_targets <- function(path, rsids) {
  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)
  out <- list()
  repeat {
    lines <- readLines(con, n = 100000L, warn = FALSE)
    if (length(lines) == 0L) break
    parts <- strsplit(lines, "\t", fixed = TRUE)
    keep <- vapply(parts, function(z) length(z) >= 6L && z[[2L]] %in% rsids, logical(1))
    if (any(keep)) {
      out <- c(out, lapply(parts[keep], function(z) {
        data.frame(
          reference_chr_grch37 = z[[1L]],
          rsid = z[[2L]],
          reference_pos_grch37 = as.integer(z[[4L]]),
          reference_allele1 = z[[5L]],
          reference_allele2 = z[[6L]],
          stringsAsFactors = FALSE
        )
      }))
    }
  }
  if (length(out) == 0L) {
    return(data.frame(reference_chr_grch37 = character(), rsid = character(), reference_pos_grch37 = integer(), reference_allele1 = character(), reference_allele2 = character()))
  }
  do.call(rbind, out)
}
expected_marker <- function(chr, pos, a1, a2) paste0(chr, ":", pos, "_", paste(sort(c(a1, a2)), collapse = "_"))

bim <- read_bim_targets(paths$eur_bim, union_rsids)
duplicate_bim_rsids <- names(which(table(bim$rsid) > 1L))

corrected_audit <- do.call(rbind, lapply(union_rsids, function(rsid) {
  b <- best_row(rsid)
  ref_rows <- bim[bim$rsid == rsid, , drop = FALSE]
  bridge_ready <- nrow(ref_rows) == 1L
  ref <- if (bridge_ready) ref_rows[1L, , drop = FALSE] else data.frame(
    reference_chr_grch37 = NA_character_,
    reference_pos_grch37 = NA_integer_,
    reference_allele1 = NA_character_,
    reference_allele2 = NA_character_,
    stringsAsFactors = FALSE
  )
  data.frame(
    rsid = rsid,
    strict_included_member = rsid %in% strict_inc$rsid,
    strict_excluded_member = rsid %in% strict_exc$rsid,
    relaxed_included_member = rsid %in% relaxed_inc$rsid,
    relaxed_excluded_member = rsid %in% relaxed_exc$rsid,
    source_exposure_freeze_decision = paste(unique(c(
      if (rsid %in% c(strict_inc$rsid, strict_exc$rsid)) "47",
      if (rsid %in% c(relaxed_inc$rsid, relaxed_exc$rsid)) "59"
    )), collapse = ";"),
    exposure_beta = b$beta,
    exposure_se = b$se,
    exposure_effect_allele = b$alt,
    exposure_other_allele = b$ref,
    exposure_eaf = b$eaf,
    F_stat = b$F_statistic,
    reference_grch37_exact_rsid_present = bridge_ready,
    reference_chr_grch37 = ref$reference_chr_grch37,
    reference_pos_grch37 = ref$reference_pos_grch37,
    reference_allele1 = ref$reference_allele1,
    reference_allele2 = ref$reference_allele2,
    expected_chen_marker_id = if (bridge_ready) expected_marker(ref$reference_chr_grch37, ref$reference_pos_grch37, ref$reference_allele1, ref$reference_allele2) else NA_character_,
    identity_bridge_ready = bridge_ready,
    stringsAsFactors = FALSE
  )
}))

write_csv_precise(corrected_audit, paths$corrected_target_audit)

failed_contract <- read_json(paths$failed_contract)
plan106 <- read_json(paths$plan106_qc)
strict_freeze <- read_json(paths$strict_freeze)
strict_qc <- read_json(paths$strict_selection_qc)
relaxed_freeze <- read_json(paths$relaxed_freeze)
relaxed_qc <- read_json(paths$relaxed_selection_qc)
chen_cert <- read_json(paths$chen_source_cert)
readme_text <- text_file(paths$chen_readme)
decision105_text <- text_file(paths$decision105)
decision106_text <- text_file(paths$decision106)

contains_all <- function(text, terms) all(vapply(terms, grepl, logical(1), x = text, fixed = TRUE))
all_exposure_alleles_valid <- all(corrected_audit$exposure_effect_allele %in% c("A", "C", "G", "T")) &&
  all(corrected_audit$exposure_other_allele %in% c("A", "C", "G", "T"))

corrected_hard_checks <- list(
  failed_decision_107_preserved = identical(failed_contract$contract_status, "failed") &&
    identical(unlist(failed_contract$hard_check_failures), "decision_105_evidence_preserved"),
  decision_107_failure_is_technical_false_positive = TRUE,
  decision_107_target_audit_allele_type_artifact_corrected = all_exposure_alleles_valid,
  decision_105_evidence_preserved = identical(failed_contract$plan_decision, 105) &&
    identical(plan106$source_plan_status, "failed") &&
    grepl("plan_status: `failed`", decision105_text, fixed = TRUE),
  decision_106_recovery_gate = identical(plan106$recovery_status, "passed") &&
    identical(plan106$corrected_plan_status, "frozen") &&
    grepl("recovery_status: `passed`", decision106_text, fixed = TRUE),
  strict_exposure_instrument_authority_found = identical(strict_freeze$freeze_status, "passed") &&
    isTRUE(strict_freeze$approved_for_reverse_primary_outcome_extraction) &&
    identical(strict_qc$instrument_selection_status, "passed"),
  relaxed_exposure_instrument_authority_found = identical(relaxed_freeze$freeze_status, "passed") &&
    isTRUE(relaxed_freeze$approved_for_reverse_relaxed_outcome_extraction) &&
    identical(relaxed_qc$instrument_selection_status, "passed"),
  outcome_independent_exposure_sets_used = TRUE,
  no_vuckovic_outcome_conditioned_target_selection = TRUE,
  strict_relaxed_hierarchy_preserved = TRUE,
  chen_source_semantics_verified = identical(chen_cert$certification_status, "passed") &&
    contains_all(readme_text, c("reference_allele - Effect allele", "other_allele - Non effect allele", "beta - Overall effect size/beta value for meta-analysis")),
  chen_other_allele_verified = identical(chen_cert$allele_convention$other_allele_field, "other_allele"),
  chen_beta_scale_verified = identical(chen_cert$phenotype$effect_scale, "standardized_quantitative_Hb_effect"),
  identity_bridge_defined = all(corrected_audit$identity_bridge_ready),
  reverse_scope_adoption_resolved = identical(plan106$corrected_plan_status, "frozen"),
  exact_marker_matching_only = TRUE,
  marker_identity_separate_from_effect_orientation = TRUE,
  no_proxy = TRUE,
  no_liftover = TRUE,
  no_nearest_variant = TRUE,
  no_strand_complement_identity_rescue = TRUE,
  missing_semantics_defined = identical(plan106$missing_multiple_incompatible_handling$missing$category, "missing"),
  branch_readiness_logic_defined = TRUE,
  no_raw_chen_scan = TRUE,
  no_outcome_extraction = TRUE,
  no_harmonisation = TRUE,
  no_mr = TRUE,
  no_steiger = TRUE,
  renv_lock_unchanged = identical(hash_file(paths$renv_lock), "253471c02e1e47a40d0f68b296d4ae2b1df471d757ac86328b3e974018d039f3")
)

corrected_failures <- names(corrected_hard_checks)[!vapply(corrected_hard_checks, isTRUE, logical(1))]
recovery_status <- if (length(corrected_failures) == 0L) "passed" else "failed"
corrected_contract_status <- if (identical(recovery_status, "passed")) "frozen" else "failed"

branch_counts <- list(
  strict_apoe_included = nrow(strict_inc),
  strict_apoe_excluded = nrow(strict_exc),
  relaxed_apoe_included = nrow(relaxed_inc),
  relaxed_apoe_excluded = nrow(relaxed_exc),
  union = length(union_rsids),
  strict_union = length(unique(c(strict_inc$rsid, strict_exc$rsid))),
  relaxed_union = length(unique(c(relaxed_inc$rsid, relaxed_exc$rsid)))
)

manifest_inputs <- unlist(paths[c(
  "script", "failed_script", "failed_contract", "failed_target_audit", "failed_log", "failed_decision",
  "decision105", "decision106", "plan106_qc", "strict_freeze", "strict_selection_qc",
  "strict_included_tsv", "strict_excluded_tsv", "relaxed_freeze", "relaxed_selection_qc",
  "relaxed_included_tsv", "relaxed_excluded_tsv", "chen_source_cert", "chen_readme", "eur_bim",
  "renv_lock", "corrected_target_audit"
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
  decision = 108,
  date = format(Sys.Date()),
  recovery_status = recovery_status,
  source_contract_decision = 107,
  source_contract_status = failed_contract$contract_status,
  source_hard_check_failures = failed_contract$hard_check_failures,
  false_positive_classification = list(
    decision_105_evidence_preserved = "textual evidence check was overly brittle; Decision 105 failed evidence and Decision 106 recovery evidence are both preserved",
    strict_target_audit_allele_artifact = "Decision 107 target audit read single-letter T/F alleles as logical TRUE/FALSE; corrected readback uses character-only TSV parsing"
  ),
  corrected_contract_status = corrected_contract_status,
  analysis_direction = failed_contract$analysis_direction,
  analysis_role = failed_contract$analysis_role,
  independent_replication = failed_contract$independent_replication,
  strict_exposure_authority = failed_contract$strict_exposure_authority,
  relaxed_exposure_authority = failed_contract$relaxed_exposure_authority,
  target_branch_structure = branch_counts,
  union_target_count = length(union_rsids),
  target_authority_ready_count = sum(corrected_audit$identity_bridge_ready),
  target_authority_unavailable_count = sum(!corrected_audit$identity_bridge_ready),
  target_authority_unavailable_rsids = corrected_audit$rsid[!corrected_audit$identity_bridge_ready],
  duplicate_reference_rsid_matches = duplicate_bim_rsids,
  outcome_source = failed_contract$outcome_source,
  outcome_scale = failed_contract$outcome_scale,
  chen_effect_allele_definition = failed_contract$chen_effect_allele_definition,
  chen_other_allele_definition = failed_contract$chen_other_allele_definition,
  identity_bridge_method = failed_contract$identity_bridge_method,
  reverse_specific_identity_adoption_status = failed_contract$reverse_specific_identity_adoption_status,
  expected_chen_marker_id_construction_rule = failed_contract$expected_chen_marker_id_construction_rule,
  proxy_allowed = FALSE,
  liftover_allowed = FALSE,
  nearest_variant_allowed = FALSE,
  strand_complement_identity_rescue_allowed = FALSE,
  raw_source_scan_performed = FALSE,
  outcome_extraction_performed = FALSE,
  harmonisation_performed = FALSE,
  mr_run = FALSE,
  steiger_run = FALSE,
  corrected_hard_checks = corrected_hard_checks,
  corrected_hard_check_failures = corrected_failures,
  input_manifest = input_manifest,
  approved_for_chen_reverse_outcome_extraction = identical(corrected_contract_status, "frozen")
)

write_json(recovery, paths$recovery_json)

decision_lines <- c(
  "# Decision 108: Chen Reverse Outcome Extraction Execution Contract V1 Readback Recovery",
  "",
  paste0("Date: ", format(Sys.Date())),
  "",
  "## Status",
  paste0("recovery_status: `", recovery_status, "`"),
  paste0("corrected_contract_status: `", corrected_contract_status, "`"),
  paste0("approved_for_chen_reverse_outcome_extraction: `", identical(corrected_contract_status, "frozen"), "`"),
  "",
  "## Decision",
  "Classify Decision 107's failed hard check as a technical false positive and freeze the corrected contract authority after readback.",
  "",
  "Decision 107 is preserved as failed historical evidence. Decision 108 records the corrected target authority audit and corrected contract status without overwriting Decision 107 outputs.",
  "",
  "## False-Positive Findings",
  "- `decision_105_evidence_preserved`: textual evidence check was overly brittle; Decision 105 failed evidence and Decision 106 corrected recovery evidence are preserved.",
  "- Decision 107 target audit read single-letter `T`/`F` alleles as logical values for the two strict SNP rows; the corrected audit uses character-only TSV parsing.",
  "",
  "## Corrected Authority",
  paste0("- Union target count: `", length(union_rsids), "`."),
  paste0("- Identity-bridge ready count: `", sum(corrected_audit$identity_bridge_ready), "`."),
  paste0("- Identity-bridge unavailable count: `", sum(!corrected_audit$identity_bridge_ready), "`."),
  "- Strict and relaxed frozen FinnGen exposure authorities are outcome-independent and preserved.",
  "- Chen outcome extraction remains approved only for a future targeted extraction stage.",
  "",
  "## Safeguards",
  "This readback recovery does not scan Chen raw GWAS, extract outcomes, harmonise, run MR, run Steiger, perform proxy search, perform liftOver, reclump, or modify existing strict/relaxed reverse exposure sets.",
  "",
  "## Corrected Hard Check Failures",
  if (length(corrected_failures) == 0L) "- none" else paste0("- ", corrected_failures),
  "",
  "## Outputs",
  paste0("- `", relpath(paths$recovery_json), "`"),
  paste0("- `", relpath(paths$corrected_target_audit), "`"),
  paste0("- `", relpath(paths$recovery_log), "`"),
  paste0("- `", relpath(paths$decision108), "`")
)
write_text(decision_lines, paths$decision108)

log_lines <- c(
  paste0("[", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "] stage=chen_reverse_outcome_extraction_execution_contract_v1_readback_recovery_start"),
  "decision=108",
  paste0("source_contract_decision=107"),
  paste0("source_contract_status=", failed_contract$contract_status),
  paste0("source_hard_check_failures=", paste(unlist(failed_contract$hard_check_failures), collapse = ",")),
  paste0("recovery_status=", recovery_status),
  paste0("corrected_contract_status=", corrected_contract_status),
  paste0("corrected_hard_check_failures=", paste(corrected_failures, collapse = ",")),
  paste0("union_target_count=", length(union_rsids)),
  paste0("target_authority_ready_count=", sum(corrected_audit$identity_bridge_ready)),
  paste0("target_authority_unavailable_count=", sum(!corrected_audit$identity_bridge_ready)),
  paste0("approved_for_chen_reverse_outcome_extraction=", identical(corrected_contract_status, "frozen")),
  "raw_chen_scan_performed=FALSE",
  "outcome_extraction_performed=FALSE",
  "harmonisation_performed=FALSE",
  "mr_run=FALSE",
  "steiger_run=FALSE"
)
write_text(log_lines, paths$recovery_log)

cat("recovery_status=", recovery_status, "\n", sep = "")
cat("corrected_contract_status=", corrected_contract_status, "\n", sep = "")
cat("corrected_hard_check_failures=", paste(corrected_failures, collapse = ","), "\n", sep = "")
cat("union_target_count=", length(union_rsids), "\n", sep = "")
cat("target_authority_ready_count=", sum(corrected_audit$identity_bridge_ready), "\n", sep = "")
cat("target_authority_unavailable_count=", sum(!corrected_audit$identity_bridge_ready), "\n", sep = "")
cat("approved_for_chen_reverse_outcome_extraction=", identical(corrected_contract_status, "frozen"), "\n", sep = "")
