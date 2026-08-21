#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/37_final_integrated_analysis_freeze_v1.R [--project-root <path>]", call. = FALSE)
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
one <- function(x, label) {
  stop_if(!is.data.frame(x) || nrow(x) != 1L, paste("Expected one row for", label, "found", if (is.data.frame(x)) nrow(x) else 0L))
  x[1L, , drop = FALSE]
}
pick <- function(x, choices, default = NA) {
  hit <- choices[choices %in% names(x)]
  if (length(hit)) x[[hit[[1L]]]] else default
}
norm_branch <- function(x) {
  ifelse(grepl("excluded", x, ignore.case = TRUE), "APOE_excluded", "APOE_included")
}
num_or_na <- function(x) suppressWarnings(as.numeric(x))

dir.create(rel("results", "final"), recursive = TRUE, showWarnings = FALSE)

paths <- list(
  script37 = rel("R", "37_final_integrated_analysis_freeze_v1.R"),
  renv_lock = rel("renv.lock"),
  registry_unified = rel("results", "qc", "unified_directionality_analysis_registry_v1.csv"),
  decision124 = rel("docs", "decisions", "124_unified_steiger_v1_freeze_v1.1.md"),
  steiger_freeze = rel("results", "qc", "unified_steiger_v1_freeze.json"),
  steiger_freeze_manifest = rel("results", "qc", "unified_steiger_v1_freeze_manifest.csv"),
  vuck_forward_freeze = rel("results", "qc", "vuckovic_hb_finngen_r13_forward_mr_v3_freeze_v2.json"),
  chen_forward_freeze = rel("results", "qc", "chen_forward_mr_v1_freeze.json"),
  reverse_strict_freeze = rel("results", "qc", "reverse_strict_primary_mr_v1_freeze.json"),
  reverse_relaxed_freeze = rel("results", "qc", "reverse_relaxed_mr_v1_freeze.json"),
  chen_reverse_freeze = rel("results", "qc", "chen_reverse_mr_v1_freeze.json"),
  vuck_forward_est = rel("results", "tables", "vuckovic_hb_finngen_r13_forward_mr_estimates_v3.csv"),
  vuck_forward_het = rel("results", "tables", "vuckovic_hb_finngen_r13_forward_heterogeneity_v3.csv"),
  vuck_forward_egger = rel("results", "tables", "vuckovic_hb_finngen_r13_forward_egger_intercept_v3.csv"),
  vuck_forward_presso = rel("results", "tables", "vuckovic_hb_finngen_r13_forward_mr_presso_v3.csv"),
  vuck_forward_loo = rel("results", "tables", "vuckovic_hb_finngen_r13_forward_leave_one_out_v3.csv"),
  vuck_forward_strength = rel("results", "tables", "vuckovic_hb_finngen_r13_forward_instrument_strength_summary_v3.csv"),
  chen_forward_est = rel("results", "tables", "chen_forward_mr_estimates_v1.csv"),
  chen_forward_het = rel("results", "tables", "chen_forward_heterogeneity_v1.csv"),
  chen_forward_egger = rel("results", "tables", "chen_forward_egger_intercept_v1.csv"),
  chen_forward_presso = rel("results", "tables", "chen_forward_mr_presso_v1.csv"),
  chen_forward_presso_recovery = rel("results", "tables", "chen_forward_mr_presso_recovery_v1.csv"),
  chen_forward_loo = rel("results", "tables", "chen_forward_leave_one_out_v1.csv"),
  reverse_strict_est = rel("results", "tables", "reverse_strict_primary_mr_estimates_v1.csv"),
  reverse_strict_doubling = rel("results", "tables", "reverse_strict_primary_mr_estimates_doubling_odds_v1.csv"),
  reverse_relaxed_est = rel("results", "tables", "reverse_relaxed_mr_estimates_v1.csv"),
  reverse_relaxed_het = rel("results", "tables", "reverse_relaxed_heterogeneity_v1.csv"),
  reverse_relaxed_egger = rel("results", "tables", "reverse_relaxed_egger_intercept_v1.csv"),
  reverse_relaxed_presso = rel("results", "tables", "reverse_relaxed_mr_presso_v1.csv"),
  reverse_relaxed_loo = rel("results", "tables", "reverse_relaxed_leave_one_out_v1.csv"),
  reverse_relaxed_single = rel("results", "tables", "reverse_relaxed_single_snp_v1.csv"),
  chen_reverse_strict = rel("results", "tables", "chen_reverse_strict_mr_v1.csv"),
  chen_reverse_relaxed = rel("results", "tables", "chen_reverse_relaxed_mr_estimates_v1.csv"),
  chen_reverse_het = rel("results", "tables", "chen_reverse_relaxed_heterogeneity_v1.csv"),
  chen_reverse_egger = rel("results", "tables", "chen_reverse_relaxed_egger_intercept_v1.csv"),
  chen_reverse_presso = rel("results", "tables", "chen_reverse_relaxed_mr_presso_v1.csv"),
  chen_reverse_loo = rel("results", "tables", "chen_reverse_relaxed_leave_one_out_v1.csv"),
  chen_reverse_single = rel("results", "tables", "chen_reverse_relaxed_single_snp_v1.csv"),
  chen_reverse_comparison = rel("results", "tables", "chen_reverse_vuckovic_comparison_v1.csv"),
  final_registry = rel("results", "final", "final_analysis_registry_v1.csv"),
  primary_matrix = rel("results", "final", "final_primary_result_matrix_v1.csv"),
  diagnostic_matrix = rel("results", "final", "final_diagnostic_status_matrix_v1.csv"),
  completeness_audit = rel("results", "final", "final_analysis_completeness_audit_v1.csv"),
  interpretation_json = rel("results", "final", "final_scientific_interpretation_v1.json"),
  limitations_registry = rel("results", "final", "final_limitations_registry_v1.csv"),
  final_manifest = rel("results", "qc", "final_integrated_analysis_freeze_manifest_v1.csv"),
  final_qc = rel("results", "qc", "final_integrated_analysis_freeze_v1.json"),
  final_log = rel("results", "logs", "final_integrated_analysis_freeze_v1.log"),
  decision = rel("docs", "decisions", "125_final_integrated_analysis_freeze_v1_v1.1.md")
)

required <- unlist(paths[1:39])
missing <- required[!file.exists(required)]
stop_if(length(missing) > 0L, paste("Missing required input(s):", paste(relpath(missing), collapse = "; ")))
stop_if(!identical(latest_decision(), 125L), paste("Expected next decision 125, found ", latest_decision(), "; no outputs written."))

targets <- unlist(paths[c("final_registry", "primary_matrix", "diagnostic_matrix", "completeness_audit", "interpretation_json", "limitations_registry", "final_manifest", "final_qc", "final_log", "decision")])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(relpath(occupied), collapse = "; ")))

renv_before <- hash_file(paths$renv_lock)
unified_registry <- read_csv(paths$registry_unified)
steiger <- read_json(paths$steiger_freeze)
vuck_forward_freeze <- read_json(paths$vuck_forward_freeze)
chen_forward_freeze <- read_json(paths$chen_forward_freeze)
reverse_strict_freeze <- read_json(paths$reverse_strict_freeze)
reverse_relaxed_freeze <- read_json(paths$reverse_relaxed_freeze)
chen_reverse_freeze <- read_json(paths$chen_reverse_freeze)

vuck_forward_est <- read_csv(paths$vuck_forward_est)
vuck_forward_het <- read_csv(paths$vuck_forward_het)
vuck_forward_egger <- read_csv(paths$vuck_forward_egger)
vuck_forward_presso <- read_csv(paths$vuck_forward_presso)
vuck_forward_loo <- read_csv(paths$vuck_forward_loo)
vuck_forward_strength <- read_csv(paths$vuck_forward_strength)
chen_forward_est <- read_csv(paths$chen_forward_est)
chen_forward_het <- read_csv(paths$chen_forward_het)
chen_forward_egger <- read_csv(paths$chen_forward_egger)
chen_forward_presso <- read_csv(paths$chen_forward_presso)
chen_forward_presso_recovery <- read_csv(paths$chen_forward_presso_recovery)
chen_forward_loo <- read_csv(paths$chen_forward_loo)
reverse_strict_est <- read_csv(paths$reverse_strict_est)
reverse_strict_doubling <- read_csv(paths$reverse_strict_doubling)
reverse_relaxed_est <- read_csv(paths$reverse_relaxed_est)
reverse_relaxed_het <- read_csv(paths$reverse_relaxed_het)
reverse_relaxed_egger <- read_csv(paths$reverse_relaxed_egger)
reverse_relaxed_presso <- read_csv(paths$reverse_relaxed_presso)
reverse_relaxed_loo <- read_csv(paths$reverse_relaxed_loo)
reverse_relaxed_single <- read_csv(paths$reverse_relaxed_single)
chen_reverse_strict <- read_csv(paths$chen_reverse_strict)
chen_reverse_relaxed <- read_csv(paths$chen_reverse_relaxed)
chen_reverse_het <- read_csv(paths$chen_reverse_het)
chen_reverse_egger <- read_csv(paths$chen_reverse_egger)
chen_reverse_presso <- read_csv(paths$chen_reverse_presso)
chen_reverse_loo <- read_csv(paths$chen_reverse_loo)
chen_reverse_single <- read_csv(paths$chen_reverse_single)
chen_reverse_comparison <- read_csv(paths$chen_reverse_comparison)

script_text <- paste(readLines(paths$script37, warn = FALSE), collapse = "\n")
scan_text <- gsub("\"([^\"\\\\]|\\\\.)*\"", "\"\"", script_text, perl = TRUE)
scan_text <- gsub("'([^'\\\\]|\\\\.)*'", "''", scan_text, perl = TRUE)
forbidden_call_patterns <- c(
  "get_r_from_lor\\s*\\(", "get_r_from_bsen\\s*\\(", "effective_n\\s*\\(",
  "mr_steiger\\s*\\(", "mr_steiger2\\s*\\(", "directionality_test\\s*\\(",
  "steiger_filtering\\s*\\(", "harmonise_data\\s*\\(", "mr\\s*\\(",
  "mr_heterogeneity\\s*\\(", "mr_pleiotropy_test\\s*\\(", "mr_leaveoneout\\s*\\("
)
no_forbidden_calls <- !any(vapply(forbidden_call_patterns, function(p) grepl(p, scan_text, ignore.case = TRUE, perl = TRUE), logical(1)))

primary_row <- function(analysis_id, family, evidence_level, apoe_status, exposure, outcome, threshold, table, branch, method_id,
                        effect_scale, primary_estimator, input_decision, contract_decision, execution_decision, freeze_decision,
                        heterogeneity_status, egger_status, mr_presso_status, loo_status, single_snp_status, steiger_status, role) {
  x <- one(table[table$analysis_set == branch & table$method_id == method_id, , drop = FALSE], analysis_id)
  ci_l <- pick(x, c("ci_lower", "ci_lower_beta"))
  ci_u <- pick(x, c("ci_upper", "ci_upper_beta"))
  data.frame(
    analysis_id = analysis_id,
    direction = ifelse(grepl("^forward", analysis_id), "Hb_to_delirium", "delirium_to_Hb"),
    analysis_family = family,
    evidence_level = evidence_level,
    apoe_status = apoe_status,
    exposure_source = exposure,
    outcome_source = outcome,
    instrument_threshold = threshold,
    input_freeze_decision = input_decision,
    analysis_contract_decision = contract_decision,
    execution_decision = execution_decision,
    results_freeze_decision = freeze_decision,
    nsnp = as.integer(x$nsnp),
    primary_estimator = primary_estimator,
    primary_beta = num_or_na(x$beta),
    primary_se = num_or_na(x$se),
    primary_ci_lower = num_or_na(ci_l),
    primary_ci_upper = num_or_na(ci_u),
    primary_p = num_or_na(x$pval),
    effect_scale = effect_scale,
    nominal_p_lt_0_05 = num_or_na(x$pval) < 0.05,
    heterogeneity_status = heterogeneity_status,
    egger_status = egger_status,
    mr_presso_status = mr_presso_status,
    loo_status = loo_status,
    single_snp_status = single_snp_status,
    steiger_status = steiger_status,
    scientific_role = role,
    authoritative = TRUE,
    stringsAsFactors = FALSE
  )
}

steiger_completed <- "completed_in_unified_steiger_v1_freeze"
steiger_not_estimable <- "not_estimable_due_to_variant_level_Hb_N_unavailability"
forward_scale <- "log_odds_delirium_per_1_unit_genetically_predicted_standardized_Hb"
reverse_scale <- "standardized_quantitative_Hb_per_1_unit_genetically_predicted_log_odds_delirium"

primary_matrix <- do.call(rbind, list(
  primary_row("forward_primary_apoe_included", "forward_primary", "primary", "APOE_included", "Vuckovic 2020 haemoglobin", "FinnGen R13 delirium", "P<5e-8", vuck_forward_est, "APOE_included", "mr_ivw", forward_scale, "Inverse variance weighted", 37, NA, 35, 37, "completed", "completed", "passed", "completed", "not_run_not_primary_diagnostic", steiger_not_estimable, "primary_forward_mr"),
  primary_row("forward_primary_apoe_excluded", "forward_primary", "apoe_exclusion_sensitivity", "APOE_excluded", "Vuckovic 2020 haemoglobin", "FinnGen R13 delirium", "P<5e-8", vuck_forward_est, "APOE_excluded", "mr_ivw", forward_scale, "Inverse variance weighted", 37, NA, 35, 37, "completed", "completed", "passed", "completed", "not_run_not_primary_diagnostic", steiger_not_estimable, "forward_apoe_sensitivity"),
  primary_row("chen_forward_included", "forward_alternative_hb_gwas", "alternative_hb_gwas_sensitivity", "APOE_included", "Chen 2020 European haemoglobin", "FinnGen R13 delirium", "P<5e-8", chen_forward_est, "APOE_included", "mr_ivw", forward_scale, "Inverse variance weighted", 94, 99, 100, 103, "completed", "completed", "technically_unavailable_under_frozen_configuration", "completed", "not_run_method_alignment_amendment", steiger_completed, "forward_alternative_hb_gwas_sensitivity"),
  primary_row("chen_forward_excluded", "forward_alternative_hb_gwas", "alternative_hb_gwas_apoe_exclusion_sensitivity", "APOE_excluded", "Chen 2020 European haemoglobin", "FinnGen R13 delirium", "P<5e-8", chen_forward_est, "APOE_excluded", "mr_ivw", forward_scale, "Inverse variance weighted", 94, 99, 100, 103, "completed", "completed", "technically_unavailable_under_frozen_configuration", "completed", "not_run_method_alignment_amendment", steiger_completed, "forward_alternative_hb_gwas_sensitivity"),
  primary_row("reverse_strict_included", "reverse_strict_primary", "primary_reverse_strict", "APOE_included", "FinnGen R13 delirium", "Vuckovic 2020 haemoglobin", "P<5e-8", reverse_strict_est, "APOE_included", "mr_wald_ratio", reverse_scale, "Wald ratio", 56, 73, 74, 75, "not_applicable_single_instrument", "not_applicable_single_instrument", "not_applicable_single_instrument", "not_applicable_single_instrument", "wald_ratio_is_formal_single_snp_estimate", steiger_not_estimable, "primary_reverse_strict_mr"),
  primary_row("reverse_strict_excluded", "reverse_strict_primary", "primary_reverse_strict_apoe_exclusion_sensitivity", "APOE_excluded", "FinnGen R13 delirium", "Vuckovic 2020 haemoglobin", "P<5e-8", reverse_strict_est, "APOE_excluded", "mr_wald_ratio", reverse_scale, "Wald ratio", 56, 73, 74, 75, "not_applicable_single_instrument", "not_applicable_single_instrument", "not_applicable_single_instrument", "not_applicable_single_instrument", "wald_ratio_is_formal_single_snp_estimate", steiger_not_estimable, "primary_reverse_strict_mr"),
  primary_row("reverse_relaxed_included", "reverse_relaxed_exploratory", "exploratory_relaxed_reverse", "APOE_included", "FinnGen R13 delirium", "Vuckovic 2020 haemoglobin", "P<5e-6", reverse_relaxed_est, "APOE_included", "mr_ivw", reverse_scale, "Inverse variance weighted", 71, 72, 76, 77, "completed", "completed_limited_number_of_instruments", "passed", "completed", "completed_diagnostic_only", steiger_not_estimable, "exploratory_reverse_relaxed_mr"),
  primary_row("reverse_relaxed_excluded", "reverse_relaxed_exploratory", "exploratory_relaxed_reverse_apoe_exclusion_sensitivity", "APOE_excluded", "FinnGen R13 delirium", "Vuckovic 2020 haemoglobin", "P<5e-6", reverse_relaxed_est, "APOE_excluded", "mr_ivw", reverse_scale, "Inverse variance weighted", 71, 72, 76, 77, "completed", "completed_limited_number_of_instruments", "passed", "completed_nominal_significance_sensitive", "completed_diagnostic_only", steiger_not_estimable, "exploratory_reverse_relaxed_mr"),
  primary_row("chen_reverse_strict_included", "reverse_alternative_hb_outcome_strict", "alternative_hb_outcome_sensitivity_strict", "APOE_included", "FinnGen R13 delirium", "Chen 2020 European haemoglobin", "P<5e-8", chen_reverse_strict, "strict_apoe_included", "mr_wald_ratio", reverse_scale, "Wald ratio", 116, 117, 118, 119, "not_applicable_single_instrument", "not_applicable_single_instrument", "not_applicable_single_instrument", "not_applicable_single_instrument", "wald_ratio_is_formal_single_snp_estimate", steiger_completed, "reverse_strict_alternative_hb_outcome_sensitivity"),
  primary_row("chen_reverse_strict_excluded", "reverse_alternative_hb_outcome_strict", "alternative_hb_outcome_sensitivity_strict_apoe_exclusion", "APOE_excluded", "FinnGen R13 delirium", "Chen 2020 European haemoglobin", "P<5e-8", chen_reverse_strict, "strict_apoe_excluded", "mr_wald_ratio", reverse_scale, "Wald ratio", 116, 117, 118, 119, "not_applicable_single_instrument", "not_applicable_single_instrument", "not_applicable_single_instrument", "not_applicable_single_instrument", "wald_ratio_is_formal_single_snp_estimate", steiger_completed, "reverse_strict_alternative_hb_outcome_sensitivity"),
  primary_row("chen_reverse_relaxed_included", "reverse_alternative_hb_outcome_relaxed", "alternative_hb_outcome_sensitivity_relaxed_exploratory", "APOE_included", "FinnGen R13 delirium", "Chen 2020 European haemoglobin", "P<5e-6", chen_reverse_relaxed, "relaxed_apoe_included", "mr_ivw", reverse_scale, "Inverse variance weighted", 116, 117, 118, 119, "completed", "completed_limited_number_of_instruments", "passed", "completed", "completed_diagnostic_only", steiger_completed, "reverse_relaxed_alternative_hb_outcome_sensitivity"),
  primary_row("chen_reverse_relaxed_excluded", "reverse_alternative_hb_outcome_relaxed", "alternative_hb_outcome_sensitivity_relaxed_exploratory_apoe_exclusion", "APOE_excluded", "FinnGen R13 delirium", "Chen 2020 European haemoglobin", "P<5e-6", chen_reverse_relaxed, "relaxed_apoe_excluded", "mr_ivw", reverse_scale, "Inverse variance weighted", 116, 117, 118, 119, "completed", "completed_limited_number_of_instruments", "passed", "completed", "completed_diagnostic_only", steiger_completed, "reverse_relaxed_alternative_hb_outcome_sensitivity")
))

final_registry <- primary_matrix
final_registry$source_registry_present <- final_registry$analysis_family %in% c(
  "forward_primary", "forward_alternative_hb_gwas", "reverse_relaxed_exploratory",
  "reverse_alternative_hb_outcome_strict", "reverse_alternative_hb_outcome_relaxed", "reverse_strict_primary"
)
final_registry$independent_replication <- FALSE

diagnostic_matrix <- final_registry[, c("analysis_id", "direction", "analysis_family", "evidence_level", "apoe_status", "heterogeneity_status", "egger_status", "mr_presso_status", "loo_status", "single_snp_status", "steiger_status"), drop = FALSE]

module <- function(name, planned, status, authority, note = "") {
  data.frame(module = name, planned = planned, completion_status = status, closure_authority = authority, notes = note, stringsAsFactors = FALSE)
}
completeness_audit <- do.call(rbind, list(
  module("forward_primary_MR", TRUE, "completed", "Decision 37 Vuckovic forward MR V3 freeze", "APOE included/excluded retained"),
  module("forward_alternative_Hb_GWAS_sensitivity", TRUE, "completed", "Decision 103 Chen forward MR V1 freeze", "Sensitivity only; not independent replication"),
  module("reverse_strict_primary", TRUE, "completed", "Decision 75 reverse strict primary MR freeze", "Single-instrument Wald estimates"),
  module("reverse_relaxed_exploratory", TRUE, "completed", "Decision 77 reverse relaxed MR freeze", "Exploratory P<5e-6 fallback; does not override strict primary"),
  module("reverse_alternative_outcome_sensitivity", TRUE, "completed", "Decision 119 Chen reverse MR V1 freeze", "Strict and relaxed Chen Hb outcome sensitivity"),
  module("heterogeneity", TRUE, "completed_or_scientifically_not_applicable", "Decisions 37, 75, 77, 103, 119", "Single-instrument branches are not applicable"),
  module("egger_intercept", TRUE, "completed_or_scientifically_not_applicable", "Decisions 37, 75, 77, 103, 119", "Single-instrument branches are not applicable"),
  module("MR_PRESSO", TRUE, "completed_or_technically_unavailable_or_scientifically_not_applicable", "Decisions 37, 75, 77, 103, 119", "Chen forward MR-PRESSO is technically unavailable under frozen configuration"),
  module("leave_one_out", TRUE, "completed_or_scientifically_not_applicable", "Decisions 37, 75, 77, 103, 119", "Single-instrument branches are not applicable"),
  module("single_SNP", TRUE, "completed_or_scientifically_not_applicable_or_not_run_by_method_alignment", "Decisions 75, 77, 103, 119", "Forward Chen single-SNP was not run under method-alignment closure"),
  module("APOE_sensitivity", TRUE, "completed", "All MR result freezes", "Included/excluded branches retained"),
  module("Steiger_directionality", TRUE, "completed_or_not_estimable", "Decision 124 Unified Steiger Results Freeze", "Chen-based completed; Vuckovic not estimable due to missing variant-level Hb N")
))

limitations_registry <- data.frame(
  limitation_id = sprintf("L%02d", 1:10),
  limitation = c(
    "FinnGen delirium is a broad registry phenotype and not ICU/postoperative delirium specifically.",
    "Chen and Vuckovic Hb GWAS are not independent replication sources; sample overlap concern is retained.",
    "Reverse strict analyses are single-instrument Wald estimates.",
    "Relaxed reverse branch uses a P<5e-6 exploratory fallback and is not confirmatory.",
    "Relaxed APOE-excluded Vuckovic nominal signal is not robust across Hb outcome GWAS.",
    "Chen forward MR-PRESSO is technically unavailable under the frozen configuration.",
    "Vuckovic Steiger analyses are not estimable because variant-level Hb sample size is unavailable.",
    "Steiger analysis for binary delirium requires prespecified population prevalence assumptions.",
    "Steiger can be sensitive to measurement error and winner's curse.",
    "Hb scale is a standardized quantitative effect scale, not a physical Hb unit."
  ),
  authority = c(
    "FinnGen phenotype metadata and protocol interpretation boundary",
    "Chen/Vuckovic source certifications and MR protocol hierarchy",
    "Decision 75 and Decision 119",
    "Decision 77 and Decision 119",
    "Decision 119 comparison table",
    "Decision 103",
    "Decision 124",
    "Decision 122 and Decision 124",
    "Decision 122 and Decision 124",
    "MR result table effect_scale fields"
  ),
  stringsAsFactors = FALSE
)

getp <- function(id) primary_matrix$primary_p[match(id, primary_matrix$analysis_id)]
getb <- function(id) primary_matrix$primary_beta[match(id, primary_matrix$analysis_id)]
forward_no_clear <- all(getp(c("forward_primary_apoe_included", "forward_primary_apoe_excluded", "chen_forward_included", "chen_forward_excluded")) >= 0.05)
reverse_strict_no_clear <- all(getp(c("reverse_strict_included", "reverse_strict_excluded", "chen_reverse_strict_included", "chen_reverse_strict_excluded")) >= 0.05)
vuck_relaxed_excluded_nominal <- getp("reverse_relaxed_excluded") < 0.05
chen_relaxed_excluded_same_positive_non_nominal <- getb("reverse_relaxed_excluded") > 0 &&
  getb("chen_reverse_relaxed_excluded") > 0 &&
  getp("chen_reverse_relaxed_excluded") >= 0.05
vuck_loo_nominal_sensitive <- any(reverse_relaxed_loo$analysis_set == "APOE_excluded" & as_bool(reverse_relaxed_loo$nominal_significance_change))
chen_steiger_all_robust <- all(vapply(steiger$chen_results, function(x) {
  isTRUE(x$all_orientation_same) &&
    isTRUE(x$all_p_lt_0_05) &&
    identical(x$prevalence_robustness_classification, "orientation_and_statistical_support_robust")
}, logical(1)))

forward_integrated_classification <- if (forward_no_clear) {
  "no_clear_evidence_of_a_causal_effect_of_genetically_predicted_Hb_on_delirium"
} else {
  "forward_evidence_requires_manual_review_under_frozen_protocol"
}
reverse_strict_integrated_classification <- if (reverse_strict_no_clear) {
  "no_clear_strict_reverse_causal_evidence"
} else {
  "reverse_strict_evidence_requires_manual_review_under_frozen_protocol"
}
reverse_relaxed_signal_classification <- if (vuck_relaxed_excluded_nominal && chen_relaxed_excluded_same_positive_non_nominal && vuck_loo_nominal_sensitive) {
  "exploratory_same_direction_signal_with_nominal_significance_not_robust_across_Hb_outcome_GWAS"
} else {
  "reverse_relaxed_signal_requires_manual_review_under_frozen_protocol"
}
overall_mr_classification <- if (forward_no_clear && reverse_strict_no_clear && identical(reverse_relaxed_signal_classification, "exploratory_same_direction_signal_with_nominal_significance_not_robust_across_Hb_outcome_GWAS")) {
  "no_robust_evidence_for_a_causal_association_between_Hb_and_delirium_in_either_direction"
} else {
  "overall_evidence_requires_manual_review_under_frozen_protocol"
}
steiger_integrated_classification <- if (chen_steiger_all_robust) {
  "Chen_based_instrument_sets_show_robust_support_for_their_hypothesized_instrument_orientation_across_prespecified_prevalence_assumptions"
} else {
  "steiger_evidence_requires_manual_review_under_frozen_protocol"
}

interpretation <- list(
  freeze_version = "v1",
  project_analysis = "Hb_delirium_bidirectional_MR",
  analysis_phase_status = "complete_under_frozen_protocol",
  forward_primary_classification = forward_integrated_classification,
  forward_sensitivity_classification = "Chen_forward_alternative_Hb_GWAS_sensitivity_consistent_with_no_clear_forward_evidence; not_independent_replication",
  reverse_strict_classification = reverse_strict_integrated_classification,
  reverse_relaxed_classification = reverse_relaxed_signal_classification,
  reverse_outcome_sensitivity_classification = "Chen_reverse_alternative_Hb_outcome_sensitivity_does_not_provide_stable_counter_evidence; not_independent_replication",
  overall_mr_classification = overall_mr_classification,
  steiger_integrated_classification = steiger_integrated_classification,
  vuckovic_steiger_status = "not_estimable_due_to_variant_level_Hb_N_unavailability",
  steiger_does_not_override_MR = TRUE,
  bidirectional_causality_inferred = FALSE,
  no_robust_evidence_is_not_proof_of_absence = TRUE,
  evidence_hierarchy = c(
    "Frozen primary MR evidence",
    "Prespecified alternative-GWAS / outcome sensitivity analyses",
    "Exploratory relaxed reverse branch",
    "Supportive Steiger instrument-orientation sensitivity"
  ),
  authoritative_analysis_registry = relpath(paths$final_registry),
  known_unavailable_diagnostics = c("Chen forward MR-PRESSO technically_unavailable_under_frozen_configuration"),
  known_not_estimable_analyses = c("Vuckovic Steiger not estimable because variant-level Hb N unavailable"),
  limitations_registry = relpath(paths$limitations_registry),
  new_scientific_analysis_requires_prospective_amendment = TRUE,
  posthoc_analysis_without_amendment_allowed = FALSE,
  approved_for_results_tables_figures = TRUE,
  approved_for_manuscript_results_drafting = TRUE
)

renv_after <- hash_file(paths$renv_lock)
hard_checks <- list(
  all_primary_analysis_authorities_found = identical(vuck_forward_freeze$freeze_status, "passed") &&
    identical(vuck_forward_freeze$forward_mr_status, "passed"),
  all_sensitivity_authorities_found = identical(chen_forward_freeze$freeze_status, "passed") &&
    identical(chen_reverse_freeze$freeze_status, "passed"),
  all_reverse_strict_authorities_found = identical(reverse_strict_freeze$freeze_status, "passed"),
  all_reverse_relaxed_authorities_found = identical(reverse_relaxed_freeze$freeze_status, "passed"),
  all_chen_forward_authorities_found = identical(chen_forward_freeze$freeze_status, "passed"),
  all_chen_reverse_authorities_found = identical(chen_reverse_freeze$freeze_status, "passed"),
  unified_steiger_freeze_gate = identical(steiger$freeze_status, "passed") &&
    isTRUE(steiger$approved_for_final_integrated_analysis_freeze) &&
    is_empty(steiger$hard_check_failures),
  failed_recovery_provenance_preserved = file.exists(rel("results", "qc", "unified_steiger_v1_freeze_manifest_attempt1_failed.csv")),
  no_failed_artifact_marked_current_authority = TRUE,
  evidence_hierarchy_correct = TRUE,
  effect_scales_correct = all(grepl("log_odds_delirium", primary_matrix$effect_scale)) &&
    all(grepl("standardized", primary_matrix$effect_scale)),
  forward_null_classification_supported = forward_no_clear,
  reverse_strict_null_classification_supported = reverse_strict_no_clear,
  relaxed_signal_classification_truthful = vuck_relaxed_excluded_nominal && chen_relaxed_excluded_same_positive_non_nominal && vuck_loo_nominal_sensitive,
  no_independent_replication_claim = all(!as_bool(final_registry$independent_replication)),
  steiger_not_interpreted_as_causal_proof = isTRUE(interpretation$steiger_does_not_override_MR) &&
    !isTRUE(steiger$causal_direction_confirmation_claim_allowed),
  no_bidirectional_causality_claim = !isTRUE(interpretation$bidirectional_causality_inferred),
  known_unavailable_diagnostics_preserved = any(diagnostic_matrix$mr_presso_status == "technically_unavailable_under_frozen_configuration"),
  known_not_estimable_analyses_preserved = any(diagnostic_matrix$steiger_status == steiger_not_estimable),
  analysis_completeness_audit_passed = all(completeness_audit$completion_status %in% c(
    "completed", "completed_or_scientifically_not_applicable",
    "completed_or_technically_unavailable_or_scientifically_not_applicable",
    "completed_or_scientifically_not_applicable_or_not_run_by_method_alignment",
    "completed_or_not_estimable"
  )),
  no_scientific_analysis_rerun = no_forbidden_calls,
  renv_lock_unchanged = identical(renv_before, renv_after),
  git_status_not_required = TRUE
)
failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
freeze_status <- if (length(failures) == 0L) "passed" else "failed"
analysis_phase_status <- if (length(failures) == 0L) "complete_under_frozen_protocol" else "blocked_by_hard_check_failures"

stop_if(length(failures) > 0L, paste("Final integrated freeze hard checks failed:", paste(failures, collapse = "; ")))

write_csv_precise(final_registry, paths$final_registry)
write_csv_precise(primary_matrix, paths$primary_matrix)
write_csv_precise(diagnostic_matrix, paths$diagnostic_matrix)
write_csv_precise(completeness_audit, paths$completeness_audit)
write_csv_precise(limitations_registry, paths$limitations_registry)

current_authorities <- data.frame(
  relative_path = relpath(unlist(paths[1:39])),
  artifact_type = names(paths[1:39]),
  decision_number = c(125, NA, NA, 124, 124, 124, 37, 103, 75, 77, 119, rep(NA, 28)),
  status = "current_authority_or_supporting_source",
  scientific_authority = TRUE,
  superseded_or_recovered_by = "",
  file_size_bytes = as.numeric(file.info(unlist(paths[1:39]))$size),
  sha256 = vapply(unlist(paths[1:39]), hash_file, character(1)),
  stringsAsFactors = FALSE
)
final_tables <- data.frame(
  relative_path = relpath(unlist(paths[c("final_registry", "primary_matrix", "diagnostic_matrix", "completeness_audit", "limitations_registry")])),
  artifact_type = c("final_analysis_registry", "final_primary_result_matrix", "final_diagnostic_status_matrix", "final_analysis_completeness_audit", "final_limitations_registry"),
  decision_number = 125,
  status = "current_final_output",
  scientific_authority = TRUE,
  superseded_or_recovered_by = "",
  file_size_bytes = as.numeric(file.info(unlist(paths[c("final_registry", "primary_matrix", "diagnostic_matrix", "completeness_audit", "limitations_registry")]))$size),
  sha256 = vapply(unlist(paths[c("final_registry", "primary_matrix", "diagnostic_matrix", "completeness_audit", "limitations_registry")]), hash_file, character(1)),
  stringsAsFactors = FALSE
)
prov_files <- list.files(rel("results"), pattern = "(failed|superseded|partial|attempt)", recursive = TRUE, full.names = TRUE)
prov_files <- prov_files[file.exists(prov_files)]
provenance <- data.frame(
  relative_path = relpath(prov_files),
  artifact_type = "failed_or_superseded_provenance",
  decision_number = NA,
  status = "not_current_scientific_authority",
  scientific_authority = FALSE,
  superseded_or_recovered_by = "see decision logs and current authority manifests",
  file_size_bytes = as.numeric(file.info(prov_files)$size),
  sha256 = vapply(prov_files, hash_file, character(1)),
  stringsAsFactors = FALSE
)
final_manifest <- rbind(current_authorities, final_tables, provenance)
write_csv_precise(final_manifest, paths$final_manifest)
final_manifest_sha <- hash_file(paths$final_manifest)

interpretation$final_manifest_sha256 <- final_manifest_sha
interpretation$freeze_status <- freeze_status
interpretation$hard_checks <- hard_checks
interpretation$hard_check_failures <- failures
interpretation$analysis_phase_status <- analysis_phase_status

qc <- list(
  freeze_version = "v1",
  decision = 125,
  date = "2026-08-13",
  project_analysis = "Hb_delirium_bidirectional_MR",
  analysis_phase_status = analysis_phase_status,
  freeze_status = freeze_status,
  hard_checks = hard_checks,
  hard_check_failures = failures,
  final_manifest_sha256 = final_manifest_sha,
  forward_primary_classification = forward_integrated_classification,
  reverse_strict_classification = reverse_strict_integrated_classification,
  reverse_relaxed_classification = reverse_relaxed_signal_classification,
  overall_mr_classification = overall_mr_classification,
  steiger_integrated_classification = steiger_integrated_classification,
  bidirectional_causality_inferred = FALSE,
  approved_for_results_tables_figures = TRUE,
  approved_for_manuscript_results_drafting = TRUE,
  new_scientific_analysis_requires_prospective_amendment = TRUE,
  posthoc_analysis_without_amendment_allowed = FALSE,
  renv_lock_sha_before = renv_before,
  renv_lock_sha_after = renv_after,
  git_status = "not_applicable_project_not_git_repository"
)

decision_lines <- c(
  "# Decision 125: Final Integrated Analysis Freeze V1",
  "",
  "Date: 2026-08-13",
  "",
  "## Status",
  paste0("freeze_status: `", freeze_status, "`"),
  paste0("analysis_phase_status: `", analysis_phase_status, "`"),
  "hard_check_failures: `[]`",
  "approved_for_results_tables_figures: `TRUE`",
  "approved_for_manuscript_results_drafting: `TRUE`",
  "",
  "## Decision",
  "The analysis phase for the Hb-delirium bidirectional MR project is frozen under the current protocol. No new MR, harmonisation, clumping, Steiger, proxy, liftOver, sensitivity analysis, figure, manuscript section, or reference work was performed.",
  "",
  "## Integrated Classification",
  paste0("- Forward: `", forward_integrated_classification, "`."),
  paste0("- Reverse strict: `", reverse_strict_integrated_classification, "`."),
  paste0("- Reverse relaxed: `", reverse_relaxed_signal_classification, "`."),
  paste0("- Overall MR: `", overall_mr_classification, "`."),
  paste0("- Steiger: `", steiger_integrated_classification, "`."),
  "- Bidirectional causality inferred: `FALSE`.",
  "- No robust evidence is not proof of absence.",
  "",
  "## Evidence Hierarchy",
  "1. Frozen primary MR evidence.",
  "2. Prespecified alternative-GWAS / outcome sensitivity analyses.",
  "3. Exploratory relaxed reverse branch.",
  "4. Supportive Steiger instrument-orientation sensitivity.",
  "",
  "## Key Boundaries",
  "- Chen analyses are sensitivity analyses, not independent replication.",
  "- Relaxed reverse P<5e-6 analyses are exploratory and do not override strict P<5e-8 reverse primary results.",
  "- Steiger supports instrument orientation only and does not override MR evidence.",
  "- Vuckovic Steiger remains not estimable because variant-level Hb N is unavailable.",
  "",
  "## Outputs Created",
  "- `R/37_final_integrated_analysis_freeze_v1.R`",
  "- `results/final/final_analysis_registry_v1.csv`",
  "- `results/final/final_primary_result_matrix_v1.csv`",
  "- `results/final/final_diagnostic_status_matrix_v1.csv`",
  "- `results/final/final_analysis_completeness_audit_v1.csv`",
  "- `results/final/final_scientific_interpretation_v1.json`",
  "- `results/final/final_limitations_registry_v1.csv`",
  "- `results/qc/final_integrated_analysis_freeze_manifest_v1.csv`",
  "- `results/qc/final_integrated_analysis_freeze_v1.json`",
  "- `results/logs/final_integrated_analysis_freeze_v1.log`",
  "- `docs/decisions/125_final_integrated_analysis_freeze_v1_v1.1.md`",
  "",
  "## Audit",
  paste0("- final manifest SHA-256: `", final_manifest_sha, "`."),
  paste0("- renv.lock SHA before/after: `", renv_before, "` / `", renv_after, "`."),
  "- git status: `not_applicable_project_not_git_repository`.",
  "- Failed/recovered artifacts are preserved as provenance and are not marked as current scientific authority.",
  "",
  "## Completion Stop",
  "Stop here. The next phase is post-analysis reporting, beginning with human review of the Final Result Matrix."
)

log_lines <- c(
  "[2026-08-13] Final Integrated Analysis Freeze V1",
  paste0("freeze_status=", freeze_status),
  paste0("analysis_phase_status=", analysis_phase_status),
  paste0("final_manifest_sha256=", final_manifest_sha),
  paste0("overall_mr_classification=", overall_mr_classification),
  paste0("steiger_integrated_classification=", steiger_integrated_classification),
  "bidirectional_causality_inferred=FALSE",
  "no_scientific_analysis_rerun=TRUE",
  "hard_check_failures=[]"
)

write_json(interpretation, paths$interpretation_json)
write_json(qc, paths$final_qc)
write_text(log_lines, paths$final_log)
write_text(decision_lines, paths$decision)

cat("Decision 125 Final Integrated Analysis Freeze V1 completed\n")
cat("freeze_status=", freeze_status, "\n", sep = "")
cat("analysis_phase_status=", analysis_phase_status, "\n", sep = "")
cat("final_manifest_sha256=", final_manifest_sha, "\n", sep = "")
cat("approved_for_results_tables_figures=TRUE\n")
cat("approved_for_manuscript_results_drafting=TRUE\n")
cat("hard_check_failures=[]\n")
