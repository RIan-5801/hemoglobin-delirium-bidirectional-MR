#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/33a_chen_reverse_mr_v1_freeze_manifest.R [--project-root <path>]", call. = FALSE)
}
setwd(root)

for (pkg in c("jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

rel <- function(...) file.path(root, ...)
norm <- function(x) gsub("\\\\", "/", x)
relpath <- function(x) norm(sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", norm(normalizePath(x, winslash = "/", mustWork = FALSE))))
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
is_empty <- function(x) is.null(x) || length(x) == 0L
records <- function(x) if (!is.data.frame(x) || nrow(x) == 0L) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
num_equal <- function(a, b, tol = 1e-10) {
  a <- as.numeric(a); b <- as.numeric(b)
  length(a) == length(b) && all((is.na(a) & is.na(b)) | (is.finite(a) & is.finite(b) & abs(a - b) <= tol))
}
bool_col <- function(x) {
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
read_table <- function(path) utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
near_zero_pair <- function(a, b) abs(as.numeric(a)) < 0.001 && abs(as.numeric(b)) < 0.001

classify_pair <- function(v_beta, c_beta, v_p, c_p) {
  v_beta <- as.numeric(v_beta); c_beta <- as.numeric(c_beta)
  v_p <- as.numeric(v_p); c_p <- as.numeric(c_p)
  if (is.finite(v_beta) && is.finite(c_beta) && sign(v_beta) == sign(c_beta) && v_p >= 0.05 && c_p >= 0.05) {
    return("same_direction_both_non_nominal")
  }
  if (is.finite(v_beta) && is.finite(c_beta) && sign(v_beta) != sign(c_beta) && near_zero_pair(v_beta, c_beta) && v_p >= 0.05 && c_p >= 0.05) {
    return("near_zero_estimates_with_trivial_sign_difference")
  }
  if (is.finite(v_beta) && is.finite(c_beta) && sign(v_beta) == sign(c_beta) && v_p < 0.05 && c_p >= 0.05) {
    return("same_direction_nominal_significance_not_reproduced_in_alternative_outcome_sensitivity")
  }
  if (is.finite(v_beta) && is.finite(c_beta) && sign(v_beta) != sign(c_beta) && !near_zero_pair(v_beta, c_beta)) {
    return("material_directional_discordance")
  }
  "other_predefined_neutral_classification"
}

paths <- list(
  decision114 = rel("docs", "decisions", "114_chen_reverse_formal_harmonisation_v1_v1.1.md"),
  decision115 = rel("docs", "decisions", "115_chen_reverse_formal_harmonisation_v1_readback_recovery_v1.1.md"),
  decision116 = rel("docs", "decisions", "116_chen_reverse_harmonised_mr_inputs_v1_freeze_v1.1.md"),
  decision117 = rel("docs", "decisions", "117_chen_reverse_mr_analysis_contract_v1_v1.1.md"),
  decision118 = rel("docs", "decisions", "118_chen_reverse_mr_v1_v1.1.md"),
  decision72 = rel("docs", "decisions", "72_reverse_relaxed_mr_analysis_contract_v1_v1.1.md"),
  decision75 = rel("docs", "decisions", "75_reverse_strict_primary_mr_v1_freeze_v1.1.md"),
  decision77 = rel("docs", "decisions", "77_reverse_relaxed_mr_v1_freeze_v1.1.md"),
  script33 = rel("R", "33_chen_reverse_mr_v1.R"),
  script33a = rel("R", "33a_chen_reverse_mr_v1_freeze_manifest.R"),
  input_freeze_json = rel("results", "qc", "chen_reverse_harmonised_mr_inputs_v1_freeze.json"),
  input_freeze_manifest = rel("results", "qc", "chen_reverse_harmonised_mr_inputs_v1_freeze_manifest.csv"),
  contract_json = rel("results", "qc", "chen_reverse_mr_analysis_contract_v1.json"),
  mr_qc = rel("results", "qc", "chen_reverse_mr_v1.json"),
  mr_log = rel("results", "logs", "chen_reverse_mr_v1.log"),
  strict_table = rel("results", "tables", "chen_reverse_strict_mr_v1.csv"),
  relaxed_table = rel("results", "tables", "chen_reverse_relaxed_mr_estimates_v1.csv"),
  heterogeneity_table = rel("results", "tables", "chen_reverse_relaxed_heterogeneity_v1.csv"),
  egger_table = rel("results", "tables", "chen_reverse_relaxed_egger_intercept_v1.csv"),
  presso_table = rel("results", "tables", "chen_reverse_relaxed_mr_presso_v1.csv"),
  loo_table = rel("results", "tables", "chen_reverse_relaxed_leave_one_out_v1.csv"),
  single_table = rel("results", "tables", "chen_reverse_relaxed_single_snp_v1.csv"),
  comparison_table = rel("results", "tables", "chen_reverse_vuckovic_comparison_v1.csv"),
  vuckovic_strict_freeze = rel("results", "qc", "reverse_strict_primary_mr_v1_freeze.json"),
  vuckovic_relaxed_freeze = rel("results", "qc", "reverse_relaxed_mr_v1_freeze.json"),
  vuckovic_relaxed_contract = rel("results", "qc", "reverse_relaxed_mr_analysis_contract_v1.json"),
  renv_lock = rel("renv.lock"),
  manifest = rel("results", "qc", "chen_reverse_mr_v1_freeze_manifest.csv"),
  freeze_json = rel("results", "qc", "chen_reverse_mr_v1_freeze.json"),
  freeze_log = rel("results", "logs", "chen_reverse_mr_v1_freeze.log"),
  decision = rel("docs", "decisions", "119_chen_reverse_mr_v1_freeze_v1.1.md")
)

required <- unlist(paths[setdiff(names(paths), c("manifest", "freeze_json", "freeze_log", "decision"))])
missing <- required[!file.exists(required)]
stop_if(length(missing) > 0L, paste("Missing required input(s):", paste(missing, collapse = "; ")))

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 119L), paste("Expected next decision 119, found ", next_decision, "; no outputs written."))

targets <- unlist(paths[c("manifest", "freeze_json", "freeze_log", "decision")])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(occupied, collapse = "; ")))

renv_before <- hash_file(paths$renv_lock)
input_freeze <- read_json(paths$input_freeze_json)
contract <- read_json(paths$contract_json)
mr_qc <- read_json(paths$mr_qc)

strict <- read_table(paths$strict_table)
relaxed <- read_table(paths$relaxed_table)
heterogeneity <- read_table(paths$heterogeneity_table)
egger <- read_table(paths$egger_table)
presso <- read_table(paths$presso_table)
loo <- read_table(paths$loo_table)
single <- read_table(paths$single_table)
comparison <- read_table(paths$comparison_table)

manifest_inputs <- setdiff(names(paths), c("manifest", "freeze_json", "freeze_log", "decision"))
manifest_records <- data.frame(
  relative_path = relpath(unlist(paths[manifest_inputs])),
  file_role = manifest_inputs,
  scientific_authority = TRUE,
  file_size_bytes = as.numeric(file.info(unlist(paths[manifest_inputs]))$size),
  sha256 = vapply(unlist(paths[manifest_inputs]), hash_file, character(1)),
  stringsAsFactors = FALSE
)
partial_manifest <- paste0(paths$manifest, ".partial")
on.exit(unlink(partial_manifest, force = TRUE), add = TRUE)
old <- options(digits = 17, scipen = 999)
utils::write.csv(manifest_records, partial_manifest, row.names = FALSE, na = "")
options(old)
manifest_sha <- hash_file(partial_manifest)

script_text <- paste(readLines(paths$script33a, warn = FALSE), collapse = "\n")
scan_text <- gsub("\"([^\"\\\\]|\\\\.)*\"", "\"\"", script_text, perl = TRUE)
scan_text <- gsub("'([^'\\\\]|\\\\.)*'", "''", scan_text, perl = TRUE)
forbidden_call_patterns <- c(
  "TwoSampleMR::", "MRPRESSO::", "mr_wald_ratio\\s*\\(", "mr\\s*\\(",
  "mr_heterogeneity\\s*\\(", "mr_pleiotropy_test\\s*\\(",
  "mr_leaveoneout\\s*\\(", "harmonise_data\\s*\\(",
  "directionality_test\\s*\\(", "mr_steiger\\s*\\(", "steiger_filtering\\s*\\("
)
no_forbidden_calls <- !any(vapply(forbidden_call_patterns, function(p) grepl(p, scan_text, ignore.case = TRUE, perl = TRUE), logical(1)))

strict_included <- strict[strict$analysis_set == "strict_apoe_included", , drop = FALSE]
strict_excluded <- strict[strict$analysis_set == "strict_apoe_excluded", , drop = FALSE]
relaxed_ivw_included <- relaxed[relaxed$analysis_set == "relaxed_apoe_included" & relaxed$method_id == "mr_ivw", , drop = FALSE]
relaxed_ivw_excluded <- relaxed[relaxed$analysis_set == "relaxed_apoe_excluded" & relaxed$method_id == "mr_ivw", , drop = FALSE]

loo_summary <- do.call(rbind, lapply(split(loo, loo$analysis_set), function(x) {
  i <- which.max(as.numeric(x$absolute_shift))
  data.frame(
    analysis_set = x$analysis_set[[1L]],
    full_ivw_beta = unique(as.numeric(x$full_ivw_beta))[1],
    max_absolute_shift = as.numeric(x$absolute_shift[[i]]),
    max_shift_rsid = x$removed_rsid[[i]],
    any_sign_change = any(bool_col(x$sign_change)),
    any_nominal_significance_change = any(bool_col(x$nominal_significance_change)),
    stringsAsFactors = FALSE
  )
}))
single_summary <- do.call(rbind, lapply(split(single, single$analysis_set), function(x) {
  i <- which.max(abs(as.numeric(x$beta)))
  data.frame(
    analysis_set = x$analysis_set[[1L]],
    max_abs_single_snp_beta = abs(as.numeric(x$beta[[i]])),
    max_abs_single_snp_rsid = x$rsid[[i]],
    nominal_p_lt_0_05_count = sum(as.numeric(x$pval) < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

comp_row <- function(domain, set, method = NULL) {
  x <- comparison[comparison$comparison_domain == domain & comparison$analysis_set == set, , drop = FALSE]
  if (!is.null(method)) x <- x[x$method == method, , drop = FALSE]
  stop_if(nrow(x) != 1L, paste("Comparison row not found exactly once:", domain, set, method))
  x
}
strict_inc_comp <- comp_row("strict_included", "included", "Wald ratio")
strict_exc_comp <- comp_row("strict_excluded", "excluded", "Wald ratio")
relaxed_inc_ivw_comp <- comp_row("relaxed_estimator", "relaxed_apoe_included", "Inverse variance weighted")
relaxed_exc_ivw_comp <- comp_row("relaxed_estimator", "relaxed_apoe_excluded", "Inverse variance weighted")

strict_included_robustness <- classify_pair(strict_inc_comp$vuckovic_value, strict_inc_comp$chen_value, strict_inc_comp$vuckovic_pval, strict_inc_comp$chen_pval)
strict_excluded_robustness <- classify_pair(strict_exc_comp$vuckovic_value, strict_exc_comp$chen_value, strict_exc_comp$vuckovic_pval, strict_exc_comp$chen_pval)
relaxed_included_robustness <- classify_pair(relaxed_inc_ivw_comp$vuckovic_value, relaxed_inc_ivw_comp$chen_value, relaxed_inc_ivw_comp$vuckovic_pval, relaxed_inc_ivw_comp$chen_pval)
relaxed_excluded_robustness <- classify_pair(relaxed_exc_ivw_comp$vuckovic_value, relaxed_exc_ivw_comp$chen_value, relaxed_exc_ivw_comp$vuckovic_pval, relaxed_exc_ivw_comp$chen_pval)
overall_pattern <- if (
  strict_included_robustness == "same_direction_both_non_nominal" &&
    strict_excluded_robustness == "same_direction_both_non_nominal" &&
    relaxed_excluded_robustness == "same_direction_nominal_significance_not_reproduced_in_alternative_outcome_sensitivity"
) {
  "broadly_consistent_with_no_clear_strict_reverse_causal_evidence; relaxed_APOE_excluded_nominal_signal_not_statistically_robust_across_Hb_outcome_GWAS"
} else {
  "predefined_neutral_outcome_sensitivity_pattern"
}

doubling_strict_ok <- num_equal(strict$beta_doubling, strict$beta * log(2), tol = 1e-12) &&
  num_equal(strict$se_doubling, strict$se * log(2), tol = 1e-12) &&
  num_equal(strict$ci_lower_doubling, strict$ci_lower * log(2), tol = 1e-12) &&
  num_equal(strict$ci_upper_doubling, strict$ci_upper * log(2), tol = 1e-12) &&
  num_equal(strict$pval_doubling, strict$pval, tol = 0)
doubling_relaxed_ok <- num_equal(relaxed$beta_doubling, relaxed$beta * log(2), tol = 1e-12) &&
  num_equal(relaxed$se_doubling, relaxed$se * log(2), tol = 1e-12) &&
  num_equal(relaxed$ci_lower_doubling, relaxed$ci_lower * log(2), tol = 1e-12) &&
  num_equal(relaxed$ci_upper_doubling, relaxed$ci_upper * log(2), tol = 1e-12) &&
  num_equal(relaxed$pval_doubling, relaxed$pval, tol = 0)
no_or_transform <- !any(grepl("^OR$|^OR_|_OR$|odds_ratio|exp\\(|g/dL|g/L", c(names(strict), names(relaxed), unlist(strict), unlist(relaxed)), ignore.case = TRUE))

renv_after <- hash_file(paths$renv_lock)
hard_checks <- list(
  decision_116_gate = identical(input_freeze$freeze_status, "passed") &&
    isTRUE(input_freeze$approved_for_chen_reverse_mr_design) &&
    is_empty(input_freeze$hard_check_failures) &&
    identical(input_freeze$manifest_sha256, hash_file(paths$input_freeze_manifest)),
  decision_117_gate = identical(contract$contract_status, "frozen") &&
    isTRUE(contract$approved_for_chen_reverse_mr_execution) &&
    is_empty(contract$hard_check_failures),
  decision_118_gate = identical(mr_qc$mr_status, "passed") &&
    isTRUE(mr_qc$approved_for_chen_reverse_results_interpretation) &&
    is_empty(mr_qc$hard_check_failures),
  no_mr_rerun = no_forbidden_calls,
  strict_results_reverified = nrow(strict) == 2L &&
    all(strict$method == "Wald ratio") &&
    all(strict$method_id == "mr_wald_ratio") &&
    all(strict$nsnp == 1L) &&
    all(strict$effect_scale == "standardized_quantitative_Hb_per_1_unit_genetically_predicted_log_odds_delirium") &&
    all(strict$heterogeneity_status == "not_applicable_single_instrument") &&
    all(strict$single_snp_diagnostic_status == "not_separately_run_wald_is_formal_estimate"),
  relaxed_results_reverified = nrow(relaxed) == 10L &&
    setequal(unique(relaxed$method_id), c("mr_ivw", "mr_weighted_median", "mr_egger_regression", "mr_weighted_mode", "mr_simple_mode")) &&
    all(relaxed$nsnp[relaxed$analysis_set == "relaxed_apoe_included"] == 10L) &&
    all(relaxed$nsnp[relaxed$analysis_set == "relaxed_apoe_excluded"] == 9L),
  heterogeneity_reverified = nrow(heterogeneity) == 4L &&
    all(c("Q", "df", "pval") %in% names(heterogeneity)) &&
    all(!bool_col(heterogeneity$automatic_snp_removal_allowed)),
  egger_intercept_reverified = nrow(egger) == 2L &&
    all(egger$egger_precision_limitation == "limited_number_of_instruments"),
  mr_presso_reverified = nrow(presso) >= 4L &&
    all(c("relaxed_apoe_included", "relaxed_apoe_excluded") %in% unique(presso$analysis_set)) &&
    all(presso$status[presso$test_type == "Global Test"] == "passed"),
  loo_reverified = nrow(loo) == 19L &&
    nrow(loo_summary) == 2L &&
    all(c("max_absolute_shift", "max_shift_rsid", "any_sign_change", "any_nominal_significance_change") %in% names(loo_summary)),
  single_snp_reverified = nrow(single) == 19L &&
    nrow(single_summary) == 2L &&
    all(single$diagnostic_role == "single_snp_diagnostic_only"),
  vuckovic_comparison_reverified = nrow(comparison) >= 18L &&
    all(!bool_col(comparison$independent_replication)) &&
    all(!bool_col(comparison$formal_meta_analysis)) &&
    all(!bool_col(comparison$difference_test_performed)),
  strict_relaxed_hierarchy_preserved = identical(mr_qc$strict_primary_superseded_by_relaxed, FALSE) &&
    identical(mr_qc$relaxed_confirmatory, FALSE),
  independent_replication_false = identical(mr_qc$independent_replication, FALSE),
  effect_scale_preserved = identical(mr_qc$analysis_direction, "delirium_to_Hb") &&
    all(grepl("standardized_quantitative_Hb", c(strict$effect_scale, relaxed$effect_scale), fixed = TRUE)),
  doubling_odds_transform_preserved = doubling_strict_ok && doubling_relaxed_ok,
  no_or_transform = no_or_transform,
  no_posthoc_filtering = isTRUE(mr_qc$hard_checks$no_posthoc_filtering),
  no_proxy = isTRUE(mr_qc$hard_checks$no_proxy),
  no_liftover = isTRUE(mr_qc$hard_checks$no_liftover),
  no_steiger = isFALSE(mr_qc$steiger_run),
  renv_lock_unchanged = identical(renv_before, renv_after) && isTRUE(mr_qc$renv_lock_unchanged),
  git_status_not_required = identical(mr_qc$git_repository_present, FALSE) &&
    identical(mr_qc$git_status, "not_applicable_project_not_git_repository")
)
failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
freeze_status <- if (length(failures) == 0L) "passed" else "failed"
approved_unified <- identical(freeze_status, "passed")

freeze <- list(
  freeze_version = "v1",
  date = "2026-08-13",
  authoritative_chen_reverse_mr_version = "v1",
  analysis_direction = "delirium_to_Hb",
  analysis_role = "reverse_alternative_hb_outcome_sensitivity",
  independent_replication = FALSE,
  source_mr_input_freeze_decision = 116,
  mr_contract_decision = 117,
  mr_execution_decision = 118,
  strict_threshold = 5e-8,
  relaxed_threshold = 5e-6,
  strict_results = records(strict),
  relaxed_results = records(relaxed),
  heterogeneity_results = records(heterogeneity),
  egger_intercept_results = records(egger),
  mr_presso_results = records(presso),
  leave_one_out_summary = records(loo_summary),
  single_snp_summary = records(single_summary),
  vuckovic_comparison = records(comparison),
  robustness_summary = list(
    strict_included_robustness = strict_included_robustness,
    strict_excluded_robustness = strict_excluded_robustness,
    relaxed_included_robustness = relaxed_included_robustness,
    relaxed_excluded_robustness = relaxed_excluded_robustness,
    overall_reverse_outcome_sensitivity_pattern = overall_pattern,
    independent_replication = FALSE,
    formal_meta_analysis = FALSE,
    difference_test_performed = FALSE
  ),
  effect_scale = list(
    exposure_scale = "log_odds_delirium",
    outcome_scale = "standardized_quantitative_Hb_effect",
    raw_interpretation = "standardized quantitative Hb change per 1-unit genetically predicted log odds of delirium",
    secondary_interpretation = "standardized quantitative Hb change per doubling genetically predicted odds of delirium",
    no_or_transform = TRUE,
    no_physical_unit_claim = TRUE
  ),
  strict_primary_superseded_by_relaxed = FALSE,
  relaxed_confirmatory = FALSE,
  steiger_run = FALSE,
  steiger_status = "deferred_to_unified_directionality_sensitivity_stage",
  software_environment = mr_qc$software_environment,
  seed = mr_qc$seed,
  git_repository_present = mr_qc$git_repository_present,
  git_status = mr_qc$git_status,
  manifest_path = relpath(paths$manifest),
  manifest_sha256 = manifest_sha,
  freeze_status = freeze_status,
  approved_for_unified_directionality_design = approved_unified,
  hard_checks = hard_checks,
  hard_check_failures = failures,
  informational_findings = list(
    decision_114_status_preserved = "failed_due_readback_consistency_technical_issue",
    decision_115_authority_preserved = "technical_readback_recovery_not_scientific_rerun",
    no_mr_rerun = TRUE,
    no_results_reinterpretation_as_independent_replication = TRUE,
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after,
    out_of_sync_message = "informational_only"
  )
)

decision_lines <- c(
  "# Decision 119: Chen Reverse MR V1 Results Freeze",
  "",
  "Date: 2026-08-13",
  "",
  "## Status",
  paste0("freeze_status: `", freeze_status, "`"),
  paste0("approved_for_unified_directionality_design: `", if (approved_unified) "TRUE" else "FALSE", "`"),
  "",
  "## Decision",
  "Freeze Decision 118 Chen Reverse MR V1 results as the authoritative Chen reverse alternative-Hb-outcome sensitivity MR result set.",
  "",
  "This freeze did not rerun Wald ratio, IVW, MR-Egger, weighted median/mode, heterogeneity, Egger intercept, MR-PRESSO, leave-one-out, single-SNP, Steiger, harmonisation, clumping, proxy, or liftOver.",
  "",
  "## Authority Gates",
  "- Decision 116 input freeze: passed.",
  "- Decision 117 MR contract: frozen and approved.",
  "- Decision 118 MR execution: passed and approved for interpretation.",
  paste0("- Decision 116 manifest SHA-256: `", input_freeze$manifest_sha256, "`."),
  "",
  "## Strict Results",
  paste0("- included `", strict_included$rsid, "`: beta=`", signif(strict_included$beta, 6), "`, 95% CI `", signif(strict_included$ci_lower, 6), " to ", signif(strict_included$ci_upper, 6), "`, P=`", signif(strict_included$pval, 6), "`; doubling beta=`", signif(strict_included$beta_doubling, 6), "`."),
  paste0("- excluded `", strict_excluded$rsid, "`: beta=`", signif(strict_excluded$beta, 6), "`, 95% CI `", signif(strict_excluded$ci_lower, 6), " to ", signif(strict_excluded$ci_upper, 6), "`, P=`", signif(strict_excluded$pval, 6), "`; doubling beta=`", signif(strict_excluded$beta_doubling, 6), "`."),
  "",
  "## Relaxed IVW Results",
  paste0("- included: beta=`", signif(relaxed_ivw_included$beta, 6), "`, 95% CI `", signif(relaxed_ivw_included$ci_lower, 6), " to ", signif(relaxed_ivw_included$ci_upper, 6), "`, P=`", signif(relaxed_ivw_included$pval, 6), "`; doubling beta=`", signif(relaxed_ivw_included$beta_doubling, 6), "`."),
  paste0("- excluded: beta=`", signif(relaxed_ivw_excluded$beta, 6), "`, 95% CI `", signif(relaxed_ivw_excluded$ci_lower, 6), " to ", signif(relaxed_ivw_excluded$ci_upper, 6), "`, P=`", signif(relaxed_ivw_excluded$pval, 6), "`; doubling beta=`", signif(relaxed_ivw_excluded$beta_doubling, 6), "`."),
  "",
  "## Diagnostics",
  paste0("- Heterogeneity rows: `", nrow(heterogeneity), "`."),
  paste0("- Egger intercept rows: `", nrow(egger), "`."),
  paste0("- MR-PRESSO Global P values: `", paste(presso$pval[presso$test_type == "Global Test"], collapse = "; "), "`; outlier counts `", paste(presso$outlier_count[presso$test_type == "Outlier Test"], collapse = "; "), "`."),
  paste0("- LOO max shift rsIDs: `", paste(loo_summary$analysis_set, loo_summary$max_shift_rsid, sep = "=", collapse = "; "), "`."),
  paste0("- Single-SNP max |beta| rsIDs: `", paste(single_summary$analysis_set, single_summary$max_abs_single_snp_rsid, sep = "=", collapse = "; "), "`."),
  "",
  "## Vuckovic vs Chen Robustness",
  paste0("- strict included: `", strict_included_robustness, "`."),
  paste0("- strict excluded: `", strict_excluded_robustness, "`."),
  paste0("- relaxed included: `", relaxed_included_robustness, "`."),
  paste0("- relaxed excluded: `", relaxed_excluded_robustness, "`."),
  paste0("- overall pattern: `", overall_pattern, "`."),
  "- independent_replication: `FALSE`.",
  "",
  "## Audit",
  paste0("- manifest SHA-256: `", manifest_sha, "`."),
  paste0("- renv.lock SHA before/after: `", renv_before, "` / `", renv_after, "`."),
  paste0("- git status: `", mr_qc$git_status, "`."),
  "- Steiger: `deferred_to_unified_directionality_sensitivity_stage`.",
  "",
  "## Hard Check Failures",
  if (length(failures) == 0L) "- none" else paste0("- `", failures, "`"),
  "",
  "## Outputs Created",
  "- `R/33a_chen_reverse_mr_v1_freeze_manifest.R`",
  "- `results/qc/chen_reverse_mr_v1_freeze_manifest.csv`",
  "- `results/qc/chen_reverse_mr_v1_freeze.json`",
  "- `results/logs/chen_reverse_mr_v1_freeze.log`",
  "- `docs/decisions/119_chen_reverse_mr_v1_freeze_v1.1.md`",
  "",
  "## Next Gate",
  "Stop here. The next separate stage is Unified Directionality / Steiger Sensitivity Framework design."
)

log_lines <- c(
  "[2026-08-13] Chen Reverse MR V1 results freeze",
  paste0("freeze_status=", freeze_status),
  paste0("manifest_sha256=", manifest_sha),
  paste0("approved_for_unified_directionality_design=", approved_unified),
  paste0("hard_check_failures=", if (length(failures) == 0L) "[]" else paste(failures, collapse = ";")),
  paste0("overall_reverse_outcome_sensitivity_pattern=", overall_pattern),
  "mr_rerun=FALSE",
  "steiger_run=FALSE"
)

stop_if(length(failures) > 0L, paste("Freeze hard checks failed:", paste(failures, collapse = "; ")))
if (!file.rename(partial_manifest, paths$manifest)) stop("Atomic rename failed: ", paths$manifest, call. = FALSE)
write_json(freeze, paths$freeze_json)
write_text(log_lines, paths$freeze_log)
write_text(decision_lines, paths$decision)

cat("Decision 119 Chen Reverse MR V1 Results Freeze completed\n")
cat("freeze_status=", freeze_status, "\n", sep = "")
cat("manifest_sha256=", manifest_sha, "\n", sep = "")
cat("approved_for_unified_directionality_design=", approved_unified, "\n", sep = "")
cat("hard_check_failures=[]\n")
