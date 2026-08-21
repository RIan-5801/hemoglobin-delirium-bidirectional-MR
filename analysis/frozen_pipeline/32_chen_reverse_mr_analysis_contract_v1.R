#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

root <- "E:/Research/hb_delirium_bidir_mr"
setwd(root)

for (pkg in c("jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing package: ", pkg, call. = FALSE)
  }
}

rel <- function(...) file.path(...)
norm_path <- function(x) gsub("\\\\", "/", x)
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
is_empty <- function(x) is.null(x) || length(x) == 0L
`%||%` <- function(x, y) if (is.null(x)) y else x

write_json_atomic <- function(x, path) {
  partial <- paste0(path, ".partial")
  jsonlite::write_json(x, partial, pretty = TRUE, auto_unbox = TRUE, null = "null")
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}

write_csv_atomic <- function(x, path) {
  partial <- paste0(path, ".partial")
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.csv(x, partial, row.names = FALSE, na = "")
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}

write_text_atomic <- function(lines, path) {
  partial <- paste0(path, ".partial")
  writeLines(lines, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}

as_branch <- function(freeze, branch_name) {
  hits <- Filter(function(x) identical(x$branch, branch_name), freeze$branch_results)
  if (length(hits) != 1L) stop("Branch not found exactly once: ", branch_name, call. = FALSE)
  hits[[1]]
}

method_ids <- function(estimator_hierarchy) {
  ids <- vapply(
    estimator_hierarchy,
    function(x) if (!is.null(x$method_id)) x$method_id else NA_character_,
    character(1)
  )
  ids[!is.na(ids)]
}

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
if (!identical(next_decision, 117L)) {
  stop("Expected next decision 117, found ", next_decision, "; no outputs written.", call. = FALSE)
}

paths <- list(
  freeze = rel("results", "qc", "chen_reverse_harmonised_mr_inputs_v1_freeze.json"),
  freeze_manifest = rel("results", "qc", "chen_reverse_harmonised_mr_inputs_v1_freeze_manifest.csv"),
  strict_contract = rel("results", "qc", "reverse_strict_primary_mr_analysis_contract_v1.json"),
  strict_freeze = rel("results", "qc", "reverse_strict_primary_mr_v1_freeze.json"),
  relaxed_contract = rel("results", "qc", "reverse_relaxed_mr_analysis_contract_v1.json"),
  relaxed_freeze = rel("results", "qc", "reverse_relaxed_mr_v1_freeze.json"),
  renv_lock = rel("renv.lock"),
  contract = rel("results", "qc", "chen_reverse_mr_analysis_contract_v1.json"),
  drift = rel("results", "qc", "chen_reverse_mr_method_drift_audit_v1.csv"),
  log = rel("results", "logs", "chen_reverse_mr_analysis_contract_v1.log"),
  decision = rel("docs", "decisions", "117_chen_reverse_mr_analysis_contract_v1_v1.1.md")
)

inputs <- unlist(paths[c(
  "freeze", "freeze_manifest", "strict_contract", "strict_freeze",
  "relaxed_contract", "relaxed_freeze", "renv_lock"
)])
missing <- inputs[!file.exists(inputs)]
if (length(missing) > 0L) {
  stop("Missing required input(s): ", paste(missing, collapse = "; "), call. = FALSE)
}

targets <- unlist(paths[c("contract", "drift", "log", "decision")])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
if (length(occupied) > 0L) {
  stop("Target or partial exists: ", paste(occupied, collapse = "; "), call. = FALSE)
}

freeze <- read_json(paths$freeze)
strict_contract <- read_json(paths$strict_contract)
strict_freeze <- read_json(paths$strict_freeze)
relaxed_contract <- read_json(paths$relaxed_contract)
relaxed_freeze <- read_json(paths$relaxed_freeze)

manifest_sha <- hash_file(paths$freeze_manifest)
renv_before <- hash_file(paths$renv_lock)
renv_after <- hash_file(paths$renv_lock)

branches <- list(
  strict_apoe_included = as_branch(freeze, "strict_apoe_included"),
  strict_apoe_excluded = as_branch(freeze, "strict_apoe_excluded"),
  relaxed_apoe_included = as_branch(freeze, "relaxed_apoe_included"),
  relaxed_apoe_excluded = as_branch(freeze, "relaxed_apoe_excluded")
)

strict_branch_plan <- lapply(branches[c("strict_apoe_included", "strict_apoe_excluded")], function(b) {
  n <- b$final_valid_count
  if (identical(n, 0L) || identical(n, 0)) {
    estimator <- "not_estimable"
    mr_estimable <- FALSE
    reason <- "no_final_valid_instrument"
    stop_for_review <- FALSE
  } else if (identical(n, 1L) || identical(n, 1)) {
    estimator <- "Wald_ratio"
    mr_estimable <- TRUE
    reason <- NULL
    stop_for_review <- FALSE
  } else {
    estimator <- "human_review_required"
    mr_estimable <- NA
    reason <- "strict_multi_iv_design_not_authorized"
    stop_for_review <- TRUE
  }
  list(
    branch = b$branch,
    role = b$analysis_role,
    final_valid_count = n,
    final_rsids_source = "Decision_116_freeze_json_dynamic_read",
    final_rsids = b$final_rsids,
    mr_estimable = mr_estimable,
    estimator_class = estimator,
    reason = reason,
    strict_multi_iv_stop_for_human_review = stop_for_review,
    diagnostics = list(
      heterogeneity = "not_estimable_single_instrument_when_n1",
      MR_Egger = "not_estimable_single_instrument_when_n1",
      Egger_intercept = "not_estimable_single_instrument_when_n1",
      weighted_median = "not_applicable_single_instrument_when_n1",
      mode = "not_applicable_single_instrument_when_n1",
      MR_PRESSO = "not_applicable_single_instrument_when_n1",
      leave_one_out = "not_run_single_instrument_when_n1",
      single_snp_diagnostic = "not_separately_counted_for_strict_n1"
    )
  )
})

relaxed_branch_plan <- lapply(branches[c("relaxed_apoe_included", "relaxed_apoe_excluded")], function(b) {
  list(
    branch = b$branch,
    role = b$analysis_role,
    final_valid_count = b$final_valid_count,
    final_rsids_source = "Decision_116_freeze_json_dynamic_read",
    final_rsids = b$final_rsids,
    estimator_class = "multi_IV_exploratory_MR",
    primary_estimator_within_exploratory_branch = "IVW",
    sensitivity_estimators = list("Weighted median", "MR-Egger"),
    supportive_estimators = list("Weighted mode", "Simple mode"),
    diagnostics = list("heterogeneity", "Egger intercept", "MR-PRESSO", "leave-one-out", "single-SNP Wald diagnostics")
  )
})

relaxed_ids <- method_ids(relaxed_contract$estimator_hierarchy)
expected_relaxed_ids <- c(
  "mr_ivw", "mr_weighted_median", "mr_egger_regression",
  "mr_weighted_mode", "mr_simple_mode"
)

hard_checks <- list(
  decision_116_input_freeze_gate = identical(freeze$freeze_status, "passed") &&
    isTRUE(freeze$approved_for_chen_reverse_mr_design) &&
    is_empty(freeze$hard_check_failures) &&
    identical(tolower(freeze$manifest_sha256), tolower(manifest_sha)) &&
    identical(freeze$authoritative_chen_reverse_harmonisation_version, "v1"),
  decision_114_failure_provenance_preserved = TRUE,
  decision_115_recovery_authority_preserved = TRUE,
  analysis_role_sensitivity_not_primary = identical(freeze$analysis_role, "reverse_alternative_hb_outcome_sensitivity"),
  independent_replication_false = identical(freeze$independent_replication, FALSE),
  strict_relaxed_hierarchy_preserved = identical(freeze$strict_threshold, 5e-8) &&
    identical(freeze$relaxed_threshold, 5e-6) &&
    identical(freeze$strict_primary_superseded_by_relaxed, FALSE) &&
    identical(freeze$relaxed_confirmatory, FALSE),
  strict_method_authority_found = identical(strict_contract$contract_status, "frozen") &&
    isTRUE(strict_contract$approved_for_reverse_strict_primary_mr_execution) &&
    is_empty(strict_contract$hard_check_failures) &&
    identical(strict_contract$estimator, "Wald_ratio") &&
    identical(strict_contract$wald_implementation$authoritative_method, "TwoSampleMR::mr_wald_ratio") &&
    identical(strict_freeze$freeze_status, "passed") &&
    is_empty(strict_freeze$hard_check_failures),
  relaxed_method_authority_found = identical(relaxed_contract$contract_status, "frozen") &&
    isTRUE(relaxed_contract$approved_for_reverse_relaxed_mr_execution) &&
    is_empty(relaxed_contract$hard_check_failures) &&
    identical(relaxed_freeze$freeze_status, "passed") &&
    is_empty(relaxed_freeze$hard_check_failures),
  strict_n1_wald_rule_defined = all(vapply(strict_branch_plan, function(x) identical(x$estimator_class, "Wald_ratio"), logical(1))),
  strict_no_multi_iv_diagnostics_when_n1 = all(vapply(
    strict_branch_plan,
    function(x) identical(x$diagnostics$MR_PRESSO, "not_applicable_single_instrument_when_n1") &&
      identical(x$diagnostics$single_snp_diagnostic, "not_separately_counted_for_strict_n1"),
    logical(1)
  )),
  relaxed_ivw_primary_within_exploratory_branch = "mr_ivw" %in% relaxed_ids &&
    any(vapply(
      relaxed_contract$estimator_hierarchy,
      function(x) identical(x$method_id, "mr_ivw") &&
        identical(x$role, "exploratory_branch_primary_estimator"),
      logical(1)
    )),
  relaxed_estimator_hierarchy_matches_decision72 = identical(relaxed_ids, expected_relaxed_ids),
  heterogeneity_plan_matches_decision72 = length(relaxed_contract$heterogeneity_tests) >= 2L &&
    all(vapply(relaxed_contract$heterogeneity_tests, function(x) isFALSE(x$automatic_snp_removal_allowed), logical(1))),
  egger_plan_matches_decision72 = length(relaxed_contract$pleiotropy_tests) == 1L &&
    identical(relaxed_contract$pleiotropy_tests[[1]]$test, "MR_Egger_intercept") &&
    isFALSE(relaxed_contract$pleiotropy_tests[[1]]$automatic_snp_removal_allowed),
  mr_presso_plan_matches_decision72 = isTRUE(relaxed_contract$mr_presso_plan$planned) &&
    identical(relaxed_contract$mr_presso_plan$NbDistribution, 10000L) &&
    isTRUE(all.equal(relaxed_contract$mr_presso_plan$SignifThreshold, 0.05)) &&
    identical(relaxed_contract$mr_presso_plan$seed, 2026L) &&
    isFALSE(relaxed_contract$mr_presso_plan$automatic_main_input_redefinition_allowed),
  loo_plan_matches_decision72 = isTRUE(relaxed_contract$leave_one_out_plan$planned) &&
    isFALSE(relaxed_contract$leave_one_out_plan$automatic_main_analysis_redefinition_allowed),
  relaxed_single_snp_plan_matches_decision72 = isTRUE(relaxed_contract$single_snp_plan$planned) &&
    identical(relaxed_contract$single_snp_plan$role, "influence_diagnostic_visualization_only") &&
    isFALSE(relaxed_contract$single_snp_plan$main_result_use_allowed),
  effect_scale_defined = identical(freeze$outcome_scale, "standardized_quantitative_Hb_effect"),
  binary_exposure_interpretation_defined = TRUE,
  doubling_odds_scale_defined = TRUE,
  no_or_transform = TRUE,
  vuckovic_comparison_noninferential = TRUE,
  no_posthoc_filtering = TRUE,
  steiger_deferred = TRUE,
  software_environment_read_only = isTRUE(relaxed_contract$software_environment$read_only_probe_only) &&
    identical(strict_contract$software_environment$TwoSampleMR_version, "0.7.9") &&
    identical(relaxed_contract$software_environment$TwoSampleMR$version, "0.7.9") &&
    identical(relaxed_contract$software_environment$MRPRESSO$version, "1.0"),
  no_mr_executed = TRUE,
  renv_lock_unchanged = identical(renv_before, renv_after)
)

hard_check_failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
contract_status <- if (length(hard_check_failures) == 0L) "frozen" else "failed"
approved_for_execution <- identical(contract_status, "frozen")

contract <- list(
  contract_version = "v1",
  decision = 117,
  date = "2026-08-13",
  analysis_direction = "delirium_to_Hb",
  analysis_role = "reverse_alternative_hb_outcome_sensitivity",
  exposure = "FinnGen_R13_F5_DELIRIUM",
  outcome = "Chen_2020_Hb_BCX2",
  independent_replication = FALSE,
  source_mr_input_freeze_decision = 116,
  source_mr_input_freeze_json = norm_path(paths$freeze),
  source_mr_input_freeze_manifest = norm_path(paths$freeze_manifest),
  source_mr_input_freeze_manifest_sha256_expected = freeze$manifest_sha256,
  source_mr_input_freeze_manifest_sha256_observed = manifest_sha,
  harmonisation_provenance = list(
    decision_114 = "scientific harmonisation execution provenance; original run recorded failed due readback consistency technical issue",
    decision_115 = "technical readback recovery authority; did not rerun scientific harmonisation",
    decision_116 = "authoritative MR-input freeze"
  ),
  thresholds = list(
    strict_threshold = freeze$strict_threshold,
    relaxed_threshold = freeze$relaxed_threshold
  ),
  evidence_hierarchy = list(
    strict = list(
      instrument_threshold = 5e-8,
      role = "reverse_strict_primary_alternative_hb_outcome_sensitivity"
    ),
    relaxed = list(
      instrument_threshold = 5e-6,
      role = "reverse_relaxed_exploratory_alternative_hb_outcome_sensitivity",
      branch_type = "protocol_triggered_exploratory_fallback",
      confirmatory_status = "exploratory_not_confirmatory"
    ),
    strict_primary_superseded_by_relaxed = FALSE,
    relaxed_confirmatory = FALSE
  ),
  branches = list(
    strict = strict_branch_plan,
    relaxed = relaxed_branch_plan
  ),
  apoe_policy = list(
    included_excluded_from_independent_clumping_logic = TRUE,
    not_same_model_removing_one_apoe_snp = TRUE,
    not_independent_replication = TRUE
  ),
  exposure_scale = "log_odds_delirium",
  outcome_scale = "standardized_quantitative_Hb_effect",
  raw_effect_interpretation =
    "change in standardized quantitative Hb per 1-unit increase in genetically predicted log odds of delirium",
  forbidden_effect_interpretations = list(
    "clinical delirium causes Hb to change",
    "change in Hb in g/dL",
    "change in Hb in g/L",
    "odds ratio transformation for continuous Hb outcome"
  ),
  doubling_odds_rescaling_enabled = TRUE,
  doubling_odds_formula = list(
    beta_doubling = "beta_raw * log(2)",
    se_doubling = "se_raw * log(2)",
    ci_lower_doubling = "ci_lower_raw * log(2)",
    ci_upper_doubling = "ci_upper_raw * log(2)",
    p_value = "unchanged"
  ),
  doubling_odds_interpretation =
    "standardized quantitative Hb change per doubling of genetically predicted odds of delirium",
  strict_estimator_plan = list(
    source_contract_decision = 73,
    source_results_freeze_decision = 75,
    estimator = "Wald_ratio",
    method = "Wald ratio",
    method_id = strict_contract$wald_implementation$method_id,
    authoritative_implementation = strict_contract$wald_implementation$authoritative_method,
    function_body_sha256 = strict_contract$wald_implementation$function_body_sha256,
    software_formula = strict_contract$wald_implementation$software_formula,
    diagnostics_not_estimable_rules = strict_contract[c(
      "heterogeneity_status", "egger_status", "egger_intercept_status",
      "weighted_median_status", "mode_status", "mr_presso_status",
      "leave_one_out_status"
    )],
    n_rules = list(
      n0 = list(mr_estimable = FALSE, reason = "no_final_valid_instrument"),
      n1 = list(estimator_class = "Wald_ratio"),
      n_ge_2 = list(status = "stop_for_human_review", reason = "strict_multi_iv_design_not_authorized")
    )
  ),
  relaxed_primary_estimator = "IVW",
  relaxed_estimator_hierarchy = relaxed_contract$estimator_hierarchy,
  relaxed_sensitivity_estimators = relaxed_contract$sensitivity_estimators,
  heterogeneity_plan = relaxed_contract$heterogeneity_tests,
  egger_intercept_plan = c(
    relaxed_contract$pleiotropy_tests[[1]],
    list(egger_precision_limitation = "limited_number_of_instruments")
  ),
  mr_presso_plan = relaxed_contract$mr_presso_plan,
  mr_presso_failure_semantics = list(
    status_if_failed = "not_estimable_or_failed",
    exact_error_required = TRUE,
    fake_global_p_forbidden = TRUE,
    fake_no_outliers_forbidden = TRUE,
    reduce_NbDistribution_forbidden = TRUE
  ),
  leave_one_out_plan = relaxed_contract$leave_one_out_plan,
  single_snp_plan = list(
    relaxed = relaxed_contract$single_snp_plan,
    strict_n1 = "formal Wald estimator only; no separate diagnostic single-SNP analysis"
  ),
  strict_relaxed_interpretation_plan = list(
    relaxed_may_not_override_strict = TRUE,
    combined_or_meta_analysis_allowed = FALSE,
    strict_null_relaxed_p_lt_0_05_classification = "exploratory_discordance",
    broadly_consistent_wording_allowed_if_same_direction_or_similar = TRUE,
    proof_wording_allowed = FALSE
  ),
  vuckovic_comparison_plan = list(
    strict_vuckovic_authority_decision = 75,
    relaxed_vuckovic_authority_decision = 77,
    comparison_role = "robustness_to_alternative_Hb_outcome_GWAS",
    same_finngen_exposure_authority = TRUE,
    outcome_gwas_changed = TRUE,
    vuckovic_chen_independent_replication = FALSE,
    meta_analysis_allowed = FALSE,
    independent_replication_test_allowed = FALSE,
    formal_difference_test_allowed = FALSE
  ),
  posthoc_instrument_filtering_allowed = FALSE,
  forbidden_posthoc_filters = list(
    "F>=30 filtering",
    "MR P-based removal",
    "outcome P-based removal",
    "effect-direction removal",
    "heterogeneity-driven removal",
    "Egger-driven removal",
    "MR-PRESSO automatic main-input redefinition",
    "LOO-driven removal"
  ),
  steiger_run = FALSE,
  steiger_status = "deferred_to_unified_directionality_sensitivity_stage",
  seed = 2026,
  software_environment = list(
    read_only_probe_only = TRUE,
    install_update_restore_snapshot_performed = FALSE,
    R_version = relaxed_contract$software_environment$R_version,
    frozen_mr_library = relaxed_contract$software_environment$frozen_mr_library,
    TwoSampleMR = relaxed_contract$software_environment$TwoSampleMR,
    MRPRESSO = relaxed_contract$software_environment$MRPRESSO,
    seed = 2026,
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after,
    renv_out_of_sync_message = "informational_only"
  ),
  method_drift_audit = norm_path(paths$drift),
  mr_run = FALSE,
  contract_status = contract_status,
  approved_for_chen_reverse_mr_execution = approved_for_execution,
  hard_checks = hard_checks,
  hard_check_failures = hard_check_failures
)

drift <- data.frame(
  parameter_domain = c(
    "analysis_direction", "exposure_authority", "outcome_gwas", "analysis_role",
    "independent_replication", "strict_threshold", "relaxed_threshold",
    "strict_relaxed_hierarchy", "apoe_sets", "strict_estimator",
    "strict_diagnostics", "relaxed_estimator_hierarchy", "ivw_implementation",
    "heterogeneity", "egger_intercept", "mr_presso", "leave_one_out",
    "single_snp", "effect_scale", "doubling_odds", "or_transform",
    "steiger", "posthoc_filtering", "software_environment", "seed",
    "result_values", "output_paths"
  ),
  vuckovic_reverse_authority = c(
    "delirium_to_Hb", "FinnGen_R13_F5_DELIRIUM", "Vuckovic_2020_Hb",
    "secondary reverse; strict primary plus relaxed exploratory fallback",
    "FALSE", "5e-8", "5e-6", "strict not superseded by relaxed; relaxed not confirmatory",
    "APOE included and APOE excluded from independent clumping logic",
    "Decision 73/75 TwoSampleMR Wald ratio for strict n=1",
    "not estimable for strict n=1", "Decision 72/77 IVW; weighted median; MR-Egger; weighted mode; simple mode",
    "Decision 72/77 mr_ivw", "Decision 72/77 Cochran Q for IVW and MR-Egger",
    "Decision 72/77 MR-Egger intercept", "Decision 72/77 MRPRESSO NbDistribution=10000 SignifThreshold=0.05",
    "Decision 72/77 IVW leave-one-out", "Decision 72/77 relaxed single-SNP Wald diagnostic only",
    "standardized inverse-normal-transformed haemoglobin outcome", "beta/se/CI multiplied by log(2); P unchanged",
    "not allowed for continuous Hb outcome", "deferred", "not allowed",
    "R 4.6.1; renv/mr-v1-library; TwoSampleMR 0.7.9; MRPRESSO 1.0", "2026",
    "Vuckovic frozen result values", "Vuckovic reverse outputs"
  ),
  chen_reverse_contract_v1 = c(
    "delirium_to_Hb", "FinnGen_R13_F5_DELIRIUM", "Chen_2020_Hb_BCX2",
    "reverse_alternative_hb_outcome_sensitivity",
    "FALSE", "5e-8", "5e-6", "strict not superseded by relaxed; relaxed not confirmatory",
    "APOE included and APOE excluded from independent clumping logic",
    "reuse Decision 73/75 TwoSampleMR Wald ratio for strict n=1",
    "not estimable for strict n=1", "reuse Decision 72/77 IVW; weighted median; MR-Egger; weighted mode; simple mode",
    "reuse Decision 72/77 mr_ivw", "reuse Decision 72/77 Cochran Q for IVW and MR-Egger",
    "reuse Decision 72/77 MR-Egger intercept", "reuse Decision 72/77 MRPRESSO NbDistribution=10000 SignifThreshold=0.05",
    "reuse Decision 72/77 IVW leave-one-out", "reuse Decision 72/77 relaxed single-SNP Wald diagnostic only",
    "standardized quantitative Chen Hb effect outcome", "beta/se/CI multiplied by log(2); P unchanged",
    "not allowed for continuous Hb outcome", "deferred_to_unified_directionality_sensitivity_stage", "not allowed",
    "R 4.6.1; renv/mr-v1-library; TwoSampleMR 0.7.9; MRPRESSO 1.0", "2026",
    "future Chen result values; no MR executed in contract", "Chen reverse contract outputs"
  ),
  drift_status = c(
    "preserved", "preserved", "allowed_difference", "allowed_difference",
    "preserved", "preserved", "preserved", "preserved",
    "preserved", "preserved", "preserved", "preserved", "preserved",
    "preserved", "preserved", "preserved", "preserved",
    "preserved", "allowed_scale_wording_difference", "preserved",
    "preserved", "preserved", "preserved", "preserved", "preserved",
    "allowed_difference", "allowed_difference"
  ),
  stringsAsFactors = FALSE
)

branch_summary <- data.frame(
  branch = names(branches),
  final_valid_count = vapply(branches, function(x) x$final_valid_count, numeric(1)),
  final_rsids = vapply(branches, function(x) x$final_rsids, character(1)),
  analysis_role = vapply(branches, function(x) x$analysis_role, character(1)),
  estimator_plan = c("Wald ratio", "Wald ratio", "IVW primary within exploratory branch", "IVW primary within exploratory branch"),
  stringsAsFactors = FALSE
)

decision_lines <- c(
  "# Decision 117: Chen Reverse MR Analysis Contract V1",
  "",
  "Date: 2026-08-13",
  "",
  "## Status",
  sprintf("contract_status: `%s`", contract_status),
  sprintf("approved_for_chen_reverse_mr_execution: `%s`", if (approved_for_execution) "TRUE" else "FALSE"),
  "",
  "## Decision",
  "Freeze the Chen reverse alternative-Hb-outcome sensitivity MR analysis contract using Decision 116 as the authoritative MR-input freeze.",
  "",
  "This contract does not run Wald ratio, IVW, MR-Egger, weighted median, weighted/simple mode, heterogeneity, Egger intercept, MR-PRESSO, leave-one-out, single-SNP analysis, MR, or Steiger.",
  "",
  "## Authority",
  "- Decision 114: scientific harmonisation execution provenance; the original run remains recorded as failed due to a readback consistency technical issue.",
  "- Decision 115: technical readback recovery authority; it did not rerun scientific harmonisation.",
  "- Decision 116: authoritative Chen reverse harmonised MR-input freeze.",
  sprintf("- Decision 116 manifest SHA-256: `%s`.", manifest_sha),
  "- Strict method authority: Decision 73 contract and Decision 75 freeze.",
  "- Relaxed method authority: Decision 72 contract and Decision 77 freeze.",
  "",
  "## Analysis Identity",
  "- Direction: `delirium_to_Hb`.",
  "- Exposure: `FinnGen_R13_F5_DELIRIUM`.",
  "- Outcome: `Chen_2020_Hb_BCX2`.",
  "- Role: `reverse_alternative_hb_outcome_sensitivity`.",
  "- Independent replication: `FALSE`.",
  "",
  "## Dynamic Branch Counts From Decision 116",
  paste0("- `", branch_summary$branch, "`: `", branch_summary$final_valid_count, "`; plan `", branch_summary$estimator_plan, "`."),
  "",
  "## Evidence Hierarchy",
  "- Strict branch: `P < 5e-8`, `reverse_strict_primary_alternative_hb_outcome_sensitivity`.",
  "- Relaxed branch: `P < 5e-6`, protocol-triggered exploratory fallback.",
  "- `strict_primary_superseded_by_relaxed=FALSE`.",
  "- `relaxed_confirmatory=FALSE`.",
  "",
  "## Strict Estimator Plan",
  "- If strict branch `n=0`: not estimable, reason `no_final_valid_instrument`.",
  "- If strict branch `n=1`: formal estimator is TwoSampleMR Wald ratio, method ID `mr_wald_ratio`.",
  "- If strict branch `n>=2`: stop for human review because strict multi-IV design is not authorized in this contract.",
  "- For strict `n=1`, IVW, MR-Egger, weighted median/mode, heterogeneity, Egger intercept, MR-PRESSO, leave-one-out, and separate single-SNP diagnostic analysis are not run.",
  "",
  "## Relaxed Estimator Hierarchy",
  "- IVW is the primary estimator within the relaxed exploratory branch only.",
  "- Weighted median and MR-Egger are sensitivity estimators.",
  "- Weighted mode and simple mode are supportive estimators.",
  "- Heterogeneity, Egger intercept, MR-PRESSO, leave-one-out, and single-SNP Wald analyses are diagnostics only and must not redefine the frozen input.",
  "",
  "## Effect Scale",
  "- Raw beta: standardized quantitative Hb change per 1-unit increase in genetically predicted log odds of delirium.",
  "- Doubling-odds scale: beta, SE, and CI are multiplied by `log(2)`; P values are unchanged.",
  "- No g/dL or g/L claim is allowed.",
  "- No OR transform is allowed because Chen Hb is a continuous outcome.",
  "",
  "## Vuckovic Comparison Plan",
  "Future Chen results may be compared against frozen Vuckovic reverse strict and relaxed results as a noninferential robustness comparison for alternative Hb outcome GWAS source only. This is not an independent replication, meta-analysis, or formal difference test.",
  "",
  "## Steiger",
  "`steiger_run=FALSE`; `steiger_status=deferred_to_unified_directionality_sensitivity_stage`.",
  "",
  "## Software",
  sprintf("- R: `%s`.", relaxed_contract$software_environment$R_version),
  sprintf("- MR library: `%s`.", relaxed_contract$software_environment$frozen_mr_library),
  sprintf("- TwoSampleMR: `%s`, RemoteSha `%s`.", relaxed_contract$software_environment$TwoSampleMR$version, relaxed_contract$software_environment$TwoSampleMR$RemoteSha),
  sprintf("- MRPRESSO: `%s`, RemoteSha `%s`.", relaxed_contract$software_environment$MRPRESSO$version, relaxed_contract$software_environment$MRPRESSO$RemoteSha),
  "- Seed: `2026`.",
  sprintf("- renv.lock SHA before/after: `%s` / `%s`.", renv_before, renv_after),
  "",
  "## Hard Check Failures",
  if (length(hard_check_failures) == 0L) "- none" else paste0("- `", hard_check_failures, "`"),
  "",
  "## Outputs Created",
  "- `results/qc/chen_reverse_mr_analysis_contract_v1.json`",
  "- `results/qc/chen_reverse_mr_method_drift_audit_v1.csv`",
  "- `results/logs/chen_reverse_mr_analysis_contract_v1.log`",
  "- `docs/decisions/117_chen_reverse_mr_analysis_contract_v1_v1.1.md`",
  "",
  "## Next Gate",
  "Only if this contract remains frozen with `hard_check_failures=[]` may a later approved Decision execute Chen Reverse MR V1."
)

log_lines <- c(
  "Chen Reverse MR Analysis Contract V1",
  "Decision: 117",
  "Date: 2026-08-13",
  sprintf("Decision 116 freeze_status: %s", freeze$freeze_status),
  sprintf("Decision 116 approved_for_chen_reverse_mr_design: %s", freeze$approved_for_chen_reverse_mr_design),
  sprintf("Decision 116 manifest SHA observed: %s", manifest_sha),
  "Decision 114 provenance preserved as failed due readback technical issue.",
  "Decision 115 preserved as technical readback recovery authority.",
  "Strict authority recovered from Decision 73/75.",
  "Relaxed authority recovered from Decision 72/77.",
  paste(sprintf("%s n=%s", branch_summary$branch, branch_summary$final_valid_count), collapse = "; "),
  sprintf("contract_status: %s", contract_status),
  sprintf("approved_for_chen_reverse_mr_execution: %s", approved_for_execution),
  sprintf("hard_check_failures: %s", if (length(hard_check_failures) == 0L) "none" else paste(hard_check_failures, collapse = "; ")),
  "MR executed: FALSE",
  "Steiger executed: FALSE",
  sprintf("renv.lock SHA before: %s", renv_before),
  sprintf("renv.lock SHA after: %s", renv_after)
)

write_json_atomic(contract, paths$contract)
write_csv_atomic(drift, paths$drift)
write_text_atomic(log_lines, paths$log)
write_text_atomic(decision_lines, paths$decision)

cat("Decision 117 Chen Reverse MR Analysis Contract V1 completed\n")
cat("contract_status=", contract_status, "\n", sep = "")
cat("approved_for_chen_reverse_mr_execution=", approved_for_execution, "\n", sep = "")
cat("hard_check_failures=", if (length(hard_check_failures) == 0L) "none" else paste(hard_check_failures, collapse = ";"), "\n", sep = "")
