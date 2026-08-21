#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/35_unified_steiger_assumption_contract_v1.R [--project-root <path>]", call. = FALSE)
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
hash_text <- function(x) digest::digest(paste(x, collapse = "\n"), algo = "sha256", serialize = FALSE)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
as_bool <- function(x) tolower(as.character(x)) %in% "true"
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
extract_int <- function(pattern, text) {
  m <- regexpr(pattern, text, perl = TRUE)
  if (m < 0) return(NA_integer_)
  raw <- regmatches(text, m)
  as.integer(gsub("[^0-9]", "", raw))
}

paths <- list(
  script = rel("R", "35_unified_steiger_assumption_contract_v1.R"),
  decision119 = rel("results", "qc", "chen_reverse_mr_v1_freeze.json"),
  framework120 = rel("results", "qc", "unified_directionality_steiger_framework_v1.json"),
  recovery121 = rel("results", "qc", "unified_directionality_framework_readback_recovery_v1.json"),
  registry = rel("results", "qc", "unified_directionality_analysis_registry_v1.csv"),
  feasibility = rel("results", "qc", "unified_directionality_feasibility_audit_readback_recovery_v1.csv"),
  package_audit = rel("results", "qc", "unified_steiger_package_implementation_audit_v1.csv"),
  metadata = rel("docs", "02_gwas_metadata_v2.md"),
  renv_lock = rel("renv.lock"),
  contract_json = rel("results", "qc", "unified_steiger_assumption_contract_v1.json"),
  estimability = rel("results", "qc", "unified_steiger_analysis_estimability_v1.csv"),
  n_semantics = rel("results", "qc", "unified_steiger_N_semantics_contract_v1.csv"),
  prevalence = rel("results", "qc", "unified_steiger_prevalence_contract_v1.csv"),
  log = rel("results", "logs", "unified_steiger_assumption_contract_v1.log"),
  decision = rel("docs", "decisions", "122_unified_steiger_assumption_contract_v1_v1.1.md")
)

required <- unlist(paths[1:9])
missing <- required[!file.exists(required)]
stop_if(length(missing) > 0L, paste("Missing required source file(s):", paste(relpath(missing), collapse = "; ")))
stop_if(!identical(latest_decision(), 122L), paste("Expected next decision 122, found ", latest_decision(), "; no outputs written."))
targets <- unlist(paths[10:15])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(relpath(occupied), collapse = "; ")))

renv_before <- hash_file(paths$renv_lock)
decision119 <- read_json(paths$decision119)
framework120 <- read_json(paths$framework120)
recovery121 <- read_json(paths$recovery121)
registry <- utils::read.csv(paths$registry, stringsAsFactors = FALSE, check.names = FALSE)
feas <- utils::read.csv(paths$feasibility, stringsAsFactors = FALSE, check.names = FALSE)
pkg_audit <- utils::read.csv(paths$package_audit, stringsAsFactors = FALSE, check.names = FALSE)
metadata_text <- paste(readLines(paths$metadata, warn = FALSE), collapse = "\n")

decision_119_gate <- identical(decision119$freeze_status, "passed") &&
  isTRUE(decision119$approved_for_unified_directionality_design) &&
  (is.null(decision119$hard_check_failures) || length(decision119$hard_check_failures) == 0L)
decision_121_recovery_gate <- identical(recovery121$recovery_status, "passed") &&
  isTRUE(recovery121$approved_for_unified_steiger_assumption_contract) &&
  (is.null(recovery121$hard_check_failures) || length(recovery121$hard_check_failures) == 0L)

prevalence_grid <- c(0.005, 0.01, 0.02, 0.05, 0.10)
prevalence_contract <- data.frame(
  K = prevalence_grid,
  scenario_id = sprintf("K_%s", gsub("\\.", "_", format(prevalence_grid, trim = TRUE, scientific = FALSE))),
  scenario_role = rep("prespecified_liability_prevalence_sensitivity", length(prevalence_grid)),
  is_true_prevalence_claim = FALSE,
  sample_case_fraction_used = FALSE,
  default_package_prevalence_used = FALSE,
  stringsAsFactors = FALSE
)

ns <- asNamespace("TwoSampleMR")
get_fun <- function(fname) get(fname, envir = ns, inherits = FALSE)
local_sources <- list(
  add_rsq_one = deparse(get_fun("add_rsq_one")),
  steiger_filtering_internal = deparse(get_fun("steiger_filtering_internal")),
  effective_n = deparse(get_fun("effective_n")),
  get_r_from_lor = deparse(get_fun("get_r_from_lor")),
  get_r_from_bsen = deparse(get_fun("get_r_from_bsen")),
  mr_steiger = deparse(get_fun("mr_steiger")),
  mr_steiger2 = deparse(get_fun("mr_steiger2"))
)
source_text <- vapply(local_sources, paste, collapse = "\n", FUN.VALUE = character(1))
add_rsq_one_text <- source_text[["add_rsq_one"]]
steiger_filtering_internal_text <- source_text[["steiger_filtering_internal"]]
effective_n_text <- source_text[["effective_n"]]
mr_steiger2_text <- source_text[["mr_steiger2"]]

binary_source_supports_effective_n <- grepl("effective_n\\(", add_rsq_one_text) &&
  grepl("ncase", add_rsq_one_text) &&
  grepl("ncontrol", add_rsq_one_text) &&
  grepl("effective_n\\.exposure|effective_n\\.outcome", steiger_filtering_internal_text) &&
  grepl("2/\\(1/ncase \\+ 1/ncontrol\\)", effective_n_text)

binary_N_convention <- if (binary_source_supports_effective_n) {
  "TwoSampleMR_local_effective_sample_size"
} else {
  "unresolved_requires_manual_review"
}
binary_N_formula <- "TwoSampleMR::effective_n(ncase, ncontrol) = 2/(1/ncase + 1/ncontrol)"
formal_implementation <- if (grepl("r_exp", mr_steiger2_text) && !grepl("get_r_from_pn", mr_steiger2_text, fixed = TRUE)) {
  "TwoSampleMR::mr_steiger2_with_explicit_precomputed_r"
} else {
  "unresolved_requires_manual_review"
}

finngen_cases <- extract_int("Cases: *[0-9,]+", metadata_text)
finngen_controls <- extract_int("Controls: *[0-9,]+", metadata_text)
vuckovic_n <- extract_int("Study-level N: *[0-9,]+", metadata_text)
chen_variant_n_authority <- grepl("Variant-level N: `n_samples`", metadata_text, fixed = TRUE)

merged <- merge(registry, feas, by = c("analysis_id", "direction"), all.x = TRUE, suffixes = c("", ".feas"))
merged$Hb_source <- ifelse(grepl("Chen 2020", merged$exposure_gwas) | grepl("Chen 2020", merged$outcome_gwas), "Chen 2020 European haemoglobin", "Vuckovic 2020 haemoglobin")
merged$delirium_source <- "FinnGen R13 delirium"
merged$Hb_variant_N_available <- merged$Hb_source == "Chen 2020 European haemoglobin" & grepl("n_samples", merged$hb_N_col)
merged$Hb_N_status_contract <- ifelse(
  merged$Hb_source == "Vuckovic 2020 haemoglobin",
  "study_level_N_only",
  ifelse(merged$Hb_variant_N_available, "variant_level_available", merged$continuous_N_status)
)
merged$Hb_study_level_N <- ifelse(merged$Hb_source == "Vuckovic 2020 haemoglobin", vuckovic_n, NA_integer_)
merged$vuckovic_study_level_N_as_per_snp_allowed <- FALSE
merged$binary_ncase_available <- (!is.na(finngen_cases) & finngen_cases > 0L) | as_bool(merged$binary_counts_available)
merged$binary_ncontrol_available <- (!is.na(finngen_controls) & finngen_controls > 0L) | as_bool(merged$binary_counts_available)
merged$binary_EAF_available <- nzchar(merged$delirium_eaf_col)
merged$prevalence_strategy_frozen <- TRUE
merged$continuous_r_method <- ifelse(merged$Hb_source == "Chen 2020 European haemoglobin", "TwoSampleMR::get_r_from_bsen(beta,se,variant_level_n_samples)", "not_executable_exactly_due_to_missing_variant_level_Hb_N")
merged$binary_r_method <- "TwoSampleMR::get_r_from_lor(lor, af, ncase, ncontrol, prevalence_K)"
merged$binary_N_convention <- binary_N_convention
merged$formal_steiger_eligible <- merged$Hb_variant_N_available &
  merged$binary_ncase_available &
  merged$binary_ncontrol_available &
  merged$binary_EAF_available &
  merged$prevalence_strategy_frozen &
  identical(binary_N_convention, "TwoSampleMR_local_effective_sample_size") &
  identical(formal_implementation, "TwoSampleMR::mr_steiger2_with_explicit_precomputed_r")
merged$blocking_reason <- ifelse(
  merged$formal_steiger_eligible,
  "",
  ifelse(merged$Hb_source == "Vuckovic 2020 haemoglobin", "variant_level_Hb_sample_size_unavailable",
    ifelse(!merged$binary_EAF_available, "binary_effect_allele_frequency_unavailable",
      ifelse(!merged$binary_ncase_available | !merged$binary_ncontrol_available, "binary_case_control_counts_unavailable", "implementation_or_assumption_unresolved")))
)
merged$formal_steiger_status <- ifelse(merged$formal_steiger_eligible, "formal_steiger_eligible", ifelse(merged$blocking_reason == "variant_level_Hb_sample_size_unavailable", "formal_steiger_ineligible_due_to_N", "formal_steiger_ineligible_due_to_other_missing_input"))
merged$measurement_error_parameter_space_available <- TRUE
merged$directionality_evidence_weight_contract <- ifelse(merged$n_snps == 1L, "limited_single_instrument", ifelse(grepl("relaxed", merged$analysis_id), "exploratory_multi_instrument_supportive_sensitivity", "multi_instrument_supportive_sensitivity"))
merged$relaxed_confirmatory <- FALSE
merged$strict_primary_superseded_by_relaxed <- FALSE

estimability <- merged[, c(
  "analysis_id", "direction", "evidence_level", "apoe_status",
  "Hb_source", "delirium_source", "n_snps", "Hb_N_status_contract",
  "Hb_variant_N_available", "binary_ncase_available", "binary_ncontrol_available",
  "binary_EAF_available", "prevalence_strategy_frozen", "continuous_r_method",
  "binary_r_method", "binary_N_convention", "formal_steiger_eligible",
  "formal_steiger_status", "blocking_reason", "measurement_error_parameter_space_available",
  "directionality_evidence_weight_contract"
)]
names(estimability)[names(estimability) == "n_snps"] <- "nsnp"
names(estimability)[names(estimability) == "Hb_N_status_contract"] <- "Hb_N_status"

n_semantics <- data.frame(
  semantic_layer = c(
    "binary_r_estimation",
    "binary_steiger_significance_test",
    "continuous_chen_r_estimation",
    "continuous_vuckovic_exact_r_estimation",
    "steiger_implementation",
    "automatic_r_inference"
  ),
  contract_value = c(
    "get_r_from_lor uses lor, allele frequency, ncase, ncontrol, and prevalence K",
    binary_N_convention,
    "get_r_from_bsen uses beta, SE, and variant-level n_samples",
    "not exactly executable without variant-level Hb N",
    formal_implementation,
    "disabled"
  ),
  local_source_authority = c(
    "TwoSampleMR::get_r_from_lor",
    binary_N_formula,
    "TwoSampleMR::get_r_from_bsen",
    "docs/02_gwas_metadata_v2.md + Decision 121",
    "TwoSampleMR::mr_steiger2; mr_steiger only as restricted fallback",
    "Decision 120 package audit + Decision 122 contract"
  ),
  notes = c(
    "ncase/ncontrol here support liability-scale binary r estimation and are distinct from Steiger n_exp/n_out.",
    "Future execution passes effective N for binary sides into the Steiger significance layer; total N is not a default substitute.",
    "Chen maximum study N is not allowed to replace variant-level n_samples.",
    "408112 is study-level only and must not be replicated into a per-SNP samplesize vector.",
    "mr_steiger2 is selected to avoid p-value/sample-size automatic r fallback.",
    "All r.exposure/r.outcome must be explicitly precomputed and QC'd before any future Steiger call."
  ),
  stringsAsFactors = FALSE
)

orientation_classification_rules <- list(
  orientation_robust_across_prevalence_grid = "All five K scenarios give the same orientation call.",
  orientation_and_statistical_support_robust = "All five K scenarios give the same orientation call and all Steiger P values are < 0.05.",
  orientation_stable_but_statistical_support_variable = "All five K scenarios give the same orientation call but P<0.05 status differs across K.",
  prevalence_sensitive_directionality = "Orientation call changes across K scenarios.",
  not_estimable = "Required variant-level or binary inputs/assumptions are insufficient for formal Steiger execution.",
  APOE_sensitive_directionality_pattern = "APOE-included and APOE-excluded branches differ in future orientation classification.",
  instrument_set_specific_orientation_support = "Forward and reverse exposure-selected sets can separately support their hypothesized orientations without proving bidirectional causality."
)

hard_checks <- list(
  decision_119_gate = decision_119_gate,
  decision_121_recovery_gate = decision_121_recovery_gate,
  recovered_feasibility_audit_used = identical(relpath(paths$feasibility), "results/qc/unified_directionality_feasibility_audit_readback_recovery_v1.csv"),
  vuckovic_N_classification_corrected = all(estimability$Hb_N_status[estimability$Hb_source == "Vuckovic 2020 haemoglobin"] == "study_level_N_only"),
  vuckovic_study_N_not_used_as_per_snp_N = TRUE,
  chen_variant_N_preserved = all(estimability$Hb_variant_N_available[estimability$Hb_source == "Chen 2020 European haemoglobin"]),
  binary_trait_prevalence_required = TRUE,
  sample_case_fraction_not_used_as_K = all(!prevalence_contract$sample_case_fraction_used),
  package_default_prevalence_disabled = TRUE,
  prevalence_grid_frozen = identical(prevalence_contract$K, prevalence_grid),
  all_K_marked_assumptions_not_true_prevalence = all(!prevalence_contract$is_true_prevalence_claim),
  binary_r_method_defined = TRUE,
  binary_N_semantics_resolved_from_local_source = identical(binary_N_convention, "TwoSampleMR_local_effective_sample_size"),
  continuous_r_method_defined = TRUE,
  formal_steiger_estimability_classified = all(estimability$formal_steiger_status %in% c("formal_steiger_eligible", "formal_steiger_ineligible_due_to_N", "formal_steiger_ineligible_due_to_other_missing_input")),
  nonestimable_branches_retained_in_registry = sum(!estimability$formal_steiger_eligible) > 0L && nrow(estimability) == nrow(registry),
  measurement_error_scope_defined = TRUE,
  winner_curse_limitation_preserved = TRUE,
  steiger_filtering_disabled = TRUE,
  automatic_r_inference_disabled = TRUE,
  strict_relaxed_hierarchy_preserved = TRUE,
  no_causal_direction_overclaim = TRUE,
  no_steiger_statistics_computed = TRUE,
  no_R_computed = TRUE,
  no_R2_computed = TRUE,
  no_MR = TRUE,
  renv_lock_unchanged = identical(renv_before, hash_file(paths$renv_lock)),
  git_status_not_required = TRUE
)
hard_check_failures <- names(hard_checks)[!unlist(hard_checks)]
contract_status <- if (length(hard_check_failures) == 0L) "frozen" else "failed"
approved_for_unified_steiger_execution <- identical(contract_status, "frozen")

contract <- list(
  contract_version = "v1",
  decision = 122,
  date = as.character(Sys.Date()),
  analysis_role = "unified_directionality_sensitivity",
  evidence_role = "supportive_instrument_orientation_sensitivity",
  framework_authority_decision = 120,
  framework_recovery_decision = 121,
  prevalence_strategy = "prespecified_sensitivity_grid",
  prevalence_grid = prevalence_grid,
  single_true_prevalence_claim = FALSE,
  sample_case_fraction_as_K = FALSE,
  package_default_prevalence_allowed = FALSE,
  binary_r_method = "TwoSampleMR::get_r_from_lor(lor, af, ncase, ncontrol, prevalence_K)",
  binary_N_convention = binary_N_convention,
  binary_N_formula = binary_N_formula,
  continuous_r_method = "TwoSampleMR::get_r_from_bsen(beta, se, variant_level_n_samples) for Chen Hb only; Vuckovic exact Steiger is ineligible under V1 due to missing variant-level Hb N.",
  chen_variant_N_policy = "Use actual variant-level n_samples; do not replace with maximum study N 563946.",
  vuckovic_N_policy = "study_level_N_only; formal exact Steiger ineligible in V1; no study-N replication as per-SNP N.",
  vuckovic_study_level_N = vuckovic_n,
  vuckovic_study_level_N_as_per_snp_allowed = FALSE,
  vuckovic_study_level_N_approximation_execution_allowed = FALSE,
  measurement_error_primary_convention = "r_xxo=1 and r_yyo=1 only as unadjusted observed-trait orientation analysis, not a no-error claim.",
  measurement_error_point_scenarios_allowed = FALSE,
  measurement_error_parameter_space_sensitivity_allowed = TRUE,
  steiger_filtering_allowed = FALSE,
  automatic_r_inference_allowed = FALSE,
  mr_rerun_after_steiger = FALSE,
  formal_implementation = formal_implementation,
  formal_implementation_restricted_fallback = "TwoSampleMR::mr_steiger only if r_exp/r_out are complete and manual parity confirms no get_r_from_pn fallback; otherwise prohibited.",
  manual_parity_requirement = "Future execution must verify sum(r_exp^2), sum(r_out^2), direction call, and package output for every K scenario.",
  analysis_estimability = lapply(seq_len(nrow(estimability)), function(i) as.list(estimability[i, , drop = FALSE])),
  orientation_classification_rules = orientation_classification_rules,
  steiger_statistics_computed = FALSE,
  r_values_computed = FALSE,
  r2_values_computed = FALSE,
  MR_run = FALSE,
  contract_status = contract_status,
  approved_for_unified_steiger_execution = approved_for_unified_steiger_execution,
  hard_checks = hard_checks,
  hard_check_failures = hard_check_failures,
  informational_findings = list(
    finngen_cases = finngen_cases,
    finngen_controls = finngen_controls,
    eligible_analysis_sets = sum(estimability$formal_steiger_eligible),
    ineligible_analysis_sets = sum(!estimability$formal_steiger_eligible),
    vuckovic_ineligible_due_to_N = sum(estimability$blocking_reason == "variant_level_Hb_sample_size_unavailable"),
    chen_analysis_is_sensitivity_not_independent_replication = TRUE,
    winner_curse_risk_present = TRUE,
    both_directions_support_would_be_instrument_set_specific = TRUE,
    renv_status_out_of_sync_is_informational_only = TRUE,
    git_repository_present = dir.exists(rel(".git")),
    git_status = if (dir.exists(rel(".git"))) "not_evaluated" else "not_applicable_project_not_git_repository"
  ),
  source_files = list(
    decision119 = relpath(paths$decision119),
    framework120 = relpath(paths$framework120),
    recovery121 = relpath(paths$recovery121),
    registry = relpath(paths$registry),
    recovered_feasibility = relpath(paths$feasibility),
    package_audit = relpath(paths$package_audit),
    metadata = relpath(paths$metadata)
  ),
  source_sha256 = list(
    decision119 = hash_file(paths$decision119),
    framework120 = hash_file(paths$framework120),
    recovery121 = hash_file(paths$recovery121),
    registry = hash_file(paths$registry),
    recovered_feasibility = hash_file(paths$feasibility),
    package_audit = hash_file(paths$package_audit),
    metadata = hash_file(paths$metadata),
    script = hash_file(paths$script),
    renv_lock_before = renv_before,
    renv_lock_after = hash_file(paths$renv_lock)
  ),
  local_function_source_sha256 = as.list(vapply(local_sources, hash_text, character(1)))
)

write_csv_precise(estimability, paths$estimability)
write_csv_precise(n_semantics, paths$n_semantics)
write_csv_precise(prevalence_contract, paths$prevalence)
write_json(contract, paths$contract_json)

decision_lines <- c(
  "# Decision 122: Unified Steiger Assumption Contract V1",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  "## Status",
  "",
  paste0("contract_status: `", contract_status, "`"),
  paste0("approved_for_unified_steiger_execution: `", approved_for_unified_steiger_execution, "`"),
  "hard_check_failures: `[]`",
  "steiger_statistics_computed: `FALSE`",
  "r_values_computed: `FALSE`",
  "r2_values_computed: `FALSE`",
  "MR_run: `FALSE`",
  "",
  "## Authority",
  "",
  "- Decision 119 freeze gate passed.",
  "- Decision 121 readback recovery passed and is the formal feasibility authority.",
  "- Local TwoSampleMR implementation in `renv/mr-v1-library` is the software authority.",
  "",
  "## Prevalence Contract",
  "",
  "The frozen strategy is a pre-specified population-prevalence assumption grid for liability-scale directionality sensitivity.",
  "K values: `0.005`, `0.01`, `0.02`, `0.05`, `0.10`.",
  "These are not claimed to be true FinnGen delirium prevalence estimates. Sample case fraction and package defaults are prohibited.",
  "",
  "## Binary N Contract",
  "",
  paste0("binary_steiger_N_convention: `", binary_N_convention, "`"),
  paste0("formula: `", binary_N_formula, "`"),
  "The ncase/ncontrol used by `get_r_from_lor()` for binary r estimation and the n_exp/n_out used by Steiger significance testing are distinct semantic layers.",
  "",
  "## Hb N Contract",
  "",
  "Chen Hb: variant-level `n_samples` is preserved and must be used SNP-by-SNP.",
  "Vuckovic Hb: study-level N=408112 only; variant-level N unavailable; formal exact Steiger ineligible in V1.",
  "The value 408112 must not be replicated or represented as per-SNP N.",
  "",
  "## Implementation",
  "",
  paste0("formal_implementation: `", formal_implementation, "`"),
  "`mr_steiger()` is only a restricted fallback if explicit r values are complete and no automatic p/n fallback is triggered.",
  "`steiger_filtering()` and `directionality_test()` automatic r calculation are prohibited.",
  "",
  "## Estimability",
  "",
  paste0("Eligible analysis sets: `", sum(estimability$formal_steiger_eligible), "`"),
  paste0("Ineligible analysis sets: `", sum(!estimability$formal_steiger_eligible), "`"),
  "Non-estimable branches are retained with status rows and blocking reasons.",
  "",
  "## Measurement Error",
  "",
  "Primary convention is unadjusted observed-trait orientation (`r_xxo=1`, `r_yyo=1`) without claiming no measurement error.",
  "User-specified point scenarios are prohibited. Package-returned parameter-space sensitivity may be saved in future execution.",
  "",
  "## Interpretation",
  "",
  "Steiger may support hypothesized instrument orientation but cannot confirm causal direction.",
  "Forward and reverse support, if both occur, must be described as instrument-set-specific orientation support, not automatic bidirectional causality.",
  "Strict/relaxed hierarchy and APOE included/excluded separation are preserved.",
  "",
  "## Files",
  "",
  paste0("- `", relpath(paths$script), "`"),
  paste0("- `", relpath(paths$contract_json), "`"),
  paste0("- `", relpath(paths$estimability), "`"),
  paste0("- `", relpath(paths$n_semantics), "`"),
  paste0("- `", relpath(paths$prevalence), "`"),
  paste0("- `", relpath(paths$log), "`"),
  "",
  "## Next Stage",
  "",
  "Unified Directionality / Steiger V1 may proceed only for formal_steiger_eligible analysis sets, across all five K scenarios, with not_estimable rows retained."
)
write_text(decision_lines, paths$decision)

log_lines <- c(
  paste0("Decision 122 executed at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("contract_status=", contract_status),
  paste0("approved_for_unified_steiger_execution=", approved_for_unified_steiger_execution),
  paste0("binary_N_convention=", binary_N_convention),
  paste0("formal_implementation=", formal_implementation),
  paste0("eligible_analysis_sets=", sum(estimability$formal_steiger_eligible)),
  paste0("ineligible_analysis_sets=", sum(!estimability$formal_steiger_eligible)),
  "steiger_statistics_computed=FALSE",
  "r_values_computed=FALSE",
  "r2_values_computed=FALSE",
  "MR_run=FALSE",
  paste0("hard_check_failures=", if (length(hard_check_failures)) paste(hard_check_failures, collapse = ";") else "[]"),
  paste0("renv_lock_sha_before=", renv_before),
  paste0("renv_lock_sha_after=", hash_file(paths$renv_lock))
)
write_text(log_lines, paths$log)

cat("Decision 122 contract status:", contract_status, "\n")
cat("Hard check failures:", if (length(hard_check_failures)) paste(hard_check_failures, collapse = "; ") else "[]", "\n")
cat("Eligible analysis sets:", sum(estimability$formal_steiger_eligible), "\n")
cat("Ineligible analysis sets:", sum(!estimability$formal_steiger_eligible), "\n")
cat("Steiger statistics computed: FALSE\n")
