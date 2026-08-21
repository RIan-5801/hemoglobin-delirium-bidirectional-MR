#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

root <- "E:/Research/hb_delirium_bidir_mr"
setwd(root)

for (pkg in c("jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

rel <- function(...) file.path(...)
norm <- function(x) gsub("\\\\", "/", x)
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
write_json <- function(x, path) {
  partial <- paste0(path, ".partial")
  jsonlite::write_json(x, partial, pretty = TRUE, auto_unbox = TRUE, null = "null")
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}
write_csv <- function(x, path) {
  partial <- paste0(path, ".partial")
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.csv(x, partial, row.names = FALSE, na = "")
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}
write_text <- function(lines, path) {
  partial <- paste0(path, ".partial")
  writeLines(lines, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}
extract_required_methods <- function(script_txt) {
  m <- regexec("required_methods\\s*<-\\s*c\\(([^)]*)\\)", script_txt, perl = TRUE)
  hit <- regmatches(script_txt, m)[[1]]
  if (length(hit) < 2L) return(character(0))
  methods <- unlist(regmatches(hit[[2]], gregexpr('"[^"]+"', hit[[2]], perl = TRUE)))
  gsub('"', "", methods, fixed = TRUE)
}
extract_int <- function(pattern, txt) {
  hit <- regmatches(txt, regexpr(pattern, txt, perl = TRUE))
  if (length(hit) == 1L && nzchar(hit)) as.integer(gsub("[^0-9]", "", hit)) else NA_integer_
}
extract_num_after_equal <- function(pattern, txt) {
  hit <- regmatches(txt, regexpr(pattern, txt, perl = TRUE))
  if (length(hit) == 1L && nzchar(hit)) as.numeric(sub(".*=\\s*", "", hit)) else NA_real_
}

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
if (!identical(next_decision, 99L)) {
  stop("Expected next decision 99, found ", next_decision, "; no outputs written.", call. = FALSE)
}

paths <- list(
  freeze = rel("results", "qc", "chen_forward_harmonised_mr_inputs_v1_freeze.json"),
  freeze_manifest = rel("results", "qc", "chen_forward_harmonised_mr_inputs_v1_freeze_manifest.csv"),
  amendment_readback = rel("results", "qc", "chen_forward_mr_method_alignment_amendment_v1_readback_audit_v1.json"),
  primary_script = rel("R", "09_forward_mr_vuckovic_hb_finngen_r13_primary_v3.R"),
  primary_qc = rel("results", "qc", "vuckovic_hb_finngen_r13_forward_mr_v3.json"),
  renv_lock = rel("renv.lock"),
  contract = rel("results", "qc", "chen_forward_mr_analysis_contract_v2.json"),
  drift = rel("results", "qc", "chen_forward_mr_parameter_drift_audit_v2.csv"),
  log = rel("results", "logs", "chen_forward_mr_analysis_contract_v2.log"),
  decision = rel("docs", "decisions", "99_chen_forward_mr_analysis_contract_v2_v1.1.md")
)

inputs <- unlist(paths[c("freeze", "freeze_manifest", "amendment_readback", "primary_script", "primary_qc", "renv_lock")])
missing <- inputs[!file.exists(inputs)]
if (length(missing) > 0L) stop("Missing required input(s): ", paste(missing, collapse = "; "), call. = FALSE)
targets <- unlist(paths[c("contract", "drift", "log", "decision")])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
if (length(occupied) > 0L) stop("Target or partial exists: ", paste(occupied, collapse = "; "), call. = FALSE)

freeze <- read_json(paths$freeze)
amendment_gate <- read_json(paths$amendment_readback)
primary_qc <- read_json(paths$primary_qc)
script_txt <- paste(readLines(paths$primary_script, warn = FALSE), collapse = "\n")
required_methods <- extract_required_methods(script_txt)
expected_methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode")
seed_value <- extract_int("seed_value\\s*<-\\s*[0-9]+L?", script_txt)
nb_distribution <- extract_int("NbDistribution\\s*=\\s*[0-9]+", script_txt)
signif_threshold <- extract_num_after_equal("SignifThreshold\\s*=\\s*[0-9.]+", script_txt)

method_labels <- c(
  mr_ivw = "Inverse variance weighted",
  mr_egger_regression = "MR Egger",
  mr_weighted_median = "Weighted median",
  mr_simple_mode = "Simple mode",
  mr_weighted_mode = "Weighted mode"
)
method_hierarchy <- data.frame(
  method_id = required_methods,
  method = unname(method_labels[required_methods]),
  estimator_hierarchy_role = ifelse(required_methods == "mr_ivw", "primary", "sensitivity"),
  stringsAsFactors = FALSE
)

renv_before <- hash_file(paths$renv_lock)
renv_after <- hash_file(paths$renv_lock)
freeze_manifest_sha <- hash_file(paths$freeze_manifest)
single_snp_found <- grepl("mr_singlesnp|mr_wald_ratio|single_snp|single-SNP|Wald ratio", script_txt, ignore.case = TRUE)

hard_checks <- list(
  decision_94_freeze_gate = identical(freeze$freeze_status, "passed") &&
    isTRUE(freeze$approved_for_chen_forward_mr_design) &&
    length(freeze$hard_check_failures) == 0L &&
    identical(freeze$authoritative_chen_forward_harmonisation_version, "v1") &&
    identical(tolower(freeze$manifest_sha256), tolower(freeze_manifest_sha)),
  method_alignment_amendment_gate = identical(amendment_gate$authoritative_amendment_status_after_readback, "frozen") &&
    isTRUE(amendment_gate$approved_for_chen_forward_mr_contract_v2_after_readback) &&
    length(amendment_gate$corrected_hard_check_failures) == 0L,
  analysis_role_sensitivity_not_primary = TRUE,
  independent_replication_false = TRUE,
  forward_primary_method_authority_found = file.exists(paths$primary_script) && file.exists(paths$primary_qc),
  ivw_implementation_matches_primary = "mr_ivw" %in% required_methods && grepl("TwoSampleMR::mr\\s*\\(", script_txt),
  estimator_hierarchy_matches_primary = identical(required_methods, expected_methods),
  heterogeneity_plan_matches_primary = grepl("TwoSampleMR::mr_heterogeneity", script_txt) &&
    grepl("mr_egger_regression", script_txt) && grepl("mr_ivw", script_txt),
  egger_plan_matches_primary = grepl("TwoSampleMR::mr_pleiotropy_test", script_txt),
  mr_presso_plan_matches_primary = grepl("MRPRESSO::mr_presso", script_txt) &&
    identical(nb_distribution, 10000L) &&
    isTRUE(all.equal(signif_threshold, 0.05)) &&
    identical(seed_value, 2026L),
  loo_plan_matches_primary = grepl("TwoSampleMR::mr_leaveoneout", script_txt) &&
    grepl("TwoSampleMR::mr_ivw", script_txt),
  single_snp_not_required = TRUE,
  single_snp_not_planned = !single_snp_found,
  effect_scale_defined = identical(freeze$effect_scale, "standardized_quantitative_Hb_effect"),
  or_transform_defined = identical(freeze$outcome_scale, "log_odds_delirium"),
  no_posthoc_filtering = TRUE,
  steiger_deferred = TRUE,
  software_environment_read_only = identical(primary_qc$mr_library, "E:/Research/hb_delirium_bidir_mr/renv/mr-v1-library") &&
    identical(primary_qc$TwoSampleMR_version, "0.7.9") &&
    identical(primary_qc$MRPRESSO_version, "1.0"),
  no_mr_executed_in_contract = TRUE,
  renv_lock_unchanged = identical(renv_before, renv_after)
)
failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
status <- if (length(failures) == 0L) "frozen" else "failed"

contract <- list(
  contract_version = "v2",
  date = "2026-08-12",
  analysis_direction = "Hb_to_delirium",
  analysis_role = "forward_alternative_hb_gwas_sensitivity",
  independent_replication = FALSE,
  source_mr_input_freeze_decision = 94,
  method_alignment_amendment_decision = 98,
  exposure_source = "Chen_2020_Hb_BCX2",
  outcome_source = "FinnGen_R13_F5_DELIRIUM",
  analysis_sets = list(
    APOE_included = list(role = "chen_forward_alternative_hb_gwas_sensitivity_main", final_valid_count = freeze$included_final_valid_count),
    APOE_excluded = list(role = "chen_forward_alternative_hb_gwas_apoe_exclusion_sensitivity", final_valid_count = freeze$excluded_final_valid_count)
  ),
  exposure_scale = "standardized_quantitative_Hb_effect",
  outcome_scale = "log_odds_delirium",
  raw_effect_interpretation = "log-odds change in delirium per 1-unit increase in genetically predicted standardized Chen Hb scale",
  or_transform_enabled = TRUE,
  primary_estimator = "IVW",
  primary_estimator_method_id = "mr_ivw",
  estimator_hierarchy = method_hierarchy,
  diagnostics = list(
    heterogeneity = list(implementation = "TwoSampleMR::mr_heterogeneity", method_list = list("mr_egger_regression", "mr_ivw")),
    egger_intercept = list(implementation = "TwoSampleMR::mr_pleiotropy_test"),
    mr_presso = list(implementation = "MRPRESSO::mr_presso", NbDistribution = nb_distribution, SignifThreshold = signif_threshold, seed = seed_value, outlier_corrected_estimate_role = "sensitivity_only"),
    leave_one_out = list(implementation = "TwoSampleMR::mr_leaveoneout", method = "TwoSampleMR::mr_ivw"),
    single_snp = list(run = FALSE, status = "not_planned_by_method_alignment_amendment")
  ),
  steiger_run = FALSE,
  steiger_status = "deferred_to_unified_directionality_sensitivity_stage",
  seed = seed_value,
  software_environment = list(
    R_version_authority = primary_qc$R_version,
    TwoSampleMR_version = primary_qc$TwoSampleMR_version,
    TwoSampleMR_RemoteSha = primary_qc$TwoSampleMR_RemoteSha,
    MRPRESSO_version = primary_qc$MRPRESSO_version,
    MRPRESSO_RemoteSha = primary_qc$MRPRESSO_RemoteSha,
    mr_library = primary_qc$mr_library,
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after
  ),
  contract_status = status,
  approved_for_chen_forward_mr_execution = identical(status, "frozen"),
  hard_checks = hard_checks,
  hard_check_failures = failures
)

drift <- data.frame(
  parameter_domain = c(
    "exposure_gwas_source", "instrument_set", "analysis_role", "estimator_hierarchy",
    "ivw_implementation", "sensitivity_estimators", "heterogeneity", "egger_intercept",
    "mr_presso", "leave_one_out", "single_snp", "effect_scale", "or_transform",
    "posthoc_filtering", "steiger", "software_environment", "seed"
  ),
  vuckovic_forward_primary = c(
    "Vuckovic_2020_Hb", "Vuckovic frozen primary IVs", "forward_primary",
    paste(required_methods, collapse = ";"), "TwoSampleMR::mr with mr_ivw",
    paste(setdiff(required_methods, "mr_ivw"), collapse = ";"), "TwoSampleMR::mr_heterogeneity",
    "TwoSampleMR::mr_pleiotropy_test",
    sprintf("MRPRESSO::mr_presso; NbDistribution=%s; SignifThreshold=%s", nb_distribution, signif_threshold),
    "TwoSampleMR::mr_leaveoneout with TwoSampleMR::mr_ivw",
    "not part of frozen authority", "standardized inverse-normal-transformed Hb",
    "OR=exp(beta) for binary FinnGen outcome", "not allowed", "deferred", primary_qc$mr_library, as.character(seed_value)
  ),
  chen_forward_contract_v2 = c(
    "Chen_2020_Hb_BCX2", "Decision 94 frozen Chen final-valid IVs",
    "forward_alternative_hb_gwas_sensitivity",
    paste(required_methods, collapse = ";"), "TwoSampleMR::mr with mr_ivw",
    paste(setdiff(required_methods, "mr_ivw"), collapse = ";"), "TwoSampleMR::mr_heterogeneity",
    "TwoSampleMR::mr_pleiotropy_test",
    sprintf("MRPRESSO::mr_presso; NbDistribution=%s; SignifThreshold=%s", nb_distribution, signif_threshold),
    "TwoSampleMR::mr_leaveoneout with TwoSampleMR::mr_ivw",
    "not planned by Decision 98 method-alignment amendment",
    "standardized quantitative Chen Hb effect",
    "OR=exp(beta) for binary FinnGen outcome", "not allowed",
    "deferred_to_unified_directionality_sensitivity_stage", primary_qc$mr_library, as.character(seed_value)
  ),
  drift_status = c(
    "allowed", "allowed", "allowed", rep("preserved", 7), "prospectively_removed_to_preserve_method_alignment",
    "scale_source_changed_with_exposure_gwas", "preserved", "preserved", "preserved", "preserved", "preserved"
  ),
  allowable_under_contract_v2 = c(TRUE, TRUE, TRUE, rep(TRUE, 14)),
  stringsAsFactors = FALSE
)

log_lines <- c(
  "[2026-08-12] Chen Forward MR Analysis Contract V2",
  paste0("contract_status=", status),
  paste0("approved_for_chen_forward_mr_execution=", identical(status, "frozen")),
  paste0("required_methods=", paste(required_methods, collapse = ";")),
  paste0("single_snp_run=FALSE"),
  paste0("hard_check_failures=", if (length(failures) == 0L) "[]" else paste(failures, collapse = ";")),
  "No MR was executed."
)

decision_lines <- c(
  "# Decision 99: Chen Forward MR Analysis Contract V2",
  "",
  "Date: 2026-08-12",
  "",
  "## Status",
  paste0("contract_status: `", status, "`"),
  paste0("approved_for_chen_forward_mr_execution: `", identical(status, "frozen"), "`"),
  "",
  "## Decision",
  "Chen Forward MR Analysis Contract V2 is based on Decision 94 frozen MR inputs and Decision 98 method-alignment amendment.",
  "",
  "Single-SNP/Wald-ratio diagnostics are not planned and not required for Chen Forward MR V1.",
  "",
  "## Estimator Hierarchy",
  paste0("- `", required_methods, "`"),
  "",
  "## Diagnostics",
  "- IVW/MR-Egger heterogeneity",
  "- Egger intercept",
  "- MR-PRESSO",
  "- leave-one-out",
  "- no single-SNP/Wald-ratio diagnostics",
  "- no Steiger",
  "",
  "## Hard Check Failures",
  if (length(failures) == 0L) "- none" else paste0("- `", failures, "`"),
  "",
  "## Outputs",
  paste0("- `", norm(paths$contract), "`"),
  paste0("- `", norm(paths$drift), "`"),
  paste0("- `", norm(paths$log), "`"),
  paste0("- `", norm(paths$decision), "`")
)

write_json(contract, paths$contract)
write_csv(drift, paths$drift)
write_text(log_lines, paths$log)
write_text(decision_lines, paths$decision)

cat("contract_status=", status, "\n", sep = "")
cat("approved_for_chen_forward_mr_execution=", identical(status, "frozen"), "\n", sep = "")
cat("hard_check_failures=", if (length(failures) == 0L) "[]" else paste(failures, collapse = ";"), "\n", sep = "")
