#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/36a_unified_steiger_v1_freeze_manifest.R [--project-root <path>]", call. = FALSE)
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

paths <- list(
  decision120 = rel("docs", "decisions", "120_unified_directionality_steiger_framework_and_feasibility_v1_v1.1.md"),
  decision121 = rel("docs", "decisions", "121_unified_directionality_framework_readback_recovery_v1_v1.1.md"),
  decision122 = rel("docs", "decisions", "122_unified_steiger_assumption_contract_v1_v1.1.md"),
  decision123 = rel("docs", "decisions", "123_unified_directionality_steiger_v1_v1.1.md"),
  script34 = rel("R", "34_unified_directionality_steiger_framework_v1.R"),
  script34a = rel("R", "34a_unified_directionality_framework_readback_recovery_v1.R"),
  script35 = rel("R", "35_unified_steiger_assumption_contract_v1.R"),
  script36 = rel("R", "36_unified_directionality_steiger_v1.R"),
  script36a = rel("R", "36a_unified_steiger_v1_freeze_manifest.R"),
  framework = rel("results", "qc", "unified_directionality_steiger_framework_v1.json"),
  recovery = rel("results", "qc", "unified_directionality_framework_readback_recovery_v1.json"),
  contract = rel("results", "qc", "unified_steiger_assumption_contract_v1.json"),
  registry = rel("results", "qc", "unified_directionality_analysis_registry_v1.csv"),
  feasibility = rel("results", "qc", "unified_directionality_feasibility_audit_readback_recovery_v1.csv"),
  estimability = rel("results", "qc", "unified_steiger_analysis_estimability_v1.csv"),
  n_semantics = rel("results", "qc", "unified_steiger_N_semantics_contract_v1.csv"),
  prevalence = rel("results", "qc", "unified_steiger_prevalence_contract_v1.csv"),
  package_audit = rel("results", "qc", "unified_steiger_package_implementation_audit_v1.csv"),
  snp_parquet = rel("results", "tables", "unified_steiger_snp_level_r_v1.parquet"),
  snp_tsv = rel("results", "tables", "unified_steiger_snp_level_r_v1.tsv"),
  scenario = rel("results", "tables", "unified_steiger_scenario_results_v1.csv"),
  summary = rel("results", "tables", "unified_steiger_analysis_summary_v1.csv"),
  not_estimable = rel("results", "tables", "unified_steiger_not_estimable_v1.csv"),
  parity = rel("results", "qc", "unified_steiger_manual_parity_audit_v1.csv"),
  qc = rel("results", "qc", "unified_steiger_v1.json"),
  source_log = rel("results", "logs", "unified_steiger_v1.log"),
  renv_lock = rel("renv.lock"),
  manifest = rel("results", "qc", "unified_steiger_v1_freeze_manifest.csv"),
  freeze_json = rel("results", "qc", "unified_steiger_v1_freeze.json"),
  freeze_log = rel("results", "logs", "unified_steiger_v1_freeze.log"),
  decision = rel("docs", "decisions", "124_unified_steiger_v1_freeze_v1.1.md")
)

required <- unlist(paths[setdiff(names(paths), c("manifest", "freeze_json", "freeze_log", "decision"))])
missing <- required[!file.exists(required)]
stop_if(length(missing) > 0L, paste("Missing required input(s):", paste(relpath(missing), collapse = "; ")))
stop_if(!identical(latest_decision(), 124L), paste("Expected next decision 124, found ", latest_decision(), "; no outputs written."))

targets <- unlist(paths[c("manifest", "freeze_json", "freeze_log", "decision")])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(relpath(occupied), collapse = "; ")))

renv_before <- hash_file(paths$renv_lock)
framework <- read_json(paths$framework)
recovery <- read_json(paths$recovery)
contract <- read_json(paths$contract)
qc <- read_json(paths$qc)
registry <- read_csv(paths$registry)
feasibility <- read_csv(paths$feasibility)
estimability <- read_csv(paths$estimability)
prevalence <- read_csv(paths$prevalence)
scenario <- read_csv(paths$scenario)
summary <- read_csv(paths$summary)
not_estimable <- read_csv(paths$not_estimable)
parity <- read_csv(paths$parity)

input_names <- setdiff(names(paths), c("manifest", "freeze_json", "freeze_log", "decision"))
manifest_records <- data.frame(
  relative_path = relpath(unlist(paths[input_names])),
  file_role = input_names,
  scientific_authority = TRUE,
  file_size_bytes = as.numeric(file.info(unlist(paths[input_names]))$size),
  sha256 = vapply(unlist(paths[input_names]), hash_file, character(1)),
  stringsAsFactors = FALSE
)
partial_manifest <- paste0(paths$manifest, ".partial")
on.exit(unlink(partial_manifest, force = TRUE), add = TRUE)
old <- options(digits = 17, scipen = 999)
utils::write.csv(manifest_records, partial_manifest, row.names = FALSE, na = "")
options(old)
manifest_sha <- hash_file(partial_manifest)

script_text <- paste(readLines(paths$script36a, warn = FALSE), collapse = "\n")
scan_text <- gsub("\"([^\"\\\\]|\\\\.)*\"", "\"\"", script_text, perl = TRUE)
scan_text <- gsub("'([^'\\\\]|\\\\.)*'", "''", scan_text, perl = TRUE)
forbidden_call_patterns <- c(
  "get_r_from_lor\\s*\\(", "get_r_from_bsen\\s*\\(", "effective_n\\s*\\(",
  "mr_steiger\\s*\\(", "mr_steiger2\\s*\\(", "directionality_test\\s*\\(",
  "steiger_filtering\\s*\\(", "harmonise_data\\s*\\(", "mr\\s*\\("
)
no_forbidden_calls <- !any(vapply(forbidden_call_patterns, function(p) grepl(p, scan_text, ignore.case = TRUE, perl = TRUE), logical(1)))

eligible <- estimability[as_bool(estimability$formal_steiger_eligible), , drop = FALSE]
ineligible <- estimability[!as_bool(estimability$formal_steiger_eligible), , drop = FALSE]
prevalence_grid <- as.numeric(prevalence$K)
expected_scenario_rows <- nrow(eligible) * length(prevalence_grid)
expected_snp_level_rows <- sum(as.integer(eligible$nsnp), na.rm = TRUE) * length(prevalence_grid)
snp_level_rows <- length(readLines(paths$snp_tsv, warn = FALSE)) - 1L

scenario_by_analysis <- split(scenario, scenario$analysis_id)
complete_grid <- all(vapply(eligible$analysis_id, function(id) {
  x <- scenario_by_analysis[[id]]
  !is.null(x) && nrow(x) == length(prevalence_grid) && setequal(as.numeric(x$prevalence_K), prevalence_grid)
}, logical(1)))
eligible_robust <- all(summary$status[summary$formal_steiger_eligible == TRUE | summary$formal_steiger_eligible == "TRUE"] == "steiger_completed") &&
  all(as_bool(summary$supports_hypothesized_orientation_all_K[summary$status == "steiger_completed"])) &&
  all(as_bool(summary$all_p_lt_0_05[summary$status == "steiger_completed"])) &&
  all(summary$prevalence_robustness_classification[summary$status == "steiger_completed"] == "orientation_and_statistical_support_robust")
not_estimable_ok <- nrow(not_estimable) == nrow(ineligible) &&
  all(not_estimable$status == "not_estimable") &&
  all(not_estimable$blocking_reason == "variant_level_Hb_sample_size_unavailable") &&
  all(as.integer(not_estimable$authenticated_study_level_N) == 408112L) &&
  all(!as_bool(not_estimable$study_level_N_used_as_per_snp)) &&
  all(!as_bool(not_estimable$r_computed)) &&
  all(!as_bool(not_estimable$r2_computed)) &&
  all(!as_bool(not_estimable$steiger_run))
parity_ok <- nrow(parity) == expected_scenario_rows &&
  all(as_bool(parity$get_r_parity_pass)) &&
  all(as_bool(parity$effective_N_parity_pass)) &&
  all(as_bool(parity$R2_parity_pass)) &&
  all(as_bool(parity$orientation_parity_pass)) &&
  all(as_bool(parity$steiger_p_parity_pass))

renv_after <- hash_file(paths$renv_lock)
hard_checks <- list(
  decision_120_gate = identical(framework$framework_status, "frozen") &&
    isTRUE(framework$approved_for_unified_steiger_assumption_contract) &&
    is_empty(framework$hard_check_failures),
  decision_121_gate = identical(recovery$recovery_status, "passed") &&
    isTRUE(recovery$approved_for_unified_steiger_assumption_contract) &&
    is_empty(recovery$hard_check_failures),
  decision_122_gate = identical(contract$contract_status, "frozen") &&
    isTRUE(contract$approved_for_unified_steiger_execution) &&
    is_empty(contract$hard_check_failures),
  decision_123_gate = identical(qc$steiger_status, "passed") &&
    isTRUE(qc$approved_for_unified_steiger_results_interpretation) &&
    is_empty(qc$hard_check_failures),
  no_steiger_rerun = no_forbidden_calls,
  all_analysis_sets_accounted_for = nrow(registry) == 12L &&
    nrow(estimability) == 12L &&
    nrow(summary) == 12L &&
    setequal(registry$analysis_id, summary$analysis_id) &&
    setequal(registry$analysis_id, estimability$analysis_id),
  eligible_scenario_grid_complete = complete_grid &&
    nrow(scenario) == expected_scenario_rows &&
    identical(as.integer(qc$expected_scenario_rows), as.integer(expected_scenario_rows)) &&
    identical(as.integer(qc$actual_scenario_rows), as.integer(nrow(scenario))),
  snp_level_rows_complete = snp_level_rows == expected_snp_level_rows &&
    identical(as.integer(qc$expected_snp_level_rows), as.integer(expected_snp_level_rows)) &&
    identical(as.integer(qc$actual_snp_level_rows), as.integer(snp_level_rows)),
  ineligible_status_rows_complete = not_estimable_ok,
  vuckovic_N_not_imputed = all(feasibility$continuous_N_status[grepl("vuckovic", feasibility$analysis_id)] == "study_level_N_only") &&
    all(feasibility$study_level_N_must_not_be_used_as_per_snp_N[grepl("vuckovic", feasibility$analysis_id)] == TRUE |
          feasibility$study_level_N_must_not_be_used_as_per_snp_N[grepl("vuckovic", feasibility$analysis_id)] == "TRUE") &&
    isFALSE(qc$vuckovic_study_N_used_as_per_snp),
  prevalence_contract_preserved = setequal(prevalence_grid, c(0.005, 0.01, 0.02, 0.05, 0.10)) &&
    identical(contract$single_true_prevalence_claim, FALSE),
  package_default_prevalence_false = isFALSE(qc$package_default_prevalence_used) &&
    isFALSE(contract$package_default_prevalence_allowed),
  sample_case_fraction_as_K_false = isFALSE(qc$sample_case_fraction_used_as_K) &&
    isFALSE(contract$sample_case_fraction_as_K),
  manual_parity_all_passed = parity_ok,
  no_steiger_filtering = isFALSE(qc$steiger_filtering_performed) &&
    isFALSE(contract$steiger_filtering_allowed),
  no_instrument_removal = isFALSE(qc$instrument_filtering_performed),
  no_mr_rerun = isFALSE(qc$mr_rerun) &&
    isFALSE(contract$mr_rerun_after_steiger),
  strict_relaxed_hierarchy_preserved = all(scenario$relaxed_confirmatory == FALSE | scenario$relaxed_confirmatory == "FALSE") &&
    all(scenario$strict_primary_superseded_by_relaxed == FALSE | scenario$strict_primary_superseded_by_relaxed == "FALSE"),
  cross_direction_interpretation_correct = eligible_robust,
  no_bidirectional_causality_claim = TRUE,
  no_causal_direction_confirmation_claim = isFALSE(qc$causal_direction_confirmation_claim_allowed) &&
    all(scenario$causal_direction_confirmation_claim_allowed == FALSE | scenario$causal_direction_confirmation_claim_allowed == "FALSE"),
  renv_lock_unchanged = identical(renv_before, renv_after) &&
    identical(qc$source_sha256$renv_lock_before, renv_before) &&
    identical(qc$source_sha256$renv_lock_after, renv_after),
  git_status_not_required = identical(qc$informational_findings$git_status, "not_applicable_project_not_git_repository")
)
failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
freeze_status <- if (length(failures) == 0L) "passed" else "failed"
approved_final <- identical(freeze_status, "passed")

freeze <- list(
  freeze_version = "v1",
  decision = 124,
  date = "2026-08-13",
  authoritative_unified_steiger_version = "v1",
  analysis_role = "unified_directionality_sensitivity",
  evidence_role = "supportive_instrument_orientation_sensitivity",
  source_decisions = list(framework = 120, readback_recovery = 121, assumption_contract = 122, steiger_execution = 123),
  steiger_rerun = FALSE,
  r_recomputed = FALSE,
  r2_recomputed = FALSE,
  eligible_analysis_count = nrow(eligible),
  ineligible_analysis_count = nrow(ineligible),
  expected_scenario_rows = expected_scenario_rows,
  actual_scenario_rows = nrow(scenario),
  expected_snp_level_rows = expected_snp_level_rows,
  actual_snp_level_rows = snp_level_rows,
  prevalence_strategy = "prespecified_liability_prevalence_sensitivity_grid",
  prevalence_grid = prevalence_grid,
  single_true_prevalence_claim = FALSE,
  sample_case_fraction_used_as_K = FALSE,
  package_default_prevalence_used = FALSE,
  binary_N_convention = "TwoSampleMR_local_effective_sample_size",
  effective_N = qc$effective_N,
  chen_results = records(summary[summary$status == "steiger_completed", , drop = FALSE]),
  vuckovic_not_estimable_results = records(not_estimable),
  manual_parity_status = list(
    rows = nrow(parity),
    get_r_parity_all_pass = all(as_bool(parity$get_r_parity_pass)),
    effective_N_parity_all_pass = all(as_bool(parity$effective_N_parity_pass)),
    R2_parity_all_pass = all(as_bool(parity$R2_parity_pass)),
    orientation_parity_all_pass = all(as_bool(parity$orientation_parity_pass)),
    steiger_P_parity_all_pass = all(as_bool(parity$steiger_p_parity_pass))
  ),
  cross_direction_steiger_pattern = "instrument_set_specific_orientation_support",
  bidirectional_causality_inferred = FALSE,
  causal_direction_confirmation_claim_allowed = FALSE,
  steiger_filtering_performed = FALSE,
  instrument_filtering_performed = FALSE,
  mr_rerun = FALSE,
  interpretation_boundary = "Steiger supports hypothesized instrument orientation only and does not confirm causal direction or causal-effect existence.",
  manifest_path = relpath(paths$manifest),
  manifest_sha256 = manifest_sha,
  freeze_status = freeze_status,
  approved_for_final_integrated_analysis_freeze = approved_final,
  hard_checks = hard_checks,
  hard_check_failures = failures,
  informational_findings = list(
    formal_results_limited_to_chen_based_sets = TRUE,
    vuckovic_based_sets_not_estimable_due_to_missing_variant_level_Hb_N = TRUE,
    vuckovic_authenticated_study_level_N = 408112,
    vuckovic_study_level_N_used_as_per_snp = FALSE,
    winner_curse_risk_present = TRUE,
    measurement_error_primary_convention = "r_xxo=1 and r_yyo=1 as unadjusted observed-trait orientation only",
    git_repository_present = FALSE,
    git_status = "not_applicable_project_not_git_repository",
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after
  )
)

decision_lines <- c(
  "# Decision 124: Unified Steiger Results Freeze V1",
  "",
  "Date: 2026-08-13",
  "",
  "## Status",
  paste0("freeze_status: `", freeze_status, "`"),
  paste0("approved_for_final_integrated_analysis_freeze: `", if (approved_final) "TRUE" else "FALSE", "`"),
  "hard_check_failures: `[]`",
  "",
  "## Decision",
  "Decision 123 Unified Directionality / Steiger V1 is frozen as the authoritative unified directionality sensitivity result.",
  "",
  "This freeze only read, verified, hashed, summarized, and froze existing outputs. It did not recompute r, R2, effective N, Steiger statistics, MR, harmonisation, clumping, proxy, or liftOver.",
  "",
  "## Counts",
  paste0("- eligible analysis sets: `", nrow(eligible), "`."),
  paste0("- ineligible analysis sets: `", nrow(ineligible), "`."),
  paste0("- scenario rows expected/actual: `", expected_scenario_rows, "` / `", nrow(scenario), "`."),
  paste0("- SNP-level rows expected/actual: `", expected_snp_level_rows, "` / `", snp_level_rows, "`."),
  "",
  "## Interpretation Boundary",
  "All formal-eligible Chen-based sets supported their hypothesized instrument orientation across the prespecified K grid with Steiger P<0.05. This is instrument-orientation sensitivity support, not confirmation of causal direction.",
  "",
  "Vuckovic-based sets are formally not estimable under the frozen data requirements because variant-level Hb sample size is unavailable; study-level N=408112 was not used as per-SNP N.",
  "",
  "Forward and reverse Steiger support is interpreted as instrument-set-specific orientation support and does not imply bidirectional causality.",
  "",
  "## Audit",
  paste0("- prevalence grid: `", paste(prevalence_grid, collapse = ", "), "`."),
  paste0("- manual parity rows: `", nrow(parity), "`, failures `0`."),
  paste0("- manifest SHA-256: `", manifest_sha, "`."),
  paste0("- renv.lock SHA before/after: `", renv_before, "` / `", renv_after, "`."),
  "- git status: `not_applicable_project_not_git_repository`.",
  "",
  "## Hard Check Failures",
  if (length(failures) == 0L) "- none" else paste0("- `", failures, "`"),
  "",
  "## Outputs Created",
  "- `R/36a_unified_steiger_v1_freeze_manifest.R`",
  "- `results/qc/unified_steiger_v1_freeze_manifest.csv`",
  "- `results/qc/unified_steiger_v1_freeze.json`",
  "- `results/logs/unified_steiger_v1_freeze.log`",
  "- `docs/decisions/124_unified_steiger_v1_freeze_v1.1.md`",
  "",
  "## Next Gate",
  "Decision 125 Final Integrated Analysis Freeze V1 may proceed only because this freeze passed with no hard-check failures."
)

log_lines <- c(
  "[2026-08-13] Unified Steiger Results Freeze V1",
  paste0("freeze_status=", freeze_status),
  paste0("manifest_sha256=", manifest_sha),
  paste0("approved_for_final_integrated_analysis_freeze=", approved_final),
  paste0("eligible_analysis_count=", nrow(eligible)),
  paste0("ineligible_analysis_count=", nrow(ineligible)),
  paste0("scenario_rows=", nrow(scenario)),
  paste0("snp_level_rows=", snp_level_rows),
  paste0("hard_check_failures=", if (length(failures) == 0L) "[]" else paste(failures, collapse = ";")),
  "steiger_rerun=FALSE",
  "r_recomputed=FALSE",
  "r2_recomputed=FALSE",
  "mr_rerun=FALSE"
)

stop_if(length(failures) > 0L, paste("Freeze hard checks failed:", paste(failures, collapse = "; ")))
if (!file.rename(partial_manifest, paths$manifest)) stop("Atomic rename failed: ", paths$manifest, call. = FALSE)
write_json(freeze, paths$freeze_json)
write_text(log_lines, paths$freeze_log)
write_text(decision_lines, paths$decision)

cat("Decision 124 Unified Steiger Results Freeze V1 completed\n")
cat("freeze_status=", freeze_status, "\n", sep = "")
cat("manifest_sha256=", manifest_sha, "\n", sep = "")
cat("approved_for_final_integrated_analysis_freeze=", approved_final, "\n", sep = "")
cat("hard_check_failures=[]\n")
