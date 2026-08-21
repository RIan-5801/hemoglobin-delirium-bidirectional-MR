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
write_text <- function(lines, path) {
  partial <- paste0(path, ".partial")
  writeLines(lines, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
if (!identical(next_decision, 97L)) {
  stop("Expected next decision 97, found ", next_decision, "; no outputs written.", call. = FALSE)
}

paths <- list(
  decision95 = rel("docs", "decisions", "95_chen_forward_mr_analysis_contract_v1_v1.1.md"),
  decision96 = rel("docs", "decisions", "96_chen_forward_mr_analysis_contract_v1_readback_audit_v1.1.md"),
  contract_v1 = rel("results", "qc", "chen_forward_mr_analysis_contract_v1.json"),
  readback = rel("results", "qc", "chen_forward_mr_analysis_contract_v1_readback_audit_v1.json"),
  primary_script = rel("R", "09_forward_mr_vuckovic_hb_finngen_r13_primary_v3.R"),
  primary_qc = rel("results", "qc", "vuckovic_hb_finngen_r13_forward_mr_v3.json"),
  renv_lock = rel("renv.lock"),
  amendment_json = rel("results", "qc", "chen_forward_mr_method_alignment_amendment_v1.json"),
  log = rel("results", "logs", "chen_forward_mr_method_alignment_amendment_v1.log"),
  decision = rel("docs", "decisions", "97_chen_forward_mr_method_alignment_amendment_v1_v1.1.md")
)

inputs <- unlist(paths[c("decision95", "decision96", "contract_v1", "readback", "primary_script", "primary_qc", "renv_lock")])
missing <- inputs[!file.exists(inputs)]
if (length(missing) > 0L) stop("Missing required input(s): ", paste(missing, collapse = "; "), call. = FALSE)

targets <- unlist(paths[c("amendment_json", "log", "decision")])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
if (length(occupied) > 0L) stop("Target or partial exists: ", paste(occupied, collapse = "; "), call. = FALSE)

contract_v1 <- read_json(paths$contract_v1)
readback <- read_json(paths$readback)
primary_qc <- read_json(paths$primary_qc)
primary_txt <- paste(readLines(paths$primary_script, warn = FALSE), collapse = "\n")

single_snp_absent <- !grepl("mr_singlesnp|mr_wald_ratio|single_snp|single-SNP|Wald ratio", primary_txt, ignore.case = TRUE)
core_methods <- readback$recovered_required_methods
expected_methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode")
renv_before <- hash_file(paths$renv_lock)
renv_after <- hash_file(paths$renv_lock)

hard_checks <- list(
  decision_95_preserved = file.exists(paths$decision95) && file.exists(paths$contract_v1),
  decision_96_preserved = file.exists(paths$decision96) && file.exists(paths$readback),
  no_mr_executed_before_amendment = isTRUE(readback$no_mr_executed) && isFALSE(readback$approved_for_chen_forward_mr_execution_after_readback),
  no_results_available_before_amendment = !file.exists(rel("results", "qc", "chen_forward_mr_v1.json")) &&
    !file.exists(rel("results", "tables", "chen_forward_mr_estimates_v1.csv")),
  single_snp_absent_from_primary_authority = single_snp_absent && isFALSE(readback$corrected_hard_checks$single_snp_plan_matches_primary),
  single_snp_removed_from_chen_plan = TRUE,
  core_estimators_unchanged = identical(core_methods, expected_methods),
  heterogeneity_plan_unchanged = isTRUE(readback$corrected_hard_checks$heterogeneity_plan_matches_primary),
  egger_plan_unchanged = isTRUE(readback$corrected_hard_checks$egger_plan_matches_primary),
  mr_presso_plan_unchanged = isTRUE(readback$corrected_hard_checks$mr_presso_configuration_matches_primary),
  loo_plan_unchanged = isTRUE(readback$corrected_hard_checks$loo_plan_matches_primary),
  seed_unchanged = identical(as.integer(contract_v1$seed), 2026L) && identical(as.integer(primary_qc$seed), 2026L),
  no_steiger = TRUE,
  result_agnostic_amendment = isFALSE(file.exists(rel("results", "qc", "chen_forward_mr_v1.json"))) &&
    isFALSE(file.exists(rel("results", "tables", "chen_forward_mr_estimates_v1.csv"))),
  renv_lock_unchanged = identical(renv_before, renv_after)
)
failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
status <- if (length(failures) == 0L) "frozen" else "failed"

amendment <- list(
  amendment_version = "v1",
  date = "2026-08-12",
  amendment_type = "prospective_result_agnostic_method_alignment_amendment",
  trigger_contract_decision = 95,
  supporting_readback_audit_decision = 96,
  mr_executed_before_amendment = FALSE,
  results_available_before_amendment = FALSE,
  forward_primary_method_authority = list(
    script = norm(paths$primary_script),
    qc = norm(paths$primary_qc),
    decisions = list(31, 35, 37),
    recovered_methods = core_methods,
    heterogeneity = "TwoSampleMR::mr_heterogeneity",
    egger_intercept = "TwoSampleMR::mr_pleiotropy_test",
    mr_presso = list(
      implementation = "MRPRESSO::mr_presso",
      NbDistribution = contract_v1$mr_presso_plan$NbDistribution,
      SignifThreshold = contract_v1$mr_presso_plan$SignifThreshold,
      seed = contract_v1$mr_presso_plan$seed
    ),
    leave_one_out = "TwoSampleMR::mr_leaveoneout with TwoSampleMR::mr_ivw",
    single_snp_or_wald_ratio_found = !single_snp_absent
  ),
  single_snp_required_before = TRUE,
  single_snp_required_after = FALSE,
  single_snp_allowed_after = FALSE,
  single_snp_status = "not_planned_to_preserve_method_alignment_with_forward_primary",
  single_snp_reason = "single-SNP/Wald-ratio diagnostics were not part of the frozen Vuckovic Forward Primary MR V3 authority",
  core_estimator_hierarchy_unchanged = TRUE,
  diagnostic_framework_unchanged_except_single_snp_removal = TRUE,
  steiger_status = "deferred_to_unified_directionality_sensitivity_stage",
  scope = list(
    analysis_direction = "Hb_to_delirium",
    analysis_role = "forward_alternative_hb_gwas_sensitivity",
    target_analysis = "Chen Forward MR V1",
    does_not_modify = list("Vuckovic Forward Primary", "Reverse analyses", "future projects")
  ),
  preserved_evidence = list(
    decision95_sha256 = hash_file(paths$decision95),
    decision96_sha256 = hash_file(paths$decision96),
    contract_v1_sha256 = hash_file(paths$contract_v1),
    readback_audit_sha256 = hash_file(paths$readback)
  ),
  amendment_status = status,
  approved_for_chen_forward_mr_contract_v2 = identical(status, "frozen"),
  hard_checks = hard_checks,
  hard_check_failures = failures,
  renv_lock_sha_before = renv_before,
  renv_lock_sha_after = renv_after
)

log_lines <- c(
  "[2026-08-12] Chen Forward MR Method-Alignment Amendment V1",
  paste0("amendment_status=", status),
  paste0("approved_for_chen_forward_mr_contract_v2=", identical(status, "frozen")),
  paste0("single_snp_required_after=", FALSE),
  paste0("single_snp_allowed_after=", FALSE),
  paste0("hard_check_failures=", if (length(failures) == 0L) "[]" else paste(failures, collapse = ";")),
  "No MR was executed."
)

decision_lines <- c(
  "# Decision 97: Chen Forward MR Method-Alignment Amendment V1",
  "",
  "Date: 2026-08-12",
  "",
  "## Status",
  paste0("amendment_status: `", status, "`"),
  paste0("approved_for_chen_forward_mr_contract_v2: `", identical(status, "frozen"), "`"),
  "",
  "## Decision",
  "This is a prospective, result-agnostic method-alignment amendment for Chen Forward MR V1.",
  "",
  "The single-SNP/Wald-ratio requirement is removed from Chen Forward MR V1 because it was not part of the frozen Vuckovic Forward Primary MR V3 authority. Chen forward MR had not been executed and no Chen MR results were available before this amendment.",
  "",
  "## Frozen Rule",
  "- `single_snp_analysis_required=FALSE`",
  "- `single_snp_analysis_allowed=FALSE` for Chen Forward MR V1",
  "- `single_snp_status=not_planned_to_preserve_method_alignment_with_forward_primary`",
  "",
  "## Unchanged Methods",
  "- IVW",
  "- MR-Egger",
  "- weighted median",
  "- weighted mode",
  "- simple mode",
  "- IVW/MR-Egger heterogeneity",
  "- Egger intercept",
  "- MR-PRESSO",
  "- leave-one-out",
  "- OR transformation",
  "- Vuckovic interpretive comparison",
  "- seed 2026",
  "",
  "## Hard Check Failures",
  if (length(failures) == 0L) "- none" else paste0("- `", failures, "`"),
  "",
  "## Outputs",
  paste0("- `", norm(paths$amendment_json), "`"),
  paste0("- `", norm(paths$log), "`"),
  paste0("- `", norm(paths$decision), "`")
)

write_json(amendment, paths$amendment_json)
write_text(log_lines, paths$log)
write_text(decision_lines, paths$decision)

cat("amendment_status=", status, "\n", sep = "")
cat("approved_for_chen_forward_mr_contract_v2=", identical(status, "frozen"), "\n", sep = "")
cat("hard_check_failures=", if (length(failures) == 0L) "[]" else paste(failures, collapse = ";"), "\n", sep = "")
