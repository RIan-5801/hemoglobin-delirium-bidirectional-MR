#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/19_reverse_strict_primary_mr_v1.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))
mr_library <- file.path(root, "renv", "mr-v1-library")
.libPaths(c(mr_library, .libPaths()))

for (pkg in c("DBI", "duckdb", "jsonlite", "digest", "TwoSampleMR")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
hash_text <- function(x) digest::digest(x, algo = "sha256")
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
sql_string <- function(path, must_work = TRUE) {
  paste0("'", gsub("'", "''", normalizePath(path, winslash = "/", mustWork = must_work), fixed = TRUE), "'")
}
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
num_equal <- function(a, b, tol = 1e-12) {
  a <- as.numeric(a); b <- as.numeric(b)
  is.finite(a) && is.finite(b) && abs(a - b) <= tol
}
write_csv_precise <- function(x, path) {
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  write.csv(x, path, row.names = FALSE, na = "")
}

contract_path <- file.path(root, "results", "qc", "reverse_strict_primary_mr_analysis_contract_v1.json")
freeze_path <- file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_formal_harmonisation_v4_freeze.json")
freeze_manifest_path <- file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_formal_harmonisation_v4_freeze_manifest.csv")
renv_lock <- file.path(root, "renv.lock")
included_path <- file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_primary_harmonised_apoe_included_v4.parquet")
excluded_path <- file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_primary_harmonised_apoe_excluded_v4.parquet")

out <- c(
  estimates = file.path(root, "results", "tables", "reverse_strict_primary_mr_estimates_v1.csv"),
  doubling = file.path(root, "results", "tables", "reverse_strict_primary_mr_estimates_doubling_odds_v1.csv"),
  qc = file.path(root, "results", "qc", "reverse_strict_primary_mr_v1.json"),
  log = file.path(root, "results", "logs", "reverse_strict_primary_mr_v1.log")
)
stop_if(any(file.exists(c(out, paste0(out, ".partial")))), "Strict primary MR V1 target or partial exists; refusing to overwrite.")
dir.create(dirname(out[["estimates"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["qc"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["log"]]), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = out[["log"]], append = TRUE)

read_one <- function(con, path, analysis_set, analysis_role, expected_rsid) {
  x <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(path)))
  stop_if(nrow(x) != 1L, paste("Expected one row for", analysis_set))
  required <- c("rsid", "exposure_beta", "exposure_se", "outcome_beta_harmonised",
                "outcome_se_harmonised", "final_valid_instrument")
  stop_if(!all(required %in% names(x)), paste("Required columns missing for", analysis_set))
  stop_if(!identical(as.character(x$rsid[[1L]]), expected_rsid), paste("Unexpected rsID for", analysis_set))
  stop_if(!isTRUE(x$final_valid_instrument[[1L]]), paste("Final-valid flag is not TRUE for", analysis_set))
  x$analysis_set_contract <- analysis_set
  x$analysis_role_contract <- analysis_role
  x
}

wald_one <- function(x, analysis_set, analysis_role) {
  b_exp <- as.numeric(x$exposure_beta[[1L]])
  se_exp <- as.numeric(x$exposure_se[[1L]])
  b_out <- as.numeric(x$outcome_beta_harmonised[[1L]])
  se_out <- as.numeric(x$outcome_se_harmonised[[1L]])
  stop_if(!all(is.finite(c(b_exp, se_exp, b_out, se_out))) || se_exp <= 0 || se_out <= 0 || b_exp == 0,
          paste("Invalid Wald inputs for", analysis_set))
  wr <- TwoSampleMR::mr_wald_ratio(b_exp = b_exp, b_out = b_out, se_exp = se_exp, se_out = se_out)
  beta <- as.numeric(wr$b)
  se <- as.numeric(wr$se)
  pval <- as.numeric(wr$pval)
  ci_lower <- beta - 1.96 * se
  ci_upper <- beta + 1.96 * se
  beta_delta <- b_out / b_exp
  se_delta <- sqrt(se_out^2 / b_exp^2 + b_out^2 * se_exp^2 / b_exp^4)
  p_delta <- stats::pnorm(abs(beta_delta) / se_delta, lower.tail = FALSE) * 2
  ci_delta_lower <- beta_delta - 1.96 * se_delta
  ci_delta_upper <- beta_delta + 1.96 * se_delta
  data.frame(
    analysis_set = analysis_set,
    analysis_role = analysis_role,
    rsid = as.character(x$rsid[[1L]]),
    nsnp = 1L,
    method = "Wald ratio",
    method_id = "mr_wald_ratio",
    authoritative_implementation = "TwoSampleMR::mr_wald_ratio",
    exposure_beta = b_exp,
    exposure_se = se_exp,
    outcome_beta = b_out,
    outcome_se = se_out,
    beta = beta,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    pval = pval,
    effect_scale = "standardized_Hb_per_1_unit_genetically_predicted_log_odds_delirium",
    manual_delta_beta = beta_delta,
    manual_delta_se = se_delta,
    manual_delta_ci_lower = ci_delta_lower,
    manual_delta_ci_upper = ci_delta_upper,
    manual_delta_pval = p_delta,
    manual_vs_software_beta_abs_diff = abs(beta_delta - beta),
    manual_vs_software_se_abs_diff = abs(se_delta - se),
    manual_vs_software_p_abs_diff = abs(p_delta - pval),
    formula_difference_expected = TRUE,
    stringsAsFactors = FALSE
  )
}

main <- function() {
  log_line("stage=reverse_strict_primary_mr_v1")
  set.seed(2026L)
  renv_before <- hash_file(renv_lock)
  contract <- jsonlite::fromJSON(contract_path, simplifyVector = FALSE)
  freeze <- jsonlite::fromJSON(freeze_path, simplifyVector = FALSE)
  manifest_sha <- hash_file(freeze_manifest_path)
  ts_desc <- utils::packageDescription("TwoSampleMR")
  wald_fun <- get("mr_wald_ratio", envir = asNamespace("TwoSampleMR"))
  wald_body_sha <- hash_text(paste(deparse(body(wald_fun), width.cutoff = 500L), collapse = "\n"))

  pre_checks <- list(
    contract_frozen = identical(contract$contract_status, "frozen"),
    contract_approved = isTRUE(contract$approved_for_reverse_strict_primary_mr_execution),
    strict_harmonisation_freeze_gate = identical(freeze$freeze_status, "passed") &&
      identical(freeze$harmonisation_status, "passed") && length(freeze$hard_check_failures) == 0L,
    freeze_manifest_sha_verified = identical(tolower(manifest_sha), tolower(freeze$manifest_sha256)),
    strict_primary_threshold_exact = num_equal(contract$instrument_threshold, 5e-8),
    wald_method_confirmed = exists("mr_wald_ratio", envir = asNamespace("TwoSampleMR"), inherits = FALSE),
    wald_body_sha_matches_contract = identical(wald_body_sha, contract$wald_implementation$function_body_sha256),
    twosamplemr_version_matches_contract = identical(as.character(utils::packageVersion("TwoSampleMR")), contract$wald_implementation$version),
    twosamplemr_sha_matches_contract = identical(as.character(ts_desc[["RemoteSha"]]), contract$wald_implementation$RemoteSha)
  )
  stop_if(!all(unlist(pre_checks)), paste("Pre-execution gate failed:", paste(names(pre_checks)[!unlist(pre_checks)], collapse = "; ")))

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  inc <- read_one(con, included_path, "APOE_included", "strict_reverse_primary_analysis", "rs429358")
  exc <- read_one(con, excluded_path, "APOE_excluded", "strict_reverse_apoe_region_exclusion_sensitivity", "rs58537897")
  results <- rbind(
    wald_one(inc, "APOE_included", "strict_reverse_primary_analysis"),
    wald_one(exc, "APOE_excluded", "strict_reverse_apoe_region_exclusion_sensitivity")
  )

  doubling <- results
  doubling$beta_raw <- doubling$beta
  doubling$se_raw <- doubling$se
  doubling$ci_lower_raw <- doubling$ci_lower
  doubling$ci_upper_raw <- doubling$ci_upper
  doubling$beta_per_doubling_odds <- doubling$beta * log(2)
  doubling$se_per_doubling_odds <- doubling$se * log(2)
  doubling$ci_lower_per_doubling <- doubling$ci_lower * log(2)
  doubling$ci_upper_per_doubling <- doubling$ci_upper * log(2)
  doubling$effect_scale_secondary <- "SD_Hb_per_doubling_of_genetically_predicted_odds_of_delirium"
  doubling <- doubling[, c("analysis_set", "analysis_role", "rsid", "nsnp", "method", "beta_raw", "se_raw",
                           "ci_lower_raw", "ci_upper_raw", "pval", "beta_per_doubling_odds",
                           "se_per_doubling_odds", "ci_lower_per_doubling", "ci_upper_per_doubling",
                           "effect_scale_secondary")]

  hard_checks <- c(pre_checks, list(
    included_nsnp_equals_one = nrow(inc) == 1L && results$nsnp[results$analysis_set == "APOE_included"] == 1L,
    excluded_nsnp_equals_one = nrow(exc) == 1L && results$nsnp[results$analysis_set == "APOE_excluded"] == 1L,
    no_ivw_single_snp_mislabel = !any(results$method == "IVW" | grepl("Inverse variance", results$method, ignore.case = TRUE)),
    effect_scale_correct = all(results$effect_scale == "standardized_Hb_per_1_unit_genetically_predicted_log_odds_delirium"),
    no_hb_physical_units = !any(grepl("g/dL|g/L", results$effect_scale, ignore.case = TRUE)),
    package_implementation_audited = TRUE,
    manual_formula_audited = all(is.finite(results$manual_delta_se)),
    included_estimate_completed = all(is.finite(unlist(results[results$analysis_set == "APOE_included", c("beta", "se", "pval")]))),
    excluded_estimate_completed = all(is.finite(unlist(results[results$analysis_set == "APOE_excluded", c("beta", "se", "pval")]))),
    doubling_odds_rescaling_correct = all(abs(doubling$beta_per_doubling_odds - results$beta * log(2)) <= 1e-15) &&
      all(abs(doubling$se_per_doubling_odds - results$se * log(2)) <= 1e-15) &&
      identical(doubling$pval, results$pval),
    diagnostics_not_estimable = TRUE,
    no_heterogeneity_run = TRUE,
    no_egger_run = TRUE,
    no_presso_run = TRUE,
    no_loo_run = TRUE,
    no_steiger_run = TRUE,
    no_relaxed_mr_run = TRUE,
    no_reharmonisation = TRUE,
    no_proxy = TRUE,
    no_liftover = TRUE,
    no_rescan_or_clumping = TRUE
  ))
  renv_after <- hash_file(renv_lock)
  hard_checks <- c(hard_checks, list(renv_lock_unchanged = identical(renv_before, renv_after)))
  hard_check_failures <- names(hard_checks)[!unlist(hard_checks)]
  status <- if (length(hard_check_failures) == 0L) "passed" else "failed"

  write_csv_precise(results, paste0(out[["estimates"]], ".partial"))
  write_csv_precise(doubling, paste0(out[["doubling"]], ".partial"))
  qc <- list(
    mr_version = "v1",
    decision = 74,
    source_contract_decision = 73,
    analysis_direction = "delirium_to_Hb",
    analysis_role = "secondary_reverse_primary",
    instrument_threshold = 5e-8,
    included_nsnp = 1L,
    excluded_nsnp = 1L,
    included_rsid = results$rsid[results$analysis_set == "APOE_included"],
    excluded_rsid = results$rsid[results$analysis_set == "APOE_excluded"],
    authoritative_method = "Wald ratio",
    wald_implementation = list(
      authoritative_implementation = "TwoSampleMR::mr_wald_ratio",
      TwoSampleMR_version = as.character(utils::packageVersion("TwoSampleMR")),
      TwoSampleMR_RemoteSha = as.character(ts_desc[["RemoteSha"]]),
      function_body_sha256 = wald_body_sha,
      software_se_formula = "se_outcome / abs(beta_exposure)",
      manual_delta_se_formula = "sqrt(se_outcome^2 / beta_exposure^2 + beta_outcome^2 * se_exposure^2 / beta_exposure^4)",
      manual_delta_audit_role = "implementation_difference_audit_not_authoritative_replacement"
    ),
    raw_results = records(results),
    doubling_odds_results = records(doubling),
    diagnostics_not_estimable = list(
      heterogeneity = "not_estimable_single_instrument",
      MR_Egger = "not_estimable_single_instrument",
      Egger_intercept = "not_estimable_single_instrument",
      weighted_median = "not_applicable_single_instrument",
      mode = "not_applicable_single_instrument",
      MR_PRESSO = "not_applicable_single_instrument",
      leave_one_out = "not_run_single_instrument",
      Steiger = "deferred_to_separate_directionality_sensitivity_contract"
    ),
    input_sha256 = list(
      included = hash_file(included_path),
      excluded = hash_file(excluded_path),
      freeze_manifest = manifest_sha
    ),
    R_version = R.version.string,
    seed = 2026L,
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after,
    renv_lock_unchanged = identical(renv_before, renv_after),
    strict_primary_mr_status = status,
    approved_for_reverse_strict_primary_results_freeze = identical(status, "passed"),
    hard_checks = hard_checks,
    hard_check_failures = hard_check_failures
  )
  jsonlite::write_json(qc, paste0(out[["qc"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  if (!identical(status, "passed")) stop("Strict primary MR V1 did not pass; partial outputs retained.", call. = FALSE)
  for (path in out[c("estimates", "doubling", "qc")]) {
    stop_if(file.exists(path), paste("Output appeared during run:", path))
    stop_if(!file.rename(paste0(path, ".partial"), path), paste("Atomic rename failed:", path))
  }
  log_line("strict_primary_mr_status=passed approved_for_reverse_strict_primary_results_freeze=TRUE hard_check_failures=[]")
}

tryCatch(main(), error = function(e) {
  log_line("status=failed error=", conditionMessage(e))
  quit(status = 1L)
})
