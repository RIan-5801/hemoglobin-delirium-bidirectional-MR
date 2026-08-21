#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/20_reverse_relaxed_mr_v1_freeze_manifest.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

for (pkg in c("digest", "jsonlite")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
num_equal <- function(a, b, tol = 1e-8) {
  a <- as.numeric(a); b <- as.numeric(b)
  length(a) == length(b) && all((is.na(a) & is.na(b)) | (is.finite(a) & is.finite(b) & abs(a - b) <= tol))
}
bool_col <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% "true"
}
as_df <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.list(x) && length(x) > 0L) {
    fields <- unique(unlist(lapply(x, names), use.names = FALSE))
    rows <- lapply(x, function(rec) {
      vals <- lapply(fields, function(k) {
        v <- rec[[k]]
        if (is.null(v) || length(v) == 0L) return(NA)
        if (is.list(v)) return(as.character(jsonlite::toJSON(v, auto_unbox = TRUE, na = "null")))
        v[[1L]]
      })
      names(vals) <- fields
      as.data.frame(vals, stringsAsFactors = FALSE, check.names = FALSE)
    })
    return(do.call(rbind, rows))
  }
  data.frame()
}
read_table <- function(rel) read.csv(file.path(root, rel), stringsAsFactors = FALSE, check.names = FALSE)

out <- c(
  manifest = file.path(root, "results", "qc", "reverse_relaxed_mr_v1_freeze_manifest.csv"),
  json = file.path(root, "results", "qc", "reverse_relaxed_mr_v1_freeze.json"),
  log = file.path(root, "results", "logs", "reverse_relaxed_mr_v1_freeze.log")
)
decision_rel <- "docs/decisions/77_reverse_relaxed_mr_v1_freeze_v1.1.md"
script_rel <- "R/20_reverse_relaxed_mr_v1_freeze_manifest.R"

stop_if(any(file.exists(c(out, paste0(out, ".partial")))), "A relaxed MR V1 freeze final or partial target exists; refusing to overwrite.")
dir.create(dirname(out[["manifest"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["log"]]), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = out[["log"]], append = TRUE)

summarise_loo <- function(loo) {
  sets <- sort(unique(loo$analysis_set))
  out <- lapply(sets, function(s) {
    z <- loo[loo$analysis_set == s, , drop = FALSE]
    i <- which.max(abs(as.numeric(z$absolute_shift)))
    data.frame(
      analysis_set = s,
      full_ivw_beta = unique(as.numeric(z$full_ivw_beta))[1],
      max_absolute_shift = as.numeric(z$absolute_shift[i]),
      max_shift_rsid = as.character(z$removed_rsid[i]),
      any_sign_change = any(bool_col(z$sign_change)),
      any_nominal_significance_change = any(bool_col(z$nominal_significance_change)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

summarise_single <- function(single) {
  sets <- sort(unique(single$analysis_set))
  out <- lapply(sets, function(s) {
    z <- single[single$analysis_set == s, , drop = FALSE]
    i <- which.max(abs(as.numeric(z$beta)))
    data.frame(
      analysis_set = s,
      max_abs_single_snp_beta = abs(as.numeric(z$beta[i])),
      max_abs_single_snp_rsid = as.character(z$rsid[i]),
      nominal_p_lt_0_05_count = sum(as.numeric(z$pval) < 0.05, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

compare_records <- function(csv_df, json_records, key_cols, num_cols, tol = 1e-4) {
  j <- as_df(json_records)
  if (nrow(csv_df) != nrow(j)) return(FALSE)
  csv_key <- do.call(paste, c(csv_df[key_cols], sep = "\r"))
  json_key <- do.call(paste, c(j[key_cols], sep = "\r"))
  csv_df <- csv_df[order(csv_key), , drop = FALSE]
  j <- j[order(json_key), , drop = FALSE]
  rownames(csv_df) <- NULL; rownames(j) <- NULL
  key_ok <- all(vapply(key_cols, function(k) identical(as.character(csv_df[[k]]), as.character(j[[k]])), logical(1)))
  num_ok <- all(vapply(num_cols, function(k) num_equal(csv_df[[k]], j[[k]], tol = tol), logical(1)))
  key_ok && num_ok
}

main <- function() {
  log_line("stage=reverse_relaxed_mr_v1_results_freeze")

  rel <- c(
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_freeze.json",
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_freeze_manifest.csv",
    "docs/decisions/72_reverse_relaxed_mr_analysis_contract_v1_v1.1.md",
    "results/qc/reverse_relaxed_mr_analysis_contract_v1.json",
    "results/qc/reverse_strict_primary_mr_v1_freeze.json",
    "results/qc/reverse_strict_primary_mr_v1_freeze_manifest.csv",
    "docs/decisions/76_reverse_relaxed_mr_v1_v1.1.md",
    "R/20_reverse_relaxed_mr_v1.R",
    "results/qc/reverse_relaxed_mr_v1.json",
    "results/logs/reverse_relaxed_mr_v1.log",
    "results/tables/reverse_relaxed_mr_estimates_v1.csv",
    "results/tables/reverse_relaxed_mr_estimates_doubling_odds_v1.csv",
    "results/tables/reverse_relaxed_heterogeneity_v1.csv",
    "results/tables/reverse_relaxed_egger_intercept_v1.csv",
    "results/tables/reverse_relaxed_mr_presso_v1.csv",
    "results/tables/reverse_relaxed_leave_one_out_v1.csv",
    "results/tables/reverse_relaxed_single_snp_v1.csv",
    "results/tables/reverse_relaxed_strict_comparison_v1.csv",
    decision_rel,
    script_rel,
    "renv.lock"
  )
  roles <- c(
    "decision_71_harmonisation_freeze_json",
    "decision_71_harmonisation_freeze_manifest",
    "decision_72_mr_contract_decision",
    "decision_72_mr_contract_json",
    "decision_75_strict_mr_freeze_json",
    "decision_75_strict_mr_freeze_manifest",
    "decision_76_reverse_relaxed_mr_decision",
    "reverse_relaxed_mr_v1_script",
    "reverse_relaxed_mr_v1_qc",
    "reverse_relaxed_mr_v1_log",
    "reverse_relaxed_mr_estimates",
    "reverse_relaxed_doubling_odds",
    "reverse_relaxed_heterogeneity",
    "reverse_relaxed_egger_intercept",
    "reverse_relaxed_mr_presso",
    "reverse_relaxed_leave_one_out",
    "reverse_relaxed_single_snp",
    "reverse_relaxed_strict_comparison",
    "decision_77_freeze_decision",
    "reverse_relaxed_mr_v1_freeze_script",
    "renv_lock"
  )
  paths <- file.path(root, rel)
  stop_if(any(!file.exists(paths)), paste("Missing relaxed MR freeze input(s):", paste(rel[!file.exists(paths)], collapse = "; ")))
  log_line("checkpoint=inputs_exist")

  harm <- jsonlite::fromJSON(file.path(root, rel[1]), simplifyVector = FALSE)
  contract <- jsonlite::fromJSON(file.path(root, rel[4]), simplifyVector = FALSE)
  strict <- jsonlite::fromJSON(file.path(root, rel[5]), simplifyVector = FALSE)
  mr_qc <- jsonlite::fromJSON(file.path(root, rel[9]), simplifyVector = FALSE)
  log_line("checkpoint=json_loaded")

  estimates <- read_table(rel[11])
  doubling <- read_table(rel[12])
  heterogeneity <- read_table(rel[13])
  egger <- read_table(rel[14])
  presso <- read_table(rel[15])
  loo <- read_table(rel[16])
  single <- read_table(rel[17])
  comparison <- read_table(rel[18])

  loo_summary <- summarise_loo(loo)
  single_summary <- summarise_single(single)
  renv_sha <- hash_file(file.path(root, "renv.lock"))
  log_line("checkpoint=tables_loaded_and_summarised")
  script_text <- paste(readLines(file.path(root, script_rel), warn = FALSE), collapse = "\n")
  forbidden_patterns <- c(
    "TwoSampleMR::mr\\s*\\(",
    "mr\\s*\\(",
    "mr_heterogeneity\\s*\\(",
    "mr_pleiotropy_test\\s*\\(",
    "mr_leaveoneout\\s*\\(",
    "mr_singlesnp\\s*\\(",
    "MRPRESSO::mr_presso\\s*\\(",
    "harmonise_data\\s*\\(",
    "mr_steiger\\s*\\(",
    "steiger_filtering\\s*\\(",
    "directionality_test\\s*\\("
  )
  no_forbidden_calls <- !any(vapply(forbidden_patterns, function(p) grepl(p, script_text, ignore.case = TRUE, perl = TRUE), logical(1)))
  log_line("checkpoint=static_scan_done no_forbidden_calls=", no_forbidden_calls)

  log_line("checkpoint=compare main_estimates start")
  main_estimates_match_qc <- compare_records(
    estimates, mr_qc$raw_results,
    key_cols = c("analysis_set", "method_id", "nsnp", "effect_scale", "estimator_hierarchy_role"),
    num_cols = c("beta", "se", "pval", "ci_lower", "ci_upper"),
    tol = 1e-4
  )
  log_line("checkpoint=compare main_estimates done value=", main_estimates_match_qc)
  log_line("checkpoint=compare heterogeneity start")
  heterogeneity_preserved <- compare_records(
    heterogeneity, mr_qc$heterogeneity_results,
    key_cols = c("analysis_set", "method"),
    num_cols = c("Q", "df", "pval"),
    tol = 1e-4
  )
  log_line("checkpoint=compare heterogeneity done value=", heterogeneity_preserved)
  log_line("checkpoint=compare egger start")
  egger_preserved <- compare_records(
    egger, mr_qc$egger_intercept_results,
    key_cols = c("analysis_set", "analysis_role"),
    num_cols = c("intercept", "se", "pval"),
    tol = 1e-4
  )
  log_line("checkpoint=compare egger done value=", egger_preserved)
  log_line("checkpoint=compare presso start")
  presso_preserved <- compare_records(
    presso, mr_qc$mr_presso_results,
    key_cols = c("analysis_set", "test_type", "metric", "value", "notes"),
    num_cols = c("pval"),
    tol = 1e-4
  )
  log_line("checkpoint=compare presso done value=", presso_preserved)
  log_line("checkpoint=compare loo start")
  loo_preserved <- compare_records(
    loo_summary, mr_qc$leave_one_out_summary,
    key_cols = c("analysis_set", "max_shift_rsid"),
    num_cols = c("max_absolute_shift"),
    tol = 1e-4
  )
  log_line("checkpoint=compare loo done value=", loo_preserved)
  log_line("checkpoint=compare single_snp start")
  single_preserved <- compare_records(
    single_summary, mr_qc$single_snp_summary,
    key_cols = c("analysis_set", "max_abs_single_snp_rsid", "nominal_p_lt_0_05_count"),
    num_cols = c("max_abs_single_snp_beta"),
    tol = 1e-4
  )
  log_line("checkpoint=compare single_snp done value=", single_preserved)
  log_line("checkpoint=qc_comparisons_done")

  raw_scale <- "standardized_Hb_per_1_unit_genetically_predicted_log_odds_delirium"
  secondary_scale <- "SD_Hb_per_doubling_of_genetically_predicted_odds_of_delirium"
  doubling_ok <- nrow(estimates) == nrow(doubling) &&
    all(estimates$analysis_set == doubling$analysis_set) &&
    all(estimates$method_id == doubling$method_id) &&
    num_equal(doubling$beta_per_doubling_odds, estimates$beta * log(2), tol = 1e-12) &&
    num_equal(doubling$se_per_doubling_odds, estimates$se * log(2), tol = 1e-12) &&
    num_equal(doubling$ci_lower_per_doubling, estimates$ci_lower * log(2), tol = 1e-12) &&
    num_equal(doubling$ci_upper_per_doubling, estimates$ci_upper * log(2), tol = 1e-12) &&
    num_equal(doubling$pval, estimates$pval, tol = 0)

  or_pattern <- "(^OR$|^OR_|_OR$|odds_ratio|exp\\(|exp_beta|beta_exponentiated)"
  no_or_transform <- !any(grepl(or_pattern, names(estimates), ignore.case = TRUE)) &&
    !any(grepl("g/dL|g/L|odds_ratio|exp\\(|exp_beta|beta_exponentiated", unlist(estimates), ignore.case = TRUE)) &&
    !any(grepl("g/dL|g/L|odds_ratio|exp\\(|exp_beta|beta_exponentiated", unlist(doubling), ignore.case = TRUE))

  strict_relaxed_hierarchy_preserved <- all(bool_col(comparison$strict_relaxed_direction_consistent)) &&
    all(comparison$interpretation_role == "interpretive_audit_only_no_meta_analysis_no_override")
  apoe_sets_not_independent_replication <- isTRUE(contract$apoe_sensitivity_interpretation$not_independent_replication)

  hard_checks <- list(
    relaxed_mr_qc_gate = identical(mr_qc$mr_status, "passed") &&
      isTRUE(mr_qc$approved_for_reverse_relaxed_results_interpretation) &&
      length(mr_qc$hard_check_failures) == 0L,
    harmonisation_freeze_gate = identical(harm$freeze_status, "passed") &&
      isTRUE(harm$approved_for_reverse_relaxed_mr_design) &&
      length(harm$hard_check_failures) == 0L &&
      identical(hash_file(file.path(root, rel[2])), harm$manifest_sha256),
    mr_contract_gate = identical(contract$contract_status, "frozen") &&
      isTRUE(contract$approved_for_reverse_relaxed_mr_execution) &&
      identical(contract$analysis_role, "secondary_reverse_exploratory_relaxed") &&
      identical(as.numeric(contract$p_threshold), 5e-6),
    strict_mr_freeze_gate = identical(strict$freeze_status, "passed") &&
      isTRUE(strict$approved_for_reverse_relaxed_mr_execution) &&
      length(strict$hard_check_failures) == 0L &&
      identical(hash_file(file.path(root, rel[6])), strict$manifest_sha256),
    main_estimates_match_qc = main_estimates_match_qc,
    doubling_odds_rescaling_verified = doubling_ok,
    continuous_outcome_scale_preserved = all(estimates$effect_scale == raw_scale) &&
      all(doubling$effect_scale_secondary == secondary_scale),
    no_or_transform = no_or_transform,
    heterogeneity_results_preserved = heterogeneity_preserved,
    egger_intercept_results_preserved = egger_preserved,
    mr_presso_results_preserved = presso_preserved,
    loo_results_preserved = loo_preserved,
    single_snp_results_preserved = single_preserved,
    no_posthoc_snp_filtering = as.integer(mr_qc$included_nsnp) == 10L &&
      as.integer(mr_qc$excluded_nsnp) == 9L &&
      sum(estimates$analysis_set == "APOE_included" & estimates$method_id == "mr_ivw" & estimates$nsnp == 10L) == 1L &&
      sum(estimates$analysis_set == "APOE_excluded" & estimates$method_id == "mr_ivw" & estimates$nsnp == 9L) == 1L,
    strict_relaxed_hierarchy_preserved = strict_relaxed_hierarchy_preserved,
    apoe_sets_not_independent_replication = apoe_sets_not_independent_replication,
    no_steiger = identical(mr_qc$steiger_run, FALSE) &&
      identical(mr_qc$steiger_status, "deferred_to_separate_directionality_sensitivity_contract"),
    software_provenance_complete = !is.null(mr_qc$software_environment$R_version) &&
      !is.null(mr_qc$software_environment$TwoSampleMR_version) &&
      nchar(mr_qc$software_environment$TwoSampleMR_RemoteSha) >= 40L &&
      !is.null(mr_qc$software_environment$MRPRESSO_version) &&
      nchar(mr_qc$software_environment$MRPRESSO_RemoteSha) >= 40L &&
      identical(as.integer(mr_qc$seed), 2026L) &&
      !is.null(mr_qc$input_sha256$decision_71_manifest) &&
      !is.null(mr_qc$input_sha256$decision_72_contract) &&
      !is.null(mr_qc$input_sha256$strict_mr_freeze_manifest),
    renv_lock_unchanged = isTRUE(mr_qc$renv_lock_unchanged) &&
      identical(mr_qc$renv_lock_sha_before, mr_qc$renv_lock_sha_after) &&
      identical(renv_sha, mr_qc$renv_lock_sha_after),
    no_mr_rerun = no_forbidden_calls
  )
  log_line("checkpoint=hard_checks_built")
  failures <- names(hard_checks)[!unlist(hard_checks)]
  status <- if (length(failures) == 0L) "passed" else "failed"

  classify_set <- function(set_name) {
    z <- estimates[estimates$analysis_set == set_name, , drop = FALSE]
    ivw <- z[z$method_id == "mr_ivw", , drop = FALSE]
    alt <- z[z$method_id != "mr_ivw", , drop = FALSE]
    hz <- heterogeneity[heterogeneity$analysis_set == set_name, , drop = FALSE]
    ez <- egger[egger$analysis_set == set_name, , drop = FALSE]
    pz <- presso[presso$analysis_set == set_name & presso$test_type == "Global Test", , drop = FALSE]
    lz <- loo_summary[loo_summary$analysis_set == set_name, , drop = FALSE]
    same_direction <- sign(as.numeric(alt$beta)) == sign(as.numeric(ivw$beta[1]))
    list(
      analysis_set = set_name,
      ivw_nominal_p_lt_0_05 = isTRUE(as.numeric(ivw$pval[1]) < 0.05),
      alternative_estimator_direction_consistency = sprintf("%d/%d alternative estimators match IVW direction", sum(same_direction), length(same_direction)),
      alternative_estimator_nominal_support_count = sum(as.numeric(alt$pval) < 0.05, na.rm = TRUE),
      loo_nominal_significance_change = isTRUE(lz$any_nominal_significance_change[1]),
      heterogeneity_evidence_detected = any(as.numeric(hz$pval) < 0.05, na.rm = TRUE),
      egger_intercept_evidence_detected = any(as.numeric(ez$pval) < 0.05, na.rm = TRUE),
      presso_global_evidence_detected = any(as.numeric(pz$pval) < 0.05, na.rm = TRUE),
      classification = if (as.numeric(ivw$pval[1]) < 0.05) "nominal_exploratory_signal" else "no_statistical_evidence_detected_in_exploratory_relaxed_analysis"
    )
  }
  exploratory_signal_classification <- list(
    classify_set("APOE_included"),
    classify_set("APOE_excluded")
  )

  manifest <- data.frame(
    file_role = roles,
    relative_path = rel,
    file_size_bytes = as.numeric(file.info(paths)$size),
    sha256 = vapply(paths, hash_file, character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(manifest, paste0(out[["manifest"]], ".partial"), row.names = FALSE)
  manifest_sha <- hash_file(paste0(out[["manifest"]], ".partial"))

  result <- list(
    freeze_version = "v1",
    decision = 77,
    authoritative_reverse_relaxed_mr_version = "v1",
    analysis_direction = "delirium_to_Hb",
    analysis_role = "secondary_reverse_exploratory_relaxed",
    branch_type = "protocol_triggered_exploratory_fallback",
    p_threshold = 5e-6,
    strict_primary_threshold = 5e-8,
    source_harmonisation_version = "v3",
    source_harmonisation_freeze_decision = 71,
    mr_contract_decision = 72,
    strict_primary_mr_freeze_decision = 75,
    source_mr_decision = 76,
    included_nsnp = as.integer(mr_qc$included_nsnp),
    excluded_nsnp = as.integer(mr_qc$excluded_nsnp),
    included_results = records(estimates[estimates$analysis_set == "APOE_included", , drop = FALSE]),
    excluded_results = records(estimates[estimates$analysis_set == "APOE_excluded", , drop = FALSE]),
    included_ivw = records(estimates[estimates$analysis_set == "APOE_included" & estimates$method_id == "mr_ivw", , drop = FALSE])[[1]],
    excluded_ivw = records(estimates[estimates$analysis_set == "APOE_excluded" & estimates$method_id == "mr_ivw", , drop = FALSE])[[1]],
    doubling_odds_results = records(doubling),
    heterogeneity_results = records(heterogeneity),
    egger_intercept_results = records(egger),
    mr_presso_results = records(presso),
    leave_one_out_summary = records(loo_summary),
    single_snp_summary = records(single_summary),
    strict_relaxed_comparison = records(comparison),
    exploratory_signal_classification = exploratory_signal_classification,
    apoe_interpretation = list(
      included_set_role = "reverse_relaxed_exploratory_main",
      excluded_set_role = "APOE-region exclusion sensitivity",
      independent_replication = FALSE,
      causal_contribution_of_APOE_not_inferred_from_set_difference = TRUE
    ),
    effect_scale = list(
      raw = raw_scale,
      raw_effect_interpretation = "standardized Hb change per 1-unit genetically predicted log odds of delirium",
      secondary = secondary_scale,
      no_or_transformation = TRUE,
      physical_unit_claim_allowed = FALSE
    ),
    steiger_run = FALSE,
    steiger_status = "deferred_to_separate_directionality_sensitivity_contract",
    strict_primary_superseded_by_relaxed = FALSE,
    relaxed_confirmatory = FALSE,
    software_environment = mr_qc$software_environment,
    seed = as.integer(mr_qc$seed),
    input_sha256 = list(
      decision_71_manifest = hash_file(file.path(root, rel[2])),
      decision_72_contract = hash_file(file.path(root, rel[4])),
      decision_75_strict_mr_freeze_manifest = hash_file(file.path(root, rel[6])),
      reverse_relaxed_mr_qc = hash_file(file.path(root, rel[9])),
      renv_lock = renv_sha
    ),
    manifest_sha256 = manifest_sha,
    freeze_status = status,
    approved_for_reverse_directionality_sensitivity_design = identical(status, "passed"),
    hard_checks = hard_checks,
    hard_check_failures = failures,
    informational_findings = list(
      confirmatory_status = "exploratory_not_confirmatory",
      relaxed_ivw_p_lt_0_05_may_only_be_reported_as_nominal_exploratory_signal = TRUE,
      strict_primary_evidence_tier_preserved = TRUE,
      apoe_sets_are_sensitivity_comparison_not_independent_replication = TRUE,
      mr_presso_no_outliers_does_not_trigger_snp_filtering = TRUE,
      steiger_deferred_to_separate_contract = TRUE,
      freeze_read_only_no_mr_functions_called = TRUE
    )
  )
  jsonlite::write_json(result, paste0(out[["json"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  if (!identical(status, "passed")) stop("Reverse relaxed MR V1 results freeze failed; partial outputs retained.", call. = FALSE)
  for (path in out[c("manifest", "json")]) {
    stop_if(file.exists(path), paste("Output appeared during run:", path))
    stop_if(!file.rename(paste0(path, ".partial"), path), paste("Atomic rename failed:", path))
  }
  log_line("freeze_status=passed approved_for_reverse_directionality_sensitivity_design=TRUE manifest_sha256=", manifest_sha)
}

tryCatch(main(), error = function(e) {
  log_line("status=failed error=", conditionMessage(e))
  quit(status = 1L)
})
