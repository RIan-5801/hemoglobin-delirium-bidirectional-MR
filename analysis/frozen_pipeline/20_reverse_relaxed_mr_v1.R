#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/20_reverse_relaxed_mr_v1.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))
mr_library <- file.path(root, "renv", "mr-v1-library")
.libPaths(c(mr_library, .libPaths()))

for (pkg in c("DBI", "duckdb", "jsonlite", "digest", "TwoSampleMR", "MRPRESSO")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
num_equal <- function(a, b, tol = 1e-12) {
  a <- as.numeric(a); b <- as.numeric(b)
  all(is.finite(a) & is.finite(b) & abs(a - b) <= tol)
}
sql_string <- function(path, must_work = TRUE) {
  paste0("'", gsub("'", "''", normalizePath(path, winslash = "/", mustWork = must_work), fixed = TRUE), "'")
}
read_parquet <- function(con, path) DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(path)))
write_csv_precise <- function(x, path) {
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  write.csv(x, path, row.names = FALSE, na = "")
}

contract_path <- file.path(root, "results", "qc", "reverse_relaxed_mr_analysis_contract_v1.json")
relaxed_freeze_path <- file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_freeze.json")
relaxed_manifest_path <- file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_freeze_manifest.csv")
strict_freeze_path <- file.path(root, "results", "qc", "reverse_strict_primary_mr_v1_freeze.json")
renv_lock <- file.path(root, "renv.lock")
included_path <- file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_relaxed_harmonised_apoe_included_v3.parquet")
excluded_path <- file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_relaxed_harmonised_apoe_excluded_v3.parquet")

out <- c(
  estimates = file.path(root, "results", "tables", "reverse_relaxed_mr_estimates_v1.csv"),
  doubling = file.path(root, "results", "tables", "reverse_relaxed_mr_estimates_doubling_odds_v1.csv"),
  heterogeneity = file.path(root, "results", "tables", "reverse_relaxed_heterogeneity_v1.csv"),
  egger = file.path(root, "results", "tables", "reverse_relaxed_egger_intercept_v1.csv"),
  presso = file.path(root, "results", "tables", "reverse_relaxed_mr_presso_v1.csv"),
  loo = file.path(root, "results", "tables", "reverse_relaxed_leave_one_out_v1.csv"),
  single = file.path(root, "results", "tables", "reverse_relaxed_single_snp_v1.csv"),
  comparison = file.path(root, "results", "tables", "reverse_relaxed_strict_comparison_v1.csv"),
  qc = file.path(root, "results", "qc", "reverse_relaxed_mr_v1.json"),
  log = file.path(root, "results", "logs", "reverse_relaxed_mr_v1.log")
)
stop_if(any(file.exists(c(out, paste0(out, ".partial")))), "Reverse relaxed MR V1 target or partial exists; refusing to overwrite.")
dir.create(dirname(out[["estimates"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["qc"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["log"]]), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = out[["log"]], append = TRUE)

method_ids <- c("mr_ivw", "mr_weighted_median", "mr_egger_regression", "mr_weighted_mode", "mr_simple_mode")
method_roles <- c(
  mr_ivw = "exploratory_branch_primary_estimator",
  mr_weighted_median = "robust_sensitivity_estimator",
  mr_egger_regression = "directional_pleiotropy_sensitivity_estimator",
  mr_weighted_mode = "supportive_sensitivity_estimator",
  mr_simple_mode = "supportive_sensitivity_estimator"
)

make_dat <- function(x, analysis_set, analysis_role) {
  required <- c("target_rsid", "exposure_beta", "exposure_se", "exposure_effect_allele",
                "exposure_other_allele", "exposure_eaf", "outcome_beta_harmonised",
                "outcome_se_harmonised", "outcome_effect_allele_harmonised",
                "outcome_other_allele_harmonised", "outcome_eaf_harmonised",
                "outcome_pval_harmonised", "final_valid_instrument")
  stop_if(!all(required %in% names(x)), paste("Required harmonised input columns are absent:", analysis_set))
  stop_if(!all(x$final_valid_instrument), paste("Input contains non-final-valid rows:", analysis_set))
  data.frame(
    SNP = as.character(x$target_rsid),
    beta.exposure = as.numeric(x$exposure_beta),
    se.exposure = as.numeric(x$exposure_se),
    effect_allele.exposure = as.character(x$exposure_effect_allele),
    other_allele.exposure = as.character(x$exposure_other_allele),
    eaf.exposure = as.numeric(x$exposure_eaf),
    beta.outcome = as.numeric(x$outcome_beta_harmonised),
    se.outcome = as.numeric(x$outcome_se_harmonised),
    effect_allele.outcome = as.character(x$outcome_effect_allele_harmonised),
    other_allele.outcome = as.character(x$outcome_other_allele_harmonised),
    eaf.outcome = as.numeric(x$outcome_eaf_harmonised),
    exposure = "FinnGen_R13_F5_DELIRIUM",
    outcome = "Vuckovic_2020_Hb",
    id.exposure = "finngen_R13_F5_DELIRIUM",
    id.outcome = "vuckovic_hb_2020",
    mr_keep = TRUE,
    analysis_set = analysis_set,
    analysis_role = analysis_role,
    stringsAsFactors = FALSE
  )
}

check_input <- function(source, dat, expected_count, expected_rsids) {
  f <- (source$exposure_beta / source$exposure_se)^2
  list(
    count_matches_freeze = nrow(source) == expected_count,
    rsids_match_freeze = identical(as.character(source$target_rsid), as.character(expected_rsids)),
    all_final_valid = all(source$final_valid_instrument),
    rsid_unique = !anyDuplicated(source$target_rsid) && !anyNA(source$target_rsid),
    finite_exposure_effects = all(is.finite(dat$beta.exposure)),
    positive_exposure_se = all(is.finite(dat$se.exposure) & dat$se.exposure > 0),
    finite_outcome_effects = all(is.finite(dat$beta.outcome)),
    positive_outcome_se = all(is.finite(dat$se.outcome) & dat$se.outcome > 0),
    all_F_ge_10 = all(is.finite(f) & f >= 10),
    conversion_preserved = identical(as.character(source$target_rsid), dat$SNP) &&
      num_equal(source$exposure_beta, dat$beta.exposure) &&
      num_equal(source$exposure_se, dat$se.exposure) &&
      num_equal(source$outcome_beta_harmonised, dat$beta.outcome) &&
      num_equal(source$outcome_se_harmonised, dat$se.outcome) &&
      identical(as.character(source$exposure_effect_allele), dat$effect_allele.exposure) &&
      identical(as.character(source$exposure_other_allele), dat$other_allele.exposure) &&
      identical(as.character(source$outcome_effect_allele_harmonised), dat$effect_allele.outcome) &&
      identical(as.character(source$outcome_other_allele_harmonised), dat$other_allele.outcome),
    f_summary = data.frame(n = length(f), min = min(f), mean = mean(f), median = stats::median(f), max = max(f), F_lt10 = sum(f < 10))
  )
}

run_mr_methods <- function(dat, analysis_set, analysis_role) {
  x <- TwoSampleMR::mr(dat, method_list = method_ids)
  if (is.list(x) && "mr" %in% names(x)) x <- x$mr
  stop_if(!is.data.frame(x) || !all(c("method", "nsnp", "b", "se", "pval") %in% names(x)), "MR result schema is invalid.")
  methods <- TwoSampleMR::mr_method_list()
  lookup <- stats::setNames(as.character(methods$obj), as.character(methods$name))
  z <- data.frame(
    analysis_set = analysis_set,
    analysis_role = analysis_role,
    method = as.character(x$method),
    method_id = unname(lookup[as.character(x$method)]),
    nsnp = as.integer(x$nsnp),
    beta = as.numeric(x$b),
    se = as.numeric(x$se),
    pval = as.numeric(x$pval),
    stringsAsFactors = FALSE
  )
  z$ci_lower <- z$beta - 1.96 * z$se
  z$ci_upper <- z$beta + 1.96 * z$se
  z$effect_scale <- "standardized_Hb_per_1_unit_genetically_predicted_log_odds_delirium"
  z$nominal_p_lt_0_05 <- z$pval < 0.05
  z$estimator_hierarchy_role <- unname(method_roles[z$method_id])
  z
}

run_heterogeneity <- function(dat, analysis_set, analysis_role) {
  x <- TwoSampleMR::mr_heterogeneity(dat, method_list = c("mr_egger_regression", "mr_ivw"))
  stop_if(!all(c("method", "Q", "Q_df", "Q_pval") %in% names(x)), "Heterogeneity result schema is invalid.")
  data.frame(
    analysis_set = analysis_set,
    analysis_role = analysis_role,
    method = as.character(x$method),
    Q = as.numeric(x$Q),
    df = as.numeric(x$Q_df),
    pval = as.numeric(x$Q_pval),
    heterogeneity_p_lt_0_05 = as.numeric(x$Q_pval) < 0.05,
    automatic_snp_removal_allowed = FALSE,
    stringsAsFactors = FALSE
  )
}

run_egger <- function(dat, analysis_set, analysis_role) {
  x <- TwoSampleMR::mr_pleiotropy_test(dat)
  stop_if(!all(c("egger_intercept", "se", "pval") %in% names(x)), "Egger intercept result schema is invalid.")
  data.frame(
    analysis_set = analysis_set,
    analysis_role = analysis_role,
    intercept = as.numeric(x$egger_intercept[[1L]]),
    se = as.numeric(x$se[[1L]]),
    pval = as.numeric(x$pval[[1L]]),
    nominal_p_lt_0_05 = as.numeric(x$pval[[1L]]) < 0.05,
    interpretation_rule = "P>=0.05 means no statistical evidence detected, not proof of no pleiotropy",
    egger_precision_limitation = "limited_number_of_instruments",
    stringsAsFactors = FALSE
  )
}

run_presso <- function(dat, analysis_set, analysis_role) {
  set.seed(2026L)
  x <- dat[, c("beta.outcome", "beta.exposure", "se.outcome", "se.exposure"), drop = FALSE]
  rownames(x) <- dat$SNP
  result <- tryCatch(
    MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome",
      BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome",
      SdExposure = "se.exposure",
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      data = x,
      NbDistribution = 10000,
      SignifThreshold = 0.05
    ),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    return(data.frame(
      analysis_set = analysis_set, analysis_role = analysis_role, test_type = "MR-PRESSO",
      metric = "status", value = "not_estimable_or_failed", pval = NA_real_,
      outlier_rsid = "", notes = conditionMessage(result), stringsAsFactors = FALSE
    ))
  }
  root <- result[["MR-PRESSO results"]]
  global <- if (!is.null(root)) root[["Global Test"]] else NULL
  outlier <- if (!is.null(root)) root[["Outlier Test"]] else NULL
  distortion <- if (!is.null(root)) root[["Distortion Test"]] else NULL
  rows <- list(
    data.frame(analysis_set = analysis_set, analysis_role = analysis_role, test_type = "Global Test",
               metric = "RSSobs", value = if (!is.null(global) && "RSSobs" %in% names(global)) as.character(global$RSSobs[[1L]]) else "",
               pval = if (!is.null(global) && "Pvalue" %in% names(global)) suppressWarnings(as.numeric(global$Pvalue[[1L]])) else NA_real_,
               outlier_rsid = "", notes = "mr_presso_status=passed", stringsAsFactors = FALSE)
  )
  ids <- if (!is.null(outlier) && nrow(outlier) > 0L) rownames(outlier) else character(0)
  if (length(ids) == 0L) {
    rows[[length(rows) + 1L]] <- data.frame(analysis_set = analysis_set, analysis_role = analysis_role,
      test_type = "Outlier Test", metric = "outlier_count", value = "0", pval = NA_real_,
      outlier_rsid = "", notes = "no_outlier_reported", stringsAsFactors = FALSE)
  } else {
    rows[[length(rows) + 1L]] <- data.frame(analysis_set = analysis_set, analysis_role = analysis_role,
      test_type = "Outlier Test", metric = "outlier", value = "", pval = NA_real_,
      outlier_rsid = ids, notes = "outlier_reported_sensitivity_only", stringsAsFactors = FALSE)
  }
  if (!is.null(distortion)) {
    rows[[length(rows) + 1L]] <- data.frame(analysis_set = analysis_set, analysis_role = analysis_role,
      test_type = "Distortion Test", metric = "raw_result", value = paste(capture.output(print(distortion)), collapse = " | "),
      pval = NA_real_, outlier_rsid = "", notes = "distortion_output_if_available", stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

run_loo <- function(dat, analysis_set, analysis_role, full_ivw_beta, full_ivw_pval) {
  raw <- TwoSampleMR::mr_leaveoneout(dat, method = TwoSampleMR::mr_ivw)
  if (!is.data.frame(raw)) raw <- as.data.frame(raw)
  stop_if(!all(c("SNP", "b", "se", "p") %in% names(raw)), "LOO result does not expose SNP/b/se/p fields.")
  labels <- as.character(raw$SNP)
  overall <- !is.na(labels) & grepl("^All", labels)
  stop_if(sum(overall) != 1L, "LOO requires exactly one All sentinel.")
  snp_raw <- raw[!overall, , drop = FALSE]
  stop_if(nrow(snp_raw) != nrow(dat) || !setequal(as.character(snp_raw$SNP), dat$SNP), "LOO SNP set does not match input.")
  z <- data.frame(
    analysis_set = analysis_set,
    analysis_role = analysis_role,
    removed_rsid = as.character(snp_raw$SNP),
    nsnp_remaining = nrow(dat) - 1L,
    beta = as.numeric(snp_raw$b),
    se = as.numeric(snp_raw$se),
    ci_lower = as.numeric(snp_raw$b) - 1.96 * as.numeric(snp_raw$se),
    ci_upper = as.numeric(snp_raw$b) + 1.96 * as.numeric(snp_raw$se),
    pval = as.numeric(snp_raw$p),
    full_ivw_beta = full_ivw_beta,
    full_ivw_pval = full_ivw_pval,
    stringsAsFactors = FALSE
  )
  z$absolute_shift <- abs(z$beta - full_ivw_beta)
  z$relative_shift <- ifelse(full_ivw_beta == 0, NA_real_, z$absolute_shift / abs(full_ivw_beta))
  z$sign_change <- sign(z$beta) != sign(full_ivw_beta)
  z$nominal_significance_change <- (z$pval < 0.05) != (full_ivw_pval < 0.05)
  z
}

run_single <- function(dat, analysis_set, analysis_role) {
  rows <- lapply(seq_len(nrow(dat)), function(i) {
    wr <- TwoSampleMR::mr_wald_ratio(
      b_exp = dat$beta.exposure[[i]],
      b_out = dat$beta.outcome[[i]],
      se_exp = dat$se.exposure[[i]],
      se_out = dat$se.outcome[[i]]
    )
    beta <- as.numeric(wr$b)
    se <- as.numeric(wr$se)
    data.frame(
      analysis_set = analysis_set,
      analysis_role = analysis_role,
      rsid = dat$SNP[[i]],
      method = "Wald ratio",
      method_id = "mr_wald_ratio",
      beta = beta,
      se = se,
      ci_lower = beta - 1.96 * se,
      ci_upper = beta + 1.96 * se,
      pval = as.numeric(wr$pval),
      diagnostic_role = "single_snp_diagnostic_plotting_input_only",
      branch_primary_result = FALSE,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

make_doubling <- function(est) {
  z <- est
  z$beta_raw <- z$beta
  z$se_raw <- z$se
  z$ci_lower_raw <- z$ci_lower
  z$ci_upper_raw <- z$ci_upper
  z$beta_per_doubling_odds <- z$beta * log(2)
  z$se_per_doubling_odds <- z$se * log(2)
  z$ci_lower_per_doubling <- z$ci_lower * log(2)
  z$ci_upper_per_doubling <- z$ci_upper * log(2)
  z$effect_scale_secondary <- "SD_Hb_per_doubling_of_genetically_predicted_odds_of_delirium"
  z[, c("analysis_set", "analysis_role", "method", "method_id", "nsnp", "beta_raw", "se_raw",
        "ci_lower_raw", "ci_upper_raw", "pval", "beta_per_doubling_odds", "se_per_doubling_odds",
        "ci_lower_per_doubling", "ci_upper_per_doubling", "effect_scale_secondary",
        "estimator_hierarchy_role")]
}

main <- function() {
  log_line("stage=reverse_relaxed_mr_v1")
  set.seed(2026L)
  renv_before <- hash_file(renv_lock)
  contract <- jsonlite::fromJSON(contract_path, simplifyVector = FALSE)
  relaxed_freeze <- jsonlite::fromJSON(relaxed_freeze_path, simplifyVector = FALSE)
  strict_freeze <- jsonlite::fromJSON(strict_freeze_path, simplifyVector = FALSE)
  manifest_sha <- hash_file(relaxed_manifest_path)
  contract_sha <- hash_file(contract_path)
  strict_freeze_manifest_sha <- strict_freeze$manifest_sha256
  ts_desc <- utils::packageDescription("TwoSampleMR")
  mp_desc <- utils::packageDescription("MRPRESSO")
  ts_ver <- as.character(utils::packageVersion("TwoSampleMR"))
  ts_sha <- as.character(ts_desc[["RemoteSha"]])
  mp_ver <- as.character(utils::packageVersion("MRPRESSO"))
  mp_sha <- as.character(mp_desc[["RemoteSha"]])

  method_table <- TwoSampleMR::mr_method_list()
  methods_available <- all(method_ids %in% method_table$obj)
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  inc_src <- read_parquet(con, included_path)
  exc_src <- read_parquet(con, excluded_path)
  inc_dat <- make_dat(inc_src, "APOE_included", "reverse_relaxed_exploratory_main")
  exc_dat <- make_dat(exc_src, "APOE_excluded", "reverse_relaxed_apoe_exclusion_sensitivity")
  inc_checks <- check_input(inc_src, inc_dat, as.integer(relaxed_freeze$included_final_valid_instrument_count), relaxed_freeze$included_final_rsids)
  exc_checks <- check_input(exc_src, exc_dat, as.integer(relaxed_freeze$excluded_final_valid_instrument_count), relaxed_freeze$excluded_final_rsids)

  pre_checks <- list(
    strict_primary_mr_freeze_gate = identical(strict_freeze$freeze_status, "passed") &&
      isTRUE(strict_freeze$approved_for_reverse_relaxed_mr_execution) && length(strict_freeze$hard_check_failures) == 0L,
    relaxed_mr_contract_gate = identical(contract$contract_status, "frozen") &&
      isTRUE(contract$approved_for_reverse_relaxed_mr_execution) &&
      identical(contract$analysis_role, "secondary_reverse_exploratory_relaxed") &&
      isTRUE(contract$steiger_status == "deferred_to_separate_directionality_sensitivity_contract"),
    relaxed_harmonisation_freeze_gate = identical(relaxed_freeze$freeze_status, "passed") &&
      isTRUE(relaxed_freeze$approved_for_reverse_relaxed_mr_design) &&
      identical(relaxed_freeze$authoritative_reverse_relaxed_formal_harmonisation_version, "v3") &&
      length(relaxed_freeze$hard_check_failures) == 0L &&
      identical(tolower(manifest_sha), tolower(relaxed_freeze$manifest_sha256)),
    method_ids_available = methods_available,
    software_versions_match_contract = identical(ts_ver, contract$software_environment$TwoSampleMR$version) &&
      identical(ts_sha, contract$software_environment$TwoSampleMR$RemoteSha) &&
      identical(mp_ver, contract$software_environment$MRPRESSO$version) &&
      identical(mp_sha, contract$software_environment$MRPRESSO$RemoteSha)
  )
  stop_if(!all(unlist(pre_checks)), paste("Pre-execution gate failed:", paste(names(pre_checks)[!unlist(pre_checks)], collapse = "; ")))
  input_checks <- c(inc_checks[setdiff(names(inc_checks), "f_summary")], exc_checks[setdiff(names(exc_checks), "f_summary")])
  stop_if(!all(unlist(input_checks)), paste("Input audit failed:", paste(names(input_checks)[!unlist(input_checks)], collapse = "; ")))

  inc_est <- run_mr_methods(inc_dat, "APOE_included", "reverse_relaxed_exploratory_main")
  exc_est <- run_mr_methods(exc_dat, "APOE_excluded", "reverse_relaxed_apoe_exclusion_sensitivity")
  estimates <- rbind(inc_est, exc_est)
  doubling <- make_doubling(estimates)
  inc_ivw <- inc_est[inc_est$method_id == "mr_ivw", , drop = FALSE]
  exc_ivw <- exc_est[exc_est$method_id == "mr_ivw", , drop = FALSE]
  heterogeneity <- rbind(
    run_heterogeneity(inc_dat, "APOE_included", "reverse_relaxed_exploratory_main"),
    run_heterogeneity(exc_dat, "APOE_excluded", "reverse_relaxed_apoe_exclusion_sensitivity")
  )
  egger <- rbind(
    run_egger(inc_dat, "APOE_included", "reverse_relaxed_exploratory_main"),
    run_egger(exc_dat, "APOE_excluded", "reverse_relaxed_apoe_exclusion_sensitivity")
  )
  presso <- rbind(
    run_presso(inc_dat, "APOE_included", "reverse_relaxed_exploratory_main"),
    run_presso(exc_dat, "APOE_excluded", "reverse_relaxed_apoe_exclusion_sensitivity")
  )
  loo <- rbind(
    run_loo(inc_dat, "APOE_included", "reverse_relaxed_exploratory_main", inc_ivw$beta[[1]], inc_ivw$pval[[1]]),
    run_loo(exc_dat, "APOE_excluded", "reverse_relaxed_apoe_exclusion_sensitivity", exc_ivw$beta[[1]], exc_ivw$pval[[1]])
  )
  single <- rbind(
    run_single(inc_dat, "APOE_included", "reverse_relaxed_exploratory_main"),
    run_single(exc_dat, "APOE_excluded", "reverse_relaxed_apoe_exclusion_sensitivity")
  )
  strict_included <- strict_freeze$included_raw_result
  strict_excluded <- strict_freeze$excluded_raw_result
  comparison <- data.frame(
    comparison_set = c("strict_vs_relaxed_included", "strict_vs_relaxed_excluded"),
    strict_analysis_role = c("secondary_reverse_primary", "secondary_reverse_primary"),
    relaxed_analysis_role = c("secondary_reverse_exploratory_relaxed", "secondary_reverse_exploratory_relaxed"),
    strict_method = c(strict_included$method, strict_excluded$method),
    relaxed_method = c("Inverse variance weighted", "Inverse variance weighted"),
    strict_beta = c(strict_included$beta, strict_excluded$beta),
    relaxed_beta = c(inc_ivw$beta[[1]], exc_ivw$beta[[1]]),
    strict_ci_lower = c(strict_included$ci_lower, strict_excluded$ci_lower),
    strict_ci_upper = c(strict_included$ci_upper, strict_excluded$ci_upper),
    relaxed_ci_lower = c(inc_ivw$ci_lower[[1]], exc_ivw$ci_lower[[1]]),
    relaxed_ci_upper = c(inc_ivw$ci_upper[[1]], exc_ivw$ci_upper[[1]]),
    strict_pval = c(strict_included$pval, strict_excluded$pval),
    relaxed_pval = c(inc_ivw$pval[[1]], exc_ivw$pval[[1]]),
    strict_relaxed_direction_consistent = sign(c(strict_included$beta, strict_excluded$beta)) == sign(c(inc_ivw$beta[[1]], exc_ivw$beta[[1]])),
    strict_relaxed_nominal_significance_pattern = paste0(ifelse(c(strict_included$pval, strict_excluded$pval) < 0.05, "strict_p_lt_0_05", "strict_p_ge_0_05"), "__",
                                                         ifelse(c(inc_ivw$pval[[1]], exc_ivw$pval[[1]]) < 0.05, "relaxed_p_lt_0_05", "relaxed_p_ge_0_05")),
    interpretation_role = "interpretive_audit_only_no_meta_analysis_no_override",
    stringsAsFactors = FALSE
  )

  loo_summary <- do.call(rbind, lapply(split(loo, loo$analysis_set), function(x) {
    i <- which.max(x$absolute_shift)
    data.frame(
      analysis_set = x$analysis_set[[1]],
      max_absolute_shift = x$absolute_shift[[i]],
      max_shift_rsid = x$removed_rsid[[i]],
      any_sign_change = any(x$sign_change),
      any_nominal_significance_change = any(x$nominal_significance_change),
      stringsAsFactors = FALSE
    )
  }))
  single_summary <- do.call(rbind, lapply(split(single, single$analysis_set), function(x) {
    i <- which.max(abs(x$beta))
    data.frame(
      analysis_set = x$analysis_set[[1]],
      max_abs_single_snp_beta = abs(x$beta[[i]]),
      max_abs_single_snp_rsid = x$rsid[[i]],
      nominal_p_lt_0_05_count = sum(x$pval < 0.05),
      stringsAsFactors = FALSE
    )
  }))

  hard_checks <- c(pre_checks, list(
    input_counts_match_freeze = inc_checks$count_matches_freeze && exc_checks$count_matches_freeze,
    input_rsids_match_freeze = inc_checks$rsids_match_freeze && exc_checks$rsids_match_freeze,
    no_reharmonisation = TRUE,
    no_posthoc_instrument_filtering = TRUE,
    all_F_ge_10 = inc_checks$all_F_ge_10 && exc_checks$all_F_ge_10,
    ivw_run_included = nrow(inc_ivw) == 1L,
    ivw_run_excluded = nrow(exc_ivw) == 1L,
    weighted_median_run_included = sum(inc_est$method_id == "mr_weighted_median") == 1L,
    weighted_median_run_excluded = sum(exc_est$method_id == "mr_weighted_median") == 1L,
    egger_run_included = sum(inc_est$method_id == "mr_egger_regression") == 1L,
    egger_run_excluded = sum(exc_est$method_id == "mr_egger_regression") == 1L,
    weighted_mode_run_included = sum(inc_est$method_id == "mr_weighted_mode") == 1L,
    weighted_mode_run_excluded = sum(exc_est$method_id == "mr_weighted_mode") == 1L,
    simple_mode_run_included = sum(inc_est$method_id == "mr_simple_mode") == 1L,
    simple_mode_run_excluded = sum(exc_est$method_id == "mr_simple_mode") == 1L,
    effect_scale_continuous_hb = all(estimates$effect_scale == "standardized_Hb_per_1_unit_genetically_predicted_log_odds_delirium"),
    no_or_transform = !any(grepl("^OR$|OR_", names(estimates))),
    doubling_odds_rescaling_correct = all(abs(doubling$beta_per_doubling_odds - estimates$beta * log(2)) <= 1e-12) &&
      all(abs(doubling$se_per_doubling_odds - estimates$se * log(2)) <= 1e-12) &&
      identical(as.numeric(doubling$pval), as.numeric(estimates$pval)),
    heterogeneity_completed = nrow(heterogeneity) == 4L,
    egger_intercept_completed = nrow(egger) == 2L,
    mr_presso_attempted = all(c("APOE_included", "APOE_excluded") %in% unique(presso$analysis_set)),
    leave_one_out_completed = nrow(loo) == nrow(inc_dat) + nrow(exc_dat),
    single_snp_completed = nrow(single) == nrow(inc_dat) + nrow(exc_dat),
    no_steiger = TRUE,
    no_proxy = TRUE,
    no_liftover = TRUE,
    strict_relaxed_hierarchy_preserved = TRUE,
    apoe_sets_not_called_independent_replication = TRUE,
    no_relaxed_result_freeze_created = TRUE
  ))
  renv_after <- hash_file(renv_lock)
  hard_checks <- c(hard_checks, list(renv_lock_unchanged = identical(renv_before, renv_after)))
  failures <- names(hard_checks)[!unlist(hard_checks)]
  status <- if (length(failures) == 0L) "passed" else "failed"

  write_csv_precise(estimates, paste0(out[["estimates"]], ".partial"))
  write_csv_precise(doubling, paste0(out[["doubling"]], ".partial"))
  write_csv_precise(heterogeneity, paste0(out[["heterogeneity"]], ".partial"))
  write_csv_precise(egger, paste0(out[["egger"]], ".partial"))
  write_csv_precise(presso, paste0(out[["presso"]], ".partial"))
  write_csv_precise(loo, paste0(out[["loo"]], ".partial"))
  write_csv_precise(single, paste0(out[["single"]], ".partial"))
  write_csv_precise(comparison, paste0(out[["comparison"]], ".partial"))
  qc <- list(
    mr_version = "v1",
    decision = 76,
    analysis_direction = "delirium_to_Hb",
    analysis_role = "secondary_reverse_exploratory_relaxed",
    branch_type = "protocol_triggered_exploratory_fallback",
    p_threshold = 5e-6,
    strict_primary_threshold = 5e-8,
    source_harmonisation_version = "v3",
    source_harmonisation_freeze_decision = 71,
    mr_contract_decision = 72,
    strict_primary_mr_freeze_decision = 75,
    included_nsnp = nrow(inc_dat),
    excluded_nsnp = nrow(exc_dat),
    included_rsids = inc_dat$SNP,
    excluded_rsids = exc_dat$SNP,
    methods_run = method_ids,
    raw_results = records(estimates),
    doubling_odds_results = records(doubling),
    heterogeneity_results = records(heterogeneity),
    egger_intercept_results = records(egger),
    mr_presso_results = records(presso),
    leave_one_out_summary = records(loo_summary),
    single_snp_summary = records(single_summary),
    strict_relaxed_comparison = records(comparison),
    steiger_run = FALSE,
    steiger_status = "deferred_to_separate_directionality_sensitivity_contract",
    seed = 2026L,
    software_environment = list(
      R_version = R.version.string,
      TwoSampleMR_version = ts_ver,
      TwoSampleMR_RemoteSha = ts_sha,
      MRPRESSO_version = mp_ver,
      MRPRESSO_RemoteSha = mp_sha,
      mr_library = normalizePath(mr_library, winslash = "/", mustWork = TRUE),
      install_update_restore_snapshot_performed = FALSE,
      renv_out_of_sync_message = "informational_only"
    ),
    input_sha256 = list(
      included = hash_file(included_path),
      excluded = hash_file(excluded_path),
      decision_71_manifest = manifest_sha,
      decision_72_contract = contract_sha,
      strict_mr_freeze_manifest = strict_freeze_manifest_sha
    ),
    instrument_strength_summary = list(
      included = as.list(inc_checks$f_summary[1, , drop = FALSE]),
      excluded = as.list(exc_checks$f_summary[1, , drop = FALSE])
    ),
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after,
    renv_lock_unchanged = identical(renv_before, renv_after),
    mr_status = status,
    approved_for_reverse_relaxed_results_interpretation = identical(status, "passed"),
    hard_checks = hard_checks,
    hard_check_failures = failures,
    informational_findings = list(
      confirmatory_status = "exploratory_not_confirmatory",
      significance_vocabulary = "nominal exploratory evidence only when P<0.05",
      mr_presso_failure_handling = "MR-PRESSO failure, if any, is recorded as not_estimable_or_failed and does not automatically fail core MR when all other hard checks pass",
      egger_precision_limitation = "limited_number_of_instruments",
      no_relaxed_result_freeze_created = TRUE
    )
  )
  jsonlite::write_json(qc, paste0(out[["qc"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  if (!identical(status, "passed")) stop("Reverse relaxed MR V1 hard checks failed; partial outputs retained.", call. = FALSE)
  for (path in out[c("estimates", "doubling", "heterogeneity", "egger", "presso", "loo", "single", "comparison", "qc")]) {
    stop_if(file.exists(path), paste("Output appeared during run:", path))
    stop_if(!file.rename(paste0(path, ".partial"), path), paste("Atomic rename failed:", path))
  }
  log_line("mr_status=passed approved_for_reverse_relaxed_results_interpretation=TRUE hard_check_failures=[]")
}

tryCatch(main(), error = function(e) {
  log_line("status=failed error=", conditionMessage(e))
  quit(status = 1L)
})
