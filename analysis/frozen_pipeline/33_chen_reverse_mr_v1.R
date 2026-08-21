#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/33_chen_reverse_mr_v1.R [--project-root <path>]", call. = FALSE)
}
setwd(root)
source(file.path(root, "renv", "activate.R"))
mr_library <- file.path(root, "renv", "mr-v1-library")
.libPaths(c(mr_library, .libPaths()))

for (pkg in c("DBI", "duckdb", "jsonlite", "digest", "TwoSampleMR", "MRPRESSO")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

rel <- function(...) file.path(root, ...)
norm <- function(x) gsub("\\\\", "/", x)
relpath <- function(x) norm(sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", norm(normalizePath(x, winslash = "/", mustWork = FALSE))))
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
hash_text <- function(x) digest::digest(x, algo = "sha256")
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
is_empty <- function(x) is.null(x) || length(x) == 0L
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
records <- function(x) if (!is.data.frame(x) || nrow(x) == 0L) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
num_equal <- function(a, b, tol = 1e-12) {
  a <- as.numeric(a); b <- as.numeric(b)
  length(a) == length(b) && all(is.finite(a) & is.finite(b) & abs(a - b) <= tol)
}
sql_string <- function(path, must_work = TRUE) {
  paste0("'", gsub("'", "''", normalizePath(path, winslash = "/", mustWork = must_work), fixed = TRUE), "'")
}
read_parquet <- function(con, path) DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(path)))

atomic_write <- function(path, writer) {
  partial <- paste0(path, ".partial")
  stop_if(file.exists(path) || file.exists(partial), paste("Output occupied:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writer(partial)
  stop_if(!file.exists(partial), paste("Writer did not create partial:", partial))
  stop_if(!file.rename(partial, path), paste("Atomic rename failed:", path))
}
write_csv_precise <- function(x, path) atomic_write(path, function(p) {
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.csv(x, p, row.names = FALSE, na = "")
})
write_json <- function(x, path) atomic_write(path, function(p) {
  jsonlite::write_json(x, p, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
})
write_text <- function(lines, path) atomic_write(path, function(p) writeLines(lines, p, useBytes = TRUE))

split_rsids <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(x)) character(0) else strsplit(as.character(x), ";", fixed = TRUE)[[1L]]
}
branch_from_freeze <- function(freeze, branch_name) {
  hits <- Filter(function(x) identical(x$branch, branch_name), freeze$branch_results)
  stop_if(length(hits) != 1L, paste("Branch not found exactly once in Decision 116:", branch_name))
  hits[[1L]]
}
method_ids_from_contract <- function(contract) {
  ids <- vapply(contract$relaxed_estimator_hierarchy, function(x) if (!is.null(x$method_id)) x$method_id else NA_character_, character(1))
  ids[!is.na(ids)]
}
method_roles_from_contract <- function(contract) {
  ids <- method_ids_from_contract(contract)
  roles <- vapply(contract$relaxed_estimator_hierarchy, function(x) if (!is.null(x$role)) x$role else NA_character_, character(1))
  stats::setNames(roles[!is.na(ids)], ids)
}
branch_role <- function(branch) {
  if (identical(branch, "strict_apoe_included")) "reverse_strict_primary_alternative_hb_outcome_sensitivity"
  else if (identical(branch, "strict_apoe_excluded")) "reverse_strict_primary_alternative_hb_outcome_sensitivity"
  else if (identical(branch, "relaxed_apoe_included")) "reverse_relaxed_exploratory_alternative_hb_outcome_sensitivity"
  else if (identical(branch, "relaxed_apoe_excluded")) "reverse_relaxed_exploratory_alternative_hb_outcome_sensitivity"
  else NA_character_
}
analysis_set_label <- function(branch) {
  sub("^strict_|^relaxed_", "", branch)
}

paths <- list(
  decision116 = rel("results", "qc", "chen_reverse_harmonised_mr_inputs_v1_freeze.json"),
  decision116_manifest = rel("results", "qc", "chen_reverse_harmonised_mr_inputs_v1_freeze_manifest.csv"),
  decision117 = rel("results", "qc", "chen_reverse_mr_analysis_contract_v1.json"),
  strict_contract = rel("results", "qc", "reverse_strict_primary_mr_analysis_contract_v1.json"),
  relaxed_contract = rel("results", "qc", "reverse_relaxed_mr_analysis_contract_v1.json"),
  strict_vuckovic_freeze = rel("results", "qc", "reverse_strict_primary_mr_v1_freeze.json"),
  relaxed_vuckovic_freeze = rel("results", "qc", "reverse_relaxed_mr_v1_freeze.json"),
  strict_inc = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_strict_apoe_included_v1.parquet"),
  strict_exc = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_strict_apoe_excluded_v1.parquet"),
  relaxed_inc = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_relaxed_apoe_included_v1.parquet"),
  relaxed_exc = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_relaxed_apoe_excluded_v1.parquet"),
  renv_lock = rel("renv.lock"),
  script = rel("R", "33_chen_reverse_mr_v1.R"),
  strict_table = rel("results", "tables", "chen_reverse_strict_mr_v1.csv"),
  relaxed_table = rel("results", "tables", "chen_reverse_relaxed_mr_estimates_v1.csv"),
  heterogeneity = rel("results", "tables", "chen_reverse_relaxed_heterogeneity_v1.csv"),
  egger = rel("results", "tables", "chen_reverse_relaxed_egger_intercept_v1.csv"),
  presso = rel("results", "tables", "chen_reverse_relaxed_mr_presso_v1.csv"),
  loo = rel("results", "tables", "chen_reverse_relaxed_leave_one_out_v1.csv"),
  single = rel("results", "tables", "chen_reverse_relaxed_single_snp_v1.csv"),
  comparison = rel("results", "tables", "chen_reverse_vuckovic_comparison_v1.csv"),
  qc = rel("results", "qc", "chen_reverse_mr_v1.json"),
  log = rel("results", "logs", "chen_reverse_mr_v1.log"),
  decision = rel("docs", "decisions", "118_chen_reverse_mr_v1_v1.1.md")
)

inputs <- unlist(paths[c(
  "decision116", "decision116_manifest", "decision117", "strict_contract",
  "relaxed_contract", "strict_vuckovic_freeze", "relaxed_vuckovic_freeze", "strict_inc", "strict_exc",
  "relaxed_inc", "relaxed_exc", "renv_lock"
)])
missing <- inputs[!file.exists(inputs)]
stop_if(length(missing) > 0L, paste("Required input missing:", paste(missing, collapse = "; ")))
targets <- unlist(paths[c(
  "strict_table", "relaxed_table", "heterogeneity", "egger", "presso",
  "loo", "single", "comparison", "qc", "log", "decision"
)])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Output or partial exists:", paste(occupied, collapse = "; ")))

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 118L), paste("Expected next decision 118, found ", next_decision, "; no outputs written."))

dir.create(dirname(paths$log), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = paths$log, append = TRUE)
log_line("stage=chen_reverse_mr_v1_start")

make_dat <- function(source, branch) {
  required <- c(
    "rsid", "exposure_beta", "exposure_se", "exposure_effect_allele",
    "exposure_other_allele", "outcome_beta_harmonised", "outcome_se_harmonised",
    "outcome_effect_allele_harmonised", "outcome_other_allele_harmonised",
    "final_valid", "F_stat"
  )
  missing_cols <- setdiff(required, names(source))
  stop_if(length(missing_cols) > 0L, paste("Required harmonised input columns absent for", branch, ":", paste(missing_cols, collapse = "; ")))
  if (!"outcome_eaf_harmonised" %in% names(source)) source$outcome_eaf_harmonised <- NA_real_
  if (!"exposure_eaf" %in% names(source)) source$exposure_eaf <- NA_real_
  stop_if(!all(source$final_valid), paste("Non-final-valid row found in MR input:", branch))
  stop_if(!identical(as.character(source$outcome_effect_allele_harmonised), as.character(source$exposure_effect_allele)) ||
            !identical(as.character(source$outcome_other_allele_harmonised), as.character(source$exposure_other_allele)),
          paste("Final alleles are not aligned for", branch))
  data.frame(
    SNP = as.character(source$rsid),
    beta.exposure = as.numeric(source$exposure_beta),
    se.exposure = as.numeric(source$exposure_se),
    effect_allele.exposure = as.character(source$exposure_effect_allele),
    other_allele.exposure = as.character(source$exposure_other_allele),
    eaf.exposure = as.numeric(source$exposure_eaf),
    beta.outcome = as.numeric(source$outcome_beta_harmonised),
    se.outcome = as.numeric(source$outcome_se_harmonised),
    effect_allele.outcome = as.character(source$outcome_effect_allele_harmonised),
    other_allele.outcome = as.character(source$outcome_other_allele_harmonised),
    eaf.outcome = as.numeric(source$outcome_eaf_harmonised),
    exposure = "FinnGen_R13_F5_DELIRIUM",
    outcome = "Chen_2020_Hb_BCX2",
    id.exposure = "finngen_R13_F5_DELIRIUM",
    id.outcome = "chen_2020_hb_bcx2",
    mr_keep = TRUE,
    stringsAsFactors = FALSE
  )
}

check_input <- function(source, dat, branch_info, branch) {
  expected_rsids <- split_rsids(branch_info$final_rsids)
  f <- (source$exposure_beta / source$exposure_se)^2
  list(
    row_count_matches_freeze = identical(nrow(source), as.integer(branch_info$final_valid_count)),
    rsids_match_freeze = identical(as.character(source$rsid), expected_rsids),
    all_final_valid = all(source$final_valid),
    rsid_unique = !anyDuplicated(source$rsid) && !anyNA(source$rsid),
    finite_exposure_effects = all(is.finite(dat$beta.exposure)),
    positive_exposure_se = all(is.finite(dat$se.exposure) & dat$se.exposure > 0),
    finite_outcome_effects = all(is.finite(dat$beta.outcome)),
    positive_outcome_se = all(is.finite(dat$se.outcome) & dat$se.outcome > 0),
    all_F_ge_10 = all(is.finite(f) & f >= 10),
    F_matches_freeze = all(abs(as.numeric(source$F_stat) - f) <= 1e-6),
    alleles_aligned = identical(as.character(source$outcome_effect_allele_harmonised), as.character(source$exposure_effect_allele)) &&
      identical(as.character(source$outcome_other_allele_harmonised), as.character(source$exposure_other_allele)),
    conversion_preserved = identical(as.character(source$rsid), dat$SNP) &&
      num_equal(source$exposure_beta, dat$beta.exposure) &&
      num_equal(source$exposure_se, dat$se.exposure) &&
      num_equal(source$outcome_beta_harmonised, dat$beta.outcome) &&
      num_equal(source$outcome_se_harmonised, dat$se.outcome) &&
      identical(as.character(source$exposure_effect_allele), dat$effect_allele.exposure) &&
      identical(as.character(source$exposure_other_allele), dat$other_allele.exposure) &&
      identical(as.character(source$outcome_effect_allele_harmonised), dat$effect_allele.outcome) &&
      identical(as.character(source$outcome_other_allele_harmonised), dat$other_allele.outcome),
    f_summary = data.frame(
      branch = branch, n = length(f), min = min(f), mean = mean(f),
      median = stats::median(f), max = max(f), F_lt10 = sum(f < 10),
      stringsAsFactors = FALSE
    )
  )
}

wald_one <- function(dat, branch, role, implementation, wald_sha) {
  stop_if(nrow(dat) != 1L, paste("Wald branch must have exactly one instrument:", branch))
  wr <- TwoSampleMR::mr_wald_ratio(
    b_exp = dat$beta.exposure[[1L]],
    b_out = dat$beta.outcome[[1L]],
    se_exp = dat$se.exposure[[1L]],
    se_out = dat$se.outcome[[1L]]
  )
  beta <- as.numeric(wr$b)
  se <- as.numeric(wr$se)
  pval <- as.numeric(wr$pval)
  ci_lower <- beta - 1.96 * se
  ci_upper <- beta + 1.96 * se
  data.frame(
    analysis_set = branch,
    analysis_role = role,
    rsid = dat$SNP[[1L]],
    nsnp = 1L,
    method = "Wald ratio",
    method_id = "mr_wald_ratio",
    implementation = implementation,
    function_body_sha256 = wald_sha,
    beta_exposure = dat$beta.exposure[[1L]],
    se_exposure = dat$se.exposure[[1L]],
    beta_outcome = dat$beta.outcome[[1L]],
    se_outcome = dat$se.outcome[[1L]],
    beta = beta,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    pval = pval,
    beta_doubling = beta * log(2),
    se_doubling = se * log(2),
    ci_lower_doubling = ci_lower * log(2),
    ci_upper_doubling = ci_upper * log(2),
    pval_doubling = pval,
    effect_scale = "standardized_quantitative_Hb_per_1_unit_genetically_predicted_log_odds_delirium",
    mr_estimable = TRUE,
    heterogeneity_status = "not_applicable_single_instrument",
    egger_intercept_status = "not_applicable_single_instrument",
    mr_presso_status = "not_applicable_single_instrument",
    loo_status = "not_applicable_single_instrument",
    single_snp_diagnostic_status = "not_separately_run_wald_is_formal_estimate",
    stringsAsFactors = FALSE
  )
}

method_lookup <- function() {
  m <- TwoSampleMR::mr_method_list()
  stats::setNames(as.character(m$obj), as.character(m$name))
}

run_mr_methods <- function(dat, branch, role, method_ids, method_roles) {
  x <- TwoSampleMR::mr(dat, method_list = method_ids)
  if (is.list(x) && "mr" %in% names(x)) x <- x$mr
  stop_if(!is.data.frame(x) || !all(c("method", "nsnp", "b", "se", "pval") %in% names(x)), "MR result schema is invalid.")
  lookup <- method_lookup()
  z <- data.frame(
    analysis_set = branch,
    analysis_role = role,
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
  z$beta_doubling <- z$beta * log(2)
  z$se_doubling <- z$se * log(2)
  z$ci_lower_doubling <- z$ci_lower * log(2)
  z$ci_upper_doubling <- z$ci_upper * log(2)
  z$pval_doubling <- z$pval
  z$effect_scale <- "standardized_quantitative_Hb_per_1_unit_genetically_predicted_log_odds_delirium"
  z$estimator_hierarchy_role <- unname(method_roles[z$method_id])
  z$nominal_p_lt_0_05 <- z$pval < 0.05
  z
}

run_heterogeneity <- function(dat, branch, role, contract) {
  ids <- vapply(contract$heterogeneity_plan, function(x) x$method_id, character(1))
  x <- TwoSampleMR::mr_heterogeneity(dat, method_list = ids)
  stop_if(!all(c("method", "Q", "Q_df", "Q_pval") %in% names(x)), "Heterogeneity result schema is invalid.")
  data.frame(
    analysis_set = branch,
    analysis_role = role,
    method = as.character(x$method),
    Q = as.numeric(x$Q),
    df = as.numeric(x$Q_df),
    pval = as.numeric(x$Q_pval),
    heterogeneity_p_lt_0_05 = as.numeric(x$Q_pval) < 0.05,
    automatic_snp_removal_allowed = FALSE,
    stringsAsFactors = FALSE
  )
}

run_egger <- function(dat, branch, role) {
  x <- TwoSampleMR::mr_pleiotropy_test(dat)
  stop_if(!all(c("egger_intercept", "se", "pval") %in% names(x)), "Egger intercept result schema is invalid.")
  pval <- as.numeric(x$pval[[1L]])
  data.frame(
    analysis_set = branch,
    analysis_role = role,
    instrument_count = nrow(dat),
    intercept = as.numeric(x$egger_intercept[[1L]]),
    se = as.numeric(x$se[[1L]]),
    pval = pval,
    nominal_p_lt_0_05 = pval < 0.05,
    interpretation_rule = if (pval < 0.05) {
      "statistical evidence consistent with directional pleiotropy"
    } else {
      "no statistical evidence of directional pleiotropy detected"
    },
    egger_precision_limitation = "limited_number_of_instruments",
    stringsAsFactors = FALSE
  )
}

run_presso <- function(dat, branch, role, config) {
  started <- Sys.time()
  set.seed(as.integer(config$seed))
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
      NbDistribution = as.integer(config$NbDistribution),
      SignifThreshold = as.numeric(config$SignifThreshold)
    ),
    error = function(e) e
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (inherits(result, "error")) {
    return(data.frame(
      analysis_set = branch, analysis_role = role, test_type = "MR-PRESSO",
      metric = "status", value = "not_estimable_or_failed", pval = NA_real_,
      outlier_count = NA_integer_, outlier_rsid = "", elapsed_seconds = elapsed,
      status = "not_estimable_or_failed", notes = conditionMessage(result),
      stringsAsFactors = FALSE
    ))
  }
  root_res <- result[["MR-PRESSO results"]]
  global <- if (!is.null(root_res)) root_res[["Global Test"]] else NULL
  outlier <- if (!is.null(root_res)) root_res[["Outlier Test"]] else NULL
  distortion <- if (!is.null(root_res)) root_res[["Distortion Test"]] else NULL
  rows <- list(data.frame(
    analysis_set = branch, analysis_role = role, test_type = "Global Test",
    metric = "RSSobs",
    value = if (!is.null(global) && "RSSobs" %in% names(global)) as.character(global$RSSobs[[1L]]) else "",
    pval = if (!is.null(global) && "Pvalue" %in% names(global)) suppressWarnings(as.numeric(global$Pvalue[[1L]])) else NA_real_,
    outlier_count = NA_integer_, outlier_rsid = "", elapsed_seconds = elapsed,
    status = "passed", notes = "mr_presso_status=passed",
    stringsAsFactors = FALSE
  ))
  ids <- if (!is.null(outlier) && is.data.frame(outlier) && nrow(outlier) > 0L) rownames(outlier) else character(0)
  if (length(ids) == 0L) {
    rows[[length(rows) + 1L]] <- data.frame(
      analysis_set = branch, analysis_role = role, test_type = "Outlier Test",
      metric = "outlier_count", value = "0", pval = NA_real_, outlier_count = 0L,
      outlier_rsid = "", elapsed_seconds = elapsed, status = "passed",
      notes = "no_outlier_reported", stringsAsFactors = FALSE
    )
  } else {
    rows[[length(rows) + 1L]] <- data.frame(
      analysis_set = branch, analysis_role = role, test_type = "Outlier Test",
      metric = "outlier", value = "", pval = NA_real_, outlier_count = length(ids),
      outlier_rsid = ids, elapsed_seconds = elapsed, status = "passed",
      notes = "outlier_reported_sensitivity_only", stringsAsFactors = FALSE
    )
  }
  if (!is.null(distortion)) {
    rows[[length(rows) + 1L]] <- data.frame(
      analysis_set = branch, analysis_role = role, test_type = "Distortion Test",
      metric = "raw_result", value = paste(capture.output(print(distortion)), collapse = " | "),
      pval = NA_real_, outlier_count = NA_integer_, outlier_rsid = "",
      elapsed_seconds = elapsed, status = "passed", notes = "distortion_output_if_available",
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

run_loo <- function(dat, branch, role, full_ivw_beta, full_ivw_pval) {
  raw <- TwoSampleMR::mr_leaveoneout(dat, method = TwoSampleMR::mr_ivw)
  if (!is.data.frame(raw)) raw <- as.data.frame(raw)
  stop_if(!all(c("SNP", "b", "se", "p") %in% names(raw)), "LOO result does not expose SNP/b/se/p fields.")
  labels <- as.character(raw$SNP)
  overall <- !is.na(labels) & grepl("^All", labels)
  stop_if(sum(overall) != 1L, "LOO requires exactly one All sentinel.")
  z0 <- raw[!overall, , drop = FALSE]
  stop_if(nrow(z0) != nrow(dat) || !setequal(as.character(z0$SNP), dat$SNP), "LOO SNP set does not match input.")
  z <- data.frame(
    analysis_set = branch, analysis_role = role, removed_rsid = as.character(z0$SNP),
    beta = as.numeric(z0$b), se = as.numeric(z0$se),
    ci_lower = as.numeric(z0$b) - 1.96 * as.numeric(z0$se),
    ci_upper = as.numeric(z0$b) + 1.96 * as.numeric(z0$se),
    pval = as.numeric(z0$p), full_ivw_beta = full_ivw_beta,
    stringsAsFactors = FALSE
  )
  z$beta_doubling <- z$beta * log(2)
  z$se_doubling <- z$se * log(2)
  z$ci_lower_doubling <- z$ci_lower * log(2)
  z$ci_upper_doubling <- z$ci_upper * log(2)
  z$absolute_shift <- abs(z$beta - full_ivw_beta)
  z$sign_change <- sign(z$beta) != sign(full_ivw_beta)
  z$nominal_significance_change <- (z$pval < 0.05) != (full_ivw_pval < 0.05)
  z
}

run_single <- function(dat, branch, role) {
  rows <- lapply(seq_len(nrow(dat)), function(i) {
    wr <- TwoSampleMR::mr_wald_ratio(
      b_exp = dat$beta.exposure[[i]], b_out = dat$beta.outcome[[i]],
      se_exp = dat$se.exposure[[i]], se_out = dat$se.outcome[[i]]
    )
    beta <- as.numeric(wr$b)
    se <- as.numeric(wr$se)
    ci_lower <- beta - 1.96 * se
    ci_upper <- beta + 1.96 * se
    data.frame(
      analysis_set = branch, analysis_role = role, rsid = dat$SNP[[i]],
      method = "Wald ratio", method_id = "mr_wald_ratio",
      beta = beta, se = se, ci_lower = ci_lower, ci_upper = ci_upper,
      pval = as.numeric(wr$pval), beta_doubling = beta * log(2),
      ci_lower_doubling = ci_lower * log(2), ci_upper_doubling = ci_upper * log(2),
      diagnostic_role = "single_snp_diagnostic_only",
      main_result_use_allowed = FALSE, stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

loo_summary <- function(loo) {
  do.call(rbind, lapply(split(loo, loo$analysis_set), function(x) {
    i <- which.max(x$absolute_shift)
    data.frame(
      analysis_set = x$analysis_set[[1L]],
      full_ivw_beta = x$full_ivw_beta[[1L]],
      max_absolute_shift = x$absolute_shift[[i]],
      max_shift_rsid = x$removed_rsid[[i]],
      any_sign_change = any(x$sign_change),
      any_nominal_significance_change = any(x$nominal_significance_change),
      stringsAsFactors = FALSE
    )
  }))
}
single_summary <- function(single) {
  do.call(rbind, lapply(split(single, single$analysis_set), function(x) {
    i <- which.max(abs(x$beta))
    data.frame(
      analysis_set = x$analysis_set[[1L]],
      max_abs_single_snp_beta = abs(x$beta[[i]]),
      max_abs_single_snp_rsid = x$rsid[[i]],
      nominal_p_lt_0_05_count = sum(x$pval < 0.05),
      stringsAsFactors = FALSE
    )
  }))
}

comparison_rows <- function(strict_results, relaxed_results, heterogeneity, egger, presso, loo_sum, single_sum, vuck_strict, vuck_relaxed, chen_rsids, vuck_rsids) {
  rows <- list()
  add <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
  si <- strict_results[strict_results$analysis_set == "strict_apoe_included", , drop = FALSE]
  se <- strict_results[strict_results$analysis_set == "strict_apoe_excluded", , drop = FALSE]
  add(
    comparison_domain = "strict_included", analysis_set = "included", method = "Wald ratio",
    vuckovic_value = vuck_strict$included_raw_result$beta, chen_value = si$beta,
    vuckovic_ci = paste(vuck_strict$included_raw_result$ci_lower, vuck_strict$included_raw_result$ci_upper, sep = ";"),
    chen_ci = paste(si$ci_lower, si$ci_upper, sep = ";"),
    vuckovic_pval = vuck_strict$included_raw_result$pval, chen_pval = si$pval,
    direction_consistent = sign(vuck_strict$included_raw_result$beta) == sign(si$beta),
    independent_replication = FALSE, formal_meta_analysis = FALSE, difference_test_performed = FALSE,
    notes = paste0("same_exposure_rsid=", identical(vuck_strict$included_rsid, si$rsid))
  )
  add(
    comparison_domain = "strict_excluded", analysis_set = "excluded", method = "Wald ratio",
    vuckovic_value = vuck_strict$excluded_raw_result$beta, chen_value = se$beta,
    vuckovic_ci = paste(vuck_strict$excluded_raw_result$ci_lower, vuck_strict$excluded_raw_result$ci_upper, sep = ";"),
    chen_ci = paste(se$ci_lower, se$ci_upper, sep = ";"),
    vuckovic_pval = vuck_strict$excluded_raw_result$pval, chen_pval = se$pval,
    direction_consistent = sign(vuck_strict$excluded_raw_result$beta) == sign(se$beta),
    independent_replication = FALSE, formal_meta_analysis = FALSE, difference_test_performed = FALSE,
    notes = paste0("same_exposure_rsid=", identical(vuck_strict$excluded_rsid, se$rsid))
  )
  for (set in c("relaxed_apoe_included", "relaxed_apoe_excluded")) {
    vset <- if (set == "relaxed_apoe_included") "APOE_included" else "APOE_excluded"
    vres <- if (set == "relaxed_apoe_included") vuck_relaxed$included_results else vuck_relaxed$excluded_results
    for (i in seq_along(vres)) {
      ch <- relaxed_results[relaxed_results$analysis_set == set & relaxed_results$method_id == vres[[i]]$method_id, , drop = FALSE]
      if (nrow(ch) != 1L) next
      add(
        comparison_domain = "relaxed_estimator", analysis_set = set, method = ch$method,
        vuckovic_value = vres[[i]]$beta, chen_value = ch$beta,
        vuckovic_ci = paste(vres[[i]]$ci_lower, vres[[i]]$ci_upper, sep = ";"),
        chen_ci = paste(ch$ci_lower, ch$ci_upper, sep = ";"),
        vuckovic_pval = vres[[i]]$pval, chen_pval = ch$pval,
        direction_consistent = sign(vres[[i]]$beta) == sign(ch$beta),
        independent_replication = FALSE, formal_meta_analysis = FALSE, difference_test_performed = FALSE,
        notes = paste0("vuckovic_set=", vset, "; hierarchy_role=", ch$estimator_hierarchy_role)
      )
    }
  }
  for (set in names(chen_rsids)) {
    v <- vuck_rsids[[set]]
    c <- chen_rsids[[set]]
    add(
      comparison_domain = "instrument_set_overlap", analysis_set = set, method = "rsid_set",
      vuckovic_value = length(v), chen_value = length(c), vuckovic_ci = "", chen_ci = "",
      vuckovic_pval = NA_real_, chen_pval = NA_real_, direction_consistent = NA,
      independent_replication = FALSE, formal_meta_analysis = FALSE, difference_test_performed = FALSE,
      notes = paste0(
        "same_final_rsid_set=", setequal(v, c), "; shared_n=", length(intersect(v, c)),
        "; vuckovic_only_n=", length(setdiff(v, c)), "; chen_only_n=", length(setdiff(c, v))
      )
    )
  }
  for (set in c("relaxed_apoe_included", "relaxed_apoe_excluded")) {
    h <- heterogeneity[heterogeneity$analysis_set == set, , drop = FALSE]
    e <- egger[egger$analysis_set == set, , drop = FALSE]
    p <- presso[presso$analysis_set == set, , drop = FALSE]
    l <- loo_sum[loo_sum$analysis_set == set, , drop = FALSE]
    s <- single_sum[single_sum$analysis_set == set, , drop = FALSE]
    add(
      comparison_domain = "diagnostic_pattern", analysis_set = set, method = "diagnostics",
      vuckovic_value = NA_real_, chen_value = NA_real_, vuckovic_ci = "", chen_ci = "",
      vuckovic_pval = NA_real_, chen_pval = NA_real_, direction_consistent = NA,
      independent_replication = FALSE, formal_meta_analysis = FALSE, difference_test_performed = FALSE,
      notes = paste0(
        "chen_heterogeneity_any_p_lt_0.05=", any(h$pval < 0.05, na.rm = TRUE),
        "; chen_egger_intercept_p=", paste(e$pval, collapse = ";"),
        "; chen_presso_status=", paste(unique(p$status), collapse = ";"),
        "; chen_loo_max_shift_rsid=", paste(l$max_shift_rsid, collapse = ";"),
        "; chen_single_max_abs_rsid=", paste(s$max_abs_single_snp_rsid, collapse = ";")
      )
    )
  }
  do.call(rbind, rows)
}

main <- function() {
  set.seed(2026L)
  renv_before <- hash_file(paths$renv_lock)
  freeze <- read_json(paths$decision116)
  contract <- read_json(paths$decision117)
  strict_contract <- read_json(paths$strict_contract)
  relaxed_contract <- read_json(paths$relaxed_contract)
  vuck_strict <- read_json(paths$strict_vuckovic_freeze)
  vuck_relaxed <- read_json(paths$relaxed_vuckovic_freeze)
  manifest_sha <- hash_file(paths$decision116_manifest)
  ts_desc <- utils::packageDescription("TwoSampleMR")
  mp_desc <- utils::packageDescription("MRPRESSO")
  ts_ver <- as.character(utils::packageVersion("TwoSampleMR"))
  ts_sha <- as.character(ts_desc[["RemoteSha"]])
  mp_ver <- as.character(utils::packageVersion("MRPRESSO"))
  mp_sha <- as.character(mp_desc[["RemoteSha"]])
  wald_fun <- get("mr_wald_ratio", envir = asNamespace("TwoSampleMR"))
  wald_body_sha <- hash_text(paste(deparse(body(wald_fun), width.cutoff = 500L), collapse = "\n"))
  git_present <- dir.exists(file.path(root, ".git"))
  git_status <- if (git_present) "present_not_queried" else "not_applicable_project_not_git_repository"

  strict_inc_info <- branch_from_freeze(freeze, "strict_apoe_included")
  strict_exc_info <- branch_from_freeze(freeze, "strict_apoe_excluded")
  relaxed_inc_info <- branch_from_freeze(freeze, "relaxed_apoe_included")
  relaxed_exc_info <- branch_from_freeze(freeze, "relaxed_apoe_excluded")
  method_ids <- method_ids_from_contract(contract)
  method_roles <- method_roles_from_contract(contract)

  pre_checks <- list(
    decision_116_input_freeze_gate = identical(freeze$freeze_status, "passed") &&
      isTRUE(freeze$approved_for_chen_reverse_mr_design) &&
      is_empty(freeze$hard_check_failures) &&
      identical(tolower(freeze$manifest_sha256), tolower(manifest_sha)),
    decision_117_contract_gate = identical(contract$contract_status, "frozen") &&
      isTRUE(contract$approved_for_chen_reverse_mr_execution) &&
      is_empty(contract$hard_check_failures) &&
      identical(contract$analysis_direction, "delirium_to_Hb") &&
      identical(contract$analysis_role, "reverse_alternative_hb_outcome_sensitivity"),
    strict_method_authority_gate = identical(strict_contract$wald_implementation$authoritative_method, "TwoSampleMR::mr_wald_ratio") &&
      identical(strict_contract$wald_implementation$function_body_sha256, wald_body_sha),
    relaxed_method_authority_gate = identical(method_ids, c("mr_ivw", "mr_weighted_median", "mr_egger_regression", "mr_weighted_mode", "mr_simple_mode")),
    software_versions_match_contract = identical(ts_ver, contract$software_environment$TwoSampleMR$version) &&
      identical(ts_sha, contract$software_environment$TwoSampleMR$RemoteSha) &&
      identical(mp_ver, contract$software_environment$MRPRESSO$version) &&
      identical(mp_sha, contract$software_environment$MRPRESSO$RemoteSha),
    vuckovic_relaxed_contract_gate = identical(relaxed_contract$contract_status, "frozen") &&
      isTRUE(relaxed_contract$approved_for_reverse_relaxed_mr_execution) &&
      is_empty(relaxed_contract$hard_check_failures)
  )
  stop_if(!all(unlist(pre_checks)), paste("Pre-execution gate failed:", paste(names(pre_checks)[!unlist(pre_checks)], collapse = "; ")))

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  strict_inc_src <- read_parquet(con, paths$strict_inc)
  strict_exc_src <- read_parquet(con, paths$strict_exc)
  relaxed_inc_src <- read_parquet(con, paths$relaxed_inc)
  relaxed_exc_src <- read_parquet(con, paths$relaxed_exc)
  strict_inc_dat <- make_dat(strict_inc_src, "strict_apoe_included")
  strict_exc_dat <- make_dat(strict_exc_src, "strict_apoe_excluded")
  relaxed_inc_dat <- make_dat(relaxed_inc_src, "relaxed_apoe_included")
  relaxed_exc_dat <- make_dat(relaxed_exc_src, "relaxed_apoe_excluded")

  input_audits <- list(
    strict_apoe_included = check_input(strict_inc_src, strict_inc_dat, strict_inc_info, "strict_apoe_included"),
    strict_apoe_excluded = check_input(strict_exc_src, strict_exc_dat, strict_exc_info, "strict_apoe_excluded"),
    relaxed_apoe_included = check_input(relaxed_inc_src, relaxed_inc_dat, relaxed_inc_info, "relaxed_apoe_included"),
    relaxed_apoe_excluded = check_input(relaxed_exc_src, relaxed_exc_dat, relaxed_exc_info, "relaxed_apoe_excluded")
  )
  input_flags <- unlist(lapply(input_audits, function(x) x[setdiff(names(x), "f_summary")]))
  stop_if(!all(input_flags), paste("Input identity gate failed:", paste(names(input_flags)[!input_flags], collapse = "; ")))

  strict_counts <- c(nrow(strict_inc_dat), nrow(strict_exc_dat))
  stop_if(any(strict_counts > 1L), "Strict n>=2 encountered; strict multi-IV design is not authorized.")
  strict_results <- rbind(
    wald_one(strict_inc_dat, "strict_apoe_included", branch_role("strict_apoe_included"), contract$strict_estimator_plan$authoritative_implementation, wald_body_sha),
    wald_one(strict_exc_dat, "strict_apoe_excluded", branch_role("strict_apoe_excluded"), contract$strict_estimator_plan$authoritative_implementation, wald_body_sha)
  )

  relaxed_inc_est <- run_mr_methods(relaxed_inc_dat, "relaxed_apoe_included", branch_role("relaxed_apoe_included"), method_ids, method_roles)
  relaxed_exc_est <- run_mr_methods(relaxed_exc_dat, "relaxed_apoe_excluded", branch_role("relaxed_apoe_excluded"), method_ids, method_roles)
  relaxed_results <- rbind(relaxed_inc_est, relaxed_exc_est)
  inc_ivw <- relaxed_inc_est[relaxed_inc_est$method_id == "mr_ivw", , drop = FALSE]
  exc_ivw <- relaxed_exc_est[relaxed_exc_est$method_id == "mr_ivw", , drop = FALSE]
  heterogeneity <- rbind(
    run_heterogeneity(relaxed_inc_dat, "relaxed_apoe_included", branch_role("relaxed_apoe_included"), contract),
    run_heterogeneity(relaxed_exc_dat, "relaxed_apoe_excluded", branch_role("relaxed_apoe_excluded"), contract)
  )
  egger <- rbind(
    run_egger(relaxed_inc_dat, "relaxed_apoe_included", branch_role("relaxed_apoe_included")),
    run_egger(relaxed_exc_dat, "relaxed_apoe_excluded", branch_role("relaxed_apoe_excluded"))
  )
  presso <- rbind(
    run_presso(relaxed_inc_dat, "relaxed_apoe_included", branch_role("relaxed_apoe_included"), contract$mr_presso_plan),
    run_presso(relaxed_exc_dat, "relaxed_apoe_excluded", branch_role("relaxed_apoe_excluded"), contract$mr_presso_plan)
  )
  loo <- rbind(
    run_loo(relaxed_inc_dat, "relaxed_apoe_included", branch_role("relaxed_apoe_included"), inc_ivw$beta[[1L]], inc_ivw$pval[[1L]]),
    run_loo(relaxed_exc_dat, "relaxed_apoe_excluded", branch_role("relaxed_apoe_excluded"), exc_ivw$beta[[1L]], exc_ivw$pval[[1L]])
  )
  single <- rbind(
    run_single(relaxed_inc_dat, "relaxed_apoe_included", branch_role("relaxed_apoe_included")),
    run_single(relaxed_exc_dat, "relaxed_apoe_excluded", branch_role("relaxed_apoe_excluded"))
  )
  loo_sum <- loo_summary(loo)
  single_sum <- single_summary(single)

  chen_rsids <- list(
    strict_included = strict_inc_dat$SNP,
    strict_excluded = strict_exc_dat$SNP,
    relaxed_included = relaxed_inc_dat$SNP,
    relaxed_excluded = relaxed_exc_dat$SNP
  )
  vuck_rsids <- list(
    strict_included = vuck_strict$included_rsid,
    strict_excluded = vuck_strict$excluded_rsid,
    relaxed_included = relaxed_contract$analysis_sets[[1L]]$final_rsids,
    relaxed_excluded = relaxed_contract$analysis_sets[[2L]]$final_rsids
  )

  comparison <- comparison_rows(
    strict_results, relaxed_results, heterogeneity, egger, presso, loo_sum, single_sum,
    vuck_strict, vuck_relaxed, chen_rsids, vuck_rsids
  )

  renv_after <- hash_file(paths$renv_lock)
  hard_checks <- c(pre_checks, list(
    strict_branch_counts_match_freeze = nrow(strict_inc_dat) == strict_inc_info$final_valid_count &&
      nrow(strict_exc_dat) == strict_exc_info$final_valid_count,
    relaxed_branch_counts_match_freeze = nrow(relaxed_inc_dat) == relaxed_inc_info$final_valid_count &&
      nrow(relaxed_exc_dat) == relaxed_exc_info$final_valid_count,
    all_branch_rsids_match_freeze = all(vapply(input_audits, function(x) x$rsids_match_freeze, logical(1))),
    all_final_alleles_aligned = all(vapply(input_audits, function(x) x$alleles_aligned, logical(1))),
    no_reharmonisation = TRUE,
    strict_n1_wald_execution_correct = nrow(strict_results) == 2L && all(strict_results$nsnp == 1L) && all(strict_results$method_id == "mr_wald_ratio"),
    strict_no_multi_iv_diagnostics = all(strict_results$heterogeneity_status == "not_applicable_single_instrument") &&
      all(strict_results$single_snp_diagnostic_status == "not_separately_run_wald_is_formal_estimate"),
    relaxed_ivw_run = sum(relaxed_results$method_id == "mr_ivw") == 2L,
    relaxed_egger_run = sum(relaxed_results$method_id == "mr_egger_regression") == 2L,
    relaxed_weighted_median_run = sum(relaxed_results$method_id == "mr_weighted_median") == 2L,
    relaxed_weighted_mode_run = sum(relaxed_results$method_id == "mr_weighted_mode") == 2L,
    relaxed_simple_mode_run = sum(relaxed_results$method_id == "mr_simple_mode") == 2L,
    relaxed_method_hierarchy_matches_contract = setequal(relaxed_results$method_id, method_ids),
    heterogeneity_completed = nrow(heterogeneity) == 4L,
    egger_intercept_completed = nrow(egger) == 2L,
    mr_presso_attempted = all(c("relaxed_apoe_included", "relaxed_apoe_excluded") %in% unique(presso$analysis_set)),
    loo_completed = nrow(loo) == nrow(relaxed_inc_dat) + nrow(relaxed_exc_dat),
    relaxed_single_snp_completed = nrow(single) == nrow(relaxed_inc_dat) + nrow(relaxed_exc_dat),
    doubling_odds_transform_correct = all(abs(strict_results$beta_doubling - strict_results$beta * log(2)) <= 1e-12) &&
      all(abs(relaxed_results$beta_doubling - relaxed_results$beta * log(2)) <= 1e-12) &&
      identical(as.numeric(strict_results$pval_doubling), as.numeric(strict_results$pval)) &&
      identical(as.numeric(relaxed_results$pval_doubling), as.numeric(relaxed_results$pval)),
    no_or_transform = !any(grepl("^OR$|OR_", c(names(strict_results), names(relaxed_results)))),
    vuckovic_comparison_noninferential = all(comparison$independent_replication == FALSE) &&
      all(comparison$formal_meta_analysis == FALSE) &&
      all(comparison$difference_test_performed == FALSE),
    independent_replication_false = identical(contract$independent_replication, FALSE),
    strict_relaxed_hierarchy_preserved = identical(contract$evidence_hierarchy$strict_primary_superseded_by_relaxed, FALSE) &&
      identical(contract$evidence_hierarchy$relaxed_confirmatory, FALSE),
    no_posthoc_filtering = TRUE,
    no_proxy = TRUE,
    no_liftover = TRUE,
    no_steiger = TRUE,
    renv_lock_unchanged = identical(renv_before, renv_after),
    git_status_not_required = !git_present
  ))
  failures <- names(hard_checks)[!unlist(hard_checks)]
  status <- if (length(failures) == 0L) "passed" else "failed"
  approved <- identical(status, "passed")

  qc <- list(
    mr_version = "v1",
    decision = 118,
    date = "2026-08-13",
    analysis_direction = "delirium_to_Hb",
    analysis_role = "reverse_alternative_hb_outcome_sensitivity",
    independent_replication = FALSE,
    source_mr_input_freeze_decision = 116,
    mr_contract_decision = 117,
    strict_threshold = 5e-8,
    relaxed_threshold = 5e-6,
    strict_results = records(strict_results),
    relaxed_results = records(relaxed_results),
    heterogeneity_results = records(heterogeneity),
    egger_intercept_results = records(egger),
    mr_presso_results = records(presso),
    leave_one_out_summary = records(loo_sum),
    single_snp_summary = records(single_sum),
    vuckovic_comparison = records(comparison),
    strict_primary_superseded_by_relaxed = FALSE,
    relaxed_confirmatory = FALSE,
    steiger_run = FALSE,
    steiger_status = "deferred_to_unified_directionality_sensitivity_stage",
    git_repository_present = git_present,
    git_status = git_status,
    git_status_hard_failure = FALSE,
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
      decision116_manifest = manifest_sha,
      decision117_contract = hash_file(paths$decision117),
      vuckovic_relaxed_contract = hash_file(paths$relaxed_contract),
      strict_apoe_included = hash_file(paths$strict_inc),
      strict_apoe_excluded = hash_file(paths$strict_exc),
      relaxed_apoe_included = hash_file(paths$relaxed_inc),
      relaxed_apoe_excluded = hash_file(paths$relaxed_exc),
      vuckovic_strict_freeze = hash_file(paths$strict_vuckovic_freeze),
      vuckovic_relaxed_freeze = hash_file(paths$relaxed_vuckovic_freeze)
    ),
    instrument_strength_summary = lapply(input_audits, function(x) as.list(x$f_summary[1, , drop = FALSE])),
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after,
    renv_lock_unchanged = identical(renv_before, renv_after),
    mr_status = status,
    approved_for_chen_reverse_results_interpretation = approved,
    hard_checks = hard_checks,
    hard_check_failures = failures,
    informational_findings = list(
      confirmatory_status = "Chen reverse is alternative-Hb-outcome sensitivity; relaxed branch is exploratory_not_confirmatory",
      decision_114_status_preserved = "failed_due_readback_consistency_technical_issue",
      decision_115_authority_preserved = "technical_readback_recovery_not_scientific_rerun",
      mr_presso_failure_handling = "if failed, exact error is reported and no fake Global P or no-outlier statement is generated",
      vuckovic_chen_independent_replication = FALSE,
      freeze_created = FALSE
    )
  )

  decision_lines <- c(
    "# Decision 118: Chen Reverse MR V1",
    "",
    "Date: 2026-08-13",
    "",
    "## Status",
    sprintf("mr_status: `%s`", status),
    sprintf("approved_for_chen_reverse_results_interpretation: `%s`", if (approved) "TRUE" else "FALSE"),
    "",
    "## Decision",
    "Execute Chen reverse alternative-Hb-outcome sensitivity MR V1 according to Decision 116 MR-input freeze and Decision 117 MR analysis contract.",
    "",
    "This run did not create a Chen Reverse MR Results Freeze, did not run Steiger, did not re-harmonise, did not reclump, did not extract outcomes, did not use proxies/liftOver, and did not apply post-hoc SNP filtering.",
    "",
    "## Authority Gates",
    sprintf("- Decision 116 manifest SHA-256: `%s`.", manifest_sha),
    "- Decision 116 freeze gate: passed.",
    "- Decision 117 contract gate: frozen and approved.",
    "",
    "## Strict Results",
    paste0("- `", strict_results$analysis_set, "` `", strict_results$rsid, "`: beta=`", signif(strict_results$beta, 6), "`, SE=`", signif(strict_results$se, 6), "`, 95% CI `", signif(strict_results$ci_lower, 6), " to ", signif(strict_results$ci_upper, 6), "`, P=`", signif(strict_results$pval, 6), "`; doubling beta=`", signif(strict_results$beta_doubling, 6), "`, doubling CI `", signif(strict_results$ci_lower_doubling, 6), " to ", signif(strict_results$ci_upper_doubling, 6), "`."),
    "",
    "## Relaxed IVW Results",
    paste0("- `", relaxed_results$analysis_set[relaxed_results$method_id == "mr_ivw"], "`: beta=`", signif(relaxed_results$beta[relaxed_results$method_id == "mr_ivw"], 6), "`, SE=`", signif(relaxed_results$se[relaxed_results$method_id == "mr_ivw"], 6), "`, 95% CI `", signif(relaxed_results$ci_lower[relaxed_results$method_id == "mr_ivw"], 6), " to ", signif(relaxed_results$ci_upper[relaxed_results$method_id == "mr_ivw"], 6), "`, P=`", signif(relaxed_results$pval[relaxed_results$method_id == "mr_ivw"], 6), "`; doubling beta=`", signif(relaxed_results$beta_doubling[relaxed_results$method_id == "mr_ivw"], 6), "`."),
    "",
    "## Diagnostics",
    sprintf("- Heterogeneity rows: `%s`.", nrow(heterogeneity)),
    sprintf("- Egger intercept rows: `%s`.", nrow(egger)),
    sprintf("- MR-PRESSO statuses: `%s`.", paste(unique(presso$status), collapse = ";")),
    sprintf("- LOO summary rows: `%s`.", nrow(loo_sum)),
    sprintf("- Single-SNP summary rows: `%s`.", nrow(single_sum)),
    "",
    "## Vuckovic Comparison",
    "Vuckovic vs Chen comparison is noninferential robustness assessment only: `independent_replication=FALSE`, `formal_meta_analysis=FALSE`, `difference_test_performed=FALSE`.",
    "",
    "## Software And Audit",
    sprintf("- R: `%s`.", R.version.string),
    sprintf("- TwoSampleMR: `%s`, RemoteSha `%s`.", ts_ver, ts_sha),
    sprintf("- MRPRESSO: `%s`, RemoteSha `%s`.", mp_ver, mp_sha),
    "- Seed: `2026`.",
    sprintf("- renv.lock SHA before/after: `%s` / `%s`.", renv_before, renv_after),
    sprintf("- git repository present: `%s`; git status: `%s`.", git_present, git_status),
    "- Steiger: `deferred_to_unified_directionality_sensitivity_stage`.",
    "",
    "## Hard Check Failures",
    if (length(failures) == 0L) "- none" else paste0("- `", failures, "`"),
    "",
    "## Outputs Created",
    "- `R/33_chen_reverse_mr_v1.R`",
    "- `results/tables/chen_reverse_strict_mr_v1.csv`",
    "- `results/tables/chen_reverse_relaxed_mr_estimates_v1.csv`",
    "- `results/tables/chen_reverse_relaxed_heterogeneity_v1.csv`",
    "- `results/tables/chen_reverse_relaxed_egger_intercept_v1.csv`",
    "- `results/tables/chen_reverse_relaxed_mr_presso_v1.csv`",
    "- `results/tables/chen_reverse_relaxed_leave_one_out_v1.csv`",
    "- `results/tables/chen_reverse_relaxed_single_snp_v1.csv`",
    "- `results/tables/chen_reverse_vuckovic_comparison_v1.csv`",
    "- `results/qc/chen_reverse_mr_v1.json`",
    "- `results/logs/chen_reverse_mr_v1.log`",
    "- `docs/decisions/118_chen_reverse_mr_v1_v1.1.md`",
    "",
    "## Next Gate",
    "Stop here for human review. Do not create Chen Reverse MR Results Freeze until separately approved."
  )

  log_tail <- c(
    sprintf("decision116_manifest_sha=%s", manifest_sha),
    sprintf("strict_apoe_included=%s", paste(signif(unlist(strict_results[strict_results$analysis_set == "strict_apoe_included", c("beta", "se", "ci_lower", "ci_upper", "pval")]), 8), collapse = ";")),
    sprintf("strict_apoe_excluded=%s", paste(signif(unlist(strict_results[strict_results$analysis_set == "strict_apoe_excluded", c("beta", "se", "ci_lower", "ci_upper", "pval")]), 8), collapse = ";")),
    sprintf("relaxed_ivw_included=%s", paste(signif(unlist(inc_ivw[, c("beta", "se", "ci_lower", "ci_upper", "pval")]), 8), collapse = ";")),
    sprintf("relaxed_ivw_excluded=%s", paste(signif(unlist(exc_ivw[, c("beta", "se", "ci_lower", "ci_upper", "pval")]), 8), collapse = ";")),
    sprintf("mr_status=%s", status),
    sprintf("approved_for_chen_reverse_results_interpretation=%s", approved),
    sprintf("hard_check_failures=%s", if (length(failures) == 0L) "none" else paste(failures, collapse = ";")),
    "freeze_created=FALSE",
    "steiger_run=FALSE"
  )
  for (line in log_tail) log_line(line)

  write_csv_precise(strict_results, paths$strict_table)
  write_csv_precise(relaxed_results, paths$relaxed_table)
  write_csv_precise(heterogeneity, paths$heterogeneity)
  write_csv_precise(egger, paths$egger)
  write_csv_precise(presso, paths$presso)
  write_csv_precise(loo, paths$loo)
  write_csv_precise(single, paths$single)
  write_csv_precise(comparison, paths$comparison)
  write_json(qc, paths$qc)
  write_text(decision_lines, paths$decision)

  if (!approved) stop("Chen reverse MR V1 hard checks failed; outputs retained for audit.", call. = FALSE)
  cat("Decision 118 Chen Reverse MR V1 completed\n")
  cat("mr_status=", status, "\n", sep = "")
  cat("approved_for_chen_reverse_results_interpretation=", approved, "\n", sep = "")
  cat("hard_check_failures=", if (length(failures) == 0L) "none" else paste(failures, collapse = ";"), "\n", sep = "")
}

tryCatch(main(), error = function(e) {
  log_line("status=failed error=", conditionMessage(e))
  quit(status = 1L)
})
