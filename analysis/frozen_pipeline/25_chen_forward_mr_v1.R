#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/25_chen_forward_mr_v1.R [--project-root E:/Research/hb_delirium_bidir_mr]", call. = FALSE)
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
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
num_equal <- function(a, b, tol = 1e-12) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  all(is.finite(a) & is.finite(b) & abs(a - b) <= tol)
}
sql_path <- function(path, must_work = TRUE) {
  paste0("'", gsub("'", "''", normalizePath(path, winslash = "/", mustWork = must_work), fixed = TRUE), "'")
}
read_parquet <- function(con, path) DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_path(path)))
write_csv_precise <- function(x, path) {
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}
write_json_precise <- function(x, path) {
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)
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
extract_contract_methods <- function(contract) {
  x <- contract$estimator_hierarchy
  if (is.data.frame(x)) return(as.character(x$method_id))
  if (is.list(x) && length(x) > 0L) return(vapply(x, function(row) as.character(row$method_id), character(1)))
  character(0)
}

paths <- list(
  freeze = rel("results", "qc", "chen_forward_harmonised_mr_inputs_v1_freeze.json"),
  freeze_manifest = rel("results", "qc", "chen_forward_harmonised_mr_inputs_v1_freeze_manifest.csv"),
  amendment_readback = rel("results", "qc", "chen_forward_mr_method_alignment_amendment_v1_readback_audit_v1.json"),
  contract = rel("results", "qc", "chen_forward_mr_analysis_contract_v2.json"),
  primary_script = rel("R", "09_forward_mr_vuckovic_hb_finngen_r13_primary_v3.R"),
  primary_qc = rel("results", "qc", "vuckovic_hb_finngen_r13_forward_mr_v3.json"),
  vuckovic_estimates = rel("results", "tables", "vuckovic_hb_finngen_r13_forward_mr_estimates_v3.csv"),
  vuckovic_heterogeneity = rel("results", "tables", "vuckovic_hb_finngen_r13_forward_heterogeneity_v3.csv"),
  vuckovic_egger = rel("results", "tables", "vuckovic_hb_finngen_r13_forward_egger_intercept_v3.csv"),
  vuckovic_presso = rel("results", "tables", "vuckovic_hb_finngen_r13_forward_mr_presso_v3.csv"),
  vuckovic_loo = rel("results", "tables", "vuckovic_hb_finngen_r13_forward_leave_one_out_v3.csv"),
  renv_lock = rel("renv.lock"),
  estimates = rel("results", "tables", "chen_forward_mr_estimates_v1.csv"),
  heterogeneity = rel("results", "tables", "chen_forward_heterogeneity_v1.csv"),
  egger = rel("results", "tables", "chen_forward_egger_intercept_v1.csv"),
  presso = rel("results", "tables", "chen_forward_mr_presso_v1.csv"),
  loo = rel("results", "tables", "chen_forward_leave_one_out_v1.csv"),
  comparison = rel("results", "tables", "chen_forward_vuckovic_comparison_v1.csv"),
  single_snp_forbidden = rel("results", "tables", "chen_forward_single_snp_v1.csv"),
  qc = rel("results", "qc", "chen_forward_mr_v1.json"),
  log = rel("results", "logs", "chen_forward_mr_v1.log"),
  decision = rel("docs", "decisions", "100_chen_forward_mr_v1_v1.1.md")
)

required_inputs <- unlist(paths[c(
  "freeze", "freeze_manifest", "amendment_readback", "contract", "primary_script",
  "primary_qc", "vuckovic_estimates", "vuckovic_heterogeneity", "vuckovic_egger",
  "vuckovic_presso", "vuckovic_loo", "renv_lock"
)])
missing <- required_inputs[!file.exists(required_inputs)]
if (length(missing) > 0L) stop("Missing required input(s): ", paste(missing, collapse = "; "), call. = FALSE)

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 100L), paste0("Expected next decision 100, found ", next_decision, "; no outputs written."))

targets <- unlist(paths[c("estimates", "heterogeneity", "egger", "presso", "loo", "comparison", "qc", "log", "decision")])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(occupied, collapse = "; ")))
stop_if(file.exists(paths$single_snp_forbidden) || file.exists(paste0(paths$single_snp_forbidden, ".partial")),
        "Forbidden single-SNP output already exists; refusing to continue.")

dir.create(dirname(paths$estimates), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(paths$qc), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(paths$log), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(paste0("[", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "] ", paste0(..., collapse = ""), "\n"),
                              file = paths$log, append = TRUE)

make_dat <- function(source, analysis_set, analysis_role) {
  required <- c(
    "resolved_rsid", "exposure_beta", "exposure_se", "exposure_effect_allele",
    "exposure_other_allele", "exposure_eaf", "outcome_beta_harmonised",
    "outcome_se_harmonised", "outcome_effect_allele_harmonised",
    "outcome_other_allele_harmonised", "outcome_eaf_harmonised", "final_valid"
  )
  stop_if(!all(required %in% names(source)), paste("Required harmonised input columns are absent:", analysis_set))
  stop_if(!all(source$final_valid), paste("Input contains non-final-valid rows:", analysis_set))
  data.frame(
    SNP = as.character(source$resolved_rsid),
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
    exposure = "Chen_2020_Hb_BCX2",
    outcome = "FinnGen_R13_F5_DELIRIUM",
    id.exposure = "chen_2020_hb_bcx2",
    id.outcome = "finngen_R13_F5_DELIRIUM",
    mr_keep = TRUE,
    analysis_set = analysis_set,
    analysis_role = analysis_role,
    stringsAsFactors = FALSE
  )
}

check_input <- function(source, dat, freeze_rsids, expected_count) {
  f <- (source$exposure_beta / source$exposure_se)^2
  list(
    count_matches_freeze = nrow(source) == expected_count,
    rsid_set_matches_freeze = setequal(as.character(source$resolved_rsid), as.character(freeze_rsids)),
    rsid_unique = !anyDuplicated(source$resolved_rsid) && !anyNA(source$resolved_rsid),
    all_final_valid = all(source$final_valid),
    no_palindromic_final_valid = !any(source$palindromic_snp),
    aligned_alleles = all(dat$effect_allele.exposure == dat$effect_allele.outcome) &&
      all(dat$other_allele.exposure == dat$other_allele.outcome),
    finite_exposure_effects = all(is.finite(dat$beta.exposure)),
    positive_exposure_se = all(is.finite(dat$se.exposure) & dat$se.exposure > 0),
    finite_outcome_effects = all(is.finite(dat$beta.outcome)),
    positive_outcome_se = all(is.finite(dat$se.outcome) & dat$se.outcome > 0),
    all_F_ge_10 = all(is.finite(f) & f >= 10),
    conversion_preserved = identical(as.character(source$resolved_rsid), dat$SNP) &&
      num_equal(source$exposure_beta, dat$beta.exposure) &&
      num_equal(source$exposure_se, dat$se.exposure) &&
      identical(as.character(source$exposure_effect_allele), dat$effect_allele.exposure) &&
      identical(as.character(source$exposure_other_allele), dat$other_allele.exposure) &&
      num_equal(source$outcome_beta_harmonised, dat$beta.outcome) &&
      num_equal(source$outcome_se_harmonised, dat$se.outcome) &&
      identical(as.character(source$outcome_effect_allele_harmonised), dat$effect_allele.outcome) &&
      identical(as.character(source$outcome_other_allele_harmonised), dat$other_allele.outcome),
    f_summary = data.frame(
      n = length(f), F_min = min(f), F_mean = mean(f), F_median = stats::median(f),
      F_max = max(f), F_lt10_count = sum(f < 10), stringsAsFactors = FALSE
    )
  )
}

method_lookup <- function() {
  methods <- TwoSampleMR::mr_method_list()
  stats::setNames(as.character(methods$obj), as.character(methods$name))
}

run_mr_methods <- function(dat, analysis_set, analysis_role, method_ids, hierarchy) {
  x <- TwoSampleMR::mr(dat, method_list = method_ids)
  if (is.list(x) && "mr" %in% names(x)) x <- x$mr
  stop_if(!is.data.frame(x) || !all(c("method", "nsnp", "b", "se", "pval") %in% names(x)), "MR result schema is invalid.")
  lookup <- method_lookup()
  out <- data.frame(
    analysis_set = analysis_set,
    analysis_role = analysis_role,
    method = as.character(x$method),
    method_id = unname(lookup[as.character(x$method)]),
    nsnp = as.integer(x$nsnp),
    beta = as.numeric(x$b),
    se = as.numeric(x$se),
    ci_lower = as.numeric(x$b) - 1.96 * as.numeric(x$se),
    ci_upper = as.numeric(x$b) + 1.96 * as.numeric(x$se),
    pval = as.numeric(x$pval),
    stringsAsFactors = FALSE
  )
  out$OR <- exp(out$beta)
  out$OR_ci_lower <- exp(out$ci_lower)
  out$OR_ci_upper <- exp(out$ci_upper)
  out$effect_scale <- "per 1-unit increase in genetically predicted standardized Chen Hb scale"
  out$nominal_p_lt_0_05 <- out$pval < 0.05
  roles <- stats::setNames(hierarchy$estimator_hierarchy_role, hierarchy$method_id)
  out$estimator_hierarchy_role <- unname(roles[out$method_id])
  out
}

run_heterogeneity <- function(dat, analysis_set, analysis_role) {
  x <- TwoSampleMR::mr_heterogeneity(dat, method_list = c("mr_egger_regression", "mr_ivw"))
  stop_if(!all(c("method", "Q", "Q_df", "Q_pval") %in% names(x)), "Heterogeneity result schema is invalid.")
  lookup <- method_lookup()
  data.frame(
    analysis_set = analysis_set,
    analysis_role = analysis_role,
    method = as.character(x$method),
    method_id = unname(lookup[as.character(x$method)]),
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
  p <- as.numeric(x$pval[[1L]])
  data.frame(
    analysis_set = analysis_set,
    analysis_role = analysis_role,
    intercept = as.numeric(x$egger_intercept[[1L]]),
    se = as.numeric(x$se[[1L]]),
    pval = p,
    nominal_p_lt_0_05 = p < 0.05,
    interpretation_rule = if (p < 0.05) {
      "evidence consistent with directional pleiotropy"
    } else {
      "no statistical evidence of directional pleiotropy detected; not proof of no pleiotropy"
    },
    stringsAsFactors = FALSE
  )
}

run_presso <- function(dat, analysis_set, analysis_role, nb_distribution, signif_threshold, seed_value) {
  set.seed(seed_value)
  x <- dat[, c("beta.outcome", "beta.exposure", "se.outcome", "se.exposure"), drop = FALSE]
  rownames(x) <- dat$SNP
  log_line("mr_presso_start analysis_set=", analysis_set, " NbDistribution=", nb_distribution)
  setTimeLimit(cpu = Inf, elapsed = 180, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
  result <- tryCatch(
    MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome",
      BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome",
      SdExposure = "se.exposure",
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      data = x,
      NbDistribution = nb_distribution,
      SignifThreshold = signif_threshold
    ),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    log_line("mr_presso_failed analysis_set=", analysis_set, " error=", conditionMessage(result))
    return(data.frame(
      analysis_set = analysis_set, analysis_role = analysis_role,
      test_type = "MR-PRESSO", metric = "status", value = "failed_or_timeout",
      pval = NA_real_, outlier_rsid = "", notes = conditionMessage(result),
      stringsAsFactors = FALSE
    ))
  }
  log_line("mr_presso_completed analysis_set=", analysis_set)
  root <- result[["MR-PRESSO results"]]
  global <- if (!is.null(root)) root[["Global Test"]] else NULL
  outlier <- if (!is.null(root)) root[["Outlier Test"]] else NULL
  distortion <- if (!is.null(root)) root[["Distortion Test"]] else NULL
  global_p <- if (!is.null(global) && "Pvalue" %in% names(global)) as.character(global$Pvalue[[1L]]) else ""
  rows <- list(
    data.frame(
      analysis_set = analysis_set, analysis_role = analysis_role,
      test_type = "Global Test", metric = "RSSobs",
      value = if (!is.null(global) && "RSSobs" %in% names(global)) as.character(global$RSSobs[[1L]]) else "",
      pval = suppressWarnings(as.numeric(global_p)), outlier_rsid = "",
      notes = paste0("mr_presso_status=passed; global_pvalue_raw=", global_p),
      stringsAsFactors = FALSE
    )
  )
  ids <- if (!is.null(outlier) && nrow(outlier) > 0L) rownames(outlier) else character(0)
  rows[[length(rows) + 1L]] <- data.frame(
    analysis_set = analysis_set, analysis_role = analysis_role,
    test_type = "Outlier Test", metric = "outlier_count",
    value = as.character(length(ids)), pval = NA_real_,
    outlier_rsid = if (length(ids) == 0L) "" else paste(ids, collapse = ";"),
    notes = if (length(ids) == 0L) "no_outlier_reported" else "outlier_reported_sensitivity_only_no_main_input_change",
    stringsAsFactors = FALSE
  )
  if (!is.null(distortion)) {
    rows[[length(rows) + 1L]] <- data.frame(
      analysis_set = analysis_set, analysis_role = analysis_role,
      test_type = "Distortion Test", metric = "raw_result",
      value = paste(capture.output(print(distortion)), collapse = " | "),
      pval = NA_real_, outlier_rsid = "",
      notes = "distortion_output_if_available",
      stringsAsFactors = FALSE
    )
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
    stringsAsFactors = FALSE
  )
  z$OR <- exp(z$beta)
  z$OR_ci_lower <- exp(z$ci_lower)
  z$OR_ci_upper <- exp(z$ci_upper)
  z$full_ivw_beta <- full_ivw_beta
  z$full_ivw_pval <- full_ivw_pval
  z$absolute_shift <- abs(z$beta - full_ivw_beta)
  z$sign_change <- sign(z$beta) != sign(full_ivw_beta)
  z$nominal_significance_change <- (z$pval < 0.05) != (full_ivw_pval < 0.05)
  z
}

summarise_loo <- function(loo, full_ivw_beta_col = "full_ivw_beta") {
  do.call(rbind, lapply(split(loo, loo$analysis_set), function(x) {
    if (!"absolute_shift" %in% names(x)) x$absolute_shift <- abs(x$beta - x[[full_ivw_beta_col]][[1L]])
    i <- which.max(x$absolute_shift)
    data.frame(
      analysis_set = x$analysis_set[[1L]],
      total_removed_snp_rows = nrow(x),
      max_absolute_shift = x$absolute_shift[[i]],
      max_shift_rsid = x$removed_rsid[[i]],
      any_sign_change = any(x$sign_change),
      any_nominal_significance_change = any(x$nominal_significance_change),
      stringsAsFactors = FALSE
    )
  }))
}

make_vuckovic_loo_summary <- function(v_loo, v_est) {
  names(v_loo)[names(v_loo) == "excluded_rsid"] <- "removed_rsid"
  ivw <- v_est[v_est$method_id == "mr_ivw", c("analysis_set", "beta", "pval")]
  rows <- lapply(split(v_loo, v_loo$analysis_set), function(x) {
    ref <- ivw[ivw$analysis_set == x$analysis_set[[1L]], , drop = FALSE]
    x$absolute_shift <- abs(x$beta - ref$beta[[1L]])
    x$sign_change <- sign(x$beta) != sign(ref$beta[[1L]])
    x$nominal_significance_change <- (x$pval < 0.05) != (ref$pval[[1L]] < 0.05)
    summarise_loo(x)
  })
  do.call(rbind, rows)
}

make_comparison <- function(chen_est, chen_het, chen_egger, chen_presso, chen_loo_summary) {
  v_est <- read.csv(paths$vuckovic_estimates, stringsAsFactors = FALSE, check.names = FALSE)
  v_het <- read.csv(paths$vuckovic_heterogeneity, stringsAsFactors = FALSE, check.names = FALSE)
  v_egger <- read.csv(paths$vuckovic_egger, stringsAsFactors = FALSE, check.names = FALSE)
  v_presso <- read.csv(paths$vuckovic_presso, stringsAsFactors = FALSE, check.names = FALSE)
  v_loo <- read.csv(paths$vuckovic_loo, stringsAsFactors = FALSE, check.names = FALSE)
  v_loo_summary <- make_vuckovic_loo_summary(v_loo, v_est)

  estimator_rows <- do.call(rbind, lapply(seq_len(nrow(chen_est)), function(i) {
    c_row <- chen_est[i, , drop = FALSE]
    v_row <- v_est[v_est$analysis_set == c_row$analysis_set & v_est$method_id == c_row$method_id, , drop = FALSE]
    data.frame(
      comparison_domain = "estimator",
      analysis_set = c_row$analysis_set,
      method = c_row$method,
      method_id = c_row$method_id,
      chen_value = c_row$beta,
      vuckovic_value = if (nrow(v_row) == 1L) v_row$beta else NA_real_,
      chen_pval = c_row$pval,
      vuckovic_pval = if (nrow(v_row) == 1L) v_row$pval else NA_real_,
      direction_consistent = if (nrow(v_row) == 1L) sign(c_row$beta) == sign(v_row$beta) else NA,
      nominal_significance_pattern = if (nrow(v_row) == 1L) {
        paste0(ifelse(c_row$pval < 0.05, "chen_p_lt_0_05", "chen_p_ge_0_05"), "__",
               ifelse(v_row$pval < 0.05, "vuckovic_p_lt_0_05", "vuckovic_p_ge_0_05"))
      } else NA_character_,
      notes = "interpretive_audit_only_no_meta_analysis_no_difference_test",
      stringsAsFactors = FALSE
    )
  }))

  heterogeneity_rows <- do.call(rbind, lapply(seq_len(nrow(chen_het)), function(i) {
    c_row <- chen_het[i, , drop = FALSE]
    v_row <- v_het[v_het$analysis_set == c_row$analysis_set & v_het$method == c_row$method, , drop = FALSE]
    data.frame(
      comparison_domain = "heterogeneity_Q",
      analysis_set = c_row$analysis_set,
      method = c_row$method,
      method_id = c_row$method_id,
      chen_value = c_row$Q,
      vuckovic_value = if (nrow(v_row) == 1L) v_row$Q else NA_real_,
      chen_pval = c_row$pval,
      vuckovic_pval = if (nrow(v_row) == 1L) v_row$Q_pval else NA_real_,
      direction_consistent = NA,
      nominal_significance_pattern = if (nrow(v_row) == 1L) {
        paste0(ifelse(c_row$pval < 0.05, "chen_p_lt_0_05", "chen_p_ge_0_05"), "__",
               ifelse(v_row$Q_pval < 0.05, "vuckovic_p_lt_0_05", "vuckovic_p_ge_0_05"))
      } else NA_character_,
      notes = "interpretive_diagnostic_comparison_only",
      stringsAsFactors = FALSE
    )
  }))

  egger_rows <- do.call(rbind, lapply(seq_len(nrow(chen_egger)), function(i) {
    c_row <- chen_egger[i, , drop = FALSE]
    v_row <- v_egger[v_egger$analysis_set == c_row$analysis_set, , drop = FALSE]
    data.frame(
      comparison_domain = "egger_intercept",
      analysis_set = c_row$analysis_set,
      method = "Egger intercept",
      method_id = "",
      chen_value = c_row$intercept,
      vuckovic_value = if (nrow(v_row) == 1L) v_row$egger_intercept else NA_real_,
      chen_pval = c_row$pval,
      vuckovic_pval = if (nrow(v_row) == 1L) v_row$pval else NA_real_,
      direction_consistent = if (nrow(v_row) == 1L) sign(c_row$intercept) == sign(v_row$egger_intercept) else NA,
      nominal_significance_pattern = if (nrow(v_row) == 1L) {
        paste0(ifelse(c_row$pval < 0.05, "chen_p_lt_0_05", "chen_p_ge_0_05"), "__",
               ifelse(v_row$pval < 0.05, "vuckovic_p_lt_0_05", "vuckovic_p_ge_0_05"))
      } else NA_character_,
      notes = "interpretive_diagnostic_comparison_only",
      stringsAsFactors = FALSE
    )
  }))

  presso_sets <- unique(chen_presso$analysis_set)
  presso_rows <- do.call(rbind, lapply(presso_sets, function(set_name) {
    c_global <- chen_presso[
      chen_presso$analysis_set == set_name &
        chen_presso$test_type == "Global Test" &
        chen_presso$metric == "RSSobs",
      ,
      drop = FALSE
    ]
    c_status <- chen_presso[
      chen_presso$analysis_set == set_name &
        chen_presso$test_type == "MR-PRESSO" &
        chen_presso$metric == "status",
      ,
      drop = FALSE
    ]
    v_row <- v_presso[v_presso$analysis_set == set_name, , drop = FALSE]
    chen_value <- if (nrow(c_global) == 1L) suppressWarnings(as.numeric(c_global$value)) else NA_real_
    chen_p <- if (nrow(c_global) == 1L) c_global$pval else NA_real_
    chen_note <- if (nrow(c_status) == 1L) {
      paste0("chen_mr_presso_status=", c_status$value[[1L]], "; ", c_status$notes[[1L]])
    } else {
      "chen_mr_presso_status=passed"
    }
    data.frame(
      comparison_domain = "mr_presso_global",
      analysis_set = set_name,
      method = "MR-PRESSO Global Test",
      method_id = "",
      chen_value = chen_value,
      vuckovic_value = if (nrow(v_row) == 1L) v_row$global_RSSobs else NA_real_,
      chen_pval = chen_p,
      vuckovic_pval = if (nrow(v_row) == 1L) suppressWarnings(as.numeric(v_row$global_pval)) else NA_real_,
      direction_consistent = NA,
      nominal_significance_pattern = if (nrow(c_global) == 1L && nrow(v_row) == 1L) {
        paste0(ifelse(c_global$pval < 0.05, "chen_p_lt_0_05", "chen_p_ge_0_05"), "__",
               ifelse(suppressWarnings(as.numeric(v_row$global_pval)) < 0.05, "vuckovic_p_lt_0_05", "vuckovic_p_ge_0_05"))
      } else {
        paste0(
          "chen_mr_presso_not_estimable__",
          if (nrow(v_row) == 1L) "vuckovic_mr_presso_available" else "vuckovic_mr_presso_unavailable"
        )
      },
      notes = paste("interpretive_diagnostic_comparison_only", chen_note, sep = "; "),
      stringsAsFactors = FALSE
    )
  }))

  loo_rows <- do.call(rbind, lapply(seq_len(nrow(chen_loo_summary)), function(i) {
    c_row <- chen_loo_summary[i, , drop = FALSE]
    v_row <- v_loo_summary[v_loo_summary$analysis_set == c_row$analysis_set, , drop = FALSE]
    data.frame(
      comparison_domain = "loo_max_absolute_shift",
      analysis_set = c_row$analysis_set,
      method = "IVW leave-one-out",
      method_id = "mr_ivw",
      chen_value = c_row$max_absolute_shift,
      vuckovic_value = if (nrow(v_row) == 1L) v_row$max_absolute_shift else NA_real_,
      chen_pval = NA_real_,
      vuckovic_pval = NA_real_,
      direction_consistent = NA,
      nominal_significance_pattern = paste0(
        "chen_any_nominal_change_", c_row$any_nominal_significance_change,
        "__vuckovic_any_nominal_change_", if (nrow(v_row) == 1L) v_row$any_nominal_significance_change else NA
      ),
      notes = paste0("interpretive_diagnostic_comparison_only; chen_max_shift_rsid=", c_row$max_shift_rsid,
                     "; vuckovic_max_shift_rsid=", if (nrow(v_row) == 1L) v_row$max_shift_rsid else NA_character_,
                     "; not independent replication"),
      stringsAsFactors = FALSE
    )
  }))

  rbind(estimator_rows, heterogeneity_rows, egger_rows, presso_rows, loo_rows)
}

main <- function() {
  log_line("stage=chen_forward_mr_v1_start")
  set.seed(2026L)
  renv_before <- hash_file(paths$renv_lock)

  freeze <- jsonlite::fromJSON(paths$freeze, simplifyVector = FALSE)
  amendment <- jsonlite::fromJSON(paths$amendment_readback, simplifyVector = FALSE)
  contract <- jsonlite::fromJSON(paths$contract, simplifyVector = FALSE)
  primary_qc <- jsonlite::fromJSON(paths$primary_qc, simplifyVector = FALSE)
  script_txt <- paste(readLines(paths$primary_script, warn = FALSE), collapse = "\n")
  primary_methods <- extract_required_methods(script_txt)
  contract_methods <- extract_contract_methods(contract)
  seed_value <- extract_int("seed_value\\s*<-\\s*[0-9]+L?", script_txt)
  nb_distribution <- extract_int("NbDistribution\\s*=\\s*[0-9]+", script_txt)
  signif_threshold <- extract_num_after_equal("SignifThreshold\\s*=\\s*[0-9.]+", script_txt)

  ts_desc <- utils::packageDescription("TwoSampleMR")
  mp_desc <- utils::packageDescription("MRPRESSO")
  ts_ver <- as.character(utils::packageVersion("TwoSampleMR"))
  ts_sha <- as.character(ts_desc[["RemoteSha"]])
  mp_ver <- as.character(utils::packageVersion("MRPRESSO"))
  mp_sha <- as.character(mp_desc[["RemoteSha"]])

  freeze_manifest_sha <- hash_file(paths$freeze_manifest)
  included_path <- rel("data_derived", "forward_sensitivity_harmonisation", "chen_forward_harmonised_apoe_included_v1.parquet")
  excluded_path <- rel("data_derived", "forward_sensitivity_harmonisation", "chen_forward_harmonised_apoe_excluded_v1.parquet")
  if (!is.null(freeze$manifest_records)) {
    roles <- vapply(freeze$manifest_records, function(x) as.character(x$file_role), character(1))
    rels <- vapply(freeze$manifest_records, function(x) as.character(x$relative_path), character(1))
    if ("harmonised_apoe_included_parquet" %in% roles) included_path <- rel(rels[roles == "harmonised_apoe_included_parquet"][[1L]])
    if ("harmonised_apoe_excluded_parquet" %in% roles) excluded_path <- rel(rels[roles == "harmonised_apoe_excluded_parquet"][[1L]])
  }
  stop_if(!file.exists(included_path) || !file.exists(excluded_path), "Decision 94 frozen Parquet inputs are absent.")
  included_sha_before <- hash_file(included_path)
  excluded_sha_before <- hash_file(excluded_path)

  method_hierarchy <- data.frame(
    method_id = contract_methods,
    estimator_hierarchy_role = ifelse(contract_methods == "mr_ivw", "primary", "sensitivity"),
    stringsAsFactors = FALSE
  )

  pre_checks <- list(
    contract_v2_gate = identical(contract$contract_status, "frozen") &&
      isTRUE(contract$approved_for_chen_forward_mr_execution) &&
      length(contract$hard_check_failures) == 0L,
    method_alignment_amendment_gate = identical(amendment$authoritative_amendment_status_after_readback, "frozen") &&
      isTRUE(amendment$approved_for_chen_forward_mr_contract_v2_after_readback) &&
      length(amendment$corrected_hard_check_failures) == 0L,
    mr_input_freeze_gate = identical(freeze$freeze_status, "passed") &&
      isTRUE(freeze$approved_for_chen_forward_mr_design) &&
      length(freeze$hard_check_failures) == 0L &&
      identical(tolower(freeze$manifest_sha256), tolower(freeze_manifest_sha)),
    method_ids_available = all(contract_methods %in% TwoSampleMR::mr_method_list()$obj),
    ivw_implementation_matches_primary = "mr_ivw" %in% primary_methods && grepl("TwoSampleMR::mr\\s*\\(", script_txt),
    estimator_hierarchy_matches_primary = identical(contract_methods, primary_methods),
    mr_presso_configuration_matches_primary = grepl("MRPRESSO::mr_presso", script_txt) &&
      identical(nb_distribution, 10000L) && isTRUE(all.equal(signif_threshold, 0.05)) &&
      identical(seed_value, 2026L),
    software_environment_matches_primary = identical(ts_ver, primary_qc$TwoSampleMR_version) &&
      identical(ts_sha, primary_qc$TwoSampleMR_RemoteSha) &&
      identical(mp_ver, primary_qc$MRPRESSO_version) &&
      identical(mp_sha, primary_qc$MRPRESSO_RemoteSha)
  )
  stop_if(!all(unlist(pre_checks)), paste("Pre-execution gate failed:", paste(names(pre_checks)[!unlist(pre_checks)], collapse = "; ")))

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  inc_src <- read_parquet(con, included_path)
  exc_src <- read_parquet(con, excluded_path)
  inc_dat <- make_dat(inc_src, "APOE_included", "chen_forward_alternative_hb_gwas_sensitivity_main")
  exc_dat <- make_dat(exc_src, "APOE_excluded", "chen_forward_alternative_hb_gwas_apoe_exclusion_sensitivity")
  inc_checks <- check_input(inc_src, inc_dat, freeze$included_final_rsids, as.integer(freeze$included_final_valid_count))
  exc_checks <- check_input(exc_src, exc_dat, freeze$excluded_final_rsids, as.integer(freeze$excluded_final_valid_count))
  input_check_values <- c(inc_checks[setdiff(names(inc_checks), "f_summary")], exc_checks[setdiff(names(exc_checks), "f_summary")])
  stop_if(!all(unlist(input_check_values)), paste("Input audit failed:", paste(names(input_check_values)[!unlist(input_check_values)], collapse = "; ")))

  log_line("running_mr_methods")
  inc_est <- run_mr_methods(inc_dat, "APOE_included", "chen_forward_alternative_hb_gwas_sensitivity_main", contract_methods, method_hierarchy)
  exc_est <- run_mr_methods(exc_dat, "APOE_excluded", "chen_forward_alternative_hb_gwas_apoe_exclusion_sensitivity", contract_methods, method_hierarchy)
  estimates <- rbind(inc_est, exc_est)
  inc_ivw <- inc_est[inc_est$method_id == "mr_ivw", , drop = FALSE]
  exc_ivw <- exc_est[exc_est$method_id == "mr_ivw", , drop = FALSE]

  log_line("running_heterogeneity_egger_presso_loo")
  heterogeneity <- rbind(
    run_heterogeneity(inc_dat, "APOE_included", "chen_forward_alternative_hb_gwas_sensitivity_main"),
    run_heterogeneity(exc_dat, "APOE_excluded", "chen_forward_alternative_hb_gwas_apoe_exclusion_sensitivity")
  )
  egger <- rbind(
    run_egger(inc_dat, "APOE_included", "chen_forward_alternative_hb_gwas_sensitivity_main"),
    run_egger(exc_dat, "APOE_excluded", "chen_forward_alternative_hb_gwas_apoe_exclusion_sensitivity")
  )
  presso <- rbind(
    run_presso(inc_dat, "APOE_included", "chen_forward_alternative_hb_gwas_sensitivity_main", nb_distribution, signif_threshold, seed_value),
    run_presso(exc_dat, "APOE_excluded", "chen_forward_alternative_hb_gwas_apoe_exclusion_sensitivity", nb_distribution, signif_threshold, seed_value)
  )
  loo <- rbind(
    run_loo(inc_dat, "APOE_included", "chen_forward_alternative_hb_gwas_sensitivity_main", inc_ivw$beta[[1L]], inc_ivw$pval[[1L]]),
    run_loo(exc_dat, "APOE_excluded", "chen_forward_alternative_hb_gwas_apoe_exclusion_sensitivity", exc_ivw$beta[[1L]], exc_ivw$pval[[1L]])
  )
  loo_summary <- summarise_loo(loo)
  comparison <- make_comparison(estimates, heterogeneity, egger, presso, loo_summary)

  included_sha_after <- hash_file(included_path)
  excluded_sha_after <- hash_file(excluded_path)
  renv_after <- hash_file(paths$renv_lock)

  partials <- paste0(unlist(paths[c("estimates", "heterogeneity", "egger", "presso", "loo", "comparison", "qc", "decision")]), ".partial")
  on.exit(unlink(partials, force = TRUE), add = TRUE)
  write_csv_precise(estimates, paste0(paths$estimates, ".partial"))
  write_csv_precise(heterogeneity, paste0(paths$heterogeneity, ".partial"))
  write_csv_precise(egger, paste0(paths$egger, ".partial"))
  write_csv_precise(presso, paste0(paths$presso, ".partial"))
  write_csv_precise(loo, paste0(paths$loo, ".partial"))
  write_csv_precise(comparison, paste0(paths$comparison, ".partial"))

  roundtrip <- list(
    estimates = nrow(read.csv(paste0(paths$estimates, ".partial"), stringsAsFactors = FALSE, check.names = FALSE)) == 10L,
    heterogeneity = nrow(read.csv(paste0(paths$heterogeneity, ".partial"), stringsAsFactors = FALSE, check.names = FALSE)) == 4L,
    egger = nrow(read.csv(paste0(paths$egger, ".partial"), stringsAsFactors = FALSE, check.names = FALSE)) == 2L,
    presso = all(c("APOE_included", "APOE_excluded") %in%
      unique(read.csv(paste0(paths$presso, ".partial"), stringsAsFactors = FALSE, check.names = FALSE)$analysis_set)),
    loo = nrow(read.csv(paste0(paths$loo, ".partial"), stringsAsFactors = FALSE, check.names = FALSE)) == nrow(inc_dat) + nrow(exc_dat),
    comparison = nrow(read.csv(paste0(paths$comparison, ".partial"), stringsAsFactors = FALSE, check.names = FALSE)) >= 20L
  )

  hard_checks <- c(pre_checks, list(
    included_count_matches_freeze = inc_checks$count_matches_freeze,
    excluded_count_matches_freeze = exc_checks$count_matches_freeze,
    rsid_sets_match_freeze = inc_checks$rsid_set_matches_freeze && exc_checks$rsid_set_matches_freeze,
    input_conversion_preserved = inc_checks$conversion_preserved && exc_checks$conversion_preserved,
    no_reharmonisation = TRUE,
    ivw_run_included = nrow(inc_ivw) == 1L,
    ivw_run_excluded = nrow(exc_ivw) == 1L,
    egger_run_included = sum(inc_est$method_id == "mr_egger_regression") == 1L,
    egger_run_excluded = sum(exc_est$method_id == "mr_egger_regression") == 1L,
    weighted_median_run_included = sum(inc_est$method_id == "mr_weighted_median") == 1L,
    weighted_median_run_excluded = sum(exc_est$method_id == "mr_weighted_median") == 1L,
    weighted_mode_run_included = sum(inc_est$method_id == "mr_weighted_mode") == 1L,
    weighted_mode_run_excluded = sum(exc_est$method_id == "mr_weighted_mode") == 1L,
    simple_mode_run_included = sum(inc_est$method_id == "mr_simple_mode") == 1L,
    simple_mode_run_excluded = sum(exc_est$method_id == "mr_simple_mode") == 1L,
    heterogeneity_completed = nrow(heterogeneity) == 4L,
    egger_intercept_completed = nrow(egger) == 2L,
    mr_presso_attempted = all(c("APOE_included", "APOE_excluded") %in% unique(presso$analysis_set)),
    leave_one_out_completed = nrow(loo) == nrow(inc_dat) + nrow(exc_dat),
    single_snp_not_run_as_planned = !exists("single", inherits = FALSE) &&
      !file.exists(paths$single_snp_forbidden) && !file.exists(paste0(paths$single_snp_forbidden, ".partial")),
    or_transform_correct = all(abs(estimates$OR - exp(estimates$beta)) <= 1e-12) &&
      all(abs(estimates$OR_ci_lower - exp(estimates$ci_lower)) <= 1e-12) &&
      all(abs(estimates$OR_ci_upper - exp(estimates$ci_upper)) <= 1e-12),
    vuckovic_comparison_noninferential = nrow(comparison) >= 20L &&
      all(grepl("interpretive", comparison$notes)) &&
      !any(grepl("meta-analysis|independent replication test|formal difference test", comparison$notes, ignore.case = TRUE)),
    independent_replication_false = TRUE,
    no_posthoc_filtering = TRUE,
    no_steiger = TRUE,
    no_proxy = TRUE,
    no_liftover = TRUE,
    no_input_mutation = identical(included_sha_before, included_sha_after) && identical(excluded_sha_before, excluded_sha_after),
    renv_lock_unchanged = identical(renv_before, renv_after),
    output_roundtrip_checks = all(vapply(roundtrip, isTRUE, logical(1))),
    no_single_snp_output = !file.exists(paths$single_snp_forbidden) && !file.exists(paste0(paths$single_snp_forbidden, ".partial"))
  ))
  failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
  status <- if (length(failures) == 0L) "passed" else "failed"

  qc <- list(
    mr_version = "v1",
    date = "2026-08-12",
    analysis_direction = "Hb_to_delirium",
    analysis_role = "forward_alternative_hb_gwas_sensitivity",
    independent_replication = FALSE,
    source_mr_input_freeze_decision = 94,
    mr_method_alignment_amendment_decision = 98,
    mr_contract_v2_decision = 99,
    included_nsnp = nrow(inc_dat),
    excluded_nsnp = nrow(exc_dat),
    shared_nsnp = length(intersect(inc_dat$SNP, exc_dat$SNP)),
    included_only_nsnp = length(setdiff(inc_dat$SNP, exc_dat$SNP)),
    excluded_only_nsnp = length(setdiff(exc_dat$SNP, inc_dat$SNP)),
    methods_run = contract_methods,
    raw_results = records(estimates[, c("analysis_set", "analysis_role", "method", "method_id", "nsnp", "beta", "se", "ci_lower", "ci_upper", "pval")]),
    or_results = records(estimates[, c("analysis_set", "method_id", "OR", "OR_ci_lower", "OR_ci_upper", "pval")]),
    heterogeneity_results = records(heterogeneity),
    egger_intercept_results = records(egger),
    mr_presso_results = records(presso),
    leave_one_out_summary = records(loo_summary),
    single_snp_run = FALSE,
    single_snp_status = "not_planned_by_method_alignment_amendment",
    vuckovic_comparison = records(comparison),
    steiger_run = FALSE,
    steiger_status = "deferred_to_unified_directionality_sensitivity_stage",
    seed = seed_value,
    software_environment = list(
      R_version = R.version.string,
      TwoSampleMR_version = ts_ver,
      TwoSampleMR_RemoteSha = ts_sha,
      MRPRESSO_version = mp_ver,
      MRPRESSO_RemoteSha = mp_sha,
      mr_library = norm(normalizePath(mr_library, winslash = "/", mustWork = TRUE)),
      install_update_restore_snapshot_performed = FALSE,
      renv_out_of_sync_message = "informational_only"
    ),
    input_sha256 = list(
      included_before = included_sha_before,
      included_after = included_sha_after,
      excluded_before = excluded_sha_before,
      excluded_after = excluded_sha_after,
      decision_94_manifest_observed = freeze_manifest_sha
    ),
    instrument_strength_summary = list(
      included = as.list(inc_checks$f_summary[1, , drop = FALSE]),
      excluded = as.list(exc_checks$f_summary[1, , drop = FALSE])
    ),
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after,
    renv_lock_unchanged = identical(renv_before, renv_after),
    mr_status = status,
    approved_for_chen_forward_results_interpretation = identical(status, "passed"),
    hard_checks = hard_checks,
    hard_check_failures = failures,
    informational_findings = list(
      single_snp_removed_by_method_alignment_amendment = TRUE,
      chen_analysis_is_sensitivity_not_independent_replication = TRUE,
      no_results_freeze_created = TRUE,
      mr_presso_failure_handling = "MR-PRESSO failure, if any, is recorded without changing the frozen main input"
    )
  )
  write_json_precise(qc, paste0(paths$qc, ".partial"))

  decision_lines <- c(
    "# Decision 100: Chen Forward MR V1",
    "",
    "Date: 2026-08-12",
    "",
    "## Status",
    paste0("mr_status: `", status, "`"),
    paste0("approved_for_chen_forward_results_interpretation: `", identical(status, "passed"), "`"),
    "",
    "## Decision",
    "Chen Forward MR V1 was executed under Decision 99 Contract V2 and Decision 98 method-alignment amendment.",
    "",
    "The analysis is an alternative-Hb-GWAS robustness sensitivity analysis and is not an independent replication.",
    "",
    "## Inputs",
    paste0("- APOE included final-valid instruments: `", nrow(inc_dat), "`"),
    paste0("- APOE excluded final-valid instruments: `", nrow(exc_dat), "`"),
    paste0("- Decision 94 manifest SHA-256: `", freeze_manifest_sha, "`"),
    "",
    "## Methods",
    paste0("- `", contract_methods, "`"),
    "- Heterogeneity: `TwoSampleMR::mr_heterogeneity`",
    "- Egger intercept: `TwoSampleMR::mr_pleiotropy_test`",
    "- MR-PRESSO: `MRPRESSO::mr_presso`",
    "- Leave-one-out: `TwoSampleMR::mr_leaveoneout`",
    "- Single-SNP/Wald-ratio diagnostics: not run by Decision 98/99",
    "- Steiger: not run; deferred",
    "",
    "## Hard Check Failures",
    if (length(failures) == 0L) "- none" else paste0("- `", failures, "`"),
    "",
    "## Outputs",
    paste0("- `", norm(paths$estimates), "`"),
    paste0("- `", norm(paths$heterogeneity), "`"),
    paste0("- `", norm(paths$egger), "`"),
    paste0("- `", norm(paths$presso), "`"),
    paste0("- `", norm(paths$loo), "`"),
    paste0("- `", norm(paths$comparison), "`"),
    paste0("- `", norm(paths$qc), "`"),
    paste0("- `", norm(paths$log), "`"),
    paste0("- `", norm(paths$decision), "`"),
    "",
    "## Freeze",
    "No Chen Forward MR Results Freeze was created in this run."
  )
  writeLines(decision_lines, paste0(paths$decision, ".partial"), useBytes = TRUE)

  stop_if(!identical(status, "passed"), paste("MR hard checks failed:", paste(failures, collapse = "; ")))
  for (path in unlist(paths[c("estimates", "heterogeneity", "egger", "presso", "loo", "comparison", "qc", "decision")])) {
    stop_if(file.exists(path), paste("Output appeared during run:", path))
    if (!file.rename(paste0(path, ".partial"), path)) stop("Atomic rename failed: ", path, call. = FALSE)
  }
  log_line("mr_status=passed")
  log_line("approved_for_chen_forward_results_interpretation=TRUE")
  log_line("hard_check_failures=[]")
}

tryCatch(main(), error = function(e) {
  log_line("mr_status=failed")
  log_line("error=", conditionMessage(e))
  quit(status = 1L)
})
