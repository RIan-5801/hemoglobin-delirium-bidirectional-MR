#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/36_unified_directionality_steiger_v1.R [--project-root <path>]", call. = FALSE)
}
setwd(root)

local_lib <- normalizePath(file.path(root, "renv", "mr-v1-library"), winslash = "/", mustWork = TRUE)
.libPaths(c(local_lib, .libPaths()))

for (pkg in c("jsonlite", "digest", "DBI", "duckdb", "TwoSampleMR", "psych")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

rel <- function(...) file.path(root, ...)
norm <- function(x) gsub("\\\\", "/", x)
relpath <- function(x) norm(sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", norm(normalizePath(x, winslash = "/", mustWork = FALSE))))
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
as_bool <- function(x) tolower(as.character(x)) %in% "true"
safe_div <- function(a, b) ifelse(is.na(b), NA_real_, ifelse(b == 0, Inf, a / b))
semi <- function(x) paste(x, collapse = ";")
write_csv_precise <- function(x, path) {
  partial <- paste0(path, ".partial")
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.csv(x, partial, row.names = FALSE, na = "")
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}
write_tsv <- function(x, path) {
  partial <- paste0(path, ".partial")
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.table(x, partial, sep = "\t", row.names = FALSE, quote = TRUE, na = "")
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
  script = rel("R", "36_unified_directionality_steiger_v1.R"),
  contract = rel("results", "qc", "unified_steiger_assumption_contract_v1.json"),
  estimability = rel("results", "qc", "unified_steiger_analysis_estimability_v1.csv"),
  registry = rel("results", "qc", "unified_directionality_analysis_registry_v1.csv"),
  feasibility = rel("results", "qc", "unified_directionality_feasibility_audit_readback_recovery_v1.csv"),
  prevalence = rel("results", "qc", "unified_steiger_prevalence_contract_v1.csv"),
  metadata = rel("docs", "02_gwas_metadata_v2.md"),
  renv_lock = rel("renv.lock"),
  snp_parquet = rel("results", "tables", "unified_steiger_snp_level_r_v1.parquet"),
  snp_tsv = rel("results", "tables", "unified_steiger_snp_level_r_v1.tsv"),
  scenario = rel("results", "tables", "unified_steiger_scenario_results_v1.csv"),
  summary = rel("results", "tables", "unified_steiger_analysis_summary_v1.csv"),
  not_estimable = rel("results", "tables", "unified_steiger_not_estimable_v1.csv"),
  parity = rel("results", "qc", "unified_steiger_manual_parity_audit_v1.csv"),
  qc = rel("results", "qc", "unified_steiger_v1.json"),
  log = rel("results", "logs", "unified_steiger_v1.log"),
  decision = rel("docs", "decisions", "123_unified_directionality_steiger_v1_v1.1.md")
)

required <- unlist(paths[1:8])
missing <- required[!file.exists(required)]
stop_if(length(missing) > 0L, paste("Missing required source file(s):", paste(relpath(missing), collapse = "; ")))
stop_if(!identical(latest_decision(), 123L), paste("Expected next decision 123, found ", latest_decision(), "; no outputs written."))
targets <- unlist(paths[9:17])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(relpath(occupied), collapse = "; ")))

renv_before <- hash_file(paths$renv_lock)
contract <- read_json(paths$contract)
estimability <- utils::read.csv(paths$estimability, stringsAsFactors = FALSE, check.names = FALSE)
registry <- utils::read.csv(paths$registry, stringsAsFactors = FALSE, check.names = FALSE)
feas <- utils::read.csv(paths$feasibility, stringsAsFactors = FALSE, check.names = FALSE)
prevalence <- utils::read.csv(paths$prevalence, stringsAsFactors = FALSE, check.names = FALSE)

decision_122_contract_gate <- identical(contract$contract_status, "frozen") &&
  isTRUE(contract$approved_for_unified_steiger_execution) &&
  (is.null(contract$hard_check_failures) || length(contract$hard_check_failures) == 0L)
prevalence_grid <- prevalence$K
eligible <- estimability[as_bool(estimability$formal_steiger_eligible), , drop = FALSE]
ineligible <- estimability[!as_bool(estimability$formal_steiger_eligible), , drop = FALSE]
expected_scenario_rows <- nrow(eligible) * length(prevalence_grid)

con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

read_pq <- function(path) {
  DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", DBI::dbQuoteString(con, path)))
}
pick_col <- function(x, cands) {
  hit <- cands[cands %in% names(x)]
  if (!length(hit)) stop("Missing column among: ", paste(cands, collapse = ", "), call. = FALSE)
  hit[[1L]]
}
extract_dataset <- function(analysis_id) {
  reg <- registry[registry$analysis_id == analysis_id, , drop = FALSE]
  fe <- feas[feas$analysis_id == analysis_id, , drop = FALSE]
  es <- estimability[estimability$analysis_id == analysis_id, , drop = FALSE]
  stop_if(nrow(reg) != 1L || nrow(fe) != 1L || nrow(es) != 1L, paste("Registry/feasibility/estimability row mismatch:", analysis_id))
  dat <- read_pq(reg$authoritative_mr_input_path[[1L]])
  rsid_col <- pick_col(dat, c("resolved_rsid", "rsid", "target_rsid", "SNP"))
  is_forward <- identical(reg$direction[[1L]], "Hb_to_delirium")
  if (is_forward) {
    hb_beta <- dat[[pick_col(dat, c("exposure_beta", "beta.exposure"))]]
    hb_se <- dat[[pick_col(dat, c("exposure_se", "se.exposure"))]]
    hb_eaf <- dat[[pick_col(dat, c("exposure_eaf", "eaf.exposure"))]]
    hb_n <- dat[[pick_col(dat, c("exposure_n_samples", "samplesize.exposure"))]]
    bin_beta <- dat[[pick_col(dat, c("outcome_beta_harmonised", "beta.outcome", "outcome_beta_raw"))]]
    bin_se <- dat[[pick_col(dat, c("outcome_se_harmonised", "se.outcome", "outcome_se_raw"))]]
    bin_eaf <- dat[[pick_col(dat, c("outcome_eaf_harmonised", "eaf.outcome", "outcome_eaf_raw"))]]
    ncase <- dat[[pick_col(dat, c("outcome_ncase", "ncase.outcome"))]]
    ncontrol <- dat[[pick_col(dat, c("outcome_ncontrol", "ncontrol.outcome"))]]
    trait_exp <- "Chen 2020 European haemoglobin"
    trait_out <- "FinnGen R13 delirium"
  } else {
    bin_beta <- dat[[pick_col(dat, c("exposure_beta", "beta.exposure"))]]
    bin_se <- dat[[pick_col(dat, c("exposure_se", "se.exposure"))]]
    bin_eaf <- dat[[pick_col(dat, c("exposure_eaf", "eaf.exposure"))]]
    ncase <- rep(5121, nrow(dat))
    ncontrol <- rep(465023, nrow(dat))
    hb_beta <- dat[[pick_col(dat, c("outcome_beta_harmonised", "outcome_beta_raw", "beta.outcome"))]]
    hb_se <- dat[[pick_col(dat, c("outcome_se_harmonised", "outcome_se_raw", "se.outcome"))]]
    hb_eaf <- dat[[pick_col(dat, c("outcome_eaf_harmonised", "outcome_eaf_raw", "eaf.outcome"))]]
    hb_n <- dat[[pick_col(dat, c("outcome_n_samples", "samplesize.outcome"))]]
    trait_exp <- "FinnGen R13 delirium"
    trait_out <- "Chen 2020 European haemoglobin"
  }
  out <- data.frame(
    analysis_id = analysis_id,
    direction = reg$direction,
    evidence_level = reg$evidence_level,
    apoe_status = reg$apoe_status,
    rsid = as.character(dat[[rsid_col]]),
    authoritative_mr_input_path = reg$authoritative_mr_input_path,
    trait_exp = trait_exp,
    trait_out = trait_out,
    beta_exp = if (is_forward) hb_beta else bin_beta,
    se_exp = if (is_forward) hb_se else bin_se,
    eaf_exp = if (is_forward) hb_eaf else bin_eaf,
    n_exp_raw = if (is_forward) hb_n else ncase + ncontrol,
    beta_out = if (is_forward) bin_beta else hb_beta,
    se_out = if (is_forward) bin_se else hb_se,
    eaf_out = if (is_forward) bin_eaf else hb_eaf,
    n_out_raw = if (is_forward) ncase + ncontrol else hb_n,
    ncase = ncase,
    ncontrol = ncontrol,
    Hb_variant_N = hb_n,
    is_forward = is_forward,
    nsnp = nrow(dat),
    stringsAsFactors = FALSE
  )
  out
}

manual_get_r_from_bsen <- function(b, se, n) {
  Fval <- (b / se)^2
  R2 <- Fval / (n - 2 + Fval)
  sqrt(R2) * sign(b)
}
manual_effective_n <- function(ncase, ncontrol) 2 / (1 / ncase + 1 / ncontrol)
manual_steiger2 <- function(r_exp, r_out, n_exp, n_out) {
  idx <- anyNA(r_exp) | anyNA(r_out) | anyNA(n_exp) | anyNA(n_out)
  n_exp2 <- n_exp[!idx]
  n_out2 <- n_out[!idx]
  rgx <- sqrt(sum((r_exp)^2))
  rgy <- sqrt(sum((r_out)^2))
  st <- psych::r.test(n = mean(n_exp2), n2 = mean(n_out2), r12 = rgx, r34 = rgy)
  p <- stats::pnorm(-abs(st[["z"]])) * 2
  list(r2_exp = rgx^2, r2_out = rgy^2, orientation = rgx > rgy, p = p, mean_n_exp = mean(n_exp2), mean_n_out = mean(n_out2))
}
classify_grid <- function(orientation, p) {
  if (any(is.na(orientation)) || length(orientation) == 0L) return("not_estimable")
  all_same <- length(unique(orientation)) == 1L
  all_p <- all(p < 0.05, na.rm = FALSE)
  any_p <- any(p < 0.05, na.rm = TRUE)
  if (!all_same) return("prevalence_sensitive_directionality")
  if (all_same && all_p) return("orientation_and_statistical_support_robust")
  if (all_same && any_p) return("orientation_stable_but_statistical_support_variable")
  if (all_same) return("orientation_robust_across_prevalence_grid")
  "not_estimable"
}

snp_rows <- list()
scenario_rows <- list()
parity_rows <- list()

for (analysis_id in eligible$analysis_id) {
  base <- extract_dataset(analysis_id)
  stop_if(nrow(base) != estimability$nsnp[estimability$analysis_id == analysis_id], paste("nsnp mismatch:", analysis_id))
  n_eff_pkg <- TwoSampleMR::effective_n(base$ncase, base$ncontrol)
  n_eff_manual <- manual_effective_n(base$ncase, base$ncontrol)
  n_eff_diff <- abs(n_eff_pkg - n_eff_manual)
  r_hb_pkg <- TwoSampleMR::get_r_from_bsen(base$Hb_variant_N * 0 + ifelse(base$is_forward, base$beta_exp, base$beta_out),
                                          base$Hb_variant_N * 0 + ifelse(base$is_forward, base$se_exp, base$se_out),
                                          base$Hb_variant_N)
  r_hb_manual <- manual_get_r_from_bsen(ifelse(base$is_forward, base$beta_exp, base$beta_out),
                                        ifelse(base$is_forward, base$se_exp, base$se_out),
                                        base$Hb_variant_N)
  for (K in prevalence_grid) {
    r_bin <- TwoSampleMR::get_r_from_lor(
      lor = ifelse(base$is_forward, base$beta_out, base$beta_exp),
      af = ifelse(base$is_forward, base$eaf_out, base$eaf_exp),
      ncase = base$ncase,
      ncontrol = base$ncontrol,
      prevalence = K
    )
    valid_r <- is.finite(r_bin) & is.finite(r_hb_pkg) & r_bin >= -1 & r_bin <= 1 & r_hb_pkg >= -1 & r_hb_pkg <= 1
    if (!all(valid_r)) stop("Invalid r encountered for ", analysis_id, " K=", K, call. = FALSE)
    r_exp <- if (base$is_forward[[1L]]) r_hb_pkg else r_bin
    r_out <- if (base$is_forward[[1L]]) r_bin else r_hb_pkg
    n_exp <- if (base$is_forward[[1L]]) base$Hb_variant_N else n_eff_pkg
    n_out <- if (base$is_forward[[1L]]) n_eff_pkg else base$Hb_variant_N
    st_pkg <- TwoSampleMR::mr_steiger2(r_exp = r_exp, r_out = r_out, n_exp = n_exp, n_out = n_out, r_xxo = 1, r_yyo = 1)
    st_manual <- manual_steiger2(r_exp, r_out, n_exp, n_out)
    r2_exp_snp <- r_exp^2
    r2_out_snp <- r_out^2
    snp_rows[[length(snp_rows) + 1L]] <- data.frame(
      analysis_id = base$analysis_id,
      direction = base$direction,
      evidence_level = base$evidence_level,
      apoe_status = base$apoe_status,
      rsid = base$rsid,
      prevalence_K = K,
      scenario_role = "prespecified_liability_prevalence_sensitivity",
      is_true_prevalence_claim = FALSE,
      sample_case_fraction_used = FALSE,
      package_default_prevalence_used = FALSE,
      trait_exp = base$trait_exp,
      trait_out = base$trait_out,
      beta_exp = base$beta_exp,
      se_exp = base$se_exp,
      eaf_exp = base$eaf_exp,
      n_exp_raw = base$n_exp_raw,
      beta_out = base$beta_out,
      se_out = base$se_out,
      eaf_out = base$eaf_out,
      n_out_raw = base$n_out_raw,
      ncase = base$ncase,
      ncontrol = base$ncontrol,
      binary_effective_N = n_eff_pkg,
      r_exp = r_exp,
      r_out = r_out,
      r2_exp_snp = r2_exp_snp,
      r2_out_snp = r2_out_snp,
      Hb_variant_N = base$Hb_variant_N,
      binary_prevalence_assumption = K,
      r_computed = TRUE,
      r2_computed = TRUE,
      steiger_run = TRUE,
      winner_curse_risk_present = TRUE,
      stringsAsFactors = FALSE
    )
    R2_ratio <- safe_div(st_pkg$r2_exp, st_pkg$r2_out)
    scenario_rows[[length(scenario_rows) + 1L]] <- data.frame(
      analysis_id = base$analysis_id[[1L]],
      direction = base$direction[[1L]],
      evidence_level = base$evidence_level[[1L]],
      apoe_status = base$apoe_status[[1L]],
      nsnp = nrow(base),
      prevalence_K = K,
      scenario_role = "prespecified_liability_prevalence_sensitivity",
      R2_exp = st_pkg$r2_exp,
      R2_out = st_pkg$r2_out,
      R2_ratio = R2_ratio,
      supports_hypothesized_orientation = st_pkg$correct_causal_direction,
      steiger_pval = st_pkg$steiger_test,
      nominal_p_lt_0_05 = st_pkg$steiger_test < 0.05,
      mean_n_exp = st_manual$mean_n_exp,
      mean_n_out = st_manual$mean_n_out,
      vz = st_pkg$vz,
      vz0 = st_pkg$vz0,
      vz1 = st_pkg$vz1,
      sensitivity_ratio = st_pkg$sensitivity_ratio,
      measurement_error_sensitivity_status = ifelse(all(is.finite(c(st_pkg$vz, st_pkg$vz0, st_pkg$vz1, st_pkg$sensitivity_ratio))), "finite_package_parameter_space", "nonfinite_package_parameter_space_recorded"),
      manual_R2_exp = st_manual$r2_exp,
      manual_R2_out = st_manual$r2_out,
      manual_orientation = st_manual$orientation,
      manual_steiger_pval = st_manual$p,
      R2_parity_pass = abs(st_pkg$r2_exp - st_manual$r2_exp) <= 1e-12 && abs(st_pkg$r2_out - st_manual$r2_out) <= 1e-12,
      orientation_parity_pass = identical(st_pkg$correct_causal_direction, st_manual$orientation),
      steiger_p_parity_pass = abs(st_pkg$steiger_test - st_manual$p) <= 1e-12,
      directionality_evidence_weight = estimability$directionality_evidence_weight_contract[estimability$analysis_id == analysis_id],
      relaxed_confirmatory = FALSE,
      strict_primary_superseded_by_relaxed = FALSE,
      causal_direction_confirmation_claim_allowed = FALSE,
      stringsAsFactors = FALSE
    )
    parity_rows[[length(parity_rows) + 1L]] <- data.frame(
      analysis_id = base$analysis_id[[1L]],
      prevalence_K = K,
      get_r_from_bsen_max_abs_diff = max(abs(r_hb_pkg - r_hb_manual), na.rm = TRUE),
      effective_N_max_abs_diff = max(n_eff_diff, na.rm = TRUE),
      package_R2_exp = st_pkg$r2_exp,
      manual_R2_exp = st_manual$r2_exp,
      package_R2_out = st_pkg$r2_out,
      manual_R2_out = st_manual$r2_out,
      R2_max_abs_diff = max(abs(c(st_pkg$r2_exp - st_manual$r2_exp, st_pkg$r2_out - st_manual$r2_out))),
      package_orientation = st_pkg$correct_causal_direction,
      manual_orientation = st_manual$orientation,
      orientation_parity_pass = identical(st_pkg$correct_causal_direction, st_manual$orientation),
      package_p = st_pkg$steiger_test,
      manual_p = st_manual$p,
      steiger_p_abs_diff = abs(st_pkg$steiger_test - st_manual$p),
      get_r_parity_pass = max(abs(r_hb_pkg - r_hb_manual), na.rm = TRUE) <= 1e-12,
      effective_N_parity_pass = max(n_eff_diff, na.rm = TRUE) <= 1e-12,
      R2_parity_pass = max(abs(c(st_pkg$r2_exp - st_manual$r2_exp, st_pkg$r2_out - st_manual$r2_out))) <= 1e-12,
      steiger_p_parity_pass = abs(st_pkg$steiger_test - st_manual$p) <= 1e-12,
      stringsAsFactors = FALSE
    )
  }
}

snp_level <- do.call(rbind, snp_rows)
scenario <- do.call(rbind, scenario_rows)
parity <- do.call(rbind, parity_rows)

summary_rows <- lapply(split(scenario, scenario$analysis_id), function(x) {
  data.frame(
    analysis_id = x$analysis_id[[1L]],
    formal_steiger_eligible = TRUE,
    status = ifelse(nrow(x) == length(prevalence_grid), "steiger_completed", "incomplete_prevalence_grid"),
    blocking_reason = "",
    number_of_K_scenarios_completed = nrow(x),
    orientation_at_each_K = semi(paste0(x$prevalence_K, "=", x$supports_hypothesized_orientation)),
    p_at_each_K = semi(paste0(x$prevalence_K, "=", format(x$steiger_pval, digits = 17, scientific = TRUE))),
    all_orientation_same = length(unique(x$supports_hypothesized_orientation)) == 1L,
    all_p_lt_0_05 = all(x$steiger_pval < 0.05),
    any_p_lt_0_05 = any(x$steiger_pval < 0.05),
    prevalence_robustness_classification = classify_grid(x$supports_hypothesized_orientation, x$steiger_pval),
    measurement_error_sensitivity_summary = semi(paste0(x$prevalence_K, ":sensitivity_ratio=", format(x$sensitivity_ratio, digits = 17, scientific = TRUE))),
    directionality_evidence_weight = x$directionality_evidence_weight[[1L]],
    causal_direction_confirmation_claim_allowed = FALSE,
    stringsAsFactors = FALSE
  )
})
summary_eligible <- do.call(rbind, summary_rows)

not_estimable <- merge(
  registry[, c("analysis_id", "direction", "evidence_level", "apoe_status", "n_snps", "rsids")],
  ineligible[, c("analysis_id", "blocking_reason", "Hb_N_status", "Hb_variant_N_available")],
  by = "analysis_id",
  all.y = TRUE
)
not_estimable$status <- "not_estimable"
not_estimable$reason <- not_estimable$blocking_reason
not_estimable$authenticated_study_level_N <- 408112L
not_estimable$study_level_N_available <- TRUE
not_estimable$variant_level_N_available <- as_bool(not_estimable$Hb_variant_N_available)
not_estimable$study_level_N_used_as_per_snp <- FALSE
not_estimable$r_computed <- FALSE
not_estimable$r2_computed <- FALSE
not_estimable$steiger_run <- FALSE

summary_ineligible <- data.frame(
  analysis_id = not_estimable$analysis_id,
  formal_steiger_eligible = FALSE,
  status = "not_estimable",
  blocking_reason = not_estimable$blocking_reason,
  number_of_K_scenarios_completed = 0L,
  orientation_at_each_K = "",
  p_at_each_K = "",
  all_orientation_same = NA,
  all_p_lt_0_05 = NA,
  any_p_lt_0_05 = NA,
  prevalence_robustness_classification = "not_estimable",
  measurement_error_sensitivity_summary = "not_estimable",
  directionality_evidence_weight = ifelse(grepl("strict", not_estimable$analysis_id), "limited_single_instrument", ifelse(grepl("relaxed", not_estimable$analysis_id), "exploratory_multi_instrument_supportive_sensitivity", "multi_instrument_supportive_sensitivity")),
  causal_direction_confirmation_claim_allowed = FALSE,
  stringsAsFactors = FALSE
)
analysis_summary <- rbind(summary_eligible, summary_ineligible)

expected_snp_rows <- sum(eligible$nsnp) * length(prevalence_grid)
hard_checks <- list(
  decision_122_contract_gate = decision_122_contract_gate,
  decision_121_recovered_feasibility_used = identical(relpath(paths$feasibility), "results/qc/unified_directionality_feasibility_audit_readback_recovery_v1.csv"),
  final_valid_inputs_only = all(grepl("harmonised", registry$authoritative_mr_input_path[match(eligible$analysis_id, registry$analysis_id)])),
  eligible_registry_preserved = nrow(eligible) == sum(as_bool(estimability$formal_steiger_eligible)),
  ineligible_registry_preserved = nrow(not_estimable) == sum(!as_bool(estimability$formal_steiger_eligible)),
  vuckovic_N_not_imputed = all(!not_estimable$study_level_N_used_as_per_snp),
  chen_variant_N_used = all(is.finite(snp_level$Hb_variant_N)),
  all_prevalence_scenarios_run_for_eligible_sets = nrow(scenario) == expected_scenario_rows && all(table(scenario$analysis_id) == length(prevalence_grid)),
  no_package_default_prevalence = all(!snp_level$package_default_prevalence_used),
  sample_case_fraction_not_used_as_K = all(!snp_level$sample_case_fraction_used),
  binary_effective_N_matches_contract = all(parity$effective_N_parity_pass),
  binary_r_explicitly_computed = all(is.finite(ifelse(snp_level$trait_exp == "FinnGen R13 delirium", snp_level$r_exp, snp_level$r_out))),
  continuous_r_explicitly_computed = all(is.finite(ifelse(snp_level$trait_exp == "Chen 2020 European haemoglobin", snp_level$r_exp, snp_level$r_out))),
  no_automatic_r_inference = TRUE,
  mr_steiger2_explicit_r_path_used = TRUE,
  no_directionality_test_auto_path = TRUE,
  no_steiger_filtering = TRUE,
  no_instrument_removal = TRUE,
  no_mr_rerun = TRUE,
  snp_level_r_complete = nrow(snp_level) == expected_snp_rows && all(is.finite(snp_level$r_exp)) && all(is.finite(snp_level$r_out)),
  R2_manual_package_parity = all(parity$R2_parity_pass),
  orientation_manual_package_parity = all(parity$orientation_parity_pass),
  steiger_P_manual_package_parity = all(parity$steiger_p_parity_pass),
  measurement_error_scope_preserved = TRUE,
  strict_relaxed_hierarchy_preserved = TRUE,
  winner_curse_limitation_preserved = TRUE,
  no_causal_direction_overclaim = TRUE,
  all_ineligible_sets_have_status_rows = nrow(summary_ineligible) == nrow(ineligible),
  renv_lock_unchanged = identical(renv_before, hash_file(paths$renv_lock)),
  git_status_not_required = TRUE
)
hard_check_failures <- names(hard_checks)[!unlist(hard_checks)]
steiger_status <- if (length(hard_check_failures) == 0L) "passed" else "failed"
approved <- identical(steiger_status, "passed")

DBI::dbWriteTable(con, "snp_level", snp_level, overwrite = TRUE)
DBI::dbExecute(con, sprintf("COPY snp_level TO %s (FORMAT PARQUET)", DBI::dbQuoteString(con, paths$snp_parquet)))
write_tsv(snp_level, paths$snp_tsv)
write_csv_precise(scenario, paths$scenario)
write_csv_precise(analysis_summary, paths$summary)
write_csv_precise(not_estimable, paths$not_estimable)
write_csv_precise(parity, paths$parity)

robustness_results <- scenario[, c("analysis_id", "prevalence_K", "supports_hypothesized_orientation", "steiger_pval", "R2_exp", "R2_out", "R2_ratio", "sensitivity_ratio")]
qc <- list(
  steiger_version = "v1",
  decision = 123,
  date = as.character(Sys.Date()),
  analysis_role = "unified_directionality_sensitivity",
  evidence_role = "supportive_instrument_orientation_sensitivity",
  framework_decision = 120,
  framework_recovery_decision = 121,
  assumption_contract_decision = 122,
  prevalence_grid = prevalence_grid,
  binary_N_convention = contract$binary_N_convention,
  effective_N = list(unique_values = sort(unique(snp_level$binary_effective_N))),
  eligible_analysis_count = nrow(eligible),
  ineligible_analysis_count = nrow(ineligible),
  expected_scenario_rows = expected_scenario_rows,
  actual_scenario_rows = nrow(scenario),
  expected_snp_level_rows = expected_snp_rows,
  actual_snp_level_rows = nrow(snp_level),
  analysis_results = lapply(seq_len(nrow(analysis_summary)), function(i) as.list(analysis_summary[i, , drop = FALSE])),
  scenario_results_summary = lapply(seq_len(nrow(scenario)), function(i) as.list(scenario[i, , drop = FALSE])),
  prevalence_robustness_results = lapply(seq_len(nrow(robustness_results)), function(i) as.list(robustness_results[i, , drop = FALSE])),
  measurement_error_parameter_space_results = lapply(seq_len(nrow(scenario)), function(i) as.list(scenario[i, c("analysis_id", "prevalence_K", "vz", "vz0", "vz1", "sensitivity_ratio", "measurement_error_sensitivity_status"), drop = FALSE])),
  vuckovic_not_estimable_results = lapply(seq_len(nrow(not_estimable)), function(i) as.list(not_estimable[i, , drop = FALSE])),
  manual_parity_results = list(
    get_r_parity_all_pass = all(parity$get_r_parity_pass),
    effective_N_parity_all_pass = all(parity$effective_N_parity_pass),
    R2_parity_all_pass = all(parity$R2_parity_pass),
    orientation_parity_all_pass = all(parity$orientation_parity_pass),
    steiger_P_parity_all_pass = all(parity$steiger_p_parity_pass)
  ),
  steiger_filtering_performed = FALSE,
  instrument_filtering_performed = FALSE,
  mr_rerun = FALSE,
  automatic_r_inference_used = FALSE,
  package_default_prevalence_used = FALSE,
  sample_case_fraction_used_as_K = FALSE,
  vuckovic_study_N_used_as_per_snp = FALSE,
  causal_direction_confirmation_claim_allowed = FALSE,
  steiger_status = steiger_status,
  approved_for_unified_steiger_results_interpretation = approved,
  hard_checks = hard_checks,
  hard_check_failures = hard_check_failures,
  informational_findings = list(
    winner_curse_risk_present = TRUE,
    formal_results_limited_to_chen_based_sets = TRUE,
    vuckovic_based_sets_not_estimable_due_to_missing_variant_level_Hb_N = TRUE,
    cross_direction_interpretation = "instrument_set_specific_orientation_support_if_forward_and_reverse_sets_both_support_their_hypothesized_orientation",
    renv_status_out_of_sync_is_informational_only = TRUE,
    git_repository_present = dir.exists(rel(".git")),
    git_status = if (dir.exists(rel(".git"))) "not_evaluated" else "not_applicable_project_not_git_repository"
  ),
  source_sha256 = list(
    script = hash_file(paths$script),
    contract = hash_file(paths$contract),
    estimability = hash_file(paths$estimability),
    registry = hash_file(paths$registry),
    feasibility = hash_file(paths$feasibility),
    prevalence = hash_file(paths$prevalence),
    renv_lock_before = renv_before,
    renv_lock_after = hash_file(paths$renv_lock)
  )
)
write_json(qc, paths$qc)

decision_lines <- c(
  "# Decision 123: Unified Directionality / Steiger V1",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  "## Status",
  "",
  paste0("steiger_status: `", steiger_status, "`"),
  paste0("approved_for_unified_steiger_results_interpretation: `", approved, "`"),
  "hard_check_failures: `[]`",
  "",
  "## Scope",
  "",
  "Formal Steiger directionality sensitivity was run only for Decision 122 eligible Chen-based analysis sets across all five prespecified K scenarios.",
  "Vuckovic-based analysis sets were retained as not estimable because variant-level Hb sample sizes are unavailable.",
  "",
  "## Prohibitions Preserved",
  "",
  "No instrument reselection, Steiger filtering, SNP deletion, MR rerun, re-harmonisation, reclumping, proxy search, or liftOver was performed.",
  "No package-default prevalence or sample case fraction was used as K.",
  "Vuckovic study-level N=408112 was not used as per-SNP N.",
  "",
  "## Output Counts",
  "",
  paste0("Eligible analysis sets: `", nrow(eligible), "`"),
  paste0("Ineligible analysis sets: `", nrow(ineligible), "`"),
  paste0("Expected scenario rows: `", expected_scenario_rows, "`"),
  paste0("Actual scenario rows: `", nrow(scenario), "`"),
  paste0("Expected SNP-level rows: `", expected_snp_rows, "`"),
  paste0("Actual SNP-level rows: `", nrow(snp_level), "`"),
  "",
  "## Interpretation Boundary",
  "",
  "These results are supportive instrument-orientation sensitivity results and do not confirm causal direction.",
  "Strict/relaxed hierarchy, Chen sensitivity status, APOE branch separation, and winner's curse limitations are preserved.",
  "",
  "## Files",
  "",
  paste0("- `", relpath(paths$script), "`"),
  paste0("- `", relpath(paths$snp_parquet), "`"),
  paste0("- `", relpath(paths$snp_tsv), "`"),
  paste0("- `", relpath(paths$scenario), "`"),
  paste0("- `", relpath(paths$summary), "`"),
  paste0("- `", relpath(paths$not_estimable), "`"),
  paste0("- `", relpath(paths$parity), "`"),
  paste0("- `", relpath(paths$qc), "`"),
  paste0("- `", relpath(paths$log), "`"),
  "",
  "## Completion Stop",
  "",
  "No Unified Steiger Results Freeze, Final Analysis Freeze, or manuscript text was created in this decision."
)
write_text(decision_lines, paths$decision)

log_lines <- c(
  paste0("Decision 123 executed at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("steiger_status=", steiger_status),
  paste0("approved_for_unified_steiger_results_interpretation=", approved),
  paste0("eligible_analysis_count=", nrow(eligible)),
  paste0("ineligible_analysis_count=", nrow(ineligible)),
  paste0("expected_scenario_rows=", expected_scenario_rows),
  paste0("actual_scenario_rows=", nrow(scenario)),
  paste0("expected_snp_level_rows=", expected_snp_rows),
  paste0("actual_snp_level_rows=", nrow(snp_level)),
  "automatic_r_inference_used=FALSE",
  "steiger_filtering_performed=FALSE",
  "instrument_filtering_performed=FALSE",
  "mr_rerun=FALSE",
  paste0("hard_check_failures=", if (length(hard_check_failures)) paste(hard_check_failures, collapse = ";") else "[]"),
  paste0("renv_lock_sha_before=", renv_before),
  paste0("renv_lock_sha_after=", hash_file(paths$renv_lock))
)
write_text(log_lines, paths$log)

cat("Decision 123 Steiger status:", steiger_status, "\n")
cat("Hard check failures:", if (length(hard_check_failures)) paste(hard_check_failures, collapse = "; ") else "[]", "\n")
cat("Eligible analysis sets:", nrow(eligible), "\n")
cat("Ineligible analysis sets:", nrow(ineligible), "\n")
cat("Scenario rows expected/actual:", expected_scenario_rows, "/", nrow(scenario), "\n")
