args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") stop("Usage: Rscript 09_forward_mr_vuckovic_hb_finngen_r13_primary_v3.R --project-root <path>", call. = FALSE)

project_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
stage <- "startup"
mr_library <- file.path(project_root, "renv", "mr-v1-library")
renv_activate <- file.path(project_root, "renv", "activate.R")
renv_lock <- file.path(project_root, "renv.lock")
harmonisation_qc_path <- file.path(project_root, "results", "qc", "vuckovic_hb_finngen_r13_primary_harmonisation_v4.json")
loo_audit_qc_path <- file.path(project_root, "results", "qc", "vuckovic_hb_finngen_r13_forward_loo_label_audit_v4.json")
included_input_path <- file.path(project_root, "data_derived", "harmonised", "vuckovic_hb_finngen_r13_primary_apoe_included_v4.parquet")
excluded_input_path <- file.path(project_root, "data_derived", "harmonised", "vuckovic_hb_finngen_r13_primary_apoe_excluded_v4.parquet")
mr_input_dir <- file.path(project_root, "data_derived", "mr_inputs")
tables_dir <- file.path(project_root, "results", "tables")
qc_dir <- file.path(project_root, "results", "qc")
logs_dir <- file.path(project_root, "results", "logs")

mr_input_included_path <- file.path(mr_input_dir, "vuckovic_hb_finngen_r13_forward_primary_apoe_included_v3.parquet")
mr_input_excluded_path <- file.path(mr_input_dir, "vuckovic_hb_finngen_r13_forward_primary_apoe_excluded_v3.parquet")
estimates_path <- file.path(tables_dir, "vuckovic_hb_finngen_r13_forward_mr_estimates_v3.csv")
heterogeneity_path <- file.path(tables_dir, "vuckovic_hb_finngen_r13_forward_heterogeneity_v3.csv")
egger_path <- file.path(tables_dir, "vuckovic_hb_finngen_r13_forward_egger_intercept_v3.csv")
strength_path <- file.path(tables_dir, "vuckovic_hb_finngen_r13_forward_instrument_strength_v3.csv")
strength_summary_path <- file.path(tables_dir, "vuckovic_hb_finngen_r13_forward_instrument_strength_summary_v3.csv")
loo_path <- file.path(tables_dir, "vuckovic_hb_finngen_r13_forward_leave_one_out_v3.csv")
loo_full_path <- file.path(tables_dir, "vuckovic_hb_finngen_r13_forward_leave_one_out_full_ivw_v3.csv")
presso_path <- file.path(tables_dir, "vuckovic_hb_finngen_r13_forward_mr_presso_v3.csv")
qc_path <- file.path(qc_dir, "vuckovic_hb_finngen_r13_forward_mr_v3.json")
log_path <- file.path(logs_dir, "vuckovic_hb_finngen_r13_forward_mr_v3.log")

expected_twosamplemr_version <- "0.7.9"
expected_twosamplemr_sha <- "3d119f20d6fc164b0c7f710f5590fee9580f2c7b"
expected_mrpresso_version <- "1.0"
expected_mrpresso_sha <- "3e3c92d7eda6dce0d1d66077373ec0f7ff4f7e87"
expected_included_sha <- "882b0e9f1c1567f0d4265e2c9b95c7ac9cbe327ce7256a8f0c9b7b5b12995aea"
expected_excluded_sha <- "183034c0b5d10035fc671b06ee964b9fe2042a9e8cebe36f45af695fc5b8c0c1"
required_methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode")
seed_value <- 2026L
effect_interpretation <- "per genetically predicted 1-unit increase in standardized inverse-normal-transformed haemoglobin"

stop_if <- function(condition, message) if (isTRUE(condition)) stop(message, call. = FALSE)
safe_log <- function(...) cat(paste0(...), "\n", file = log_path, append = TRUE)
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
numeric_equal <- function(x, y, tolerance = 1e-12) is.finite(x) && is.finite(y) && abs(as.numeric(x) - as.numeric(y)) <= tolerance
read_parquet <- function(con, path) DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", normalizePath(path, winslash = "/", mustWork = TRUE)))
write_parquet <- function(con, data, path) {
  table_name <- "forward_mr_v3_write"
  DBI::dbWriteTable(con, table_name, data, temporary = TRUE, overwrite = TRUE)
  DBI::dbExecute(con, sprintf("COPY %s TO '%s' (FORMAT PARQUET)", table_name, normalizePath(path, winslash = "/", mustWork = FALSE)))
  DBI::dbRemoveTable(con, table_name)
}
write_csv_precise <- function(data, path) {
  old <- options(digits = 17, scipen = 999); on.exit(options(old), add = TRUE)
  write.csv(data, path, row.names = FALSE, na = "")
}
write_partial <- function(writer, final_path) writer(paste0(final_path, ".partial"))
publish_partial <- function(final_path) {
  if (!file.rename(paste0(final_path, ".partial"), final_path)) stop(paste0("Atomic rename failed: ", final_path), call. = FALSE)
}

make_mr_data <- function(source, analysis_set, analysis_role) {
  required <- c("rsid", "exposure_beta", "exposure_se", "effect_allele", "other_allele", "exposure_eaf", "outcome_beta", "outcome_se", "outcome_effect_allele_aligned", "outcome_other_allele_aligned", "outcome_eaf", "mr_keep_primary", "palindromic_snp", "exposure_n_study", "exposure_effect_scale", "outcome_ncase", "outcome_ncontrol", "outcome_n_study", "outcome_effect_scale", "outcome_pval")
  stop_if(!all(required %in% names(source)), "Required harmonised input columns are absent.")
  data.frame(SNP = as.character(source$rsid), beta.exposure = as.numeric(source$exposure_beta), se.exposure = as.numeric(source$exposure_se), effect_allele.exposure = as.character(source$effect_allele), other_allele.exposure = as.character(source$other_allele), eaf.exposure = as.numeric(source$exposure_eaf), beta.outcome = as.numeric(source$outcome_beta), se.outcome = as.numeric(source$outcome_se), effect_allele.outcome = as.character(source$outcome_effect_allele_aligned), other_allele.outcome = as.character(source$outcome_other_allele_aligned), eaf.outcome = as.numeric(source$outcome_eaf), exposure = "Vuckovic_2020_Hb", outcome = "FinnGen_R13_F5_DELIRIUM", id.exposure = "vuckovic_hb_2020", id.outcome = "finngen_R13_F5_DELIRIUM", mr_keep = TRUE, analysis_set = analysis_set, analysis_role = analysis_role, stringsAsFactors = FALSE)
}

check_input <- function(source, dat) list(
  input_count = nrow(source) == 312L, rsid_unique = !anyNA(source$rsid) && !anyDuplicated(source$rsid),
  all_mr_keep = all(source$mr_keep_primary) && all(dat$mr_keep), no_palindromic_primary = !any(source$palindromic_snp),
  aligned_alleles = all(dat$effect_allele.exposure == dat$effect_allele.outcome) && all(dat$other_allele.exposure == dat$other_allele.outcome),
  finite_exposure_effects = all(is.finite(dat$beta.exposure)), positive_exposure_se = all(is.finite(dat$se.exposure) & dat$se.exposure > 0),
  finite_outcome_effects = all(is.finite(dat$beta.outcome)), positive_outcome_se = all(is.finite(dat$se.outcome) & dat$se.outcome > 0),
  valid_outcome_p = all(is.finite(source$outcome_pval) & source$outcome_pval >= 0 & source$outcome_pval <= 1),
  exposure_metadata = all(source$exposure_n_study == 408112) && all(source$exposure_effect_scale == "standardized_inverse_normal_transformed_Hb"),
  outcome_metadata = all(source$outcome_ncase == 5121) && all(source$outcome_ncontrol == 465023) && all(source$outcome_n_study == 470144) && all(source$outcome_effect_scale == "log_odds")
)

run_mr_methods <- function(dat, analysis_set, analysis_role) {
  x <- TwoSampleMR::mr(dat, method_list = required_methods)
  if (is.list(x) && "mr" %in% names(x)) x <- x$mr
  stop_if(!is.data.frame(x) || !all(c("method", "nsnp", "b", "se", "pval") %in% names(x)), "MR result schema is invalid.")
  methods <- TwoSampleMR::mr_method_list(); lookup <- stats::setNames(as.character(methods$obj), as.character(methods$name))
  out <- data.frame(analysis_set = analysis_set, analysis_role = analysis_role, method_id = unname(lookup[as.character(x$method)]), method = as.character(x$method), nsnp = as.integer(x$nsnp), beta = as.numeric(x$b), se = as.numeric(x$se), pval = as.numeric(x$pval), stringsAsFactors = FALSE)
  out$ci_lower_beta <- out$beta - 1.96 * out$se; out$ci_upper_beta <- out$beta + 1.96 * out$se
  out$OR <- exp(out$beta); out$OR_lci <- exp(out$ci_lower_beta); out$OR_uci <- exp(out$ci_upper_beta)
  out
}
run_heterogeneity <- function(dat, analysis_set, analysis_role) {
  x <- TwoSampleMR::mr_heterogeneity(dat, method_list = c("mr_egger_regression", "mr_ivw"))
  stop_if(!all(c("method", "Q", "Q_df", "Q_pval") %in% names(x)), "Heterogeneity result schema is invalid.")
  data.frame(analysis_set = analysis_set, analysis_role = analysis_role, method = as.character(x$method), Q = as.numeric(x$Q), Q_df = as.numeric(x$Q_df), Q_pval = as.numeric(x$Q_pval), stringsAsFactors = FALSE)
}
run_egger <- function(dat, analysis_set, analysis_role) {
  x <- TwoSampleMR::mr_pleiotropy_test(dat)
  stop_if(!all(c("egger_intercept", "se", "pval") %in% names(x)), "Egger-intercept result schema is invalid.")
  data.frame(analysis_set = analysis_set, analysis_role = analysis_role, egger_intercept = as.numeric(x$egger_intercept[[1L]]), se = as.numeric(x$se[[1L]]), pval = as.numeric(x$pval[[1L]]), stringsAsFactors = FALSE)
}
run_loo <- function(dat, analysis_set, analysis_role, full_ivw) {
  raw <- TwoSampleMR::mr_leaveoneout(dat, method = TwoSampleMR::mr_ivw)
  if (!is.data.frame(raw)) raw <- as.data.frame(raw)
  stop_if(!all(c("SNP", "b", "se", "p") %in% names(raw)), "LOO does not expose the V4-verified SNP/b/se/p interface.")
  labels <- as.character(raw$SNP); overall <- !is.na(labels) & grepl("^All", labels)
  stop_if(sum(overall) != 1L, "LOO requires exactly one All sentinel.")
  loo_snps <- labels[!overall]; input_snps <- as.character(dat$SNP)
  stop_if(length(loo_snps) != 312L || length(unique(loo_snps)) != 312L || !setequal(input_snps, loo_snps) || length(setdiff(input_snps, loo_snps)) != 0L || length(setdiff(loo_snps, input_snps)) != 0L, "LOO real-SNP set gate failed.")
  input_freq <- table(input_snps); loo_freq <- table(loo_snps)
  stop_if(!all(input_freq == 1L) || !all(loo_freq == 1L) || !identical(sort(names(input_freq)), sort(names(loo_freq))), "LOO real-SNP frequency gate failed.")
  snp_raw <- raw[!overall, , drop = FALSE]; all_raw <- raw[overall, , drop = FALSE]
  loo <- data.frame(analysis_set = analysis_set, analysis_role = analysis_role, excluded_rsid = as.character(snp_raw$SNP), nsnp_remaining = 311L, beta = as.numeric(snp_raw$b), se = as.numeric(snp_raw$se), pval = as.numeric(snp_raw$p), stringsAsFactors = FALSE)
  loo$ci_lower_beta <- loo$beta - 1.96 * loo$se; loo$ci_upper_beta <- loo$beta + 1.96 * loo$se
  loo$OR <- exp(loo$beta); loo$OR_lci <- exp(loo$ci_lower_beta); loo$OR_uci <- exp(loo$ci_upper_beta)
  overall_match <- numeric_equal(all_raw$b[[1L]], full_ivw$beta[[1L]]) && numeric_equal(all_raw$se[[1L]], full_ivw$se[[1L]]) && numeric_equal(all_raw$p[[1L]], full_ivw$pval[[1L]])
  full <- data.frame(analysis_set = analysis_set, analysis_role = analysis_role, overall_label = as.character(all_raw$SNP[[1L]]), beta = as.numeric(all_raw$b[[1L]]), se = as.numeric(all_raw$se[[1L]]), pval = as.numeric(all_raw$p[[1L]]), stringsAsFactors = FALSE)
  full$OR <- exp(full$beta); full$OR_lci <- exp(full$beta - 1.96 * full$se); full$OR_uci <- exp(full$beta + 1.96 * full$se)
  list(loo = loo, full = full, overall_match = overall_match, raw_rows = nrow(raw), real_snp_rows = length(loo_snps), real_snp_unique = length(unique(loo_snps)), exact_set = setequal(input_snps, loo_snps), frequency_match = identical(sort(names(input_freq)), sort(names(loo_freq))) && all(input_freq == 1L) && all(loo_freq == 1L))
}
summarise_strength <- function(dat, analysis_set, analysis_role) {
  f <- (dat$beta.exposure / dat$se.exposure)^2
  per_snp <- data.frame(analysis_set = analysis_set, analysis_role = analysis_role, rsid = dat$SNP, exposure_beta = dat$beta.exposure, exposure_se = dat$se.exposure, F_statistic = f, weak_instrument_F_lt_10 = f < 10, stringsAsFactors = FALSE)
  summary <- data.frame(analysis_set = analysis_set, analysis_role = analysis_role, nsnp = nrow(dat), F_min = min(f), F_mean = mean(f), F_median = stats::median(f), F_max = max(f), F_lt_10_count = sum(f < 10), stringsAsFactors = FALSE)
  list(per_snp = per_snp, summary = summary)
}
run_presso <- function(dat, analysis_set, analysis_role) {
  set.seed(seed_value); x <- dat[, c("beta.outcome", "beta.exposure", "se.outcome", "se.exposure"), drop = FALSE]; rownames(x) <- dat$SNP
  result <- tryCatch(MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = x, NbDistribution = 10000, SignifThreshold = 0.05), error = function(e) e)
  if (inherits(result, "error")) return(data.frame(analysis_set = analysis_set, analysis_role = analysis_role, mr_presso_status = "failed", global_RSSobs = NA_real_, global_pval = NA_character_, outlier_count = NA_integer_, outlier_rsids = NA_character_, mr_presso_error = conditionMessage(result), stringsAsFactors = FALSE))
  root <- result[["MR-PRESSO results"]]; global <- if (!is.null(root)) root[["Global Test"]] else NULL; outlier <- if (!is.null(root)) root[["Outlier Test"]] else NULL
  ids <- if (!is.null(outlier) && nrow(outlier) > 0L) rownames(outlier) else character(0)
  data.frame(analysis_set = analysis_set, analysis_role = analysis_role, mr_presso_status = "passed", global_RSSobs = if (!is.null(global) && "RSSobs" %in% names(global)) as.numeric(global$RSSobs[[1L]]) else NA_real_, global_pval = if (!is.null(global) && "Pvalue" %in% names(global)) as.character(global$Pvalue[[1L]]) else NA_character_, outlier_count = length(ids), outlier_rsids = paste(ids, collapse = ";"), mr_presso_error = NA_character_, stringsAsFactors = FALSE)
}

main <- function() {
  output_targets <- c(mr_input_included_path, mr_input_excluded_path, estimates_path, heterogeneity_path, egger_path, strength_path, strength_summary_path, loo_path, loo_full_path, presso_path, qc_path, log_path)
  stop_if(any(file.exists(c(output_targets, paste0(output_targets, ".partial")))), "Forward MR V3 target or partial exists; refusing to overwrite.")
  stop_if(!file.exists(renv_activate) || !dir.exists(mr_library) || !file.exists(renv_lock), "Frozen MR runtime paths are absent.")
  dir.create(mr_input_dir, recursive = TRUE, showWarnings = FALSE); dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE); dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE); dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  source(renv_activate); .libPaths(c(mr_library, .libPaths()))
  deps <- list(dependency_library_present = dir.exists(mr_library), twosamplemr_available = requireNamespace("TwoSampleMR", quietly = TRUE), mrpresso_available = requireNamespace("MRPRESSO", quietly = TRUE))
  stop_if(!all(vapply(deps, isTRUE, logical(1))), "Frozen MR packages are unavailable.")
  ts_desc <- utils::packageDescription("TwoSampleMR"); mp_desc <- utils::packageDescription("MRPRESSO")
  ts_ver <- as.character(utils::packageVersion("TwoSampleMR")); mp_ver <- as.character(utils::packageVersion("MRPRESSO")); ts_sha <- as.character(ts_desc$RemoteSha); mp_sha <- as.character(mp_desc$RemoteSha)
  deps$twosamplemr_version_exact <- identical(ts_ver, expected_twosamplemr_version); deps$twosamplemr_sha_exact <- identical(ts_sha, expected_twosamplemr_sha); deps$mrpresso_version_exact <- identical(mp_ver, expected_mrpresso_version); deps$mrpresso_sha_exact <- identical(mp_sha, expected_mrpresso_sha)
  stop_if(!all(vapply(deps, isTRUE, logical(1))), "Frozen MR package version or RemoteSha gate failed.")
  set.seed(seed_value); renv_lock_sha_before <- hash_file(renv_lock)
  safe_log("R_version=", R.version.string); safe_log("TwoSampleMR_version=", ts_ver); safe_log("TwoSampleMR_RemoteSha=", ts_sha); safe_log("MRPRESSO_version=", mp_ver); safe_log("MRPRESSO_RemoteSha=", mp_sha); safe_log("seed=", seed_value)

  stage <<- "authoritative_gates"
  stop_if(!all(file.exists(c(harmonisation_qc_path, loo_audit_qc_path, included_input_path, excluded_input_path))), "Authoritative V4 inputs are absent.")
  hqc <- jsonlite::fromJSON(harmonisation_qc_path, simplifyVector = FALSE); lqc <- jsonlite::fromJSON(loo_audit_qc_path, simplifyVector = FALSE)
  harmonisation_gate <- identical(hqc$harmonisation_status, "passed") && isTRUE(hqc$approved_for_forward_primary_mr) && length(hqc$hard_check_failures) == 0L
  loo_gate_set <- function(x) isTRUE(x$source_vs_tsmr_exact_set_match) && isTRUE(x$source_vs_tsmr_exact_order_match) && identical(as.integer(x$analysis_stratum_count), 1L) && identical(as.integer(x$loo_total_row_count), 313L) && identical(as.integer(x$loo_overall_row_count), 1L) && identical(as.integer(x$loo_non_overall_row_count), 312L) && identical(as.integer(x$loo_non_overall_distinct_count), 312L) && isTRUE(x$exact_set_match) && identical(as.integer(x$input_not_in_loo_count), 0L) && identical(as.integer(x$loo_not_in_input_count), 0L) && identical(as.integer(x$loo_duplicate_count), 0L)
  loo_v4_audit_gate <- identical(lqc$diagnostic_status, "completed") && isTRUE(lqc$approved_for_loo_mapping_fix) && identical(lqc$explicit_diagnostic_reason, "v2_loo_set_check_implementation_error") && loo_gate_set(lqc$included) && loo_gate_set(lqc$excluded) && length(lqc$hard_check_failures) == 0L
  stop_if(!harmonisation_gate || !loo_v4_audit_gate, "Harmonisation or LOO V4 gate failed.")
  included_sha_before <- hash_file(included_input_path); excluded_sha_before <- hash_file(excluded_input_path)
  input_sha_matches_v4_audit <- identical(tolower(included_sha_before), expected_included_sha) && identical(tolower(excluded_sha_before), expected_excluded_sha)
  stop_if(!input_sha_matches_v4_audit, "Harmonised input SHA-256 differs from the V4 LOO audit.")

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE)); on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  included_source <- read_parquet(con, included_input_path); excluded_source <- read_parquet(con, excluded_input_path)
  included_mr <- make_mr_data(included_source, "APOE_included", "forward_primary"); excluded_mr <- make_mr_data(excluded_source, "APOE_excluded", "forward_apoe_exclusion_sensitivity")
  inc_checks <- check_input(included_source, included_mr); exc_checks <- check_input(excluded_source, excluded_mr)
  shared_nsnp <- length(intersect(included_mr$SNP, excluded_mr$SNP)); included_only_nsnp <- length(setdiff(included_mr$SNP, excluded_mr$SNP)); excluded_only_nsnp <- length(setdiff(excluded_mr$SNP, included_mr$SNP))
  membership <- list(shared_snp_count = shared_nsnp == 311L, included_only_count = included_only_nsnp == 1L, excluded_only_count = excluded_only_nsnp == 1L)
  stop_if(!all(vapply(c(inc_checks, exc_checks, membership), isTRUE, logical(1))), "Formal MR input QC gate failed.")

  stage <<- "formal_mr"
  methods_available <- all(required_methods %in% TwoSampleMR::mr_method_list()$obj); stop_if(!methods_available, "Prespecified MR method is unavailable.")
  inc_strength <- summarise_strength(included_mr, "APOE_included", "forward_primary"); exc_strength <- summarise_strength(excluded_mr, "APOE_excluded", "forward_apoe_exclusion_sensitivity")
  inc_est <- run_mr_methods(included_mr, "APOE_included", "forward_primary"); exc_est <- run_mr_methods(excluded_mr, "APOE_excluded", "forward_apoe_exclusion_sensitivity"); estimates <- rbind(inc_est, exc_est)
  inc_ivw <- inc_est[inc_est$method_id == "mr_ivw", , drop = FALSE]; exc_ivw <- exc_est[exc_est$method_id == "mr_ivw", , drop = FALSE]
  heterogeneity <- rbind(run_heterogeneity(included_mr, "APOE_included", "forward_primary"), run_heterogeneity(excluded_mr, "APOE_excluded", "forward_apoe_exclusion_sensitivity"))
  egger <- rbind(run_egger(included_mr, "APOE_included", "forward_primary"), run_egger(excluded_mr, "APOE_excluded", "forward_apoe_exclusion_sensitivity"))
  inc_loo <- run_loo(included_mr, "APOE_included", "forward_primary", inc_ivw); exc_loo <- run_loo(excluded_mr, "APOE_excluded", "forward_apoe_exclusion_sensitivity", exc_ivw)
  leave_one_out <- rbind(inc_loo$loo, exc_loo$loo); full_ivw <- rbind(inc_loo$full, exc_loo$full)
  inc_presso <- run_presso(included_mr, "APOE_included", "forward_primary"); exc_presso <- run_presso(excluded_mr, "APOE_excluded", "forward_apoe_exclusion_sensitivity"); presso <- rbind(inc_presso, exc_presso)

  five_inc <- nrow(inc_est) == 5L && setequal(inc_est$method_id, required_methods) && all(is.finite(unlist(inc_est[c("beta", "se", "pval")])))
  five_exc <- nrow(exc_est) == 5L && setequal(exc_est$method_id, required_methods) && all(is.finite(unlist(exc_est[c("beta", "se", "pval")])))
  loo_checks <- list(loo_overall_count_included = inc_loo$raw_rows == 313L, loo_overall_count_excluded = exc_loo$raw_rows == 313L, loo_real_snp_count_included = inc_loo$real_snp_rows == 312L, loo_real_snp_count_excluded = exc_loo$real_snp_rows == 312L, loo_real_snp_unique_included = inc_loo$real_snp_unique == 312L, loo_real_snp_unique_excluded = exc_loo$real_snp_unique == 312L, loo_excluded_snp_set_match_included = inc_loo$exact_set, loo_excluded_snp_set_match_excluded = exc_loo$exact_set, loo_frequency_match_included = inc_loo$frequency_match, loo_frequency_match_excluded = exc_loo$frequency_match, loo_no_overall_in_snp_table_included = !any(grepl("^All", inc_loo$loo$excluded_rsid)), loo_no_overall_in_snp_table_excluded = !any(grepl("^All", exc_loo$loo$excluded_rsid)), loo_overall_matches_primary_ivw_included = inc_loo$overall_match, loo_overall_matches_primary_ivw_excluded = exc_loo$overall_match)

  stage <<- "atomic_write_and_readback"
  partial_targets <- setNames(c(mr_input_included_path, mr_input_excluded_path, estimates_path, heterogeneity_path, egger_path, strength_path, strength_summary_path, loo_path, loo_full_path, presso_path, qc_path), c("inc_input", "exc_input", "estimates", "heterogeneity", "egger", "strength", "strength_summary", "loo", "loo_full", "presso", "qc"))
  on.exit(unlink(paste0(unname(partial_targets), ".partial"), force = TRUE), add = TRUE)
  write_partial(function(p) write_parquet(con, included_mr, p), partial_targets[["inc_input"]]); write_partial(function(p) write_parquet(con, excluded_mr, p), partial_targets[["exc_input"]])
  write_partial(function(p) write_csv_precise(estimates, p), partial_targets[["estimates"]]); write_partial(function(p) write_csv_precise(heterogeneity, p), partial_targets[["heterogeneity"]]); write_partial(function(p) write_csv_precise(egger, p), partial_targets[["egger"]]); write_partial(function(p) write_csv_precise(rbind(inc_strength$per_snp, exc_strength$per_snp), p), partial_targets[["strength"]]); write_partial(function(p) write_csv_precise(rbind(inc_strength$summary, exc_strength$summary), p), partial_targets[["strength_summary"]]); write_partial(function(p) write_csv_precise(leave_one_out, p), partial_targets[["loo"]]); write_partial(function(p) write_csv_precise(full_ivw, p), partial_targets[["loo_full"]]); write_partial(function(p) write_csv_precise(presso, p), partial_targets[["presso"]])
  rb <- list(inc_input = nrow(read_parquet(con, paste0(partial_targets[["inc_input"]], ".partial"))) == 312L, exc_input = nrow(read_parquet(con, paste0(partial_targets[["exc_input"]], ".partial"))) == 312L, estimates = nrow(read.csv(paste0(partial_targets[["estimates"]], ".partial"), stringsAsFactors = FALSE, check.names = FALSE)) == 10L, heterogeneity = nrow(read.csv(paste0(partial_targets[["heterogeneity"]], ".partial"), stringsAsFactors = FALSE, check.names = FALSE)) == 4L, egger = nrow(read.csv(paste0(partial_targets[["egger"]], ".partial"), stringsAsFactors = FALSE, check.names = FALSE)) == 2L, strength = nrow(read.csv(paste0(partial_targets[["strength"]], ".partial"), stringsAsFactors = FALSE, check.names = FALSE)) == 624L, loo = nrow(read.csv(paste0(partial_targets[["loo"]], ".partial"), stringsAsFactors = FALSE, check.names = FALSE)) == 624L, loo_full = nrow(read.csv(paste0(partial_targets[["loo_full"]], ".partial"), stringsAsFactors = FALSE, check.names = FALSE)) == 2L, presso = nrow(read.csv(paste0(partial_targets[["presso"]], ".partial"), stringsAsFactors = FALSE, check.names = FALSE)) == 2L)
  output_roundtrip_checks <- all(vapply(rb, isTRUE, logical(1)))
  included_sha_after <- hash_file(included_input_path); excluded_sha_after <- hash_file(excluded_input_path); renv_lock_sha_after <- hash_file(renv_lock)
  no_input_mutation <- identical(included_sha_before, included_sha_after) && identical(excluded_sha_before, excluded_sha_after); renv_lock_unchanged <- identical(tolower(renv_lock_sha_before), tolower(renv_lock_sha_after))
  hard_checks <- c(deps, list(harmonisation_gate = harmonisation_gate, loo_v4_audit_gate = loo_v4_audit_gate, input_sha_matches_v4_audit = input_sha_matches_v4_audit, included_input_count = inc_checks$input_count, excluded_input_count = exc_checks$input_count, included_rsid_unique = inc_checks$rsid_unique, excluded_rsid_unique = exc_checks$rsid_unique, shared_snp_count = membership$shared_snp_count, included_only_count = membership$included_only_count, excluded_only_count = membership$excluded_only_count, all_mr_keep = inc_checks$all_mr_keep && exc_checks$all_mr_keep, no_palindromic_primary = inc_checks$no_palindromic_primary && exc_checks$no_palindromic_primary, no_multiple_reintroduced = !isTRUE(hqc$multiple_target_reintroduced), no_missing_reintroduced = isTRUE(hqc$hard_checks$no_missing_target_reintroduced), aligned_alleles = inc_checks$aligned_alleles && exc_checks$aligned_alleles, finite_exposure_effects = inc_checks$finite_exposure_effects && exc_checks$finite_exposure_effects, positive_exposure_se = inc_checks$positive_exposure_se && exc_checks$positive_exposure_se, finite_outcome_effects = inc_checks$finite_outcome_effects && exc_checks$finite_outcome_effects, positive_outcome_se = inc_checks$positive_outcome_se && exc_checks$positive_outcome_se, valid_outcome_p = inc_checks$valid_outcome_p && exc_checks$valid_outcome_p, required_mr_methods_present = methods_available, instrument_strength_completed_included = nrow(inc_strength$per_snp) == 312L, instrument_strength_completed_excluded = nrow(exc_strength$per_snp) == 312L, five_methods_completed_included = five_inc, five_methods_completed_excluded = five_exc, ivw_completed_included = nrow(inc_ivw) == 1L, ivw_completed_excluded = nrow(exc_ivw) == 1L, heterogeneity_completed_included = nrow(heterogeneity[heterogeneity$analysis_set == "APOE_included", ]) == 2L, heterogeneity_completed_excluded = nrow(heterogeneity[heterogeneity$analysis_set == "APOE_excluded", ]) == 2L, egger_intercept_completed_included = nrow(egger[egger$analysis_set == "APOE_included", ]) == 1L, egger_intercept_completed_excluded = nrow(egger[egger$analysis_set == "APOE_excluded", ]) == 1L, leave_one_out_completed_included = nrow(inc_loo$loo) == 312L, leave_one_out_completed_excluded = nrow(exc_loo$loo) == 312L, mr_presso_attempted_included = TRUE, mr_presso_attempted_excluded = TRUE, output_roundtrip_checks = output_roundtrip_checks, no_input_mutation = no_input_mutation, renv_lock_unchanged = renv_lock_unchanged), loo_checks)
  failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]; status <- if (length(failures) == 0L) "passed" else "failed"
  qc <- list(mr_version = "v3", analysis_direction = "Hb_to_delirium", primary_analysis_set = "APOE_included", sensitivity_analysis_set = "APOE_excluded", source_harmonisation_version = "v4", source_loo_label_audit_version = "v4", source_loo_diagnostic_reason = "v2_loo_set_check_implementation_error", mr_status = status, approved_for_forward_results_interpretation = identical(status, "passed"), effect_interpretation = effect_interpretation, included_nsnp = nrow(included_mr), excluded_nsnp = nrow(excluded_mr), shared_nsnp = shared_nsnp, included_only_nsnp = included_only_nsnp, excluded_only_nsnp = excluded_only_nsnp, primary_ivw = inc_ivw, apoe_excluded_ivw = exc_ivw, instrument_strength_summary = rbind(inc_strength$summary, exc_strength$summary), heterogeneity_summary = heterogeneity, egger_intercept_summary = egger, leave_one_out_summary = list(included = list(total_rows = inc_loo$raw_rows, real_snp_rows = inc_loo$real_snp_rows, overall_rows = 1L), excluded = list(total_rows = exc_loo$raw_rows, real_snp_rows = exc_loo$real_snp_rows, overall_rows = 1L)), mr_presso_summary = presso, loo_overall_matches_primary_ivw_included = inc_loo$overall_match, loo_overall_matches_primary_ivw_excluded = exc_loo$overall_match, steiger_status = "deferred_to_directionality_sensitivity_stage", R_version = R.version.string, TwoSampleMR_version = ts_ver, TwoSampleMR_RemoteSha = ts_sha, MRPRESSO_version = mp_ver, MRPRESSO_RemoteSha = mp_sha, seed = seed_value, mr_library = normalizePath(mr_library, winslash = "/", mustWork = TRUE), harmonised_input_sha256 = list(included_before = included_sha_before, excluded_before = excluded_sha_before, included_after = included_sha_after, excluded_after = excluded_sha_after), renv_lock_sha_before = renv_lock_sha_before, renv_lock_sha_after = renv_lock_sha_after, renv_lock_unchanged = renv_lock_unchanged, hard_checks = hard_checks, hard_check_failures = failures, informational_findings = list(mr_presso_failed_sets = presso$analysis_set[presso$mr_presso_status == "failed"], mr_presso_errors = presso$mr_presso_error[presso$mr_presso_status == "failed"], steiger = "deferred; this is not a hard failure"))
  write_partial(function(p) writeLines(jsonlite::toJSON(qc, auto_unbox = TRUE, pretty = TRUE, digits = NA), p), partial_targets[["qc"]])
  stop_if(status != "passed", paste0("MR hard checks failed: ", paste(failures, collapse = "; ")))
  for (path in unname(partial_targets)) publish_partial(path)
  safe_log("renv_lock_sha_before=", renv_lock_sha_before); safe_log("renv_lock_sha_after=", renv_lock_sha_after); safe_log("mr_status=passed"); safe_log("hard_check_failures="); safe_log("approved_for_forward_results_interpretation=TRUE")
  0L
}

exit_status <- tryCatch(main(), error = function(e) { safe_log("mr_status=failed"); safe_log("stage=", stage); safe_log("error=", conditionMessage(e)); 1L })
quit(status = exit_status)
