args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") stop("Usage: Rscript 09_forward_mr_vuckovic_hb_finngen_r13_loo_label_audit_v4.R --project-root <path>", call. = FALSE)

root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
stage <- "startup"
mr_library <- file.path(root, "renv", "mr-v1-library")
renv_activate <- file.path(root, "renv", "activate.R")
harmonisation_qc_path <- file.path(root, "results", "qc", "vuckovic_hb_finngen_r13_primary_harmonisation_v4.json")
included_path <- file.path(root, "data_derived", "harmonised", "vuckovic_hb_finngen_r13_primary_apoe_included_v4.parquet")
excluded_path <- file.path(root, "data_derived", "harmonised", "vuckovic_hb_finngen_r13_primary_apoe_excluded_v4.parquet")
inc_audit_path <- file.path(root, "results", "qc", "vuckovic_hb_finngen_r13_forward_loo_label_audit_included_v4.csv")
exc_audit_path <- file.path(root, "results", "qc", "vuckovic_hb_finngen_r13_forward_loo_label_audit_excluded_v4.csv")
mismatch_path <- file.path(root, "results", "qc", "vuckovic_hb_finngen_r13_forward_loo_label_mismatches_v4.csv")
strata_path <- file.path(root, "results", "qc", "vuckovic_hb_finngen_r13_forward_loo_strata_v4.csv")
qc_path <- file.path(root, "results", "qc", "vuckovic_hb_finngen_r13_forward_loo_label_audit_v4.json")
log_path <- file.path(root, "results", "logs", "vuckovic_hb_finngen_r13_forward_loo_label_audit_v4.log")
expected_ts_version <- "0.7.9"; expected_ts_sha <- "3d119f20d6fc164b0c7f710f5590fee9580f2c7b"; expected_mp_version <- "1.0"; expected_mp_sha <- "3e3c92d7eda6dce0d1d66077373ec0f7ff4f7e87"

safe_log <- function(...) cat(paste0(...), "\n", file = log_path, append = TRUE)
stop_if <- function(condition, message) if (isTRUE(condition)) stop(message, call. = FALSE)
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
read_parquet <- function(con, path) DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", normalizePath(path, winslash = "/", mustWork = TRUE)))
write_csv_precise <- function(data, path) { old <- options(digits = 17, scipen = 999); on.exit(options(old), add = TRUE); write.csv(data, path, row.names = FALSE, na = "") }
write_partial <- function(writer, final_path) writer(paste0(final_path, ".partial"))
publish_partial <- function(final_path) if (!file.rename(paste0(final_path, ".partial"), final_path)) stop(paste0("Atomic rename failed: ", final_path), call. = FALSE)

summarise_loo_strata <- function(raw) {
  required <- c("id.exposure", "id.outcome", "exposure", "outcome", "SNP")
  if (!all(required %in% names(raw))) stop("Required LOO stratum fields are missing.", call. = FALSE)
  key <- interaction(as.character(raw$id.exposure), as.character(raw$id.outcome), as.character(raw$exposure), as.character(raw$outcome), drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(raw)), key)
  result <- do.call(rbind, lapply(groups, function(idx) data.frame(
    id.exposure = as.character(raw$id.exposure[idx[[1L]]]), id.outcome = as.character(raw$id.outcome[idx[[1L]]]),
    exposure = as.character(raw$exposure[idx[[1L]]]), outcome = as.character(raw$outcome[idx[[1L]]]),
    row_count = length(idx), unique_SNP_count = length(unique(as.character(raw$SNP[idx]))), stringsAsFactors = FALSE
  )))
  rownames(result) <- NULL
  result
}

run_strata_unit_test <- function() {
  fixture <- data.frame(id.exposure = c("e1", "e1", "e1", "e2"), id.outcome = c("o1", "o1", "o1", "o2"), exposure = c("E", "E", "E", "E2"), outcome = c("O", "O", "O", "O2"), SNP = c("rs1", "rs2", "All", "rs3"), stringsAsFactors = FALSE)
  x <- summarise_loo_strata(fixture)
  is.data.frame(x) && identical(names(x), c("id.exposure", "id.outcome", "exposure", "outcome", "row_count", "unique_SNP_count")) && nrow(x) == 2L &&
    x$row_count[[1L]] == 3L && x$unique_SNP_count[[1L]] == 3L && x$row_count[[2L]] == 1L && x$unique_SNP_count[[2L]] == 1L &&
    (is.integer(x$row_count) || is.numeric(x$row_count)) && (is.integer(x$unique_SNP_count) || is.numeric(x$unique_SNP_count)) &&
    !any(vapply(x, is.list, logical(1))) && !any(vapply(x, is.matrix, logical(1)))
}

make_dat <- function(x, set_name) data.frame(
  SNP = as.character(x$rsid), beta.exposure = as.numeric(x$exposure_beta), se.exposure = as.numeric(x$exposure_se), effect_allele.exposure = as.character(x$effect_allele), other_allele.exposure = as.character(x$other_allele), eaf.exposure = as.numeric(x$exposure_eaf),
  beta.outcome = as.numeric(x$outcome_beta), se.outcome = as.numeric(x$outcome_se), effect_allele.outcome = as.character(x$outcome_effect_allele_aligned), other_allele.outcome = as.character(x$outcome_other_allele_aligned), eaf.outcome = as.numeric(x$outcome_eaf),
  exposure = "Vuckovic_2020_Hb", outcome = "FinnGen_R13_F5_DELIRIUM", id.exposure = "vuckovic_hb_2020", id.outcome = "finngen_R13_F5_DELIRIUM", mr_keep = TRUE, analysis_set = set_name, stringsAsFactors = FALSE
)

diagnose_labels <- function(x) data.frame(raw_label = as.character(x), trimmed_label = trimws(as.character(x)), lowercase_label = tolower(as.character(x)), trimmed_lowercase_label = tolower(trimws(as.character(x))), character_count = nchar(as.character(x), type = "chars"), canonical_rsid = grepl("^rs[0-9]+$", as.character(x)), stringsAsFactors = FALSE)

frequency_table <- function(input, loo, analysis_set) {
  labs <- sort(unique(c(input, loo))); a <- table(input); b <- table(loo)
  ia <- as.integer(a[labs]); ib <- as.integer(b[labs]); ia[is.na(ia)] <- 0L; ib[is.na(ib)] <- 0L
  data.frame(analysis_set = analysis_set, label = labs, input_count = ia, loo_count = ib, count_match = ia == ib, stringsAsFactors = FALSE)
}

audit_one <- function(source, dat, analysis_set) {
  raw <- TwoSampleMR::mr_leaveoneout(dat, method = TwoSampleMR::mr_ivw)
  if (!is.data.frame(raw)) raw <- as.data.frame(raw)
  stop_if(!("SNP" %in% names(raw)), "LOO result does not contain SNP.")
  safe_log("LOO analysis_set=", analysis_set); safe_log("LOO class=", paste(class(raw), collapse = ",")); safe_log("LOO nrow=", nrow(raw)); safe_log("LOO ncol=", ncol(raw)); safe_log("LOO names=", paste(names(raw), collapse = "|")); safe_log("LOO SNP head=", paste(head(as.character(raw$SNP), 10), collapse = "|")); safe_log("LOO SNP tail=", paste(tail(as.character(raw$SNP), 10), collapse = "|"))
  labels <- as.character(raw$SNP); overall <- !is.na(labels) & grepl("^All", labels); loo_snps <- labels[!overall]; input_snps <- as.character(dat$SNP); source_rsids <- as.character(source$rsid)
  strata <- summarise_loo_strata(raw); strata$analysis_set <- analysis_set; strata <- strata[, c("analysis_set", "id.exposure", "id.outcome", "exposure", "outcome", "row_count", "unique_SNP_count")]
  audit <- cbind(data.frame(analysis_set = analysis_set, loo_row_number = seq_len(nrow(raw)), exposure = as.character(raw$exposure), outcome = as.character(raw$outcome), id.exposure = as.character(raw$id.exposure), id.outcome = as.character(raw$id.outcome), raw_loo_SNP = labels, is_overall = overall, stringsAsFactors = FALSE), diagnose_labels(labels))
  audit$exactly_in_input_set <- !overall & audit$raw_loo_SNP %in% input_snps
  input_not_in_loo <- setdiff(input_snps, loo_snps); loo_not_in_input <- setdiff(loo_snps, input_snps)
  mismatch <- rbind(if (length(input_not_in_loo)) cbind(data.frame(analysis_set = analysis_set, mismatch_direction = "input_not_in_loo", stringsAsFactors = FALSE), diagnose_labels(input_not_in_loo)) else NULL, if (length(loo_not_in_input)) cbind(data.frame(analysis_set = analysis_set, mismatch_direction = "loo_not_in_input", stringsAsFactors = FALSE), diagnose_labels(loo_not_in_input)) else NULL)
  if (is.null(mismatch)) mismatch <- data.frame(analysis_set = character(), mismatch_direction = character(), raw_label = character(), trimmed_label = character(), lowercase_label = character(), trimmed_lowercase_label = character(), character_count = integer(), canonical_rsid = logical(), stringsAsFactors = FALSE)
  position <- if (length(input_snps) == length(loo_snps)) data.frame(position = seq_along(input_snps), input_SNP = input_snps, loo_SNP = loo_snps, exact_match = input_snps == loo_snps, stringsAsFactors = FALSE) else data.frame(position = integer(), input_SNP = character(), loo_SNP = character(), exact_match = logical(), stringsAsFactors = FALSE)
  list(audit = audit, mismatch = mismatch, strata = strata, frequency = frequency_table(input_snps, loo_snps, analysis_set), summary = list(
    input_count = length(source_rsids), tsmr_pre_loo_count = length(input_snps), source_rsid_missing_count = sum(is.na(source_rsids) | source_rsids == ""), tsmr_snp_missing_count = sum(is.na(input_snps) | input_snps == ""), source_rsid_duplicate_count = sum(duplicated(source_rsids)), tsmr_snp_duplicate_count = sum(duplicated(input_snps)), source_vs_tsmr_exact_set_match = setequal(source_rsids, input_snps), source_vs_tsmr_exact_order_match = identical(source_rsids, input_snps), source_only_rsids = setdiff(source_rsids, input_snps), tsmr_only_snps = setdiff(input_snps, source_rsids), analysis_stratum_count = nrow(strata), loo_total_row_count = nrow(raw), loo_overall_row_count = sum(overall), loo_overall_labels = unique(labels[overall]), loo_non_overall_row_count = length(loo_snps), loo_non_overall_distinct_count = length(unique(loo_snps)), loo_duplicate_count = sum(duplicated(loo_snps)), exact_set_match = setequal(input_snps, loo_snps), input_not_in_loo_count = length(input_not_in_loo), loo_not_in_input_count = length(loo_not_in_input), trimmed_set_match = setequal(trimws(input_snps), trimws(loo_snps)), lowercase_set_match = setequal(tolower(input_snps), tolower(loo_snps)), trimmed_lowercase_set_match = setequal(tolower(trimws(input_snps)), tolower(trimws(loo_snps))), canonical_input_count = sum(diagnose_labels(input_snps)$canonical_rsid), noncanonical_input_count = sum(!diagnose_labels(input_snps)$canonical_rsid), canonical_loo_count = sum(diagnose_labels(loo_snps)$canonical_rsid), noncanonical_loo_count = sum(!diagnose_labels(loo_snps)$canonical_rsid), position_exact_match_count = sum(position$exact_match), input_is_sorted = identical(input_snps, sort(input_snps)), loo_is_sorted = identical(loo_snps, sort(loo_snps)), sorted_vector_exact_match = identical(sort(input_snps), sort(loo_snps)), position_table = position
  ))
}

main <- function() {
  targets <- c(inc_audit_path, exc_audit_path, mismatch_path, strata_path, qc_path, log_path)
  stop_if(any(file.exists(c(targets, paste0(targets, ".partial")))), "LOO label audit V4 target or partial exists; refusing to overwrite.")
  source(renv_activate); .libPaths(c(mr_library, .libPaths()))
  dependency_gate <- requireNamespace("TwoSampleMR", quietly = TRUE) && requireNamespace("MRPRESSO", quietly = TRUE) && identical(as.character(utils::packageVersion("TwoSampleMR")), expected_ts_version) && identical(as.character(utils::packageDescription("TwoSampleMR")$RemoteSha), expected_ts_sha) && identical(as.character(utils::packageVersion("MRPRESSO")), expected_mp_version) && identical(as.character(utils::packageDescription("MRPRESSO")$RemoteSha), expected_mp_sha)
  stop_if(!dependency_gate, "Frozen dependency gate failed.")
  strata_summary_unit_test <- run_strata_unit_test(); stop_if(!strata_summary_unit_test, "Strata summary unit test failed."); safe_log("strata_summary_unit_test=passed")
  hqc <- jsonlite::fromJSON(harmonisation_qc_path, simplifyVector = FALSE); harmonisation_gate <- identical(hqc$harmonisation_status, "passed") && isTRUE(hqc$approved_for_forward_primary_mr) && length(hqc$hard_check_failures) == 0L; stop_if(!harmonisation_gate, "Harmonisation V4 gate failed.")
  stage <<- "loo_label_diagnostic"; inc_sha_before <- hash_file(included_path); exc_sha_before <- hash_file(excluded_path); con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE)); on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  included <- read_parquet(con, included_path); excluded <- read_parquet(con, excluded_path); inc <- audit_one(included, make_dat(included, "APOE_included"), "APOE_included"); exc <- audit_one(excluded, make_dat(excluded, "APOE_excluded"), "APOE_excluded")
  mismatch <- rbind(inc$mismatch, exc$mismatch); strata <- rbind(inc$strata, exc$strata); inc_sha_after <- hash_file(included_path); exc_sha_after <- hash_file(excluded_path); no_input_mutation <- identical(inc_sha_before, inc_sha_after) && identical(exc_sha_before, exc_sha_after)
  hard_checks <- list(dependency_gate = dependency_gate, harmonisation_gate = harmonisation_gate, strata_summary_unit_test = strata_summary_unit_test, included_input_count = nrow(included) == 312L, excluded_input_count = nrow(excluded) == 312L, source_to_tsmr_mapping_audited_included = TRUE, source_to_tsmr_mapping_audited_excluded = TRUE, analysis_strata_audited_included = TRUE, analysis_strata_audited_excluded = TRUE, loo_executed_included = inc$summary$loo_total_row_count > 0L, loo_executed_excluded = exc$summary$loo_total_row_count > 0L, overall_rows_identified_included = inc$summary$loo_overall_row_count >= 0L, overall_rows_identified_excluded = exc$summary$loo_overall_row_count >= 0L, exact_set_difference_computed_included = TRUE, exact_set_difference_computed_excluded = TRUE, frequency_audit_completed_included = nrow(inc$frequency) > 0L, frequency_audit_completed_excluded = nrow(exc$frequency) > 0L, mismatch_rows_fully_exported = nrow(mismatch) == inc$summary$input_not_in_loo_count + inc$summary$loo_not_in_input_count + exc$summary$input_not_in_loo_count + exc$summary$loo_not_in_input_count, no_input_mutation = no_input_mutation)
  failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
  fully_consistent <- inc$summary$source_vs_tsmr_exact_set_match && exc$summary$source_vs_tsmr_exact_set_match && inc$summary$analysis_stratum_count == 1L && exc$summary$analysis_stratum_count == 1L && inc$summary$loo_overall_row_count == 1L && exc$summary$loo_overall_row_count == 1L && inc$summary$loo_non_overall_row_count == 312L && exc$summary$loo_non_overall_row_count == 312L && inc$summary$loo_non_overall_distinct_count == 312L && exc$summary$loo_non_overall_distinct_count == 312L && inc$summary$exact_set_match && exc$summary$exact_set_match
  explicit_reason <- if (fully_consistent) "v2_loo_set_check_implementation_error" else "not_uniquely_confirmed"
  qc <- list(audit_version = "v4", audit_role = "forward_mr_leave_one_out_label_diagnostic", diagnostic_status = if (length(failures)) "failed" else "completed", approved_for_forward_results_interpretation = FALSE, approved_for_loo_mapping_fix = fully_consistent, explicit_diagnostic_reason = explicit_reason, TwoSampleMR_version = as.character(utils::packageVersion("TwoSampleMR")), TwoSampleMR_RemoteSha = as.character(utils::packageDescription("TwoSampleMR")$RemoteSha), included = inc$summary, excluded = exc$summary, frequency_audit = rbind(inc$frequency, exc$frequency), hard_checks = hard_checks, hard_check_failures = failures, input_sha256_before = list(included = inc_sha_before, excluded = exc_sha_before), input_sha256_after = list(included = inc_sha_after, excluded = exc_sha_after))
  stage <<- "atomic_write"; partial_targets <- c(inc_audit_path, exc_audit_path, mismatch_path, strata_path, qc_path); on.exit(unlink(paste0(partial_targets, ".partial"), force = TRUE), add = TRUE)
  write_partial(function(p) write_csv_precise(inc$audit, p), inc_audit_path); write_partial(function(p) write_csv_precise(exc$audit, p), exc_audit_path); write_partial(function(p) write_csv_precise(mismatch, p), mismatch_path); write_partial(function(p) write_csv_precise(strata, p), strata_path); write_partial(function(p) writeLines(jsonlite::toJSON(qc, auto_unbox = TRUE, pretty = TRUE, digits = NA), p), qc_path)
  for (p in partial_targets) publish_partial(p)
  safe_log("diagnostic_status=completed"); safe_log("explicit_diagnostic_reason=", explicit_reason); safe_log("approved_for_loo_mapping_fix=", fully_consistent); safe_log("hard_check_failures=", paste(failures, collapse = ";")); 0L
}

status <- tryCatch(main(), error = function(e) { safe_log("diagnostic_status=failed"); safe_log("stage=", stage); safe_log("error=", conditionMessage(e)); 1L })
quit(status = status)
