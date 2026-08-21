#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

pkgs <- c("digest", "DBI", "duckdb", "jsonlite")
if (any(!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))) {
  stop("Required installed package missing; no automatic installation.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript.exe R/13_reverse_primary_harmonisation_preflight_v1.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}

root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

stop_if <- function(x, m) if (isTRUE(x)) stop(m, call. = FALSE)
hash_file <- function(p) digest::digest(file = p, algo = "sha256", serialize = FALSE)
sql_string <- function(p, must_work = TRUE) paste0("'", gsub("'", "''", normalizePath(p, winslash = "/", mustWork = must_work), fixed = TRUE), "'")
sql_ident <- function(con, x) as.character(DBI::dbQuoteIdentifier(con, x))
comp <- function(x) unname(c(A = "T", T = "A", C = "G", G = "C")[x])
setkey <- function(a, b) paste(sort(c(a, b)), collapse = "/")
records <- function(x) lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))

freeze_json <- file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_outcome_extraction_v1_freeze.json")
unique_pq <- file.path(root, "data_derived", "reverse_outcome_extraction", "vuckovic_hb_reverse_primary_unique_matches_v1.parquet")
source_path <- file.path(root, "data_raw", "gwas", "ebi-a-GCST90002384.vcf.gz")

out <- c(
  pq = file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_primary_harmonisation_preflight_v1.parquet"),
  tsv = file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_primary_harmonisation_preflight_v1.tsv"),
  json = file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_harmonisation_preflight_v1.json"),
  log = file.path(root, "results", "logs", "vuckovic_hb_reverse_primary_harmonisation_preflight_v1.log")
)
partials <- paste0(out, ".partial")

stop_if(any(!file.exists(c(freeze_json, unique_pq, source_path))), "Required freeze or extraction input is missing.")
stop_if(any(file.exists(c(out, partials))), "A preflight final or partial target exists; refusing to overwrite.")
dir.create(dirname(out[["pq"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["json"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["log"]]), recursive = TRUE, showWarnings = FALSE)
log_line <- function(x) cat(sprintf("%s %s\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), x), file = out[["log"]], append = TRUE)

main <- function() {
  log_line("stage=reverse_primary_harmonisation_preflight_v1")
  freeze <- jsonlite::fromJSON(freeze_json, simplifyVector = FALSE)
  stop_if(!identical(freeze$freeze_status, "passed"), "Outcome extraction freeze gate is not passed.")
  stop_if(!isTRUE(freeze$approved_for_reverse_primary_harmonisation_preflight), "Freeze did not approve harmonisation preflight.")

  source_sha_before <- hash_file(source_path)
  stop_if(!identical(source_sha_before, freeze$vuckovic_source_sha256), "Vuckovic source SHA differs from extraction freeze.")

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  x <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(unique_pq)))
  stop_if(nrow(x) != 2L, "Expected exactly two unique extraction rows.")

  x <- x[order(match(x$membership, c("apoe_included", "apoe_excluded"))), , drop = FALSE]
  valid_exp <- x$exposure_effect_allele_raw %in% c("A", "C", "G", "T") & x$exposure_other_allele_raw %in% c("A", "C", "G", "T")
  valid_out <- x$outcome_effect_allele_raw %in% c("A", "C", "G", "T") & x$outcome_other_allele_raw %in% c("A", "C", "G", "T")
  exact <- valid_exp & valid_out & x$exposure_effect_allele_raw == x$outcome_effect_allele_raw & x$exposure_other_allele_raw == x$outcome_other_allele_raw
  swapped <- valid_exp & valid_out & x$exposure_effect_allele_raw == x$outcome_other_allele_raw & x$exposure_other_allele_raw == x$outcome_effect_allele_raw
  strand_exact <- valid_exp & valid_out & x$exposure_effect_allele_raw == comp(x$outcome_effect_allele_raw) & x$exposure_other_allele_raw == comp(x$outcome_other_allele_raw)
  strand_swapped <- valid_exp & valid_out & x$exposure_effect_allele_raw == comp(x$outcome_other_allele_raw) & x$exposure_other_allele_raw == comp(x$outcome_effect_allele_raw)

  preflight <- data.frame(
    analysis_set = ifelse(x$membership == "apoe_included", "APOE included", "APOE excluded"),
    target_rsid = x$target_rsid,
    exposure_effect_allele = x$exposure_effect_allele_raw,
    exposure_other_allele = x$exposure_other_allele_raw,
    exposure_beta = x$exposure_beta_raw,
    exposure_se = x$exposure_se_raw,
    exposure_pval = x$exposure_pval_raw,
    exposure_eaf = x$exposure_eaf_raw,
    outcome_effect_allele_raw = x$outcome_effect_allele_raw,
    outcome_other_allele_raw = x$outcome_other_allele_raw,
    outcome_beta_raw = x$outcome_beta_raw,
    outcome_se_raw = x$outcome_se_raw,
    outcome_pval_raw = x$outcome_pval_raw,
    outcome_eaf_raw = x$outcome_eaf_raw,
    valid_exposure_alleles = valid_exp,
    valid_outcome_alleles = valid_out,
    palindromic_snp = valid_exp & mapply(setkey, x$exposure_effect_allele_raw, x$exposure_other_allele_raw) %in% c("A/T", "C/G"),
    raw_orientation_class = ifelse(!valid_exp | !valid_out, "invalid",
      ifelse(exact, "exact_match",
      ifelse(swapped, "swapped_match",
      ifelse(strand_exact, "strand_exact_match",
      ifelse(strand_swapped, "strand_swapped_match", "incompatible"))))),
    beta_flip_required = swapped | strand_swapped,
    eaf_flip_required = swapped | strand_swapped,
    strand_flip_required = strand_exact | strand_swapped,
    beta_flip_performed = FALSE,
    eaf_flip_performed = FALSE,
    harmonisation_performed = FALSE,
    source_match_unique = x$match_class == "unique",
    stringsAsFactors = FALSE
  )
  preflight$eligible_for_formal_harmonisation <- preflight$source_match_unique &
    preflight$valid_exposure_alleles & preflight$valid_outcome_alleles &
    preflight$raw_orientation_class %in% c("exact_match", "swapped_match", "strand_exact_match", "strand_swapped_match") &
    !preflight$raw_orientation_class %in% c("invalid", "incompatible")

  hard_check_failures <- character()
  add_fail <- function(x) hard_check_failures <<- unique(c(hard_check_failures, x))
  if (!identical(freeze$freeze_status, "passed")) add_fail("outcome_extraction_freeze_gate_failed")
  if (!any(preflight$analysis_set == "APOE included" & preflight$target_rsid == freeze$included_target_rsid)) add_fail("included_target_present_failed")
  if (!any(preflight$analysis_set == "APOE excluded" & preflight$target_rsid == freeze$excluded_target_rsid)) add_fail("excluded_target_present_failed")
  if (!all(preflight$source_match_unique)) add_fail("all_source_matches_unique_failed")
  if (!all(preflight$valid_exposure_alleles)) add_fail("exposure_alleles_valid_failed")
  if (!all(preflight$valid_outcome_alleles)) add_fail("outcome_alleles_valid_failed")
  if (any(is.na(preflight$raw_orientation_class) | preflight$raw_orientation_class == "")) add_fail("orientation_classified_failed")
  if (any(is.na(preflight$palindromic_snp))) add_fail("palindromic_status_classified_failed")
  if (any(preflight$beta_flip_performed)) add_fail("beta_flip_performed")
  if (any(preflight$eaf_flip_performed)) add_fail("eaf_flip_performed")
  if (any(preflight$harmonisation_performed)) add_fail("harmonisation_performed")
  if (sum(preflight$analysis_set == "APOE included") != 1L) add_fail("included_preflight_complete_failed")
  if (sum(preflight$analysis_set == "APOE excluded") != 1L) add_fail("excluded_preflight_complete_failed")
  if (!all(preflight$target_rsid %in% c(freeze$included_target_rsid, freeze$excluded_target_rsid))) add_fail("target_identity_failed")

  write.table(preflight, paste0(out[["tsv"]], ".partial"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  tmp <- "reverse_primary_harmonisation_preflight_v1"
  DBI::dbWriteTable(con, tmp, preflight, temporary = TRUE, overwrite = TRUE)
  DBI::dbExecute(con, sprintf("COPY %s TO %s (FORMAT PARQUET, COMPRESSION ZSTD)", sql_ident(con, tmp), sql_string(paste0(out[["pq"]], ".partial"), must_work = FALSE)))
  staged_tsv <- read.delim(paste0(out[["tsv"]], ".partial"), check.names = FALSE)
  staged_pq <- DBI::dbGetQuery(con, sprintf("SELECT target_rsid, analysis_set FROM read_parquet(%s)", sql_string(paste0(out[["pq"]], ".partial"))))
  parquet_tsv_consistency <- nrow(staged_tsv) == nrow(staged_pq) && setequal(staged_tsv$target_rsid, staged_pq$target_rsid)
  if (!parquet_tsv_consistency) add_fail("parquet_tsv_consistency_failed")

  source_sha_after <- hash_file(source_path)
  if (!identical(source_sha_after, source_sha_before)) add_fail("no_source_mutation_failed")

  hard_checks <- list(
    outcome_extraction_freeze_gate = identical(freeze$freeze_status, "passed"),
    included_target_present = any(preflight$analysis_set == "APOE included" & preflight$target_rsid == freeze$included_target_rsid),
    excluded_target_present = any(preflight$analysis_set == "APOE excluded" & preflight$target_rsid == freeze$excluded_target_rsid),
    all_source_matches_unique = all(preflight$source_match_unique),
    exposure_alleles_valid = all(preflight$valid_exposure_alleles),
    outcome_alleles_valid = all(preflight$valid_outcome_alleles),
    orientation_classified = all(!is.na(preflight$raw_orientation_class) & preflight$raw_orientation_class != ""),
    palindromic_status_classified = all(!is.na(preflight$palindromic_snp)),
    no_beta_flip_performed = !any(preflight$beta_flip_performed),
    no_eaf_flip_performed = !any(preflight$eaf_flip_performed),
    no_harmonisation_performed = !any(preflight$harmonisation_performed),
    included_preflight_complete = sum(preflight$analysis_set == "APOE included") == 1L,
    excluded_preflight_complete = sum(preflight$analysis_set == "APOE excluded") == 1L,
    parquet_tsv_consistency = parquet_tsv_consistency,
    no_source_mutation = identical(source_sha_after, source_sha_before),
    no_reverse_mr = TRUE,
    relaxed_threshold_not_started = TRUE
  )
  preflight_status <- if (length(hard_check_failures) == 0L &&
    all(preflight$eligible_for_formal_harmonisation) &&
    all(!preflight$palindromic_snp) &&
    all(preflight$raw_orientation_class %in% c("exact_match", "swapped_match", "strand_exact_match", "strand_swapped_match"))) "passed" else "failed"
  approved <- identical(preflight_status, "passed")

  count_class <- function(z) sum(preflight$raw_orientation_class == z)
  qc <- list(
    preflight_version = "v1",
    analysis_direction = "delirium_to_Hb",
    analysis_role = "secondary_reverse_primary",
    source_outcome_extraction_version = "v1",
    included_preflight_count = sum(preflight$analysis_set == "APOE included"),
    excluded_preflight_count = sum(preflight$analysis_set == "APOE excluded"),
    included_eligible_for_formal_harmonisation_count = sum(preflight$analysis_set == "APOE included" & preflight$eligible_for_formal_harmonisation),
    excluded_eligible_for_formal_harmonisation_count = sum(preflight$analysis_set == "APOE excluded" & preflight$eligible_for_formal_harmonisation),
    palindromic_count = sum(preflight$palindromic_snp),
    exact_match_count = count_class("exact_match"),
    swapped_match_count = count_class("swapped_match"),
    strand_exact_match_count = count_class("strand_exact_match"),
    strand_swapped_match_count = count_class("strand_swapped_match"),
    incompatible_count = count_class("incompatible"),
    invalid_count = count_class("invalid"),
    beta_flip_required_count = sum(preflight$beta_flip_required),
    eaf_flip_required_count = sum(preflight$eaf_flip_required),
    beta_flip_performed_count = sum(preflight$beta_flip_performed),
    eaf_flip_performed_count = sum(preflight$eaf_flip_performed),
    harmonisation_performed = FALSE,
    relaxed_threshold = 5e-6,
    relaxed_threshold_status = "not_triggered",
    source_sha = list(before = source_sha_before, after = source_sha_after),
    preflight_rows = records(preflight),
    preflight_status = preflight_status,
    approved_for_reverse_primary_formal_harmonisation = approved,
    hard_checks = hard_checks,
    hard_check_failures = hard_check_failures
  )
  jsonlite::write_json(qc, paste0(out[["json"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")

  for (p in out[c("pq", "tsv", "json")]) {
    stop_if(file.exists(p), paste("Output appeared during run:", p))
    stop_if(!file.rename(paste0(p, ".partial"), p), paste("Atomic rename failed:", p))
  }
  log_line(sprintf("SUCCESS preflight_status=%s approved_for_reverse_primary_formal_harmonisation=%s hard_check_failures=%s", preflight_status, approved, paste(hard_check_failures, collapse = ";")))
  if (!approved) stop("Preflight did not pass; formal harmonisation requires audit.", call. = FALSE)
}

tryCatch(main(), error = function(e) {
  log_line(paste0("TERMINATED_PRIMARY: ", conditionMessage(e)))
  quit(status = 1L)
})

