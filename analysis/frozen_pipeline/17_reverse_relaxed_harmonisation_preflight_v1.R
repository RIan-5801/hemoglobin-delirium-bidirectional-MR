#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/17_reverse_relaxed_harmonisation_preflight_v1.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

for (pkg in c("DBI", "duckdb", "jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
hash_file <- function(path) digest::digest(file = path, algo = "sha256")
sql_string <- function(path, must_work = TRUE) paste0("'", gsub("'", "''", normalizePath(path, winslash = "/", mustWork = must_work), fixed = TRUE), "'")
sql_ident <- function(con, x) as.character(DBI::dbQuoteIdentifier(con, x))
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
comp <- function(x) unname(c(A = "T", T = "A", C = "G", G = "C")[x])
setkey <- function(a, b) paste(sort(c(a, b)), collapse = "/")
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))

out <- c(
  pq = file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1.parquet"),
  tsv = file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1.tsv"),
  json = file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1.json"),
  log = file.path(root, "results", "logs", "vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1.log")
)
stop_if(any(file.exists(c(out, paste0(out, ".partial")))), "A relaxed preflight final or partial target exists; refusing to overwrite.")
dir.create(dirname(out[["pq"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["json"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["log"]]), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = out[["log"]], append = TRUE)

audit_palindrome_rule_scope <- function() {
  forward_qc_path <- file.path(root, "results", "qc", "vuckovic_hb_finngen_r13_primary_harmonisation_v4.json")
  forward_decision_path <- file.path(root, "docs", "decisions", "30_vuckovic_hb_finngen_r13_primary_harmonisation_v4_v1.1.md")
  forward_rule <- NA_character_
  if (file.exists(forward_qc_path)) {
    fq <- jsonlite::fromJSON(forward_qc_path, simplifyVector = FALSE)
    forward_rule <- fq$informational_findings$palindromic_rule
  }
  decision_files <- list.files(file.path(root, "docs", "decisions"), pattern = "\\.md$", full.names = TRUE)
  protocol_files <- list.files(file.path(root, "docs", "protocol"), pattern = "\\.md$", full.names = TRUE)
  all_files <- c(decision_files, protocol_files)
  txt <- unlist(lapply(all_files, function(p) paste(readLines(p, warn = FALSE, encoding = "UTF-8"), collapse = "\n")), use.names = FALSE)
  project_wide_hit <- any(grepl("palindrom", txt, ignore.case = TRUE) & grepl("project-wide|bidirectional", txt, ignore.case = TRUE))
  list(
    project_wide_or_bidirectional_palindrome_rule_found = project_wide_hit,
    project_wide_rule_applied = FALSE,
    forward_specific_rule_found = !is.na(forward_rule) && nzchar(forward_rule),
    forward_specific_rule_decision = if (file.exists(forward_decision_path)) 30L else NA_integer_,
    forward_specific_rule_text = forward_rule,
    forward_specific_rule_not_automatically_applicable_to_reverse_relaxed_branch = !project_wide_hit && !is.na(forward_rule) && nzchar(forward_rule)
  )
}

main <- function() {
  log_line("stage=reverse_relaxed_harmonisation_preflight_v1")
  freeze_path <- file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_outcome_extraction_v1_freeze.json")
  unique_path <- file.path(root, "data_derived", "reverse_outcome_extraction", "vuckovic_hb_reverse_relaxed_unique_matches_v1.parquet")
  source_path <- file.path(root, "data_raw", "gwas", "ebi-a-GCST90002384.vcf.gz")
  renv_lock <- file.path(root, "renv.lock")
  decision_path <- file.path(root, "docs", "decisions", "62_vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1_v1.1.md")
  stop_if(any(!file.exists(c(freeze_path, unique_path, source_path, renv_lock, decision_path))), "Required preflight input is missing.")
  freeze <- jsonlite::fromJSON(freeze_path, simplifyVector = FALSE)
  stop_if(!identical(freeze$freeze_status, "passed"), "Relaxed outcome extraction freeze gate failed.")
  stop_if(!isTRUE(freeze$approved_for_reverse_relaxed_harmonisation_preflight), "Freeze did not approve relaxed harmonisation preflight.")
  stop_if(length(freeze$hard_check_failures) != 0L, "Freeze hard_check_failures is not empty.")
  source_sha_before <- hash_file(source_path)
  renv_before <- hash_file(renv_lock)

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  x <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(unique_path)))
  stop_if(nrow(x) != as.integer(freeze$unique_match_count), "Unique-match row count differs from freeze.")
  valid_exp <- x$exposure_effect_allele %in% c("A", "C", "G", "T") & x$exposure_other_allele %in% c("A", "C", "G", "T")
  valid_out <- x$outcome_effect_allele_raw %in% c("A", "C", "G", "T") & x$outcome_other_allele_raw %in% c("A", "C", "G", "T")
  exact <- valid_exp & valid_out & x$exposure_effect_allele == x$outcome_effect_allele_raw & x$exposure_other_allele == x$outcome_other_allele_raw
  swapped <- valid_exp & valid_out & x$exposure_effect_allele == x$outcome_other_allele_raw & x$exposure_other_allele == x$outcome_effect_allele_raw
  strand_exact <- valid_exp & valid_out & x$exposure_effect_allele == comp(x$outcome_effect_allele_raw) & x$exposure_other_allele == comp(x$outcome_other_allele_raw)
  strand_swapped <- valid_exp & valid_out & x$exposure_effect_allele == comp(x$outcome_other_allele_raw) & x$exposure_other_allele == comp(x$outcome_effect_allele_raw)
  orientation <- ifelse(!valid_exp | !valid_out, "invalid",
    ifelse(exact, "exact_match",
    ifelse(swapped, "swapped_match",
    ifelse(strand_exact, "strand_exact_match",
    ifelse(strand_swapped, "strand_swapped_match", "incompatible")))))
  pal <- valid_exp & mapply(setkey, x$exposure_effect_allele, x$exposure_other_allele) %in% c("A/T", "C/G")
  resolvable <- orientation %in% c("exact_match", "swapped_match", "strand_exact_match", "strand_swapped_match")
  preflight <- data.frame(
    target_rsid = x$target_rsid,
    included_member = as.logical(x$included_member),
    excluded_member = as.logical(x$excluded_member),
    exposure_effect_allele = x$exposure_effect_allele,
    exposure_other_allele = x$exposure_other_allele,
    exposure_beta = as.numeric(x$exposure_beta),
    exposure_se = as.numeric(x$exposure_se),
    exposure_pval = as.numeric(x$exposure_pval),
    exposure_eaf = as.numeric(x$exposure_eaf),
    exposure_chr_grch38 = as.integer(x$exposure_chr_grch38),
    exposure_pos_grch38 = as.integer(x$exposure_pos_grch38),
    outcome_effect_allele_raw = x$outcome_effect_allele_raw,
    outcome_other_allele_raw = x$outcome_other_allele_raw,
    outcome_beta_raw = as.numeric(x$outcome_beta_raw),
    outcome_se_raw = as.numeric(x$outcome_se_raw),
    outcome_pval_raw = as.numeric(x$outcome_pval_raw),
    outcome_eaf_raw = as.numeric(x$outcome_eaf_raw),
    outcome_chr_grch37 = x$outcome_chr_raw,
    outcome_pos_grch37 = x$outcome_pos_raw,
    source_match_count = 1L,
    source_match_unique = TRUE,
    valid_exposure_alleles = valid_exp,
    valid_outcome_alleles = valid_out,
    palindromic_snp = pal,
    raw_orientation_class = orientation,
    beta_flip_required = (swapped | strand_swapped) & !pal,
    eaf_flip_required = (swapped | strand_swapped) & !pal & !is.na(x$outcome_eaf_raw),
    strand_flip_required = (strand_exact | strand_swapped) & !pal,
    beta_flip_performed = FALSE,
    eaf_flip_performed = FALSE,
    strand_flip_performed = FALSE,
    requires_palindromic_rule_adjudication = pal,
    eligible_for_formal_harmonisation_without_rule = valid_exp & valid_out & resolvable & !pal,
    potentially_harmonisable_after_rule_adjudication = valid_exp & valid_out & resolvable & pal,
    harmonisation_performed = FALSE,
    stringsAsFactors = FALSE
  )

  rule_audit <- audit_palindrome_rule_scope()
  write.table(preflight, paste0(out[["tsv"]], ".partial"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  nm <- "reverse_relaxed_harmonisation_preflight_v1"
  DBI::dbWriteTable(con, nm, preflight, temporary = TRUE, overwrite = TRUE)
  DBI::dbExecute(con, sprintf("COPY %s TO %s (FORMAT PARQUET, COMPRESSION ZSTD)", sql_ident(con, nm), sql_string(paste0(out[["pq"]], ".partial"), must_work = FALSE)))
  staged_tsv <- read.delim(paste0(out[["tsv"]], ".partial"), check.names = FALSE)
  staged_pq <- DBI::dbGetQuery(con, sprintf("SELECT target_rsid FROM read_parquet(%s)", sql_string(paste0(out[["pq"]], ".partial"))))
  parquet_tsv_consistency <- nrow(staged_tsv) == nrow(staged_pq) && setequal(staged_tsv$target_rsid, staged_pq$target_rsid) && !anyDuplicated(names(staged_tsv))
  source_sha_after <- hash_file(source_path)
  renv_after <- hash_file(renv_lock)

  hard_check_failures <- character()
  if (!identical(freeze$freeze_status, "passed")) hard_check_failures <- c(hard_check_failures, "relaxed_outcome_extraction_freeze_gate_failed")
  if (nrow(preflight) != as.integer(freeze$union_target_count)) hard_check_failures <- c(hard_check_failures, "all_targets_present_failed")
  if (!all(preflight$source_match_unique)) hard_check_failures <- c(hard_check_failures, "all_matches_unique_failed")
  if (!all(preflight$valid_exposure_alleles)) hard_check_failures <- c(hard_check_failures, "exposure_alleles_valid_failed")
  if (!all(preflight$valid_outcome_alleles)) hard_check_failures <- c(hard_check_failures, "outcome_alleles_valid_failed")
  if (any(is.na(preflight$raw_orientation_class) | preflight$raw_orientation_class == "")) hard_check_failures <- c(hard_check_failures, "orientation_recomputed_failed")
  if (any(is.na(preflight$palindromic_snp))) hard_check_failures <- c(hard_check_failures, "palindromic_status_recomputed_failed")
  if (rule_audit$project_wide_rule_applied) hard_check_failures <- c(hard_check_failures, "palindromic_decision_applied")
  if (any(preflight$beta_flip_performed)) hard_check_failures <- c(hard_check_failures, "beta_flip_performed")
  if (any(preflight$eaf_flip_performed)) hard_check_failures <- c(hard_check_failures, "eaf_flip_performed")
  if (any(preflight$strand_flip_performed)) hard_check_failures <- c(hard_check_failures, "strand_flip_performed")
  if (any(preflight$harmonisation_performed)) hard_check_failures <- c(hard_check_failures, "harmonisation_performed")
  if (sum(preflight$included_member) != as.integer(freeze$included_target_count)) hard_check_failures <- c(hard_check_failures, "included_preflight_complete_failed")
  if (sum(preflight$excluded_member) != as.integer(freeze$excluded_target_count)) hard_check_failures <- c(hard_check_failures, "excluded_preflight_complete_failed")
  if (!parquet_tsv_consistency) hard_check_failures <- c(hard_check_failures, "parquet_tsv_consistency_failed")
  if (!identical(source_sha_before, source_sha_after)) hard_check_failures <- c(hard_check_failures, "source_mutation_detected")
  if (!identical(renv_before, renv_after)) hard_check_failures <- c(hard_check_failures, "renv_lock_changed")

  count_class <- function(z) sum(preflight$raw_orientation_class == z)
  requires_count <- sum(preflight$requires_palindromic_rule_adjudication)
  preflight_status <- if (length(hard_check_failures) == 0L &&
    all(preflight$source_match_unique) &&
    all(preflight$valid_exposure_alleles) &&
    all(preflight$valid_outcome_alleles) &&
    count_class("incompatible") == 0L &&
    count_class("invalid") == 0L) "passed" else "failed"
  pal_rule_required <- requires_count > 0L
  approved_formal <- identical(preflight_status, "passed") && !pal_rule_required && isTRUE(rule_audit$project_wide_or_bidirectional_palindrome_rule_found)
  hard_checks <- list(
    relaxed_outcome_extraction_freeze_gate = identical(freeze$freeze_status, "passed"),
    all_targets_present = nrow(preflight) == as.integer(freeze$union_target_count),
    all_matches_unique = all(preflight$source_match_unique),
    exposure_alleles_valid = all(preflight$valid_exposure_alleles),
    outcome_alleles_valid = all(preflight$valid_outcome_alleles),
    orientation_recomputed = all(!is.na(preflight$raw_orientation_class) & preflight$raw_orientation_class != ""),
    palindromic_status_recomputed = all(!is.na(preflight$palindromic_snp)),
    palindromic_rule_scope_audited = TRUE,
    no_palindromic_decision_applied = !rule_audit$project_wide_rule_applied,
    no_beta_flip_performed = !any(preflight$beta_flip_performed),
    no_eaf_flip_performed = !any(preflight$eaf_flip_performed),
    no_strand_flip_performed = !any(preflight$strand_flip_performed),
    no_harmonisation_performed = !any(preflight$harmonisation_performed),
    included_preflight_complete = sum(preflight$included_member) == as.integer(freeze$included_target_count),
    excluded_preflight_complete = sum(preflight$excluded_member) == as.integer(freeze$excluded_target_count),
    parquet_tsv_consistency = parquet_tsv_consistency,
    no_source_mutation = identical(source_sha_before, source_sha_after),
    no_reverse_mr = TRUE
  )
  qc <- list(
    preflight_version = "v1",
    decision = 62,
    analysis_direction = "delirium_to_Hb",
    analysis_role = "secondary_reverse_exploratory_relaxed",
    source_instrument_version = "v2",
    source_outcome_extraction_version = "v1",
    source_outcome_extraction_freeze_version = "v1",
    p_threshold = 5e-6,
    union_preflight_count = nrow(preflight),
    included_preflight_count = sum(preflight$included_member),
    excluded_preflight_count = sum(preflight$excluded_member),
    palindromic_count = sum(preflight$palindromic_snp),
    requires_palindromic_rule_adjudication_count = requires_count,
    included_without_rule_eligible_count = sum(preflight$included_member & preflight$eligible_for_formal_harmonisation_without_rule),
    excluded_without_rule_eligible_count = sum(preflight$excluded_member & preflight$eligible_for_formal_harmonisation_without_rule),
    included_potential_after_rule_count = sum(preflight$included_member & (preflight$eligible_for_formal_harmonisation_without_rule | preflight$potentially_harmonisable_after_rule_adjudication)),
    excluded_potential_after_rule_count = sum(preflight$excluded_member & (preflight$eligible_for_formal_harmonisation_without_rule | preflight$potentially_harmonisable_after_rule_adjudication)),
    exact_match_count = count_class("exact_match"),
    swapped_match_count = count_class("swapped_match"),
    strand_exact_match_count = count_class("strand_exact_match"),
    strand_swapped_match_count = count_class("strand_swapped_match"),
    incompatible_count = count_class("incompatible"),
    invalid_count = count_class("invalid"),
    beta_flip_required_count = sum(preflight$beta_flip_required),
    eaf_flip_required_count = sum(preflight$eaf_flip_required),
    strand_flip_required_count = sum(preflight$strand_flip_required),
    beta_flip_performed_count = sum(preflight$beta_flip_performed),
    eaf_flip_performed_count = sum(preflight$eaf_flip_performed),
    strand_flip_performed_count = sum(preflight$strand_flip_performed),
    harmonisation_performed = FALSE,
    palindromic_rule_adjudication_required = pal_rule_required,
    preflight_status = preflight_status,
    approved_for_reverse_relaxed_formal_harmonisation = approved_formal,
    preflight_rows = records(preflight),
    palindrome_rule_scope_audit = rule_audit,
    source_sha = list(before = source_sha_before, after = source_sha_after),
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after,
    hard_checks = hard_checks,
    hard_check_failures = hard_check_failures,
    informational_findings = list(
      formal_harmonisation_not_performed = TRUE,
      palindromic_exact_match_still_requires_rule = pal_rule_required,
      forward_specific_rule_not_applied = isTRUE(rule_audit$forward_specific_rule_not_automatically_applicable_to_reverse_relaxed_branch)
    )
  )
  jsonlite::write_json(qc, paste0(out[["json"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  if (!identical(preflight_status, "passed")) stop("Relaxed harmonisation preflight failed; partial outputs retained.", call. = FALSE)
  for (path in out[c("pq", "tsv", "json")]) {
    stop_if(file.exists(path), paste("Output appeared during run:", path))
    stop_if(!file.rename(paste0(path, ".partial"), path), paste("Atomic rename failed:", path))
  }
  log_line("preflight_status=passed palindromic_rule_adjudication_required=", pal_rule_required, " approved_for_reverse_relaxed_formal_harmonisation=", approved_formal)
}

tryCatch(main(), error = function(e) {
  log_line("status=failed error=", conditionMessage(e))
  quit(status = 1L)
})
