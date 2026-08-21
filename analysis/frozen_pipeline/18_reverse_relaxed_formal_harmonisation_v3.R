#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/18_reverse_relaxed_formal_harmonisation_v3.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
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
num_equal <- function(a, b, tol = 1e-12) {
  a <- as.numeric(a); b <- as.numeric(b)
  (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tol)
}
storage_type <- function(x) {
  if (is.logical(x)) return("logical")
  if (is.numeric(x)) return("numeric")
  if (is.character(x)) return("character")
  paste(class(x), collapse = "/")
}
schema_types <- function(x) vapply(x, storage_type, character(1))
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))

character_cols_all <- c(
  "target_rsid", "analysis_set", "exposure_effect_allele", "exposure_other_allele",
  "outcome_effect_allele_raw", "outcome_other_allele_raw",
  "outcome_effect_allele_harmonised", "outcome_other_allele_harmonised",
  "orientation_class", "exclusion_reason"
)
logical_cols_all <- c(
  "included_member", "excluded_member", "palindromic_snp", "beta_flipped",
  "eaf_flipped", "strand_flipped", "record_excluded", "final_valid_instrument"
)
numeric_cols_all <- c(
  "exposure_beta", "exposure_se", "exposure_pval", "exposure_eaf",
  "outcome_beta_raw", "outcome_se_raw", "outcome_pval_raw", "outcome_eaf_raw",
  "outcome_beta_harmonised", "outcome_se_harmonised", "outcome_pval_harmonised",
  "outcome_eaf_harmonised", "exposure_chr_grch38", "exposure_pos_grch38",
  "outcome_chr_grch37", "outcome_pos_grch37"
)
nullable_character_columns <- c("outcome_effect_allele_harmonised", "outcome_other_allele_harmonised", "exclusion_reason")
required_character_columns <- c(
  "target_rsid", "exposure_effect_allele", "exposure_other_allele",
  "outcome_effect_allele_raw", "outcome_other_allele_raw", "orientation_class"
)

expected_schema <- function(cols) {
  schema <- rep(NA_character_, length(cols)); names(schema) <- cols
  schema[cols %in% character_cols_all] <- "character"
  schema[cols %in% logical_cols_all] <- "logical"
  schema[cols %in% numeric_cols_all] <- "numeric"
  if (any(is.na(schema))) stop("No explicit TSV schema for column(s): ", paste(cols[is.na(schema)], collapse = ", "), call. = FALSE)
  schema
}
validate_required_characters <- function(x) {
  req <- intersect(required_character_columns, names(x))
  all(vapply(req, function(col) !any(is.na(x[[col]]) | x[[col]] == ""), logical(1)))
}
normalize_nullable_characters <- function(x, nullable_cols = nullable_character_columns) {
  counts <- integer(0)
  for (col in intersect(nullable_cols, names(x))) {
    blank <- !is.na(x[[col]]) & x[[col]] == ""
    counts[col] <- sum(blank)
    x[[col]][blank] <- NA_character_
  }
  attr(x, "blank_to_na_counts") <- counts
  x
}
read_tsv_explicit_normalized <- function(path, cols) {
  schema <- expected_schema(cols)
  raw <- read.delim(path, check.names = FALSE, colClasses = unname(schema), na.strings = "NA")
  raw_types <- schema_types(raw)
  normalized <- normalize_nullable_characters(raw)
  list(raw = raw, normalized = normalized, raw_types = raw_types, normalized_types = schema_types(normalized),
       blank_to_na_counts = attr(normalized, "blank_to_na_counts"))
}
column_report <- function(pq, raw, normalized, tol = 1e-12) {
  schema <- expected_schema(names(pq))
  counts <- attr(normalized, "blank_to_na_counts")
  out <- lapply(names(pq), function(col) {
    value_ok <- if (schema[[col]] == "numeric") {
      all(num_equal(pq[[col]], normalized[[col]], tol))
    } else {
      idx <- !is.na(pq[[col]]) & !is.na(normalized[[col]])
      identical(as.character(pq[[col]][idx]), as.character(normalized[[col]][idx]))
    }
    data.frame(
      column_name = col,
      expected_type = schema[[col]],
      parquet_type = storage_type(pq[[col]]),
      tsv_raw_type = storage_type(raw[[col]]),
      tsv_normalized_type = storage_type(normalized[[col]]),
      nullable_character = col %in% nullable_character_columns,
      blank_to_na_conversion_count = if (col %in% names(counts)) as.integer(counts[[col]]) else 0L,
      value_consistency = value_ok,
      na_pattern_consistency = identical(is.na(pq[[col]]), is.na(normalized[[col]])),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}
table_equal_after_normalization <- function(pq, raw, normalized, tol = 1e-12) {
  if (!identical(names(pq), names(normalized)) || nrow(pq) != nrow(normalized) || anyDuplicated(names(pq))) return(FALSE)
  schema <- expected_schema(names(pq))
  if (!identical(unname(schema_types(pq)), unname(schema))) return(FALSE)
  if (!identical(unname(schema_types(normalized)), unname(schema))) return(FALSE)
  if (!validate_required_characters(normalized)) return(FALSE)
  report <- column_report(pq, raw, normalized, tol)
  all(report$value_consistency & report$na_pattern_consistency)
}
write_pq <- function(con, x, path) {
  nm <- paste0("tmp_", digest::digest(path, algo = "xxhash32", serialize = FALSE))
  DBI::dbWriteTable(con, nm, x, temporary = TRUE, overwrite = TRUE)
  DBI::dbExecute(con, sprintf("COPY %s TO %s (FORMAT PARQUET, COMPRESSION ZSTD)", sql_ident(con, nm), sql_string(path, must_work = FALSE)))
  DBI::dbRemoveTable(con, nm)
}
check_pq_tsv <- function(con, pq_path, tsv_path) {
  pq <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(pq_path)))
  tsv <- read_tsv_explicit_normalized(tsv_path, names(pq))
  report <- column_report(pq, tsv$raw, tsv$normalized)
  list(
    consistent = table_equal_after_normalization(pq, tsv$raw, tsv$normalized),
    expected = as.list(expected_schema(names(pq))),
    parquet = as.list(schema_types(pq)),
    tsv_raw = as.list(tsv$raw_types),
    tsv_normalized = as.list(tsv$normalized_types),
    blank_to_na_counts = as.list(tsv$blank_to_na_counts),
    required_character_fields_nonmissing = validate_required_characters(tsv$normalized),
    field_report = report
  )
}
fixture_check <- function(con) {
  fixture <- data.frame(
    target_rsid = c("rs_fixture_1", "rs_fixture_2"),
    included_member = c(TRUE, TRUE),
    excluded_member = c(FALSE, FALSE),
    exposure_effect_allele = c("A", "C"),
    exposure_other_allele = c("G", "T"),
    exposure_beta = c(0.1, 0.2),
    exposure_se = c(0.01, 0.02),
    exposure_pval = c(1e-6, 2e-6),
    exposure_eaf = c(0.2, 0.3),
    outcome_effect_allele_raw = c("A", "C"),
    outcome_other_allele_raw = c("G", "T"),
    outcome_beta_raw = c(0.01, 0.02),
    outcome_se_raw = c(0.001, 0.002),
    outcome_pval_raw = c(0.7, 0.8),
    outcome_eaf_raw = c(0.21, 0.31),
    outcome_effect_allele_harmonised = c("A", NA_character_),
    outcome_other_allele_harmonised = c("G", NA_character_),
    outcome_beta_harmonised = c(0.01, NA_real_),
    outcome_se_harmonised = c(0.001, NA_real_),
    outcome_pval_harmonised = c(0.7, NA_real_),
    outcome_eaf_harmonised = c(0.21, NA_real_),
    palindromic_snp = c(FALSE, TRUE),
    orientation_class = c("exact_match", "exact_match"),
    beta_flipped = c(FALSE, FALSE),
    eaf_flipped = c(FALSE, FALSE),
    strand_flipped = c(FALSE, FALSE),
    record_excluded = c(FALSE, TRUE),
    exclusion_reason = c(NA_character_, "palindromic_snp_excluded_by_reverse_relaxed_rule_v1"),
    final_valid_instrument = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  pq <- tempfile(fileext = ".parquet")
  tsv <- tempfile(fileext = ".tsv")
  write_pq(con, fixture, pq)
  write.table(fixture, tsv, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  chk <- check_pq_tsv(con, pq, tsv)
  harmonised_blank <- chk$field_report$blank_to_na_conversion_count[match("outcome_effect_allele_harmonised", chk$field_report$column_name)] == 1L &&
    chk$field_report$blank_to_na_conversion_count[match("outcome_other_allele_harmonised", chk$field_report$column_name)] == 1L
  exclusion_blank <- chk$field_report$blank_to_na_conversion_count[match("exclusion_reason", chk$field_report$column_name)] == 1L
  list(
    status = if (chk$consistent && harmonised_blank && exclusion_blank && chk$required_character_fields_nonmissing) "passed" else "failed",
    check = chk
  )
}

failure_artifact_rel <- c(
  "R/18_reverse_relaxed_formal_harmonisation_v1.R",
  "docs/decisions/65_vuckovic_hb_reverse_relaxed_formal_harmonisation_v1_v1.1.md",
  "results/logs/vuckovic_hb_reverse_relaxed_formal_harmonisation_v1.log",
  "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_master_v1.parquet.partial",
  "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_master_v1.tsv.partial",
  "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_apoe_included_v1.parquet.partial",
  "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_apoe_included_v1.tsv.partial",
  "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_apoe_excluded_v1.parquet.partial",
  "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_apoe_excluded_v1.tsv.partial",
  "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_counts_v1.csv.partial",
  "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v1.json.partial",
  "R/18_reverse_relaxed_formal_harmonisation_v2.R",
  "docs/decisions/66_vuckovic_hb_reverse_relaxed_formal_harmonisation_v2_readback_fix_v1.1.md",
  "results/logs/vuckovic_hb_reverse_relaxed_formal_harmonisation_v2.log",
  "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_master_v2.parquet.partial",
  "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_master_v2.tsv.partial",
  "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_apoe_included_v2.parquet.partial",
  "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_apoe_included_v2.tsv.partial",
  "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_apoe_excluded_v2.parquet.partial",
  "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_apoe_excluded_v2.tsv.partial",
  "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_counts_v2.csv.partial",
  "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v2.json.partial"
)

out <- c(
  master_pq = file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_relaxed_harmonised_master_v3.parquet"),
  master_tsv = file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_relaxed_harmonised_master_v3.tsv"),
  inc_pq = file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_relaxed_harmonised_apoe_included_v3.parquet"),
  inc_tsv = file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_relaxed_harmonised_apoe_included_v3.tsv"),
  exc_pq = file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_relaxed_harmonised_apoe_excluded_v3.parquet"),
  exc_tsv = file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_relaxed_harmonised_apoe_excluded_v3.tsv"),
  counts = file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_formal_harmonisation_counts_v3.csv"),
  json = file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_formal_harmonisation_v3.json"),
  log = file.path(root, "results", "logs", "vuckovic_hb_reverse_relaxed_formal_harmonisation_v3.log")
)
stop_if(any(file.exists(c(out, paste0(out, ".partial")))), "A reverse relaxed formal harmonisation V3 final or partial target exists; refusing to overwrite.")
dir.create(dirname(out[["master_pq"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["json"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["log"]]), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = out[["log"]], append = TRUE)

main <- function() {
  log_line("stage=reverse_relaxed_formal_harmonisation_v3_tsv_na_normalization_fix")
  renv_lock <- file.path(root, "renv.lock")
  renv_before <- hash_file(renv_lock)
  failure_paths <- file.path(root, failure_artifact_rel)
  stop_if(any(!file.exists(failure_paths)), paste("Missing failed-version artifact(s):", paste(failure_artifact_rel[!file.exists(failure_paths)], collapse = "; ")))
  sha_before <- vapply(failure_paths, hash_file, character(1)); names(sha_before) <- failure_artifact_rel

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  fixture <- fixture_check(con)
  stop_if(!identical(fixture$status, "passed"), "V3 fixture failed.")

  preflight_freeze_path <- file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1_freeze.json")
  rule_path <- file.path(root, "results", "qc", "reverse_relaxed_palindromic_handling_rule_v1.json")
  preflight_pq <- file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1.parquet")
  decision_path <- file.path(root, "docs", "decisions", "67_vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_tsv_na_normalization_fix_v1.1.md")
  stop_if(any(!file.exists(c(preflight_freeze_path, rule_path, preflight_pq, decision_path, renv_lock))), "Required V3 formal harmonisation input is missing.")
  preflight_freeze <- jsonlite::fromJSON(preflight_freeze_path, simplifyVector = FALSE)
  rule <- jsonlite::fromJSON(rule_path, simplifyVector = FALSE)
  pre <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(preflight_pq)))

  h <- data.frame(
    target_rsid = pre$target_rsid,
    included_member = pre$included_member,
    excluded_member = pre$excluded_member,
    exposure_effect_allele = pre$exposure_effect_allele,
    exposure_other_allele = pre$exposure_other_allele,
    exposure_beta = pre$exposure_beta,
    exposure_se = pre$exposure_se,
    exposure_pval = pre$exposure_pval,
    exposure_eaf = pre$exposure_eaf,
    outcome_effect_allele_raw = pre$outcome_effect_allele_raw,
    outcome_other_allele_raw = pre$outcome_other_allele_raw,
    outcome_beta_raw = pre$outcome_beta_raw,
    outcome_se_raw = pre$outcome_se_raw,
    outcome_pval_raw = pre$outcome_pval_raw,
    outcome_eaf_raw = pre$outcome_eaf_raw,
    outcome_effect_allele_harmonised = NA_character_,
    outcome_other_allele_harmonised = NA_character_,
    outcome_beta_harmonised = NA_real_,
    outcome_se_harmonised = NA_real_,
    outcome_pval_harmonised = NA_real_,
    outcome_eaf_harmonised = NA_real_,
    palindromic_snp = pre$palindromic_snp,
    orientation_class = pre$raw_orientation_class,
    beta_flipped = FALSE,
    eaf_flipped = FALSE,
    strand_flipped = FALSE,
    record_excluded = FALSE,
    exclusion_reason = NA_character_,
    final_valid_instrument = FALSE,
    stringsAsFactors = FALSE
  )

  pal <- h$palindromic_snp
  exact <- !pal & h$orientation_class == "exact_match"
  swapped <- !pal & h$orientation_class == "swapped_match"
  strand_exact <- !pal & h$orientation_class == "strand_exact_match"
  strand_swapped <- !pal & h$orientation_class == "strand_swapped_match"
  incompatible <- !pal & h$orientation_class == "incompatible"
  invalid <- !pal & h$orientation_class == "invalid"

  h$record_excluded[pal] <- TRUE
  h$exclusion_reason[pal] <- "palindromic_snp_excluded_by_reverse_relaxed_rule_v1"

  direct <- exact | strand_exact
  h$outcome_effect_allele_harmonised[direct] <- h$exposure_effect_allele[direct]
  h$outcome_other_allele_harmonised[direct] <- h$exposure_other_allele[direct]
  h$outcome_beta_harmonised[direct] <- h$outcome_beta_raw[direct]
  h$outcome_se_harmonised[direct] <- h$outcome_se_raw[direct]
  h$outcome_pval_harmonised[direct] <- h$outcome_pval_raw[direct]
  h$outcome_eaf_harmonised[direct] <- h$outcome_eaf_raw[direct]
  h$strand_flipped[strand_exact] <- TRUE

  reverse <- swapped | strand_swapped
  h$outcome_effect_allele_harmonised[reverse] <- h$exposure_effect_allele[reverse]
  h$outcome_other_allele_harmonised[reverse] <- h$exposure_other_allele[reverse]
  h$outcome_beta_harmonised[reverse] <- -h$outcome_beta_raw[reverse]
  h$outcome_se_harmonised[reverse] <- h$outcome_se_raw[reverse]
  h$outcome_pval_harmonised[reverse] <- h$outcome_pval_raw[reverse]
  h$outcome_eaf_harmonised[reverse] <- ifelse(is.na(h$outcome_eaf_raw[reverse]), NA_real_, 1 - h$outcome_eaf_raw[reverse])
  h$beta_flipped[reverse] <- TRUE
  h$eaf_flipped[reverse] <- TRUE
  h$strand_flipped[strand_swapped] <- TRUE

  h$record_excluded[incompatible | invalid] <- TRUE
  h$exclusion_reason[incompatible] <- "incompatible_alleles"
  h$exclusion_reason[invalid] <- "invalid_alleles"
  h$final_valid_instrument <- !h$record_excluded & h$orientation_class %in% c("exact_match", "swapped_match", "strand_exact_match", "strand_swapped_match")
  included_subset <- h[h$included_member & h$final_valid_instrument, , drop = FALSE]
  excluded_subset <- h[h$excluded_member & h$final_valid_instrument, , drop = FALSE]

  counts <- data.frame(
    metric = c("union_input_count", "included_input_count", "excluded_input_count", "palindromic_count_union",
               "included_palindromic_excluded_count", "excluded_palindromic_excluded_count",
               "included_final_valid_instrument_count", "excluded_final_valid_instrument_count",
               "shared_final_valid_instrument_count", "exact_match_count", "swapped_match_count",
               "strand_exact_match_count", "strand_swapped_match_count", "incompatible_count", "invalid_count",
               "beta_flip_count", "eaf_flip_count", "strand_flip_count"),
    value = c(nrow(h), sum(h$included_member), sum(h$excluded_member), sum(h$palindromic_snp),
              sum(h$included_member & h$record_excluded & h$palindromic_snp),
              sum(h$excluded_member & h$record_excluded & h$palindromic_snp),
              nrow(included_subset), nrow(excluded_subset),
              sum(h$included_member & h$excluded_member & h$final_valid_instrument),
              sum(h$orientation_class == "exact_match"), sum(h$orientation_class == "swapped_match"),
              sum(h$orientation_class == "strand_exact_match"), sum(h$orientation_class == "strand_swapped_match"),
              sum(h$orientation_class == "incompatible"), sum(h$orientation_class == "invalid"),
              sum(h$beta_flipped), sum(h$eaf_flipped), sum(h$strand_flipped)),
    stringsAsFactors = FALSE
  )

  hard_check_failures <- character()
  add_fail <- function(x) hard_check_failures <<- unique(c(hard_check_failures, x))
  if (!identical(preflight_freeze$freeze_status, "passed")) add_fail("preflight_freeze_gate_failed")
  if (!identical(preflight_freeze$approved_for_reverse_relaxed_formal_harmonisation, FALSE)) add_fail("preflight_freeze_unexpected_formal_approval")
  if (!identical(rule$rule_status, "frozen") || !isTRUE(rule$approved_for_reverse_relaxed_formal_harmonisation)) add_fail("palindrome_rule_freeze_gate_failed")
  if (!identical(rule$analysis_direction, "delirium_to_Hb") || !identical(rule$analysis_role, "secondary_reverse_exploratory_relaxed") ||
      !identical(rule$rule_scope, "reverse_relaxed_formal_harmonisation") ||
      !identical(rule$palindromic_definition, "A/T_or_C/G_unordered_allele_pair") ||
      !identical(rule$palindromic_action, "exclude") || isTRUE(rule$eaf_threshold_used) || isTRUE(rule$eaf_based_reinclusion)) add_fail("scientific_parameter_drift_detected")
  if (nrow(h) != as.integer(preflight_freeze$union_preflight_count)) add_fail("all_preflight_records_accounted_for_failed")
  if (any(h$palindromic_snp & !h$record_excluded)) add_fail("palindromic_records_not_excluded_by_rule")
  if (any(h$palindromic_snp & (h$beta_flipped | h$eaf_flipped | h$strand_flipped | h$final_valid_instrument))) add_fail("palindromic_flip_or_reinclusion_detected")
  if (any(!h$palindromic_snp & h$orientation_class %in% c("exact_match", "swapped_match", "strand_exact_match", "strand_swapped_match") & h$record_excluded)) add_fail("nonpalindromic_resolvable_record_excluded")

  exact_np <- !h$palindromic_snp & h$orientation_class == "exact_match"
  if (any(exact_np)) {
    if (!all(h$outcome_effect_allele_harmonised[exact_np] == h$exposure_effect_allele[exact_np] &
             h$outcome_other_allele_harmonised[exact_np] == h$exposure_other_allele[exact_np])) add_fail("exact_match_allele_invariance_failed")
    if (!all(num_equal(h$outcome_beta_harmonised[exact_np], h$outcome_beta_raw[exact_np]))) add_fail("exact_match_beta_invariance_failed")
    if (!all(num_equal(h$outcome_se_harmonised[exact_np], h$outcome_se_raw[exact_np]))) add_fail("exact_match_se_invariance_failed")
    if (!all(num_equal(h$outcome_pval_harmonised[exact_np], h$outcome_pval_raw[exact_np]))) add_fail("exact_match_pval_invariance_failed")
    if (!all(num_equal(h$outcome_eaf_harmonised[exact_np], h$outcome_eaf_raw[exact_np]))) add_fail("exact_match_eaf_invariance_failed")
  }
  if (!all(num_equal(h$exposure_beta, pre$exposure_beta)) || !all(num_equal(h$exposure_se, pre$exposure_se)) ||
      !all(num_equal(h$exposure_pval, pre$exposure_pval)) || !all(num_equal(h$exposure_eaf, pre$exposure_eaf)) ||
      !identical(h$exposure_effect_allele, pre$exposure_effect_allele) ||
      !identical(h$exposure_other_allele, pre$exposure_other_allele)) add_fail("exposure_fields_unchanged_failed")

  write_pq(con, h, paste0(out[["master_pq"]], ".partial"))
  write.table(h, paste0(out[["master_tsv"]], ".partial"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  write_pq(con, included_subset, paste0(out[["inc_pq"]], ".partial"))
  write.table(included_subset, paste0(out[["inc_tsv"]], ".partial"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  write_pq(con, excluded_subset, paste0(out[["exc_pq"]], ".partial"))
  write.table(excluded_subset, paste0(out[["exc_tsv"]], ".partial"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  write.csv(counts, paste0(out[["counts"]], ".partial"), row.names = FALSE)

  master_check <- check_pq_tsv(con, paste0(out[["master_pq"]], ".partial"), paste0(out[["master_tsv"]], ".partial"))
  included_check <- check_pq_tsv(con, paste0(out[["inc_pq"]], ".partial"), paste0(out[["inc_tsv"]], ".partial"))
  excluded_check <- check_pq_tsv(con, paste0(out[["exc_pq"]], ".partial"), paste0(out[["exc_tsv"]], ".partial"))
  if (!master_check$consistent) add_fail("master_parquet_tsv_consistency_failed")
  if (!included_check$consistent) add_fail("included_parquet_tsv_consistency_failed")
  if (!excluded_check$consistent) add_fail("excluded_parquet_tsv_consistency_failed")
  if (!master_check$required_character_fields_nonmissing || !included_check$required_character_fields_nonmissing || !excluded_check$required_character_fields_nonmissing) add_fail("required_character_fields_nonmissing_failed")
  nullable_na_exact <- function(chk) {
    r <- chk$field_report[chk$field_report$column_name %in% nullable_character_columns, , drop = FALSE]
    all(r$na_pattern_consistency)
  }
  if (!nullable_na_exact(master_check)) add_fail("master_nullable_character_na_pattern_exact_failed")
  if (!nullable_na_exact(included_check)) add_fail("included_nullable_character_na_pattern_exact_failed")
  if (!nullable_na_exact(excluded_check)) add_fail("excluded_nullable_character_na_pattern_exact_failed")

  sha_after <- vapply(failure_paths, hash_file, character(1)); names(sha_after) <- failure_artifact_rel
  artifacts_preserved <- identical(sha_before, sha_after)
  if (!artifacts_preserved) add_fail("failed_version_artifacts_changed")
  renv_after <- hash_file(renv_lock)
  if (!identical(renv_before, renv_after)) add_fail("renv_lock_changed")
  v1_idx <- grepl("_v1\\.|Decision 65|65_|v1\\.R|_v1_", names(sha_before))
  v2_idx <- grepl("_v2\\.|Decision 66|66_|v2\\.R|_v2_", names(sha_before))

  hard_checks <- list(
    v1_failure_artifacts_preserved = identical(sha_before[v1_idx], sha_after[v1_idx]),
    v2_failure_artifacts_preserved = identical(sha_before[v2_idx], sha_after[v2_idx]),
    v3_not_using_v1_partial = TRUE,
    v3_not_using_v2_partial = TRUE,
    preflight_freeze_gate = identical(preflight_freeze$freeze_status, "passed"),
    palindrome_rule_freeze_gate = identical(rule$rule_status, "frozen") && isTRUE(rule$approved_for_reverse_relaxed_formal_harmonisation),
    scientific_parameter_drift_none = !("scientific_parameter_drift_detected" %in% hard_check_failures),
    fixture_v3_passed = identical(fixture$status, "passed"),
    explicit_tsv_schema_used = TRUE,
    nullable_character_policy_explicit = identical(nullable_character_columns, c("outcome_effect_allele_harmonised", "outcome_other_allele_harmonised", "exclusion_reason")),
    no_global_blank_to_na_conversion = TRUE,
    required_character_fields_nonmissing = !("required_character_fields_nonmissing_failed" %in% hard_check_failures),
    master_nullable_character_na_pattern_exact = nullable_na_exact(master_check),
    included_nullable_character_na_pattern_exact = nullable_na_exact(included_check),
    excluded_nullable_character_na_pattern_exact = nullable_na_exact(excluded_check),
    all_preflight_records_accounted_for = nrow(h) == as.integer(preflight_freeze$union_preflight_count),
    palindromic_rule_applied = !any(h$palindromic_snp & !h$record_excluded),
    no_eaf_reinclusion = !any(h$palindromic_snp & (h$eaf_flipped | h$final_valid_instrument)),
    exact_match_invariance = !any(grepl("^exact_match_", hard_check_failures)),
    exposure_fields_unchanged = !("exposure_fields_unchanged_failed" %in% hard_check_failures),
    master_parquet_tsv_consistency = master_check$consistent,
    included_parquet_tsv_consistency = included_check$consistent,
    excluded_parquet_tsv_consistency = excluded_check$consistent,
    renv_lock_unchanged = identical(renv_before, renv_after),
    no_reverse_mr = TRUE,
    no_proxy = TRUE,
    no_liftover = TRUE
  )
  harmonisation_status <- if (length(hard_check_failures) == 0L && all(unlist(hard_checks))) "passed" else "failed"
  field_report <- rbind(
    cbind(table_name = "master", master_check$field_report),
    cbind(table_name = "included", included_check$field_report),
    cbind(table_name = "excluded", excluded_check$field_report)
  )
  qc <- list(
    harmonisation_version = "v3",
    decision = 67,
    supersedes_failed_versions = c("v1", "v2"),
    v1_failure_type = "technical_default_tsv_type_inference",
    v2_failure_type = "technical_nullable_character_blank_vs_na_readback",
    scientific_rules_changed = FALSE,
    readback_fix = "explicit_schema_plus_field_specific_nullable_character_blank_to_na_normalization",
    global_blank_to_na_used = FALSE,
    nullable_character_columns = nullable_character_columns,
    required_character_columns = required_character_columns,
    fixture_status = fixture$status,
    fixture_field_report = records(fixture$check$field_report),
    analysis_direction = "delirium_to_Hb",
    analysis_role = "secondary_reverse_exploratory_relaxed",
    p_threshold = 5e-6,
    source_instrument_version = "v2",
    source_outcome_extraction_version = "v1",
    source_preflight_version = "v1",
    source_preflight_freeze_version = "v1",
    palindromic_rule_version = "v1",
    palindromic_rule = "exclude_all_palindromic_snps_without_eaf_based_reinclusion",
    included_input_count = sum(h$included_member),
    excluded_input_count = sum(h$excluded_member),
    union_input_count = nrow(h),
    palindromic_count_union = sum(h$palindromic_snp),
    included_palindromic_excluded_count = sum(h$included_member & h$record_excluded & h$palindromic_snp),
    excluded_palindromic_excluded_count = sum(h$excluded_member & h$record_excluded & h$palindromic_snp),
    included_final_valid_instrument_count = nrow(included_subset),
    excluded_final_valid_instrument_count = nrow(excluded_subset),
    shared_final_valid_instrument_count = sum(h$included_member & h$excluded_member & h$final_valid_instrument),
    exact_match_count = sum(h$orientation_class == "exact_match"),
    swapped_match_count = sum(h$orientation_class == "swapped_match"),
    strand_exact_match_count = sum(h$orientation_class == "strand_exact_match"),
    strand_swapped_match_count = sum(h$orientation_class == "strand_swapped_match"),
    incompatible_count = sum(h$orientation_class == "incompatible"),
    invalid_count = sum(h$orientation_class == "invalid"),
    beta_flip_count = sum(h$beta_flipped),
    eaf_flip_count = sum(h$eaf_flipped),
    strand_flip_count = sum(h$strand_flipped),
    harmonisation_rows = records(h),
    field_consistency_report = records(field_report),
    blank_to_na_conversion_counts = list(
      master = master_check$blank_to_na_counts,
      included = included_check$blank_to_na_counts,
      excluded = excluded_check$blank_to_na_counts
    ),
    master_parquet_tsv_consistency = master_check$consistent,
    included_parquet_tsv_consistency = included_check$consistent,
    excluded_parquet_tsv_consistency = excluded_check$consistent,
    failed_version_artifacts = records(data.frame(
      relative_path = failure_artifact_rel,
      sha256_before = unname(sha_before),
      sha256_after = unname(sha_after),
      unchanged = unname(sha_before) == unname(sha_after),
      stringsAsFactors = FALSE
    )),
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after,
    renv_lock_unchanged = identical(renv_before, renv_after),
    harmonisation_status = harmonisation_status,
    approved_for_reverse_relaxed_mr_input_freeze = identical(harmonisation_status, "passed"),
    hard_checks = hard_checks,
    hard_check_failures = hard_check_failures
  )
  jsonlite::write_json(qc, paste0(out[["json"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  if (!identical(harmonisation_status, "passed")) stop("Reverse relaxed formal harmonisation V3 failed; partial outputs retained.", call. = FALSE)
  for (path in out) {
    stop_if(file.exists(path), paste("Output appeared during run:", path))
    stop_if(!file.rename(paste0(path, ".partial"), path), paste("Atomic rename failed:", path))
  }
  log_line("harmonisation_status=passed included_final_valid=", nrow(included_subset),
           " excluded_final_valid=", nrow(excluded_subset), " hard_check_failures=[]")
}

tryCatch(main(), error = function(e) {
  log_line("status=failed error=", conditionMessage(e))
  quit(status = 1L)
})
