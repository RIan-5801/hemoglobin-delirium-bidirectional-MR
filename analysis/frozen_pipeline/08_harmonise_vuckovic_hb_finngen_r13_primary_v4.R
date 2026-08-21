args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript 08_harmonise_vuckovic_hb_finngen_r13_primary_v4.R --project-root <path>", call. = FALSE)
}

root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
stage <- "startup"
input_path <- file.path(root, "data_derived", "harmonisation_preflight", "vuckovic_hb_finngen_r13_unique_matches_preflight_v2.parquet")
preflight_qc_path <- file.path(root, "results", "qc", "vuckovic_hb_finngen_r13_harmonisation_preflight_v2.json")
out_dir <- file.path(root, "data_derived", "harmonised")
qc_dir <- file.path(root, "results", "qc")
all_audit_path <- file.path(out_dir, "vuckovic_hb_finngen_r13_primary_all_audit_v4.parquet")
inc_pq_path <- file.path(out_dir, "vuckovic_hb_finngen_r13_primary_apoe_included_v4.parquet")
inc_tsv_path <- file.path(out_dir, "vuckovic_hb_finngen_r13_primary_apoe_included_v4.tsv")
exc_pq_path <- file.path(out_dir, "vuckovic_hb_finngen_r13_primary_apoe_excluded_v4.parquet")
exc_tsv_path <- file.path(out_dir, "vuckovic_hb_finngen_r13_primary_apoe_excluded_v4.tsv")
excluded_tsv_path <- file.path(out_dir, "vuckovic_hb_finngen_r13_primary_excluded_variants_v4.tsv")
qc_path <- file.path(qc_dir, "vuckovic_hb_finngen_r13_primary_harmonisation_v4.json")
counts_path <- file.path(qc_dir, "vuckovic_hb_finngen_r13_primary_harmonisation_counts_v4.csv")
log_path <- file.path(root, "results", "logs", "vuckovic_hb_finngen_r13_primary_harmonisation_v4.log")

safe_log <- function(...) cat(paste0(...), "\n", file = log_path, append = TRUE)

preflight_placeholder_columns <- c(
  "outcome_beta_flipped",
  "harmonisation_performed",
  "record_excluded"
)
preflight_placeholder_rename_map <- c(
  outcome_beta_flipped = "preflight_outcome_beta_flipped",
  harmonisation_performed = "preflight_harmonisation_performed",
  record_excluded = "preflight_record_excluded"
)

assert_unique_column_names <- function(data, object_name) {
  duplicate_names <- unique(names(data)[duplicated(names(data))])
  if (length(duplicate_names) > 0L) {
    stop(paste0(object_name, " contains duplicate columns: ", paste(duplicate_names, collapse = ", ")), call. = FALSE)
  }
  TRUE
}

validate_preflight_contract <- function(data) {
  columns_present <- all(preflight_placeholder_columns %in% names(data))
  values_valid <- columns_present && all(vapply(
    preflight_placeholder_columns,
    function(column_name) {
      x <- data[[column_name]]
      is.logical(x) && !anyNA(x) && all(x == FALSE)
    },
    logical(1)
  ))
  list(
    columns_present = columns_present,
    values_valid = values_valid,
    optional_eaf_flip_absent = !("outcome_eaf_flipped" %in% names(data))
  )
}

rename_preflight_placeholders <- function(data) {
  missing_columns <- setdiff(names(preflight_placeholder_rename_map), names(data))
  if (length(missing_columns) > 0L) {
    stop(paste0("Missing preflight placeholder columns: ", paste(missing_columns, collapse = ", ")), call. = FALSE)
  }
  new_names <- unname(preflight_placeholder_rename_map)
  if (length(intersect(new_names, names(data))) > 0L) {
    stop("Preflight placeholder destination column already exists.", call. = FALSE)
  }
  for (old_name in names(preflight_placeholder_rename_map)) {
    names(data)[names(data) == old_name] <- preflight_placeholder_rename_map[[old_name]]
  }
  assert_unique_column_names(data, "renamed_preflight")
  data
}

complement_allele <- function(x) {
  unname(c(A = "T", C = "G", G = "C", T = "A")[[toupper(as.character(x))]])
}

transform_row <- function(row) {
  orientation <- as.character(row$raw_allele_orientation)
  beta <- as.numeric(row$outcome_beta_raw)
  eaf <- as.numeric(row$outcome_eaf_raw)
  outcome_effect <- as.character(row$outcome_effect_allele_raw)
  outcome_other <- as.character(row$outcome_other_allele_raw)
  beta_flipped <- FALSE
  eaf_flipped <- FALSE
  strand_complemented <- FALSE
  harmonisation_performed <- FALSE
  record_excluded <- TRUE
  mr_keep_primary <- FALSE
  harmonisation_action <- "exclude_incompatible"
  harmonisation_status <- "excluded_primary"
  exclusion_reason <- "incompatible_alleles"

  if (orientation == "exact_match") {
    harmonisation_action <- "exact_keep"
    harmonisation_status <- "kept_primary"
    exclusion_reason <- NA_character_
    harmonisation_performed <- TRUE
    record_excluded <- FALSE
    mr_keep_primary <- TRUE
  } else if (orientation == "swapped_match") {
    harmonisation_action <- "swapped_keep"
    harmonisation_status <- "kept_primary"
    exclusion_reason <- NA_character_
    beta <- -beta
    beta_flipped <- TRUE
    if (!is.na(eaf)) {
      eaf <- 1 - eaf
      eaf_flipped <- TRUE
    }
    outcome_effect <- as.character(row$outcome_other_allele_raw)
    outcome_other <- as.character(row$outcome_effect_allele_raw)
    harmonisation_performed <- TRUE
    record_excluded <- FALSE
    mr_keep_primary <- TRUE
  } else if (orientation == "strand_exact_match") {
    harmonisation_action <- "strand_exact_keep"
    harmonisation_status <- "kept_primary"
    exclusion_reason <- NA_character_
    outcome_effect <- complement_allele(row$outcome_effect_allele_raw)
    outcome_other <- complement_allele(row$outcome_other_allele_raw)
    strand_complemented <- TRUE
    harmonisation_performed <- TRUE
    record_excluded <- FALSE
    mr_keep_primary <- TRUE
  } else if (orientation == "strand_swapped_match") {
    harmonisation_action <- "strand_swapped_keep"
    harmonisation_status <- "kept_primary"
    exclusion_reason <- NA_character_
    beta <- -beta
    beta_flipped <- TRUE
    if (!is.na(eaf)) {
      eaf <- 1 - eaf
      eaf_flipped <- TRUE
    }
    outcome_effect <- complement_allele(row$outcome_other_allele_raw)
    outcome_other <- complement_allele(row$outcome_effect_allele_raw)
    strand_complemented <- TRUE
    harmonisation_performed <- TRUE
    record_excluded <- FALSE
    mr_keep_primary <- TRUE
  } else if (orientation == "palindromic_snp") {
    harmonisation_action <- "exclude_palindromic"
    exclusion_reason <- "palindromic_excluded_from_primary"
  } else if (orientation == "allele_missing_or_invalid") {
    harmonisation_action <- "exclude_invalid"
    exclusion_reason <- "allele_missing_or_invalid"
  }

  data.frame(
    outcome_effect_allele_aligned = outcome_effect,
    outcome_other_allele_aligned = outcome_other,
    outcome_beta = beta,
    outcome_se = as.numeric(row$outcome_se_raw),
    outcome_pval = as.numeric(row$outcome_pval_raw),
    outcome_eaf = eaf,
    outcome_beta_flipped = beta_flipped,
    outcome_eaf_flipped = eaf_flipped,
    strand_complemented = strand_complemented,
    harmonisation_performed = harmonisation_performed,
    record_excluded = record_excluded,
    harmonisation_action = harmonisation_action,
    harmonisation_status = harmonisation_status,
    exclusion_reason = exclusion_reason,
    mr_keep_primary = mr_keep_primary,
    stringsAsFactors = FALSE
  )
}

make_contract_fixture <- function() {
  data.frame(
    rsid = "rs_fixture",
    raw_allele_orientation = "exact_match",
    outcome_effect_allele_raw = "A",
    outcome_other_allele_raw = "G",
    outcome_beta_raw = 0.2,
    outcome_se_raw = 0.05,
    outcome_pval_raw = 0.01,
    outcome_eaf_raw = 0.3,
    outcome_beta_flipped = FALSE,
    harmonisation_performed = FALSE,
    record_excluded = FALSE,
    stringsAsFactors = FALSE
  )
}

run_schema_unit_tests <- function() {
  fixture <- make_contract_fixture()
  contract <- validate_preflight_contract(fixture)
  renamed <- rename_preflight_placeholders(fixture)
  renamed_columns_present <- all(unname(preflight_placeholder_rename_map) %in% names(renamed)) &&
    !any(names(preflight_placeholder_rename_map) %in% names(renamed))
  formal <- cbind(renamed, transform_row(renamed))
  formal_fields_present <- all(c(
    "outcome_beta_flipped", "outcome_eaf_flipped",
    "harmonisation_performed", "record_excluded"
  ) %in% names(formal))
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(con, "v4_schema_fixture", formal, temporary = TRUE)
  registered <- DBI::dbGetQuery(con, "SELECT * FROM v4_schema_fixture")
  list(
    preflight_contract_unit_test = isTRUE(contract$columns_present) &&
      isTRUE(contract$values_valid) && isTRUE(contract$optional_eaf_flip_absent),
    placeholder_rename_unit_test = renamed_columns_present,
    formal_field_creation_unit_test = formal_fields_present && !anyDuplicated(names(formal)),
    duplicate_column_guard_unit_test = !anyDuplicated(names(renamed)) && !anyDuplicated(names(formal)),
    duckdb_registration_unit_test = nrow(registered) == 1L && !anyDuplicated(names(registered))
  )
}

run_harmonisation_unit_tests <- function() {
  exact <- transform_row(rename_preflight_placeholders(make_contract_fixture()))
  swapped_fixture <- make_contract_fixture()
  swapped_fixture$raw_allele_orientation <- "swapped_match"
  swapped_fixture$outcome_effect_allele_raw <- "G"
  swapped_fixture$outcome_other_allele_raw <- "A"
  swapped <- transform_row(rename_preflight_placeholders(swapped_fixture))
  pal_fixture <- make_contract_fixture()
  pal_fixture$raw_allele_orientation <- "palindromic_snp"
  pal <- transform_row(rename_preflight_placeholders(pal_fixture))
  isTRUE(exact$mr_keep_primary) && !isTRUE(exact$outcome_beta_flipped) &&
    !isTRUE(exact$outcome_eaf_flipped) &&
    isTRUE(swapped$mr_keep_primary) && isTRUE(swapped$outcome_beta_flipped) &&
    isTRUE(swapped$outcome_eaf_flipped) && identical(swapped$outcome_beta, -0.2) &&
    identical(swapped$outcome_eaf, 0.7) && identical(swapped$outcome_se, 0.05) &&
    identical(swapped$outcome_pval, 0.01) && !isTRUE(pal$mr_keep_primary) &&
    isTRUE(pal$record_excluded) && !isTRUE(pal$harmonisation_performed)
}

read_parquet <- function(con, path) {
  DBI::dbGetQuery(con, sprintf(
    "SELECT * FROM read_parquet('%s')",
    normalizePath(path, winslash = "/", mustWork = TRUE)
  ))
}

write_parquet <- function(con, data, path) {
  table_name <- "v4_harmonisation_write"
  DBI::dbWriteTable(con, table_name, data, temporary = TRUE, overwrite = TRUE)
  DBI::dbExecute(con, sprintf(
    "COPY %s TO '%s' (FORMAT PARQUET)", table_name,
    normalizePath(path, winslash = "/", mustWork = FALSE)
  ))
  DBI::dbRemoveTable(con, table_name)
}

numeric_equal <- function(x, y, tolerance = 1e-12) {
  if (length(x) != length(y)) return(FALSE)
  x <- as.numeric(x)
  y <- as.numeric(y)
  na_both <- is.na(x) & is.na(y)
  both <- !is.na(x) & !is.na(y)
  checks <- rep(FALSE, length(x))
  checks[na_both] <- TRUE
  checks[both] <- is.finite(x[both]) & is.finite(y[both]) & abs(x[both] - y[both]) <= tolerance
  all(checks)
}

value_equal <- function(x, y) {
  if (is.numeric(x) || is.integer(x)) return(numeric_equal(x, y))
  if (length(x) != length(y)) return(FALSE)
  na_both <- is.na(x) & is.na(y)
  both <- !is.na(x) & !is.na(y)
  all(na_both | (both & as.character(x) == as.character(y)))
}

compare_parquet_tsv <- function(parquet_data, tsv_data) {
  if (nrow(parquet_data) != nrow(tsv_data) || !identical(names(parquet_data), names(tsv_data))) return(FALSE)
  if (anyDuplicated(parquet_data$rsid) || anyDuplicated(tsv_data$rsid) || !setequal(parquet_data$rsid, tsv_data$rsid)) return(FALSE)
  parquet_data <- parquet_data[order(parquet_data$rsid), , drop = FALSE]
  tsv_data <- tsv_data[order(tsv_data$rsid), , drop = FALSE]
  all(vapply(names(parquet_data), function(column_name) value_equal(parquet_data[[column_name]], tsv_data[[column_name]]), logical(1)))
}

write_partial <- function(writer, final_path) {
  writer(paste0(final_path, ".partial"))
}

rename_partial <- function(final_path) {
  partial_path <- paste0(final_path, ".partial")
  if (!file.rename(partial_path, final_path)) stop(paste0("Atomic rename failed: ", final_path), call. = FALSE)
}

main <- function() {
  targets <- c(
    all_audit_path, inc_pq_path, inc_tsv_path, exc_pq_path, exc_tsv_path,
    excluded_tsv_path, qc_path, counts_path, log_path
  )
  if (any(file.exists(c(targets, paste0(targets, ".partial"))))) {
    stop("Primary harmonisation V4 target or partial exists; refusing to overwrite.", call. = FALSE)
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (!file.exists(input_path) || !file.exists(preflight_qc_path)) stop("Required preflight V2 input is absent.", call. = FALSE)

  stage <<- "preflight_gate"
  preflight_qc <- jsonlite::fromJSON(preflight_qc_path, simplifyVector = FALSE)
  if (preflight_qc$preflight_status != "passed" || !isTRUE(preflight_qc$approved_for_harmonisation_rule_freeze)) {
    stop("Preflight V2 QC approval gate failed.", call. = FALSE)
  }

  stage <<- "unit_tests"
  schema_tests <- run_schema_unit_tests()
  if (!all(vapply(schema_tests, isTRUE, logical(1)))) stop("V4 schema unit tests failed.", call. = FALSE)
  harmonisation_unit_tests <- run_harmonisation_unit_tests()
  if (!isTRUE(harmonisation_unit_tests)) stop("V4 harmonisation unit tests failed.", call. = FALSE)
  for (test_name in names(schema_tests)) safe_log(test_name, "=passed")
  safe_log("harmonisation_unit_tests=passed")

  stage <<- "input_read_and_contract_validation"
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  preflight <- read_parquet(con, input_path)
  contract <- validate_preflight_contract(preflight)
  safe_log("preflight_contract_columns_present=", contract$columns_present)
  safe_log("preflight_contract_values_valid=", contract$values_valid)
  safe_log("preflight_optional_eaf_flip_absent=", contract$optional_eaf_flip_absent)
  if (!all(vapply(contract, isTRUE, logical(1)))) stop("Preflight V2 contract gate failed.", call. = FALSE)
  if (nrow(preflight) != 362L || anyNA(preflight$rsid) || anyDuplicated(preflight$rsid) || "rs8070692" %in% preflight$rsid) {
    stop("Preflight V2 input integrity gate failed.", call. = FALSE)
  }
  preflight <- rename_preflight_placeholders(preflight)

  stage <<- "formal_harmonisation"
  source_columns <- c(
    "exposure_effect_allele", "exposure_other_allele", "exposure_beta", "exposure_se", "exposure_eaf",
    "outcome_effect_allele", "outcome_other_allele", "outcome_beta", "outcome_se", "outcome_pval", "outcome_eaf"
  )
  if (!all(source_columns %in% names(preflight))) stop("Required raw source columns absent.", call. = FALSE)
  for (column_name in source_columns) preflight[[paste0(column_name, "_raw")]] <- preflight[[column_name]]
  preflight <- preflight[, setdiff(names(preflight), c("outcome_beta", "outcome_se", "outcome_pval", "outcome_eaf")), drop = FALSE]
  formal_rows <- do.call(rbind, lapply(seq_len(nrow(preflight)), function(i) transform_row(preflight[i, , drop = FALSE])))
  all_audit <- cbind(preflight, formal_rows)
  all_audit$effect_allele <- all_audit$exposure_effect_allele_raw
  all_audit$other_allele <- all_audit$exposure_other_allele_raw
  all_audit$exposure_beta <- all_audit$exposure_beta_raw
  all_audit$exposure_se <- all_audit$exposure_se_raw
  all_audit$exposure_eaf <- all_audit$exposure_eaf_raw
  all_audit$exposure_source <- "Vuckovic_2020_Hb"
  all_audit$exposure_trait <- "haemoglobin_concentration"
  all_audit$exposure_n_study <- 408112
  all_audit$exposure_n_variant <- NA_real_
  all_audit$exposure_effect_scale <- "standardized_inverse_normal_transformed_Hb"
  all_audit$outcome_source <- "FinnGen_R13_F5_DELIRIUM"
  all_audit$outcome_trait <- "delirium"
  all_audit$outcome_ncase <- 5121
  all_audit$outcome_ncontrol <- 465023
  all_audit$outcome_n_study <- 470144
  all_audit$outcome_effect_scale <- "log_odds"
  assert_unique_column_names(all_audit, "all_audit")

  stage <<- "subsets_and_hard_checks"
  exact <- all_audit$raw_allele_orientation == "exact_match"
  swapped <- all_audit$raw_allele_orientation == "swapped_match"
  palindromic <- all_audit$raw_allele_orientation == "palindromic_snp"
  keep <- all_audit$mr_keep_primary
  included <- all_audit[keep & all_audit$instrument_membership %in% c("shared", "included_only"), , drop = FALSE]
  excluded <- all_audit[keep & all_audit$instrument_membership %in% c("shared", "excluded_only"), , drop = FALSE]
  excluded_variants <- all_audit[!keep, , drop = FALSE]
  assert_unique_column_names(included, "included")
  assert_unique_column_names(excluded, "excluded")
  assert_unique_column_names(excluded_variants, "excluded_variants")

  hard_checks <- list(
    preflight_gate = preflight_qc$preflight_status == "passed" && isTRUE(preflight_qc$approved_for_harmonisation_rule_freeze),
    preflight_contract_columns_present = contract$columns_present,
    preflight_contract_values_valid = contract$values_valid,
    preflight_optional_eaf_flip_absent = contract$optional_eaf_flip_absent,
    preflight_contract_unit_test = schema_tests$preflight_contract_unit_test,
    placeholder_rename_unit_test = schema_tests$placeholder_rename_unit_test,
    formal_field_creation_unit_test = schema_tests$formal_field_creation_unit_test,
    duplicate_column_guard_unit_test = schema_tests$duplicate_column_guard_unit_test,
    duckdb_registration_unit_test = schema_tests$duckdb_registration_unit_test,
    fixture_schema_complete = TRUE,
    harmonisation_unit_tests = harmonisation_unit_tests,
    input_integrity = nrow(all_audit) == 362L && !anyNA(all_audit$rsid) && !anyDuplicated(all_audit$rsid),
    unique_input_count = nrow(all_audit) == 362L,
    exact_match_count = sum(exact) == 311L,
    swapped_match_count = sum(swapped) == 2L,
    palindromic_count = sum(palindromic) == 49L,
    no_incompatible_variants = !any(all_audit$raw_allele_orientation %in% c("incompatible_alleles", "allele_missing_or_invalid")),
    all_exact_unflipped = all(!all_audit$outcome_beta_flipped[exact]) && all(!all_audit$outcome_eaf_flipped[exact]),
    all_swapped_flipped = all(all_audit$outcome_beta_flipped[swapped]) && all(all_audit$outcome_eaf_flipped[swapped]),
    flip_count_equals_swapped_count = sum(all_audit$outcome_beta_flipped) == sum(swapped),
    se_unchanged = numeric_equal(all_audit$outcome_se, all_audit$outcome_se_raw),
    pval_unchanged = numeric_equal(all_audit$outcome_pval, all_audit$outcome_pval_raw),
    aligned_alleles_match_exposure = all(all_audit$outcome_effect_allele_aligned[keep] == all_audit$effect_allele[keep]) && all(all_audit$outcome_other_allele_aligned[keep] == all_audit$other_allele[keep]),
    no_palindromic_in_primary = !any(included$palindromic_snp) && !any(excluded$palindromic_snp),
    all_palindromic_in_exclusion_audit = nrow(excluded_variants) == 49L && all(excluded_variants$palindromic_snp) && all(excluded_variants$exclusion_reason == "palindromic_excluded_from_primary"),
    included_membership_correct = all(included$instrument_membership %in% c("shared", "included_only")),
    excluded_membership_correct = all(excluded$instrument_membership %in% c("shared", "excluded_only")),
    included_rsid_unique = !anyDuplicated(included$rsid),
    excluded_rsid_unique = !anyDuplicated(excluded$rsid),
    all_audit_column_names_unique = !anyDuplicated(names(all_audit)),
    included_column_names_unique = !anyDuplicated(names(included)),
    excluded_column_names_unique = !anyDuplicated(names(excluded)),
    excluded_variants_column_names_unique = !anyDuplicated(names(excluded_variants)),
    no_multiple_target_reintroduced = !any(all_audit$rsid == "rs8070692"),
    no_missing_target_reintroduced = nrow(all_audit) == 362L,
    no_proxy = TRUE,
    no_liftover = TRUE
  )
  prewrite_failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
  if (length(prewrite_failures) > 0L) stop(paste0("Pre-write hard checks failed: ", paste(prewrite_failures, collapse = "; ")), call. = FALSE)

  stage <<- "partial_write_and_readback"
  partial_targets <- c(all_audit_path, inc_pq_path, inc_tsv_path, exc_pq_path, exc_tsv_path, excluded_tsv_path, counts_path, qc_path)
  on.exit(unlink(paste0(partial_targets, ".partial"), force = TRUE), add = TRUE)
  write_partial(function(path) write_parquet(con, all_audit, path), all_audit_path)
  write_partial(function(path) write_parquet(con, included, path), inc_pq_path)
  write_partial(function(path) write.table(included, path, sep = "\t", quote = FALSE, row.names = FALSE, na = ""), inc_tsv_path)
  write_partial(function(path) write_parquet(con, excluded, path), exc_pq_path)
  write_partial(function(path) write.table(excluded, path, sep = "\t", quote = FALSE, row.names = FALSE, na = ""), exc_tsv_path)
  excluded_audit_columns <- c("rsid", "raw_allele_orientation", "palindromic_snp", "instrument_membership", "exposure_eaf_raw", "outcome_eaf_raw", "exposure_maf", "outcome_maf", "preferred_frequency_orientation", "frequency_orientation_margin", "exclusion_reason")
  if (!all(excluded_audit_columns %in% names(excluded_variants))) stop("Excluded-variant audit columns absent.", call. = FALSE)
  write_partial(function(path) write.table(excluded_variants[, excluded_audit_columns, drop = FALSE], path, sep = "\t", quote = FALSE, row.names = FALSE, na = ""), excluded_tsv_path)

  hard_checks$included_parquet_tsv_consistency <- compare_parquet_tsv(read_parquet(con, paste0(inc_pq_path, ".partial")), read.delim(paste0(inc_tsv_path, ".partial"), check.names = FALSE, stringsAsFactors = FALSE, na.strings = ""))
  hard_checks$excluded_parquet_tsv_consistency <- compare_parquet_tsv(read_parquet(con, paste0(exc_pq_path, ".partial")), read.delim(paste0(exc_tsv_path, ".partial"), check.names = FALSE, stringsAsFactors = FALSE, na.strings = ""))
  failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
  if (length(failures) > 0L) stop(paste0("Readback hard checks failed: ", paste(failures, collapse = "; ")), call. = FALSE)

  counts <- data.frame(
    metric = c("input_unique", "exact_match", "swapped_match", "palindromic_excluded", "primary_apoe_included", "primary_apoe_excluded", "primary_shared", "primary_included_only", "primary_excluded_only"),
    count = c(nrow(all_audit), sum(exact), sum(swapped), sum(palindromic), nrow(included), nrow(excluded), sum(keep & all_audit$instrument_membership == "shared"), sum(keep & all_audit$instrument_membership == "included_only"), sum(keep & all_audit$instrument_membership == "excluded_only")),
    stringsAsFactors = FALSE
  )
  write_partial(function(path) write.csv(counts, path, row.names = FALSE), counts_path)
  qc <- list(
    harmonisation_version = "v4",
    analysis_role = "forward_primary_non_palindromic",
    source_extraction_version = "v9",
    source_final_audit_version = "v7",
    source_preflight_version = "v2",
    harmonisation_status = "passed",
    approved_for_forward_primary_mr = TRUE,
    preflight_contract_columns = preflight_placeholder_columns,
    preflight_contract_columns_present = contract$columns_present,
    preflight_contract_values_valid = contract$values_valid,
    preflight_optional_eaf_flip_absent = contract$optional_eaf_flip_absent,
    preflight_placeholder_rename_map = preflight_placeholder_rename_map,
    input_unique_count = nrow(all_audit),
    exact_match_count = sum(exact),
    swapped_match_count = sum(swapped),
    palindromic_excluded_count = sum(palindromic),
    incompatible_excluded_count = 0L,
    invalid_excluded_count = 0L,
    outcome_beta_flip_count = sum(all_audit$outcome_beta_flipped),
    outcome_eaf_flip_count = sum(all_audit$outcome_eaf_flipped),
    included_primary_count = nrow(included),
    excluded_primary_count = nrow(excluded),
    shared_primary_count = sum(keep & all_audit$instrument_membership == "shared"),
    included_only_primary_count = sum(keep & all_audit$instrument_membership == "included_only"),
    excluded_only_primary_count = sum(keep & all_audit$instrument_membership == "excluded_only"),
    all_audit_duplicate_column_count = sum(duplicated(names(all_audit))),
    included_duplicate_column_count = sum(duplicated(names(included))),
    excluded_duplicate_column_count = sum(duplicated(names(excluded))),
    excluded_variants_duplicate_column_count = sum(duplicated(names(excluded_variants))),
    multiple_target_reintroduced = FALSE,
    proxy_used = FALSE,
    liftover_used = FALSE,
    hard_checks = hard_checks,
    hard_check_failures = character(0),
    informational_findings = list(palindromic_rule = "all palindromic SNPs excluded from primary without EAF-based re-inclusion")
  )
  write_partial(function(path) writeLines(jsonlite::toJSON(qc, auto_unbox = TRUE, pretty = TRUE), path), qc_path)

  stage <<- "atomic_publish"
  for (final_path in c(all_audit_path, inc_pq_path, inc_tsv_path, exc_pq_path, exc_tsv_path, excluded_tsv_path, counts_path, qc_path)) rename_partial(final_path)
  safe_log("harmonisation_status=passed")
  safe_log("hard_check_failures=")
  safe_log("approved_for_forward_primary_mr=TRUE")
  0L
}

exit_status <- tryCatch(
  main(),
  error = function(error) {
    safe_log("harmonisation_status=failed")
    safe_log("stage=", stage)
    safe_log("error=", conditionMessage(error))
    1L
  }
)
quit(status = exit_status)
