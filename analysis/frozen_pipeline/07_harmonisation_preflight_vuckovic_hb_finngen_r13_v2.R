args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript 07_harmonisation_preflight_vuckovic_hb_finngen_r13_v2.R --project-root <path>", call. = FALSE)
}

root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
stage <- "startup"
preflight_dir <- file.path(root, "data_derived", "harmonisation_preflight")
qc_dir <- file.path(root, "results", "qc")
main_output <- file.path(preflight_dir, "vuckovic_hb_finngen_r13_unique_matches_preflight_v2.parquet")
pal_output <- file.path(preflight_dir, "vuckovic_hb_finngen_r13_palindromic_preflight_v2.tsv")
qc_output <- file.path(qc_dir, "vuckovic_hb_finngen_r13_harmonisation_preflight_v2.json")
orientation_output <- file.path(qc_dir, "vuckovic_hb_finngen_r13_harmonisation_orientation_counts_v2.csv")
log_output <- file.path(root, "results", "logs", "vuckovic_hb_finngen_r13_harmonisation_preflight_v2.log")

safe_log <- function(...) {
  cat(paste0(...), "\n", file = log_output, append = TRUE)
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

is_valid_allele <- function(x) {
  length(x) == 1L && !is.na(x) && toupper(as.character(x)) %in% c("A", "C", "G", "T")
}

complement_allele <- function(x) {
  mapping <- c(A = "T", C = "G", G = "C", T = "A")
  unname(mapping[[toupper(as.character(x))]])
}

is_palindromic_pair <- function(effect_allele, other_allele) {
  if (!is_valid_allele(effect_allele) || !is_valid_allele(other_allele)) return(FALSE)
  allele_set <- sort(c(toupper(effect_allele), toupper(other_allele)))
  identical(allele_set, c("A", "T")) || identical(allele_set, c("C", "G"))
}

classify_orientation <- function(ea, oa, oe, oo) {
  if (!all(vapply(c(ea, oa, oe, oo), is_valid_allele, logical(1)))) return("allele_missing_or_invalid")
  ea <- toupper(ea)
  oa <- toupper(oa)
  oe <- toupper(oe)
  oo <- toupper(oo)
  if (is_palindromic_pair(ea, oa)) return("palindromic_snp")
  if (identical(ea, oe) && identical(oa, oo)) return("exact_match")
  if (identical(ea, oo) && identical(oa, oe)) return("swapped_match")
  if (identical(ea, complement_allele(oe)) && identical(oa, complement_allele(oo))) return("strand_exact_match")
  if (identical(ea, complement_allele(oo)) && identical(oa, complement_allele(oe))) return("strand_swapped_match")
  "incompatible_alleles"
}

frequency_diagnostic <- function(exposure_eaf, outcome_eaf) {
  exposure_eaf <- suppressWarnings(as.numeric(exposure_eaf))
  outcome_eaf <- suppressWarnings(as.numeric(outcome_eaf))
  if (length(exposure_eaf) != 1L || length(outcome_eaf) != 1L || is.na(exposure_eaf) || is.na(outcome_eaf) || !is.finite(exposure_eaf) || !is.finite(outcome_eaf)) {
    return(list(exposure_maf = NA_real_, outcome_maf = NA_real_, frequency_difference_same = NA_real_, frequency_difference_flipped = NA_real_, frequency_orientation_margin = NA_real_, preferred_frequency_orientation = "missing_frequency"))
  }
  same_difference <- abs(exposure_eaf - outcome_eaf)
  flipped_difference <- abs(exposure_eaf - (1 - outcome_eaf))
  preference <- "tie"
  if (same_difference < flipped_difference) preference <- "same"
  if (flipped_difference < same_difference) preference <- "flipped"
  list(exposure_maf = min(exposure_eaf, 1 - exposure_eaf), outcome_maf = min(outcome_eaf, 1 - outcome_eaf), frequency_difference_same = same_difference, frequency_difference_flipped = flipped_difference, frequency_orientation_margin = abs(same_difference - flipped_difference), preferred_frequency_orientation = preference)
}

classify_instrument_membership <- function(
  in_apoe_included_input,
  in_apoe_excluded_input
) {
  if (length(in_apoe_included_input) != length(in_apoe_excluded_input)) {
    stop("APOE membership vectors have different lengths.", call. = FALSE)
  }
  included <- as.logical(in_apoe_included_input)
  excluded <- as.logical(in_apoe_excluded_input)
  membership <- rep("invalid_membership", length(included))
  valid_flags <- !is.na(included) & !is.na(excluded)
  membership[valid_flags & included & excluded] <- "shared"
  membership[valid_flags & included & !excluded] <- "included_only"
  membership[valid_flags & !included & excluded] <- "excluded_only"
  membership
}

read_parquet_file <- function(con, path) {
  input_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", input_path))
}

write_parquet_file <- function(con, data, path) {
  table_name <- "preflight_v2_write"
  DBI::dbWriteTable(con, table_name, data, temporary = TRUE, overwrite = TRUE)
  output_path <- paste0(normalizePath(dirname(path), winslash = "/", mustWork = TRUE), "/", basename(path))
  DBI::dbExecute(con, sprintf("COPY %s TO '%s' (FORMAT PARQUET)", table_name, output_path))
  DBI::dbRemoveTable(con, table_name)
}

write_atomic <- function(writer, final_path) {
  partial_path <- paste0(final_path, ".partial")
  on.exit(if (file.exists(partial_path)) unlink(partial_path), add = TRUE)
  writer(partial_path)
  if (!file.rename(partial_path, final_path)) stop(paste("Atomic rename failed:", final_path), call. = FALSE)
}

run_unit_tests <- function() {
  orientation_pass <- all(c(
    classify_orientation("A", "G", "A", "G") == "exact_match",
    classify_orientation("A", "G", "G", "A") == "swapped_match",
    classify_orientation("A", "G", "T", "C") == "strand_exact_match",
    classify_orientation("A", "G", "C", "T") == "strand_swapped_match",
    classify_orientation("A", "T", "A", "T") == "palindromic_snp",
    classify_orientation("C", "G", "G", "C") == "palindromic_snp",
    classify_orientation("A", "G", "A", "C") == "incompatible_alleles",
    classify_orientation("N", "G", "A", "C") == "allele_missing_or_invalid"
  ))
  frequency_pass <- all(c(
    frequency_diagnostic(0.20, 0.22)$preferred_frequency_orientation == "same",
    frequency_diagnostic(0.20, 0.78)$preferred_frequency_orientation == "flipped",
    frequency_diagnostic(NA_real_, 0.20)$preferred_frequency_orientation == "missing_frequency"
  ))
  observed <- classify_instrument_membership(c(TRUE, TRUE, FALSE, FALSE, NA), c(TRUE, FALSE, TRUE, FALSE, TRUE))
  membership_pass <- identical(observed, c("shared", "included_only", "excluded_only", "invalid_membership", "invalid_membership"))
  mixed <- classify_instrument_membership(c(TRUE, TRUE, FALSE), c(TRUE, FALSE, TRUE))
  mixed_pass <- identical(mixed, c("shared", "included_only", "excluded_only"))
  list(allele_orientation_unit_tests = orientation_pass, frequency_diagnostic_unit_tests = frequency_pass, membership_unit_test = membership_pass, mixed_membership_unit_test = mixed_pass)
}

main <- function() {
  stage <<- "preflight_paths"
  outputs <- c(main_output, pal_output, qc_output, orientation_output, log_output)
  if (any(file.exists(c(outputs, paste0(outputs, ".partial"))))) stop("Preflight V2 target or partial target exists; refusing to overwrite.", call. = FALSE)
  if (!dir.exists(preflight_dir)) dir.create(preflight_dir, recursive = TRUE)

  stage <<- "unit_tests"
  unit_tests <- run_unit_tests()
  if (!all(vapply(unit_tests, isTRUE, logical(1)))) stop("V2 unit tests failed before V9 artifacts were read.", call. = FALSE)
  for (name in names(unit_tests)) safe_log(name, "=passed")

  stage <<- "v7_audit_gate"
  v7_qc <- jsonlite::fromJSON(file.path(qc_dir, "vuckovic_hb_finngen_r13_extraction_v9_final_audit_v7.json"), simplifyVector = FALSE)
  if (v7_qc$audit_status != "passed" || !isTRUE(v7_qc$approved_for_harmonisation_development)) stop("V7 final-audit gate failed.", call. = FALSE)
  manifest <- read.csv(file.path(qc_dir, "vuckovic_hb_finngen_r13_extraction_v9_output_sha_manifest_v4.csv"), stringsAsFactors = FALSE)

  stage <<- "input_sha_gate"
  paths <- c(
    allele = file.path(root, "data_derived", "outcome_extracts", "audit", "vuckovic_hb_finngen_r13_allele_audit_v9.parquet"),
    target = file.path(root, "data_derived", "outcome_extracts", "vuckovic_hb_target_union_v9.parquet"),
    included = file.path(root, "data_derived", "outcome_extracts", "vuckovic_hb_finngen_r13_apoe_included_v9.parquet"),
    excluded = file.path(root, "data_derived", "outcome_extracts", "vuckovic_hb_finngen_r13_apoe_excluded_v9.parquet")
  )
  manifest_keys <- c(allele = "allele", target = "target_union", included = "inc_pq", excluded = "exc_pq")
  observed_sha <- vapply(paths, sha256_file, character(1))
  expected_sha <- vapply(names(paths), function(name) manifest$observed_sha256[manifest$key == manifest_keys[[name]]], character(1))
  if (!all(tolower(observed_sha) == tolower(expected_sha))) stop("V9 input SHA gate failed.", call. = FALSE)

  stage <<- "small_input_read"
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  preflight <- read_parquet_file(con, paths[["allele"]])
  target <- read_parquet_file(con, paths[["target"]])
  included_output <- read_parquet_file(con, paths[["included"]])
  excluded_output <- read_parquet_file(con, paths[["excluded"]])

  stage <<- "rowwise_membership"
  if (nrow(preflight) != 362L || anyNA(preflight$rsid) || anyDuplicated(preflight$rsid)) stop("Unique-match input integrity gate failed.", call. = FALSE)
  if (nrow(included_output) != 361L || nrow(excluded_output) != 361L) stop("Formal output count gate failed.", call. = FALSE)
  multiple_target_rsid <- "rs8070692"
  if (multiple_target_rsid %in% c(preflight$rsid, included_output$rsid, excluded_output$rsid)) stop("Multiple target entered unique output.", call. = FALSE)
  target_index <- match(preflight$rsid, target$rsid)
  if (anyNA(target_index)) stop("Preflight rsID absent from target union.", call. = FALSE)
  preflight$in_apoe_included_input <- target$in_apoe_included_input[target_index]
  preflight$in_apoe_excluded_input <- target$in_apoe_excluded_input[target_index]
  preflight$instrument_membership <- classify_instrument_membership(preflight$in_apoe_included_input, preflight$in_apoe_excluded_input)
  expected_membership <- classify_instrument_membership(preflight$in_apoe_included_input, preflight$in_apoe_excluded_input)
  membership_rowwise_match <- identical(as.character(preflight$instrument_membership), as.character(expected_membership))
  invalid_membership_count <- sum(preflight$instrument_membership == "invalid_membership")
  shared_rows <- preflight$instrument_membership == "shared"
  included_only_rows <- preflight$instrument_membership == "included_only"
  excluded_only_rows <- preflight$instrument_membership == "excluded_only"
  shared_unique_count <- sum(shared_rows)
  included_only_unique_count <- sum(included_only_rows)
  excluded_only_unique_count <- sum(excluded_only_rows)
  included_unique_count <- sum(preflight$instrument_membership %in% c("shared", "included_only"))
  excluded_unique_count <- sum(preflight$instrument_membership %in% c("shared", "excluded_only"))

  stage <<- "orientation_and_frequency_diagnostics"
  preflight$raw_allele_orientation <- mapply(classify_orientation, preflight$exposure_effect_allele, preflight$exposure_other_allele, preflight$outcome_effect_allele, preflight$outcome_other_allele, USE.NAMES = FALSE)
  preflight$palindromic_snp <- preflight$raw_allele_orientation == "palindromic_snp"
  frequency_rows <- lapply(seq_len(nrow(preflight)), function(i) frequency_diagnostic(preflight$exposure_eaf[[i]], preflight$outcome_eaf[[i]]))
  frequency_df <- as.data.frame(do.call(rbind, lapply(frequency_rows, as.data.frame)), stringsAsFactors = FALSE)
  preflight <- cbind(preflight, frequency_df)
  preflight$harmonisation_performed <- FALSE
  preflight$outcome_beta_flipped <- FALSE
  preflight$record_excluded <- FALSE
  orientation_counts <- as.data.frame(table(preflight$raw_allele_orientation), stringsAsFactors = FALSE)
  names(orientation_counts) <- c("raw_allele_orientation", "count")
  orientation_counts$preflight_version <- "v2"
  palindromic <- preflight[preflight$palindromic_snp, , drop = FALSE]

  stage <<- "hard_checks"
  hard_checks <- list(
    v7_audit_gate = v7_qc$audit_status == "passed" && isTRUE(v7_qc$approved_for_harmonisation_development),
    input_sha_consistency = all(tolower(observed_sha) == tolower(expected_sha)),
    allele_orientation_unit_tests = isTRUE(unit_tests$allele_orientation_unit_tests),
    frequency_diagnostic_unit_tests = isTRUE(unit_tests$frequency_diagnostic_unit_tests),
    membership_unit_test = isTRUE(unit_tests$membership_unit_test),
    mixed_membership_unit_test = isTRUE(unit_tests$mixed_membership_unit_test),
    unique_match_count = nrow(preflight) == 362L,
    unique_rsid_integrity = !anyNA(preflight$rsid) && !anyDuplicated(preflight$rsid),
    multiple_target_excluded = !any(multiple_target_rsid %in% c(preflight$rsid, included_output$rsid, excluded_output$rsid)),
    membership_rowwise_match = membership_rowwise_match,
    no_invalid_membership = invalid_membership_count == 0L,
    shared_unique_count = shared_unique_count == 360L,
    included_only_unique_count = included_only_unique_count == 1L,
    excluded_only_unique_count = excluded_only_unique_count == 1L,
    membership_total = shared_unique_count + included_only_unique_count + excluded_only_unique_count == 362L,
    shared_membership_logic = all(preflight$in_apoe_included_input[shared_rows] & preflight$in_apoe_excluded_input[shared_rows]),
    included_only_membership_logic = all(preflight$in_apoe_included_input[included_only_rows] & !preflight$in_apoe_excluded_input[included_only_rows]),
    excluded_only_membership_logic = all(!preflight$in_apoe_included_input[excluded_only_rows] & preflight$in_apoe_excluded_input[excluded_only_rows]),
    included_unique_count = included_unique_count == 361L,
    excluded_unique_count = excluded_unique_count == 361L,
    included_output_rsid_match = setequal(included_output$rsid, preflight$rsid[preflight$instrument_membership %in% c("shared", "included_only")]),
    excluded_output_rsid_match = setequal(excluded_output$rsid, preflight$rsid[preflight$instrument_membership %in% c("shared", "excluded_only")]),
    palindromic_count = nrow(palindromic) == 49L && !anyDuplicated(palindromic$rsid) && all(palindromic$raw_allele_orientation == "palindromic_snp"),
    preflight_output_integrity = all(c(preflight$harmonisation_performed == FALSE, preflight$outcome_beta_flipped == FALSE, preflight$record_excluded == FALSE)),
    no_harmonisation_performed = !any(preflight$harmonisation_performed),
    no_beta_flips = sum(preflight$outcome_beta_flipped) == 0L,
    no_record_exclusions = sum(preflight$record_excluded) == 0L
  )
  hard_check_failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
  preflight_status <- if (length(hard_check_failures) == 0L) "passed" else "failed"

  stage <<- "write_outputs"
  write_atomic(function(path) write_parquet_file(con, preflight, path), main_output)
  pal_columns <- c("rsid", "exposure_effect_allele", "exposure_other_allele", "outcome_effect_allele", "outcome_other_allele", "exposure_eaf", "outcome_eaf", "exposure_maf", "outcome_maf", "frequency_difference_same", "frequency_difference_flipped", "frequency_orientation_margin", "preferred_frequency_orientation", "instrument_membership")
  pal_export <- palindromic[, pal_columns, drop = FALSE]
  pal_export$alleles <- paste0(pal_export$exposure_effect_allele, "/", pal_export$exposure_other_allele)
  pal_export <- pal_export[, c("rsid", "alleles", pal_columns[6:14]), drop = FALSE]
  write_atomic(function(path) write.table(pal_export, path, sep = "\t", quote = FALSE, row.names = FALSE, na = ""), pal_output)
  write_atomic(function(path) write.csv(orientation_counts, path, row.names = FALSE), orientation_output)
  qc <- list(
    preflight_version = "v2", source_extraction_version = "v9", source_final_audit_version = "v7",
    preflight_status = preflight_status, approved_for_harmonisation_rule_freeze = identical(preflight_status, "passed"),
    unique_match_count = nrow(preflight), included_unique_count = included_unique_count, excluded_unique_count = excluded_unique_count,
    shared_unique_count = shared_unique_count, included_only_unique_count = included_only_unique_count, excluded_only_unique_count = excluded_only_unique_count,
    invalid_membership_count = invalid_membership_count, membership_rowwise_match = membership_rowwise_match,
    orientation_counts = orientation_counts, palindromic_count = nrow(palindromic),
    palindromic_frequency_complete_count = sum(palindromic$preferred_frequency_orientation != "missing_frequency"),
    palindromic_frequency_missing_count = sum(palindromic$preferred_frequency_orientation == "missing_frequency"),
    preferred_same_count = sum(palindromic$preferred_frequency_orientation == "same"),
    preferred_flipped_count = sum(palindromic$preferred_frequency_orientation == "flipped"),
    preferred_tie_count = sum(palindromic$preferred_frequency_orientation == "tie"),
    multiple_target_excluded = !any(multiple_target_rsid %in% c(preflight$rsid, included_output$rsid, excluded_output$rsid)), multiple_target_rsid = multiple_target_rsid,
    harmonisation_performed = FALSE, outcome_beta_flip_count = 0, record_exclusion_count = 0,
    hard_checks = hard_checks, hard_check_failures = hard_check_failures,
    informational_findings = list(input_sha256 = observed_sha, expected_input_sha256 = expected_sha, no_frequency_threshold_applied = TRUE, no_proxy_or_liftover = TRUE)
  )
  write_atomic(function(path) writeLines(jsonlite::toJSON(qc, auto_unbox = TRUE, pretty = TRUE), path), qc_output)
  safe_log("preflight_status=", preflight_status)
  safe_log("hard_check_failures=", paste(hard_check_failures, collapse = ";"))
  safe_log("approved_for_harmonisation_rule_freeze=", qc$approved_for_harmonisation_rule_freeze)
  if (preflight_status == "passed") return(0L)
  1L
}

exit_status <- tryCatch(main(), error = function(error) {
  safe_log("preflight_status=failed")
  safe_log("stage=", stage)
  safe_log("error=", conditionMessage(error))
  1L
})

quit(status = exit_status)
