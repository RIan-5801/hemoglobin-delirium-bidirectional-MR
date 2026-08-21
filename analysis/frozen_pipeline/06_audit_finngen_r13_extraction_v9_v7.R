args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2L) {
  stop("Usage: Rscript 06_audit_finngen_r13_extraction_v9_v7.R --project-root <path>", call. = FALSE)
}

if (args[[1L]] != "--project-root") {
  stop("First argument must be --project-root.", call. = FALSE)
}

root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
stage <- "startup"
qc_dir <- file.path(root, "results", "qc")
out_dir <- file.path(root, "data_derived", "outcome_extracts")
audit_json_path <- file.path(qc_dir, "vuckovic_hb_finngen_r13_extraction_v9_final_audit_v7.json")
manifest_path <- file.path(qc_dir, "vuckovic_hb_finngen_r13_extraction_v9_output_sha_manifest_v4.csv")
audit_log_path <- file.path(root, "results", "logs", "vuckovic_hb_finngen_r13_extraction_v9_final_audit_v7.log")
v9_script_path <- file.path(root, "R", "05_extract_finngen_r13_outcomes_v9.R")

numeric_equal_scalar <- function(
  x,
  y
) {
  valid <- c(
    length(x) == 1L,
    length(y) == 1L,
    !is.na(x),
    !is.na(y),
    is.finite(as.numeric(x)),
    is.finite(as.numeric(y))
  )

  isTRUE(all(c(
    all(valid),
    as.numeric(x) == as.numeric(y)
  )))
}

numeric_vector_equal <- function(
  x,
  y,
  tolerance = 1e-12
) {
  if (length(x) != length(y)) {
    return(FALSE)
  }

  x_numeric <- suppressWarnings(as.numeric(x))
  y_numeric <- suppressWarnings(as.numeric(y))

  same_missing <-
    is.na(x_numeric) &
    is.na(y_numeric)

  both_present <-
    !is.na(x_numeric) &
    !is.na(y_numeric)

  comparisons <- rep(FALSE, length(x_numeric))
  comparisons[same_missing] <- TRUE

  comparisons[both_present] <-
    is.finite(x_numeric[both_present]) &
    is.finite(y_numeric[both_present]) &
    abs(
      x_numeric[both_present] -
        y_numeric[both_present]
    ) <= tolerance

  all(comparisons)
}

character_vector_equal <- function(
  x,
  y
) {
  if (length(x) != length(y)) {
    return(FALSE)
  }

  x_character <- as.character(x)
  y_character <- as.character(y)

  same_missing <-
    is.na(x_character) &
    is.na(y_character)

  both_present <-
    !is.na(x_character) &
    !is.na(y_character)

  comparisons <- rep(FALSE, length(x_character))
  comparisons[same_missing] <- TRUE

  comparisons[both_present] <-
    x_character[both_present] ==
    y_character[both_present]

  all(comparisons)
}

compare_parquet_tsv <- function(
  parquet_df,
  tsv_df
) {
  parquet_df <- parquet_df[order(parquet_df$rsid), , drop = FALSE]
  tsv_df <- tsv_df[order(tsv_df$rsid), , drop = FALSE]

  character_fields <- c(
    "rsid",
    "exposure_effect_allele",
    "exposure_other_allele",
    "outcome_effect_allele",
    "outcome_other_allele",
    "outcome_source",
    "outcome_trait",
    "outcome_build",
    "outcome_effect_scale"
  )
  numeric_fields <- c(
    "exposure_beta",
    "exposure_se",
    "outcome_beta",
    "outcome_se",
    "outcome_pval",
    "outcome_eaf"
  )
  candidate_fields <- union(character_fields, numeric_fields)
  parquet_has_field <- candidate_fields %in% names(parquet_df)
  tsv_has_field <- candidate_fields %in% names(tsv_df)
  one_sided_missing_fields <- candidate_fields[xor(parquet_has_field, tsv_has_field)]
  shared_character_fields <- character_fields[character_fields %in% names(parquet_df)]
  shared_character_fields <- shared_character_fields[shared_character_fields %in% names(tsv_df)]
  shared_numeric_fields <- numeric_fields[numeric_fields %in% names(parquet_df)]
  shared_numeric_fields <- shared_numeric_fields[shared_numeric_fields %in% names(tsv_df)]
  character_mismatch_fields <- shared_character_fields[!vapply(
    shared_character_fields,
    function(field) character_vector_equal(parquet_df[[field]], tsv_df[[field]]),
    logical(1)
  )]
  numeric_mismatch_fields <- shared_numeric_fields[!vapply(
    shared_numeric_fields,
    function(field) numeric_vector_equal(parquet_df[[field]], tsv_df[[field]]),
    logical(1)
  )]
  parquet_duplicate_rsid_count <- sum(duplicated(parquet_df$rsid))
  tsv_duplicate_rsid_count <- sum(duplicated(tsv_df$rsid))
  consistent <- all(c(
    nrow(parquet_df) == nrow(tsv_df),
    setequal(names(parquet_df), names(tsv_df)),
    setequal(parquet_df$rsid, tsv_df$rsid),
    !anyNA(parquet_df$rsid),
    !anyNA(tsv_df$rsid),
    parquet_duplicate_rsid_count == 0L,
    tsv_duplicate_rsid_count == 0L,
    length(character_mismatch_fields) == 0L,
    length(numeric_mismatch_fields) == 0L,
    length(one_sided_missing_fields) == 0L
  ))

  list(
    parquet_row_count = nrow(parquet_df),
    tsv_row_count = nrow(tsv_df),
    column_set_match = setequal(names(parquet_df), names(tsv_df)),
    rsid_set_match = setequal(parquet_df$rsid, tsv_df$rsid),
    parquet_duplicate_rsid_count = parquet_duplicate_rsid_count,
    tsv_duplicate_rsid_count = tsv_duplicate_rsid_count,
    character_mismatch_count = length(character_mismatch_fields),
    character_mismatch_fields = character_mismatch_fields,
    numeric_mismatch_count = length(numeric_mismatch_fields),
    numeric_mismatch_fields = numeric_mismatch_fields,
    one_sided_missing_fields = one_sided_missing_fields,
    consistency = consistent
  )
}

normalize_v9_sha_sequence <- function(
  x,
  verified_key_order
) {
  values <- as.character(unlist(x, use.names = FALSE))
  valid_length <- length(values) == length(verified_key_order)
  valid_values <- all(c(
    !anyNA(values),
    !any(values == ""),
    all(grepl("^[0-9a-fA-F]{64}$", values))
  ))

  if (!isTRUE(valid_length)) {
    stop("V9 QC SHA sequence length does not match verified key order.", call. = FALSE)
  }

  if (!isTRUE(valid_values)) {
    stop("V9 QC SHA sequence contains invalid SHA-256 values.", call. = FALSE)
  }

  setNames(tolower(values), verified_key_order)
}

read_parquet_file <- function(
  con,
  path
) {
  normalized_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", normalized_path))
}

sha256_file <- function(
  path
) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

add_hard_check <- function(
  checks,
  name,
  value
) {
  if (length(value) != 1L) {
    stop(paste("Hard check did not resolve to one value:", name), call. = FALSE)
  }

  checks[[name]] <- isTRUE(value)
  checks
}

safe_log <- function(
  ...
) {
  cat(paste0(...), "\n", file = audit_log_path, append = TRUE)
}

find_assignments_recursive <- function(
  expression_list,
  target_name
) {
  found <- list()

  walk <- function(node) {
    is_assignment <- FALSE
    if (is.call(node)) {
      if (length(node) >= 3L) {
        is_assignment <- isTRUE(all(c(
          identical(node[[1L]], as.name("<-")),
          is.symbol(node[[2L]]),
          identical(as.character(node[[2L]]), target_name)
        )))
      }
    }

    if (isTRUE(is_assignment)) {
      found[[length(found) + 1L]] <<- node
    }

    is_composite <- any(c(is.call(node), is.expression(node), is.pairlist(node)))
    if (isTRUE(is_composite)) {
      for (index in seq_along(node)) {
        walk(node[[index]])
      }
    }

    invisible(NULL)
  }

  walk(expression_list)
  found
}

write_atomic_text <- function(
  text,
  final_path
) {
  partial_path <- paste0(final_path, ".partial")
  on.exit({
    if (file.exists(partial_path)) {
      unlink(partial_path)
    }
  }, add = TRUE)
  writeLines(text, partial_path, useBytes = TRUE)
  moved <- file.rename(partial_path, final_path)
  if (!isTRUE(moved)) {
    stop(paste("Atomic rename failed:", final_path), call. = FALSE)
  }
}

run_comparison_unit_tests <- function() {
  n <- 361L
  tolerance <- 1e-12
  base_numeric <- rep(0, n)
  base_numeric[c(7L, 222L)] <- NA_real_
  within_numeric <- base_numeric
  within_numeric[[100L]] <- tolerance / 2
  outside_numeric <- base_numeric
  outside_numeric[[100L]] <- tolerance * 2
  actual_within_difference <- abs(within_numeric[[100L]] - base_numeric[[100L]])
  actual_outside_difference <- abs(outside_numeric[[100L]] - base_numeric[[100L]])
  character_a <- paste0("rs", seq_len(n))
  character_a[c(5L, 333L)] <- NA_character_
  character_b <- character_a
  character_mismatch <- character_a
  character_mismatch[[12L]] <- "different"
  character_one_sided_na <- character_a
  character_one_sided_na[[5L]] <- "rs5"
  data_a <- data.frame(
    rsid = paste0("rs", seq_len(n)),
    exposure_effect_allele = rep("A", n),
    exposure_other_allele = rep("G", n),
    exposure_beta = base_numeric,
    exposure_se = rep(0.01, n),
    outcome_effect_allele = rep("C", n),
    outcome_other_allele = rep("T", n),
    outcome_beta = base_numeric / 10,
    outcome_se = rep(0.02, n),
    outcome_pval = rep(0.05, n),
    outcome_eaf = rep(0.2, n),
    stringsAsFactors = FALSE
  )
  data_numeric_mismatch <- data_a
  data_numeric_mismatch$outcome_beta[[101L]] <- tolerance * 2
  data_character_mismatch <- data_a
  data_character_mismatch$outcome_effect_allele[[102L]] <- "G"
  numeric_mismatch <- compare_parquet_tsv(data_a, data_numeric_mismatch)
  character_mismatch_result <- compare_parquet_tsv(data_a, data_character_mismatch)

  checks <- list(
    comparison_unit_test_length = n,
    comparison_tolerance = tolerance,
    actual_within_difference = actual_within_difference,
    actual_outside_difference = actual_outside_difference,
    within_difference_less_than_tolerance = isTRUE(actual_within_difference < tolerance),
    outside_difference_greater_than_tolerance = isTRUE(actual_outside_difference > tolerance),
    numeric_equal_test = numeric_vector_equal(base_numeric, base_numeric, tolerance),
    numeric_tolerance_within_test = numeric_vector_equal(base_numeric, within_numeric, tolerance),
    numeric_tolerance_outside_test = !numeric_vector_equal(base_numeric, outside_numeric, tolerance),
    numeric_single_sided_na_test = !numeric_vector_equal(base_numeric, replace(base_numeric, 7L, 1), tolerance),
    numeric_length_mismatch_test = !numeric_vector_equal(base_numeric, base_numeric[-1L], tolerance),
    numeric_infinite_test = !numeric_vector_equal(c(Inf), c(Inf), tolerance),
    numeric_nan_both_missing_test = numeric_vector_equal(c(NaN), c(NaN), tolerance),
    character_equal_test = character_vector_equal(character_a, character_b),
    character_mismatch_test = !character_vector_equal(character_a, character_mismatch),
    character_single_sided_na_test = !character_vector_equal(character_a, character_one_sided_na),
    character_length_mismatch_test = !character_vector_equal(character_a, character_a[-1L]),
    data_frame_equal_test = compare_parquet_tsv(data_a, data_a)$consistency,
    data_frame_numeric_mismatch_test = isTRUE(all(c(
      numeric_mismatch$numeric_mismatch_count == 1L,
      "outcome_beta" %in% numeric_mismatch$numeric_mismatch_fields,
      !numeric_mismatch$consistency
    ))),
    data_frame_character_mismatch_test = isTRUE(all(c(
      character_mismatch_result$character_mismatch_count == 1L,
      "outcome_effect_allele" %in% character_mismatch_result$character_mismatch_fields,
      !character_mismatch_result$consistency
    )))
  )
  boolean_check_names <- names(checks)[vapply(checks, is.logical, logical(1))]
  checks$vector_comparison_unit_tests <- all(vapply(checks[boolean_check_names], isTRUE, logical(1)))
  checks
}

main <- function() {
  stage <<- "preflight"
  target_paths <- c(audit_json_path, manifest_path, audit_log_path)
  target_partials <- paste0(target_paths, ".partial")
  occupied <- c(target_paths, target_partials)[file.exists(c(target_paths, target_partials))]
  if (length(occupied) > 0L) {
    stop("V7 audit target or partial target already exists; refusing to overwrite.", call. = FALSE)
  }

  stage <<- "comparison_unit_tests"
  comparison_unit_tests <- run_comparison_unit_tests()
  if (!isTRUE(comparison_unit_tests$vector_comparison_unit_tests)) {
    stop("Vector comparison unit tests failed before V9 artifacts were read.", call. = FALSE)
  }
  for (name in names(comparison_unit_tests)) {
    value <- comparison_unit_tests[[name]]
    rendered <- if (is.logical(value)) if (isTRUE(value)) "passed" else "failed" else as.character(value)
    safe_log(name, "=", rendered)
  }

  stage <<- "v9_output_preflight"
  paths <- c(
    target_union = file.path(out_dir, "vuckovic_hb_target_union_v9.parquet"),
    union_matches = file.path(out_dir, "vuckovic_hb_finngen_r13_union_matches_v9.parquet"),
    inc_pq = file.path(out_dir, "vuckovic_hb_finngen_r13_apoe_included_v9.parquet"),
    inc_tsv = file.path(out_dir, "vuckovic_hb_finngen_r13_apoe_included_v9.tsv"),
    exc_pq = file.path(out_dir, "vuckovic_hb_finngen_r13_apoe_excluded_v9.parquet"),
    exc_tsv = file.path(out_dir, "vuckovic_hb_finngen_r13_apoe_excluded_v9.tsv"),
    missing = file.path(out_dir, "audit", "vuckovic_hb_finngen_r13_missing_v9.parquet"),
    multiple = file.path(out_dir, "audit", "vuckovic_hb_finngen_r13_multiple_matches_v9.parquet"),
    allele = file.path(out_dir, "audit", "vuckovic_hb_finngen_r13_allele_audit_v9.parquet"),
    qc = file.path(qc_dir, "vuckovic_hb_finngen_r13_extraction_v9.json"),
    match_csv = file.path(qc_dir, "vuckovic_hb_finngen_r13_match_counts_v9.csv"),
    allele_csv = file.path(qc_dir, "vuckovic_hb_finngen_r13_allele_counts_v9.csv"),
    log = file.path(root, "results", "logs", "vuckovic_hb_finngen_r13_extraction_v9.log")
  )
  missing_output_files <- unname(paths[!file.exists(paths)])
  existing_paths <- paths[file.exists(paths)]
  zero_byte_output_files <- unname(existing_paths[file.info(existing_paths)$size <= 0])
  partial_files <- unique(c(
    list.files(out_dir, pattern = "\\.partial$", recursive = TRUE, full.names = TRUE),
    list.files(qc_dir, pattern = "v9.*\\.partial$", full.names = TRUE)
  ))
  if (length(missing_output_files) > 0L) {
    stop("One or more required V9 artifacts are missing.", call. = FALSE)
  }
  if (length(zero_byte_output_files) > 0L) {
    stop("One or more required V9 artifacts are zero bytes.", call. = FALSE)
  }

  stage <<- "v9_output_key_order"
  parsed_v9 <- parse(file = v9_script_path)
  out_assignments <- find_assignments_recursive(parsed_v9, "out")
  v9_out_definition_count <- length(out_assignments)
  if (v9_out_definition_count != 1L) {
    stop("Expected exactly one out assignment in V9 source.", call. = FALSE)
  }
  audit_env <- new.env(parent = globalenv())
  audit_env$root <- "/AUDIT_DUMMY_ROOT"
  eval(out_assignments[[1L]], envir = audit_env)
  v9_out_key_order <- names(audit_env$out)
  expected_v9_out_key_order <- c(
    "target_union", "union_matches", "inc_pq", "inc_tsv", "exc_pq", "exc_tsv",
    "missing", "multiple", "allele", "qc", "match_csv", "allele_csv", "log"
  )
  v9_out_key_order_verified <- identical(v9_out_key_order, expected_v9_out_key_order)
  if (!isTRUE(v9_out_key_order_verified)) {
    stop("V9 output key order was not verified from source.", call. = FALSE)
  }
  expected_sha_key_set <- c(
    "target_union", "union_matches", "inc_pq", "inc_tsv", "exc_pq", "exc_tsv",
    "missing", "multiple", "allele", "match_csv", "allele_csv"
  )
  v9_sha_key_order <- v9_out_key_order[v9_out_key_order %in% expected_sha_key_set]
  sha_key_order_verified <- isTRUE(all(c(
    length(v9_sha_key_order) == 11L,
    setequal(v9_sha_key_order, expected_sha_key_set),
    identical(v9_sha_key_order, expected_sha_key_set)
  )))
  if (!isTRUE(sha_key_order_verified)) {
    stop("V9 SHA key order was not verified from V9 output order.", call. = FALSE)
  }

  stage <<- "qc_sha_mapping"
  qc <- jsonlite::fromJSON(paths[["qc"]], simplifyVector = FALSE)
  raw_sha_sequence <- qc$output_sha256_excluding_qc_json_and_log
  qc_sha_raw_class <- paste(class(raw_sha_sequence), collapse = ";")
  qc_sha_raw_type <- typeof(raw_sha_sequence)
  qc_sha_raw_length <- length(raw_sha_sequence)
  qc_sha_raw_names_present <- FALSE
  if (!is.null(names(raw_sha_sequence))) {
    qc_sha_raw_names_present <- any(nzchar(names(raw_sha_sequence)))
  }
  sha_map <- normalize_v9_sha_sequence(raw_sha_sequence, v9_sha_key_order)

  stage <<- "output_sha_verification"
  observed_sha_map <- vapply(
    v9_sha_key_order,
    function(key) sha256_file(paths[[key]]),
    character(1)
  )
  names_match <- identical(names(observed_sha_map), names(sha_map))
  positional_sha_match <- identical(tolower(unname(observed_sha_map)), tolower(unname(sha_map)))
  unordered_sha_set_match <- setequal(tolower(unname(observed_sha_map)), tolower(unname(sha_map)))
  sha_match <- tolower(observed_sha_map) == tolower(sha_map)
  sha_details <- data.frame(
    key = v9_sha_key_order,
    file = unname(paths[v9_sha_key_order]),
    sequence_position = seq_along(v9_sha_key_order),
    expected_sha256 = unname(sha_map),
    observed_sha256 = unname(observed_sha_map),
    sha_match = unname(sha_match),
    stringsAsFactors = FALSE
  )

  stage <<- "small_output_read"
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  target <- read_parquet_file(con, paths[["target_union"]])
  union_matches <- read_parquet_file(con, paths[["union_matches"]])
  included <- read_parquet_file(con, paths[["inc_pq"]])
  excluded <- read_parquet_file(con, paths[["exc_pq"]])
  missing <- read_parquet_file(con, paths[["missing"]])
  multiple <- read_parquet_file(con, paths[["multiple"]])
  allele <- read_parquet_file(con, paths[["allele"]])
  included_tsv <- read.delim(paths[["inc_tsv"]], check.names = FALSE, stringsAsFactors = FALSE)
  excluded_tsv <- read.delim(paths[["exc_tsv"]], check.names = FALSE, stringsAsFactors = FALSE)
  allele_counts_csv <- read.csv(paths[["allele_csv"]], stringsAsFactors = FALSE)
  v9_log <- readLines(paths[["log"]], warn = FALSE)

  stage <<- "small_output_checks"
  included_comparison <- compare_parquet_tsv(included, included_tsv)
  excluded_comparison <- compare_parquet_tsv(excluded, excluded_tsv)
  status_counts <- table(target$finngen_match_status)
  count_status <- function(status) {
    if (status %in% names(status_counts)) {
      return(as.numeric(status_counts[[status]]))
    }
    0
  }
  unique_match_count <- count_status("unique_match")
  missing_match_count <- count_status("missing_in_finngen")
  multiple_match_count <- count_status("multiple_finngen_matches")
  included_target_rsids <- target$rsid[target$in_apoe_included_input]
  excluded_target_rsids <- target$rsid[target$in_apoe_excluded_input]
  shared_target_rsids <- intersect(included_target_rsids, excluded_target_rsids)
  unique_target_rsids <- target$rsid[target$finngen_match_status == "unique_match"]
  included_unique_match_count <- sum(included_target_rsids %in% unique_target_rsids)
  excluded_unique_match_count <- sum(excluded_target_rsids %in% unique_target_rsids)
  shared_unique_match_count <- sum(shared_target_rsids %in% unique_target_rsids)
  included_only_unique_match_count <- sum(setdiff(included_target_rsids, excluded_target_rsids) %in% unique_target_rsids)
  excluded_only_unique_match_count <- sum(setdiff(excluded_target_rsids, included_target_rsids) %in% unique_target_rsids)
  multiple_rsids <- sort(unique(multiple$rsid))
  missing_rsids <- sort(unique(missing$rsid))
  union_match_rsid_column <- if ("rsid" %in% names(union_matches)) "rsid" else "matched_rsid_token"
  union_matches_are_targets <- FALSE
  if (union_match_rsid_column %in% names(union_matches)) {
    union_matches_are_targets <- all(union_matches[[union_match_rsid_column]] %in% target$rsid)
  }
  actual_allele_counts <- as.data.frame(table(allele$allele_audit_status), stringsAsFactors = FALSE)
  names(actual_allele_counts) <- c("allele_audit_status", "count")
  actual_allele_counts <- actual_allele_counts[order(actual_allele_counts$allele_audit_status), , drop = FALSE]
  allele_counts_csv <- allele_counts_csv[order(allele_counts_csv$allele_audit_status), , drop = FALSE]
  allele_counts_consistency <- isTRUE(all(c(
    identical(as.character(actual_allele_counts$allele_audit_status), as.character(allele_counts_csv$allele_audit_status)),
    numeric_vector_equal(actual_allele_counts$count, allele_counts_csv$count)
  )))
  all_outcomes <- rbind(included, excluded)
  eaf_valid <- is.na(as.numeric(all_outcomes$outcome_eaf)) |
    (is.finite(as.numeric(all_outcomes$outcome_eaf)) & as.numeric(all_outcomes$outcome_eaf) >= 0 & as.numeric(all_outcomes$outcome_eaf) <= 1)
  metadata_valid <- all(c(
    all(all_outcomes$outcome_source == "FinnGen_R13_F5_DELIRIUM"),
    all(all_outcomes$outcome_trait == "delirium"),
    all(all_outcomes$outcome_build == "GRCh38"),
    all(as.numeric(all_outcomes$outcome_ncase) == 5121),
    all(as.numeric(all_outcomes$outcome_ncontrol) == 465023),
    all(as.numeric(all_outcomes$outcome_n_study) == 470144)
  ))
  outcome_field_quality <- all(c(
    metadata_valid,
    all(is.finite(as.numeric(all_outcomes$outcome_beta))),
    all(is.finite(as.numeric(all_outcomes$outcome_se))),
    all(as.numeric(all_outcomes$outcome_se) > 0),
    all(is.finite(as.numeric(all_outcomes$outcome_pval))),
    all(as.numeric(all_outcomes$outcome_pval) >= 0),
    all(as.numeric(all_outcomes$outcome_pval) <= 1),
    all(eaf_valid)
  ))
  outcome_directory_names <- list.files(out_dir, recursive = TRUE, full.names = FALSE)
  no_harmonisation_proxy_or_flip_results <- !any(grepl("harmoni[sz]|proxy|effect[_ ]?flip|beta[_ ]?flip", outcome_directory_names, ignore.case = TRUE))
  frozen_finngen_sha <- "85637f0f3358807964d4f8a3e500293168a706f1c08c65f3fc5512b65df40ed8"
  input_sha_consistency <- all(tolower(c(
    qc$input_sha256$finngen_before,
    qc$input_sha256$finngen_after,
    qc$certified_input_sha256
  )) == frozen_finngen_sha)
  log_success <- any(grepl("SUCCESS:", v9_log, fixed = TRUE))
  no_v9_terminated_primary <- !any(grepl("TERMINATED_PRIMARY", v9_log, fixed = TRUE))
  no_v9_cleanup_warning <- !any(grepl("CLEANUP_WARNING", v9_log, fixed = TRUE))

  stage <<- "hard_checks"
  hard_checks <- list()
  hard_checks <- add_hard_check(hard_checks, "comparison_unit_tests", comparison_unit_tests$vector_comparison_unit_tests)
  hard_checks <- add_hard_check(hard_checks, "output_files_complete", length(missing_output_files) == 0L)
  hard_checks <- add_hard_check(hard_checks, "no_zero_byte_outputs", length(zero_byte_output_files) == 0L)
  hard_checks <- add_hard_check(hard_checks, "no_partial_files", length(partial_files) == 0L)
  hard_checks <- add_hard_check(hard_checks, "v9_log_success", log_success)
  hard_checks <- add_hard_check(hard_checks, "no_v9_terminated_primary", no_v9_terminated_primary)
  hard_checks <- add_hard_check(hard_checks, "no_v9_cleanup_warning", no_v9_cleanup_warning)
  hard_checks <- add_hard_check(hard_checks, "certified_row_count", numeric_equal_scalar(qc$certified_data_rows, 21326687))
  hard_checks <- add_hard_check(hard_checks, "observed_row_count", numeric_equal_scalar(qc$observed_data_rows, 21326687))
  hard_checks <- add_hard_check(hard_checks, "row_counts_equal", numeric_equal_scalar(qc$certified_data_rows, qc$observed_data_rows))
  hard_checks <- add_hard_check(hard_checks, "input_sha_consistency", input_sha_consistency)
  hard_checks <- add_hard_check(hard_checks, "target_union_count", nrow(target) == 395L)
  hard_checks <- add_hard_check(hard_checks, "included_target_count", length(included_target_rsids) == 393L)
  hard_checks <- add_hard_check(hard_checks, "excluded_target_count", length(excluded_target_rsids) == 394L)
  hard_checks <- add_hard_check(hard_checks, "shared_target_count", length(shared_target_rsids) == 392L)
  hard_checks <- add_hard_check(hard_checks, "match_classification_total", unique_match_count + missing_match_count + multiple_match_count == 395L)
  hard_checks <- add_hard_check(hard_checks, "unique_match_count", unique_match_count == 362L)
  hard_checks <- add_hard_check(hard_checks, "missing_match_count", missing_match_count == 32L)
  hard_checks <- add_hard_check(hard_checks, "multiple_match_count", multiple_match_count == 1L)
  hard_checks <- add_hard_check(hard_checks, "missing_audit_count", length(missing_rsids) == 32L)
  hard_checks <- add_hard_check(hard_checks, "multiple_audit_count", length(multiple_rsids) == 1L)
  hard_checks <- add_hard_check(hard_checks, "allele_audit_count", all(c(nrow(allele) == 362L, !anyDuplicated(allele$rsid), !any(allele$rsid %in% c(missing$rsid, multiple$rsid)))))
  hard_checks <- add_hard_check(hard_checks, "union_matches_target_membership", union_matches_are_targets)
  hard_checks <- add_hard_check(hard_checks, "included_output_membership", all(c(all(included$rsid %in% intersect(included_target_rsids, unique_target_rsids)), !anyDuplicated(included$rsid))))
  hard_checks <- add_hard_check(hard_checks, "excluded_output_membership", all(c(all(excluded$rsid %in% intersect(excluded_target_rsids, unique_target_rsids)), !anyDuplicated(excluded$rsid))))
  hard_checks <- add_hard_check(hard_checks, "included_parquet_tsv", included_comparison$consistency)
  hard_checks <- add_hard_check(hard_checks, "excluded_parquet_tsv", excluded_comparison$consistency)
  hard_checks <- add_hard_check(hard_checks, "allele_counts_consistency", allele_counts_consistency)
  hard_checks <- add_hard_check(hard_checks, "outcome_field_quality", outcome_field_quality)
  hard_checks <- add_hard_check(hard_checks, "returned_non_target_rows", numeric_equal_scalar(qc$returned_non_target_detail_rows, 0))
  hard_checks <- add_hard_check(hard_checks, "v9_out_key_order", v9_out_key_order_verified)
  hard_checks <- add_hard_check(hard_checks, "output_sha_keys", all(c(names_match, sha_key_order_verified)))
  hard_checks <- add_hard_check(hard_checks, "output_sha_values", all(c(all(sha_match), positional_sha_match, unordered_sha_set_match)))
  hard_checks <- add_hard_check(hard_checks, "multiple_not_selected", all(c(
    !any(multiple_rsids %in% included$rsid),
    !any(multiple_rsids %in% excluded$rsid)
  )))
  hard_checks <- add_hard_check(hard_checks, "no_harmonisation_proxy_or_flip_results", no_harmonisation_proxy_or_flip_results)
  hard_check_failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
  audit_status <- if (length(hard_check_failures) == 0L) "passed" else "failed"

  stage <<- "write_outputs"
  manifest_partial <- paste0(manifest_path, ".partial")
  write.csv(sha_details, manifest_partial, row.names = FALSE, quote = TRUE)
  if (!file.rename(manifest_partial, manifest_path)) {
    stop("Could not atomically publish SHA manifest.", call. = FALSE)
  }
  multiple_match_diagnostics <- list(
    rsids = multiple_rsids,
    candidate_rows = nrow(multiple),
    distinct_rsids = length(multiple_rsids),
    present_in_allele_audit = sum(multiple_rsids %in% allele$rsid),
    present_in_included_output = sum(multiple_rsids %in% included$rsid),
    present_in_excluded_output = sum(multiple_rsids %in% excluded$rsid),
    automatic_selection_performed = FALSE
  )
  audit <- list(
    audit_version = "v7",
    production_version = "v9",
    audit_status = audit_status,
    approved_for_harmonisation_development = identical(audit_status, "passed"),
    comparison_unit_tests = comparison_unit_tests,
    output_file_count = length(paths),
    missing_output_files = missing_output_files,
    zero_byte_output_files = zero_byte_output_files,
    partial_file_count = length(partial_files),
    partial_files = partial_files,
    log_success = log_success,
    certified_data_rows = as.numeric(qc$certified_data_rows),
    observed_data_rows = as.numeric(qc$observed_data_rows),
    row_count_numeric_match = numeric_equal_scalar(qc$certified_data_rows, qc$observed_data_rows),
    input_sha_consistency = input_sha_consistency,
    target_union_count = nrow(target),
    included_target_count = length(included_target_rsids),
    excluded_target_count = length(excluded_target_rsids),
    shared_target_count = length(shared_target_rsids),
    match_counts = list(unique_match = unique_match_count, missing_in_finngen = missing_match_count, multiple_finngen_matches = multiple_match_count),
    included_match_counts = list(unique_match = included_unique_match_count, missing = sum(included_target_rsids %in% missing_rsids), multiple = sum(included_target_rsids %in% multiple_rsids)),
    excluded_match_counts = list(unique_match = excluded_unique_match_count, missing = sum(excluded_target_rsids %in% missing_rsids), multiple = sum(excluded_target_rsids %in% multiple_rsids)),
    shared_match_counts = list(unique_match = shared_unique_match_count, included_only_unique_match = included_only_unique_match_count, excluded_only_unique_match = excluded_only_unique_match_count),
    multiple_match_diagnostics = multiple_match_diagnostics,
    allele_audit_counts = actual_allele_counts,
    returned_non_target_detail_rows = as.numeric(qc$returned_non_target_detail_rows),
    included_parquet_tsv_comparison = included_comparison,
    excluded_parquet_tsv_comparison = excluded_comparison,
    v9_out_definition_count = v9_out_definition_count,
    v9_out_key_order = v9_out_key_order,
    v9_out_key_order_verified = v9_out_key_order_verified,
    qc_sha_raw_class = qc_sha_raw_class,
    qc_sha_raw_type = qc_sha_raw_type,
    qc_sha_raw_length = qc_sha_raw_length,
    qc_sha_raw_names_present = qc_sha_raw_names_present,
    v9_qc_sha_names_preserved = FALSE,
    v9_qc_sha_serialization_limitation = "output SHA values were stored as an unnamed JSON array",
    sha_mapping_method = "positional_mapping_from_verified_v9_out_key_order",
    output_sha_expected_key_count = length(sha_map),
    output_sha_verified_count = sum(sha_match),
    output_sha_mismatch_count = sum(!sha_match),
    output_sha_missing_key_count = sum(!(v9_sha_key_order %in% names(sha_map))),
    output_sha_extra_key_count = sum(!(names(sha_map) %in% v9_sha_key_order)),
    positional_sha_match = positional_sha_match,
    unordered_sha_set_match = unordered_sha_set_match,
    output_sha_details = sha_details,
    hard_checks = hard_checks,
    hard_check_failures = hard_check_failures,
    informational_findings = list(
      noncanonical_rsids_rows = as.numeric(qc$rsids_audit$noncanonical_rsids_rows),
      palindromic_snp = sum(allele$allele_audit_status == "palindromic_snp"),
      future_qc_sha_serialization_recommendation = "serialize output hashes as a named list or explicit key-sha table"
    )
  )
  write_atomic_text(jsonlite::toJSON(audit, auto_unbox = TRUE, pretty = TRUE, null = "null"), audit_json_path)
  safe_log("audit_status=", audit_status)
  safe_log("hard_check_failures=", paste(hard_check_failures, collapse = ";"))
  safe_log("approved_for_harmonisation_development=", audit$approved_for_harmonisation_development)

  if (audit_status == "passed") {
    return(0L)
  }

  1L
}

exit_status <- tryCatch(
  main(),
  error = function(error) {
    safe_log("audit_status=failed")
    safe_log("stage=", if (exists("stage", inherits = TRUE)) stage else "unavailable")
    safe_log("error=", conditionMessage(error))
    1L
  }
)

quit(status = exit_status)
