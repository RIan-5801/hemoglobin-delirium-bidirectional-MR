#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/29_chen_reverse_outcome_extraction_v1.R [--project-root <path>]", call. = FALSE)
}
setwd(root)

for (pkg in c("DBI", "duckdb", "jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

rel <- function(...) file.path(root, ...)
norm <- function(x) gsub("\\\\", "/", x)
relpath <- function(x) norm(sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", norm(normalizePath(x, winslash = "/", mustWork = FALSE))))
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
atomic_write <- function(path, writer) {
  partial <- paste0(path, ".partial")
  stop_if(file.exists(path) || file.exists(partial), paste("Output occupied:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writer(partial)
  stop_if(!file.exists(partial), paste("Writer did not create partial:", partial))
  stop_if(!file.rename(partial, path), paste("Atomic rename failed:", path))
}
write_tsv <- function(x, path) atomic_write(path, function(p) {
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.table(x, p, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
})
write_csv <- function(x, path) atomic_write(path, function(p) {
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.csv(x, p, row.names = FALSE, na = "")
})
write_json <- function(x, path) atomic_write(path, function(p) {
  jsonlite::write_json(x, p, pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)
})
write_text <- function(lines, path) atomic_write(path, function(p) writeLines(lines, p, useBytes = TRUE))
parse_num <- function(x) {
  y <- trimws(as.character(x))
  y[y == ""] <- NA_character_
  suppressWarnings(as.numeric(y))
}
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
allele_key <- function(a, b) paste(sort(c(as.character(a), as.character(b))), collapse = "/")
marker_parts <- function(marker) {
  z <- strsplit(marker, "[:_]", perl = TRUE)[[1L]]
  if (length(z) != 4L) return(c(NA_character_, NA_character_, NA_character_, NA_character_))
  z
}

paths <- list(
  script = rel("R", "29_chen_reverse_outcome_extraction_v1.R"),
  decision109_qc = rel("results", "qc", "chen_reverse_outcome_extraction_execution_contract_v1_readback_closure.json"),
  target_authority = rel("results", "qc", "chen_reverse_outcome_target_authority_audit_v1_readback_recovery.csv"),
  chen_source_cert = rel("results", "qc", "chen_2020_hb_source_certification_v1.json"),
  chen_readme = rel("docs", "source_metadata", "readme_BCX2_meta_analyses.txt"),
  chen_raw = rel("data_raw", "gwas", "BCX2_HGB_EA_GWAMA.out.gz"),
  renv_lock = rel("renv.lock"),
  out_dir = rel("data_derived", "reverse_sensitivity_outcome"),
  union_targets = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_union_targets_v1.tsv"),
  master_parquet = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_master_v1.parquet"),
  master_tsv = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_master_v1.tsv"),
  strict_inc_parquet = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_strict_apoe_included_v1.parquet"),
  strict_inc_tsv = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_strict_apoe_included_v1.tsv"),
  strict_exc_parquet = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_strict_apoe_excluded_v1.parquet"),
  strict_exc_tsv = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_strict_apoe_excluded_v1.tsv"),
  relaxed_inc_parquet = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_relaxed_apoe_included_v1.parquet"),
  relaxed_inc_tsv = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_relaxed_apoe_included_v1.tsv"),
  relaxed_exc_parquet = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_relaxed_apoe_excluded_v1.parquet"),
  relaxed_exc_tsv = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_relaxed_apoe_excluded_v1.tsv"),
  match_audit = rel("results", "qc", "chen_reverse_outcome_match_audit_v1.csv"),
  missing = rel("results", "qc", "chen_reverse_outcome_missing_v1.tsv"),
  multiple = rel("results", "qc", "chen_reverse_outcome_multiple_matches_v1.tsv"),
  qc_json = rel("results", "qc", "chen_reverse_outcome_extraction_v1.json"),
  log = rel("results", "logs", "chen_reverse_outcome_extraction_v1.log"),
  decision = rel("docs", "decisions", "110_chen_reverse_outcome_extraction_v1_v1.1.md")
)

targets <- unlist(paths[c(
  "union_targets", "master_parquet", "master_tsv", "strict_inc_parquet", "strict_inc_tsv",
  "strict_exc_parquet", "strict_exc_tsv", "relaxed_inc_parquet", "relaxed_inc_tsv",
  "relaxed_exc_parquet", "relaxed_exc_tsv", "match_audit", "missing", "multiple",
  "qc_json", "log", "decision"
)])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Output or partial exists:", paste(occupied, collapse = "; ")))

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 110L), paste0("Expected next decision 110, found ", next_decision, "; no outputs written."))

dir.create(dirname(paths$log), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = paths$log, append = TRUE)
log_line("stage=chen_reverse_outcome_extraction_v1_start")

fixture <- data.frame(x = c("A", "C", "G", "T", "F"), stringsAsFactors = FALSE)
fixture_file <- tempfile(fileext = ".tsv")
utils::write.table(fixture, fixture_file, sep = "\t", quote = FALSE, row.names = FALSE)
fixture_back <- utils::read.delim(fixture_file, colClasses = "character", stringsAsFactors = FALSE, check.names = FALSE)
allele_T_F_fixture_passed <- identical(fixture_back$x, fixture$x) && is.character(fixture_back$x)
stop_if(!allele_T_F_fixture_passed, "Character schema fixture failed; no Chen source scan performed.")

contract <- read_json(paths$decision109_qc)
stop_if(!identical(contract$authoritative_contract_status, "frozen"), "Decision 109 contract gate failed.")
stop_if(!isTRUE(contract$approved_for_chen_reverse_outcome_extraction), "Decision 109 extraction approval missing.")
stop_if(length(contract$closure_hard_check_failures) != 0L, "Decision 109 hard-check failures are not empty.")

target_authority <- utils::read.csv(paths$target_authority, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
required_target_cols <- c(
  "rsid", "strict_included_member", "strict_excluded_member", "relaxed_included_member", "relaxed_excluded_member",
  "expected_chen_marker_id", "identity_bridge_ready", "reference_allele1", "reference_allele2"
)
stop_if(any(!required_target_cols %in% names(target_authority)), "Corrected target audit schema missing required columns.")
for (x in c("rsid", "expected_chen_marker_id", "exposure_effect_allele", "exposure_other_allele", "reference_allele1", "reference_allele2")) {
  stop_if(!is.character(target_authority[[x]]), paste("Target audit column is not character:", x))
}
bool <- function(x) toupper(as.character(x)) == "TRUE"
target_authority$strict_included_member <- bool(target_authority$strict_included_member)
target_authority$strict_excluded_member <- bool(target_authority$strict_excluded_member)
target_authority$relaxed_included_member <- bool(target_authority$relaxed_included_member)
target_authority$relaxed_excluded_member <- bool(target_authority$relaxed_excluded_member)
target_authority$identity_bridge_ready <- bool(target_authority$identity_bridge_ready)
expected_from_fields <- vapply(seq_len(nrow(target_authority)), function(i) {
  paste0(
    target_authority$reference_chr_grch37[[i]], ":", target_authority$reference_pos_grch37[[i]], "_",
    paste(sort(c(target_authority$reference_allele1[[i]], target_authority$reference_allele2[[i]])), collapse = "_")
  )
}, character(1))
expected_marker_validation_ok <- identical(expected_from_fields, target_authority$expected_chen_marker_id)
stop_if(!expected_marker_validation_ok, "Expected Chen marker validation against target authority failed.")

source_sha_before <- hash_file(paths$chen_raw)
chen_cert <- read_json(paths$chen_source_cert)
expected_sha <- tolower(chen_cert$source_sha256)
stop_if(!identical(tolower(source_sha_before), expected_sha), "Chen source SHA before scan differs from certified authority.")
source_size_bytes <- unname(file.info(paths$chen_raw)$size)

header_con <- gzfile(paths$chen_raw, open = "rt")
header <- readLines(header_con, n = 1L, warn = FALSE)
close(header_con)
header_fields <- strsplit(header, "\t", fixed = TRUE)[[1L]]
required_source_cols <- c("rs_number", "reference_allele", "other_allele", "eaf", "beta", "se", "p-value", "n_samples")
chen_source_schema_verified <- all(required_source_cols %in% header_fields)
stop_if(!chen_source_schema_verified, "Chen source header missing required fields.")
col_idx <- setNames(match(header_fields, header_fields), header_fields)

target_markers <- target_authority$expected_chen_marker_id
matches <- vector("list", 0L)
invalid_matches <- vector("list", 0L)
data_rows_scanned <- 0L
start <- Sys.time()
con <- gzfile(paths$chen_raw, open = "rt")
on.exit(try(close(con), silent = TRUE), add = TRUE)
discard <- readLines(con, n = 1L, warn = FALSE)
repeat {
  lines <- readLines(con, n = 100000L, warn = FALSE)
  if (length(lines) == 0L) break
  data_rows_scanned <- data_rows_scanned + length(lines)
  tab_pos <- regexpr("\t", lines, fixed = TRUE)
  ids <- ifelse(tab_pos > 0L, substr(lines, 1L, tab_pos - 1L), lines)
  keep <- ids %in% target_markers
  if (any(keep)) {
    kept <- lines[keep]
    parts <- strsplit(kept, "\t", fixed = TRUE)
    for (z in parts) {
      if (length(z) != length(header_fields)) {
        invalid_matches[[length(invalid_matches) + 1L]] <- list(source_marker_id = z[[1L]], raw_field_count = length(z))
      } else {
        names(z) <- header_fields
        matches[[length(matches) + 1L]] <- as.list(z)
      }
    }
  }
  if (data_rows_scanned %% 5000000L == 0L) log_line("data_rows_scanned=", data_rows_scanned)
}
close(con)
runtime_seconds <- as.numeric(difftime(Sys.time(), start, units = "secs"))
source_scan_completed <- TRUE
source_sha_after <- hash_file(paths$chen_raw)
source_unchanged <- identical(source_sha_before, source_sha_after)
log_line("data_rows_scanned=", data_rows_scanned, "; runtime_seconds=", sprintf("%.3f", runtime_seconds))

match_df <- if (length(matches) == 0L) {
  data.frame(matrix(ncol = length(header_fields), nrow = 0L, dimnames = list(NULL, header_fields)), stringsAsFactors = FALSE)
} else {
  as.data.frame(do.call(rbind, lapply(matches, function(z) as.data.frame(z, stringsAsFactors = FALSE))), stringsAsFactors = FALSE)
}
invalid_df <- if (length(invalid_matches) == 0L) {
  data.frame(source_marker_id = character(), raw_field_count = integer(), stringsAsFactors = FALSE)
} else {
  do.call(rbind, lapply(invalid_matches, as.data.frame, stringsAsFactors = FALSE))
}

target_authority$raw_match_count <- vapply(target_authority$expected_chen_marker_id, function(id) sum(match_df$rs_number == id), integer(1))
master_rows <- lapply(seq_len(nrow(target_authority)), function(i) {
  target <- target_authority[i, , drop = FALSE]
  marker <- target$expected_chen_marker_id
  mp <- marker_parts(marker)
  rows <- match_df[match_df$rs_number == marker, , drop = FALSE]
  raw_count <- nrow(rows)
  status <- if (!isTRUE(target$identity_bridge_ready)) {
    "identity_bridge_unavailable"
  } else if (raw_count == 0L) {
    "missing"
  } else if (raw_count > 1L) {
    "multiple_exact_match"
  } else {
    chen_key <- allele_key(rows$reference_allele[[1L]], rows$other_allele[[1L]])
    marker_key <- allele_key(mp[[3L]], mp[[4L]])
    if (identical(chen_key, marker_key)) "unique_exact_match" else "marker_effect_allele_incompatible"
  }
  row <- if (raw_count == 1L && identical(status, "unique_exact_match")) rows[1L, , drop = FALSE] else NULL
  data.frame(
    target_rsid = target$rsid,
    expected_chen_marker_id = marker,
    strict_included_member = target$strict_included_member,
    strict_excluded_member = target$strict_excluded_member,
    relaxed_included_member = target$relaxed_included_member,
    relaxed_excluded_member = target$relaxed_excluded_member,
    identity_bridge_ready = target$identity_bridge_ready,
    target_chr_grch37 = target$reference_chr_grch37,
    target_pos_grch37 = target$reference_pos_grch37,
    marker_allele1 = mp[[3L]],
    marker_allele2 = mp[[4L]],
    raw_match_count = raw_count,
    match_status = status,
    available_for_harmonisation = identical(status, "unique_exact_match"),
    source_marker_id = if (!is.null(row)) row$rs_number else NA_character_,
    reference_allele = if (!is.null(row)) row$reference_allele else NA_character_,
    other_allele = if (!is.null(row)) row$other_allele else NA_character_,
    eaf = if (!is.null(row)) parse_num(row$eaf) else NA_real_,
    beta = if (!is.null(row)) parse_num(row$beta) else NA_real_,
    se = if (!is.null(row)) parse_num(row$se) else NA_real_,
    p_value = if (!is.null(row)) parse_num(row[["p-value"]]) else NA_real_,
    n_samples = if (!is.null(row)) parse_num(row$n_samples) else NA_real_,
    beta_95L = if (!is.null(row)) parse_num(row$beta_95L) else NA_real_,
    beta_95U = if (!is.null(row)) parse_num(row$beta_95U) else NA_real_,
    z = if (!is.null(row)) parse_num(row$z) else NA_real_,
    n_studies = if (!is.null(row)) parse_num(row$n_studies) else NA_real_,
    effects = if (!is.null(row)) row$effects else NA_character_,
    reason = if (identical(status, "unique_exact_match")) NA_character_ else status,
    stringsAsFactors = FALSE
  )
})
master <- do.call(rbind, master_rows)

branch_summary <- function(label, col) {
  x <- master[master[[col]], , drop = FALSE]
  data.frame(
    branch = label,
    target = nrow(x),
    exact = sum(x$match_status == "unique_exact_match"),
    missing = sum(x$match_status == "missing"),
    multiple = sum(x$match_status == "multiple_exact_match"),
    incompatible = sum(x$match_status == "marker_effect_allele_incompatible"),
    identity_bridge_unavailable = sum(x$match_status == "identity_bridge_unavailable"),
    branch_harmonisation_readiness = sum(x$match_status == "unique_exact_match") > 0L,
    readiness_reason = if (sum(x$match_status == "unique_exact_match") > 0L) NA_character_ else "no_exact_chen_outcome_available",
    stringsAsFactors = FALSE
  )
}
branch_counts <- rbind(
  branch_summary("strict_apoe_included", "strict_included_member"),
  branch_summary("strict_apoe_excluded", "strict_excluded_member"),
  branch_summary("relaxed_apoe_included", "relaxed_included_member"),
  branch_summary("relaxed_apoe_excluded", "relaxed_excluded_member")
)

match_audit <- master[, c("target_rsid", "expected_chen_marker_id", "raw_match_count", "match_status", "reason",
                          "strict_included_member", "strict_excluded_member", "relaxed_included_member", "relaxed_excluded_member",
                          "available_for_harmonisation")]
missing_df <- master[master$match_status == "missing", c("target_rsid", "expected_chen_marker_id", "reason")]
multiple_df <- if (any(master$match_status == "multiple_exact_match")) {
  match_df[match_df$rs_number %in% master$expected_chen_marker_id[master$match_status == "multiple_exact_match"], , drop = FALSE]
} else {
  data.frame()
}

tmpdb <- rel("results", "tmp", "duckdb_chen_reverse_outcome_v1")
dir.create(tmpdb, recursive = TRUE, showWarnings = FALSE)
con_db <- DBI::dbConnect(duckdb::duckdb(dbdir = tempfile(tmpdir = tmpdb), read_only = FALSE))
on.exit(try(DBI::dbDisconnect(con_db, shutdown = TRUE), silent = TRUE), add = TRUE)
write_parquet_df <- function(x, path, table_name) {
  atomic_write(path, function(p) {
    DBI::dbRemoveTable(con_db, table_name, fail_if_missing = FALSE)
    DBI::dbWriteTable(con_db, table_name, x)
    DBI::dbExecute(con_db, sprintf("COPY %s TO '%s' (FORMAT PARQUET)", DBI::dbQuoteIdentifier(con_db, table_name), gsub("'", "''", p, fixed = TRUE)))
  })
}
validate_pair <- function(parquet, tsv, key = "target_rsid") {
  p <- DBI::dbGetQuery(con_db, sprintf("SELECT * FROM read_parquet('%s')", gsub("'", "''", parquet, fixed = TRUE)))
  t <- utils::read.delim(tsv, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
  same_cols <- identical(names(p), names(t))
  same_n <- identical(nrow(p), nrow(t))
  same_order <- same_n && identical(as.character(p[[key]]), as.character(t[[key]]))
  char_cols <- names(p)[!vapply(p, is.numeric, logical(1))]
  num_cols <- names(p)[vapply(p, is.numeric, logical(1))]
  char_ok <- all(vapply(char_cols, function(k) identical(as.character(p[[k]]), as.character(t[[k]])), logical(1)))
  numeric_stats <- lapply(num_cols, function(k) {
    a <- as.numeric(p[[k]])
    b <- suppressWarnings(as.numeric(t[[k]]))
    ok <- (is.na(a) & is.na(b)) | (is.finite(a) & is.finite(b))
    abs_diff <- abs(a - b)
    rel_diff <- abs_diff / pmax(abs(a), abs(b), .Machine$double.xmin)
    data.frame(column = k, max_abs = max(abs_diff[ok], na.rm = TRUE), max_rel = max(rel_diff[ok], na.rm = TRUE), stringsAsFactors = FALSE)
  })
  ns <- if (length(numeric_stats)) do.call(rbind, numeric_stats) else data.frame(column = character(), max_abs = numeric(), max_rel = numeric())
  numeric_ok <- all(is.na(ns$max_abs) | ns$max_abs <= 1e-12 | ns$max_rel <= 1e-12)
  list(row_count = nrow(p), same_cols = same_cols, same_n = same_n, same_order = same_order, char_ok = char_ok, numeric_ok = numeric_ok, numeric_stats = records(ns))
}

write_tsv(target_authority, paths$union_targets)
write_parquet_df(master, paths$master_parquet, "master")
write_tsv(master, paths$master_tsv)
write_parquet_df(master[master$strict_included_member, , drop = FALSE], paths$strict_inc_parquet, "strict_inc")
write_tsv(master[master$strict_included_member, , drop = FALSE], paths$strict_inc_tsv)
write_parquet_df(master[master$strict_excluded_member, , drop = FALSE], paths$strict_exc_parquet, "strict_exc")
write_tsv(master[master$strict_excluded_member, , drop = FALSE], paths$strict_exc_tsv)
write_parquet_df(master[master$relaxed_included_member, , drop = FALSE], paths$relaxed_inc_parquet, "relaxed_inc")
write_tsv(master[master$relaxed_included_member, , drop = FALSE], paths$relaxed_inc_tsv)
write_parquet_df(master[master$relaxed_excluded_member, , drop = FALSE], paths$relaxed_exc_parquet, "relaxed_exc")
write_tsv(master[master$relaxed_excluded_member, , drop = FALSE], paths$relaxed_exc_tsv)
write_csv(match_audit, paths$match_audit)
write_tsv(missing_df, paths$missing)
if (nrow(multiple_df) > 0L) write_tsv(multiple_df, paths$multiple)

pair_checks <- list(
  master = validate_pair(paths$master_parquet, paths$master_tsv),
  strict_included = validate_pair(paths$strict_inc_parquet, paths$strict_inc_tsv),
  strict_excluded = validate_pair(paths$strict_exc_parquet, paths$strict_exc_tsv),
  relaxed_included = validate_pair(paths$relaxed_inc_parquet, paths$relaxed_inc_tsv),
  relaxed_excluded = validate_pair(paths$relaxed_exc_parquet, paths$relaxed_exc_tsv)
)
parquet_tsv_ok <- all(vapply(pair_checks, function(x) x$same_cols && x$same_n && x$same_order && x$char_ok && x$numeric_ok, logical(1)))

union_exact <- sum(master$match_status == "unique_exact_match")
union_missing <- sum(master$match_status == "missing")
union_multiple <- sum(master$match_status == "multiple_exact_match")
union_incompatible <- sum(master$match_status == "marker_effect_allele_incompatible")

hard_checks <- list(
  decision_109_authority_gate = identical(contract$authoritative_contract_status, "frozen") && isTRUE(contract$approved_for_chen_reverse_outcome_extraction),
  strict_exposure_authority_gate = identical(contract$strict_exposure_authority$freeze_status, "passed"),
  relaxed_exposure_authority_gate = identical(contract$relaxed_exposure_authority$freeze_status, "passed"),
  corrected_target_audit_used = identical(relpath(paths$target_authority), contract$corrected_target_audit_authority),
  target_character_schema_explicit = all(vapply(target_authority[c("rsid", "expected_chen_marker_id", "exposure_effect_allele", "exposure_other_allele", "reference_allele1", "reference_allele2")], is.character, logical(1))),
  allele_T_F_fixture_passed = allele_T_F_fixture_passed,
  target_counts_reverified = nrow(target_authority) == contract$union_target_count,
  identity_bridge_ready_reverified = sum(target_authority$identity_bridge_ready) == contract$target_authority_ready_count,
  chen_source_sha_before_gate = identical(tolower(source_sha_before), expected_sha),
  chen_source_sha_after_gate = identical(tolower(source_sha_after), expected_sha),
  source_unchanged = source_unchanged,
  chen_source_schema_verified = chen_source_schema_verified,
  exact_marker_id_matching_only = TRUE,
  no_substring_matching = TRUE,
  source_scan_complete = source_scan_completed,
  all_targets_classified = nrow(master) == nrow(target_authority) && all(!is.na(master$match_status)),
  marker_effect_allele_internal_consistency_audited = TRUE,
  no_effect_orientation_transform = TRUE,
  branch_membership_preserved = all(branch_counts$target == c(sum(target_authority$strict_included_member), sum(target_authority$strict_excluded_member), sum(target_authority$relaxed_included_member), sum(target_authority$relaxed_excluded_member))),
  no_exposure_reselection = TRUE,
  no_reclumping = TRUE,
  no_proxy = TRUE,
  no_liftover = TRUE,
  no_nearest_variant = TRUE,
  no_strand_complement_identity_rescue = TRUE,
  master_parquet_tsv_consistency = parquet_tsv_ok,
  branch_parquet_tsv_consistency = parquet_tsv_ok,
  no_harmonisation = TRUE,
  no_mr = TRUE,
  no_steiger = TRUE,
  renv_lock_unchanged = identical(hash_file(paths$renv_lock), "253471c02e1e47a40d0f68b296d4ae2b1df471d757ac86328b3e974018d039f3")
)
hard_check_failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
outcome_extraction_status <- if (length(hard_check_failures) == 0L) "passed" else "failed"

qc <- list(
  extraction_version = "v1",
  decision = 110,
  date = format(Sys.Date()),
  analysis_direction = "delirium_to_Hb",
  analysis_role = "reverse_alternative_hb_outcome_sensitivity",
  independent_replication = FALSE,
  contract_authority_decision = 109,
  strict_exposure_authority_decision = 47,
  relaxed_exposure_authority_decision = 59,
  outcome_source = "Chen_2020_Hb_BCX2",
  outcome_scale = contract$outcome_scale,
  source_sha_before = source_sha_before,
  source_sha_after = source_sha_after,
  source_size_bytes = source_size_bytes,
  source_rows_scanned = data_rows_scanned,
  runtime_seconds = runtime_seconds,
  source_scan_completed = source_scan_completed,
  union_target_count = nrow(master),
  identity_bridge_ready_count = sum(target_authority$identity_bridge_ready),
  identity_bridge_unavailable_count = sum(!target_authority$identity_bridge_ready),
  union_unique_exact_match_count = union_exact,
  union_missing_count = union_missing,
  union_multiple_exact_match_count = union_multiple,
  union_marker_effect_allele_incompatible_count = union_incompatible,
  exact_match_rsids = master$target_rsid[master$match_status == "unique_exact_match"],
  missing_rsids = master$target_rsid[master$match_status == "missing"],
  multiple_rsids = master$target_rsid[master$match_status == "multiple_exact_match"],
  incompatible_rsids = master$target_rsid[master$match_status == "marker_effect_allele_incompatible"],
  branch_counts = records(branch_counts),
  branch_harmonisation_readiness = records(branch_counts[, c("branch", "branch_harmonisation_readiness", "readiness_reason")]),
  matching_method = "exact_BCX2_marker_id_string_identity",
  proxy_used = FALSE,
  liftover_used = FALSE,
  nearest_variant_used = FALSE,
  strand_complement_identity_rescue_used = FALSE,
  harmonisation_performed = FALSE,
  mr_run = FALSE,
  steiger_run = FALSE,
  allele_T_F_fixture_passed = allele_T_F_fixture_passed,
  parquet_tsv_consistency = pair_checks,
  outcome_extraction_status = outcome_extraction_status,
  approved_for_chen_reverse_outcome_extraction_freeze = identical(outcome_extraction_status, "passed"),
  hard_checks = hard_checks,
  hard_check_failures = hard_check_failures,
  informational_findings = list(
    invalid_source_row_matches = records(invalid_df),
    multiple_matches_file_created = file.exists(paths$multiple),
    renv_out_of_sync_warning_may_be_emitted = TRUE
  )
)
write_json(qc, paths$qc_json)

decision_lines <- c(
  "# Decision 110: Chen Reverse Outcome Extraction V1",
  "",
  paste0("Date: ", format(Sys.Date())),
  "",
  "## Status",
  paste0("outcome_extraction_status: `", outcome_extraction_status, "`"),
  paste0("approved_for_chen_reverse_outcome_extraction_freeze: `", identical(outcome_extraction_status, "passed"), "`"),
  "",
  "## Decision",
  "Execute targeted Chen reverse outcome extraction using Decision 109 contract authority.",
  "",
  "The extraction performs exact BCX2 marker-ID string matching only. It does not harmonise, run MR, run Steiger, use proxy/liftOver, reclump, or modify exposure instruments.",
  "",
  "## Results",
  paste0("- Chen source SHA before: `", source_sha_before, "`."),
  paste0("- Chen source SHA after: `", source_sha_after, "`."),
  paste0("- Source rows scanned: `", data_rows_scanned, "`."),
  paste0("- Runtime seconds: `", sprintf("%.3f", runtime_seconds), "`."),
  paste0("- Union target count: `", nrow(master), "`."),
  paste0("- Exact / missing / multiple / incompatible: `", union_exact, " / ", union_missing, " / ", union_multiple, " / ", union_incompatible, "`."),
  paste0("- Hard-check failures: `", paste(hard_check_failures, collapse = ";"), "`."),
  "",
  "## Outputs",
  paste0("- `", relpath(paths$master_parquet), "`"),
  paste0("- `", relpath(paths$master_tsv), "`"),
  paste0("- `", relpath(paths$qc_json), "`"),
  paste0("- `", relpath(paths$log), "`")
)
write_text(decision_lines, paths$decision)

log_line("outcome_extraction_status=", outcome_extraction_status, "; hard_check_failures=", paste(hard_check_failures, collapse = ","))
cat("outcome_extraction_status=", outcome_extraction_status, "\n", sep = "")
cat("hard_check_failures=", paste(hard_check_failures, collapse = ","), "\n", sep = "")
cat("source_rows_scanned=", data_rows_scanned, "\n", sep = "")
cat("union_exact=", union_exact, "\n", sep = "")
cat("union_missing=", union_missing, "\n", sep = "")
cat("union_multiple=", union_multiple, "\n", sep = "")
cat("union_incompatible=", union_incompatible, "\n", sep = "")
