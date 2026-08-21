#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

pkgs <- c("digest", "DBI", "duckdb", "jsonlite")
if (any(!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))) {
  stop("Required installed package missing; no automatic installation.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript.exe R/12_reverse_primary_outcome_extraction_v1_freeze_manifest.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}

root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

stop_if <- function(x, m) if (isTRUE(x)) stop(m, call. = FALSE)
hash_file <- function(p) digest::digest(file = p, algo = "sha256", serialize = FALSE)
sql_string <- function(p, must_work = TRUE) paste0("'", gsub("'", "''", normalizePath(p, winslash = "/", mustWork = must_work), fixed = TRUE), "'")
records <- function(x) lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
num_equal <- function(a, b, tol = 1e-12) {
  a <- as.numeric(a); b <- as.numeric(b)
  (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tol)
}

paths <- list(
  script = file.path(root, "R", "12_extract_vuckovic_hb_reverse_primary_outcomes_v1.R"),
  targets = file.path(root, "data_derived", "reverse_outcome_extraction", "finngen_r13_delirium_reverse_primary_targets_v1.tsv"),
  all_pq = file.path(root, "data_derived", "reverse_outcome_extraction", "vuckovic_hb_reverse_primary_all_matches_v1.parquet"),
  all_tsv = file.path(root, "data_derived", "reverse_outcome_extraction", "vuckovic_hb_reverse_primary_all_matches_v1.tsv"),
  unique_pq = file.path(root, "data_derived", "reverse_outcome_extraction", "vuckovic_hb_reverse_primary_unique_matches_v1.parquet"),
  unique_tsv = file.path(root, "data_derived", "reverse_outcome_extraction", "vuckovic_hb_reverse_primary_unique_matches_v1.tsv"),
  status_csv = file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_target_match_status_v1.csv"),
  extraction_json = file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_outcome_extraction_v1.json"),
  extraction_log = file.path(root, "results", "logs", "vuckovic_hb_reverse_primary_outcome_extraction_v1.log"),
  decision48 = file.path(root, "docs", "decisions", "48_vuckovic_hb_reverse_primary_outcome_extraction_v1_v1.1.md"),
  upstream_freeze_json = file.path(root, "results", "qc", "finngen_r13_delirium_reverse_primary_instruments_v4_freeze_v3.json"),
  upstream_freeze_manifest = file.path(root, "results", "qc", "finngen_r13_delirium_reverse_primary_instruments_v4_freeze_manifest_v3.csv"),
  source = file.path(root, "data_raw", "gwas", "ebi-a-GCST90002384.vcf.gz")
)

out <- c(
  manifest = file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_outcome_extraction_v1_freeze_manifest.csv"),
  json = file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_outcome_extraction_v1_freeze.json"),
  log = file.path(root, "results", "logs", "vuckovic_hb_reverse_primary_outcome_extraction_v1_freeze.log")
)
partials <- paste0(out, ".partial")

stop_if(any(!file.exists(unlist(paths))), "Required extraction or upstream authority file is missing.")
stop_if(any(file.exists(c(out, partials))), "A freeze final or partial target exists; refusing to overwrite.")
dir.create(dirname(out[["manifest"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["log"]]), recursive = TRUE, showWarnings = FALSE)

log_line <- function(x) cat(sprintf("%s %s\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), x), file = out[["log"]], append = TRUE)

main <- function() {
  log_line("stage=reverse_primary_outcome_extraction_v1_freeze")
  qc <- jsonlite::fromJSON(paths$extraction_json, simplifyVector = FALSE)
  upstream <- jsonlite::fromJSON(paths$upstream_freeze_json, simplifyVector = FALSE)
  hard_check_failures <- character()
  add_fail <- function(x) hard_check_failures <<- unique(c(hard_check_failures, x))

  if (!identical(qc$extraction_status, "passed")) add_fail("extraction_status_not_passed")
  if (!isTRUE(qc$approved_for_reverse_primary_harmonisation_preflight)) add_fail("not_approved_for_preflight")
  if (length(qc$hard_check_failures) != 0L) add_fail("extraction_hard_check_failures_not_empty")
  if (!identical(upstream$freeze_status, "passed")) add_fail("upstream_reverse_instrument_freeze_not_passed")

  expected_source_sha <- "c68c98c7800c59d9d64cb88739e2e245f8aad6e4e886455cbf7424661afc3d41"
  source_sha_before <- hash_file(paths$source)
  if (!identical(source_sha_before, expected_source_sha)) add_fail("vuckovic_source_sha_mismatch")
  if (!identical(tolower(qc$vuckovic_source$sha256), expected_source_sha)) add_fail("qc_vuckovic_source_sha_mismatch")

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  all_pq <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(paths$all_pq)))
  unique_pq <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(paths$unique_pq)))
  all_tsv <- read.delim(paths$all_tsv, check.names = FALSE)
  unique_tsv <- read.delim(paths$unique_tsv, check.names = FALSE)
  targets <- read.delim(paths$targets, check.names = FALSE)
  status <- read.csv(paths$status_csv, check.names = FALSE)

  compare_cols <- c("target_rsid", "membership", "match_class", "exposure_effect_allele_raw", "exposure_other_allele_raw",
                    "outcome_effect_allele_raw", "outcome_other_allele_raw", "outcome_beta_raw", "outcome_se_raw",
                    "outcome_lp_raw", "outcome_pval_raw", "outcome_eaf_raw")
  if (!identical(names(all_pq), names(all_tsv))) add_fail("all_matches_parquet_tsv_schema_mismatch")
  if (!identical(names(unique_pq), names(unique_tsv))) add_fail("unique_matches_parquet_tsv_schema_mismatch")
  if (nrow(all_pq) != nrow(all_tsv) || nrow(unique_pq) != nrow(unique_tsv)) add_fail("parquet_tsv_row_count_mismatch")
  for (col in intersect(compare_cols, names(all_pq))) {
    if (is.numeric(all_pq[[col]]) || is.numeric(all_tsv[[col]])) {
      if (!all(num_equal(all_pq[[col]], all_tsv[[col]]))) add_fail(paste0("all_matches_numeric_mismatch_", col))
    } else if (!identical(as.character(all_pq[[col]]), as.character(all_tsv[[col]]))) {
      add_fail(paste0("all_matches_character_mismatch_", col))
    }
  }
  for (col in intersect(compare_cols, names(unique_pq))) {
    if (is.numeric(unique_pq[[col]]) || is.numeric(unique_tsv[[col]])) {
      if (!all(num_equal(unique_pq[[col]], unique_tsv[[col]]))) add_fail(paste0("unique_matches_numeric_mismatch_", col))
    } else if (!identical(as.character(unique_pq[[col]]), as.character(unique_tsv[[col]]))) {
      add_fail(paste0("unique_matches_character_mismatch_", col))
    }
  }

  if (!setequal(targets$target_rsid, status$target_rsid) || !setequal(targets$target_rsid, unique_pq$target_rsid)) add_fail("target_rsid_set_mismatch")
  if (!identical(sort(targets$target_rsid), sort(c("rs429358", "rs58537897")))) add_fail("target_rsids_not_expected_from_formal_files")
  if (!identical(qc$targets$included, "rs429358") || !identical(qc$targets$excluded, "rs58537897")) add_fail("qc_target_identity_mismatch")
  if (any(status$match_class != "unique")) add_fail("non_unique_status_present")
  if (nrow(unique_pq) != 2L || nrow(all_pq) != 2L) add_fail("match_row_count_not_two")

  counts <- list(
    included_target_count = sum(targets$membership == "apoe_included"),
    excluded_target_count = sum(targets$membership == "apoe_excluded"),
    union_target_count = length(unique(targets$target_rsid)),
    unique_match_count = sum(status$match_class == "unique"),
    missing_match_count = sum(status$match_class == "missing"),
    multiple_match_count = sum(status$match_class == "multiple"),
    included_unique_match_count = sum(status$membership == "apoe_included" & status$match_class == "unique"),
    excluded_unique_match_count = sum(status$membership == "apoe_excluded" & status$match_class == "unique")
  )

  hard_checks <- list(
    extraction_status_passed = identical(qc$extraction_status, "passed"),
    approved_for_reverse_primary_harmonisation_preflight = isTRUE(qc$approved_for_reverse_primary_harmonisation_preflight),
    extraction_hard_check_failures_empty = length(qc$hard_check_failures) == 0L,
    target_identity_from_formal_files = identical(sort(targets$target_rsid), sort(c("rs429358", "rs58537897"))),
    all_source_matches_unique = all(status$match_class == "unique"),
    parquet_tsv_consistency = length(grep("parquet_tsv|numeric_mismatch|character_mismatch", hard_check_failures)) == 0L,
    no_source_mutation = identical(source_sha_before, expected_source_sha),
    no_vuckovic_rescan = TRUE,
    no_harmonisation = TRUE,
    no_reverse_mr = TRUE,
    relaxed_threshold_not_started = TRUE
  )
  freeze_status <- if (length(hard_check_failures) == 0L && all(unlist(hard_checks))) "passed" else "failed"
  approved <- identical(freeze_status, "passed")

  rel <- c(
    "R/12_extract_vuckovic_hb_reverse_primary_outcomes_v1.R",
    "data_derived/reverse_outcome_extraction/finngen_r13_delirium_reverse_primary_targets_v1.tsv",
    "data_derived/reverse_outcome_extraction/vuckovic_hb_reverse_primary_all_matches_v1.parquet",
    "data_derived/reverse_outcome_extraction/vuckovic_hb_reverse_primary_all_matches_v1.tsv",
    "data_derived/reverse_outcome_extraction/vuckovic_hb_reverse_primary_unique_matches_v1.parquet",
    "data_derived/reverse_outcome_extraction/vuckovic_hb_reverse_primary_unique_matches_v1.tsv",
    "results/qc/vuckovic_hb_reverse_primary_target_match_status_v1.csv",
    "results/qc/vuckovic_hb_reverse_primary_outcome_extraction_v1.json",
    "results/logs/vuckovic_hb_reverse_primary_outcome_extraction_v1.log",
    "docs/decisions/48_vuckovic_hb_reverse_primary_outcome_extraction_v1_v1.1.md",
    "results/qc/finngen_r13_delirium_reverse_primary_instruments_v4_freeze_v3.json",
    "results/qc/finngen_r13_delirium_reverse_primary_instruments_v4_freeze_manifest_v3.csv"
  )
  manifest <- data.frame(
    file_role = "authority",
    relative_path = rel,
    file_size_bytes = file.info(file.path(root, rel))$size,
    sha256 = vapply(file.path(root, rel), hash_file, character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(manifest, paste0(out[["manifest"]], ".partial"), row.names = FALSE)
  manifest_sha256 <- hash_file(paste0(out[["manifest"]], ".partial"))

  js <- c(list(
    freeze_version = "v1",
    authoritative_reverse_primary_outcome_extraction_version = "v1",
    analysis_direction = "delirium_to_Hb",
    included_target_rsid = unique(status$target_rsid[status$membership == "apoe_included"]),
    excluded_target_rsid = unique(status$target_rsid[status$membership == "apoe_excluded"]),
    vuckovic_source_path = qc$vuckovic_source$path,
    vuckovic_source_sha256 = expected_source_sha,
    manifest_sha256 = manifest_sha256,
    freeze_status = freeze_status,
    approved_for_reverse_primary_harmonisation_preflight = approved,
    hard_checks = hard_checks,
    hard_check_failures = hard_check_failures
  ), counts)
  jsonlite::write_json(js, paste0(out[["json"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")

  for (p in out[c("manifest", "json")]) {
    stop_if(file.exists(p), paste("Output appeared during run:", p))
    stop_if(!file.rename(paste0(p, ".partial"), p), paste("Atomic rename failed:", p))
  }
  log_line(sprintf("SUCCESS freeze_status=%s manifest_sha256=%s hard_check_failures=%s", freeze_status, manifest_sha256, paste(hard_check_failures, collapse = ";")))
  if (!approved) stop("Freeze did not pass; harmonisation preflight must not proceed.", call. = FALSE)
}

tryCatch(main(), error = function(e) {
  log_line(paste0("TERMINATED_PRIMARY: ", conditionMessage(e)))
  quit(status = 1L)
})

