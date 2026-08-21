#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/17_reverse_relaxed_harmonisation_preflight_v1_freeze_manifest.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

for (pkg in c("DBI", "duckdb", "jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
hash_file <- function(path) digest::digest(file = path, algo = "sha256")
sql_string <- function(path, must_work = TRUE) {
  paste0("'", gsub("'", "''", normalizePath(path, winslash = "/", mustWork = must_work), fixed = TRUE), "'")
}
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
num_equal <- function(a, b, tol = 1e-12) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tol)
}
canon <- function(x) {
  if (is.logical(x)) return(toupper(as.character(x)))
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}
table_equal <- function(a, b, tol = 1e-12) {
  if (!identical(names(a), names(b)) || nrow(a) != nrow(b) || anyDuplicated(names(a))) return(FALSE)
  for (col in names(a)) {
    if (is.numeric(a[[col]]) || is.numeric(b[[col]])) {
      if (!all(num_equal(a[[col]], b[[col]], tol))) return(FALSE)
    } else if (!identical(canon(a[[col]]), canon(b[[col]]))) {
      return(FALSE)
    }
  }
  TRUE
}

out <- c(
  manifest = file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1_freeze_manifest.csv"),
  json = file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1_freeze.json"),
  log = file.path(root, "results", "logs", "vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1_freeze.log")
)
stop_if(any(file.exists(c(out, paste0(out, ".partial")))), "A reverse relaxed preflight freeze final or partial target exists; refusing to overwrite.")
dir.create(dirname(out[["manifest"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["log"]]), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = out[["log"]], append = TRUE)

main <- function() {
  log_line("stage=reverse_relaxed_harmonisation_preflight_v1_freeze")
  renv_lock <- file.path(root, "renv.lock")
  renv_before <- hash_file(renv_lock)
  rel <- c(
    "docs/decisions/59_finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze_v1.1.md",
    "docs/decisions/61_vuckovic_hb_reverse_relaxed_outcome_extraction_v1_freeze_v1.1.md",
    "docs/decisions/62_vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1_v1.1.md",
    "docs/decisions/63_vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1_freeze_v1.1.md",
    "R/17_reverse_relaxed_harmonisation_preflight_v1.R",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1.parquet",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1.tsv",
    "results/qc/vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1.json",
    "results/logs/vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1.log",
    "results/qc/vuckovic_hb_reverse_relaxed_outcome_extraction_v1_freeze.json",
    "results/qc/finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze.json",
    "renv.lock"
  )
  paths <- file.path(root, rel)
  stop_if(any(!file.exists(paths)), paste("Missing freeze input(s):", paste(rel[!file.exists(paths)], collapse = "; ")))

  qc <- jsonlite::fromJSON(file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1.json"), simplifyVector = FALSE)
  outcome_freeze <- jsonlite::fromJSON(file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_outcome_extraction_v1_freeze.json"), simplifyVector = FALSE)
  instrument_freeze <- jsonlite::fromJSON(file.path(root, "results", "qc", "finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze.json"), simplifyVector = FALSE)

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  pq <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(file.path(root, rel[6]))))
  tsv <- read.delim(file.path(root, rel[7]), check.names = FALSE)

  failures <- character()
  add_fail <- function(x) failures <<- unique(c(failures, x))
  if (!identical(qc$preflight_status, "passed")) add_fail("preflight_status_not_passed")
  if (length(qc$hard_check_failures) != 0L) add_fail("preflight_hard_check_failures_not_empty")
  if (!identical(qc$harmonisation_performed, FALSE)) add_fail("harmonisation_performed_not_false")
  if (!isTRUE(qc$palindromic_rule_adjudication_required)) add_fail("palindromic_rule_adjudication_required_not_true")
  if (!identical(qc$approved_for_reverse_relaxed_formal_harmonisation, FALSE)) add_fail("formal_harmonisation_already_approved")
  if (as.integer(qc$beta_flip_performed_count) != 0L) add_fail("beta_flip_performed_count_not_zero")
  if (as.integer(qc$eaf_flip_performed_count) != 0L) add_fail("eaf_flip_performed_count_not_zero")
  if (as.integer(qc$strand_flip_performed_count) != 0L) add_fail("strand_flip_performed_count_not_zero")
  if (!identical(outcome_freeze$freeze_status, "passed")) add_fail("outcome_extraction_freeze_not_passed")
  if (!identical(instrument_freeze$freeze_status, "passed")) add_fail("instrument_freeze_not_passed")
  if (!table_equal(pq, tsv)) add_fail("parquet_tsv_consistency_failed")

  required_cols <- c(
    "target_rsid", "included_member", "excluded_member",
    "exposure_effect_allele", "exposure_other_allele", "exposure_beta", "exposure_se", "exposure_pval", "exposure_eaf",
    "outcome_effect_allele_raw", "outcome_other_allele_raw", "outcome_beta_raw", "outcome_se_raw", "outcome_pval_raw", "outcome_eaf_raw",
    "palindromic_snp", "raw_orientation_class", "requires_palindromic_rule_adjudication",
    "eligible_for_formal_harmonisation_without_rule", "potentially_harmonisable_after_rule_adjudication"
  )
  missing_cols <- setdiff(required_cols, names(pq))
  if (length(missing_cols) > 0L) add_fail(paste0("missing_required_columns:", paste(missing_cols, collapse = ",")))
  if (anyDuplicated(pq$target_rsid)) add_fail("duplicate_target_rsids")
  if (nrow(pq) != as.integer(qc$union_preflight_count)) add_fail("union_preflight_count_mismatch")
  if (sum(pq$included_member) != as.integer(qc$included_preflight_count)) add_fail("included_preflight_count_mismatch")
  if (sum(pq$excluded_member) != as.integer(qc$excluded_preflight_count)) add_fail("excluded_preflight_count_mismatch")
  if (sum(pq$palindromic_snp) != as.integer(qc$palindromic_count)) add_fail("palindromic_count_mismatch")
  if (sum(pq$requires_palindromic_rule_adjudication) != as.integer(qc$requires_palindromic_rule_adjudication_count)) add_fail("requires_rule_count_mismatch")
  if (sum(pq$raw_orientation_class == "exact_match") != as.integer(qc$exact_match_count)) add_fail("exact_match_count_mismatch")
  if (sum(pq$raw_orientation_class == "swapped_match") != as.integer(qc$swapped_match_count)) add_fail("swapped_match_count_mismatch")
  if (sum(pq$raw_orientation_class == "strand_exact_match") != as.integer(qc$strand_exact_match_count)) add_fail("strand_exact_match_count_mismatch")
  if (sum(pq$raw_orientation_class == "strand_swapped_match") != as.integer(qc$strand_swapped_match_count)) add_fail("strand_swapped_match_count_mismatch")
  if (sum(pq$raw_orientation_class == "incompatible") != as.integer(qc$incompatible_count)) add_fail("incompatible_count_mismatch")
  if (sum(pq$raw_orientation_class == "invalid") != as.integer(qc$invalid_count)) add_fail("invalid_count_mismatch")

  manifest <- data.frame(
    file_role = c("instrument_freeze_decision", "outcome_extraction_freeze_decision", "preflight_decision",
                  "preflight_freeze_decision", "preflight_script", "preflight_parquet", "preflight_tsv",
                  "preflight_qc_json", "preflight_log", "outcome_extraction_freeze_json",
                  "instrument_freeze_json", "renv_lock"),
    relative_path = rel,
    file_size_bytes = as.numeric(file.info(paths)$size),
    sha256 = vapply(paths, hash_file, character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(manifest, paste0(out[["manifest"]], ".partial"), row.names = FALSE)
  manifest_sha <- hash_file(paste0(out[["manifest"]], ".partial"))
  renv_after <- hash_file(renv_lock)
  if (!identical(renv_before, renv_after)) add_fail("renv_lock_changed")

  freeze_status <- if (length(failures) == 0L) "passed" else "failed"
  result <- list(
    freeze_version = "v1",
    decision = 63,
    authoritative_reverse_relaxed_harmonisation_preflight_version = "v1",
    source_relaxed_instrument_version = "v2",
    source_outcome_extraction_version = "v1",
    source_outcome_extraction_freeze_version = "v1",
    analysis_direction = "delirium_to_Hb",
    analysis_role = "secondary_reverse_exploratory_relaxed",
    branch_type = "protocol_prespecified_fallback",
    p_threshold = 5e-6,
    union_preflight_count = nrow(pq),
    included_preflight_count = sum(pq$included_member),
    excluded_preflight_count = sum(pq$excluded_member),
    palindromic_count = sum(pq$palindromic_snp),
    requires_palindromic_rule_adjudication_count = sum(pq$requires_palindromic_rule_adjudication),
    exact_match_count = sum(pq$raw_orientation_class == "exact_match"),
    swapped_match_count = sum(pq$raw_orientation_class == "swapped_match"),
    strand_exact_match_count = sum(pq$raw_orientation_class == "strand_exact_match"),
    strand_swapped_match_count = sum(pq$raw_orientation_class == "strand_swapped_match"),
    incompatible_count = sum(pq$raw_orientation_class == "incompatible"),
    invalid_count = sum(pq$raw_orientation_class == "invalid"),
    beta_flip_performed_count = as.integer(qc$beta_flip_performed_count),
    eaf_flip_performed_count = as.integer(qc$eaf_flip_performed_count),
    strand_flip_performed_count = as.integer(qc$strand_flip_performed_count),
    harmonisation_performed = FALSE,
    palindromic_rule_adjudication_required = TRUE,
    approved_for_reverse_relaxed_formal_harmonisation = FALSE,
    manifest_sha256 = manifest_sha,
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after,
    freeze_status = freeze_status,
    hard_checks = list(
      preflight_status_passed = identical(qc$preflight_status, "passed"),
      preflight_hard_check_failures_empty = length(qc$hard_check_failures) == 0L,
      no_harmonisation_performed = identical(qc$harmonisation_performed, FALSE),
      palindromic_rule_adjudication_required = isTRUE(qc$palindromic_rule_adjudication_required),
      formal_harmonisation_not_approved = identical(qc$approved_for_reverse_relaxed_formal_harmonisation, FALSE),
      no_beta_flip_performed = as.integer(qc$beta_flip_performed_count) == 0L,
      no_eaf_flip_performed = as.integer(qc$eaf_flip_performed_count) == 0L,
      no_strand_flip_performed = as.integer(qc$strand_flip_performed_count) == 0L,
      parquet_tsv_consistency = table_equal(pq, tsv),
      all_target_rsids_unique = !anyDuplicated(pq$target_rsid),
      counts_match_qc = !any(grepl("_count_mismatch$", failures)),
      renv_lock_unchanged = identical(renv_before, renv_after),
      no_reextraction = TRUE,
      no_reverse_mr = TRUE
    ),
    hard_check_failures = failures
  )
  jsonlite::write_json(result, paste0(out[["json"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  if (!identical(freeze_status, "passed")) stop("Reverse relaxed preflight freeze failed; partial outputs retained.", call. = FALSE)
  for (path in out[c("manifest", "json")]) {
    stop_if(file.exists(path), paste("Output appeared during run:", path))
    stop_if(!file.rename(paste0(path, ".partial"), path), paste("Atomic rename failed:", path))
  }
  log_line("freeze_status=passed manifest_sha256=", manifest_sha)
}

tryCatch(main(), error = function(e) {
  log_line("status=failed error=", conditionMessage(e))
  quit(status = 1L)
})
