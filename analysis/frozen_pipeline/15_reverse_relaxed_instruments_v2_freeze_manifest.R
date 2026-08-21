#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/15_reverse_relaxed_instruments_v2_freeze_manifest.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

for (pkg in c("DBI", "duckdb", "jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
hash_file <- function(path) digest::digest(file = path, algo = "sha256")
sql_string <- function(path, must_work = TRUE) paste0("'", gsub("'", "''", normalizePath(path, winslash = "/", mustWork = must_work), fixed = TRUE), "'")
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

out <- c(
  manifest = file.path(root, "results", "qc", "finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze_manifest.csv"),
  json = file.path(root, "results", "qc", "finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze.json"),
  log = file.path(root, "results", "logs", "finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze.log")
)
stop_if(any(file.exists(c(out, paste0(out, ".partial")))), "A relaxed V2 freeze target or partial exists; refusing to overwrite.")
dir.create(dirname(out[["manifest"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["log"]]), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = out[["log"]], append = TRUE)

main <- function() {
  log_line("stage=reverse_relaxed_instruments_v2_freeze")
  rel <- c(
    "docs/decisions/56_vuckovic_hb_reverse_primary_formal_harmonisation_v4_freeze_v1.1.md",
    "docs/decisions/58_finngen_r13_delirium_reverse_exploratory_relaxed_instrument_selection_v2_technical_recovery_v1.1.md",
    "docs/decisions/59_finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze_v1.1.md",
    "R/15_reverse_relaxed_instruments_finngen_r13_delirium_p5e-6_v2.R",
    "data_derived/reverse_instruments/finngen_r13_delirium_p5e-6_exploratory_candidates_v2.parquet",
    "data_derived/reverse_instruments/finngen_r13_delirium_p5e-6_exploratory_eligible_candidates_v2.parquet",
    "data_derived/reverse_instruments/finngen_r13_delirium_p5e-6_exploratory_apoe_included_clumped_v2.parquet",
    "data_derived/reverse_instruments/finngen_r13_delirium_p5e-6_exploratory_apoe_included_clumped_v2.tsv",
    "data_derived/reverse_instruments/finngen_r13_delirium_p5e-6_exploratory_apoe_excluded_clumped_v2.parquet",
    "data_derived/reverse_instruments/finngen_r13_delirium_p5e-6_exploratory_apoe_excluded_clumped_v2.tsv",
    "results/qc/finngen_r13_delirium_p5e-6_exploratory_instrument_selection_counts_v2.csv",
    "results/qc/finngen_r13_delirium_p5e-6_exploratory_not_in_reference_v2.tsv",
    "results/qc/finngen_r13_delirium_p5e-6_exploratory_schema_audit_v2.csv",
    "results/qc/finngen_r13_delirium_p5e-6_exploratory_primary_overlap_audit_v2.csv",
    "results/qc/finngen_r13_delirium_p5e-6_exploratory_instrument_selection_v2.json",
    "results/logs/finngen_r13_delirium_p5e-6_exploratory_instrument_selection_v2.log",
    "results/qc/vuckovic_hb_reverse_primary_formal_harmonisation_v4_freeze.json",
    "results/qc/vuckovic_hb_reverse_primary_formal_harmonisation_v4_freeze_manifest.csv",
    "results/qc/finngen_R13_F5_DELIRIUM_input_certification_v1.json",
    "results/qc/ld_reference_manifest_v1.csv",
    "tools/plink2/plink2.exe",
    "renv.lock"
  )
  paths <- file.path(root, rel)
  stop_if(any(!file.exists(paths)), paste("Missing freeze input(s):", paste(rel[!file.exists(paths)], collapse = "; ")))

  qc_path <- file.path(root, "results", "qc", "finngen_r13_delirium_p5e-6_exploratory_instrument_selection_v2.json")
  qc <- jsonlite::fromJSON(qc_path, simplifyVector = FALSE)
  trigger <- jsonlite::fromJSON(file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_formal_harmonisation_v4_freeze.json"), simplifyVector = FALSE)

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  read_pq <- function(path) DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(path)))
  included <- read_pq(file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_p5e-6_exploratory_apoe_included_clumped_v2.parquet"))
  excluded <- read_pq(file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_p5e-6_exploratory_apoe_excluded_clumped_v2.parquet"))

  included_rsids <- as.character(included$rsid)
  excluded_rsids <- as.character(excluded$rsid)
  shared_rsids <- intersect(included_rsids, excluded_rsids)
  included_only <- setdiff(included_rsids, excluded_rsids)
  excluded_only <- setdiff(excluded_rsids, included_rsids)
  union_rsids <- union(included_rsids, excluded_rsids)
  calc_strength <- function(x) {
    f <- (as.numeric(x$beta) / as.numeric(x$se))^2
    list(nSNP = length(f), F_min = min(f), F_mean = mean(f), F_median = median(f), F_max = max(f), F_lt_10_count = sum(f < 10))
  }
  fi <- calc_strength(included)
  fe <- calc_strength(excluded)

  failures <- character()
  add_fail <- function(x) failures <<- unique(c(failures, x))
  if (!identical(qc$instrument_selection_version, "v2")) add_fail("instrument_selection_version_not_v2")
  if (!identical(qc$instrument_selection_status, "passed")) add_fail("instrument_selection_status_not_passed")
  if (!isTRUE(qc$approved_for_reverse_relaxed_outcome_extraction)) add_fail("not_approved_for_outcome_extraction")
  if (length(qc$hard_check_failures) != 0L) add_fail("source_hard_check_failures_present")
  if (!identical(qc$analysis_direction, "delirium_to_Hb")) add_fail("analysis_direction_mismatch")
  if (!identical(qc$analysis_role, "secondary_reverse_exploratory_relaxed")) add_fail("analysis_role_mismatch")
  if (abs(as.numeric(qc$p_threshold) - 5e-6) > 0) add_fail("p_threshold_mismatch")
  if (abs(as.numeric(qc$ld_r2) - 0.001) > 0 || as.integer(qc$ld_window_kb) != 10000L) add_fail("ld_parameter_mismatch")
  if (!identical(qc$proxy_used, FALSE) || !identical(qc$liftover_used, FALSE)) add_fail("proxy_or_liftover_not_false")
  if (!identical(trigger$freeze_status, "passed") || !isTRUE(trigger$approved_for_reverse_relaxed_threshold_branch)) add_fail("protocol_trigger_freeze_gate_failed")
  if (nrow(included) != as.integer(qc$included_nsnp)) add_fail("included_nsnp_mismatch")
  if (nrow(excluded) != as.integer(qc$excluded_nsnp)) add_fail("excluded_nsnp_mismatch")
  if (length(shared_rsids) != as.integer(qc$shared_nsnp)) add_fail("shared_nsnp_mismatch")
  if (length(included_only) != as.integer(qc$included_only_nsnp)) add_fail("included_only_nsnp_mismatch")
  if (length(excluded_only) != as.integer(qc$excluded_only_nsnp)) add_fail("excluded_only_nsnp_mismatch")
  if (fi$F_lt_10_count != 0L || fe$F_lt_10_count != 0L) add_fail("F_lt_10_present")
  if (fi$F_lt_10_count != as.integer(qc$instrument_strength_included$F_lt_10_count)) add_fail("included_F_lt_10_qc_mismatch")
  if (fe$F_lt_10_count != as.integer(qc$instrument_strength_excluded$F_lt_10_count)) add_fail("excluded_F_lt_10_qc_mismatch")
  if (!isTRUE(qc$included_parquet_tsv_consistency) || !isTRUE(qc$excluded_parquet_tsv_consistency)) add_fail("parquet_tsv_consistency_not_true")

  manifest <- data.frame(
    file_role = c(
      "trigger_freeze_decision", "v2_decision", "v2_freeze_decision", "v2_script",
      "v2_candidates_parquet", "v2_eligible_candidates_parquet", "v2_included_parquet", "v2_included_tsv",
      "v2_excluded_parquet", "v2_excluded_tsv", "v2_counts_csv", "v2_not_in_reference_audit",
      "v2_schema_audit", "v2_strict_primary_overlap_audit", "v2_qc_json", "v2_log",
      "formal_harmonisation_v4_freeze_json", "formal_harmonisation_v4_freeze_manifest",
      "finngen_source_certification", "ld_reference_manifest", "plink2_binary", "renv_lock"
    ),
    relative_path = rel,
    file_size_bytes = as.numeric(file.info(paths)$size),
    sha256 = vapply(paths, hash_file, character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(manifest, paste0(out[["manifest"]], ".partial"), row.names = FALSE)
  manifest_sha <- hash_file(paste0(out[["manifest"]], ".partial"))
  freeze_status <- if (length(failures) == 0L) "passed" else "failed"
  freeze <- list(
    freeze_version = "v1",
    decision = 59,
    authoritative_reverse_relaxed_instrument_version = "v2",
    source_instrument_selection_decision = 58,
    analysis_direction = "delirium_to_Hb",
    analysis_role = "secondary_reverse_exploratory_relaxed",
    branch_type = "protocol_prespecified_fallback",
    primary_threshold = 5e-8,
    p_threshold = 5e-6,
    ld_r2 = 0.001,
    ld_window_kb = 10000,
    proxy_used = FALSE,
    liftover_used = FALSE,
    included_nsnp = nrow(included),
    excluded_nsnp = nrow(excluded),
    shared_nsnp = length(shared_rsids),
    included_only_nsnp = length(included_only),
    excluded_only_nsnp = length(excluded_only),
    union_nsnp = length(union_rsids),
    included_rsids = included_rsids,
    excluded_rsids = excluded_rsids,
    shared_rsids = shared_rsids,
    included_only_rsids = included_only,
    excluded_only_rsids = excluded_only,
    included_F_summary_recomputed = fi,
    excluded_F_summary_recomputed = fe,
    manifest_sha256 = manifest_sha,
    freeze_status = freeze_status,
    approved_for_reverse_relaxed_outcome_extraction = identical(freeze_status, "passed"),
    hard_checks = list(
      source_qc_passed = identical(qc$instrument_selection_status, "passed"),
      source_qc_hard_check_failures_empty = length(qc$hard_check_failures) == 0L,
      source_qc_approved_for_extraction = isTRUE(qc$approved_for_reverse_relaxed_outcome_extraction),
      protocol_trigger_freeze_gate = !("protocol_trigger_freeze_gate_failed" %in% failures),
      p_threshold_matches = abs(as.numeric(qc$p_threshold) - 5e-6) == 0,
      ld_parameters_match = abs(as.numeric(qc$ld_r2) - 0.001) == 0 && as.integer(qc$ld_window_kb) == 10000L,
      no_proxy = identical(qc$proxy_used, FALSE),
      no_liftover = identical(qc$liftover_used, FALSE),
      dynamic_instrument_counts_match_qc = !any(c("included_nsnp_mismatch", "excluded_nsnp_mismatch", "shared_nsnp_mismatch", "included_only_nsnp_mismatch", "excluded_only_nsnp_mismatch") %in% failures),
      recomputed_F_lt_10_absent = fi$F_lt_10_count == 0L && fe$F_lt_10_count == 0L,
      parquet_tsv_consistency = isTRUE(qc$included_parquet_tsv_consistency) && isTRUE(qc$excluded_parquet_tsv_consistency)
    ),
    hard_check_failures = failures
  )
  jsonlite::write_json(freeze, paste0(out[["json"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  if (!identical(freeze_status, "passed")) stop("Relaxed V2 freeze failed; partial outputs retained.", call. = FALSE)
  for (path in out[c("manifest", "json")]) {
    stop_if(file.exists(path), paste("Output appeared during run:", path))
    stop_if(!file.rename(paste0(path, ".partial"), path), paste("Atomic rename failed:", path))
  }
  log_line("freeze_status=passed manifest_sha256=", manifest_sha, " union_nsnp=", length(union_rsids))
}

tryCatch(main(), error = function(e) {
  log_line("status=failed error=", conditionMessage(e))
  quit(status = 1L)
})
