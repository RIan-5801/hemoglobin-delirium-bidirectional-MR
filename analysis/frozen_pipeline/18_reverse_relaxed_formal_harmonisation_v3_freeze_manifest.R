#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/18_reverse_relaxed_formal_harmonisation_v3_freeze_manifest.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
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
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
f_summary <- function(x) {
  f <- (as.numeric(x$exposure_beta) / as.numeric(x$exposure_se))^2
  list(
    n = length(f),
    min = min(f),
    mean = mean(f),
    median = stats::median(f),
    max = max(f),
    F_lt10 = sum(f < 10),
    values_by_rsid = records(data.frame(target_rsid = x$target_rsid, F = f, stringsAsFactors = FALSE))
  )
}

out <- c(
  manifest = file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_freeze_manifest.csv"),
  json = file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_freeze.json"),
  log = file.path(root, "results", "logs", "vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_freeze.log")
)
stop_if(any(file.exists(c(out, paste0(out, ".partial")))), "A V3 freeze final or partial target exists; refusing to overwrite.")
dir.create(dirname(out[["manifest"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["log"]]), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = out[["log"]], append = TRUE)

main <- function() {
  log_line("stage=reverse_relaxed_formal_harmonisation_v3_freeze")
  rel <- c(
    "docs/decisions/59_finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze_v1.1.md",
    "docs/decisions/61_vuckovic_hb_reverse_relaxed_outcome_extraction_v1_freeze_v1.1.md",
    "docs/decisions/63_vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1_freeze_v1.1.md",
    "docs/decisions/64_reverse_relaxed_palindromic_handling_rule_v1_v1.1.md",
    "docs/decisions/65_vuckovic_hb_reverse_relaxed_formal_harmonisation_v1_v1.1.md",
    "docs/decisions/66_vuckovic_hb_reverse_relaxed_formal_harmonisation_v2_readback_fix_v1.1.md",
    "docs/decisions/67_vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_tsv_na_normalization_fix_v1.1.md",
    "docs/decisions/68_vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_execution_closure_v1.1.md",
    "docs/decisions/69_vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_execution_closure_v2_v1.1.md",
    "docs/decisions/70_vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_execution_closure_v3_v1.1.md",
    "docs/decisions/71_vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_freeze_v1.1.md",
    "R/18_reverse_relaxed_formal_harmonisation_v3.R",
    "R/18_reverse_relaxed_formal_harmonisation_v3_postrun_closure_audit_v3.R",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_master_v3.parquet",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_master_v3.tsv",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_apoe_included_v3.parquet",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_apoe_included_v3.tsv",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_apoe_excluded_v3.parquet",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_apoe_excluded_v3.tsv",
    "results/qc/finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze.json",
    "results/qc/vuckovic_hb_reverse_relaxed_outcome_extraction_v1_freeze.json",
    "results/qc/vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1_freeze.json",
    "results/qc/reverse_relaxed_palindromic_handling_rule_v1.json",
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_counts_v3.csv",
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3.json",
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_closure_audit_v3.json",
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_output_manifest_v3.csv",
    "results/logs/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3.log",
    "results/logs/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_closure_audit_v3.log",
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_closure_audit.json.partial",
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_output_manifest.csv.partial",
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_closure_audit_v2.json.partial",
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_output_manifest_v2.csv.partial",
    "renv.lock"
  )
  roles <- c(
    "decision_59_instrument_freeze", "decision_61_outcome_freeze", "decision_63_preflight_freeze",
    "decision_64_palindrome_rule", "decision_65_failed_v1", "decision_66_failed_v2",
    "decision_67_formal_harmonisation_v3", "decision_68_failed_closure_v1",
    "decision_69_failed_closure_v2", "decision_70_passed_closure_v3", "decision_71_freeze",
    "v3_script", "closure_v3_script", "master_parquet", "master_tsv", "included_parquet",
    "included_tsv", "excluded_parquet", "excluded_tsv", "instrument_freeze_json",
    "outcome_freeze_json", "preflight_freeze_json", "palindrome_rule_json", "counts_v3",
    "qc_v3", "closure_v3_qc", "closure_v3_manifest", "v3_log", "closure_v3_log",
    "failed_closure_v1_qc_partial", "failed_closure_v1_manifest_partial",
    "failed_closure_v2_qc_partial", "failed_closure_v2_manifest_partial", "renv_lock"
  )
  paths <- file.path(root, rel)
  stop_if(any(!file.exists(paths)), paste("Missing freeze input(s):", paste(rel[!file.exists(paths)], collapse = "; ")))
  qc <- jsonlite::fromJSON(file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_formal_harmonisation_v3.json"), simplifyVector = FALSE)
  closure <- jsonlite::fromJSON(file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_closure_audit_v3.json"), simplifyVector = FALSE)
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  master <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_relaxed_harmonised_master_v3.parquet"))))
  inc <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_relaxed_harmonised_apoe_included_v3.parquet"))))
  exc <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(file.path(root, "data_derived", "reverse_harmonisation", "vuckovic_hb_reverse_relaxed_harmonised_apoe_excluded_v3.parquet"))))
  included_final_rsids <- inc$target_rsid
  excluded_final_rsids <- exc$target_rsid
  shared_final <- sort(intersect(included_final_rsids, excluded_final_rsids))
  strength_included <- f_summary(inc)
  strength_excluded <- f_summary(exc)
  manifest <- data.frame(
    file_role = roles,
    relative_path = rel,
    file_size_bytes = as.numeric(file.info(paths)$size),
    sha256 = vapply(paths, hash_file, character(1)),
    authority_class = ifelse(grepl("failed|partial", roles), "technical_provenance_non_authoritative", "authoritative_or_provenance"),
    stringsAsFactors = FALSE
  )
  write.csv(manifest, paste0(out[["manifest"]], ".partial"), row.names = FALSE)
  manifest_sha <- hash_file(paste0(out[["manifest"]], ".partial"))

  recomputed_counts <- list(
    included_final_valid_instrument_count = length(included_final_rsids),
    excluded_final_valid_instrument_count = length(excluded_final_rsids),
    shared_final_valid_instrument_count = length(shared_final)
  )
  hard_checks <- list(
    closure_v3_passed = identical(closure$closure_status, "passed") && isTRUE(closure$approved_for_v3_freeze),
    harmonisation_v3_partial_count_zero = as.integer(closure$harmonisation_v3_partial_count) == 0L,
    v3_qc_passed = identical(qc$harmonisation_status, "passed") && length(qc$hard_check_failures) == 0L,
    v3_outputs_exist = all(file.exists(paths[roles %in% c("master_parquet", "master_tsv", "included_parquet", "included_tsv", "excluded_parquet", "excluded_tsv", "counts_v3", "qc_v3", "v3_log")])),
    counts_match_qc = identical(as.integer(recomputed_counts$included_final_valid_instrument_count), as.integer(qc$included_final_valid_instrument_count)) &&
      identical(as.integer(recomputed_counts$excluded_final_valid_instrument_count), as.integer(qc$excluded_final_valid_instrument_count)) &&
      identical(as.integer(recomputed_counts$shared_final_valid_instrument_count), as.integer(qc$shared_final_valid_instrument_count)),
    final_subset_membership = all(inc$included_member & inc$final_valid_instrument) && all(exc$excluded_member & exc$final_valid_instrument),
    final_rsids_read_successfully = length(included_final_rsids) > 0L && length(excluded_final_rsids) > 0L,
    instrument_strength_completed = !is.na(strength_included$min) && !is.na(strength_excluded$min),
    included_F_lt10_zero = as.integer(strength_included$F_lt10) == 0L,
    excluded_F_lt10_zero = as.integer(strength_excluded$F_lt10) == 0L,
    manifest_complete = nrow(manifest) == length(rel),
    renv_lock_unchanged = identical(qc$renv_lock_sha_before, qc$renv_lock_sha_after),
    no_reverse_mr = TRUE,
    no_proxy = TRUE,
    no_liftover = TRUE
  )
  hard_check_failures <- names(hard_checks)[!unlist(hard_checks)]
  freeze_status <- if (length(hard_check_failures) == 0L) "passed" else "failed"
  result <- list(
    freeze_version = "v1",
    decision = 71,
    authoritative_reverse_relaxed_formal_harmonisation_version = "v3",
    analysis_direction = "delirium_to_Hb",
    analysis_role = "secondary_reverse_exploratory_relaxed",
    p_threshold = 5e-6,
    source_instrument_version = "v2",
    source_outcome_extraction_version = "v1",
    source_preflight_version = "v1",
    palindromic_rule_version = "v1",
    postrun_closure_audit_version = "v3",
    included_final_valid_instrument_count = length(included_final_rsids),
    excluded_final_valid_instrument_count = length(excluded_final_rsids),
    shared_final_valid_instrument_count = length(shared_final),
    included_final_rsids = included_final_rsids,
    excluded_final_rsids = excluded_final_rsids,
    shared_final_rsids = shared_final,
    instrument_strength_included = strength_included,
    instrument_strength_excluded = strength_excluded,
    manifest_sha256 = manifest_sha,
    freeze_status = freeze_status,
    approved_for_reverse_relaxed_mr_design = identical(freeze_status, "passed"),
    hard_checks = hard_checks,
    hard_check_failures = hard_check_failures,
    informational_findings = list(
      original_v3_process_nonzero_exit = TRUE,
      original_v3_nonzero_exit_note = "Original Formal Harmonisation V3 R process returned a nonzero exit code due solely to post-scientific-output log bookkeeping.",
      independent_closure_audit_v3_confirmation = c(
        "no unexecuted scientific steps",
        "no unexecuted QC/approval steps",
        "all scientific outputs complete",
        "all scientific hard checks passed",
        "only bookkeeping/logging statements remained",
        "no Formal Harmonisation V3 scientific partials existed"
      ),
      failed_closure_partials_retained_as_technical_provenance = TRUE,
      no_reverse_mr_run = TRUE
    )
  )
  jsonlite::write_json(result, paste0(out[["json"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  if (!identical(freeze_status, "passed")) stop("V3 freeze failed; partial outputs retained.", call. = FALSE)
  for (path in out[c("manifest", "json")]) {
    stop_if(file.exists(path), paste("Output appeared during run:", path))
    stop_if(!file.rename(paste0(path, ".partial"), path), paste("Atomic rename failed:", path))
  }
  log_line("freeze_status=passed approved_for_reverse_relaxed_mr_design=TRUE manifest_sha256=", manifest_sha)
}

tryCatch(main(), error = function(e) {
  log_line("status=failed error=", conditionMessage(e))
  quit(status = 1L)
})
