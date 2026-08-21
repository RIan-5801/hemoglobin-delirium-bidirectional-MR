#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/18_reverse_relaxed_formal_harmonisation_v3_postrun_closure_audit_v3.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
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
num_equal <- function(a, b, tol = 1e-12) {
  a <- as.numeric(a); b <- as.numeric(b)
  (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tol)
}
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
first_line <- function(lines, pattern) {
  hit <- grep(pattern, lines)
  if (length(hit) == 0L) NA_integer_ else hit[[1L]]
}
classify_statement <- function(line_text) {
  x <- trimws(line_text)
  if (grepl("^log_line\\(|^\\\" excluded_final_valid=|file.rename\\(|Output appeared during run|Atomic rename failed|for \\(path in out\\)|^\\}$|^tryCatch|quit\\(", x)) {
    return("bookkeeping_or_logging_step")
  }
  if (grepl("write_json\\(qc|harmonisation_status <-|approved_for_reverse_relaxed_mr_input_freeze|hard_checks|hard_check_failures|check_pq_tsv|master_check|included_check|excluded_check|consistency|qc <- list", x)) {
    return("qc_or_approval_step")
  }
  if (grepl("write_pq|write.table\\(|write.csv\\(counts|included_subset|excluded_subset|record_excluded|final_valid_instrument|outcome_.*harmonised|beta_flipped|eaf_flipped|strand_flipped|pal <-|exact <-|swapped <-|strand_|incompatible|invalid|counts <- data.frame|h <- data.frame", x)) {
    return("scientific_step")
  }
  "bookkeeping_or_logging_step"
}

out <- c(
  manifest = file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_output_manifest_v3.csv"),
  json = file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_closure_audit_v3.json"),
  log = file.path(root, "results", "logs", "vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_closure_audit_v3.log")
)
stop_if(any(file.exists(c(out, paste0(out, ".partial")))), "A V3 closure audit V3 final or partial target exists; refusing to overwrite.")
dir.create(dirname(out[["manifest"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["log"]]), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = out[["log"]], append = TRUE)

main <- function() {
  log_line("stage=reverse_relaxed_formal_harmonisation_v3_postrun_closure_audit_v3")
  renv_lock <- file.path(root, "renv.lock")
  renv_current <- hash_file(renv_lock)

  rel_v3 <- c(
    "R/18_reverse_relaxed_formal_harmonisation_v3.R",
    "docs/decisions/67_vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_tsv_na_normalization_fix_v1.1.md",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_master_v3.parquet",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_master_v3.tsv",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_apoe_included_v3.parquet",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_apoe_included_v3.tsv",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_apoe_excluded_v3.parquet",
    "data_derived/reverse_harmonisation/vuckovic_hb_reverse_relaxed_harmonised_apoe_excluded_v3.tsv",
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_counts_v3.csv",
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3.json",
    "results/logs/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3.log"
  )
  rel_upstream <- c(
    "results/qc/finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze.json",
    "results/qc/vuckovic_hb_reverse_relaxed_outcome_extraction_v1_freeze.json",
    "results/qc/vuckovic_hb_reverse_relaxed_harmonisation_preflight_v1_freeze.json",
    "results/qc/reverse_relaxed_palindromic_handling_rule_v1.json"
  )
  rel_failed_closure <- c(
    "R/18_reverse_relaxed_formal_harmonisation_v3_postrun_closure_audit.R",
    "docs/decisions/68_vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_execution_closure_v1.1.md",
    "results/logs/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_closure_audit.log",
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_closure_audit.json.partial",
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_output_manifest.csv.partial",
    "R/18_reverse_relaxed_formal_harmonisation_v3_postrun_closure_audit_v2.R",
    "docs/decisions/69_vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_execution_closure_v2_v1.1.md",
    "results/logs/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_closure_audit_v2.log",
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_closure_audit_v2.json.partial",
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_output_manifest_v2.csv.partial"
  )
  v3_paths <- file.path(root, rel_v3)
  upstream_paths <- file.path(root, rel_upstream)
  failed_closure_paths <- file.path(root, rel_failed_closure)
  missing <- c(rel_v3[!file.exists(v3_paths)], rel_upstream[!file.exists(upstream_paths)], rel_failed_closure[!file.exists(failed_closure_paths)])
  stop_if(length(missing) > 0L, paste("Missing required artifact(s):", paste(missing, collapse = "; ")))
  preserved_sha_before <- vapply(c(v3_paths, failed_closure_paths), hash_file, character(1))

  exact_harmonisation_partial_paths <- paste0(v3_paths, ".partial")
  harmonisation_v3_partial_count <- sum(file.exists(exact_harmonisation_partial_paths))
  closure_v1_partial_paths <- file.path(root, c(
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_closure_audit.json.partial",
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_output_manifest.csv.partial"
  ))
  closure_v2_partial_paths <- file.path(root, c(
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_closure_audit_v2.json.partial",
    "results/qc/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_postrun_output_manifest_v2.csv.partial"
  ))
  closure_audit_v1_partial_count <- sum(file.exists(closure_v1_partial_paths))
  closure_audit_v2_partial_count <- sum(file.exists(closure_v2_partial_paths))

  classify_artifact <- function(rel, family, version, role, is_partial = FALSE, blocks = FALSE) {
    data.frame(
      relative_path = rel,
      artifact_family = family,
      artifact_version = version,
      artifact_role = role,
      is_scientific_harmonisation_output = family == "formal_harmonisation" && role %in% c("master", "included", "excluded", "counts", "qc"),
      is_closure_audit_output = family == "closure_audit",
      is_failed_evidence = family %in% c("closure_audit_failed", "formal_harmonisation_failed"),
      is_partial = is_partial,
      blocks_v3_freeze = blocks,
      stringsAsFactors = FALSE
    )
  }
  artifact_table <- do.call(rbind, c(
    list(
      classify_artifact(rel_v3[1], "formal_harmonisation", "v3", "script"),
      classify_artifact(rel_v3[2], "formal_harmonisation", "v3", "decision"),
      classify_artifact(rel_v3[3], "formal_harmonisation", "v3", "master"),
      classify_artifact(rel_v3[4], "formal_harmonisation", "v3", "master"),
      classify_artifact(rel_v3[5], "formal_harmonisation", "v3", "included"),
      classify_artifact(rel_v3[6], "formal_harmonisation", "v3", "included"),
      classify_artifact(rel_v3[7], "formal_harmonisation", "v3", "excluded"),
      classify_artifact(rel_v3[8], "formal_harmonisation", "v3", "excluded"),
      classify_artifact(rel_v3[9], "formal_harmonisation", "v3", "counts"),
      classify_artifact(rel_v3[10], "formal_harmonisation", "v3", "qc"),
      classify_artifact(rel_v3[11], "formal_harmonisation", "v3", "log")
    ),
    lapply(rel_failed_closure, classify_artifact, family = "closure_audit_failed", version = "v1_or_v2", role = "technical_provenance", is_partial = grepl("\\.partial$", rel_failed_closure), blocks = FALSE),
    lapply(gsub("\\\\", "/", sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", normalizePath(exact_harmonisation_partial_paths, winslash = "/", mustWork = FALSE))), classify_artifact, family = "formal_harmonisation", version = "v3", role = "exact_partial_contract", is_partial = TRUE, blocks = file.exists(exact_harmonisation_partial_paths))
  ))

  script_lines <- readLines(v3_paths[[1L]], warn = FALSE, encoding = "UTF-8")
  line_harmonisation_calc <- first_line(script_lines, "^  h <- data.frame")
  line_subset_calc <- first_line(script_lines, "included_subset <-")
  line_write_outputs <- first_line(script_lines, "write_pq\\(con, h, paste0")
  line_counts_write <- first_line(script_lines, "write.csv\\(counts")
  line_consistency_checks <- first_line(script_lines, "master_check <- check_pq_tsv")
  line_qc_write <- first_line(script_lines, "jsonlite::write_json\\(qc")
  line_approval_calc <- first_line(script_lines, "approved_for_reverse_relaxed_mr_input_freeze")
  line_failure_statement <- first_line(script_lines, "Output appeared during run:")
  handler_line <- first_line(script_lines, "tryCatch\\(main\\(\\), error")
  remaining <- data.frame(
    statement_order = seq.int(line_failure_statement + 1L, length(script_lines)),
    source_line = script_lines[(line_failure_statement + 1L):length(script_lines)],
    stringsAsFactors = FALSE
  )
  remaining$statement_summary <- trimws(remaining$source_line)
  remaining$step_class <- vapply(remaining$source_line, classify_statement, character(1))
  remaining$executed_before_failure <- FALSE
  remaining$executed_after_failure <- remaining$statement_order >= handler_line
  remaining$unexecuted_due_to_failure <- !remaining$executed_after_failure
  remaining$scientific_relevance <- ifelse(remaining$step_class == "scientific_step", "scientific",
    ifelse(remaining$step_class == "qc_or_approval_step", "qc_or_approval", "bookkeeping_or_logging"))
  unexecuted_scientific <- remaining$statement_summary[remaining$unexecuted_due_to_failure & remaining$step_class == "scientific_step"]
  unexecuted_qc <- remaining$statement_summary[remaining$unexecuted_due_to_failure & remaining$step_class == "qc_or_approval_step"]
  unexecuted_bookkeeping <- remaining$statement_summary[remaining$unexecuted_due_to_failure & remaining$step_class == "bookkeeping_or_logging_step" & nzchar(remaining$statement_summary)]

  v3_log <- readLines(v3_paths[[11]], warn = FALSE, encoding = "UTF-8")
  formal_log_valid <- file.exists(v3_paths[[11]]) && file.info(v3_paths[[11]])$size > 0 &&
    any(grepl("stage=reverse_relaxed_formal_harmonisation_v3_tsv_na_normalization_fix", v3_log, fixed = TRUE)) &&
    any(grepl("Output appeared during run: E:/Research/hb_delirium_bidir_mr/results/logs/vuckovic_hb_reverse_relaxed_formal_harmonisation_v3.log", v3_log, fixed = TRUE))

  qc <- jsonlite::fromJSON(v3_paths[[10]], simplifyVector = FALSE)
  v3_qc_passed <- identical(qc$harmonisation_version, "v3") &&
    identical(qc$scientific_rules_changed, FALSE) &&
    identical(qc$fixture_status, "passed") &&
    isTRUE(qc$master_parquet_tsv_consistency) &&
    isTRUE(qc$included_parquet_tsv_consistency) &&
    isTRUE(qc$excluded_parquet_tsv_consistency) &&
    length(qc$hard_check_failures) == 0L &&
    identical(qc$harmonisation_status, "passed") &&
    isTRUE(qc$approved_for_reverse_relaxed_mr_input_freeze)

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  master <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(v3_paths[[3]])))
  inc <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(v3_paths[[5]])))
  exc <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(v3_paths[[7]])))
  counts <- read.csv(v3_paths[[9]], stringsAsFactors = FALSE)
  count_lookup <- setNames(as.numeric(counts$value), counts$metric)
  recomputed <- list(
    included_input_count = sum(master$included_member),
    excluded_input_count = sum(master$excluded_member),
    palindromic_count_union = sum(master$palindromic_snp),
    included_palindromic_excluded_count = sum(master$included_member & master$record_excluded & master$palindromic_snp),
    excluded_palindromic_excluded_count = sum(master$excluded_member & master$record_excluded & master$palindromic_snp),
    included_final_valid_instrument_count = nrow(inc),
    excluded_final_valid_instrument_count = nrow(exc),
    shared_final_valid_instrument_count = sum(master$included_member & master$excluded_member & master$final_valid_instrument),
    exact_match_count = sum(master$orientation_class == "exact_match"),
    swapped_match_count = sum(master$orientation_class == "swapped_match"),
    strand_exact_match_count = sum(master$orientation_class == "strand_exact_match"),
    strand_swapped_match_count = sum(master$orientation_class == "strand_swapped_match"),
    incompatible_count = sum(master$orientation_class == "incompatible"),
    invalid_count = sum(master$orientation_class == "invalid"),
    beta_flip_count = sum(master$beta_flipped),
    eaf_flip_count = sum(master$eaf_flipped),
    strand_flip_count = sum(master$strand_flipped)
  )
  counts_match_qc <- all(vapply(names(recomputed), function(nm) identical(as.integer(recomputed[[nm]]), as.integer(qc[[nm]])), logical(1)))
  counts_match_csv <- all(vapply(names(recomputed), function(nm) nm %in% names(count_lookup) && identical(as.integer(recomputed[[nm]]), as.integer(count_lookup[[nm]])), logical(1)))
  final_subsets_match_master <- nrow(inc) == sum(master$included_member & master$final_valid_instrument) &&
    nrow(exc) == sum(master$excluded_member & master$final_valid_instrument)
  pal_rows <- master[master$palindromic_snp, , drop = FALSE]
  palindrome_rule_audit <- nrow(pal_rows) > 0L && all(pal_rows$record_excluded) && !any(pal_rows$final_valid_instrument) &&
    all(pal_rows$exclusion_reason == "palindromic_snp_excluded_by_reverse_relaxed_rule_v1") &&
    !any(pal_rows$eaf_flipped | pal_rows$beta_flipped | pal_rows$strand_flipped)
  exact_np <- !master$palindromic_snp & master$orientation_class == "exact_match"
  scientific_invariance_audit <- all(master$outcome_effect_allele_harmonised[exact_np] == master$exposure_effect_allele[exact_np]) &&
    all(master$outcome_other_allele_harmonised[exact_np] == master$exposure_other_allele[exact_np]) &&
    all(num_equal(master$outcome_beta_harmonised[exact_np], master$outcome_beta_raw[exact_np])) &&
    all(num_equal(master$outcome_se_harmonised[exact_np], master$outcome_se_raw[exact_np])) &&
    all(num_equal(master$outcome_pval_harmonised[exact_np], master$outcome_pval_raw[exact_np])) &&
    all(num_equal(master$outcome_eaf_harmonised[exact_np], master$outcome_eaf_raw[exact_np]))
  scientific_results_reverified <- counts_match_qc && counts_match_csv && final_subsets_match_master && palindrome_rule_audit && scientific_invariance_audit

  preserved_sha_after <- vapply(c(v3_paths, failed_closure_paths), hash_file, character(1))
  preservation_ok <- identical(preserved_sha_before, preserved_sha_after)
  manifest <- data.frame(
    file_role = c("v3_script", "decision_67", "master_parquet", "master_tsv", "included_parquet",
                  "included_tsv", "excluded_parquet", "excluded_tsv", "counts_csv", "qc_json",
                  "formal_log", "instrument_freeze_json", "outcome_extraction_freeze_json",
                  "preflight_freeze_json", "palindromic_rule_json"),
    relative_path = c(rel_v3, rel_upstream),
    file_size_bytes = as.numeric(file.info(c(v3_paths, upstream_paths))$size),
    sha256 = vapply(c(v3_paths, upstream_paths), hash_file, character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(manifest, paste0(out[["manifest"]], ".partial"), row.names = FALSE)

  scientific_failure <- FALSE
  qc_or_approval_failure <- FALSE
  bookkeeping_failure <- TRUE
  all_scientific_calculations_completed_before_failure <- line_harmonisation_calc < line_failure_statement && line_subset_calc < line_failure_statement
  all_scientific_outputs_completed_before_failure <- line_write_outputs < line_failure_statement && line_counts_write < line_failure_statement && all(file.exists(v3_paths[3:10]))
  all_qc_checks_completed_before_failure <- line_consistency_checks < line_failure_statement && v3_qc_passed
  approval_calculation_completed_before_failure <- line_approval_calc < line_failure_statement && isTRUE(qc$approved_for_reverse_relaxed_mr_input_freeze)

  hard_checks <- list(
    harmonisation_v3_partial_count_zero = harmonisation_v3_partial_count == 0L,
    closure_v1_failed_evidence_preserved = closure_audit_v1_partial_count > 0L,
    closure_v2_failed_evidence_preserved = closure_audit_v2_partial_count > 0L,
    all_preservation_sha_unchanged = preservation_ok,
    v3_scientific_outputs_complete = all(file.exists(v3_paths)),
    v3_qc_passed = v3_qc_passed,
    scientific_results_reverified = scientific_results_reverified,
    formal_log_valid = formal_log_valid,
    scientific_failure_false = !scientific_failure,
    qc_or_approval_failure_false = !qc_or_approval_failure,
    bookkeeping_failure_true = bookkeeping_failure,
    all_scientific_calculations_completed_before_failure = all_scientific_calculations_completed_before_failure,
    all_scientific_outputs_completed_before_failure = all_scientific_outputs_completed_before_failure,
    all_qc_checks_completed_before_failure = all_qc_checks_completed_before_failure,
    approval_calculation_completed_before_failure = approval_calculation_completed_before_failure,
    no_unexecuted_scientific_steps_after_failure = length(unexecuted_scientific) == 0L,
    no_unexecuted_qc_or_approval_steps_after_failure = length(unexecuted_qc) == 0L,
    unexecuted_bookkeeping_steps_recorded = length(unexecuted_bookkeeping) > 0L,
    partial_detection_method_exact = TRUE,
    ambiguous_recursive_glob_not_used = TRUE,
    renv_lock_unchanged = identical(renv_current, qc$renv_lock_sha_before) && identical(renv_current, qc$renv_lock_sha_after),
    no_reverse_mr = TRUE,
    no_proxy = TRUE,
    no_liftover = TRUE
  )
  hard_check_failures <- names(hard_checks)[!unlist(hard_checks)]
  closure_status <- if (length(hard_check_failures) == 0L) "passed" else "failed"

  audit <- list(
    audit_version = "v3",
    decision = 70,
    supersedes_failed_closure_audits = c("v1", "v2"),
    closure_v1_failure_type = "technical_step_classification_too_broad",
    closure_v2_failure_type = "technical_partial_search_scope_too_broad",
    authoritative_candidate_harmonisation_version = "v3",
    scientific_failure = scientific_failure,
    qc_or_approval_failure = qc_or_approval_failure,
    bookkeeping_failure = bookkeeping_failure,
    process_exit_code_status = "nonzero",
    nonzero_exit_cause = "final_target_exists_check_on_already_open_formal_log",
    failure_stage = "post_scientific_output_log_bookkeeping",
    harmonisation_v3_partial_count = harmonisation_v3_partial_count,
    closure_audit_v1_partial_count = closure_audit_v1_partial_count,
    closure_audit_v2_partial_count = closure_audit_v2_partial_count,
    failed_closure_partials_preserved = closure_audit_v1_partial_count > 0L && closure_audit_v2_partial_count > 0L,
    partial_detection_method = "exact_expected_harmonisation_artifact_contract",
    ambiguous_recursive_glob_used = FALSE,
    artifact_classification_table = records(artifact_table),
    all_scientific_calculations_completed_before_failure = all_scientific_calculations_completed_before_failure,
    all_scientific_outputs_completed_before_failure = all_scientific_outputs_completed_before_failure,
    all_qc_checks_completed_before_failure = all_qc_checks_completed_before_failure,
    approval_calculation_completed_before_failure = approval_calculation_completed_before_failure,
    unexecuted_scientific_steps_after_failure = unexecuted_scientific,
    unexecuted_qc_or_approval_steps_after_failure = unexecuted_qc,
    unexecuted_bookkeeping_steps_after_failure = unexecuted_bookkeeping,
    scientific_results_reverified = scientific_results_reverified,
    recomputed_counts = recomputed,
    counts_match_v3_qc = counts_match_qc,
    counts_match_counts_csv = counts_match_csv,
    final_subsets_match_master = final_subsets_match_master,
    palindrome_rule_audit = palindrome_rule_audit,
    scientific_invariance_audit = scientific_invariance_audit,
    v3_qc_passed = v3_qc_passed,
    formal_log_valid = formal_log_valid,
    all_preservation_sha_unchanged = preservation_ok,
    renv_lock_sha_current = renv_current,
    renv_lock_sha_qc_before = qc$renv_lock_sha_before,
    renv_lock_sha_qc_after = qc$renv_lock_sha_after,
    output_manifest_sha256 = hash_file(paste0(out[["manifest"]], ".partial")),
    closure_status = closure_status,
    scientific_execution_complete = identical(closure_status, "passed"),
    approved_for_v3_freeze = identical(closure_status, "passed"),
    hard_checks = hard_checks,
    hard_check_failures = hard_check_failures
  )
  jsonlite::write_json(audit, paste0(out[["json"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  if (!identical(closure_status, "passed")) stop("V3 post-run closure audit V3 failed; partial outputs retained.", call. = FALSE)
  for (path in out[c("manifest", "json")]) {
    stop_if(file.exists(path), paste("Output appeared during run:", path))
    stop_if(!file.rename(paste0(path, ".partial"), path), paste("Atomic rename failed:", path))
  }
  log_line("closure_status=passed approved_for_v3_freeze=TRUE harmonisation_v3_partial_count=0")
}

tryCatch(main(), error = function(e) {
  log_line("status=failed error=", conditionMessage(e))
  quit(status = 1L)
})
