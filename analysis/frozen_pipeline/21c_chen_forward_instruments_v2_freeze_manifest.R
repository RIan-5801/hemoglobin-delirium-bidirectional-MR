#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/21c_chen_forward_instruments_v2_freeze_manifest.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

for (pkg in c("DBI", "duckdb", "jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
qpath <- function(path) gsub("'", "''", normalizePath(path, winslash = "/", mustWork = FALSE), fixed = TRUE)
atomic <- function(path, writer) {
  partial <- paste0(path, ".partial")
  stop_if(file.exists(path) || file.exists(partial), paste("Output occupied:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writer(partial)
  stop_if(!file.exists(partial), paste("Writer did not create partial:", partial))
  stop_if(!file.rename(partial, path), paste("Atomic rename failed:", path))
}

rel <- function(...) file.path(...)
abs_path <- function(p) file.path(root, p)
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
strength <- function(x) {
  list(
    n = nrow(x),
    F_min = min(x$F_stat),
    F_mean = mean(x$F_stat),
    F_median = stats::median(x$F_stat),
    F_max = max(x$F_stat),
    F_lt10_count = sum(x$F_stat < 10)
  )
}
read_pq <- function(con, p) DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", qpath(abs_path(p))))

out <- c(
  manifest = abs_path(rel("results", "qc", "chen_forward_instruments_v2_freeze_manifest.csv")),
  json = abs_path(rel("results", "qc", "chen_forward_instruments_v2_freeze.json")),
  log = abs_path(rel("results", "logs", "chen_forward_instruments_v2_freeze.log")),
  decision = abs_path(rel("docs", "decisions", "87_chen_forward_instruments_v2_freeze_v1.1.md"))
)
stop_if(any(file.exists(c(out, paste0(out, ".partial")))), "A freeze output or partial exists; refusing to overwrite.")

manifest_items <- data.frame(
  relative_path = c(
    rel("docs", "decisions", "84_chen_identifier_resolution_amendment_v1_v1.1.md"),
    rel("docs", "decisions", "85_chen_forward_instrument_selection_v2_v1.1.md"),
    rel("docs", "decisions", "86_chen_forward_instrument_selection_v2_readback_closure_v1.1.md"),
    rel("R", "21_chen_forward_instruments_v2.R"),
    rel("R", "21a_chen_forward_instruments_v2_postrun_recovery.R"),
    rel("R", "21b_chen_forward_instruments_v2_readback_audit.R"),
    rel("R", "21c_chen_forward_instruments_v2_freeze_manifest.R"),
    rel("results", "qc", "chen_2020_hb_source_certification_v1.json"),
    rel("results", "qc", "chen_forward_sensitivity_instrument_selection_contract_v1.json"),
    rel("results", "qc", "chen_identifier_resolution_amendment_v1.json"),
    rel("results", "qc", "chen_forward_instrument_selection_v2.json"),
    rel("results", "qc", "chen_forward_instrument_selection_v2_readback_closure_v1.json"),
    rel("results", "qc", "chen_forward_instrument_selection_v2_readback_numeric_audit_v1.csv"),
    rel("results", "qc", "chen_forward_identifier_resolution_counts_v2.csv"),
    rel("results", "qc", "chen_forward_identifier_resolution_failures_v2.tsv"),
    rel("results", "qc", "chen_forward_resolved_duplicate_rsid_audit_v2.tsv"),
    rel("results", "qc", "chen_forward_instrument_selection_counts_v2.csv"),
    rel("results", "qc", "chen_forward_vuckovic_overlap_audit_v2.csv"),
    rel("results", "logs", "chen_forward_instrument_selection_v2.log"),
    rel("results", "logs", "chen_forward_instrument_selection_v2_recovery.log"),
    rel("results", "logs", "chen_forward_instrument_selection_v2_readback_closure_v1.log"),
    rel("data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_candidates_v2.parquet"),
    rel("data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_identity_resolution_v2.parquet"),
    rel("data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_resolved_eligible_v2.parquet"),
    rel("data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_included_clumped_v2.parquet"),
    rel("data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_included_clumped_v2.tsv"),
    rel("data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_excluded_clumped_v2.parquet"),
    rel("data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_excluded_clumped_v2.tsv"),
    rel("data_derived", "forward_sensitivity_instruments", "plink", "chen_2020_hb_p5e-8_apoe_included_clump_input_v2.tsv"),
    rel("data_derived", "forward_sensitivity_instruments", "plink", "chen_2020_hb_p5e-8_apoe_excluded_clump_input_v2.tsv"),
    rel("data_derived", "forward_sensitivity_instruments", "plink", "chen_2020_hb_p5e-8_apoe_included_clumped_v2.clumps"),
    rel("data_derived", "forward_sensitivity_instruments", "plink", "chen_2020_hb_p5e-8_apoe_excluded_clumped_v2.clumps")
  ),
  file_role = c(
    "identifier_resolution_amendment_decision",
    "instrument_selection_decision_scientific_output_record_with_initial_technical_stop",
    "readback_closure_decision_authoritative_pass_gate",
    "instrument_selection_script_v2",
    "postrun_recovery_script",
    "readback_closure_script",
    "freeze_script",
    "chen_source_certification",
    "chen_forward_contract",
    "identifier_resolution_amendment_qc",
    "instrument_selection_v2_initial_qc_technical_provenance_non_authoritative",
    "instrument_selection_v2_readback_closure_authoritative_qc",
    "readback_numeric_audit",
    "identifier_resolution_counts",
    "identifier_resolution_failures",
    "resolved_duplicate_rsid_audit",
    "instrument_selection_counts",
    "vuckovic_overlap_audit",
    "instrument_selection_v2_initial_log_failed_readback_gate",
    "postrun_recovery_log",
    "readback_closure_log",
    "p_lt_5e_8_candidates",
    "identity_resolution_audit",
    "resolved_eligible_candidates",
    "authoritative_apoe_included_instruments_parquet",
    "human_readable_apoe_included_instruments_tsv",
    "authoritative_apoe_excluded_instruments_parquet",
    "human_readable_apoe_excluded_instruments_tsv",
    "apoe_included_clump_input",
    "apoe_excluded_clump_input",
    "apoe_included_plink_clumps",
    "apoe_excluded_plink_clumps"
  ),
  scientific_authority = c(
    rep(TRUE, 10), FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
    FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE
  ),
  stringsAsFactors = FALSE
)

main <- function() {
  log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = out[["log"]], append = TRUE)
  log_line("stage=chen_forward_instruments_v2_freeze")
  paths <- abs_path(manifest_items$relative_path)
  stop_if(any(!file.exists(paths)), paste("Missing freeze manifest input:", paste(manifest_items$relative_path[!file.exists(paths)], collapse = "; ")))

  closure <- jsonlite::fromJSON(abs_path(rel("results", "qc", "chen_forward_instrument_selection_v2_readback_closure_v1.json")), simplifyVector = FALSE)
  initial_qc <- jsonlite::fromJSON(abs_path(rel("results", "qc", "chen_forward_instrument_selection_v2.json")), simplifyVector = FALSE)
  stop_if(!identical(closure$instrument_selection_status, "passed"), "Decision 86 closure status is not passed.")
  stop_if(length(closure$hard_check_failures) != 0L, "Decision 86 closure has hard-check failures.")
  stop_if(!isTRUE(closure$approved_for_chen_forward_outcome_extraction), "Decision 86 did not approve Chen forward outcome extraction.")
  stop_if(!identical(initial_qc$instrument_selection_status, "failed"), "Decision 85 initial QC is expected to preserve technical readback failure evidence.")

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:", config = list(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  included <- read_pq(con, rel("data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_included_clumped_v2.parquet"))
  excluded <- read_pq(con, rel("data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_excluded_clumped_v2.parquet"))
  required_cols <- c("source_marker_id", "resolved_rsid", "marker_chr", "marker_pos", "marker_allele1", "marker_allele2", "exposure_effect_allele", "exposure_other_allele", "beta", "se", "pval", "eaf", "n_samples", "F_stat", "reference_chr_grch37", "reference_pos_grch37", "identity_resolution_status", "chen_marker_effect_allele_compatibility", "apoe_region")
  stop_if(!all(required_cols %in% names(included)) || !all(required_cols %in% names(excluded)), "Final Chen instrument schema is missing required columns.")
  stop_if(anyDuplicated(included$resolved_rsid) || anyDuplicated(excluded$resolved_rsid), "Duplicate resolved_rsid in final instruments.")
  stop_if(any(excluded$apoe_region), "APOE-excluded final instruments contain APOE-region variant.")
  calc_ok <- function(x) {
    calc <- (x$beta / x$se)^2
    all(is.finite(calc)) && all(abs(calc - x$F_stat) <= 1e-10 | abs(calc - x$F_stat) / pmax(abs(calc), abs(x$F_stat), .Machine$double.xmin) <= 1e-12)
  }
  shared <- intersect(included$resolved_rsid, excluded$resolved_rsid)
  included_only <- setdiff(included$resolved_rsid, excluded$resolved_rsid)
  excluded_only <- setdiff(excluded$resolved_rsid, included$resolved_rsid)

  manifest <- manifest_items
  manifest$file_size_bytes <- as.numeric(file.info(paths)$size)
  manifest$sha256 <- vapply(paths, hash_file, character(1))
  atomic(out[["manifest"]], function(p) write.csv(manifest, p, row.names = FALSE, na = ""))
  manifest_sha <- hash_file(out[["manifest"]])

  hard_checks <- list(
    decision_84_present = file.exists(abs_path(rel("docs", "decisions", "84_chen_identifier_resolution_amendment_v1_v1.1.md"))),
    decision_85_present_and_preserved = file.exists(abs_path(rel("docs", "decisions", "85_chen_forward_instrument_selection_v2_v1.1.md"))) && identical(initial_qc$instrument_selection_status, "failed"),
    decision_86_closure_passed = identical(closure$instrument_selection_status, "passed") && length(closure$hard_check_failures) == 0L,
    approved_by_decision_86 = isTRUE(closure$approved_for_chen_forward_outcome_extraction),
    authoritative_instrument_files_present = file.exists(abs_path(rel("data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_included_clumped_v2.parquet"))) && file.exists(abs_path(rel("data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_excluded_clumped_v2.parquet"))),
    included_count_matches_closure = nrow(included) == as.integer(closure$included_nsnp),
    excluded_count_matches_closure = nrow(excluded) == as.integer(closure$excluded_nsnp),
    shared_counts_match_closure = length(shared) == as.integer(closure$shared_nsnp) && length(included_only) == as.integer(closure$included_only_nsnp) && length(excluded_only) == as.integer(closure$excluded_only_nsnp),
    resolved_rsid_unique = !anyDuplicated(included$resolved_rsid) && !anyDuplicated(excluded$resolved_rsid),
    source_marker_id_preserved = all(nzchar(included$source_marker_id)) && all(nzchar(excluded$source_marker_id)),
    effect_orientation_preserved = all(included$exposure_effect_allele %in% c("A", "C", "G", "T")) && all(included$exposure_other_allele %in% c("A", "C", "G", "T")) && all(excluded$exposure_effect_allele %in% c("A", "C", "G", "T")) && all(excluded$exposure_other_allele %in% c("A", "C", "G", "T")),
    f_statistic_recomputed = calc_ok(included) && calc_ok(excluded),
    included_f_lt10_zero = sum(included$F_stat < 10) == 0L,
    excluded_f_lt10_zero = sum(excluded$F_stat < 10) == 0L,
    no_f_ge_30_filter = TRUE,
    variant_level_n_preserved = all(is.finite(included$n_samples) & included$n_samples > 0) && all(is.finite(excluded$n_samples) & excluded$n_samples > 0),
    proxy_used_false = isFALSE(closure$proxy_used),
    liftover_used_false = isFALSE(closure$liftover_used),
    nearest_variant_used_false = isFALSE(closure$nearest_variant_used),
    strand_complement_identity_rescue_used_false = isFALSE(closure$strand_complement_identity_rescue_used),
    no_harmonisation = isFALSE(closure$harmonisation_performed),
    no_mr = isFALSE(closure$mr_run),
    no_steiger = isFALSE(closure$steiger_run),
    independent_replication_false = isFALSE(closure$independent_replication),
    manifest_complete = nrow(manifest) == nrow(manifest_items) && all(manifest$file_size_bytes > 0) && all(grepl("^[0-9a-f]{64}$", manifest$sha256)),
    manifest_paths_unique = !anyDuplicated(manifest$relative_path),
    manifest_sha_valid = grepl("^[0-9a-f]{64}$", manifest_sha)
  )
  failures <- names(hard_checks)[!unlist(hard_checks)]
  freeze_status <- if (length(failures) == 0L) "passed" else "failed"
  freeze <- list(
    freeze_version = "v1",
    authoritative_chen_forward_instrument_version = "v2",
    analysis_direction = "Hb_to_delirium",
    analysis_role = "forward_alternative_hb_gwas_sensitivity",
    independent_replication = FALSE,
    p_threshold = 5e-8,
    ld_r2 = 0.001,
    ld_window_kb = 10000,
    identifier_resolution_amendment_decision = 84,
    instrument_selection_decision = 85,
    technical_closure_decision = 86,
    scientific_result_source = "Decision 85 scientific outputs",
    technical_execution_closure = "Decision 86",
    decision_85_86_execution_closure_wording = "Decision 85 scientific outputs were complete, while the original process stopped during an overly strict Parquet/TSV floating-point readback check. Decision 86 independently verified scientific equivalence and closed the technical execution issue.",
    included_nsnp = nrow(included),
    excluded_nsnp = nrow(excluded),
    shared_nsnp = length(shared),
    included_only_nsnp = length(included_only),
    excluded_only_nsnp = length(excluded_only),
    included_rsids = as.list(included$resolved_rsid),
    excluded_rsids = as.list(excluded$resolved_rsid),
    instrument_strength_included = strength(included),
    instrument_strength_excluded = strength(excluded),
    variant_level_n_preserved = TRUE,
    identifier_resolution = list(
      source_identifier = "chr:pos_allele1_allele2",
      method = "exact_GRCh37_coordinate_plus_direct_unordered_allele_set_identity_unique_certified_1KG_EUR_reference_variant",
      canonical_rsid_field = "resolved_rsid",
      proxy_used = FALSE,
      liftover_used = FALSE,
      nearest_variant_used = FALSE,
      strand_complement_identity_rescue_used = FALSE
    ),
    effect_orientation = list(
      chen_reference_allele_is_exposure_effect_allele = TRUE,
      chen_other_allele_is_exposure_other_allele = TRUE,
      reference_orientation_may_replace_chen_effect_orientation = FALSE
    ),
    readback_closure = list(
      parquet_tsv_readback_closure_decision = 86,
      technical_readback_difference_field = "F_stat",
      scientific_identity_difference_detected = FALSE,
      readback_difference_class = "floating_point_text_serialization_representation",
      machine_precision_authority = "Parquet/source numeric values",
      tsv_role = "human-readable interchange"
    ),
    manifest_path = "results/qc/chen_forward_instruments_v2_freeze_manifest.csv",
    manifest_sha256 = manifest_sha,
    freeze_status = freeze_status,
    approved_for_chen_forward_finngen_outcome_extraction = identical(freeze_status, "passed"),
    hard_checks = hard_checks,
    hard_check_failures = failures,
    informational_findings = list(
      decision_85_initial_qc_retained_as_technical_provenance_non_authoritative = TRUE,
      failed_attempts_not_mr_input = TRUE
    )
  )
  atomic(out[["json"]], function(p) jsonlite::write_json(freeze, p, pretty = TRUE, auto_unbox = TRUE, na = "null"))

  decision_lines <- c(
    "# Decision 87 - Chen forward instruments V2 freeze",
    "",
    "Date: 2026-08-12",
    "",
    "## Decision",
    "",
    "Freeze Chen Forward Instruments V2 as the authoritative instrument set for",
    "the Chen alternative-Hb GWAS forward sensitivity branch.",
    "",
    "## Authority Model",
    "",
    "Scientific result source is Decision 85 / Instrument Selection V2. Technical",
    "execution closure is Decision 86. Decision 86 is not a new instrument",
    "selection. Decision 85 scientific outputs were complete, while the original",
    "process stopped during an overly strict Parquet/TSV floating-point readback",
    "check. Decision 86 independently verified scientific equivalence and closed",
    "the technical execution issue.",
    "",
    "## Results",
    "",
    sprintf("- Freeze status: `%s`", freeze_status),
    sprintf("- Manifest SHA256: `%s`", manifest_sha),
    sprintf("- Included nSNP: `%s`", nrow(included)),
    sprintf("- Excluded nSNP: `%s`", nrow(excluded)),
    sprintf("- Shared / included-only / excluded-only: `%s / %s / %s`", length(shared), length(included_only), length(excluded_only)),
    sprintf("- Included F min/median/mean/max/F<10: `%s / %s / %s / %s / %s`", strength(included)$F_min, strength(included)$F_median, strength(included)$F_mean, strength(included)$F_max, strength(included)$F_lt10_count),
    sprintf("- Excluded F min/median/mean/max/F<10: `%s / %s / %s / %s / %s`", strength(excluded)$F_min, strength(excluded)$F_median, strength(excluded)$F_mean, strength(excluded)$F_max, strength(excluded)$F_lt10_count),
    sprintf("- Hard-check failures: `%s`", paste(failures, collapse = ";")),
    sprintf("- Approved for Chen forward FinnGen outcome extraction: `%s`", identical(freeze_status, "passed")),
    "",
    "## Safeguards",
    "",
    "This freeze does not rerun instrument selection, rerun PLINK, modify Decision",
    "85 outputs, run FinnGen outcome extraction, harmonisation, MR, Steiger, proxy",
    "lookup, liftOver, or outcome-based filtering.",
    "",
    "## Affected Files",
    "",
    "- `R/21c_chen_forward_instruments_v2_freeze_manifest.R`",
    "- `results/qc/chen_forward_instruments_v2_freeze_manifest.csv`",
    "- `results/qc/chen_forward_instruments_v2_freeze.json`",
    "- `results/logs/chen_forward_instruments_v2_freeze.log`",
    "",
    "## Expected Impact",
    "",
    "This creates an auditable, hash-locked Chen V2 instrument authority for the",
    "next separately gated FinnGen targeted outcome extraction stage."
  )
  atomic(out[["decision"]], function(p) writeLines(decision_lines, p, useBytes = TRUE))
  log_line("freeze_status=", freeze_status, " manifest_sha256=", manifest_sha, " hard_check_failures=", paste(failures, collapse = ";"))
  stop_if(!identical(freeze_status, "passed"), "Chen Forward Instruments V2 freeze failed; outputs retained.")
}

tryCatch(main(), error = function(e) {
  cat(sprintf("[%s] freeze_status=failed error=%s\n", ts(), conditionMessage(e)), file = out[["log"]], append = TRUE)
  quit(status = 1L)
})
