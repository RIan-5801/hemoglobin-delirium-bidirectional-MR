options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/22a_chen_forward_finngen_outcome_extraction_v2_freeze_manifest.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

for (pkg in c("jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing package: ", pkg, call. = FALSE)
  }
}

stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
rel <- function(path) {
  sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", root), "/?"), "", normalizePath(path, winslash = "/", mustWork = FALSE))
}
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
atomic <- function(path, writer) {
  partial <- paste0(path, ".partial")
  stop_if(file.exists(path) || file.exists(partial), paste("Output occupied:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writer(partial)
  stop_if(!file.exists(partial), paste("Writer did not create partial:", partial))
  stop_if(!file.rename(partial, path), paste("Atomic rename failed:", path))
}
empty_chr <- function(x) if (length(x) == 0L) character() else as.character(x)
same_set <- function(a, b) identical(sort(unique(as.character(a))), sort(unique(as.character(b))))
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))

paths <- c(
  script = file.path(root, "R", "22a_chen_forward_finngen_outcome_extraction_v2_freeze_manifest.R"),
  renv_lock = file.path(root, "renv.lock"),
  finngen_source = file.path(root, "data_raw", "gwas", "finngen_R13_F5_DELIRIUM.gz"),
  metadata = file.path(root, "docs", "02_gwas_metadata_v2.md"),
  decision_87 = file.path(root, "docs", "decisions", "87_chen_forward_instruments_v2_freeze_v1.1.md"),
  decision_89 = file.path(root, "docs", "decisions", "89_chen_forward_finngen_outcome_extraction_v2_technical_recovery_v1.1.md"),
  decision_90 = file.path(root, "docs", "decisions", "90_chen_forward_finngen_outcome_extraction_v2_readback_closure_v1.1.md"),
  v1_failed_log = file.path(root, "results", "logs", "chen_forward_finngen_outcome_extraction_v1.log"),
  instrument_freeze_qc = file.path(root, "results", "qc", "chen_forward_instruments_v2_freeze.json"),
  outcome_v2_script = file.path(root, "R", "22b_chen_forward_finngen_outcome_extraction_v2_technical_recovery.R"),
  outcome_closure_script = file.path(root, "R", "22c_chen_forward_finngen_outcome_extraction_v2_readback_closure.R"),
  outcome_v2_qc = file.path(root, "results", "qc", "chen_forward_finngen_outcome_extraction_v2.json"),
  outcome_closure_qc = file.path(root, "results", "qc", "chen_forward_finngen_outcome_extraction_v2_readback_closure_v1.json"),
  outcome_v2_log = file.path(root, "results", "logs", "chen_forward_finngen_outcome_extraction_v2.log"),
  outcome_closure_log = file.path(root, "results", "logs", "chen_forward_finngen_outcome_extraction_v2_readback_closure_v1.log"),
  union_targets = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_union_targets_v2.tsv"),
  match_audit = file.path(root, "results", "qc", "chen_forward_finngen_outcome_match_audit_v2.csv"),
  missing = file.path(root, "results", "qc", "chen_forward_finngen_outcome_missing_v2.tsv"),
  master_parquet = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_master_v2.parquet"),
  master_tsv = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_master_v2.tsv"),
  included_parquet = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_apoe_included_v2.parquet"),
  included_tsv = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_apoe_included_v2.tsv"),
  excluded_parquet = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_apoe_excluded_v2.parquet"),
  excluded_tsv = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_apoe_excluded_v2.tsv"),
  included_instruments = file.path(root, "data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_included_clumped_v2.tsv"),
  excluded_instruments = file.path(root, "data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_excluded_clumped_v2.tsv"),
  freeze_manifest = file.path(root, "results", "qc", "chen_forward_finngen_outcome_extraction_v2_freeze_manifest.csv"),
  freeze_json = file.path(root, "results", "qc", "chen_forward_finngen_outcome_extraction_v2_freeze.json"),
  freeze_log = file.path(root, "results", "logs", "chen_forward_finngen_outcome_extraction_v2_freeze.log"),
  freeze_decision = file.path(root, "docs", "decisions", "91_chen_forward_finngen_outcome_extraction_v2_freeze_v1.1.md")
)

for (p in paths[c("script", "renv_lock", "finngen_source", "metadata", "decision_87", "decision_89", "decision_90", "v1_failed_log", "instrument_freeze_qc", "outcome_v2_script", "outcome_closure_script", "outcome_v2_qc", "outcome_closure_qc", "outcome_v2_log", "outcome_closure_log", "union_targets", "match_audit", "missing", "master_parquet", "master_tsv", "included_parquet", "included_tsv", "excluded_parquet", "excluded_tsv", "included_instruments", "excluded_instruments")]) {
  stop_if(!file.exists(p), paste("Missing required input:", p))
}
for (p in paths[c("freeze_manifest", "freeze_json", "freeze_log", "freeze_decision")]) {
  stop_if(file.exists(p) || file.exists(paste0(p, ".partial")), paste("Output occupied:", p))
}

renv_before <- hash_file(paths[["renv_lock"]])
source_sha <- tolower(hash_file(paths[["finngen_source"]]))
expected_source_sha <- "85637f0f3358807964d4f8a3e500293168a706f1c08c65f3fc5512b65df40ed8"

instrument_freeze <- jsonlite::fromJSON(paths[["instrument_freeze_qc"]], simplifyVector = FALSE)
outcome_closure <- jsonlite::fromJSON(paths[["outcome_closure_qc"]], simplifyVector = FALSE)
included_inst <- read.delim(paths[["included_instruments"]], sep = "\t", check.names = FALSE)
excluded_inst <- read.delim(paths[["excluded_instruments"]], sep = "\t", check.names = FALSE)
targets <- read.delim(paths[["union_targets"]], sep = "\t", check.names = FALSE)
master <- read.delim(paths[["master_tsv"]], sep = "\t", check.names = FALSE, na.strings = c(""))
included_out <- read.delim(paths[["included_tsv"]], sep = "\t", check.names = FALSE, na.strings = c(""))
excluded_out <- read.delim(paths[["excluded_tsv"]], sep = "\t", check.names = FALSE, na.strings = c(""))
match_audit <- read.csv(paths[["match_audit"]], check.names = FALSE)
missing <- read.delim(paths[["missing"]], sep = "\t", check.names = FALSE)

included_ids <- unique(as.character(included_inst$resolved_rsid))
excluded_ids <- unique(as.character(excluded_inst$resolved_rsid))
union_ids <- sort(unique(c(included_ids, excluded_ids)))
shared_ids <- intersect(included_ids, excluded_ids)
included_only_ids <- setdiff(included_ids, excluded_ids)
excluded_only_ids <- setdiff(excluded_ids, included_ids)

target_included_ids <- as.character(targets$resolved_rsid[targets$in_apoe_included_input])
target_excluded_ids <- as.character(targets$resolved_rsid[targets$in_apoe_excluded_input])
target_union_ids <- as.character(targets$resolved_rsid)

exact_ids <- as.character(master$resolved_rsid[master$outcome_match_status == "unique_exact_match"])
missing_ids <- as.character(missing$resolved_rsid)
included_exact_ids <- as.character(included_out$resolved_rsid[included_out$outcome_match_status == "unique_exact_match"])
excluded_exact_ids <- as.character(excluded_out$resolved_rsid[excluded_out$outcome_match_status == "unique_exact_match"])

matched_shared_ids <- intersect(exact_ids, shared_ids)
matched_included_only_ids <- intersect(exact_ids, included_only_ids)
matched_excluded_only_ids <- intersect(exact_ids, excluded_only_ids)

manifest_inputs <- data.frame(
  relative_path = rel(paths[c("metadata", "decision_87", "decision_89", "decision_90", "v1_failed_log", "instrument_freeze_qc", "outcome_v2_script", "outcome_closure_script", "outcome_v2_qc", "outcome_closure_qc", "outcome_v2_log", "outcome_closure_log", "union_targets", "match_audit", "missing", "master_parquet", "master_tsv", "included_parquet", "included_tsv", "excluded_parquet", "excluded_tsv", "included_instruments", "excluded_instruments", "renv_lock")]),
  file_role = c(
    "authenticated_gwas_metadata",
    "source_instrument_freeze_decision",
    "technical_recovery_decision",
    "technical_readback_closure_decision",
    "v1_failed_attempt_log_non_authoritative",
    "source_instrument_freeze_qc",
    "outcome_extraction_v2_script",
    "outcome_readback_closure_script",
    "outcome_extraction_v2_qc_technical_provenance",
    "outcome_extraction_v2_readback_closure_authoritative_qc",
    "outcome_extraction_v2_log",
    "outcome_readback_closure_log",
    "outcome_union_targets_v2",
    "outcome_match_audit_v2",
    "outcome_missing_list_v2",
    "outcome_master_parquet_authority",
    "outcome_master_tsv_human_readable",
    "outcome_apoe_included_parquet_authority",
    "outcome_apoe_included_tsv_human_readable",
    "outcome_apoe_excluded_parquet_authority",
    "outcome_apoe_excluded_tsv_human_readable",
    "decision_87_included_instrument_tsv",
    "decision_87_excluded_instrument_tsv",
    "renv_lock"
  ),
  scientific_authority = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  file_size_bytes = as.numeric(file.info(paths[c("metadata", "decision_87", "decision_89", "decision_90", "v1_failed_log", "instrument_freeze_qc", "outcome_v2_script", "outcome_closure_script", "outcome_v2_qc", "outcome_closure_qc", "outcome_v2_log", "outcome_closure_log", "union_targets", "match_audit", "missing", "master_parquet", "master_tsv", "included_parquet", "included_tsv", "excluded_parquet", "excluded_tsv", "included_instruments", "excluded_instruments", "renv_lock")])$size),
  sha256 = vapply(paths[c("metadata", "decision_87", "decision_89", "decision_90", "v1_failed_log", "instrument_freeze_qc", "outcome_v2_script", "outcome_closure_script", "outcome_v2_qc", "outcome_closure_qc", "outcome_v2_log", "outcome_closure_log", "union_targets", "match_audit", "missing", "master_parquet", "master_tsv", "included_parquet", "included_tsv", "excluded_parquet", "excluded_tsv", "included_instruments", "excluded_instruments", "renv_lock")], hash_file, character(1)),
  stringsAsFactors = FALSE
)

hard_checks <- list(
  decision_87_freeze_passed = identical(instrument_freeze$freeze_status, "passed"),
  decision_87_approved_for_outcome_extraction = isTRUE(instrument_freeze$approved_for_chen_forward_finngen_outcome_extraction),
  decision_89_present = file.exists(paths[["decision_89"]]),
  decision_90_closure_passed = identical(outcome_closure$outcome_extraction_status, "passed"),
  decision_90_approved_for_harmonisation_design = isTRUE(outcome_closure$approved_for_chen_forward_harmonisation_design),
  decision_90_hard_failures_empty = length(outcome_closure$hard_check_failures) == 0L,
  v1_failure_evidence_preserved_non_authoritative = file.exists(paths[["v1_failed_log"]]),
  source_sha_matches_certification = identical(source_sha, expected_source_sha),
  source_sha_matches_decision_90_after = identical(source_sha, outcome_closure$source_sha_after_readback_closure),
  included_target_set_exact = same_set(included_ids, target_included_ids),
  excluded_target_set_exact = same_set(excluded_ids, target_excluded_ids),
  union_target_set_exact = same_set(union_ids, target_union_ids),
  included_target_count_recomputed = length(included_ids) == outcome_closure$included_target_count,
  excluded_target_count_recomputed = length(excluded_ids) == outcome_closure$excluded_target_count,
  shared_target_count_recomputed = length(shared_ids) == outcome_closure$shared_target_count,
  union_target_count_recomputed = length(union_ids) == outcome_closure$union_target_count,
  union_exact_match_count_recomputed = length(exact_ids) == outcome_closure$union_unique_exact_match_count,
  union_missing_count_recomputed = length(missing_ids) == outcome_closure$union_missing_count,
  included_exact_match_count_recomputed = length(included_exact_ids) == outcome_closure$included_exact_match_count,
  included_missing_count_recomputed = length(setdiff(included_ids, included_exact_ids)) == outcome_closure$included_missing_count,
  excluded_exact_match_count_recomputed = length(excluded_exact_ids) == outcome_closure$excluded_exact_match_count,
  excluded_missing_count_recomputed = length(setdiff(excluded_ids, excluded_exact_ids)) == outcome_closure$excluded_missing_count,
  no_unresolved_multiple_matches = identical(outcome_closure$union_multiple_match_count, 0L),
  missing_classification_complete = nrow(missing) == length(missing_ids) && all(missing$finngen_match_status == "missing"),
  exact_rsid_matching_only = identical(outcome_closure$matching_method, "exact_canonical_rsid_token"),
  no_proxy = !isTRUE(outcome_closure$proxy_used),
  no_liftover = !isTRUE(outcome_closure$liftover_used),
  no_coordinate_matching = !isTRUE(outcome_closure$coordinate_matching_used),
  cross_build_coordinate_matching_used_false = !isTRUE(outcome_closure$cross_build_coordinate_comparison_used),
  no_harmonisation = !isTRUE(outcome_closure$harmonisation_performed),
  no_mr = !isTRUE(outcome_closure$mr_run),
  no_steiger = !isTRUE(outcome_closure$steiger_run),
  manifest_complete = nrow(manifest_inputs) == 24L,
  renv_lock_unchanged = identical(renv_before, hash_file(paths[["renv_lock"]]))
)
failures <- names(hard_checks)[!unlist(hard_checks)]
status <- if (length(failures) == 0L) "passed" else "failed"
approved <- identical(status, "passed")

atomic(paths[["freeze_manifest"]], function(p) write.csv(manifest_inputs, p, row.names = FALSE, na = ""))
manifest_sha <- hash_file(paths[["freeze_manifest"]])

freeze <- list(
  freeze_version = "v1",
  authoritative_chen_forward_outcome_extraction_version = "v2",
  analysis_direction = "Hb_to_delirium",
  analysis_role = "forward_alternative_hb_gwas_sensitivity",
  independent_replication = FALSE,
  source_instrument_freeze_decision = 87,
  technical_recovery_decision = 89,
  technical_readback_closure_decision = 90,
  scientific_result_source = "FinnGen Outcome Extraction V2",
  technical_recovery_provenance = "Decision 89",
  technical_readback_closure = "Decision 90",
  v1_failed_attempt_role = "technical_provenance_non_authoritative",
  outcome_source = "FinnGen_R13_F5_DELIRIUM",
  source_path = "data_raw/gwas/finngen_R13_F5_DELIRIUM.gz",
  source_sha256 = source_sha,
  source_sha_matches_certification = identical(source_sha, expected_source_sha),
  outcome_build = "GRCh38",
  outcome_metadata = outcome_closure$outcome_metadata,
  chen_exposure_build = "GRCh37",
  cross_build_coordinate_matching_used = FALSE,
  cross_build_position_equality_required = FALSE,
  matching_method = "exact_canonical_rsid",
  included_target_count = length(included_ids),
  excluded_target_count = length(excluded_ids),
  shared_target_count = length(shared_ids),
  union_target_count = length(union_ids),
  union_exact_match_count = length(exact_ids),
  union_missing_count = length(missing_ids),
  union_multiple_match_count = outcome_closure$union_multiple_match_count,
  included_exact_match_count = length(included_exact_ids),
  included_missing_count = length(setdiff(included_ids, included_exact_ids)),
  excluded_exact_match_count = length(excluded_exact_ids),
  excluded_missing_count = length(setdiff(excluded_ids, excluded_exact_ids)),
  matched_shared_count = length(matched_shared_ids),
  matched_included_only_count = length(matched_included_only_ids),
  matched_excluded_only_count = length(matched_excluded_only_ids),
  missing_semantics = "exposure_instrument_not_available_by_exact_rsid_in_FinnGen_outcome",
  missing_rsids = sort(missing_ids),
  matched_shared_rsids = sort(matched_shared_ids),
  matched_included_only_rsids = sort(matched_included_only_ids),
  matched_excluded_only_rsids = sort(matched_excluded_only_ids),
  proxy_used = FALSE,
  liftover_used = FALSE,
  coordinate_matching_used = FALSE,
  nearest_variant_used = FALSE,
  outcome_p_filtering_used = FALSE,
  readback_semantics = list(
    parquet_authority = "machine-precision numeric authority",
    tsv_role = "human-readable interchange",
    numeric_tolerance = "accepted absolute + relative tolerance from Decision 90",
    character_blank_na_equivalence_closed_by_decision_90 = TRUE
  ),
  manifest_path = "results/qc/chen_forward_finngen_outcome_extraction_v2_freeze_manifest.csv",
  manifest_sha256 = manifest_sha,
  freeze_status = status,
  approved_for_chen_forward_harmonisation_preflight = approved,
  hard_checks = hard_checks,
  hard_check_failures = failures,
  informational_findings = list(
    v1_failed_attempt_retained_non_authoritative = TRUE,
    missing_outcome_targets_allowed_without_proxy = TRUE,
    missing_outcomes_not_harmonisation_exclusions = TRUE,
    no_outcome_based_filtering_performed = TRUE
  ),
  manifest_records = records(manifest_inputs),
  renv_lock_sha_before = renv_before,
  renv_lock_sha_after = hash_file(paths[["renv_lock"]])
)

atomic(paths[["freeze_json"]], function(p) jsonlite::write_json(freeze, p, pretty = TRUE, auto_unbox = TRUE, na = "null"))

decision_lines <- c(
  "# Decision 91 - Chen forward FinnGen outcome extraction V2 freeze",
  "",
  "Date: 2026-08-12",
  "Status: outcome extraction freeze",
  "",
  "## Decision",
  "",
  "Freeze FinnGen Outcome Extraction V2 as the authoritative outcome extraction",
  "set for the Chen alternative-Hb-GWAS forward sensitivity analysis.",
  "",
  "## Authority Model",
  "",
  "- scientific_result_source: `FinnGen Outcome Extraction V2`",
  "- technical_recovery_provenance: `Decision 89`",
  "- technical_readback_closure: `Decision 90`",
  "- V1 failed attempt role: `technical_provenance_non_authoritative`",
  "",
  "## Results",
  "",
  sprintf("- freeze_status: `%s`", status),
  sprintf("- approved_for_chen_forward_harmonisation_preflight: `%s`", approved),
  sprintf("- manifest SHA-256: `%s`", manifest_sha),
  sprintf("- included target/exact/missing: `%d/%d/%d`", length(included_ids), length(included_exact_ids), length(setdiff(included_ids, included_exact_ids))),
  sprintf("- excluded target/exact/missing: `%d/%d/%d`", length(excluded_ids), length(excluded_exact_ids), length(setdiff(excluded_ids, excluded_exact_ids))),
  sprintf("- union target/exact/missing/multiple: `%d/%d/%d/%d`", length(union_ids), length(exact_ids), length(missing_ids), outcome_closure$union_multiple_match_count),
  sprintf("- matched shared/included-only/excluded-only: `%d/%d/%d`", length(matched_shared_ids), length(matched_included_only_ids), length(matched_excluded_only_ids)),
  sprintf("- FinnGen source SHA-256: `%s`", source_sha),
  sprintf("- hard_check_failures: `%s`", paste(failures, collapse = ";")),
  "",
  "## Missing Semantics",
  "",
  "Missing instruments are classified as",
  "`exposure_instrument_not_available_by_exact_rsid_in_FinnGen_outcome`.",
  "They are not harmonisation exclusions because formal harmonisation has not yet begun.",
  "",
  "## Safeguards",
  "",
  "No proxy lookup, liftOver, coordinate matching, nearest-variant matching,",
  "outcome-P filtering, harmonisation, MR, or Steiger analysis was performed.",
  "",
  "## Expected Impact",
  "",
  "If passed, this freeze authorizes only the next separately gated Chen forward",
  "harmonisation contract and preflight stage. It does not authorize formal",
  "harmonisation or MR."
)
atomic(paths[["freeze_decision"]], function(p) writeLines(decision_lines, p, useBytes = TRUE))

atomic(paths[["freeze_log"]], function(p) {
  writeLines(c(
    sprintf("[%s] stage=chen_forward_finngen_outcome_extraction_v2_freeze", ts()),
    sprintf("[%s] freeze_status=%s approved_for_preflight=%s", ts(), status, approved),
    sprintf("[%s] union_targets=%d union_exact=%d union_missing=%d multiple=%d", ts(), length(union_ids), length(exact_ids), length(missing_ids), outcome_closure$union_multiple_match_count),
    sprintf("[%s] hard_check_failures=%s", ts(), paste(failures, collapse = ";"))
  ), p, useBytes = TRUE)
})

stop_if(!identical(status, "passed"), "Outcome extraction freeze failed; QC retained.")
message("Outcome extraction V2 freeze completed: ", status)
