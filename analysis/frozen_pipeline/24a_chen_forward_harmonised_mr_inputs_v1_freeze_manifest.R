options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/24a_chen_forward_harmonised_mr_inputs_v1_freeze_manifest.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

for (pkg in c("DBI", "duckdb", "jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing package: ", pkg, call. = FALSE)
  }
}

stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
sql_string <- function(path, must_work = TRUE) paste0("'", gsub("'", "''", normalizePath(path, winslash = "/", mustWork = must_work), fixed = TRUE), "'")
rel <- function(path) sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", root), "/?"), "", normalizePath(path, winslash = "/", mustWork = FALSE))
same_set <- function(a, b) identical(sort(unique(as.character(a))), sort(unique(as.character(b))))
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
atomic <- function(path, writer) {
  partial <- paste0(path, ".partial")
  stop_if(file.exists(path) || file.exists(partial), paste("Output occupied:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writer(partial)
  stop_if(!file.exists(partial), paste("Writer did not create partial:", partial))
  stop_if(!file.rename(partial, path), paste("Atomic rename failed:", path))
}
f_summary <- function(x) {
  list(
    n = length(x),
    F_min = if (length(x)) min(x, na.rm = TRUE) else NA_real_,
    F_mean = if (length(x)) mean(x, na.rm = TRUE) else NA_real_,
    F_median = if (length(x)) median(x, na.rm = TRUE) else NA_real_,
    F_max = if (length(x)) max(x, na.rm = TRUE) else NA_real_,
    F_lt10_count = sum(x < 10, na.rm = TRUE)
  )
}

paths <- c(
  script = file.path(root, "R", "24a_chen_forward_harmonised_mr_inputs_v1_freeze_manifest.R"),
  renv_lock = file.path(root, "renv.lock"),
  decision_78 = file.path(root, "docs", "decisions", "78_chen_2020_hb_source_certification_v1_v1.1.md"),
  decision_83 = file.path(root, "docs", "decisions", "83_chen_2020_hb_official_source_dictionary_audit_v1_v1.1.md"),
  decision_84 = file.path(root, "docs", "decisions", "84_chen_identifier_resolution_amendment_v1_v1.1.md"),
  decision_87 = file.path(root, "docs", "decisions", "87_chen_forward_instruments_v2_freeze_v1.1.md"),
  decision_91 = file.path(root, "docs", "decisions", "91_chen_forward_finngen_outcome_extraction_v2_freeze_v1.1.md"),
  decision_92 = file.path(root, "docs", "decisions", "92_chen_forward_harmonisation_contract_and_preflight_v1_v1.1.md"),
  decision_93 = file.path(root, "docs", "decisions", "93_chen_forward_formal_harmonisation_v1_v1.1.md"),
  formal_script = file.path(root, "R", "24_chen_forward_formal_harmonisation_v1.R"),
  formal_qc = file.path(root, "results", "qc", "chen_forward_formal_harmonisation_v1.json"),
  formal_counts = file.path(root, "results", "qc", "chen_forward_formal_harmonisation_counts_v1.csv"),
  transform_audit = file.path(root, "results", "qc", "chen_forward_formal_harmonisation_transform_audit_v1.csv"),
  excluded_snps = file.path(root, "results", "qc", "chen_forward_formal_harmonisation_excluded_snps_v1.tsv"),
  formal_log = file.path(root, "results", "logs", "chen_forward_formal_harmonisation_v1.log"),
  master_parquet = file.path(root, "data_derived", "forward_sensitivity_harmonisation", "chen_forward_harmonised_master_v1.parquet"),
  master_tsv = file.path(root, "data_derived", "forward_sensitivity_harmonisation", "chen_forward_harmonised_master_v1.tsv"),
  included_parquet = file.path(root, "data_derived", "forward_sensitivity_harmonisation", "chen_forward_harmonised_apoe_included_v1.parquet"),
  included_tsv = file.path(root, "data_derived", "forward_sensitivity_harmonisation", "chen_forward_harmonised_apoe_included_v1.tsv"),
  excluded_parquet = file.path(root, "data_derived", "forward_sensitivity_harmonisation", "chen_forward_harmonised_apoe_excluded_v1.parquet"),
  excluded_tsv = file.path(root, "data_derived", "forward_sensitivity_harmonisation", "chen_forward_harmonised_apoe_excluded_v1.tsv"),
  freeze_manifest = file.path(root, "results", "qc", "chen_forward_harmonised_mr_inputs_v1_freeze_manifest.csv"),
  freeze_json = file.path(root, "results", "qc", "chen_forward_harmonised_mr_inputs_v1_freeze.json"),
  freeze_log = file.path(root, "results", "logs", "chen_forward_harmonised_mr_inputs_v1_freeze.log"),
  freeze_decision = file.path(root, "docs", "decisions", "94_chen_forward_harmonised_mr_inputs_v1_freeze_v1.1.md")
)
for (p in paths[c("script", "renv_lock", "decision_78", "decision_83", "decision_84", "decision_87", "decision_91", "decision_92", "decision_93", "formal_script", "formal_qc", "formal_counts", "transform_audit", "excluded_snps", "formal_log", "master_parquet", "master_tsv", "included_parquet", "included_tsv", "excluded_parquet", "excluded_tsv")]) {
  stop_if(!file.exists(p), paste("Missing required input:", p))
}
for (p in paths[c("freeze_manifest", "freeze_json", "freeze_log", "freeze_decision")]) {
  stop_if(file.exists(p) || file.exists(paste0(p, ".partial")), paste("Output occupied:", p))
}

renv_before <- hash_file(paths[["renv_lock"]])
formal_qc <- jsonlite::fromJSON(paths[["formal_qc"]], simplifyVector = FALSE)

con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
read_pq <- function(path) DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(path)))
master <- read_pq(paths[["master_parquet"]])
included <- read_pq(paths[["included_parquet"]])
excluded <- read_pq(paths[["excluded_parquet"]])

included_ids <- as.character(included$resolved_rsid)
excluded_ids <- as.character(excluded$resolved_rsid)
shared <- intersect(included_ids, excluded_ids)
included_only <- setdiff(included_ids, excluded_ids)
excluded_only <- setdiff(excluded_ids, included_ids)

manifest_inputs <- data.frame(
  relative_path = rel(paths[c("decision_78", "decision_83", "decision_84", "decision_87", "decision_91", "decision_92", "decision_93", "formal_script", "formal_qc", "formal_counts", "transform_audit", "excluded_snps", "formal_log", "master_parquet", "master_tsv", "included_parquet", "included_tsv", "excluded_parquet", "excluded_tsv", "renv_lock")]),
  file_role = c(
    "chen_source_certification_decision",
    "chen_official_source_dictionary_decision",
    "chen_identifier_resolution_amendment_decision",
    "chen_forward_instrument_freeze_decision",
    "chen_forward_outcome_extraction_freeze_decision",
    "chen_forward_harmonisation_contract_preflight_decision",
    "chen_forward_formal_harmonisation_decision",
    "formal_harmonisation_script",
    "formal_harmonisation_qc",
    "formal_harmonisation_counts",
    "formal_harmonisation_transform_audit",
    "formal_harmonisation_excluded_snp_audit",
    "formal_harmonisation_log",
    "harmonised_master_parquet",
    "harmonised_master_tsv",
    "harmonised_apoe_included_parquet",
    "harmonised_apoe_included_tsv",
    "harmonised_apoe_excluded_parquet",
    "harmonised_apoe_excluded_tsv",
    "renv_lock"
  ),
  scientific_authority = TRUE,
  file_size_bytes = as.numeric(file.info(paths[c("decision_78", "decision_83", "decision_84", "decision_87", "decision_91", "decision_92", "decision_93", "formal_script", "formal_qc", "formal_counts", "transform_audit", "excluded_snps", "formal_log", "master_parquet", "master_tsv", "included_parquet", "included_tsv", "excluded_parquet", "excluded_tsv", "renv_lock")])$size),
  sha256 = vapply(paths[c("decision_78", "decision_83", "decision_84", "decision_87", "decision_91", "decision_92", "decision_93", "formal_script", "formal_qc", "formal_counts", "transform_audit", "excluded_snps", "formal_log", "master_parquet", "master_tsv", "included_parquet", "included_tsv", "excluded_parquet", "excluded_tsv", "renv_lock")], hash_file, character(1)),
  stringsAsFactors = FALSE
)
atomic(paths[["freeze_manifest"]], function(p) write.csv(manifest_inputs, p, row.names = FALSE, na = ""))
manifest_sha <- hash_file(paths[["freeze_manifest"]])

hard_checks <- list(
  formal_harmonisation_gate = identical(formal_qc$harmonisation_status, "passed") && isTRUE(formal_qc$approved_for_chen_forward_mr_input_freeze) && length(formal_qc$hard_check_failures) == 0L,
  formal_outputs_complete = all(file.exists(paths[c("master_parquet", "master_tsv", "included_parquet", "included_tsv", "excluded_parquet", "excluded_tsv", "formal_counts", "transform_audit", "excluded_snps", "formal_qc", "formal_log")])),
  final_counts_reverified = nrow(included) == formal_qc$final_valid_counts$included &&
    nrow(excluded) == formal_qc$final_valid_counts$excluded &&
    length(shared) == formal_qc$final_valid_counts$shared &&
    length(included_only) == formal_qc$final_valid_counts$included_only &&
    length(excluded_only) == formal_qc$final_valid_counts$excluded_only,
  final_rsids_reverified = same_set(included_ids, unlist(formal_qc$final_valid_rsid_sets$included)) &&
    same_set(excluded_ids, unlist(formal_qc$final_valid_rsid_sets$excluded)),
  no_palindromic_final_valid = !any(included$palindromic_snp) && !any(excluded$palindromic_snp),
  alleles_fully_aligned = all(included$outcome_effect_allele_harmonised == included$exposure_effect_allele) &&
    all(included$outcome_other_allele_harmonised == included$exposure_other_allele) &&
    all(excluded$outcome_effect_allele_harmonised == excluded$exposure_effect_allele) &&
    all(excluded$outcome_other_allele_harmonised == excluded$exposure_other_allele),
  all_se_positive = all(included$exposure_se > 0) && all(included$outcome_se_harmonised > 0) &&
    all(excluded$exposure_se > 0) && all(excluded$outcome_se_harmonised > 0),
  all_F_ge_10 = all(included$exposure_F_stat >= 10) && all(excluded$exposure_F_stat >= 10),
  no_posthoc_filtering = nrow(master) == formal_qc$final_valid_counts$union + formal_qc$palindromic_exclusion_counts$union +
    formal_qc$incompatible_exclusion_counts$union + formal_qc$invalid_exclusion_counts$union,
  analysis_role_sensitivity = identical(formal_qc$analysis_role, "forward_alternative_hb_gwas_sensitivity"),
  independent_replication_false = TRUE,
  software_provenance_complete = all(file.exists(paths[c("formal_script", "script", "renv_lock")])),
  renv_lock_unchanged = identical(renv_before, hash_file(paths[["renv_lock"]]))
)
failures <- names(hard_checks)[!unlist(hard_checks)]
status <- if (length(failures) == 0L) "passed" else "failed"
approved <- identical(status, "passed")

freeze <- list(
  freeze_version = "v1",
  authoritative_chen_forward_harmonisation_version = "v1",
  analysis_direction = "Hb_to_delirium",
  analysis_role = "forward_alternative_hb_gwas_sensitivity",
  independent_replication = FALSE,
  source_decisions = list(chen_certification = 78, chen_dictionary = 83, identifier_resolution = 84, instrument_freeze = 87, outcome_freeze = 91, harmonisation_contract_preflight = 92, formal_harmonisation = 93),
  included_final_valid_count = nrow(included),
  excluded_final_valid_count = nrow(excluded),
  shared_final_valid_count = length(shared),
  included_only_final_valid_count = length(included_only),
  excluded_only_final_valid_count = length(excluded_only),
  included_final_rsids = sort(included_ids),
  excluded_final_rsids = sort(excluded_ids),
  shared_final_rsids = sort(shared),
  included_only_final_rsids = sort(included_only),
  excluded_only_final_rsids = sort(excluded_only),
  instrument_strength_included = f_summary(included$exposure_F_stat),
  instrument_strength_excluded = f_summary(excluded$exposure_F_stat),
  palindromic_rule = formal_qc$palindromic_rule,
  effect_scale = "standardized_quantitative_Hb_effect",
  outcome_scale = "log_odds_delirium",
  manifest_path = "results/qc/chen_forward_harmonised_mr_inputs_v1_freeze_manifest.csv",
  manifest_sha256 = manifest_sha,
  freeze_status = status,
  approved_for_chen_forward_mr_design = approved,
  hard_checks = hard_checks,
  hard_check_failures = failures,
  manifest_records = records(manifest_inputs),
  informational_findings = list(
    freeze_did_not_reharmonise = TRUE,
    no_mr_run = TRUE,
    no_posthoc_filtering = TRUE,
    chen_analysis_is_sensitivity_not_replication = TRUE
  ),
  renv_lock_sha_before = renv_before,
  renv_lock_sha_after = hash_file(paths[["renv_lock"]])
)
atomic(paths[["freeze_json"]], function(p) jsonlite::write_json(freeze, p, pretty = TRUE, auto_unbox = TRUE, na = "null"))

decision_lines <- c(
  "# Decision 94 - Chen forward harmonised MR inputs V1 freeze",
  "",
  "Date: 2026-08-12",
  "Status: harmonised MR-input freeze",
  "",
  "## Decision",
  "",
  "Freeze the Decision 93 Chen forward harmonised final-valid APOE included and",
  "excluded datasets as the authoritative MR-input candidates for a future Chen",
  "forward MR design stage.",
  "",
  "## Results",
  "",
  sprintf("- freeze_status: `%s`", status),
  sprintf("- manifest SHA-256: `%s`", manifest_sha),
  sprintf("- included/excluded final-valid: `%d/%d`", nrow(included), nrow(excluded)),
  sprintf("- shared/included-only/excluded-only: `%d/%d/%d`", length(shared), length(included_only), length(excluded_only)),
  sprintf("- hard_check_failures: `%s`", paste(failures, collapse = ";")),
  sprintf("- approved_for_chen_forward_mr_design: `%s`", approved),
  "",
  "## Safeguards",
  "",
  "This freeze did not re-harmonise, run MR, run Steiger, reselect instruments,",
  "re-clump, use proxy lookup, use liftOver, use outcome-based filtering, or",
  "apply an F>=30 filter.",
  "",
  "## Expected Impact",
  "",
  "If passed, this freeze authorizes only a future separately approved Chen",
  "forward MR analysis contract. It does not authorize MR."
)
atomic(paths[["freeze_decision"]], function(p) writeLines(decision_lines, p, useBytes = TRUE))

atomic(paths[["freeze_log"]], function(p) {
  writeLines(c(
    sprintf("[%s] stage=chen_forward_harmonised_mr_inputs_v1_freeze", ts()),
    sprintf("[%s] freeze_status=%s approved_for_mr_design=%s", ts(), status, approved),
    sprintf("[%s] included_final=%d excluded_final=%d shared=%d included_only=%d excluded_only=%d", ts(), nrow(included), nrow(excluded), length(shared), length(included_only), length(excluded_only)),
    sprintf("[%s] hard_check_failures=%s", ts(), paste(failures, collapse = ";"))
  ), p, useBytes = TRUE)
})

stop_if(!identical(status, "passed"), "MR input freeze failed; QC retained.")
message("Chen forward harmonised MR input freeze completed: ", status)
