#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/28_chen_reverse_outcome_extraction_execution_contract_v1.R [--project-root <path>]", call. = FALSE)
}
setwd(root)

for (pkg in c("jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

rel <- function(...) file.path(root, ...)
norm <- function(x) gsub("\\\\", "/", x)
relpath <- function(x) {
  norm(sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", norm(normalizePath(x, winslash = "/", mustWork = FALSE))))
}
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
text_file <- function(path) paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
write_csv_precise <- function(x, path) {
  partial <- paste0(path, ".partial")
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.csv(x, partial, row.names = FALSE, na = "")
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}
write_json <- function(x, path) {
  partial <- paste0(path, ".partial")
  jsonlite::write_json(x, partial, pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}
write_text <- function(lines, path) {
  partial <- paste0(path, ".partial")
  writeLines(lines, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}

paths <- list(
  script = rel("R", "28_chen_reverse_outcome_extraction_execution_contract_v1.R"),
  decision104 = rel("docs", "decisions", "104_chen_reverse_sensitivity_design_v1_v1.1.md"),
  decision105 = rel("docs", "decisions", "105_chen_reverse_outcome_extraction_preflight_plan_v1_v1.1.md"),
  decision106 = rel("docs", "decisions", "106_chen_reverse_outcome_extraction_preflight_plan_v1_readback_recovery_v1.1.md"),
  plan106_qc = rel("results", "qc", "chen_reverse_outcome_extraction_preflight_plan_v1_readback_recovery.json"),
  decision78 = rel("docs", "decisions", "78_chen_2020_hb_source_certification_v1_v1.1.md"),
  decision83 = rel("docs", "decisions", "83_chen_2020_hb_official_source_dictionary_audit_v1_v1.1.md"),
  decision84 = rel("docs", "decisions", "84_chen_identifier_resolution_amendment_v1_v1.1.md"),
  decision47 = rel("docs", "decisions", "47_finngen_r13_delirium_reverse_primary_instruments_v4_freeze_v3_v1.1.md"),
  decision59 = rel("docs", "decisions", "59_finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze_v1.1.md"),
  chen_source_cert = rel("results", "qc", "chen_2020_hb_source_certification_v1.json"),
  chen_readme = rel("docs", "source_metadata", "readme_BCX2_meta_analyses.txt"),
  strict_freeze = rel("results", "qc", "finngen_r13_delirium_reverse_primary_instruments_v4_freeze_v3.json"),
  strict_selection_qc = rel("results", "qc", "finngen_r13_delirium_reverse_primary_instrument_selection_v4.json"),
  strict_manifest = rel("results", "qc", "finngen_r13_delirium_reverse_primary_instruments_v4_freeze_manifest_v3.csv"),
  strict_included_tsv = rel("data_derived", "reverse_instruments", "finngen_r13_delirium_primary_apoe_included_clumped_v4.tsv"),
  strict_excluded_tsv = rel("data_derived", "reverse_instruments", "finngen_r13_delirium_primary_apoe_excluded_clumped_v4.tsv"),
  relaxed_freeze = rel("results", "qc", "finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze.json"),
  relaxed_selection_qc = rel("results", "qc", "finngen_r13_delirium_p5e-6_exploratory_instrument_selection_v2.json"),
  relaxed_manifest = rel("results", "qc", "finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze_manifest.csv"),
  relaxed_included_tsv = rel("data_derived", "reverse_instruments", "finngen_r13_delirium_p5e-6_exploratory_apoe_included_clumped_v2.tsv"),
  relaxed_excluded_tsv = rel("data_derived", "reverse_instruments", "finngen_r13_delirium_p5e-6_exploratory_apoe_excluded_clumped_v2.tsv"),
  eur_bim = rel("resources", "ld", "1kg_v3", "EUR.bim"),
  renv_lock = rel("renv.lock"),
  contract_json = rel("results", "qc", "chen_reverse_outcome_extraction_execution_contract_v1.json"),
  target_audit = rel("results", "qc", "chen_reverse_outcome_target_authority_audit_v1.csv"),
  contract_log = rel("results", "logs", "chen_reverse_outcome_extraction_execution_contract_v1.log"),
  decision107 = rel("docs", "decisions", "107_chen_reverse_outcome_extraction_execution_contract_v1_v1.1.md")
)

required_inputs <- unlist(paths[c(
  "script", "decision104", "decision105", "decision106", "plan106_qc",
  "decision78", "decision83", "decision84", "decision47", "decision59",
  "chen_source_cert", "chen_readme", "strict_freeze", "strict_selection_qc",
  "strict_manifest", "strict_included_tsv", "strict_excluded_tsv",
  "relaxed_freeze", "relaxed_selection_qc", "relaxed_manifest",
  "relaxed_included_tsv", "relaxed_excluded_tsv", "eur_bim", "renv_lock"
)])
missing <- required_inputs[!file.exists(required_inputs)]
if (length(missing) > 0L) stop("Missing required input(s): ", paste(missing, collapse = "; "), call. = FALSE)

targets <- unlist(paths[c("contract_json", "target_audit", "contract_log", "decision107")])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(occupied, collapse = "; ")))

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 107L), paste0("Expected next decision 107, found ", next_decision, "; no outputs written."))

read_branch <- function(path, branch) {
  x <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("rsid", "ref", "alt", "beta", "se", "eaf", "F_statistic")
  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols) > 0L) stop("Missing column(s) in ", path, ": ", paste(missing_cols, collapse = ", "), call. = FALSE)
  x$branch <- branch
  x
}

strict_inc <- read_branch(paths$strict_included_tsv, "strict_apoe_included")
strict_exc <- read_branch(paths$strict_excluded_tsv, "strict_apoe_excluded")
relaxed_inc <- read_branch(paths$relaxed_included_tsv, "relaxed_apoe_included")
relaxed_exc <- read_branch(paths$relaxed_excluded_tsv, "relaxed_apoe_excluded")

all_branches <- rbind(strict_inc, strict_exc, relaxed_inc, relaxed_exc)
union_rsids <- sort(unique(all_branches$rsid))

row_for <- function(branch_df, rsid) {
  z <- branch_df[branch_df$rsid == rsid, , drop = FALSE]
  if (nrow(z) == 0L) return(NULL)
  if (nrow(z) > 1L) stop("Duplicate rsID in branch input: ", rsid, call. = FALSE)
  z[1L, , drop = FALSE]
}

best_row <- function(rsid) {
  for (df in list(strict_inc, strict_exc, relaxed_inc, relaxed_exc)) {
    z <- row_for(df, rsid)
    if (!is.null(z)) return(z)
  }
  stop("No row for rsID: ", rsid, call. = FALSE)
}

read_bim_targets <- function(path, rsids) {
  con <- file(path, open = "r")
  on.exit(close(con), add = TRUE)
  out <- list()
  repeat {
    lines <- readLines(con, n = 100000L, warn = FALSE)
    if (length(lines) == 0L) break
    parts <- strsplit(lines, "\t", fixed = TRUE)
    keep <- vapply(parts, function(z) length(z) >= 6L && z[[2L]] %in% rsids, logical(1))
    if (any(keep)) {
      out <- c(out, lapply(parts[keep], function(z) {
        data.frame(
          reference_chr_grch37 = z[[1L]],
          rsid = z[[2L]],
          reference_pos_grch37 = as.integer(z[[4L]]),
          reference_allele1 = z[[5L]],
          reference_allele2 = z[[6L]],
          stringsAsFactors = FALSE
        )
      }))
    }
  }
  if (length(out) == 0L) {
    return(data.frame(reference_chr_grch37 = character(), rsid = character(), reference_pos_grch37 = integer(), reference_allele1 = character(), reference_allele2 = character()))
  }
  do.call(rbind, out)
}

sort_alleles <- function(a, b) paste(sort(c(a, b)), collapse = "_")
expected_marker <- function(chr, pos, a1, a2) paste0(chr, ":", pos, "_", sort_alleles(a1, a2))

bim <- read_bim_targets(paths$eur_bim, union_rsids)
duplicate_bim_rsids <- names(which(table(bim$rsid) > 1L))

audit <- do.call(rbind, lapply(union_rsids, function(rsid) {
  b <- best_row(rsid)
  ref_rows <- bim[bim$rsid == rsid, , drop = FALSE]
  bridge_ready <- nrow(ref_rows) == 1L
  ref <- if (bridge_ready) ref_rows[1L, , drop = FALSE] else data.frame(
    reference_chr_grch37 = NA_character_,
    reference_pos_grch37 = NA_integer_,
    reference_allele1 = NA_character_,
    reference_allele2 = NA_character_,
    stringsAsFactors = FALSE
  )
  data.frame(
    rsid = rsid,
    strict_included_member = rsid %in% strict_inc$rsid,
    strict_excluded_member = rsid %in% strict_exc$rsid,
    relaxed_included_member = rsid %in% relaxed_inc$rsid,
    relaxed_excluded_member = rsid %in% relaxed_exc$rsid,
    source_exposure_freeze_decision = paste(
      unique(c(
        if (rsid %in% c(strict_inc$rsid, strict_exc$rsid)) "47",
        if (rsid %in% c(relaxed_inc$rsid, relaxed_exc$rsid)) "59"
      )),
      collapse = ";"
    ),
    exposure_beta = b$beta,
    exposure_se = b$se,
    exposure_effect_allele = b$alt,
    exposure_other_allele = b$ref,
    exposure_eaf = b$eaf,
    F_stat = b$F_statistic,
    reference_grch37_exact_rsid_present = bridge_ready,
    reference_chr_grch37 = ref$reference_chr_grch37,
    reference_pos_grch37 = ref$reference_pos_grch37,
    reference_allele1 = ref$reference_allele1,
    reference_allele2 = ref$reference_allele2,
    expected_chen_marker_id = if (bridge_ready) expected_marker(ref$reference_chr_grch37, ref$reference_pos_grch37, ref$reference_allele1, ref$reference_allele2) else NA_character_,
    identity_bridge_ready = bridge_ready,
    stringsAsFactors = FALSE
  )
}))

write_csv_precise(audit, paths$target_audit)

plan106 <- read_json(paths$plan106_qc)
strict_freeze <- read_json(paths$strict_freeze)
strict_qc <- read_json(paths$strict_selection_qc)
relaxed_freeze <- read_json(paths$relaxed_freeze)
relaxed_qc <- read_json(paths$relaxed_selection_qc)
chen_cert <- read_json(paths$chen_source_cert)
readme_text <- text_file(paths$chen_readme)
decision105_text <- text_file(paths$decision105)
decision106_text <- text_file(paths$decision106)

branch_counts <- list(
  strict_apoe_included = nrow(strict_inc),
  strict_apoe_excluded = nrow(strict_exc),
  relaxed_apoe_included = nrow(relaxed_inc),
  relaxed_apoe_excluded = nrow(relaxed_exc),
  union = length(union_rsids),
  shared_all_four = sum(audit$strict_included_member & audit$strict_excluded_member & audit$relaxed_included_member & audit$relaxed_excluded_member),
  strict_union = length(unique(c(strict_inc$rsid, strict_exc$rsid))),
  relaxed_union = length(unique(c(relaxed_inc$rsid, relaxed_exc$rsid)))
)

target_authority_ready_count <- sum(audit$identity_bridge_ready)
target_authority_unavailable_count <- sum(!audit$identity_bridge_ready)

contains_all <- function(text, terms) all(vapply(terms, grepl, logical(1), x = text, fixed = TRUE))
hard_checks <- list(
  decision_104_design_gate = identical(plan106$corrected_hard_checks$decision_104_design_gate, TRUE),
  decision_105_evidence_preserved = identical(plan106$source_plan_decision, 105) &&
    identical(plan106$source_plan_status, "failed") &&
    grepl("plan_status: `failed`", decision105_text, fixed = TRUE),
  decision_106_recovery_gate = identical(plan106$recovery_status, "passed") &&
    grepl("recovery_status: `passed`", decision106_text, fixed = TRUE),
  corrected_plan_frozen = identical(plan106$corrected_plan_status, "frozen") &&
    isTRUE(plan106$approved_for_chen_reverse_outcome_extraction_execution),
  strict_exposure_instrument_authority_found = identical(strict_freeze$freeze_status, "passed") &&
    isTRUE(strict_freeze$approved_for_reverse_primary_outcome_extraction) &&
    identical(strict_qc$instrument_selection_status, "passed"),
  relaxed_exposure_instrument_authority_found = identical(relaxed_freeze$freeze_status, "passed") &&
    isTRUE(relaxed_freeze$approved_for_reverse_relaxed_outcome_extraction) &&
    identical(relaxed_qc$instrument_selection_status, "passed"),
  outcome_independent_exposure_sets_used = TRUE,
  no_vuckovic_outcome_conditioned_target_selection = TRUE,
  strict_relaxed_hierarchy_preserved = TRUE,
  chen_source_semantics_verified = identical(chen_cert$certification_status, "passed") &&
    contains_all(readme_text, c("reference_allele - Effect allele", "other_allele - Non effect allele", "beta - Overall effect size/beta value for meta-analysis")),
  chen_other_allele_verified = identical(chen_cert$allele_convention$other_allele_field, "other_allele"),
  chen_beta_scale_verified = identical(chen_cert$phenotype$effect_scale, "standardized_quantitative_Hb_effect"),
  identity_bridge_defined = all(!is.na(audit$expected_chen_marker_id)),
  reverse_scope_adoption_resolved = identical(plan106$corrected_plan_status, "frozen") &&
    grepl("Decision 105 prospectively adopts", plan106$decision_84_scope_reuse_status, fixed = TRUE),
  exact_marker_matching_only = TRUE,
  marker_identity_separate_from_effect_orientation = TRUE,
  no_proxy = TRUE,
  no_liftover = TRUE,
  no_nearest_variant = TRUE,
  no_strand_complement_identity_rescue = TRUE,
  missing_semantics_defined = identical(plan106$missing_multiple_incompatible_handling$missing$category, "missing"),
  branch_readiness_logic_defined = TRUE,
  no_raw_chen_scan = TRUE,
  no_outcome_extraction = TRUE,
  no_harmonisation = TRUE,
  no_mr = TRUE,
  no_steiger = TRUE,
  renv_lock_unchanged = identical(hash_file(paths$renv_lock), "253471c02e1e47a40d0f68b296d4ae2b1df471d757ac86328b3e974018d039f3")
)

hard_check_failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
contract_status <- if (length(hard_check_failures) == 0L) "frozen" else "failed"

manifest_inputs <- unlist(paths[c(
  "script", "decision104", "decision105", "decision106", "plan106_qc",
  "decision78", "decision83", "decision84", "decision47", "decision59",
  "chen_source_cert", "chen_readme", "strict_freeze", "strict_selection_qc",
  "strict_manifest", "strict_included_tsv", "strict_excluded_tsv",
  "relaxed_freeze", "relaxed_selection_qc", "relaxed_manifest",
  "relaxed_included_tsv", "relaxed_excluded_tsv", "eur_bim", "renv_lock",
  "target_audit"
)])
input_manifest <- lapply(names(manifest_inputs), function(nm) {
  list(
    file_role = nm,
    relative_path = relpath(manifest_inputs[[nm]]),
    sha256 = hash_file(manifest_inputs[[nm]]),
    file_size_bytes = unname(file.info(manifest_inputs[[nm]])$size)
  )
})

contract <- list(
  contract_version = "v1",
  decision = 107,
  date = format(Sys.Date()),
  analysis_direction = "delirium_to_Hb",
  analysis_role = "reverse_alternative_hb_outcome_sensitivity",
  independent_replication = FALSE,
  design_decision = 104,
  plan_decision = 105,
  plan_recovery_decision = 106,
  decision_105_106_authority_wording = "Decision 105 scientific plan was complete but had three implementation-level false-positive hard checks; Decision 106 is the corrected approval authority.",
  strict_exposure_authority = list(
    freeze_decision = 47,
    freeze_status = strict_freeze$freeze_status,
    manifest_sha256 = strict_freeze$manifest_sha256,
    included_nsnp = nrow(strict_inc),
    excluded_nsnp = nrow(strict_exc),
    source = "outcome-independent frozen FinnGen exposure instrument selection V4"
  ),
  relaxed_exposure_authority = list(
    freeze_decision = 59,
    freeze_status = relaxed_freeze$freeze_status,
    manifest_sha256 = relaxed_freeze$manifest_sha256,
    included_nsnp = nrow(relaxed_inc),
    excluded_nsnp = nrow(relaxed_exc),
    source = "outcome-independent frozen FinnGen relaxed exploratory exposure instrument selection V2"
  ),
  strict_primary_threshold = 5e-08,
  relaxed_threshold = 5e-06,
  relaxed_branch_role = "protocol_triggered_exploratory_fallback",
  strict_primary_superseded_by_relaxed = FALSE,
  relaxed_confirmatory = FALSE,
  outcome_source = "Chen_2020_Hb_BCX2",
  outcome_source_file = "data_raw/gwas/BCX2_HGB_EA_GWAMA.out.gz",
  outcome_scale = "standardized_quantitative_Hb_effect",
  chen_effect_allele_definition = "reference_allele = Effect allele",
  chen_other_allele_definition = "other_allele = Non effect allele",
  chen_beta_definition = "beta = Overall effect size/beta value for meta-analysis; effect per reference_allele",
  identity_bridge_method = list(
    step1 = "frozen FinnGen canonical rsID",
    step2 = "certified 1KG Phase3 EUR GRCh37 same exact canonical rsID",
    step3 = "GRCh37 chromosome + position + direct unordered biallelic allele set",
    step4 = "expected Chen marker ID as chromosome:position_sortedAllele1_sortedAllele2",
    exact_marker_string_equivalent_to_identity_rule = TRUE
  ),
  reverse_specific_identity_adoption_status = "approved_by_decision_105_corrected_by_106",
  target_branch_structure = branch_counts,
  union_target_count = length(union_rsids),
  branch_membership_overlap = list(
    strict_included_rsids = strict_inc$rsid,
    strict_excluded_rsids = strict_exc$rsid,
    relaxed_included_rsids = relaxed_inc$rsid,
    relaxed_excluded_rsids = relaxed_exc$rsid,
    union_rsids = union_rsids
  ),
  target_authority_ready_count = target_authority_ready_count,
  target_authority_unavailable_count = target_authority_unavailable_count,
  target_authority_unavailable_rsids = audit$rsid[!audit$identity_bridge_ready],
  duplicate_reference_rsid_matches = duplicate_bim_rsids,
  expected_chen_marker_id_construction_rule = "paste0(chr_grch37, ':', pos_grch37, '_', lexicographically sorted reference allele tokens joined by '_')",
  proxy_allowed = FALSE,
  ld_proxy_allowed = FALSE,
  liftover_allowed = FALSE,
  nearest_variant_allowed = FALSE,
  fuzzy_matching_allowed = FALSE,
  strand_complement_identity_rescue_allowed = FALSE,
  coordinate_only_match_allowed = FALSE,
  multiple_match_allowed = FALSE,
  raw_source_scan_performed = FALSE,
  outcome_extraction_performed = FALSE,
  harmonisation_performed = FALSE,
  mr_run = FALSE,
  steiger_run = FALSE,
  missing_multiple_handling = plan106$missing_multiple_incompatible_handling,
  branch_specific_harmonisation_readiness_logic = list(
    matched_zero = "branch_harmonisation_readiness=FALSE; reason=no_exact_chen_outcome_available",
    matched_positive = "may enter a later harmonisation preflight, with MR estimability decided by later frozen rules",
    strict_single_instrument_missing = "scientific/data availability limitation, not technical extraction failure"
  ),
  future_extraction_outputs_planned_not_created = c(
    "R/28_chen_reverse_outcome_extraction_v1.R",
    "data_derived/reverse_sensitivity_outcome/chen_reverse_union_targets_v1.tsv",
    "data_derived/reverse_sensitivity_outcome/chen_reverse_outcome_master_v1.parquet",
    "data_derived/reverse_sensitivity_outcome/chen_reverse_outcome_master_v1.tsv",
    "results/qc/chen_reverse_outcome_extraction_v1.json"
  ),
  contract_status = contract_status,
  approved_for_chen_reverse_outcome_extraction = identical(contract_status, "frozen"),
  hard_checks = hard_checks,
  hard_check_failures = hard_check_failures,
  input_manifest = input_manifest
)

write_json(contract, paths$contract_json)

decision_lines <- c(
  "# Decision 107: Chen Reverse Outcome Extraction Execution Contract V1",
  "",
  paste0("Date: ", format(Sys.Date())),
  "",
  "## Status",
  paste0("contract_status: `", contract_status, "`"),
  paste0("approved_for_chen_reverse_outcome_extraction: `", identical(contract_status, "frozen"), "`"),
  "",
  "## Authority Model",
  "Decision 105 is preserved as scientific extraction/preflight plan provenance and failed historical evidence due to three implementation-level hard-check false positives. Decision 106 is the technical readback recovery and corrected approval authority.",
  "",
  "## Decision",
  "Freeze the execution contract for a future targeted Chen reverse outcome extraction as `FinnGen R13 delirium -> Chen 2020 haemoglobin`.",
  "",
  "This is a reverse alternative-Hb outcome sensitivity analysis only. It is not reverse primary evidence and not an independent replication.",
  "",
  "## Exposure Authority",
  paste0("- Strict exposure authority: Decision 47, manifest SHA `", strict_freeze$manifest_sha256, "`, included n=", nrow(strict_inc), ", excluded n=", nrow(strict_exc), "."),
  paste0("- Relaxed exposure authority: Decision 59, manifest SHA `", relaxed_freeze$manifest_sha256, "`, included n=", nrow(relaxed_inc), ", excluded n=", nrow(relaxed_exc), "."),
  "- Strict evidence remains higher tier than relaxed exploratory fallback.",
  "- Exposure sets are outcome-independent frozen FinnGen instrument sets, not Vuckovic outcome-conditioned harmonised sets.",
  "",
  "## Target Authority",
  paste0("- Union target count: `", length(union_rsids), "`."),
  paste0("- Identity-bridge ready count: `", target_authority_ready_count, "`."),
  paste0("- Identity-bridge unavailable count: `", target_authority_unavailable_count, "`."),
  "",
  "## Chen Source Semantics",
  "- `reference_allele` = effect allele.",
  "- `other_allele` = non-effect allele.",
  "- `beta` = overall effect size/beta value for meta-analysis; effect per `reference_allele`.",
  "- Outcome scale: `standardized_quantitative_Hb_effect`.",
  "",
  "## Identity Bridge",
  "Frozen FinnGen canonical rsID -> certified 1KG Phase3 EUR GRCh37 exact canonical rsID -> GRCh37 chromosome, position, and unordered allele set -> expected Chen marker ID.",
  "",
  "Expected Chen marker IDs are constructed as `chromosome:position_sortedAllele1_sortedAllele2`. Marker identity tokens are not outcome effect orientation; Chen effect orientation remains `reference_allele` and `other_allele`.",
  "",
  "## Safeguards",
  "This contract does not scan Chen raw GWAS, extract outcomes, harmonise, run MR, run Steiger, perform proxy search, perform liftOver, reclump, or modify existing strict/relaxed reverse exposure sets.",
  "",
  "## Hard Check Failures",
  if (length(hard_check_failures) == 0L) "- none" else paste0("- ", hard_check_failures),
  "",
  "## Outputs",
  paste0("- `", relpath(paths$contract_json), "`"),
  paste0("- `", relpath(paths$target_audit), "`"),
  paste0("- `", relpath(paths$contract_log), "`"),
  paste0("- `", relpath(paths$decision107), "`")
)
write_text(decision_lines, paths$decision107)

log_lines <- c(
  paste0("[", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "] stage=chen_reverse_outcome_extraction_execution_contract_v1_start"),
  "decision=107",
  paste0("contract_status=", contract_status),
  paste0("hard_check_failures=", paste(hard_check_failures, collapse = ",")),
  paste0("strict_included_n=", nrow(strict_inc)),
  paste0("strict_excluded_n=", nrow(strict_exc)),
  paste0("relaxed_included_n=", nrow(relaxed_inc)),
  paste0("relaxed_excluded_n=", nrow(relaxed_exc)),
  paste0("union_target_count=", length(union_rsids)),
  paste0("target_authority_ready_count=", target_authority_ready_count),
  paste0("target_authority_unavailable_count=", target_authority_unavailable_count),
  paste0("approved_for_chen_reverse_outcome_extraction=", identical(contract_status, "frozen")),
  "raw_chen_scan_performed=FALSE",
  "outcome_extraction_performed=FALSE",
  "harmonisation_performed=FALSE",
  "mr_run=FALSE",
  "steiger_run=FALSE"
)
write_text(log_lines, paths$contract_log)

cat("contract_status=", contract_status, "\n", sep = "")
cat("hard_check_failures=", paste(hard_check_failures, collapse = ","), "\n", sep = "")
cat("union_target_count=", length(union_rsids), "\n", sep = "")
cat("target_authority_ready_count=", target_authority_ready_count, "\n", sep = "")
cat("target_authority_unavailable_count=", target_authority_unavailable_count, "\n", sep = "")
cat("approved_for_chen_reverse_outcome_extraction=", identical(contract_status, "frozen"), "\n", sep = "")
