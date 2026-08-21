#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/27_chen_reverse_outcome_extraction_preflight_plan_v1.R [--project-root <path>]", call. = FALSE)
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
stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
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
  script = rel("R", "27_chen_reverse_outcome_extraction_preflight_plan_v1.R"),
  decision104 = rel("docs", "decisions", "104_chen_reverse_sensitivity_design_v1_v1.1.md"),
  design104_qc = rel("results", "qc", "chen_reverse_sensitivity_design_v1.json"),
  decision78 = rel("docs", "decisions", "78_chen_2020_hb_source_certification_v1_v1.1.md"),
  decision83 = rel("docs", "decisions", "83_chen_2020_hb_official_source_dictionary_audit_v1_v1.1.md"),
  decision84 = rel("docs", "decisions", "84_chen_identifier_resolution_amendment_v1_v1.1.md"),
  chen_source_cert = rel("results", "qc", "chen_2020_hb_source_certification_v1.json"),
  chen_readme = rel("docs", "source_metadata", "readme_BCX2_meta_analyses.txt"),
  reverse_primary_instrument_freeze = rel("results", "qc", "finngen_r13_delirium_reverse_primary_instruments_v4_freeze_v3.json"),
  reverse_strict_mr_freeze = rel("results", "qc", "reverse_strict_primary_mr_v1_freeze.json"),
  reverse_relaxed_mr_freeze = rel("results", "qc", "reverse_relaxed_mr_v1_freeze.json"),
  relaxed_harmonisation_freeze = rel("results", "qc", "vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_freeze.json"),
  reverse_relaxed_palindrome_rule = rel("results", "qc", "reverse_relaxed_palindromic_handling_rule_v1.json"),
  renv_lock = rel("renv.lock"),
  plan_qc = rel("results", "qc", "chen_reverse_outcome_extraction_preflight_plan_v1.json"),
  plan_log = rel("results", "logs", "chen_reverse_outcome_extraction_preflight_plan_v1.log"),
  decision105 = rel("docs", "decisions", "105_chen_reverse_outcome_extraction_preflight_plan_v1_v1.1.md")
)

required_inputs <- unlist(paths[c(
  "script", "decision104", "design104_qc", "decision78", "decision83", "decision84",
  "chen_source_cert", "chen_readme", "reverse_primary_instrument_freeze",
  "reverse_strict_mr_freeze", "reverse_relaxed_mr_freeze",
  "relaxed_harmonisation_freeze", "reverse_relaxed_palindrome_rule", "renv_lock"
)])
missing <- required_inputs[!file.exists(required_inputs)]
if (length(missing) > 0L) stop("Missing required input(s): ", paste(missing, collapse = "; "), call. = FALSE)

targets <- unlist(paths[c("plan_qc", "plan_log", "decision105")])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(occupied, collapse = "; ")))

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 105L), paste0("Expected next decision 105, found ", next_decision, "; no outputs written."))

text_file <- function(path) paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
decision78_text <- text_file(paths$decision78)
decision83_text <- text_file(paths$decision83)
decision84_text <- text_file(paths$decision84)
readme_text <- text_file(paths$chen_readme)
decision104_text <- text_file(paths$decision104)

design104 <- read_json(paths$design104_qc)
chen_cert <- read_json(paths$chen_source_cert)
reverse_primary <- read_json(paths$reverse_primary_instrument_freeze)
reverse_strict <- read_json(paths$reverse_strict_mr_freeze)
reverse_relaxed <- read_json(paths$reverse_relaxed_mr_freeze)
relaxed_harmonisation_freeze <- read_json(paths$relaxed_harmonisation_freeze)
palindrome_rule <- read_json(paths$reverse_relaxed_palindrome_rule)

contains_all <- function(text, terms) all(vapply(terms, grepl, logical(1), x = text, fixed = TRUE))

chen_effect_allele_definition <- "reference_allele = Effect allele"
chen_other_allele_definition <- "other_allele = Non effect allele"
chen_beta_definition <- "beta = Overall effect size/beta value for meta-analysis"
chen_beta_scale <- "standardized_quantitative_Hb_effect"
chen_eaf_definition <- "eaf = effect allele frequency"
chen_variant_n_definition <- "n_samples = Number of samples with marker present"

decision84_scope_reuse_status <- paste(
  "Decision 84 is not reused as a project-wide authority.",
  "Decision 105 prospectively adopts the same exact GRCh37 coordinate plus unordered allele-set identity philosophy",
  "for Chen as reverse outcome only, while preserving Decision 84's original forward-only scope."
)

identity_bridge_plan <- list(
  exposure_identity = "frozen FinnGen canonical rsID from existing strict/relaxed reverse exposure sets",
  bridge_reference = "certified 1000 Genomes Phase 3 EUR GRCh37 exact canonical-rsID identity",
  target_identity = "GRCh37 chromosome + position + unordered biallelic allele set",
  chen_match_rule = "exact chromosome + exact GRCh37 position + exact unordered marker allele set in Chen rs_number",
  forbidden = c("direct GRCh38-vs-GRCh37 coordinate comparison", "proxy", "LD proxy", "nearest variant", "fuzzy matching", "strand-complement rescue", "liftOver")
)

missing_handling <- list(
  category = "missing",
  definition = "frozen FinnGen exposure target unavailable by approved exact identity in Chen outcome source",
  action = "retain exposure-set provenance and report missing; do not modify exposure instruments"
)

multiple_handling <- list(
  category = "multiple_exact_match",
  definition = "more than one Chen exact coordinate plus unordered allele-set match or unresolved reference identity",
  action = "do not silently choose; block that target from later harmonisation unless separately approved"
)

incompatible_handling <- list(
  category = "marker_effect_allele_incompatible",
  definition = "Chen reference_allele/other_allele unordered set does not match the parsed marker unordered allele set",
  action = "classify and hold for audit; do not strand-complement rescue or flip beta"
)

hard_checks <- list(
  decision_104_design_gate = identical(design104$design_status, "passed") &&
    isTRUE(design104$approved_for_chen_reverse_sensitivity_outcome_extraction_plan) &&
    grepl("design_status: `passed`", decision104_text, fixed = TRUE),
  chen_source_certification_passed = identical(chen_cert$certification_status, "passed") &&
    identical(chen_cert$decision, 78),
  chen_official_dictionary_present = contains_all(readme_text, c(
    "reference_allele - Effect allele",
    "other_allele - Non effect allele",
    "eaf - effect allele frequency",
    "beta - Overall effect size/beta value for meta-analysis",
    "n_samples - Number of samples with marker present"
  )),
  chen_other_allele_reverified = identical(chen_cert$allele_convention$other_allele_field, "other_allele") &&
    grepl("other_allele - Non effect allele", readme_text, fixed = TRUE),
  chen_beta_scale_reverified = identical(chen_cert$phenotype$effect_scale, chen_beta_scale) &&
    grepl("Effect scale: standardized quantitative Hb effect", decision78_text, fixed = TRUE),
  documented_vs_inferred_scale_separated = isFALSE(chen_cert$phenotype$physical_unit_claim_allowed) &&
    identical(chen_cert$phenotype$forbidden_interpretations, c("g/dL", "g/L", "clinical haemoglobin concentration units")),
  strict_reverse_instruments_frozen_and_reused = identical(reverse_primary$freeze_status, "passed") &&
    isTRUE(reverse_primary$approved_for_reverse_primary_outcome_extraction),
  strict_reverse_mr_frozen = identical(reverse_strict$freeze_status, "passed") &&
    identical(reverse_strict$instrument_threshold, 5e-08),
  relaxed_reverse_instruments_frozen_and_reused = identical(relaxed_harmonisation_freeze$freeze_status, "passed"),
  relaxed_reverse_exploratory_frozen = identical(reverse_relaxed$freeze_status, "passed") &&
    identical(reverse_relaxed$p_threshold, 5e-06) &&
    identical(reverse_relaxed$branch_type, "protocol_triggered_exploratory_fallback"),
  no_exposure_reselection = TRUE,
  cross_build_identity_plan_defined = length(identity_bridge_plan) > 0L,
  no_direct_cross_build_coordinate_matching = TRUE,
  decision_84_scope_checked = grepl("It does not apply to Vuckovic, FinnGen, reverse analyses", decision84_text, fixed = TRUE),
  decision_84_reverse_adoption_is_new_decision = TRUE,
  no_proxy = TRUE,
  no_liftover = TRUE,
  no_nearest_variant = TRUE,
  no_fuzzy_matching = TRUE,
  no_strand_complement_rescue = TRUE,
  missing_handling_defined = identical(missing_handling$category, "missing"),
  multiple_handling_defined = identical(multiple_handling$category, "multiple_exact_match"),
  incompatible_handling_defined = identical(incompatible_handling$category, "marker_effect_allele_incompatible"),
  no_chen_raw_scan = TRUE,
  outcome_extraction_performed = FALSE,
  no_harmonisation = TRUE,
  no_mr = TRUE,
  no_steiger = TRUE
)

hard_check_failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
plan_status <- if (length(hard_check_failures) == 0L) "frozen" else "failed"

manifest_inputs <- unlist(paths[c(
  "script", "decision104", "design104_qc", "decision78", "decision83", "decision84",
  "chen_source_cert", "chen_readme", "reverse_primary_instrument_freeze",
  "reverse_strict_mr_freeze", "reverse_relaxed_mr_freeze",
  "relaxed_harmonisation_freeze", "reverse_relaxed_palindrome_rule", "renv_lock"
)])
input_manifest <- lapply(names(manifest_inputs), function(nm) {
  list(
    file_role = nm,
    relative_path = relpath(manifest_inputs[[nm]]),
    sha256 = hash_file(manifest_inputs[[nm]]),
    file_size_bytes = unname(file.info(manifest_inputs[[nm]])$size)
  )
})

plan <- list(
  plan_version = "v1",
  decision = 105,
  date = format(Sys.Date()),
  plan_status = plan_status,
  analysis_direction = "delirium_to_Hb",
  analysis_role = "reverse_alternative_hb_outcome_sensitivity",
  independent_replication = FALSE,
  source_authorities = list(
    design_gate_decision = 104,
    chen_source_certification_decision = 78,
    chen_official_dictionary_audit_decision = 83,
    chen_forward_identifier_resolution_decision = 84,
    reverse_strict_mr_freeze_decision = reverse_strict$decision,
    reverse_relaxed_mr_freeze_decision = reverse_relaxed$decision
  ),
  chen_outcome_semantics_reverification = list(
    source_path = chen_cert$source_path,
    source_sha256 = chen_cert$source_sha256,
    genome_build = chen_cert$genome_build,
    effect_allele_definition = chen_effect_allele_definition,
    other_allele_definition = chen_other_allele_definition,
    beta_definition = chen_beta_definition,
    beta_scale = chen_beta_scale,
    eaf_definition = chen_eaf_definition,
    variant_level_N_definition = chen_variant_n_definition,
    physical_unit_claim_allowed = FALSE,
    forbidden_interpretations = chen_cert$phenotype$forbidden_interpretations,
    source_documentation_evidence = c(
      relpath(paths$decision78),
      relpath(paths$decision83),
      relpath(paths$chen_readme),
      relpath(paths$chen_source_cert)
    )
  ),
  exposure_instrument_authority_plan = list(
    strict_branch_preserved = TRUE,
    relaxed_branch_preserved = TRUE,
    strict_threshold = 5e-08,
    relaxed_exploratory_threshold = 5e-06,
    exposure_reselection_allowed = FALSE,
    reclumping_allowed = FALSE,
    exposure_set_change_due_to_chen_outcome_availability_allowed = FALSE,
    evidence_hierarchy = "strict reverse primary > relaxed exploratory fallback"
  ),
  identity_bridge_plan = identity_bridge_plan,
  decision_84_scope_reuse_status = decision84_scope_reuse_status,
  targeted_outcome_extraction_future_plan = list(
    extraction_performed_now = FALSE,
    future_target_union = c("strict_APOE_included", "strict_APOE_excluded", "relaxed_APOE_included", "relaxed_APOE_excluded"),
    future_single_source_scan_preferred = TRUE,
    chen_match_status_categories = c(
      "unique_exact_match",
      "missing",
      "multiple_exact_match",
      "marker_effect_allele_incompatible",
      "invalid"
    ),
    branch_membership_must_be_retained = TRUE
  ),
  outcome_allele_value_qc_plan = list(
    retain_fields = c(
      "source_marker_id", "resolved_exposure_rsid", "reference_allele", "other_allele",
      "beta", "se", "eaf", "p-value", "n_samples", "GRCh37_coordinate_provenance"
    ),
    unordered_marker_vs_effect_allele_set_check_required = TRUE,
    beta_flipping_performed_in_extraction_stage = FALSE
  ),
  missing_multiple_incompatible_handling = list(
    missing = missing_handling,
    multiple = multiple_handling,
    incompatible = incompatible_handling
  ),
  harmonisation_preflight_readiness_plan = list(
    harmonisation_performed_now = FALSE,
    future_exposure_orientation_source = "frozen FinnGen reverse instruments",
    future_outcome_orientation_source = "Chen reference_allele effect allele and other_allele non-effect allele after this re-verification",
    palindrome_rule_status = "deferred; must be separately justified for Chen reverse sensitivity before harmonisation",
    steiger_status = "deferred"
  ),
  prohibited_now = list(
    proxy_allowed = FALSE,
    liftover_allowed = FALSE,
    nearest_variant_allowed = FALSE,
    fuzzy_matching_allowed = FALSE,
    strand_complement_rescue_allowed = FALSE,
    outcome_extraction_performed = FALSE,
    harmonisation_performed = FALSE,
    mr_run = FALSE,
    steiger_run = FALSE,
    chen_raw_scan_performed = FALSE
  ),
  hard_checks = hard_checks,
  hard_check_failures = hard_check_failures,
  input_manifest = input_manifest,
  raw_file_hash_recomputation = "not_recomputed_in_plan_stage; no Chen raw GWAS scan was performed",
  approved_for_chen_reverse_outcome_extraction_execution = identical(plan_status, "frozen")
)

write_json(plan, paths$plan_qc)

decision_lines <- c(
  "# Decision 105: Chen Reverse Outcome Extraction Preflight Plan V1",
  "",
  paste0("Date: ", format(Sys.Date())),
  "",
  "## Status",
  paste0("plan_status: `", plan_status, "`"),
  paste0("approved_for_chen_reverse_outcome_extraction_execution: `", identical(plan_status, "frozen"), "`"),
  "",
  "## Decision",
  "Freeze a plan for future Chen reverse sensitivity outcome extraction as `FinnGen R13 delirium -> Chen 2020 haemoglobin`.",
  "",
  "This is a reverse alternative-Hb outcome sensitivity analysis only. It is not reverse primary evidence and not an independent replication.",
  "",
  "This decision records source semantics, exposure-set authority, exact identity matching rules, missing/multiple/incompatible handling, and next-stage gates. It does not execute Chen outcome extraction.",
  "",
  "## Source Semantics",
  "- `reference_allele` = effect allele.",
  "- `other_allele` = non-effect allele.",
  "- `beta` = overall effect size/beta value for meta-analysis.",
  "- `eaf` = effect allele frequency.",
  "- `n_samples` = number of samples with marker present.",
  "- Beta scale is `standardized_quantitative_Hb_effect`; g/dL, g/L, and clinical-unit interpretations are not allowed.",
  "",
  "## Exposure Instrument Authority",
  "- Strict reverse instruments must be reused from the frozen `P < 5e-8` reverse exposure authority.",
  "- Relaxed instruments must be reused only as the frozen protocol-triggered exploratory `P < 5e-6` fallback.",
  "- Exposure reselection, reclumping, and modification due to Chen outcome availability are not allowed.",
  "",
  "## Identity Bridge",
  "- Start from frozen FinnGen canonical rsIDs.",
  "- Resolve through certified 1000 Genomes Phase 3 EUR GRCh37 exact canonical-rsID identity.",
  "- Match Chen markers only by exact GRCh37 chromosome, exact position, and exact unordered biallelic marker allele set.",
  "- Direct GRCh38-vs-GRCh37 coordinate comparison is not allowed.",
  "- Proxy, LD proxy, nearest-variant matching, fuzzy matching, strand-complement rescue, and liftOver are not allowed.",
  "",
  "## Decision 84 Scope",
  "Decision 84 remains forward-only. Decision 105 prospectively adopts the same exact GRCh37 coordinate plus unordered allele-set identity philosophy for Chen as a reverse outcome source, without changing Decision 84's original scope.",
  "",
  "## Missing, Multiple, And Incompatible Handling",
  "- Missing Chen outcome records must be reported and must not change the exposure instrument sets.",
  "- Multiple exact matches must not be silently selected.",
  "- Marker/effect allele incompatibility must be classified and held for audit; no strand-complement rescue or beta flip is allowed in extraction.",
  "",
  "## Safeguards",
  "This plan does not scan Chen raw GWAS, extract outcomes, harmonise, run MR, run Steiger, perform proxy search, perform liftOver, reclump, or modify existing strict/relaxed reverse exposure sets.",
  "",
  "## Hard Check Failures",
  if (length(hard_check_failures) == 0L) "- none" else paste0("- ", hard_check_failures),
  "",
  "## Outputs",
  paste0("- `", relpath(paths$plan_qc), "`"),
  paste0("- `", relpath(paths$plan_log), "`"),
  paste0("- `", relpath(paths$decision105), "`")
)
write_text(decision_lines, paths$decision105)

log_lines <- c(
  paste0("[", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "] stage=chen_reverse_outcome_extraction_preflight_plan_v1_start"),
  "decision=105",
  paste0("plan_status=", plan_status),
  paste0("hard_check_failures=", paste(hard_check_failures, collapse = ",")),
  paste0("approved_for_chen_reverse_outcome_extraction_execution=", identical(plan_status, "frozen")),
  "chen_raw_scan_performed=FALSE",
  "outcome_extraction_performed=FALSE",
  "harmonisation_performed=FALSE",
  "mr_run=FALSE",
  "steiger_run=FALSE"
)
write_text(log_lines, paths$plan_log)

cat("plan_status=", plan_status, "\n", sep = "")
cat("hard_check_failures=", paste(hard_check_failures, collapse = ","), "\n", sep = "")
cat("approved_for_chen_reverse_outcome_extraction_execution=", identical(plan_status, "frozen"), "\n", sep = "")
