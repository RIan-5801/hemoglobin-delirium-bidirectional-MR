#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/26_chen_reverse_sensitivity_design_v1.R [--project-root <path>]", call. = FALSE)
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
  script = rel("R", "26_chen_reverse_sensitivity_design_v1.R"),
  metadata = rel("docs", "02_gwas_metadata_v2.md"),
  protocol_freeze = rel("docs", "decisions", "04_protocol_freeze_v1.1.md"),
  chen_forward_freeze_decision = rel("docs", "decisions", "103_chen_forward_mr_v1_freeze_v1.1.md"),
  chen_forward_freeze_qc = rel("results", "qc", "chen_forward_mr_v1_freeze.json"),
  reverse_primary_instrument_freeze_qc = rel("results", "qc", "finngen_r13_delirium_reverse_primary_instruments_v4_freeze_v3.json"),
  reverse_strict_mr_freeze_qc = rel("results", "qc", "reverse_strict_primary_mr_v1_freeze.json"),
  reverse_relaxed_mr_freeze_qc = rel("results", "qc", "reverse_relaxed_mr_v1_freeze.json"),
  renv_lock = rel("renv.lock"),
  design_qc = rel("results", "qc", "chen_reverse_sensitivity_design_v1.json"),
  design_log = rel("results", "logs", "chen_reverse_sensitivity_design_v1.log"),
  decision = rel("docs", "decisions", "104_chen_reverse_sensitivity_design_v1_v1.1.md")
)

required_inputs <- unlist(paths[c(
  "script", "metadata", "protocol_freeze", "chen_forward_freeze_decision",
  "chen_forward_freeze_qc", "reverse_primary_instrument_freeze_qc",
  "reverse_strict_mr_freeze_qc", "reverse_relaxed_mr_freeze_qc", "renv_lock"
)])
missing <- required_inputs[!file.exists(required_inputs)]
if (length(missing) > 0L) stop("Missing required input(s): ", paste(missing, collapse = "; "), call. = FALSE)

targets <- unlist(paths[c("design_qc", "design_log", "decision")])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(occupied, collapse = "; ")))

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 104L), paste0("Expected next decision 104, found ", next_decision, "; no outputs written."))

metadata_lines <- readLines(paths$metadata, warn = FALSE, encoding = "UTF-8")
protocol_lines <- readLines(paths$protocol_freeze, warn = FALSE, encoding = "UTF-8")
decision103_lines <- readLines(paths$chen_forward_freeze_decision, warn = FALSE, encoding = "UTF-8")
chen_forward_freeze <- read_json(paths$chen_forward_freeze_qc)
reverse_primary_instruments <- read_json(paths$reverse_primary_instrument_freeze_qc)
reverse_strict_mr <- read_json(paths$reverse_strict_mr_freeze_qc)
reverse_relaxed_mr <- read_json(paths$reverse_relaxed_mr_freeze_qc)

expected <- list(
  finngen_raw = "data_raw/gwas/finngen_R13_F5_DELIRIUM.gz",
  finngen_sha = "85637F0F3358807964D4F8A3E500293168A706F1C08C65F3FC5512B65DF40ED8",
  chen_raw = "data_raw/gwas/BCX2_HGB_EA_GWAMA.out.gz",
  chen_sha = "F1DFEA8897CB29F39D891B7922BED2EA95A869BB864D39C56DB73E5D69F8ABF8"
)

metadata_text <- paste(metadata_lines, collapse = "\n")
protocol_text <- paste(protocol_lines, collapse = "\n")
decision103_text <- paste(decision103_lines, collapse = "\n")

hard_checks <- list(
  decision_103_passed = identical(chen_forward_freeze$freeze_status, "passed"),
  decision_103_approved_design = isTRUE(chen_forward_freeze$approved_for_chen_reverse_sensitivity_design) &&
    grepl("approved_for_chen_reverse_sensitivity_design: `TRUE`", decision103_text, fixed = TRUE),
  reverse_primary_instruments_frozen = identical(reverse_primary_instruments$freeze_status, "passed") &&
    isTRUE(reverse_primary_instruments$approved_for_reverse_primary_outcome_extraction),
  reverse_strict_mr_frozen = identical(reverse_strict_mr$freeze_status, "passed"),
  reverse_relaxed_mr_frozen = identical(reverse_relaxed_mr$freeze_status, "passed"),
  chen_forward_not_independent_replication = identical(chen_forward_freeze$independent_replication, FALSE),
  metadata_finngen_hash_present = grepl(expected$finngen_raw, metadata_text, fixed = TRUE) &&
    grepl(expected$finngen_sha, metadata_text, fixed = TRUE),
  metadata_chen_hash_present = grepl(expected$chen_raw, metadata_text, fixed = TRUE) &&
    grepl(expected$chen_sha, metadata_text, fixed = TRUE),
  finngen_effect_allele_documented = grepl("Effect allele: `alt`", metadata_text, fixed = TRUE),
  chen_effect_allele_documented = grepl("Effect allele: `reference_allele`", metadata_text, fixed = TRUE),
  chen_other_allele_requires_recheck = grepl("Other-allele field definition: must be rechecked", metadata_text, fixed = TRUE),
  reverse_sensitivity_role_frozen = grepl("FinnGen", protocol_text, fixed = TRUE) &&
    grepl("Chen is reverse sensitivity", protocol_text, fixed = TRUE),
  strict_reverse_threshold_preserved = identical(reverse_strict_mr$instrument_threshold, 5e-08),
  relaxed_reverse_exploratory_preserved = identical(reverse_relaxed_mr$p_threshold, 5e-06) &&
    identical(reverse_relaxed_mr$branch_type, "protocol_triggered_exploratory_fallback"),
  no_raw_data_modification = TRUE,
  no_gwas_scan = TRUE,
  no_instrument_reselection = TRUE,
  no_ld_clumping = TRUE,
  no_outcome_extraction = TRUE,
  no_proxy_or_liftover = TRUE,
  no_harmonisation = TRUE,
  no_mr = TRUE,
  no_steiger = TRUE,
  output_paths_versioned_and_new = TRUE
)

hard_check_failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
design_status <- if (length(hard_check_failures) == 0L) "passed" else "failed"

manifest_inputs <- unlist(paths[c(
  "script", "metadata", "protocol_freeze", "chen_forward_freeze_decision",
  "chen_forward_freeze_qc", "reverse_primary_instrument_freeze_qc",
  "reverse_strict_mr_freeze_qc", "reverse_relaxed_mr_freeze_qc", "renv_lock"
)])
input_manifest <- lapply(names(manifest_inputs), function(nm) {
  list(
    file_role = nm,
    relative_path = relpath(manifest_inputs[[nm]]),
    sha256 = hash_file(manifest_inputs[[nm]]),
    file_size_bytes = unname(file.info(manifest_inputs[[nm]])$size)
  )
})

qc <- list(
  design_version = "v1",
  decision = 104,
  date = format(Sys.Date()),
  design_status = design_status,
  analysis_direction = "delirium_to_Hb",
  analysis_role = "reverse_sensitivity",
  exposure = list(
    dataset = "FinnGen R13 F5_DELIRIUM",
    raw_file = expected$finngen_raw,
    build = "GRCh38",
    effect_allele = "alt",
    phenotype_scale = "log odds of delirium",
    cases = 5121,
    controls = 465023,
    study_n = 470144,
    authenticated_sha256 = expected$finngen_sha
  ),
  outcome = list(
    dataset = "Chen 2020 BCX2 European haemoglobin",
    raw_file = expected$chen_raw,
    build = "GRCh37",
    effect_allele = "reference_allele",
    other_allele_status = "must_be_rechecked_against_authenticated_source_field_dictionary_before_harmonisation",
    phenotype_scale = "haemoglobin concentration",
    sample_size_rule = "use variant-level n_samples where downstream implementation requires variant-level N",
    authenticated_sha256 = expected$chen_sha
  ),
  source_authorities = list(
    protocol_freeze_decision = 4,
    chen_forward_freeze_decision = 103,
    reverse_primary_instrument_freeze = "finngen_r13_delirium_reverse_primary_instruments_v4_freeze_v3",
    reverse_strict_mr_freeze_decision = reverse_strict_mr$decision,
    reverse_relaxed_mr_freeze_decision = reverse_relaxed_mr$decision
  ),
  frozen_parameters = list(
    strict_reverse_threshold = 5e-08,
    relaxed_exploratory_threshold = 5e-06,
    relaxed_threshold_label = "protocol_triggered_exploratory_fallback_only",
    ld_clumping_r2 = 0.001,
    ld_clumping_window_kb = 10000,
    ld_panel_ancestry = "EUR",
    apoe_policy = "maintain APOE-included and APOE-excluded branches as prespecified non-result-driven branches",
    main_instrument_universe = "autosomal biallelic A/C/G/T SNPs with authenticated build and alleles",
    weak_instrument_policy = "report F < 10 rather than silently remove unless separately approved"
  ),
  stage_boundary = list(
    current_stage = "design_preflight_only",
    approved_actions = c("record design decision", "lock source authorities", "define next-stage gates"),
    prohibited_in_this_stage = c(
      "full GWAS scan",
      "P-value selection",
      "LD clumping",
      "proxy search",
      "liftOver",
      "outcome extraction",
      "harmonisation",
      "MR",
      "Steiger",
      "manuscript conclusion change"
    ),
    next_stage_requires_new_A_to_H_plan_and_user_approval = TRUE
  ),
  next_stage_gate = list(
    proposed_next_stage = "Chen reverse sensitivity outcome extraction/preflight",
    must_verify_before_execution = c(
      "Chen source field dictionary for other allele and beta scale",
      "rsID/coordinate identity strategy using authenticated metadata only",
      "GRCh38 exposure versus GRCh37 outcome handling without undocumented liftOver",
      "target paths are new and versioned",
      "strict and exploratory branches remain separately labelled",
      "APOE included/excluded branches remain non-result-driven"
    )
  ),
  input_manifest = input_manifest,
  raw_file_hash_recomputation = "not_recomputed_in_design_stage; authenticated SHA-256 values were read from frozen metadata, and no raw GWAS processing was performed",
  hard_checks = hard_checks,
  hard_check_failures = hard_check_failures,
  approved_for_chen_reverse_sensitivity_outcome_extraction_plan = identical(design_status, "passed")
)

write_json(qc, paths$design_qc)

decision_lines <- c(
  "# Decision 104: Chen Reverse Sensitivity Design V1",
  "",
  paste0("Date: ", format(Sys.Date())),
  "",
  "## Status",
  paste0("design_status: `", design_status, "`"),
  paste0("approved_for_chen_reverse_sensitivity_outcome_extraction_plan: `", identical(design_status, "passed"), "`"),
  "",
  "## Decision",
  "Open the Chen reverse sensitivity pathway as `FinnGen R13 delirium -> Chen 2020 haemoglobin`.",
  "",
  "This is a reverse sensitivity analysis only. It is not primary evidence and not an independent replication.",
  "",
  "This decision records design/preflight boundaries only. It does not authorize outcome extraction, harmonisation, MR, Steiger, proxy search, liftOver, reclumping, or instrument reselection.",
  "",
  "## Rationale",
  "- Decision 103 passed and approved Chen reverse sensitivity design.",
  "- Chen 2020 is the prespecified alternative Hb GWAS for reverse sensitivity.",
  "- Existing reverse strict and relaxed Vuckovic analyses remain frozen and unchanged.",
  "",
  "## Frozen Parameters",
  "- Strict reverse threshold: `P < 5e-8`.",
  "- Exploratory relaxed threshold: `P < 5e-6`, labelled as protocol-triggered exploratory fallback only.",
  "- LD clumping, if later approved: `r2 < 0.001`, `10000 kb`, EUR LD panel.",
  "- APOE branches must remain APOE-included and APOE-excluded, with no result-driven exclusion.",
  "",
  "## Required Next-Stage Gates",
  "- Recheck Chen source field dictionary for other allele and beta scale before harmonisation.",
  "- Keep FinnGen GRCh38 and Chen GRCh37 build handling explicit; no undocumented liftOver or coordinate substitution.",
  "- Use only authenticated metadata for effect allele, other allele, rsID/coordinate identity, trait definition, and phenotype scale.",
  "- Prepare a new A-H plan and obtain user approval before any extraction or harmonisation.",
  "",
  "## Hard Check Failures",
  if (length(hard_check_failures) == 0L) "- none" else paste0("- ", hard_check_failures),
  "",
  "## Outputs",
  paste0("- `", relpath(paths$design_qc), "`"),
  paste0("- `", relpath(paths$design_log), "`"),
  paste0("- `", relpath(paths$decision), "`")
)
write_text(decision_lines, paths$decision)

log_lines <- c(
  paste0("[", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "] stage=chen_reverse_sensitivity_design_v1_start"),
  paste0("decision=104"),
  paste0("design_status=", design_status),
  paste0("hard_check_failures=", paste(hard_check_failures, collapse = ",")),
  paste0("approved_for_chen_reverse_sensitivity_outcome_extraction_plan=", identical(design_status, "passed")),
  paste0("no_raw_data_modification=TRUE"),
  paste0("no_gwas_scan=TRUE"),
  paste0("no_harmonisation=TRUE"),
  paste0("no_mr=TRUE")
)
write_text(log_lines, paths$design_log)

cat("design_status=", design_status, "\n", sep = "")
cat("hard_check_failures=", paste(hard_check_failures, collapse = ","), "\n", sep = "")
cat("approved_for_chen_reverse_sensitivity_outcome_extraction_plan=", identical(design_status, "passed"), "\n", sep = "")
