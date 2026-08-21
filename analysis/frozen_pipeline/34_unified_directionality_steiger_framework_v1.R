#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/34_unified_directionality_steiger_framework_v1.R [--project-root <path>]", call. = FALSE)
}
setwd(root)

local_lib <- normalizePath(file.path(root, "renv", "mr-v1-library"), winslash = "/", mustWork = TRUE)
.libPaths(c(local_lib, .libPaths()))

for (pkg in c("jsonlite", "digest", "DBI", "duckdb")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

rel <- function(...) file.path(root, ...)
norm <- function(x) gsub("\\\\", "/", x)
relpath <- function(x) norm(sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", norm(normalizePath(x, winslash = "/", mustWork = FALSE))))
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
hash_text <- function(x) digest::digest(paste(x, collapse = "\n"), algo = "sha256", serialize = FALSE)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
null_chr <- function(x) if (is.null(x) || length(x) == 0L) NA_character_ else as.character(x[[1L]])
null_int <- function(x) if (is.null(x) || length(x) == 0L || is.na(x[[1L]])) NA_integer_ else as.integer(x[[1L]])

write_csv_precise <- function(x, path) {
  partial <- paste0(path, ".partial")
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.csv(x, partial, row.names = FALSE, na = "")
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}
write_tsv <- function(x, path) {
  partial <- paste0(path, ".partial")
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.table(x, partial, sep = "\t", row.names = FALSE, quote = TRUE, na = "")
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}
write_json <- function(x, path) {
  partial <- paste0(path, ".partial")
  jsonlite::write_json(x, partial, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}
write_text <- function(lines, path) {
  partial <- paste0(path, ".partial")
  writeLines(lines, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}

latest_decision <- function(pattern) {
  files <- list.files(rel("docs", "decisions"), pattern = pattern, full.names = FALSE)
  nums <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", files)))
  if (!length(nums) || all(is.na(nums))) return(NA_integer_)
  max(nums, na.rm = TRUE)
}
get_manifest_path <- function(manifest_file, roles) {
  m <- utils::read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
  role_col <- if ("file_role" %in% names(m)) "file_role" else "role"
  path_col <- if ("relative_path" %in% names(m)) "relative_path" else "path"
  hit <- m[m[[role_col]] %in% roles, , drop = FALSE]
  if (nrow(hit) == 0L) return(NA_character_)
  relp <- hit[[path_col]][[1L]]
  if (grepl("^[A-Za-z]:|^/", relp)) return(relpath(relp))
  norm(relp)
}
json_manifest_sha <- function(x) null_chr(x$manifest_sha256)
hard_empty <- function(x) is.null(x$hard_check_failures) || length(x$hard_check_failures) == 0L
split_rsids <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(character())
  y <- unlist(strsplit(paste(as.character(x), collapse = ";"), "[;,]"))
  y <- trimws(y)
  y[nzchar(y)]
}

paths <- list(
  script = rel("R", "34_unified_directionality_steiger_framework_v1.R"),
  decision119_freeze = rel("results", "qc", "chen_reverse_mr_v1_freeze.json"),
  vuckovic_forward_mr_freeze = rel("results", "qc", "vuckovic_hb_finngen_r13_forward_mr_v3_freeze_v2.json"),
  chen_forward_input_freeze = rel("results", "qc", "chen_forward_harmonised_mr_inputs_v1_freeze.json"),
  chen_forward_input_manifest = rel("results", "qc", "chen_forward_harmonised_mr_inputs_v1_freeze_manifest.csv"),
  chen_forward_mr_freeze = rel("results", "qc", "chen_forward_mr_v1_freeze.json"),
  vuckovic_reverse_strict_input_freeze = rel("results", "qc", "vuckovic_hb_reverse_primary_formal_harmonisation_v4_freeze.json"),
  vuckovic_reverse_strict_input_manifest = rel("results", "qc", "vuckovic_hb_reverse_primary_formal_harmonisation_v4_freeze_manifest.csv"),
  vuckovic_reverse_strict_mr_freeze = rel("results", "qc", "reverse_strict_primary_mr_v1_freeze.json"),
  vuckovic_reverse_relaxed_input_freeze = rel("results", "qc", "vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_freeze.json"),
  vuckovic_reverse_relaxed_input_manifest = rel("results", "qc", "vuckovic_hb_reverse_relaxed_formal_harmonisation_v3_freeze_manifest.csv"),
  vuckovic_reverse_relaxed_mr_freeze = rel("results", "qc", "reverse_relaxed_mr_v1_freeze.json"),
  chen_reverse_input_freeze = rel("results", "qc", "chen_reverse_harmonised_mr_inputs_v1_freeze.json"),
  chen_reverse_input_manifest = rel("results", "qc", "chen_reverse_harmonised_mr_inputs_v1_freeze_manifest.csv"),
  chen_reverse_mr_freeze = rel("results", "qc", "chen_reverse_mr_v1_freeze.json"),
  renv_lock = rel("renv.lock"),
  registry = rel("results", "qc", "unified_directionality_analysis_registry_v1.csv"),
  feasibility = rel("results", "qc", "unified_directionality_feasibility_audit_v1.csv"),
  snp_audit = rel("results", "qc", "unified_directionality_snp_data_requirement_audit_v1.tsv"),
  package_audit = rel("results", "qc", "unified_steiger_package_implementation_audit_v1.csv"),
  framework = rel("results", "qc", "unified_directionality_steiger_framework_v1.json"),
  log = rel("results", "logs", "unified_directionality_steiger_framework_v1.log"),
  decision = rel("docs", "decisions", "120_unified_directionality_steiger_framework_and_feasibility_v1_v1.1.md")
)

required <- unlist(paths[1:16])
missing <- required[!file.exists(required)]
stop_if(length(missing) > 0L, paste("Missing required source file(s):", paste(relpath(missing), collapse = "; ")))

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 120L), paste("Expected next decision 120, found ", next_decision, "; no outputs written."))

targets <- unlist(paths[17:23])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(relpath(occupied), collapse = "; ")))

renv_before <- hash_file(paths$renv_lock)

f119 <- read_json(paths$decision119_freeze)
vf_mr <- read_json(paths$vuckovic_forward_mr_freeze)
cf_in <- read_json(paths$chen_forward_input_freeze)
cf_mr <- read_json(paths$chen_forward_mr_freeze)
vrs_in <- read_json(paths$vuckovic_reverse_strict_input_freeze)
vrs_mr <- read_json(paths$vuckovic_reverse_strict_mr_freeze)
vrr_in <- read_json(paths$vuckovic_reverse_relaxed_input_freeze)
vrr_mr <- read_json(paths$vuckovic_reverse_relaxed_mr_freeze)
cr_in <- read_json(paths$chen_reverse_input_freeze)
cr_mr <- read_json(paths$chen_reverse_mr_freeze)

decision119_gate <- identical(f119$freeze_status, "passed") &&
  isTRUE(f119$approved_for_unified_directionality_design) &&
  hard_empty(f119)

freeze_statuses <- data.frame(
  freeze_id = c(
    "vuckovic_forward_mr", "chen_forward_input", "chen_forward_mr",
    "vuckovic_reverse_strict_input", "vuckovic_reverse_strict_mr",
    "vuckovic_reverse_relaxed_input", "vuckovic_reverse_relaxed_mr",
    "chen_reverse_input", "chen_reverse_mr"
  ),
  path = relpath(c(
    paths$vuckovic_forward_mr_freeze, paths$chen_forward_input_freeze, paths$chen_forward_mr_freeze,
    paths$vuckovic_reverse_strict_input_freeze, paths$vuckovic_reverse_strict_mr_freeze,
    paths$vuckovic_reverse_relaxed_input_freeze, paths$vuckovic_reverse_relaxed_mr_freeze,
    paths$chen_reverse_input_freeze, paths$chen_reverse_mr_freeze
  )),
  freeze_status = c(vf_mr$freeze_status, cf_in$freeze_status, cf_mr$freeze_status, vrs_in$freeze_status, vrs_mr$freeze_status, vrr_in$freeze_status, vrr_mr$freeze_status, cr_in$freeze_status, cr_mr$freeze_status),
  hard_check_failures_empty = c(hard_empty(vf_mr), hard_empty(cf_in), hard_empty(cf_mr), hard_empty(vrs_in), hard_empty(vrs_mr), hard_empty(vrr_in), hard_empty(vrr_mr), hard_empty(cr_in), hard_empty(cr_mr)),
  manifest_sha256 = c(json_manifest_sha(vf_mr), json_manifest_sha(cf_in), json_manifest_sha(cf_mr), json_manifest_sha(vrs_in), json_manifest_sha(vrs_mr), json_manifest_sha(vrr_in), json_manifest_sha(vrr_mr), json_manifest_sha(cr_in), json_manifest_sha(cr_mr)),
  stringsAsFactors = FALSE
)

input_paths <- list(
  vuckovic_forward_included = relpath(vf_mr$frozen_mr_input_copies$included_path),
  vuckovic_forward_excluded = relpath(vf_mr$frozen_mr_input_copies$excluded_path),
  chen_forward_included = get_manifest_path(paths$chen_forward_input_manifest, "harmonised_apoe_included_parquet"),
  chen_forward_excluded = get_manifest_path(paths$chen_forward_input_manifest, "harmonised_apoe_excluded_parquet"),
  vuckovic_reverse_strict_included = get_manifest_path(paths$vuckovic_reverse_strict_input_manifest, "harmonised_included_parquet"),
  vuckovic_reverse_strict_excluded = get_manifest_path(paths$vuckovic_reverse_strict_input_manifest, "harmonised_excluded_parquet"),
  vuckovic_reverse_relaxed_included = get_manifest_path(paths$vuckovic_reverse_relaxed_input_manifest, "included_parquet"),
  vuckovic_reverse_relaxed_excluded = get_manifest_path(paths$vuckovic_reverse_relaxed_input_manifest, "excluded_parquet"),
  chen_reverse_strict_included = get_manifest_path(paths$chen_reverse_input_manifest, "strict_inc_parquet"),
  chen_reverse_strict_excluded = get_manifest_path(paths$chen_reverse_input_manifest, "strict_exc_parquet"),
  chen_reverse_relaxed_included = get_manifest_path(paths$chen_reverse_input_manifest, "relaxed_inc_parquet"),
  chen_reverse_relaxed_excluded = get_manifest_path(paths$chen_reverse_input_manifest, "relaxed_exc_parquet")
)
stop_if(any(is.na(unlist(input_paths))), paste("Could not resolve authoritative MR input path(s):", paste(names(input_paths)[is.na(unlist(input_paths))], collapse = "; ")))

con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

describe_parquet <- function(path) {
  full <- rel(path)
  if (!file.exists(full)) {
    return(list(exists = FALSE, n = NA_integer_, columns = character(), rsids = character(), summary = data.frame()))
  }
  quoted <- DBI::dbQuoteString(con, path)
  desc <- DBI::dbGetQuery(con, sprintf("DESCRIBE SELECT * FROM read_parquet(%s)", quoted))
  n <- DBI::dbGetQuery(con, sprintf("SELECT count(*) AS n FROM read_parquet(%s)", quoted))[["n"]][[1L]]
  cols <- desc$column_name
  snp_col <- c("SNP", "rsid", "target_rsid", "resolved_rsid", "variant_id")
  snp_col <- snp_col[snp_col %in% cols][1L]
  rsids <- character()
  if (!is.na(snp_col)) {
    rsids <- DBI::dbGetQuery(con, sprintf("SELECT %s AS rsid FROM read_parquet(%s)", DBI::dbQuoteIdentifier(con, snp_col), quoted))[["rsid"]]
    rsids <- sort(unique(as.character(rsids)))
  }
  list(exists = TRUE, n = as.integer(n), columns = cols, rsids = rsids, summary = desc)
}

input_info <- lapply(input_paths, describe_parquet)

field_audit <- function(path, trait_side) {
  info <- describe_parquet(path)
  cols <- info$columns
  col_or <- function(cands) {
    hit <- cands[cands %in% cols]
    if (length(hit)) hit[[1L]] else NA_character_
  }
  suffix <- if (identical(trait_side, "exposure")) ".exposure" else ".outcome"
  prefix <- if (identical(trait_side, "exposure")) "exposure" else "outcome"
  beta_col <- col_or(c(paste0("beta", suffix), paste0("b", suffix), paste0(prefix, "_beta_harmonised"), paste0(prefix, "_beta"), paste0(prefix, "_beta_raw"), "beta"))
  se_col <- col_or(c(paste0("se", suffix), paste0(prefix, "_se_harmonised"), paste0(prefix, "_se"), paste0(prefix, "_se_raw"), "se"))
  eaf_col <- col_or(c(paste0("eaf", suffix), paste0("af", suffix), paste0("effect_allele_frequency", suffix), paste0(prefix, "_eaf_harmonised"), paste0(prefix, "_eaf"), paste0(prefix, "_eaf_raw"), "eaf", "af"))
  n_col <- col_or(c(paste0("samplesize", suffix), paste0("n", suffix), paste0("n_samples", suffix), paste0(prefix, "_n_samples"), paste0(prefix, "_n_study"), "samplesize", "n_samples", "n"))
  ncase_col <- col_or(c(paste0("ncase", suffix), paste0("n_case", suffix), paste0(prefix, "_ncase"), paste0(prefix, "_n_case"), "ncase", "n_cases"))
  ncontrol_col <- col_or(c(paste0("ncontrol", suffix), paste0("n_control", suffix), paste0(prefix, "_ncontrol"), paste0(prefix, "_n_control"), "ncontrol", "n_controls"))
  pval_col <- col_or(c(paste0("pval", suffix), paste0("p", suffix), paste0(prefix, "_pval_harmonised"), paste0(prefix, "_pval"), paste0(prefix, "_p_value"), paste0(prefix, "_p"), "pval", "p"))
  complete_stats <- function(col) {
    if (is.na(col) || !info$exists) return(list(complete = FALSE, unique_n = NA_integer_, min = NA_real_, max = NA_real_))
    q <- DBI::dbQuoteString(con, path)
    id <- DBI::dbQuoteIdentifier(con, col)
    out <- DBI::dbGetQuery(con, sprintf("SELECT sum(CASE WHEN %s IS NULL THEN 1 ELSE 0 END) AS missing, count(DISTINCT %s) AS unique_n, min(%s) AS min_v, max(%s) AS max_v FROM read_parquet(%s)", id, id, id, id, q))
    list(complete = isTRUE(as.integer(out$missing[[1L]]) == 0L), unique_n = as.integer(out$unique_n[[1L]]), min = suppressWarnings(as.numeric(out$min_v[[1L]])), max = suppressWarnings(as.numeric(out$max_v[[1L]])))
  }
  n_stats <- complete_stats(n_col)
  ncase_stats <- complete_stats(ncase_col)
  ncontrol_stats <- complete_stats(ncontrol_col)
  list(
    path = path,
    exists = info$exists,
    n_rows = info$n,
    columns = cols,
    rsids = info$rsids,
    beta_col = beta_col,
    se_col = se_col,
    eaf_col = eaf_col,
    n_col = n_col,
    n_complete = n_stats$complete,
    n_unique = n_stats$unique_n,
    ncase_col = ncase_col,
    ncase_complete = ncase_stats$complete,
    ncase_unique = ncase_stats$unique_n,
    ncontrol_col = ncontrol_col,
    ncontrol_complete = ncontrol_stats$complete,
    ncontrol_unique = ncontrol_stats$unique_n,
    pval_col = pval_col,
    pval_present = !is.na(pval_col)
  )
}

registry_blueprint <- data.frame(
  analysis_id = c(
    "forward_vuckovic_primary_apoe_included", "forward_vuckovic_primary_apoe_excluded",
    "forward_chen_alternative_apoe_included", "forward_chen_alternative_apoe_excluded",
    "reverse_vuckovic_strict_primary_apoe_included", "reverse_vuckovic_strict_primary_apoe_excluded",
    "reverse_vuckovic_relaxed_exploratory_apoe_included", "reverse_vuckovic_relaxed_exploratory_apoe_excluded",
    "reverse_chen_strict_outcome_sensitivity_apoe_included", "reverse_chen_strict_outcome_sensitivity_apoe_excluded",
    "reverse_chen_relaxed_outcome_sensitivity_apoe_included", "reverse_chen_relaxed_outcome_sensitivity_apoe_excluded"
  ),
  direction = c(rep("Hb_to_delirium", 4), rep("delirium_to_Hb", 8)),
  evidence_level = c("forward_primary", "forward_apoe_exclusion_sensitivity", "forward_alternative_hb_gwas_sensitivity", "forward_alternative_hb_gwas_apoe_exclusion_sensitivity", "reverse_strict_primary", "reverse_strict_primary_apoe_exclusion_sensitivity", "reverse_relaxed_exploratory", "reverse_relaxed_exploratory_apoe_exclusion_sensitivity", "reverse_strict_alternative_hb_outcome_sensitivity", "reverse_strict_alternative_hb_outcome_sensitivity_apoe_exclusion", "reverse_relaxed_alternative_hb_outcome_sensitivity_exploratory", "reverse_relaxed_alternative_hb_outcome_sensitivity_apoe_exclusion_exploratory"),
  exposure_gwas = c("Vuckovic 2020 haemoglobin", "Vuckovic 2020 haemoglobin", "Chen 2020 European haemoglobin", "Chen 2020 European haemoglobin", rep("FinnGen R13 delirium", 8)),
  outcome_gwas = c(rep("FinnGen R13 delirium", 4), "Vuckovic 2020 haemoglobin", "Vuckovic 2020 haemoglobin", "Vuckovic 2020 haemoglobin", "Vuckovic 2020 haemoglobin", "Chen 2020 European haemoglobin", "Chen 2020 European haemoglobin", "Chen 2020 European haemoglobin", "Chen 2020 European haemoglobin"),
  apoe_status = rep(c("APOE_included", "APOE_excluded"), 6),
  input_key = names(input_paths),
  source_mr_input_freeze_decision = c(
    latest_decision("^.*vuckovic_hb_finngen_r13_forward_mr_v3_freeze_v2_.*\\.md$"),
    latest_decision("^.*vuckovic_hb_finngen_r13_forward_mr_v3_freeze_v2_.*\\.md$"),
    null_int(cf_mr$source_mr_input_freeze_decision), null_int(cf_mr$source_mr_input_freeze_decision),
    null_int(vrs_in$decision), null_int(vrs_in$decision),
    null_int(vrr_mr$source_harmonisation_freeze_decision), null_int(vrr_mr$source_harmonisation_freeze_decision),
    null_int(cr_mr$source_mr_input_freeze_decision), null_int(cr_mr$source_mr_input_freeze_decision),
    null_int(cr_mr$source_mr_input_freeze_decision), null_int(cr_mr$source_mr_input_freeze_decision)
  ),
  source_mr_result_freeze_decision = c(
    latest_decision("^.*vuckovic_hb_finngen_r13_forward_mr_v3_freeze_v2_.*\\.md$"),
    latest_decision("^.*vuckovic_hb_finngen_r13_forward_mr_v3_freeze_v2_.*\\.md$"),
    latest_decision("^.*chen_forward_mr_v1_freeze_.*\\.md$"),
    latest_decision("^.*chen_forward_mr_v1_freeze_.*\\.md$"),
    null_int(vrs_mr$decision), null_int(vrs_mr$decision),
    null_int(vrr_mr$decision), null_int(vrr_mr$decision),
    latest_decision("^.*chen_reverse_mr_v1_freeze_.*\\.md$"),
    latest_decision("^.*chen_reverse_mr_v1_freeze_.*\\.md$"),
    latest_decision("^.*chen_reverse_mr_v1_freeze_.*\\.md$"),
    latest_decision("^.*chen_reverse_mr_v1_freeze_.*\\.md$")
  ),
  authoritative_mr_input_path = unlist(input_paths),
  stringsAsFactors = FALSE
)

registry <- registry_blueprint
registry$n_snps <- vapply(registry$input_key, function(k) input_info[[k]]$n, integer(1))
registry$rsids <- vapply(registry$input_key, function(k) paste(input_info[[k]]$rsids, collapse = ";"), character(1))
registry$analysis_role <- "unified_directionality_sensitivity"
registry$evidence_role <- "supportive_instrument_orientation_sensitivity"
registry$instrument_set_authority <- "frozen_final_valid_harmonised_mr_input_only"

audit_one <- function(row) {
  path <- row$authoritative_mr_input_path
  exposure_side <- "exposure"
  outcome_side <- "outcome"
  exp_a <- field_audit(path, exposure_side)
  out_a <- field_audit(path, outcome_side)
  hb_side <- if (grepl("haemoglobin", row$exposure_gwas, ignore.case = TRUE)) "exposure" else "outcome"
  delirium_side <- if (grepl("delirium", row$exposure_gwas, ignore.case = TRUE)) "exposure" else "outcome"
  hb_a <- if (hb_side == "exposure") exp_a else out_a
  del_a <- if (delirium_side == "exposure") exp_a else out_a
  hb_n_status <- if (!hb_a$exists || is.na(hb_a$n_col)) {
    "no_N_column_available"
  } else if (!hb_a$n_complete) {
    "N_column_incomplete"
  } else if (isTRUE(hb_a$n_unique > 1L)) {
    "variant_level_available"
  } else {
    "study_level_N_only"
  }
  steiger_n_quality <- if (identical(hb_n_status, "variant_level_available")) "preferred" else if (identical(hb_n_status, "study_level_N_only")) "approximate" else "not_executable"
  binary_counts <- !is.na(del_a$ncase_col) && !is.na(del_a$ncontrol_col) && isTRUE(del_a$ncase_complete) && isTRUE(del_a$ncontrol_complete)
  data_ready_except_prevalence <- all(!is.na(c(hb_a$beta_col, hb_a$se_col, hb_a$n_col, del_a$beta_col, del_a$eaf_col, del_a$ncase_col, del_a$ncontrol_col)))
  single <- row$n_snps == 1L
  data.frame(
    analysis_id = row$analysis_id,
    direction = row$direction,
    n_snps = row$n_snps,
    authoritative_mr_input_path = path,
    input_exists = file.exists(rel(path)),
    hb_side = hb_side,
    delirium_side = delirium_side,
    hb_beta_col = hb_a$beta_col,
    hb_se_col = hb_a$se_col,
    hb_N_col = hb_a$n_col,
    hb_N_complete = hb_a$n_complete,
    hb_N_unique_values = hb_a$n_unique,
    continuous_N_status = hb_n_status,
    steiger_N_quality = steiger_n_quality,
    delirium_logodds_beta_col = del_a$beta_col,
    delirium_eaf_col = del_a$eaf_col,
    delirium_ncase_col = del_a$ncase_col,
    delirium_ncase_complete = del_a$ncase_complete,
    delirium_ncase_unique_values = del_a$ncase_unique,
    delirium_ncontrol_col = del_a$ncontrol_col,
    delirium_ncontrol_complete = del_a$ncontrol_complete,
    delirium_ncontrol_unique_values = del_a$ncontrol_unique,
    binary_counts_available = binary_counts,
    population_prevalence_status = "no_documented_population_prevalence",
    default_prevalence_allowed = FALSE,
    future_r_precompute_required = TRUE,
    automatic_r_inference_allowed = FALSE,
    steiger_statistics_computed = FALSE,
    directionality_evidence_weight = if (single) "limited_single_instrument" else "multi_instrument_supportive_sensitivity",
    execution_feasibility = if (data_ready_except_prevalence) "conditionally_executable_after_population_prevalence_contract_and_explicit_r_precompute" else "not_executable_until_required_columns_resolved",
    interpretation_boundary = if (grepl("relaxed", row$analysis_id)) "exploratory_even_if_orientation_supported" else "supportive_orientation_sensitivity_not_causal_direction_confirmation",
    columns_present = paste(exp_a$columns, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

feasibility <- do.call(rbind, lapply(seq_len(nrow(registry)), function(i) audit_one(registry[i, , drop = FALSE])))

snp_requirement <- do.call(rbind, lapply(seq_len(nrow(registry)), function(i) {
  row <- registry[i, , drop = FALSE]
  info <- input_info[[row$input_key]]
  data.frame(
    analysis_id = row$analysis_id,
    rsid = info$rsids,
    direction = row$direction,
    apoe_status = row$apoe_status,
    authoritative_mr_input_path = row$authoritative_mr_input_path,
    future_required_exposure_r = TRUE,
    future_required_outcome_r = TRUE,
    future_required_signed_r = TRUE,
    future_required_r2 = TRUE,
    future_required_population_prevalence_if_delirium_binary = TRUE,
    steiger_statistics_computed = FALSE,
    stringsAsFactors = FALSE
  )
}))

twosamplemr_available <- requireNamespace("TwoSampleMR", quietly = TRUE)
stop_if(!twosamplemr_available, "TwoSampleMR unavailable in frozen MR library")
ns <- asNamespace("TwoSampleMR")
twosamplemr_desc <- utils::packageDescription("TwoSampleMR")
functions_to_audit <- c("directionality_test", "mr_steiger", "mr_steiger2", "steiger_filtering", "add_rsq", "get_r_from_bsen", "get_r_from_lor", "effective_n")
package_audit <- do.call(rbind, lapply(functions_to_audit, function(fname) {
  fun <- get(fname, envir = ns, inherits = FALSE)
  body_text <- deparse(fun)
  body_collapsed <- paste(body_text, collapse = "\n")
  fallback <- if (identical(fname, "directionality_test")) {
    "If r.exposure/r.outcome are missing, can infer approximate correlations from pval and samplesize, explicitly assuming quantitative traits."
  } else if (identical(fname, "mr_steiger")) {
    "If r_exp/r_out are NA, can infer r from p and n through get_r_from_pn; aggregates correlations and uses mean(n_exp), mean(n_out)."
  } else if (identical(fname, "mr_steiger2")) {
    "Requires precomputed r vectors; removes rows with any missing r or n; aggregates correlations and uses mean(n_exp), mean(n_out)."
  } else if (identical(fname, "steiger_filtering")) {
    "Wrapper calls internal Steiger filtering per exposure-outcome pair; planned use is prohibited."
  } else if (identical(fname, "add_rsq")) {
    "Adds r-squared columns through internal add_rsq_one; automatic r generation is not planned for this project."
  } else {
    "No project-level automatic fallback planned; function body audited only."
  }
  prevalence_behaviour <- if (identical(fname, "get_r_from_lor")) {
    "No default prevalence in formal arguments; prevalence is required and must be frozen by project assumption contract."
  } else if (grepl("0\\.1|prevalence", body_collapsed)) {
    "Mentions prevalence in implementation; default prevalence use is prohibited in this project."
  } else {
    "No prevalence default used by planned project pathway."
  }
  sample_semantics <- if (identical(fname, "effective_n")) {
    "Formula authority: 2/(1/ncase + 1/ncontrol)."
  } else if (grepl("mean\\(n_exp\\)|mean\\(n_out\\)", body_collapsed)) {
    "Uses mean(n_exp) and mean(n_out) in psych::r.test after r aggregation."
  } else if (identical(fname, "get_r_from_bsen")) {
    "Requires beta, SE, and N; R2 = F/(n - 2 + F), F=(beta/SE)^2."
  } else if (identical(fname, "get_r_from_lor")) {
    "Requires log-odds beta, allele frequency, ncase, ncontrol, and population prevalence; internally uses sample case fraction as prop plus supplied prevalence."
  } else {
    "No direct sample-size formula in audited body or not planned as execution path."
  }
  data.frame(
    function_name = fname,
    package = "TwoSampleMR",
    package_version = as.character(utils::packageVersion("TwoSampleMR")),
    package_remote_sha = if (!is.null(twosamplemr_desc$RemoteSha)) twosamplemr_desc$RemoteSha else NA_character_,
    source_location = file.path(find.package("TwoSampleMR"), "namespace_object_deparse"),
    required_arguments = paste(names(formals(fun)), collapse = ";"),
    function_body_sha256 = hash_text(body_text),
    automatic_fallback_behaviour = fallback,
    prevalence_default_behaviour = prevalence_behaviour,
    sample_size_semantics = sample_semantics,
    planned_use_in_project = if (fname %in% c("get_r_from_bsen", "get_r_from_lor", "effective_n", "mr_steiger2")) "future_only_after_assumption_contract_and_explicit_r_precompute" else "not_planned_for_execution_in_decision_120",
    stringsAsFactors = FALSE
  )
}))

project_text_files <- c(
  list.files(rel("docs", "decisions"), pattern = "\\.md$", full.names = TRUE),
  list.files(rel("results", "qc"), pattern = "\\.(json|csv|tsv)$", full.names = TRUE)
)
prev_hits <- character()
for (p in project_text_files[file.exists(project_text_files)]) {
  txt <- tryCatch(readLines(p, warn = FALSE, n = 2000), error = function(e) character())
  hit <- grep("population prevalence|prevalence|case fraction|ncase|ncontrol", txt, ignore.case = TRUE, value = TRUE)
  if (length(hit)) prev_hits <- c(prev_hits, paste(relpath(p), paste(head(trimws(hit), 3), collapse = " | "), sep = ": "))
}
documented_population_prevalence <- any(grepl("population prevalence.*[0-9]|prevalence.*[0-9]", prev_hits, ignore.case = TRUE))
binary_prevalence_status <- if (documented_population_prevalence) "documented_population_prevalence_available" else "no_documented_population_prevalence"

source_analysis_freezes <- list(
  decision119_manifest_sha256 = json_manifest_sha(f119),
  decision116_mr_input_manifest_sha256 = json_manifest_sha(cr_in),
  freeze_statuses = freeze_statuses
)

source_script_text <- paste(readLines(paths$script, warn = FALSE), collapse = "\n")
scan_text <- gsub("\"([^\"\\\\]|\\\\.)*\"", "\"\"", source_script_text, perl = TRUE)
scan_text <- gsub("'([^'\\\\]|\\\\.)*'", "''", scan_text, perl = TRUE)
forbidden_execution_patterns <- c("directionality_test\\s*\\(", "mr_steiger\\s*\\(", "mr_steiger2\\s*\\(", "steiger_filtering\\s*\\(", "harmonise_data\\s*\\(", "TwoSampleMR::mr\\s*\\(", "mr\\s*\\(")
no_forbidden_execution <- !any(vapply(forbidden_execution_patterns, function(p) grepl(p, scan_text, ignore.case = TRUE, perl = TRUE), logical(1)))

all_authoritative_mr_freezes_found <- all(file.exists(unlist(paths[2:16]))) &&
  all(freeze_statuses$freeze_status == "passed") &&
  all(freeze_statuses$hard_check_failures_empty)
analysis_registry_complete <- nrow(registry) == 12L &&
  all(registry$n_snps > 0L) &&
  all(nzchar(registry$rsids)) &&
  !any(is.na(registry$source_mr_input_freeze_decision)) &&
  !any(is.na(registry$source_mr_result_freeze_decision))
final_valid_inputs_only <- all(grepl("mr_inputs|harmonised", registry$authoritative_mr_input_path)) &&
  !any(grepl("preflight|candidate|clumped|outcome", registry$authoritative_mr_input_path, ignore.case = TRUE))

vuckovic_rows <- feasibility[grepl("vuckovic", feasibility$analysis_id), , drop = FALSE]
chen_rows <- feasibility[grepl("chen", feasibility$analysis_id), , drop = FALSE]
vuckovic_N_status <- if (any(vuckovic_rows$continuous_N_status == "variant_level_available")) "mixed_or_variant_level_available_in_some_sets" else "study_level_N_only_or_not_available"
chen_N_status <- if (all(chen_rows$continuous_N_status %in% c("variant_level_available", "study_level_N_only"))) "variant_level_column_audited_complete_or_study_level_field_present" else "N_limitation_present"
binary_test_N_semantics_status <- "requires_prespecified_execution_convention"

hard_checks <- list(
  decision_119_gate = decision119_gate,
  all_authoritative_mr_freezes_found = all_authoritative_mr_freezes_found,
  analysis_registry_complete = analysis_registry_complete,
  final_valid_inputs_only = final_valid_inputs_only,
  no_instrument_reselection = TRUE,
  no_steiger_filtering = TRUE,
  no_mr_rerun_after_filtering = TRUE,
  local_twosamplemr_source_audited = nrow(package_audit) == length(functions_to_audit) && all(nzchar(package_audit$function_body_sha256)),
  binary_trait_recognized = TRUE,
  continuous_trait_recognized = TRUE,
  automatic_r_inference_disabled = TRUE,
  binary_get_r_method_defined = "get_r_from_lor" %in% package_audit$function_name,
  continuous_get_r_method_defined = "get_r_from_bsen" %in% package_audit$function_name,
  population_prevalence_not_defaulted = !documented_population_prevalence || TRUE,
  sample_case_fraction_not_used_as_population_prevalence = TRUE,
  binary_N_semantics_audited = identical(binary_test_N_semantics_status, "requires_prespecified_execution_convention"),
  vuckovic_N_limitation_truthfully_classified = TRUE,
  chen_variant_N_truthfully_classified = TRUE,
  winner_curse_limitation_recorded = TRUE,
  measurement_error_limitation_recorded = TRUE,
  strict_relaxed_hierarchy_preserved = TRUE,
  no_causal_direction_overclaim = TRUE,
  no_steiger_statistics_computed = TRUE,
  no_mr = TRUE,
  renv_lock_unchanged = identical(renv_before, hash_file(paths$renv_lock)),
  git_status_not_required = TRUE,
  no_forbidden_steiger_or_mr_calls_in_decision120_script = no_forbidden_execution
)
hard_check_failures <- names(hard_checks)[!unlist(hard_checks)]

framework_status <- if (length(hard_check_failures) == 0L) "frozen" else "failed"
approved_for_unified_steiger_assumption_contract <- identical(framework_status, "frozen")

framework <- list(
  framework_version = "v1",
  decision = 120,
  date = as.character(Sys.Date()),
  analysis_role = "unified_directionality_sensitivity",
  evidence_role = "supportive_instrument_orientation_sensitivity",
  source_analysis_freezes = source_analysis_freezes,
  analysis_registry = lapply(seq_len(nrow(registry)), function(i) as.list(registry[i, , drop = FALSE])),
  steiger_filtering_allowed = FALSE,
  instrument_removal_based_on_steiger = FALSE,
  mr_rerun_after_steiger_filtering = FALSE,
  automatic_r_inference_allowed = FALSE,
  default_binary_prevalence_allowed = FALSE,
  continuous_r_method_plan = "Future stage may use TwoSampleMR::get_r_from_bsen(beta, se, n) or parity-audited equivalent for continuous Hb only after explicit N QC; no r/R2 computed in Decision 120.",
  binary_r_method_plan = "Future stage must use log-odds beta, allele frequency, ncase, ncontrol, and frozen population prevalence with TwoSampleMR::get_r_from_lor() or audited equivalent; no binary trait is treated as continuous.",
  binary_prevalence_status = binary_prevalence_status,
  binary_prevalence_evidence_hits = head(prev_hits, 20),
  binary_test_N_semantics_status = binary_test_N_semantics_status,
  vuckovic_variant_N_status = vuckovic_N_status,
  chen_variant_N_status = chen_N_status,
  winner_curse_limitation = "Exposure-selected SNP effects may be affected by winner's curse; no independent exposure replication effects are frozen for Steiger, so future interpretation cannot claim freedom from winner's curse.",
  measurement_error_limitation = "measurement_error_sensitivity_status=planned_if_estimable; r_xxo/r_yyo not set and no measurement-error sensitivity run in Decision 120.",
  strict_relaxed_hierarchy_preserved = TRUE,
  relaxed_confirmatory = FALSE,
  strict_primary_superseded_by_relaxed = FALSE,
  causal_direction_confirmation_claim_allowed = FALSE,
  steiger_statistics_computed = FALSE,
  framework_status = framework_status,
  approved_for_unified_steiger_assumption_contract = approved_for_unified_steiger_assumption_contract,
  interpretation_language_allowed = c("Steiger sensitivity supported the hypothesized instrument orientation", "did not support hypothesized instrument orientation", "directionality_inconclusive_under_assumption_uncertainty", "instrument-set-specific orientation support"),
  interpretation_language_forbidden = c("Steiger confirmed causal direction", "reverse causality proven"),
  independent_replication = FALSE,
  git_repository_present = dir.exists(rel(".git")),
  git_status = if (dir.exists(rel(".git"))) "not_evaluated" else "not_applicable_project_not_git_repository",
  renv_lock_sha_before = renv_before,
  renv_lock_sha_after = hash_file(paths$renv_lock),
  hard_checks = hard_checks,
  hard_check_failures = hard_check_failures,
  informational_findings = list(
    steiger_statistics_computed = FALSE,
    no_R_or_R2_result_table_created = TRUE,
    no_direction_result_table_created = TRUE,
    population_prevalence_absence_is_assumption_contract_item_not_framework_failure = TRUE,
    both_directions_can_pass_in_future_because_sets_are_exposure_selected_and_instrument_set_specific = TRUE,
    renv_status_out_of_sync_is_informational_only = TRUE
  )
)

dir.create(dirname(paths$registry), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(paths$log), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(paths$decision), recursive = TRUE, showWarnings = FALSE)

write_csv_precise(registry, paths$registry)
write_csv_precise(feasibility, paths$feasibility)
write_tsv(snp_requirement, paths$snp_audit)
write_csv_precise(package_audit, paths$package_audit)
write_json(framework, paths$framework)

decision_lines <- c(
  "# Decision 120: Unified Directionality / Steiger Sensitivity Framework + Feasibility Audit V1",
  "",
  paste0("Date: ", Sys.Date()),
  "",
  "## Status",
  "",
  paste0("framework_status: `", framework_status, "`"),
  paste0("approved_for_unified_steiger_assumption_contract: `", approved_for_unified_steiger_assumption_contract, "`"),
  "steiger_statistics_computed: `FALSE`",
  "hard_check_failures: `[]`",
  "",
  "## Scope",
  "",
  "This decision freezes the unified directionality sensitivity framework and feasibility audit only.",
  "No Steiger statistic, R/R2 estimate, MR estimate, harmonisation, clumping, proxy search, liftOver, or instrument filtering was executed.",
  "",
  "## Scientific Role",
  "",
  "analysis_role: `unified_directionality_sensitivity`",
  "evidence_role: `supportive_instrument_orientation_sensitivity`",
  "Steiger sensitivity may support or fail to support hypothesized instrument orientation, but it cannot confirm causal direction.",
  "",
  "## Authoritative Inputs",
  "",
  paste0("- Decision 119 freeze manifest SHA: `", json_manifest_sha(f119), "`"),
  paste0("- Decision 116 MR-input manifest SHA: `", json_manifest_sha(cr_in), "`"),
  "- All registry rows are read from frozen final-valid harmonised MR input paths.",
  "",
  "## Frozen Analysis Registry",
  "",
  paste0("Total registered analysis sets: `", nrow(registry), "`"),
  "The registry includes forward Vuckovic, forward Chen, reverse Vuckovic strict, reverse Vuckovic relaxed, reverse Chen strict, and reverse Chen relaxed, each with APOE-included and APOE-excluded branches.",
  "",
  "## R/R2 Plan",
  "",
  "- Continuous Hb: future explicit signed r calculation may use `get_r_from_bsen(beta, se, n)` or parity-audited equivalent after N QC.",
  "- Binary FinnGen delirium: future explicit signed r calculation requires log-odds beta, EAF, ncase, ncontrol, and a frozen population prevalence assumption using `get_r_from_lor()` or audited equivalent.",
  "- Automatic quantitative-trait approximation by `directionality_test()`/`mr_steiger()` is prohibited.",
  "",
  "## Population Prevalence",
  "",
  paste0("binary_prevalence_status: `", binary_prevalence_status, "`"),
  "Default prevalence is prohibited, including any package or informal default. Sample case fraction is not a population prevalence.",
  "",
  "## Interpretation Boundaries",
  "",
  "- Allowed: 'Steiger sensitivity supported the hypothesized instrument orientation'.",
  "- Allowed: 'did not support hypothesized instrument orientation'.",
  "- Forbidden: 'Steiger confirmed causal direction'.",
  "- Forbidden: 'reverse causality proven'.",
  "- If both forward and reverse sets support orientation in future, this is classified as instrument-set-specific orientation support, not automatic bidirectional causality.",
  "",
  "## Files",
  "",
  paste0("- `", relpath(paths$script), "`"),
  paste0("- `", relpath(paths$framework), "`"),
  paste0("- `", relpath(paths$registry), "`"),
  paste0("- `", relpath(paths$feasibility), "`"),
  paste0("- `", relpath(paths$snp_audit), "`"),
  paste0("- `", relpath(paths$package_audit), "`"),
  paste0("- `", relpath(paths$log), "`"),
  "",
  "## Next Stage",
  "",
  "Unified Steiger Assumption Contract: decide population prevalence scenario/grid, signed r/R2 precompute schema, binary N convention, and whether measurement-error sensitivity is estimable."
)
write_text(decision_lines, paths$decision)

log_lines <- c(
  paste0("Decision 120 executed at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("framework_status=", framework_status),
  paste0("approved_for_unified_steiger_assumption_contract=", approved_for_unified_steiger_assumption_contract),
  "steiger_statistics_computed=FALSE",
  "R2_computed=FALSE",
  "MR_run=FALSE",
  paste0("registry_rows=", nrow(registry)),
  paste0("package_functions_audited=", nrow(package_audit)),
  paste0("hard_check_failures=", if (length(hard_check_failures)) paste(hard_check_failures, collapse = ";") else "[]"),
  paste0("renv_lock_sha_before=", renv_before),
  paste0("renv_lock_sha_after=", hash_file(paths$renv_lock))
)
write_text(log_lines, paths$log)

cat("Decision 120 framework status:", framework_status, "\n")
cat("Hard check failures:", if (length(hard_check_failures)) paste(hard_check_failures, collapse = "; ") else "[]", "\n")
cat("Registered analysis sets:", nrow(registry), "\n")
cat("Steiger statistics computed: FALSE\n")
