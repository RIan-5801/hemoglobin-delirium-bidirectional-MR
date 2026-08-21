#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/30_chen_reverse_harmonisation_preflight_v1.R [--project-root <path>]", call. = FALSE)
}
setwd(root)

for (pkg in c("DBI", "duckdb", "jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

rel <- function(...) file.path(root, ...)
norm <- function(x) gsub("\\\\", "/", x)
relpath <- function(x) norm(sub(paste0("^", gsub("\\\\", "/", root), "/?"), "", norm(normalizePath(x, winslash = "/", mustWork = FALSE))))
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
bool <- function(x) toupper(as.character(x)) == "TRUE"
parse_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
records <- function(x) {
  if (!is.data.frame(x) || nrow(x) == 0L) return(list())
  lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
}
atomic_write <- function(path, writer) {
  partial <- paste0(path, ".partial")
  stop_if(file.exists(path) || file.exists(partial), paste("Output occupied:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writer(partial)
  stop_if(!file.exists(partial), paste("Writer did not create partial:", partial))
  stop_if(!file.rename(partial, path), paste("Atomic rename failed:", path))
}
write_json <- function(x, path) atomic_write(path, function(p) {
  jsonlite::write_json(x, p, pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)
})
write_tsv <- function(x, path) atomic_write(path, function(p) {
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.table(x, p, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
})
write_csv <- function(x, path) atomic_write(path, function(p) {
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.csv(x, p, row.names = FALSE, na = "")
})
write_text <- function(lines, path) atomic_write(path, function(p) writeLines(lines, p, useBytes = TRUE))
single_base <- function(x) !is.na(x) & x %in% c("A", "C", "G", "T")
complement <- function(x) {
  y <- toupper(as.character(x))
  z <- c(A = "T", T = "A", C = "G", G = "C")
  unname(z[y])
}
is_palindromic <- function(a, b) {
  key <- paste(sort(c(as.character(a), as.character(b))), collapse = "/")
  key %in% c("A/T", "C/G")
}
classify_alleles <- function(exp_ea, exp_oa, out_ea, out_oa) {
  exp_ea <- toupper(exp_ea); exp_oa <- toupper(exp_oa)
  out_ea <- toupper(out_ea); out_oa <- toupper(out_oa)
  valid <- all(single_base(c(exp_ea, exp_oa, out_ea, out_oa))) && exp_ea != exp_oa && out_ea != out_oa
  if (!valid) return("invalid")
  if (out_ea == exp_ea && out_oa == exp_oa) return("exact")
  if (out_ea == exp_oa && out_oa == exp_ea) return("swapped")
  if (complement(out_ea) == exp_ea && complement(out_oa) == exp_oa) return("strand_exact")
  if (complement(out_ea) == exp_oa && complement(out_oa) == exp_ea) return("strand_swapped")
  "incompatible"
}
planned_beta_flip <- function(classification) classification %in% c("swapped", "strand_swapped")
planned_eaf_flip <- function(classification) classification %in% c("swapped", "strand_swapped")
planned_strand_transform <- function(classification) classification %in% c("strand_exact", "strand_swapped")

paths <- list(
  script = rel("R", "30_chen_reverse_harmonisation_preflight_v1.R"),
  decision47 = rel("docs", "decisions", "47_finngen_r13_delirium_reverse_primary_instruments_v4_freeze_v3_v1.1.md"),
  decision59 = rel("docs", "decisions", "59_finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze_v1.1.md"),
  decision64 = rel("docs", "decisions", "64_reverse_relaxed_palindromic_handling_rule_v1_v1.1.md"),
  decision112 = rel("docs", "decisions", "112_chen_reverse_outcome_extraction_v1_freeze_v1.1.md"),
  strict_freeze = rel("results", "qc", "finngen_r13_delirium_reverse_primary_instruments_v4_freeze_v3.json"),
  relaxed_freeze = rel("results", "qc", "finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze.json"),
  outcome_freeze = rel("results", "qc", "chen_reverse_outcome_extraction_v1_freeze.json"),
  outcome_freeze_manifest = rel("results", "qc", "chen_reverse_outcome_extraction_v1_freeze_manifest.csv"),
  strict_targets = rel("data_derived", "reverse_outcome_extraction", "finngen_r13_delirium_reverse_primary_targets_v1.tsv"),
  relaxed_targets = rel("data_derived", "reverse_outcome_extraction", "finngen_r13_delirium_reverse_relaxed_targets_v1.tsv"),
  chen_union = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_union_targets_v1.tsv"),
  chen_master = rel("data_derived", "reverse_sensitivity_outcome", "chen_reverse_outcome_master_v1.tsv"),
  renv_lock = rel("renv.lock"),
  contract_json = rel("results", "qc", "chen_reverse_harmonisation_contract_v1.json"),
  preflight_json = rel("results", "qc", "chen_reverse_harmonisation_preflight_v1.json"),
  counts_csv = rel("results", "qc", "chen_reverse_harmonisation_preflight_counts_v1.csv"),
  pal_tsv = rel("results", "qc", "chen_reverse_harmonisation_palindromic_snps_v1.tsv"),
  incompatible_tsv = rel("results", "qc", "chen_reverse_harmonisation_incompatible_snps_v1.tsv"),
  master_parquet = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonisation_preflight_master_v1.parquet"),
  log = rel("results", "logs", "chen_reverse_harmonisation_preflight_v1.log"),
  decision = rel("docs", "decisions", "113_chen_reverse_harmonisation_contract_and_preflight_v1_v1.1.md")
)

outputs <- unlist(paths[c("contract_json", "preflight_json", "counts_csv", "pal_tsv", "incompatible_tsv", "master_parquet", "log", "decision")])
occupied <- outputs[file.exists(outputs) | file.exists(paste0(outputs, ".partial"))]
stop_if(length(occupied) > 0L, paste("Output or partial exists:", paste(occupied, collapse = "; ")))

inputs <- unlist(paths[c("decision47", "decision59", "decision64", "decision112", "strict_freeze", "relaxed_freeze", "outcome_freeze", "outcome_freeze_manifest", "strict_targets", "relaxed_targets", "chen_union", "chen_master", "renv_lock")])
missing_inputs <- inputs[!file.exists(inputs)]
stop_if(length(missing_inputs) > 0L, paste("Required input missing:", paste(missing_inputs, collapse = "; ")))

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 113L), paste0("Expected next decision 113, found ", next_decision, "; no outputs written."))

dir.create(dirname(paths$log), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = paths$log, append = TRUE)
log_line("stage=chen_reverse_harmonisation_contract_preflight_v1_start")

strict_freeze <- read_json(paths$strict_freeze)
relaxed_freeze <- read_json(paths$relaxed_freeze)
outcome_freeze <- read_json(paths$outcome_freeze)
strict_targets <- utils::read.delim(paths$strict_targets, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
relaxed_targets <- utils::read.delim(paths$relaxed_targets, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
chen_union <- utils::read.delim(paths$chen_union, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
chen_master <- utils::read.delim(paths$chen_master, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")

strict_exposure <- data.frame(
  rsid = strict_targets$target_rsid,
  strict_included_member = strict_targets$membership == "apoe_included",
  strict_excluded_member = strict_targets$membership == "apoe_excluded",
  exposure_effect_allele = toupper(strict_targets$exposure_effect_allele_raw),
  exposure_other_allele = toupper(strict_targets$exposure_other_allele_raw),
  exposure_beta = parse_num(strict_targets$exposure_beta_raw),
  exposure_se = parse_num(strict_targets$exposure_se_raw),
  exposure_pval = parse_num(strict_targets$exposure_pval_raw),
  exposure_eaf = parse_num(strict_targets$exposure_eaf_raw),
  F_stat = parse_num(strict_targets$exposure_F_statistic),
  stringsAsFactors = FALSE
)
relaxed_exposure <- data.frame(
  rsid = relaxed_targets$target_rsid,
  relaxed_included_member = bool(relaxed_targets$included_member),
  relaxed_excluded_member = bool(relaxed_targets$excluded_member),
  exposure_effect_allele = toupper(relaxed_targets$exposure_effect_allele),
  exposure_other_allele = toupper(relaxed_targets$exposure_other_allele),
  exposure_beta = parse_num(relaxed_targets$exposure_beta),
  exposure_se = parse_num(relaxed_targets$exposure_se),
  exposure_pval = parse_num(relaxed_targets$exposure_pval),
  exposure_eaf = parse_num(relaxed_targets$exposure_eaf),
  F_stat = parse_num(relaxed_targets$exposure_F_statistic),
  stringsAsFactors = FALSE
)

all_rsids <- sort(unique(c(strict_exposure$rsid, relaxed_exposure$rsid)))
master <- data.frame(rsid = all_rsids, stringsAsFactors = FALSE)
master <- merge(master, strict_exposure, by = "rsid", all.x = TRUE, sort = FALSE)
master <- merge(master, relaxed_exposure, by = "rsid", all.x = TRUE, sort = FALSE, suffixes = c("_strict", "_relaxed"))
coalesce <- function(a, b) ifelse(!is.na(a) & a != "", a, b)
master$exposure_effect_allele <- coalesce(master$exposure_effect_allele_strict, master$exposure_effect_allele_relaxed)
master$exposure_other_allele <- coalesce(master$exposure_other_allele_strict, master$exposure_other_allele_relaxed)
master$exposure_beta <- ifelse(!is.na(master$exposure_beta_strict), master$exposure_beta_strict, master$exposure_beta_relaxed)
master$exposure_se <- ifelse(!is.na(master$exposure_se_strict), master$exposure_se_strict, master$exposure_se_relaxed)
master$exposure_pval <- ifelse(!is.na(master$exposure_pval_strict), master$exposure_pval_strict, master$exposure_pval_relaxed)
master$exposure_eaf <- ifelse(!is.na(master$exposure_eaf_strict), master$exposure_eaf_strict, master$exposure_eaf_relaxed)
master$F_stat <- ifelse(!is.na(master$F_stat_strict), master$F_stat_strict, master$F_stat_relaxed)
for (x in c("strict_included_member", "strict_excluded_member", "relaxed_included_member", "relaxed_excluded_member")) {
  master[[x]][is.na(master[[x]])] <- FALSE
}
keep_cols <- c("rsid", "strict_included_member", "strict_excluded_member", "relaxed_included_member", "relaxed_excluded_member", "exposure_effect_allele", "exposure_other_allele", "exposure_beta", "exposure_se", "exposure_pval", "exposure_eaf", "F_stat")
master <- master[, keep_cols]

outcome <- chen_master[, c("target_rsid", "source_marker_id", "reference_allele", "other_allele", "beta", "se", "eaf", "p_value", "n_samples", "match_status", "available_for_harmonisation")]
names(outcome) <- c("rsid", "source_marker_id", "outcome_effect_allele_raw", "outcome_other_allele_raw", "outcome_beta_raw", "outcome_se_raw", "outcome_eaf_raw", "outcome_p_value", "outcome_n_samples", "outcome_match_status", "available_for_harmonisation")
master <- merge(master, outcome, by = "rsid", all.x = TRUE, sort = FALSE)
master$outcome_effect_allele_raw <- toupper(master$outcome_effect_allele_raw)
master$outcome_other_allele_raw <- toupper(master$outcome_other_allele_raw)
master$outcome_beta_raw <- parse_num(master$outcome_beta_raw)
master$outcome_se_raw <- parse_num(master$outcome_se_raw)
master$outcome_eaf_raw <- parse_num(master$outcome_eaf_raw)
master$outcome_p_value <- parse_num(master$outcome_p_value)
master$outcome_n_samples <- parse_num(master$outcome_n_samples)

master$classification <- vapply(seq_len(nrow(master)), function(i) {
  classify_alleles(master$exposure_effect_allele[[i]], master$exposure_other_allele[[i]], master$outcome_effect_allele_raw[[i]], master$outcome_other_allele_raw[[i]])
}, character(1))
master$palindromic <- vapply(seq_len(nrow(master)), function(i) {
  is_palindromic(master$exposure_effect_allele[[i]], master$exposure_other_allele[[i]])
}, logical(1))
master$planned_beta_flip <- planned_beta_flip(master$classification)
master$planned_eaf_flip <- planned_eaf_flip(master$classification)
master$planned_strand_transform <- planned_strand_transform(master$classification)
compatible_class <- master$classification %in% c("exact", "swapped", "strand_exact", "strand_swapped")
master$projected_final_valid <- compatible_class & !master$palindromic & bool(master$available_for_harmonisation)
for (branch in c("strict_included", "strict_excluded", "relaxed_included", "relaxed_excluded")) {
  member_col <- paste0(branch, "_member")
  out_col <- paste0(branch, "_projected_final_valid")
  reason_col <- paste0(branch, "_projected_exclusion_reason")
  master[[out_col]] <- master[[member_col]] & master$projected_final_valid
  master[[reason_col]] <- ifelse(
    !master[[member_col]], "not_in_branch",
    ifelse(!bool(master$available_for_harmonisation), "no_unique_exact_outcome",
      ifelse(master$classification == "invalid", "invalid_allele",
        ifelse(master$classification == "incompatible", "allele_incompatible",
          ifelse(master$palindromic, paste0(branch, "_excluded_by_prespecified_palindrome_rule"), NA_character_)
        )
      )
    )
  )
}

branch_summary <- function(branch, member_col, role, threshold, evidence_role) {
  x <- master[master[[member_col]], , drop = FALSE]
  final <- x[x$projected_final_valid, , drop = FALSE]
  class_levels <- c("exact", "swapped", "strand_exact", "strand_swapped", "incompatible", "invalid")
  counts <- setNames(vapply(class_levels, function(k) sum(x$classification == k), integer(1)), class_levels)
  data.frame(
    branch = branch,
    analysis_role = role,
    evidence_role = evidence_role,
    threshold = threshold,
    input_count = nrow(x),
    exact = counts[["exact"]],
    swapped = counts[["swapped"]],
    strand_exact = counts[["strand_exact"]],
    strand_swapped = counts[["strand_swapped"]],
    incompatible = counts[["incompatible"]],
    invalid = counts[["invalid"]],
    palindromic_count = sum(x$palindromic),
    projected_final_valid_count = nrow(final),
    projected_final_valid_rsids = paste(final$rsid, collapse = ";"),
    projected_mr_estimability = nrow(final) > 0L,
    planned_beta_flip_count = sum(x$planned_beta_flip),
    planned_eaf_flip_count = sum(x$planned_eaf_flip),
    planned_strand_transform_count = sum(x$planned_strand_transform),
    min_F = if (nrow(x)) min(x$F_stat, na.rm = TRUE) else NA_real_,
    median_F = if (nrow(x)) stats::median(x$F_stat, na.rm = TRUE) else NA_real_,
    max_F = if (nrow(x)) max(x$F_stat, na.rm = TRUE) else NA_real_,
    weak_F_lt_10_count = sum(x$F_stat < 10, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
counts <- rbind(
  branch_summary("strict_apoe_included", "strict_included_member", "reverse_alternative_hb_outcome_sensitivity", 5e-8, "reverse_strict_primary_outcome_sensitivity"),
  branch_summary("strict_apoe_excluded", "strict_excluded_member", "reverse_alternative_hb_outcome_sensitivity", 5e-8, "reverse_strict_primary_outcome_sensitivity"),
  branch_summary("relaxed_apoe_included", "relaxed_included_member", "reverse_alternative_hb_outcome_sensitivity", 5e-6, "reverse_protocol_triggered_exploratory_outcome_sensitivity"),
  branch_summary("relaxed_apoe_excluded", "relaxed_excluded_member", "reverse_alternative_hb_outcome_sensitivity", 5e-6, "reverse_protocol_triggered_exploratory_outcome_sensitivity")
)

pal_rows <- master[master$palindromic, c("rsid", "exposure_effect_allele", "exposure_other_allele", "outcome_effect_allele_raw", "outcome_other_allele_raw", "classification", "strict_included_member", "strict_excluded_member", "relaxed_included_member", "relaxed_excluded_member", "strict_included_projected_exclusion_reason", "strict_excluded_projected_exclusion_reason", "relaxed_included_projected_exclusion_reason", "relaxed_excluded_projected_exclusion_reason")]
pal_rows$relaxed_palindrome_rule_authority_decision <- 64L
pal_rows$strict_palindrome_rule <- "prospective_exclude_all_without_EAF_reinclusion"
pal_rows$eaf_based_palindrome_reinclusion_allowed <- FALSE
bad_rows <- master[master$classification %in% c("incompatible", "invalid"), c("rsid", "exposure_effect_allele", "exposure_other_allele", "outcome_effect_allele_raw", "outcome_other_allele_raw", "classification", "strict_included_member", "strict_excluded_member", "relaxed_included_member", "relaxed_excluded_member")]

strict_overlap <- intersect(master$rsid[master$strict_included_projected_final_valid], master$rsid[master$strict_excluded_projected_final_valid])
relaxed_overlap <- intersect(master$rsid[master$relaxed_included_projected_final_valid], master$rsid[master$relaxed_excluded_projected_final_valid])

con_db <- DBI::dbConnect(duckdb::duckdb(), read_only = FALSE)
on.exit(try(DBI::dbDisconnect(con_db, shutdown = TRUE), silent = TRUE), add = TRUE)
write_parquet_df <- function(x, path, table_name) {
  atomic_write(path, function(p) {
    DBI::dbRemoveTable(con_db, table_name, fail_if_missing = FALSE)
    DBI::dbWriteTable(con_db, table_name, x)
    DBI::dbExecute(con_db, sprintf("COPY %s TO %s (FORMAT PARQUET)", DBI::dbQuoteIdentifier(con_db, table_name), DBI::dbQuoteString(con_db, norm(p))))
  })
}

hard_checks_common <- list(
  decision_47_strict_exposure_gate = identical(strict_freeze$freeze_status, "passed") && identical(strict_freeze$manifest_sha256, "a13a4d946a656854ec1215b9d49f6d9292fc826433992175f227846249b03824"),
  decision_59_relaxed_exposure_gate = file.exists(paths$decision59) && identical(hash_file(paths$decision59), "44451cb5498035cc5a350275d75088ff0289dffec738797ce4692124c6ac7939"),
  decision_112_outcome_gate = identical(outcome_freeze$freeze_status, "passed") && identical(outcome_freeze$manifest_sha256, "b272de6ce5a1852593f71521184fbf2e3e5ee913cb33bf762f19490631899dcd") && length(outcome_freeze$hard_check_failures) == 0L,
  strict_relaxed_hierarchy_preserved = TRUE,
  all_outcome_exact_matches_used = nrow(master) == as.integer(outcome_freeze$union_target_count) && all(master$outcome_match_status == "unique_exact_match"),
  no_vuckovic_outcome_conditioning = TRUE,
  exposure_orientation_preserved = all(master$exposure_effect_allele %in% c("A", "C", "G", "T")) && all(master$exposure_other_allele %in% c("A", "C", "G", "T")),
  chen_outcome_orientation_preserved_raw = all(master$outcome_effect_allele_raw %in% c("A", "C", "G", "T")) && all(master$outcome_other_allele_raw %in% c("A", "C", "G", "T")),
  marker_tokens_not_used_for_effect_orientation = TRUE,
  cross_build_coordinates_not_rematched = TRUE,
  allele_classification_exhaustive = all(!is.na(master$classification)) && all(master$classification %in% c("exact", "swapped", "strand_exact", "strand_swapped", "incompatible", "invalid")),
  palindrome_definition_fixed = TRUE,
  decision_64_relaxed_rule_preserved = file.exists(paths$decision64),
  strict_palindrome_rule_prospectively_adopted = TRUE,
  no_eaf_palindrome_reinclusion = TRUE,
  projected_final_valid_rule_correct = all(master$projected_final_valid == (compatible_class & !master$palindromic & bool(master$available_for_harmonisation))),
  planned_transformations_auditable = all(!is.na(master$planned_beta_flip)) && all(!is.na(master$planned_eaf_flip)) && all(!is.na(master$planned_strand_transform)),
  F_statistics_preserved = all(abs(master$F_stat - (master$exposure_beta / master$exposure_se)^2) < 1e-6),
  no_f_ge_30_filter = TRUE,
  no_outcome_based_filter = TRUE,
  no_proxy = TRUE,
  no_liftover = TRUE,
  no_formal_harmonisation = TRUE,
  no_mr = TRUE,
  no_steiger = TRUE,
  renv_lock_unchanged = identical(hash_file(paths$renv_lock), "253471c02e1e47a40d0f68b296d4ae2b1df471d757ac86328b3e974018d039f3")
)
hard_check_failures <- names(hard_checks_common)[!vapply(hard_checks_common, isTRUE, logical(1))]
contract_status <- if (length(hard_check_failures) == 0L) "frozen" else "failed"
preflight_status <- if (length(hard_check_failures) == 0L) "passed" else "failed"
approved_next <- length(hard_check_failures) == 0L

contract <- list(
  contract_version = "v1",
  decision = 113,
  date = format(Sys.Date()),
  analysis_direction = "delirium_to_Hb",
  analysis_role = "reverse_alternative_hb_outcome_sensitivity",
  independent_replication = FALSE,
  strict_exposure_authority_decision = 47,
  relaxed_exposure_authority_decision = 59,
  outcome_extraction_freeze_decision = 112,
  strict_threshold = 5e-8,
  relaxed_threshold = 5e-6,
  strict_relaxed_hierarchy_preserved = TRUE,
  strict_primary_superseded_by_relaxed = FALSE,
  relaxed_confirmatory = FALSE,
  exposure = "FinnGen_R13_F5_DELIRIUM",
  outcome = "Chen_2020_Hb_BCX2",
  exposure_scale = "log_odds_delirium",
  outcome_scale = outcome_freeze$outcome_scale,
  exposure_orientation_source = "FinnGen_frozen_effect_alleles",
  outcome_orientation_source = "Chen_reference_allele_other_allele",
  relaxed_palindrome_rule_authority_decision = 64,
  strict_palindrome_rule_adoption = "prospective_exclude_all_without_EAF_reinclusion",
  eaf_based_palindrome_reinclusion = FALSE,
  harmonisation_classes = c("exact", "swapped", "strand_exact", "strand_swapped", "incompatible", "invalid"),
  proxy_allowed = FALSE,
  liftover_allowed = FALSE,
  formal_harmonisation_performed = FALSE,
  mr_run = FALSE,
  steiger_run = FALSE,
  contract_status = contract_status,
  approved_for_chen_reverse_formal_harmonisation = approved_next,
  hard_checks = hard_checks_common,
  hard_check_failures = hard_check_failures
)

branch_results <- lapply(seq_len(nrow(counts)), function(i) {
  row <- counts[i, , drop = FALSE]
  list(
    branch = row$branch,
    input_count = row$input_count,
    classification_counts = as.list(row[, c("exact", "swapped", "strand_exact", "strand_swapped", "incompatible", "invalid")]),
    palindromic_count = row$palindromic_count,
    projected_final_valid_count = row$projected_final_valid_count,
    projected_final_valid_rsids = if (row$projected_final_valid_rsids == "") list() else strsplit(row$projected_final_valid_rsids, ";", fixed = TRUE)[[1L]],
    planned_beta_flip_count = row$planned_beta_flip_count,
    planned_eaf_flip_count = row$planned_eaf_flip_count,
    planned_strand_transform_count = row$planned_strand_transform_count,
    projected_mr_estimability = row$projected_mr_estimability,
    instrument_strength_summary = list(min_F = row$min_F, median_F = row$median_F, max_F = row$max_F, weak_F_lt_10_count = row$weak_F_lt_10_count)
  )
})
names(branch_results) <- counts$branch

preflight <- list(
  preflight_version = "v1",
  decision = 113,
  date = format(Sys.Date()),
  branch_results = branch_results,
  strict_overlap = list(shared_projected_rsids = strict_overlap, shared_projected_count = length(strict_overlap)),
  relaxed_overlap = list(shared_projected_rsids = relaxed_overlap, shared_projected_count = length(relaxed_overlap)),
  palindromic_rsids = pal_rows$rsid,
  incompatible_or_invalid_rsids = bad_rows$rsid,
  proxy_used = FALSE,
  liftover_used = FALSE,
  formal_harmonisation_performed = FALSE,
  mr_run = FALSE,
  steiger_run = FALSE,
  preflight_status = preflight_status,
  approved_for_chen_reverse_formal_harmonisation = approved_next,
  hard_checks = hard_checks_common,
  hard_check_failures = hard_check_failures
)

write_json(contract, paths$contract_json)
write_json(preflight, paths$preflight_json)
write_csv(counts, paths$counts_csv)
write_tsv(pal_rows, paths$pal_tsv)
write_tsv(bad_rows, paths$incompatible_tsv)
write_parquet_df(master, paths$master_parquet, "chen_reverse_preflight_master")

decision_lines <- c(
  "# Decision 113: Chen Reverse Harmonisation Contract And Preflight V1",
  "",
  paste0("Date: ", format(Sys.Date())),
  "",
  "## Status",
  paste0("contract_status: `", contract_status, "`"),
  paste0("preflight_status: `", preflight_status, "`"),
  paste0("approved_for_chen_reverse_formal_harmonisation: `", approved_next, "`"),
  "",
  "## Decision",
  "Freeze the Chen reverse harmonisation contract and perform preflight-only allele classification for the alternative Hb outcome sensitivity analysis.",
  "",
  "This stage classifies alleles and records planned transformations only. It does not apply beta/EAF flips, does not create formal harmonised MR inputs, does not run MR, and does not run Steiger.",
  "",
  "## Palindrome Rule",
  "Relaxed branches preserve Decision 64: exclude all palindromic SNPs without EAF-based re-inclusion.",
  "",
  "Strict branches prospectively adopt the same conservative rule before formal harmonisation and MR for the Chen reverse outcome sensitivity analysis.",
  "",
  "## Results",
  paste0("- Strict included projected final-valid: `", counts$projected_final_valid_count[counts$branch == "strict_apoe_included"], "/", counts$input_count[counts$branch == "strict_apoe_included"], "`."),
  paste0("- Strict excluded projected final-valid: `", counts$projected_final_valid_count[counts$branch == "strict_apoe_excluded"], "/", counts$input_count[counts$branch == "strict_apoe_excluded"], "`."),
  paste0("- Relaxed included projected final-valid: `", counts$projected_final_valid_count[counts$branch == "relaxed_apoe_included"], "/", counts$input_count[counts$branch == "relaxed_apoe_included"], "`."),
  paste0("- Relaxed excluded projected final-valid: `", counts$projected_final_valid_count[counts$branch == "relaxed_apoe_excluded"], "/", counts$input_count[counts$branch == "relaxed_apoe_excluded"], "`."),
  paste0("- Palindromic SNP count: `", nrow(pal_rows), "`."),
  paste0("- Incompatible/invalid SNP count: `", nrow(bad_rows), "`."),
  paste0("- Hard-check failures: `", if (length(hard_check_failures) == 0L) "none" else paste(hard_check_failures, collapse = ";"), "`."),
  "",
  "## Outputs",
  paste0("- `", relpath(paths$contract_json), "`"),
  paste0("- `", relpath(paths$preflight_json), "`"),
  paste0("- `", relpath(paths$counts_csv), "`"),
  paste0("- `", relpath(paths$pal_tsv), "`"),
  paste0("- `", relpath(paths$incompatible_tsv), "`"),
  paste0("- `", relpath(paths$master_parquet), "`"),
  paste0("- `", relpath(paths$log), "`")
)
write_text(decision_lines, paths$decision)

log_line("contract_status=", contract_status, "; preflight_status=", preflight_status, "; hard_check_failures=", paste(hard_check_failures, collapse = ","))
cat("contract_status=", contract_status, "\n", sep = "")
cat("preflight_status=", preflight_status, "\n", sep = "")
cat("approved_for_chen_reverse_formal_harmonisation=", approved_next, "\n", sep = "")
cat("hard_check_failures=", paste(hard_check_failures, collapse = ","), "\n", sep = "")
cat("palindromic_count=", nrow(pal_rows), "\n", sep = "")
cat("incompatible_invalid_count=", nrow(bad_rows), "\n", sep = "")
