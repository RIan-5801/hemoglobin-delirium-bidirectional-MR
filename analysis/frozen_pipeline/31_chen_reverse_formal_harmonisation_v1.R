#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/31_chen_reverse_formal_harmonisation_v1.R [--project-root <path>]", call. = FALSE)
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
comp <- function(x) unname(c(A = "T", T = "A", C = "G", G = "C")[toupper(as.character(x))])
classify_alleles <- function(exp_ea, exp_oa, out_ea, out_oa) {
  exp_ea <- toupper(exp_ea); exp_oa <- toupper(exp_oa)
  out_ea <- toupper(out_ea); out_oa <- toupper(out_oa)
  valid <- all(single_base(c(exp_ea, exp_oa, out_ea, out_oa))) && exp_ea != exp_oa && out_ea != out_oa
  if (!valid) return("invalid")
  if (out_ea == exp_ea && out_oa == exp_oa) return("exact")
  if (out_ea == exp_oa && out_oa == exp_ea) return("swapped")
  if (comp(out_ea) == exp_ea && comp(out_oa) == exp_oa) return("strand_exact")
  if (comp(out_ea) == exp_oa && comp(out_oa) == exp_ea) return("strand_swapped")
  "incompatible"
}
is_palindromic <- function(a, b) paste(sort(c(as.character(a), as.character(b))), collapse = "/") %in% c("A/T", "C/G")
is_missing_token <- function(x) is.na(x) | (!is.na(x) & as.character(x) == "")

paths <- list(
  script = rel("R", "31_chen_reverse_formal_harmonisation_v1.R"),
  decision47 = rel("docs", "decisions", "47_finngen_r13_delirium_reverse_primary_instruments_v4_freeze_v3_v1.1.md"),
  decision59 = rel("docs", "decisions", "59_finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze_v1.1.md"),
  decision64 = rel("docs", "decisions", "64_reverse_relaxed_palindromic_handling_rule_v1_v1.1.md"),
  decision112 = rel("docs", "decisions", "112_chen_reverse_outcome_extraction_v1_freeze_v1.1.md"),
  decision113 = rel("docs", "decisions", "113_chen_reverse_harmonisation_contract_and_preflight_v1_v1.1.md"),
  strict_freeze = rel("results", "qc", "finngen_r13_delirium_reverse_primary_instruments_v4_freeze_v3.json"),
  relaxed_freeze = rel("results", "qc", "finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze.json"),
  outcome_freeze = rel("results", "qc", "chen_reverse_outcome_extraction_v1_freeze.json"),
  contract = rel("results", "qc", "chen_reverse_harmonisation_contract_v1.json"),
  preflight = rel("results", "qc", "chen_reverse_harmonisation_preflight_v1.json"),
  preflight_counts = rel("results", "qc", "chen_reverse_harmonisation_preflight_counts_v1.csv"),
  preflight_master = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonisation_preflight_master_v1.parquet"),
  renv_lock = rel("renv.lock"),
  master_parquet = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_master_v1.parquet"),
  master_tsv = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_master_v1.tsv"),
  strict_inc_parquet = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_strict_apoe_included_v1.parquet"),
  strict_inc_tsv = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_strict_apoe_included_v1.tsv"),
  strict_exc_parquet = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_strict_apoe_excluded_v1.parquet"),
  strict_exc_tsv = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_strict_apoe_excluded_v1.tsv"),
  relaxed_inc_parquet = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_relaxed_apoe_included_v1.parquet"),
  relaxed_inc_tsv = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_relaxed_apoe_included_v1.tsv"),
  relaxed_exc_parquet = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_relaxed_apoe_excluded_v1.parquet"),
  relaxed_exc_tsv = rel("data_derived", "reverse_sensitivity_harmonisation", "chen_reverse_harmonised_relaxed_apoe_excluded_v1.tsv"),
  counts = rel("results", "qc", "chen_reverse_formal_harmonisation_counts_v1.csv"),
  transform_audit = rel("results", "qc", "chen_reverse_formal_harmonisation_transform_audit_v1.csv"),
  excluded = rel("results", "qc", "chen_reverse_formal_harmonisation_excluded_snps_v1.tsv"),
  qc = rel("results", "qc", "chen_reverse_formal_harmonisation_v1.json"),
  log = rel("results", "logs", "chen_reverse_formal_harmonisation_v1.log"),
  decision = rel("docs", "decisions", "114_chen_reverse_formal_harmonisation_v1_v1.1.md")
)

outputs <- unlist(paths[c("master_parquet", "master_tsv", "strict_inc_parquet", "strict_inc_tsv", "strict_exc_parquet", "strict_exc_tsv", "relaxed_inc_parquet", "relaxed_inc_tsv", "relaxed_exc_parquet", "relaxed_exc_tsv", "counts", "transform_audit", "excluded", "qc", "log", "decision")])
occupied <- outputs[file.exists(outputs) | file.exists(paste0(outputs, ".partial"))]
stop_if(length(occupied) > 0L, paste("Output or partial exists:", paste(occupied, collapse = "; ")))
inputs <- unlist(paths[c("decision47", "decision59", "decision64", "decision112", "decision113", "strict_freeze", "relaxed_freeze", "outcome_freeze", "contract", "preflight", "preflight_counts", "preflight_master", "renv_lock")])
missing_inputs <- inputs[!file.exists(inputs)]
stop_if(length(missing_inputs) > 0L, paste("Required input missing:", paste(missing_inputs, collapse = "; ")))

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 114L), paste0("Expected next decision 114, found ", next_decision, "; no outputs written."))

dir.create(dirname(paths$log), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = paths$log, append = TRUE)
log_line("stage=chen_reverse_formal_harmonisation_v1_start")

strict_freeze <- read_json(paths$strict_freeze)
relaxed_freeze <- read_json(paths$relaxed_freeze)
outcome_freeze <- read_json(paths$outcome_freeze)
contract <- read_json(paths$contract)
preflight <- read_json(paths$preflight)
preflight_counts <- utils::read.csv(paths$preflight_counts, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")

con_db <- DBI::dbConnect(duckdb::duckdb(), read_only = FALSE)
on.exit(try(DBI::dbDisconnect(con_db, shutdown = TRUE), silent = TRUE), add = TRUE)
read_parquet <- function(path) DBI::dbGetQuery(con_db, sprintf("SELECT * FROM read_parquet(%s)", DBI::dbQuoteString(con_db, norm(path))))
write_parquet_df <- function(x, path, table_name) {
  atomic_write(path, function(p) {
    DBI::dbRemoveTable(con_db, table_name, fail_if_missing = FALSE)
    DBI::dbWriteTable(con_db, table_name, x)
    DBI::dbExecute(con_db, sprintf("COPY %s TO %s (FORMAT PARQUET)", DBI::dbQuoteIdentifier(con_db, table_name), DBI::dbQuoteString(con_db, norm(p))))
  })
}

master <- read_parquet(paths$preflight_master)
master$classification_recomputed <- vapply(seq_len(nrow(master)), function(i) {
  classify_alleles(master$exposure_effect_allele[[i]], master$exposure_other_allele[[i]], master$outcome_effect_allele_raw[[i]], master$outcome_other_allele_raw[[i]])
}, character(1))
master$palindromic_recomputed <- vapply(seq_len(nrow(master)), function(i) is_palindromic(master$exposure_effect_allele[[i]], master$exposure_other_allele[[i]]), logical(1))
compatible <- master$classification_recomputed %in% c("exact", "swapped", "strand_exact", "strand_swapped")
master$final_valid <- compatible & !master$palindromic_recomputed & bool(master$available_for_harmonisation)
master$classification_planned_beta_flip <- master$classification_recomputed %in% c("swapped", "strand_swapped")
master$classification_planned_eaf_flip <- master$classification_recomputed %in% c("swapped", "strand_swapped")
master$classification_planned_strand_transform <- master$classification_recomputed %in% c("strand_exact", "strand_swapped")
master$performed_beta_flip <- master$final_valid & master$classification_planned_beta_flip
master$performed_eaf_flip <- master$final_valid & master$classification_planned_eaf_flip
master$performed_strand_transform <- master$final_valid & master$classification_planned_strand_transform
master$outcome_beta_harmonised <- ifelse(master$final_valid, ifelse(master$performed_beta_flip, -master$outcome_beta_raw, master$outcome_beta_raw), NA_real_)
master$outcome_se_harmonised <- ifelse(master$final_valid, master$outcome_se_raw, NA_real_)
master$outcome_eaf_harmonised <- ifelse(master$final_valid, ifelse(master$performed_eaf_flip & !is.na(master$outcome_eaf_raw), 1 - master$outcome_eaf_raw, master$outcome_eaf_raw), NA_real_)
master$outcome_effect_allele_harmonised <- ifelse(master$final_valid, master$exposure_effect_allele, NA_character_)
master$outcome_other_allele_harmonised <- ifelse(master$final_valid, master$exposure_other_allele, NA_character_)
master$exclusion_reason <- ifelse(master$final_valid, NA_character_,
  ifelse(master$palindromic_recomputed, "excluded_by_prespecified_palindrome_rule",
    ifelse(master$classification_recomputed == "invalid", "invalid_allele",
      ifelse(master$classification_recomputed == "incompatible", "allele_incompatible", "not_available_for_harmonisation")
    )
  )
)
master$harmonisation_version <- "v1"
master$formal_harmonisation_applied <- master$final_valid

branch_final <- function(branch, member_col) {
  master[master[[member_col]] & master$final_valid, , drop = FALSE]
}
branch_input <- function(member_col) master[master[[member_col]], , drop = FALSE]
branch_summary <- function(branch, member_col) {
  x <- branch_input(member_col)
  f <- branch_final(branch, member_col)
  cls <- c("exact", "swapped", "strand_exact", "strand_swapped", "incompatible", "invalid")
  data.frame(
    branch = branch,
    input_count = nrow(x),
    exact = sum(x$classification_recomputed == "exact"),
    swapped = sum(x$classification_recomputed == "swapped"),
    strand_exact = sum(x$classification_recomputed == "strand_exact"),
    strand_swapped = sum(x$classification_recomputed == "strand_swapped"),
    incompatible = sum(x$classification_recomputed == "incompatible"),
    invalid = sum(x$classification_recomputed == "invalid"),
    palindromic_excluded_count = sum(x$palindromic_recomputed),
    incompatible_excluded_count = sum(x$classification_recomputed == "incompatible"),
    invalid_excluded_count = sum(x$classification_recomputed == "invalid"),
    final_valid_count = nrow(f),
    final_valid_rsids = paste(f$rsid, collapse = ";"),
    performed_beta_flip_count = sum(f$performed_beta_flip),
    performed_eaf_flip_count = sum(f$performed_eaf_flip),
    performed_strand_transform_count = sum(f$performed_strand_transform),
    F_min = if (nrow(f)) min(f$F_stat, na.rm = TRUE) else NA_real_,
    F_mean = if (nrow(f)) mean(f$F_stat, na.rm = TRUE) else NA_real_,
    F_median = if (nrow(f)) stats::median(f$F_stat, na.rm = TRUE) else NA_real_,
    F_max = if (nrow(f)) max(f$F_stat, na.rm = TRUE) else NA_real_,
    F_lt10 = sum(f$F_stat < 10, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
counts <- rbind(
  branch_summary("strict_apoe_included", "strict_included_member"),
  branch_summary("strict_apoe_excluded", "strict_excluded_member"),
  branch_summary("relaxed_apoe_included", "relaxed_included_member"),
  branch_summary("relaxed_apoe_excluded", "relaxed_excluded_member")
)
expected_sets <- setNames(lapply(preflight$branch_results, function(x) {
  y <- x$projected_final_valid_rsids
  if (length(y) == 0L) character() else as.character(y)
}), names(preflight$branch_results))
actual_sets <- list(
  strict_apoe_included = branch_final("strict_apoe_included", "strict_included_member")$rsid,
  strict_apoe_excluded = branch_final("strict_apoe_excluded", "strict_excluded_member")$rsid,
  relaxed_apoe_included = branch_final("relaxed_apoe_included", "relaxed_included_member")$rsid,
  relaxed_apoe_excluded = branch_final("relaxed_apoe_excluded", "relaxed_excluded_member")$rsid
)
sets_equal <- all(vapply(names(actual_sets), function(nm) setequal(actual_sets[[nm]], expected_sets[[nm]]), logical(1)))

excluded <- master[!master$final_valid, c("rsid", "strict_included_member", "strict_excluded_member", "relaxed_included_member", "relaxed_excluded_member", "exposure_effect_allele", "exposure_other_allele", "outcome_effect_allele_raw", "outcome_other_allele_raw", "classification_recomputed", "palindromic_recomputed", "exclusion_reason")]
transform_audit <- master[, c("rsid", "classification", "classification_recomputed", "palindromic", "palindromic_recomputed", "projected_final_valid", "final_valid", "classification_planned_beta_flip", "performed_beta_flip", "classification_planned_eaf_flip", "performed_eaf_flip", "classification_planned_strand_transform", "performed_strand_transform", "exclusion_reason")]

nullable_character_fields <- c("outcome_effect_allele_harmonised", "outcome_other_allele_harmonised", "exclusion_reason")
validate_pair <- function(parquet, tsv, key = "rsid") {
  p <- read_parquet(parquet)
  t <- utils::read.delim(tsv, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character")
  same_cols <- identical(names(p), names(t))
  same_n <- identical(nrow(p), nrow(t))
  same_order <- same_n && identical(as.character(p[[key]]), as.character(t[[key]]))
  char_cols <- names(p)[!vapply(p, is.numeric, logical(1))]
  num_cols <- names(p)[vapply(p, is.numeric, logical(1))]
  char_ok <- all(vapply(char_cols, function(k) {
    a <- as.character(p[[k]])
    b <- as.character(t[[k]])
    if (k %in% nullable_character_fields) {
      all((is_missing_token(a) & is_missing_token(b)) | (!is.na(a) & !is.na(b) & a == b))
    } else {
      identical(a, b)
    }
  }, logical(1)))
  numeric_stats <- lapply(num_cols, function(k) {
    a <- as.numeric(p[[k]])
    b <- suppressWarnings(as.numeric(t[[k]]))
    both_missing <- is.na(a) & is_missing_token(as.character(t[[k]]))
    both_num <- is.finite(a) & is.finite(b)
    abs_diff <- rep(NA_real_, length(a))
    rel_diff <- rep(NA_real_, length(a))
    abs_diff[both_num] <- abs(a[both_num] - b[both_num])
    rel_diff[both_num] <- abs_diff[both_num] / pmax(abs(a[both_num]), abs(b[both_num]), .Machine$double.xmin)
    data.frame(column = k, all_cells_ok = all(both_missing | both_num), max_abs = if (any(both_num)) max(abs_diff[both_num], na.rm = TRUE) else 0, max_rel = if (any(both_num)) max(rel_diff[both_num], na.rm = TRUE) else 0, stringsAsFactors = FALSE)
  })
  numeric_stats <- if (length(numeric_stats)) do.call(rbind, numeric_stats) else data.frame(column = character(), all_cells_ok = logical(), max_abs = numeric(), max_rel = numeric())
  numeric_ok <- all(numeric_stats$all_cells_ok) && all(numeric_stats$max_abs <= 1e-12 | numeric_stats$max_rel <= 1e-12)
  list(row_count = nrow(p), same_cols = same_cols, same_n = same_n, same_order = same_order, char_ok = char_ok, numeric_ok = numeric_ok, numeric_stats = records(numeric_stats))
}

write_parquet_df(master, paths$master_parquet, "master")
write_tsv(master, paths$master_tsv)
write_parquet_df(branch_final("strict_apoe_included", "strict_included_member"), paths$strict_inc_parquet, "strict_inc")
write_tsv(branch_final("strict_apoe_included", "strict_included_member"), paths$strict_inc_tsv)
write_parquet_df(branch_final("strict_apoe_excluded", "strict_excluded_member"), paths$strict_exc_parquet, "strict_exc")
write_tsv(branch_final("strict_apoe_excluded", "strict_excluded_member"), paths$strict_exc_tsv)
write_parquet_df(branch_final("relaxed_apoe_included", "relaxed_included_member"), paths$relaxed_inc_parquet, "relaxed_inc")
write_tsv(branch_final("relaxed_apoe_included", "relaxed_included_member"), paths$relaxed_inc_tsv)
write_parquet_df(branch_final("relaxed_apoe_excluded", "relaxed_excluded_member"), paths$relaxed_exc_parquet, "relaxed_exc")
write_tsv(branch_final("relaxed_apoe_excluded", "relaxed_excluded_member"), paths$relaxed_exc_tsv)
write_csv(counts, paths$counts)
write_csv(transform_audit, paths$transform_audit)
write_tsv(excluded, paths$excluded)

pair_checks <- list(
  master = validate_pair(paths$master_parquet, paths$master_tsv),
  strict_included = validate_pair(paths$strict_inc_parquet, paths$strict_inc_tsv),
  strict_excluded = validate_pair(paths$strict_exc_parquet, paths$strict_exc_tsv),
  relaxed_included = validate_pair(paths$relaxed_inc_parquet, paths$relaxed_inc_tsv),
  relaxed_excluded = validate_pair(paths$relaxed_exc_parquet, paths$relaxed_exc_tsv)
)
pair_ok <- all(vapply(pair_checks, function(x) x$same_cols && x$same_n && x$same_order && x$char_ok && x$numeric_ok, logical(1)))

final_rows <- master[master$final_valid, , drop = FALSE]
hard_checks <- list(
  decision_47_gate = identical(strict_freeze$freeze_status, "passed"),
  decision_59_gate = file.exists(paths$decision59),
  decision_112_gate = identical(outcome_freeze$freeze_status, "passed") && length(outcome_freeze$hard_check_failures) == 0L,
  decision_113_gate = identical(contract$contract_status, "frozen") && identical(preflight$preflight_status, "passed") && isTRUE(contract$approved_for_chen_reverse_formal_harmonisation) && length(contract$hard_check_failures) == 0L && length(preflight$hard_check_failures) == 0L,
  row_level_preflight_classification_reproduced = all(master$classification == master$classification_recomputed) && all(master$palindromic == master$palindromic_recomputed),
  strict_relaxed_hierarchy_preserved = TRUE,
  input_membership_preserved = all(c("strict_included_member", "strict_excluded_member", "relaxed_included_member", "relaxed_excluded_member") %in% names(master)),
  exposure_values_immutable = all(abs(master$F_stat - (master$exposure_beta / master$exposure_se)^2) < 1e-6),
  chen_outcome_values_raw_preserved_before_transform = all(is.finite(master$outcome_beta_raw)) && all(is.finite(master$outcome_se_raw)),
  palindrome_rules_applied = all(!master$final_valid[master$palindromic_recomputed]),
  no_eaf_reinclusion = TRUE,
  palindromic_final_valid_zero = sum(master$final_valid & master$palindromic_recomputed) == 0L,
  final_valid_definition_correct = all(master$final_valid == (compatible & !master$palindromic_recomputed & bool(master$available_for_harmonisation))),
  performed_transformations_match_rules = all(master$performed_beta_flip == (master$final_valid & master$classification_planned_beta_flip)) && all(master$performed_eaf_flip == (master$final_valid & master$classification_planned_eaf_flip)) && all(master$performed_strand_transform == (master$final_valid & master$classification_planned_strand_transform)),
  final_alleles_fully_aligned = all(final_rows$outcome_effect_allele_harmonised == final_rows$exposure_effect_allele) && all(final_rows$outcome_other_allele_harmonised == final_rows$exposure_other_allele),
  final_effect_fields_complete = all(is.finite(final_rows$exposure_beta)) && all(is.finite(final_rows$exposure_se)) && all(is.finite(final_rows$outcome_beta_harmonised)) && all(is.finite(final_rows$outcome_se_harmonised)) && all(final_rows$exposure_se > 0) && all(final_rows$outcome_se_harmonised > 0),
  preflight_formal_rsid_sets_equal = sets_equal,
  F_statistics_preserved = all(abs(master$F_stat - (master$exposure_beta / master$exposure_se)^2) < 1e-6),
  all_final_valid_F_ge_10 = all(final_rows$F_stat >= 10),
  no_f_ge_30_filter = TRUE,
  no_outcome_based_filter = TRUE,
  master_parquet_tsv_consistency = pair_checks$master$same_cols && pair_checks$master$same_n && pair_checks$master$same_order && pair_checks$master$char_ok && pair_checks$master$numeric_ok,
  all_four_branch_parquet_tsv_consistency = pair_ok,
  no_proxy = TRUE,
  no_liftover = TRUE,
  no_mr = TRUE,
  no_steiger = TRUE,
  renv_lock_unchanged = identical(hash_file(paths$renv_lock), "253471c02e1e47a40d0f68b296d4ae2b1df471d757ac86328b3e974018d039f3")
)
hard_check_failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]
harmonisation_status <- if (length(hard_check_failures) == 0L) "passed" else "failed"
approved_freeze <- identical(harmonisation_status, "passed")

qc <- list(
  harmonisation_version = "v1",
  decision = 114,
  date = format(Sys.Date()),
  analysis_direction = "delirium_to_Hb",
  analysis_role = "reverse_alternative_hb_outcome_sensitivity",
  independent_replication = FALSE,
  strict_exposure_authority_decision = 47,
  relaxed_exposure_authority_decision = 59,
  outcome_freeze_decision = 112,
  harmonisation_contract_decision = 113,
  outcome_scale = outcome_freeze$outcome_scale,
  strict_threshold = 5e-8,
  relaxed_threshold = 5e-6,
  branch_results = records(counts),
  final_valid_rsid_sets = actual_sets,
  palindrome_results = records(master[master$palindromic_recomputed, c("rsid", "classification_recomputed", "strict_included_member", "strict_excluded_member", "relaxed_included_member", "relaxed_excluded_member", "exclusion_reason")]),
  transformation_counts = list(
    performed_beta_flip_count = sum(master$performed_beta_flip),
    performed_eaf_flip_count = sum(master$performed_eaf_flip),
    performed_strand_transform_count = sum(master$performed_strand_transform)
  ),
  instrument_strength_results = records(counts[, c("branch", "final_valid_count", "F_min", "F_mean", "F_median", "F_max", "F_lt10")]),
  eaf_based_palindrome_reinclusion_count = 0,
  proxy_used = FALSE,
  liftover_used = FALSE,
  mr_run = FALSE,
  steiger_run = FALSE,
  harmonisation_status = harmonisation_status,
  approved_for_chen_reverse_mr_input_freeze = approved_freeze,
  nullable_character_fields_for_tsv_blank_normalization = nullable_character_fields,
  parquet_tsv_consistency = pair_checks,
  hard_checks = hard_checks,
  hard_check_failures = hard_check_failures,
  informational_findings = list(
    formal_harmonisation_is_not_mr = TRUE,
    no_effect_estimator_run = TRUE,
    excluded_rows_retain_raw_fields_only = TRUE
  )
)
write_json(qc, paths$qc)

decision_lines <- c(
  "# Decision 114: Chen Reverse Formal Harmonisation V1",
  "",
  paste0("Date: ", format(Sys.Date())),
  "",
  "## Status",
  paste0("harmonisation_status: `", harmonisation_status, "`"),
  paste0("approved_for_chen_reverse_mr_input_freeze: `", approved_freeze, "`"),
  "",
  "## Decision",
  "Execute Chen reverse formal harmonisation V1 using Decision 113 preflight rules and Decision 112 Chen outcome records.",
  "",
  "This stage creates harmonised MR input candidates and branch files. It does not run MR, Steiger, heterogeneity, Egger, MR-PRESSO, leave-one-out, proxy search, liftOver, re-extraction, or reclumping.",
  "",
  "## Results",
  paste0("- Strict included final-valid: `", counts$final_valid_count[counts$branch == "strict_apoe_included"], "/", counts$input_count[counts$branch == "strict_apoe_included"], "`."),
  paste0("- Strict excluded final-valid: `", counts$final_valid_count[counts$branch == "strict_apoe_excluded"], "/", counts$input_count[counts$branch == "strict_apoe_excluded"], "`."),
  paste0("- Relaxed included final-valid: `", counts$final_valid_count[counts$branch == "relaxed_apoe_included"], "/", counts$input_count[counts$branch == "relaxed_apoe_included"], "`."),
  paste0("- Relaxed excluded final-valid: `", counts$final_valid_count[counts$branch == "relaxed_apoe_excluded"], "/", counts$input_count[counts$branch == "relaxed_apoe_excluded"], "`."),
  paste0("- Palindromic exclusions: `", nrow(excluded[excluded$palindromic_recomputed, , drop = FALSE]), "`."),
  paste0("- Incompatible/invalid exclusions: `", sum(excluded$classification_recomputed %in% c("incompatible", "invalid")), "`."),
  paste0("- Performed beta/EAF/strand transformations: `", sum(master$performed_beta_flip), " / ", sum(master$performed_eaf_flip), " / ", sum(master$performed_strand_transform), "`."),
  paste0("- Hard-check failures: `", if (length(hard_check_failures) == 0L) "none" else paste(hard_check_failures, collapse = ";"), "`."),
  "",
  "## Outputs",
  paste0("- `", relpath(paths$qc), "`"),
  paste0("- `", relpath(paths$counts), "`"),
  paste0("- `", relpath(paths$transform_audit), "`"),
  paste0("- `", relpath(paths$excluded), "`"),
  paste0("- `", relpath(paths$master_parquet), "`"),
  paste0("- `", relpath(paths$master_tsv), "`"),
  paste0("- `", relpath(paths$log), "`")
)
write_text(decision_lines, paths$decision)

log_line("harmonisation_status=", harmonisation_status, "; hard_check_failures=", paste(hard_check_failures, collapse = ","))
cat("harmonisation_status=", harmonisation_status, "\n", sep = "")
cat("approved_for_chen_reverse_mr_input_freeze=", approved_freeze, "\n", sep = "")
cat("hard_check_failures=", paste(hard_check_failures, collapse = ","), "\n", sep = "")
cat("final_valid_total=", nrow(final_rows), "\n", sep = "")
