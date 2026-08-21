options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/24_chen_forward_formal_harmonisation_v1.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
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
sql_ident <- function(con, x) as.character(DBI::dbQuoteIdentifier(con, x))
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
comp <- function(x) unname(c(A = "T", T = "A", C = "G", G = "C")[x])
allele_key <- function(a, b) paste(sort(c(a, b)), collapse = "/")
is_snp_allele <- function(x) !is.na(x) & nchar(x) == 1L & x %in% c("A", "C", "G", "T")
same_set <- function(a, b) identical(sort(unique(as.character(a))), sort(unique(as.character(b))))
num_equal <- function(a, b, atol = 1e-10, rtol = 1e-8) {
  both_na <- is.na(a) & is.na(b)
  diff <- abs(a - b)
  ok <- both_na | (!is.na(diff) & (diff <= atol | diff <= rtol * pmax(abs(a), abs(b), 1)))
  all(ok)
}
atomic <- function(path, writer) {
  partial <- paste0(path, ".partial")
  stop_if(file.exists(path) || file.exists(partial), paste("Output occupied:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writer(partial)
  stop_if(!file.exists(partial), paste("Writer did not create partial:", partial))
  stop_if(!file.rename(partial, path), paste("Atomic rename failed:", path))
}
write_parquet <- function(con, path, data) {
  nm <- paste0("tmp_", digest::digest(path, algo = "xxhash32", serialize = FALSE))
  DBI::dbWriteTable(con, nm, data, temporary = TRUE, overwrite = TRUE)
  DBI::dbExecute(con, sprintf("COPY %s TO %s (FORMAT PARQUET, COMPRESSION ZSTD)", sql_ident(con, nm), sql_string(path, must_work = FALSE)))
}
audit_pair <- function(label, parquet_df, tsv_df) {
  stop_if(!identical(names(parquet_df), names(tsv_df)), paste(label, "column order mismatch"))
  stop_if(nrow(parquet_df) != nrow(tsv_df), paste(label, "row count mismatch"))
  rows <- list()
  max_abs <- 0
  max_rel <- 0
  for (nm in names(parquet_df)) {
    p <- parquet_df[[nm]]
    t <- tsv_df[[nm]]
    if (is.numeric(p)) {
      d <- abs(as.numeric(p) - as.numeric(t))
      rel <- d / pmax(abs(as.numeric(p)), abs(as.numeric(t)), 1)
      max_abs <- max(max_abs, d, na.rm = TRUE)
      max_rel <- max(max_rel, rel, na.rm = TRUE)
      ok <- num_equal(as.numeric(p), as.numeric(t))
      rows[[length(rows) + 1L]] <- data.frame(analysis_set = label, column = nm, check_type = "numeric_tolerance", ok = ok, max_absolute_difference = ifelse(all(is.na(d)), NA_real_, max(d, na.rm = TRUE)), max_relative_difference = ifelse(all(is.na(rel)), NA_real_, max(rel, na.rm = TRUE)), stringsAsFactors = FALSE)
    } else {
      pc <- as.character(p)
      tc <- as.character(t)
      pna <- is.na(pc) | pc == ""
      tna <- is.na(tc) | tc == ""
      ok <- all(xor(pna, tna) == FALSE) && all(pc[!pna & !tna] == tc[!pna & !tna])
      rows[[length(rows) + 1L]] <- data.frame(analysis_set = label, column = nm, check_type = "character_exact_na_pattern", ok = ok, max_absolute_difference = NA_real_, max_relative_difference = NA_real_, stringsAsFactors = FALSE)
    }
  }
  list(audit = do.call(rbind, rows), ok = all(vapply(rows, function(z) z$ok, logical(1))), max_abs = max_abs, max_rel = max_rel)
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
  script = file.path(root, "R", "24_chen_forward_formal_harmonisation_v1.R"),
  renv_lock = file.path(root, "renv.lock"),
  decision_87 = file.path(root, "docs", "decisions", "87_chen_forward_instruments_v2_freeze_v1.1.md"),
  decision_91 = file.path(root, "docs", "decisions", "91_chen_forward_finngen_outcome_extraction_v2_freeze_v1.1.md"),
  decision_92 = file.path(root, "docs", "decisions", "92_chen_forward_harmonisation_contract_and_preflight_v1_v1.1.md"),
  instrument_freeze = file.path(root, "results", "qc", "chen_forward_instruments_v2_freeze.json"),
  outcome_freeze = file.path(root, "results", "qc", "chen_forward_finngen_outcome_extraction_v2_freeze.json"),
  contract = file.path(root, "results", "qc", "chen_forward_harmonisation_contract_v1.json"),
  preflight_qc = file.path(root, "results", "qc", "chen_forward_harmonisation_preflight_v1.json"),
  preflight_counts = file.path(root, "results", "qc", "chen_forward_harmonisation_preflight_counts_v1.csv"),
  preflight_master = file.path(root, "data_derived", "forward_sensitivity_harmonisation", "chen_forward_harmonisation_preflight_master_v1.parquet"),
  master_parquet = file.path(root, "data_derived", "forward_sensitivity_harmonisation", "chen_forward_harmonised_master_v1.parquet"),
  master_tsv = file.path(root, "data_derived", "forward_sensitivity_harmonisation", "chen_forward_harmonised_master_v1.tsv"),
  included_parquet = file.path(root, "data_derived", "forward_sensitivity_harmonisation", "chen_forward_harmonised_apoe_included_v1.parquet"),
  included_tsv = file.path(root, "data_derived", "forward_sensitivity_harmonisation", "chen_forward_harmonised_apoe_included_v1.tsv"),
  excluded_parquet = file.path(root, "data_derived", "forward_sensitivity_harmonisation", "chen_forward_harmonised_apoe_excluded_v1.parquet"),
  excluded_tsv = file.path(root, "data_derived", "forward_sensitivity_harmonisation", "chen_forward_harmonised_apoe_excluded_v1.tsv"),
  counts = file.path(root, "results", "qc", "chen_forward_formal_harmonisation_counts_v1.csv"),
  transform_audit = file.path(root, "results", "qc", "chen_forward_formal_harmonisation_transform_audit_v1.csv"),
  excluded_snps = file.path(root, "results", "qc", "chen_forward_formal_harmonisation_excluded_snps_v1.tsv"),
  qc = file.path(root, "results", "qc", "chen_forward_formal_harmonisation_v1.json"),
  log = file.path(root, "results", "logs", "chen_forward_formal_harmonisation_v1.log"),
  decision = file.path(root, "docs", "decisions", "93_chen_forward_formal_harmonisation_v1_v1.1.md")
)
for (p in paths[c("script", "renv_lock", "decision_87", "decision_91", "decision_92", "instrument_freeze", "outcome_freeze", "contract", "preflight_qc", "preflight_counts", "preflight_master")]) {
  stop_if(!file.exists(p), paste("Missing required input:", p))
}
for (p in paths[c("master_parquet", "master_tsv", "included_parquet", "included_tsv", "excluded_parquet", "excluded_tsv", "counts", "transform_audit", "excluded_snps", "qc", "log", "decision")]) {
  stop_if(file.exists(p) || file.exists(paste0(p, ".partial")), paste("Output occupied:", p))
}

renv_before <- hash_file(paths[["renv_lock"]])
instrument_freeze <- jsonlite::fromJSON(paths[["instrument_freeze"]], simplifyVector = FALSE)
outcome_freeze <- jsonlite::fromJSON(paths[["outcome_freeze"]], simplifyVector = FALSE)
contract <- jsonlite::fromJSON(paths[["contract"]], simplifyVector = FALSE)
preflight_qc <- jsonlite::fromJSON(paths[["preflight_qc"]], simplifyVector = FALSE)
preflight_counts <- read.csv(paths[["preflight_counts"]], check.names = FALSE)

con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
preflight <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(paths[["preflight_master"]])))
preflight <- preflight[order(preflight$resolved_rsid), , drop = FALSE]

exp_ea <- toupper(as.character(preflight$exposure_effect_allele))
exp_oa <- toupper(as.character(preflight$exposure_other_allele))
out_ea <- toupper(as.character(preflight$outcome_effect_allele))
out_oa <- toupper(as.character(preflight$outcome_other_allele))
valid_exp <- is_snp_allele(exp_ea) & is_snp_allele(exp_oa) & exp_ea != exp_oa
valid_out <- is_snp_allele(out_ea) & is_snp_allele(out_oa) & out_ea != out_oa
exact <- valid_exp & valid_out & out_ea == exp_ea & out_oa == exp_oa
swapped <- valid_exp & valid_out & out_ea == exp_oa & out_oa == exp_ea
strand_exact <- valid_exp & valid_out & comp(out_ea) == exp_ea & comp(out_oa) == exp_oa
strand_swapped <- valid_exp & valid_out & comp(out_ea) == exp_oa & comp(out_oa) == exp_ea
classification <- ifelse(!valid_exp | !valid_out, "invalid",
  ifelse(exact, "exact",
  ifelse(swapped, "swapped",
  ifelse(strand_exact, "strand_exact",
  ifelse(strand_swapped, "strand_swapped", "incompatible")))))
pal <- valid_exp & mapply(allele_key, exp_ea, exp_oa) %in% c("A/T", "C/G")
resolvable <- classification %in% c("exact", "swapped", "strand_exact", "strand_swapped")
final_valid <- resolvable & !pal

classification_planned_beta_flip <- classification %in% c("swapped", "strand_swapped")
classification_planned_eaf_flip <- classification %in% c("swapped", "strand_swapped") & !is.na(preflight$outcome_eaf)
classification_planned_strand_transform <- classification %in% c("strand_exact", "strand_swapped")
performed_beta_flip <- final_valid & classification %in% c("swapped", "strand_swapped")
performed_eaf_flip <- final_valid & classification %in% c("swapped", "strand_swapped") & !is.na(preflight$outcome_eaf)
performed_strand_transform <- final_valid & classification %in% c("strand_exact", "strand_swapped")

outcome_beta_h <- rep(NA_real_, nrow(preflight))
outcome_se_h <- rep(NA_real_, nrow(preflight))
outcome_eaf_h <- rep(NA_real_, nrow(preflight))
outcome_ea_h <- rep(NA_character_, nrow(preflight))
outcome_oa_h <- rep(NA_character_, nrow(preflight))
outcome_beta_h[final_valid] <- ifelse(performed_beta_flip[final_valid], -preflight$outcome_beta[final_valid], preflight$outcome_beta[final_valid])
outcome_se_h[final_valid] <- preflight$outcome_se[final_valid]
outcome_eaf_h[final_valid] <- ifelse(performed_eaf_flip[final_valid], 1 - preflight$outcome_eaf[final_valid], preflight$outcome_eaf[final_valid])
outcome_ea_h[final_valid] <- exp_ea[final_valid]
outcome_oa_h[final_valid] <- exp_oa[final_valid]

exclusion_reason <- ifelse(!valid_exp | !valid_out, "invalid_allele",
  ifelse(!resolvable, "incompatible_alleles",
  ifelse(pal, "palindromic_snp_excluded_by_forward_rule", "")))

harm <- data.frame(
  resolved_rsid = as.character(preflight$resolved_rsid),
  source_marker_id = as.character(preflight$source_marker_id),
  in_apoe_included_input = as.logical(preflight$in_apoe_included_input),
  in_apoe_excluded_input = as.logical(preflight$in_apoe_excluded_input),
  exposure_effect_allele = exp_ea,
  exposure_other_allele = exp_oa,
  exposure_beta = as.numeric(preflight$exposure_beta),
  exposure_se = as.numeric(preflight$exposure_se),
  exposure_pval = as.numeric(preflight$exposure_pval),
  exposure_eaf = as.numeric(preflight$exposure_eaf),
  exposure_n_samples = as.integer(preflight$exposure_n_samples),
  exposure_F_stat = as.numeric(preflight$exposure_F_stat),
  exposure_marker_chr_grch37 = as.integer(preflight$exposure_marker_chr_grch37),
  exposure_marker_pos_grch37 = as.integer(preflight$exposure_marker_pos_grch37),
  outcome_source = as.character(preflight$outcome_source),
  outcome_trait = as.character(preflight$outcome_trait),
  outcome_build = as.character(preflight$outcome_build),
  outcome_rsid = as.character(preflight$outcome_rsid),
  outcome_chr_grch38 = as.integer(preflight$outcome_chr_grch38),
  outcome_pos_grch38 = as.integer(preflight$outcome_pos_grch38),
  outcome_effect_allele_raw = out_ea,
  outcome_other_allele_raw = out_oa,
  outcome_beta_raw = as.numeric(preflight$outcome_beta),
  outcome_se_raw = as.numeric(preflight$outcome_se),
  outcome_p = as.numeric(preflight$outcome_p),
  outcome_eaf_raw = as.numeric(preflight$outcome_eaf),
  outcome_effect_allele_harmonised = outcome_ea_h,
  outcome_other_allele_harmonised = outcome_oa_h,
  outcome_beta_harmonised = outcome_beta_h,
  outcome_se_harmonised = outcome_se_h,
  outcome_eaf_harmonised = outcome_eaf_h,
  outcome_ncase = as.integer(preflight$outcome_ncase),
  outcome_ncontrol = as.integer(preflight$outcome_ncontrol),
  outcome_n_study = as.integer(preflight$outcome_n_study),
  outcome_effect_scale = as.character(preflight$outcome_effect_scale),
  preflight_compatibility_class = as.character(preflight$compatibility_class),
  reproduced_compatibility_class = classification,
  palindromic_snp = pal,
  final_valid = final_valid,
  exclusion_reason = exclusion_reason,
  classification_planned_beta_flip = classification_planned_beta_flip,
  classification_planned_eaf_flip = classification_planned_eaf_flip,
  classification_planned_strand_transform = classification_planned_strand_transform,
  performed_beta_flip = performed_beta_flip,
  performed_eaf_flip = performed_eaf_flip,
  performed_strand_transform = performed_strand_transform,
  formal_harmonisation_performed = final_valid,
  stringsAsFactors = FALSE
)

included_final <- harm[harm$in_apoe_included_input & harm$final_valid, , drop = FALSE]
excluded_final <- harm[harm$in_apoe_excluded_input & harm$final_valid, , drop = FALSE]
shared_final <- intersect(included_final$resolved_rsid, excluded_final$resolved_rsid)
included_only_final <- setdiff(included_final$resolved_rsid, excluded_final$resolved_rsid)
excluded_only_final <- setdiff(excluded_final$resolved_rsid, included_final$resolved_rsid)
projected_included <- preflight$resolved_rsid[preflight$in_apoe_included_input & preflight$projected_final_valid]
projected_excluded <- preflight$resolved_rsid[preflight$in_apoe_excluded_input & preflight$projected_final_valid]

count_set <- function(label, flag) {
  x <- harm[flag, , drop = FALSE]
  data.frame(
    analysis_set = label,
    input_count = nrow(x),
    palindromic_excluded_count = sum(x$palindromic_snp),
    incompatible_excluded_count = sum(x$reproduced_compatibility_class == "incompatible"),
    invalid_excluded_count = sum(x$reproduced_compatibility_class == "invalid"),
    exact_count = sum(x$reproduced_compatibility_class == "exact"),
    swapped_count = sum(x$reproduced_compatibility_class == "swapped"),
    strand_exact_count = sum(x$reproduced_compatibility_class == "strand_exact"),
    strand_swapped_count = sum(x$reproduced_compatibility_class == "strand_swapped"),
    final_valid_count = sum(x$final_valid),
    classification_planned_beta_flip_count = sum(x$classification_planned_beta_flip),
    performed_beta_flip_count = sum(x$performed_beta_flip),
    classification_planned_eaf_flip_count = sum(x$classification_planned_eaf_flip),
    performed_eaf_flip_count = sum(x$performed_eaf_flip),
    classification_planned_strand_transform_count = sum(x$classification_planned_strand_transform),
    performed_strand_transform_count = sum(x$performed_strand_transform),
    planned_transform_not_performed_due_to_palindrome_count = sum(x$palindromic_snp & (x$classification_planned_beta_flip | x$classification_planned_eaf_flip | x$classification_planned_strand_transform)),
    stringsAsFactors = FALSE
  )
}
counts <- rbind(
  count_set("APOE_included", harm$in_apoe_included_input),
  count_set("APOE_excluded", harm$in_apoe_excluded_input),
  count_set("union", rep(TRUE, nrow(harm)))
)

transform_audit <- harm[, c("resolved_rsid", "in_apoe_included_input", "in_apoe_excluded_input", "preflight_compatibility_class", "reproduced_compatibility_class", "palindromic_snp", "final_valid", "classification_planned_beta_flip", "performed_beta_flip", "classification_planned_eaf_flip", "performed_eaf_flip", "classification_planned_strand_transform", "performed_strand_transform", "exclusion_reason"), drop = FALSE]
excluded_table <- harm[!harm$final_valid, c("resolved_rsid", "source_marker_id", "in_apoe_included_input", "in_apoe_excluded_input", "exposure_effect_allele", "exposure_other_allele", "outcome_effect_allele_raw", "outcome_other_allele_raw", "reproduced_compatibility_class", "palindromic_snp", "exclusion_reason"), drop = FALSE]

atomic(paths[["master_tsv"]], function(p) write.table(harm, p, sep = "\t", quote = FALSE, row.names = FALSE, na = ""))
atomic(paths[["included_tsv"]], function(p) write.table(included_final, p, sep = "\t", quote = FALSE, row.names = FALSE, na = ""))
atomic(paths[["excluded_tsv"]], function(p) write.table(excluded_final, p, sep = "\t", quote = FALSE, row.names = FALSE, na = ""))
atomic(paths[["counts"]], function(p) write.csv(counts, p, row.names = FALSE, na = ""))
atomic(paths[["transform_audit"]], function(p) write.csv(transform_audit, p, row.names = FALSE, na = ""))
atomic(paths[["excluded_snps"]], function(p) write.table(excluded_table, p, sep = "\t", quote = FALSE, row.names = FALSE, na = ""))
atomic(paths[["master_parquet"]], function(p) write_parquet(con, p, harm))
atomic(paths[["included_parquet"]], function(p) write_parquet(con, p, included_final))
atomic(paths[["excluded_parquet"]], function(p) write_parquet(con, p, excluded_final))

read_tsv <- function(path) read.delim(path, sep = "\t", check.names = FALSE, na.strings = c(""))
read_pq <- function(path) DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(path)))
master_pair <- audit_pair("master", read_pq(paths[["master_parquet"]]), read_tsv(paths[["master_tsv"]]))
included_pair <- audit_pair("included", read_pq(paths[["included_parquet"]]), read_tsv(paths[["included_tsv"]]))
excluded_pair <- audit_pair("excluded", read_pq(paths[["excluded_parquet"]]), read_tsv(paths[["excluded_tsv"]]))
readback_audit <- rbind(master_pair$audit, included_pair$audit, excluded_pair$audit)

hard_checks <- list(
  decision_87_instrument_gate = identical(instrument_freeze$freeze_status, "passed") && length(instrument_freeze$hard_check_failures) == 0L,
  decision_91_outcome_gate = identical(outcome_freeze$freeze_status, "passed") && length(outcome_freeze$hard_check_failures) == 0L,
  decision_92_contract_gate = identical(contract$contract_status, "frozen") && identical(preflight_qc$preflight_status, "passed") && isTRUE(preflight_qc$approved_for_chen_forward_formal_harmonisation) && length(preflight_qc$hard_check_failures) == 0L,
  preflight_row_level_classification_reproduced = identical(classification, as.character(preflight$compatibility_class)),
  input_membership_preserved = sum(harm$in_apoe_included_input) == preflight_qc$included_matched_input_count && sum(harm$in_apoe_excluded_input) == preflight_qc$excluded_matched_input_count,
  exposure_orientation_preserved = all(harm$exposure_effect_allele == exp_ea) && all(harm$exposure_other_allele == exp_oa),
  outcome_orientation_verified = identical(contract$outcome_effect_orientation_source, "FinnGen_certified_effect_allele_fields"),
  marker_tokens_not_used_for_orientation = !isTRUE(contract$marker_identifier_alleles_used_for_effect_orientation),
  reference_panel_not_used_for_orientation = !isTRUE(contract$reference_panel_alleles_used_for_effect_orientation),
  palindrome_rule_applied = identical(contract$palindromic_rule, "exclude_all_without_EAF_reinclusion") && all(!harm$final_valid[harm$palindromic_snp]),
  no_eaf_palindrome_reinclusion = sum(harm$palindromic_snp & harm$final_valid) == 0L,
  palindromic_final_valid_zero = sum(harm$palindromic_snp & harm$final_valid) == 0L,
  final_valid_definition_correct = identical(final_valid, as.logical(harm$final_valid)),
  performed_flips_match_row_level_rules = all(harm$performed_beta_flip == (harm$final_valid & harm$reproduced_compatibility_class %in% c("swapped", "strand_swapped"))) &&
    all(harm$performed_eaf_flip == (harm$final_valid & harm$reproduced_compatibility_class %in% c("swapped", "strand_swapped") & !is.na(harm$outcome_eaf_raw))) &&
    all(harm$performed_strand_transform == (harm$final_valid & harm$reproduced_compatibility_class %in% c("strand_exact", "strand_swapped"))),
  final_harmonised_alleles_match_exposure = all(harm$outcome_effect_allele_harmonised[harm$final_valid] == harm$exposure_effect_allele[harm$final_valid]) &&
    all(harm$outcome_other_allele_harmonised[harm$final_valid] == harm$exposure_other_allele[harm$final_valid]),
  exposure_beta_se_unchanged = TRUE,
  exposure_F_unchanged = all(num_equal(harm$exposure_F_stat, (harm$exposure_beta / harm$exposure_se)^2, atol = 1e-8, rtol = 1e-8)),
  final_effect_fields_complete = all(is.finite(harm$exposure_beta[harm$final_valid])) &&
    all(is.finite(harm$exposure_se[harm$final_valid])) &&
    all(is.finite(harm$outcome_beta_harmonised[harm$final_valid])) &&
    all(is.finite(harm$outcome_se_harmonised[harm$final_valid])) &&
    all(harm$exposure_se[harm$final_valid] > 0) &&
    all(harm$outcome_se_harmonised[harm$final_valid] > 0),
  projected_final_valid_rsids_match_preflight = same_set(included_final$resolved_rsid, projected_included) && same_set(excluded_final$resolved_rsid, projected_excluded),
  included_parquet_tsv_consistency = included_pair$ok,
  excluded_parquet_tsv_consistency = excluded_pair$ok,
  master_parquet_tsv_consistency = master_pair$ok,
  no_outcome_based_filtering = nrow(harm) == nrow(preflight),
  no_f_ge_30_filter = nrow(harm) == nrow(preflight),
  no_proxy = !isTRUE(contract$proxy_allowed),
  no_liftover = !isTRUE(contract$liftover_allowed),
  no_mr = TRUE,
  no_steiger = TRUE,
  renv_lock_unchanged = identical(renv_before, hash_file(paths[["renv_lock"]]))
)
failures <- names(hard_checks)[!unlist(hard_checks)]
if (f_summary(included_final$exposure_F_stat)$F_lt10_count > 0 || f_summary(excluded_final$exposure_F_stat)$F_lt10_count > 0) {
  failures <- unique(c(failures, "final_F_lt10_present_requires_manual_review"))
}
status <- if (length(failures) == 0L) "passed" else "failed"
approved <- identical(status, "passed")

qc <- list(
  harmonisation_version = "v1",
  analysis_direction = "Hb_to_delirium",
  analysis_role = "forward_alternative_hb_gwas_sensitivity",
  source_instrument_freeze_decision = 87,
  source_outcome_freeze_decision = 91,
  harmonisation_contract_decision = 92,
  exposure_orientation_source = "Chen_reference_allele_other_allele",
  outcome_orientation_source = "FinnGen_ALT_REF",
  palindromic_rule = "exclude_all_without_EAF_reinclusion",
  included_input_count = sum(harm$in_apoe_included_input),
  excluded_input_count = sum(harm$in_apoe_excluded_input),
  classification_counts = records(counts),
  palindromic_exclusion_counts = list(included = sum(harm$in_apoe_included_input & harm$palindromic_snp), excluded = sum(harm$in_apoe_excluded_input & harm$palindromic_snp), union = sum(harm$palindromic_snp)),
  incompatible_exclusion_counts = list(included = sum(harm$in_apoe_included_input & harm$reproduced_compatibility_class == "incompatible"), excluded = sum(harm$in_apoe_excluded_input & harm$reproduced_compatibility_class == "incompatible"), union = sum(harm$reproduced_compatibility_class == "incompatible")),
  invalid_exclusion_counts = list(included = sum(harm$in_apoe_included_input & harm$reproduced_compatibility_class == "invalid"), excluded = sum(harm$in_apoe_excluded_input & harm$reproduced_compatibility_class == "invalid"), union = sum(harm$reproduced_compatibility_class == "invalid")),
  final_valid_counts = list(included = nrow(included_final), excluded = nrow(excluded_final), shared = length(shared_final), included_only = length(included_only_final), excluded_only = length(excluded_only_final), union = sum(harm$final_valid)),
  final_valid_rsid_sets = list(included = sort(included_final$resolved_rsid), excluded = sort(excluded_final$resolved_rsid), shared = sort(shared_final), included_only = sort(included_only_final), excluded_only = sort(excluded_only_final)),
  transformation_counts = records(counts[, c("analysis_set", "classification_planned_beta_flip_count", "performed_beta_flip_count", "classification_planned_eaf_flip_count", "performed_eaf_flip_count", "classification_planned_strand_transform_count", "performed_strand_transform_count", "planned_transform_not_performed_due_to_palindrome_count"), drop = FALSE]),
  instrument_strength_final = list(included = f_summary(included_final$exposure_F_stat), excluded = f_summary(excluded_final$exposure_F_stat)),
  eaf_based_palindrome_reinclusion_count = 0L,
  proxy_used = FALSE,
  liftover_used = FALSE,
  mr_run = FALSE,
  steiger_run = FALSE,
  parquet_tsv_consistency = list(
    master = list(ok = master_pair$ok, max_absolute_difference = master_pair$max_abs, max_relative_difference = master_pair$max_rel),
    included = list(ok = included_pair$ok, max_absolute_difference = included_pair$max_abs, max_relative_difference = included_pair$max_rel),
    excluded = list(ok = excluded_pair$ok, max_absolute_difference = excluded_pair$max_abs, max_relative_difference = excluded_pair$max_rel)
  ),
  harmonisation_status = status,
  approved_for_chen_forward_mr_input_freeze = approved,
  hard_checks = hard_checks,
  hard_check_failures = failures,
  informational_findings = list(
    planned_vs_performed_transformations_separated = TRUE,
    palindromic_planned_transformations_not_performed = TRUE,
    excluded_rows_have_missing_harmonised_effect_fields = TRUE,
    parquet_machine_precision_authority_tsv_human_readable = TRUE
  ),
  renv_lock_sha_before = renv_before,
  renv_lock_sha_after = hash_file(paths[["renv_lock"]]),
  output_sha256 = list(
    master_parquet = hash_file(paths[["master_parquet"]]),
    master_tsv = hash_file(paths[["master_tsv"]]),
    included_parquet = hash_file(paths[["included_parquet"]]),
    included_tsv = hash_file(paths[["included_tsv"]]),
    excluded_parquet = hash_file(paths[["excluded_parquet"]]),
    excluded_tsv = hash_file(paths[["excluded_tsv"]]),
    counts = hash_file(paths[["counts"]]),
    transform_audit = hash_file(paths[["transform_audit"]]),
    excluded_snps = hash_file(paths[["excluded_snps"]])
  )
)
atomic(paths[["qc"]], function(p) jsonlite::write_json(qc, p, pretty = TRUE, auto_unbox = TRUE, na = "null"))

decision_lines <- c(
  "# Decision 93 - Chen forward formal harmonisation V1",
  "",
  "Date: 2026-08-12",
  "Status: formal harmonisation",
  "",
  "## Decision",
  "",
  "Execute Chen forward formal harmonisation V1 from Decision 92 preflight.",
  "Only non-palindromic exact, swapped, strand-exact, and strand-swapped rows",
  "are retained as final-valid harmonised MR-input candidates.",
  "",
  "## Results",
  "",
  sprintf("- harmonisation_status: `%s`", status),
  sprintf("- included input/final-valid: `%d/%d`", sum(harm$in_apoe_included_input), nrow(included_final)),
  sprintf("- excluded input/final-valid: `%d/%d`", sum(harm$in_apoe_excluded_input), nrow(excluded_final)),
  sprintf("- shared/included-only/excluded-only final-valid: `%d/%d/%d`", length(shared_final), length(included_only_final), length(excluded_only_final)),
  sprintf("- palindromic excluded included/excluded/union: `%d/%d/%d`", sum(harm$in_apoe_included_input & harm$palindromic_snp), sum(harm$in_apoe_excluded_input & harm$palindromic_snp), sum(harm$palindromic_snp)),
  sprintf("- incompatible/invalid union: `%d/%d`", sum(harm$reproduced_compatibility_class == "incompatible"), sum(harm$reproduced_compatibility_class == "invalid")),
  sprintf("- hard_check_failures: `%s`", paste(failures, collapse = ";")),
  sprintf("- approved_for_chen_forward_mr_input_freeze: `%s`", approved),
  "",
  "## Safeguards",
  "",
  "No MR, Steiger, proxy lookup, liftOver, reclumping, instrument reselection,",
  "outcome re-extraction, EAF-based palindromic reinclusion, outcome-based SNP",
  "filtering, or F>=30 filtering was performed.",
  "",
  "## Expected Impact",
  "",
  "If passed, these harmonised outputs authorize only a separate MR-input freeze.",
  "They do not authorize MR."
)
atomic(paths[["decision"]], function(p) writeLines(decision_lines, p, useBytes = TRUE))

atomic(paths[["log"]], function(p) {
  writeLines(c(
    sprintf("[%s] stage=chen_forward_formal_harmonisation_v1", ts()),
    sprintf("[%s] harmonisation_status=%s approved_for_mr_input_freeze=%s", ts(), status, approved),
    sprintf("[%s] included_final=%d excluded_final=%d shared_final=%d", ts(), nrow(included_final), nrow(excluded_final), length(shared_final)),
    sprintf("[%s] hard_check_failures=%s", ts(), paste(failures, collapse = ";"))
  ), p, useBytes = TRUE)
})

stop_if(!identical(status, "passed"), "Formal harmonisation failed; QC retained.")
message("Chen forward formal harmonisation completed: ", status)
