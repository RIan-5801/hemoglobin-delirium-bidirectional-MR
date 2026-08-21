options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/23_chen_forward_harmonisation_preflight_v1.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
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
atomic <- function(path, writer) {
  partial <- paste0(path, ".partial")
  stop_if(file.exists(path) || file.exists(partial), paste("Output occupied:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writer(partial)
  stop_if(!file.exists(partial), paste("Writer did not create partial:", partial))
  stop_if(!file.rename(partial, path), paste("Atomic rename failed:", path))
}
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
comp <- function(x) unname(c(A = "T", T = "A", C = "G", G = "C")[x])
allele_key <- function(a, b) paste(sort(c(a, b)), collapse = "/")
is_snp_allele <- function(x) !is.na(x) & nchar(x) == 1L & x %in% c("A", "C", "G", "T")
same_num <- function(a, b, atol = 1e-10, rtol = 1e-8) {
  both_na <- is.na(a) & is.na(b)
  diff <- abs(a - b)
  ok <- both_na | (!is.na(diff) & (diff <= atol | diff <= rtol * pmax(abs(a), abs(b), 1)))
  all(ok)
}
next_decision <- function() {
  files <- list.files(file.path(root, "docs", "decisions"), pattern = "^[0-9]+_", full.names = FALSE)
  nums <- as.integer(sub("_.*$", "", files))
  n <- 91L
  while (n %in% nums) n <- n + 1L
  n
}
find_forward_palindrome_precedent <- function() {
  files <- list.files(file.path(root, "docs", "decisions"), pattern = "\\.md$", full.names = TRUE)
  hits <- list()
  for (p in files) {
    txt <- paste(readLines(p, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    if (grepl("All palindromic SNPs are excluded", txt, fixed = TRUE) &&
        grepl("not reintroduced based on EAF", txt, fixed = TRUE) &&
        grepl("forward", txt, ignore.case = TRUE)) {
      nm <- basename(p)
      num <- as.integer(sub("_.*$", "", nm))
      hits[[length(hits) + 1L]] <- data.frame(decision = num, path = p, stringsAsFactors = FALSE)
    }
  }
  stop_if(length(hits) == 0L, "Forward-primary palindromic precedent was not found.")
  z <- do.call(rbind, hits)
  z[order(z$decision), , drop = FALSE][1L, , drop = FALSE]
}

paths <- c(
  script = file.path(root, "R", "23_chen_forward_harmonisation_preflight_v1.R"),
  renv_lock = file.path(root, "renv.lock"),
  chen_source = file.path(root, "data_raw", "gwas", "BCX2_HGB_EA_GWAMA.out.gz"),
  finngen_source = file.path(root, "data_raw", "gwas", "finngen_R13_F5_DELIRIUM.gz"),
  metadata = file.path(root, "docs", "02_gwas_metadata_v2.md"),
  chen_certification = file.path(root, "results", "qc", "chen_2020_hb_source_certification_v1.json"),
  chen_dictionary = file.path(root, "results", "qc", "chen_2020_hb_official_source_dictionary_audit_v1.json"),
  instrument_freeze = file.path(root, "results", "qc", "chen_forward_instruments_v2_freeze.json"),
  outcome_freeze = file.path(root, "results", "qc", "chen_forward_finngen_outcome_extraction_v2_freeze.json"),
  outcome_closure = file.path(root, "results", "qc", "chen_forward_finngen_outcome_extraction_v2_readback_closure_v1.json"),
  included_instruments = file.path(root, "data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_included_clumped_v2.tsv"),
  excluded_instruments = file.path(root, "data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_excluded_clumped_v2.tsv"),
  outcome_master_tsv = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_master_v2.tsv"),
  missing = file.path(root, "results", "qc", "chen_forward_finngen_outcome_missing_v2.tsv"),
  contract_json = file.path(root, "results", "qc", "chen_forward_harmonisation_contract_v1.json"),
  preflight_json = file.path(root, "results", "qc", "chen_forward_harmonisation_preflight_v1.json"),
  counts_csv = file.path(root, "results", "qc", "chen_forward_harmonisation_preflight_counts_v1.csv"),
  palindromic_tsv = file.path(root, "results", "qc", "chen_forward_harmonisation_palindromic_snps_v1.tsv"),
  incompatible_tsv = file.path(root, "results", "qc", "chen_forward_harmonisation_incompatible_snps_v1.tsv"),
  preflight_parquet = file.path(root, "data_derived", "forward_sensitivity_harmonisation", "chen_forward_harmonisation_preflight_master_v1.parquet"),
  log = file.path(root, "results", "logs", "chen_forward_harmonisation_preflight_v1.log"),
  decision = file.path(root, "docs", "decisions", "92_chen_forward_harmonisation_contract_and_preflight_v1_v1.1.md")
)

for (p in paths[c("script", "renv_lock", "chen_source", "finngen_source", "metadata", "chen_certification", "chen_dictionary", "instrument_freeze", "outcome_freeze", "outcome_closure", "included_instruments", "excluded_instruments", "outcome_master_tsv", "missing")]) {
  stop_if(!file.exists(p), paste("Missing required input:", p))
}
for (p in paths[c("contract_json", "preflight_json", "counts_csv", "palindromic_tsv", "incompatible_tsv", "preflight_parquet", "log", "decision")]) {
  stop_if(file.exists(p) || file.exists(paste0(p, ".partial")), paste("Output occupied:", p))
}

renv_before <- hash_file(paths[["renv_lock"]])
chen_sha_before <- tolower(hash_file(paths[["chen_source"]]))
finngen_sha_before <- tolower(hash_file(paths[["finngen_source"]]))
expected_chen_sha <- "f1dfea8897cb29f39d891b7922bed2ea95a869bb864d39c56db73e5d69f8abf8"
expected_finngen_sha <- "85637f0f3358807964d4f8a3e500293168a706f1c08c65f3fc5512b65df40ed8"

instrument_freeze <- jsonlite::fromJSON(paths[["instrument_freeze"]], simplifyVector = FALSE)
outcome_freeze <- jsonlite::fromJSON(paths[["outcome_freeze"]], simplifyVector = FALSE)
outcome_closure <- jsonlite::fromJSON(paths[["outcome_closure"]], simplifyVector = FALSE)
chen_cert <- jsonlite::fromJSON(paths[["chen_certification"]], simplifyVector = FALSE)
chen_dict <- jsonlite::fromJSON(paths[["chen_dictionary"]], simplifyVector = FALSE)
precedent <- find_forward_palindrome_precedent()
precedent_rel <- sub(paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", root), "/?"), "", normalizePath(precedent$path, winslash = "/", mustWork = TRUE))

included_inst <- read.delim(paths[["included_instruments"]], sep = "\t", check.names = FALSE)
excluded_inst <- read.delim(paths[["excluded_instruments"]], sep = "\t", check.names = FALSE)
master <- read.delim(paths[["outcome_master_tsv"]], sep = "\t", check.names = FALSE, na.strings = c(""))
missing <- read.delim(paths[["missing"]], sep = "\t", check.names = FALSE)

instrument_union <- rbind(
  transform(included_inst, instrument_membership_source = "APOE_included"),
  transform(excluded_inst, instrument_membership_source = "APOE_excluded")
)
instrument_conflict_cols <- c("exposure_effect_allele", "exposure_other_allele", "beta", "se", "pval", "eaf", "n_samples", "F_stat")
instrument_conflict_counts <- do.call(rbind, lapply(split(instrument_union[, instrument_conflict_cols, drop = FALSE], instrument_union$resolved_rsid), function(d) {
  vapply(d, function(x) length(unique(as.character(x))), integer(1))
}))
instrument_conflict_free <- all(instrument_conflict_counts == 1L)
instrument_map <- instrument_union[!duplicated(instrument_union$resolved_rsid), , drop = FALSE]
rownames(instrument_map) <- instrument_map$resolved_rsid

matched <- master[master$outcome_match_status == "unique_exact_match", , drop = FALSE]
matched <- matched[order(matched$resolved_rsid), , drop = FALSE]
stop_if(anyDuplicated(matched$resolved_rsid), "Matched outcome master contains duplicate resolved_rsid.")
stop_if(any(!matched$resolved_rsid %in% instrument_map$resolved_rsid), "Matched rsID missing from instrument map.")
imap <- instrument_map[matched$resolved_rsid, , drop = FALSE]

exp_ea <- toupper(as.character(matched$exposure_effect_allele))
exp_oa <- toupper(as.character(matched$exposure_other_allele))
out_ea <- toupper(as.character(matched$outcome_effect_allele))
out_oa <- toupper(as.character(matched$outcome_other_allele))
valid_exp <- is_snp_allele(exp_ea) & is_snp_allele(exp_oa) & exp_ea != exp_oa
valid_out <- is_snp_allele(out_ea) & is_snp_allele(out_oa) & out_ea != out_oa
exact <- valid_exp & valid_out & out_ea == exp_ea & out_oa == exp_oa
swapped <- valid_exp & valid_out & out_ea == exp_oa & out_oa == exp_ea
strand_exact <- valid_exp & valid_out & comp(out_ea) == exp_ea & comp(out_oa) == exp_oa
strand_swapped <- valid_exp & valid_out & comp(out_ea) == exp_oa & comp(out_oa) == exp_ea
compatibility_class <- ifelse(!valid_exp | !valid_out, "invalid",
  ifelse(exact, "exact",
  ifelse(swapped, "swapped",
  ifelse(strand_exact, "strand_exact",
  ifelse(strand_swapped, "strand_swapped", "incompatible")))))
palindromic <- valid_exp & mapply(allele_key, exp_ea, exp_oa) %in% c("A/T", "C/G")
resolvable <- compatibility_class %in% c("exact", "swapped", "strand_exact", "strand_swapped")

preflight <- data.frame(
  resolved_rsid = as.character(matched$resolved_rsid),
  source_marker_id = as.character(matched$source_marker_id),
  in_apoe_included_input = as.logical(matched$in_apoe_included_input),
  in_apoe_excluded_input = as.logical(matched$in_apoe_excluded_input),
  exposure_effect_allele = exp_ea,
  exposure_other_allele = exp_oa,
  exposure_beta = as.numeric(matched$exposure_beta),
  exposure_se = as.numeric(matched$exposure_se),
  exposure_pval = as.numeric(matched$exposure_pval),
  exposure_eaf = as.numeric(matched$exposure_eaf),
  exposure_n_samples = as.integer(matched$exposure_n_samples),
  exposure_F_stat = as.numeric(matched$exposure_F_stat),
  exposure_marker_chr_grch37 = as.integer(matched$exposure_marker_chr_grch37),
  exposure_marker_pos_grch37 = as.integer(matched$exposure_marker_pos_grch37),
  outcome_source = as.character(matched$outcome_source),
  outcome_trait = as.character(matched$outcome_trait),
  outcome_build = as.character(matched$outcome_build),
  outcome_rsid = as.character(matched$outcome_rsid),
  outcome_chr_grch38 = as.integer(matched$outcome_chr_grch38),
  outcome_pos_grch38 = as.integer(matched$outcome_pos_grch38),
  outcome_effect_allele = out_ea,
  outcome_other_allele = out_oa,
  outcome_beta = as.numeric(matched$outcome_beta),
  outcome_se = as.numeric(matched$outcome_se),
  outcome_p = as.numeric(matched$outcome_p),
  outcome_eaf = as.numeric(matched$outcome_eaf),
  outcome_ncase = as.integer(matched$outcome_ncase),
  outcome_ncontrol = as.integer(matched$outcome_ncontrol),
  outcome_n_study = as.integer(matched$outcome_n_study),
  outcome_effect_scale = as.character(matched$outcome_effect_scale),
  compatibility_class = compatibility_class,
  valid_exposure_alleles = valid_exp,
  valid_outcome_alleles = valid_out,
  palindromic_snp = palindromic,
  planned_beta_flip = compatibility_class %in% c("swapped", "strand_swapped"),
  planned_eaf_flip = compatibility_class %in% c("swapped", "strand_swapped") & !is.na(as.numeric(matched$outcome_eaf)),
  planned_strand_transform = compatibility_class %in% c("strand_exact", "strand_swapped"),
  beta_flip_performed = FALSE,
  eaf_flip_performed = FALSE,
  strand_transform_performed = FALSE,
  formal_harmonisation_performed = FALSE,
  projected_final_valid = resolvable & !palindromic,
  exclusion_reason = ifelse(!valid_exp | !valid_out, "invalid_allele",
    ifelse(!resolvable, "incompatible_alleles",
    ifelse(palindromic, "palindromic_snp_excluded_by_forward_rule", ""))),
  stringsAsFactors = FALSE
)

included_ids <- unique(as.character(included_inst$resolved_rsid))
excluded_ids <- unique(as.character(excluded_inst$resolved_rsid))
missing_ids <- unique(as.character(missing$resolved_rsid))
union_ids <- sort(unique(c(included_ids, excluded_ids)))

count_set <- function(label, flag) {
  x <- preflight[flag, , drop = FALSE]
  data.frame(
    analysis_set = label,
    matched_input_count = nrow(x),
    exact_count = sum(x$compatibility_class == "exact"),
    swapped_count = sum(x$compatibility_class == "swapped"),
    strand_exact_count = sum(x$compatibility_class == "strand_exact"),
    strand_swapped_count = sum(x$compatibility_class == "strand_swapped"),
    palindromic_count = sum(x$palindromic_snp),
    incompatible_count = sum(x$compatibility_class == "incompatible"),
    invalid_count = sum(x$compatibility_class == "invalid"),
    projected_final_valid_count = sum(x$projected_final_valid),
    planned_beta_flip_count = sum(x$planned_beta_flip),
    planned_eaf_flip_count = sum(x$planned_eaf_flip),
    planned_strand_transform_count = sum(x$planned_strand_transform),
    stringsAsFactors = FALSE
  )
}
counts <- rbind(
  count_set("APOE_included", preflight$in_apoe_included_input),
  count_set("APOE_excluded", preflight$in_apoe_excluded_input),
  count_set("union", rep(TRUE, nrow(preflight)))
)

pal_table <- preflight[preflight$palindromic_snp, c("resolved_rsid", "source_marker_id", "in_apoe_included_input", "in_apoe_excluded_input", "exposure_effect_allele", "exposure_other_allele", "outcome_effect_allele", "outcome_other_allele", "compatibility_class", "projected_final_valid", "exclusion_reason"), drop = FALSE]
bad_table <- preflight[preflight$compatibility_class %in% c("incompatible", "invalid"), c("resolved_rsid", "source_marker_id", "in_apoe_included_input", "in_apoe_excluded_input", "exposure_effect_allele", "exposure_other_allele", "outcome_effect_allele", "outcome_other_allele", "compatibility_class", "projected_final_valid", "exclusion_reason"), drop = FALSE]

hard_checks <- list(
  instrument_freeze_gate = identical(instrument_freeze$freeze_status, "passed") && length(instrument_freeze$hard_check_failures) == 0L,
  outcome_extraction_freeze_gate = identical(outcome_freeze$freeze_status, "passed") && isTRUE(outcome_freeze$approved_for_chen_forward_harmonisation_preflight) && length(outcome_freeze$hard_check_failures) == 0L,
  decision_90_closure_gate = identical(outcome_closure$outcome_extraction_status, "passed") && length(outcome_closure$hard_check_failures) == 0L,
  analysis_role_sensitivity = identical(outcome_freeze$analysis_role, "forward_alternative_hb_gwas_sensitivity"),
  independent_replication_false = !isTRUE(instrument_freeze$independent_replication),
  matched_input_membership_valid = nrow(preflight) == outcome_freeze$union_exact_match_count &&
    sum(preflight$in_apoe_included_input) == outcome_freeze$included_exact_match_count &&
    sum(preflight$in_apoe_excluded_input) == outcome_freeze$excluded_exact_match_count &&
    all(preflight$resolved_rsid %in% union_ids),
  missing_outcomes_not_rescued = length(missing_ids) == outcome_freeze$union_missing_count &&
    !any(missing_ids %in% preflight$resolved_rsid),
  chen_effect_orientation_preserved = identical(chen_cert$allele_convention$effect_allele_field, "reference_allele") &&
    identical(chen_cert$allele_convention$other_allele_field, "other_allele") &&
    identical(chen_dict$reference_allele_documented_definition, "reference_allele - Effect allele") &&
    identical(chen_dict$other_allele_documented_definition, "other_allele - Non effect allele") &&
    all(exp_ea == toupper(as.character(imap$exposure_effect_allele))) &&
    all(exp_oa == toupper(as.character(imap$exposure_other_allele))),
  marker_tokens_not_used_as_effect_orientation = TRUE,
  reference_panel_orientation_not_used_as_effect_orientation = TRUE,
  finngen_effect_orientation_verified = identical(outcome_freeze$outcome_metadata$effect_allele, "alt") &&
    identical(outcome_freeze$outcome_metadata$other_allele, "ref") &&
    identical(outcome_freeze$outcome_metadata$beta, "beta") &&
    identical(outcome_freeze$outcome_metadata$eaf, "af_alt"),
  allele_classification_exhaustive = all(compatibility_class %in% c("exact", "swapped", "strand_exact", "strand_swapped", "incompatible", "invalid")) &&
    length(compatibility_class) == nrow(preflight),
  palindrome_definition_fixed = TRUE,
  forward_palindrome_rule_explicitly_adopted = is.finite(precedent$decision) && nrow(precedent) == 1L,
  no_eaf_palindrome_reinclusion = all(!preflight$projected_final_valid[preflight$palindromic_snp]),
  nonpalindromic_strand_logic_preserved = TRUE,
  no_frequency_based_alignment = TRUE,
  no_outcome_based_filtering = nrow(preflight) == outcome_freeze$union_exact_match_count,
  no_f_ge_30_filtering = nrow(preflight) == outcome_freeze$union_exact_match_count && instrument_conflict_free,
  f_statistic_preserved = same_num(preflight$exposure_F_stat, as.numeric(imap$F_stat)),
  planned_flips_auditable = all(!is.na(preflight$planned_beta_flip)) && all(!is.na(preflight$planned_eaf_flip)) && all(!is.na(preflight$planned_strand_transform)),
  no_proxy = !isTRUE(outcome_freeze$proxy_used),
  no_liftover = !isTRUE(outcome_freeze$liftover_used),
  no_coordinate_matching = !isTRUE(outcome_freeze$coordinate_matching_used),
  no_formal_harmonisation = !any(preflight$formal_harmonisation_performed),
  no_mr = TRUE,
  no_steiger = TRUE,
  chen_source_sha_valid = identical(chen_sha_before, expected_chen_sha),
  finngen_source_sha_valid = identical(finngen_sha_before, expected_finngen_sha),
  renv_lock_unchanged = identical(renv_before, hash_file(paths[["renv_lock"]]))
)
failures <- names(hard_checks)[!unlist(hard_checks)]
contract_status <- if (length(failures) == 0L) "frozen" else "failed"
preflight_status <- if (length(failures) == 0L) "passed" else "failed"
approved <- identical(preflight_status, "passed")

contract <- list(
  contract_version = "v1",
  analysis_direction = "Hb_to_delirium",
  analysis_role = "forward_alternative_hb_gwas_sensitivity",
  exposure = "Chen 2020 Hb",
  outcome = "FinnGen R13 F5_DELIRIUM",
  independent_replication = FALSE,
  source_instrument_freeze_decision = 87,
  source_outcome_extraction_freeze_decision = 91,
  technical_recovery_decision = 89,
  technical_readback_closure_decision = 90,
  exposure_effect_orientation_source = "Chen_reference_allele_other_allele",
  outcome_effect_orientation_source = "FinnGen_certified_effect_allele_fields",
  marker_identifier_alleles_used_for_effect_orientation = FALSE,
  reference_panel_alleles_used_for_effect_orientation = FALSE,
  harmonisation_classes = c("exact", "swapped", "strand_exact", "strand_swapped", "incompatible", "invalid"),
  palindromic_definition = "unordered allele set A/T or C/G",
  palindromic_rule = "exclude_all_without_EAF_reinclusion",
  palindromic_rule_precedent_decision = precedent$decision,
  palindromic_rule_precedent_path = precedent_rel,
  palindromic_rule_scope = "chen_forward_alternative_hb_gwas_sensitivity",
  palindromic_rule_timing = "adopted_before_formal_harmonisation_and_MR",
  palindromic_rule_result_agnostic = TRUE,
  eaf_based_palindrome_reinclusion = FALSE,
  proxy_allowed = FALSE,
  liftover_allowed = FALSE,
  coordinate_matching_allowed = FALSE,
  posthoc_instrument_filtering_allowed = FALSE,
  formal_harmonisation_performed = FALSE,
  mr_run = FALSE,
  steiger_run = FALSE,
  contract_status = contract_status,
  approved_for_formal_harmonisation = approved,
  hard_checks = hard_checks,
  hard_check_failures = failures
)

preflight_json <- list(
  preflight_version = "v1",
  analysis_direction = "Hb_to_delirium",
  analysis_role = "forward_alternative_hb_gwas_sensitivity",
  source_instrument_freeze_decision = 87,
  source_outcome_extraction_freeze_decision = 91,
  included_matched_input_count = counts$matched_input_count[counts$analysis_set == "APOE_included"],
  excluded_matched_input_count = counts$matched_input_count[counts$analysis_set == "APOE_excluded"],
  included_classification_counts = as.list(counts[counts$analysis_set == "APOE_included", c("exact_count", "swapped_count", "strand_exact_count", "strand_swapped_count", "incompatible_count", "invalid_count")]),
  excluded_classification_counts = as.list(counts[counts$analysis_set == "APOE_excluded", c("exact_count", "swapped_count", "strand_exact_count", "strand_swapped_count", "incompatible_count", "invalid_count")]),
  included_palindromic_count = counts$palindromic_count[counts$analysis_set == "APOE_included"],
  excluded_palindromic_count = counts$palindromic_count[counts$analysis_set == "APOE_excluded"],
  included_projected_final_valid_count = counts$projected_final_valid_count[counts$analysis_set == "APOE_included"],
  excluded_projected_final_valid_count = counts$projected_final_valid_count[counts$analysis_set == "APOE_excluded"],
  shared_projected_final_valid_count = sum(preflight$in_apoe_included_input & preflight$in_apoe_excluded_input & preflight$projected_final_valid),
  planned_beta_flip_count = counts$planned_beta_flip_count[counts$analysis_set == "union"],
  planned_eaf_flip_count = counts$planned_eaf_flip_count[counts$analysis_set == "union"],
  planned_strand_transform_count = counts$planned_strand_transform_count[counts$analysis_set == "union"],
  palindromic_rsids = sort(preflight$resolved_rsid[preflight$palindromic_snp]),
  incompatible_or_invalid_rsids = sort(preflight$resolved_rsid[preflight$compatibility_class %in% c("incompatible", "invalid")]),
  proxy_used = FALSE,
  liftover_used = FALSE,
  coordinate_matching_used = FALSE,
  formal_harmonisation_performed = FALSE,
  mr_run = FALSE,
  steiger_run = FALSE,
  preflight_status = preflight_status,
  approved_for_chen_forward_formal_harmonisation = approved,
  hard_checks = hard_checks,
  hard_check_failures = failures,
  informational_findings = list(
    missing_outcome_targets_not_in_harmonisation_preflight = TRUE,
    palindromic_variants_classified_and_projected_excluded = TRUE,
    incompatible_and_invalid_variants_retained_in_audit = TRUE,
    no_beta_or_eaf_values_rewritten = TRUE
  ),
  source_sha = list(chen_before = chen_sha_before, finngen_before = finngen_sha_before, chen_after = tolower(hash_file(paths[["chen_source"]])), finngen_after = tolower(hash_file(paths[["finngen_source"]]))),
  renv_lock_sha_before = renv_before,
  renv_lock_sha_after = hash_file(paths[["renv_lock"]])
)

atomic(paths[["counts_csv"]], function(p) write.csv(counts, p, row.names = FALSE, na = ""))
atomic(paths[["palindromic_tsv"]], function(p) write.table(pal_table, p, sep = "\t", quote = FALSE, row.names = FALSE, na = ""))
atomic(paths[["incompatible_tsv"]], function(p) write.table(bad_table, p, sep = "\t", quote = FALSE, row.names = FALSE, na = ""))

con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
atomic(paths[["preflight_parquet"]], function(p) {
  DBI::dbWriteTable(con, "chen_forward_harmonisation_preflight_master_v1", preflight, temporary = TRUE, overwrite = TRUE)
  DBI::dbExecute(con, sprintf("COPY %s TO %s (FORMAT PARQUET, COMPRESSION ZSTD)", sql_ident(con, "chen_forward_harmonisation_preflight_master_v1"), sql_string(p, must_work = FALSE)))
})
pq_ids <- DBI::dbGetQuery(con, sprintf("SELECT resolved_rsid FROM read_parquet(%s) ORDER BY resolved_rsid", sql_string(paths[["preflight_parquet"]])))
if (!identical(as.character(pq_ids$resolved_rsid), sort(as.character(preflight$resolved_rsid)))) {
  failures <- unique(c(failures, "preflight_parquet_readback_failed"))
  contract$hard_check_failures <- failures
  preflight_json$hard_check_failures <- failures
  contract$contract_status <- "failed"
  preflight_json$preflight_status <- "failed"
  contract$approved_for_formal_harmonisation <- FALSE
  preflight_json$approved_for_chen_forward_formal_harmonisation <- FALSE
}

atomic(paths[["contract_json"]], function(p) jsonlite::write_json(contract, p, pretty = TRUE, auto_unbox = TRUE, na = "null"))
atomic(paths[["preflight_json"]], function(p) jsonlite::write_json(preflight_json, p, pretty = TRUE, auto_unbox = TRUE, na = "null"))

decision_lines <- c(
  "# Decision 92 - Chen forward harmonisation contract and preflight V1",
  "",
  "Date: 2026-08-12",
  "Status: contract and preflight",
  "",
  "## Decision",
  "",
  "Freeze the Chen forward alternative-Hb-GWAS sensitivity harmonisation contract",
  "and execute a non-harmonising preflight classification using Decision 87 frozen",
  "Chen instruments and Decision 91 frozen FinnGen exact-rsID outcome extraction.",
  "",
  "## Orientation Authorities",
  "",
  "- Exposure effect allele: Chen `reference_allele`.",
  "- Exposure other allele: Chen `other_allele`.",
  "- Outcome effect allele: FinnGen `alt`.",
  "- Outcome other allele: FinnGen `ref`.",
  "- Marker identifier allele tokens and reference-panel alleles are not used for effect orientation.",
  "",
  "## Palindromic Rule",
  "",
  sprintf("- Rule: `exclude_all_without_EAF_reinclusion`."),
  sprintf("- Precedent Decision: `%s`.", precedent$decision),
  "- Scope: `chen_forward_alternative_hb_gwas_sensitivity`.",
  "- Timing: adopted before formal harmonisation and MR.",
  "",
  "## Preflight Results",
  "",
  sprintf("- contract_status: `%s`", contract$contract_status),
  sprintf("- preflight_status: `%s`", preflight_json$preflight_status),
  sprintf("- included matched input: `%d`", preflight_json$included_matched_input_count),
  sprintf("- excluded matched input: `%d`", preflight_json$excluded_matched_input_count),
  sprintf("- included projected final valid: `%d`", preflight_json$included_projected_final_valid_count),
  sprintf("- excluded projected final valid: `%d`", preflight_json$excluded_projected_final_valid_count),
  sprintf("- shared projected final valid: `%d`", preflight_json$shared_projected_final_valid_count),
  sprintf("- planned beta/eaf/strand-transform counts: `%d/%d/%d`", preflight_json$planned_beta_flip_count, preflight_json$planned_eaf_flip_count, preflight_json$planned_strand_transform_count),
  sprintf("- hard_check_failures: `%s`", paste(failures, collapse = ";")),
  "",
  "## Safeguards",
  "",
  "This stage did not run TwoSampleMR harmonise_data, formal harmonisation, MR,",
  "Steiger, proxy lookup, liftOver, coordinate matching, outcome-P filtering,",
  "frequency-based alignment, EAF-based palindromic reinclusion, re-clumping,",
  "or instrument reselection.",
  "",
  "## Expected Impact",
  "",
  "If passed, these outputs authorize only a future separately approved Chen",
  "forward formal harmonisation stage. They do not authorize MR."
)
atomic(paths[["decision"]], function(p) writeLines(decision_lines, p, useBytes = TRUE))

atomic(paths[["log"]], function(p) {
  writeLines(c(
    sprintf("[%s] stage=chen_forward_harmonisation_contract_preflight_v1", ts()),
    sprintf("[%s] contract_status=%s preflight_status=%s approved_for_formal_harmonisation=%s", ts(), contract$contract_status, preflight_json$preflight_status, preflight_json$approved_for_chen_forward_formal_harmonisation),
    sprintf("[%s] included_matched=%d excluded_matched=%d palindromic_union=%d incompatible_or_invalid_union=%d", ts(), preflight_json$included_matched_input_count, preflight_json$excluded_matched_input_count, length(preflight_json$palindromic_rsids), length(preflight_json$incompatible_or_invalid_rsids)),
    sprintf("[%s] hard_check_failures=%s", ts(), paste(failures, collapse = ";"))
  ), p, useBytes = TRUE)
})

stop_if(length(failures) != 0L, "Chen forward harmonisation preflight failed; QC retained.")
message("Chen forward harmonisation contract and preflight completed: ", preflight_status)
