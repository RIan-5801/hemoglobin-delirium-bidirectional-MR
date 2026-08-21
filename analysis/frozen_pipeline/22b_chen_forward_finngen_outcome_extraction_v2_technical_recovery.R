#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/22b_chen_forward_finngen_outcome_extraction_v2_technical_recovery.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

for (pkg in c("DBI", "duckdb", "jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
qpath <- function(path) gsub("'", "''", normalizePath(path, winslash = "/", mustWork = FALSE), fixed = TRUE)
strict_num <- function(x, label) {
  raw <- trimws(as.character(x))
  miss <- is.na(raw) | raw %in% c("", "NA", ".")
  y <- suppressWarnings(as.numeric(ifelse(miss, NA_character_, raw)))
  stop_if(any(!miss & is.na(y)), paste("Numeric conversion failed:", label))
  y
}
atomic <- function(path, writer) {
  partial <- paste0(path, ".partial")
  stop_if(file.exists(path) || file.exists(partial), paste("Output occupied:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writer(partial)
  stop_if(!file.exists(partial), paste("Writer did not create partial:", partial))
  stop_if(!file.rename(partial, path), paste("Atomic rename failed:", path))
}
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
write_tsv_precise <- function(x, path) {
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
write_csv_precise <- function(x, path) {
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  write.csv(x, path, row.names = FALSE, na = "")
}
validate_pair <- function(con, parquet_path, tsv_path, key = "resolved_rsid") {
  p <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", qpath(parquet_path)))
  t <- read.delim(tsv_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  same_cols <- identical(names(p), names(t))
  same_n <- nrow(p) == nrow(t)
  same_key_order <- same_n && identical(as.character(p[[key]]), as.character(t[[key]]))
  char_cols <- names(p)[vapply(p, function(z) is.character(z) || is.logical(z), logical(1))]
  num_cols <- names(p)[vapply(p, is.numeric, logical(1))]
  char_ok <- all(vapply(char_cols, function(k) identical(as.character(p[[k]]), as.character(t[[k]])), logical(1)))
  num_ok <- all(vapply(num_cols, function(k) {
    a <- as.numeric(p[[k]])
    b <- suppressWarnings(as.numeric(t[[k]]))
    same_na <- is.na(a) & is.na(b)
    finite_pair <- is.finite(a) & is.finite(b)
    abs_diff <- abs(a - b)
    rel_diff <- abs_diff / pmax(abs(a), abs(b), .Machine$double.xmin)
    all(same_na | (finite_pair & (abs_diff <= 1e-12 | rel_diff <= 1e-14)))
  }, logical(1)))
  list(row_count = nrow(p), same_cols = same_cols, same_n = same_n, same_key_order = same_key_order, char_ok = char_ok, numeric_with_tolerance_ok = num_ok)
}

build_union_query <- function(source_sql, target_rsids) {
  stop_if(length(target_rsids) < 1L || any(!grepl("^rs[0-9]+$", target_rsids)), "Target rsID list invalid for exact-token query.")
  vals <- paste(sprintf("'%s'", target_rsids), collapse = ", ")
  sprintf("WITH raw AS (
      SELECT row_number() OVER() AS scan_row_number, \"#chrom\" AS fg_chrom, pos AS fg_pos, ref AS fg_ref, alt AS fg_alt,
        rsids AS fg_rsids_raw, nearest_genes AS fg_nearest_genes, pval AS fg_pval_raw, mlogp AS fg_mlogp_raw,
        beta AS fg_beta_raw, sebeta AS fg_se_raw, af_alt AS fg_af_alt_raw, af_alt_cases AS fg_af_alt_cases_raw,
        af_alt_controls AS fg_af_alt_controls_raw, COUNT(*) OVER() AS observed_data_rows,
        SUM(CASE WHEN trim(coalesce(rsids,''))='' THEN 1 ELSE 0 END) OVER() AS missing_rsids_rows,
        SUM(CASE WHEN instr(coalesce(rsids,''),',')>0 AND instr(coalesce(rsids,''),';')=0 THEN 1 ELSE 0 END) OVER() AS comma_delimited_rows,
        SUM(CASE WHEN instr(coalesce(rsids,''),';')>0 AND instr(coalesce(rsids,''),',')=0 THEN 1 ELSE 0 END) OVER() AS semicolon_delimited_rows,
        SUM(CASE WHEN instr(coalesce(rsids,''),';')>0 AND instr(coalesce(rsids,''),',')>0 THEN 1 ELSE 0 END) OVER() AS mixed_delimiter_rows,
        SUM(CASE WHEN trim(coalesce(rsids,''))<>'' AND NOT regexp_full_match(trim(coalesce(rsids,'')),'^rs[0-9]+([[:space:]]*[,;][[:space:]]*rs[0-9]+)*$') THEN 1 ELSE 0 END) OVER() AS noncanonical_rsids_rows
      FROM %s
    ),
    tok AS (
      SELECT raw.*, trim(t.token) AS matched_rsid_token
      FROM raw
      CROSS JOIN UNNEST(string_split(replace(coalesce(fg_rsids_raw,''),';',','), ',')) AS t(token)
    )
    SELECT *
    FROM tok
    WHERE scan_row_number=1 OR matched_rsid_token IN (%s)", source_sql, vals)
}

fg <- file.path(root, "data_raw", "gwas", "finngen_R13_F5_DELIRIUM.gz")
cert_path <- file.path(root, "results", "qc", "finngen_R13_F5_DELIRIUM_input_certification_v1.json")
freeze_path <- file.path(root, "results", "qc", "chen_forward_instruments_v2_freeze.json")
first_attempt_union_targets <- file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_union_targets_v1.tsv")
first_attempt_log <- file.path(root, "results", "logs", "chen_forward_finngen_outcome_extraction_v1.log")
included_inst <- file.path(root, "data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_included_clumped_v2.parquet")
excluded_inst <- file.path(root, "data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_excluded_clumped_v2.parquet")
renv_lock <- file.path(root, "renv.lock")
out <- c(
  union_targets = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_union_targets_v2.tsv"),
  master_parquet = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_master_v2.parquet"),
  master_tsv = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_master_v2.tsv"),
  included_parquet = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_apoe_included_v2.parquet"),
  included_tsv = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_apoe_included_v2.tsv"),
  excluded_parquet = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_apoe_excluded_v2.parquet"),
  excluded_tsv = file.path(root, "data_derived", "forward_sensitivity_outcome", "chen_forward_finngen_outcome_apoe_excluded_v2.tsv"),
  match_audit = file.path(root, "results", "qc", "chen_forward_finngen_outcome_match_audit_v2.csv"),
  missing = file.path(root, "results", "qc", "chen_forward_finngen_outcome_missing_v2.tsv"),
  qc_json = file.path(root, "results", "qc", "chen_forward_finngen_outcome_extraction_v2.json"),
  log = file.path(root, "results", "logs", "chen_forward_finngen_outcome_extraction_v2.log"),
  decision = file.path(root, "docs", "decisions", "89_chen_forward_finngen_outcome_extraction_v2_technical_recovery_v1.1.md")
)
stop_if(any(file.exists(c(out, paste0(out, ".partial")))), "An outcome extraction V2 output or partial exists; refusing to overwrite.")
for (p in out) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = out[["log"]], append = TRUE)

write_decision89 <- function(summary) {
  lines <- c(
    "# Decision 89 - Chen forward FinnGen outcome extraction V2 technical recovery",
    "",
    "Date: 2026-08-12",
    "",
    "## Decision",
    "",
    "Execute technical recovery V2 for Chen forward FinnGen targeted outcome",
    "extraction after V1 stopped before publishing outcome results because DuckDB",
    "could not create a temporary directory under the R session temp path.",
    "",
    "## Rationale",
    "",
    "V2 preserves the Decision 88 scientific extraction design, starts again from",
    "the Decision 87 frozen Chen instrument authority and the certified FinnGen raw",
    "source, and uses a stable project-local DuckDB temp directory. V1's union",
    "target file and failure log are retained as technical provenance and are not",
    "used as authoritative MR inputs.",
    "",
    "## Results",
    "",
    sprintf("- FinnGen source SHA before: `%s`", summary$source_sha_before),
    sprintf("- FinnGen source SHA after: `%s`", summary$source_sha_after),
    sprintf("- Source rows scanned: `%s`", summary$source_rows_scanned),
    sprintf("- Runtime seconds: `%s`", summary$runtime_seconds),
    sprintf("- Included targets: `%s`", summary$included_target_count),
    sprintf("- Excluded targets: `%s`", summary$excluded_target_count),
    sprintf("- Union targets: `%s`", summary$union_target_count),
    sprintf("- Union exact / missing / multiple: `%s / %s / %s`", summary$union_unique_exact_match_count, summary$union_missing_count, summary$union_multiple_match_count),
    sprintf("- Included exact / missing: `%s / %s`", summary$included_exact_match_count, summary$included_missing_count),
    sprintf("- Excluded exact / missing: `%s / %s`", summary$excluded_exact_match_count, summary$excluded_missing_count),
    sprintf("- Outcome extraction status: `%s`", summary$outcome_extraction_status),
    sprintf("- Approved for Chen forward harmonisation design: `%s`", summary$approved_for_chen_forward_harmonisation_design),
    sprintf("- Hard-check failures: `%s`", paste(summary$hard_check_failures, collapse = ";")),
    "",
    "## Safeguards",
    "",
    "This stage does not perform harmonisation, MR, Steiger, proxy lookup, LD proxy",
    "lookup, coordinate matching, liftOver, nearest-variant matching, substring",
    "matching, instrument reselection, PLINK clumping, or outcome-based filtering.",
    "",
    "## Affected Files",
    "",
    "- `R/22b_chen_forward_finngen_outcome_extraction_v2_technical_recovery.R`",
    "- `data_derived/forward_sensitivity_outcome/chen_forward_finngen_union_targets_v2.tsv`",
    "- `data_derived/forward_sensitivity_outcome/chen_forward_finngen_outcome_master_v2.parquet`",
    "- `data_derived/forward_sensitivity_outcome/chen_forward_finngen_outcome_master_v2.tsv`",
    "- `data_derived/forward_sensitivity_outcome/chen_forward_finngen_outcome_apoe_included_v2.parquet`",
    "- `data_derived/forward_sensitivity_outcome/chen_forward_finngen_outcome_apoe_included_v2.tsv`",
    "- `data_derived/forward_sensitivity_outcome/chen_forward_finngen_outcome_apoe_excluded_v2.parquet`",
    "- `data_derived/forward_sensitivity_outcome/chen_forward_finngen_outcome_apoe_excluded_v2.tsv`",
    "- `results/qc/chen_forward_finngen_outcome_match_audit_v2.csv`",
    "- `results/qc/chen_forward_finngen_outcome_missing_v2.tsv`",
    "- `results/qc/chen_forward_finngen_outcome_extraction_v2.json`",
    "- `results/logs/chen_forward_finngen_outcome_extraction_v2.log`",
    "",
    "## Expected Impact",
    "",
    "This creates Chen forward sensitivity outcome files for a later, separately",
    "approved harmonisation design/preflight. It does not alter existing Vuckovic,",
    "reverse, instrument-selection, or MR results."
  )
  atomic(out[["decision"]], function(p) writeLines(lines, p, useBytes = TRUE))
}

main <- function() {
  started <- Sys.time()
  log_line("stage=chen_forward_finngen_outcome_extraction_v2_technical_recovery")
  for (p in c(fg, cert_path, freeze_path, first_attempt_union_targets, first_attempt_log, included_inst, excluded_inst, renv_lock)) {
    stop_if(!file.exists(p), paste("Missing required input:", p))
  }
  cert <- jsonlite::fromJSON(cert_path, simplifyVector = FALSE)
  freeze <- jsonlite::fromJSON(freeze_path, simplifyVector = FALSE)
  stop_if(!identical(cert$certification_status, "passed") || !isTRUE(cert$header_exact_match), "FinnGen certification gate failed.")
  expected_fg_sha <- tolower(cert$certified_input_sha256)
  source_sha_before <- tolower(hash_file(fg))
  stop_if(!identical(source_sha_before, expected_fg_sha), "FinnGen source SHA before scan differs from certification.")
  stop_if(!identical(freeze$freeze_status, "passed") || !isTRUE(freeze$approved_for_chen_forward_finngen_outcome_extraction) || length(freeze$hard_check_failures) != 0L, "Chen instrument freeze gate failed.")
  renv_before <- hash_file(renv_lock)

  duckdb_tmp <- file.path(root, "results", "tmp", "duckdb_chen_forward_finngen_outcome_v2")
  dir.create(duckdb_tmp, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(DUCKDB_TEMP_DIRECTORY = duckdb_tmp)
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:", config = list(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(con, "PRAGMA threads=8")
  DBI::dbExecute(con, "PRAGMA memory_limit='8GB'")
  DBI::dbExecute(con, sprintf("PRAGMA temp_directory='%s'", qpath(duckdb_tmp)))
  included <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", qpath(included_inst)))
  excluded <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", qpath(excluded_inst)))
  stop_if(anyDuplicated(included$resolved_rsid) || anyDuplicated(excluded$resolved_rsid), "Duplicate resolved_rsid in frozen instruments.")

  mk_targets <- function(x, in_inc, in_exc) {
    data.frame(
      resolved_rsid = x$resolved_rsid,
      source_marker_id = x$source_marker_id,
      in_apoe_included_input = in_inc,
      in_apoe_excluded_input = in_exc,
      exposure_effect_allele = x$exposure_effect_allele,
      exposure_other_allele = x$exposure_other_allele,
      exposure_beta = x$beta,
      exposure_se = x$se,
      exposure_pval = x$pval,
      exposure_eaf = x$eaf,
      exposure_n_samples = x$n_samples,
      exposure_F_stat = x$F_stat,
      exposure_marker_chr_grch37 = x$marker_chr,
      exposure_marker_pos_grch37 = x$marker_pos,
      exposure_reference_chr_grch37 = x$reference_chr_grch37,
      exposure_reference_pos_grch37 = x$reference_pos_grch37,
      exposure_apoe_region = x$apoe_region,
      stringsAsFactors = FALSE
    )
  }
  targets <- mk_targets(included, TRUE, FALSE)
  targets$in_apoe_excluded_input <- targets$resolved_rsid %in% excluded$resolved_rsid
  exc_only <- excluded[!excluded$resolved_rsid %in% included$resolved_rsid, , drop = FALSE]
  if (nrow(exc_only) > 0L) targets <- rbind(targets, mk_targets(exc_only, FALSE, TRUE))
  targets <- targets[order(targets$resolved_rsid), , drop = FALSE]
  rownames(targets) <- NULL
  stop_if(anyDuplicated(targets$resolved_rsid), "Union target set contains duplicate rsIDs.")
  included_target_count <- nrow(included)
  excluded_target_count <- nrow(excluded)
  shared_target_count <- length(intersect(included$resolved_rsid, excluded$resolved_rsid))
  union_target_count <- nrow(targets)
  stop_if(included_target_count != as.integer(freeze$included_nsnp) || excluded_target_count != as.integer(freeze$excluded_nsnp), "Target counts differ from freeze.")
  atomic(out[["union_targets"]], function(p) write_tsv_precise(targets, p))

  header <- c("#chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls")
  hcon <- gzfile(fg, open = "rt", encoding = "UTF-8")
  header_line <- tryCatch(readLines(hcon, n = 1L, warn = FALSE), finally = close(hcon))
  observed_header <- strsplit(header_line, "\t", fixed = TRUE)[[1L]]
  stop_if(!identical(observed_header, header), "Raw FinnGen header gate failed.")
  schema <- paste(sprintf("'%s': 'VARCHAR'", header), collapse = ", ")
  source_sql <- sprintf("read_csv('%s', header=true, delim='\\t', columns={%s}, auto_detect=false, strict_mode=true)", qpath(fg), schema)
  query <- build_union_query(source_sql, targets$resolved_rsid)
  plan <- DBI::dbGetQuery(con, paste("EXPLAIN", query))
  explain_plan_sha256 <- digest::digest(paste(unlist(plan), collapse = "\n"), algo = "sha256", serialize = FALSE)
  log_line("DuckDB_EXPLAIN_STATUS=passed explain_plan_sha256=", explain_plan_sha256)
  scan <- DBI::dbGetQuery(con, query)
  sentinel <- scan[scan$scan_row_number == 1L, , drop = FALSE]
  stop_if(nrow(sentinel) != 1L, "Scan sentinel row count is not one.")
  observed_rows <- as.numeric(sentinel$observed_data_rows[[1L]])
  stop_if(observed_rows != as.numeric(cert$certified_data_rows), "FinnGen observed data rows differ from certification.")
  detail <- scan[scan$scan_row_number != 1L & scan$matched_rsid_token %in% targets$resolved_rsid, , drop = FALSE]
  stop_if(any(!(detail$matched_rsid_token %in% targets$resolved_rsid)), "Returned non-target detail row detected.")

  target_counts <- as.data.frame(table(detail$matched_rsid_token), stringsAsFactors = FALSE)
  names(target_counts) <- c("resolved_rsid", "finngen_match_count")
  targets2 <- merge(targets, target_counts, by = "resolved_rsid", all.x = TRUE, sort = FALSE)
  targets2$finngen_match_count[is.na(targets2$finngen_match_count)] <- 0L
  targets2$finngen_match_status <- ifelse(targets2$finngen_match_count == 0L, "missing", ifelse(targets2$finngen_match_count == 1L, "unique_exact_match", "multiple_source_matches"))
  multiple <- targets2[targets2$finngen_match_status == "multiple_source_matches", , drop = FALSE]
  stop_if(nrow(multiple) != 0L, "Multiple FinnGen source matches detected; stopping.")
  matched <- merge(targets2[targets2$finngen_match_status == "unique_exact_match", , drop = FALSE], detail, by.x = "resolved_rsid", by.y = "matched_rsid_token", all.x = TRUE, sort = FALSE)
  stop_if(nrow(matched) != sum(targets2$finngen_match_status == "unique_exact_match"), "Matched master row count differs from target classification.")

  master <- data.frame(
    resolved_rsid = matched$resolved_rsid,
    source_marker_id = matched$source_marker_id,
    in_apoe_included_input = matched$in_apoe_included_input,
    in_apoe_excluded_input = matched$in_apoe_excluded_input,
    exposure_effect_allele = matched$exposure_effect_allele,
    exposure_other_allele = matched$exposure_other_allele,
    exposure_beta = matched$exposure_beta,
    exposure_se = matched$exposure_se,
    exposure_pval = matched$exposure_pval,
    exposure_eaf = matched$exposure_eaf,
    exposure_n_samples = matched$exposure_n_samples,
    exposure_F_stat = matched$exposure_F_stat,
    exposure_marker_chr_grch37 = matched$exposure_marker_chr_grch37,
    exposure_marker_pos_grch37 = matched$exposure_marker_pos_grch37,
    outcome_source = "FinnGen_R13_F5_DELIRIUM",
    outcome_trait = "delirium",
    outcome_build = "GRCh38",
    outcome_rsid = matched$resolved_rsid,
    outcome_chr_grch38 = strict_num(matched$fg_chrom, "outcome_chr_grch38"),
    outcome_pos_grch38 = strict_num(matched$fg_pos, "outcome_pos_grch38"),
    outcome_effect_allele = toupper(matched$fg_alt),
    outcome_other_allele = toupper(matched$fg_ref),
    outcome_beta = strict_num(matched$fg_beta_raw, "outcome_beta"),
    outcome_se = strict_num(matched$fg_se_raw, "outcome_se"),
    outcome_p = strict_num(matched$fg_pval_raw, "outcome_p"),
    outcome_log10p = strict_num(matched$fg_mlogp_raw, "outcome_log10p"),
    outcome_eaf = strict_num(matched$fg_af_alt_raw, "outcome_eaf"),
    outcome_af_alt_cases = strict_num(matched$fg_af_alt_cases_raw, "outcome_af_alt_cases"),
    outcome_af_alt_controls = strict_num(matched$fg_af_alt_controls_raw, "outcome_af_alt_controls"),
    outcome_rsids_raw = matched$fg_rsids_raw,
    outcome_nearest_genes = matched$fg_nearest_genes,
    outcome_ncase = 5121L,
    outcome_ncontrol = 465023L,
    outcome_n_study = 470144L,
    outcome_effect_scale = "log(OR) for delirium per ALT allele",
    outcome_match_status = "unique_exact_match",
    stringsAsFactors = FALSE
  )
  master <- master[order(master$resolved_rsid), , drop = FALSE]
  included_out <- master[master$in_apoe_included_input, , drop = FALSE]
  excluded_out <- master[master$in_apoe_excluded_input, , drop = FALSE]
  missing <- targets2[targets2$finngen_match_status == "missing", , drop = FALSE]
  match_audit <- data.frame(
    analysis_set = c("union", "included", "excluded"),
    target_count = c(union_target_count, included_target_count, excluded_target_count),
    exact_match_count = c(nrow(master), nrow(included_out), nrow(excluded_out)),
    missing_count = c(nrow(missing), sum(targets2$resolved_rsid %in% included$resolved_rsid & targets2$finngen_match_status == "missing"), sum(targets2$resolved_rsid %in% excluded$resolved_rsid & targets2$finngen_match_status == "missing")),
    multiple_match_count = c(0L, 0L, 0L),
    stringsAsFactors = FALSE
  )

  write_pair <- function(df, parquet_path, tsv_path) {
    atomic(parquet_path, function(p) {
      DBI::dbWriteTable(con, "write_df", df, overwrite = TRUE)
      DBI::dbExecute(con, sprintf("COPY write_df TO '%s' (FORMAT PARQUET)", qpath(p)))
      DBI::dbExecute(con, "DROP TABLE write_df")
    })
    atomic(tsv_path, function(p) write_tsv_precise(df, p))
  }
  write_pair(master, out[["master_parquet"]], out[["master_tsv"]])
  write_pair(included_out, out[["included_parquet"]], out[["included_tsv"]])
  write_pair(excluded_out, out[["excluded_parquet"]], out[["excluded_tsv"]])
  atomic(out[["match_audit"]], function(p) write_csv_precise(match_audit, p))
  atomic(out[["missing"]], function(p) write_tsv_precise(missing, p))
  rt_master <- validate_pair(con, out[["master_parquet"]], out[["master_tsv"]])
  rt_included <- validate_pair(con, out[["included_parquet"]], out[["included_tsv"]])
  rt_excluded <- validate_pair(con, out[["excluded_parquet"]], out[["excluded_tsv"]])

  source_sha_after <- tolower(hash_file(fg))
  renv_after <- hash_file(renv_lock)
  hard_checks <- list(
    chen_instrument_freeze_gate = identical(freeze$freeze_status, "passed") && isTRUE(freeze$approved_for_chen_forward_finngen_outcome_extraction) && length(freeze$hard_check_failures) == 0L,
    finngen_source_certification_gate = identical(cert$certification_status, "passed") && isTRUE(cert$header_exact_match),
    source_sha_before_gate = identical(source_sha_before, expected_fg_sha),
    source_sha_after_gate = identical(source_sha_after, expected_fg_sha),
    source_unchanged = identical(source_sha_before, source_sha_after),
    v1_failure_evidence_preserved = file.exists(first_attempt_union_targets) && file.exists(first_attempt_log),
    v2_started_from_freeze_not_v1_targets = TRUE,
    union_target_set_valid = union_target_count == length(unique(targets$resolved_rsid)) && included_target_count > 0L && excluded_target_count > 0L,
    exact_rsid_matching_only = TRUE,
    no_substring_matching = TRUE,
    no_proxy = TRUE,
    no_liftover = TRUE,
    no_coordinate_matching = TRUE,
    cross_build_coordinate_comparison_used_false = TRUE,
    source_scan_complete = observed_rows == as.numeric(cert$certified_data_rows),
    no_unresolved_multiple_matches = nrow(multiple) == 0L,
    at_least_one_exact_match = nrow(master) > 0L,
    included_membership_preserved = setequal(included_out$resolved_rsid, intersect(master$resolved_rsid, included$resolved_rsid)),
    excluded_membership_preserved = setequal(excluded_out$resolved_rsid, intersect(master$resolved_rsid, excluded$resolved_rsid)),
    no_outcome_based_filtering = TRUE,
    master_parquet_tsv_consistency = all(unlist(rt_master[c("same_cols", "same_n", "same_key_order", "char_ok", "numeric_with_tolerance_ok")])),
    included_parquet_tsv_consistency = all(unlist(rt_included[c("same_cols", "same_n", "same_key_order", "char_ok", "numeric_with_tolerance_ok")])),
    excluded_parquet_tsv_consistency = all(unlist(rt_excluded[c("same_cols", "same_n", "same_key_order", "char_ok", "numeric_with_tolerance_ok")])),
    no_harmonisation = TRUE,
    no_mr = TRUE,
    no_steiger = TRUE,
    renv_lock_unchanged = identical(renv_before, renv_after)
  )
  failures <- names(hard_checks)[!unlist(hard_checks)]
  status <- if (length(failures) == 0L) "passed" else "failed"
  approved <- identical(status, "passed")
  qc <- list(
    extraction_version = "v2",
    decision = 89,
    technical_recovery_for_decision = 88,
    analysis_direction = "Hb_to_delirium",
    analysis_role = "forward_alternative_hb_gwas_sensitivity",
    source_instrument_freeze_decision = 87,
    outcome_source = "FinnGen_R13_F5_DELIRIUM",
    outcome_build = "GRCh38",
    outcome_metadata = list(effect_allele = "alt", other_allele = "ref", beta = "beta", se = "sebeta", p = "pval", eaf = "af_alt", ncase = 5121L, ncontrol = 465023L, n_study = 470144L, effect_scale = "log(OR) for delirium per ALT allele"),
    matching_method = "exact_canonical_rsid_token",
    proxy_used = FALSE,
    liftover_used = FALSE,
    coordinate_matching_used = FALSE,
    cross_build_coordinate_comparison_used = FALSE,
    included_target_count = included_target_count,
    excluded_target_count = excluded_target_count,
    shared_target_count = shared_target_count,
    union_target_count = union_target_count,
    union_unique_exact_match_count = nrow(master),
    union_missing_count = nrow(missing),
    union_multiple_match_count = nrow(multiple),
    included_exact_match_count = nrow(included_out),
    included_missing_count = sum(targets2$resolved_rsid %in% included$resolved_rsid & targets2$finngen_match_status == "missing"),
    excluded_exact_match_count = nrow(excluded_out),
    excluded_missing_count = sum(targets2$resolved_rsid %in% excluded$resolved_rsid & targets2$finngen_match_status == "missing"),
    source_scan_completed = TRUE,
    source_rows_scanned = as.integer(observed_rows),
    runtime_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    source_sha_before = source_sha_before,
    source_sha_after = source_sha_after,
    harmonisation_performed = FALSE,
    mr_run = FALSE,
    steiger_run = FALSE,
    parquet_tsv_consistency = list(master = rt_master, included = rt_included, excluded = rt_excluded),
    rsids_audit = list(
      missing_rsids_rows = as.numeric(sentinel$missing_rsids_rows[[1L]]),
      comma_delimited_rows = as.numeric(sentinel$comma_delimited_rows[[1L]]),
      semicolon_delimited_rows = as.numeric(sentinel$semicolon_delimited_rows[[1L]]),
      mixed_delimiter_rows = as.numeric(sentinel$mixed_delimiter_rows[[1L]]),
      noncanonical_rsids_rows = as.numeric(sentinel$noncanonical_rsids_rows[[1L]])
    ),
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after,
    renv_lock_unchanged = identical(renv_before, renv_after),
    output_sha256 = as.list(vapply(out[names(out) %in% c("union_targets", "master_parquet", "master_tsv", "included_parquet", "included_tsv", "excluded_parquet", "excluded_tsv", "match_audit", "missing")], hash_file, character(1))),
    recovery = list(
      v1_failure_reason = "DuckDB temporary directory creation failed before outcome outputs were published",
      stable_project_temp_directory_used = duckdb_tmp,
      scientific_design_changed = FALSE,
      raw_gwas_modified = FALSE
    ),
    outcome_extraction_status = status,
    approved_for_chen_forward_harmonisation_design = approved,
    hard_checks = hard_checks,
    hard_check_failures = failures,
    informational_findings = list(
      missing_outcome_targets_allowed_without_proxy = TRUE,
      outcome_based_filtering_not_performed = TRUE,
      subsets_derived_from_same_union_master = TRUE,
      v1_union_targets_retained_as_technical_provenance_non_authoritative = TRUE
    )
  )
  atomic(out[["qc_json"]], function(p) jsonlite::write_json(qc, p, pretty = TRUE, auto_unbox = TRUE, na = "null"))
  write_decision89(qc)
  log_line("outcome_extraction_status=", status, " union_targets=", union_target_count, " union_exact=", nrow(master), " union_missing=", nrow(missing), " hard_check_failures=", paste(failures, collapse = ";"))
  stop_if(!identical(status, "passed"), "Chen forward FinnGen outcome extraction V2 failed; QC retained.")
}

tryCatch(main(), error = function(e) {
  log_line("status=failed error=", conditionMessage(e))
  quit(status = 1L)
})
