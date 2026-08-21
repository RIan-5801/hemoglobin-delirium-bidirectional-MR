#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "renv", "activate.R"))
required <- c("DBI", "duckdb", "jsonlite", "digest")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required project packages: ", paste(missing_packages, collapse = ", "))

input <- file.path(project_root, "data_raw", "gwas", "ebi-a-GCST90002384.vcf.gz")
instrument_dir <- file.path(project_root, "data_derived", "instruments")
quarantine_dir <- file.path(instrument_dir, "quarantine")
qc_dir <- file.path(project_root, "results", "qc")
log_dir <- file.path(project_root, "results", "logs")
dir.create(instrument_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(quarantine_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

candidate_parquet <- file.path(instrument_dir, "vuckovic_hb_p5e8_candidates_v3.parquet")
candidate_tsv <- file.path(instrument_dir, "vuckovic_hb_p5e8_candidates_v3.tsv.gz")
quarantine_parquet <- file.path(quarantine_dir, "vuckovic_hb_p5e8_unresolved_v3.parquet")
qc_json <- file.path(qc_dir, "vuckovic_hb_fullscan_v3_qc.json")
format_csv <- file.path(qc_dir, "vuckovic_hb_format_patterns_v3.csv")
duplicates_csv <- file.path(qc_dir, "vuckovic_hb_candidate_duplicates_v3.csv")
by_chr_csv <- file.path(qc_dir, "vuckovic_hb_candidates_by_chr_v3.csv")
log_out <- file.path(log_dir, "vuckovic_hb_candidate_selection_v3.log")
targets <- c(candidate_parquet, candidate_tsv, quarantine_parquet, qc_json, format_csv, duplicates_csv, by_chr_csv, log_out)

if (!file.exists(input)) stop("Input VCF not found: ", input)
partial_targets <- paste0(targets, ".partial")
if (any(file.exists(targets))) stop("One or more final V3 target files already exist; refusing to overwrite.")
if (any(file.exists(partial_targets))) stop("One or more V3 partial files already exist; refusing to overwrite.")

log_line <- function(...) {
  text <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), " | ", paste(..., collapse = ""))
  cat(text, "\n", file = log_out, append = TRUE)
  message(text)
}

stage_start <- function(name) {
  started <- Sys.time()
  log_line("STAGE_START=", name, "; time=", format(started, "%Y-%m-%dT%H:%M:%S%z"))
  started
}

stage_end <- function(name, started) {
  finished <- Sys.time()
  log_line("STAGE_END=", name, "; time=", format(finished, "%Y-%m-%dT%H:%M:%S%z"),
           "; duration_seconds=", as.numeric(difftime(finished, started, units = "secs")))
  invisible(finished)
}

log_final <- log_out
log_out <- paste0(log_final, ".partial")
options(error = function() {
  err <- geterrmessage()
  try(log_line("ERROR stage termination: ", err), silent = TRUE)
  quit(save = "no", status = 1, runLast = FALSE)
})

find_vcf_header <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  i <- 0L
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("#CHROM header was not found in VCF input.")
    i <- i + 1L
    if (startsWith(line, "#CHROM")) return(list(line_number = i, fields = strsplit(line, "\t", fixed = TRUE)[[1]]))
  }
}

sql_path <- function(conn, path) as.character(DBI::dbQuoteString(conn, normalizePath(path, winslash = "/", mustWork = FALSE)))
sql_ident <- function(conn, name) as.character(DBI::dbQuoteIdentifier(conn, name))
df_to_records <- function(x) {
  if (!is.data.frame(x)) stop("JSON conversion expected a data.frame.")
  lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
}
range_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(list(min = NA_real_, max = NA_real_))
  list(min = min(x), max = max(x))
}
elapsed_seconds <- function(start, end) as.numeric(difftime(end, start, units = "secs"))

start_time <- Sys.time()
sha_stage <- stage_start("sha_validation")
expected_sha256 <- "C68C98C7800C59D9D64CB88739E2E245F8AAD6E4E886455CBF7424661AFC3D41"
calculated_sha256_before <- toupper(digest::digest(input, algo = "sha256", file = TRUE))
sha256_before_matches <- identical(calculated_sha256_before, expected_sha256)
log_line("SHA-256 before scan calculated=", calculated_sha256_before, "; expected=", expected_sha256, "; match=", sha256_before_matches)
if (!sha256_before_matches) stop("Input SHA-256 mismatch. No VCF parsing was performed.")
stage_end("sha_validation", sha_stage)

file_read_stage <- stage_start("file_read")
header <- find_vcf_header(input)
required_vcf_columns <- c("#CHROM", "POS", "ID", "REF", "ALT", "FORMAT")
if (!all(required_vcf_columns %in% header$fields)) stop("VCF #CHROM header lacks required columns.")
sample_columns_header <- header$fields[(match("FORMAT", header$fields) + 1L):length(header$fields)]
if (length(sample_columns_header) != 1L) stop("Expected exactly one VCF sample column in header; found: ", paste(sample_columns_header, collapse = ", "))
log_line("Detected #CHROM header at physical line ", header$line_number, "; sample column=", sample_columns_header)

conn <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:", read_only = FALSE)
on.exit(DBI::dbDisconnect(conn, shutdown = TRUE), add = TRUE)
DBI::dbExecute(conn, "PRAGMA threads=8;")
DBI::dbExecute(conn, "PRAGMA memory_limit='8GB';")
log_line("DuckDB resource limits set: threads=8; memory_limit=8GB.")
input_sql <- sql_path(conn, input)
scan_sql <- paste0(
  "read_csv_auto(", input_sql,
  ", delim = '\\t', header = true, skip = ", header$line_number - 1L,
  ", compression = 'gzip', all_varchar = true, sample_size = 10000, ignore_errors = false)"
)
probe <- tryCatch(DBI::dbGetQuery(conn, paste0("SELECT * FROM ", scan_sql, " LIMIT 10")), error = function(e) e)
if (inherits(probe, "error")) {
  log_line("DuckDB direct BGZF/gzip probe FAILED: ", conditionMessage(probe))
  stop("DuckDB cannot directly scan this BGZF/gzip VCF. No fallback was used.")
}
if (nrow(probe) != 10L || !all(required_vcf_columns %in% names(probe))) {
  log_line("DuckDB direct BGZF/gzip probe schema FAILED. Rows=", nrow(probe), "; columns=", paste(names(probe), collapse = ","))
  stop("DuckDB scan did not return the expected VCF schema. No fallback was used.")
}
sample_columns <- names(probe)[(match("FORMAT", names(probe)) + 1L):ncol(probe)]
if (length(sample_columns) != 1L) stop("DuckDB scan did not return exactly one FORMAT sample column.")
sample_col_sql <- sql_ident(conn, sample_columns[[1]])
log_line("DuckDB direct BGZF/gzip probe succeeded.")
stage_end("file_read", file_read_stage)

format_unit_stage <- stage_start("format_audit_100k_unit_test")
format_summary_sql <- function(limit_clause = "") {
  paste0(
    "WITH format_base AS (",
    " SELECT COALESCE(NULLIF(\"FORMAT\", ''), '<MISSING>') AS format_value",
    " FROM ", scan_sql, limit_clause,
    " )",
    " SELECT format_value AS format_pattern, COUNT(*)::BIGINT AS n_records,",
    " strpos(':' || format_value || ':', ':ES:') > 0 AS has_es,",
    " strpos(':' || format_value || ':', ':SE:') > 0 AS has_se,",
    " strpos(':' || format_value || ':', ':LP:') > 0 AS has_lp,",
    " strpos(':' || format_value || ':', ':AF:') > 0 AS has_af",
    " FROM format_base",
    " GROUP BY format_value",
    " ORDER BY n_records DESC"
  )
}
format_unit_sql <- format_summary_sql(" LIMIT 100000")
format_unit <- tryCatch(DBI::dbGetQuery(conn, format_unit_sql), error = function(e) e)
if (inherits(format_unit, "error")) {
  log_line("FORMAT SQL unit test FAILED. SQL=", format_unit_sql, "; error=", conditionMessage(format_unit))
  stop("FORMAT SQL unit test failed; no full scan was started.")
}
unit_count <- sum(format_unit$n_records)
unit_expected <- format_unit$format_pattern == "ES:SE:LP:AF:ID"
unit_core <- format_unit$has_es & format_unit$has_se & format_unit$has_lp & format_unit$has_af
if (!identical(as.numeric(unit_count), 100000) || !any(unit_expected) || !all(unit_core)) {
  log_line("FORMAT SQL unit test FAILED validation. SQL=", format_unit_sql,
           "; rows=", nrow(format_unit), "; n_records_sum=", unit_count,
           "; expected_pattern_present=", any(unit_expected), "; all_patterns_have_core_keys=", all(unit_core))
  stop("FORMAT SQL unit test validation failed; no full scan was started.")
}
log_line("FORMAT SQL unit test PASSED. rows=", nrow(format_unit), "; n_records_sum=", unit_count,
         "; expected_pattern_present=", any(unit_expected), "; all_patterns_have_core_keys=", all(unit_core))


stage_end("format_audit_100k_unit_test", format_unit_stage)

full_read_stage <- stage_start("file_read_full_count")
record_count <- DBI::dbGetQuery(conn, paste0("SELECT count(*) AS n FROM ", scan_sql))$n[[1]]
expected_records <- 40262327
log_line("Full scan record count=", record_count, "; expected=", expected_records)
if (!identical(as.numeric(record_count), as.numeric(expected_records))) {
  stop("Full VCF record count differs from the approved expected value; no candidate output was created.")
}
stage_end("file_read_full_count", full_read_stage)

common_cte <- paste0(
"WITH raw AS (
  SELECT row_number() OVER () AS source_row_number,
         \"#CHROM\" AS chr_raw, \"POS\" AS pos_raw, \"ID\" AS rsid_raw,
         \"REF\" AS ref_raw, \"ALT\" AS alt_raw, \"FORMAT\" AS format_pattern,
         ", sample_col_sql, " AS sample_value
  FROM ", scan_sql, "
), tokens AS (
  SELECT *, string_split(format_pattern, ':') AS format_keys, string_split(sample_value, ':') AS format_values
  FROM raw
), parsed AS (
  SELECT *,
    list_position(format_keys, 'ES') AS es_key_position,
    list_position(format_keys, 'SE') AS se_key_position,
    list_position(format_keys, 'LP') AS lp_key_position,
    list_position(format_keys, 'AF') AS af_key_position,
    try_cast(list_extract(format_values, list_position(format_keys, 'ES')) AS DOUBLE) AS beta,
    try_cast(list_extract(format_values, list_position(format_keys, 'SE')) AS DOUBLE) AS se,
    try_cast(list_extract(format_values, list_position(format_keys, 'LP')) AS DOUBLE) AS lp,
    try_cast(list_extract(format_values, list_position(format_keys, 'AF')) AS DOUBLE) AS eaf,
    try_cast(pos_raw AS BIGINT) AS pos,
    upper(ref_raw) AS ref_upper,
    upper(alt_raw) AS alt_upper
  FROM tokens
), base AS (
  SELECT *,
    (pos IS NOT NULL AND pos > 0) AS position_valid,
    (rsid_raw IS NULL OR rsid_raw = '' OR rsid_raw = '.') AS rsid_missing,
    (ref_raw IS NULL OR alt_raw IS NULL OR ref_raw = '' OR alt_raw = '' OR ref_raw = '.' OR alt_raw = '.') AS allele_missing,
    (alt_raw IS NOT NULL AND strpos(alt_raw, ',') > 0) AS multi_allelic,
    (ref_raw IS NOT NULL AND alt_raw IS NOT NULL AND ref_raw <> '' AND alt_raw <> '' AND ref_raw <> '.' AND alt_raw <> '.' AND
      (NOT regexp_full_match(ref_upper, '^[ACGT]+$') OR NOT regexp_full_match(alt_upper, '^[ACGT,]+$'))) AS allele_invalid,
    (ref_raw IS NOT NULL AND alt_raw IS NOT NULL AND strpos(alt_raw, ',') = 0 AND length(ref_upper) <> length(alt_upper)) AS indel,
    coalesce(NOT (format_pattern IN ('ES:SE:LP:AF:ID', 'ES:SE:LP:AF')), true) AS format_pattern_unexpected,
    (beta IS NULL OR NOT isfinite(beta)) AS beta_nonfinite,
    (se IS NULL OR NOT isfinite(se) OR se <= 0) AS se_nonpositive,
    (lp IS NULL OR NOT isfinite(lp) OR lp < 0) AS lp_invalid,
    (eaf IS NOT NULL AND (NOT isfinite(eaf) OR eaf < 0 OR eaf > 1)) AS eaf_out_of_range,
    (chr_raw NOT IN ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22')) AS non_autosomal_contig,
    CASE WHEN pos IS NOT NULL AND pos > 0 AND chr_raw IS NOT NULL AND ref_raw IS NOT NULL AND alt_raw IS NOT NULL
         THEN concat(chr_raw, ':', pos, ':', ref_raw, ':', alt_raw) ELSE NULL END AS variant_id_source,
    CASE WHEN pos IS NOT NULL AND pos > 0 AND chr_raw IS NOT NULL AND ref_upper IS NOT NULL AND alt_upper IS NOT NULL
         THEN CASE WHEN ref_upper <= alt_upper THEN concat(chr_raw, ':', pos, ':', ref_upper, ':', alt_upper)
                   ELSE concat(chr_raw, ':', pos, ':', alt_upper, ':', ref_upper) END
         ELSE NULL END AS allele_set_key
  FROM parsed
) ")

format_audit_stage <- stage_start("format_audit")
format_df <- DBI::dbGetQuery(conn, format_summary_sql())
utils::write.csv(format_df, paste0(format_csv, ".partial"), row.names = FALSE)
log_line("FORMAT pattern audit completed; patterns=", nrow(format_df))
stage_end("format_audit", format_audit_stage)

global_qc_stage <- stage_start("global_qc_aggregation")
full_qc <- DBI::dbGetQuery(conn, paste0(common_cte, "
SELECT
  count(*) AS records,
  sum(CASE WHEN NOT position_valid THEN 1 ELSE 0 END) AS position_invalid,
  sum(CASE WHEN rsid_missing THEN 1 ELSE 0 END) AS rsid_missing,
  sum(CASE WHEN allele_missing THEN 1 ELSE 0 END) AS allele_missing,
  sum(CASE WHEN allele_invalid THEN 1 ELSE 0 END) AS allele_invalid,
  sum(CASE WHEN indel THEN 1 ELSE 0 END) AS indel,
  sum(CASE WHEN multi_allelic THEN 1 ELSE 0 END) AS multi_allelic,
  sum(CASE WHEN format_pattern_unexpected THEN 1 ELSE 0 END) AS format_pattern_unexpected,
  sum(CASE WHEN es_key_position IS NULL THEN 1 ELSE 0 END) AS format_missing_ES,
  sum(CASE WHEN se_key_position IS NULL THEN 1 ELSE 0 END) AS format_missing_SE,
  sum(CASE WHEN lp_key_position IS NULL THEN 1 ELSE 0 END) AS format_missing_LP,
  sum(CASE WHEN af_key_position IS NULL THEN 1 ELSE 0 END) AS format_missing_AF,
  sum(CASE WHEN beta_nonfinite THEN 1 ELSE 0 END) AS beta_nonfinite,
  sum(CASE WHEN se_nonpositive THEN 1 ELSE 0 END) AS se_nonpositive,
  sum(CASE WHEN lp IS NULL THEN 1 ELSE 0 END) AS lp_missing,
  sum(CASE WHEN lp IS NOT NULL AND NOT isfinite(lp) THEN 1 ELSE 0 END) AS lp_nonfinite,
  sum(CASE WHEN lp IS NOT NULL AND lp < 0 THEN 1 ELSE 0 END) AS lp_negative,
  sum(CASE WHEN eaf_out_of_range THEN 1 ELSE 0 END) AS eaf_out_of_range,
  sum(CASE WHEN non_autosomal_contig THEN 1 ELSE 0 END) AS non_autosomal_contig,
  sum(CASE WHEN lp IS NOT NULL AND isfinite(lp) AND lp > 7.301029995663981 THEN 1 ELSE 0 END) AS significant_lp_records,
  min(beta) FILTER (WHERE isfinite(beta)) AS beta_min,
  max(beta) FILTER (WHERE isfinite(beta)) AS beta_max,
  min(se) FILTER (WHERE isfinite(se)) AS se_min,
  max(se) FILTER (WHERE isfinite(se)) AS se_max,
  min(lp) FILTER (WHERE isfinite(lp)) AS lp_min,
  max(lp) FILTER (WHERE isfinite(lp)) AS lp_max,
  min(eaf) FILTER (WHERE isfinite(eaf)) AS eaf_min,
  max(eaf) FILTER (WHERE isfinite(eaf)) AS eaf_max
FROM base
"))
full_duplicates <- DBI::dbGetQuery(conn, paste0(common_cte, "
SELECT coalesce(sum(n), 0) AS duplicate_records, count(*) AS duplicate_groups
FROM (
  SELECT variant_id_source, count(*) AS n
  FROM base
  WHERE variant_id_source IS NOT NULL
  GROUP BY variant_id_source
  HAVING count(*) > 1
) d
"))
log_line("Full QC and duplicate audit completed.")
stage_end("global_qc_aggregation", global_qc_stage)

candidate_stage_timer <- stage_start("candidate_extraction")
DBI::dbExecute(conn, paste0("CREATE TEMP TABLE significant_stage AS ", common_cte, "
SELECT *,
  (es_key_position IS NULL OR se_key_position IS NULL OR lp_key_position IS NULL OR af_key_position IS NULL OR
   beta IS NULL OR se IS NULL OR eaf IS NULL) AS core_field_missing,
  rtrim(concat(
    CASE WHEN rsid_missing THEN 'rsid_missing;' ELSE '' END,
    CASE WHEN NOT position_valid THEN 'position_invalid;' ELSE '' END,
    CASE WHEN allele_missing THEN 'allele_missing;' ELSE '' END,
    CASE WHEN allele_invalid THEN 'allele_invalid;' ELSE '' END,
    CASE WHEN indel THEN 'indel;' ELSE '' END,
    CASE WHEN multi_allelic THEN 'multi_allelic;' ELSE '' END,
    CASE WHEN format_pattern_unexpected THEN 'format_pattern_unexpected;' ELSE '' END,
    CASE WHEN beta_nonfinite THEN 'beta_nonfinite;' ELSE '' END,
    CASE WHEN se_nonpositive THEN 'se_nonpositive;' ELSE '' END,
    CASE WHEN lp_invalid THEN 'lp_invalid;' ELSE '' END,
    CASE WHEN lp > 300 THEN 'pval_capped;' ELSE '' END,
    CASE WHEN eaf_out_of_range THEN 'eaf_out_of_range;' ELSE '' END,
    CASE WHEN non_autosomal_contig THEN 'non_autosomal_contig;' ELSE '' END
  ), ';') AS qc_flag_base
FROM base
WHERE lp IS NOT NULL AND isfinite(lp) AND lp > 7.301029995663981
"))
significant_count <- DBI::dbGetQuery(conn, "SELECT count(*) AS n FROM significant_stage")$n[[1]]
DBI::dbExecute(conn, "
CREATE TEMP TABLE candidate_stage AS
SELECT *
FROM significant_stage
WHERE NOT core_field_missing
")
candidate_count <- DBI::dbGetQuery(conn, "SELECT count(*) AS n FROM candidate_stage")$n[[1]]
quarantine_count <- DBI::dbGetQuery(conn, "SELECT count(*) AS n FROM significant_stage WHERE core_field_missing")$n[[1]]
log_line("Significant stage=", significant_count, "; main candidates=", candidate_count, "; quarantine=", quarantine_count)

DBI::dbExecute(conn, "
CREATE TEMP TABLE candidate_final AS
SELECT
  source_row_number,
  'Vuckovic2020_GCST90002384' AS source,
  'hemoglobin' AS trait,
  'GRCh37' AS build,
  chr_raw AS chr,
  pos,
  variant_id_source,
  allele_set_key,
  rsid_raw AS rsid,
  alt_raw AS effect_allele,
  ref_raw AS other_allele,
  beta,
  se,
  CASE WHEN lp <= 300 THEN power(10.0, -lp) ELSE 1e-300 END AS pval,
  lp AS log10p,
  lp AS lp_raw,
  eaf,
  CAST(NULL AS BIGINT) AS n_variant,
  CAST(408112 AS BIGINT) AS n_study,
  'standardized inverse-normal transformed hemoglobin; beta per ALT allele' AS effect_scale,
  format_pattern,
  power(beta / se, 2) AS F_stat,
  rtrim(concat(
    qc_flag_base,
    CASE WHEN count(*) OVER (PARTITION BY variant_id_source) > 1 THEN 'duplicate_variant_id_source;' ELSE '' END,
    CASE WHEN count(*) OVER (PARTITION BY allele_set_key) > 1 THEN 'duplicate_allele_set_key;' ELSE '' END
  ), ';') AS qc_flag,
  NOT (NOT position_valid OR allele_missing OR allele_invalid OR format_pattern_unexpected OR beta_nonfinite OR se_nonpositive OR lp_invalid OR eaf_out_of_range) AS qc_pass_core
FROM candidate_stage
")
quarantine_stage_timer <- stage_start("quarantine_extraction")
quarantine_path_sql <- sql_path(conn, paste0(quarantine_parquet, ".partial"))
DBI::dbExecute(conn, paste0("
COPY (
  SELECT source_row_number, chr_raw AS raw_chr, pos_raw AS raw_pos, rsid_raw AS raw_rsid,
         ref_raw AS raw_ref, alt_raw AS raw_alt, format_pattern, sample_value,
         es_key_position, se_key_position, lp_key_position, af_key_position,
         beta, se, lp, eaf, core_field_missing, qc_flag_base AS failure_reasons
  FROM significant_stage
  WHERE core_field_missing
) TO ", quarantine_path_sql, " (FORMAT PARQUET, COMPRESSION ZSTD)
"))
stage_end("quarantine_extraction", quarantine_stage_timer)
candidate_path_sql <- sql_path(conn, paste0(candidate_parquet, ".partial"))
candidate_tsv_sql <- sql_path(conn, paste0(candidate_tsv, ".partial"))
DBI::dbExecute(conn, paste0("COPY candidate_final TO ", candidate_path_sql, " (FORMAT PARQUET, COMPRESSION ZSTD)"))
DBI::dbExecute(conn, paste0("COPY candidate_final TO ", candidate_tsv_sql, " (FORMAT CSV, HEADER, DELIMITER '\\t', COMPRESSION GZIP)"))
stage_end("candidate_extraction", candidate_stage_timer)

candidate_duplicates <- DBI::dbGetQuery(conn, "
SELECT variant_id_source, count(*) AS records
FROM candidate_final
WHERE variant_id_source IS NOT NULL
GROUP BY variant_id_source
HAVING count(*) > 1
ORDER BY records DESC, variant_id_source
")
utils::write.csv(candidate_duplicates, paste0(duplicates_csv, ".partial"), row.names = FALSE)
candidates_by_chr <- DBI::dbGetQuery(conn, "
SELECT chr, count(*) AS records
FROM candidate_final
GROUP BY chr
ORDER BY CASE WHEN try_cast(chr AS INTEGER) IS NULL THEN 999 ELSE try_cast(chr AS INTEGER) END, chr
")
utils::write.csv(candidates_by_chr, paste0(by_chr_csv, ".partial"), row.names = FALSE)

candidate_qc_flags <- DBI::dbGetQuery(conn, "
SELECT flag AS qc_flag, count(*) AS records
FROM (
  SELECT unnest(string_split(CASE WHEN qc_flag = '' THEN 'qc_none' ELSE qc_flag END, ';')) AS flag
  FROM candidate_final
) x
GROUP BY flag
ORDER BY records DESC, qc_flag
")
candidate_f <- DBI::dbGetQuery(conn, "
SELECT min(F_stat) AS min, median(F_stat) AS median, avg(F_stat) AS mean, max(F_stat) AS max,
       sum(CASE WHEN F_stat < 10 THEN 1 ELSE 0 END) AS below_10
FROM candidate_final
")
candidate_dup_summary <- DBI::dbGetQuery(conn, "
SELECT coalesce(sum(n), 0) AS duplicate_records, count(*) AS duplicate_groups
FROM (
  SELECT variant_id_source, count(*) AS n
  FROM candidate_final
  WHERE variant_id_source IS NOT NULL
  GROUP BY variant_id_source
  HAVING count(*) > 1
) d
")
candidate_allele_dup_summary <- DBI::dbGetQuery(conn, "
SELECT coalesce(sum(n), 0) AS duplicate_records, count(*) AS duplicate_groups
FROM (
  SELECT allele_set_key, count(*) AS n
  FROM candidate_final
  WHERE allele_set_key IS NOT NULL
  GROUP BY allele_set_key
  HAVING count(*) > 1
) d
")

output_validation_stage <- stage_start("output_readback_validation")
parquet_sql <- paste0("read_parquet(", candidate_path_sql, ")")
tsv_sql <- paste0("read_csv_auto(", candidate_tsv_sql, ", delim = '\\t', header = true, compression = 'gzip', all_varchar = true)")
parquet_count <- DBI::dbGetQuery(conn, paste0("SELECT count(*) AS n FROM ", parquet_sql))$n[[1]]
tsv_count <- DBI::dbGetQuery(conn, paste0("SELECT count(*) AS n FROM ", tsv_sql))$n[[1]]
variant_set_difference <- DBI::dbGetQuery(conn, paste0("
SELECT count(*) AS n FROM (
  (SELECT DISTINCT variant_id_source FROM ", parquet_sql, "
   EXCEPT
   SELECT DISTINCT variant_id_source FROM ", tsv_sql, ")
  UNION ALL
  (SELECT DISTINCT variant_id_source FROM ", tsv_sql, "
   EXCEPT
   SELECT DISTINCT variant_id_source FROM ", parquet_sql, ")
) d
"))$n[[1]]
schema <- DBI::dbGetQuery(conn, paste0("DESCRIBE SELECT * FROM ", parquet_sql))
first10 <- DBI::dbGetQuery(conn, paste0("SELECT * FROM ", parquet_sql, " LIMIT 10"))
last10 <- DBI::dbGetQuery(conn, paste0("SELECT * FROM ", parquet_sql, " LIMIT 10 OFFSET ", max(0L, as.integer(candidate_count) - 10L)))
outside_threshold <- DBI::dbGetQuery(conn, paste0("SELECT count(*) AS n FROM ", parquet_sql, " WHERE log10p <= 7.301029995663981"))$n[[1]]
calculated_sha256_after <- toupper(digest::digest(input, algo = "sha256", file = TRUE))
sha256_after_matches <- identical(calculated_sha256_after, expected_sha256)
output_field_names <- schema$column_name
required_output_fields <- c("source_row_number", "source", "trait", "build", "chr", "pos", "variant_id_source", "allele_set_key", "rsid", "effect_allele", "other_allele", "beta", "se", "pval", "log10p", "lp_raw", "eaf", "n_variant", "n_study", "effect_scale", "format_pattern", "F_stat", "qc_flag", "qc_pass_core")
schema_matches <- identical(output_field_names, required_output_fields)
validation_complete <- identical(as.numeric(parquet_count), as.numeric(candidate_count)) &&
  identical(as.numeric(tsv_count), as.numeric(candidate_count)) &&
  identical(as.numeric(variant_set_difference), 0) &&
  identical(as.numeric(outside_threshold), 0) &&
  sha256_after_matches && schema_matches &&
  nrow(first10) == min(10L, as.integer(candidate_count)) &&
  nrow(last10) == min(10L, as.integer(candidate_count))
if (!validation_complete) {
  log_line("Output validation FAILED: parquet=", parquet_count, "; tsv=", tsv_count, "; candidate=", candidate_count,
           "; set_difference=", variant_set_difference, "; outside_threshold=", outside_threshold,
           "; sha_after_match=", sha256_after_matches, "; schema_matches=", schema_matches)
  stop("One or more candidate-output validations failed; run is not marked successful.")
}

end_time <- Sys.time()
qc <- list(
  run = list(
    start_time = format(start_time, "%Y-%m-%dT%H:%M:%S%z"),
    end_time = format(end_time, "%Y-%m-%dT%H:%M:%S%z"),
    duration_seconds = elapsed_seconds(start_time, end_time),
    input = normalizePath(input, winslash = "/", mustWork = TRUE),
    input_sha256_before = list(calculated = calculated_sha256_before, expected = expected_sha256, matches_expected = sha256_before_matches),
    input_sha256_after = list(calculated = calculated_sha256_after, expected = expected_sha256, matches_expected = sha256_after_matches),
    vcf_metadata_lines = header$line_number - 1L,
    data_records = record_count,
    expected_data_records = expected_records,
    record_count_matches = identical(as.numeric(record_count), as.numeric(expected_records)),
    r_version = R.version.string,
    packages = as.list(vapply(c("renv", required), function(x) as.character(utils::packageVersion(x)), character(1))),
    validation_complete = validation_complete
  ),
  format_patterns = df_to_records(format_df),
  full_qc = df_to_records(full_qc),
  full_variant_id_source_duplicates = df_to_records(full_duplicates),
  candidates = list(
    lp_threshold = 7.301029995663981,
    raw_significant_records = significant_count,
    main_records = candidate_count,
    quarantine_records = quarantine_count,
    by_chromosome = df_to_records(candidates_by_chr),
    qc_flag_counts = df_to_records(candidate_qc_flags),
    F_stat = df_to_records(candidate_f),
    variant_id_source_duplicates = df_to_records(candidate_dup_summary),
    allele_set_key_duplicates = df_to_records(candidate_allele_dup_summary)
  ),
  output_validation = list(
    parquet_records = parquet_count,
    tsv_records = tsv_count,
    row_count_matches = identical(as.numeric(parquet_count), as.numeric(candidate_count)) && identical(as.numeric(tsv_count), as.numeric(candidate_count)),
    variant_id_source_set_difference = variant_set_difference,
    variant_id_source_sets_match = identical(as.numeric(variant_set_difference), 0),
    schema_matches_design = schema_matches,
    schema = df_to_records(schema),
    all_candidates_above_lp_threshold = identical(as.numeric(outside_threshold), 0),
    records_at_or_below_threshold = outside_threshold,
    first10 = df_to_records(first10),
    last10 = df_to_records(last10)
  )
)
tryCatch(
  jsonlite::write_json(qc, paste0(qc_json, ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null"),
  error = function(e) {
    log_line("JSON write FAILED; object=qc; error=", conditionMessage(e))
    stop("QC JSON write failed; run is not marked successful.")
  }
)
stage_end("output_readback_validation", output_validation_stage)
log_line("Completed successfully. Full records=", record_count, "; candidates=", candidate_count, "; quarantine=", quarantine_count,
         "; duration_seconds=", elapsed_seconds(start_time, end_time))

final_map <- setNames(
  c(candidate_parquet, candidate_tsv, quarantine_parquet, format_csv, duplicates_csv, by_chr_csv, qc_json, log_final),
  c(paste0(candidate_parquet, ".partial"), paste0(candidate_tsv, ".partial"), paste0(quarantine_parquet, ".partial"),
    paste0(format_csv, ".partial"), paste0(duplicates_csv, ".partial"), paste0(by_chr_csv, ".partial"),
    paste0(qc_json, ".partial"), log_out)
)
for (partial_path in names(final_map)) {
  final_path <- unname(final_map[[partial_path]])
  if (!file.rename(partial_path, final_path)) {
    stop("Failed to promote validated partial output: ", partial_path, " -> ", final_path)
  }
}



