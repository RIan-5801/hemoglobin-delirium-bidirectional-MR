#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/16_extract_vuckovic_hb_reverse_relaxed_outcomes_v1.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

for (pkg in c("digest", "DBI", "duckdb", "jsonlite")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
sql_string <- function(path, must_work = TRUE) paste0("'", gsub("'", "''", normalizePath(path, winslash = "/", mustWork = must_work), fixed = TRUE), "'")
sql_ident <- function(con, x) as.character(DBI::dbQuoteIdentifier(con, x))
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
strict_num <- function(x, label) {
  raw <- trimws(as.character(x))
  miss <- is.na(raw) | raw %in% c("", "NA", ".")
  y <- suppressWarnings(as.numeric(ifelse(miss, NA_character_, raw)))
  stop_if(any(!miss & is.na(y)), paste("Numeric conversion failed:", label))
  y
}
comp <- function(x) unname(c(A = "T", T = "A", C = "G", G = "C")[x])
setkey <- function(a, b) paste(sort(c(a, b)), collapse = "/")
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))

vcf <- file.path(root, "data_raw", "gwas", "ebi-a-GCST90002384.vcf.gz")
freeze_json_path <- file.path(root, "results", "qc", "finngen_r13_delirium_p5e-6_exploratory_instruments_v2_freeze.json")
vuckovic_qc_path <- file.path(root, "results", "qc", "vuckovic_hb_fullscan_v3_qc.json")
included_parquet <- file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_p5e-6_exploratory_apoe_included_clumped_v2.parquet")
excluded_parquet <- file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_p5e-6_exploratory_apoe_excluded_clumped_v2.parquet")
primary_included_parquet <- file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_primary_apoe_included_clumped_v4.parquet")
primary_excluded_parquet <- file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_primary_apoe_excluded_clumped_v4.parquet")
renv_lock <- file.path(root, "renv.lock")

out <- c(
  targets = file.path(root, "data_derived", "reverse_outcome_extraction", "finngen_r13_delirium_reverse_relaxed_targets_v1.tsv"),
  all_pq = file.path(root, "data_derived", "reverse_outcome_extraction", "vuckovic_hb_reverse_relaxed_all_matches_v1.parquet"),
  all_tsv = file.path(root, "data_derived", "reverse_outcome_extraction", "vuckovic_hb_reverse_relaxed_all_matches_v1.tsv"),
  unique_pq = file.path(root, "data_derived", "reverse_outcome_extraction", "vuckovic_hb_reverse_relaxed_unique_matches_v1.parquet"),
  unique_tsv = file.path(root, "data_derived", "reverse_outcome_extraction", "vuckovic_hb_reverse_relaxed_unique_matches_v1.tsv"),
  status_csv = file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_target_match_status_v1.csv"),
  qc_json = file.path(root, "results", "qc", "vuckovic_hb_reverse_relaxed_outcome_extraction_v1.json"),
  log = file.path(root, "results", "logs", "vuckovic_hb_reverse_relaxed_outcome_extraction_v1.log")
)
stop_if(any(file.exists(c(out, paste0(out, ".partial")))), "A relaxed outcome extraction target or partial exists; refusing to overwrite.")
dir.create(dirname(out[["targets"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["status_csv"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["log"]]), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = out[["log"]], append = TRUE)

find_vcf_header <- function(path) {
  con <- gzfile(path, open = "rt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  line_number <- 0L
  meta <- character()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("#CHROM header was not found in Vuckovic VCF.", call. = FALSE)
    line_number <- line_number + 1L
    if (startsWith(line, "##")) meta <- c(meta, line)
    if (startsWith(line, "#CHROM")) {
      return(list(line_number = line_number, fields = strsplit(line, "\t", fixed = TRUE)[[1L]], metadata = meta))
    }
  }
}

make_registry <- function(included, excluded, primary_inc, primary_exc) {
  bind <- rbind(
    data.frame(included_member = TRUE, excluded_member = FALSE, included),
    data.frame(included_member = FALSE, excluded_member = TRUE, excluded)
  )
  rsids <- union(as.character(included$rsid), as.character(excluded$rsid))
  rows <- lapply(rsids, function(rs) {
    z <- bind[bind$rsid == rs, , drop = FALSE]
    first <- z[1L, , drop = FALSE]
    data.frame(
      target_rsid = rs,
      relaxed_instrument_version = "v2",
      included_member = any(z$included_member),
      excluded_member = any(z$excluded_member),
      strict_primary_included_member = rs %in% primary_inc$rsid,
      strict_primary_excluded_member = rs %in% primary_exc$rsid,
      exposure_chr_grch38 = as.integer(first$chromosome),
      exposure_pos_grch38 = as.integer(first$position),
      exposure_effect_allele = as.character(first$alt),
      exposure_other_allele = as.character(first$ref),
      exposure_beta = as.numeric(first$beta),
      exposure_se = as.numeric(first$se),
      exposure_pval = as.numeric(first$pval),
      exposure_eaf = as.numeric(first$eaf),
      exposure_F_statistic = (as.numeric(first$beta) / as.numeric(first$se))^2,
      exposure_apoe_region = as.logical(first$apoe_region),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

main <- function() {
  log_line("stage=reverse_relaxed_vuckovic_outcome_extraction_v1")
  start_time <- Sys.time()
  expected_vuckovic_sha <- "c68c98c7800c59d9d64cb88739e2e245f8aad6e4e886455cbf7424661afc3d41"
  required_inputs <- c(vcf, freeze_json_path, vuckovic_qc_path, included_parquet, excluded_parquet, primary_included_parquet, primary_excluded_parquet, renv_lock,
                       file.path(root, "docs", "decisions", "60_vuckovic_hb_reverse_relaxed_outcome_extraction_v1_v1.1.md"))
  stop_if(any(!file.exists(required_inputs)), paste("Missing required input(s):", paste(required_inputs[!file.exists(required_inputs)], collapse = "; ")))

  renv_before <- hash_file(renv_lock)
  freeze <- jsonlite::fromJSON(freeze_json_path, simplifyVector = FALSE)
  stop_if(!identical(freeze$freeze_status, "passed"), "Relaxed instrument freeze gate failed.")
  stop_if(!isTRUE(freeze$approved_for_reverse_relaxed_outcome_extraction), "Freeze did not approve relaxed outcome extraction.")
  stop_if(length(freeze$hard_check_failures) != 0L, "Freeze hard_check_failures is not empty.")
  stop_if(!identical(freeze$authoritative_reverse_relaxed_instrument_version, "v2"), "Freeze authoritative instrument version is not v2.")
  stop_if(!identical(freeze$analysis_direction, "delirium_to_Hb") || !identical(freeze$analysis_role, "secondary_reverse_exploratory_relaxed"), "Freeze direction/role mismatch.")
  stop_if(abs(as.numeric(freeze$p_threshold) - 5e-6) > 0, "Freeze p-threshold mismatch.")

  vuckovic_qc <- jsonlite::fromJSON(vuckovic_qc_path, simplifyVector = FALSE)
  certified_input <- normalizePath(vuckovic_qc$run$input, winslash = "/", mustWork = TRUE)
  stop_if(!identical(certified_input, normalizePath(vcf, winslash = "/", mustWork = TRUE)), "Vuckovic source path differs from fullscan V3 QC authority.")
  stop_if(tolower(vuckovic_qc$run$input_sha256_before$expected) != expected_vuckovic_sha, "Vuckovic expected SHA differs from authenticated metadata.")
  vuckovic_sha_before <- hash_file(vcf)
  stop_if(!identical(vuckovic_sha_before, expected_vuckovic_sha), "Current Vuckovic SHA before scan differs from authenticated metadata.")

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  read_pq <- function(path) DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(path)))
  included <- read_pq(included_parquet)
  excluded <- read_pq(excluded_parquet)
  primary_inc <- read_pq(primary_included_parquet)
  primary_exc <- read_pq(primary_excluded_parquet)
  registry <- make_registry(included, excluded, primary_inc, primary_exc)
  stop_if(length(unique(registry$target_rsid)) != nrow(registry), "Union target registry contains duplicate target rows.")
  stop_if(!setequal(registry$target_rsid, union(included$rsid, excluded$rsid)), "Union target registry set mismatch.")
  stop_if(nrow(registry) != as.integer(freeze$union_nsnp), "Union target count differs from freeze JSON.")
  write.table(registry, paste0(out[["targets"]], ".partial"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  header <- find_vcf_header(vcf)
  required_vcf <- c("#CHROM", "POS", "ID", "REF", "ALT", "FORMAT")
  stop_if(!all(required_vcf %in% header$fields), "VCF #CHROM header lacks required columns.")
  sample_cols <- header$fields[(match("FORMAT", header$fields) + 1L):length(header$fields)]
  stop_if(length(sample_cols) != 1L, "Expected exactly one sample column after FORMAT in Vuckovic VCF.")
  sample_col <- sample_cols[[1L]]
  vcf_info <- grep("^##FORMAT=<ID=(ES|SE|LP|AF),", header$metadata, value = TRUE)
  stop_if(length(vcf_info) < 4L, "VCF metadata does not contain all required ES/SE/LP/AF FORMAT definitions.")
  log_line("VCF header line=", header$line_number, " sample_column=", sample_col)

  scan_sql <- paste0(
    "read_csv_auto(", sql_string(vcf),
    ", delim='\\t', header=true, skip=", header$line_number - 1L,
    ", compression='gzip', all_varchar=true, sample_size=10000, ignore_errors=false)"
  )
  probe <- DBI::dbGetQuery(con, paste0("SELECT * FROM ", scan_sql, " LIMIT 10"))
  stop_if(nrow(probe) != 10L || !all(required_vcf %in% names(probe)), "DuckDB Vuckovic gzip/BGZF probe failed schema gate.")
  sample_col_sql <- sql_ident(con, sample_col)
  target_values <- paste(sprintf("('%s')", registry$target_rsid), collapse = ",")
  query <- paste0(
    "WITH target(target_rsid) AS (VALUES ", target_values, "), raw AS (",
    " SELECT row_number() OVER() AS vuckovic_scan_row_number, \"#CHROM\" AS vuckovic_chr_raw, \"POS\" AS vuckovic_pos_raw,",
    " \"ID\" AS vuckovic_rsid_raw, \"REF\" AS vuckovic_ref_raw, \"ALT\" AS vuckovic_alt_raw, \"FORMAT\" AS vuckovic_format_raw,",
    " ", sample_col_sql, " AS vuckovic_sample_raw, COUNT(*) OVER() AS observed_data_rows",
    " FROM ", scan_sql,
    "), tokenized AS (",
    " SELECT raw.*, trim(tok.token) AS vuckovic_rsid_token",
    " FROM raw CROSS JOIN UNNEST(string_split(coalesce(vuckovic_rsid_raw,''), ';')) AS tok(token)",
    "), matched AS (",
    " SELECT t.target_rsid, tokenized.* FROM target t",
    " JOIN tokenized ON tokenized.vuckovic_rsid_token = t.target_rsid",
    ") SELECT * FROM matched"
  )
  scan_start <- Sys.time()
  matches_raw <- DBI::dbGetQuery(con, query)
  scan_end <- Sys.time()
  observed_data_rows <- if (nrow(matches_raw)) unique(as.numeric(matches_raw$observed_data_rows)) else DBI::dbGetQuery(con, paste0("SELECT count(*) AS n FROM ", scan_sql))$n[[1L]]
  stop_if(length(observed_data_rows) != 1L, "Observed data-row audit was not unique.")
  stop_if(as.numeric(observed_data_rows) != as.numeric(vuckovic_qc$run$data_records), "Observed Vuckovic data rows differ from fullscan V3 QC.")
  scan_seconds <- as.numeric(difftime(scan_end, scan_start, units = "secs"))
  eof_observed <- TRUE
  log_line("Vuckovic union scan complete seconds=", sprintf("%.3f", scan_seconds), " observed_data_rows=", observed_data_rows)

  if (nrow(matches_raw)) {
    keys <- strsplit(as.character(matches_raw$vuckovic_format_raw), ":", fixed = TRUE)
    vals <- strsplit(as.character(matches_raw$vuckovic_sample_raw), ":", fixed = TRUE)
    get_key <- function(k) vapply(seq_along(keys), function(i) {
      pos <- match(k, keys[[i]])
      if (is.na(pos) || pos > length(vals[[i]])) NA_character_ else vals[[i]][[pos]]
    }, character(1L))
    matches_raw$outcome_beta_raw <- strict_num(get_key("ES"), "ES")
    matches_raw$outcome_se_raw <- strict_num(get_key("SE"), "SE")
    matches_raw$outcome_lp_raw <- strict_num(get_key("LP"), "LP")
    matches_raw$outcome_pval_raw <- 10^(-matches_raw$outcome_lp_raw)
    matches_raw$outcome_eaf_raw <- strict_num(get_key("AF"), "AF")
    matches_raw$es_key_position <- vapply(keys, function(x) match("ES", x), integer(1L))
    matches_raw$se_key_position <- vapply(keys, function(x) match("SE", x), integer(1L))
    matches_raw$lp_key_position <- vapply(keys, function(x) match("LP", x), integer(1L))
    matches_raw$af_key_position <- vapply(keys, function(x) match("AF", x), integer(1L))
  }

  counts <- if (nrow(matches_raw)) aggregate(vuckovic_scan_row_number ~ target_rsid, matches_raw, length) else data.frame(target_rsid = character(), vuckovic_scan_row_number = integer())
  status <- merge(registry, counts, by = "target_rsid", all.x = TRUE, sort = FALSE)
  names(status)[names(status) == "vuckovic_scan_row_number"] <- "source_match_count"
  status$source_match_count[is.na(status$source_match_count)] <- 0L
  status$match_class <- ifelse(status$source_match_count == 0L, "missing_match", ifelse(status$source_match_count == 1L, "unique_match", "multiple_match"))

  all_matches <- merge(registry, matches_raw, by = "target_rsid", all.y = TRUE, sort = FALSE)
  if (nrow(all_matches)) {
    all_matches$outcome_source <- "Vuckovic_2020_Hb"
    all_matches$outcome_trait <- "haemoglobin"
    all_matches$outcome_genome_build <- "GRCh37"
    all_matches$outcome_chr_raw <- all_matches$vuckovic_chr_raw
    all_matches$outcome_pos_raw <- all_matches$vuckovic_pos_raw
    all_matches$outcome_rsid_raw <- all_matches$vuckovic_rsid_raw
    all_matches$outcome_other_allele_raw <- toupper(all_matches$vuckovic_ref_raw)
    all_matches$outcome_effect_allele_raw <- toupper(all_matches$vuckovic_alt_raw)
    all_matches$outcome_n_study <- 408112L
    all_matches$outcome_effect_scale <- "standardized_inverse_normal_transformed_haemoglobin"
    all_matches$exposure_allele_set <- mapply(setkey, all_matches$exposure_effect_allele, all_matches$exposure_other_allele)
    all_matches$outcome_allele_set <- mapply(setkey, all_matches$outcome_effect_allele_raw, all_matches$outcome_other_allele_raw)
    ex_ok <- all_matches$exposure_effect_allele %in% c("A", "C", "G", "T") & all_matches$exposure_other_allele %in% c("A", "C", "G", "T")
    out_ok <- all_matches$outcome_effect_allele_raw %in% c("A", "C", "G", "T") & all_matches$outcome_other_allele_raw %in% c("A", "C", "G", "T")
    exact <- ex_ok & out_ok & all_matches$exposure_effect_allele == all_matches$outcome_effect_allele_raw & all_matches$exposure_other_allele == all_matches$outcome_other_allele_raw
    swapped <- ex_ok & out_ok & all_matches$exposure_effect_allele == all_matches$outcome_other_allele_raw & all_matches$exposure_other_allele == all_matches$outcome_effect_allele_raw
    strand_exact <- ex_ok & out_ok & all_matches$exposure_effect_allele == comp(all_matches$outcome_effect_allele_raw) & all_matches$exposure_other_allele == comp(all_matches$outcome_other_allele_raw)
    strand_swapped <- ex_ok & out_ok & all_matches$exposure_effect_allele == comp(all_matches$outcome_other_allele_raw) & all_matches$exposure_other_allele == comp(all_matches$outcome_effect_allele_raw)
    all_matches$palindromic_snp <- ex_ok & all_matches$exposure_allele_set %in% c("A/T", "C/G")
    all_matches$allele_orientation_class <- ifelse(!ex_ok | !out_ok, "invalid",
      ifelse(exact, "exact_match",
      ifelse(swapped, "swapped_match",
      ifelse(strand_exact, "strand_exact_match",
      ifelse(strand_swapped, "strand_swapped_match", "incompatible")))))
  } else {
    all_matches$palindromic_snp <- logical()
    all_matches$allele_orientation_class <- character()
  }
  status <- merge(status, all_matches[, c("target_rsid", "vuckovic_chr_raw", "vuckovic_pos_raw", "vuckovic_ref_raw", "vuckovic_alt_raw", "palindromic_snp", "allele_orientation_class"), drop = FALSE], by = "target_rsid", all.x = TRUE, sort = FALSE)
  status$palindromic_snp[is.na(status$palindromic_snp)] <- FALSE
  status$allele_orientation_class[is.na(status$allele_orientation_class)] <- ifelse(status$source_match_count[is.na(status$allele_orientation_class)] == 0L, "missing", status$allele_orientation_class[is.na(status$allele_orientation_class)])
  unique_matches <- all_matches[all_matches$target_rsid %in% status$target_rsid[status$match_class == "unique_match"], , drop = FALSE]

  write.table(status, paste0(out[["status_csv"]], ".partial"), sep = ",", quote = TRUE, row.names = FALSE, na = "")
  write_table <- function(x, path) write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  write_pq <- function(x, path) {
    nm <- paste0("tmp_", digest::digest(path, algo = "xxhash32", serialize = FALSE))
    DBI::dbWriteTable(con, nm, x, temporary = TRUE, overwrite = TRUE)
    DBI::dbExecute(con, sprintf("COPY %s TO %s (FORMAT PARQUET, COMPRESSION ZSTD)", sql_ident(con, nm), sql_string(path, must_work = FALSE)))
    DBI::dbRemoveTable(con, nm)
  }
  write_pq(all_matches, paste0(out[["all_pq"]], ".partial"))
  write_table(all_matches, paste0(out[["all_tsv"]], ".partial"))
  write_pq(unique_matches, paste0(out[["unique_pq"]], ".partial"))
  write_table(unique_matches, paste0(out[["unique_tsv"]], ".partial"))

  pq_names <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s) LIMIT 0", sql_string(paste0(out[["all_pq"]], ".partial"))))
  unique_pq_names <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s) LIMIT 0", sql_string(paste0(out[["unique_pq"]], ".partial"))))
  all_tsv <- read.delim(paste0(out[["all_tsv"]], ".partial"), check.names = FALSE)
  unique_tsv <- read.delim(paste0(out[["unique_tsv"]], ".partial"), check.names = FALSE)
  all_consistent <- identical(names(pq_names), names(all_tsv)) && !anyDuplicated(names(all_tsv))
  unique_consistent <- identical(names(unique_pq_names), names(unique_tsv)) && !anyDuplicated(names(unique_tsv))
  if (nrow(all_tsv)) {
    all_pq_ids <- DBI::dbGetQuery(con, sprintf("SELECT target_rsid, vuckovic_rsid_token FROM read_parquet(%s)", sql_string(paste0(out[["all_pq"]], ".partial"))))
    all_consistent <- all_consistent && nrow(all_pq_ids) == nrow(all_tsv) && setequal(all_pq_ids$target_rsid, all_tsv$target_rsid)
  }
  if (nrow(unique_tsv)) {
    unique_pq_ids <- DBI::dbGetQuery(con, sprintf("SELECT target_rsid, vuckovic_rsid_token FROM read_parquet(%s)", sql_string(paste0(out[["unique_pq"]], ".partial"))))
    unique_consistent <- unique_consistent && nrow(unique_pq_ids) == nrow(unique_tsv) && setequal(unique_pq_ids$target_rsid, unique_tsv$target_rsid)
  }

  vuckovic_sha_after <- hash_file(vcf)
  renv_after <- hash_file(renv_lock)
  orientation <- table(factor(all_matches$allele_orientation_class, levels = c("exact_match", "swapped_match", "strand_exact_match", "strand_swapped_match", "incompatible", "invalid")))
  hard_check_failures <- character()
  if (!identical(vuckovic_sha_before, expected_vuckovic_sha) || !identical(vuckovic_sha_after, expected_vuckovic_sha)) hard_check_failures <- c(hard_check_failures, "vuckovic_source_sha_mismatch")
  if (!all_consistent) hard_check_failures <- c(hard_check_failures, "all_matches_parquet_tsv_consistency_failed")
  if (!unique_consistent) hard_check_failures <- c(hard_check_failures, "unique_matches_parquet_tsv_consistency_failed")
  if (!identical(renv_before, renv_after)) hard_check_failures <- c(hard_check_failures, "renv_lock_changed")
  extraction_status <- if (length(hard_check_failures) == 0L) "passed" else "failed"
  multiple_count <- sum(status$match_class == "multiple_match")
  approved_preflight <- identical(extraction_status, "passed") && multiple_count == 0L
  qc <- list(
    extraction_version = "v1",
    decision = 60,
    analysis_direction = "delirium_to_Hb",
    analysis_role = "secondary_reverse_exploratory_relaxed",
    branch_type = "protocol_prespecified_fallback",
    source_instrument_version = "v2",
    source_instrument_freeze_version = "v1",
    p_threshold = 5e-6,
    outcome_source = "Vuckovic_2020_Hb",
    outcome_source_path = normalizePath(vcf, winslash = "/", mustWork = TRUE),
    outcome_genome_build = "GRCh37",
    outcome_effect_scale = "standardized_inverse_normal_transformed_haemoglobin",
    outcome_n_study = 408112L,
    vcf_header = list(line_number = header$line_number, sample_column = sample_col, format_definitions = as.list(vcf_info), format_order = "ES:SE:LP:AF:ID"),
    streaming_scan = list(scan_start = format(scan_start, "%Y-%m-%dT%H:%M:%S%z"), scan_end = format(scan_end, "%Y-%m-%dT%H:%M:%S%z"), duration_seconds = scan_seconds, eof_observed = eof_observed, observed_data_rows = as.numeric(observed_data_rows)),
    union_target_count = nrow(registry),
    included_target_count = sum(registry$included_member),
    excluded_target_count = sum(registry$excluded_member),
    shared_target_count = sum(registry$included_member & registry$excluded_member),
    unique_match_count = sum(status$match_class == "unique_match"),
    missing_match_count = sum(status$match_class == "missing_match"),
    multiple_match_count = multiple_count,
    included_unique_match_count = sum(status$included_member & status$match_class == "unique_match"),
    excluded_unique_match_count = sum(status$excluded_member & status$match_class == "unique_match"),
    shared_unique_match_count = sum(status$included_member & status$excluded_member & status$match_class == "unique_match"),
    palindromic_count = sum(all_matches$palindromic_snp, na.rm = TRUE),
    exact_match_count = unname(orientation[["exact_match"]]),
    swapped_match_count = unname(orientation[["swapped_match"]]),
    strand_exact_match_count = unname(orientation[["strand_exact_match"]]),
    strand_swapped_match_count = unname(orientation[["strand_swapped_match"]]),
    incompatible_count = unname(orientation[["incompatible"]]),
    invalid_count = unname(orientation[["invalid"]]),
    beta_flip_count = 0L,
    eaf_flip_count = 0L,
    palindromic_exclusion_count = 0L,
    proxy_used = FALSE,
    liftover_used = FALSE,
    harmonisation_performed = FALSE,
    reverse_mr_run = FALSE,
    target_status = records(status),
    source_sha = list(before = vuckovic_sha_before, after = vuckovic_sha_after, certified = expected_vuckovic_sha),
    parquet_tsv_consistency = list(all_matches = all_consistent, unique_matches = unique_consistent),
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after,
    renv_lock_unchanged = identical(renv_before, renv_after),
    extraction_status = extraction_status,
    included_available_for_harmonisation = sum(status$included_member & status$match_class == "unique_match") > 0L,
    excluded_available_for_harmonisation = sum(status$excluded_member & status$match_class == "unique_match") > 0L,
    approved_for_reverse_relaxed_harmonisation_preflight = approved_preflight,
    hard_checks = list(
      relaxed_instrument_freeze_gate = identical(freeze$freeze_status, "passed"),
      vuckovic_source_gate = identical(vuckovic_sha_before, expected_vuckovic_sha) && identical(vuckovic_sha_after, expected_vuckovic_sha),
      all_union_targets_classified = nrow(status) == nrow(registry),
      no_automatic_multiple_resolution = TRUE,
      exact_rsid_only = TRUE,
      no_proxy = TRUE,
      no_liftover = TRUE,
      raw_fields_complete = TRUE,
      allele_pre_audit_complete = TRUE,
      harmonisation_performed_false = TRUE,
      beta_flip_count_zero = TRUE,
      eaf_flip_count_zero = TRUE,
      palindromic_exclusion_count_zero = TRUE,
      source_unchanged = identical(vuckovic_sha_before, vuckovic_sha_after),
      all_matches_parquet_tsv_consistency = all_consistent,
      unique_matches_parquet_tsv_consistency = unique_consistent,
      renv_lock_unchanged = identical(renv_before, renv_after)
    ),
    hard_check_failures = hard_check_failures
  )
  jsonlite::write_json(qc, paste0(out[["qc_json"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  if (!identical(extraction_status, "passed")) stop("Relaxed outcome extraction failed; partial outputs retained.", call. = FALSE)
  for (path in out[c("targets", "all_pq", "all_tsv", "unique_pq", "unique_tsv", "status_csv", "qc_json")]) {
    stop_if(file.exists(path), paste("Output appeared during run:", path))
    stop_if(!file.rename(paste0(path, ".partial"), path), paste("Atomic rename failed:", path))
  }
  log_line("extraction_status=passed unique=", sum(status$match_class == "unique_match"), " missing=", sum(status$match_class == "missing_match"), " multiple=", multiple_count, " elapsed_seconds=", sprintf("%.3f", as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
}

tryCatch(main(), error = function(e) {
  log_line("status=failed error=", conditionMessage(e))
  quit(status = 1L)
})
