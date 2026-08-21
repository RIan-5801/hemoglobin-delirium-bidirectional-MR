#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

pkgs <- c("digest", "DBI", "duckdb", "jsonlite")
if (any(!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))) {
  stop("Required installed package missing; no automatic installation.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript.exe R/12_extract_vuckovic_hb_reverse_primary_outcomes_v1.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}

root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

stop_if <- function(x, m) if (isTRUE(x)) stop(m, call. = FALSE)
hash_file <- function(p) digest::digest(file = p, algo = "sha256", serialize = FALSE)
sql_string <- function(p, must_work = TRUE) {
  paste0("'", gsub("'", "''", normalizePath(p, winslash = "/", mustWork = must_work), fixed = TRUE), "'")
}
sql_ident <- function(con, x) as.character(DBI::dbQuoteIdentifier(con, x))
strict_num <- function(x, label) {
  raw <- trimws(as.character(x))
  miss <- is.na(raw) | raw %in% c("", "NA", ".")
  y <- suppressWarnings(as.numeric(ifelse(miss, NA_character_, raw)))
  stop_if(any(!miss & is.na(y)), paste("Numeric conversion failed:", label))
  y
}
comp <- function(x) unname(c(A = "T", T = "A", C = "G", G = "C")[x])
setkey <- function(a, b) paste(sort(c(a, b)), collapse = "/")
records <- function(x) {
  if (!is.data.frame(x)) return(list())
  lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
}

vcf <- file.path(root, "data_raw", "gwas", "ebi-a-GCST90002384.vcf.gz")
freeze_json_path <- file.path(root, "results", "qc", "finngen_r13_delirium_reverse_primary_instruments_v4_freeze_v3.json")
freeze_manifest_path <- file.path(root, "results", "qc", "finngen_r13_delirium_reverse_primary_instruments_v4_freeze_manifest_v3.csv")
vuckovic_qc_path <- file.path(root, "results", "qc", "vuckovic_hb_fullscan_v3_qc.json")
included_parquet <- file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_primary_apoe_included_clumped_v4.parquet")
excluded_parquet <- file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_primary_apoe_excluded_clumped_v4.parquet")
decision_path <- file.path(root, "docs", "decisions", "48_vuckovic_hb_reverse_primary_outcome_extraction_v1_v1.1.md")

out <- c(
  targets = file.path(root, "data_derived", "reverse_outcome_extraction", "finngen_r13_delirium_reverse_primary_targets_v1.tsv"),
  all_pq = file.path(root, "data_derived", "reverse_outcome_extraction", "vuckovic_hb_reverse_primary_all_matches_v1.parquet"),
  all_tsv = file.path(root, "data_derived", "reverse_outcome_extraction", "vuckovic_hb_reverse_primary_all_matches_v1.tsv"),
  unique_pq = file.path(root, "data_derived", "reverse_outcome_extraction", "vuckovic_hb_reverse_primary_unique_matches_v1.parquet"),
  unique_tsv = file.path(root, "data_derived", "reverse_outcome_extraction", "vuckovic_hb_reverse_primary_unique_matches_v1.tsv"),
  status_csv = file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_target_match_status_v1.csv"),
  qc_json = file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_outcome_extraction_v1.json"),
  log = file.path(root, "results", "logs", "vuckovic_hb_reverse_primary_outcome_extraction_v1.log")
)
partials <- paste0(out, ".partial")

if (!file.exists(vcf)) stop("Authenticated Vuckovic source missing: ", vcf, call. = FALSE)
if (any(!file.exists(c(freeze_json_path, freeze_manifest_path, vuckovic_qc_path, included_parquet, excluded_parquet, decision_path)))) {
  stop("Required authoritative context file missing.", call. = FALSE)
}
stop_if(any(file.exists(c(out, partials))), "A V1 final or partial extraction target exists; refusing to overwrite.")

dir.create(dirname(out[["targets"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["status_csv"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["log"]]), recursive = TRUE, showWarnings = FALSE)

log_created <- FALSE
log_line <- function(x) cat(sprintf("%s %s\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), x), file = out[["log"]], append = TRUE)
safe_log_line <- function(x) if (isTRUE(log_created)) try(log_line(x), silent = TRUE)
safe_disconnect <- function(con) {
  if (!is.null(con)) try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
}

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

main <- function() {
  log_created <<- TRUE
  start_time <- Sys.time()
  log_line("Protocol=docs/protocol/analysis_plan_v1.1.md; Decision=48; stage=Vuckovic Hb reverse-primary targeted outcome extraction V1.")

  expected_vuckovic_sha <- tolower("C68C98C7800C59D9D64CB88739E2E245F8AAD6E4E886455CBF7424661AFC3D41")
  freeze <- jsonlite::fromJSON(freeze_json_path, simplifyVector = FALSE)
  freeze_manifest_sha_observed <- hash_file(freeze_manifest_path)
  stop_if(!identical(freeze$freeze_status, "passed"), "Freeze gate failed: freeze_status is not passed.")
  stop_if(!identical(freeze$authoritative_reverse_primary_instrument_version, "v4"), "Freeze gate failed: authoritative instrument version is not v4.")
  stop_if(!identical(freeze$precision_audit_version, "v3"), "Freeze gate failed: precision audit is not v3.")
  stop_if(!isTRUE(freeze$approved_for_reverse_primary_outcome_extraction), "Freeze gate failed: extraction not approved.")
  stop_if(length(freeze$hard_check_failures) != 0L, "Freeze gate failed: hard_check_failures is not empty.")
  stop_if(!identical(tolower(freeze$manifest_sha256), freeze_manifest_sha_observed), "Freeze manifest SHA does not match Freeze JSON.")
  log_line(sprintf("Freeze gate passed; manifest_sha256=%s", freeze_manifest_sha_observed))

  vuckovic_qc <- jsonlite::fromJSON(vuckovic_qc_path, simplifyVector = FALSE)
  certified_input <- normalizePath(vuckovic_qc$run$input, winslash = "/", mustWork = TRUE)
  stop_if(!identical(certified_input, normalizePath(vcf, winslash = "/", mustWork = TRUE)), "Vuckovic source path differs from fullscan V3 QC authority.")
  stop_if(!isTRUE(vuckovic_qc$run$input_sha256_before$matches_expected) || !isTRUE(vuckovic_qc$run$input_sha256_after$matches_expected), "Vuckovic fullscan V3 source SHA gate failed.")
  stop_if(tolower(vuckovic_qc$run$input_sha256_before$expected) != expected_vuckovic_sha, "Vuckovic expected SHA differs from authenticated metadata.")

  vuckovic_sha_before <- hash_file(vcf)
  stop_if(!identical(vuckovic_sha_before, expected_vuckovic_sha), "Current Vuckovic SHA before scan differs from authenticated metadata.")
  log_line(sprintf("Vuckovic source SHA before=%s", vuckovic_sha_before))

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(safe_disconnect(con), add = TRUE)

  read_pq <- function(path) DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(path)))
  inc <- read_pq(included_parquet)
  exc <- read_pq(excluded_parquet)
  required_iv <- c("chromosome", "position", "rsid", "ref", "alt", "beta", "se", "pval", "eaf", "F_statistic", "finngen_log10p", "plink_log10p", "apoe_region", "instrument_membership")
  stop_if(!identical(names(inc), required_iv) || !identical(names(exc), required_iv), "Frozen V4 instrument Parquet schema mismatch.")
  stop_if(nrow(inc) != 1L || nrow(exc) != 1L, "Frozen V4 targeted extraction expects one included and one excluded target.")
  stop_if(!identical(inc$rsid, freeze$included_instrument_rsid) || !identical(exc$rsid, freeze$excluded_instrument_rsid), "Freeze JSON rsIDs differ from V4 Parquet targets.")
  shared_targets <- intersect(inc$rsid, exc$rsid)
  union_targets <- union(inc$rsid, exc$rsid)
  stop_if(length(shared_targets) != 0L || length(union_targets) != 2L, "Target count gate failed against Freeze V3.")

  mk_target <- function(x, membership) {
    data.frame(
      target_rsid = x$rsid,
      membership = membership,
      exposure_source = "FinnGen_R13_F5_DELIRIUM",
      exposure_trait = "delirium",
      exposure_build = "GRCh38",
      exposure_chr = as.character(x$chromosome),
      exposure_pos = as.character(x$position),
      exposure_other_allele_raw = toupper(as.character(x$ref)),
      exposure_effect_allele_raw = toupper(as.character(x$alt)),
      exposure_beta_raw = as.numeric(x$beta),
      exposure_se_raw = as.numeric(x$se),
      exposure_pval_raw = as.numeric(x$pval),
      exposure_eaf_raw = as.numeric(x$eaf),
      exposure_F_statistic = as.numeric(x$F_statistic),
      exposure_apoe_region = as.logical(x$apoe_region),
      stringsAsFactors = FALSE
    )
  }
  targets <- rbind(mk_target(inc, "apoe_included"), mk_target(exc, "apoe_excluded"))
  write.table(targets, paste0(out[["targets"]], ".partial"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")

  header <- find_vcf_header(vcf)
  required_vcf <- c("#CHROM", "POS", "ID", "REF", "ALT", "FORMAT")
  stop_if(!all(required_vcf %in% header$fields), "VCF #CHROM header lacks required columns.")
  sample_cols <- header$fields[(match("FORMAT", header$fields) + 1L):length(header$fields)]
  stop_if(length(sample_cols) != 1L, "Expected exactly one sample column after FORMAT in Vuckovic VCF.")
  sample_col <- sample_cols[[1L]]
  vcf_info <- grep("^##FORMAT=<ID=(ES|SE|LP|AF),", header$metadata, value = TRUE)
  stop_if(length(vcf_info) < 4L, "VCF metadata does not contain all required ES/SE/LP/AF FORMAT definitions.")
  log_line(sprintf("VCF header line=%s; sample_column=%s", header$line_number, sample_col))

  scan_sql <- paste0(
    "read_csv_auto(", sql_string(vcf),
    ", delim='\\t', header=true, skip=", header$line_number - 1L,
    ", compression='gzip', all_varchar=true, sample_size=10000, ignore_errors=false)"
  )
  probe <- DBI::dbGetQuery(con, paste0("SELECT * FROM ", scan_sql, " LIMIT 10"))
  stop_if(nrow(probe) != 10L || !all(required_vcf %in% names(probe)), "DuckDB Vuckovic gzip/BGZF probe failed schema gate.")
  stop_if(length(names(probe)[(match("FORMAT", names(probe)) + 1L):ncol(probe)]) != 1L, "DuckDB probe did not return exactly one sample column.")
  sample_col_sql <- sql_ident(con, sample_col)

  target_values <- paste(sprintf("('%s')", targets$target_rsid), collapse = ",")
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
  eof_observed <- TRUE
  observed_data_rows <- if (nrow(matches_raw)) unique(as.numeric(matches_raw$observed_data_rows)) else {
    DBI::dbGetQuery(con, paste0("SELECT count(*) AS n FROM ", scan_sql))$n[[1L]]
  }
  stop_if(length(observed_data_rows) != 1L, "Observed data-row audit was not unique.")
  stop_if(as.numeric(observed_data_rows) != as.numeric(vuckovic_qc$run$data_records), "Observed Vuckovic data rows differ from fullscan V3 QC.")
  scan_seconds <- as.numeric(difftime(scan_end, scan_start, units = "secs"))
  log_line(sprintf("Streaming scan completed; seconds=%.3f; eof_observed=%s; observed_data_rows=%s", scan_seconds, eof_observed, observed_data_rows))

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
    matches_raw$outcome_eaf_raw <- strict_num(get_key("AF"), "AF")
    matches_raw$es_key_position <- vapply(keys, function(x) match("ES", x), integer(1L))
    matches_raw$se_key_position <- vapply(keys, function(x) match("SE", x), integer(1L))
    matches_raw$lp_key_position <- vapply(keys, function(x) match("LP", x), integer(1L))
    matches_raw$af_key_position <- vapply(keys, function(x) match("AF", x), integer(1L))
  }

  status <- merge(targets, aggregate(vuckovic_scan_row_number ~ target_rsid, matches_raw, length), by = "target_rsid", all.x = TRUE)
  names(status)[names(status) == "vuckovic_scan_row_number"] <- "match_count"
  status$match_count[is.na(status$match_count)] <- 0L
  status$match_class <- ifelse(status$match_count == 0L, "missing", ifelse(status$match_count == 1L, "unique", "multiple"))
  stop_if(any(status$match_class == "multiple"), "Multiple Vuckovic matches detected; no automatic row selection permitted.")

  all_matches <- merge(status, matches_raw, by = "target_rsid", all.y = TRUE, sort = FALSE)
  if (nrow(all_matches)) {
    all_matches$outcome_source <- "Vuckovic2020_GCST90002384"
    all_matches$outcome_trait <- "hemoglobin"
    all_matches$outcome_build <- "GRCh37"
    all_matches$outcome_chr_raw <- all_matches$vuckovic_chr_raw
    all_matches$outcome_pos_raw <- all_matches$vuckovic_pos_raw
    all_matches$outcome_rsid_raw <- all_matches$vuckovic_rsid_raw
    all_matches$outcome_other_allele_raw <- toupper(all_matches$vuckovic_ref_raw)
    all_matches$outcome_effect_allele_raw <- toupper(all_matches$vuckovic_alt_raw)
    all_matches$outcome_pval_raw <- 10^(-all_matches$outcome_lp_raw)
    all_matches$outcome_n_study <- 408112L
    all_matches$outcome_effect_scale <- "standardized inverse-normal transformed haemoglobin; beta per ALT allele"
    all_matches$exposure_allele_set <- mapply(setkey, all_matches$exposure_effect_allele_raw, all_matches$exposure_other_allele_raw)
    all_matches$outcome_allele_set <- mapply(setkey, all_matches$outcome_effect_allele_raw, all_matches$outcome_other_allele_raw)
    ex_ok <- all_matches$exposure_effect_allele_raw %in% c("A", "C", "G", "T") & all_matches$exposure_other_allele_raw %in% c("A", "C", "G", "T")
    out_ok <- all_matches$outcome_effect_allele_raw %in% c("A", "C", "G", "T") & all_matches$outcome_other_allele_raw %in% c("A", "C", "G", "T")
    exact <- ex_ok & out_ok & all_matches$exposure_effect_allele_raw == all_matches$outcome_effect_allele_raw & all_matches$exposure_other_allele_raw == all_matches$outcome_other_allele_raw
    swapped <- ex_ok & out_ok & all_matches$exposure_effect_allele_raw == all_matches$outcome_other_allele_raw & all_matches$exposure_other_allele_raw == all_matches$outcome_effect_allele_raw
    strand_exact <- ex_ok & out_ok & all_matches$exposure_effect_allele_raw == comp(all_matches$outcome_effect_allele_raw) & all_matches$exposure_other_allele_raw == comp(all_matches$outcome_other_allele_raw)
    strand_swapped <- ex_ok & out_ok & all_matches$exposure_effect_allele_raw == comp(all_matches$outcome_other_allele_raw) & all_matches$exposure_other_allele_raw == comp(all_matches$outcome_effect_allele_raw)
    all_matches$palindromic_snp <- ex_ok & all_matches$exposure_allele_set %in% c("A/T", "C/G")
    all_matches$orientation_class <- ifelse(!ex_ok | !out_ok, "invalid",
      ifelse(exact, "exact_match",
      ifelse(swapped, "swapped_match",
      ifelse(strand_exact, "strand_exact_match",
      ifelse(strand_swapped, "strand_swapped_match", "incompatible")))))
  }

  unique_matches <- all_matches[all_matches$match_class == "unique", , drop = FALSE]
  write.table(status[, c("target_rsid", "membership", "match_count", "match_class", "exposure_effect_allele_raw", "exposure_other_allele_raw", "exposure_beta_raw", "exposure_se_raw", "exposure_pval_raw", "exposure_eaf_raw", "exposure_chr", "exposure_pos")], paste0(out[["status_csv"]], ".partial"), sep = ",", row.names = FALSE, quote = TRUE, na = "")

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

  staged_all_pq <- DBI::dbGetQuery(con, sprintf("SELECT target_rsid, membership FROM read_parquet(%s)", sql_string(paste0(out[["all_pq"]], ".partial"))))
  staged_unique_pq <- DBI::dbGetQuery(con, sprintf("SELECT target_rsid, membership FROM read_parquet(%s)", sql_string(paste0(out[["unique_pq"]], ".partial"))))
  staged_all_tsv <- read.delim(paste0(out[["all_tsv"]], ".partial"), check.names = FALSE)
  staged_unique_tsv <- read.delim(paste0(out[["unique_tsv"]], ".partial"), check.names = FALSE)
  all_consistent <- nrow(staged_all_pq) == nrow(staged_all_tsv) && setequal(staged_all_pq$target_rsid, staged_all_tsv$target_rsid)
  unique_consistent <- nrow(staged_unique_pq) == nrow(staged_unique_tsv) && setequal(staged_unique_pq$target_rsid, staged_unique_tsv$target_rsid)
  stop_if(!all_consistent || !unique_consistent, "Parquet/TSV consistency check failed.")

  vuckovic_sha_after <- hash_file(vcf)
  stop_if(!identical(vuckovic_sha_after, vuckovic_sha_before) || !identical(vuckovic_sha_after, expected_vuckovic_sha), "Vuckovic SHA changed during scan or differs from authenticated SHA.")

  hard_check_failures <- character()
  extraction_status <- if (length(hard_check_failures) == 0L && !any(status$match_class == "multiple") && all_consistent && unique_consistent) "passed" else "failed"
  approved_for_reverse_primary_harmonisation_preflight <- identical(extraction_status, "passed")
  included_available_for_harmonisation <- any(unique_matches$membership == "apoe_included")
  excluded_available_for_harmonisation <- any(unique_matches$membership == "apoe_excluded")

  qc <- list(
    decision_number = 48L,
    freeze_gate = list(status = freeze$freeze_status, manifest_sha256 = freeze_manifest_sha_observed),
    targets = list(included = inc$rsid, excluded = exc$rsid, included_target_count = nrow(inc), excluded_target_count = nrow(exc), shared_target_count = length(shared_targets), union_target_count = length(union_targets)),
    vuckovic_source = list(path = normalizePath(vcf, winslash = "/", mustWork = TRUE), sha256 = expected_vuckovic_sha, build = "GRCh37", effect_allele = "ALT", other_allele = "REF", beta = "ES", se = "SE", p = "LP = -log10(P)", eaf = "AF", n_study = 408112L, effect_scale = "standardized inverse-normal transformed haemoglobin"),
    vcf_header = list(line_number = header$line_number, sample_column = sample_col, format_definitions = as.list(vcf_info), ES_position = "FORMAT token position per row", SE_position = "FORMAT token position per row", LP_position = "FORMAT token position per row", AF_position = "FORMAT token position per row"),
    streaming_scan = list(duration_seconds = scan_seconds, eof_observed = eof_observed, observed_data_rows = as.numeric(observed_data_rows)),
    match_status = records(status),
    unique_matches = records(unique_matches),
    counts = list(unique = sum(status$match_class == "unique"), missing = sum(status$match_class == "missing"), multiple = sum(status$match_class == "multiple")),
    source_sha = list(before = vuckovic_sha_before, after = vuckovic_sha_after, certified = expected_vuckovic_sha),
    parquet_tsv_consistency = list(all_matches = all_consistent, unique_matches = unique_consistent),
    hard_check_failures = hard_check_failures,
    extraction_status = extraction_status,
    included_available_for_harmonisation = included_available_for_harmonisation,
    excluded_available_for_harmonisation = excluded_available_for_harmonisation,
    approved_for_reverse_primary_harmonisation_preflight = approved_for_reverse_primary_harmonisation_preflight,
    relaxed_threshold_status = "not_triggered",
    prohibited_actions = list(harmonisation = FALSE, beta_flip = FALSE, eaf_flip = FALSE, palindromic_exclusion = FALSE, reverse_mr = FALSE, steiger = FALSE, mr_presso = FALSE, proxy = FALSE, liftover = FALSE)
  )
  jsonlite::write_json(qc, paste0(out[["qc_json"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")

  finalizable <- c(out[["targets"]], out[["all_pq"]], out[["all_tsv"]], out[["unique_pq"]], out[["unique_tsv"]], out[["status_csv"]], out[["qc_json"]])
  for (fp in finalizable) {
    stop_if(file.exists(fp), paste("Output appeared during run:", fp))
    stop_if(!file.rename(paste0(fp, ".partial"), fp), paste("Atomic rename failed:", fp))
  }
  log_line(sprintf("SUCCESS: extraction_status=%s; unique=%s; missing=%s; multiple=%s; duration_seconds=%.3f", extraction_status, sum(status$match_class == "unique"), sum(status$match_class == "missing"), sum(status$match_class == "multiple"), as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
}

tryCatch(main(), error = function(e) {
  safe_log_line(paste0("TERMINATED_PRIMARY: ", conditionMessage(e)))
  quit(status = 1L)
})

