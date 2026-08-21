#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/15_reverse_relaxed_instruments_finngen_r13_delirium_p5e-6_v2.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

for (pkg in c("DBI", "duckdb", "jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
hash_file <- function(path) digest::digest(file = path, algo = "sha256")
sql_string <- function(path, must_work = TRUE) {
  paste0("'", gsub("'", "''", normalizePath(path, winslash = "/", mustWork = must_work), fixed = TRUE), "'")
}
sql_ident <- function(con, x) as.character(DBI::dbQuoteIdentifier(con, x))
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

p_threshold <- 5e-6
p_threshold_log10 <- -log10(p_threshold)
finngen_sha <- "85637f0f3358807964d4f8a3e500293168a706f1c08c65f3fc5512b65df40ed8"
plink_zip_sha <- "dafbf939c4f106e11660ab5b7c545e749e515cd0bc6622b87d83ee0cb434b5ea"
plink_sha <- "247491bfca7512e070dc99d6565e9fc56f3a52ad5afc01286016271d34c4992f"

out <- c(
  candidates = file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_p5e-6_exploratory_candidates_v2.parquet"),
  eligible = file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_p5e-6_exploratory_eligible_candidates_v2.parquet"),
  included_pq = file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_p5e-6_exploratory_apoe_included_clumped_v2.parquet"),
  included_tsv = file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_p5e-6_exploratory_apoe_included_clumped_v2.tsv"),
  excluded_pq = file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_p5e-6_exploratory_apoe_excluded_clumped_v2.parquet"),
  excluded_tsv = file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_p5e-6_exploratory_apoe_excluded_clumped_v2.tsv"),
  counts = file.path(root, "results", "qc", "finngen_r13_delirium_p5e-6_exploratory_instrument_selection_counts_v2.csv"),
  not_ref = file.path(root, "results", "qc", "finngen_r13_delirium_p5e-6_exploratory_not_in_reference_v2.tsv"),
  schema = file.path(root, "results", "qc", "finngen_r13_delirium_p5e-6_exploratory_schema_audit_v2.csv"),
  overlap = file.path(root, "results", "qc", "finngen_r13_delirium_p5e-6_exploratory_primary_overlap_audit_v2.csv"),
  dup = file.path(root, "results", "qc", "finngen_r13_delirium_p5e-6_exploratory_duplicate_rsid_audit_v2.tsv"),
  json = file.path(root, "results", "qc", "finngen_r13_delirium_p5e-6_exploratory_instrument_selection_v2.json"),
  log = file.path(root, "results", "logs", "finngen_r13_delirium_p5e-6_exploratory_instrument_selection_v2.log")
)
stop_if(any(file.exists(c(out, paste0(out, ".partial")))), "A relaxed V2 final or partial target exists; refusing to overwrite.")
for (d in unique(dirname(out))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = out[["log"]], append = TRUE)

write_pq <- function(con, x, path) {
  nm <- paste0("tmp_", digest::digest(path, algo = "xxhash32", serialize = FALSE))
  DBI::dbWriteTable(con, nm, x, temporary = TRUE, overwrite = TRUE)
  DBI::dbExecute(con, sprintf("COPY %s TO %s (FORMAT PARQUET, COMPRESSION ZSTD)", sql_ident(con, nm), sql_string(path, must_work = FALSE)))
  DBI::dbRemoveTable(con, nm)
}
read_pq <- function(con, path) DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_string(path)))
num <- function(x) suppressWarnings(as.numeric(x))

final_schema <- function(x, membership, plink_log10p) {
  data.frame(
    chromosome = as.integer(x$chromosome),
    position = as.integer(x$position),
    rsid = as.character(x$rsid),
    ref = as.character(x$ref),
    alt = as.character(x$alt),
    beta = as.numeric(x$beta),
    se = as.numeric(x$se),
    pval = as.numeric(x$pval),
    eaf = as.numeric(x$eaf),
    F_statistic = as.numeric(x$F_statistic),
    finngen_log10p = as.numeric(x$finngen_log10p),
    plink_log10p = as.numeric(plink_log10p),
    apoe_region = as.logical(x$apoe_region),
    instrument_membership = membership,
    stringsAsFactors = FALSE
  )
}
schema_counts <- function(x, analysis_set) {
  data.frame(
    analysis_set = analysis_set,
    column_count = length(names(x)),
    duplicate_column_count = sum(duplicated(names(x))),
    raw_LOG10P_count = sum(names(x) == "LOG10P"),
    plink_log10p_count = sum(names(x) == "plink_log10p"),
    finngen_log10p_count = sum(names(x) == "finngen_log10p"),
    schema_unique = !anyDuplicated(names(x)),
    duckdb_registration_ok = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
check_pq_tsv <- function(con, pq_path, tsv_path) {
  pq <- read_pq(con, pq_path)
  tsv <- read.delim(tsv_path, check.names = FALSE)
  identical(names(pq), names(tsv)) && nrow(pq) == nrow(tsv) && !anyDuplicated(names(pq)) && setequal(as.character(pq$rsid), as.character(tsv$rsid))
}

restore_plink <- function() {
  zip_path <- file.path(root, "tools", "plink2", "plink2_win64_20260504.zip")
  final_path <- file.path(root, "tools", "plink2", "plink2.exe")
  tmp_dir <- file.path(root, "tools", "plink2", ".restore_v2_tmp")
  tmp_exe <- file.path(tmp_dir, "plink2.exe")
  hard <- list(
    plink_zip_certification_gate = FALSE,
    plink_target_absent_before_restore = FALSE,
    plink_temp_extraction_completed = FALSE,
    plink_temp_sha_matches_expected = FALSE,
    plink_atomic_publish_completed = FALSE,
    plink_final_sha_matches_expected = FALSE,
    plink_version_command_passed = FALSE,
    plink_restoration_no_overwrite = FALSE
  )
  stop_if(!file.exists(zip_path), "PLINK2 zip missing.")
  zip_sha <- hash_file(zip_path)
  hard$plink_zip_certification_gate <- identical(zip_sha, plink_zip_sha)
  stop_if(!hard$plink_zip_certification_gate, "PLINK2 zip SHA mismatch.")
  binary_target_preexisted <- file.exists(final_path)
  hard$plink_target_absent_before_restore <- !binary_target_preexisted
  hard$plink_restoration_no_overwrite <- !binary_target_preexisted
  stop_if(binary_target_preexisted, "Target PLINK2 binary already exists; refusing to overwrite.")
  stop_if(file.exists(tmp_dir), "PLINK2 restoration temp directory already exists; refusing to reuse.")
  dir.create(tmp_dir, recursive = FALSE, showWarnings = FALSE)
  utils::unzip(zip_path, files = "plink2.exe", exdir = tmp_dir, overwrite = FALSE)
  hard$plink_temp_extraction_completed <- file.exists(tmp_exe)
  stop_if(!hard$plink_temp_extraction_completed, "Temporary PLINK2 extraction failed.")
  temp_sha <- hash_file(tmp_exe)
  hard$plink_temp_sha_matches_expected <- identical(temp_sha, plink_sha)
  stop_if(!hard$plink_temp_sha_matches_expected, "Temporary PLINK2 SHA mismatch.")
  hard$plink_atomic_publish_completed <- file.rename(tmp_exe, final_path)
  stop_if(!hard$plink_atomic_publish_completed, "Atomic publish of PLINK2 binary failed.")
  unlink(tmp_dir, recursive = TRUE, force = TRUE)
  final_sha <- hash_file(final_path)
  hard$plink_final_sha_matches_expected <- identical(final_sha, plink_sha)
  stop_if(!hard$plink_final_sha_matches_expected, "Final PLINK2 SHA mismatch.")
  version_output <- system2(final_path, args = "--version", stdout = TRUE, stderr = TRUE)
  version_status <- attr(version_output, "status")
  hard$plink_version_command_passed <- is.null(version_status) || identical(as.integer(version_status), 0L)
  stop_if(!hard$plink_version_command_passed, "PLINK2 version command failed.")
  list(
    zip_path = zip_path,
    zip_sha256 = zip_sha,
    archive_member = "plink2.exe",
    temporary_binary_sha256 = temp_sha,
    final_binary_path = final_path,
    final_binary_sha256 = final_sha,
    plink_version_output = paste(version_output, collapse = " | "),
    binary_target_preexisted = binary_target_preexisted,
    overwrite_performed = FALSE,
    scientific_parameters_changed = FALSE,
    plink_binary_restoration_status = "passed",
    hard_checks = hard
  )
}

main <- function() {
  log_line("stage=reverse_exploratory_relaxed_threshold_instrument_selection_v2")
  renv_lock <- file.path(root, "renv.lock")
  renv_before <- hash_file(renv_lock)
  restoration <- restore_plink()
  log_line("plink_binary_restoration_status=passed version=", restoration$plink_version_output)

  fg <- file.path(root, "data_raw", "gwas", "finngen_R13_F5_DELIRIUM.gz")
  plink <- restoration$final_binary_path
  eur <- file.path(root, "resources", "ld", "1kg_v3", "EUR")
  v1_required <- c(
    file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_p5e-6_exploratory_candidates_v1.parquet.partial"),
    file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_p5e-6_exploratory_eligible_candidates_v1.parquet.partial"),
    file.path(root, "results", "qc", "finngen_r13_delirium_p5e-6_exploratory_not_in_reference_v1.tsv.partial"),
    file.path(root, "results", "logs", "finngen_r13_delirium_p5e-6_exploratory_instrument_selection_v1.log")
  )
  inputs <- c(
    fg, plink, paste0(eur, c(".bed", ".bim", ".fam")), renv_lock,
    file.path(root, "results", "qc", "finngen_R13_F5_DELIRIUM_input_certification_v1.json"),
    file.path(root, "results", "qc", "finngen_grch38_apoe_region_certification_v1.json"),
    file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_formal_harmonisation_v4_freeze.json"),
    file.path(root, "results", "qc", "finngen_r13_delirium_reverse_primary_instrument_selection_v4.json"),
    file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_primary_apoe_included_clumped_v4.parquet"),
    file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_primary_apoe_excluded_clumped_v4.parquet"),
    file.path(root, "docs", "decisions", "58_finngen_r13_delirium_reverse_exploratory_relaxed_instrument_selection_v2_technical_recovery_v1.1.md"),
    v1_required
  )
  stop_if(any(!file.exists(inputs)), paste("Missing input(s):", paste(inputs[!file.exists(inputs)], collapse = "; ")))
  freeze <- jsonlite::fromJSON(file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_formal_harmonisation_v4_freeze.json"), simplifyVector = FALSE)
  cert <- jsonlite::fromJSON(file.path(root, "results", "qc", "finngen_R13_F5_DELIRIUM_input_certification_v1.json"), simplifyVector = FALSE)
  apoe <- jsonlite::fromJSON(file.path(root, "results", "qc", "finngen_grch38_apoe_region_certification_v1.json"), simplifyVector = FALSE)
  primary <- jsonlite::fromJSON(file.path(root, "results", "qc", "finngen_r13_delirium_reverse_primary_instrument_selection_v4.json"), simplifyVector = FALSE)
  stop_if(!identical(freeze$freeze_status, "passed") || !isTRUE(freeze$approved_for_reverse_relaxed_threshold_instrument_selection), "Decision 56 freeze gate failed.")
  stop_if(!identical(freeze$relaxed_threshold_status, "triggered_not_started") || abs(as.numeric(freeze$relaxed_threshold) - p_threshold) > 0, "Relaxed threshold trigger mismatch.")
  stop_if(!identical(cert$certification_status, "passed") || as.character(cert$certified_input_sha256) != finngen_sha, "FinnGen certification gate failed.")
  stop_if(!identical(apoe$certification_status, "passed") || !isTRUE(apoe$approved_for_reverse_apoe_exclusion), "APOE certification gate failed.")
  stop_if(!identical(primary$instrument_selection_status, "passed") || abs(as.numeric(primary$p_threshold) - 5e-8) > 0, "Strict primary instrument-selection gate failed.")
  finngen_before <- hash_file(fg)
  stop_if(finngen_before != finngen_sha, "FinnGen SHA before mismatch.")
  stop_if(hash_file(plink) != plink_sha, "PLINK2 SHA mismatch before analysis.")

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  header <- c("#chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls")
  schema <- paste(sprintf("'%s': 'VARCHAR'", header), collapse = ", ")
  src <- sprintf("read_csv(%s, header=true, delim='\\t', columns={%s}, auto_detect=false, strict_mode=true)", sql_string(fg), schema)
  scan_start <- ts()
  raw <- DBI::dbGetQuery(con, sprintf(
    "SELECT \"#chrom\" chr, pos, ref, alt, rsids, pval, beta, sebeta, af_alt FROM %s WHERE TRY_CAST(pval AS DOUBLE)>0 AND TRY_CAST(pval AS DOUBLE)<%.17g",
    src, p_threshold
  ))
  scan_end <- ts()
  log_line("single_streaming_scan raw_candidate_count=", nrow(raw))

  x <- data.frame(
    chromosome = as.integer(raw$chr),
    position = as.integer(raw$pos),
    rsid_raw = as.character(raw$rsids),
    ref = toupper(trimws(raw$ref)),
    alt = toupper(trimws(raw$alt)),
    beta = num(raw$beta),
    se = num(raw$sebeta),
    pval = num(raw$pval),
    eaf = num(raw$af_alt),
    stringsAsFactors = FALSE
  )
  x$rsid <- trimws(x$rsid_raw)
  x$finngen_log10p <- -log10(x$pval)
  x$F_statistic <- (x$beta / x$se)^2
  x$apoe_region <- x$chromosome == 19L & x$position >= 44000000L & x$position <= 46000000L
  reason <- rep("eligible", nrow(x))
  mark <- function(idx, label) reason[reason == "eligible" & idx] <<- label
  mark(is.na(x$chromosome) | x$chromosome < 1L | x$chromosome > 22L, "non_autosomal")
  mark(is.na(x$rsid) | x$rsid == "", "missing_rsid")
  mark(!is.na(x$rsid) & x$rsid != "" & !grepl("^rs[0-9]+$", x$rsid), "noncanonical_rsid")
  mark(is.na(x$ref) | is.na(x$alt) | nchar(x$ref) != 1L | nchar(x$alt) != 1L, "non_biallelic")
  mark(grepl("[^ACGT]", x$ref) | grepl("[^ACGT]", x$alt), "non_acgt")
  mark(x$ref == x$alt, "same_allele")
  mark(!is.finite(x$beta), "invalid_beta")
  mark(!is.finite(x$se) | x$se <= 0, "invalid_se")
  mark(!is.finite(x$pval) | x$pval <= 0 | x$pval >= p_threshold, "invalid_p")
  x$eligibility_status <- ifelse(reason == "eligible", "eligible_pre_reference", "ineligible")
  x$exclusion_reason <- ifelse(reason == "eligible", NA, reason)
  write_pq(con, x, paste0(out[["candidates"]], ".partial"))

  eligible <- x[reason == "eligible", , drop = FALSE]
  write_pq(con, eligible, paste0(out[["eligible"]], ".partial"))
  dup <- eligible[duplicated(eligible$rsid) | duplicated(eligible$rsid, fromLast = TRUE), , drop = FALSE]
  duplicate_rsid_count <- length(unique(dup$rsid))
  if (nrow(dup) > 0L) {
    write_tsv(dup, paste0(out[["dup"]], ".partial"))
    stop("Eligible duplicate rsID detected; no undocumented duplicate resolution applied.", call. = FALSE)
  }

  DBI::dbWriteTable(con, "eligible_ids", data.frame(rsid = eligible$rsid), temporary = TRUE, overwrite = TRUE)
  bim_src <- sprintf("read_csv(%s, header=false, delim='\\t', columns={'rch':'VARCHAR','rsid':'VARCHAR','cm':'VARCHAR','rpos':'VARCHAR','a1':'VARCHAR','a2':'VARCHAR'}, auto_detect=false, strict_mode=true)", sql_string(paste0(eur, ".bim")))
  ref <- DBI::dbGetQuery(con, sprintf("SELECT b.* FROM %s b INNER JOIN eligible_ids i ON b.rsid=i.rsid", bim_src))
  DBI::dbRemoveTable(con, "eligible_ids")
  eligible_ref <- merge(eligible, ref, by = "rsid", all.x = TRUE, sort = FALSE)
  not_ref <- eligible_ref[is.na(eligible_ref$rpos), c("rsid", "chromosome", "position", "pval", "finngen_log10p"), drop = FALSE]
  write_tsv(not_ref, paste0(out[["not_ref"]], ".partial"))
  present <- eligible_ref[!is.na(eligible_ref$rpos), , drop = FALSE]
  present$SNP <- present$rsid
  present$CHR <- as.integer(present$rch)
  present$BP <- as.integer(present$rpos)
  present$LOG10P <- present$finngen_log10p
  included_input <- present
  excluded_input <- present[!present$apoe_region, , drop = FALSE]
  apoe_split_ok <- nrow(excluded_input) == nrow(present) - sum(present$apoe_region) &&
    setequal(excluded_input$rsid, setdiff(present$rsid, present$rsid[present$apoe_region]))
  stop_if(!apoe_split_ok, "APOE split gate failed.")

  clump <- function(dat, label) {
    stop_if(nrow(dat) < 1L, paste("No SNPs available for", label, "clumping."))
    in_path <- tempfile(fileext = ".tsv")
    out_prefix <- tempfile()
    on.exit(unlink(c(in_path, paste0(out_prefix, c(".clumps", ".log"))), force = TRUE), add = TRUE)
    write_tsv(dat[, c("SNP", "CHR", "BP", "LOG10P", "alt")], in_path)
    plink_args <- c(
      "--bfile", eur,
      "--clump", in_path,
      "--clump-id-field", "SNP",
      "--clump-p-field", "LOG10P",
      "--clump-a1-field", "alt",
      "--clump-force-a1",
      "--clump-log10",
      "--clump-log10-p1", sprintf("%.15f", p_threshold_log10),
      "--clump-log10-p2", "0",
      "--clump-r2", "0.001",
      "--clump-kb", "10000",
      "--clump-unphased",
      "--threads", "8",
      "--memory", "8000",
      "--out", out_prefix
    )
    log_line("PLINK ", label, " binary=", plink)
    log_line("PLINK ", label, " binary_sha=", hash_file(plink))
    log_line("PLINK ", label, " command=", paste(c(shQuote(plink), shQuote(plink_args)), collapse = " "))
    status <- system2(plink, args = plink_args, stdout = TRUE, stderr = TRUE)
    exit_status <- attr(status, "status")
    exit_status <- if (is.null(exit_status)) 0L else as.integer(exit_status)
    log_line("PLINK ", label, " exit_status=", exit_status)
    log_line("PLINK ", label, " output=", paste(status, collapse = " | "))
    stop_if(exit_status != 0L, paste("PLINK failed for", label))
    stop_if(!file.exists(paste0(out_prefix, ".clumps")), paste("PLINK clumps output missing for", label))
    cl <- read.delim(paste0(out_prefix, ".clumps"), check.names = FALSE)
    id_col <- intersect(c("ID", "SNP"), names(cl))
    stop_if(length(id_col) != 1L, paste("PLINK clumps ID column not recognised for", label))
    ids <- as.character(cl[[id_col]])
    stop_if(anyDuplicated(ids) > 0L || !all(ids %in% dat$rsid), paste("PLINK index verification failed for", label))
    list(data = dat[match(ids, dat$rsid), , drop = FALSE], clumps = cl, ids = ids, command = paste(c(shQuote(plink), shQuote(plink_args)), collapse = " "), exit_status = exit_status)
  }

  inc_clump <- clump(included_input, "apoe_included")
  exc_clump <- clump(excluded_input, "apoe_excluded")
  clump_log10p <- function(clumped) {
    id_col <- intersect(c("ID", "SNP"), names(clumped$clumps))
    if ("LOG10P" %in% names(clumped$clumps)) as.numeric(clumped$clumps$LOG10P[match(clumped$ids, clumped$clumps[[id_col]])]) else rep(NA_real_, length(clumped$ids))
  }
  included <- final_schema(inc_clump$data, "apoe_included_exploratory_p5e-6", clump_log10p(inc_clump))
  excluded <- final_schema(exc_clump$data, "apoe_excluded_exploratory_p5e-6", clump_log10p(exc_clump))
  schema_audit <- rbind(schema_counts(included, "APOE included"), schema_counts(excluded, "APOE excluded"))
  write.csv(schema_audit, paste0(out[["schema"]], ".partial"), row.names = FALSE)
  stop_if(any(schema_audit$duplicate_column_count != 0L) || any(schema_audit$raw_LOG10P_count != 0L), "Final schema audit failed.")

  write_pq(con, included, paste0(out[["included_pq"]], ".partial"))
  write_tsv(included, paste0(out[["included_tsv"]], ".partial"))
  write_pq(con, excluded, paste0(out[["excluded_pq"]], ".partial"))
  write_tsv(excluded, paste0(out[["excluded_tsv"]], ".partial"))
  inc_consistency <- check_pq_tsv(con, paste0(out[["included_pq"]], ".partial"), paste0(out[["included_tsv"]], ".partial"))
  exc_consistency <- check_pq_tsv(con, paste0(out[["excluded_pq"]], ".partial"), paste0(out[["excluded_tsv"]], ".partial"))

  primary_inc <- read_pq(con, file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_primary_apoe_included_clumped_v4.parquet"))
  primary_exc <- read_pq(con, file.path(root, "data_derived", "reverse_instruments", "finngen_r13_delirium_primary_apoe_excluded_clumped_v4.parquet"))
  overlap <- data.frame(
    analysis_set = c("APOE included", "APOE excluded"),
    strict_primary_nsnp = c(nrow(primary_inc), nrow(primary_exc)),
    relaxed_exploratory_nsnp = c(nrow(included), nrow(excluded)),
    overlap_nsnp = c(length(intersect(primary_inc$rsid, included$rsid)), length(intersect(primary_exc$rsid, excluded$rsid))),
    relaxed_only_nsnp = c(length(setdiff(included$rsid, primary_inc$rsid)), length(setdiff(excluded$rsid, primary_exc$rsid))),
    strict_only_nsnp = c(length(setdiff(primary_inc$rsid, included$rsid)), length(setdiff(primary_exc$rsid, excluded$rsid))),
    overlap_rsids = c(paste(intersect(primary_inc$rsid, included$rsid), collapse = ";"), paste(intersect(primary_exc$rsid, excluded$rsid), collapse = ";")),
    stringsAsFactors = FALSE
  )
  write.csv(overlap, paste0(out[["overlap"]], ".partial"), row.names = FALSE)

  strength <- function(z) {
    list(
      nSNP = nrow(z),
      F_min = min(z$F_statistic),
      F_mean = mean(z$F_statistic),
      F_median = median(z$F_statistic),
      F_max = max(z$F_statistic),
      F_lt_10_count = sum(z$F_statistic < 10)
    )
  }
  inc_strength <- strength(included)
  exc_strength <- strength(excluded)
  shared <- intersect(included$rsid, excluded$rsid)
  exclusion_counts <- as.data.frame(table(reason), stringsAsFactors = FALSE)
  names(exclusion_counts) <- c("eligibility_or_exclusion_reason", "count")
  counts <- data.frame(
    metric = c(
      "raw_candidate_count", "eligible_candidate_count", "duplicate_rsid_count",
      "eligible_reference_present_count", "candidate_not_in_reference_count",
      "apoe_candidate_excluded_count", "included_clump_input_count", "excluded_clump_input_count",
      "included_plink_index_count", "excluded_plink_index_count", "included_nsnp", "excluded_nsnp",
      "shared_nsnp", "included_only_nsnp", "excluded_only_nsnp",
      "included_F_lt_10_count", "excluded_F_lt_10_count"
    ),
    value = c(
      nrow(x), nrow(eligible), duplicate_rsid_count,
      nrow(present), nrow(not_ref), sum(present$apoe_region),
      nrow(included_input), nrow(excluded_input), length(inc_clump$ids), length(exc_clump$ids),
      nrow(included), nrow(excluded), length(shared),
      length(setdiff(included$rsid, excluded$rsid)), length(setdiff(excluded$rsid, included$rsid)),
      inc_strength$F_lt_10_count, exc_strength$F_lt_10_count
    ),
    stringsAsFactors = FALSE
  )
  write.csv(counts, paste0(out[["counts"]], ".partial"), row.names = FALSE)

  finngen_after <- hash_file(fg)
  renv_after <- hash_file(renv_lock)
  hard_check_failures <- character()
  if (!inc_consistency) hard_check_failures <- c(hard_check_failures, "included_parquet_tsv_consistency_failed")
  if (!exc_consistency) hard_check_failures <- c(hard_check_failures, "excluded_parquet_tsv_consistency_failed")
  if (finngen_after != finngen_sha) hard_check_failures <- c(hard_check_failures, "finngen_sha_after_mismatch")
  if (renv_before != renv_after) hard_check_failures <- c(hard_check_failures, "renv_lock_changed")
  if (!all(unlist(restoration$hard_checks))) hard_check_failures <- c(hard_check_failures, names(restoration$hard_checks)[!unlist(restoration$hard_checks)])
  if (!all(file.exists(v1_required))) hard_check_failures <- c(hard_check_failures, "v1_failure_artifacts_preserved_failed")
  if (!inc_strength$F_lt_10_count == 0L) hard_check_failures <- c(hard_check_failures, "included_F_lt_10_present")
  if (!exc_strength$F_lt_10_count == 0L) hard_check_failures <- c(hard_check_failures, "excluded_F_lt_10_present")
  if (duplicate_rsid_count != 0L) hard_check_failures <- c(hard_check_failures, "duplicate_rsid_present")
  if (!apoe_split_ok) hard_check_failures <- c(hard_check_failures, "apoe_split_set_identity_failed")
  status <- if (length(hard_check_failures) == 0L) "passed" else "failed"
  qc <- list(
    instrument_selection_version = "v2",
    decision = 58,
    supersedes_failed_version = "v1",
    v1_failure_reason = "missing_project_local_plink2_binary_before_plink_clumping",
    v1_failure_artifacts_preserved = all(file.exists(v1_required)),
    v2_not_using_v1_partial = TRUE,
    analysis_direction = "delirium_to_Hb",
    analysis_role = "secondary_reverse_exploratory_relaxed",
    strict_primary_instrument_version = "v4",
    trigger_freeze_version = "v1",
    trigger_freeze_decision = 56,
    p_threshold = p_threshold,
    threshold_label = "P < 5e-6",
    exploratory_only = TRUE,
    ld_r2 = 0.001,
    ld_window_kb = 10000,
    ld_reference = "1000G_Phase3_EUR",
    genome_build = "FinnGen GRCh38 exposure coordinates; exact rsID match to EUR LD reference",
    apoe_policy = "APOE included in candidate universe; separate APOE-included and APOE-excluded clumping branches",
    apoe_region = "chr19:44000000-46000000",
    plink_restoration = restoration,
    plink_included_command = inc_clump$command,
    plink_excluded_command = exc_clump$command,
    plink_included_exit_status = inc_clump$exit_status,
    plink_excluded_exit_status = exc_clump$exit_status,
    source_finngen_sha256_before = finngen_before,
    source_finngen_sha256_after = finngen_after,
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after,
    renv_lock_unchanged = identical(renv_before, renv_after),
    raw_candidate_count = nrow(x),
    eligibility_exclusion_counts = exclusion_counts,
    eligible_candidate_count = nrow(eligible),
    duplicate_rsid_count = duplicate_rsid_count,
    eligible_reference_present_count = nrow(present),
    candidate_not_in_reference_count = nrow(not_ref),
    apoe_candidate_excluded_count = sum(present$apoe_region),
    apoe_split_set_identity_ok = apoe_split_ok,
    included_clump_input_count = nrow(included_input),
    excluded_clump_input_count = nrow(excluded_input),
    included_plink_index_count = length(inc_clump$ids),
    excluded_plink_index_count = length(exc_clump$ids),
    included_nsnp = nrow(included),
    excluded_nsnp = nrow(excluded),
    shared_nsnp = length(shared),
    included_only_nsnp = length(setdiff(included$rsid, excluded$rsid)),
    excluded_only_nsnp = length(setdiff(excluded$rsid, included$rsid)),
    instrument_strength_included = inc_strength,
    instrument_strength_excluded = exc_strength,
    weak_instruments_reported_not_removed = TRUE,
    included_parquet_tsv_consistency = inc_consistency,
    excluded_parquet_tsv_consistency = exc_consistency,
    proxy_used = FALSE,
    liftover_used = FALSE,
    relaxed_harmonisation_started = FALSE,
    reverse_mr_started = FALSE,
    strict_primary_results_modified = FALSE,
    hard_checks = c(restoration$hard_checks, list(
      v1_failure_artifacts_preserved = all(file.exists(v1_required)),
      v2_not_using_v1_partial = TRUE,
      finngen_sha_before_matches_expected = identical(finngen_before, finngen_sha),
      finngen_sha_after_matches_expected = identical(finngen_after, finngen_sha),
      renv_lock_unchanged = identical(renv_before, renv_after),
      duplicate_rsid_absent = duplicate_rsid_count == 0L,
      included_F_lt_10_absent = inc_strength$F_lt_10_count == 0L,
      excluded_F_lt_10_absent = exc_strength$F_lt_10_count == 0L,
      apoe_split_set_identity_ok = apoe_split_ok,
      included_parquet_tsv_consistency = inc_consistency,
      excluded_parquet_tsv_consistency = exc_consistency
    )),
    hard_check_failures = hard_check_failures,
    instrument_selection_status = status,
    approved_for_reverse_relaxed_outcome_extraction = identical(status, "passed"),
    informational_findings = list(
      scan_start = scan_start,
      scan_end = scan_end,
      renv_out_of_sync_warning_may_be_emitted_by_project_activation = TRUE,
      schema_audit = schema_audit,
      primary_overlap_audit = overlap
    )
  )
  jsonlite::write_json(qc, paste0(out[["json"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  if (!identical(status, "passed")) stop("Relaxed instrument selection V2 failed; partial outputs retained.", call. = FALSE)
  for (path in out[c("candidates", "eligible", "included_pq", "included_tsv", "excluded_pq", "excluded_tsv", "counts", "not_ref", "schema", "overlap", "json")]) {
    stop_if(file.exists(path), paste("Output appeared during run:", path))
    stop_if(!file.rename(paste0(path, ".partial"), path), paste("Atomic rename failed:", path))
  }
  log_line("instrument_selection_status=passed included_nsnp=", nrow(included), " excluded_nsnp=", nrow(excluded), " hard_check_failures=0")
}

tryCatch(main(), error = function(e) {
  log_line("status=failed error=", conditionMessage(e))
  quit(status = 1L)
})
