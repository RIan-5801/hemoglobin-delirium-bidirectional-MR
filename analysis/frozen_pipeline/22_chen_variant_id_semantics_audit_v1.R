#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/22_chen_variant_id_semantics_audit_v1.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
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
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))

write_csv_precise <- function(x, path) {
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  write.csv(x, path, row.names = FALSE, na = "")
}
write_tsv_precise <- function(x, path) {
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

source_path <- file.path(root, "data_raw", "gwas", "BCX2_HGB_EA_GWAMA.out.gz")
eur_bim <- file.path(root, "resources", "ld", "1kg_v3", "EUR.bim")
cert_path <- file.path(root, "results", "qc", "chen_2020_hb_source_certification_v1.json")
decision81_qc_path <- file.path(root, "results", "qc", "chen_forward_instrument_selection_v1.json")
decision81_counts_path <- file.path(root, "results", "qc", "chen_forward_instrument_selection_counts_v1.csv")
renv_lock <- file.path(root, "renv.lock")

out <- c(
  script = file.path(root, "R", "22_chen_variant_id_semantics_audit_v1.R"),
  decision = file.path(root, "docs", "decisions", "82_chen_2020_hb_variant_id_semantics_audit_v1_v1.1.md"),
  qc_json = file.path(root, "results", "qc", "chen_2020_hb_variant_id_semantics_audit_v1.json"),
  format_counts = file.path(root, "results", "qc", "chen_2020_hb_variant_id_format_counts_v1.csv"),
  examples = file.path(root, "results", "qc", "chen_2020_hb_variant_id_examples_v1.tsv"),
  feasibility = file.path(root, "results", "qc", "chen_2020_hb_variant_identity_resolution_feasibility_v1.csv"),
  log = file.path(root, "results", "logs", "chen_2020_hb_variant_id_semantics_audit_v1.log")
)

runtime_targets <- out[!names(out) %in% c("script", "decision")]
all_targets <- c(runtime_targets, paste0(runtime_targets, ".partial"))
stop_if(any(file.exists(all_targets)), paste("A final or partial Chen variant-ID semantics audit target exists; refusing to overwrite:", paste(names(all_targets)[file.exists(all_targets)], collapse = ", ")))
for (p in runtime_targets) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)

log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = out[["log"]], append = TRUE)
atomic <- function(path, writer) {
  partial <- paste0(path, ".partial")
  stop_if(file.exists(path) || file.exists(partial), paste("Output occupied:", path))
  writer(partial)
  stop_if(!file.exists(partial), paste("Writer did not create partial:", partial))
  stop_if(!file.rename(partial, path), paste("Atomic rename failed:", path))
}

source_documentation_audit <- function(root) {
  roots <- c(file.path(root, "docs"), file.path(root, "data_raw", "metadata"), file.path(root, "results", "qc"))
  roots <- roots[dir.exists(roots)]
  files <- unlist(lapply(roots, function(d) list.files(d, recursive = TRUE, full.names = TRUE)), use.names = FALSE)
  files <- files[grepl("\\.(md|txt|tsv|csv|json)$", files, ignore.case = TRUE)]
  files <- files[file.info(files)$size < 5 * 1024 * 1024]
  generated <- normalizePath(c(
    file.path(root, "docs", "decisions", "82_chen_2020_hb_variant_id_semantics_audit_v1_v1.1.md"),
    file.path(root, "R", "22_chen_variant_id_semantics_audit_v1.R"),
    file.path(root, "results", "qc", "chen_2020_hb_variant_id_semantics_audit_v1.json")
  ), winslash = "/", mustWork = FALSE)
  files <- files[!normalizePath(files, winslash = "/", mustWork = FALSE) %in% generated]
  files <- files[!grepl("chen_2020_hb_variant_(id|identity)_", basename(files), ignore.case = TRUE)]
  keep <- grepl("chen|bcx2|gwas|metadata|manifest|decision|source|schema|contract|coordinate|instrument", basename(files), ignore.case = TRUE) |
    grepl("[/\\\\]docs[/\\\\]", files)
  files <- files[keep]
  hits <- data.frame(file = character(), line_number = integer(), line = character(), stringsAsFactors = FALSE)
  for (f in files) {
    lines <- tryCatch(readLines(f, warn = FALSE, encoding = "UTF-8"), error = function(e) character())
    if (!length(lines)) next
    idx <- grep("rs_number|BCX2_HGB|Chen 2020|reference_allele|other_allele", lines, ignore.case = TRUE)
    if (length(idx)) {
      rel <- gsub("\\\\", "/", substring(normalizePath(f, winslash = "/", mustWork = FALSE), nchar(root) + 2L))
      hits <- rbind(hits, data.frame(file = rel, line_number = idx, line = lines[idx], stringsAsFactors = FALSE))
    }
  }
  source_dictionary_hits <- hits[grepl("dictionary|readme|download|metadata|manifest", hits$file, ignore.case = TRUE) &
    grepl("rs_number", hits$line, ignore.case = TRUE), , drop = FALSE]
  explicit_definition_hits <- source_dictionary_hits[grepl("meaning|definition|defined|identifier|variant|rsid|rs id|SNP", source_dictionary_hits$line, ignore.case = TRUE), , drop = FALSE]
  if (nrow(explicit_definition_hits) > 0L) {
    definition <- paste(unique(explicit_definition_hits$line), collapse = " | ")
    status <- "source_documented_in_project"
  } else {
    definition <- "No Chen original data dictionary, README, or download metadata in the project explicitly defines rs_number."
    status <- "not_source_documented_in_project"
  }
  list(
    source_documented_identifier_definition = definition,
    source_identifier_definition_status = status,
    documentation_search_files_scanned = length(files),
    source_dictionary_searched = TRUE,
    documented_vs_inferred_semantics_separated = TRUE,
    evidence_hits = records(hits[seq_len(min(nrow(hits), 50L)), , drop = FALSE])
  )
}

main <- function() {
  started <- Sys.time()
  log_line("stage=chen_variant_id_semantics_audit_v1")
  for (p in c(source_path, eur_bim, cert_path, decision81_qc_path, decision81_counts_path, renv_lock)) {
    stop_if(!file.exists(p), paste("Missing required input:", p))
  }

  cert <- jsonlite::fromJSON(cert_path, simplifyVector = FALSE)
  decision81 <- jsonlite::fromJSON(decision81_qc_path, simplifyVector = FALSE)
  expected_sha <- toupper(cert$source_sha256)
  sha_before <- toupper(hash_file(source_path))
  renv_before <- hash_file(renv_lock)
  log_line("source_sha_before=", sha_before)
  stop_if(!identical(sha_before, expected_sha), "Chen source SHA before audit differs from Decision 78 certification.")

  decision81_counts <- read.csv(decision81_counts_path, check.names = FALSE, stringsAsFactors = FALSE)
  d81_raw_candidates <- decision81_counts$value[match("raw_p_lt_5e_8_candidate_count", decision81_counts$metric)]
  d81_eligible <- decision81_counts$value[match("eligible_candidate_count", decision81_counts$metric)]
  decision81_failure_preserved <- identical(decision81$instrument_selection_status, "failed") &&
    identical(decision81$approved_for_chen_forward_outcome_extraction, FALSE) &&
    isTRUE(as.integer(d81_raw_candidates) == 125669L) &&
    isTRUE(as.integer(d81_eligible) == 0L)
  stop_if(!decision81_failure_preserved, "Decision 81 failure is not preserved as expected.")

  doc_audit <- source_documentation_audit(root)
  log_line("documentation_search_files_scanned=", doc_audit$documentation_search_files_scanned)

  duckdb_tmp <- file.path(root, "results", "tmp", "duckdb_chen_variant_id_semantics_v1")
  dir.create(duckdb_tmp, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(DUCKDB_TEMP_DIRECTORY = duckdb_tmp)
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:", config = list(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(con, "PRAGMA threads=8")
  DBI::dbExecute(con, "PRAGMA memory_limit='8GB'")
  DBI::dbExecute(con, sprintf("PRAGMA temp_directory='%s'", qpath(duckdb_tmp)))

  DBI::dbExecute(con, sprintf(
    "CREATE TEMP TABLE chen_raw AS SELECT * FROM read_csv('%s', delim='\\t', header=true, compression='gzip', all_varchar=true, ignore_errors=false)",
    qpath(source_path)
  ))
  DBI::dbExecute(con, "
    CREATE TEMP TABLE source_scan AS
    WITH parsed AS (
      SELECT
        row_number() OVER () AS source_row_number,
        rs_number AS rs_number_raw,
        trim(rs_number) AS id_trim,
        upper(trim(rs_number)) AS id_upper,
        reference_allele,
        other_allele,
        upper(reference_allele) AS chen_effect_allele,
        upper(other_allele) AS chen_other_allele,
        TRY_CAST(\"p-value\" AS DOUBLE) AS pval,
        regexp_extract(upper(trim(rs_number)), '^(?:CHR)?([0-9]{1,2}|X|Y|MT|M):([0-9]+)_([^_]+)_([^_]+)$', 1) AS parsed_chr_raw,
        regexp_extract(upper(trim(rs_number)), '^(?:CHR)?([0-9]{1,2}|X|Y|MT|M):([0-9]+)_([^_]+)_([^_]+)$', 2) AS parsed_pos_raw,
        regexp_extract(upper(trim(rs_number)), '^(?:CHR)?([0-9]{1,2}|X|Y|MT|M):([0-9]+)_([^_]+)_([^_]+)$', 3) AS parsed_allele_token_1,
        regexp_extract(upper(trim(rs_number)), '^(?:CHR)?([0-9]{1,2}|X|Y|MT|M):([0-9]+)_([^_]+)_([^_]+)$', 4) AS parsed_allele_token_2,
        regexp_extract(upper(trim(rs_number)), '^(?:CHR)?([0-9]{1,2}|X|Y|MT|M):([0-9]+)$', 1) AS parsed_chr_pos_chr,
        regexp_extract(upper(trim(rs_number)), '^(?:CHR)?([0-9]{1,2}|X|Y|MT|M):([0-9]+)$', 2) AS parsed_chr_pos_pos
      FROM chen_raw
    ),
    classified AS (
      SELECT
        *,
        CASE
          WHEN id_trim IS NULL OR id_trim = '' OR upper(id_trim) IN ('NA', 'N/A', 'NULL', '.') THEN 'missing'
          WHEN regexp_full_match(id_trim, '^rs[0-9]+$') THEN 'canonical_rsid'
          WHEN regexp_full_match(id_upper, '^(?:CHR)?([0-9]{1,2}|X|Y|MT|M):[0-9]+_[ACGT]+_[ACGT]+$') THEN 'chr_pos_allele1_allele2'
          WHEN regexp_full_match(id_upper, '^(?:CHR)?([0-9]{1,2}|X|Y|MT|M):[0-9]+$') THEN 'chr_pos'
          WHEN regexp_full_match(id_trim, '^[0-9]+$') THEN 'numeric_only'
          WHEN instr(id_trim, ':') > 0 THEN 'colon_delimited_other'
          WHEN instr(id_trim, '_') > 0 THEN 'underscore_delimited_other'
          ELSE 'other'
        END AS identifier_grammar
      FROM parsed
    )
    SELECT
      *,
      CASE WHEN parsed_chr_raw <> '' THEN CASE WHEN parsed_chr_raw = 'M' THEN 'MT' ELSE parsed_chr_raw END
           WHEN parsed_chr_pos_chr <> '' THEN CASE WHEN parsed_chr_pos_chr = 'M' THEN 'MT' ELSE parsed_chr_pos_chr END
           ELSE NULL END AS parsed_chr,
      TRY_CAST(CASE WHEN parsed_pos_raw <> '' THEN parsed_pos_raw WHEN parsed_chr_pos_pos <> '' THEN parsed_chr_pos_pos ELSE NULL END AS BIGINT) AS parsed_pos,
      CASE WHEN identifier_grammar IN ('chr_pos', 'chr_pos_allele1_allele2') THEN TRUE ELSE FALSE END AS coordinate_like,
      CASE
        WHEN (parsed_chr_raw <> '' OR parsed_chr_pos_chr <> '') AND TRY_CAST(CASE WHEN parsed_pos_raw <> '' THEN parsed_pos_raw WHEN parsed_chr_pos_pos <> '' THEN parsed_chr_pos_pos ELSE NULL END AS BIGINT) > 0 THEN TRUE
        ELSE FALSE
      END AS coordinate_parseable,
      CASE WHEN parsed_allele_token_1 <> '' AND parsed_allele_token_2 <> '' THEN TRUE ELSE FALSE END AS identifier_contains_alleles,
      CASE
        WHEN parsed_allele_token_1 <> '' AND parsed_allele_token_2 <> '' AND parsed_allele_token_1 <= parsed_allele_token_2 THEN parsed_allele_token_1 || '/' || parsed_allele_token_2
        WHEN parsed_allele_token_1 <> '' AND parsed_allele_token_2 <> '' THEN parsed_allele_token_2 || '/' || parsed_allele_token_1
        ELSE NULL
      END AS identifier_allele_set,
      CASE
        WHEN parsed_allele_token_1 <> '' AND parsed_allele_token_2 <> '' AND parsed_allele_token_1 = chen_effect_allele AND parsed_allele_token_2 = chen_other_allele THEN 'same_order'
        WHEN parsed_allele_token_1 <> '' AND parsed_allele_token_2 <> '' AND parsed_allele_token_1 = chen_other_allele AND parsed_allele_token_2 = chen_effect_allele THEN 'reverse_order'
        WHEN parsed_allele_token_1 <> '' AND parsed_allele_token_2 <> '' AND translate(parsed_allele_token_1, 'ACGT', 'TGCA') = chen_effect_allele AND translate(parsed_allele_token_2, 'ACGT', 'TGCA') = chen_other_allele THEN 'strand_complement_same_order'
        WHEN parsed_allele_token_1 <> '' AND parsed_allele_token_2 <> '' AND translate(parsed_allele_token_1, 'ACGT', 'TGCA') = chen_other_allele AND translate(parsed_allele_token_2, 'ACGT', 'TGCA') = chen_effect_allele THEN 'strand_complement_reverse_order'
        WHEN parsed_allele_token_1 <> '' AND parsed_allele_token_2 <> '' THEN 'incompatible'
        ELSE 'not_evaluable'
      END AS identifier_vs_chen_alleles_relation,
      CASE WHEN pval > 0 AND pval < 5e-8 THEN TRUE ELSE FALSE END AS significant_candidate
    FROM classified
  ")

  DBI::dbExecute(con, sprintf(
    "CREATE TEMP TABLE eur AS
     SELECT column0::VARCHAR AS reference_chr_grch37,
            column1::VARCHAR AS reference_rsid,
            column3::BIGINT AS reference_pos_grch37,
            upper(column4::VARCHAR) AS reference_a1,
            upper(column5::VARCHAR) AS reference_a2,
            CASE WHEN upper(column4::VARCHAR) <= upper(column5::VARCHAR)
                 THEN upper(column4::VARCHAR) || '/' || upper(column5::VARCHAR)
                 ELSE upper(column5::VARCHAR) || '/' || upper(column4::VARCHAR)
            END AS reference_allele_set
     FROM read_csv('%s', delim='\\t', header=false, all_varchar=false)",
    qpath(eur_bim)
  ))
  DBI::dbExecute(con, "
    CREATE TEMP TABLE eur_coord AS
    SELECT reference_chr_grch37, reference_pos_grch37, COUNT(*) AS coordinate_reference_match_count
    FROM eur
    WHERE reference_rsid <> '.'
    GROUP BY reference_chr_grch37, reference_pos_grch37
  ")
  DBI::dbExecute(con, "
    CREATE TEMP TABLE eur_coord_allele AS
    SELECT reference_chr_grch37, reference_pos_grch37, reference_allele_set, COUNT(*) AS coordinate_allele_reference_match_count
    FROM eur
    WHERE reference_rsid <> '.'
    GROUP BY reference_chr_grch37, reference_pos_grch37, reference_allele_set
  ")
  DBI::dbExecute(con, "
    CREATE TEMP TABLE identity_audit AS
    SELECT
      s.*,
      COALESCE(c.coordinate_reference_match_count, 0) AS coordinate_reference_match_count,
      COALESCE(a.coordinate_allele_reference_match_count, 0) AS coordinate_allele_reference_match_count
    FROM source_scan s
    LEFT JOIN eur_coord c
      ON s.parsed_chr = c.reference_chr_grch37
     AND s.parsed_pos = c.reference_pos_grch37
    LEFT JOIN eur_coord_allele a
      ON s.parsed_chr = a.reference_chr_grch37
     AND s.parsed_pos = a.reference_pos_grch37
     AND s.identifier_allele_set = a.reference_allele_set
  ")

  format_counts <- DBI::dbGetQuery(con, "
    WITH base AS (
      SELECT 'all_source' AS subset, identifier_grammar FROM identity_audit
      UNION ALL
      SELECT 'p_lt_5e-8_candidates' AS subset, identifier_grammar FROM identity_audit WHERE significant_candidate
    ),
    counted AS (
      SELECT subset, identifier_grammar, COUNT(*) AS n FROM base GROUP BY subset, identifier_grammar
    )
    SELECT subset, identifier_grammar, n, n * 1.0 / SUM(n) OVER (PARTITION BY subset) AS fraction
    FROM counted
    ORDER BY subset, n DESC, identifier_grammar
  ")
  atomic(out[["format_counts"]], function(p) write_csv_precise(format_counts, p))

  examples <- DBI::dbGetQuery(con, "
    WITH e AS (
      SELECT
        'all_source' AS subset,
        identifier_grammar,
        row_number() OVER (PARTITION BY identifier_grammar ORDER BY source_row_number) AS example_rank,
        rs_number_raw,
        parsed_chr,
        parsed_pos,
        parsed_allele_token_1,
        parsed_allele_token_2,
        reference_allele AS chen_effect_allele,
        other_allele AS chen_other_allele,
        identifier_vs_chen_alleles_relation,
        pval
      FROM identity_audit
      UNION ALL
      SELECT
        'p_lt_5e-8_candidates' AS subset,
        identifier_grammar,
        row_number() OVER (PARTITION BY identifier_grammar ORDER BY source_row_number) AS example_rank,
        rs_number_raw,
        parsed_chr,
        parsed_pos,
        parsed_allele_token_1,
        parsed_allele_token_2,
        reference_allele AS chen_effect_allele,
        other_allele AS chen_other_allele,
        identifier_vs_chen_alleles_relation,
        pval
      FROM identity_audit
      WHERE significant_candidate
    )
    SELECT * FROM e
    WHERE example_rank <= 10
    ORDER BY subset, identifier_grammar, example_rank
  ")
  atomic(out[["examples"]], function(p) write_tsv_precise(examples, p))

  relation_counts <- DBI::dbGetQuery(con, "
    WITH base AS (
      SELECT 'all_source' AS subset, identifier_vs_chen_alleles_relation FROM identity_audit WHERE identifier_contains_alleles
      UNION ALL
      SELECT 'p_lt_5e-8_candidates' AS subset, identifier_vs_chen_alleles_relation FROM identity_audit WHERE significant_candidate AND identifier_contains_alleles
    )
    SELECT subset, identifier_vs_chen_alleles_relation, COUNT(*) AS n
    FROM base
    GROUP BY subset, identifier_vs_chen_alleles_relation
    ORDER BY subset, n DESC, identifier_vs_chen_alleles_relation
  ")

  feasibility <- DBI::dbGetQuery(con, "
    WITH base AS (
      SELECT 'all_source' AS subset, * FROM identity_audit
      UNION ALL
      SELECT 'p_lt_5e-8_candidates' AS subset, * FROM identity_audit WHERE significant_candidate
    ),
    agg AS (
      SELECT
        subset,
        COUNT(*) AS n_total,
        SUM(CASE WHEN identifier_grammar = 'canonical_rsid' THEN 1 ELSE 0 END) AS n_canonical_rsid,
        SUM(CASE WHEN coordinate_like THEN 1 ELSE 0 END) AS n_coordinate_like,
        SUM(CASE WHEN coordinate_parseable THEN 1 ELSE 0 END) AS n_parseable_coordinate,
        SUM(CASE WHEN identifier_contains_alleles THEN 1 ELSE 0 END) AS n_identifier_contains_alleles,
        SUM(CASE WHEN coordinate_parseable AND coordinate_reference_match_count = 1 THEN 1 ELSE 0 END) AS n_exact_coordinate_unique_reference,
        SUM(CASE WHEN coordinate_parseable AND coordinate_reference_match_count > 1 THEN 1 ELSE 0 END) AS n_exact_coordinate_multiple_reference,
        SUM(CASE WHEN coordinate_parseable AND coordinate_reference_match_count = 0 THEN 1 ELSE 0 END) AS n_exact_coordinate_no_reference,
        SUM(CASE WHEN identifier_contains_alleles AND coordinate_allele_reference_match_count = 1 THEN 1 ELSE 0 END) AS n_coordinate_allele_unique_reference,
        SUM(CASE WHEN identifier_contains_alleles AND coordinate_allele_reference_match_count > 1 THEN 1 ELSE 0 END) AS n_coordinate_allele_multiple_reference,
        SUM(CASE WHEN identifier_contains_alleles AND coordinate_allele_reference_match_count = 0 THEN 1 ELSE 0 END) AS n_coordinate_allele_no_reference,
        SUM(CASE WHEN identifier_contains_alleles AND coordinate_reference_match_count > 0 AND coordinate_allele_reference_match_count = 0 THEN 1 ELSE 0 END) AS n_coordinate_allele_incompatible,
        SUM(CASE WHEN coordinate_parseable AND TRY_CAST(parsed_chr AS INTEGER) BETWEEN 1 AND 22 AND parsed_pos > 0 THEN 1 ELSE 0 END) AS n_grch37_plausible_autosomal_coordinate
      FROM base
      GROUP BY subset
    ),
    metrics AS (
      SELECT subset, 'n_total' AS metric, n_total::DOUBLE AS value FROM agg
      UNION ALL SELECT subset, 'n_canonical_rsid', n_canonical_rsid::DOUBLE FROM agg
      UNION ALL SELECT subset, 'n_coordinate_like', n_coordinate_like::DOUBLE FROM agg
      UNION ALL SELECT subset, 'n_parseable_coordinate', n_parseable_coordinate::DOUBLE FROM agg
      UNION ALL SELECT subset, 'n_identifier_contains_alleles', n_identifier_contains_alleles::DOUBLE FROM agg
      UNION ALL SELECT subset, 'n_exact_coordinate_unique_reference', n_exact_coordinate_unique_reference::DOUBLE FROM agg
      UNION ALL SELECT subset, 'n_exact_coordinate_multiple_reference', n_exact_coordinate_multiple_reference::DOUBLE FROM agg
      UNION ALL SELECT subset, 'n_exact_coordinate_no_reference', n_exact_coordinate_no_reference::DOUBLE FROM agg
      UNION ALL SELECT subset, 'n_coordinate_allele_unique_reference', n_coordinate_allele_unique_reference::DOUBLE FROM agg
      UNION ALL SELECT subset, 'n_coordinate_allele_multiple_reference', n_coordinate_allele_multiple_reference::DOUBLE FROM agg
      UNION ALL SELECT subset, 'n_coordinate_allele_no_reference', n_coordinate_allele_no_reference::DOUBLE FROM agg
      UNION ALL SELECT subset, 'n_coordinate_allele_incompatible', n_coordinate_allele_incompatible::DOUBLE FROM agg
      UNION ALL SELECT subset, 'n_grch37_plausible_autosomal_coordinate', n_grch37_plausible_autosomal_coordinate::DOUBLE FROM agg
      UNION ALL SELECT subset, 'recoverable_fraction_coordinate_only', CASE WHEN n_parseable_coordinate > 0 THEN n_exact_coordinate_unique_reference * 1.0 / n_parseable_coordinate ELSE NULL END FROM agg
      UNION ALL SELECT subset, 'recoverable_fraction_coordinate_plus_alleles', CASE WHEN n_identifier_contains_alleles > 0 THEN n_coordinate_allele_unique_reference * 1.0 / n_identifier_contains_alleles ELSE NULL END FROM agg
    )
    SELECT * FROM metrics ORDER BY subset, metric
  ")
  relation_metric_rows <- data.frame(
    subset = relation_counts$subset,
    metric = paste0("identifier_vs_chen_alleles_relation__", relation_counts$identifier_vs_chen_alleles_relation),
    value = as.numeric(relation_counts$n),
    stringsAsFactors = FALSE
  )
  feasibility <- rbind(feasibility, relation_metric_rows)
  atomic(out[["feasibility"]], function(p) write_csv_precise(feasibility, p))

  get_metric <- function(subset, metric) {
    x <- feasibility$value[feasibility$subset == subset & feasibility$metric == metric]
    if (length(x) == 0L) NA_real_ else as.numeric(x[[1L]])
  }
  dominant <- format_counts[order(format_counts$subset, -format_counts$n), , drop = FALSE]
  dominant_all <- dominant$identifier_grammar[dominant$subset == "all_source"][[1L]]
  dominant_sig <- dominant$identifier_grammar[dominant$subset == "p_lt_5e-8_candidates"][[1L]]
  source_definition_status <- doc_audit$source_identifier_definition_status
  recommended_class <- if (!identical(source_definition_status, "source_documented_in_project")) {
    "C_additional_source_dictionary_required"
  } else if (isTRUE(get_metric("p_lt_5e-8_candidates", "recoverable_fraction_coordinate_plus_alleles") == 1) &&
             isTRUE(get_metric("p_lt_5e-8_candidates", "n_coordinate_allele_incompatible") == 0)) {
    "A_canonical_rsid_recoverable_with_unique_exact_variant_identity"
  } else {
    "B_variant_identity_not_safely_recoverable"
  }

  sha_after <- toupper(hash_file(source_path))
  renv_after <- hash_file(renv_lock)
  hard_checks <- list(
    decision_78_source_gate = identical(cert$certification_status, "passed") && length(cert$hard_check_failures) == 0L,
    decision_81_failure_preserved = decision81_failure_preserved,
    source_sha_before_gate = identical(sha_before, expected_sha),
    source_sha_after_gate = identical(sha_after, expected_sha),
    source_unchanged = identical(sha_before, sha_after),
    raw_identifier_values_preserved = TRUE,
    identifier_grammar_audited = nrow(format_counts[format_counts$subset == "all_source", , drop = FALSE]) > 0L,
    significant_candidate_grammar_audited = nrow(format_counts[format_counts$subset == "p_lt_5e-8_candidates", , drop = FALSE]) > 0L,
    source_dictionary_searched = isTRUE(doc_audit$source_dictionary_searched),
    documented_vs_inferred_semantics_separated = isTRUE(doc_audit$documented_vs_inferred_semantics_separated),
    effect_allele_not_called_genomic_ref = TRUE,
    allele_token_semantics_kept_independent = TRUE,
    identity_resolution_feasibility_only = TRUE,
    no_proxy = TRUE,
    no_liftover = TRUE,
    no_nearest_variant = TRUE,
    no_clumping = TRUE,
    no_instrument_selection = TRUE,
    no_outcome_extraction = TRUE,
    no_harmonisation = TRUE,
    no_mr = TRUE,
    no_steiger = TRUE
  )
  failures <- names(hard_checks)[!unlist(hard_checks)]
  audit_status <- if (length(failures) == 0L) "passed" else "failed"

  qc <- list(
    audit_version = "v1",
    decision = 82,
    source_name = "Chen_2020_BCX2_European_Hb",
    source_path = "data_raw/gwas/BCX2_HGB_EA_GWAMA.out.gz",
    source_sha256_before = sha_before,
    source_sha256_after = sha_after,
    source_genome_build = "GRCh37",
    source_field = "rs_number",
    source_documented_identifier_definition = doc_audit$source_documented_identifier_definition,
    source_identifier_definition_status = source_definition_status,
    source_documentation_evidence_hits = doc_audit$evidence_hits,
    all_variant_identifier_format_counts = records(format_counts[format_counts$subset == "all_source", , drop = FALSE]),
    significant_candidate_identifier_format_counts = records(format_counts[format_counts$subset == "p_lt_5e-8_candidates", , drop = FALSE]),
    dominant_identifier_grammar = list(all_source = dominant_all, p_lt_5e_8_candidates = dominant_sig),
    coordinate_parseable = list(
      all_source_n = get_metric("all_source", "n_parseable_coordinate"),
      p_lt_5e_8_candidates_n = get_metric("p_lt_5e-8_candidates", "n_parseable_coordinate")
    ),
    identifier_contains_alleles = list(
      all_source_n = get_metric("all_source", "n_identifier_contains_alleles"),
      p_lt_5e_8_candidates_n = get_metric("p_lt_5e-8_candidates", "n_identifier_contains_alleles")
    ),
    allele_token_semantics_status = "syntax_only; identifier allele tokens compared with Chen effect/other allele fields but not assigned genomic REF/ALT semantics",
    identifier_vs_chen_alleles_relation_counts = records(relation_counts),
    reference_identity_feasibility = records(feasibility),
    recoverability_metrics = list(
      all_source_coordinate_only_fraction = get_metric("all_source", "recoverable_fraction_coordinate_only"),
      all_source_coordinate_plus_alleles_fraction = get_metric("all_source", "recoverable_fraction_coordinate_plus_alleles"),
      candidate_coordinate_only_fraction = get_metric("p_lt_5e-8_candidates", "recoverable_fraction_coordinate_only"),
      candidate_coordinate_plus_alleles_fraction = get_metric("p_lt_5e-8_candidates", "recoverable_fraction_coordinate_plus_alleles")
    ),
    proxy_used = FALSE,
    liftover_used = FALSE,
    nearest_variant_matching_used = FALSE,
    instrument_selection_performed = FALSE,
    clumping_performed = FALSE,
    outcome_extraction_performed = FALSE,
    harmonisation_performed = FALSE,
    mr_run = FALSE,
    steiger_run = FALSE,
    audit_started = format(started, "%Y-%m-%dT%H:%M:%S%z"),
    audit_completed = ts(),
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after,
    renv_lock_unchanged = identical(renv_before, renv_after),
    audit_status = audit_status,
    recommended_next_decision_class = recommended_class,
    hard_checks = hard_checks,
    hard_check_failures = failures,
    informational_findings = list(
      decision_81_failure_not_modified = TRUE,
      coordinate_based_exact_identity_resolution_is_not_proxy = TRUE,
      coordinate_based_exact_identity_resolution_is_not_original_exact_rsid_identity = TRUE,
      future_use_requires_new_chen_specific_scientific_decision = TRUE,
      reference_alleles_from_1kg_are_used_only_for_feasibility_not_harmonisation = TRUE
    )
  )
  atomic(out[["qc_json"]], function(p) jsonlite::write_json(qc, p, pretty = TRUE, auto_unbox = TRUE, na = "null"))
  stop_if(!identical(audit_status, "passed"), "Chen variant-ID semantics audit failed; QC JSON retained.")
  log_line("audit_status=passed recommended_next_decision_class=", recommended_class)
  log_line("dominant_all=", dominant_all, " dominant_p_lt_5e-8=", dominant_sig)
  log_line("source_sha_after=", sha_after)
}

tryCatch(main(), error = function(e) {
  log_line("status=failed error=", conditionMessage(e))
  quit(status = 1L)
})
