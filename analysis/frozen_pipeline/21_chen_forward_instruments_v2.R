#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/21_chen_forward_instruments_v2.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
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
allele_set <- function(a, b) mapply(function(x, y) paste(sort(c(toupper(x), toupper(y))), collapse = "/"), a, b, USE.NAMES = FALSE)
comp_base <- function(x) chartr("ACGT", "TGCA", toupper(x))
comp_set <- function(a, b) allele_set(comp_base(a), comp_base(b))
is_acgt <- function(x) grepl("^[ACGT]$", toupper(x))
num_equal <- function(a, b, tol = 1e-12) {
  a <- as.numeric(a); b <- as.numeric(b)
  length(a) == length(b) && all((is.na(a) & is.na(b)) | (is.finite(a) & is.finite(b) & abs(a - b) <= tol))
}
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
atomic <- function(path, writer) {
  partial <- paste0(path, ".partial")
  stop_if(file.exists(path) || file.exists(partial), paste("Output occupied:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writer(partial)
  stop_if(!file.exists(partial), paste("Writer did not create partial:", partial))
  stop_if(!file.rename(partial, path), paste("Atomic rename failed:", path))
}

source_path <- file.path(root, "data_raw", "gwas", "BCX2_HGB_EA_GWAMA.out.gz")
eur_prefix <- file.path(root, "resources", "ld", "1kg_v3", "EUR")
eur_bim <- paste0(eur_prefix, ".bim")
plink <- file.path(root, "tools", "plink2", "plink2.exe")
cert_path <- file.path(root, "results", "qc", "chen_2020_hb_source_certification_v1.json")
contract_path <- file.path(root, "results", "qc", "chen_forward_sensitivity_instrument_selection_contract_v1.json")
coord_path <- file.path(root, "results", "qc", "chen_2020_hb_grch37_coordinate_authority_v1.json")
official_dict_path <- file.path(root, "results", "qc", "chen_2020_hb_official_source_dictionary_audit_v1.json")
amendment_path <- file.path(root, "results", "qc", "chen_identifier_resolution_amendment_v1.json")
decision81_path <- file.path(root, "docs", "decisions", "81_chen_forward_instrument_selection_v1_v1.1.md")
decision82_format_path <- file.path(root, "results", "qc", "chen_2020_hb_variant_id_format_counts_v1.csv")
decision82_resolution_path <- file.path(root, "results", "qc", "chen_2020_hb_variant_identity_resolution_feasibility_v1.csv")
vuckovic_inst <- file.path(root, "data_derived", "instruments", "vuckovic_hb_instruments_apoe_included_v2.tsv")
renv_lock <- file.path(root, "renv.lock")

out <- c(
  script = file.path(root, "R", "21_chen_forward_instruments_v2.R"),
  candidates = file.path(root, "data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_candidates_v2.parquet"),
  identity_resolution = file.path(root, "data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_identity_resolution_v2.parquet"),
  resolved_eligible = file.path(root, "data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_resolved_eligible_v2.parquet"),
  included_input = file.path(root, "data_derived", "forward_sensitivity_instruments", "plink", "chen_2020_hb_p5e-8_apoe_included_clump_input_v2.tsv"),
  excluded_input = file.path(root, "data_derived", "forward_sensitivity_instruments", "plink", "chen_2020_hb_p5e-8_apoe_excluded_clump_input_v2.tsv"),
  included_clumps = file.path(root, "data_derived", "forward_sensitivity_instruments", "plink", "chen_2020_hb_p5e-8_apoe_included_clumped_v2.clumps"),
  excluded_clumps = file.path(root, "data_derived", "forward_sensitivity_instruments", "plink", "chen_2020_hb_p5e-8_apoe_excluded_clumped_v2.clumps"),
  included_plink_log = file.path(root, "data_derived", "forward_sensitivity_instruments", "plink", "chen_2020_hb_p5e-8_apoe_included_clumped_v2.log"),
  excluded_plink_log = file.path(root, "data_derived", "forward_sensitivity_instruments", "plink", "chen_2020_hb_p5e-8_apoe_excluded_clumped_v2.log"),
  included_missing_id = file.path(root, "data_derived", "forward_sensitivity_instruments", "plink", "chen_2020_hb_p5e-8_apoe_included_clumped_v2.clumps.missing_id"),
  excluded_missing_id = file.path(root, "data_derived", "forward_sensitivity_instruments", "plink", "chen_2020_hb_p5e-8_apoe_excluded_clumped_v2.clumps.missing_id"),
  included_missing_allele = file.path(root, "data_derived", "forward_sensitivity_instruments", "plink", "chen_2020_hb_p5e-8_apoe_included_clumped_v2.clumps.missing_allele"),
  excluded_missing_allele = file.path(root, "data_derived", "forward_sensitivity_instruments", "plink", "chen_2020_hb_p5e-8_apoe_excluded_clumped_v2.clumps.missing_allele"),
  included_inst_parquet = file.path(root, "data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_included_clumped_v2.parquet"),
  included_inst_tsv = file.path(root, "data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_included_clumped_v2.tsv"),
  excluded_inst_parquet = file.path(root, "data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_excluded_clumped_v2.parquet"),
  excluded_inst_tsv = file.path(root, "data_derived", "forward_sensitivity_instruments", "chen_2020_hb_p5e-8_apoe_excluded_clumped_v2.tsv"),
  resolution_counts = file.path(root, "results", "qc", "chen_forward_identifier_resolution_counts_v2.csv"),
  resolution_failures = file.path(root, "results", "qc", "chen_forward_identifier_resolution_failures_v2.tsv"),
  duplicate_audit = file.path(root, "results", "qc", "chen_forward_resolved_duplicate_rsid_audit_v2.tsv"),
  counts = file.path(root, "results", "qc", "chen_forward_instrument_selection_counts_v2.csv"),
  overlap = file.path(root, "results", "qc", "chen_forward_vuckovic_overlap_audit_v2.csv"),
  qc_json = file.path(root, "results", "qc", "chen_forward_instrument_selection_v2.json"),
  log = file.path(root, "results", "logs", "chen_forward_instrument_selection_v2.log"),
  decision = file.path(root, "docs", "decisions", "85_chen_forward_instrument_selection_v2_v1.1.md")
)

target_names <- setdiff(names(out), "script")
all_targets <- c(out[target_names], paste0(out[target_names], ".partial"))
stop_if(any(file.exists(all_targets)), paste("A final or partial Chen forward V2 target exists; refusing to overwrite:", paste(names(all_targets)[file.exists(all_targets)], collapse = ", ")))
for (p in out[target_names]) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)

log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = out[["log"]], append = TRUE)

read_clumps <- function(path) {
  x <- read.table(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE, comment.char = "")
  norm <- toupper(gsub("[^A-Za-z0-9]", "", names(x)))
  hit <- match(c("ID", "SNP", "INDEXSNP"), norm, nomatch = 0L)
  hit <- hit[hit > 0L]
  stop_if(length(hit) == 0L, paste("Cannot identify index SNP column:", path))
  x$index_snp <- as.character(x[[hit[[1L]]]])
  x
}

validate_roundtrip <- function(con, parquet_path, tsv_path, original) {
  p <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", qpath(parquet_path)))
  t <- read.delim(tsv_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  same_cols <- identical(names(p), names(t)) && identical(names(p), names(original))
  same_n <- nrow(p) == nrow(t) && nrow(p) == nrow(original)
  same_rsid <- identical(as.character(p$resolved_rsid), as.character(t$resolved_rsid)) && identical(as.character(p$resolved_rsid), as.character(original$resolved_rsid))
  char_cols <- names(original)[vapply(original, function(z) is.character(z) || is.logical(z), logical(1))]
  num_cols <- names(original)[vapply(original, is.numeric, logical(1))]
  char_ok <- all(vapply(char_cols, function(k) identical(as.character(p[[k]]), as.character(t[[k]])) && identical(as.character(p[[k]]), as.character(original[[k]])), logical(1)))
  num_ok <- all(vapply(num_cols, function(k) num_equal(p[[k]], t[[k]]) && num_equal(p[[k]], original[[k]]), logical(1)))
  list(row_count = nrow(p), same_cols = same_cols, same_n = same_n, same_rsid = same_rsid, char_ok = char_ok, num_ok = num_ok)
}

count_missing <- function(path) {
  if (!file.exists(path)) return(0L)
  sum(nzchar(trimws(readLines(path, warn = FALSE))))
}

write_decision85 <- function(summary) {
  lines <- c(
    "# Decision 85 - Chen forward instrument selection V2",
    "",
    "Date: 2026-08-12",
    "",
    "## Decision",
    "",
    "Execute Chen 2020 Hb to FinnGen R13 `F5_DELIRIUM` forward alternative-Hb",
    "GWAS sensitivity instrument selection V2 under Decision 84.",
    "",
    "## Rationale",
    "",
    "Decision 84 freezes a Chen-specific rule permitting official coordinate marker",
    "IDs to be resolved by exact GRCh37 coordinate plus unordered allele-set identity",
    "in the certified 1000 Genomes Phase 3 EUR reference. V2 therefore starts from",
    "the raw Chen source rather than from Decision 81's canonical-rsID failure set.",
    "",
    "## Scope and Safeguards",
    "",
    "This stage extracts Chen records with `0 < p-value < 5e-8`, parses official",
    "marker IDs, resolves only unique exact coordinate-plus-allele identities, audits",
    "duplicate resolved rsIDs, checks Chen marker/effect allele compatibility without",
    "beta flipping, creates APOE-included and APOE-excluded inputs before clumping,",
    "and runs independent PLINK clumping for both sets.",
    "",
    "This stage does not run FinnGen outcome extraction, harmonisation, MR, Steiger,",
    "proxy search, LD proxy search, liftOver, nearest-variant matching, fuzzy",
    "matching, or strand-complement identity rescue. Chen remains an alternative-Hb",
    "GWAS sensitivity analysis and not an independent replication.",
    "",
    "## Results Recorded",
    "",
    sprintf("- Source SHA before: `%s`", summary$source_sha256_before),
    sprintf("- Source SHA after: `%s`", summary$source_sha256_after),
    sprintf("- Raw P<5e-8 candidate count: `%s`", summary$raw_candidate_count),
    sprintf("- Pre-resolution eligible count: `%s`", summary$pre_resolution_eligible_count),
    sprintf("- Unique exact identity count: `%s`", summary$unique_exact_identity_count),
    sprintf("- Resolved duplicate rsID count: `%s`", summary$resolved_duplicate_rsid_count),
    sprintf("- Resolved eligible count: `%s`", summary$resolved_eligible_count),
    sprintf("- APOE-window candidate count: `%s`", summary$apoe_window_candidate_count),
    sprintf("- Included clump input count: `%s`", summary$included_clump_input_count),
    sprintf("- Excluded clump input count: `%s`", summary$excluded_clump_input_count),
    sprintf("- Included final nSNP: `%s`", summary$included_nsnp),
    sprintf("- Excluded final nSNP: `%s`", summary$excluded_nsnp),
    sprintf("- Shared / included-only / excluded-only: `%s / %s / %s`", summary$shared_nsnp, summary$included_only_nsnp, summary$excluded_only_nsnp),
    sprintf("- Instrument-selection status: `%s`", summary$instrument_selection_status),
    sprintf("- Approved for Chen forward outcome extraction: `%s`", summary$approved_for_chen_forward_outcome_extraction),
    "",
    "## Affected Files",
    "",
    "- `R/21_chen_forward_instruments_v2.R`",
    "- `data_derived/forward_sensitivity_instruments/chen_2020_hb_p5e-8_candidates_v2.parquet`",
    "- `data_derived/forward_sensitivity_instruments/chen_2020_hb_p5e-8_identity_resolution_v2.parquet`",
    "- `data_derived/forward_sensitivity_instruments/chen_2020_hb_p5e-8_resolved_eligible_v2.parquet`",
    "- `data_derived/forward_sensitivity_instruments/chen_2020_hb_p5e-8_apoe_included_clumped_v2.parquet`",
    "- `data_derived/forward_sensitivity_instruments/chen_2020_hb_p5e-8_apoe_included_clumped_v2.tsv`",
    "- `data_derived/forward_sensitivity_instruments/chen_2020_hb_p5e-8_apoe_excluded_clumped_v2.parquet`",
    "- `data_derived/forward_sensitivity_instruments/chen_2020_hb_p5e-8_apoe_excluded_clumped_v2.tsv`",
    "- `results/qc/chen_forward_identifier_resolution_counts_v2.csv`",
    "- `results/qc/chen_forward_identifier_resolution_failures_v2.tsv`",
    "- `results/qc/chen_forward_resolved_duplicate_rsid_audit_v2.tsv`",
    "- `results/qc/chen_forward_instrument_selection_counts_v2.csv`",
    "- `results/qc/chen_forward_vuckovic_overlap_audit_v2.csv`",
    "- `results/qc/chen_forward_instrument_selection_v2.json`",
    "- `results/logs/chen_forward_instrument_selection_v2.log`",
    "",
    "## Expected Impact",
    "",
    "This creates Chen forward sensitivity instruments for a later, separately",
    "approved FinnGen targeted outcome extraction. It does not alter existing",
    "Vuckovic/FinnGen results or Decision 81 failure evidence."
  )
  atomic(out[["decision"]], function(p) writeLines(lines, p, useBytes = TRUE))
}

main <- function() {
  started <- Sys.time()
  log_line("stage=chen_forward_instrument_selection_v2")

  required_inputs <- c(
    source_path, eur_bim, paste0(eur_prefix, ".bed"), paste0(eur_prefix, ".fam"),
    plink, cert_path, contract_path, coord_path, official_dict_path, amendment_path,
    decision81_path, decision82_format_path, decision82_resolution_path, renv_lock
  )
  for (p in required_inputs) stop_if(!file.exists(p), paste("Missing required input:", p))

  cert <- jsonlite::fromJSON(cert_path, simplifyVector = FALSE)
  contract <- jsonlite::fromJSON(contract_path, simplifyVector = FALSE)
  coord <- jsonlite::fromJSON(coord_path, simplifyVector = FALSE)
  official <- jsonlite::fromJSON(official_dict_path, simplifyVector = FALSE)
  amendment <- jsonlite::fromJSON(amendment_path, simplifyVector = FALSE)

  expected_sha <- toupper(cert$source_sha256)
  sha_before <- toupper(hash_file(source_path))
  renv_before <- hash_file(renv_lock)
  plink_sha <- hash_file(plink)
  plink_version <- paste(system2(plink, "--version", stdout = TRUE), collapse = " ")
  log_line("source_sha_before=", sha_before)
  log_line("plink_version=", plink_version)
  log_line("plink_sha256=", plink_sha)

  stop_if(!identical(sha_before, expected_sha), "Chen source SHA before scan differs from certification.")
  stop_if(!identical(cert$certification_status, "passed") || length(cert$hard_check_failures) != 0L, "Decision 78 gate failed.")
  stop_if(!identical(contract$contract_status, "frozen") || length(contract$hard_check_failures) != 0L, "Decision 79 gate failed.")
  stop_if(!identical(coord$coordinate_authority_status, "passed") || length(coord$hard_check_failures) != 0L, "Decision 80 gate failed.")
  stop_if(!identical(official$audit_status, "passed") || length(official$hard_check_failures) != 0L, "Decision 83 gate failed.")
  stop_if(!identical(amendment$amendment_status, "frozen") || !isTRUE(amendment$approved_for_chen_forward_instrument_selection_v2) || length(amendment$hard_check_failures) != 0L, "Decision 84 amendment gate failed.")

  duckdb_tmp <- file.path(root, "results", "tmp", "duckdb_chen_forward_v2")
  dir.create(duckdb_tmp, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(DUCKDB_TEMP_DIRECTORY = duckdb_tmp)
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:", config = list(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbExecute(con, "PRAGMA threads=8")
  DBI::dbExecute(con, "PRAGMA memory_limit='8GB'")
  DBI::dbExecute(con, sprintf("PRAGMA temp_directory='%s'", qpath(duckdb_tmp)))

  expected_header <- c("rs_number", "reference_allele", "other_allele", "eaf", "beta", "se", "beta_95L", "beta_95U", "z", "p-value", "_-log10_p-value", "q_statistic", "q-p-value", "i2", "n_studies", "n_samples", "effects")
  certified_header <- unlist(cert$raw_header, use.names = FALSE)
  header_line <- readLines(gzfile(source_path, "rt"), n = 1L, warn = FALSE)
  observed_header <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  if (!identical(observed_header, expected_header) && identical(observed_header, certified_header)) {
    expected_header <- certified_header
  }
  stop_if(!identical(observed_header, expected_header), "Chen source header differs from Decision 78.")

  DBI::dbExecute(con, sprintf(
    "CREATE TEMP TABLE chen_raw AS SELECT * FROM read_csv('%s', delim='\\t', header=true, compression='gzip', all_varchar=true, ignore_errors=false)",
    qpath(source_path)
  ))
  DBI::dbExecute(con, "CREATE TEMP TABLE chen_raw_indexed AS SELECT row_number() OVER () AS source_row_number, * FROM chen_raw")
  DBI::dbExecute(con, sprintf(
    "CREATE TEMP TABLE eur AS SELECT column0::VARCHAR AS reference_chr_grch37, column1::VARCHAR AS resolved_rsid, column3::BIGINT AS reference_pos_grch37, upper(column4::VARCHAR) AS reference_a1, upper(column5::VARCHAR) AS reference_a2 FROM read_csv('%s', delim='\\t', header=false, all_varchar=false)",
    qpath(eur_bim)
  ))
  eur_dups <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM (SELECT resolved_rsid FROM eur WHERE resolved_rsid <> '.' GROUP BY resolved_rsid HAVING COUNT(*) > 1)")$n[[1]]
  stop_if(as.integer(eur_dups) != 0L, "EUR reference contains duplicate rsIDs; stopping.")

  DBI::dbExecute(con, "
    CREATE TEMP TABLE candidates AS
    SELECT
      source_row_number,
      rs_number AS source_marker_id,
      reference_allele AS exposure_effect_allele,
      other_allele AS exposure_other_allele,
      TRY_CAST(eaf AS DOUBLE) AS eaf,
      TRY_CAST(beta AS DOUBLE) AS beta,
      TRY_CAST(se AS DOUBLE) AS se,
      TRY_CAST(\"p-value\" AS DOUBLE) AS pval,
      TRY_CAST(\"_-log10_p-value\" AS DOUBLE) AS log10p,
      TRY_CAST(n_studies AS BIGINT) AS n_studies,
      TRY_CAST(n_samples AS DOUBLE) AS n_samples,
      effects AS effects,
      TRY_CAST(beta AS DOUBLE) * TRY_CAST(beta AS DOUBLE) / (TRY_CAST(se AS DOUBLE) * TRY_CAST(se AS DOUBLE)) AS F_stat,
      regexp_extract(rs_number, '^([0-9]+):([0-9]+)_([A-Za-z]+)_([A-Za-z]+)$', 1)::VARCHAR AS marker_chr,
      TRY_CAST(regexp_extract(rs_number, '^([0-9]+):([0-9]+)_([A-Za-z]+)_([A-Za-z]+)$', 2) AS BIGINT) AS marker_pos,
      upper(regexp_extract(rs_number, '^([0-9]+):([0-9]+)_([A-Za-z]+)_([A-Za-z]+)$', 3))::VARCHAR AS marker_allele1,
      upper(regexp_extract(rs_number, '^([0-9]+):([0-9]+)_([A-Za-z]+)_([A-Za-z]+)$', 4))::VARCHAR AS marker_allele2
    FROM chen_raw_indexed
    WHERE TRY_CAST(\"p-value\" AS DOUBLE) > 0 AND TRY_CAST(\"p-value\" AS DOUBLE) < 5e-8
  ")
  DBI::dbExecute(con, "
    CREATE TEMP TABLE candidates_eligible AS
    SELECT *,
      CASE
        WHEN source_marker_id IS NULL OR NOT regexp_full_match(source_marker_id, '^[0-9]+:[0-9]+_[A-Za-z]+_[A-Za-z]+$') THEN 'marker_parse_failure'
        WHEN marker_chr NOT IN ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22') THEN 'non_autosomal'
        WHEN marker_pos IS NULL OR marker_pos <= 0 THEN 'marker_parse_failure'
        WHEN NOT regexp_full_match(marker_allele1, '^[ACGT]$') OR NOT regexp_full_match(marker_allele2, '^[ACGT]$') THEN 'non_snp_marker'
        WHEN marker_allele1 = marker_allele2 THEN 'identical_marker_alleles'
        WHEN beta IS NULL OR NOT isfinite(beta) THEN 'beta_nonfinite'
        WHEN se IS NULL OR NOT isfinite(se) OR se <= 0 THEN 'se_nonpositive'
        WHEN pval IS NULL OR NOT isfinite(pval) OR NOT (pval > 0 AND pval < 5e-8) THEN 'p_not_strictly_lt_5e-8'
        WHEN n_samples IS NULL OR NOT isfinite(n_samples) OR n_samples <= 0 THEN 'variant_level_n_missing_or_nonpositive'
        ELSE 'pre_resolution_eligible'
      END AS pre_resolution_status,
      CASE WHEN marker_allele1 < marker_allele2 THEN marker_allele1 || '/' || marker_allele2 ELSE marker_allele2 || '/' || marker_allele1 END AS marker_allele_set,
      CASE WHEN upper(exposure_effect_allele) < upper(exposure_other_allele) THEN upper(exposure_effect_allele) || '/' || upper(exposure_other_allele) ELSE upper(exposure_other_allele) || '/' || upper(exposure_effect_allele) END AS chen_effect_allele_set,
      CASE WHEN marker_chr = '19' AND marker_pos BETWEEN 44000000 AND 46000000 THEN TRUE ELSE FALSE END AS apoe_region
    FROM candidates
  ")
  DBI::dbExecute(con, sprintf("COPY candidates_eligible TO '%s' (FORMAT PARQUET)", qpath(paste0(out[["candidates"]], ".partial"))))
  stop_if(!file.rename(paste0(out[["candidates"]], ".partial"), out[["candidates"]]), "Candidate parquet promotion failed.")

  DBI::dbExecute(con, "
    CREATE TEMP TABLE eur_coord_counts AS
    SELECT reference_chr_grch37, reference_pos_grch37, COUNT(*) AS coordinate_match_count
    FROM eur
    GROUP BY reference_chr_grch37, reference_pos_grch37
  ")
  DBI::dbExecute(con, "
    CREATE TEMP TABLE exact_matches AS
    SELECT
      c.source_marker_id,
      e.resolved_rsid,
      e.reference_chr_grch37,
      e.reference_pos_grch37,
      e.reference_a1,
      e.reference_a2,
      CASE WHEN e.reference_a1 < e.reference_a2 THEN e.reference_a1 || '/' || e.reference_a2 ELSE e.reference_a2 || '/' || e.reference_a1 END AS reference_allele_set
    FROM candidates_eligible c
    JOIN eur e
      ON c.marker_chr = e.reference_chr_grch37
     AND c.marker_pos = e.reference_pos_grch37
     AND c.marker_allele_set = CASE WHEN e.reference_a1 < e.reference_a2 THEN e.reference_a1 || '/' || e.reference_a2 ELSE e.reference_a2 || '/' || e.reference_a1 END
    WHERE c.pre_resolution_status = 'pre_resolution_eligible'
  ")
  DBI::dbExecute(con, "
    CREATE TEMP TABLE exact_match_counts AS
    SELECT source_marker_id, COUNT(*) AS exact_match_count
    FROM exact_matches
    GROUP BY source_marker_id
  ")
  DBI::dbExecute(con, "
    CREATE TEMP TABLE first_exact_match AS
    SELECT *
    FROM (
      SELECT *, row_number() OVER (PARTITION BY source_marker_id ORDER BY resolved_rsid) AS rn
      FROM exact_matches
    )
    WHERE rn = 1
  ")
  DBI::dbExecute(con, "
    CREATE TEMP TABLE resolution AS
    SELECT
      c.*,
      COALESCE(cc.coordinate_match_count, 0) AS coordinate_match_count,
      COALESCE(emc.exact_match_count, 0) AS exact_match_count,
      fem.resolved_rsid,
      fem.reference_chr_grch37,
      fem.reference_pos_grch37,
      fem.reference_a1,
      fem.reference_a2,
      fem.reference_allele_set,
      CASE
        WHEN c.pre_resolution_status <> 'pre_resolution_eligible' THEN c.pre_resolution_status
        WHEN COALESCE(cc.coordinate_match_count, 0) = 0 THEN 'coordinate_missing'
        WHEN COALESCE(cc.coordinate_match_count, 0) > 1 AND COALESCE(emc.exact_match_count, 0) = 0 THEN 'coordinate_multiple'
        WHEN COALESCE(emc.exact_match_count, 0) = 0 THEN 'coordinate_match_allele_incompatible'
        WHEN COALESCE(emc.exact_match_count, 0) = 1 THEN 'unique_exact_coordinate_allele_identity'
        ELSE 'multiple_exact_allele_match'
      END AS identity_resolution_status,
      CASE
        WHEN c.pre_resolution_status <> 'pre_resolution_eligible' THEN 'not_evaluated_pre_resolution_ineligible'
        WHEN c.chen_effect_allele_set = c.marker_allele_set THEN 'same_unordered_set'
        WHEN (CASE WHEN replace(replace(replace(replace(c.exposure_effect_allele, 'A', 't'), 'T', 'a'), 'C', 'g'), 'G', 'c') < replace(replace(replace(replace(c.exposure_other_allele, 'A', 't'), 'T', 'a'), 'C', 'g'), 'G', 'c') THEN upper(replace(replace(replace(replace(c.exposure_effect_allele, 'A', 't'), 'T', 'a'), 'C', 'g'), 'G', 'c')) || '/' || upper(replace(replace(replace(replace(c.exposure_other_allele, 'A', 't'), 'T', 'a'), 'C', 'g'), 'G', 'c')) ELSE upper(replace(replace(replace(replace(c.exposure_other_allele, 'A', 't'), 'T', 'a'), 'C', 'g'), 'G', 'c')) || '/' || upper(replace(replace(replace(replace(c.exposure_effect_allele, 'A', 't'), 'T', 'a'), 'C', 'g'), 'G', 'c')) END) = c.marker_allele_set THEN 'strand_complement_set_not_rescued'
        ELSE 'incompatible'
      END AS chen_marker_effect_allele_compatibility
    FROM candidates_eligible c
    LEFT JOIN eur_coord_counts cc
      ON c.marker_chr = cc.reference_chr_grch37
     AND c.marker_pos = cc.reference_pos_grch37
    LEFT JOIN exact_match_counts emc
      ON c.source_marker_id = emc.source_marker_id
    LEFT JOIN first_exact_match fem
      ON c.source_marker_id = fem.source_marker_id
  ")
  DBI::dbExecute(con, sprintf("COPY resolution TO '%s' (FORMAT PARQUET)", qpath(paste0(out[["identity_resolution"]], ".partial"))))
  stop_if(!file.rename(paste0(out[["identity_resolution"]], ".partial"), out[["identity_resolution"]]), "Identity-resolution parquet promotion failed.")

  DBI::dbExecute(con, "
    CREATE TEMP TABLE resolved_predup AS
    SELECT *
    FROM resolution
    WHERE pre_resolution_status = 'pre_resolution_eligible'
      AND identity_resolution_status = 'unique_exact_coordinate_allele_identity'
      AND chen_marker_effect_allele_compatibility = 'same_unordered_set'
      AND resolved_rsid IS NOT NULL
  ")
  duplicate_audit <- DBI::dbGetQuery(con, "
    SELECT *
    FROM resolved_predup
    WHERE resolved_rsid IN (
      SELECT resolved_rsid FROM resolved_predup GROUP BY resolved_rsid HAVING COUNT(*) > 1
    )
    ORDER BY resolved_rsid, pval, source_marker_id
  ")
  atomic(out[["duplicate_audit"]], function(p) {
    if (nrow(duplicate_audit) == 0L) {
      write_tsv_precise(data.frame(note = "no_duplicate_resolved_rsids"), p)
    } else {
      write_tsv_precise(duplicate_audit, p)
    }
  })
  resolved_duplicate_rsid_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM (SELECT resolved_rsid FROM resolved_predup GROUP BY resolved_rsid HAVING COUNT(*) > 1)")$n[[1]]
  if (as.integer(resolved_duplicate_rsid_count) != 0L) {
    stop("Duplicate resolved rsIDs detected; stopping before clumping.", call. = FALSE)
  }

  DBI::dbExecute(con, "CREATE TEMP TABLE resolved_eligible AS SELECT * FROM resolved_predup")
  DBI::dbExecute(con, sprintf("COPY resolved_eligible TO '%s' (FORMAT PARQUET)", qpath(paste0(out[["resolved_eligible"]], ".partial"))))
  stop_if(!file.rename(paste0(out[["resolved_eligible"]], ".partial"), out[["resolved_eligible"]]), "Resolved-eligible parquet promotion failed.")

  resolution_breakdown <- DBI::dbGetQuery(con, "SELECT identity_resolution_status, COUNT(*) AS n FROM resolution GROUP BY identity_resolution_status ORDER BY identity_resolution_status")
  compatibility_breakdown <- DBI::dbGetQuery(con, "SELECT chen_marker_effect_allele_compatibility, COUNT(*) AS n FROM resolution GROUP BY chen_marker_effect_allele_compatibility ORDER BY chen_marker_effect_allele_compatibility")
  raw_candidate_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM candidates")$n[[1]]
  pre_resolution_eligible_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM candidates_eligible WHERE pre_resolution_status='pre_resolution_eligible'")$n[[1]]
  unique_exact_identity_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM resolution WHERE identity_resolution_status='unique_exact_coordinate_allele_identity'")$n[[1]]
  resolved_eligible_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM resolved_eligible")$n[[1]]
  apoe_window_candidate_count <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM resolved_eligible WHERE apoe_region")$n[[1]]
  resolution_counts <- data.frame(
    metric = c(
      "raw_candidate_count", "pre_resolution_eligible_count", "marker_parse_failure_count",
      "non_snp_marker_count", "coordinate_missing_count", "coordinate_allele_incompatible_count",
      "coordinate_multiple_count", "multiple_exact_allele_match_count", "unique_exact_identity_count",
      "chen_effect_marker_allele_incompatible_count", "chen_effect_marker_strand_complement_not_rescued_count",
      "resolved_duplicate_rsid_count", "resolved_eligible_count", "reference_resolution_fraction"
    ),
    value = c(
      as.numeric(raw_candidate_count), as.numeric(pre_resolution_eligible_count),
      DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM resolution WHERE identity_resolution_status='marker_parse_failure'")$n[[1]],
      DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM resolution WHERE identity_resolution_status='non_snp_marker'")$n[[1]],
      DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM resolution WHERE identity_resolution_status='coordinate_missing'")$n[[1]],
      DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM resolution WHERE identity_resolution_status='coordinate_match_allele_incompatible'")$n[[1]],
      DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM resolution WHERE identity_resolution_status='coordinate_multiple'")$n[[1]],
      DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM resolution WHERE identity_resolution_status='multiple_exact_allele_match'")$n[[1]],
      as.numeric(unique_exact_identity_count),
      DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM resolution WHERE chen_marker_effect_allele_compatibility='incompatible' AND identity_resolution_status='unique_exact_coordinate_allele_identity'")$n[[1]],
      DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM resolution WHERE chen_marker_effect_allele_compatibility='strand_complement_set_not_rescued' AND identity_resolution_status='unique_exact_coordinate_allele_identity'")$n[[1]],
      as.numeric(resolved_duplicate_rsid_count),
      as.numeric(resolved_eligible_count),
      if (as.numeric(pre_resolution_eligible_count) > 0) as.numeric(unique_exact_identity_count) / as.numeric(pre_resolution_eligible_count) else NA_real_
    ),
    stringsAsFactors = FALSE
  )
  atomic(out[["resolution_counts"]], function(p) write_csv_precise(resolution_counts, p))
  resolution_failures <- DBI::dbGetQuery(con, "
    SELECT *
    FROM resolution
    WHERE NOT (
      pre_resolution_status = 'pre_resolution_eligible'
      AND identity_resolution_status = 'unique_exact_coordinate_allele_identity'
      AND chen_marker_effect_allele_compatibility = 'same_unordered_set'
    )
    ORDER BY identity_resolution_status, source_marker_id
  ")
  atomic(out[["resolution_failures"]], function(p) write_tsv_precise(resolution_failures, p))

  resolved_eligible <- DBI::dbGetQuery(con, "SELECT * FROM resolved_eligible ORDER BY pval, resolved_rsid")
  stop_if(nrow(resolved_eligible) == 0L, "No resolved eligible candidates; stopping before clumping.")
  stop_if(anyDuplicated(resolved_eligible$resolved_rsid), "Duplicate resolved rsID after duplicate audit.")
  stop_if(any(resolved_eligible$F_stat < 10), "F<10 detected before clumping; stopping for manual review.")

  clump_cols <- c("resolved_rsid", "pval", "source_marker_id", "marker_chr", "marker_pos", "marker_allele1", "marker_allele2", "exposure_effect_allele", "exposure_other_allele", "beta", "se", "eaf", "n_samples", "F_stat", "apoe_region")
  included_input <- resolved_eligible[, clump_cols]
  included_input$SNP <- included_input$resolved_rsid
  included_input$P <- included_input$pval
  excluded_input <- included_input[!included_input$apoe_region, , drop = FALSE]
  stop_if(nrow(excluded_input) == 0L, "APOE-excluded clump input is empty.")
  atomic(out[["included_input"]], function(p) write_tsv_precise(included_input, p))
  atomic(out[["excluded_input"]], function(p) write_tsv_precise(excluded_input, p))

  run_plink <- function(label, input_path, out_clumps) {
    prefix <- sub("\\.clumps$", "", out_clumps)
    args <- c(
      "--bfile", eur_prefix,
      "--clump", input_path,
      "--clump-id-field", "SNP",
      "--clump-p-field", "P",
      "--clump-p1", "5e-8",
      "--clump-p2", "1",
      "--clump-r2", "0.001",
      "--clump-kb", "10000",
      "--clump-unphased",
      "--threads", "8",
      "--memory", "8000",
      "--out", prefix
    )
    log_line("plink_", label, "_command=", paste(c(plink, args), collapse = " "))
    status <- system2(plink, args = args, stdout = TRUE, stderr = TRUE)
    exit_code <- attr(status, "status")
    if (is.null(exit_code)) exit_code <- 0L
    log_line("plink_", label, "_exit_code=", exit_code)
    if (length(status)) for (line in status) log_line("plink_", label, "_stdout=", line)
    stop_if(exit_code != 0L, paste("PLINK failed for", label))
    stop_if(!file.exists(out_clumps), paste("PLINK did not create .clumps for", label))
  }
  run_plink("included", out[["included_input"]], out[["included_clumps"]])
  run_plink("excluded", out[["excluded_input"]], out[["excluded_clumps"]])

  make_instruments <- function(clumps_path, input_df, role) {
    cl <- read_clumps(clumps_path)
    stop_if(anyDuplicated(cl$index_snp), paste("Duplicate index SNP in", role))
    i <- match(cl$index_snp, input_df$resolved_rsid)
    stop_if(any(is.na(i)), paste("Clumped index SNP missing from input for", role))
    z <- input_df[i, , drop = FALSE]
    z$analysis_set <- role
    z$analysis_role <- if (role == "APOE_included") "forward_alternative_hb_gwas_sensitivity_main" else "APOE_region_exclusion_sensitivity"
    z$reference_chr_grch37 <- z$marker_chr
    z$reference_pos_grch37 <- z$marker_pos
    z$identity_resolution_status <- "unique_exact_coordinate_allele_identity"
    z$chen_marker_effect_allele_compatibility <- "same_unordered_set"
    z[, c(
      "analysis_set", "analysis_role", "source_marker_id", "resolved_rsid",
      "marker_chr", "marker_pos", "marker_allele1", "marker_allele2",
      "exposure_effect_allele", "exposure_other_allele", "beta", "se", "pval",
      "eaf", "n_samples", "F_stat", "reference_chr_grch37", "reference_pos_grch37",
      "identity_resolution_status", "chen_marker_effect_allele_compatibility", "apoe_region"
    )]
  }
  inst_included <- make_instruments(out[["included_clumps"]], included_input, "APOE_included")
  inst_excluded <- make_instruments(out[["excluded_clumps"]], excluded_input, "APOE_excluded")
  stop_if(any(inst_included$F_stat < 10) || any(inst_excluded$F_stat < 10), "F<10 detected in final Chen instruments; manual review required.")

  write_pair <- function(df, parquet_path, tsv_path) {
    atomic(parquet_path, function(p) {
      DBI::dbWriteTable(con, "write_df", df, overwrite = TRUE)
      DBI::dbExecute(con, sprintf("COPY write_df TO '%s' (FORMAT PARQUET)", qpath(p)))
      DBI::dbExecute(con, "DROP TABLE write_df")
    })
    atomic(tsv_path, function(p) write_tsv_precise(df, p))
  }
  write_pair(inst_included, out[["included_inst_parquet"]], out[["included_inst_tsv"]])
  write_pair(inst_excluded, out[["excluded_inst_parquet"]], out[["excluded_inst_tsv"]])
  rt_included <- validate_roundtrip(con, out[["included_inst_parquet"]], out[["included_inst_tsv"]], inst_included)
  rt_excluded <- validate_roundtrip(con, out[["excluded_inst_parquet"]], out[["excluded_inst_tsv"]], inst_excluded)
  stop_if(!all(unlist(rt_included[c("same_cols", "same_n", "same_rsid", "char_ok", "num_ok")])), "Included Parquet/TSV validation failed.")
  stop_if(!all(unlist(rt_excluded[c("same_cols", "same_n", "same_rsid", "char_ok", "num_ok")])), "Excluded Parquet/TSV validation failed.")

  strength <- function(x) list(n = nrow(x), F_min = min(x$F_stat), F_mean = mean(x$F_stat), F_median = stats::median(x$F_stat), F_max = max(x$F_stat), F_lt10_count = sum(x$F_stat < 10))
  shared <- intersect(inst_included$resolved_rsid, inst_excluded$resolved_rsid)
  included_only <- setdiff(inst_included$resolved_rsid, inst_excluded$resolved_rsid)
  excluded_only <- setdiff(inst_excluded$resolved_rsid, inst_included$resolved_rsid)
  vuckovic_overlap <- data.frame(
    metric = c("vuckovic_included_n", "chen_included_n", "overlap_count", "chen_only_count", "vuckovic_only_count"),
    value = NA_integer_,
    stringsAsFactors = FALSE
  )
  if (file.exists(vuckovic_inst)) {
    v <- read.delim(vuckovic_inst, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    vuckovic_overlap$value <- c(length(unique(v$rsid)), length(unique(inst_included$resolved_rsid)), length(intersect(unique(v$rsid), unique(inst_included$resolved_rsid))), length(setdiff(unique(inst_included$resolved_rsid), unique(v$rsid))), length(setdiff(unique(v$rsid), unique(inst_included$resolved_rsid))))
  }
  atomic(out[["overlap"]], function(p) write_csv_precise(vuckovic_overlap, p))

  decision82_resolution <- read.csv(decision82_resolution_path, stringsAsFactors = FALSE)
  decision82_format <- read.csv(decision82_format_path, stringsAsFactors = FALSE)
  d82_p_sig <- decision82_resolution[decision82_resolution$subset == "p_lt_5e-8_candidates", , drop = FALSE]
  d82_unique <- d82_p_sig$value[d82_p_sig$metric == "n_coordinate_allele_unique_reference"]
  d82_raw <- decision82_format$n[decision82_format$subset == "p_lt_5e-8_candidates" & decision82_format$identifier_grammar == "chr_pos_allele1_allele2"]
  d82_compare <- list(
    decision82_raw_candidate_count = if (length(d82_raw)) as.integer(d82_raw[[1]]) else NA_integer_,
    decision82_unique_coordinate_plus_allele_count = if (length(d82_unique)) as.integer(d82_unique[[1]]) else NA_integer_,
    v2_raw_candidate_count = as.integer(raw_candidate_count),
    v2_unique_exact_identity_count = as.integer(unique_exact_identity_count),
    comparison_role = "informational_not_pass_gate"
  )

  counts <- data.frame(
    metric = c(
      "raw_candidate_count", "pre_resolution_eligible_count", "unique_exact_identity_count",
      "resolved_duplicate_rsid_count", "resolved_eligible_count", "apoe_window_candidate_count",
      "included_clump_input_count", "excluded_clump_input_count", "included_nsnp", "excluded_nsnp",
      "shared_nsnp", "included_only_nsnp", "excluded_only_nsnp",
      "included_F_lt10", "excluded_F_lt10"
    ),
    value = c(
      as.integer(raw_candidate_count), as.integer(pre_resolution_eligible_count), as.integer(unique_exact_identity_count),
      as.integer(resolved_duplicate_rsid_count), as.integer(resolved_eligible_count), as.integer(apoe_window_candidate_count),
      nrow(included_input), nrow(excluded_input), nrow(inst_included), nrow(inst_excluded),
      length(shared), length(included_only), length(excluded_only),
      sum(inst_included$F_stat < 10), sum(inst_excluded$F_stat < 10)
    ),
    stringsAsFactors = FALSE
  )
  atomic(out[["counts"]], function(p) write_csv_precise(counts, p))

  sha_after <- toupper(hash_file(source_path))
  renv_after <- hash_file(renv_lock)
  hard_checks <- list(
    decision_78_source_gate = identical(cert$certification_status, "passed") && length(cert$hard_check_failures) == 0L,
    decision_79_contract_gate = identical(contract$contract_status, "frozen") && length(contract$hard_check_failures) == 0L,
    decision_80_coordinate_gate = identical(coord$coordinate_authority_status, "passed") && length(coord$hard_check_failures) == 0L,
    decision_81_failure_preserved = file.exists(decision81_path),
    decision_82_semantics_gate = file.exists(decision82_format_path) && file.exists(decision82_resolution_path),
    decision_83_official_dictionary_gate = identical(official$audit_status, "passed") && length(official$hard_check_failures) == 0L,
    identifier_resolution_amendment_gate = identical(amendment$amendment_status, "frozen") && isTRUE(amendment$approved_for_chen_forward_instrument_selection_v2),
    source_sha_before_after = identical(sha_before, expected_sha) && identical(sha_after, expected_sha) && identical(sha_before, sha_after),
    p_threshold_preserved = TRUE,
    marker_parser_valid = as.integer(pre_resolution_eligible_count) > 0L,
    effect_orientation_separate = TRUE,
    exact_coordinate_allele_identity_only = TRUE,
    no_coordinate_only_resolution = TRUE,
    no_strand_complement_identity_rescue = TRUE,
    no_proxy = TRUE,
    no_liftover = TRUE,
    no_nearest_variant = TRUE,
    resolved_rsid_uniqueness = as.integer(resolved_duplicate_rsid_count) == 0L && !anyDuplicated(resolved_eligible$resolved_rsid),
    chen_marker_effect_allele_compatibility = all(resolved_eligible$chen_marker_effect_allele_compatibility == "same_unordered_set"),
    apoe_exclusion_before_independent_clumping = all(!excluded_input$apoe_region),
    included_independent_clumping_completed = file.exists(out[["included_clumps"]]),
    excluded_independent_clumping_completed = file.exists(out[["excluded_clumps"]]),
    plink_binary_verified = identical(tolower(plink_sha), "247491bfca7512e070dc99d6565e9fc56f3a52ad5afc01286016271d34c4992f") && grepl("PLINK v2.0.0-a.7.1", plink_version, fixed = TRUE),
    variant_level_n_preserved = all(is.finite(inst_included$n_samples) & inst_included$n_samples > 0) && all(is.finite(inst_excluded$n_samples) & inst_excluded$n_samples > 0),
    instrument_strength_completed = TRUE,
    no_f_ge_30_filter = TRUE,
    included_parquet_tsv_consistency = all(unlist(rt_included[c("same_cols", "same_n", "same_rsid", "char_ok", "num_ok")])),
    excluded_parquet_tsv_consistency = all(unlist(rt_excluded[c("same_cols", "same_n", "same_rsid", "char_ok", "num_ok")])),
    independent_replication_false = identical(contract$independent_replication, FALSE),
    no_outcome_extraction = TRUE,
    no_harmonisation = TRUE,
    no_mr = TRUE,
    no_steiger = TRUE
  )
  failures <- names(hard_checks)[!unlist(hard_checks)]
  status <- if (
    length(failures) == 0L &&
      as.integer(resolved_eligible_count) > 0L &&
      as.integer(resolved_duplicate_rsid_count) == 0L &&
      nrow(inst_included) > 0L &&
      nrow(inst_excluded) > 0L &&
      sum(inst_included$F_stat < 10) == 0L &&
      sum(inst_excluded$F_stat < 10) == 0L
  ) "passed" else "failed"
  approved_outcome <- identical(status, "passed")
  qc <- list(
    instrument_selection_version = "v2",
    decision = 85,
    analysis_direction = "Hb_to_delirium",
    analysis_role = "forward_alternative_hb_gwas_sensitivity",
    exposure_source = "Chen_2020_Hb_BCX2",
    independent_replication = FALSE,
    source_sha256_before = sha_before,
    source_sha256_after = sha_after,
    source_sha256_certified = expected_sha,
    identifier_resolution_amendment_version = "v1",
    identifier_resolution_method = "exact_GRCh37_coordinate_plus_unordered_allele_set_unique_reference_match",
    p_threshold = 5e-8,
    ld_r2 = 0.001,
    ld_window_kb = 10000,
    raw_candidate_count = as.integer(raw_candidate_count),
    pre_resolution_eligible_count = as.integer(pre_resolution_eligible_count),
    unique_exact_identity_count = as.integer(unique_exact_identity_count),
    identity_resolution_failure_counts = records(resolution_breakdown),
    chen_marker_effect_allele_compatibility_counts = records(compatibility_breakdown),
    resolved_duplicate_rsid_count = as.integer(resolved_duplicate_rsid_count),
    resolved_eligible_count = as.integer(resolved_eligible_count),
    apoe_window_candidate_count = as.integer(apoe_window_candidate_count),
    included_clump_input_count = nrow(included_input),
    excluded_clump_input_count = nrow(excluded_input),
    included_nsnp = nrow(inst_included),
    excluded_nsnp = nrow(inst_excluded),
    shared_nsnp = length(shared),
    included_only_nsnp = length(included_only),
    excluded_only_nsnp = length(excluded_only),
    instrument_strength_included = strength(inst_included),
    instrument_strength_excluded = strength(inst_excluded),
    variant_level_n_preserved = TRUE,
    proxy_used = FALSE,
    liftover_used = FALSE,
    nearest_variant_used = FALSE,
    strand_complement_identity_rescue_used = FALSE,
    outcome_extraction_performed = FALSE,
    harmonisation_performed = FALSE,
    mr_run = FALSE,
    steiger_run = FALSE,
    parquet_tsv_consistency = list(included = rt_included, excluded = rt_excluded),
    decision82_resolution_comparison = d82_compare,
    vuckovic_overlap_audit = records(vuckovic_overlap),
    plink = list(
      path = "tools/plink2/plink2.exe",
      sha256 = plink_sha,
      version = plink_version,
      included_command_recorded_in_log = TRUE,
      excluded_command_recorded_in_log = TRUE,
      included_missing_id = count_missing(out[["included_missing_id"]]),
      excluded_missing_id = count_missing(out[["excluded_missing_id"]]),
      included_missing_allele = count_missing(out[["included_missing_allele"]]),
      excluded_missing_allele = count_missing(out[["excluded_missing_allele"]])
    ),
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after,
    renv_lock_unchanged = identical(renv_before, renv_after),
    instrument_selection_status = status,
    approved_for_chen_forward_outcome_extraction = approved_outcome,
    hard_checks = hard_checks,
    hard_check_failures = failures,
    informational_findings = list(
      chen_vs_vuckovic_independent_replication = FALSE,
      decision82_json_not_parsed = TRUE,
      decision82_csv_used_for_informational_comparison = TRUE,
      excluded_clumping_not_derived_by_posthoc_deletion_from_included_clumps = TRUE
    )
  )
  atomic(out[["qc_json"]], function(p) jsonlite::write_json(qc, p, pretty = TRUE, auto_unbox = TRUE, na = "null"))
  write_decision85(qc)
  log_line("instrument_selection_status=", status, " approved_for_chen_forward_outcome_extraction=", approved_outcome, " included_nsnp=", nrow(inst_included), " excluded_nsnp=", nrow(inst_excluded), " hard_check_failures=", paste(failures, collapse = ";"))
  stop_if(!identical(status, "passed"), "Chen forward instrument selection V2 failed; QC JSON retained.")
}

tryCatch(main(), error = function(e) {
  log_line("status=failed error=", conditionMessage(e))
  quit(status = 1L)
})
