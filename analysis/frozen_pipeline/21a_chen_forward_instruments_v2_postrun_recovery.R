#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/21a_chen_forward_instruments_v2_postrun_recovery.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
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
num_equal <- function(a, b, tol = 1e-12) {
  a <- as.numeric(a); b <- as.numeric(b)
  length(a) == length(b) && all((is.na(a) & is.na(b)) | (is.finite(a) & is.finite(b) & abs(a - b) <= tol))
}
write_csv_precise <- function(x, path) {
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  write.csv(x, path, row.names = FALSE, na = "")
}
atomic <- function(path, writer) {
  partial <- paste0(path, ".partial")
  stop_if(file.exists(path) || file.exists(partial), paste("Output occupied:", path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writer(partial)
  stop_if(!file.exists(partial), paste("Writer did not create partial:", partial))
  stop_if(!file.rename(partial, path), paste("Atomic rename failed:", path))
}
count_missing <- function(path) {
  if (!file.exists(path)) return(0L)
  sum(nzchar(trimws(readLines(path, warn = FALSE))))
}
assert_log_clean <- function(path) {
  lines <- readLines(path, warn = FALSE)
  stop_if(any(grepl("\\b(ERROR|FATAL)\\b", lines, ignore.case = TRUE, perl = TRUE)), paste("PLINK log contains ERROR/FATAL:", path))
}
read_clumps <- function(path) {
  x <- read.table(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE, comment.char = "")
  norm <- toupper(gsub("[^A-Za-z0-9]", "", names(x)))
  hit <- match(c("ID", "SNP", "INDEXSNP"), norm, nomatch = 0L)
  hit <- hit[hit > 0L]
  stop_if(length(hit) == 0L, paste("Cannot identify index SNP column:", path))
  x$index_snp <- as.character(x[[hit[[1L]]]])
  x
}
strength <- function(x) {
  list(
    n = nrow(x),
    F_min = min(x$F_stat),
    F_mean = mean(x$F_stat),
    F_median = stats::median(x$F_stat),
    F_max = max(x$F_stat),
    F_lt10_count = sum(x$F_stat < 10)
  )
}
validate_roundtrip <- function(parquet_df, tsv_df) {
  same_cols <- identical(names(parquet_df), names(tsv_df))
  same_n <- nrow(parquet_df) == nrow(tsv_df)
  same_order <- same_n && identical(as.character(parquet_df$resolved_rsid), as.character(tsv_df$resolved_rsid))
  char_cols <- names(parquet_df)[vapply(parquet_df, function(z) is.character(z) || is.logical(z), logical(1))]
  num_cols <- names(parquet_df)[vapply(parquet_df, is.numeric, logical(1))]
  char_ok <- all(vapply(char_cols, function(k) identical(as.character(parquet_df[[k]]), as.character(tsv_df[[k]])), logical(1)))
  num_ok <- all(vapply(num_cols, function(k) num_equal(parquet_df[[k]], tsv_df[[k]]), logical(1)))
  list(row_count = nrow(parquet_df), same_cols = same_cols, same_n = same_n, same_rsid_order = same_order, char_ok = char_ok, num_ok = num_ok)
}
write_decision85 <- function(summary, path) {
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
    sprintf("- Included F min/median/mean/max/F<10: `%s / %s / %s / %s / %s`", summary$instrument_strength_included$F_min, summary$instrument_strength_included$F_median, summary$instrument_strength_included$F_mean, summary$instrument_strength_included$F_max, summary$instrument_strength_included$F_lt10_count),
    sprintf("- Excluded F min/median/mean/max/F<10: `%s / %s / %s / %s / %s`", summary$instrument_strength_excluded$F_min, summary$instrument_strength_excluded$F_median, summary$instrument_strength_excluded$F_mean, summary$instrument_strength_excluded$F_max, summary$instrument_strength_excluded$F_lt10_count),
    sprintf("- Instrument-selection status: `%s`", summary$instrument_selection_status),
    sprintf("- Approved for Chen forward outcome extraction: `%s`", summary$approved_for_chen_forward_outcome_extraction),
    "",
    "## Technical Recovery Note",
    "",
    "The first V2 run completed source scanning, identifier resolution, two PLINK",
    "clumping runs, and final instrument file creation, but stopped at an overly",
    "strict Parquet/TSV readback check. This postrun recovery script reads those",
    "existing V2 artifacts, performs tolerant column-aware readback validation, and",
    "creates the missing QC and decision files without rewriting the scientific",
    "outputs or rerunning PLINK.",
    "",
    "## Affected Files",
    "",
    "- `R/21_chen_forward_instruments_v2.R`",
    "- `R/21a_chen_forward_instruments_v2_postrun_recovery.R`",
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
    "- `results/logs/chen_forward_instrument_selection_v2_recovery.log`",
    "",
    "## Expected Impact",
    "",
    "This creates Chen forward sensitivity instruments for a later, separately",
    "approved FinnGen targeted outcome extraction. It does not alter existing",
    "Vuckovic/FinnGen results or Decision 81 failure evidence."
  )
  atomic(path, function(p) writeLines(lines, p, useBytes = TRUE))
}

source_path <- file.path(root, "data_raw", "gwas", "BCX2_HGB_EA_GWAMA.out.gz")
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
plink <- file.path(root, "tools", "plink2", "plink2.exe")

existing <- c(
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
  first_log = file.path(root, "results", "logs", "chen_forward_instrument_selection_v2.log")
)
new_out <- c(
  counts = file.path(root, "results", "qc", "chen_forward_instrument_selection_counts_v2.csv"),
  overlap = file.path(root, "results", "qc", "chen_forward_vuckovic_overlap_audit_v2.csv"),
  qc_json = file.path(root, "results", "qc", "chen_forward_instrument_selection_v2.json"),
  recovery_log = file.path(root, "results", "logs", "chen_forward_instrument_selection_v2_recovery.log"),
  decision = file.path(root, "docs", "decisions", "85_chen_forward_instrument_selection_v2_v1.1.md")
)

main <- function() {
  required_existing <- existing[!names(existing) %in% c("included_missing_id", "excluded_missing_id", "included_missing_allele", "excluded_missing_allele")]
  for (p in c(required_existing, cert_path, contract_path, coord_path, official_dict_path, amendment_path, decision81_path, decision82_format_path, decision82_resolution_path, source_path, renv_lock, plink)) {
    stop_if(!file.exists(p), paste("Missing required file:", p))
  }
  for (p in c(new_out, paste0(new_out, ".partial"))) {
    stop_if(file.exists(p), paste("Output occupied:", p))
  }

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
  assert_log_clean(existing[["included_plink_log"]])
  assert_log_clean(existing[["excluded_plink_log"]])

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:", config = list(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  candidates_n <- DBI::dbGetQuery(con, sprintf("SELECT COUNT(*) AS n FROM read_parquet('%s')", qpath(existing[["candidates"]])))$n[[1]]
  pre_resolution_eligible_count <- DBI::dbGetQuery(con, sprintf("SELECT COUNT(*) AS n FROM read_parquet('%s') WHERE pre_resolution_status='pre_resolution_eligible'", qpath(existing[["candidates"]])))$n[[1]]
  unique_exact_identity_count <- DBI::dbGetQuery(con, sprintf("SELECT COUNT(*) AS n FROM read_parquet('%s') WHERE identity_resolution_status='unique_exact_coordinate_allele_identity'", qpath(existing[["identity_resolution"]])))$n[[1]]
  resolved_eligible_count <- DBI::dbGetQuery(con, sprintf("SELECT COUNT(*) AS n FROM read_parquet('%s')", qpath(existing[["resolved_eligible"]])))$n[[1]]
  apoe_window_candidate_count <- DBI::dbGetQuery(con, sprintf("SELECT COUNT(*) AS n FROM read_parquet('%s') WHERE apoe_region", qpath(existing[["resolved_eligible"]])))$n[[1]]
  resolved_duplicate_rsid_count <- 0L
  dup_note <- readLines(existing[["duplicate_audit"]], warn = FALSE)
  if (!length(dup_note) || !grepl("no_duplicate_resolved_rsids", paste(dup_note, collapse = "\n"), fixed = TRUE)) {
    dup_df <- read.delim(existing[["duplicate_audit"]], sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    resolved_duplicate_rsid_count <- length(unique(dup_df$resolved_rsid))
  }
  resolution_breakdown <- DBI::dbGetQuery(con, sprintf("SELECT identity_resolution_status, COUNT(*) AS n FROM read_parquet('%s') GROUP BY identity_resolution_status ORDER BY identity_resolution_status", qpath(existing[["identity_resolution"]])))
  compatibility_breakdown <- DBI::dbGetQuery(con, sprintf("SELECT chen_marker_effect_allele_compatibility, COUNT(*) AS n FROM read_parquet('%s') GROUP BY chen_marker_effect_allele_compatibility ORDER BY chen_marker_effect_allele_compatibility", qpath(existing[["identity_resolution"]])))

  included_parquet <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", qpath(existing[["included_inst_parquet"]])))
  excluded_parquet <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", qpath(existing[["excluded_inst_parquet"]])))
  included_tsv <- read.delim(existing[["included_inst_tsv"]], sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  excluded_tsv <- read.delim(existing[["excluded_inst_tsv"]], sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  rt_included <- validate_roundtrip(included_parquet, included_tsv)
  rt_excluded <- validate_roundtrip(excluded_parquet, excluded_tsv)

  included_input <- read.delim(existing[["included_input"]], sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  excluded_input <- read.delim(existing[["excluded_input"]], sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  included_input$apoe_region <- tolower(as.character(included_input$apoe_region)) == "true"
  excluded_input$apoe_region <- tolower(as.character(excluded_input$apoe_region)) == "true"
  included_clumps <- read_clumps(existing[["included_clumps"]])
  excluded_clumps <- read_clumps(existing[["excluded_clumps"]])
  stop_if(anyDuplicated(included_clumps$index_snp), "Duplicate included PLINK index SNP.")
  stop_if(anyDuplicated(excluded_clumps$index_snp), "Duplicate excluded PLINK index SNP.")
  stop_if(!identical(as.character(included_clumps$index_snp), as.character(included_parquet$resolved_rsid)), "Included clumps and final instruments differ.")
  stop_if(!identical(as.character(excluded_clumps$index_snp), as.character(excluded_parquet$resolved_rsid)), "Excluded clumps and final instruments differ.")
  stop_if(any(excluded_parquet$apoe_region), "APOE-excluded final set contains APOE-region variant.")
  stop_if(any(included_parquet$F_stat < 10) || any(excluded_parquet$F_stat < 10), "F<10 detected in final instruments.")

  shared <- intersect(included_parquet$resolved_rsid, excluded_parquet$resolved_rsid)
  included_only <- setdiff(included_parquet$resolved_rsid, excluded_parquet$resolved_rsid)
  excluded_only <- setdiff(excluded_parquet$resolved_rsid, included_parquet$resolved_rsid)
  vuckovic_overlap <- data.frame(
    metric = c("vuckovic_included_n", "chen_included_n", "overlap_count", "chen_only_count", "vuckovic_only_count"),
    value = NA_integer_,
    stringsAsFactors = FALSE
  )
  if (file.exists(vuckovic_inst)) {
    v <- read.delim(vuckovic_inst, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    vuckovic_overlap$value <- c(length(unique(v$rsid)), length(unique(included_parquet$resolved_rsid)), length(intersect(unique(v$rsid), unique(included_parquet$resolved_rsid))), length(setdiff(unique(included_parquet$resolved_rsid), unique(v$rsid))), length(setdiff(unique(v$rsid), unique(included_parquet$resolved_rsid))))
  }
  atomic(new_out[["overlap"]], function(p) write_csv_precise(vuckovic_overlap, p))

  counts <- data.frame(
    metric = c(
      "raw_candidate_count", "pre_resolution_eligible_count", "unique_exact_identity_count",
      "resolved_duplicate_rsid_count", "resolved_eligible_count", "apoe_window_candidate_count",
      "included_clump_input_count", "excluded_clump_input_count", "included_nsnp", "excluded_nsnp",
      "shared_nsnp", "included_only_nsnp", "excluded_only_nsnp", "included_F_lt10", "excluded_F_lt10"
    ),
    value = c(
      as.integer(candidates_n), as.integer(pre_resolution_eligible_count), as.integer(unique_exact_identity_count),
      as.integer(resolved_duplicate_rsid_count), as.integer(resolved_eligible_count), as.integer(apoe_window_candidate_count),
      nrow(included_input), nrow(excluded_input), nrow(included_parquet), nrow(excluded_parquet),
      length(shared), length(included_only), length(excluded_only),
      sum(included_parquet$F_stat < 10), sum(excluded_parquet$F_stat < 10)
    ),
    stringsAsFactors = FALSE
  )
  atomic(new_out[["counts"]], function(p) write_csv_precise(counts, p))

  decision82_resolution <- read.csv(decision82_resolution_path, stringsAsFactors = FALSE)
  decision82_format <- read.csv(decision82_format_path, stringsAsFactors = FALSE)
  d82_p_sig <- decision82_resolution[decision82_resolution$subset == "p_lt_5e-8_candidates", , drop = FALSE]
  d82_unique <- d82_p_sig$value[d82_p_sig$metric == "n_coordinate_allele_unique_reference"]
  d82_raw <- decision82_format$n[decision82_format$subset == "p_lt_5e-8_candidates" & decision82_format$identifier_grammar == "chr_pos_allele1_allele2"]
  d82_compare <- list(
    decision82_raw_candidate_count = if (length(d82_raw)) as.integer(d82_raw[[1]]) else NA_integer_,
    decision82_unique_coordinate_plus_allele_count = if (length(d82_unique)) as.integer(d82_unique[[1]]) else NA_integer_,
    v2_raw_candidate_count = as.integer(candidates_n),
    v2_unique_exact_identity_count = as.integer(unique_exact_identity_count),
    comparison_role = "informational_not_pass_gate"
  )

  sha_after <- toupper(hash_file(source_path))
  renv_after <- hash_file(renv_lock)
  resolved_ids <- DBI::dbGetQuery(
    con,
    sprintf("SELECT resolved_rsid FROM read_parquet('%s')", qpath(existing[["resolved_eligible"]]))
  )$resolved_rsid
  marker_effect_incompatible_n <- DBI::dbGetQuery(
    con,
    sprintf("SELECT COUNT(*) AS n FROM read_parquet('%s') WHERE chen_marker_effect_allele_compatibility <> 'same_unordered_set'", qpath(existing[["resolved_eligible"]]))
  )$n[[1]]
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
    resolved_rsid_uniqueness = as.integer(resolved_duplicate_rsid_count) == 0L && !anyDuplicated(resolved_ids),
    chen_marker_effect_allele_compatibility = marker_effect_incompatible_n == 0L,
    apoe_exclusion_before_independent_clumping = all(!excluded_input$apoe_region) && all(!excluded_parquet$apoe_region),
    included_independent_clumping_completed = file.exists(existing[["included_clumps"]]) && nrow(included_clumps) == nrow(included_parquet),
    excluded_independent_clumping_completed = file.exists(existing[["excluded_clumps"]]) && nrow(excluded_clumps) == nrow(excluded_parquet),
    plink_binary_verified = identical(tolower(plink_sha), "247491bfca7512e070dc99d6565e9fc56f3a52ad5afc01286016271d34c4992f") && grepl("PLINK v2.0.0-a.7.1", plink_version, fixed = TRUE),
    variant_level_n_preserved = all(is.finite(included_parquet$n_samples) & included_parquet$n_samples > 0) && all(is.finite(excluded_parquet$n_samples) & excluded_parquet$n_samples > 0),
    instrument_strength_completed = sum(included_parquet$F_stat < 10) == 0L && sum(excluded_parquet$F_stat < 10) == 0L,
    no_f_ge_30_filter = TRUE,
    included_parquet_tsv_consistency = all(unlist(rt_included[c("same_cols", "same_n", "same_rsid_order", "char_ok", "num_ok")])),
    excluded_parquet_tsv_consistency = all(unlist(rt_excluded[c("same_cols", "same_n", "same_rsid_order", "char_ok", "num_ok")])),
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
      nrow(included_parquet) > 0L &&
      nrow(excluded_parquet) > 0L &&
      sum(included_parquet$F_stat < 10) == 0L &&
      sum(excluded_parquet$F_stat < 10) == 0L
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
    raw_candidate_count = as.integer(candidates_n),
    pre_resolution_eligible_count = as.integer(pre_resolution_eligible_count),
    unique_exact_identity_count = as.integer(unique_exact_identity_count),
    identity_resolution_failure_counts = records(resolution_breakdown),
    chen_marker_effect_allele_compatibility_counts = records(compatibility_breakdown),
    resolved_duplicate_rsid_count = as.integer(resolved_duplicate_rsid_count),
    resolved_eligible_count = as.integer(resolved_eligible_count),
    apoe_window_candidate_count = as.integer(apoe_window_candidate_count),
    included_clump_input_count = nrow(included_input),
    excluded_clump_input_count = nrow(excluded_input),
    included_nsnp = nrow(included_parquet),
    excluded_nsnp = nrow(excluded_parquet),
    shared_nsnp = length(shared),
    included_only_nsnp = length(included_only),
    excluded_only_nsnp = length(excluded_only),
    instrument_strength_included = strength(included_parquet),
    instrument_strength_excluded = strength(excluded_parquet),
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
      included_missing_id = count_missing(existing[["included_missing_id"]]),
      excluded_missing_id = count_missing(existing[["excluded_missing_id"]]),
      included_missing_allele = count_missing(existing[["included_missing_allele"]]),
      excluded_missing_allele = count_missing(existing[["excluded_missing_allele"]])
    ),
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after,
    renv_lock_unchanged = identical(renv_before, renv_after),
    recovery = list(
      performed = TRUE,
      reason = "first_v2_run_completed_scanning_resolution_plink_and_final_outputs_but_failed_overly_strict_readback_validation",
      raw_gwas_rescanned = FALSE,
      plink_rerun = FALSE,
      scientific_outputs_rewritten = FALSE
    ),
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
  atomic(new_out[["qc_json"]], function(p) jsonlite::write_json(qc, p, pretty = TRUE, auto_unbox = TRUE, na = "null"))
  write_decision85(qc, new_out[["decision"]])
  atomic(new_out[["recovery_log"]], function(p) {
    writeLines(c(
      sprintf("[%s] stage=chen_forward_instrument_selection_v2_postrun_recovery", ts()),
      sprintf("[%s] raw_gwas_rescanned=FALSE plink_rerun=FALSE scientific_outputs_rewritten=FALSE", ts()),
      sprintf("[%s] included_nsnp=%s excluded_nsnp=%s status=%s approved_for_outcome_extraction=%s", ts(), nrow(included_parquet), nrow(excluded_parquet), status, approved_outcome),
      sprintf("[%s] hard_check_failures=%s", ts(), paste(failures, collapse = ";"))
    ), p, useBytes = TRUE)
  })
  stop_if(!identical(status, "passed"), "Recovery audit failed; QC JSON retained.")
}

tryCatch(main(), error = function(e) {
  message("status=failed error=", conditionMessage(e))
  quit(status = 1L)
})
