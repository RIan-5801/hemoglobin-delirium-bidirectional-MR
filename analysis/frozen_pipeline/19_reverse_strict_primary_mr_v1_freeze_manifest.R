#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root") {
  stop("Usage: Rscript R/19_reverse_strict_primary_mr_v1_freeze_manifest.R --project-root E:/Research/hb_delirium_bidir_mr", call. = FALSE)
}
root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
source(file.path(root, "renv", "activate.R"))

for (pkg in c("digest", "jsonlite")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
ts <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))

out <- c(
  manifest = file.path(root, "results", "qc", "reverse_strict_primary_mr_v1_freeze_manifest.csv"),
  json = file.path(root, "results", "qc", "reverse_strict_primary_mr_v1_freeze.json"),
  log = file.path(root, "results", "logs", "reverse_strict_primary_mr_v1_freeze.log")
)
stop_if(any(file.exists(c(out, paste0(out, ".partial")))), "A strict MR freeze final or partial target exists; refusing to overwrite.")
dir.create(dirname(out[["manifest"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out[["log"]]), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(..., collapse = "")), file = out[["log"]], append = TRUE)

main <- function() {
  log_line("stage=reverse_strict_primary_mr_v1_freeze")
  rel <- c(
    "docs/decisions/73_reverse_strict_primary_mr_analysis_contract_v1_v1.1.md",
    "docs/decisions/74_reverse_strict_primary_mr_v1_v1.1.md",
    "docs/decisions/75_reverse_strict_primary_mr_v1_freeze_v1.1.md",
    "R/19_reverse_strict_primary_mr_v1.R",
    "R/19_reverse_strict_primary_mr_v1_freeze_manifest.R",
    "results/qc/reverse_strict_primary_mr_analysis_contract_v1.json",
    "results/tables/reverse_strict_primary_mr_estimates_v1.csv",
    "results/tables/reverse_strict_primary_mr_estimates_doubling_odds_v1.csv",
    "results/qc/reverse_strict_primary_mr_v1.json",
    "results/logs/reverse_strict_primary_mr_v1.log",
    "docs/decisions/56_vuckovic_hb_reverse_primary_formal_harmonisation_v4_freeze_v1.1.md",
    "results/qc/vuckovic_hb_reverse_primary_formal_harmonisation_v4_freeze.json",
    "results/qc/vuckovic_hb_reverse_primary_formal_harmonisation_v4_freeze_manifest.csv",
    "renv.lock"
  )
  roles <- c(
    "decision_73_contract", "decision_74_strict_mr", "decision_75_freeze", "strict_mr_script",
    "strict_mr_freeze_script", "strict_mr_contract_json", "strict_raw_estimates",
    "strict_doubling_estimates", "strict_mr_qc", "strict_mr_log",
    "strict_harmonisation_freeze_decision", "strict_harmonisation_freeze_json",
    "strict_harmonisation_freeze_manifest", "renv_lock"
  )
  paths <- file.path(root, rel)
  stop_if(any(!file.exists(paths)), paste("Missing strict freeze input(s):", paste(rel[!file.exists(paths)], collapse = "; ")))
  strict_qc <- jsonlite::fromJSON(file.path(root, "results", "qc", "reverse_strict_primary_mr_v1.json"), simplifyVector = FALSE)
  contract <- jsonlite::fromJSON(file.path(root, "results", "qc", "reverse_strict_primary_mr_analysis_contract_v1.json"), simplifyVector = FALSE)
  freeze_manifest_sha <- hash_file(file.path(root, "results", "qc", "vuckovic_hb_reverse_primary_formal_harmonisation_v4_freeze_manifest.csv"))
  raw <- read.csv(file.path(root, "results", "tables", "reverse_strict_primary_mr_estimates_v1.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  doubling <- read.csv(file.path(root, "results", "tables", "reverse_strict_primary_mr_estimates_doubling_odds_v1.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  hard_checks <- list(
    strict_mr_status_passed = identical(strict_qc$strict_primary_mr_status, "passed"),
    strict_mr_approved_for_freeze = isTRUE(strict_qc$approved_for_reverse_strict_primary_results_freeze),
    strict_mr_hard_failures_empty = length(strict_qc$hard_check_failures) == 0L,
    strict_contract_frozen = identical(contract$contract_status, "frozen"),
    strict_contract_approved = isTRUE(contract$approved_for_reverse_strict_primary_mr_execution),
    strict_harmonisation_freeze_manifest_sha_verified = identical(freeze_manifest_sha, strict_qc$input_sha256$freeze_manifest),
    included_nsnp_equals_one = as.integer(strict_qc$included_nsnp) == 1L,
    excluded_nsnp_equals_one = as.integer(strict_qc$excluded_nsnp) == 1L,
    method_wald_ratio_only = all(raw$method == "Wald ratio"),
    no_ivw_single_snp_mislabel = !any(grepl("IVW|Inverse variance", raw$method, ignore.case = TRUE)),
    raw_and_qc_rsids_match = identical(as.character(raw$rsid[raw$analysis_set == "APOE_included"]), as.character(strict_qc$included_rsid)) &&
      identical(as.character(raw$rsid[raw$analysis_set == "APOE_excluded"]), as.character(strict_qc$excluded_rsid)),
    doubling_rescaling_correct = all(abs(doubling$beta_per_doubling_odds - raw$beta * log(2)) <= 1e-12) &&
      all(abs(doubling$se_per_doubling_odds - raw$se * log(2)) <= 1e-12) &&
      identical(as.numeric(doubling$pval), as.numeric(raw$pval)),
    strict_relaxed_hierarchy_preserved = TRUE,
    no_relaxed_mr_run_by_freeze = TRUE,
    no_reharmonisation = TRUE,
    no_proxy = TRUE,
    no_liftover = TRUE,
    no_rescan_or_clumping = TRUE
  )
  failures <- names(hard_checks)[!unlist(hard_checks)]
  status <- if (length(failures) == 0L) "passed" else "failed"
  manifest <- data.frame(
    file_role = roles,
    relative_path = rel,
    file_size_bytes = as.numeric(file.info(paths)$size),
    sha256 = vapply(paths, hash_file, character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(manifest, paste0(out[["manifest"]], ".partial"), row.names = FALSE)
  manifest_sha <- hash_file(paste0(out[["manifest"]], ".partial"))
  result <- list(
    freeze_version = "v1",
    decision = 75,
    authoritative_reverse_strict_primary_mr_version = "v1",
    analysis_direction = "delirium_to_Hb",
    analysis_role = "secondary_reverse_primary",
    instrument_threshold = 5e-8,
    included_nsnp = strict_qc$included_nsnp,
    excluded_nsnp = strict_qc$excluded_nsnp,
    included_rsid = strict_qc$included_rsid,
    excluded_rsid = strict_qc$excluded_rsid,
    authoritative_method = "Wald ratio",
    included_raw_result = records(raw[raw$analysis_set == "APOE_included", , drop = FALSE])[[1]],
    excluded_raw_result = records(raw[raw$analysis_set == "APOE_excluded", , drop = FALSE])[[1]],
    included_doubling_odds_result = records(doubling[doubling$analysis_set == "APOE_included", , drop = FALSE])[[1]],
    excluded_doubling_odds_result = records(doubling[doubling$analysis_set == "APOE_excluded", , drop = FALSE])[[1]],
    diagnostics_not_estimable = strict_qc$diagnostics_not_estimable,
    evidence_hierarchy = list(
      strict_primary = "P<5e-8 secondary_reverse_primary",
      relaxed_branch = "P<5e-6 protocol_triggered_exploratory_fallback",
      relaxed_may_not_override_strict_primary = TRUE
    ),
    manifest_sha256 = manifest_sha,
    freeze_status = status,
    approved_for_reverse_relaxed_mr_execution = identical(status, "passed"),
    hard_checks = hard_checks,
    hard_check_failures = failures
  )
  jsonlite::write_json(result, paste0(out[["json"]], ".partial"), pretty = TRUE, auto_unbox = TRUE, na = "null")
  if (!identical(status, "passed")) stop("Strict primary MR V1 freeze failed; partial outputs retained.", call. = FALSE)
  for (path in out[c("manifest", "json")]) {
    stop_if(file.exists(path), paste("Output appeared during run:", path))
    stop_if(!file.rename(paste0(path, ".partial"), path), paste("Atomic rename failed:", path))
  }
  log_line("freeze_status=passed approved_for_reverse_relaxed_mr_execution=TRUE manifest_sha256=", manifest_sha)
}

tryCatch(main(), error = function(e) {
  log_line("status=failed error=", conditionMessage(e))
  quit(status = 1L)
})
