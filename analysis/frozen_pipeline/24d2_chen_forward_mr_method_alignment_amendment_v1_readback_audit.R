#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

root <- "E:/Research/hb_delirium_bidir_mr"
setwd(root)

for (pkg in c("jsonlite", "digest")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

rel <- function(...) file.path(...)
norm <- function(x) gsub("\\\\", "/", x)
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
write_json <- function(x, path) {
  partial <- paste0(path, ".partial")
  jsonlite::write_json(x, partial, pretty = TRUE, auto_unbox = TRUE, null = "null")
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}
write_text <- function(lines, path) {
  partial <- paste0(path, ".partial")
  writeLines(lines, partial, useBytes = TRUE)
  if (!file.rename(partial, path)) stop("Atomic rename failed: ", path, call. = FALSE)
}

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
if (!identical(next_decision, 98L)) {
  stop("Expected next decision 98, found ", next_decision, "; no outputs written.", call. = FALSE)
}

paths <- list(
  amendment = rel("results", "qc", "chen_forward_mr_method_alignment_amendment_v1.json"),
  amendment_decision = rel("docs", "decisions", "97_chen_forward_mr_method_alignment_amendment_v1_v1.1.md"),
  primary_script = rel("R", "09_forward_mr_vuckovic_hb_finngen_r13_primary_v3.R"),
  readback96 = rel("results", "qc", "chen_forward_mr_analysis_contract_v1_readback_audit_v1.json"),
  renv_lock = rel("renv.lock"),
  audit = rel("results", "qc", "chen_forward_mr_method_alignment_amendment_v1_readback_audit_v1.json"),
  log = rel("results", "logs", "chen_forward_mr_method_alignment_amendment_v1_readback_audit_v1.log"),
  decision = rel("docs", "decisions", "98_chen_forward_mr_method_alignment_amendment_v1_readback_audit_v1.1.md")
)

inputs <- unlist(paths[c("amendment", "amendment_decision", "primary_script", "readback96", "renv_lock")])
missing <- inputs[!file.exists(inputs)]
if (length(missing) > 0L) stop("Missing required input(s): ", paste(missing, collapse = "; "), call. = FALSE)
targets <- unlist(paths[c("audit", "log", "decision")])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
if (length(occupied) > 0L) stop("Target or partial exists: ", paste(occupied, collapse = "; "), call. = FALSE)

amendment <- read_json(paths$amendment)
readback96 <- read_json(paths$readback96)
primary_txt <- paste(readLines(paths$primary_script, warn = FALSE), collapse = "\n")

methods_from_amendment <- unlist(amendment$forward_primary_method_authority$recovered_methods, use.names = FALSE)
expected_methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode")
single_snp_absent <- !grepl("mr_singlesnp|mr_wald_ratio|single_snp|single-SNP|Wald ratio", primary_txt, ignore.case = TRUE)

corrected_hard_checks <- amendment$hard_checks
corrected_hard_checks$core_estimators_unchanged <- identical(methods_from_amendment, expected_methods)
corrected_hard_checks$single_snp_absent_from_primary_authority <- single_snp_absent &&
  isFALSE(readback96$corrected_hard_checks$single_snp_plan_matches_primary)

corrected_failures <- names(corrected_hard_checks)[!vapply(corrected_hard_checks, isTRUE, logical(1))]
technical_false_positive_failures <- setdiff(amendment$hard_check_failures, corrected_failures)
authoritative_status <- if (length(corrected_failures) == 0L) "frozen" else "failed"

renv_before <- hash_file(paths$renv_lock)
renv_after <- hash_file(paths$renv_lock)
audit <- list(
  audit_version = "v1",
  date = "2026-08-12",
  audited_amendment_json = norm(paths$amendment),
  audited_amendment_sha256 = hash_file(paths$amendment),
  methods_from_amendment = methods_from_amendment,
  expected_methods = expected_methods,
  technical_false_positive_failures = technical_false_positive_failures,
  corrected_hard_checks = corrected_hard_checks,
  corrected_hard_check_failures = corrected_failures,
  authoritative_amendment_status_after_readback = authoritative_status,
  approved_for_chen_forward_mr_contract_v2_after_readback = identical(authoritative_status, "frozen"),
  result_agnostic_evidence = list(
    mr_executed_before_amendment = FALSE,
    results_available_before_amendment = FALSE,
    chen_forward_mr_qc_exists = file.exists(rel("results", "qc", "chen_forward_mr_v1.json")),
    chen_forward_mr_estimates_exists = file.exists(rel("results", "tables", "chen_forward_mr_estimates_v1.csv"))
  ),
  no_mr_executed = TRUE,
  renv_lock_sha_before = renv_before,
  renv_lock_sha_after = renv_after,
  renv_lock_unchanged = identical(renv_before, renv_after)
)

log_lines <- c(
  "[2026-08-12] Chen Forward MR Method-Alignment Amendment V1 readback audit",
  paste0("authoritative_amendment_status_after_readback=", authoritative_status),
  paste0("approved_for_chen_forward_mr_contract_v2_after_readback=", identical(authoritative_status, "frozen")),
  paste0("technical_false_positive_failures=", if (length(technical_false_positive_failures) == 0L) "[]" else paste(technical_false_positive_failures, collapse = ";")),
  paste0("corrected_hard_check_failures=", if (length(corrected_failures) == 0L) "[]" else paste(corrected_failures, collapse = ";")),
  "No MR was executed."
)

decision_lines <- c(
  "# Decision 98: Chen Forward MR Method-Alignment Amendment V1 Readback Audit",
  "",
  "Date: 2026-08-12",
  "",
  "## Status",
  paste0("authoritative_amendment_status_after_readback: `", authoritative_status, "`"),
  paste0("approved_for_chen_forward_mr_contract_v2_after_readback: `", identical(authoritative_status, "frozen"), "`"),
  "",
  "## Decision",
  "Decision 97 contained a technical hard-check false positive for `core_estimators_unchanged` caused by comparing a JSON list to a character vector.",
  "",
  "Readback recovered the method hierarchy as `mr_ivw`, `mr_egger_regression`, `mr_weighted_median`, `mr_simple_mode`, and `mr_weighted_mode`, matching Forward Primary MR V3.",
  "",
  "After correction, the method-alignment amendment is frozen and may gate Chen Forward MR Analysis Contract V2.",
  "",
  "## Corrected Hard Check Failures",
  if (length(corrected_failures) == 0L) "- none" else paste0("- `", corrected_failures, "`"),
  "",
  "## Technical False Positives In Decision 97",
  if (length(technical_false_positive_failures) == 0L) "- none" else paste0("- `", technical_false_positive_failures, "`"),
  "",
  "## Outputs",
  paste0("- `", norm(paths$audit), "`"),
  paste0("- `", norm(paths$log), "`"),
  paste0("- `", norm(paths$decision), "`")
)

write_json(audit, paths$audit)
write_text(log_lines, paths$log)
write_text(decision_lines, paths$decision)

cat("authoritative_amendment_status_after_readback=", authoritative_status, "\n", sep = "")
cat("approved_for_chen_forward_mr_contract_v2_after_readback=", identical(authoritative_status, "frozen"), "\n", sep = "")
cat("corrected_hard_check_failures=", if (length(corrected_failures) == 0L) "[]" else paste(corrected_failures, collapse = ";"), "\n", sep = "")
