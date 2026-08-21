#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--project-root")) {
  root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (length(args) == 0L) {
  root <- "E:/Research/hb_delirium_bidir_mr"
} else {
  stop("Usage: Rscript R/25a_chen_forward_mr_presso_recovery_v1.R [--project-root E:/Research/hb_delirium_bidir_mr]", call. = FALSE)
}
setwd(root)

source(file.path(root, "renv", "activate.R"))
mr_library <- file.path(root, "renv", "mr-v1-library")
.libPaths(c(mr_library, .libPaths()))

for (pkg in c("DBI", "duckdb", "jsonlite", "digest", "MRPRESSO")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg, call. = FALSE)
}

rel <- function(...) file.path(root, ...)
norm <- function(x) gsub("\\\\", "/", x)
hash_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
read_json <- function(path) jsonlite::fromJSON(path, simplifyVector = FALSE)
stop_if <- function(x, msg) if (isTRUE(x)) stop(msg, call. = FALSE)
records <- function(x) if (!is.data.frame(x)) list() else lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
sql_path <- function(path, must_work = TRUE) {
  paste0("'", gsub("'", "''", normalizePath(path, winslash = "/", mustWork = must_work), fixed = TRUE), "'")
}
read_parquet <- function(con, path) DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_path(path)))
write_csv_precise <- function(x, path) {
  old <- options(digits = 17, scipen = 999)
  on.exit(options(old), add = TRUE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}
write_json_precise <- function(x, path) {
  jsonlite::write_json(x, path, pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)
}
extract_int <- function(pattern, txt) {
  hit <- regmatches(txt, regexpr(pattern, txt, perl = TRUE))
  if (length(hit) == 1L && nzchar(hit)) as.integer(gsub("[^0-9]", "", hit)) else NA_integer_
}
extract_num_after_equal <- function(pattern, txt) {
  hit <- regmatches(txt, regexpr(pattern, txt, perl = TRUE))
  if (length(hit) == 1L && nzchar(hit)) as.numeric(sub(".*=\\s*", "", hit)) else NA_real_
}

paths <- list(
  freeze = rel("results", "qc", "chen_forward_harmonised_mr_inputs_v1_freeze.json"),
  contract = rel("results", "qc", "chen_forward_mr_analysis_contract_v2.json"),
  mr_qc = rel("results", "qc", "chen_forward_mr_v1.json"),
  mr_presso_original = rel("results", "tables", "chen_forward_mr_presso_v1.csv"),
  mr_script = rel("R", "25_chen_forward_mr_v1.R"),
  primary_script = rel("R", "09_forward_mr_vuckovic_hb_finngen_r13_primary_v3.R"),
  primary_qc = rel("results", "qc", "vuckovic_hb_finngen_r13_forward_mr_v3.json"),
  attempt1 = rel("results", "logs", "chen_forward_mr_v1_attempt1_hung.log"),
  attempt2 = rel("results", "logs", "chen_forward_mr_v1_attempt2_mrpresso_timeout_checks_failed.log"),
  attempt3 = rel("results", "logs", "chen_forward_mr_v1_attempt3_mrpresso_timeout_comparison_rows_failed.log"),
  attempt4 = rel("results", "logs", "chen_forward_mr_v1_attempt4_mrpresso_timeout_comparison_note_failed.log"),
  mr_log = rel("results", "logs", "chen_forward_mr_v1.log"),
  renv_lock = rel("renv.lock"),
  recovery_table = rel("results", "tables", "chen_forward_mr_presso_recovery_v1.csv"),
  recovery_qc = rel("results", "qc", "chen_forward_mr_presso_recovery_v1.json"),
  recovery_log = rel("results", "logs", "chen_forward_mr_presso_recovery_v1.log"),
  decision = rel("docs", "decisions", "101_chen_forward_mr_presso_technical_recovery_v1_v1.1.md")
)

required_inputs <- unlist(paths[c("freeze", "contract", "mr_qc", "mr_presso_original", "mr_script", "primary_script", "primary_qc", "mr_log", "renv_lock")])
missing <- required_inputs[!file.exists(required_inputs)]
if (length(missing) > 0L) stop("Missing required input(s): ", paste(missing, collapse = "; "), call. = FALSE)

decision_files <- list.files(rel("docs", "decisions"), pattern = "^[0-9]+_.*\\.md$", full.names = FALSE)
decision_numbers <- suppressWarnings(as.integer(sub("^([0-9]+)_.*$", "\\1", decision_files)))
next_decision <- max(decision_numbers, na.rm = TRUE) + 1L
stop_if(!identical(next_decision, 101L), paste0("Expected next decision 101, found ", next_decision, "; no outputs written."))

targets <- unlist(paths[c("recovery_table", "recovery_qc", "recovery_log", "decision")])
occupied <- targets[file.exists(targets) | file.exists(paste0(targets, ".partial"))]
stop_if(length(occupied) > 0L, paste("Target or partial exists:", paste(occupied, collapse = "; ")))
dir.create(dirname(paths$recovery_table), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(paths$recovery_qc), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(paths$recovery_log), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(paste0("[", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "] ", paste0(..., collapse = ""), "\n"),
                              file = paths$recovery_log, append = TRUE)

make_dat <- function(source, analysis_set, analysis_role) {
  required <- c("resolved_rsid", "exposure_beta", "exposure_se", "outcome_beta_harmonised", "outcome_se_harmonised", "final_valid")
  stop_if(!all(required %in% names(source)), paste("Required columns absent:", analysis_set))
  stop_if(!all(source$final_valid), paste("Non-final-valid row found:", analysis_set))
  data.frame(
    SNP = as.character(source$resolved_rsid),
    beta.exposure = as.numeric(source$exposure_beta),
    se.exposure = as.numeric(source$exposure_se),
    beta.outcome = as.numeric(source$outcome_beta_harmonised),
    se.outcome = as.numeric(source$outcome_se_harmonised),
    analysis_set = analysis_set,
    analysis_role = analysis_role,
    stringsAsFactors = FALSE
  )
}

run_presso_set <- function(dat, analysis_set, analysis_role, nb_distribution, signif_threshold, seed_value) {
  log_line("mr_presso_recovery_start analysis_set=", analysis_set, " NbDistribution=", nb_distribution)
  set.seed(seed_value)
  x <- dat[, c("beta.outcome", "beta.exposure", "se.outcome", "se.exposure"), drop = FALSE]
  rownames(x) <- dat$SNP
  start_time <- Sys.time()
  result <- tryCatch(
    MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome",
      BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome",
      SdExposure = "se.exposure",
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      data = x,
      NbDistribution = nb_distribution,
      SignifThreshold = signif_threshold
    ),
    error = function(e) e
  )
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if (inherits(result, "error")) {
    log_line("mr_presso_recovery_failed analysis_set=", analysis_set, " elapsed_seconds=", round(elapsed, 3),
             " error=", conditionMessage(result))
    return(list(
      status = "failed",
      rows = data.frame(
        analysis_set = analysis_set, analysis_role = analysis_role, test_type = "MR-PRESSO",
        metric = "status", value = "failed", pval = NA_real_, outlier_rsid = "",
        notes = conditionMessage(result), elapsed_seconds = elapsed, stringsAsFactors = FALSE
      ),
      raw_result = NULL,
      error = conditionMessage(result),
      elapsed_seconds = elapsed
    ))
  }
  log_line("mr_presso_recovery_completed analysis_set=", analysis_set, " elapsed_seconds=", round(elapsed, 3))
  root <- result[["MR-PRESSO results"]]
  global <- if (!is.null(root)) root[["Global Test"]] else NULL
  outlier <- if (!is.null(root)) root[["Outlier Test"]] else NULL
  distortion <- if (!is.null(root)) root[["Distortion Test"]] else NULL
  main <- result[["Main MR results"]]
  global_p <- if (!is.null(global) && "Pvalue" %in% names(global)) as.character(global$Pvalue[[1L]]) else NA_character_
  rows <- list(
    data.frame(
      analysis_set = analysis_set, analysis_role = analysis_role, test_type = "Global Test",
      metric = "RSSobs", value = if (!is.null(global) && "RSSobs" %in% names(global)) as.character(global$RSSobs[[1L]]) else "",
      pval = suppressWarnings(as.numeric(global_p)), outlier_rsid = "",
      notes = paste0("mr_presso_status=completed; global_pvalue_raw=", global_p),
      elapsed_seconds = elapsed, stringsAsFactors = FALSE
    )
  )
  ids <- if (!is.null(outlier) && nrow(outlier) > 0L) rownames(outlier) else character(0)
  rows[[length(rows) + 1L]] <- data.frame(
    analysis_set = analysis_set, analysis_role = analysis_role, test_type = "Outlier Test",
    metric = "outlier_count", value = as.character(length(ids)), pval = NA_real_,
    outlier_rsid = if (length(ids) == 0L) "" else paste(ids, collapse = ";"),
    notes = if (length(ids) == 0L) "no_outlier_reported" else "outlier_reported_sensitivity_only_no_main_input_change",
    elapsed_seconds = elapsed, stringsAsFactors = FALSE
  )
  if (!is.null(distortion)) {
    rows[[length(rows) + 1L]] <- data.frame(
      analysis_set = analysis_set, analysis_role = analysis_role, test_type = "Distortion Test",
      metric = "raw_result", value = paste(capture.output(print(distortion)), collapse = " | "),
      pval = NA_real_, outlier_rsid = "", notes = "distortion_output_if_available",
      elapsed_seconds = elapsed, stringsAsFactors = FALSE
    )
  }
  if (!is.null(main)) {
    rows[[length(rows) + 1L]] <- data.frame(
      analysis_set = analysis_set, analysis_role = analysis_role, test_type = "Main MR results",
      metric = "raw_result", value = paste(capture.output(print(main)), collapse = " | "),
      pval = NA_real_, outlier_rsid = "",
      notes = "MR-PRESSO main/corrected estimates are sensitivity diagnostics only and do not replace Decision 100 full-input IVW",
      elapsed_seconds = elapsed, stringsAsFactors = FALSE
    )
  }
  list(status = "completed", rows = do.call(rbind, rows), raw_result = result, error = NA_character_, elapsed_seconds = elapsed)
}

main <- function() {
  log_line("stage=chen_forward_mr_presso_recovery_v1_start")
  renv_before <- hash_file(paths$renv_lock)
  freeze <- read_json(paths$freeze)
  contract <- read_json(paths$contract)
  mr_qc <- read_json(paths$mr_qc)
  primary_qc <- read_json(paths$primary_qc)
  mr_script_txt <- paste(readLines(paths$mr_script, warn = FALSE), collapse = "\n")
  primary_script_txt <- paste(readLines(paths$primary_script, warn = FALSE), collapse = "\n")
  original_presso <- read.csv(paths$mr_presso_original, stringsAsFactors = FALSE, check.names = FALSE)

  nb_distribution <- contract$diagnostics$mr_presso$NbDistribution
  signif_threshold <- contract$diagnostics$mr_presso$SignifThreshold
  seed_value <- contract$diagnostics$mr_presso$seed
  primary_nb <- extract_int("NbDistribution\\s*=\\s*[0-9]+", primary_script_txt)
  primary_signif <- extract_num_after_equal("SignifThreshold\\s*=\\s*[0-9.]+", primary_script_txt)
  original_timeout_match <- regmatches(mr_script_txt, regexpr("setTimeLimit\\(cpu\\s*=\\s*Inf,\\s*elapsed\\s*=\\s*[0-9]+", mr_script_txt, perl = TRUE))
  original_timeout_value <- if (length(original_timeout_match) == 1L && nzchar(original_timeout_match)) {
    as.integer(gsub("[^0-9]", "", original_timeout_match))
  } else {
    NA_integer_
  }

  included_path <- rel("data_derived", "forward_sensitivity_harmonisation", "chen_forward_harmonised_apoe_included_v1.parquet")
  excluded_path <- rel("data_derived", "forward_sensitivity_harmonisation", "chen_forward_harmonised_apoe_excluded_v1.parquet")
  if (!is.null(freeze$manifest_records)) {
    roles <- vapply(freeze$manifest_records, function(x) as.character(x$file_role), character(1))
    rels <- vapply(freeze$manifest_records, function(x) as.character(x$relative_path), character(1))
    if ("harmonised_apoe_included_parquet" %in% roles) included_path <- rel(rels[roles == "harmonised_apoe_included_parquet"][[1L]])
    if ("harmonised_apoe_excluded_parquet" %in% roles) excluded_path <- rel(rels[roles == "harmonised_apoe_excluded_parquet"][[1L]])
  }
  included_sha_before <- hash_file(included_path)
  excluded_sha_before <- hash_file(excluded_path)

  pre_checks <- list(
    decision_100_core_mr_gate = identical(mr_qc$mr_status, "passed") &&
      isTRUE(mr_qc$approved_for_chen_forward_results_interpretation) &&
      length(mr_qc$hard_check_failures) == 0L,
    original_mr_presso_failed_due_to_elapsed_limit = all(original_presso$value == "failed_or_timeout") &&
      all(grepl("elapsed time limit", original_presso$notes, ignore.case = TRUE)),
    no_input_or_scientific_error_in_decision_100 = isTRUE(mr_qc$hard_checks$mr_presso_configuration_matches_primary) &&
      isTRUE(mr_qc$hard_checks$input_conversion_preserved) &&
      isTRUE(mr_qc$hard_checks$software_environment_matches_primary),
    scientific_parameters_match_contract = identical(as.integer(nb_distribution), 10000L) &&
      isTRUE(all.equal(as.numeric(signif_threshold), 0.05)) &&
      identical(as.integer(seed_value), 2026L),
    scientific_parameters_match_primary = identical(as.integer(primary_nb), 10000L) &&
      isTRUE(all.equal(as.numeric(primary_signif), 0.05)),
    mrpresso_package_matches_primary = identical(as.character(utils::packageVersion("MRPRESSO")), primary_qc$MRPRESSO_version) &&
      identical(as.character(utils::packageDescription("MRPRESSO")[["RemoteSha"]]), primary_qc$MRPRESSO_RemoteSha),
    original_timeout_source_identified = grepl("setTimeLimit", mr_script_txt, fixed = TRUE) &&
      identical(original_timeout_value, 180L),
    no_single_snp = isFALSE(mr_qc$single_snp_run),
    no_steiger = isFALSE(mr_qc$steiger_run)
  )
  stop_if(!all(unlist(pre_checks)), paste("Recovery authority gate failed:", paste(names(pre_checks)[!unlist(pre_checks)], collapse = "; ")))

  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  inc_src <- read_parquet(con, included_path)
  exc_src <- read_parquet(con, excluded_path)
  inc_dat <- make_dat(inc_src, "APOE_included", "chen_forward_alternative_hb_gwas_sensitivity_main")
  exc_dat <- make_dat(exc_src, "APOE_excluded", "chen_forward_alternative_hb_gwas_apoe_exclusion_sensitivity")
  input_checks <- list(
    included_count_matches_freeze = nrow(inc_dat) == as.integer(freeze$included_final_valid_count),
    excluded_count_matches_freeze = nrow(exc_dat) == as.integer(freeze$excluded_final_valid_count),
    included_rsids_match_freeze = setequal(inc_dat$SNP, freeze$included_final_rsids),
    excluded_rsids_match_freeze = setequal(exc_dat$SNP, freeze$excluded_final_rsids),
    included_numeric_fields_finite = all(is.finite(inc_dat$beta.exposure), is.finite(inc_dat$se.exposure), is.finite(inc_dat$beta.outcome), is.finite(inc_dat$se.outcome)),
    excluded_numeric_fields_finite = all(is.finite(exc_dat$beta.exposure), is.finite(exc_dat$se.exposure), is.finite(exc_dat$beta.outcome), is.finite(exc_dat$se.outcome)),
    included_se_positive = all(inc_dat$se.exposure > 0, inc_dat$se.outcome > 0),
    excluded_se_positive = all(exc_dat$se.exposure > 0, exc_dat$se.outcome > 0)
  )
  stop_if(!all(unlist(input_checks)), paste("Recovery input gate failed:", paste(names(input_checks)[!unlist(input_checks)], collapse = "; ")))

  included_result <- run_presso_set(inc_dat, "APOE_included", "chen_forward_alternative_hb_gwas_sensitivity_main", nb_distribution, signif_threshold, seed_value)
  gc(verbose = FALSE)
  excluded_result <- run_presso_set(exc_dat, "APOE_excluded", "chen_forward_alternative_hb_gwas_apoe_exclusion_sensitivity", nb_distribution, signif_threshold, seed_value)
  recovery_table <- rbind(included_result$rows, excluded_result$rows)

  included_sha_after <- hash_file(included_path)
  excluded_sha_after <- hash_file(excluded_path)
  renv_after <- hash_file(paths$renv_lock)
  recovery_status <- if (identical(included_result$status, "completed") && identical(excluded_result$status, "completed")) {
    "completed"
  } else if (grepl("time|elapsed|timeout", paste(c(included_result$error, excluded_result$error), collapse = " "), ignore.case = TRUE)) {
    "not_completed_under_frozen_configuration_due_to_computational_timeout"
  } else {
    "failed_non_timeout"
  }

  hard_checks <- c(pre_checks, input_checks, list(
    scientific_parameters_changed = FALSE,
    input_sets_changed = FALSE,
    execution_architecture_only_changed = TRUE,
    seed_policy_matches_decision_100 = TRUE,
    included_attempted = TRUE,
    excluded_attempted = TRUE,
    exact_errors_preserved_if_failed = all(!is.na(c(included_result$error, excluded_result$error)[c(included_result$status, excluded_result$status) != "completed"])),
    no_outlier_deletion = TRUE,
    no_main_input_change = identical(included_sha_before, included_sha_after) && identical(excluded_sha_before, excluded_sha_after),
    renv_lock_unchanged = identical(renv_before, renv_after)
  ))
  failures <- names(hard_checks)[!vapply(hard_checks, isTRUE, logical(1))]

  qc <- list(
    recovery_version = "v1",
    recovery_type = "technical_execution_recovery",
    source_mr_decision = 100,
    scientific_parameters_changed = FALSE,
    input_sets_changed = FALSE,
    mr_presso_version = as.character(utils::packageVersion("MRPRESSO")),
    mr_presso_sha = as.character(utils::packageDescription("MRPRESSO")[["RemoteSha"]]),
    NbDistribution = nb_distribution,
    SignifThreshold = signif_threshold,
    seed_policy = "set.seed(2026) immediately before each analysis set, matching Decision 100 run_presso semantics",
    original_timeout_source = "R setTimeLimit in R/25_chen_forward_mr_v1.R",
    original_timeout_value = original_timeout_value,
    recovery_timeout_policy = "no R setTimeLimit applied; bounded by external execution monitoring only",
    technical_execution_changes = list(
      included_excluded_run_sequentially = TRUE,
      removed_r_setTimeLimit = TRUE,
      separate_recovery_log = TRUE,
      garbage_collection_between_sets = TRUE,
      scientific_parameters_changed = FALSE
    ),
    included_status = included_result$status,
    excluded_status = excluded_result$status,
    included_result = records(included_result$rows),
    excluded_result = records(excluded_result$rows),
    recovery_status = recovery_status,
    hard_checks = hard_checks,
    hard_check_failures = failures,
    input_sha256 = list(
      included_before = included_sha_before,
      included_after = included_sha_after,
      excluded_before = excluded_sha_before,
      excluded_after = excluded_sha_after
    ),
    renv_lock_sha_before = renv_before,
    renv_lock_sha_after = renv_after
  )

  partial_table <- paste0(paths$recovery_table, ".partial")
  partial_qc <- paste0(paths$recovery_qc, ".partial")
  partial_decision <- paste0(paths$decision, ".partial")
  on.exit(unlink(c(partial_table, partial_qc, partial_decision), force = TRUE), add = TRUE)
  write_csv_precise(recovery_table, partial_table)
  write_json_precise(qc, partial_qc)

  decision_lines <- c(
    "# Decision 101: Chen Forward MR-PRESSO Technical Recovery V1",
    "",
    "Date: 2026-08-12",
    "",
    "## Status",
    paste0("recovery_status: `", recovery_status, "`"),
    paste0("included_status: `", included_result$status, "`"),
    paste0("excluded_status: `", excluded_result$status, "`"),
    "",
    "## Decision",
    "MR-PRESSO was retried as a technical execution recovery after Decision 100 elapsed-time failures.",
    "",
    "The recovery changed only execution architecture: included and excluded sets were run sequentially, the artificial R `setTimeLimit(elapsed=180)` was removed, a separate log was written, and garbage collection was performed between sets.",
    "",
    "Scientific parameters were unchanged: `NbDistribution=10000`, `SignifThreshold=0.05`, seed `2026`, MRPRESSO version/SHA, and Decision 94 frozen input sets.",
    "",
    "## Hard Check Failures",
    if (length(failures) == 0L) "- none" else paste0("- `", failures, "`"),
    "",
    "## Outputs",
    paste0("- `", norm(paths$recovery_table), "`"),
    paste0("- `", norm(paths$recovery_qc), "`"),
    paste0("- `", norm(paths$recovery_log), "`"),
    paste0("- `", norm(paths$decision), "`")
  )
  writeLines(decision_lines, partial_decision, useBytes = TRUE)

  stop_if(length(failures) > 0L, paste("Recovery hard checks failed:", paste(failures, collapse = "; ")))
  for (p in c(paths$recovery_table, paths$recovery_qc, paths$decision)) {
    stop_if(file.exists(p), paste("Output appeared during run:", p))
    if (!file.rename(paste0(p, ".partial"), p)) stop("Atomic rename failed: ", p, call. = FALSE)
  }
  log_line("recovery_status=", recovery_status)
  log_line("hard_check_failures=[]")
}

tryCatch(main(), error = function(e) {
  log_line("recovery_status=failed")
  log_line("error=", conditionMessage(e))
  quit(status = 1L)
})
