#!/usr/bin/env Rscript

# Independent, read-only certification of the compressed FinnGen R13 input.
# This script does not select variants or create outcome-association datasets.

options(stringsAsFactors = FALSE, warn = 1)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Required installed package jsonlite is unavailable. Do not install packages automatically.", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || args[[1L]] != "--project-root" || !nzchar(args[[2L]])) {
  stop("Usage: Rscript.exe R/05a_certify_finngen_r13_input_v1.R --project-root E:/Research/hb_delirium_bidir_mr",
       call. = FALSE)
}

project_root <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
input_path <- file.path(project_root, "data_raw", "gwas", "finngen_R13_F5_DELIRIUM.gz")
expected_sha256 <- "85637f0f3358807964d4f8a3e500293168a706f1c08c65f3fc5512b65df40ed8"
expected_columns <- c("#chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls")
expected_header_raw <- paste(expected_columns, collapse = "\t")

qc_path <- file.path(project_root, "results", "qc", "finngen_R13_F5_DELIRIUM_input_certification_v1.json")
log_path <- file.path(project_root, "results", "logs", "finngen_R13_F5_DELIRIUM_input_certification_v1.log")
qc_partial_path <- paste0(qc_path, ".partial")
log_partial_path <- paste0(log_path, ".partial")
log_created_by_this_run <- FALSE

stop_if <- function(condition, message) {
  if (isTRUE(condition)) stop(message, call. = FALSE)
}

sha256_file <- function(path) {
  output <- suppressWarnings(system2("certutil.exe",
                                     args = c("-hashfile", shQuote(normalizePath(path, winslash = "\\")), "SHA256"),
                                     stdout = TRUE, stderr = TRUE))
  value_lines <- trimws(output[grepl("^[0-9A-Fa-f ]+$", output)])
  value <- tolower(gsub(" ", "", value_lines))
  value <- value[nchar(value) == 64L]
  stop_if(length(value) != 1L, sprintf("Unable to obtain one SHA-256 digest for: %s", path))
  value[[1L]]
}

log_line <- function(message) {
  cat(sprintf("%s %s\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), message), file = log_path, append = TRUE)
}

abort_with_log <- function(message) {
  log_line(paste0("FAILED: ", message))
  stop(message, call. = FALSE)
}

read_lines_checked <- function(connection, n) {
  read_warnings <- character()
  value <- withCallingHandlers(
    readLines(connection, n = n, warn = TRUE, encoding = "UTF-8"),
    warning = function(warning) {
      read_warnings <<- c(read_warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  if (length(read_warnings) > 0L) {
    abort_with_log(paste("Read warning(s) encountered during gzip scan:", paste(read_warnings, collapse = " | ")))
  }
  value
}

main <- function() {
  runtime_targets <- c(qc_path, log_path, qc_partial_path, log_partial_path)
  stop_if(any(file.exists(runtime_targets)),
          paste("A certification target or partial file already exists; refusing to overwrite:",
                paste(runtime_targets[file.exists(runtime_targets)], collapse = "; ")))
  stop_if(!file.exists(input_path), sprintf("Missing raw input: %s", input_path))

  input_info_before <- file.info(input_path)
  stop_if(is.na(input_info_before$size) || input_info_before$size <= 0, "Raw input file is empty or file size is unavailable.")

  dir.create(dirname(qc_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  log_created_by_this_run <<- TRUE
  log_line("Protocol=docs/protocol/analysis_plan_v1.1.md")
  log_line("Stage=independent read-only input certification; no target-rsID extraction or association processing.")
  on.exit({
    if (file.exists(qc_partial_path)) unlink(qc_partial_path, force = FALSE)
  }, add = TRUE)

  sha256_before <- sha256_file(input_path)
  stop_if(!identical(sha256_before, expected_sha256),
          sprintf("Input SHA-256 mismatch before scan: observed=%s expected=%s", sha256_before, expected_sha256))
  log_line(sprintf("Input pre-scan SHA-256 verified: %s", sha256_before))

  scan_start <- Sys.time()
  connection <- gzfile(input_path, open = "rt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)

  header_lines <- read_lines_checked(connection, n = 1L)
  header_present <- length(header_lines) == 1L
  if (!header_present) abort_with_log("Header is missing; no certification JSON will be published.")
  header_raw <- header_lines[[1L]]
  observed_column_names <- strsplit(header_raw, "\t", fixed = TRUE)[[1L]]
  observed_column_count <- length(observed_column_names)
  header_exact_match <- identical(header_raw, expected_header_raw) &&
    identical(observed_column_names, expected_columns) && observed_column_count == 13L
  if (!header_exact_match) {
    abort_with_log(sprintf("Header mismatch: observed_column_count=%s; observed_header_raw=%s", observed_column_count, header_raw))
  }

  observed_data_rows <- 0
  blank_data_line_count <- 0
  malformed_column_count <- 0
  read_or_parse_error_count <- 0
  chunk_size <- 100000L

  repeat {
    lines <- read_lines_checked(connection, n = chunk_size)
    if (length(lines) == 0L) break
    observed_data_rows <- observed_data_rows + length(lines)
    blank <- is.na(lines) | trimws(lines) == ""
    blank_data_line_count <- blank_data_line_count + sum(blank)
    tab_count <- nchar(lines, type = "chars") - nchar(gsub("\t", "", lines, fixed = TRUE), type = "chars")
    malformed <- is.na(tab_count) | tab_count != 12L
    malformed_column_count <- malformed_column_count + sum(malformed)
    if (any(blank) || any(malformed)) {
      read_or_parse_error_count <- read_or_parse_error_count + sum(blank | malformed)
      abort_with_log(sprintf("Data-line validation failed during complete gzip scan: blank=%s malformed_column_count=%s.",
                             blank_data_line_count, malformed_column_count))
    }
  }
  scan_end <- Sys.time()
  scan_duration_seconds <- as.numeric(difftime(scan_end, scan_start, units = "secs"))

  observed_physical_lines <- observed_data_rows + 1L
  stop_if(!header_exact_match, "Internal consistency error: physical-line count cannot be certified without exact header.")
  stop_if(observed_physical_lines != observed_data_rows + 1L, "Physical-line arithmetic validation failed.")

  sha256_after <- sha256_file(input_path)
  sha256_match_before <- identical(sha256_before, expected_sha256)
  sha256_match_after <- identical(sha256_after, expected_sha256) && identical(sha256_after, sha256_before)
  stop_if(!sha256_match_after,
          sprintf("Input SHA-256 mismatch after scan: observed=%s expected=%s pre_scan=%s", sha256_after, expected_sha256, sha256_before))
  input_info_after <- file.info(input_path)
  stop_if(is.na(input_info_after$size) || input_info_after$size != input_info_before$size,
          "Input file size changed during certification.")
  stop_if(!identical(input_info_after$mtime, input_info_before$mtime),
          "Input file modification time changed during certification.")

  certification <- list(
    certification_status = "passed",
    input_path = input_path,
    expected_sha256 = expected_sha256,
    sha256_before = sha256_before,
    sha256_after = sha256_after,
    sha256_match_before = sha256_match_before,
    sha256_match_after = sha256_match_after,
    file_size_bytes = as.numeric(input_info_before$size),
    last_write_time = format(input_info_before$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    scan_start_time = format(scan_start, "%Y-%m-%dT%H:%M:%S%z"),
    scan_end_time = format(scan_end, "%Y-%m-%dT%H:%M:%S%z"),
    scan_duration_seconds = scan_duration_seconds,
    gzip_complete_read = TRUE,
    read_or_parse_error_count = read_or_parse_error_count,
    header_present = header_present,
    observed_column_count = observed_column_count,
    observed_column_names = observed_column_names,
    header_exact_match = header_exact_match,
    header_raw = header_raw,
    observed_data_rows = observed_data_rows,
    observed_physical_lines = observed_physical_lines,
    certified_data_rows = observed_data_rows,
    certified_header_rows = 1L,
    certified_physical_lines = observed_data_rows + 1L,
    certified_input_sha256 = expected_sha256,
    blank_data_line_count = blank_data_line_count,
    malformed_column_count = malformed_column_count
  )

  writeLines(jsonlite::toJSON(certification, auto_unbox = TRUE, pretty = TRUE, null = "null"), qc_partial_path)
  stop_if(file.exists(qc_path), sprintf("Refusing to overwrite output that appeared during execution: %s", qc_path))
  stop_if(!file.rename(qc_partial_path, qc_path), "Atomic certification JSON rename failed.")
  stop_if(file.exists(qc_partial_path), "Partial certification JSON remains after successful rename.")
  log_line(sprintf("SUCCESS: certification_status=passed; observed_data_rows=%s; observed_physical_lines=%s; header_exact_match=%s; gzip_complete_read=TRUE.",
                   observed_data_rows, observed_physical_lines, header_exact_match))
  log_line(sprintf("Input post-scan SHA-256 verified: %s", sha256_after))
}

tryCatch(
  main(),
  error = function(error) {
    if (isTRUE(log_created_by_this_run) && file.exists(log_path)) log_line(paste0("TERMINATED: ", conditionMessage(error)))
    quit(status = 1L)
  }
)
