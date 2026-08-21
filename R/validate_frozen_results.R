
#!/usr/bin/env Rscript

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
path <- file.path(root, "results", "final", "final_primary_result_matrix.csv")
if (!file.exists(path)) stop("Missing final result matrix: ", path, call. = FALSE)

x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
expected_ids <- c(
  "forward_primary_apoe_included", "forward_primary_apoe_excluded",
  "chen_forward_included", "chen_forward_excluded",
  "reverse_strict_included", "reverse_strict_excluded",
  "reverse_relaxed_included", "reverse_relaxed_excluded",
  "chen_reverse_strict_included", "chen_reverse_strict_excluded",
  "chen_reverse_relaxed_included", "chen_reverse_relaxed_excluded"
)
if (!identical(as.character(x$analysis_id), expected_ids)) stop("Unexpected analysis IDs or order.", call. = FALSE)

tol <- 1e-12
lower <- x$primary_beta - 1.96 * x$primary_se
upper <- x$primary_beta + 1.96 * x$primary_se
pcalc <- 2 * stats::pnorm(abs(x$primary_beta / x$primary_se), lower.tail = FALSE)

checks <- c(
  row_count = nrow(x) == 12L,
  ci_lower = max(abs(lower - x$primary_ci_lower)) < tol,
  ci_upper = max(abs(upper - x$primary_ci_upper)) < tol,
  primary_p = max(abs(pcalc - x$primary_p)) < tol,
  nominal_flags = identical(as.logical(x$nominal_p_lt_0_05), x$primary_p < 0.05),
  forward_direction = all(x$direction[x$analysis_id %in% c(
    "forward_primary_apoe_included", "forward_primary_apoe_excluded",
    "chen_forward_included", "chen_forward_excluded"
  )] == "Hb_to_delirium"),
  reverse_direction = all(x$direction[grepl("^reverse_|^chen_reverse_", x$analysis_id)] == "delirium_to_Hb"),
  snp_counts = identical(as.integer(x$nsnp), c(312L, 312L, 380L, 380L, 1L, 1L, 10L, 9L, 1L, 1L, 10L, 9L)),
  chen_forward_presso = all(x$mr_presso_status[x$analysis_id %in% c(
    "chen_forward_included", "chen_forward_excluded"
  )] == "technically_unavailable_under_frozen_configuration")
)

print(data.frame(check = names(checks), passed = unname(checks)), row.names = FALSE)
if (!all(checks)) stop("Frozen-result validation failed.", call. = FALSE)
cat("All frozen-result checks passed.
")
