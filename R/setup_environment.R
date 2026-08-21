
#!/usr/bin/env Rscript

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
lockfile <- file.path(root, "renv.lock")
lib <- file.path(root, "renv", "public-v1-library")

if (!file.exists(lockfile)) stop("Missing lockfile: ", lockfile, call. = FALSE)
dir.create(lib, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

renv::restore(
  project = root,
  library = lib,
  lockfile = lockfile,
  prompt = FALSE
)

.libPaths(c(lib, .libPaths()))
packages <- c("DBI", "digest", "duckdb", "jsonlite", "renv", "psych", "TwoSampleMR", "MRPRESSO")
versions <- vapply(packages, function(x) as.character(utils::packageVersion(x)), character(1))
print(data.frame(package = packages, version = versions), row.names = FALSE)
cat("Locked environment restored in:", lib, "\n")
