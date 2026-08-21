#!/usr/bin/env Rscript

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
required <- c("renv", "DBI", "duckdb", "jsonlite", "digest")
input <- file.path(project_root, "data_raw", "gwas", "ebi-a-GCST90002384.vcf.gz")

if (!file.exists(input)) stop("Input VCF not found: ", input)
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!file.exists(file.path(project_root, "renv", "activate.R"))) {
  r_minor <- strsplit(R.version$minor, "\\.")[[1]][1]
  bootstrap_lib <- file.path(project_root, "renv", "library", paste0("R-", R.version$major, ".", r_minor), R.version$platform)
  dir.create(bootstrap_lib, recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("renv", lib.loc = bootstrap_lib, quietly = TRUE)) {
    install.packages("renv", lib = bootstrap_lib, dependencies = FALSE)
  }
  library(renv, lib.loc = bootstrap_lib)
  renv::init(project = project_root, bare = TRUE, restart = FALSE)
}

source(file.path(project_root, "renv", "activate.R"))
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) renv::install(missing, prompt = FALSE)
renv::snapshot(project = project_root, prompt = FALSE, force = TRUE)

installed <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
if (!all(installed)) stop("Missing required project packages: ", paste(required[!installed], collapse = ", "))

report <- data.frame(
  item = c("R", required),
  version = c(R.version.string, vapply(required, function(x) as.character(utils::packageVersion(x)), character(1))),
  stringsAsFactors = FALSE
)
print(report, row.names = FALSE)
