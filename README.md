[![DOI](https://zenodo.org/badge/1341875967.svg)](https://doi.org/10.5281/zenodo.22047484)

Archived release v1.0.0: https://doi.org/10.5281/zenodo.22047485

# Haemoglobin and delirium: bidirectional Mendelian randomization

This repository candidate contains the frozen analysis code and validated summary outputs for a bidirectional two-sample Mendelian randomization study of haemoglobin concentration and delirium.

> **Release status:** `v1.0.0` upload package. The public environment is locked in `renv.lock`, and the frozen-result validator passed locally in R 4.6.1. Complete the private-repository inspection in `docs/release_steps.md` before changing the repository visibility to public.

## Study design

- Forward primary analysis: Vuckovic 2020 haemoglobin GWAS → FinnGen R13 delirium.
- Forward sensitivity analysis: Chen 2020 BCX2 haemoglobin GWAS → FinnGen R13 delirium.
- Reverse strict primary analysis: genome-wide significant FinnGen delirium instruments → haemoglobin.
- Reverse relaxed exploratory analysis: `P < 5×10⁻⁶` delirium instruments → haemoglobin.
- APOE-included and APOE-excluded branches were retained.
- Steiger analyses were supportive instrument-orientation sensitivity analyses and were not interpreted as causal proof.

## Repository contents

- `analysis/frozen_pipeline/`: selected final frozen scripts with original filenames.
- `renv.lock`: complete public R environment lockfile (R 4.6.1; 89 packages).
- `R/setup_environment.R`: restores the locked project-local R library.
- `R/validate_frozen_results.R`: independently validates the authoritative result matrix.
- `results/final/`: recovered final result and diagnostic matrices.
- `results/tables/`: frozen MR and sensitivity-analysis summary tables.
- `figures/`: final main and supplementary figures.
- `docs/`: data-source, code-map, validation and release documentation.
- `provenance/`: original incomplete lockfile and the documented reporting-direction correction.

## Data availability

GWAS files are not redistributed. See `docs/data_sources.md` for source identifiers, expected filenames and placement. The repository excludes raw and processed GWAS files, LD-reference files, local executables, internal manuscript materials and logs.

## Environment

The analysis used R 4.6.1, TwoSampleMR 0.7.9 at commit `3d119f20d6fc164b0c7f710f5590fee9580f2c7b`, MRPRESSO 1.0 at commit `3e3c92d7eda6dce0d1d66077373ec0f7ff4f7e87`, and psych 2.6.5. The public lockfile also records DBI 1.3.0, duckdb 1.5.5, renv 1.2.3 and all declared dependencies. See `environment/known_software_versions.csv`.

From the repository root:

```r
source("R/setup_environment.R")
source("R/validate_frozen_results.R")
```

The original analysis used staged manual clumping and project-specific freeze gates. The selected scripts therefore preserve the transparent frozen audit trail rather than claiming an unverified one-command raw-data rerun.

The frozen scripts also retain the original project-root literals used during the audited analysis. They contain no Windows username or credential, but some are not portable when run without `--project-root`. This provenance limitation is documented rather than silently altering the frozen source; the repository does not claim a verified one-command raw-data rerun.

## Key interpretation boundary

The frozen analysis found no robust evidence for a causal association between haemoglobin and delirium in either direction. This is not proof that an association is absent. One APOE-excluded relaxed reverse estimate was nominally significant but exploratory and not robust across haemoglobin outcome GWASs.

## License

The code is released under the MIT License. GWAS data remain subject to their original source terms.
