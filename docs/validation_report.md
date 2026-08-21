
# Validation report

## Overall assessment: Validated upload package

The frozen numerical results are internally consistent and suitable for inclusion in a code archive. A complete public `renv.lock` has been generated for R 4.6.1. It contains 89 packages, including DBI 1.3.0, duckdb 1.5.5, renv 1.2.3, psych 2.6.5, TwoSampleMR 0.7.9 and MRPRESSO 1.0. The lockfile parsed successfully and all required package records were present. The environment was restored locally, `renv` 1.2.3 was confirmed after restarting R, and all nine checks in `R/validate_frozen_results.R` passed on 2026-08-21.

## Verified checks

- Twelve authoritative analysis branches were present.
- All primary 95% confidence intervals matched `beta ± 1.96 × SE` within floating-point tolerance.
- All IVW and Wald-ratio primary P values matched independent two-sided normal calculations.
- Forward ORs and their confidence intervals matched exponentiation of the log-odds estimates.
- Heterogeneity P values matched chi-square calculations.
- Primary SNP counts were 312/312 for Vuckovic forward, 380/380 for Chen forward, 1/1 for strict reverse branches, and 10/9 for relaxed reverse branches.
- The nominally significant Vuckovic relaxed APOE-excluded result remained explicitly exploratory (`beta=0.0145566`, `P=0.0403466`) and was not robust in the Chen outcome sensitivity.
- Chen forward MR-PRESSO remained `technically_unavailable_under_frozen_configuration`; no Global Test P value or outlier count was invented.
- Seventy-six manifest-listed core artifacts available across the supplied code and result packages matched their recorded SHA-256 hashes; no mismatch was found among those core artifacts.
- The public lockfile SHA-256 is `d765ff790c4bd8edf8cd5aec854dd9cfbc2613aa6541af4e7fc47840078841c9`.
- TwoSampleMR was locked to commit `3d119f20d6fc164b0c7f710f5590fee9580f2c7b`; MRPRESSO was locked to commit `3e3c92d7eda6dce0d1d66077373ec0f7ff4f7e87`.

## Corrected reporting metadata

The early final matrices labelled the two Chen forward sensitivity rows as `delirium_to_Hb`. The recorded readback recovery corrected only those direction cells to `Hb_to_delirium`; numerical results and interpretation were unchanged. This public candidate uses the recovered matrices.

## Remaining repository checks before public visibility

1. The archived frozen scripts retain their original `E:/Research/hb_delirium_bidir_mr` path literals as provenance. These do not expose a Windows username, but they are not portable defaults. Run parameterized scripts with `--project-root`, and do not claim a verified one-command raw-data rerun until the path handling has been refactored and retested.
2. Upload to a private GitHub repository first and inspect the tracked files before changing visibility to public.
3. Connect the public repository to Zenodo and archive GitHub Release `v1.0.0`.
