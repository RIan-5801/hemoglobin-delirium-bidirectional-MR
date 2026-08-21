
# Frozen pipeline map

The scripts in `analysis/frozen_pipeline/` preserve the final project-specific analysis and audit trail. They are intentionally retained with their original filenames so that the methods can be traced to the frozen outputs.

| Stage | Main authoritative scripts |
|---|---|
| Vuckovic Hb preparation and clumping | `02_*_v3.R`, `03_*_v1.R`, `04_*_v2.R`, `04_run_*_v2.ps1` |
| FinnGen extraction and harmonisation | `05_*_v9.R`, `06_*_v9_v7.R`, `07_*_v2.R`, `08_*_v4.R` |
| Forward primary MR | `09_*primary_v3.R` and its V2 freeze manifest |
| Reverse strict and relaxed MR | `10_*` through `20_*`, using the final version named in each prefix |
| Chen forward sensitivity | `21_*` through `25_*` |
| Chen reverse sensitivity | `26_*` through `33_*` |
| Steiger orientation sensitivity | `34_*` through `36_*` |
| Integrated result freeze and direction correction | `37_*`, `37a_*`, and `37b_*` |

Older superseded or failed versions, manuscript-writing scripts, literature-acquisition scripts, local binaries and empty history files were excluded from this public candidate.
