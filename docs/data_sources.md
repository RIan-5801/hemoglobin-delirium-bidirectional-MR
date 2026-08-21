
# Data sources

The repository does not redistribute GWAS summary-statistics files. Obtain the files from their original sources and place them under `data_raw/gwas/` using the filenames below.

| Role | Study/source | Local filename used by the frozen scripts |
|---|---|---|
| Haemoglobin primary GWAS | Vuckovic et al. 2020; GWAS Catalog GCST90002384 | `ebi-a-GCST90002384.vcf.gz` |
| Delirium GWAS | FinnGen Release 13; endpoint F5_DELIRIUM | `finngen_R13_F5_DELIRIUM.gz` |
| Alternative haemoglobin GWAS | Chen et al. 2020; BCX2 European HGB GWAMA | `BCX2_HGB_EA_GWAMA.out.gz` |

LD clumping also requires a compatible European ancestry reference panel. The frozen analysis used the 1000 Genomes Phase 3 European reference under `resources/ld/1kg_v3/EUR`. Reference files and PLINK2 binaries are not redistributed.

Users must comply with the access and reuse terms of each original data source.
