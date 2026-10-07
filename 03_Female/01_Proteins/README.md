# 01_Proteins

Which proteins PAS changed in females.

| | |
|---|---|
| Reads | `02_Statistics/01_Fit/c_data/fit.rds`, `02_Statistics/02_Enrichment/c_data/enrichment.rds` |
| Writes | `c_data/01_proteins.xlsx`, `c_data/proteins.rds`, `b_reports/01_proteins.pdf` |
| Run | `Rscript 03_Female/01_Proteins/a_script/01_proteins.R` |

- **Method:** `report_proteins()` in `R/comparisons.R` keeps `PAS_vs_VEH_female` from `02_Statistics/01_Fit` (limpa `dpcDE`, BH within universe) and the STRING connectivity from `02_Statistics/02_Enrichment`. Tier: FDR ≤ 0.05, else Π ≤ 0.05, else p ≤ 0.05.
- **Report:** volcanoes for both universes, counts per tier and direction, a heatmap of every protein at p ≤ 0.05, then a per-mouse grid of every one of them with log2 FC, r with PAS, p, FDR and Π.
- **Result:** full proteome 152 proteins at p ≤ 0.05 (92 expected by chance), 109 at Π ≤ 0.05, 4 at FDR ≤ 0.05 (Dusp3 and Cand2 down; Arl6ip5 and Bcap31 up); mitochondrial universe 47, 18, 0. The Π lists are more connected in STRING than chance (up: 49 proteins, 65 edges, 2.3 times expected; down: 60 proteins, 65 edges, 3.4 times).
- **Note:** Π ranks on p and fold change together and controls no error rate. Connected lists are expected when subunits of one complex move together, so the STRING result is descriptive.
