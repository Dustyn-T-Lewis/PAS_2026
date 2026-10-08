# 01_Proteins

Which proteins PAS changed in males.

| | |
|---|---|
| Reads | `02_Statistics/01_Fit/c_data/fit.rds`, `02_Statistics/02_Enrichment/c_data/enrichment.rds` |
| Writes | `c_data/01_proteins.xlsx`, `c_data/proteins.rds`, `b_reports/01_proteins.pdf` |
| Run | `Rscript 04_Male/01_Proteins/a_script/01_proteins.R` |

- **Method:** `report_proteins()` in `R/comparisons.R` keeps `PAS_vs_VEH_male` from `02_Statistics/01_Fit` (limpa `dpcDE`, BH within universe) and the STRING connectivity from `02_Statistics/02_Enrichment`. Tier: FDR ≤ 0.05, else Π ≤ 0.05, else p ≤ 0.05.
- **Report:** volcanoes for both universes, counts per tier and direction, a heatmap of every protein at p ≤ 0.05, then a grid of each one's values by mouse, with log2 FC, r with PAS, p, FDR and Π.
- **Result:** full proteome 133 proteins at p ≤ 0.05 (92 expected), 78 at Π ≤ 0.05, none at FDR ≤ 0.05; mitochondrial universe 80, 27, 0. The Π lists are more connected in STRING than chance (up: 29 proteins, 9 edges, 2.3 times expected; down: 49 proteins, 103 edges, 4.9 times).
- **Note:** Π ranks on p and fold change together and controls no error rate. Connected lists are expected when subunits of one complex move together, so the STRING result is descriptive.
