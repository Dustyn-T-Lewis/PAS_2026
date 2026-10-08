# 01_Proteins

Which proteins respond to PAS differently in females and males.

| | |
|---|---|
| Reads | `02_Statistics/01_Fit/c_data/fit.rds`, `02_Statistics/02_Enrichment/c_data/enrichment.rds` |
| Writes | `c_data/01_proteins.xlsx`, `c_data/proteins.rds`, `b_reports/01_proteins.pdf` |
| Run | `Rscript 05_Sex_Difference/01_Proteins/a_script/01_proteins.R` |

- **Method:** `report_proteins()` keeps `treatment_by_sex` from `02_Statistics/01_Fit`. log2 FC is the female PAS response minus the male; there is no r column, because the contrast is not a two-group comparison.
- **Report:** volcanoes, counts per tier, a heatmap of every protein at p ≤ 0.05 across all four groups, then a grid of each one's values by mouse.
- **Result:** full proteome 156 at p ≤ 0.05 (92 expected), 163 at Π ≤ 0.05, none at FDR ≤ 0.05; mitochondrial universe 94, 69, 1 (Htra2, larger response in females). The Π lists are connected in STRING (up: 80 proteins, 264 edges, 2.4 times expected; down: 83 proteins, 92 edges, 1.4 times).
