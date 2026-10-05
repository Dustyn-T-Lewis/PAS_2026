# 01_Signal

One table of every primary test family, so a reader can see how many nominal hits exceed chance
before reading any single hit. Refits nothing.

- **Method:** per family (level, question, universe, term): tests, p ≤ 0.05, chance (0.05 × tests), excess, Π ≤ 0.05 (proteins), FDR ≤ 0.10 and ≤ 0.05, and share of true effects as 1 − `limma::propTrueNull(p)` for families of 100 tests or more.
- **Inputs:** `02_contrasts.xlsx` and `04_associate.xlsx` from `02_Differential_Expression`, `03_Pathway_Enrichment` and `04_Network` (proteins, pathways by collection, modules)
- **Outputs:** `c_data/01_signal.xlsx` (sheets: signal, plot_index; 166 families), `b_reports/01_signal.pdf`
- **Result:** full universe, PAS_vs_VEH_female 152 at p ≤ 0.05 against 92 by chance, 4 at FDR ≤ 0.05, estimated 8.5% true effects; PAS_vs_VEH_male 133, none at FDR, 3.4%; treatment_by_sex 156, none at FDR, 10.8%. Sex contrasts carry most of the signal: 720 (VEH) and 966 (PAS) at p ≤ 0.05, 524 and 835 at FDR ≤ 0.05, 48% and 60% true effects. Protein-phenotype families: CSA pooled slope 466 at p ≤ 0.05, 169 at FDR ≤ 0.05.
- **Note:** the counts use strict p ≤ 0.05 and FDR ≤ 0.05 (`tier_counts()` in `R/helpers.R`).
