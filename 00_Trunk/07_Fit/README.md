# 02_Contrasts

Tests which proteins each contrast moved, in the full proteome and in the MitoCarta universe.
Normalising and fitting each universe separately lets the mito results borrow information only
from mitochondrial proteins.

- **Method:** `limma::normalizeBetweenArrays(method = "quantile")` per universe; `limpa::dpcDE(sample.weights = TRUE)` (precision weights from the quantification standard errors, plus sample weights); `contrasts.fit()` and `eBayes()`; BH within each contrast and universe; Π beside p and FDR.
- **Inputs:** `01_Preprocess/02_Quantification/c_data/proteins.rds`, `02_Differential_Expression/01_Design/c_data/design.rds`
- **Outputs:** `c_data/02_contrasts.xlsx` (sheets: summary, results, pi_lists, plot_index), `c_data/contrasts.rds` (normalised matrices, fits, weights, results, Π lists; read by `03_Classify`, `04_Associate` and most later stages), `c_data/study_full/` and `c_data/study_mito/` (enrichVolcano study tables read by `03_Pathway_Enrichment` and `04_Network`), `b_reports/02_contrasts_<contrast>.pdf` (volcano and p-value histogram, one page per universe)
- **Result (full, 1,842 tests, 92 expected at p ≤ 0.05):** Female PAS: 152 at p ≤ 0.05, 4 at FDR ≤ 0.05, 8 at FDR ≤ 0.10. Male PAS: 133, none at FDR ≤ 0.10. Interaction: 156, none at FDR ≤ 0.10. Sex in VEH: 720, 524 at FDR ≤ 0.05. Sex in PAS: 966, 835 at FDR ≤ 0.05.
- **Result (mito, 568 tests, 28 expected):** Female PAS: 47, 3 at FDR ≤ 0.10. Male PAS: 80, none. Interaction: 94, 1 at FDR ≤ 0.05, 3 at FDR ≤ 0.10. Sex in VEH: 233, 173 at FDR ≤ 0.05. Sex in PAS: 334, 288 at FDR ≤ 0.05.
- **Note:** The script stops unless the full-universe FDR ≤ 0.10 counts are 8 (female PAS), 661 (sex in VEH) and 972 (sex in PAS). `pi_lists` (Π ≤ 0.05, full universe: 109 female, 78 male, 163 interaction) feed `04_Network/06_Mechanism`. To filter, use `results` with `universe`, `contrast` and `adj.P.Val`.
