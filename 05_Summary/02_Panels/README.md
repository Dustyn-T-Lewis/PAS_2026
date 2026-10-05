# 02_Panels

Asks whether a combination of features separates the two groups of each comparison in mice the
model has not seen, for proteins, set scores, module eigengenes, all three combined, each with and
without the phenotypes, and the phenotypes alone.

- **Method:** `mixOmics::splsda` / `block.splsda` (DIABLO; 1 component, keepX fixed in advance: 20 proteins, 10 sets, 3 modules, 3 phenotypes; design weight 0.1), scored by `perf` leave-one-out, centroids distance. `nestedcv::nestcv.glmnet` (alpha 0.5, lambda.1se, inner leave-one-out), scored by leave-pair-out (25 outer folds, one mouse per group). Balanced error rate (BER). Null: all 126 five-and-five relabellings refitted the same way, so p is exact (smallest 1/126 = 0.008).
- **Inputs:** `02_Differential_Expression/02_Contrasts/c_data/contrasts.rds`, `03_Pathway_Enrichment/01_Scores/c_data/set_scores.rds`, `04_Network/01_Build_Modules/c_data/modules.rds`, `00_Input/phenotypes.xlsx`
- **Outputs:** `c_data/02_panels.xlsx` (sheets: performance, selected_features, plot_index, settings), `c_data/panels_cache.rds`, `b_reports/02_panels_<comparison>.pdf` (4)
- **Result:** PAS_vs_VEH_female: 6 of 18 panels at p ≤ 0.05; phenotypes alone BER 0.10 with both methods (p = 0.024), DIABLO on combined proteome levels BER 0.10 (p = 0.024), elastic net on modules plus phenotypes BER 0.12 (p = 0.016). Proteins alone do not reach p ≤ 0.05 (DIABLO BER 0.20, p = 0.10). PAS_vs_VEH_male: 1 of 18, DIABLO on modules BER 0.10 (p = 0.032). Sex in VEH and in PAS: all 36 panels at p ≤ 0.05, most with BER 0.
- **Note:** leave-pair-out keeps elastic-net training balanced at 4 against 4 (Airola et al. 2011); leave-one-out unbalances it and pushes the null BER towards 1. Many male elastic nets select nothing (`empty_models`). In `selected_features` the elastic-net weight is the share of the 25 folds that selected the feature. Results are cached per method and refit when inputs or settings change.
