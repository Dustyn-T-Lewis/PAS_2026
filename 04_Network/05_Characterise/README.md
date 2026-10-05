# 05_Characterise

Describes what each module is: which gene sets it is enriched for, its hub proteins, and whether
its members interact in STRING more than chance.

- **Method:** `fgsea::fora` per module and collection (GO:BP, Hallmark, KEGG, Reactome; detected proteome as background, set size 15 to 500, BH within collection); STRINGdb v12.0 (mouse, score ≥ 400) `get_ppi_enrichment` with the detected proteome as background, BH across modules.
- **Inputs:** `../01_Build_Modules/c_data/modules.rds`, `03_Pathway_Enrichment/00_Gene_Sets/c_data/gene_sets.rds`, `02_Differential_Expression/02_Contrasts/c_data/study_full/da_results.csv`
- **Outputs:** `c_data/05_characterise.xlsx` (sheets: overview, annotation, hubs, string, plot_index), `c_data/string_cache.rds`, `b_reports/05_characterise.pdf`
- **Result:** six modules have a set at FDR ≤ 0.05: black is contractile (striated muscle contraction, FDR 1e-13; STRING 10.6-fold), yellow OXPHOS (Hallmark, FDR 2e-10; 1.5-fold), turquoise translation and mitochondrial translation (FDR 2e-7), red vesicle transport (FDR 5e-5; 2.8-fold), magenta platelet response and haemostasis (FDR 2e-6; 3.1-fold), blue carbohydrate biosynthesis (FDR 0.013). Brown, green and pink have nominal top terms only (green amino-acid catabolism, p = 0.006). Eight of nine modules are STRING-enriched at FDR ≤ 0.05; pink is not (1.17-fold, FDR 0.36).
- **Note:** STRING answers are cached and refetched only when `modules.rds` changes (md5).
