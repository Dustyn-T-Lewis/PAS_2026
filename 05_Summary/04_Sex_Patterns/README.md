# 04_Sex_Patterns

Labels each protein, pathway and module response to PAS as shared, detected in one sex, or
sex-differential. Only the direct interaction test claims the sexes differ. Refits nothing.

- **Method:** `sex_label()` in `R/helpers.R` from the female and male PAS effects and the treatment_by_sex test: sex-differential (interaction p ≤ 0.05); shared (both sexes p ≤ 0.05, same sign); detected in females or in males (one sex p ≤ 0.05); opposite (both p ≤ 0.05, opposite sign, interaction not significant); neither. Every row carries all three FDRs.
- **Inputs:** `02_contrasts.xlsx` from `02_Differential_Expression`, `03_Pathway_Enrichment` (fgsea sheet) and `04_Network`
- **Outputs:** `c_data/04_sex_patterns.xlsx` (sheets: summary, patterns, plot_index), `b_reports/04_sex_patterns.pdf` (female against male effect scatter per level)
- **Result:** proteins, full universe: 156 sex-differential (92 expected by chance), 7 shared, 75 detected in females, 80 in males, 1,524 neither. Mito universe: 94, 2, 14, 29. Pathways (full): GO:BP 92 sex-differential, Reactome 23, Hallmark 4, KEGG 1; MitoCarta 9. Modules: all 9 neither. No response is labelled opposite.
- **Note:** labels use nominal p; filter on the FDR columns for a stricter list. "Detected in one sex" does not mean the sexes differ.
