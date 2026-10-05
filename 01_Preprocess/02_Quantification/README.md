# 02_Quantification

Rolls precursors up to proteins with limpa's detection probability curve, so each protein gets an
abundance, a standard error and a precursor count per sample. Adds MitoCarta annotation and
records how mitochondrial enrichment and sample medians differ before any normalisation.

- **Method:** `limpa::dpcQuant(y, "Protein.Group", dpc.slope = 0.7)`, then `filterByDetection(n.samples = 5)` (the smallest group size). `dpcCN()` fits slope 0.455 and `dpcON(robust = TRUE)` 0.344, both below the 0.7 to 0.9 limpa's FAQ calls typical for DIA-NN, so the preset 0.7 is used; the fitted 0.455 is quantified beside it for comparison and no reported result uses it. MitoCarta 3.0 matched on gene symbol or UniProt accession, because symbol alone misses renamed genes (`Atp5a1` is now `Atp5f1a`). Mito share compared by Welch `t.test()`.
- **Inputs:** `01_Preprocess/01_Filtering/c_data/precursors_filtered.rds`, `00_Input/Mouse.MitoCarta3.0.xls`
- **Outputs:** `c_data/02_quantify.xlsx` (sheets: dpc_parameters, mito_share, enrichment_tests, normalization_evidence), `c_data/proteins.rds` (read by every later stage), `c_data/dpcQuant/proteins_slope_0.7.rds` and `proteins_slope_0.455.rds` (checkpoints), `b_reports/02_quantify_figures.pdf` (detection curve, missingness against abundance, mito share by group, abundance boxplots, `plotMDSUsingSEs()`)
- **Result:** 1,842 proteins by 20 samples, 568 in MitoCarta. MitoCarta proteins carry a median 47% of signal in females and 58 to 59% in males: sex p = 7.5e-09, treatment p = 0.78. Sample medians span 16.5 to 17.3 log2 before normalisation.
- **Note:** `02_quantify.R` sources `02_quantify_run.R` and loads the checkpoints, so it runs in seconds. Rebuild them (about 2.5 minutes per slope) with `Rscript 01_Preprocess/02_Quantification/a_script/02_quantify_run.R` or `REFIT_DPCQUANT=TRUE` when the precursor matrix changes; a checkpoint from a different filtering run stops the script. MitoCarta is an annotation column, not a filter. The sex difference in enrichment means every sex contrast partly reflects pellet composition.
