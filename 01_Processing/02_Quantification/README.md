# 02_Quantification

Rolls precursors up to proteins with limpa's detection probability curve, so every protein has a value and a standard error in every mouse. Adds MitoCarta and normalises the full proteome and the mitochondrial universe separately, by the method a replicate-spread benchmark picks.

| | |
|---|---|
| Reads | `01_Processing/01_Filtering/c_data/precursors_filtered.rds`, `00_Setup/Input/Mouse.MitoCarta3.0.xls` |
| Writes | `c_data/02_quantification.xlsx`, `c_data/proteins.rds` (read by 03_Design to 06_Modules and 02_Statistics/01_Fit), `c_data/dpcquant_slope_0.7.rds` (checkpoint), `b_reports/02_quantification.pdf` |
| Run | `Rscript 01_Processing/02_Quantification/a_script/02_quantification.R` |

- **Method:** `limpa::dpcQuant(dpc.slope = 0.7)`, then `filterByDetection(n.samples = 5)`. `dpcCN()` fits slope 0.476 and robust `dpcON()` 0.349, below the 0.7 to 1 limpa calls representative of real MS data, so the preset 0.7 is used. MitoCarta 3.0 matched on symbol or UniProt accession (symbol alone misses renamed genes such as `Atp5f1a`).
- **Normalisation:** four candidates per universe (none, median centring, quantile, cyclic loess), scored on replicate spread within the four groups: PMAD, PEV and PCV, the NormalyzerDE measures. Only none and median centring are eligible, because they shift a mouse's values by one constant and so leave limpa's standard errors valid. The eligible method with the lowest PMAD is used; the rule was written before the first run and no contrast enters it. Sheet `normalization` holds every candidate's scores.
- **Result:** 1,962 proteins in 20 mice, 595 in MitoCarta. Normalisation: none for the full proteome (PMAD 0.409 against 0.442 for median centring; PEV and PCV slightly favour median centring), median centring for MitoCarta (PMAD 0.387 against 0.388). MitoCarta proteins carry 45 to 52% of signal in females and 52 to 61% in males (Welch sex p = 1.9e-08; treatment p = 0.83).
- **Basis:** limpa vignette, "The detection probability curve" (slopes of 0.7 to 1 are typical, and the estimate "may be under-estimated if the variability of each peptide is large") and "Differential expression" (dpcDE weights by the quantification standard errors); Willforss et al. 2019 (NormalyzerDE) for the spread measures; limma's `normalizeMedianValues` divides by medians, which suits unlogged intensities, so median centring is written out.
- **Note:** The script loads the checkpoint and runs in seconds. Refit it (about 2.5 min) with `REFIT_DPCQUANT=TRUE Rscript 01_Processing/02_Quantification/a_script/02_quantification.R`; a checkpoint from another filtering run stops the script. Because pellet composition differs by sex, a male-female comparison partly reflects mitochondrial enrichment; within-sex PAS contrasts do not, since PAS and vehicle pellets do not differ.
