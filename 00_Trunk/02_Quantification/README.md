# 02_Quantification

Rolls precursors up to proteins with limpa's detection probability curve, so every protein has a value and a standard error in every mouse. Adds MitoCarta and normalises the full proteome and the mitochondrial universe separately.

| | |
|---|---|
| Reads | `01_Filtering/c_data/precursors_filtered.rds`, `00_Input/Mouse.MitoCarta3.0.xls` |
| Writes | `c_data/02_quantification.xlsx`, `c_data/proteins.rds` (read by 03 to 07), `c_data/dpcQuant/proteins_slope_0.7.rds` (checkpoint), `b_reports/02_quantification.pdf` |
| Run | `Rscript 00_Trunk/02_Quantification/a_script/02_quantification.R` |

- **Method:** `limpa::dpcQuant(dpc.slope = 0.7)`, then `filterByDetection(n.samples = 5)`. `dpcCN()` fits slope 0.455 and robust `dpcON()` 0.344, below the 0.7 to 0.9 limpa's FAQ calls typical for DIA-NN, so the preset 0.7 is used. MitoCarta 3.0 matched on symbol or UniProt accession (symbol alone misses renamed genes such as `Atp5f1a`). Quantile normalisation (`normalizeBetweenArrays`) per universe.
- **Result:** 1,842 proteins in 20 mice, 568 in MitoCarta. MitoCarta proteins carry 46 to 53% of signal in females and more in males (Welch sex p = 7.5e-09; treatment p = 0.78).
- **Note:** The script loads the checkpoint and runs in seconds. Rebuild it (about 2.5 min) with `Rscript 00_Trunk/02_Quantification/a_script/02_quantification_run.R` or `REFIT_DPCQUANT=TRUE`; a checkpoint from another filtering run stops the script. Because pellet composition differs by sex, a male-female comparison partly reflects mitochondrial enrichment; within-sex PAS contrasts do not, since PAS and vehicle pellets do not differ.
