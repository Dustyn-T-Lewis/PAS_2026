# 02_Quantification

Rolls precursors up to proteins with limpa's detection probability curve, so every protein has a value and a standard error in every mouse. Adds MitoCarta annotation and normalises the full proteome and the mitochondrial universe separately.

| | |
|---|---|
| Reads | `01_Processing/01_Filtering/c_data/precursors_filtered.rds`, `00_Setup/Input/Mouse.MitoCarta3.0.xls` |
| Writes | `c_data/02_quantification.xlsx`, `c_data/proteins.rds`, `c_data/dpcquant_slope_0.7.rds` (checkpoint), `b_reports/02_quantification.pdf` |
| Run | `Rscript 01_Processing/02_Quantification/a_script/02_quantification.R`; loads the checkpoint. `REFIT_DPCQUANT=TRUE` refits (about 2.5 min), as does any change to the filtered precursors. |

## Method

- `limpa::dpcQuant(dpc.slope = 0.7)`, then `filterByDetection(n.samples = 5)`. The fitted slopes (`dpcCN()` 0.476, `dpcON()` 0.349) fall below the 0.7 to 1 limpa calls typical, so the preset 0.7 is used.
- MitoCarta 3.0 matched on symbol or UniProt accession, since symbol alone misses renamed genes such as `Atp5f1a`.
- Normalisation: none, median centring, quantile and cyclic loess are scored on within-group spread (PMAD, PEV, PCV from NormalyzerDE). Only none and median centring are eligible, because they keep limpa's standard errors valid. The eligible method with the lower PMAD is used. No contrast enters the choice.

## Results

- 1,962 proteins in 20 mice, 595 in MitoCarta.
- Full proteome: no normalisation (PMAD 0.409 against 0.442 for median centring). MitoCarta: median centring (0.387 against 0.388).
- MitoCarta proteins carry 45 to 52% of signal in female pellets and 52 to 61% in male pellets (Welch sex p = 1.9e-08; treatment p = 0.83). A male against female comparison partly reflects mitochondrial enrichment; a within-sex PAS contrast does not.
