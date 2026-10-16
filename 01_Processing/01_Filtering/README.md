# 01_Filtering

Turns DIA-NN precursors into a filtered, contaminant-free precursor matrix. No mouse is dropped and nothing is filtered on partial missingness, because limpa's detection probability curve (02_Quantification) models the missing values.

| | |
|---|---|
| Reads | `00_Setup/Input/report.parquet`, `00_Setup/Input/metadata.csv`, `00_Setup/Input/report_stats.xlsx` |
| Writes | `c_data/01_filtering.xlsx`, `c_data/precursors_filtered.rds` (read by 02_Quantification, 03_Design, 06_Modules, 02_Statistics/04_PCA, 02_Statistics/07_Sensitivity and `04_Figures/S1_Processing`), `b_reports/01_filtering.pdf` |
| Run | `Rscript 01_Processing/01_Filtering/a_script/01_filtering.R` |

- **Method:** `limpa::readDIANN()` with the q-value columns for a match-between-runs (MBR) search (`Q.Value`, `Lib.Q.Value`, `Lib.PG.Q.Value`) at 0.01; all-missing precursors removed; cRAP-only groups removed; `filterNonProteotypicPeptides()`, then `filterCompoundProteins()`; groups whose symbols all match keratin, immunoglobulin or haemoglobin removed. Singleton peptides are kept.
- **Result:** 26,133 precursors in 2,297 groups after cRAP removal, down to 21,600 precursors in 2,094 groups holding 74.2% of signal. Removing shared peptides accounts for most of the loss (98.3% to 75.0% of signal).
- **Basis:** limpa vignette, "Data input": the q-value set for MBR searches (after Nordmann) and the two peptide filters, which are optional and "help with downstream interpretation". MBR is inferred, not recorded: every `Q.Value` and `Lib.Q.Value` in the report is already below 0.01, so DIA-NN filtered on the library pass. To confirm with the facility.
- **Note:** limpa 1.4.0's `removeNARows()` fails on an EList once any row is all missing (it indexes with an extra empty subscript), so the script does the same subset itself. This removes one precursor, seen only in the dropped run `AI_PAS6`. The script's printed summary reports whether a cRAP group shares an accession with a mouse group (`crap_duplicates_proteome`) and which of 12 named muscle proteins (`Ckm`, `Cycs`, `Mb`, `Sod1` and others) are lost (`must_keep_missing`). Albumin, fibrinogen and complement are kept; `contaminant_share` reports their share per group. All 64 groups that `filterCompoundProteins()` removes list several UniProt accessions for a single gene, not a mixture of genes. Of those 63 genes, 55 (e.g. `Dysf`, `Bin1`, `Myom2`, `Utrn`, `Sptan1`, `Tln2`) are absent from the final 1,962-protein table; the other 8 survive through a single-accession group.
