# 01_Filtering

Turns DIA-NN precursors into a filtered, contaminant-free precursor matrix. No mouse is dropped and nothing is filtered on partial missingness, because limpa's detection curve (02) models those gaps.

| | |
|---|---|
| Reads | `00_Input/report.parquet`, `00_Input/metadata.csv`, `00_Input/report_stats.xlsx` |
| Writes | `c_data/01_filtering.xlsx`, `c_data/precursors_filtered.rds` (read by 02 and 10), `b_reports/01_filtering.pdf` |
| Run | `Rscript 00_Trunk/01_Filtering/a_script/01_filtering.R` |

- **Method:** `limpa::readDIANN()` with five q-value columns (`Q.Value`, `Lib.Q.Value`, `Lib.PG.Q.Value`, `Global.Q.Value`, `Global.PG.Q.Value`) at 0.01; all-missing precursors removed; cRAP-only groups removed; `filterNonProteotypicPeptides()`, then `filterCompoundProteins()`; groups whose symbols all match keratin, immunoglobulin or haemoglobin removed. Singleton peptides are kept.
- **Result:** 24,634 precursors in 2,104 groups after cRAP removal, down to 20,355 precursors in 1,923 groups holding 74.2% of signal. Most of the loss is shared peptides (98.3% to 74.9% of signal).
- **Note:** The script stops if a cRAP group shares an accession with a mouse group, or if any of 12 named muscle proteins (`Ckm`, `Cycs`, `Mb`, `Sod1` and others) is lost. Albumin, fibrinogen and complement are kept; `contaminant_share` reports their share per group.
