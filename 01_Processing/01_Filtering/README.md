# 01_Filtering

Turns DIA-NN precursors into a filtered, contaminant-free precursor matrix. No mouse is dropped and nothing is filtered on missingness, because limpa models missing values in 02_Quantification.

| | |
|---|---|
| Reads | `00_Setup/Input/report.parquet`, `00_Setup/Input/metadata.csv`, `00_Setup/Input/report_stats.xlsx` |
| Writes | `c_data/01_filtering.xlsx`, `c_data/precursors_filtered.rds`, `b_reports/01_filtering.pdf` |
| Run | `Rscript 01_Processing/01_Filtering/a_script/01_filtering.R` |

## Method

- `limpa::readDIANN()` at q ≤ 0.01 on `Q.Value`, `Lib.Q.Value` and `Lib.PG.Q.Value`, the set the limpa vignette gives for match-between-runs searches. The report does not record whether MBR was used; every `Q.Value` and `Lib.Q.Value` in it is already below 0.01.
- Removed in order: all-missing precursors, cRAP-only groups, shared peptides (`filterNonProteotypicPeptides()`), compound groups (`filterCompoundProteins()`), and groups whose symbols are all keratin, immunoglobulin or haemoglobin.
- Singleton peptides are kept. Albumin, fibrinogen and complement are kept; `contaminant_share` gives their share per group.

## Results

- 26,310 precursors in 2,309 groups pass the q-value filters; 21,600 in 2,094 groups remain.
- Removing shared peptides costs the most: retained signal falls from 98.3% to 75.0%, because muscle's abundant isoform families share peptides. Final retained signal is 74.2%.
- The 64 compound groups removed each list several accessions for one gene, not a mix of genes.
