# 00 · Input

Study data. Nothing here runs.

| File | What one row is | Used by |
|---|---|---|
| `report.parquet` | one precursor in one MS run | `01_Filtering` |
| `metadata.csv` | one MS run, 21 rows | `01_Filtering` |
| `report_stats.xlsx` | one MS run, DIA-NN's own QC table | `01_Filtering` |
| `Mouse.MitoCarta3.0.xls` | one mouse gene, 1,140 mitochondrial | `02_Quantification`, `03_Pathway_Enrichment/00_Gene_Sets` |
| `phenotypes.xlsx` | one mouse, 41 physiology measures, with a `dictionary` and `notes` sheet | the `04_Associate` steps and `05_Summary/02_Panels`, through `R/helpers.R` |
| `VEH_PAS_HLI_Master Data.xlsx` | one measure per row, one mouse per column, as delivered | nothing; the source of `phenotypes.xlsx` |

`report.parquet` is the delivered `Mus.report.parquet`, renamed. This DIA-NN main report holds 297,290 precursor-by-run rows
(26,800 precursors, 21 runs), decoys stripped, label-free single channel. `Precursor.Quantity`
has no zeros or NAs, so a missing precursor is an absent row, as limpa's detection-probability
model wants. `01_Filtering` checks and removes the `cRAP-` contaminants the search appended.

Two delivered files are unused and not kept here. `Mus.report.pr_matrix.tsv` has the same
precursors in wide form but no q-value columns, so limpa could not filter it.
`AI_Proteins_20samples.xlsx` is DIA-NN's MaxLFQ protein table, and a pre-summarised matrix
discards the precursor-level modelling that is limpa's point.

Twenty mice gave twenty-one runs. `AI_PAS6` returned 1,055 protein groups against a cohort
median near 1,700 and was re-injected as `AI_PAS6r`; both are in the report. `metadata.csv` sets
`include = FALSE` on the original and `TRUE` on the re-run, with the same `mouse` on both rows,
so swapping the decision means flipping two flags.

Sex is recorded only in `metadata.csv`, entered by hand from the sample sheet: `VEH1-5` and
`PAS1-5` are female, `VEH6-10` and `PAS6-10` male. `sample_id` is rewritten to carry it, so
`AI_PAS6r` becomes `PAS_M1`.

`phenotypes.xlsx` is the delivered master sheet turned to one row per mouse by hand, values
unchanged, keyed on `mouse` so it joins `metadata.csv` and every proteomic table. It holds three
histology measures, 36 Oroboros oxygen-consumption readings and four respiratory conductances
from soleus and plantaris on two fuels, and two AoxBC values measured in nine of the ten males.
The `dictionary` sheet says what each column is and which units are assumed; the `notes` sheet
records the gaps, the two mistyped mouse numbers in the source and what is still to confirm with
Dr. Ismaeel, chief among them that respirometry was run on soleus and plantaris while the
proteome comes from gastrocnemius.
