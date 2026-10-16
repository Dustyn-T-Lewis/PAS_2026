# Input

Study data and the reference files the pipeline reads, so it runs without a network connection. Nothing here runs.

| File | What one row is | Used by |
|---|---|---|
| `report.parquet` | one precursor in one MS run | `01_Processing/01_Filtering`, `02_Statistics/05_Composition` (haemoglobin) |
| `metadata.csv` | one MS run, 21 rows | `01_Processing/01_Filtering` |
| `report_stats.xlsx` | one MS run, DIA-NN's own QC table | `01_Processing/01_Filtering` |
| `Mouse.MitoCarta3.0.xls` | one mouse gene, 1,140 mitochondrial | `01_Processing/02_Quantification`, `01_Processing/04_Gene_Sets` |
| `corum_complete.txt` | one CORUM 5.3 complex (release 2026-04-14), every organism | `01_Processing/04_Gene_Sets` |
| `STRING_v12/` | STRING v12 mouse aliases, protein info and links at score ≥ 400, as STRINGdb downloads them | `01_Processing/04_Gene_Sets`, `02_Statistics/02_Enrichment`, `02_Statistics/06_Candidates` |
| `phenotypes.xlsx` | one mouse, 41 physiology measures, with a `dictionary` and `notes` sheet | `02_Statistics/01_Fit`; six are analysed |
| `VEH_PAS_HLI_Master Data.xlsx` | one measure per row, one mouse per column, as delivered | nothing; the source of `phenotypes.xlsx` |

`report.parquet` is the delivered `Mus.report.parquet`, renamed. This DIA-NN main report holds 297,290 precursor-by-run rows (26,800 precursors, 21 runs), decoys stripped, label-free single channel. `Precursor.Quantity` has no zeros or NAs, so a missing precursor is an absent row, as limpa's detection-probability model wants. `01_Processing/01_Filtering` checks and removes the `cRAP-` contaminants the search appended.

Two delivered files are unused and not kept here. `Mus.report.pr_matrix.tsv` has the same precursors in wide form but no q-value columns, so limpa could not filter it. `AI_Proteins_20samples.xlsx` is DIA-NN's MaxLFQ protein table, and a pre-summarised matrix discards the precursor-level modelling that is limpa's point.

Twenty mice gave twenty-one runs. `AI_PAS6` identified 703 protein groups against a median of 1,430 for the included runs (`report_stats.xlsx`, `Proteins.Identified`) and was re-injected as `AI_PAS6r` (994); both are in the report. `metadata.csv` sets `include = FALSE` on the original and `TRUE` on the re-run, with the same `mouse` on both rows, so swapping the decision means flipping two flags.

Sex is recorded only in `metadata.csv`, entered by hand from the sample sheet: `VEH1-5` and `PAS1-5` are female, `VEH6-10` and `PAS6-10` male. `sample_id` is rewritten to carry it, so `AI_PAS6r` becomes `PAS_M1`.

`phenotypes.xlsx` is the delivered master sheet turned to one row per mouse by hand, values unchanged, keyed on `mouse` so it joins `metadata.csv` and every proteomic table. It holds three histology measures, 32 Oroboros oxygen-consumption readings, four respiratory conductances from soleus and plantaris on two fuels, and two AoxBC values measured in nine of the ten males. The `dictionary` sheet says what each column is and which units are assumed. The `notes` sheet records the gaps, the two mistyped mouse numbers in the source, and what is still to confirm with Dr. Ismaeel, including that respirometry was run on soleus and plantaris while the proteome comes from gastrocnemius.
