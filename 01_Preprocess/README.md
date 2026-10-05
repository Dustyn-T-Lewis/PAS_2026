# 01 · Preprocess

How much of each protein each mouse has, with an uncertainty attached. The stage turns the
DIA-NN precursor report into a protein table of 1,842 proteins by 20 samples, 568 of them in
MitoCarta 3.0. Every later stage reads `02_Quantification/c_data/proteins.rds`; keep it, not
abundances alone, because the statistics downstream use its standard errors.

| Sub-stage | Does | Writes |
|---|---|---|
| `01_Filtering` | q-value, proteotypic, compound-group, cRAP and contaminant filters; run QC | `precursors_filtered.rds`, `01_filter.xlsx`, `01_filter_figures.pdf` |
| `02_Quantification` | `dpcQuant()` roll-up to proteins, detection filter, MitoCarta annotation, enrichment and normalisation evidence | `proteins.rds`, `dpcQuant/` checkpoints, `02_quantify.xlsx`, `02_quantify_figures.pdf` |

Each sub-stage holds `a_script/` (the script), `b_reports/` (PDFs) and `c_data/` (workbook and
`.rds`). The two workbooks here hold only their result sheets: they have no `read_me`,
`input_manifest` or `package_versions` sheet, unlike the workbooks written through
`R/helpers.R::write_workbook()`.

Run order, from the repository root:

```sh
Rscript 01_Preprocess/01_Filtering/a_script/01_filter.R
Rscript 01_Preprocess/02_Quantification/a_script/02_quantify.R
```

Nothing is normalised in this stage. Quantile normalisation runs per protein universe in
`02_Differential_Expression/02_Contrasts`.
