# Input

Study data and reference files. The pipeline reads them from here; only STRINGdb contacts a server at run time.

| File | One row is | Read by |
|---|---|---|
| `report.parquet` | one precursor in one MS run | `01_Processing/01_Filtering`, `02_Statistics/05_Composition` |
| `metadata.csv` | one MS run (21 rows) | `01_Processing/01_Filtering` |
| `report_stats.xlsx` | one MS run, DIA-NN's QC table | `01_Processing/01_Filtering` |
| `Mouse.MitoCarta3.0.xls` | one mouse gene (1,140 mitochondrial) | `01_Processing/02_Quantification`, `01_Processing/04_Gene_Sets` |
| `corum_complete.txt` | one CORUM 5.3 complex (release 2026-04-14), all organisms | `01_Processing/04_Gene_Sets` |
| `STRING_v12/` | STRING v12 mouse aliases, protein info and links at score ≥ 400 | `01_Processing/04_Gene_Sets`, `02_Statistics/02_Enrichment`, `02_Statistics/06_Candidates` |
| `phenotypes.xlsx` | one mouse, 41 measures, plus `dictionary` and `notes` sheets | `02_Statistics/01_Fit` |
| `VEH_PAS_HLI_Master Data.xlsx` | one measure, one column per mouse, as delivered | nothing; source of `phenotypes.xlsx` |

## Proteomics

`report.parquet` is the delivered DIA-NN main report `Mus.report.parquet`, renamed: 297,290 rows, 26,800 precursors, 21 runs, decoys removed. `Precursor.Quantity` has no zeros or NAs, so a missing precursor is an absent row, which is what limpa expects. The delivered wide matrix (`Mus.report.pr_matrix.tsv`, no q-values) and MaxLFQ protein table (`AI_Proteins_20samples.xlsx`) are not used and not kept.

Twenty mice gave 21 runs. Run `AI_PAS6` identified 703 protein groups against a median of 1,430 for the other runs and was re-injected as `AI_PAS6r` (994). `metadata.csv` sets `include = FALSE` on the original and `TRUE` on the re-run.

Sex is recorded only in `metadata.csv`, entered by hand from the sample sheet: `VEH1-5` and `PAS1-5` are female, `VEH6-10` and `PAS6-10` male. `sample_id` carries it, so `AI_PAS6r` is `PAS_M1`.

## Phenotypes

`phenotypes.xlsx` is the master sheet turned to one row per mouse by hand, values unchanged, keyed on `mouse`. It holds three histology measures (CSA, capillary density, satellite cells per fibre), 32 Oroboros oxygen-consumption readings, four respiratory conductances (soleus and plantaris, pyruvate and octanoylcarnitine) and two AoxBC values. Six are analysed: the three histology measures and three conductances. `dictionary` gives each column's meaning and assumed unit; `notes` lists gaps and the two mistyped mouse numbers in the source.

What the delivered data do not say:

- The proteome comes from gastrocnemius; respirometry was run on soleus and plantaris.
- The muscle used for histology is not stated, nor the satellite-cell marker.
- CSA has no unit; µm² is assumed.
- Capillary density has no definition (per fibre or per area).
- Respiratory conductance has no unit; it is taken as the slope of JO2 against ΔG_ATP.
- JO2 normalisation is not stated; wet weight is assumed.
- AoxBC is not defined, and it was measured in nine of the ten males and no females.
- The sheet does not say whether the measured limb is the ligated one.
