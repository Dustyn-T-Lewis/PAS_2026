# 03_Comparisons

The three questions, answered by slicing the `02_Statistics` results. Nothing here refits a model, and no comparison reads another.

![One workflow of four steps, run for the female, male and sex-difference comparisons into their own folders](schematic.png)

| Folder | Holds |
|---|---|
| [`Female`](Female/README.md) | PAS_F − VEH_F (`PAS_vs_VEH_female`), steps 01 to 04 |
| [`Male`](Male/README.md) | PAS_M − VEH_M (`PAS_vs_VEH_male`), steps 01 to 04 |
| [`Sex_Difference`](Sex_Difference/README.md) | (PAS_F − VEH_F) − (PAS_M − VEH_M) (`treatment_by_sex`), steps 01 and 02 |
| `code/` | the four steps, written once and called by `<comparison>/a_script/run_comparison.R` |

| Script in `code/` | Writes |
|---|---|
| `01_proteins.R` | the contrast for every protein, tiered FDR ≤ 0.05, else Π ≤ 0.05, else p ≤ 0.05; STRING connectivity of the Π lists |
| `02_pathways.R` | the contrast's fgsea rows, each beside the set's singscore contrast |
| `03_phenotypes.R` | phenotype effects, and phenotype associations under the pooled and within-treatment models |
| `04_classification.R` | the sex's single-score ROC rows and the perfect-separator test |

Each step writes `c_data/<step>.xlsx` and a preview `b_reports/<step>.pdf` into its comparison folder.

## Run

```sh
Rscript 03_Comparisons/Female/a_script/run_comparison.R
```

The same for `Male` and `Sex_Difference`. Each needs only `02_Statistics` steps 01 to 03.
