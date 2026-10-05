# 02 · Differential Expression

Which single proteins moved with each contrast, which alone separate two groups, and which track
the six phenotypes. Proteins are tested in two universes, all 1,842 and the 568 MitoCarta
proteins, with BH correction within each contrast and universe. A protein counts as significant at
p ≤ 0.05 and FDR ≤ 0.05; protein results are also reported at FDR ≤ 0.10. Π = p^|log2FC| ranks
proteins and controls no error rate.

| Sub-stage | Does | Writes |
|---|---|---|
| `01_Design` | cell-means design, five contrasts, protein universes, sign checks | `design.rds`, `01_design.xlsx` |
| `02_Contrasts` | quantile normalisation per universe, `dpcDE()` fit, moderated t per contrast | `contrasts.rds`, `02_contrasts.xlsx`, `study_full/`, `study_mito/`, one PDF per contrast |
| `03_Classify` | AUC and exact Wilcoxon p per protein for four two-group comparisons | `03_classify.xlsx`, `03_classify.pdf` |
| `04_Associate` | protein against phenotype slopes under three models, DGCA correlation classes | `protein_association.rds`, `04_associate.xlsx`, one PDF per phenotype |

Contrasts: `PAS_vs_VEH_female`, `PAS_vs_VEH_male`, `treatment_by_sex` (female PAS effect minus
male PAS effect), `sex_in_VEH`, `sex_in_PAS` (both male minus female). The interaction equals
`sex_in_VEH - sex_in_PAS`. The first three are primary, the sex contrasts supplementary.

Each sub-stage holds `a_script/`, `b_reports/` and `c_data/`. The workbooks of `02_Contrasts`,
`03_Classify` and `04_Associate` open with a `read_me` sheet and end with `input_manifest` (inputs
with md5) and `package_versions`; `plot_index` maps results to PDF pages. `01_design.xlsx` has none
of these sheets.

Run order, from the repository root. `03_Classify` and `04_Associate` read only `02_Contrasts`
output, so either can follow it; the repository README runs them after the pathway and network
stages.

```sh
Rscript 02_Differential_Expression/01_Design/a_script/01_design.R
Rscript 02_Differential_Expression/02_Contrasts/a_script/02_contrasts.R
Rscript 02_Differential_Expression/04_Associate/a_script/04_associate.R
Rscript 02_Differential_Expression/03_Classify/a_script/03_classify.R
```
