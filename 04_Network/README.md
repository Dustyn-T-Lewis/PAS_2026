# 04 · Network

Co-abundance modules from WGCNA on the full universe (1,842 proteins): which module eigengenes move
with treatment or sex, which alone separate two groups, which track a phenotype, what each module
contains, and whether the proteins PAS moved are linked in STRING or the INDRA literature graph.
Every test is reported nominally with its statistics so results can be filtered; the working rule is
p ≤ 0.05 and FDR ≤ 0.05. Modules move with sex; none moves with PAS.

| Sub-stage | Does | Writes |
|---|---|---|
| `01_Build_Modules` | signed bicor WGCNA: 9 modules plus grey, eigengenes, kME, hubs | `modules.rds`, workbook, PDF |
| `02_Contrasts` | eigengene contrasts (emmeans) and module-trait correlations | `module_contrasts.rds`, workbook, 6 PDFs |
| `03_Classify` | AUC of each eigengene alone for the four two-group comparisons | workbook, PDF |
| `04_Associate` | eigengene against each phenotype under three models; DGCA | workbook, PDF |
| `05_Characterise` | gene-set over-representation, hubs and STRING enrichment per module | `string_cache.rds`, workbook, PDF |
| `06_Mechanism` | STRING and INDRA on the female FDR list and the Π lists | `mechanism_cache.rds`, workbook, PDF |

Every workbook opens with a `read_me` sheet (what each sheet holds) and ends with `input_manifest`
(files read, with md5) and `package_versions`.

Run order (`run_all.R`): `01`, `02`, `05`, `06`; then `04` after `02_Differential_Expression/04_Associate`,
and `03` after `02_Differential_Expression/03_Classify`. `03` and `04` only need `01`'s `modules.rds`;
the later slot keeps each level's classify and associate steps together. `05` and `06` share
STRING's downloaded network files in `stringdb_cache/` (local, untracked).
