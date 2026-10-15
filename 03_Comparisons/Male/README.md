# Male

Does PAS change the proteome and the phenotypes of males? Every number here comes from the `PAS_vs_VEH_male` contrast (PAS_M − VEH_M) of the shared fit. Run with `Rscript 03_Comparisons/Male/a_script/run_comparison.R`; the methods are in [`03_Comparisons`](../README.md).

![The male contrast and its four steps with their main results](schematic.png)

| Step | Result | Files |
|---|---|---|
| 01 Proteins | Full proteome: 133 proteins at p ≤ 0.05 (92 expected), 78 at Π ≤ 0.05, none at FDR ≤ 0.05. Mitochondrial universe: 80, 27, 0. Both Π lists are more connected in STRING than chance (up: 29 proteins, 9 edges, 2.3 times expected; down: 49 proteins, 103 edges, 4.9 times). | `01_proteins` |
| 02 Pathways | GO:BP 93 at p ≤ 0.05, 4 at FDR ≤ 0.05: lipid homeostasis, plasma membrane organisation, actin filament regulation and supramolecular fibre organisation down. Reactome 16 and 1 (striated muscle contraction down). MitoCarta 7 and 3 (central dogma, translation and mitochondrial ribosome up), but all three have singscore p ≥ 0.89, so the mouse-level test does not support them. | `02_pathways` |
| 03 Phenotypes | Phenotype effects: none at p ≤ 0.05; the closest is satellite cells per fibre, +43% (p 0.086). Phenotype associations at p ≤ 0.05: 2,185 of 11,052 protein tests (553 expected), 488 at FDR ≤ 0.05; pathway scores 167 of 5,076, fewer than the 254 expected; complexes 21 of 558; modules 4 of 108; clusters 5 of 78. With no phenotype effect, these associations follow variation between mice, not a PAS response. | `03_phenotypes` |
| 04 Classification | Pathway scores: 20 of 889 at p ≤ 0.05 (44 expected). 3 separate PAS from vehicle perfectly, against a median of 2 over all 252 relabellings (p 0.41). Modules 2 of 18, complexes 1 of 93, clusters 0 of 13. None at FDR ≤ 0.05. | `04_classification` |

Each step's workbook is `c_data/<step>.xlsx` and its report `b_reports/<step>.pdf`.
