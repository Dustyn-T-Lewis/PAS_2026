# Male

Does PAS change the proteome and the phenotypes of males? Every number here comes from the `PAS_vs_VEH_male` contrast (PAS_M − VEH_M) of the shared fit. Run with `Rscript 03_Comparisons/Male/a_script/run_comparison.R`; the methods are in [`03_Comparisons`](../README.md).

![The male contrast and its four steps with their main results](schematic.png)

| Step | Result | Files |
|---|---|---|
| 01 Proteins | Full proteome: 193 proteins at p ≤ 0.05 (98 expected), 106 at Π ≤ 0.05, none at FDR ≤ 0.05. Mitochondrial universe: 37, 9, 0. Both Π lists are more connected in STRING than chance (up: 27 proteins, 5 edges, 5.0 times expected; down: 79 proteins, 203 edges, 3.9 times). | `01_proteins` |
| 02 Pathways | GO:BP 161 at p ≤ 0.05 (38 expected), 35 at FDR ≤ 0.05, 34 of them down: actin filament organisation and its regulation, cytoskeleton and supramolecular fibre organisation, cell and cell-matrix adhesion, plasma membrane organisation. 17 of the 35 have a singscore contrast at p ≤ 0.05, so the mouse-level test supports about half. Reactome 16 and 0; MitoCarta 6 and 0. | `02_pathways` |
| 03 Phenotypes | Phenotype effects: none at p ≤ 0.05; the closest is satellite cells per fibre, +43% (p 0.086). Phenotype associations at p ≤ 0.05, pooled and within treatment: proteins 756 and 706 of 11,772 tests (589 expected), 9 and 3 at FDR ≤ 0.05; pathway scores 138 and 82 of 5,364, fewer than the 268 expected; complexes 19 and 16 of 588; modules 3 and 1 of 84; clusters 1 and 0 of 42. With no phenotype effect, these associations follow variation between mice, not a PAS response. | `03_phenotypes` |
| 04 Classification | Pathway scores: 22 of 938 at p ≤ 0.05 (47 expected). 3 separate PAS from vehicle perfectly, against a median of 2 over all 252 relabellings (p 0.45). Modules 2 of 14, complexes 3 of 98, clusters 0 of 7. None at FDR ≤ 0.05. | `04_classification` |

Each step's workbook is `c_data/<step>.xlsx` and its preview `b_reports/<step>.pdf`.
