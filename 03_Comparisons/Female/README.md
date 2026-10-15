# Female

Does PAS change the proteome and the phenotypes of females? Every number here comes from the `PAS_vs_VEH_female` contrast (PAS_F − VEH_F) of the shared fit. Run with `Rscript 03_Comparisons/Female/a_script/run_comparison.R`; the methods are in [`03_Comparisons`](../README.md).

![The female contrast and its four steps with their main results](schematic.png)

| Step | Result | Files |
|---|---|---|
| 01 Proteins | Full proteome: 152 proteins at p ≤ 0.05 (92 expected by chance), 109 at Π ≤ 0.05, 4 at FDR ≤ 0.05 (Dusp3 and Cand2 down; Arl6ip5 and Bcap31 up). Mitochondrial universe: 47, 18, 0. Both Π lists are more connected in STRING than chance (up: 49 proteins, 65 edges, 2.3 times expected; down: 60 proteins, 65 edges, 3.4 times). | `01_proteins` |
| 02 Pathways | GO:BP 102 at p ≤ 0.05 (36 expected), 7 at FDR ≤ 0.05: muscle system process, muscle contraction, myofibril assembly, sarcomere organisation, multicellular organismal movement and carbohydrate catabolism up, cell cycle down. Reactome 17 and 0; MitoCarta 5 and 0. Of the 7, only cell cycle has a singscore contrast at p ≤ 0.05 (p 0.039). | `02_pathways` |
| 03 Phenotypes | Phenotype effects: capillary density +138% (p 0.005), CSA +22% (p 0.002), satellite cells per fibre +47% (p 0.009); the three respiratory conductances do not change (p ≥ 0.30). Phenotype associations at p ≤ 0.05: 803 of 11,052 protein tests (553 expected), 20 at FDR ≤ 0.05, 14 of them with satellite cells; pathway scores 248 of 5,076 (254 expected), none at FDR; complexes 44 of 558, one at FDR; modules 7 of 108; clusters 5 of 78. | `03_phenotypes` |
| 04 Classification | Pathway scores: 47 of 889 at p ≤ 0.05 (44 expected). 17 separate PAS from vehicle perfectly, against a median of 3 over all 252 relabellings (p 0.095). Modules 1 of 18, complexes 6 of 93, clusters 1 of 13. None at FDR ≤ 0.05. | `04_classification` |

Each step's workbook is `c_data/<step>.xlsx` and its report `b_reports/<step>.pdf`.
