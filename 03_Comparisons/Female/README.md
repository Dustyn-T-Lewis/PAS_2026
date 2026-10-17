# Female

Does PAS change the proteome and the phenotypes of females? Every number here comes from the `PAS_vs_VEH_female` contrast (PAS_F − VEH_F) of the shared fit. Run with `Rscript 03_Comparisons/Female/a_script/run_comparison.R`; the methods are in [`03_Comparisons`](../README.md).

![The female contrast and its four steps with their main results](schematic.png)

| Step | Result | Files |
|---|---|---|
| 01 Proteins | Full proteome: 177 proteins at p ≤ 0.05, 107 at Π ≤ 0.05, 3 at FDR ≤ 0.05 (Cops3, Cops8 and Dusp3, all down). Mitochondrial universe: 42, 17, 0. Both Π lists are more connected in STRING than chance (up: 44 proteins, 26 edges, 2.6 times expected; down: 63 proteins, 60 edges, 2.4 times). | `01_proteins` |
| 02 Pathways | GO:BP 95 at p ≤ 0.05, 10 at FDR ≤ 0.05: muscle system process, muscle, skeletal and striated muscle contraction, neuromuscular process, multicellular organismal movement and small-molecule catabolism up; cell cycle, mitotic cell cycle and proteasome-mediated protein catabolism down. Only the three down sets have a singscore contrast at p ≤ 0.05. Reactome 24 and 0. Hallmark 2 and 1 (myogenesis up). KEGG 6 and 0. MitoCarta 11 and 2 (metabolism and carbohydrate metabolism up, singscore p 0.20 and 0.69). | `02_pathways` |
| 03 Phenotypes | Phenotype effects: capillary density +138% (p 0.005), CSA +22% (p 0.002), satellite cells per fibre +47% (p 0.009); the three respiratory conductances do not change (p ≥ 0.30). Phenotype associations at p ≤ 0.05, pooled and within treatment: proteins 764 and 779 of 11,772 tests, 5 and 2 at FDR ≤ 0.05; full-proteome pathway scores 383 and 274 of 5,766, 0 and 16 at FDR; complexes 71 and 46 of 588; modules 2 and 1 of 48; clusters 2 of 42 in both. | `03_phenotypes` |
| 04 Classification | Pathway scores: 47 of 1,005 at p ≤ 0.05. 16 separate PAS from vehicle perfectly, against a median of 4 over all 252 relabellings (p 0.13). Modules 1 of 8 at p ≤ 0.05, and it separates perfectly (null median 0, p 0.063); complexes 6 of 98, 1 perfect (p 0.35); clusters 0 of 7. No level has more perfect separators than the relabelling null at p ≤ 0.05, and none reaches FDR ≤ 0.05. | `04_classification` |

Each step's workbook is `c_data/<step>.xlsx` and its preview `b_reports/<step>.pdf`.
