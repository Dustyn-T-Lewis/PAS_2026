# Sex_Difference

Does the PAS response differ between the sexes? The test is the interaction contrast `treatment_by_sex`, (PAS_F − VEH_F) − (PAS_M − VEH_M); a positive value means PAS raised a protein or set more in females than in males. A feature that changes in one sex and not the other is not, on that basis, a sex difference; only this contrast tests it. Run with `Rscript 03_Comparisons/Sex_Difference/a_script/run_comparison.R`; the methods are in [`03_Comparisons`](../README.md).

![The interaction contrast and its two steps with their main results](schematic.png)

| Step | Result | Files |
|---|---|---|
| 01 Proteins | Full proteome: 156 at p ≤ 0.05 (92 expected), 163 at Π ≤ 0.05, none at FDR ≤ 0.05. Mitochondrial universe: 94, 69, 1 (Htra2, a larger response in females). The Π lists are connected in STRING (up: 80 proteins, 264 edges, 2.4 times expected; down: 83 proteins, 92 edges, 1.4 times). | `01_proteins` |
| 02 Pathways | GO:BP 92 at p ≤ 0.05, 0 at FDR ≤ 0.05. Reactome 19 and 11: higher in the female response are muscle contraction, striated muscle contraction, extracellular matrix organisation, ECM degradation, ECM proteoglycans, non-integrin membrane–ECM interactions, collagen formation, scavenger receptors and cellular response to chemical stress; lower are cell cycle and mitotic cell cycle. Seven of the 11 have singscore p ≤ 0.05. MitoCarta 9 and 6: central dogma, mitochondrial ribosome and translation lower in the female response, metabolism, carbohydrate and lipid metabolism higher; none has singscore p below 0.14, so the mouse-level test does not support them. | `02_pathways` |

There are no phenotype or classification steps. The phenotype interaction (capillary density p 0.014; the other five p ≥ 0.21) is in `02_Statistics/01_Fit` (sheet `phenotype_interaction`) and in Figure 1 panel e. Classification needs a label each mouse carries, and the interaction is a difference of differences between groups.

Each step's workbook is `c_data/<step>.xlsx` and its report `b_reports/<step>.pdf`.
