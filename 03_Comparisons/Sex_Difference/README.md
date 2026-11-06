# Sex_Difference

Contrast: `treatment_by_sex`, (PAS_F − VEH_F) − (PAS_M − VEH_M). A positive value means PAS raised a protein or set more in females than in males. A change in one sex only is not, on that basis, a sex difference; only this contrast tests it. Run with `Rscript 03_Comparisons/Sex_Difference/a_script/run_comparison.R`.

![The interaction contrast and its two steps with their main results](schematic.png)

## Proteins

- Full proteome: 162 at p ≤ 0.05, 142 at Π ≤ 0.05, none at FDR ≤ 0.05.
- Mitochondrial universe: 46, 26 and 0.
- 40 of the 142 Π ≤ 0.05 proteins were never detected in at least one group, 33 of them in PAS males, so that group's value is limpa's estimate; the `detected_*` columns filter them.
- Both Π lists have more STRING edges than chance: higher in females 294 edges among 97 proteins (3.2 times expected), lower 23 among 45 (3.3 times).

## Pathways

- GO:BP 93 sets at p ≤ 0.05 and 20 at FDR ≤ 0.05. 19 are higher in the female response (muscle contraction, myofibril assembly, detoxification, redox homeostasis); proteasome-mediated catabolism is lower. 7 of the 20 have singscore p ≤ 0.05.
- Reactome 19 and 6, all higher in females (muscle contraction, collagen and extracellular matrix, response to chemical stress); 2 have singscore p ≤ 0.05. Hallmark 3 and 3; KEGG 8 and 2.
- MitoCarta 12 and 11: the broad metabolism set, ROS and glutathione, detoxification, lysine, amino acid, lipid and carbohydrate metabolism higher; complex I, the mitochondrial ribosome and translation lower. 3 have singscore p ≤ 0.05.

## Phenotypes

No phenotype step runs here. The phenotype interaction is in `02_Statistics/01_Fit` (sheet `phenotype_interaction`): only capillary density differs (p 0.014; the other five p 0.22 to 0.95).

## Classification

Not run. Classification needs a label each mouse carries, and the interaction is a difference of group differences.
