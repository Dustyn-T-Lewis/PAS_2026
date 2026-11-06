# Male

Contrast: `PAS_vs_VEH_male`, PAS_M − VEH_M. Run with `Rscript 03_Comparisons/Male/a_script/run_comparison.R`.

![The male contrast and its four steps with their main results](schematic.png)

## Proteins

- Full proteome: 193 at p ≤ 0.05, 106 at Π ≤ 0.05, none at FDR ≤ 0.05.
- Mitochondrial universe: 37, 9 and 0.
- 8 of the 106 Π ≤ 0.05 proteins were never detected in one male group (6 in PAS males, 2 in vehicle males), so that group's value is limpa's estimate; `detected_PAS_M` and `detected_VEH_M` filter them.
- Both Π lists have more STRING edges than chance: up 5 edges among 27 proteins (5.0 times expected), down 203 among 79 (3.9 times).

## Pathways

- GO:BP 161 sets at p ≤ 0.05 and 35 at FDR ≤ 0.05, 34 of them down: actin and cytoskeleton organisation, cell and cell-matrix adhesion, plasma membrane organisation. 17 of the 35 have singscore p ≤ 0.05.
- KEGG 7 and 5, all down (ECM-receptor interaction, cardiac muscle contraction, three cardiomyopathy or myocarditis sets). Hallmark 6 and 1 (EMT down). Reactome 16 and 0. MitoCarta 6 and 0.
- These falls match the lower vessel, nerve, red-cell and myofibril scores of male PAS pellets (`02_Statistics/05_Composition`).

## Phenotypes

- No phenotype effect at p ≤ 0.05; the closest is satellite cells per fibre, +43% (p 0.086).
- Protein associations at p ≤ 0.05, pooled and within treatment: 756 and 706 of 11,772 tests; 9 and 3 at FDR ≤ 0.05. With no phenotype effect, these follow variation between mice, not a PAS response.

## Classification

- 28 of 1,005 pathway scores at p ≤ 0.05. 5 separate perfectly, against a median of 2 over the 252 relabellings (p 0.33).
- One of 8 modules at p ≤ 0.05, none perfect; 3 of 98 complexes at p ≤ 0.05.
- Nothing reaches FDR ≤ 0.05.
