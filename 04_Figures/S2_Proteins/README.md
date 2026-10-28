# S2_Proteins

Supplementary Figure 2: each comparison's own protein list, cut into co-expression clusters, and the sets each cluster over-represents.

![Supplementary Figure 2](b_reports/S2_composite.png)

| | |
|---|---|
| Reads | `02_Statistics/01_Fit/c_data/fit.rds`, `01_Processing/02_Quantification/c_data/proteins.rds`, `02_Statistics/02_Enrichment/c_data/enrichment.rds` |
| Writes | `b_reports/panels/<panel>.pdf/.png`, `b_reports/S2_composite.pdf/.png` (180 × 242 mm), `c_data/S2.xlsx` |
| Run | `Rscript 04_Figures/S2_Proteins/a_script/S2_composite.R` (each panel script also runs alone) |

Each panel is its own script in `a_script/panels/`; the composite runs them, adds the letters, stitches them over one shared key and writes one workbook. Panels read saved outputs and compute nothing beyond counts, filtering and scaling.

## Panels

### a. Female PAS against vehicle

![Panel a](b_reports/panels/a_pi_female.png)

Each row is one of the 107 full-proteome proteins at Π ≤ 0.05 in the female PAS against vehicle contrast (`fit.rds`: `protein_contrasts`), and each column one of the 20 mice, grouped as Female VEH, Female PAS, Male VEH and Male PAS under the group colour bar. A cell is the protein's normalised log2 abundance (`E_norm`) minus its mean over the vehicle mice of the same sex, so green is above own-sex vehicle and magenta below, with the colour full at ±1.5 (the key's end breaks); a grey cell marks a mouse in which the protein was not detected (no precursor observed). Rows are split by the sign of the female fold change, 44 up (log2 FC 0.55 to 1.81) and 63 down (−2.08 to −0.57), and each block is cut into co-expression clusters in 02_Enrichment (`enrichment.rds`: `heatmap_clusters`). The tree is average-linkage clustering on 1 − Pearson r across the 20 mice of the drawn values, each pair of proteins over the mice where both were detected (`coexpression_tree()` in `helpers.R`), and `dynamicTreeCut::cutreeDynamic(method = "hybrid", deepSplit = 2, minClusterSize = 10)` cuts it. A pair sharing fewer than 5 such mice counts as uncorrelated (17 of 946 pairs in the up block and 121 of 1,953 in the down block here; 26 of 351 and 54 of 3,081 in panel b; 119 of 4,656 and 94 of 990 in panel c); a protein detected in fewer than 5 mice would stay off the tree, and none does. Clusters are numbered down the panel and labelled with their size, and rows run in tree order within each. A protein the cut leaves unassigned is drawn at the end of its block, under the block's last cluster label but outside its count and its ORA (here Siglec1 in the up block, Tubb2b and Rps15 in the down block). Up is one cluster, C1 (43); down splits into C2 (51) and C3 (10). S6 a orders its rows on the same tree (`coexpression_order()`), with no cut.

Beside each cluster are the three sets with the smallest p (ties broken by FDR) among those its proteins over-represent at p ≤ 0.05 against the detected proteome (`enrichment.rds`: `cluster_ora`, `fgsea::fora` on GO:BP, Reactome, Hallmark and KEGG, BH within cluster and collection), fewer where fewer qualify, and "no set at p ≤ 0.05" where none does. Overlapping sets are not merged, so near-duplicates such as Reactome and GO:BP muscle contraction can both appear. A cluster's bars sit in a frame centred on it, moved up or down only where the frame would overlap the next one or run past its block; bar length is −log10 p on one axis (0 to 10.2) shared by panels a to c, so a bar of a given length is the same p in every panel. Fill gives the database, a black outline marks FDR ≤ 0.05 and the dotted tick p = 0.05. Each label gives the set's proteins in the cluster over its size, in sentence case with a few long words shortened (Reg., Dev., Proc., Org., Pos., Neg., Dep., Loc., Biosynth.), on at most two lines of 24 characters; a name two lines cannot hold ends at its last whole word with an ellipsis (four here; the workbook has them in full), and where two sets of one cluster share a name a tag gives the database (GO, R, H or K). No female set reaches FDR ≤ 0.05. C1: ER to Golgi vesicle-mediated transport (3/17, p = 0.005), adaptive immune system (5/67, 0.014), regulation of homotypic cell-cell adhesion (3/26, 0.018). C2: Reactome cell cycle (6/47, p = 0.001, FDR 0.17), proteasome-mediated ubiquitin-dependent protein catabolic process (4/22, 0.002), modification-dependent macromolecule catabolic process (4/30, 0.007). C3 has none at p ≤ 0.05.

The clusters are cut from proteins chosen for responding, and ORA on a cluster of 10 to 51 proteins has little power, with BH run within that cluster only, so read the bars as a description of each cluster and not as a test.

The female PAS block separates by construction, since the female contrast chose the rows; the male columns are the informative part, and only 6 of the 107 are also on the male list. Of the 107, 20 have no precursor observed in any male run, so their male cells are grey.

### b. Male PAS against vehicle

![Panel b](b_reports/panels/b_pi_male.png)

The same layout for the 106 proteins at Π ≤ 0.05 in the male PAS against vehicle contrast. Down outnumbers up, 79 (log2 FC −2.95 to −0.67) against 27 (0.60 to 3.14), and the male PAS block carries that pattern because the male contrast chose the rows. Up splits into C1 (12) and C2 (15), down into C3 (26) and C4 (46), with 7 down proteins unassigned (Srsf5, Srrm2, Nufip2, Ddah1, Eif4b, Faf1, Vcp). C3 holds the muscle sets, all three at FDR ≤ 0.05: Reactome striated muscle contraction (7/21, FDR 5e-7), GO:BP muscle contraction (10/80, 2e-5) and Reactome muscle contraction (8/43, 2e-6), out of 24 C3 sets at FDR ≤ 0.05. C4's three also pass: sensory system development (6/24, 0.007), plasma membrane organization (5/19, 0.012) and sensory organ development (6/31, 0.012), out of 5. The up clusters' sets, C1 cell redox homeostasis (2/15), regulation of neuron projection development (2/21) and cellular oxidant detoxification (2/25), and C2 muscle cell differentiation (3/58), muscle structure development (3/86) and muscle cell development (2/35), have p between 0.003 and 0.029 and FDR 1. Both up frames are taller than their 12 and 15 rows, so C1's frame starts level with the group colour bar.

The female columns are the informative part; of the 6 proteins shared with the female list, 3 fall in both sexes, 2 rise in females and fall in males, and 1 does the reverse. None of the 106 reaches FDR ≤ 0.05; the smallest FDR is 0.094. With five mice per group, a protein at Π ≤ 0.05 is a ranked candidate for follow-up and not a confirmed change.

The solid green stripe in one Female VEH and one Female PAS column, here and in panel c, is two mice, VEH_F5 and PAS_F2: 14 of the male list's proteins (fast-muscle sarcomere proteins such as Myh4, Tnnt3, Tnnc2 and Mybpc2, plus Acta1 and Col1a1) sit more than 1.5 log2 above the median of both mice's group-mates (20 of the 142 in panel c). Their limpa sample weights are 0.49 (VEH_F5, lowest of the females) and 0.92 (PAS_F2), their run depths 13,377 and 18,352 precursors, their haemoglobin shares 0.42% and 0.34%, and they score highest in their groups on the 05_Composition myofibril marker set.

### c. Treatment × sex

![Panel c](b_reports/panels/c_pi_interaction.png)

The same layout for the 142 proteins at Π ≤ 0.05 in the treatment × sex contrast, the female PAS effect minus the male PAS effect. Rows are split by the sign of that difference: F > M (97 proteins, log2 FC 0.75 to 3.60) where the female effect is the more positive, and F < M (45, −2.76 to −0.91) where it is the more negative. F > M splits into C1 (14), C2 (30) and C3 (50), with 3 unassigned (Siglec1, Abcg2, Sec11a), and F < M into C4 (24) and C5 (21). C3 holds the same muscle sets as the male C3: Hallmark myogenesis (16/90, FDR 6e-9), Reactome striated muscle contraction (9/21, 8e-8) and GO:BP muscle contraction (13/80, 3e-5), out of 36 at FDR ≤ 0.05. C2 has Hallmark EMT (6/49, FDR 0.002), endothelium development (3/15, p = 0.001, FDR 0.54) and Reactome extracellular matrix organization (5/57, p = 0.001, FDR 0.10); C1 has no set at p ≤ 0.05. The F < M sets, C4 microtubule-based movement (3/33), establishment of organelle localization (3/40) and regulation of lipid biosynthetic process (2/22), and C5 cell-cell signaling (4/83), KEGG neurotrophin signaling (2/15) and KEGG oocyte meiosis (2/15), have p between 0.007 and 0.029 and FDR 0.18 to 1.

A positive difference can come from a rise in females, a fall in males or both, and the heatmap shows which: many F > M rows are green in female PAS and magenta in male PAS. Of the 142, 60 are also on the female list and 45 on the male list, and 15 have no precursor observed in any male run, so their male cells are grey and the difference rests on the females alone. No protein reaches FDR ≤ 0.05 (smallest 0.21).

### d. Volcano plots

![Panel d](b_reports/panels/d_volcanoes.png)

One point per full-proteome protein, 1,962 per contrast, with log2 fold change on the x axis and −log10 p on the y axis; facets are female PAS against vehicle, male PAS against vehicle and the female minus male difference. Points at Π ≤ 0.05 are green when up and magenta when down, the rest grey; the dashed line marks p = 0.05, and a black ring and the gene name mark FDR ≤ 0.05 (BH within contrast). Π is a curved boundary, so a fold change above 2 can pass it with p above 0.05 (7 female, 3 male and 25 interaction proteins do) and a small fold change can fail it with p below 0.05.

At p ≤ 0.05 the contrasts give 177 (female), 193 (male) and 162 (interaction) proteins. Only three proteins pass FDR ≤ 0.05, all female and all down: Cops3 (log2 FC −1.25, FDR 0.022), Dusp3 (−1.39, 0.022) and Cops8 (−2.08, 0.025).

Π is p^|log2 FC| (Xiao et al. 2014); it ranks on p and fold change together and controls no error rate. The heatmaps' rows are chosen on the contrast they show, so the selected sex's PAS columns will always look separated; the other sex's columns and the volcano's FDR rings are the fairer read.

## Workbook

`c_data/S2.xlsx` holds the numbers each panel draws: for a to c, one sheet of the heatmap's proteins in drawn order with their direction and cluster (`<panel>_heatmap`, cluster empty where unassigned) and one of its ORA bars (`<panel>_ora`, the drawn `cluster_ora` rows as stored); for d, one sheet. `read_me`, `input_manifest` and `package_versions` follow.
