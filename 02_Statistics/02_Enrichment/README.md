# 02_Enrichment

Runs every pathway test once, per contrast and per phenotype, so the comparisons and figures only slice the results.

| | |
|---|---|
| Reads | `01_Processing/02_Quantification/c_data/proteins.rds`, `01_Processing/03_Design/c_data/design.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `02_Statistics/01_Fit/c_data/fit.rds`, `00_Setup/Input/STRING_v12/` |
| Writes | `c_data/02_enrichment.xlsx`, `c_data/enrichment.rds`, `b_reports/02_enrichment.pdf` |
| Run | `Rscript 02_Statistics/02_Enrichment/a_script/02_enrichment.R` |

## Method

- Contrasts: `enrichVolcano::run_enrichment(tests = "fgsea")` on each contrast's moderated t; GO:BP, Reactome, Hallmark and KEGG LEGACY in the full proteome, MitoCarta in the mitochondrial universe; BH within contrast and collection.
- Each fgsea row carries two checks: the set's singscore contrast from 01_Fit, tested across mice, and `limma::camera`, which estimates the set's inter-protein correlation. fgsea permutes genes and runs optimistic; where they disagree, trust the checks.
- Shown sets: within each test, sets at FDR ≤ 0.05 are grouped by `representative_sets()` (Jaccard tree cut at 0.70); each group is shown by its most specific set within 10-fold of the group's best p. F03 and F04 draw these.
- Phenotypes: fgsea on each protein's association t, run and seeded once per sex, model and phenotype, as each contrast is.
- Also: GO Slim compartments and the 19 literature and ROS sets (`fgsea(nPermSimple = 10000)`, BH within each list; a compartment holding more than half the ranked genes is skipped, which drops organelle); cross-sex signatures (each sex's Π ≤ 0.05 lists tested in the other sex by `limma::fry`); ORA (`fgsea::fora`) of the S3 heatmap clusters, the sheet keeping sets at p ≤ 0.05; STRING connectivity of the Π lists.
- S3 clusters: each direction's proteins clustered on their own, average linkage on 1 − Pearson r over mice where both proteins were detected, cut by `dynamicTreeCut` (hybrid, deepSplit 2, minClusterSize 10).

## Results

- Sets at p ≤ 0.05 / FDR ≤ 0.05, GO:BP, Reactome, Hallmark, KEGG, MitoCarta: female 95/10, 24/0, 2/1, 6/0, 11/2; male 161/35, 16/0, 6/1, 7/5, 6/0; interaction 93/20, 19/6, 3/3, 8/2, 12/11.
- Sets shown after overlap grouping: 8 of 13 at FDR ≤ 0.05 (female), 14 of 41 (male), 24 of 42 (interaction).
- No cross-sex signature reaches p ≤ 0.05 in the other sex by fry (p 0.074 to 0.98), though 50 to 78% of each list moves the same way.
- Cluster ORA reaches FDR ≤ 0.05 in four clusters (male C3 and C4, interaction C2 and C3); male C3 and interaction C3 hold the muscle contraction and myogenesis sets. No female cluster has one. The clusters are built from responding proteins, so this ORA describes them and tests little.
