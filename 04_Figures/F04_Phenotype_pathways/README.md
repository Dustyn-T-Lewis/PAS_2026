# F04_Phenotype_pathways

Figure 4: within each sex, which pathways track each phenotype across the ten mice.

![Figure 4](b_reports/F04_composite.png)

| | |
|---|---|
| Reads | `02_Statistics/02_Enrichment/c_data/enrichment.rds` |
| Writes | `b_reports/panels/a_trait_pathways.pdf/.png`, `b_reports/F04_composite.pdf/.png` (180 × 126 mm), `c_data/F04.xlsx` |
| Run | `Rscript 04_Figures/F04_Phenotype_pathways/a_script/F04_composite.R` (the panel script also runs alone) |

The figure is one panel, its own script in `a_script/panels/`; the composite runs it, adds the letter and writes one workbook. The panel reads saved outputs and computes nothing beyond filtering and ordering: its rows are the sets 02_Enrichment flags as `shown` under the within-treatment model. For each sex and trait, the sets at fgsea FDR ≤ 0.05 (both universes together) are grouped by an average-linkage Jaccard tree on their tested members, cut at 0.70, and each group is shown by its most specific set whose p is within 10-fold of the group's best; `representative_sets()` in `helpers.R` gives the rule and its sources (EnrichmentMap, DAVID, HRvLR). The pooled model is not drawn; both models are in the `association_fgsea` sheet of `02_enrichment.xlsx`.

## Panels

### a. Pathways per phenotype

![Panel a](b_reports/panels/a_trait_pathways.png)

Two grids side by side, female traits on the left and male traits on the right, sharing one key. One dot per gene set and trait: fgsea (`association_fgsea` in `enrichment.rds`) on each protein's association t with that trait across the sex's ten mice, under the within-treatment model. Fill is the NES, green when the set's proteins rise with the trait and magenta when they fall, on a scale set by the largest |NES| drawn (±2.7); dot size is −log10 FDR, capped at 4; the square beside each name gives its collection. Rows come in one block per trait: the sets that trait shows, by FDR, up to five, each listed once per sex, and every row is drawn against all six traits. Females fill all six blocks (30 rows); males show six sets for satellite cells, but four are already listed under CSA or capillary density, so that block has two rows (27 rows). Of 586 female set–trait pairs at FDR ≤ 0.05, 220 are shown and 30 drawn; for males, 165 of 410 are shown and 27 drawn. Of the drawn cells, 96 of 180 in females and 58 of 162 in males reach FDR ≤ 0.05.

In females, complex I biogenesis, respiratory chain assembly and mitochondrion organisation rise with CSA (NES 1.9 to 2.3), and aerobic respiration and OXPHOS rise with satellite cells (2.0 to 2.2) and with CSA. Muscle contraction sets fall with capillary density (NES −2.4 to −2.7, FDR down to 1.9e-12) and rise with satellite cells and both pyruvate conductances. Extracellular matrix, EMT and cardiomyopathy sets fall with soleus octanoylcarnitine conductance (NES −1.9 to −2.3), and focal adhesion, ECM-receptor and cardiomyopathy sets rise with plantaris pyruvate conductance (1.8 to 2.1). The sets selected for soleus pyruvate conductance (haemopoiesis, vesicle transport, membrane organisation, blood vessel morphogenesis, immune regulation) fall with it (NES −1.8 to −2.1) and move little with the other traits.

In males, OXPHOS rises with CSA (NES 2.1) and more strongly with capillaries (2.5), while UV response down, the mitochondrial ribosome and RNA metabolism fall with CSA. The citric acid cycle and ion transport rise with capillaries (2.3 and 2.1). The sets selected for soleus pyruvate conductance (fatty acid metabolism, mitochondrial protein degradation, TP53-regulated metabolic genes, cardiac muscle contraction, organic acid catabolism) fall with it (NES −1.8 to −2.0), and OXPHOS falls with it too (−2.2); several of the same sets rise with capillaries and with soleus octanoylcarnitine conductance. A falling OXPHOS signal with rising respiratory conductance runs against expectation; respirometry was done on soleus and the proteome comes from gastrocnemius, so treat this pairing as an unexplained association. Cell morphogenesis, extracellular matrix, endocytosis, myogenesis and EMT sets rise with plantaris pyruvate conductance (NES 1.8 to 2.1).

These associations are co-movement across mice, not PAS effects. The gene-permutation FDR runs optimistic when genes are correlated (Wu & Smyth 2012), so two checks sit beside every set in the `association_fgsea` sheet: camera (`camera_p`, `camera_fdr`), which estimates the set's own inter-protein correlation, and the set's singscore slope on the trait, tested across mice (`score_p`, `score_fdr`). For the 30 female rows at their selecting trait, camera gives p ≤ 0.05 for 4 and FDR ≤ 0.05 for none. The singscore slope gives p ≤ 0.05 for 7 and FDR ≤ 0.05 for 2, dilated cardiomyopathy (0.031) and tight junction (0.050), both falling with soleus octanoylcarnitine conductance. For the 27 male rows, camera gives p ≤ 0.05 for 4 and the singscore slope for 4, with neither reaching FDR ≤ 0.05. Read the grids as rankings of which biology moves with each trait, not as confirmed set-level associations.

## Workbook

`c_data/F04.xlsx` holds one sheet, `a_trait_pathways`, with every drawn cell (sex, trait block, set, collection, universe, phenotype, NES, fgsea p and FDR, whether that phenotype shows the set, camera p and FDR), plus `read_me`, `input_manifest` and `package_versions`.
