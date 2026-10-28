# F03_Pathways

Figure 3: which pathways each comparison moves.

![Figure 3](b_reports/F03_composite.png)

| | |
|---|---|
| Reads | `02_Statistics/02_Enrichment/c_data/enrichment.rds` |
| Writes | `b_reports/panels/a_pathways.pdf/.png`, `b_reports/F03_composite.pdf/.png` (212 × 108 mm), `c_data/F03.xlsx` |
| Run | `Rscript 04_Figures/F03_Pathways/a_script/F03_composite.R` (the panel script also runs alone) |

The figure is one panel, its own script in `a_script/panels/`; the composite runs it and writes one workbook, with no panel letter. The panel reads saved outputs and computes nothing beyond filtering, ordering and counting: its rows are the sets 02_Enrichment flags as `shown`, and each row is drawn in all three contrasts and at every trait in both sexes, the same columns as F04, and bars to the right count each contrast's distinct significant sets. Within one contrast, the sets at fgsea FDR ≤ 0.05 are grouped by an average-linkage Jaccard tree on their tested members, cut at 0.70, and each group is shown by its most specific set whose p is within 10-fold of the group's best; `representative_sets()` in `helpers.R` gives the rule and its sources (EnrichmentMap, DAVID, HRvLR).

## Panels

### a. Pathways per contrast

![Panel a](b_reports/panels/a_pathways.png)

One row per gene set. The columns come in three blocks: Contrast (PAS against vehicle in females, the same in males, and the interaction, Treatment × sex, the female response minus the male), then Female traits and Male traits, the within-treatment association fgsea at each of the six traits in that sex, as in F04. Only the contrast columns choose the rows; the trait columns show how the same sets track the phenotypes. Each dot's fill is the fgsea normalised enrichment score, magenta for sets whose proteins fall and green for sets that rise, on a scale set by the largest |NES| drawn (±2.7); its area grows with −log10 FDR from a small minimum, so a cell near FDR 1 still shows its colour, and stops growing at FDR 0.0001; the key reads in FDR (0.1, 0.05, 0.01, ≤ 0.0001). The square before each row gives its collection: GO:BP teal, Reactome indigo, Hallmark sand, KEGG legacy light blue, MitoCarta ochre (Paul Tol's colours, away from the direction and sex colours). The rows sit in four blocks. Female, Male and Treatment × sex hold the full-proteome sets that contrast shows, by FDR then p, up to 6 each; MitoCarta holds the MitoCarta pathways (mitochondrial universe, 595 proteins) shown in any contrast, by their smallest p, up to 6. A set is listed once, in the first block that reaches it. That gives 24 rows, 6 in each block. A thin dark outline marks the cell that put each row in its block ("Chosen in this contrast" in the key): the block's own contrast in the first three blocks, and in the MitoCarta block every contrast that shows the set (female and the interaction for metabolism and carbohydrate metabolism, the interaction alone for the other four). Each row sits on a light band in its block's colour, the same as the block's strip, carried across all 15 columns.

The bars to the right count distinct significant sets: sets at fgsea FDR ≤ 0.05 that 02_Enrichment flags as `shown`, so each group of overlapping sets counts once. Each contrast block has an up bar (▲) above a down bar (▼), split by database inside one outline, with the total at its end; the MitoCarta block has one pair per contrast. Females have 4 up and 2 down, males 1 up and 13 down, and Treatment × sex 15 up and 1 down. In MitoCarta, females have 2 up, the interaction 6 up and 2 down, and males none. The bars count every distinct set, so they run past the six rows a block draws; the counts before overlap grouping are in the workbook.

Females show six full-proteome sets, four up (myogenesis at FDR 0.00029, muscle system process, skeletal muscle contraction, small molecule catabolism) and two down (mitotic cell cycle, proteasome-mediated catabolism). Males show 14, and the top 6 drawn are all down: supramolecular fibre organisation (FDR 1.5e-05), cell adhesion, EMT, actin filament bundle organisation, ECM-receptor interaction and viral myocarditis. The interaction shows 16; three are already listed in the blocks above (EMT, proteasome-mediated catabolism, viral myocarditis), and the top 6 of the other 13 are drawn, all higher in the female response and led by striated muscle contraction (FDR 9.2e-05).

Read across the columns and the sexes mostly disagree: female and male NES differ in sign in 14 of the 18 full-proteome rows, and the interaction takes the female sign in 17. Because the interaction subtracts the male response, a green interaction dot can come from a male fall as much as a female rise, and all 6 male rows pair a male fall with a small female NES (|NES| ≤ 1.22). One male FDR-level set, actin filament based process (male FDR 0.00014), sits under Treatment × sex because another set represents its group in males. The male falls in cytoskeletal, adhesion and extracellular matrix sets line up with the lower vessel, nerve, red-cell and myofibril content of male PAS pellets (S1), so they may reflect pellet composition as much as a change in the fibres.

In the MitoCarta block, females gain the broad metabolism set (NES 1.65) and carbohydrate metabolism (1.97), both at FDR 0.0033 and both also shown in the interaction. The other four rows come from the interaction alone: ROS and glutathione metabolism (NES 2.02, FDR 0.0024), lysine metabolism (1.86, 0.023), complex I subunits (−1.70, 0.031) and the mitochondrial ribosome (−1.60, 0.048). Amino acid and lipid metabolism (both FDR 0.050 in the interaction) fall below the cap. No MitoCarta set reaches FDR ≤ 0.05 in males. The female and male NES differ in sign in all six rows, so each interaction result combines a move one way in females with a move the other way in males. For the four interaction-only sets, neither sex reaches FDR ≤ 0.05 alone (female 0.10 to 0.30, male 0.12 to 0.34).

A shift in the pellet's MitoCarta share between PAS and vehicle was not detected within either sex (S1 e), but several other marker scores do differ (S1 f: sarcoplasmic reticulum in females; vessel, nerve, red cell and myofibril in males), so pellet composition cannot be ruled out as a contributor to the within-sex results here. ROS and glutathione metabolism in the interaction is one of the three rows camera supports at p ≤ 0.05 (p 0.013), and 02_Enrichment notes that several FDR-level MitoCarta terms have singscore p > 0.2.

## Notes

fgsea permutes genes, so its FDR runs optimistic when genes are correlated (Wu & Smyth 2012). Two checks sit beside every set in the 02_Enrichment workbook: camera (`camera_p`, `camera_fdr`), which estimates the set's own inter-protein correlation, and the set's singscore contrast tested across mice. Of the 24 rows, camera gives p ≤ 0.05 in the contrast that shows the set for 3 and FDR ≤ 0.05 for none. The 3 are proteasome-mediated catabolism, lower in female PAS (camera p 0.003), and nervous system process (p 0.043) and MitoCarta ROS and glutathione metabolism (p 0.013) in the interaction. Read the rest as rankings, not confirmed set shifts.

## Workbook

`c_data/F03.xlsx` holds `a_pathways`, with every drawn cell, 360 in all (row block, set, collection, universe, column block, contrast or sex and phenotype, NES, fgsea p, FDR, whether that column shows the set, whether the cell is outlined (`selected`), camera p and FDR); `counts`, the sets at FDR ≤ 0.05 per block, contrast, direction and collection, distinct (`n_distinct`, what the bars draw) and all (`n_all`); plus `read_me`, `input_manifest` and `package_versions`.
