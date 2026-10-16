# F03_Pathways

Figure 3: which pathways each comparison moves.

![Figure 3](b_reports/F03_composite.png)

| | |
|---|---|
| Reads | `02_Statistics/02_Enrichment/c_data/enrichment.rds` |
| Writes | `b_reports/panels/<panel>.pdf/.png`, `b_reports/F03_composite.pdf/.png` (180 × 150 mm), `c_data/F03.xlsx` |
| Run | `Rscript 04_Figures/F03_Pathways/a_script/F03_composite.R` (each panel script also runs alone) |

Each panel is its own script in `a_script/panels/`; the composite runs them, adds the letters, stitches them and checks that the rows are exactly the sets at FDR ≤ 0.05 that collapsePathways keeps, each drawn in all three contrasts. Panels read saved outputs and compute nothing beyond counts, filtering and scaling.

## Panels

| Panel | Shows | Source |
|---|---|---|
| a | GO:BP, Reactome, Hallmark and KEGG legacy sets in the full proteome at fgsea FDR ≤ 0.05 in any contrast that collapsePathways keeps, listed once under the first contrast that selects them and drawn in all three. Dot fill: NES; size: −log10 FDR; square: collection. | `enrichment.rds`: `results` |
| b | MitoCarta pathways in the mitochondrial universe, selected and drawn the same way. | `enrichment.rds`: `results` |

fgsea permutes genes, so its FDR runs optimistic when genes are correlated (Wu & Smyth 2012). Two checks sit beside every set in the 02_Enrichment workbook: camera (`camera_p`, `camera_fdr`), which estimates the set's own inter-protein correlation, and the set's singscore contrast tested across mice. Of the 42 kept sets at fgsea FDR ≤ 0.05 in the three contrasts, camera gives p ≤ 0.05 for 7 and FDR ≤ 0.05 for none. The 7 include proteasome-mediated catabolism, lower in female PAS (camera p 0.001) and in the interaction (p 0.008), and MitoCarta ROS and glutathione metabolism in the interaction (p 0.009). Read the rest as rankings, not confirmed set shifts.

## Workbook

`c_data/F03.xlsx` holds one sheet per panel with the numbers it draws, plus `read_me`, `input_manifest` and `package_versions`.
