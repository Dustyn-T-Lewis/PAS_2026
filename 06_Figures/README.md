# 06 · Figures

Manuscript figures, drawn only from saved outputs of stages 01 to 05. Landscape, 297 mm wide.

```
F0x/
  a_script/main/F0x.R         stitcher: runs the panel scripts, assembles the figure
  a_script/main/panels/*.R    one script per panel; each also runs alone
  a_script/supp/S0y.R         stitcher for each supplementary figure behind F0x
  a_script/supp/panels/*.R
  b_reports/main/             F0x.pdf, F0x.png, panels/*.pdf
  b_reports/supp/             S0y.pdf, S0y.png, panels/*.pdf
  c_data/                     F0x_data.xlsx, S0y_data.xlsx: one sheet per panel, read_me first
```

| Figure | Shows | Supplements |
|---|---|---|
| `F01` | global proteome: PCA with PERMANOVA and MitoCarta share, sex differences, proteins and pathways PAS moves | S1 QC, S2 further global views, S3 compartments |
| `F02` | pathways PAS moves in each sex: volcano rings, term tree | S4 PAS-responsive proteins, literature-informed pathways |
| `F03` | whether the sexes respond differently: female vs male effects, labels, signature carry-over, interaction rings | |
| `F04` | whether the signature tracks phenotype: PAS effect on phenotypes, proteome links, signature score, classification | |
| `F05` | the biology behind the sex difference: phenotype-link pathways, modules, OXPHOS, leading-edge proteins | |

Shared code: `R/figures.R` (theme, palettes, `contrast_bands()`, `pathway_ring()`,
`pathway_matrix()`, `save_figure()`) and `R/panels.R` (panel saving, `run_panels()`,
`write_panel_data()`). Colours: female PAS vermillion, male PAS blue, interaction purple, sex
contrasts grey; groups light (vehicle) and dark (PAS). Significant means p ≤ 0.05 or FDR ≤ 0.05;
proteins are also shown at FDR ≤ 0.10.
