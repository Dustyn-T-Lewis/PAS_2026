# F01 · Global proteome

- **a** PCA of the 20 mice, groups shaded. PERMANOVA: sex R² 0.38, p = 0.001; treatment R² 0.04,
  p = 0.27; interaction p = 0.26. MitoCarta share fitted as a vector (envfit R² 0.89, p = 0.001).
- **b** Sex differences at FDR ≤ 0.05: 258 down, 266 up in vehicle; 406 down, 429 up in PAS.
- **c** Proteins PAS moves at p ≤ 0.05 and Π ≤ 0.05 (side by side; the tiers overlap). FDR ≤ 0.05:
  4 in females, none in males or for the interaction.
- **d** fgsea terms at FDR ≤ 0.05 per contrast, up and down, by database.

Supplements:

- **S1** QC per mouse: data depth, contamination, normalisation, MitoCarta share.
- **S2** Variance by component, response magnitude (RRPP, null), effect sizes, mitochondrial MDS,
  female vs male effects (ρ −0.07 full, −0.35 mito), p-value histograms, baseline sex pathways.
- **S3** Compartment shifts by fry for each sex and the interaction; none at FDR ≤ 0.05.

Sources: `05_Summary/05_Global`, `05_Summary/06_Magnitude`, `05_Summary/01_Signal`,
`02_Differential_Expression/02_Contrasts`, `03_Pathway_Enrichment/02_Contrasts`, `07_Compartments`,
`08_Display`, `01_Preprocess`.
