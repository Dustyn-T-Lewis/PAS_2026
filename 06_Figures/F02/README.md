# F02 · Pathways by sex

- **a, b** enrichVolcano rings, full proteome (Hallmark, Reactome without disease terms, KEGG
  reference, GO:BP without GO Slim) and mitochondrial proteome (MitoCarta), female above male; the
  12 strongest collapsed fgsea terms each, * FDR ≤ 0.05. Females: myogenesis, muscle contraction,
  myofibril assembly up; cell cycle down. Males: OXPHOS, EMT, actin and lipid terms down.
- **c** The 15 strongest female and male terms, clustered on shared genes (EnrichmentMap
  similarity, cut at 0.375), with NES and FDR in each sex and for the interaction.

Supplement:

- **S4** (a) the 318 proteins at p ≤ 0.05 in either sex or the interaction, z-scored within sex;
  (b) literature-informed pathways (Burke 2026, Wang 2019) and the ROS sets, BH across the list.
  Females: striated muscle contraction up as predicted (list FDR 0.008); energy and translation
  sets flat. Chosen after results, all shown.

Sources: `03_Pathway_Enrichment/08_Display`, `02_Differential_Expression/02_Contrasts`.
