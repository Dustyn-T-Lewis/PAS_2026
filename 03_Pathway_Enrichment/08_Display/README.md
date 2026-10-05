# 08_Display

Chooses once which terms the pathway figures draw, so every figure agrees. Nothing is retested and
no p-value changes, except the BH across the literature-informed list.

- **Method:** `display_pool()` in `R/helpers.R` keeps Hallmark, Reactome without disease and
  infection terms, KEGG MEDICUS reference pathways, GO:BP without the GO Slim terms, and all
  MitoCarta. `enrichVolcano::dedup_terms(method = "enrichmentmap", term_threshold = 1)` then
  collapses every term across databases (pooled under one name for the call).
  - Ring terms: 12 smallest-p collapsed terms per universe for female PAS, male PAS and
    treatment_by_sex, either direction.
  - Term tree: 15 smallest-p uncollapsed full-proteome terms per sex, average-linkage clustering on
    the EnrichmentMap coefficient (mean of Jaccard and overlap; Merico et al. 2010), cut at 0.375.
    Each cluster takes the name of its smallest-p member.
  - Baseline sex terms: 10 smallest-p collapsed full-proteome terms for sex in VEH and sex in PAS
    (male minus female), for Figure 1.
  - Literature-informed sets from `00_Gene_Sets`, BH across the list within each primary contrast
    (`list_fdr`).
- **Inputs:** `00_Gene_Sets/c_data/gene_sets.rds`, `02_Contrasts/c_data/enrichment.rds`
- **Outputs:** `c_data/08_display.xlsx` (sheets: ring_terms, tree_terms, tree_dots, cohesion,
  sex_dots, literature), `c_data/display.rds` (read by the figure scripts); no PDFs
- **Result:** Ring terms at FDR ≤ 0.05: full proteome 5, 5 and 7 of 12 (female, male,
  interaction); MitoCarta 0, 1 and 3. Tree: 29 terms in 15 clusters, 6 with more than one term.
  Baseline sex: 18 terms. Literature list at `list_fdr` ≤ 0.05: striated muscle contraction up in
  female PAS (NES 1.98) and the interaction (2.14), down in male PAS (-2.06); OXPHOS down in male
  PAS (-1.66), up in the interaction (1.51); MitoCarta mitochondrial central dogma up in male PAS
  (2.09), down in the interaction (-1.81).
- **Note:** The literature sheet carries 17 of the 19 named sets. Reactome ROS detoxification is too
  small to test; GOBP_CYTOPLASMIC_TRANSLATION was tested in 02 (32 proteins) but `display_pool()`
  removes it as a GO Slim term.
