# 04_Gene_Sets

Builds every protein grouping the later steps test, before any result is read, and a table of every detected protein with the groups it belongs to.

| | |
|---|---|
| Reads | `00_Setup/Input/Mouse.MitoCarta3.0.xls`, `00_Setup/Input/corum_complete.txt`, `00_Setup/Input/STRING_v12/`, `01_Processing/02_Quantification/c_data/proteins.rds` |
| Writes | `c_data/04_gene_sets.xlsx`, `c_data/gene_sets.rds` (read by 05_Scores, 06_Modules, 02_Statistics/02_Enrichment and `03_Comparisons` step 02). No PDF. |
| Run | `Rscript 01_Processing/04_Gene_Sets/a_script/04_gene_sets.R` |

- **Tested pathways:** GO:BP and Reactome (disease and infection branch excluded) through `enrichVolcano::load_gene_sets()` in the full proteome; MitoCarta 3.0 pathways in the mitochondrial universe. A set is tested when 15 to 500 of its genes are detected (MitoCarta 10 to 500), counted on the symbols enrichVolcano maps, so fgsea (02_Statistics/02_Enrichment) and singscore (05_Scores) use the same list: 727 GO:BP, 119 Reactome, 43 MitoCarta.
- **Other groupings:** GO Slim cellular compartments (genes annotated to a slim term or any descendant); CORUM 5.3 (2026-04-14), all mammalian complexes mapped to mouse with babelgene, kept at ≥ 5 detected subunits (93); STRING v12 (score ≥ 400) walktrap clusters of ≥ 5 proteins (13), named by their best Hallmark or Reactome term at FDR ≤ 0.05. Each GO:BP term carries its GO Slim BP families, for display only. Hallmark names clusters and modules and supplies four of the literature and ROS sets.
- **Workbook:** `membership` has one row per detected protein with the count and names of its tested groups per collection; `membership_long` has one row per protein and group.
- **Note:** The 12 literature-informed sets were named on 2026-10-05, after the first contrast results had been seen (Burke et al. 2026; Wang et al. 2019); 7 ROS sets were named in advance. All are reported in 02_Statistics/02_Enrichment whatever they show.
