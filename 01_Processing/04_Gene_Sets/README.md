# 04_Gene_Sets

Builds every protein grouping the later steps test, from annotation alone without reading the samples, and a table of every detected protein with the groups it belongs to.

| | |
|---|---|
| Reads | `00_Setup/Input/Mouse.MitoCarta3.0.xls`, `00_Setup/Input/corum_complete.txt`, `00_Setup/Input/STRING_v12/`, `01_Processing/02_Quantification/c_data/proteins.rds` |
| Writes | `c_data/04_gene_sets.xlsx`, `c_data/gene_sets.rds` (read by 05_Scores, 06_Modules, 02_Statistics/02_Enrichment, 02_Statistics/06_Candidates and `04_Figures/F05_Complexes`). No PDF. |
| Run | `Rscript 01_Processing/04_Gene_Sets/a_script/04_gene_sets.R` |

- **Tested pathways:** GO:BP, Reactome (disease and infection branch excluded) and Hallmark through `enrichVolcano::load_gene_sets()`, and KEGG LEGACY through `msigdbr` (human MSigDB only, mapped to mouse orthologs), in the full proteome; MitoCarta 3.0 pathways in the mitochondrial universe. A set is tested when 15 to 500 of its genes are detected (MitoCarta 10 to 500). Genes are counted on the symbols enrichVolcano maps, so fgsea (02_Statistics/02_Enrichment) and singscore (05_Scores) test the same list: 755 GO:BP, 139 Reactome, 34 Hallmark, 33 KEGG, 44 MitoCarta.
- **Other groupings:** GO Slim cellular compartments (genes annotated to a slim term or any descendant); CORUM 5.3 (2026-04-14), all mammalian complexes mapped to mouse with babelgene, kept at ≥ 5 detected subunits (98); STRING v12 (score ≥ 400) walktrap clusters of ≥ 5 proteins (7), named by their best Hallmark or Reactome term at FDR ≤ 0.05. Each GO:BP term carries its GO Slim BP families, for display only. Hallmark names clusters and modules and supplies four of the literature and ROS sets.
- **Basis:** the 15 to 500 window is the GSEA default (Subramanian et al. 2005) and enrichVolcano's; STRING clusters come from `STRINGdb::get_clusters()`; GO Slim generic is the GO Consortium's subset.
- **Workbook:** `membership` has one row per detected protein with a gene symbol (1,961 of 1,962), giving the count and names of its tested groups per collection; `membership_long` has one row per protein and group. In `catalog`, `tested` also marks the literature and ROS sets, so Reactome shows 140: the ROS set detoxification of reactive oxygen species has 13 detected genes, below the window, and is tested only in that list.
- **Note:** The 12 literature-informed sets were named on 2026-10-05, after the first contrast results had been seen (Burke et al. 2026; Wang et al. 2019); the 7 ROS sets were named in advance. All are reported in 02_Statistics/02_Enrichment whatever they show.
