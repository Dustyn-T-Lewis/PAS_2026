# 04_Gene_Sets

Builds every protein grouping the later steps test, from annotation alone, and a table of each detected protein's groups.

| | |
|---|---|
| Reads | `00_Setup/Input/Mouse.MitoCarta3.0.xls`, `00_Setup/Input/corum_complete.txt`, `00_Setup/Input/STRING_v12/`, `01_Processing/02_Quantification/c_data/proteins.rds` |
| Writes | `c_data/04_gene_sets.xlsx`, `c_data/gene_sets.rds`; no PDF |
| Run | `Rscript 01_Processing/04_Gene_Sets/a_script/04_gene_sets.R` |

## Method

- Full proteome: GO:BP, Reactome (disease and infection branch excluded) and Hallmark through `enrichVolcano::load_gene_sets()`; KEGG LEGACY through `msigdbr`, mapped to mouse orthologs. Mitochondrial universe: MitoCarta 3.0 pathways.
- A set is tested when 15 to 500 of its genes are detected (MitoCarta 10 to 500), counted on enrichVolcano's symbols so fgsea and singscore test the same sets.
- Other groupings: GO Slim cellular compartments, tested in 02_Enrichment when 15 genes up to half the ranked genes are detected (`side_list_size()`); CORUM 5.3 mammalian complexes mapped to mouse with babelgene, kept at ≥ 5 detected subunits; STRING v12 (score ≥ 400) walktrap clusters of ≥ 5 proteins, named by their best Hallmark or Reactome term.
- 12 literature-informed sets were named on 2026-10-05, after the first contrast results had been seen; 7 ROS sets were named in advance.

## Results

- Tested: 755 GO:BP, 139 Reactome, 34 Hallmark, 33 KEGG, 44 MitoCarta, 98 complexes, 7 STRING clusters.
- `membership` has one row per detected protein with a gene symbol (1,961 of 1,962) and its tested groups per collection.
