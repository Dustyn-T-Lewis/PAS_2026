# 00_Gene_Sets

Builds every gene-set collection the later steps test, and names the ROS and literature-informed
sets that are reported whatever they show. Tests nothing.

- **Method:** `enrichVolcano::load_gene_sets(c("Hallmark", "Reactome", "KEGG", "GO:BP"), species = "Mus musculus")`
  (msigdbr 26.1.0, org.Mm.eg.db 3.23.0; KEGG is KEGG MEDICUS). MitoCarta 3.0 pathways from sheet 4
  of the MitoCarta file. GO Slim terms from the `goslim_generic.obo` that enrichVolcano ships; each
  GO:BP term carries its GO Slim BP ancestors (`GO.db::GOBPANCESTOR`). Compartment sets hold genes
  annotated to a GO Slim cellular component or any descendant (`GOALL` in org.Mm.eg.db).
- **Inputs:** `00_Input/Mouse.MitoCarta3.0.xls`
- **Outputs:** `c_data/00_gene_sets.xlsx` (sheets: catalog, ros_sets, literature_sets, versions),
  `c_data/gene_sets.rds` (read by 01, 02, 04, 07, 08)
- **Result:** 7,402 GO:BP, 1,817 Reactome, 642 KEGG MEDICUS, 50 Hallmark, 149 MitoCarta and 25
  compartment sets before trimming to the data. 7 ROS sets named in advance; 12 more
  literature-informed sets (energy status, translation, oxidative fibre shift, atrophy defence,
  pipecolate handling), all present.
- **Note:** The 12 literature-informed sets were named on 2026-10-05, after the contrasts had been
  seen (Burke et al. 2026; Wang et al. 2019). They are literature-informed, not a priori.
