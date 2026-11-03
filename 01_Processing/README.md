# 01_Processing

From the DIA-NN report to protein values per mouse, plus every protein grouping later steps test. Nothing here compares groups. Steps run in number order.

![The six processing steps, their methods and the files they write](schematic.png)

| Step | Holds |
|---|---|
| [`01_Filtering`](01_Filtering/README.md) | the filtered precursor matrix |
| [`02_Quantification`](02_Quantification/README.md) | protein values and standard errors per mouse, MitoCarta annotation, normalisation |
| [`03_Design`](03_Design/README.md) | the design matrix, surrogate variable and five contrasts |
| [`04_Gene_Sets`](04_Gene_Sets/README.md) | pathways, compartments, complexes, STRING clusters and the protein membership table |
| [`05_Scores`](05_Scores/README.md) | one singscore per mouse per set |
| [`06_Modules`](06_Modules/README.md) | WGCNA co-abundance modules and their eigengenes |
