# 02_Statistics

Every model and test, run once on all 20 mice after `01_Processing`. The comparisons and figures read these results and refit nothing. Steps run in number order.

![The seven statistics steps, the five contrasts and who reads the results](schematic.png)

| Step | Holds |
|---|---|
| [`01_Fit`](01_Fit/README.md) | protein, score and eigengene contrasts; phenotype effects and associations |
| [`02_Enrichment`](02_Enrichment/README.md) | fgsea per contrast and per phenotype, with camera, ORA, cross-sex signatures and STRING connectivity |
| [`03_Classification`](03_Classification/README.md) | how well each single score separates PAS from vehicle |
| [`04_PCA`](04_PCA/README.md) | PCA and PERMANOVA of the proteome and the phenotypes |
| [`05_Composition`](05_Composition/README.md) | pellet composition scores and their PAS tests |
| [`06_Candidates`](06_Candidates/README.md) | follow-up candidate proteins by six routes |
| [`07_Sensitivity`](07_Sensitivity/README.md) | protein results refitted under other choices |
