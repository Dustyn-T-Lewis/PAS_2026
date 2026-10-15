# 02_Statistics

Every model and test, run once on all 20 mice. The three comparisons read these results; nothing here is computed again downstream.

![The four statistics steps, the five contrasts and who reads the results](schematic.png)

| Step | Asks | Method |
|---|---|---|
| [`01_Fit`](01_Fit/README.md) | contrasts, phenotype effects and phenotype associations | limpa `dpcDE`, limma |
| [`02_Enrichment`](02_Enrichment/README.md) | which pathways shift, per contrast and per phenotype | fgsea (enrichVolcano), redundancy flags, cross-sex signatures, STRING |
| [`03_Classification`](03_Classification/README.md) | does one score separate PAS from vehicle | pROC AUC, exact Wilcoxon, 252-relabelling null |
| [`04_PCA`](04_PCA/README.md) | what shapes the whole proteome and the phenotypes | PCA, PERMANOVA |

Steps run in number order after `01_Processing`. Each holds `a_script/`, `b_reports/`, `c_data/` and a README naming what it reads and writes.
