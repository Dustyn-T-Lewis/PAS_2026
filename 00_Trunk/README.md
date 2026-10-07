# 00 · Trunk

Everything the arcs read, computed once on all 20 mice. No step here answers an arc's question by itself; the arcs slice these results.

| Step | Asks | Method |
|---|---|---|
| [`00_Input`](00_Input/README.md) | what was delivered | data only |
| [`01_Filtering`](01_Filtering/README.md) | which precursors to keep | limpa `readDIANN`, q ≤ 0.01, contaminant and shared-peptide filters |
| [`02_Quantification`](02_Quantification/README.md) | protein abundance per mouse | limpa `dpcQuant`, MitoCarta, quantile normalisation per universe |
| [`03_Design`](03_Design/README.md) | which comparisons | four cell means, five contrasts |
| [`04_Gene_Sets`](04_Gene_Sets/README.md) | which protein groups to test | GO:BP, Reactome, MitoCarta, compartments, CORUM, STRING; membership table |
| [`05_Scores`](05_Scores/README.md) | set activity per mouse | singscore |
| [`06_Modules`](06_Modules/README.md) | co-abundance modules | WGCNA after removing sex |
| [`07_Fit`](07_Fit/README.md) | contrasts, phenotype effects and associations | limpa `dpcDE`, limma |
| [`08_Enrichment`](08_Enrichment/README.md) | pathway shifts per contrast and per phenotype | fgsea (enrichVolcano), redundancy flags, cross-sex signatures, STRING |
| [`09_Classification`](09_Classification/README.md) | does one score separate PAS from vehicle | pROC AUC, exact Wilcoxon, 252-relabelling null |
| [`10_Global`](10_Global/README.md) | what shapes the whole proteome | PCA, PERMANOVA |

Steps run in number order (`run_all.R`). Each holds `a_script/` (one script), `b_reports/` (one PDF, where there is something to draw), `c_data/` (one workbook plus the `.rds` later steps read) and a README naming what it reads and writes.
