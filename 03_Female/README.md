# 01 · Female

Does PAS change the proteome and the phenotypes of females? Every number comes from the `PAS_vs_VEH_female` contrast (PAS_F − VEH_F) of the shared fit and from the tests 02_Statistics runs on it; this comparison only slices and reports.

| Substage | Asks | Method |
|---|---|---|
| [`01_Proteins`](01_Proteins/README.md) | which proteins PAS changed | limpa `dpcDE` contrast; p, Π and FDR tiers; STRING connectivity of the Π lists |
| [`02_Pathways`](02_Pathways/README.md) | which pathways PAS shifted | fgsea on GO:BP, Reactome and MitoCarta, beside the singscore contrast; compartments; literature sets |
| [`03_Phenotypes`](03_Phenotypes/README.md) | what PAS did to each phenotype, and which features move with each phenotype | phenotype effects (% change); phenotype associations across the females' ten mice (limma), fgsea on association t |
| [`04_Classification`](04_Classification/README.md) | does one score separate PAS from vehicle | pROC AUC, exact Wilcoxon, 252-relabelling null |

Each substage holds `a_script/` (two lines calling `R/comparisons.R`), `b_reports/` (one PDF: overview pages, then a grid of every feature at p ≤ 0.05), `c_data/` (one workbook and one `.rds`) and a README. To read a result by hand, open the workbook's `read_me` sheet first; every row carries its feature, statistic, p, FDR and r.
