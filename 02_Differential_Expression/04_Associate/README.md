# 04_Associate

Finds proteins whose abundance tracks each of six phenotypes across mice, after the four group
means are accounted for. Three models separate a shared relationship from a female-only one and
from a sex difference in slope.

- **Method:** `limpa::dpcDE(sample.weights = TRUE)` then `eBayes()`, phenotype scaled to SD 1. Models: `pooled slope` (one slope beside the four cell means, 20 mice), `females only` (10 mice, two group means plus slope), `slope difference, male minus female` (the direct test that the relationship differs by sex). BH within phenotype and model. Supplementary: `DGCA::ddcorAll()` Pearson r per group of each two-group comparison, Fisher z, permutation q (200 permutations, seed 20260910).
- **Inputs:** `01_Preprocess/02_Quantification/c_data/proteins.rds`, `02_Differential_Expression/02_Contrasts/c_data/contrasts.rds` (quantile-normalised full universe), `00_Input/phenotypes.xlsx`
- **Outputs:** `c_data/04_associate.xlsx` (sheets: summary, association, dgca, plot_index), `c_data/protein_association.rds` (read by `03_Pathway_Enrichment/04_Associate`), `b_reports/04_associate_<phenotype>.pdf` (page 1 volcanoes for the three models, page 2 DGCA r by group)
- **Result (1,842 tests, 92 expected at p ≤ 0.05):** CSA, pooled slope: 466 at p ≤ 0.05, 169 at FDR ≤ 0.05, 327 at FDR ≤ 0.10; females only: 84, none. Soleus octanoyl conductance, females only: 245, 1 at FDR ≤ 0.05, 4 at FDR ≤ 0.10; slope difference: 149, 2 at FDR ≤ 0.05. Plantaris pyruvate, pooled: 288, 4 at FDR ≤ 0.10. Capillary density, slope difference: 303, 1 at FDR ≤ 0.10. Satellite cells and soleus pyruvate: none at FDR ≤ 0.10 in any model.
- **Note:** Full universe only. The script stops unless CSA pooled gives 327 at FDR ≤ 0.10. Soleus and plantaris pyruvate conductance give fewer nominal hits than chance in several models (34 to 52).
