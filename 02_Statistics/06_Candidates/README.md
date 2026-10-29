# 06_Candidates

Candidate proteins for follow-up, chosen by six routes and scored by how many each passes. These are discovery-stage candidates for a western blot or targeted assay in new mice, not validated biomarkers.

| | |
|---|---|
| Reads | `02_Statistics/01_Fit/c_data/fit.rds`, `02_Statistics/02_Enrichment/c_data/enrichment.rds`, `01_Processing/06_Modules/c_data/modules.rds`, `01_Processing/02_Quantification/c_data/proteins.rds`, `01_Processing/03_Design/c_data/design.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, STRING v12 |
| Writes | `c_data/06_candidates.xlsx`, `c_data/candidates.rds` (read by `04_Figures/F02_Proteins` for F02 and S4), `b_reports/06_candidates.pdf` |
| Run | `Rscript 02_Statistics/06_Candidates/a_script/06_candidates.R` |

A protein is considered in a sex only if it is detected in at least 4 of 5 mice in each group of that sex.

| Route | Rule | What it answers |
|---|---|---|
| Statistics | p ≤ 0.01 within the sex | Is the change larger than noise? |
| Every mouse | AUC ≥ 0.92 or ≤ 0.08 on surrogate-adjusted values | Is it consistent mouse by mouse? |
| Both sexes | same direction, p ≤ 0.10 in each sex | Does it hold in independent animals? |
| Phenotype | within-treatment association p ≤ 0.01 with any trait | Is it tied to an outcome, not only to the treatment? |
| Hub | kME ≥ 0.8 in a module with within-sex PAS p ≤ 0.05 | Does it sit at the centre of a responsive module? |
| Pathway | leading edge of a kept FDR-level fgsea set, p ≤ 0.05 | Does it carry a pathway result? |

- **Grouping:** candidates passing two or more routes are linked by STRING v12 (combined score ≥ 700) and split into communities (igraph walktrap); each community with two or more members is named by its smallest-p over-represented Hallmark, Reactome or GO:BP term (fgsea::fora against the detected proteome).
- **Result:** females 255 proteins pass any route, 58 pass two or more, 14 pass three or more (Cul5, Fkbp1a and Atp2a1 pass four). Males 208, 47 and 11. No female module responds to PAS after the module rebuild, so the hub route adds no female candidate.
- **Caveats:** route thresholds are deliberately loose and the routes are not independent (a strong protein tends to pass statistics and every-mouse together). The phenotype route crosses tissues (respiration from soleus and plantaris). Classes flagged in 05_Composition (sarcoplasmic reticulum in females; nerve, red cell, vessel and myofibril in males) should be read with that step's results beside them.
