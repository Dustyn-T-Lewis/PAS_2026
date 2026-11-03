# 06_Candidates

Picks candidate proteins for follow-up by six routes and scores each by how many it passes. These are discovery-stage candidates for a western blot or targeted assay in new mice, not validated biomarkers.

| | |
|---|---|
| Reads | `02_Statistics/01_Fit/c_data/fit.rds`, `02_Statistics/02_Enrichment/c_data/enrichment.rds`, `01_Processing/02_Quantification/c_data/proteins.rds`, `01_Processing/03_Design/c_data/design.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `01_Processing/06_Modules/c_data/modules.rds`, `00_Setup/Input/STRING_v12/` |
| Writes | `c_data/06_candidates.xlsx`, `c_data/candidates.rds`, `b_reports/06_candidates.pdf` |
| Run | `Rscript 02_Statistics/06_Candidates/a_script/06_candidates.R` |

## Method

A protein is considered in a sex only if detected in at least 4 of 5 mice in each group of that sex.

| Route | Rule |
|---|---|
| Statistics | p ≤ 0.01 within the sex |
| Every mouse | AUC ≥ 0.92 or ≤ 0.08 on surrogate-adjusted values |
| Both sexes | same direction, p ≤ 0.10 in each sex |
| Phenotype | within-treatment association p ≤ 0.01 with any phenotype |
| Hub | kME ≥ 0.8 in a module with within-sex PAS p ≤ 0.05 |
| Pathway | in the leading edge of a kept FDR-level fgsea set, with p ≤ 0.05 |

- Candidates passing two or more routes are linked by STRING v12 (combined score ≥ 700) and split into walktrap communities, each named by its smallest-p over-represented Hallmark, Reactome or GO:BP term.
- The thresholds are loose and the routes overlap: a strong protein tends to pass statistics and every mouse together.

## Results

- Females: 255 proteins pass any route, 58 pass two or more, 14 pass three or more. Cul5, Fkbp1a and Atp2a1 pass four.
- Males: 208, 47 and 11.
- No female module responds to PAS, so the hub route adds no female candidate.
- Read sarcoplasmic reticulum candidates in females, and nerve, red-cell, vessel and myofibril candidates in males, beside 05_Composition.
