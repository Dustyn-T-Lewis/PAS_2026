# F02_Proteins

Figure 2: within each sex, which proteins are the strongest candidates for follow-up, how they group by known interactions, and how every mouse looks for the strongest of them.

![Figure 2](b_reports/F02_composite.png)

| | |
|---|---|
| Reads | `02_Statistics/06_Candidates/c_data/candidates.rds`, `02_Statistics/01_Fit/c_data/fit.rds`, `01_Processing/02_Quantification/c_data/proteins.rds` |
| Writes | `b_reports/panels/<panel>.pdf/.png`, `b_reports/F02_composite.pdf/.png` (180 × 240 mm), `c_data/F02.xlsx` |
| Run | `Rscript 04_Figures/F02_Proteins/a_script/F02_composite.R` (each panel script also runs alone) |

The composite checks every panel's protein count against `06_candidates.xlsx`. Candidates come from six routes (statistics, every mouse, both sexes, phenotype, module hub, pathway leading edge); 02_Statistics/06_Candidates describes each. They are discovery-stage candidates, not validated biomarkers.

## Panels

| Panel | Shows | Source |
|---|---|---|
| a | Female proteins passing two or more routes (58) as a STRING v12 network (combined score ≥ 700). Shaded hulls: STRING communities, named by each community's smallest-p over-represented Hallmark, Reactome or GO:BP term. Node size: routes passed; colour: log2 FC. Proteins without a STRING partner are listed under the panel. | `candidates.rds` |
| b | Male proteins passing two or more routes (47), drawn the same way. Several male communities are nerve (Mbp, Nefl, Nefm), red-cell (Spta1, Slc4a1, Ank1, Epb41) and extracellular matrix proteins, classes that fall in the male composition scores (S1). | `candidates.rds` |
| c | Every mouse for the 14 female proteins passing three or more routes: log2 abundance minus the own-sex vehicle mean, all four groups, bar at the group mean; open circles where the protein was not detected. | `fit.rds`, `proteins.rds` |
| d | The same for the 11 male proteins passing three or more routes. | `fit.rds`, `proteins.rds` |

The moderated F-test heatmap and its pattern enrichment, which this figure replaced, are in `S6_Protein_patterns`.

## Workbook

`c_data/F02.xlsx` holds one sheet per panel, plus `read_me`, `input_manifest` and `package_versions`.
