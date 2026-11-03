# 05_Composition

Scores what went into each pellet besides mitochondria and tests whether PAS changed it within each sex. A PAS effect on the pellet mix could pass for a change in the tissue.

| | |
|---|---|
| Reads | `01_Processing/02_Quantification/c_data/proteins.rds`, `00_Setup/Input/report.parquet` (haemoglobin) |
| Writes | `c_data/05_composition.xlsx`, `c_data/composition.rds`, `b_reports/05_composition.pdf` |
| Run | `Rscript 02_Statistics/05_Composition/a_script/05_composition.R` |

## Method

- Marker scores (sarcoplasmic reticulum, endoplasmic reticulum, myofibril, nerve, red cell, vessel): mean of z-scored log2 abundance over each group's markers detected in at least 15 of 20 mice.
- Haemoglobin, removed from the proteome as a contaminant, is scored as its share of precursor signal in the report. The MitoCarta share comes from 02_Quantification.
- Welch t, PAS against vehicle within each sex, with the difference and its 95% CI.

## Results

- Females: sarcoplasmic reticulum rises with PAS (+0.66 SD, 95% CI 0.34 to 0.98, p 0.0015); no other score moves.
- Males: vessel (p 0.0003), nerve (p 0.008), red cell (p 0.02) and myofibril (p 0.02) fall with PAS.
- Haemoglobin and the MitoCarta share do not change in either sex.
- Protein and pathway results in these classes may reflect what was pelleted or dissected. 07_Sensitivity refits with the pellet share as a covariate.
