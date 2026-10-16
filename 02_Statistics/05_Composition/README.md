# 05_Composition

What went into each pellet, scored per mouse and tested for PAS within each sex.

| | |
|---|---|
| Reads | `01_Processing/02_Quantification/c_data/proteins.rds`, `00_Setup/Input/report.parquet` (haemoglobin), `00_Setup/Input/metadata.csv` |
| Writes | `c_data/05_composition.xlsx`, `c_data/composition.rds` (read by 07_Sensitivity and `04_Figures/S1_Processing`), `b_reports/05_composition.pdf` |
| Run | `Rscript 02_Statistics/05_Composition/a_script/05_composition.R` |

- **Why:** a mitochondria-enriched pellet also carries sarcoplasmic reticulum, myofibril, nerve, blood and vessel material, and a PAS effect on that mix could pass for a change in the tissue.
- **Scores:** for each marker group, the mean over its proteins detected in at least 15 of 20 mice of each protein's log2 abundance, centred and scaled across mice. Haemoglobin, removed from the proteome as a contaminant, is scored as its share of precursor signal in the raw DIA-NN report (q ≤ 0.01). The pellet MitoCarta share comes from 02_Quantification.
- **Test:** Welch t, PAS against vehicle within each sex, with the difference and its 95% CI.
- **Result:** in females the sarcoplasmic reticulum score rises with PAS (+0.66 SD, 95% CI 0.34 to 0.98, p 0.0015); no other female score moves. In males the nerve (p 0.008), red-cell (p 0.02), vessel (p 0.0003) and myofibril (p 0.02) scores fall with PAS. Haemoglobin and the MitoCarta share do not change in either sex.
- **Reading:** part of the female SR and myofibril signal and much of the male nerve, blood and vessel signal can come from what was pelleted or dissected. Protein and pathway results in those classes should be read with this beside them; 07_Sensitivity refits the contrasts with the pellet share as a covariate.
