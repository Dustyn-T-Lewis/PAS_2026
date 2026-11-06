# 01_Fit

Fits every model the comparisons read, once on all 20 mice: protein, score and eigengene contrasts, phenotype effects and phenotype associations.

| | |
|---|---|
| Reads | `01_Processing/02_Quantification/c_data/proteins.rds`, `01_Processing/03_Design/c_data/design.rds`, `01_Processing/04_Gene_Sets/c_data/gene_sets.rds`, `01_Processing/05_Scores/c_data/scores.rds`, `01_Processing/06_Modules/c_data/modules.rds`, `00_Setup/Input/phenotypes.xlsx` |
| Writes | `c_data/01_fit.xlsx`, `c_data/fit.rds`, `b_reports/01_fit.pdf` |
| Run | `Rscript 02_Statistics/01_Fit/a_script/01_fit.R` |

## Method

- Proteins: `limpa::dpcDE(sample.weights = TRUE)` per universe on the design. Each contrast is refitted as the first coefficient of an equivalent design (`limma::contrastAsCoef()`) and moderated with `eBayes()`, because `contrasts.fit()` is not exact under limpa's per-protein weights. BH within universe and contrast; Π, TREAT (fold change > 1.2) and ashr lfsr beside each row.
- Every protein row carries `detected_VEH_F`, `detected_PAS_F`, `detected_VEH_M` and `detected_PAS_M`: the mice in each group with an observed precursor. The other values are limpa's inference, so a protein never detected in a group still has a fold change there.
- Set scores: limma with `eBayes()` on the same design. Module eigengenes: ordinary least-squares t, since eight features cannot inform an eBayes prior.
- Shown score sets: within each contrast and universe, sets at p ≤ 0.05 are grouped by `representative_sets()` (Jaccard tree cut at 0.70) and each group is shown by one set (`overlap_group`, `shown`).
- Phenotype effects: % change of PAS against the same sex's vehicle, Welch t on log values. Interaction: `lm(log(value) ~ treatment * sex)`.
- Phenotype associations: each feature regressed on one z-scored phenotype within one sex, with the surrogate variable fitted within sex. `pooled` uses all ten mice of the sex; `treatment` adds group means, so it asks for co-movement among mice given the same treatment. `pooled_mito` adds the pellet MitoCarta share as a check.
- Every limma table carries `r_equiv = t / sqrt(t^2 + df)`, so contrasts and associations share one colour scale in the figures.

## Results

- Full proteome (1,962 tests), at p ≤ 0.05 / Π ≤ 0.05 / FDR ≤ 0.05: female 177 / 107 / 3 (Cops3, Cops8, Dusp3, all down); male 193 / 106 / 0; interaction 162 / 142 / 0.
- Of the Π ≤ 0.05 proteins, 1 of 107 female, 8 of 106 male and 40 of 142 interaction proteins were never detected in at least one group of their contrast (33 of the 40 in PAS males), so the value there is limpa's estimate.
- Phenotype effects in females: capillary density +138% (p 0.005), CSA +22% (p 0.002), satellite cells +47% (p 0.009). None in males. Interaction only for capillary density (p 0.014).
- Shown score sets at p ≤ 0.05: 42 of 98 in females, 26 of 55 in males, 41 of 70 for the interaction.
