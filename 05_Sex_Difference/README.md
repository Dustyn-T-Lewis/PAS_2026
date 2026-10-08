# 05_Sex_Difference

Does the PAS response differ between the sexes? The test is the interaction contrast `treatment_by_sex`, (PAS_F − VEH_F) − (PAS_M − VEH_M); positive values mean PAS raised a protein or set more in females than in males. A feature that changes in one sex and not the other is not, on that basis, a sex difference; only this contrast tests it. Male–female differences at baseline are not a question here.

| Step | Asks | Method |
|---|---|---|
| [`01_Proteins`](01_Proteins/README.md) | which proteins respond differently | limpa `dpcDE` interaction contrast; p, Π and FDR tiers; STRING |
| [`02_Pathways`](02_Pathways/README.md) | which pathways respond differently | fgsea on GO:BP, Reactome and MitoCarta, beside the singscore contrast |

The phenotype interaction (capillary density p 0.014; the other five p ≥ 0.21) is in `02_Statistics/01_Fit` (`phenotype_interaction`) and in `06_Figures/F01_Overview` panel e. There is no classification here: the interaction is a difference of differences between groups, not a label any mouse carries.
