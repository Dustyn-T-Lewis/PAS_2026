# 03_Biclusters

Looks for groups of proteins that behave alike in a subset of mice, without using group labels.
Exploratory.

- **Method:** `biclust::biclust(method = BCPlaid())` on the protein-standardised normalised full matrix, 100 seeds; biclusters whose protein sets overlap at Jaccard ≥ 0.5 are merged (complete linkage); core proteins and mice are those found in at least half of a merged bicluster's runs.
- **Inputs:** `02_Differential_Expression/02_Contrasts/c_data/contrasts.rds`
- **Outputs:** `c_data/03_biclusters.xlsx` (sheets: biclusters, plot_index), `b_reports/03_biclusters.pdf`
- **Result:** 28 merged biclusters. The stable ones split by sex, not treatment: the most recurrent (78% of seeds, 161 proteins) holds six males (3 PAS, 3 VEH); the next (69%, 86 proteins) eight females including all five PAS females; the third (48%, 121 proteins) again eight females (5 PAS, 3 VEH).
- **Note:** recurrence across seeds measures stability of the algorithm, not a significance test.
