# 06_Magnitude

Asks how far PAS moves each mouse's whole proteome, by sex, against the spread among vehicle mice.

- **Method:** each mouse's Euclidean distance from its own sex's vehicle centroid on the normalised full universe; a vehicle mouse is measured against the other four vehicle mice of its sex. `RRPP::lm.rrpp(distance ~ sex * treatment)`, 999 permutations.
- **Inputs:** `02_Differential_Expression/02_Contrasts/c_data/contrasts.rds`
- **Outputs:** `c_data/06_magnitude.xlsx` (sheets: magnitude, tests, medians), `c_data/magnitude.rds` (read by `06_Figures/F01`)
- **Result:** PAS mice sit no further out than vehicle mice (treatment R² 0.08, p = 0.24) and the shift does not differ by sex (p = 0.68); sex p = 0.51. Median distances: VEH F 33.6, PAS F 30.1, VEH M 30.6, PAS M 27.8.
- **Note:** measuring vehicle mice against four rather than five others inflates their expected squared distance slightly (×1.25 against ×1.2 for PAS mice), which works against finding a PAS shift. Noted, not corrected.
