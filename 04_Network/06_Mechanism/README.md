# 06_Mechanism

Asks whether the proteins PAS moved are linked to each other in STRING, and which regulators the
INDRA literature graph says would move them, using each service's own documented test.

- **Method:** STRINGdb v12.0 (mouse, score ≥ 400) `get_ppi_enrichment` against the 1,816 STRING-mapped detected proteins, BH across lists. INDRA CoGEx REST API: `discrete_analysis` on the female FDR ≤ 0.10 list (upstream, kinases, TFs, downstream; BH, min evidence 2) and `signed_analysis` on the female and male Π ≤ 0.05 up and down lists (BH within contrast here). Mouse symbols mapped to human with `babelgene`; background 1,758 human orthologs.
- **Inputs:** `02_Differential_Expression/02_Contrasts/c_data/contrasts.rds`
- **Outputs:** `c_data/06_mechanism.xlsx` (sheets: string, indra_summary, indra_discrete, indra_signed, coverage, settings, plot_index), `c_data/mechanism_cache.rds`, `b_reports/06_mechanism.pdf`
- **Result:** the 8 female FDR proteins share no STRING edge (0 observed, 0 expected). Every Π list is more connected than chance at FDR ≤ 0.05: female down 3.3-fold (60 proteins), female up 2.4-fold (49), male down 5.0-fold (49), male up 2.3-fold (29, FDR 0.017), interaction down 1.5-fold (83), interaction up 2.4-fold (80). No INDRA regulator reaches q ≤ 0.05 in any analysis (smallest q 0.99).
- **Note:** connectivity shows the lists are coherent, not a mechanism. The interaction is not given to `signed_analysis` because its sign is a sex difference, not a direction in tissue. Service answers were fetched 2026-10-04 and are cached; the cache is refreshed when the lists change.
