# Shared_functions

R code sourced by more than one script, so a contrast, a colour or a workbook layout is defined once.

| File | Holds | Sourced by |
|---|---|---|
| `helpers.R` | constants (seed, α, contrasts, comparisons, phenotypes, set-size window); paths (`stage_paths()`); the cached dpcQuant fit (`dpcquant_cached()`); `normalize()`; `run_depth_log2()`; `string_db()`; Π and tier counts; the overlap rule for shown sets (`representative_sets()`); heatmap values and the co-expression tree (`coexpression_tree()`); preview PDFs; `write_workbook()` | every script |
| `figures.R` | figure theme, group and direction colours, `save_figure()`, `write_figure_workbook()`, label helpers, and the plots several panels share (`pca_panel()`, `pathway_grid()`, `protein_heatmap()`, `candidate_network()`, `candidate_dots()` and others) | `helpers.R` |
| `comparisons.R` | where a comparison writes (`comparison_paths()`) and which mice it uses (`comparison_mice()`) | the four scripts in `03_Comparisons/code/` |
