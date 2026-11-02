# Shared_functions

Code sourced by more than one script, so a contrast, a colour or a workbook layout is defined once.

| File | Holds | Sourced by |
|---|---|---|
| `helpers.R` | seed, α, the q-value cutoff, contrasts, the three comparisons, phenotypes, test window; input and output paths (`stage_paths()`); the overlap rule for shown pathway sets (`representative_sets()`); the cached dpcQuant fit (`dpcquant_cached()`); `normalize()`; `run_depth_log2()`; STRING (`string_db()`); Π and tier counts; `one_per_gene()`; `label_groups()` for naming clusters and modules; the heatmap values (`centre_on_vehicle()`, `heatmap_values()`) and the co-expression tree 02_Enrichment cuts (`coexpression_tree()`); `preview_pdf()` and `gsea_table_page()` for each step's preview; `write_workbook()` (opens on `read_me`, ends with `package_versions`) | every script |
| `figures.R` | HRvLR's figure style: group colours and a three-tone contrast palette from the design panel, green and magenta for up and down, `theme_figure()`, `contrast_strips()` (ggh4x boxes), `save_figure()` for panels and composites, into `b_reports/` or its `supp/` (`figure_dir()`); `run_panels()` and `write_figure_workbook()` (panel sheets only) for the composites; `format_p()`, `protein_label()` and `label_with_database()` for labels; the plots more than one panel script draws: `group_point_scales()`, `pca_panel()`, `pathway_grid()` with `list_once()` and `count_bars()`, `pi_heatmap_panel()` with `protein_heatmap()` and `ora_bars()`, `score_cells()`, `candidate_network()`, `candidate_dots()` | `helpers.R` |
| `comparisons.R` | where a comparison writes (`comparison_paths()`) and which mice it looks at (`comparison_mice()`) | the four steps in `03_Comparisons/code/` |
