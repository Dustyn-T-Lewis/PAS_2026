# Shared_functions

Code sourced by more than one script, so a contrast, a colour or a workbook layout is defined once.

| File | Holds | Sourced by |
|---|---|---|
| `helpers.R` | seed, α, contrasts, the three comparisons, phenotypes, test window; input and output paths; STRING; Π and tier counts; `write_workbook()` (opens on `read_me`, ends with `input_manifest` and `package_versions`) | every script |
| `figures.R` | colours, `theme_figure()`, `save_figure()`, `write_report()`, heatmap and phenotype-effect plots | `helpers.R` |
| `comparisons.R` | where a comparison writes, which mice it looks at, its colours and the report grids | the four steps in `03_Comparisons/code/` |
