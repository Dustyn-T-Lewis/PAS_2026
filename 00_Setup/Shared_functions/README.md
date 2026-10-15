# Shared_functions

Code sourced by more than one script, so a contrast, a colour or a workbook layout is defined once.

| File | Holds | Sourced by |
|---|---|---|
| `helpers.R` | seed, α, contrasts, the three comparisons, phenotypes, test window; input and output paths; STRING; Π and tier counts; `write_workbook()` (opens on `read_me`, ends with `input_manifest` and `package_versions`) | every script |
| `figures.R` | HRvLR's figure style: group colours and a three-tone contrast palette from the design panel, green and magenta for up and down, `theme_figure()`, `contrast_strips()` (ggh4x boxes), `save_figure()`; `preview_pdf()` for each step's preview; the phenotype-effect plot | `helpers.R` |
| `comparisons.R` | where a comparison writes, which mice it looks at, a protein's tier | the four steps in `03_Comparisons/code/` |
