# Shared_functions

Code sourced by more than one script, so a contrast, a colour or a workbook layout is defined once.

| File | Holds | Sourced by |
|---|---|---|
| `helpers.R` | seed, α, contrasts, the three comparisons, phenotypes, test window; input and output paths; STRING; Π and tier counts; `write_workbook()` (opens on `read_me`, ends with `input_manifest` and `package_versions`) | every script |
| `figures.R` | colours, `theme_figure()`, `save_figure()` for figures; `preview_pdf()`, which writes each step's preview from package or base plots; the phenotype-effect plot | `helpers.R` |
| `comparisons.R` | where a comparison writes, which mice it looks at, a protein's tier | the four steps in `03_Comparisons/code/` |
