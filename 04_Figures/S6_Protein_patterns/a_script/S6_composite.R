# Supplementary Figure 6, protein patterns: the proteins PAS changes in either sex (moderated F),
# by response pattern, and what each pattern holds. Each panel is its own script in panels/; this
# file runs them, stitches them and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("a_ftest_heatmap", "b_pattern_ora")
panels <- run_panels("S6_Protein_patterns", PANELS)

figure <- wrap_plots(map(panels, "plot"), widths = c(116, 25))
out <- stage_paths("04_Figures", "S6_Protein_patterns")
save_figure(figure, out$reports, "S6_composite", width_mm = 180, height_mm = 145)

write_figure_workbook(
  panels,
  c(
    paste(
      "Panel a: each protein at moderated F p ≤ 0.05 (full proteome), its pattern, F, p, FDR,",
      "both PAS log2 FCs, and its pooled r with each trait per sex."
    ),
    "Panel b: the three smallest-p ORA sets per pattern: overlap, size, p, BH FDR."
  ),
  file.path(out$data, "S6.xlsx")
)
