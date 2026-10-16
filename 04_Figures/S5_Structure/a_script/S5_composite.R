# Supplementary Figure 5, structure: whether more features separate PAS from vehicle perfectly
# than chance allows, and what the main proteome components track. Each panel is its own script in
# panels/; this file runs them, stitches them and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("a_perfect_separators", "b_component_r2")
panels <- run_panels("S5_Structure", PANELS)

figure <- wrap_plots(map(panels, "plot"), widths = c(1, 1))
out <- stage_paths("04_Figures", "S5_Structure")
save_figure(figure, out$reports, "S5_composite", width_mm = 180, height_mm = 60)

write_figure_workbook(
  panels,
  c(
    "Panel a: perfect separators per sex and level: observed, null median and 95th centile, p.",
    "Panel b: R² of sex, treatment, interaction, MitoCarta share and run depth on PC1 to PC4."
  ),
  file.path(out$data, "S5.xlsx")
)
