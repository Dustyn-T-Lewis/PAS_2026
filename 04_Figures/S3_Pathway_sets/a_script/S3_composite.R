# Supplementary Figure 3, side lists: the GO Slim compartments and the sets named in advance by
# hypothesis, per contrast. Each panel is its own script in panels/; this file runs them, checks
# them against 02_enrichment.xlsx, stitches them with one key and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("a_compartments", "b_literature")
panels <- run_panels("S3_Pathway_sets", PANELS)

workbook <- output_path("02_Statistics", "02_Enrichment", "02_enrichment.xlsx")
stopifnot(
  nrow(panels$a_compartments$data) == nrow(readxl::read_excel(workbook, "compartments")),
  nrow(panels$b_literature$data) == nrow(readxl::read_excel(workbook, "literature"))
)

figure <- wrap_plots(map(panels, "plot"), widths = c(1, 1)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
out <- stage_paths("04_Figures", "S3_Pathway_sets")
save_figure(figure, out$reports, "S3_composite", width_mm = 180, height_mm = 120)

write_figure_workbook(
  panels,
  c(
    "Panel a: GO Slim compartments per contrast: NES, BH FDR, size.",
    "Panel b: sets named in advance per contrast: hypothesis, expected direction, NES, FDR, size."
  ),
  file.path(out$data, "S3.xlsx")
)

sessionInfo()
