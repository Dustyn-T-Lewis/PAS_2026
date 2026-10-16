# Figure 4, phenotype-association pathways: within each sex, which pathway sets track each
# phenotype across the mice of that sex, pooled over treatment and within treatment. Each panel is
# its own script in panels/; this file runs them, stitches them with one shared key and writes one
# workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("a_female", "b_male")
panels <- run_panels("F04_Phenotype_pathways", PANELS)

figure <- wrap_plots(
  panels$a_female$plot, panels$b_male$plot + guides(colour = "none"),
  ncol = 2
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
out <- stage_paths("04_Figures", "F04_Phenotype_pathways")
save_figure(figure, out$reports, "F04_composite", width_mm = 180, height_mm = 150)

write_figure_workbook(
  panels,
  c(
    "Panel a: each female set shown, per phenotype and model: NES, FDR, collection, trait.",
    "Panel b: each male set shown, per phenotype and model: NES, FDR, collection, trait."
  ),
  file.path(out$data, "F04.xlsx")
)
