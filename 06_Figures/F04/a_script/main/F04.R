#!/usr/bin/env Rscript
# Figure 4, whether the PAS signature tracks muscle phenotype: (a) PAS effect on each phenotype,
# (b) how much of the proteome tracks each phenotype, (c) the female signature score against the
# phenotypes PAS changed, (d) classification of PAS against vehicle. Runs the panel scripts in
# panels/ and assembles them.

source(here::here("R", "panels.R"))

panels <- run_panels("F04", "main", c(
  "A_phenotype_effects", "B_proteome_links", "C_signature_score", "D_classification"
))
plots <- map(panels, "plot")
top_row <- plots$A_phenotype_effects + plots$B_proteome_links + plot_layout(widths = c(1, 1.15))
bottom_row <- plots$C_signature_score + plots$D_classification + plot_layout(widths = c(1.35, 1))
figure <- top_row / bottom_row + plot_annotation(tag_levels = "a")
save_figure(figure, figure_dirs("F04", "main")$reports, "F04",
  height_mm = 165, width_mm = LANDSCAPE_MM
)
write_panel_data(panels, "F04", "F04")

sessionInfo()
