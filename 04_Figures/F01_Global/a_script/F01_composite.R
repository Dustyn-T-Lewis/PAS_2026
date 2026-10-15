# Figure 1, the global proteome: the design and the phenotypes on top, the proteome below. Each
# panel is its own script in panels/, run here in its own environment; this file stitches the
# returned plots, checks the plotted counts against their source tables, and writes one workbook
# with a sheet per panel.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c(
  "a_design", "b_phenotype_effects", "c_phenotype_pca", "d_proteome_pca", "e_protein_bars",
  "f_signature_barcode"
)
panels <- run_panels("F01_Global", PANELS)

fit_workbook <- output_path("02_Statistics", "01_Fit", "01_fit.xlsx")
fit_summary <- readxl::read_excel(fit_workbook, "contrast_summary") |>
  filter(level == "protein", contrast %in% COMPARISON_CONTRASTS)
plotted <- panels$e_protein_bars$data |>
  summarise(pi = sum(pi), fdr = sum(fdr), .by = c(universe, contrast)) |>
  inner_join(fit_summary, by = c("universe", "contrast"))
stopifnot(
  nrow(plotted) == 6, plotted$pi == plotted$pi_le_0.05, plotted$fdr == plotted$fdr_le_0.05,
  nrow(panels$d_proteome_pca$data) == sum(panels$a_design$data$mice),
  nrow(panels$f_signature_barcode$data) == 4
)

# One 2 x 3 grid: patchwork lines up titles, plot areas, axes and keys along each row and plot
# edges down each column.
figure <- wrap_plots(map(panels, "plot"), ncol = 3, heights = c(60, 68))
out <- stage_paths("04_Figures", "F01_Global")
save_figure(figure, out$reports, "F01_composite", width_mm = 180, height_mm = 128)

write_figure_workbook(
  panels,
  c(
    "Panel a: mice per group.",
    "Panel b: PAS vs vehicle per phenotype and sex, % change, 95% CI, Welch p, interaction p.",
    "Panel c: phenotype PC1 and PC2 per mouse, and each trait's loading.",
    "Panel d: proteome PC1 and PC2 per mouse, with the pellet's MitoCarta share.",
    "Panel e: per universe, contrast and direction, proteins at p or Π, Π and FDR ≤ 0.05.",
    "Panel f: each sex's p ≤ 0.05 up and down lists in the other sex, fry direction and p."
  ),
  file.path(out$data, "F01.xlsx")
)

sessionInfo()
