# Figure 1, the global proteome: the design and the phenotypes on top, the proteome below. Each
# panel is its own script in panels/, run here in its own environment; this file stitches the
# returned plots and writes one workbook with a sheet per panel.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c(
  "a_design", "b_phenotype_pca", "c_proteome_pca", "d_phenotype_effects", "e_protein_bars",
  "f_signature_barcode"
)
panels <- run_panels("F01_Global", PANELS)

# One 2 x 3 grid: patchwork lines up titles, plot areas, axes and keys along each row and plot
# edges down each column. The design image is wrapped whole so it fills its cell rather than the
# short plot area the PCA keys leave.
plots <- map(panels, "plot")
plots$a_design <- wrap_elements(full = plots$a_design) + theme(plot.margin = margin(t = 0))
figure <- wrap_plots(plots, ncol = 3, widths = c(66, 57, 57), heights = c(64, 64))
out <- stage_paths("04_Figures", "F01_Global")
save_figure(figure, out$reports, "F01_composite", width_mm = 180, height_mm = 128)

write_figure_workbook(
  panels,
  c(
    "Panel a: mice per group.",
    "Panel b: phenotype PC1 and PC2 per mouse, and each trait's loading.",
    "Panel c: proteome PC1 and PC2 per mouse, with the pellet's MitoCarta share.",
    "Panel d: PAS vs vehicle per phenotype and sex, % change, 95% CI, Welch p, interaction p.",
    "Panel e: per universe, contrast and direction, proteins at p ≤ 0.05, Π ≤ 0.05, FDR ≤ 0.05.",
    "Panel f: each sex's Π ≤ 0.05 up and down lists in the other sex, fry direction and p."
  ),
  file.path(out$data, "F01.xlsx")
)
