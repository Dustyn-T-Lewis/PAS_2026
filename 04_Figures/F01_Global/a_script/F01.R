# Figure 1, the global proteome: the design and the phenotypes on top, the proteome below. Each
# panel is its own script beside this one, run here in its own environment; this file stitches
# the returned plots and writes one workbook with a sheet per panel.
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
save_figure(figure, "F01_Global", "F01", width_mm = 180, height_mm = 128)
write_figure_workbook(panels, "F01_Global", "F01")
