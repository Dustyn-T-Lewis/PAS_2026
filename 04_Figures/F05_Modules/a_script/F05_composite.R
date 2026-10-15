# Figure 5, protein modules: every WGCNA module against the contrasts and the phenotypes, the
# modules that alone separate PAS from vehicle, and the strongest module-trait pairs mouse by
# mouse. Each panel is its own script in panels/; this file runs them, checks them against the
# statistics workbooks, stitches them and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("a_module_traits", "b_module_roc", "c_eigengene_scatter")
panels <- run_panels("F05_Modules", PANELS)

roc <- readxl::read_excel(
  output_path("02_Statistics", "03_Classification", "03_classification.xlsx"), "roc"
) |>
  filter(level == "module")
cells <- panels$a_module_traits$data
stopifnot(
  nrow(cells) == 14 * (length(COMPARISON_CONTRASTS) + 2 * length(TRAIT_SHORT)),
  !anyNA(cells$r), all(abs(cells$r) <= 1),
  setequal(
    paste(panels$b_module_roc$data$comparison, panels$b_module_roc$data$feature),
    with(filter(roc, p <= ALPHA), paste(comparison, feature))
  ),
  all(count(panels$c_eigengene_scatter$data, cell)$n == 10)
)

top <- panels$a_module_traits$plot
bottom <- wrap_plots(panels$b_module_roc$plot, panels$c_eigengene_scatter$plot, widths = c(110, 70))
# The bottom row is wrapped whole, so panel a does not line up with b's axis.
figure <- top / wrap_elements(full = bottom) + plot_layout(heights = c(80, 78))
out <- stage_paths("04_Figures", "F05_Modules")
save_figure(figure, out$reports, "F05_composite", width_mm = 180, height_mm = 160)

write_figure_workbook(
  panels,
  c(
    paste(
      "Panel a: per module and column, the r-equivalent of the limma t, p, BH FDR, and for traits",
      "Pearson r; module size, MitoCarta share and top annotated set."
    ),
    "Panel b: modules at exact Wilcoxon p ≤ 0.05 in either sex: AUC, DeLong CI, p, FDR.",
    "Panel c: eigengene and trait value per mouse for the four strongest module-trait pairs."
  ),
  file.path(out$data, "F05.xlsx")
)

sessionInfo()
