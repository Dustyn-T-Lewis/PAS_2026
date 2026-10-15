# Figure 3, contrast pathways: the full-proteome collections and MitoCarta, each set at fgsea FDR
# <= 0.05 that collapsePathways keeps, under the contrast that first selects it. Each panel is its
# own script in panels/; this file runs them, checks the rows against 02_enrichment.xlsx, stitches
# them with one shared key and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("a_full_matrix", "b_mitocarta_matrix")
panels <- run_panels("F03_Pathways", PANELS)

fgsea <- readxl::read_excel(
  output_path("02_Statistics", "02_Enrichment", "02_enrichment.xlsx"), "fgsea"
)
kept <- filter(fgsea, padj <= ALPHA, dedup_status == "kept")
plotted <- \(panel) unique(panels[[panel]]$data$term)
stopifnot(
  setequal(plotted("a_full_matrix"), kept$term[kept$universe == "full"]),
  setequal(plotted("b_mitocarta_matrix"), kept$term[kept$universe == "mito"]),
  all(count(bind_rows(map(panels, "data")), term)$n == length(COMPARISON_CONTRASTS))
)

figure <- wrap_plots(
  a = panels$a_full_matrix$plot,
  b = panels$b_mitocarta_matrix$plot + guides(colour = "none"), key = guide_area(),
  design = "ab\nak", widths = c(110, 70), heights = c(11, 25)
) +
  plot_layout(guides = "collect")
out <- stage_paths("04_Figures", "F03_Pathways")
save_figure(figure, out$reports, "F03_composite", width_mm = 180, height_mm = 150)

write_figure_workbook(
  panels,
  c(
    "Panel a: each full-proteome set shown, per contrast: NES, FDR, collection, group.",
    "Panel b: each MitoCarta set shown, in every contrast: NES, FDR, selecting contrast."
  ),
  file.path(out$data, "F03.xlsx")
)

sessionInfo()
