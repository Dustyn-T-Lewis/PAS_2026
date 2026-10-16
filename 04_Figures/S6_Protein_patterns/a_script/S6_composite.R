# Supplementary Figure 6, protein patterns: the proteins PAS changes in either sex (moderated F),
# by response pattern, and what each pattern holds. Each panel is its own script in panels/; this
# file runs them, checks the plotted counts against 01_fit.xlsx and 02_enrichment.xlsx, stitches
# them and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("a_ftest_heatmap", "b_pattern_ora")
panels <- run_panels("S6_Protein_patterns", PANELS)

fit_workbook <- output_path("02_Statistics", "01_Fit", "01_fit.xlsx")
ftest <- readxl::read_excel(fit_workbook, "protein_ftest") |>
  filter(universe == "full", P.Value <= ALPHA)
ora <- readxl::read_excel(
  output_path("02_Statistics", "02_Enrichment", "02_enrichment.xlsx"), "pattern_ora"
)
heatmap_rows <- panels$a_ftest_heatmap$data
stopifnot(
  setequal(heatmap_rows$feature, ftest$feature),
  identical(
    deframe(count(heatmap_rows, pattern = as.character(pattern))),
    deframe(count(ftest, pattern))
  ),
  isTRUE(all.equal(
    panels$b_pattern_ora$data$fdr,
    ora$fdr[match(
      paste(panels$b_pattern_ora$data$pattern, panels$b_pattern_ora$data$term),
      paste(ora$pattern, ora$term)
    )]
  ))
)

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

sessionInfo()
