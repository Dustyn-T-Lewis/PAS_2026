# Supplementary Figure 8, sensitivity: how the protein results move when one processing or model
# choice changes. Each panel is its own script in panels/; this file runs them, checks them against
# 07_sensitivity.xlsx, stitches them and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("a_counts", "b_named")
panels <- run_panels("S8_Sensitivity", PANELS)

workbook <- output_path("02_Statistics", "07_Sensitivity", "07_sensitivity.xlsx")
counts <- readxl::read_excel(workbook, "counts")
named <- readxl::read_excel(workbook, "named")
stopifnot(
  all(panels$a_counts$data$p_le_0.05 == counts$p_le_0.05),
  nrow(panels$b_named$data) == nrow(named),
  isTRUE(all.equal(panels$b_named$data$adj.P.Val, named$adj.P.Val))
)

figure <- panels$a_counts$plot / panels$b_named$plot + plot_layout(heights = c(62, 100))
out <- stage_paths("04_Figures", "S8_Sensitivity")
save_figure(figure, out$reports, "S8_composite", width_mm = 180, height_mm = 162)

write_figure_workbook(
  panels,
  c(
    paste(
      "Panel a: per variant and contrast, proteins at p ≤ 0.05, FDR ≤ 0.10 and FDR ≤ 0.05, and how",
      "many of the baseline's 50 smallest p-values stay in the variant's 50 smallest."
    ),
    "Panel b: each protein at FDR ≤ 0.10 in the baseline: log2 FC, p and FDR under each variant."
  ),
  file.path(out$data, "S8.xlsx")
)

sessionInfo()
