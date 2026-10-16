# Supplementary Figure 8, sensitivity: how the protein results move when one processing or model
# choice changes. Each panel is its own script in panels/; this file runs them, stitches them and
# writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("a_counts", "b_named")
top_n <- readRDS(output_path("02_Statistics", "07_Sensitivity", "sensitivity.rds"))$top_n
panels <- run_panels("S8_Sensitivity", PANELS)

figure <- panels$a_counts$plot / panels$b_named$plot + plot_layout(heights = c(62, 100))
out <- stage_paths("04_Figures", "S8_Sensitivity")
save_figure(figure, out$reports, "S8_composite", width_mm = 180, height_mm = 162)

write_figure_workbook(
  panels,
  c(
    str_glue(
      "Panel a: per variant and contrast, proteins at p ≤ 0.05, FDR ≤ 0.10 and FDR ≤ 0.05, and ",
      "how many of the baseline's {top_n} smallest p-values stay in the variant's {top_n} smallest."
    ),
    "Panel b: each protein at FDR ≤ 0.10 in the baseline: log2 FC, p and FDR under each variant."
  ),
  file.path(out$data, "S8.xlsx")
)
