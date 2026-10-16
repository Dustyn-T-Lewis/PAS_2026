# Figure 2, candidate proteins: within each sex, the proteins passing two or more of the six
# candidate routes as a STRING network, and every mouse for those passing three or more. Each
# panel is its own script in panels/; this file runs them, checks them against
# 06_candidates.xlsx, stitches them and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("a_female_network", "b_male_network", "c_female_candidates", "d_male_candidates")
panels <- run_panels("F02_Proteins", PANELS)

candidates <- readxl::read_excel(
  output_path("02_Statistics", "06_Candidates", "06_candidates.xlsx"), "candidates"
)
expected <- \(s, k) sum(candidates$sex == s & candidates$routes >= k)
stopifnot(
  nrow(panels$a_female_network$data) == expected("female", 2),
  nrow(panels$b_male_network$data) == expected("male", 2),
  nrow(panels$c_female_candidates$data) == expected("female", 3),
  nrow(panels$d_male_candidates$data) == expected("male", 3)
)

top <- wrap_plots(panels$a_female_network$plot, panels$b_male_network$plot, ncol = 2)
figure <- wrap_elements(full = top) /
  panels$c_female_candidates$plot /
  panels$d_male_candidates$plot +
  plot_layout(heights = c(110, 62, 62))
out <- stage_paths("04_Figures", "F02_Proteins")
save_figure(figure, out$reports, "F02_composite", width_mm = 180, height_mm = 240)

write_figure_workbook(
  panels,
  c(
    "Panel a: female candidates passing two or more routes: routes, STRING community and name.",
    "Panel b: male candidates passing two or more routes: routes, STRING community and name.",
    "Panel c: female candidates passing three or more routes: effect, p, FDR, lfsr, detection.",
    "Panel d: male candidates passing three or more routes: effect, p, FDR, lfsr, detection."
  ),
  file.path(out$data, "F02.xlsx")
)

sessionInfo()
