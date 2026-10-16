# Supplementary Figure 7, candidate routes: for each sex, every protein passing two or more of the
# six candidate routes and which ones it passed. The full ranked table of every candidate, one route
# or more, is the workbook's last sheet. Each panel is its own script in panels/; this file runs
# them, checks them against 06_candidates.xlsx, stitches them and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("a_female_routes", "b_male_routes")
panels <- run_panels("S7_Candidates", PANELS)

candidates <- readxl::read_excel(
  output_path("02_Statistics", "06_Candidates", "06_candidates.xlsx"), "candidates"
)
expected <- \(s) sum(candidates$sex == s & candidates$routes >= 2)
stopifnot(
  nrow(panels$a_female_routes$data) == expected("female"),
  nrow(panels$b_male_routes$data) == expected("male"),
  all(map_lgl(panels, \(p) {
    all(p$data$routes == rowSums(select(p$data, all_of(names(CANDIDATE_ROUTES)))))
  }))
)

figure <- wrap_plots(map(panels, "plot"), ncol = 2)
out <- stage_paths("04_Figures", "S7_Candidates")
save_figure(figure, out$reports, "S7_composite", width_mm = 180, height_mm = 165)

panels$all_candidates <- list(
  data = candidates |> arrange(sex, desc(routes), P.Value),
  inputs = output_path("02_Statistics", "06_Candidates", "06_candidates.xlsx")
)
write_figure_workbook(
  panels,
  c(
    "Panel a: female candidates passing two or more routes: each route, log2 FC, p, lfsr, set.",
    "Panel b: male candidates passing two or more routes: each route, log2 FC, p, lfsr, set.",
    "Every candidate in either sex passing one route or more, ranked by routes and then p."
  ),
  file.path(out$data, "S7.xlsx")
)

sessionInfo()
