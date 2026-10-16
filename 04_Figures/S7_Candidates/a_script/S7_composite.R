# Supplementary Figure 7, candidate routes: for each sex, every candidate passing community_routes
# or more routes and which of the six routes it passed. The full ranked table of every
# candidate, one route or more, is the workbook's last sheet. Each panel is its own script in
# panels/; this file runs them, stitches them and writes one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("a_female_routes", "b_male_routes")
panels <- run_panels("S7_Candidates", PANELS)

figure <- wrap_plots(map(panels, "plot"), ncol = 2)
out <- stage_paths("04_Figures", "S7_Candidates")
save_figure(figure, out$reports, "S7_composite", width_mm = 180, height_mm = 165)

candidates <- output_path("02_Statistics", "06_Candidates", "candidates.rds")
cand <- readRDS(candidates)
panels$all_candidates <- list(
  data = arrange(cand$candidates, sex, desc(routes), P.Value), inputs = candidates
)
routes_hold <- sprintf(
  "candidates passing %d or more routes: each route, log2 FC, p, lfsr, set.",
  cand$thresholds$community_routes
)
write_figure_workbook(
  panels,
  c(
    paste("Panel a: female", routes_hold), paste("Panel b: male", routes_hold),
    "Every candidate in either sex passing one route or more, ranked by routes and then p."
  ),
  file.path(out$data, "S7.xlsx")
)
