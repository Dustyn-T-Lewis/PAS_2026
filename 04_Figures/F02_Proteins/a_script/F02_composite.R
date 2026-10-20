# Figure 2, candidate proteins: within each sex, the candidates 06_Candidates grouped into STRING
# communities as a network, and every mouse for those passing CANDIDATE_DOT_ROUTES or more routes.
# Each panel is its own script in panels/; this file runs them, stitches them and writes one
# workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("a_female_network", "b_male_network", "c_female_candidates", "d_male_candidates")
panels <- run_panels("F02_Proteins", PANELS)

top <- wrap_plots(panels$a_female_network$plot, panels$b_male_network$plot, ncol = 2)
# c and d share one group key, collected under d.
figure <- wrap_elements(full = top) /
  panels$c_female_candidates$plot /
  panels$d_male_candidates$plot +
  plot_layout(heights = c(90, 62, 62), guides = "collect") &
  theme(legend.position = "bottom")
out <- stage_paths("04_Figures", "F02_Proteins")
save_figure(figure, out$reports, "F02_composite", width_mm = 180, height_mm = 220)

thresholds <- readRDS(output_path("02_Statistics", "06_Candidates", "candidates.rds"))$thresholds
network_hold <- sprintf(
  "candidates passing %d or more routes: routes, STRING community, name and has_partner.",
  thresholds$community_routes
)
dots_hold <- sprintf(
  "candidates passing %d or more routes: effect, p, FDR, lfsr, detection.", CANDIDATE_DOT_ROUTES
)
write_figure_workbook(
  panels,
  c(
    paste("Panel a: female", network_hold), paste("Panel b: male", network_hold),
    paste("Panel c: female", dots_hold), paste("Panel d: male", dots_hold)
  ),
  file.path(out$data, "F02.xlsx")
)
