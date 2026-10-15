# S2 a: female PAS vs vehicle, the proteins at Π <= 0.05 in the full proteome across all twenty
# mice, centred on each sex's vehicle mean and split by direction (01_Fit).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(fit = output_path("02_Statistics", "01_Fit", "fit.rds"))
panel <- pi_heatmap_panel(readRDS(inputs[["fit"]]), "PAS_vs_VEH_female", "PAS vs vehicle, females",
  width_mm = 50, legend = FALSE
)

save_panel(panel$plot, "S2_Proteins", "a_pi_female", 50, 130)
list(plot = panel$plot, data = panel$data, inputs = inputs)
