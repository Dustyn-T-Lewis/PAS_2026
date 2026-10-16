# S2 b: male PAS vs vehicle, the proteins at Π <= 0.05 in the full proteome across every mouse,
# centred on each sex's vehicle mean and split by direction (01_Fit).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(fit = output_path("02_Statistics", "01_Fit", "fit.rds"))
panel <- pi_heatmap_panel(
  readRDS(inputs[["fit"]]), "PAS_vs_VEH_male",
  width_mm = 50, legend = FALSE
)

save_panel(panel$plot, "S2_Proteins", "b_pi_male", panel$width_mm, panel$height_mm)
list(plot = panel$plot, data = panel$data, inputs = inputs)
