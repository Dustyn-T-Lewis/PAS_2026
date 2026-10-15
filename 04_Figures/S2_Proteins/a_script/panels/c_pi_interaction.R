# S2 c: Treatment × sex, the proteins at Π <= 0.05 in the full proteome across all twenty mice,
# centred on each sex's vehicle mean and split by direction (01_Fit).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(fit = output_path("02_Statistics", "01_Fit", "fit.rds"))
panel <- pi_heatmap_panel(readRDS(inputs[["fit"]]), "treatment_by_sex", "Treatment × sex",
  width_mm = 74
)

save_panel(panel$plot, "S2_Proteins", "c_pi_interaction", 74, 130)
list(plot = panel$plot, data = panel$data, inputs = inputs)
