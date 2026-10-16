# F02 d: every mouse for the male candidates passing three or more routes, as log2 abundance
# against the own-sex vehicle mean in all four groups; open circles mark mice in which the protein
# was not detected.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(
  candidates = output_path("02_Statistics", "06_Candidates", "candidates.rds"),
  fit = output_path("02_Statistics", "01_Fit", "fit.rds"),
  proteins = output_path("01_Processing", "02_Quantification", "proteins.rds")
)
cand <- readRDS(inputs[["candidates"]])$candidates
fit <- readRDS(inputs[["fit"]])
quant <- readRDS(inputs[["proteins"]])
targets <- as_tibble(fit$targets)
shown <- cand |>
  filter(sex == "male", routes >= 3) |>
  arrange(desc(routes), P.Value)
values <- fit$E_norm$full[shown$feature, targets$sample_id, drop = FALSE]
detected <- quant$proteins$other$n.observations[shown$feature, targets$sample_id, drop = FALSE] > 0
labels <- sprintf("%s (%d)", shown$gene, shown$routes) |>
  make.unique(sep = " #") |>
  set_names(shown$feature)

plot <- candidate_dots(values, detected, targets, labels,
  title = "Male candidates, every mouse",
  subtitle = "Proteins passing three or more routes (routes in brackets); bar: group mean"
)

save_panel(plot, "F02_Proteins", "d_male_candidates", 180, 62)
list(
  plot = plot,
  data = select(
    shown, gene, logFC, P.Value, adj.P.Val, lfsr, auc, routes, detected_veh, detected_pas
  ),
  inputs = inputs
)
