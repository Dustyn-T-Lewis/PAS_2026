# S1 c: the label-blind normalisation benchmark from 02_Quantification. Within-group spread
# (NormalyzerDE's PCV, PMAD, PEV; lower is tighter) per method and universe; only methods that keep
# limpa's standard errors were eligible, and the chosen one is filled.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(
  quantification = output_path("01_Processing", "02_Quantification", "proteins.rds")
)
benchmark <- readRDS(inputs[["quantification"]])$normalization |>
  pivot_longer(c(pcv, pmad, pev), names_to = "measure", values_to = "spread") |>
  mutate(
    method = factor(method, rev(unique(method))),
    measure = factor(toupper(measure), c("PCV", "PMAD", "PEV")),
    universe = factor(if_else(universe == "full", "Full", "Mito"), c("Full", "Mito")),
    status = case_when(chosen ~ "chosen", eligible ~ "eligible", .default = "changes SEs")
  )

plot <- ggplot(benchmark, aes(spread, method)) +
  geom_point(aes(shape = status, colour = status), size = 1.4, stroke = 0.5) +
  facet_grid(universe ~ measure, scales = "free_x") +
  scale_shape_manual(values = c(chosen = 16, eligible = 1, `changes SEs` = 4), name = NULL) +
  scale_colour_manual(
    values = c(chosen = "black", eligible = "black", `changes SEs` = "grey60"), name = NULL
  ) +
  scale_x_continuous(breaks = scales::breaks_pretty(2), expand = expansion(mult = 0.3)) +
  labs(
    x = "within-group spread", y = NULL, title = "Normalisation",
    subtitle = "Lower is tighter; chosen among SE-preserving methods"
  ) +
  theme_figure() +
  theme(strip.text = element_text(size = 5.5), panel.spacing.x = unit(3, "mm"))

save_figure(plot, "F01_Global", "S1c_normalisation", 55, 55)
list(plot = plot, data = benchmark)
