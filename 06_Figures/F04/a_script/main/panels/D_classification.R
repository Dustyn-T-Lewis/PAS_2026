#!/usr/bin/env Rscript
# F04 d. Balanced error rate for telling PAS from vehicle in unseen mice (05_Summary/02_Panels).

source(here::here("R", "panels.R"))

PHENOTYPE_LABELS <- c(
  csa_mean = "CSA", capillary_density = "Capillary density",
  satellite_cells_per_fiber = "Satellite cells / fibre",
  conductance_sol_pyr = "Soleus pyruvate conductance",
  conductance_pla_pyr = "Plantaris pyruvate conductance",
  conductance_sol_oct = "Soleus octanoyl conductance"
)
MODEL_LABELS <- c(
  "pooled slope" = "Pooled", "females only" = "Female",
  "slope difference, male minus female" = "Δ slope"
)
INPUT_LABELS <- c(
  "phenotypes only" = "Phenotypes", "proteins features" = "Proteins",
  "sets features" = "Pathways", "modules features" = "Modules",
  "combined features" = "Proteome combined",
  "combined features + phenotypes" = "Proteome + phenotypes"
)
SEX_COLOURS_SHORT <- c(F = "#D55E00", M = "#0072B2")

inputs <- c(
  panels = here("05_Summary", "02_Panels", "c_data", "02_panels.xlsx")
)

performance <- readxl::read_excel(inputs[["panels"]], "performance") |>
  filter(comparison %in% c("PAS_vs_VEH_female", "PAS_vs_VEH_male")) |>
  mutate(input_key = if_else(
    level == "phenotypes", "phenotypes only", str_glue("{level} {input}")
  )) |>
  filter(input_key %in% names(INPUT_LABELS)) |>
  mutate(
    input = factor(INPUT_LABELS[input_key], rev(INPUT_LABELS)),
    sex = if_else(comparison == "PAS_vs_VEH_female", "Female", "Male"),
    method = if_else(str_starts(method, "DIABLO"), "DIABLO", "Elastic net"),
    label = if_else(p_value <= ALPHA, sprintf("%.2f *", ber), sprintf("%.2f", ber))
  )
panel <- ggplot(performance, aes(ber, input, fill = method)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = label),
    position = position_dodge(width = 0.75), hjust = -0.15,
    size = 1.6
  ) +
  facet_wrap(~sex) +
  scale_fill_manual(values = c(DIABLO = "grey35", "Elastic net" = "grey70"), name = NULL) +
  scale_x_continuous(limits = c(0, 0.75), breaks = c(0, 0.25, 0.5)) +
  labs(
    x = "Balanced error rate (lower is better; dashed: chance)", y = NULL,
    title = "Telling PAS from vehicle in unseen mice",
    subtitle = "Leave-one-out (DIABLO), leave-pair-out (elastic net); * exact permutation p ≤ 0.05"
  ) +
  theme_figure()

panel_data <- select(performance, comparison, input, method, ber, null_median, p_value)
panel_note <- "Balanced error rate and permutation p per input and method (05_Summary/02_Panels)."
save_panel(panel, "F04", "main", "D_classification", 120, 75)
