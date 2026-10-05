#!/usr/bin/env Rscript
# F05 b. WGCNA modules: top annotated term, PAS contrasts (emmeans) and within-group phenotype
# slopes (04_Network).

source(here::here("R", "panels.R"))

PHENOTYPE_SHORT <- c(
  conductance_sol_pyr = "Soleus pyruvate", capillary_density = "Capillary density",
  csa_mean = "CSA", conductance_sol_oct = "Soleus octanoyl"
)
MODEL_SHORT <- c("females only" = "F", "slope difference, male minus female" = "Δ")
term_label <- \(x, n = 40) str_trunc(enrichVolcano::clean_label(x, width = 200), n)

inputs <- c(
  modules = here("04_Network", "02_Contrasts", "c_data", "02_contrasts.xlsx"),
  module_traits = here("04_Network", "04_Associate", "c_data", "04_associate.xlsx"),
  characterise = here("04_Network", "05_Characterise", "c_data", "05_characterise.xlsx")
)
read_sheet <- \(file, sheet) readxl::read_excel(inputs[[file]], sheet, guess_max = 1e5)

overview <- read_sheet("characterise", "overview") |>
  mutate(
    identity = term_label(str_split_i(top_terms, ";", 1), 30),
    row = str_glue("{module} ({proteins}): {identity}")
  )
module_effects <- read_sheet("modules", "results") |>
  filter(contrast %in% PRIMARY) |>
  transmute(module, column = CONTRAST_LABELS[contrast], t, p = P.Value)
module_traits <- read_sheet("module_traits", "association") |>
  filter(model == "pooled slope", phenotype %in% names(PHENOTYPE_SHORT)) |>
  transmute(module, column = PHENOTYPE_SHORT[phenotype], t, p = P.Value)
module_data <- bind_rows(module_effects, module_traits) |>
  left_join(select(overview, module, row), by = "module") |>
  mutate(
    column = factor(column, c(CONTRAST_LABELS[PRIMARY], PHENOTYPE_SHORT)),
    block = if_else(column %in% CONTRAST_LABELS, "PAS effect", "Tracks phenotype (within group)"),
    block = factor(block, c("PAS effect", "Tracks phenotype (within group)"))
  )
panel <- ggplot(module_data, aes(column, row, fill = t)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_point(
    data = \(d) filter(d, p <= ALPHA), shape = 21, fill = NA, colour = "black",
    size = 2.6, stroke = 0.6
  ) +
  facet_grid(~block, scales = "free_x", space = "free_x") +
  scale_fill_distiller(palette = "RdBu", limits = c(-3, 3), oob = scales::squish) +
  labs(
    x = NULL, y = NULL, title = "WGCNA modules: identity, PAS effect, phenotype link",
    subtitle = paste(
      "Module (proteins): top annotated term; t from emmeans contrasts or within-group",
      "\nslope; ring: p ≤ 0.05"
    )
  ) +
  theme_figure() +
  theme(
    axis.text.x = element_text(size = 5, angle = 35, hjust = 1),
    axis.text.y = element_text(size = 5.5),
    panel.grid = element_blank(), strip.text = element_text(size = 6)
  )

panel_data <- module_data
panel_note <- "Module PAS contrasts and within-group phenotype slopes (04_Network)."
save_panel(panel, "F05", "main", "B_modules", 140, 90)
