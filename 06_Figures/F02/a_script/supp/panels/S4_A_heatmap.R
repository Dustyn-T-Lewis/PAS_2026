#!/usr/bin/env Rscript
# S4 a. Every protein at p ≤ 0.05 for PAS in either sex or for the interaction, z-scored within sex
# so the PAS response shows rather than the sex difference; rows clustered; a strip marks which
# contrast selected each protein.

source(here::here("R", "panels.R"))

SELECTORS <- c(
  PAS_vs_VEH_female = "Female PAS", PAS_vs_VEH_male = "Male PAS",
  treatment_by_sex = "Interaction"
)

inputs <- c(
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
de <- readRDS(inputs[["contrasts"]])
targets <- as_tibble(de$targets) |>
  mutate(group = factor(group, FIGURE_GROUPS)) |>
  arrange(group)

selected <- de$results |>
  filter(universe == "full", contrast %in% names(SELECTORS), P.Value <= ALPHA) |>
  distinct(protein, contrast) |>
  mutate(hit = TRUE) |>
  pivot_wider(names_from = contrast, values_from = hit, values_fill = FALSE)
m <- de$E_norm$full[selected$protein, targets$sample_id]
centred <- m
for (s in unique(targets$sex)) {
  cols <- targets$sex == s
  centred[, cols] <- t(scale(t(m[, cols])))
}
limit <- unname(quantile(abs(centred), 0.99))
row_order <- rownames(centred)[cluster_order(centred)]
tile_data <- as_tibble(centred, rownames = "protein") |>
  pivot_longer(-protein, names_to = "sample_id", values_to = "z") |>
  left_join(select(targets, sample_id, group), by = "sample_id") |>
  mutate(
    protein = factor(protein, rev(row_order)),
    sample_id = factor(sample_id, targets$sample_id),
    z = pmax(pmin(z, limit), -limit)
  )
marks <- selected |>
  pivot_longer(-protein, names_to = "contrast", values_to = "hit") |>
  filter(hit) |>
  mutate(
    protein = factor(protein, rev(row_order)),
    selector = factor(SELECTORS[contrast], SELECTORS)
  )
panel_marks <- ggplot(marks, aes(selector, protein, fill = selector)) +
  geom_tile() +
  scale_fill_manual(
    values = set_names(CONTRAST_COLOURS[CONTRAST_LABELS[names(SELECTORS)]], SELECTORS),
    guide = "none"
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 7) +
  theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5))
panel_tiles <- ggplot(tile_data, aes(sample_id, protein, fill = z)) +
  geom_raster() +
  facet_grid(~group, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-limit, limit),
    name = "z\n(within sex)"
  ) +
  labs(x = NULL, y = NULL) +
  theme_figure() +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
    panel.spacing = unit(0.8, "mm"), strip.text = element_text(size = 6)
  )
panel_tiles <- panel_tiles +
  labs(
    title = "PAS-responsive proteins",
    subtitle = str_glue(
      "{nrow(selected)} proteins at p ≤ 0.05 in female PAS, male PAS or the interaction\n",
      "(left strip); centred within sex; rows clustered (1 − Pearson r, average)"
    )
  )
panel <- panel_marks + panel_tiles + plot_layout(widths = c(0.08, 1))
panel_data <- left_join(selected, distinct(de$results, protein, gene), by = "protein")
panel_note <- "Proteins in the heatmap and the contrasts (p ≤ 0.05) that selected them."
save_panel(panel, "F02", "supp", "S4_A_heatmap", 140, 160)
