# F05 a: each WGCNA module eigengene against the three contrasts and the six phenotypes, after YvO's
# module-trait heatmap. Every cell holds the r-equivalent of its limma t, t / sqrt(t^2 + df), so a
# cell's value and its outline come from the same test: the contrast, or the eigengene's slope on
# the trait within one sex (pooled model, adjusted for the surrogate variable). Solid outline: BH
# FDR <= 0.05 over the 14 modules; dashed: p <= 0.05 only; bold: |r| >= 0.4. Left, each module's
# size and MitoCarta share. Pearson r across the ten mice is in the workbook beside each cell.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(
  fit = output_path("02_Statistics", "01_Fit", "fit.rds"),
  modules = output_path("01_Processing", "06_Modules", "modules.rds")
)
fit <- readRDS(inputs[["fit"]])
sizes <- readRDS(inputs[["modules"]])$sizes |>
  filter(module != "grey") |>
  arrange(proteins) |>
  mutate(
    name = if_else(is.na(label) | label == "unannotated", module,
      sprintf("%s · %s", module, set_label(label, 30))
    ),
    row = factor(module, module)
  )

contrast_cells <- fit$score_contrasts |>
  filter(level == "module", contrast %in% COMPARISON_CONTRASTS) |>
  transmute(
    module = feature, section = "Contrasts", column = CONTRAST_SHORT[contrast], r = r_equiv,
    P.Value, adj.P.Val
  )
trait_cells <- fit$score_associations |>
  filter(level == "module", model == "pooled") |>
  transmute(
    module = feature, section = if_else(sex == "female", "Female traits", "Male traits"),
    column = TRAIT_SHORT[phenotype], r = r_equiv, pearson_r = r, P.Value, adj.P.Val
  )
cells <- bind_rows(contrast_cells, trait_cells) |>
  mutate(
    row = factor(module, levels(sizes$row)),
    section = factor(section, c("Contrasts", "Female traits", "Male traits")),
    column = factor(column, c(CONTRAST_SHORT, TRAIT_SHORT)),
    evidence = case_when(
      adj.P.Val <= ALPHA ~ "FDR ≤ 0.05", P.Value <= ALPHA ~ "p ≤ 0.05",
      .default = "none"
    ),
    face = if_else(abs(r) >= 0.4, "bold", "plain")
  )

heatmap <- ggplot(cells, aes(column, row)) +
  geom_tile(aes(fill = r), colour = "grey85", linewidth = 0.2) +
  geom_tile(
    aes(linetype = evidence), filter(cells, evidence != "none"),
    fill = NA, colour = "black", linewidth = 0.45, width = 0.92, height = 0.88
  ) +
  geom_text(aes(label = sprintf("%.2f", r), fontface = face), size = 1.55) +
  facet_grid(~section, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = "#8C510A", high = "#01665E", limits = c(-1, 1), breaks = c(-1, 0, 1),
    name = "r-equivalent",
    guide = guide_colourbar(barwidth = unit(22, "mm"), barheight = unit(2, "mm"))
  ) +
  scale_linetype_manual(
    values = c("FDR ≤ 0.05" = "solid", "p ≤ 0.05" = "22"), name = NULL,
    guide = guide_legend(override.aes = list(fill = "white"))
  ) +
  scale_x_discrete(position = "top", expand = expansion(0)) +
  scale_y_discrete(expand = expansion(0)) +
  labs(x = NULL, y = NULL) +
  theme_figure() +
  theme(
    axis.text.y = element_blank(), panel.grid = element_blank(),
    axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 5.5),
    strip.text = element_text(face = "bold", size = 5.8), strip.placement = "outside",
    panel.spacing.x = unit(1.5, "mm")
  )

size_bar <- ggplot(sizes, aes(proteins, row)) +
  geom_col(aes(fill = I(module)), colour = "grey30", linewidth = 0.2, width = 0.75) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * mito_share)), hjust = 1.15, size = 1.55) +
  scale_x_reverse(expand = expansion(mult = c(0, 0.25))) +
  scale_y_discrete(labels = deframe(select(sizes, row, name)), position = "right") +
  coord_cartesian(clip = "off") +
  labs(x = "proteins (% MitoCarta)", y = NULL) +
  theme_figure() +
  theme(
    panel.grid.major.y = element_blank(), axis.text.y.right = element_text(size = 5.3, hjust = 0)
  )

plot <- wrap_elements(plot = (size_bar | heatmap) + plot_layout(widths = c(1, 4.2))) +
  labs(
    title = "Modules against contrasts and phenotypes",
    subtitle = "r-equivalent of each limma t; traits within sex, adjusted for SV1; bold |r| ≥ 0.4"
  ) +
  theme_figure()

save_panel(plot, "F05_Modules", "a_module_traits", 180, 80)
list(
  plot = plot, data = left_join(select(cells, -row), select(sizes, -row), by = "module"),
  inputs = inputs
)
