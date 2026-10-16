# F05 a: each WGCNA module eigengene against the three contrasts and the six phenotypes, after YvO's
# module-trait heatmap. Every cell holds the r-equivalent of its limma t, t / sqrt(t^2 + df), so a
# cell's value and its outline come from the same test: the contrast, or the eigengene's slope on
# the trait within one sex (pooled model, adjusted for the surrogate variable). Solid outline: BH
# FDR <= 0.05 over the modules; dashed: p <= 0.05 only; bold: |r| >= BOLD_R. Left, each module's
# size, MitoCarta share and stability (median best-match Jaccard over subsamples of mice, from
# 06_Modules). Pearson r across the mice of the sex is in the workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(
  fit = output_path("02_Statistics", "01_Fit", "fit.rds"),
  modules = output_path("01_Processing", "06_Modules", "modules.rds")
)
fit <- readRDS(inputs[["fit"]])
BOLD_R <- 0.4
sizes <- readRDS(inputs[["modules"]])$sizes |>
  filter(module != "grey") |>
  arrange(proteins) |>
  mutate(
    module_j = sprintf("%s (J %.2f)", module, median_jaccard),
    name = if_else(is.na(label) | label == "unannotated", module_j,
      sprintf("%s · %s", module_j, set_label(label, 30))
    ),
    row = factor(module, module)
  )

modules_only <- \(x) filter(x, level == "module")
evidence_levels <- c(fdr = paste("FDR ≤", ALPHA), p = paste("p ≤", ALPHA))
cells <- score_cells(modules_only(fit$score_contrasts), modules_only(fit$score_associations)) |>
  select(module = feature, section, column, r = r_equiv, P.Value, adj.P.Val, pearson_r) |>
  mutate(
    row = factor(module, levels(sizes$row)),
    evidence = case_when(
      adj.P.Val <= ALPHA ~ evidence_levels[["fdr"]], P.Value <= ALPHA ~ evidence_levels[["p"]],
      .default = "none"
    ),
    face = if_else(abs(r) >= BOLD_R, "bold", "plain")
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
    low = R_COLOURS[["low"]], high = R_COLOURS[["high"]], limits = c(-1, 1), breaks = c(-1, 0, 1),
    name = "r-equivalent",
    guide = guide_colourbar(barwidth = unit(22, "mm"), barheight = unit(2, "mm"))
  ) +
  scale_linetype_manual(
    values = set_names(c("solid", "22"), evidence_levels), name = NULL,
    guide = guide_legend(override.aes = list(fill = "white"))
  ) +
  scale_x_discrete(position = "top", expand = expansion(0)) +
  scale_y_discrete(expand = expansion(0)) +
  labs(x = NULL, y = NULL) +
  theme_figure() +
  theme(
    axis.text.y = element_blank(), panel.grid = element_blank(),
    axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 5.5),
    strip.text = element_text(size = 5.8), strip.placement = "outside",
    panel.spacing.x = unit(1.5, "mm")
  )

size_breaks <- keep(pretty(c(0, max(sizes$proteins)), n = 3), \(x) x <= max(sizes$proteins))
size_bar <- ggplot(sizes, aes(proteins, row)) +
  geom_col(aes(fill = I(module)), colour = "grey30", linewidth = 0.2, width = 0.75) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * mito_share)), hjust = 1.15, size = 5.3 / .pt) +
  scale_x_reverse(breaks = size_breaks, expand = expansion(mult = c(0.3, 0))) +
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
    subtitle = sprintf(
      "r-equivalent of each t; bold: absolute r ≥ %s; J: subsample stability", BOLD_R
    )
  ) +
  theme_figure()

save_panel(plot, "F05_Modules", "a_module_traits", 180, 80)
list(
  plot = plot, data = left_join(select(cells, -row), select(sizes, -row), by = "module"),
  inputs = inputs
)
