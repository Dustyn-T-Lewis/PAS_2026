# F05 (complexes version) a: every CORUM complex and MitoCarta pathway whose singscore differs at
# p <= 0.05 in at least one of the three contrasts, against the contrasts and against each trait
# within each sex (pooled model, adjusted for the surrogate variable). Dot fill: r-equivalent of
# the limma t, so contrasts and traits share one scale; size: -log10 p; ring: BH FDR <= 0.05 over
# every scored set in the universe. Labels give the measured members.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(
  fit = output_path("02_Statistics", "01_Fit", "fit.rds"),
  gene_sets = output_path("01_Processing", "04_Gene_Sets", "gene_sets.rds")
)
fit <- readRDS(inputs[["fit"]])
gs <- readRDS(inputs[["gene_sets"]])

in_scope <- \(x) filter(x, level == "complex" | collection == "MitoCarta")
scored <- score_cells(in_scope(fit$score_contrasts), in_scope(fit$score_associations))
kept <- scored |>
  filter(section == "Contrasts") |>
  filter(any(P.Value <= ALPHA), .by = feature) |>
  distinct(level, feature)
measured <- gs$membership |>
  filter(collection %in% c("CORUM", "MitoCarta")) |>
  count(term, name = "members")
rows <- scored |>
  filter(column == CONTRAST_SHORT[["PAS_vs_VEH_female"]]) |>
  semi_join(kept, by = "feature") |>
  left_join(measured, by = c("feature" = "term")) |>
  mutate(
    group = if_else(level == "complex", "CORUM complexes", "MitoCarta pathways"),
    label = sprintf("%s (%d)", set_label(str_remove(feature, "_\\d+$"), 46), members)
  ) |>
  arrange(group, r_equiv) |>
  mutate(row = factor(feature, rev(feature)))
cells <- scored |>
  inner_join(select(rows, feature, row, group, label, members), by = "feature") |>
  mutate(log_p = pmin(-log10(P.Value), LOG_CAP))

plot <- ggplot(cells, aes(column, row)) +
  geom_point(
    aes(size = log_p, fill = r_equiv, stroke = adj.P.Val <= ALPHA),
    shape = 21, colour = "grey15"
  ) +
  ggh4x::facet_grid2(group ~ section,
    scales = "free", space = "free", switch = "y", labeller = labeller(group = label_wrap_gen(12))
  ) +
  scale_y_discrete(labels = deframe(distinct(cells, row, label))) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient2(
    low = R_COLOURS[["low"]], high = R_COLOURS[["high"]],
    limits = c(-1, 1), breaks = c(-1, 0, 1), name = "r-equivalent"
  ) +
  fdr_scales(name = "-log10 p") +
  guides(
    fill = colour_bar(order = 1), size = guide_legend(order = 2),
    stroke = guide_legend(order = 3, override.aes = list(size = 2.4, fill = "white"))
  ) +
  labs(
    x = NULL, y = NULL, title = "Protein complexes and mitochondrial pathways",
    subtitle = sprintf("Singscore p ≤ %s in a contrast; measured members in brackets", ALPHA)
  ) +
  theme_figure() +
  theme(
    axis.text.y = element_text(size = 5.2),
    axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 5.5),
    strip.placement = "outside", strip.text.y.left = element_text(angle = 90),
    strip.text = element_text(size = 5.5), panel.spacing = unit(1.2, "mm")
  )

save_panel(plot, "F05_Complexes", "a_set_matrix", 180, 120)
list(
  plot = plot,
  data = cells |>
    select(group, feature, label, members, section, column, r_equiv, P.Value, adj.P.Val) |>
    mutate(across(c(section, column), as.character)),
  inputs = inputs
)
