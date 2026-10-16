# F05 (complexes version) a: every CORUM complex and MitoCarta pathway whose singscore differs at
# p <= 0.05 in at least one of the three contrasts, against the contrasts and against each trait
# within each sex (pooled model, adjusted for the surrogate variable). Dot fill: r-equivalent of
# the limma t, so contrasts and traits share one scale; size: -log10 p; ring: BH FDR <= 0.05 over
# every scored set in the universe (1,066 full, 44 MitoCarta). Labels give the measured members.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(
  fit = output_path("02_Statistics", "01_Fit", "fit.rds"),
  gene_sets = output_path("01_Processing", "04_Gene_Sets", "gene_sets.rds")
)
fit <- readRDS(inputs[["fit"]])
gs <- readRDS(inputs[["gene_sets"]])
P_CAP <- 4
SECTIONS <- c("Contrasts", "Female traits", "Male traits")

in_scope <- \(x) filter(x, level == "complex" | collection == "MitoCarta")
contrast_cells <- fit$score_contrasts |>
  in_scope() |>
  filter(contrast %in% COMPARISON_CONTRASTS) |>
  transmute(
    level, feature,
    section = "Contrasts", column = CONTRAST_SHORT[contrast], r_equiv,
    P.Value, adj.P.Val
  )
kept <- contrast_cells |>
  filter(any(P.Value <= ALPHA), .by = feature) |>
  distinct(level, feature)
trait_cells <- fit$score_associations |>
  in_scope() |>
  filter(model == "pooled") |>
  transmute(
    level, feature,
    section = if_else(sex == "female", "Female traits", "Male traits"),
    column = TRAIT_SHORT[phenotype], r_equiv, P.Value, adj.P.Val
  )

measured <- gs$membership |>
  filter(collection %in% c("CORUM", "MitoCarta")) |>
  count(term, name = "members")
rows <- contrast_cells |>
  filter(column == CONTRAST_SHORT[["PAS_vs_VEH_female"]]) |>
  semi_join(kept, by = "feature") |>
  left_join(measured, by = c("feature" = "term")) |>
  mutate(
    group = if_else(level == "complex", "CORUM complexes", "MitoCarta pathways"),
    name = enrichVolcano::clean_label(str_remove(feature, "_\\d+$"), width = 200),
    label = sprintf("%s (%d)", str_trunc(name, 46), members)
  ) |>
  arrange(group, r_equiv) |>
  mutate(row = factor(feature, rev(feature)))

cells <- bind_rows(contrast_cells, trait_cells) |>
  inner_join(select(rows, feature, row, group, label, members), by = "feature") |>
  mutate(
    section = factor(section, SECTIONS),
    column = factor(column, c(CONTRAST_SHORT, TRAIT_SHORT)),
    log_p = pmin(-log10(P.Value), P_CAP)
  )

plot <- ggplot(cells, aes(column, row)) +
  geom_point(
    aes(size = log_p, fill = r_equiv, stroke = adj.P.Val <= ALPHA),
    shape = 21, colour = "grey15"
  ) +
  ggh4x::facet_grid2(group ~ section,
    scales = "free", space = "free", switch = "y"
  ) +
  scale_y_discrete(labels = deframe(distinct(cells, row, label))) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient2(
    low = DIRECTION_COLOURS[["Down"]], high = DIRECTION_COLOURS[["Up"]],
    limits = c(-1, 1), breaks = c(-1, 0, 1), name = "r-equivalent",
    guide = guide_colourbar(order = 1, barwidth = unit(14, "mm"), barheight = unit(2, "mm"))
  ) +
  scale_size(
    range = c(0.5, 3.2), limits = c(0, P_CAP), breaks = c(1, 2, P_CAP),
    labels = c("1", "2", paste0(P_CAP, "+")), name = "-log10 p",
    guide = guide_legend(order = 2)
  ) +
  scale_discrete_manual(
    "stroke",
    values = c(`TRUE` = 0.7, `FALSE` = 0.15), breaks = "TRUE", labels = "FDR ≤ 0.05",
    name = NULL, guide = guide_legend(order = 3, override.aes = list(size = 2.4, fill = "white"))
  ) +
  labs(
    x = NULL, y = NULL, title = "Protein complexes and mitochondrial pathways",
    subtitle = paste(
      "Sets whose singscore differs at p ≤ 0.05 in a contrast; traits within sex, adjusted",
      "for SV1; measured members in brackets"
    )
  ) +
  theme_figure() +
  theme(
    axis.text.y = element_text(size = 5.2),
    axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 5.5),
    strip.placement = "outside", strip.text.y.left = element_text(angle = 90),
    strip.text = element_text(face = "bold", size = 5.5), panel.spacing = unit(1.2, "mm")
  )

save_panel(plot, "F05_Complexes", "a_set_matrix", 180, 120)
list(
  plot = plot,
  data = cells |>
    select(group, feature, label, members, section, column, r_equiv, P.Value, adj.P.Val) |>
    mutate(across(c(section, column), as.character)),
  inputs = inputs
)
