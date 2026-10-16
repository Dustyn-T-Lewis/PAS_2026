# S6 a: the proteins PAS changes in either sex (moderated F over the two PAS contrasts, p <= 0.05)
# across all twenty mice, each centred on its own sex's vehicle mean, split by response pattern
# and clustered within each split. Beside each protein, its r with each phenotype within each sex.
# Rows were chosen for differing, so the clustering describes them and tests nothing.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(fit = output_path("02_Statistics", "01_Fit", "fit.rds"))
fit <- readRDS(inputs[["fit"]])
PATTERNS <- c("up in both", "down in both", "female up, male down", "female down, male up")

selected <- fit$protein_ftest |>
  filter(universe == "full", P.Value <= ALPHA) |>
  mutate(pattern = factor(pattern, PATTERNS)) |>
  arrange(pattern)
targets <- arrange(fit$targets, factor(group, FIGURE_GROUPS))
values <- fit$E_norm$full[selected$feature, targets$sample_id]
vehicle_mean <- sapply(c(F = "VEH_F", M = "VEH_M"), \(g) {
  rowMeans(values[, targets$group == g])
})
centred <- values - vehicle_mean[, targets$sex]

r_by_trait <- fit$protein_associations |>
  filter(model == "pooled", universe == "full", feature %in% selected$feature) |>
  transmute(
    feature,
    column = paste(TRAIT_SHORT[phenotype], if_else(sex == "female", "F", "M")), r
  ) |>
  pivot_wider(names_from = column, values_from = r)
r_matrix <- as.matrix(column_to_rownames(r_by_trait, "feature"))[
  selected$feature, paste(rep(TRAIT_SHORT, each = 2), c("F", "M"))
]
r_colours <- circlize::colorRamp2(c(-1, 0, 1), c("#8C510A", "white", "#01665E"))
right <- ComplexHeatmap::rowAnnotation(
  r = r_matrix, col = list(r = r_colours), na_col = "grey90",
  annotation_name_gp = grid::gpar(fontsize = 5), annotation_name_rot = 90,
  simple_anno_size = unit(1.8, "mm"), gap = unit(1, "mm"),
  annotation_legend_param = list(
    r = list(
      title = "r with\ntrait", title_gp = grid::gpar(fontsize = 5.5, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 5.5), grid_width = unit(2, "mm")
    )
  )
)
counts <- count(selected, pattern)
split <- factor(
  sprintf("%s\n(%d)", selected$pattern, counts$n[match(selected$pattern, counts$pattern)]),
  sprintf("%s\n(%d)", counts$pattern, counts$n)
)

heatmap <- protein_heatmap(centred, targets$group,
  row_split = split, right = right, limit = 1.5, width_mm = 116, height_mm = 132
)
plot <- wrap_elements(plot = heatmap) +
  labs(
    title = "Proteins PAS changes in either sex",
    subtitle = sprintf(
      "Moderated F over both PAS contrasts, p ≤ 0.05 (%d proteins); centred on own-sex vehicle",
      nrow(selected)
    )
  ) +
  theme_figure()

save_panel(plot, "S6_Protein_patterns", "a_ftest_heatmap", 116, 145)
list(
  plot = plot,
  data = selected |>
    select(feature, gene, pattern, F_stat, P.Value, adj.P.Val, starts_with("logFC")) |>
    left_join(r_by_trait, by = "feature"),
  inputs = inputs
)
