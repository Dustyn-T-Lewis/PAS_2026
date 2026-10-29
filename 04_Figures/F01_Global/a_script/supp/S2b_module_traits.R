# S2 b: each WGCNA module eigengene against the three contrasts and the six phenotypes, after YvO's
# module-trait heatmap. Every cell holds the r-equivalent of its limma t, t / sqrt(t^2 + df), so a
# cell's value and its outline come from the same test: the contrast, or the eigengene's slope on
# the trait within one sex (pooled model, adjusted for the surrogate variable). Solid outline: BH
# FDR <= 0.05 over the modules; dashed: p <= 0.05 only; bold: |r| >= BOLD_R. Left, each module's
# size, MitoCarta share and stability (median best-match Jaccard over subsamples of mice, from
# 06_Modules). Right, each eigengene's AUC for PAS against vehicle within a sex (03_Classification),
# filled on 2 * AUC - 1 so chance is white and the colours keep their direction; outlined on the
# exact Wilcoxon p and its BH FDR over the modules. Pearson r across the mice of the sex is in the
# workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(
  fit = output_path("02_Statistics", "01_Fit", "fit.rds"),
  modules = output_path("01_Processing", "06_Modules", "modules.rds"),
  classification = output_path("02_Statistics", "03_Classification", "classification.rds")
)
fit <- readRDS(inputs[["fit"]])
AUC_SECTION <- "AUC (PAS vs vehicle)"
BOLD_R <- 0.4
sizes <- readRDS(inputs[["modules"]])$sizes |>
  filter(module != "grey") |>
  arrange(proteins) |>
  mutate(
    module_j = sprintf("%s (J %.2f)", module, median_jaccard),
    name = if_else(is.na(label) | label == "unannotated", module_j,
      sprintf("%s · %s", module_j, set_label(label, 200))
    ),
    row = factor(module, module)
  )

modules_only <- \(x) filter(x, level == "module")
evidence_levels <- c(fdr = paste("FDR ≤", ALPHA), p = paste("p ≤", ALPHA))
auc_cells <- readRDS(inputs[["classification"]])$roc |>
  modules_only() |>
  transmute(
    module = feature, section = AUC_SECTION, column = str_to_sentence(comparison), r = 2 * auc - 1,
    P.Value = p, adj.P.Val = fdr, auc, auc_low, auc_high
  )
score_sections <- score_cells(
  modules_only(fit$score_contrasts), modules_only(fit$score_associations)
)
cells <- score_sections |>
  select(module = feature, section, column, r = r_equiv, P.Value, adj.P.Val, pearson_r) |>
  mutate(across(c(section, column), as.character)) |>
  bind_rows(auc_cells) |>
  mutate(
    section = factor(section, c(levels(score_sections$section), AUC_SECTION)),
    column = factor(column, unique(c(levels(score_sections$column), "Female", "Male"))),
    row = factor(module, levels(sizes$row)),
    text = sprintf("%.2f", coalesce(auc, r)),
    evidence = case_when(
      adj.P.Val <= ALPHA ~ evidence_levels[["fdr"]], P.Value <= ALPHA ~ evidence_levels[["p"]],
      .default = "none"
    ),
    face = if_else(is.na(auc) & abs(r) >= BOLD_R, "bold", "plain")
  )

# The AUC block is two columns wide, so its title takes two lines.
section_labels <- set_names(levels(cells$section)) |>
  replace(AUC_SECTION, "AUC (PAS\nvs vehicle)")
heatmap <- ggplot(cells, aes(column, row)) +
  geom_tile(aes(fill = r), colour = "grey85", linewidth = 0.2) +
  geom_tile(
    aes(linetype = evidence), filter(cells, evidence != "none"),
    fill = NA, colour = "black", linewidth = 0.45, width = 0.92, height = 0.88
  ) +
  geom_text(aes(label = text, fontface = face), size = 1.55) +
  facet_grid(
    ~section,
    scales = "free_x", space = "free_x", labeller = as_labeller(section_labels)
  ) +
  scale_fill_distiller(
    palette = "PiYG", direction = 1, limits = c(-1, 1), breaks = c(-1, 0, 1),
    name = "r-equivalent\nor 2·AUC − 1",
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
  geom_col(aes(fill = module), colour = "grey30", linewidth = 0.2, width = 0.75) +
  scale_fill_manual(values = MODULE_COLOURS, guide = "none") +
  geom_text(aes(label = sprintf("%.0f%%", 100 * mito_share)), hjust = 1.15, size = 5.3 / .pt) +
  scale_x_reverse(breaks = size_breaks, expand = expansion(mult = c(0.3, 0))) +
  scale_y_discrete(
    labels = deframe(transmute(sizes, row, name = str_wrap(name, 46))), position = "right"
  ) +
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
      paste(
        "r-equivalent of each t, bold if absolute r ≥ %s; AUC cells show the AUC,",
        "filled by 2·AUC − 1; J: subsample stability"
      ),
      BOLD_R
    )
  ) +
  theme_figure()

save_figure(plot, "F01_Global", "S2b_module_traits", 180, 92)
list(
  plot = plot,
  data = left_join(select(cells, -row, -text), select(sizes, -row), by = "module")
)
