#!/usr/bin/env Rscript
# Figure 1: the design, then the three questions side by side (PAS vs vehicle in females, in males,
# and the treatment-by-sex interaction) before any single result. Reads saved outputs only and
# computes nothing beyond counting and scaling; the README lists each panel's source.

source(here::here("R", "helpers.R"))

out <- stage_paths("Figures", "F01_Overview")
inputs <- c(
  design = here("Figures", "F01_Overview", "assets", "design.png"),
  fit = here("00_Trunk", "07_Fit", "c_data", "fit.rds"),
  enrichment = here("00_Trunk", "08_Enrichment", "c_data", "enrichment.rds"),
  global = here("00_Trunk", "10_Global", "c_data", "global.rds"),
  fit_workbook = here("00_Trunk", "07_Fit", "c_data", "07_fit.xlsx")
)
fit <- readRDS(inputs[["fit"]])
enr <- readRDS(inputs[["enrichment"]])
global <- readRDS(inputs[["global"]])
SHORT <- c(PAS_vs_VEH_female = "Female", PAS_vs_VEH_male = "Male", treatment_by_sex = "F − M")
n_mice <- nrow(fit$targets)
per_group <- unique(count(fit$targets, group)$n)
stopifnot(length(per_group) == 1)
permanova_text <- \(tab) {
  with(filter(tab, term %in% c("sex", "treatment", "sex:treatment")), paste(
    "PERMANOVA",
    paste(sprintf("%s R² %.2f, p = %.2g", sub("sex:treatment", "interaction", term), r2, p),
      collapse = "\n"
    ),
    sep = "\n"
  ))
}

# PCA scores with a shaded 95% t-ellipse per group, a test summary in the top-left corner, and
# arrows (x, y, label) from the origin, already scaled to the score axes.
pca_panel <- function(scores, variance, title, subtitle, stats, arrows) {
  axis <- \(i) sprintf("PC%d (%.1f%%)", i, 100 * variance[i])
  scores <- mutate(scores, group = factor(group, FIGURE_GROUPS))
  ggplot(scores, aes(PC1, PC2)) +
    stat_ellipse(
      aes(colour = group, fill = group),
      geom = "polygon", type = "t", level = 0.95, linewidth = 0.3, alpha = 0.12
    ) +
    geom_point(aes(colour = group, shape = sex), size = 1.3) +
    annotate(
      "text",
      x = -Inf, y = Inf, label = stats, hjust = -0.05, vjust = 1.3, size = 1.9,
      lineheight = 0.95, colour = "grey20"
    ) +
    scale_colour_manual(values = GROUP_COLOURS, aesthetics = c("colour", "fill"), name = NULL) +
    scale_shape_manual(values = SEX_SHAPES, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.32))) +
    labs(x = axis(1), y = axis(2), title = title, subtitle = subtitle) +
    theme_figure() +
    geom_segment(
      aes(0, 0, xend = x, yend = y), arrows,
      arrow = grid::arrow(length = grid::unit(1, "mm")), linewidth = 0.3, colour = "grey30"
    ) +
    ggrepel::geom_text_repel(
      aes(x, y, label = label), arrows,
      size = 1.7, colour = "grey20", min.segment.length = Inf, box.padding = 0.15, seed = SEED
    )
}

# One sex's PAS signature (proteins at p <= 0.05, up and down) placed in the other sex's ranking,
# drawn as limma's barcode (Wu & Smyth 2012): members as ticks, the up list above the axis and
# the down list below, with a tricube moving average of membership along the rank. ranks holds
# the named moderated t of each target ranking; signatures comes from 08_Enrichment.
signature_barcode <- function(signatures, ranks, span = 0.12, min_height = 0.34) {
  lane_name <- \(source, ranked_in) {
    sprintf("%s signature\nin %s ranking", str_to_sentence(source), ranked_in)
  }
  lanes <- signatures |>
    mutate(
      lane = factor(
        lane_name(source, ranked_in), lane_name(c("female", "male"), c("male", "female"))
      ),
      side = if_else(direction == "up", 1, -1),
      ranks = map2(members, ranked_in, \(m, target) {
        position <- rank(-ranks[[target]], ties.method = "first")
        unname(position[str_split_1(m, ";")])
      })
    )
  n_total <- unique(lengths(ranks))
  stopifnot(length(n_total) == 1, all(lengths(lanes$ranks) == lanes$n_ranked))
  ticks <- select(lanes, lane, direction, side, rank = ranks) |>
    unnest(rank)
  worms <- lanes |>
    mutate(worm = map(ranks, \(r) {
      smooth <- limma::tricubeMovingAverage(as.integer(seq_len(n_total) %in% r), span = span)
      tibble(rank = seq_len(n_total), height = min_height + (1 - min_height) * smooth / max(smooth))
    })) |>
    select(lane, direction, side, worm) |>
    unnest(worm) |>
    mutate(y = side * height, y0 = side * min_height)
  labels <- lanes |>
    mutate(
      text = sprintf(
        "%s (n = %d), %.0f%% same way, p %.2g", direction, n_ranked, 100 * same_way, p
      ),
      x = if_else(direction == "up", n_total * 0.995, n_total * 0.005),
      hjust = if_else(direction == "up", 1, 0), y = side * 0.97, vjust = if_else(side > 0, 1, 0)
    )
  ggplot() +
    geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
    geom_segment(
      aes(rank, side * 0.06, xend = rank, yend = side * 0.3, colour = direction), ticks,
      linewidth = 0.2, alpha = 0.8
    ) +
    geom_ribbon(aes(rank, ymin = y0, ymax = y, fill = direction), worms, alpha = 0.25) +
    geom_line(aes(rank, y, colour = direction), worms, linewidth = 0.4) +
    geom_text(
      aes(x, y, label = text, hjust = hjust, vjust = vjust, colour = direction), labels,
      size = 1.9, fontface = "bold", show.legend = FALSE
    ) +
    facet_wrap(~lane, ncol = 1, strip.position = "left") +
    scale_colour_manual(
      values = c(up = DIRECTION_COLOURS[["Up"]], down = DIRECTION_COLOURS[["Down"]]),
      aesthetics = c("colour", "fill"), guide = "none"
    ) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(-1.05, 1.05), breaks = NULL, expand = expansion(0)) +
    theme_figure() +
    theme(
      strip.placement = "outside", strip.text.y.left = element_text(angle = 90, size = 5.5),
      panel.grid.major = element_blank(), panel.spacing.y = unit(2, "mm")
    )
}

# a: the design (assets/design.png).
panel_design <- ggplot() +
  annotation_custom(grid::rasterGrob(png::readPNG(inputs[["design"]]), interpolate = TRUE)) +
  labs(
    title = "Design",
    subtitle = str_glue(
      "{n_mice} mice, n = {per_group} per group; contrast colours carry through b to g"
    )
  ) +
  theme_void(base_size = 7) +
  theme(
    plot.title = element_text(size = 7.5, face = "bold"),
    plot.subtitle = element_text(size = 5.5, face = "italic", colour = "grey35")
  )

# b: proteins per contrast. The light bar counts p <= 0.05 or Π <= 0.05; the dark bar, drawn over
# it from zero, counts Π <= 0.05. Dashed lines: proteins at p <= 0.05 expected by chance per side.
protein_counts <- fit$protein_contrasts |>
  filter(contrast %in% ARC_CONTRASTS) |>
  mutate(direction = if_else(logFC > 0, "Up", "Down")) |>
  summarise(
    tests = n(), nominal = sum(P.Value <= ALPHA | pi_score <= ALPHA), pi = sum(pi_score <= ALPHA),
    fdr = sum(adj.P.Val <= ALPHA), .by = c(universe, contrast, direction)
  ) |>
  mutate(tests = sum(tests), .by = c(universe, contrast))
bar_fills <- c(
  "Down, p or Π" = "#9ECAE1", "Down, Π" = DIRECTION_COLOURS[["Down"]],
  "Up, p or Π" = "#FCBBA1", "Up, Π" = DIRECTION_COLOURS[["Up"]]
)
universe_size <- distinct(protein_counts, universe, tests) |> deframe()
universe_label <- c(
  full = sprintf("Full proteome (%s)", format(universe_size[["full"]], big.mark = ",")),
  mito = sprintf("Mitochondrial (%s)", format(universe_size[["mito"]], big.mark = ","))
)
fdr_text <- protein_counts |>
  summarise(fdr = sum(fdr), chance = ALPHA / 2 * tests[1], .by = c(universe, contrast)) |>
  mutate(
    band = map_chr(CONTRAST_COLOURS[CONTRAST_LABELS[contrast]], lighten, amount = 0.79),
    row = sprintf("%s\nFDR %d", SHORT[contrast], fdr), universe = universe_label[universe]
  )
bars <- protein_counts |>
  pivot_longer(c(nominal, pi), names_to = "tier", values_to = "n") |>
  mutate(
    signed = if_else(direction == "Up", n, -n),
    fill = paste0(direction, if_else(tier == "pi", ", Π", ", p or Π")),
    universe = universe_label[universe]
  ) |>
  left_join(select(fdr_text, universe, contrast, row, band), by = c("universe", "contrast")) |>
  mutate(
    row = factor(row, unique(fdr_text$row[order(match(fdr_text$contrast, rev(ARC_CONTRASTS)))]))
  ) |>
  arrange(tier == "pi")
panel_proteins <- ggplot(bars, aes(signed, row)) +
  geom_tile(
    aes(x = 0, fill = I(band)), distinct(bars, row, universe, band),
    width = Inf, height = 0.95
  ) +
  geom_col(aes(fill = fill), width = 0.6, position = "identity") +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  geom_vline(
    aes(xintercept = chance), distinct(fdr_text, universe, chance),
    linetype = "dashed", linewidth = 0.25
  ) +
  geom_vline(
    aes(xintercept = -chance), distinct(fdr_text, universe, chance),
    linetype = "dashed", linewidth = 0.25
  ) +
  facet_wrap(~universe, ncol = 1, scales = "free") +
  scale_fill_manual(values = bar_fills, breaks = names(bar_fills), name = NULL) +
  scale_x_continuous(labels = abs, expand = expansion(mult = 0.08)) +
  labs(
    x = "proteins (down | up)", y = NULL, title = "Proteins per contrast",
    subtitle = "Light: p ≤ 0.05 or Π ≤ 0.05; dark: Π ≤ 0.05; dashed: chance at p ≤ 0.05"
  ) +
  theme_figure() +
  theme(
    legend.position = "bottom", legend.key.size = unit(2.5, "mm"),
    panel.grid.major.y = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 2))

# c: proteome PCA of all twenty mice with the MitoCarta share fitted as a vector (10_Global).
scores <- global$proteome_scores
reach <- 0.8 * max(abs(c(scores$PC1, scores$PC2)))
mito_arrow <- with(global$mito_vector, tibble(
  x = PC1 * reach, y = PC2 * reach, label = sprintf("MitoCarta share, R² %.2f", r2)
))
panel_proteome <- pca_panel(
  scores, global$proteome_variance, "Proteome",
  sprintf("%d mice, %s proteins", n_mice, format(universe_size[["full"]], big.mark = ",")),
  permanova_text(global$permanova), mito_arrow
)

# d: phenotype PCA of the six traits, scaled, with their loadings.
phenotype_scores <- global$phenotype_scores
reach <- 0.8 * max(abs(c(phenotype_scores$PC1, phenotype_scores$PC2)))
loadings <- global$phenotype_loadings |>
  transmute(
    x = PC1 / max(abs(c(PC1, PC2))) * reach, y = PC2 / max(abs(c(PC1, PC2))) * reach,
    label = PHENOTYPE_LABELS[phenotype]
  )
panel_phenotype <- pca_panel(
  phenotype_scores, global$phenotype_variance, "Phenotype", "Six traits, scaled",
  permanova_text(global$phenotype_permanova), loadings
) +
  theme(legend.position = "none")

# e: phenotype effects of PAS within each sex, with the interaction p beside each row.
interaction_p <- with(fit$phenotype_interaction, set_names(P.Value, PHENOTYPE_LABELS[phenotype]))
panel_effects <- phenotype_effect_plot(fit$phenotype_effects,
  title = "Phenotype effects",
  subtitle = "PAS vs vehicle within sex, Welch t on log values; F − M: interaction p"
) +
  scale_y_discrete(labels = \(x) sprintf("%s\nF − M p %.2g", x, interaction_p[x])) +
  theme(legend.position = "bottom", legend.key.size = unit(2.5, "mm"))

# f: every set at fgsea FDR <= 0.05 in any contrast that the EnrichmentMap flag keeps there, listed
# once under the first contrast that selects it and drawn in all three. † marks rows whose
# singscore contrast, tested across mice, has p > 0.05 in the selecting contrast.
selected <- enr$results |>
  filter(padj <= ALPHA, dedup_status == "kept") |>
  mutate(selected_by = factor(contrast, ARC_CONTRASTS)) |>
  arrange(selected_by, p) |>
  distinct(term, .keep_all = TRUE) |>
  transmute(
    term, collection, selected_by,
    order = row_number(),
    label = paste0(
      enrichVolcano::clean_label(term, width = 48), if_else(score_p > ALPHA, " †", "")
    )
  )
pathways <- enr$results |>
  inner_join(select(selected, term, selected_by, order, label), by = "term") |>
  mutate(
    column = factor(SHORT[contrast], SHORT),
    facet = factor(SHORT[as.character(selected_by)], SHORT),
    log_fdr = pmin(-log10(padj), 4), row = reorder(term, -order)
  )
nes_limit <- ceiling(max(abs(pathways$NES)) * 10) / 10
databases <- distinct(pathways, row, facet, collection) |>
  mutate(column = factor(" ", c(" ", SHORT)))
panel_pathways <- ggplot(pathways, aes(column, row)) +
  geom_point(aes(size = log_fdr, fill = NES), shape = 21, stroke = 0.2) +
  geom_point(
    aes(x = 0.4, colour = collection), databases,
    shape = 15, size = 1.8
  ) +
  facet_grid(facet ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(labels = deframe(distinct(pathways, row, label))) +
  scale_x_discrete(expand = expansion(add = c(0.9, 0.5))) +
  scale_colour_manual(values = DATABASE_COLOURS, name = "Collection") +
  scale_fill_gradient2(
    low = DIRECTION_COLOURS[["Down"]], high = DIRECTION_COLOURS[["Up"]],
    limits = c(-nes_limit, nes_limit), name = "NES"
  ) +
  scale_size_area(
    max_size = 3, limits = c(0, 4), breaks = c(1, 2, 4), labels = c("1", "2", "4+"),
    name = "-log10 FDR"
  ) +
  guides(
    fill = guide_colourbar(order = 1), size = guide_legend(order = 2),
    colour = guide_legend(order = 3)
  ) +
  labs(
    x = NULL, y = NULL, title = "Pathways at FDR ≤ 0.05",
    subtitle = paste(
      "fgsea, GO:BP and Reactome (full proteome), MitoCarta (mitochondrial);",
      "grouped by the contrast that selects them; † singscore p > 0.05"
    )
  ) +
  theme_figure() +
  theme(
    axis.text.y = element_text(size = 5.2), strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, size = 5.5, face = "bold"),
    panel.spacing = unit(1, "mm"), legend.key.size = unit(2.5, "mm")
  )

# g: each sex's p <= 0.05 signature in the other sex's ranking (08_Enrichment, signatures).
target_ranks <- list(
  male = enr$ranks$PAS_vs_VEH_male$full, female = enr$ranks$PAS_vs_VEH_female$full
)
panel_barcode <- signature_barcode(enr$signatures, target_ranks) +
  labs(
    x = sprintf(
      "rank in the other sex's ranking (n = %s)",
      format(length(target_ranks$male), big.mark = ",")
    ),
    y = NULL, title = "Cross-sex signatures",
    subtitle = "Each sex's p ≤ 0.05 proteins in the other's ranking; fgsea p"
  )

# The plotted counts agree with the tables 07_Fit and 08_Enrichment wrote before anything is saved.
fit_summary <- readxl::read_excel(inputs[["fit_workbook"]], "contrast_summary") |>
  filter(level == "protein", contrast %in% ARC_CONTRASTS)
plotted <- summarise(protein_counts, pi = sum(pi), fdr = sum(fdr), .by = c(universe, contrast)) |>
  inner_join(fit_summary, by = c("universe", "contrast"))
pathway_fdr <- count(filter(pathways, padj <= ALPHA), contrast, name = "plotted") |>
  inner_join(summarise(enr$summary, tested = sum(fdr_le_0.05), .by = contrast), by = "contrast")
stopifnot(
  nrow(plotted) == 6, plotted$pi == plotted$pi_le_0.05, plotted$fdr == plotted$fdr_le_0.05,
  nrow(scores) == n_mice, nrow(phenotype_scores) == n_mice,
  all(count(pathways, term)$n == 3), pathway_fdr$plotted <= pathway_fdr$tested,
  nrow(enr$signatures) == 4
)

figure <- (panel_design | panel_proteins | panel_proteome | panel_phenotype) /
  (panel_effects | panel_pathways | panel_barcode) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))
figure[[2]] <- figure[[2]] + plot_layout(widths = c(1, 2, 1))
save_figure(figure, out$reports, "f01_overview", width_mm = 300, height_mm = 190)

write_workbook(
  list(
    protein_counts = protein_counts,
    pca_proteome = scores,
    pca_phenotype = phenotype_scores,
    permanova = bind_rows(
      mutate(global$permanova, data = "proteome", .before = 1),
      mutate(global$phenotype_permanova, data = "phenotype", .before = 1)
    ),
    phenotype_effects = fit$phenotype_effects,
    phenotype_interaction = fit$phenotype_interaction,
    pathways = select(
      pathways, term, label, collection, selected_by, contrast, NES, p, padj, score_p, r_group
    ),
    signatures = select(enr$signatures, -members)
  ),
  c(
    paste(
      "Panel b: per universe, contrast and direction, proteins at p or Π ≤ 0.05,",
      "at Π ≤ 0.05, at FDR ≤ 0.05."
    ),
    "Panel c: proteome PC1 to PC4 per mouse, with technical variables.",
    "Panel d: phenotype PC1 and PC2 per mouse.",
    "Panels c, d: PERMANOVA for the proteome and the phenotypes.",
    "Panel e: PAS vs vehicle per phenotype and sex, % change, 95% CI, Welch p.",
    "Panel e: treatment by sex per phenotype, female response relative to male.",
    "Panel f: each selected set in all three contrasts: NES, p, FDR, singscore p, r with group.",
    "Panel g: each sex's p ≤ 0.05 up and down lists in the other sex's ranking."
  ),
  out$workbook, inputs
)

sessionInfo()
