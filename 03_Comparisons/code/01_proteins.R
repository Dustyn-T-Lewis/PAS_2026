# Step 01 of every comparison: which proteins its contrast changed.
# Run by 03_Comparisons/<comparison>/a_script/run_comparison.R.

source(here::here("00_Setup", "Shared_functions", "comparisons.R"))

# Proteins: the comparison's contrast sliced from the shared fit, its tiers and Π lists, and the
# STRING connectivity of those lists from 02_Enrichment (descriptive).
report_proteins <- function(comparison) {
  out <- comparison_paths(comparison, "01_proteins")
  inputs <- c(
    fit = output_path("02_Statistics", "01_Fit", "fit.rds"),
    enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds")
  )
  fit <- readRDS(inputs[["fit"]])
  contrast <- COMPARISONS[[comparison]]$contrast
  label <- COMPARISONS[[comparison]]$label
  targets <- filter(fit$targets, sample_id %in% comparison_mice(comparison, fit$targets))

  results <- fit$protein_contrasts |>
    filter(contrast == !!contrast) |>
    mutate(
      tier = tier_of(P.Value, adj.P.Val, pi_score), direction = if_else(logFC > 0, "up", "down")
    ) |>
    arrange(universe, P.Value)
  tier_table <- results |>
    reframe(tier_counts(P.Value, adj.P.Val, pi_score), .by = c(universe, direction))
  pi_lists <- results |>
    filter(universe == "full", pi_score <= ALPHA, !is.na(gene), gene != "") |>
    select(direction, feature, gene, logFC, P.Value, adj.P.Val, pi_score, r_group) |>
    arrange(direction, pi_score)
  string <- filter(readRDS(inputs[["enrichment"]])$string, contrast == !!contrast)

  volcano <- results |>
    mutate(universe = if_else(universe == "full", "Full proteome", "Mitochondrial universe")) |>
    ggplot(aes(logFC, -log10(P.Value), colour = tier)) +
    geom_point(size = 0.6) +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed", linewidth = 0.3) +
    ggrepel::geom_text_repel(
      data = \(d) filter(d, adj.P.Val <= ALPHA | rank(P.Value) <= 10, .by = universe),
      aes(label = gene), size = 2, colour = "black", max.overlaps = Inf,
      min.segment.length = 0, seed = SEED
    ) +
    facet_wrap(~universe) +
    scale_colour_manual(values = tier_colours(comparison), name = NULL, drop = FALSE) +
    labs(
      x = "log2 fold change", y = "-log10 p",
      title = str_glue("{label}: {CONTRAST_LABELS[contrast]}, every protein"),
      subtitle = "limpa dpcDE, BH within universe; labelled: FDR ≤ 0.05 and the 10 smallest p"
    ) +
    theme_figure() +
    theme(legend.position = "bottom")
  tiers <- c("p_le_0.05", "pi_le_0.05", "fdr_le_0.05")
  counts <- tier_table |>
    pivot_longer(all_of(tiers), names_to = "tier", values_to = "n") |>
    mutate(
      tier = factor(tier, tiers, c("p ≤ 0.05", "Π ≤ 0.05", "FDR ≤ 0.05")),
      universe = if_else(universe == "full", "Full proteome", "Mitochondrial universe"),
      n = if_else(direction == "down", -n, n)
    ) |>
    ggplot(aes(n, tier, fill = direction)) +
    geom_col(width = 0.7) +
    geom_vline(
      aes(xintercept = if_else(direction == "down", -1, 1) * expected_by_chance),
      data = \(d) filter(d, tier == "p ≤ 0.05"), linetype = "dashed", linewidth = 0.3
    ) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_text(aes(label = abs(n), hjust = if_else(n < 0, 1.2, -0.2)), size = 2.2) +
    facet_wrap(~universe, scales = "free_x") +
    scale_fill_manual(
      values = c(down = DIRECTION_COLOURS[["Down"]], up = DIRECTION_COLOURS[["Up"]]), name = NULL
    ) +
    scale_x_continuous(labels = abs, expand = expansion(mult = 0.2)) +
    labs(
      x = "Proteins (down | up)", y = NULL, title = "Proteins per tier",
      subtitle = paste(
        "Dashed: proteins at p ≤ 0.05 expected by chance in each direction;",
        "Π can include p > 0.05"
      )
    ) +
    theme_figure()
  hits <- filter(results, universe == "full", P.Value <= ALPHA)
  heatmap <- group_heatmap(fit$E_norm$full[hits$feature, targets$sample_id, drop = FALSE], targets,
    title = str_glue("All {nrow(hits)} proteins at p ≤ 0.05, per mouse"),
    subtitle = "Rows z-scored across these mice, clustered on 1 - Pearson r (display only)"
  )
  is_interaction <- is.na(COMPARISONS[[comparison]]$sex)
  strips <- with(hits, sprintf(
    "%s\nlog2FC %.2f%s\np %.2g, FDR %.2g, Π %.2g", gene, logFC,
    if (is_interaction) "" else sprintf(", r %.2f", r_group), P.Value, adj.P.Val, pi_score
  ))
  grids <- mouse_grid(long_values(fit$E_norm$full, hits$feature, strips, targets),
    title = str_glue("{label}: every protein at p ≤ 0.05, in order of p"),
    subtitle = if (is_interaction) {
      paste(
        "log2 abundance per mouse, bar at the group median;",
        "log2FC: female PAS response minus male PAS response"
      )
    } else {
      "log2 abundance per mouse, bar at the group median; r: Pearson r with PAS (PAS = 1)"
    }
  )
  plot_index <- write_report(file.path(out$reports, "01_proteins.pdf"), c(
    list(
      "volcanoes, every protein" = volcano, "proteins per tier" = counts,
      "heatmap of proteins at p ≤ 0.05" = heatmap
    ),
    set_names(grids, sprintf("values by mouse, page %d", seq_along(grids)))
  ))

  saveRDS(
    list(results = results, summary = tier_table, pi_lists = pi_lists, string = string),
    file.path(out$data, "01_proteins.rds")
  )
  write_workbook(
    list(
      summary = tier_table, results = select(results, -genes), pi_lists = pi_lists, string = string,
      plot_index = plot_index
    ),
    c(
      "Per universe and direction: p ≤ 0.05, chance, Π ≤ 0.05, FDR tiers.",
      paste(
        "Every protein in both universes, by p:",
        "log2 FC with 95% CI, t, p, BH FDR, Π, r with PAS, tier."
      ),
      "Full-universe proteins at Π ≤ 0.05, by direction.",
      "STRING v12 connectivity of each list against the detected proteome (descriptive).",
      "Which PDF page shows what."
    ),
    out$workbook, inputs
  )
  tier_table
}
