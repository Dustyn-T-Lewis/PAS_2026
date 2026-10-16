# Step 01 of every comparison: which proteins its contrast changed.
# Run by 03_Comparisons/<comparison>/a_script/run_comparison.R.

source(here::here("00_Setup", "Shared_functions", "comparisons.R"))

# The strongest tier a protein reaches: FDR, then Π, then nominal p.
tier_of <- function(p, fdr, pi) {
  factor(
    case_when(
      fdr <= ALPHA ~ "FDR ≤ 0.05", pi <= ALPHA ~ "Π ≤ 0.05",
      p <= ALPHA ~ "p ≤ 0.05",
      .default = "p > 0.05"
    ),
    c("p > 0.05", "p ≤ 0.05", "Π ≤ 0.05", "FDR ≤ 0.05")
  )
}

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

  volcano <- \(u) {
    \() {
      d <- filter(results, universe == u)
      plot(d$logFC, -log10(d$P.Value),
        pch = 16, cex = 0.4, col = if_else(d$P.Value <= ALPHA, "grey20", "grey75"),
        xlab = "log2 fold change", ylab = "-log10 p",
        main = str_glue("{label}, {CONTRAST_LABELS[contrast]}, {u}")
      )
      abline(h = -log10(ALPHA), lty = 2)
      top <- slice_min(d, P.Value, n = 10)
      text(top$logFC, -log10(top$P.Value), top$gene, cex = 0.5, pos = 3)
    }
  }
  hits <- filter(results, universe == "full", P.Value <= ALPHA)
  plot_index <- preview_pdf(file.path(out$reports, "01_proteins.pdf"), list(
    "volcano, full proteome" = volcano("full"),
    "volcano, MitoCarta" = volcano("mito"),
    "proteins at p ≤ 0.05, by mouse" = \() {
      heatmap(fit$E_norm$full[hits$feature, targets$sample_id, drop = FALSE],
        scale = "row", labRow = hits$gene, cexRow = 0.3, cexCol = 0.6,
        col = hcl.colors(50, "Blue-Red 3"), ColSideColors = unname(GROUP_COLOURS[targets$group])
      )
    }
  ))

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
