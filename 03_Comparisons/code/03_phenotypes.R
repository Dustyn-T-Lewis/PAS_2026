# Step 03 of the female and male comparisons: phenotype effects and phenotype associations.
# Run by 03_Comparisons/<comparison>/a_script/run_comparison.R.

source(here::here("00_Setup", "Shared_functions", "comparisons.R"))

# The protein with the smallest association p whose slope keeps its sign and p <= 0.05 when any
# one mouse is left out (ordinary least squares within the sex), so an example is not one mouse's
# doing.
robust_example <- function(associations, E, x, mice, max_tries = 50) {
  candidates <- slice_min(associations, P.Value, n = max_tries, with_ties = FALSE)
  for (i in seq_len(nrow(candidates))) {
    y <- E[candidates$feature[i], mice]
    holds <- map_lgl(seq_along(mice), \(k) {
      f <- summary(lm(y[-k] ~ x[-k]))$coefficients
      sign(f[2, 1]) == sign(candidates$estimate[i]) && f[2, 4] <= ALPHA
    })
    if (all(holds)) {
      return(candidates[i, ])
    }
  }
  candidates[0, ]
}

# Associations, two different questions kept apart. Phenotype effects: what PAS did to each
# phenotype in this sex (% change against vehicle, a treatment effect). Phenotype associations:
# which proteins, pathway scores, modules, complexes and clusters move with each phenotype across
# this sex's ten PAS and vehicle mice (01_Fit; r pooled and within each group beside every test),
# under two models: pooled across PAS and vehicle, and within treatment (treatment in the model).
# An association is co-movement, not a PAS effect. Pathway patterns also come from fgsea on the
# protein association t (02_Enrichment); that test permutes genes, so its FDR runs optimistic
# against the mouse-level singscore association of the same set.
report_phenotypes <- function(comparison) {
  out <- comparison_paths(comparison, "03_phenotypes")
  inputs <- c(
    fit = output_path("02_Statistics", "01_Fit", "fit.rds"),
    enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds")
  )
  fit <- readRDS(inputs[["fit"]])
  enrichment <- readRDS(inputs[["enrichment"]])
  sex_code <- COMPARISONS[[comparison]]$sex
  label <- COMPARISONS[[comparison]]$label
  targets <- filter(fit$targets, sample_id %in% comparison_mice(comparison, fit$targets))
  mice <- targets$sample_id
  phenotype_of <- \(p) fit$phenotypes[[p]][match(mice, fit$targets$sample_id)]

  effects <- filter(fit$phenotype_effects, sex == sex_code)
  associations <- fit$protein_associations |>
    mutate(level = "protein") |>
    bind_rows(fit$score_associations) |>
    filter(sex == comparison) |>
    mutate(label = coalesce(gene, feature), .after = feature)
  association_summary <- associations |>
    reframe(tier_counts(P.Value, adj.P.Val), .by = c(level, universe, phenotype, model)) |>
    select(-pi_le_0.05)
  associations <- filter(associations, model != "pooled_mito") |>
    arrange(level, phenotype, model, P.Value)
  # One row per feature and phenotype: the pooled association with the within-treatment p beside it.
  pooled <- filter(associations, model == "pooled") |>
    left_join(
      associations |>
        filter(model == "treatment") |>
        select(level, universe, feature, phenotype,
          p_treatment = P.Value, fdr_treatment = adj.P.Val
        ),
      by = c("level", "universe", "feature", "phenotype"), relationship = "one-to-one",
      unmatched = "error"
    )
  r_group <- bind_rows(
    mutate(fit$protein_contrasts, level = "protein"), fit$score_contrasts
  ) |>
    filter(contrast == COMPARISONS[[comparison]]$contrast) |>
    select(level, universe, feature, r_group, group_p = P.Value, group_fdr = adj.P.Val)
  fgsea <- enrichment$association_fgsea |>
    filter(sex == comparison) |>
    select(-sex)

  examples <- map(PHENOTYPES, \(p) {
    pooled |>
      filter(level == "protein", universe == "full", phenotype == !!p) |>
      robust_example(fit$E_norm$full, phenotype_of(p), mice)
  }) |>
    list_rbind()

  protein_counts <- association_summary |>
    filter(level == "protein", universe == "full", model != "pooled_mito")
  plot_index <- preview_pdf(file.path(out$reports, "03_phenotypes.pdf"), list(
    "phenotype effects" = \() {
      par(mar = c(4, 10, 2, 1))
      y <- seq_len(nrow(effects))
      plot(effects$pct_change, y,
        xlim = range(effects$ci_low, effects$ci_high), pch = 16, yaxt = "n", ylab = "",
        xlab = "% change, PAS vs vehicle (95% CI)", main = str_glue("{label}: phenotype effects")
      )
      segments(effects$ci_low, y, effects$ci_high, y)
      abline(v = 0, lty = 2)
      axis(2, y, PHENOTYPE_LABELS[effects$phenotype], las = 1, cex.axis = 0.7)
    },
    "protein associations at p ≤ 0.05, by phenotype and model" = \() {
      par(mar = c(8, 4, 2, 1))
      counts <- xtabs(p_le_0.05 ~ model + phenotype, protein_counts)
      colnames(counts) <- PHENOTYPE_LABELS[colnames(counts)]
      barplot(counts,
        beside = TRUE, legend.text = TRUE, las = 2, cex.names = 0.7,
        ylab = "proteins at p ≤ 0.05", main = str_glue("{label}: protein associations")
      )
    },
    "strongest robust protein per phenotype" = \() {
      par(mfrow = c(2, 3))
      pwalk(examples, \(feature, phenotype, label, ...) {
        x <- phenotype_of(phenotype)
        y <- fit$E_norm$full[feature, mice]
        plot(x, y,
          pch = 16, col = GROUP_COLOURS[targets$group], xlab = PHENOTYPE_LABELS[[phenotype]],
          ylab = "log2 abundance", main = label
        )
        abline(lm(y ~ x), lty = 2)
      })
    }
  ))

  write_workbook(
    list(
      phenotype_effects = effects, association_summary = association_summary,
      protein_associations = filter(associations, level == "protein"),
      score_associations = filter(associations, level != "protein"),
      pathway_association_fgsea = fgsea, treatment_r = r_group, examples = examples,
      plot_index = plot_index
    ),
    c(
      "Phenotype effects: PAS vs vehicle per phenotype in this sex, % change, 95% CI, Welch p.",
      paste(
        "Per level, universe, phenotype and model: associations at p ≤ 0.05, FDR tiers;",
        "pooled_mito adds the pellet mito share to the pooled model."
      ),
      paste(
        "Every protein's association with each phenotype across this sex's ten mice, pooled and",
        "within treatment: slope (log2 per SD), 95% CI, t, p, BH FDR, r pooled, in VEH, in PAS."
      ),
      paste(
        "Every pathway, module, complex and cluster score's association with each phenotype,",
        "same columns."
      ),
      paste(
        "fgsea on protein association t per phenotype and model (gene permutation),",
        "beside the same set's singscore association (mouse level)."
      ),
      paste(
        "Treatment effect beside the associations: every feature's r with PAS (PAS = 1)",
        "in this sex, with the contrast's p and FDR."
      ),
      "The robust example protein per phenotype.",
      "Which PDF page shows what."
    ),
    out$workbook
  )
  association_summary
}
