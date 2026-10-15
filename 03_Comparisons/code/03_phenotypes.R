# Step 03 of the female and male comparisons: phenotype effects and phenotype associations.
# Run by 03_Comparisons/<comparison>/a_script/run_comparison.R.

source(here::here("00_Setup", "Shared_functions", "comparisons.R"))

# Scatter of feature against phenotype per facet: a line within vehicle, a line within PAS and a
# dashed line through all the sex's mice, the association that is tested.
association_grid <- function(long, title, subtitle, ncol = 5, nrow = 4) {
  paginate(long, per_page = ncol * nrow, draw = \(d) {
    ggplot(d, aes(phenotype_value, value)) +
      geom_smooth(
        aes(group = 1),
        method = "lm", formula = y ~ x, se = FALSE, colour = "grey20", linetype = "dashed",
        linewidth = 0.4
      ) +
      geom_smooth(
        aes(colour = group),
        method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.4
      ) +
      geom_point(aes(colour = group), size = 1.1) +
      facet_wrap(~strip, scales = "free", ncol = ncol, nrow = nrow, drop = FALSE) +
      scale_colour_manual(values = GROUP_COLOURS, name = NULL) +
      labs(x = "Phenotype", y = "Feature value", title = title, subtitle = subtitle) +
      theme_figure() +
      theme(
        strip.text = element_text(size = 4.6, face = "plain"), axis.text = element_text(size = 4.5),
        legend.position = "bottom"
      )
  })
}

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
# which proteins, pathway scores (GO:BP, Reactome, MitoCarta), modules, complexes and clusters move
# with each phenotype across this sex's ten PAS and vehicle mice (01_Fit; r pooled and within each
# group beside every test), under two models: pooled across PAS and vehicle, and within treatment
# (treatment in the model). An association is co-movement, not a PAS effect. Pathway patterns also
# come from fgsea on the protein association t (02_Enrichment); that test permutes genes, so its FDR
# runs optimistic against the mouse-level singscore association of the same set.
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

  effects_page <- phenotype_effect_plot(effects,
    title = str_glue("{label}: phenotype effects of PAS"),
    subtitle = "% change against vehicle, Welch t on log values, n = 5 v 5"
  )
  grid_page <- association_summary |>
    filter(model != "pooled_mito") |>
    mutate(
      model = factor(model, c("pooled", "treatment"), c("Pooled", "Within treatment")),
      level = factor(
        paste0(LEVEL_LABELS[level], if_else(universe == "mito" & level == "protein", ", mito", "")),
        c("Proteins", "Proteins, mito", LEVEL_LABELS[-1])
      ),
      trait = factor(PHENOTYPE_LABELS[phenotype], rev(PHENOTYPE_LABELS))
    ) |>
    summarise(
      across(c(tests, p_le_0.05, expected_by_chance, fdr_le_0.05), sum),
      .by = c(model, level, trait)
    ) |>
    mutate(ratio = p_le_0.05 / expected_by_chance) |>
    ggplot(aes(level, trait, fill = ratio)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = str_glue("{p_le_0.05} / {tests}\nFDR {fdr_le_0.05}")), size = 1.9) +
    facet_wrap(~model) +
    scale_fill_gradient(low = "white", high = "grey35", name = "p ≤ 0.05 /\nchance") +
    labs(
      x = NULL, y = NULL, title = str_glue("{label}: phenotype associations by level"),
      subtitle = paste(
        "This sex's ten mice, pooled across PAS and vehicle or within treatment;",
        "p ≤ 0.05 of tests; BH within level, phenotype and model"
      )
    ) +
    theme_figure()
  check_page <- association_summary |>
    filter(level == "protein", universe == "full", model %in% c("pooled", "pooled_mito")) |>
    select(phenotype, model, p_le_0.05, fdr_le_0.05) |>
    pivot_longer(c(p_le_0.05, fdr_le_0.05), names_to = "count") |>
    mutate(
      trait = factor(PHENOTYPE_LABELS[phenotype], rev(PHENOTYPE_LABELS)),
      count = if_else(count == "p_le_0.05", "p ≤ 0.05", "FDR ≤ 0.05"),
      model = if_else(model == "pooled_mito", "with pellet mito share", "without")
    ) |>
    ggplot(aes(value, trait, fill = model)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    facet_wrap(~count, scales = "free_x") +
    scale_fill_manual(
      values = c(without = "grey40", `with pellet mito share` = "grey75"), name = NULL
    ) +
    labs(
      x = "Proteins", y = NULL,
      title = "Check: protein associations with and without pellet mitochondrial share",
      subtitle = "Share centred within group, added to the association model"
    ) +
    theme_figure()
  group_rows <- pooled |>
    filter(level %in% c("module", "complex", "cluster"), universe == "full") |>
    select(level, universe, feature, phenotype, r, adj.P.Val) |>
    left_join(r_group, by = c("level", "universe", "feature"), relationship = "many-to-one")
  heat <- bind_rows(
    distinct(group_rows, level, feature, r = r_group, fdr = group_fdr) |>
      mutate(column = "PAS vs vehicle"),
    transmute(group_rows, level, feature, r, fdr = adj.P.Val, column = PHENOTYPE_LABELS[phenotype])
  ) |>
    mutate(
      column = factor(column, c("PAS vs vehicle", PHENOTYPE_LABELS)),
      panel = if_else(column == "PAS vs vehicle", "Treatment effect", "Phenotype associations"),
      level = LEVEL_LABELS[level]
    )
  heat_pages <- heat |>
    mutate(
      chunk = ceiling(match(feature, unique(feature)) / 32), .by = level
    ) |>
    mutate(level = factor(level, LEVEL_LABELS)) |>
    group_by(level, chunk)
  keys <- group_keys(heat_pages)
  r_heatmaps <- group_split(heat_pages) |>
    set_names(sprintf("%s, page %d", keys$level, keys$chunk)) |>
    map(\(d) {
      ggplot(d, aes(column, feature, fill = r)) +
        geom_tile(colour = "white") +
        geom_tile(
          data = \(x) filter(x, fdr <= ALPHA), fill = NA, colour = "black", linewidth = 0.5
        ) +
        geom_text(aes(label = sprintf("%.2f", r)), size = 1.7) +
        facet_grid(~panel, scales = "free_x", space = "free_x") +
        scale_y_discrete(labels = \(x) {
          str_trunc(enrichVolcano::clean_label(x, width = 200), 45, side = "center")
        }) +
        scale_fill_gradient2(
          low = DIRECTION_COLOURS[["Down"]], high = DIRECTION_COLOURS[["Up"]], limits = c(-1, 1)
        ) +
        labs(
          x = NULL, y = NULL,
          title = str_glue("{label}: {d$level[1]}, treatment effect and phenotype associations"),
          subtitle = paste(
            "Pearson r across this sex's ten mice (treatment: PAS = 1);",
            "box: FDR ≤ 0.05 (limma, BH within level)"
          )
        ) +
        theme_figure() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
    })
  example_long <- examples |>
    rowwise() |>
    reframe(
      strip = sprintf(
        "%s · %s\nr %.2f (VEH %.2f, PAS %.2f)\np %.2g, FDR %.2g",
        PHENOTYPE_LABELS[phenotype], label, r,
        r_veh, r_pas, P.Value, adj.P.Val
      ),
      sample_id = mice, value = fit$E_norm$full[feature, mice],
      phenotype_value = phenotype_of(phenotype)
    ) |>
    left_join(select(targets, sample_id, group), by = "sample_id") |>
    mutate(group = factor(group, FIGURE_GROUPS), strip = strip_factor(strip))
  examples_page <- association_grid(example_long,
    title = str_glue("{label}: strongest robust protein per phenotype"),
    subtitle =
      "Smallest association p that keeps its sign and p ≤ 0.05 with any one mouse left out",
    ncol = 3, nrow = 2
  )

  score_matrix <- rbind(fit$scores$full, fit$scores$mito, fit$scores$module)
  feature_matrix <- \(level, universe) {
    if (level == "protein") fit$E_norm[[universe]] else score_matrix
  }
  hits <- pooled |>
    filter(P.Value <= ALPHA | p_treatment <= ALPHA) |>
    arrange(match(level, names(LEVEL_LABELS)), match(phenotype, PHENOTYPES), P.Value)
  page_key <- paste(hits$level, hits$phenotype)
  association_pages <- split(hits, factor(page_key, unique(page_key))) |>
    unname() |>
    map(\(d) {
      level <- d$level[1]
      p <- d$phenotype[1]
      strips <- sprintf(
        "%s%s\nr %.2f (VEH %.2f, PAS %.2f)\np %.2g pooled, %.2g within", short(d$label, 32),
        if_else(level == "protein" & d$universe == "mito", " (mito)", ""),
        d$r, d$r_veh, d$r_pas, d$P.Value, d$p_treatment
      )
      long <- pmap(list(d$feature, d$universe, strips), \(f, u, s) {
        tibble(
          strip = s, sample_id = mice, value = feature_matrix(level, u)[f, mice],
          phenotype_value = phenotype_of(p)
        )
      }) |>
        list_rbind() |>
        left_join(select(targets, sample_id, group), by = "sample_id") |>
        mutate(group = factor(group, FIGURE_GROUPS), strip = strip_factor(strip))
      effect <- filter(effects, phenotype == !!p)
      pages <- association_grid(long,
        title = str_glue(
          "{label}: {LEVEL_LABELS[level]} associated with {PHENOTYPE_LABELS[p]}, ",
          "p ≤ 0.05 in either model (n = {nrow(d)})"
        ),
        subtitle = sprintf(
          paste0(
            "Phenotype effect of PAS in this sex: %+.0f%% (p %.2g)%s.\n",
            "Lines within vehicle and within PAS (5 mice each), the within-treatment model; ",
            "dashed: all 10 mice, the pooled model"
          ),
          effect$pct_change, effect$P.Value,
          if (effect$P.Value > ALPHA) {
            ", so these associations follow variation between mice, not a PAS response"
          } else {
            ""
          }
        )
      )
      set_names(
        pages,
        sprintf(
          "association grid, %s, %s, page %d", LEVEL_LABELS[level], PHENOTYPE_LABELS[p],
          seq_along(pages)
        )
      )
    }) |>
    unlist(recursive = FALSE)
  plot_index <- write_report(file.path(out$reports, "03_phenotypes.pdf"), c(
    list(
      "phenotype effects" = effects_page, "phenotype association counts by level" = grid_page,
      "mito-share check" = check_page
    ),
    set_names(r_heatmaps, paste("r heatmap,", names(r_heatmaps))),
    set_names(examples_page, "robust per-phenotype examples"),
    association_pages
  ))

  saveRDS(
    list(
      effects = effects, associations = associations, association_summary = association_summary,
      fgsea = fgsea,
      examples = examples, treatment_r = r_group
    ),
    file.path(out$data, "03_phenotypes.rds")
  )
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
        "Per level, universe, phenotype and model: associations at p ≤ 0.05, chance, FDR tiers;",
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
    out$workbook, inputs
  )
  association_summary
}
