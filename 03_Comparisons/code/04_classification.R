# Step 04 of the female and male comparisons: which single scores separate PAS from vehicle.
# Run by 03_Comparisons/<comparison>/a_script/run_comparison.R.

source(here::here("00_Setup", "Shared_functions", "comparisons.R"))

# ROC curve points for one feature, drawn on the side it separates.
roc_points <- function(x, case, auc) {
  r <- pROC::roc(
    factor(case, c(FALSE, TRUE)), x,
    direction = if (auc >= 0.5) "<" else ">", quiet = TRUE
  )
  tibble(fpr = 1 - r$specificities, tpr = r$sensitivities) |> arrange(fpr, tpr)
}

# Classification: which single pathway scores, modules, complexes and clusters separate PAS from
# vehicle mice in this sex, sliced from 03_Classification (AUC, DeLong CI, exact Wilcoxon p, BH
# within level, perfect separations against all 252 relabellings). Female and male comparisons only.
report_classification <- function(comparison) {
  sex_code <- COMPARISONS[[comparison]]$sex
  stopifnot(!is.na(sex_code))
  out <- comparison_paths(comparison, "04_classification")
  inputs <- c(
    classification = output_path("02_Statistics", "03_Classification", "classification.rds"),
    fit = output_path("02_Statistics", "01_Fit", "fit.rds")
  )
  cls <- readRDS(inputs[["classification"]])
  label <- COMPARISONS[[comparison]]$label
  targets <- readRDS(inputs[["fit"]])$targets
  keep <- targets$sex == sex_code
  case <- targets$treatment[keep] == "PAS"
  roc <- filter(cls$roc, comparison == !!comparison) |>
    select(-comparison)
  perfect <- filter(cls$perfect, comparison == !!comparison) |>
    select(-comparison)
  tier_table <- filter(cls$summary, comparison == !!comparison) |>
    select(-comparison)

  perfect_page <- perfect |>
    mutate(level = factor(LEVEL_LABELS[-1][level], rev(LEVEL_LABELS[-1]))) |>
    ggplot(aes(observed, level)) +
    geom_errorbar(
      aes(xmin = 0, xmax = null_q95),
      width = 0.3, colour = "grey60", orientation = "y"
    ) +
    geom_point(aes(x = null_median), shape = 124, size = 3, colour = "grey40") +
    geom_point(size = 2, colour = comparison_colour(comparison)) +
    geom_text(aes(label = sprintf("p %.3f", p)), nudge_y = 0.3, hjust = 0, size = 2) +
    labs(
      x = "Features separating PAS from vehicle perfectly", y = NULL,
      title = str_glue("{label}: perfect separations against chance"),
      subtitle = paste(
        "Dot: observed; grey: median and 95th percentile",
        "over all 252 relabellings of the same mice"
      )
    ) +
    theme_figure()
  auc_page <- roc |>
    mutate(level = factor(LEVEL_LABELS[-1][level], LEVEL_LABELS[-1])) |>
    ggplot(aes(auc)) +
    geom_histogram(breaks = seq(0, 1, 0.05), fill = "grey65", colour = "white", linewidth = 0.2) +
    facet_wrap(~level, scales = "free_y", ncol = 1) +
    labs(
      x = "AUC (1: higher in PAS)", y = "Features",
      title = str_glue("{label}: AUC of every feature")
    ) +
    theme_figure()
  curve_pages <- roc |>
    filter(p <= ALPHA) |>
    group_by(level) |>
    group_split() |>
    map(\(d) {
      level <- d$level[1]
      m <- cls$features[[level]][, keep, drop = FALSE]
      long <- pmap(d, \(feature, auc, p, fdr, r_group, ...) {
        roc_points(m[feature, ], case, auc) |>
          mutate(strip = sprintf(
            "%s\nAUC %.2f (%s), r %.2f\np %.3f, FDR %.2g", short(feature, 24), auc,
            if (auc >= 0.5) "higher" else "lower", r_group, p, fdr
          ))
      }) |>
        list_rbind() |>
        mutate(strip = strip_factor(strip))
      pages <- paginate(long, \(x) {
        ggplot(x, aes(fpr, tpr)) +
          geom_abline(linetype = "dashed", linewidth = 0.3, colour = "grey60") +
          geom_step(direction = "vh", colour = comparison_colour(comparison), linewidth = 0.5) +
          facet_wrap(~strip, ncol = 5, nrow = GRID_PER_PAGE / 5, drop = FALSE) +
          coord_equal() +
          labs(
            x = "1 - specificity", y = "Sensitivity",
            title = str_glue(
              "{label}: {LEVEL_LABELS[-1][level]}, PAS vs vehicle, ",
              "every feature at p ≤ 0.05 (n = {nrow(d)})"
            ),
            subtitle = paste(
              "Exact Wilcoxon p, BH within level; higher or lower in PAS;",
              "r with PAS (PAS = 1)"
            )
          ) +
          theme_figure() +
          theme(
            strip.text = element_text(size = 4.6, face = "plain"),
            axis.text = element_text(size = 4.5)
          )
      })
      set_names(
        pages,
        sprintf("ROC grid, %s, page %d", LEVEL_LABELS[-1][level], seq_along(pages))
      )
    }) |>
    unlist(recursive = FALSE)
  plot_index <- write_report(file.path(out$reports, "04_classification.pdf"), c(
    list("perfect separations against chance" = perfect_page, "AUC distributions" = auc_page),
    curve_pages
  ))

  saveRDS(
    list(roc = roc, perfect = perfect, summary = tier_table),
    file.path(out$data, "04_classification.rds")
  )
  write_workbook(
    list(summary = tier_table, roc = roc, perfect_separations = perfect, plot_index = plot_index),
    c(
      "Per level: features, at p ≤ 0.05, expected by chance, FDR ≤ 0.05, perfect separations.",
      paste(
        "Every feature, by p: AUC (PAS vs vehicle), DeLong 95% CI (empty at perfect separation),",
        "exact Wilcoxon p, BH FDR within level, r with PAS."
      ),
      "Per level: perfect separations against all 252 relabellings.",
      "Which PDF page shows what."
    ),
    out$workbook, inputs
  )
  tier_table
}
