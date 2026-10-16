# Step 04 of the female and male comparisons: which single scores separate PAS from vehicle.
# Run by 03_Comparisons/<comparison>/a_script/run_comparison.R.

source(here::here("00_Setup", "Shared_functions", "comparisons.R"))

# Classification: which single pathway scores, modules, complexes and clusters separate PAS from
# vehicle mice in this sex, sliced from 03_Classification (AUC, DeLong CI, exact Wilcoxon p, BH
# within level, perfect separations against every relabelling). Female and male comparisons only.
report_classification <- function(comparison) {
  sex_code <- COMPARISONS[[comparison]]$sex
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

  best <- slice_min(roc, p, n = 4, with_ties = FALSE)
  plot_index <- preview_pdf(file.path(out$reports, "04_classification.pdf"), list(
    "AUC of every pathway score" = \() {
      hist(roc$auc[roc$level == "pathway"],
        breaks = seq(0, 1, 0.05), xlab = "AUC (1: higher in PAS)",
        main = str_glue("{label}: pathway scores")
      )
    },
    "four best separators" = \() {
      walk2(best$feature, seq_along(best$feature), \(f, i) {
        level <- best$level[i]
        pROC::plot.roc(case, cls$features[[level]][f, keep],
          levels = c(FALSE, TRUE), quiet = TRUE, add = i > 1, col = i,
          main = str_glue("{label}: four smallest p")
        )
      })
      legend("bottomright", str_trunc(best$feature, 45),
        col = seq_along(best$feature), lwd = 2, bty = "n", cex = 0.6
      )
    }
  ))

  write_workbook(
    list(summary = tier_table, roc = roc, perfect_separations = perfect, plot_index = plot_index),
    c(
      "Per level: features, at p ≤ 0.05, expected by chance, FDR tiers, perfect separations.",
      paste(
        "Every feature, by p: AUC (PAS vs vehicle), DeLong 95% CI (empty at perfect separation),",
        "exact Wilcoxon p, BH FDR within level, r with PAS."
      ),
      str_glue(
        "Per level: perfect separations against all",
        " {perfect$relabellings[1]} relabellings."
      ),
      "Which PDF page shows what."
    ),
    out$workbook, inputs
  )
  tier_table
}
