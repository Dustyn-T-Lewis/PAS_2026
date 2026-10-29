# How well each pathway score, module eigengene, complex score and cluster score alone separates
# PAS from vehicle mice within each sex (five against five): AUC (pROC; 1 = higher in PAS) with
# DeLong's 95% interval, left empty at perfect separation where it is undefined, the exact Wilcoxon
# p, BH within sex and level, and r with group beside it. How many features separate perfectly is
# judged against every relabelling of the same ten mice. Proteins are left out: their ROC
# repeats the protein contrast. The sex-difference comparison has no classification; the interaction
# contrast is its test.

source(here::here("00_Setup", "Shared_functions", "helpers.R"))

out <- stage_paths("02_Statistics", "03_Classification")
inputs <- c(fit = output_path("02_Statistics", "01_Fit", "fit.rds"))
fit <- readRDS(inputs[["fit"]])
targets <- fit$targets
LEVELS <- c("pathway", "module", "complex", "cluster")

levels_of <- distinct(fit$score_contrasts, feature, level, collection)
all_scores <- rbind(fit$scores$full, fit$scores$mito, fit$scores$module)[, targets$sample_id]
features <- map(set_names(LEVELS), \(lvl) {
  all_scores[levels_of$feature[levels_of$level == lvl], , drop = FALSE]
})

# The mice of one sex, and which of them were given PAS.
sex_mice <- function(comparison) {
  in_sex <- targets$sex == COMPARISONS[[comparison]]$sex
  list(in_sex = in_sex, case = targets$treatment[in_sex] == "PAS")
}

# Each feature's ROC (kept for the curves) and its AUC, DeLong interval and exact Wilcoxon p.
roc_table <- function(m, case) {
  y <- factor(case, levels = c(FALSE, TRUE))
  rocs <- map(set_names(rownames(m)), \(f) {
    pROC::roc(y, m[f, ], levels = levels(y), direction = "<", quiet = TRUE)
  })
  table <- imap(rocs, \(roc, f) {
    auc <- as.numeric(roc$auc)
    ci <- if (auc %in% c(0, 1)) {
      rep(NA_real_, 3)
    } else {
      as.numeric(pROC::ci.auc(roc, method = "delong"))
    }
    tibble(
      feature = f, auc = auc, auc_low = ci[1], auc_high = ci[3],
      p = wilcox.test(m[f, case], m[f, !case], exact = TRUE)$p.value
    )
  }) |>
    list_rbind() |>
    mutate(fdr = p.adjust(p, "BH"))
  list(table = table, rocs = rocs)
}

# Features that put every case above (or below) every control, against the same count under each
# relabelling of the mice (choose(10, 5) for five against five): the exact null for that count.
perfect_null <- function(m, case) {
  keys <- \(idx) apply(idx, 1, \(i) paste(sort(i), collapse = ","))
  n <- sum(case)
  ord <- t(apply(m, 1, order))
  top <- keys(ord[, (ncol(m) - n + 1):ncol(m), drop = FALSE])
  bottom <- keys(ord[, seq_len(n), drop = FALSE])
  count <- \(labelled) sum(top == labelled | bottom == labelled)
  null <- map_int(combn(ncol(m), n, simplify = FALSE), \(i) count(paste(sort(i), collapse = ",")))
  observed <- count(paste(which(case), collapse = ","))
  tibble(
    observed = observed, null_median = median(null), null_q95 = unname(quantile(null, 0.95)),
    relabellings = length(null), p = mean(null >= observed)
  )
}

results <- map(set_names(names(SEX_COMPARISONS)), \(comparison) {
  mice <- sex_mice(comparison)
  imap(features, \(m, level) {
    m <- m[, mice$in_sex]
    tested <- roc_table(m, mice$case)
    list(
      roc = mutate(tested$table, comparison = comparison, level = level, .before = 1),
      curves = tested$rocs,
      perfect = perfect_null(m, mice$case) |>
        mutate(comparison = comparison, level = level, .before = 1)
    )
  })
})
r_group <- fit$score_contrasts |>
  filter(contrast %in% map_chr(SEX_COMPARISONS, "contrast")) |>
  transmute(comparison = COMPARISON_OF_CONTRAST[contrast], feature, r_group)
roc <- map(results, \(r) list_rbind(map(r, "roc"))) |>
  list_rbind() |>
  left_join(distinct(levels_of, feature, collection), by = "feature") |>
  left_join(r_group, by = c("comparison", "feature")) |>
  relocate(collection, .after = level) |>
  arrange(comparison, factor(level, LEVELS), p)
perfect <- list_rbind(map(results, \(r) list_rbind(map(r, "perfect"))))
tier_table <- roc |>
  reframe(tier_counts(p, fdr), perfect = sum(auc %in% c(0, 1)), .by = c(comparison, level)) |>
  select(comparison, level, features = tests, everything(), -pi_le_0.05)

# pROC's own curves for the four pathway scores with the smallest p in one sex.
roc_curves <- \(comparison) {
  \() {
    best <- roc |>
      filter(comparison == !!comparison, level == "pathway") |>
      slice_min(p, n = 4, with_ties = FALSE)
    walk2(best$feature, seq_along(best$feature), \(f, i) {
      pROC::plot.roc(results[[comparison]]$pathway$curves[[f]],
        add = i > 1, col = i, main = COMPARISONS[[comparison]]$label
      )
    })
    legend("bottomright", str_trunc(best$feature, 45),
      col = seq_along(best$feature), lwd = 2,
      bty = "n", cex = 0.6
    )
  }
}
plot_index <- preview_pdf(file.path(out$reports, "03_classification.pdf"), list(
  "AUC of every pathway score" = \() {
    par(mfrow = c(1, 2))
    walk(names(SEX_COMPARISONS), \(cmp) {
      hist(roc$auc[roc$comparison == cmp & roc$level == "pathway"],
        breaks = seq(0, 1, 0.05), main = COMPARISONS[[cmp]]$label, xlab = "AUC (1: higher in PAS)"
      )
    })
  },
  "best pathway separators, female" = roc_curves("female"),
  "best pathway separators, male" = roc_curves("male")
))

saveRDS(
  list(roc = roc, perfect = perfect, summary = tier_table, features = features),
  file.path(out$data, "classification.rds")
)
write_workbook(
  list(summary = tier_table, roc = roc, perfect_separations = perfect, plot_index = plot_index),
  c(
    paste(
      "Per comparison and level: features, at p ≤ 0.05, FDR tiers,",
      "perfect separations."
    ),
    paste(
      "Every feature, PAS vs vehicle within sex: AUC, DeLong 95% CI (empty at perfect separation),",
      "exact Wilcoxon p, BH FDR within comparison and level, r with group (PAS = 1)."
    ),
    str_glue(
      "Per comparison and level: perfect separations against all",
      " {perfect$relabellings[1]} relabellings."
    ),
    "Which PDF page shows what."
  ),
  out$workbook
)

tier_table
perfect
