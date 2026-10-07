#!/usr/bin/env Rscript
# How well each pathway score, module eigengene, complex score and cluster score alone separates
# PAS from vehicle mice within each sex (five against five): AUC (pROC; 1 = higher in PAS) with
# DeLong's 95% interval, left empty at perfect separation where it is undefined, the exact Wilcoxon
# p, BH within sex and level, and r with group beside it. How many features separate perfectly is
# judged against all 252 relabellings of the same ten mice. Proteins are left out: their ROC
# repeats the protein contrast. The sex-difference arc has no classification; the interaction
# contrast is its test.

source(here::here("R", "helpers.R"))

out <- stage_paths("00_Trunk", "09_Classification")
inputs <- c(fit = here("00_Trunk", "07_Fit", "c_data", "fit.rds"))
fit <- readRDS(inputs[["fit"]])
targets <- fit$targets
LEVELS <- c(pathway = "Pathways", module = "Modules", complex = "Complexes", cluster = "Clusters")

levels_of <- distinct(fit$score_contrasts, feature, level, collection)
all_scores <- rbind(fit$scores$full, fit$scores$mito, fit$scores$module)[, targets$sample_id]
features <- map(set_names(names(LEVELS)), \(lvl) {
  all_scores[levels_of$feature[levels_of$level == lvl], , drop = FALSE]
})

roc_table <- function(m, case) {
  y <- factor(case, levels = c(FALSE, TRUE))
  map(rownames(m), \(f) {
    x <- m[f, ]
    roc <- pROC::roc(y, x, levels = levels(y), direction = "<", quiet = TRUE)
    auc <- as.numeric(roc$auc)
    ci <- if (auc %in% c(0, 1)) {
      rep(NA_real_, 3)
    } else {
      as.numeric(pROC::ci.auc(roc, method = "delong"))
    }
    tibble(
      feature = f, auc = auc, auc_low = ci[1], auc_high = ci[3],
      p = wilcox.test(x[case], x[!case], exact = TRUE)$p.value
    )
  }) |>
    list_rbind() |>
    mutate(fdr = p.adjust(p, "BH"))
}

# Features that put every case above (or below) every control, against the same count under each
# of the choose(10, 5) = 252 relabellings: the exact null for that count.
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

arcs <- keep(ARCS, \(a) !is.na(a$sex))
results <- imap(arcs, \(a, arc) {
  in_sex <- targets$sex == a$sex
  case <- targets$treatment[in_sex] == "PAS"
  imap(features, \(m, level) {
    list(
      roc = roc_table(m[, in_sex], case) |> mutate(arc = arc, level = level, .before = 1),
      perfect = perfect_null(m[, in_sex], case) |> mutate(arc = arc, level = level, .before = 1)
    )
  })
})
r_group <- fit$score_contrasts |>
  filter(contrast %in% map_chr(arcs, "contrast")) |>
  transmute(arc = ARC_OF_CONTRAST[contrast], feature, r_group)
roc <- map(results, \(r) list_rbind(map(r, "roc"))) |>
  list_rbind() |>
  left_join(distinct(levels_of, feature, collection), by = "feature") |>
  left_join(r_group, by = c("arc", "feature")) |>
  relocate(collection, .after = level) |>
  arrange(arc, factor(level, names(LEVELS)), p)
perfect <- list_rbind(map(results, \(r) list_rbind(map(r, "perfect"))))
tier_table <- roc |>
  summarise(
    features = n(), p_le_0.05 = sum(p <= ALPHA), expected_by_chance = ALPHA * n(),
    fdr_le_0.05 = sum(fdr <= ALPHA), perfect = sum(auc %in% c(0, 1)), .by = c(arc, level)
  )
# Perfect separations counted from mouse orderings must equal the features with AUC 0 or 1.
stopifnot(identical(
  perfect$observed,
  left_join(select(perfect, arc, level), tier_table, by = c("arc", "level"))$perfect
))

arc_label <- \(x) map_chr(ARCS[x], "label")
perfect_page <- perfect |>
  mutate(level = factor(LEVELS[level], rev(LEVELS)), arc = arc_label(arc)) |>
  ggplot(aes(observed, level)) +
  geom_errorbar(aes(xmin = 0, xmax = null_q95), width = 0.3, colour = "grey60", orientation = "y") +
  geom_point(aes(x = null_median), shape = 124, size = 3, colour = "grey40") +
  geom_point(size = 2) +
  geom_text(aes(label = sprintf("%d, p %.3f", observed, p)), nudge_y = 0.3, hjust = 0, size = 2) +
  facet_wrap(~arc) +
  labs(
    x = "Features separating PAS from vehicle perfectly", y = NULL,
    title = "Perfect separations against chance, PAS vs vehicle within sex",
    subtitle = "Dot: observed; grey: median and 95th percentile over all 252 relabellings"
  ) +
  theme_figure()
auc_page <- roc |>
  mutate(level = factor(LEVELS[level], LEVELS), arc = arc_label(arc)) |>
  ggplot(aes(auc)) +
  geom_histogram(breaks = seq(0, 1, 0.05), fill = "grey65", colour = "white", linewidth = 0.2) +
  facet_grid(level ~ arc, scales = "free_y") +
  labs(
    x = "AUC (1: higher in PAS)", y = "Features", title = "AUC of every feature, PAS vs vehicle"
  ) +
  theme_figure()
plot_index <- write_report(file.path(out$reports, "09_classification.pdf"), list(
  "perfect separations against chance" = perfect_page, "AUC distributions" = auc_page
))

saveRDS(
  list(roc = roc, perfect = perfect, summary = tier_table, features = features),
  file.path(out$data, "classification.rds")
)
write_workbook(
  list(summary = tier_table, roc = roc, perfect_separations = perfect, plot_index = plot_index),
  c(
    paste(
      "Per arc and level: features, at p ≤ 0.05, expected by chance, FDR ≤ 0.05,",
      "perfect separations."
    ),
    paste(
      "Every feature, PAS vs vehicle within sex: AUC, DeLong 95% CI (empty at perfect separation),",
      "exact Wilcoxon p, BH FDR within arc and level, r with group (PAS = 1)."
    ),
    "Per arc and level: perfect separations against all 252 relabellings.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

tier_table

sessionInfo()
