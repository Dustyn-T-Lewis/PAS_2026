# How much the protein results depend on choices made once. Each row changes one choice from the
# baseline and refits the three comparison contrasts in the full proteome exactly as 01_Fit does
# (contrastAsCoef, dpcDE, eBayes):
#   detection slope ALT_SLOPES instead of 02_Quantification's (dpcQuant rerun, restricted to the
#   baseline proteins and normalised by the method the full universe chose);
#   measured run depth instead of the surrogate variable, or no adjustment;
#   no per-mouse sample weights;
#   without the re-injected male PAS run (PAS_M1);
#   pellet MitoCarta share added as a covariate.
# Reported per row and contrast: proteins at p <= 0.05, FDR <= 0.10 and FDR <= 0.05, overlap of
# the TOP_N smallest p-values with the baseline, and the p and FDR of every protein at FDR <= 0.10
# in the baseline.

source(here::here("00_Setup", "Shared_functions", "helpers.R"))
library(limpa)

ALT_SLOPES <- c(0.5, 1.0)
TOP_N <- 50

out <- stage_paths("02_Statistics", "07_Sensitivity")
inputs <- c(
  precursors = output_path("01_Processing", "01_Filtering", "precursors_filtered.rds"),
  proteins = output_path("01_Processing", "02_Quantification", "proteins.rds"),
  design = output_path("01_Processing", "03_Design", "design.rds"),
  fit = output_path("02_Statistics", "01_Fit", "fit.rds")
)
filtered <- readRDS(inputs[["precursors"]])
quant <- readRDS(inputs[["proteins"]])
d <- readRDS(inputs[["design"]])
fit <- readRDS(inputs[["fit"]])
targets <- as_tibble(quant$proteins$targets)
run_depth <- run_depth_log2(filtered$run_qc, targets$sample_id)
group_cols <- setdiff(colnames(d$design), "sv1")

baseline_y <- quant$proteins[quant$universes$full, ]
baseline_y$E <- quant$E_norm$full
full_method <- filter(quant$normalization, universe == "full", chosen)$method

requantify <- function(slope) {
  q <- dpcquant_cached(filtered$y, slope, dpcquant_cache(out$data, slope))[rownames(baseline_y), ]
  q$genes <- baseline_y$genes
  q$E <- normalize(q$E, full_method)
  q
}

with_covariate <- \(value) cbind(d$design[, group_cols], covariate = value - mean(value))
slopes <- map(ALT_SLOPES, \(slope) list(y = requantify(slope), design = d$design, weights = TRUE))
variants <- c(
  list(baseline = list(y = baseline_y, design = d$design, weights = TRUE)),
  set_names(slopes, sprintf("slope %.1f", ALT_SLOPES)),
  list(
    `run depth instead of sv1` = list(
      y = baseline_y, design = with_covariate(run_depth), weights = TRUE
    ),
    `no adjustment` = list(y = baseline_y, design = d$design[, group_cols], weights = TRUE),
    `no sample weights` = list(y = baseline_y, design = d$design, weights = FALSE),
    `without PAS_M1` = list(
      y = baseline_y[, targets$sample_id != "PAS_M1"],
      design = d$design[targets$sample_id != "PAS_M1", ], weights = TRUE
    ),
    `pellet share covariate` = list(
      y = baseline_y,
      design = cbind(d$design, mito = quant$pct_mito[targets$sample_id] - mean(quant$pct_mito)),
      weights = TRUE
    )
  )
)

fit_variant <- function(v) {
  map(set_names(COMPARISON_CONTRASTS), \(cn) {
    contrast <- set_names(rep(0, ncol(v$design)), colnames(v$design))
    contrast[group_cols] <- d$contrasts[group_cols, cn]
    reparametrised <- limma::contrastAsCoef(v$design, contrast)$design
    f <- limma::eBayes(dpcDE(v$y, reparametrised, sample.weights = v$weights, plot = FALSE))
    limma::topTable(f, coef = 1, number = Inf, sort.by = "none") |>
      as_tibble(rownames = "feature") |>
      transmute(feature, contrast = cn, logFC, P.Value, adj.P.Val)
  }) |>
    list_rbind()
}
results <- imap(variants, \(v, name) mutate(fit_variant(v), variant = name, .before = 1)) |>
  list_rbind() |>
  left_join(distinct(filter(fit$protein_contrasts, universe == "full"), feature, gene),
    by = "feature"
  )

baseline <- filter(results, variant == "baseline")
top_baseline <- baseline |>
  slice_min(P.Value, n = TOP_N, by = contrast, with_ties = FALSE) |>
  select(contrast, feature)
counts <- results |>
  reframe(tier_counts(P.Value, adj.P.Val), .by = c(variant, contrast)) |>
  left_join(
    results |>
      slice_min(P.Value, n = TOP_N, by = c(variant, contrast), with_ties = FALSE) |>
      semi_join(top_baseline, by = c("contrast", "feature")) |>
      count(variant, contrast, name = "top_shared"),
    by = c("variant", "contrast")
  ) |>
  mutate(variant = factor(variant, names(variants))) |>
  arrange(contrast, variant)
named <- baseline |>
  filter(adj.P.Val <= FDR_LOOSE) |>
  select(contrast, feature) |>
  inner_join(results, by = c("contrast", "feature")) |>
  select(variant, contrast, feature, gene, logFC, P.Value, adj.P.Val)
# The baseline row refits 01_Fit's model, so its p-values should match 01_Fit's.
baseline_vs_fit <- fit$protein_contrasts |>
  filter(universe == "full", contrast %in% COMPARISON_CONTRASTS) |>
  select(contrast, feature, fit_p = P.Value) |>
  inner_join(baseline, by = c("contrast", "feature")) |>
  summarise(max_abs_p_difference = max(abs(fit_p - P.Value)))

plot_index <- preview_pdf(file.path(out$reports, "07_sensitivity.pdf"), list(
  "proteins at p <= 0.05 per variant and contrast" = \() {
    m <- xtabs(p_le_0.05 ~ variant + contrast, counts)
    par(mar = c(10, 4, 2, 1))
    barplot(t(m), beside = TRUE, las = 2, legend.text = colnames(m), ylab = "p <= 0.05")
  }
))

saveRDS(
  list(counts = counts, named = named, top_n = TOP_N),
  file.path(out$data, "sensitivity.rds")
)
write_workbook(
  list(counts = select(counts, -tests), named = p_fdr(named), plot_index = plot_index),
  c(
    str_glue(
      "Per variant and contrast, of the {nrow(baseline_y)} full-proteome proteins: proteins at ",
      "p <= 0.05, FDR <= 0.10, FDR <= 0.05, and how many of the baseline's {TOP_N} smallest ",
      "p-values stay in the variant's {TOP_N} smallest."
    ),
    "Every protein at FDR <= 0.10 in the baseline: log2 FC, p and FDR under each variant.",
    "Which PDF page shows what."
  ),
  out$workbook
)

baseline_vs_fit
counts
