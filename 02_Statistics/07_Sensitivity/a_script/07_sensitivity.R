# How much the protein results depend on choices made once. Each row changes one choice from the
# baseline and refits the three comparison contrasts in the full proteome exactly as 01_Fit does
# (contrastAsCoef, dpcDE, eBayes):
#   detection slope 0.5 or 1.0 instead of 0.7 (dpcQuant rerun, same 1,962 proteins);
#   measured run depth instead of the surrogate variable, or no adjustment;
#   no per-mouse sample weights;
#   without the re-injected male PAS run (PAS_M1);
#   pellet MitoCarta share added as a covariate.
# Reported per row and contrast: proteins at p <= 0.05, FDR <= 0.10 and FDR <= 0.05, overlap of
# the 50 smallest p-values with the baseline, and the p and FDR of every protein at FDR <= 0.10 in
# the baseline.

source(here::here("00_Setup", "Shared_functions", "helpers.R"))
library(limma)
library(limpa)

ALT_SLOPES <- c(0.5, 1.0)
TOP_N <- 50

out <- stage_paths("02_Statistics", "07_Sensitivity")
inputs <- c(
  precursors = here("01_Processing", "01_Filtering", "c_data", "precursors_filtered.rds"),
  proteins = here("01_Processing", "02_Quantification", "c_data", "proteins.rds"),
  design = here("01_Processing", "03_Design", "c_data", "design.rds"),
  fit = here("02_Statistics", "01_Fit", "c_data", "fit.rds")
)
filtered <- readRDS(inputs[["precursors"]])
quant <- readRDS(inputs[["proteins"]])
d <- readRDS(inputs[["design"]])
fit <- readRDS(inputs[["fit"]])
targets <- as_tibble(quant$proteins$targets)
run_depth <- log2(filtered$run_qc$precursors[match(targets$sample_id, filtered$run_qc$sample_id)])
stopifnot(identical(rownames(d$design), targets$sample_id), !anyNA(run_depth))
group_cols <- setdiff(colnames(d$design), "sv1")

baseline_y <- quant$proteins[quant$universes$full, ]
baseline_y$E <- quant$E_norm$full

# dpcQuant at other slopes; each takes about two and a half minutes, so results are kept and
# reused while the filtered precursors are unchanged.
requantify <- function(slope) {
  checkpoint <- file.path(out$data, sprintf("dpcquant_slope_%s.rds", slope))
  fingerprint <- c(dim(filtered$y), sum(filtered$y$E, na.rm = TRUE))
  if (!file.exists(checkpoint) || !identical(attr(readRDS(checkpoint), "fp"), fingerprint)) {
    q <- dpcQuant(filtered$y, "Protein.Group", dpc.slope = slope, verbose = FALSE)
    attr(q, "fp") <- fingerprint
    saveRDS(q, checkpoint)
  }
  q <- readRDS(checkpoint)[rownames(baseline_y), ]
  q$genes <- baseline_y$genes
  q
}

with_covariate <- \(value) cbind(d$design[, group_cols], covariate = value - mean(value))
variants <- list(
  baseline = list(y = baseline_y, design = d$design, weights = TRUE),
  `slope 0.5` = list(y = requantify(0.5), design = d$design, weights = TRUE),
  `slope 1.0` = list(y = requantify(1.0), design = d$design, weights = TRUE),
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

fit_variant <- function(v) {
  map(set_names(COMPARISON_CONTRASTS), \(cn) {
    contrast <- set_names(rep(0, ncol(v$design)), colnames(v$design))
    contrast[group_cols] <- d$contrasts[group_cols, cn]
    reparametrised <- contrastAsCoef(v$design, contrast)$design
    f <- eBayes(dpcDE(v$y, reparametrised, sample.weights = v$weights, plot = FALSE))
    topTable(f, coef = 1, number = Inf, sort.by = "none") |>
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
  summarise(
    p_le_0.05 = sum(P.Value <= ALPHA), fdr_le_0.10 = sum(adj.P.Val <= 0.10),
    fdr_le_0.05 = sum(adj.P.Val <= ALPHA), .by = c(variant, contrast)
  ) |>
  left_join(
    results |>
      slice_min(P.Value, n = TOP_N, by = c(variant, contrast), with_ties = FALSE) |>
      semi_join(top_baseline, by = c("contrast", "feature")) |>
      count(variant, contrast, name = "top50_shared"),
    by = c("variant", "contrast")
  ) |>
  mutate(variant = factor(variant, names(variants))) |>
  arrange(contrast, variant)
named <- baseline |>
  filter(adj.P.Val <= 0.10) |>
  select(contrast, feature) |>
  inner_join(results, by = c("contrast", "feature")) |>
  select(variant, contrast, feature, gene, logFC, P.Value, adj.P.Val)
stopifnot(
  nrow(counts) == length(variants) * length(COMPARISON_CONTRASTS),
  all(counts$top50_shared[counts$variant == "baseline"] == TOP_N)
)
# The baseline row must reproduce 01_Fit.
check <- fit$protein_contrasts |>
  filter(universe == "full", contrast %in% COMPARISON_CONTRASTS) |>
  select(contrast, feature, fit_p = P.Value) |>
  inner_join(baseline, by = c("contrast", "feature"))
stopifnot(isTRUE(all.equal(check$fit_p, check$P.Value)))

plot_index <- preview_pdf(file.path(out$reports, "07_sensitivity.pdf"), list(
  "proteins at p <= 0.05 per variant and contrast" = \() {
    m <- xtabs(p_le_0.05 ~ variant + contrast, counts)
    par(mar = c(10, 4, 2, 1))
    barplot(t(m), beside = TRUE, las = 2, legend.text = colnames(m), ylab = "p <= 0.05")
    abline(h = ALPHA * nrow(baseline_y), lty = 2)
  }
))

saveRDS(
  list(counts = counts, named = named, results = results),
  file.path(out$data, "sensitivity.rds")
)
write_workbook(
  list(counts = counts, named = named, plot_index = plot_index),
  c(
    paste(
      "Per variant and contrast: proteins at p <= 0.05, FDR <= 0.10, FDR <= 0.05, and how many of",
      "the baseline's 50 smallest p-values stay in the variant's 50 smallest."
    ),
    "Every protein at FDR <= 0.10 in the baseline: log2 FC, p and FDR under each variant.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

counts

sessionInfo()
