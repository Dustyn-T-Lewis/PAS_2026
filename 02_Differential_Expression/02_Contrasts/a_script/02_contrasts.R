#!/usr/bin/env Rscript
# Which proteins each contrast moved, in the full proteome and in the mitochondrial universe:
# quantile normalisation per universe, limpa's dpcDE with precision weights, moderated t, BH within
# each contrast and universe, and Π beside p and FDR. Writes the enrichVolcano study folders the
# pathway stage reads.

source(here::here("R", "helpers.R"))
library(limma)
library(limpa)

out <- stage_paths("02_Differential_Expression", "02_Contrasts")
inputs <- c(
  proteins = here("01_Preprocess", "02_Quantification", "c_data", "proteins.rds"),
  design = here("02_Differential_Expression", "01_Design", "c_data", "design.rds")
)

quant <- readRDS(inputs[["proteins"]])
proteins <- quant$proteins
d <- readRDS(inputs[["design"]])
stopifnot(
  identical(rownames(d$design), colnames(proteins$E)),
  identical(colnames(d$contrasts), CONTRASTS)
)

E_norm <- map(d$universes, \(keep) normalizeBetweenArrays(proteins$E[keep, ], quant$normalization))
dpc_fits <- imap(d$universes, \(keep, universe) {
  p <- proteins[keep, ]
  p$E <- E_norm[[universe]]
  dpcDE(p, d$design, sample.weights = TRUE, plot = FALSE)
})
fits <- map(dpc_fits, \(f) eBayes(contrasts.fit(f, d$contrasts)))
# dpcDE keeps its precision weights on fit$EList; fit$weights is NULL.
weights <- map(dpc_fits, \(f) f$EList$weights)

results <- imap(fits, \(fit, universe) {
  map(CONTRASTS, \(cn) {
    topTable(fit, coef = cn, number = Inf, sort.by = "none") |>
      rownames_to_column("protein") |>
      transmute(
        universe = universe, contrast = cn, protein, gene = str_split_i(Genes, ";", 1),
        genes = Genes, is_mito, logFC, AveExpr, t, P.Value, adj.P.Val,
        pi_score = pi_score(P.Value, logFC)
      )
  }) |>
    list_rbind()
}) |>
  list_rbind()

summary <- results |>
  reframe(
    tier_counts(P.Value, adj.P.Val, pi_score),
    share_true_effects = 1 - propTrueNull(P.Value),
    .by = c(universe, contrast)
  )

# The headline counts every later stage is read against.
headline <- \(u, cn) summary$fdr_le_0.10[summary$universe == u & summary$contrast == cn]
stopifnot(
  headline("full", "PAS_vs_VEH_female") == 8, headline("full", "sex_in_VEH") == 661,
  headline("full", "sex_in_PAS") == 972
)

# Π ≤ 0.05 lists, split by direction, for the list-based tools (STRING, INDRA).
pi_lists <- results |>
  filter(universe == "full", pi_score <= ALPHA, !is.na(gene), gene != "") |>
  mutate(direction = if_else(logFC > 0, "up", "down")) |>
  select(contrast, direction, protein, gene, logFC, P.Value, adj.P.Val, pi_score) |>
  arrange(contrast, direction, pi_score)

# enrichVolcano reads a study as five tables; the contrast expressions are rebuilt from the
# design's contrast matrix so the two cannot disagree.
stopifnot(all(d$contrasts %in% c(-1, 0, 1)))
contrast_expression <- apply(d$contrasts, 2, \(k) {
  paste(c(paste(names(k)[k > 0], collapse = " + "), names(k)[k < 0]), collapse = " - ")
})
iwalk(d$universes, \(keep, universe) {
  dir <- file.path(out$data, paste0("study_", universe))
  dir.create(dir, showWarnings = FALSE)
  results |>
    filter(universe == !!universe) |>
    transmute(protein, gene, contrast, logFC, t, P.Value, adj.P.Val, AveExpr) |>
    readr::write_csv(file.path(dir, "da_results.csv"))
  as_tibble(E_norm[[universe]], rownames = "protein") |>
    readr::write_csv(file.path(dir, "matrix.csv"))
  as_tibble(weights[[universe]], rownames = "protein") |>
    set_names(c("protein", colnames(E_norm[[universe]]))) |>
    readr::write_csv(file.path(dir, "weights.csv"))
  tibble(sample = colnames(proteins$E), group = proteins$targets$group) |>
    readr::write_csv(file.path(dir, "samples.csv"))
  tibble(name = names(contrast_expression), expression = unname(contrast_expression)) |>
    readr::write_csv(file.path(dir, "contrasts.csv"))
})

volcano <- function(x, title) {
  plot(x$logFC, -log10(x$P.Value),
    pch = 19, cex = 0.45,
    col = case_when(
      x$adj.P.Val <= 0.10 ~ "#B2182B", x$pi_score <= ALPHA ~ "#EF8A62", TRUE ~ "grey70"
    ),
    xlab = "log2 fold change", ylab = "-log10 p", main = title, cex.main = 0.9
  )
  abline(v = 0, h = -log10(ALPHA), lty = 3)
  legend("topleft", c("FDR ≤ 0.10", "Π ≤ 0.05", "other"),
    pch = 19, col = c("#B2182B", "#EF8A62", "grey70"), bty = "n", cex = 0.7
  )
}
plot_index <- map(CONTRASTS, \(cn) {
  pages <- imap(d$universes, \(keep, universe) {
    x <- filter(results, universe == !!universe, contrast == cn)
    s <- filter(summary, universe == !!universe, contrast == cn)
    function() {
      par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
      volcano(x, str_glue("{cn}, {universe} universe"))
      hist(x$P.Value,
        breaks = seq(0, 1, 0.05), col = "grey80", border = "white", xlab = "p",
        main = sprintf(
          "%d of %d at p ≤ 0.05; true-effect share %.2f", s$p_le_0.05, s$tests,
          s$share_true_effects
        ), cex.main = 0.8
      )
      abline(h = s$tests * 0.05, lty = 2, col = "#B2182B")
    }
  })
  names(pages) <- str_glue("volcano and p-value histogram, {names(pages)} universe")
  write_pdf(file.path(out$reports, str_glue("02_contrasts_{cn}.pdf")), pages) |>
    mutate(contrast = cn, .before = 1)
}) |>
  list_rbind()

saveRDS(
  list(
    design = d$design, contrasts = d$contrasts, universes = d$universes, E_norm = E_norm,
    fits = fits, dpc_fits = dpc_fits, weights = weights, results = results, pi_lists = pi_lists,
    pct_mito = quant$pct_mito, targets = proteins$targets, genes = proteins$genes
  ),
  file.path(out$data, "contrasts.rds")
)

write_workbook(
  list(summary = summary, results = results, pi_lists = pi_lists, plot_index = plot_index),
  c(
    "Per universe and contrast: p ≤ 0.05, chance, Π ≤ 0.05, FDR tiers, true-effect share.",
    "Every protein, universe and contrast: log2 fold change, moderated t, p, BH FDR, Π.",
    "Proteins with Π ≤ 0.05 in the full universe, by contrast and direction.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

summary

sessionInfo()
