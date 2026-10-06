#!/usr/bin/env Rscript
# Precursors to a quantified, MitoCarta-annotated protein matrix, normalised once for the full
# proteome and once for the mitochondrial universe. The README gives the reasons.

# 02_quantification_run.R attaches the helpers, limpa and readxl, and carries the c_data paths,
# quantify_or_load() and the MitoCarta lookup.
source(here::here("00_Trunk", "02_Quantification", "a_script", "02_quantification_run.R"))

REFIT <- isTRUE(as.logical(Sys.getenv("REFIT_DPCQUANT", "FALSE")))
NORMALIZATION <- "quantile"

out <- stage_paths("00_Trunk", "02_Quantification")
pdf(file.path(out$reports, "02_quantification.pdf"), width = 7, height = 5)

y <- readRDS(file.path(FILTER_DATA, "precursors_filtered.rds"))$y

dpc_cn <- dpcCN(y$E)
dpc_on <- dpcON(y$E, robust = TRUE)
plotDPC(dpc_cn)
plotAveVsMis(y)

dpc_parameters <- tibble(
  estimator = c("dpcCN", "dpcON(robust)"),
  intercept = c(dpc_cn$dpc[1], dpc_on$dpc[1]),
  slope = c(dpc_cn$dpc[2], dpc_on$dpc[2])
)

# The detection bar sits at the smallest group size, derived from the targets so dropping a
# mouse cannot leave the code and the stated rule disagreeing.
proteins <- quantify_or_load(y, REFIT)
proteins <- proteins[filterByDetection(proteins, n.samples = min(table(proteins$targets$group))), ]

mitocarta <- read_mitocarta()
proteins$genes$is_mito <- seq_len(nrow(proteins)) %in%
  mitocarta_rows(proteins, mitocarta$symbol, mitocarta)
localization <- function(keys, table) {
  hit <- mitocarta$sub_localization[match(keys, table)]
  hit <- hit[!is.na(hit)]
  if (length(hit)) hit[1] else NA_character_
}
proteins$genes$sub_localization <- map2_chr(
  str_split(proteins$genes$Genes, ";"), str_split(rownames(proteins), ";"),
  \(g, a) coalesce(localization(g, mitocarta$symbol), localization(a, mitocarta$accession))
)

pct_mito <- 100 * colSums(2^proteins$E[proteins$genes$is_mito, ]) / colSums(2^proteins$E)
mito_share <- tibble(
  sample_id = colnames(proteins), treatment = proteins$targets$treatment,
  sex = proteins$targets$sex, group = proteins$targets$group, pct_mito_signal = pct_mito
)
boxplot(pct_mito_signal ~ group,
  data = mito_share, outline = FALSE, ylab = "% signal from MitoCarta proteins"
)
stripchart(pct_mito_signal ~ group, data = mito_share, vertical = TRUE, add = TRUE, pch = 16)

enrichment_tests <- tibble(
  factor = c("treatment", "sex"),
  p_value = c(
    t.test(pct_mito_signal ~ treatment, mito_share)$p.value,
    t.test(pct_mito_signal ~ sex, mito_share)$p.value
  )
)

boxplot(proteins$E,
  las = 2, cex.axis = 0.6, outline = FALSE,
  col = as.integer(factor(proteins$targets$group)), ylab = "log2 protein abundance"
)
normalization_evidence <- tibble(
  sample_id = colnames(proteins), group = proteins$targets$group,
  median = apply(proteins$E, 2, median), iqr = apply(proteins$E, 2, IQR)
)

# Each universe is normalised on its own proteins, so the mitochondrial results do not depend on
# how much non-mitochondrial protein each pellet carried.
universes <- list(full = rep(TRUE, nrow(proteins)), mito = proteins$genes$is_mito)
E_norm <- map(universes, \(keep) normalizeBetweenArrays(proteins$E[keep, ], NORMALIZATION))

group <- factor(proteins$targets$group)
plotMDSUsingSEs(proteins, pch = 16, col = as.integer(group))
legend("bottomleft", levels(group), col = seq_along(levels(group)), pch = 16, bty = "n")
invisible(dev.off())

saveRDS(
  list(
    proteins = proteins, E_norm = E_norm, universes = universes, pct_mito = pct_mito,
    mito_share = mito_share, normalization = NORMALIZATION, dpc_slope = DPC_SLOPE
  ),
  file.path(QUANT_DATA, "proteins.rds")
)

write_workbook(
  list(
    dpc_parameters = dpc_parameters, mito_share = mito_share,
    enrichment_tests = enrichment_tests, normalization_evidence = normalization_evidence
  ),
  c(
    "Detection-probability curve intercept and slope, by estimator.",
    "MitoCarta share of signal per mouse.",
    "Tests of the MitoCarta share against sex and treatment.",
    "Median and IQR of log2 abundance per mouse before normalisation."
  ),
  out$workbook,
  c(
    precursors = file.path(FILTER_DATA, "precursors_filtered.rds"),
    mitocarta = input_path("Mouse.MitoCarta3.0.xls")
  )
)

tibble(
  n_proteins = nrow(proteins), samples = ncol(proteins), mitocarta = sum(proteins$genes$is_mito),
  median_spread_log2 = diff(range(normalization_evidence$median))
)

sessionInfo()
