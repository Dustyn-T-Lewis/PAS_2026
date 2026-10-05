#!/usr/bin/env Rscript
# Precursors to a quantified, MitoCarta-annotated protein matrix, plus the evidence behind the
# normalization choice. The README gives the reasons; normalization itself runs in 02_Differential.

# 02_quantify_run.R attaches here, limpa, readxl, dplyr, purrr and stringr, and carries the
# c_data paths, quantify_or_load() and the MitoCarta lookup.
source(here::here("01_Preprocess", "02_Quantification", "a_script", "02_quantify_run.R"))
source(here::here("R", "helpers.R"))

library(tibble)
library(writexl)

# The detection bar sits at the smallest group size, derived from the targets so dropping a
# mouse cannot leave the code and the stated rule disagreeing.
filter_by_detection <- function(p) {
  p[filterByDetection(p, n.samples = min(table(p$targets$group))), ]
}

REFIT <- isTRUE(as.logical(Sys.getenv("REFIT_DPCQUANT", "FALSE")))
DPC_SLOPE <- 0.7
NORMALIZATION <- "quantile"

BASE <- here("01_Preprocess", "02_Quantification")
FIG <- file.path(BASE, "b_reports", "02_quantify_figures.pdf")
dir.create(file.path(BASE, "c_data"), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(FIG), recursive = TRUE, showWarnings = FALSE)
pdf(FIG, width = 7, height = 5)

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
FITTED_SLOPE <- round(dpc_cn$dpc[2], 3)

proteins <- filter_by_detection(quantify_or_load(y, DPC_SLOPE, REFIT))
proteins_fitted_slope <- filter_by_detection(quantify_or_load(y, FITTED_SLOPE, REFIT))

mitocarta <- read_mitocarta()

annotate_mito <- function(p) {
  p$genes$is_mito <- seq_len(nrow(p)) %in% mitocarta_rows(p, mitocarta$symbol, mitocarta)
  localization <- function(keys, table) {
    hit <- mitocarta$sub_localization[match(keys, table)]
    hit <- hit[!is.na(hit)]
    if (length(hit)) hit[1] else NA_character_
  }
  p$genes$sub_localization <- map2_chr(
    str_split(p$genes$Genes, ";"), str_split(rownames(p), ";"),
    \(g, a) {
      coalesce(localization(g, mitocarta$symbol), localization(a, mitocarta$accession))
    }
  )
  p
}
proteins <- annotate_mito(proteins)
proteins_fitted_slope <- annotate_mito(proteins_fitted_slope)

pct_mito <- 100 * colSums(2^proteins$E[proteins$genes$is_mito, ]) / colSums(2^proteins$E)

mito_share <- tibble(
  sample_id = colnames(proteins),
  treatment = proteins$targets$treatment,
  sex = proteins$targets$sex,
  group = proteins$targets$group,
  pct_mito_signal = pct_mito
)

boxplot(pct_mito_signal ~ group,
  data = mito_share, outline = FALSE, ylab = "% signal from MitoCarta proteins"
)
stripchart(pct_mito_signal ~ group,
  data = mito_share, vertical = TRUE, add = TRUE, pch = 16
)

enrichment_tests <- tibble(
  factor = c("treatment", "sex"),
  p_value = c(
    t.test(pct_mito_signal ~ treatment, mito_share)$p.value,
    t.test(pct_mito_signal ~ sex, mito_share)$p.value
  )
)

boxplot(proteins$E,
  las = 2, cex.axis = 0.6, outline = FALSE,
  col = as.integer(factor(proteins$targets$group)),
  ylab = "log2 protein abundance"
)

normalization_evidence <- tibble(
  sample_id = colnames(proteins),
  group = proteins$targets$group,
  median = apply(proteins$E, 2, median),
  iqr = apply(proteins$E, 2, IQR)
)
median_spread <- diff(range(normalization_evidence$median))

group <- factor(proteins$targets$group)
plotMDSUsingSEs(proteins, pch = 16, col = as.integer(group))
legend("bottomleft", levels(group), col = seq_along(levels(group)), pch = 16, bty = "n")
invisible(dev.off())

saveRDS(
  list(
    proteins = proteins, proteins_fitted_slope = proteins_fitted_slope,
    dpc_cn = dpc_cn, pct_mito = pct_mito, normalization = NORMALIZATION,
    dpc_slope = DPC_SLOPE, fitted_slope = FITTED_SLOPE,
    dpc_parameters = dpc_parameters, mito_share = mito_share,
    enrichment_tests = enrichment_tests, normalization_evidence = normalization_evidence
  ),
  file.path(QUANT_DATA, "proteins.rds")
)

write_workbook(
  list(
    dpc_parameters = dpc_parameters,
    mito_share = mito_share,
    enrichment_tests = enrichment_tests,
    normalization_evidence = normalization_evidence
  ),
  c(
    "Detection-probability curve intercept and slope, by estimator.",
    "MitoCarta share of signal per mouse.",
    "Tests of the MitoCarta share against sex and treatment.",
    "Median and IQR of log2 abundance per mouse before normalisation."
  ),
  file.path(QUANT_DATA, "02_quantify.xlsx"),
  c(
    precursors = file.path(FILTER_DATA, "precursors_filtered.rds"),
    mitocarta = here("00_Input", "Mouse.MitoCarta3.0.xls")
  )
)

tibble(
  fitted_slope = FITTED_SLOPE, preset_slope = DPC_SLOPE,
  protein_groups = nrow(proteins), samples = ncol(proteins),
  mitocarta = sum(proteins$genes$is_mito), pct_mitocarta = 100 * mean(proteins$genes$is_mito),
  median_spread_log2 = median_spread, normalization = NORMALIZATION
)
mutate(enrichment_tests, p_value = format.pval(p_value, digits = 2))

sessionInfo()
