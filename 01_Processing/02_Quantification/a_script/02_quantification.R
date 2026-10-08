# Precursors to a quantified, MitoCarta-annotated protein matrix, normalised once for the full
# proteome and once for the mitochondrial universe. The README gives the reasons. dpcQuant() takes
# about two and a half minutes, so its result is checkpointed in c_data/dpcQuant/ and refitted only
# with REFIT_DPCQUANT=TRUE.

source(here::here("R", "helpers.R"))
library(limpa)

DPC_SLOPE <- 0.7
NORMALIZATION <- "quantile"
REFIT <- isTRUE(as.logical(Sys.getenv("REFIT_DPCQUANT", "FALSE")))

out <- stage_paths("01_Processing", "02_Quantification")
inputs <- c(
  precursors = here("01_Processing", "01_Filtering", "c_data", "precursors_filtered.rds"),
  mitocarta = input_path("Mouse.MitoCarta3.0.xls")
)
checkpoint <- file.path(out$data, sprintf("dpcquant_slope_%s.rds", DPC_SLOPE))
y <- readRDS(inputs[["precursors"]])$y

# A checkpoint built from a different filtering run would report old numbers silently.
fingerprint <- \(y) c(dim(y), sum(y$E, na.rm = TRUE))
if (REFIT || !file.exists(checkpoint)) {
  quantified <- dpcQuant(y, "Protein.Group", dpc.slope = DPC_SLOPE, verbose = FALSE)
  attr(quantified, "source_fingerprint") <- fingerprint(y)
  saveRDS(quantified, checkpoint)
}
quantified <- readRDS(checkpoint)
stopifnot(identical(attr(quantified, "source_fingerprint"), fingerprint(y)))
inputs <- c(inputs, dpcquant_checkpoint = checkpoint)

pdf(file.path(out$reports, "02_quantification.pdf"), width = 7, height = 5)
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
proteins <- quantified[
  filterByDetection(quantified, n.samples = min(table(quantified$targets$group))),
]

# MitoCarta matched on symbol or UniProt accession: symbol alone misses renamed genes.
mitocarta <- readxl::read_excel(inputs[["mitocarta"]], sheet = 2) |>
  transmute(
    symbol = Symbol, accession = UniProt, sub_localization = MitoCarta3.0_SubMitoLocalization
  )
group_symbols <- str_split(proteins$genes$Genes, ";")
group_accessions <- str_split(rownames(proteins), ";")
proteins$genes$is_mito <- map_lgl(group_symbols, \(g) any(g %in% mitocarta$symbol)) |
  map_lgl(group_accessions, \(a) any(a %in% mitocarta$accession))
first_hit <- function(keys, table) {
  hit <- mitocarta$sub_localization[match(keys, table)]
  hit <- hit[!is.na(hit)]
  if (length(hit)) hit[1] else NA_character_
}
proteins$genes$sub_localization <- map2_chr(group_symbols, group_accessions, \(g, a) {
  coalesce(first_hit(g, mitocarta$symbol), first_hit(a, mitocarta$accession))
})

pct_mito <- 100 * colSums(2^proteins$E[proteins$genes$is_mito, ]) / colSums(2^proteins$E)
mito_share <- tibble(
  sample_id = colnames(proteins), treatment = proteins$targets$treatment,
  sex = proteins$targets$sex, group = proteins$targets$group, pct_mito_signal = pct_mito
)
boxplot(pct_mito_signal ~ group,
  data = mito_share, outline = FALSE, ylab = "% signal from MitoCarta proteins"
)
stripchart(pct_mito_signal ~ group, data = mito_share, vertical = TRUE, add = TRUE, pch = 16)

mito_share_tests <- tibble(
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
  file.path(out$data, "proteins.rds")
)
write_workbook(
  list(
    dpc_parameters = dpc_parameters, mito_share = mito_share,
    mito_share_tests = mito_share_tests, normalization_evidence = normalization_evidence
  ),
  c(
    "Detection-probability curve intercept and slope, by estimator.",
    "MitoCarta share of signal per mouse.",
    "Welch t-tests of the MitoCarta share against sex and against treatment.",
    "Median and IQR of log2 abundance per mouse before normalisation."
  ),
  out$workbook, inputs
)

tibble(
  n_proteins = nrow(proteins), samples = ncol(proteins), mitocarta = sum(proteins$genes$is_mito),
  median_spread_log2 = diff(range(normalization_evidence$median))
)

sessionInfo()
