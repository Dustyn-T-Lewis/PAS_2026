#!/usr/bin/env Rscript
# dpcQuant() is the slow call in this pipeline, about two and a half minutes per slope. It runs
# from here so 02_quantify.R runs in seconds off the checkpoint and refits only when asked.
# Rebuild both checkpoints with:
#   Rscript 01_Preprocess/02_Quantification/a_script/02_quantify_run.R

library(here)
library(limpa)
library(readxl)
library(dplyr)
library(purrr)
library(stringr)

QUANT_DATA <- here("01_Preprocess", "02_Quantification", "c_data")
FILTER_DATA <- here("01_Preprocess", "01_Filtering", "c_data")

fingerprint <- function(y) c(dim(y), sum(y$E, na.rm = TRUE))

# Symbol and UniProt accession both, for the reason given in the README.
read_mitocarta <- function() {
  read_excel(here("00_Input", "Mouse.MitoCarta3.0.xls"), sheet = 2) |>
    transmute(
      symbol = Symbol, accession = UniProt,
      sub_localization = MitoCarta3.0_SubMitoLocalization
    )
}

mitocarta_rows <- function(p, symbols, mitocarta) {
  accessions <- mitocarta$accession[mitocarta$symbol %in% symbols]
  which(
    map_lgl(str_split(p$genes$Genes, ";"), \(g) any(g %in% symbols)) |
      map_lgl(str_split(rownames(p), ";"), \(a) any(a %in% accessions))
  )
}

quantify_or_load <- function(y, slope, refit = FALSE) {
  path <- file.path(QUANT_DATA, "dpcQuant", sprintf("proteins_slope_%s.rds", slope))
  if (!refit && file.exists(path)) {
    p <- readRDS(path)
    # A checkpoint built from a different filtering run would report old numbers silently.
    stopifnot(identical(attr(p, "source_fingerprint"), fingerprint(y)))
    return(p)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  # dpc.slope, not a fitted dpc object, so limpa re-estimates the intercept at each slope and
  # the second arm tests the slope alone. Re-estimated -7.983 against dpcCN's own -7.937.
  p <- dpcQuant(y, "Protein.Group", dpc.slope = slope, verbose = FALSE)
  attr(p, "source_fingerprint") <- fingerprint(y)
  saveRDS(p, path)
  p
}

if (sys.nframe() == 0L) {
  y <- readRDS(file.path(FILTER_DATA, "precursors_filtered.rds"))$y
  dpc_cn <- dpcCN(y$E)
  for (s in c(0.7, round(dpc_cn$dpc[2], 3))) {
    quantify_or_load(y, s, refit = TRUE)
    message("rebuilt checkpoint for slope ", s)
  }
}
