#!/usr/bin/env Rscript
# dpcQuant() is the slow call in this pipeline, about two and a half minutes. It runs from here so
# 02_quantification.R runs in seconds off the checkpoint and refits only when asked. Rebuild with:
#   Rscript 00_Trunk/02_Quantification/a_script/02_quantification_run.R

source(here::here("R", "helpers.R"))
library(limpa)
library(readxl)

QUANT_DATA <- here("00_Trunk", "02_Quantification", "c_data")
FILTER_DATA <- here("00_Trunk", "01_Filtering", "c_data")
DPC_SLOPE <- 0.7

fingerprint <- function(y) c(dim(y), sum(y$E, na.rm = TRUE))

# Symbol and UniProt accession both, for the reason given in the README.
read_mitocarta <- function() {
  read_excel(input_path("Mouse.MitoCarta3.0.xls"), sheet = 2) |>
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

quantify_or_load <- function(y, refit = FALSE) {
  path <- file.path(QUANT_DATA, "dpcQuant", sprintf("proteins_slope_%s.rds", DPC_SLOPE))
  if (!refit && file.exists(path)) {
    p <- readRDS(path)
    # A checkpoint built from a different filtering run would report old numbers silently.
    stopifnot(identical(attr(p, "source_fingerprint"), fingerprint(y)))
    return(p)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  p <- dpcQuant(y, "Protein.Group", dpc.slope = DPC_SLOPE, verbose = FALSE)
  attr(p, "source_fingerprint") <- fingerprint(y)
  saveRDS(p, path)
  p
}

if (sys.nframe() == 0L) {
  quantify_or_load(readRDS(file.path(FILTER_DATA, "precursors_filtered.rds"))$y, refit = TRUE)
  message("rebuilt checkpoint for slope ", DPC_SLOPE)
}
