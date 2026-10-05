#!/usr/bin/env Rscript
# DIA-NN precursors to a filtered, contaminant-free matrix. The README gives the reason for each
# filter.

library(here)
library(limpa)
library(dplyr)
library(stringr)
library(purrr)
library(readr)
library(tibble)
library(readxl)
library(nanoparquet)
source(here::here("R", "helpers.R"))

BASE <- here("01_Preprocess", "01_Filtering")
DAT <- file.path(BASE, "c_data")
FIG <- file.path(BASE, "b_reports", "01_filter_figures.pdf")
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(FIG), recursive = TRUE, showWarnings = FALSE)

total_signal <- function(m) sum(2^m, na.rm = TRUE)

step_summary <- function(elist, step, signal_start) {
  tibble(
    step = step,
    precursors = nrow(elist),
    protein_groups = n_distinct(elist$genes$Protein.Group),
    pct_signal = 100 * total_signal(elist$E) / signal_start
  )
}

# Built from a vector so a dropped "|" at a line wrap cannot widen the pattern.
CONTAMINANT_FAMILIES <- str_c(
  "^(",
  str_c(
    c(
      "Krt[0-9][0-9a-z-]*", "Krtap[0-9][0-9a-z-]*",
      "Ig[hkl]c", "Igh[gmaed][0-9a-z-]*", "Ig[kl]v[0-9a-z-]*", "Jchain",
      "Hb[ab]-[a-z0-9]+", "Hba", "Hbb"
    ),
    collapse = "|"
  ),
  ")$"
)

MUST_KEEP <- c(
  "Ckm", "Ckmt2", "Cycs", "Mb", "Aldoa", "Cat",
  "Fabp3", "Glud1", "Prdx1", "Prdx2", "Sod1", "Txn"
)

PLASMA <- c(
  "Alb", "Fga", "Fgb", "Fgg", "Trf", "Ttr", "Ahsg", "Apoa1", "Apoa4", "Apob", "Apoe",
  "Hp", "Hpx", "C3", "C4b", "Cfb", "Serpina1a", "Serpina1b", "Serpina1c", "Serpina1d",
  "Serpina1e", "Serpina3k", "Serpinc1", "Serping1", "Plg", "Kng1", "Itih2", "Gc", "Vtn"
)

# Protein.Ids is added for the cRAP check. log = TRUE is the default, so y$E is already log2.
y <- readDIANN(
  here("00_Input", "report.parquet"),
  annotation.columns = c(
    "Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "Proteotypic"
  ),
  q.columns = c(
    "Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "Global.Q.Value", "Global.PG.Q.Value"
  ),
  q.cutoffs = rep(0.01, 5)
)

metadata <- read_csv(here("00_Input", "metadata.csv"), show_col_types = FALSE) |>
  filter(include)

i <- match(metadata$run, colnames(y))
stopifnot(!anyNA(i))
runs_dropped <- setdiff(colnames(y), colnames(y)[i])

y <- y[, i]
colnames(y) <- metadata$sample_id

# A data frame, not a tibble: limma subsets $targets with drop-style indexing.
y$targets <- as.data.frame(metadata)
rownames(y$targets) <- metadata$sample_id

all_na <- sum(rowSums(!is.na(y$E)) == 0)
y <- removeNARows(y)

signal_start <- total_signal(y$E)
precursors_start <- nrow(y)

raw_totals <- read_parquet(
  here("00_Input", "report.parquet"),
  col_select = c("Run", "Precursor.Quantity")
) |>
  summarise(raw_total = sum(Precursor.Quantity), .by = Run) |>
  rename(run = Run) |>
  inner_join(
    read_excel(here("00_Input", "report_stats.xlsx")) |>
      transmute(
        run = str_remove(str_extract(File.Name, "[^\\\\]+$"), "\\.raw$"),
        precursors = Precursors.Identified,
        proteins = Proteins.Identified,
        norm_instability = Normalisation.Instability
      ),
    by = "run"
  ) |>
  inner_join(select(metadata, run, sample_id, group), by = "run") |>
  mutate(raw_rel = raw_total / median(raw_total)) |>
  arrange(raw_rel)

members <- str_split(y$genes$Protein.Ids, ";")
crap_only <- map_lgl(members, \(m) all(str_starts(m, "cRAP-")))
crap_mixed <- str_detect(y$genes$Protein.Ids, "cRAP-") & !crap_only

stopifnot(!any(map_lgl(members, \(m) {
  any(str_remove(m[str_starts(m, "cRAP-")], "^cRAP-") %in% m[!str_starts(m, "cRAP-")])
})))

crap_detected <- tibble(
  protein = y$genes$Protein.Names[crap_only],
  signal = rowSums(2^y$E[crap_only, , drop = FALSE], na.rm = TRUE)
) |>
  summarise(
    precursors = n(), pct_signal = round(100 * sum(signal) / signal_start, 3),
    .by = protein
  ) |>
  arrange(desc(pct_signal))

crap_removed <- n_distinct(y$genes$Protein.Group[crap_only])
crap_precursors <- sum(crap_only)
y <- y[!crap_only, ]

summary_start <- step_summary(y, "cRAP groups removed", signal_start)

y <- filterNonProteotypicPeptides(y)
summary_proteotypic <- step_summary(y, "shared peptides removed", signal_start)

y <- filterCompoundProteins(y)
summary_compound <- step_summary(y, "compound groups removed", signal_start)

gene_symbols <- str_split(y$genes$Genes, ";")
contaminant <- map_lgl(gene_symbols, \(g) isTRUE(all(str_detect(g, CONTAMINANT_FAMILIES))))
is_plasma <- map_lgl(gene_symbols, \(g) isTRUE(all(g %in% PLASMA)))

share <- function(rows) {
  colSums(2^y$E[rows, , drop = FALSE], na.rm = TRUE) / colSums(2^y$E, na.rm = TRUE)
}

contaminant_per_run <- tibble(
  sample_id = colnames(y$E), group = y$targets$group,
  pct_contaminant = 100 * share(contaminant),
  pct_plasma = 100 * share(is_plasma)
)
contaminant_share <- contaminant_per_run |>
  summarise(
    n = n(),
    median_pct_contaminant = round(median(pct_contaminant), 2),
    median_pct_plasma = round(median(pct_plasma), 2),
    .by = group
  )

removed <- tibble(
  symbol = y$genes$Genes[contaminant],
  signal = rowSums(2^y$E[contaminant, , drop = FALSE], na.rm = TRUE)
) |>
  summarise(precursors = n(), pct_signal = 100 * sum(signal) / signal_start, .by = symbol) |>
  arrange(desc(pct_signal))

y <- y[!contaminant, ]
summary_contaminants <- step_summary(y, "contaminants removed", signal_start)

# Catches an over-reaching contaminant rule, and a search whose FASTA duplicates the proteome.
stopifnot(all(MUST_KEEP %in% unlist(str_split(y$genes$Genes, ";"))))

filter_log <- bind_rows(
  summary_start, summary_proteotypic, summary_compound, summary_contaminants
)

saveRDS(
  list(
    y = y, metadata = metadata, signal_start = signal_start,
    precursors_start = precursors_start,
    run_qc = raw_totals, filter_log = filter_log
  ),
  file.path(DAT, "precursors_filtered.rds")
)

write_workbook(
  list(
    filter_log = filter_log,
    contaminants_removed = removed,
    crap_detected = crap_detected,
    contaminant_share = contaminant_share,
    contaminant_per_run = contaminant_per_run,
    run_qc = raw_totals
  ),
  c(
    "Precursors, protein groups and share of signal after each filter.",
    "Contaminant symbols removed, with precursors and share of signal.",
    "cRAP entries the search appended, with precursors and share of signal.",
    "Median contaminant and plasma share of signal per group.",
    "Contaminant and plasma share of signal per mouse, before removal.",
    "Per run: raw signal, precursors and protein groups identified."
  ),
  file.path(DAT, "01_filter.xlsx"),
  c(
    report = here("00_Input", "report.parquet"), metadata = here("00_Input", "metadata.csv"),
    report_stats = here("00_Input", "report_stats.xlsx")
  )
)

pdf(FIG, width = 7, height = 5)
group <- factor(raw_totals$group)
plot(raw_totals$proteins, raw_totals$raw_rel,
  col = as.integer(group), pch = 16,
  xlab = "protein groups identified", ylab = "raw signal / cohort median"
)
legend("topleft", levels(group), col = seq_along(levels(group)), pch = 16, bty = "n")
invisible(dev.off())

tibble(
  runs_dropped = str_c(runs_dropped, collapse = ", "), samples = ncol(y),
  all_na_precursors = all_na, crap_groups = crap_removed, crap_precursors = crap_precursors,
  crap_mixed_groups = sum(crap_mixed),
  plasma_spread = diff(range(contaminant_share$median_pct_plasma))
)
filter_log

sessionInfo()
