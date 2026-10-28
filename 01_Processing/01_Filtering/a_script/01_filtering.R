# DIA-NN precursors to a filtered, contaminant-free matrix. The README gives the reason for each
# filter.

source(here::here("00_Setup", "Shared_functions", "helpers.R"))
library(limpa)

out <- stage_paths("01_Processing", "01_Filtering")
inputs <- c(
  report = input_path("report.parquet"), metadata = input_path("metadata.csv"),
  report_stats = input_path("report_stats.xlsx")
)

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
# The q-value columns are limpa's set for a search with match-between-runs (vignette, "Data
# input"): every Q.Value and Lib.Q.Value in the report is already below Q_CUTOFF, so DIA-NN
# filtered on the library pass.
y <- readDIANN(
  inputs[["report"]],
  annotation.columns = c(
    "Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "Proteotypic"
  ),
  q.columns = c("Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value"),
  q.cutoffs = Q_CUTOFF
)

metadata <- readr::read_csv(inputs[["metadata"]], show_col_types = FALSE) |>
  filter(include)
raw_quantity <- nanoparquet::read_parquet(
  inputs[["report"]],
  col_select = c("Run", "Precursor.Quantity")
)
report_stats <- readxl::read_excel(inputs[["report_stats"]])

i <- match(metadata$run, colnames(y))
runs_dropped <- setdiff(colnames(y), colnames(y)[i])

y <- y[, i]
colnames(y) <- metadata$sample_id

# A data frame, not a tibble: limma subsets $targets with drop-style indexing.
y$targets <- as.data.frame(metadata)
rownames(y$targets) <- metadata$sample_id

# The subset limpa::removeNARows() performs, written out: in limpa 1.4.0 its EList method indexes
# with an extra empty subscript, which limma's `[.EList` rejects once any row is all missing.
detected <- rowSums(!is.na(y$E)) > 0
all_na <- sum(!detected)
y <- y[detected, ]

signal_start <- total_signal(y$E)
summary_entering <- step_summary(y, "q-value filters passed", signal_start)

raw_totals <- raw_quantity |>
  summarise(raw_total = sum(Precursor.Quantity), .by = Run) |>
  rename(run = Run) |>
  inner_join(
    report_stats |>
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

crap_duplicates_proteome <- any(map_lgl(members, \(m) {
  any(str_remove(m[str_starts(m, "cRAP-")], "^cRAP-") %in% m[!str_starts(m, "cRAP-")])
}))

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

# An over-reaching contaminant rule would remove one of these muscle proteins.
must_keep_missing <- setdiff(MUST_KEEP, unlist(str_split(y$genes$Genes, ";")))

filter_log <- bind_rows(
  summary_entering, summary_start, summary_proteotypic, summary_compound, summary_contaminants
)

saveRDS(
  list(y = y, run_qc = raw_totals, filter_log = filter_log),
  file.path(out$data, "precursors_filtered.rds")
)

group <- factor(raw_totals$group, FIGURE_GROUPS)
plot_index <- preview_pdf(file.path(out$reports, "01_filtering.pdf"), list(
  "precursor intensity per run" = \() limma::plotDensities(y, legend = FALSE),
  "precursors after each filter" = \() {
    par(mar = c(4, 13, 1, 1))
    barplot(rev(filter_log$precursors),
      names.arg = rev(filter_log$step), horiz = TRUE, las = 1, xlab = "precursors"
    )
  },
  "raw signal against protein groups" = \() {
    plot(raw_totals$proteins, raw_totals$raw_rel,
      col = GROUP_COLOURS[as.character(group)], pch = 16,
      xlab = "protein groups identified", ylab = "raw signal / cohort median"
    )
    legend("topleft", levels(group), col = GROUP_COLOURS[levels(group)], pch = 16, bty = "n")
  }
))

write_workbook(
  list(
    filter_log = filter_log,
    contaminants_removed = removed,
    crap_detected = crap_detected,
    contaminant_share = contaminant_share,
    contaminant_per_run = contaminant_per_run,
    run_qc = raw_totals,
    plot_index = plot_index
  ),
  c(
    "Precursors, protein groups and share of signal entering filtering and after each filter.",
    "Contaminant symbols removed, with precursors and share of signal.",
    "cRAP entries the search appended, with precursors and share of signal.",
    "Median contaminant and plasma share of signal per group.",
    "Contaminant and plasma share of signal per mouse, before removal.",
    "Per run: raw signal, precursors and protein groups identified.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

tibble(
  runs_dropped = str_c(runs_dropped, collapse = ", "), samples = ncol(y),
  runs_unmatched = sum(is.na(i)), runs_with_qc = nrow(raw_totals),
  crap_duplicates_proteome = crap_duplicates_proteome,
  must_keep_missing = str_c(must_keep_missing, collapse = ", "),
  all_na_precursors = all_na, crap_groups = crap_removed, crap_precursors = crap_precursors,
  crap_mixed_groups = sum(crap_mixed),
  plasma_spread = diff(range(contaminant_share$median_pct_plasma))
)
filter_log
