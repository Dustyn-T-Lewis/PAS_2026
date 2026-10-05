#!/usr/bin/env Rscript
# Groups of proteins that behave alike in a subset of mice, by biclust's plaid model on the
# protein-standardised full matrix. One run says nothing about stability, so plaid runs under 100
# seeds; biclusters whose protein sets overlap at Jaccard 0.5 or more are treated as one, and each
# is reported with how many seeds found it and which groups its mice come from. Exploratory.

source(here::here("R", "helpers.R"))
library(biclust)
# biclust attaches MASS, whose select() masks dplyr's, so dplyr::select is named in full.

N_SEEDS <- 100
JACCARD <- 0.5

out <- stage_paths("05_Summary", "03_Biclusters")
inputs <- c(
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
de <- readRDS(inputs[["contrasts"]])
z <- t(scale(t(de$E_norm$full)))
group <- set_names(de$targets$group, colnames(z))
gene <- set_names(str_split_i(de$genes$Genes, ";", 1), rownames(de$genes))

found <- map(seq_len(N_SEEDS), \(s) {
  set.seed(SEED + s)
  # Plaid seeds its layers with k-means, whose step-limit notice is harmless.
  b <- muffle(biclust(z, method = BCPlaid(), verbose = FALSE), "Quick-TRANSfer")
  map(seq_len(b@Number), \(k) {
    list(seed = s, proteins = rownames(z)[b@RowxNumber[, k]], mice = colnames(z)[b@NumberxCol[k, ]])
  })
}) |>
  list_flatten()
stopifnot(length(found) > 1)

jaccard <- \(a, b) length(intersect(a, b)) / length(union(a, b))
distance <- Vectorize(\(i, j) 1 - jaccard(found[[i]]$proteins, found[[j]]$proteins))
d <- outer(seq_along(found), seq_along(found), distance)
cluster <- cutree(hclust(as.dist(d), "complete"), h = 1 - JACCARD)

biclusters <- map(sort(unique(cluster)), \(k) {
  members <- found[cluster == k]
  proteins <- table(unlist(map(members, "proteins")))
  mice <- table(unlist(map(members, "mice")))
  core_proteins <- names(proteins)[proteins >= length(members) / 2]
  core_mice <- names(mice)[mice >= length(members) / 2]
  tibble(
    bicluster = k, seeds = n_distinct(map_int(members, "seed")), recurrence = seeds / N_SEEDS,
    proteins = length(core_proteins), mice = length(core_mice),
    groups = paste(names(table(group[core_mice])), table(group[core_mice]), collapse = "; "),
    genes = paste(sort(gene[core_proteins]), collapse = ";")
  )
}) |>
  list_rbind() |>
  arrange(desc(recurrence))

plot_index <- write_pdf(file.path(out$reports, "03_biclusters.pdf"), list(
  "bicluster recurrence across seeds" = \() {
    par(mar = c(4, 4, 3, 1))
    plot(biclusters$recurrence, biclusters$proteins,
      pch = 19, log = "y",
      xlab = str_glue("share of {N_SEEDS} seeds finding it"), ylab = "core proteins",
      main = "Plaid biclusters: how often each recurs and how large it is"
    )
  }
))

write_workbook(
  list(biclusters = biclusters, plot_index = plot_index),
  c(
    "Each recurrent bicluster: seeds finding it, core proteins and mice, mice by group, genes.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)
biclusters |> dplyr::select(-genes)
sessionInfo()
