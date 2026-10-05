#!/usr/bin/env Rscript
# Which pathways track each phenotype. Lead: fgsea through enrichVolcano on proteins ranked by their
# phenotype t from 02_Differential_Expression/04_Associate, so the answer reads as NES, FDR and
# leading edge like the treatment results. Supplementary: limma on the singscores. Both under the
# three models: pooled slope, females only, male-minus-female slope difference.

source(here::here("R", "helpers.R"))
library(limma)
library(enrichVolcano)
BiocParallel::register(BiocParallel::SerialParam())

out <- stage_paths("03_Pathway_Enrichment", "04_Associate")
inputs <- c(
  association = here(
    "02_Differential_Expression", "04_Associate", "c_data", "protein_association.rds"
  ),
  gene_sets = here("03_Pathway_Enrichment", "00_Gene_Sets", "c_data", "gene_sets.rds"),
  scores = here("03_Pathway_Enrichment", "01_Scores", "c_data", "set_scores.rds"),
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
assoc <- readRDS(inputs[["association"]])$association
gs <- readRDS(inputs[["gene_sets"]])
sc <- readRDS(inputs[["scores"]])
de <- readRDS(inputs[["contrasts"]])
ph <- read_phenotypes(rownames(de$targets))

# Each phenotype-and-model ranking is one "contrast" to enrichVolcano.
ranked <- assoc |>
  transmute(protein, gene,
    contrast = str_glue("{phenotype} | {model}"), logFC = log2fc_per_sd, t,
    P.Value, adj.P.Val
  ) |>
  as_da(species = "Mus musculus")
set.seed(SEED)
all_sets <- c(gs$full, gs$mito)
gsea <- suppressMessages(run_enrichment(ranked, all_sets, tests = "fgsea", min_size = 15))$fgsea
gsea <- dedup_terms(gsea, unlist(unname(all_sets), recursive = FALSE), method = "enrichmentmap")
phenotype_gsea <- gsea@results |>
  separate_wider_delim(contrast, " | ", names = c("phenotype", "model")) |>
  mutate(leading_edge = map_chr(leading_edge, \(x) paste(x, collapse = ";"))) |>
  select(phenotype, model,
    collection = database, term, NES = score, p, padj, size, direction,
    dedup_status, leading_edge
  )

male <- as.numeric(de$targets$sex == "M")
female <- de$targets$sex == "F"
score_fit <- function(m, design, coef) {
  eBayes(lmFit(m, design)) |>
    topTable(coef = coef, number = Inf, sort.by = "none") |>
    rownames_to_column("set") |>
    transmute(set, score_per_sd = logFC, t, P.Value)
}
singscore <- expand_grid(phenotype = PHENOTYPES, universe = names(sc$scores)) |>
  pmap(\(phenotype, universe) {
    m <- sc$scores[[universe]]
    x <- ph[[phenotype]]
    bind_rows(
      mutate(score_fit(m, with_phenotype(de$design, x), "phenotype"), model = "pooled slope"),
      mutate(score_fit(
        m[, female], with_phenotype(de$design[female, c("PAS_F", "VEH_F")], x[female]),
        "phenotype"
      ), model = "females only"),
      mutate(score_fit(m, with_sex_slope(de$design, x, male), "phenotype_male_minus_female"),
        model = "slope difference, male minus female"
      )
    ) |>
      mutate(phenotype = phenotype, universe = universe, .before = 1)
  }) |>
  list_rbind() |>
  mutate(collection = sc$collection[set]) |>
  mutate(adj.P.Val = p.adjust(P.Value, "BH"), .by = c(phenotype, model, collection)) |>
  relocate(phenotype, model, universe, collection)

summary <- bind_rows(
  phenotype_gsea |>
    reframe(tier_counts(p, padj), .by = c(phenotype, model, collection)) |>
    mutate(test = "phenotype fgsea", .before = 1),
  singscore |>
    reframe(tier_counts(P.Value, adj.P.Val), .by = c(phenotype, model, collection)) |>
    mutate(test = "singscore limma", .before = 1)
)

nes_bars <- function(x, title) {
  x <- slice_min(x, p, n = 12, with_ties = FALSE) |> arrange(NES)
  barplot(x$NES,
    names.arg = str_trunc(x$term, 45), horiz = TRUE, las = 1, cex.names = 0.5,
    col = if_else(x$padj <= ALPHA, "#B2182B", "grey70"), xlab = "NES", main = title, cex.main = 0.8
  )
  abline(v = 0)
}
plot_index <- map(PHENOTYPES, \(p) {
  x <- filter(phenotype_gsea, phenotype == p, is.na(dedup_status) | dedup_status != "redundant")
  write_pdf(file.path(out$reports, str_glue("04_associate_{p}.pdf")), list(
    "phenotype fgsea, 12 smallest p per model (red FDR ≤ 0.05)" = \() {
      par(mfrow = c(1, 3), mar = c(4, 17, 3, 1))
      walk(unique(x$model), \(m) nes_bars(filter(x, model == m), str_glue("{p}: {m}")))
    }
  )) |>
    mutate(phenotype = p, .before = 1)
}) |>
  list_rbind()

write_workbook(
  list(
    summary = summary, phenotype_gsea = phenotype_gsea, singscore = singscore,
    plot_index = plot_index
  ),
  c(
    "Per test, phenotype, model and collection: tests, p ≤ 0.05, chance, FDR ≤ 0.10 and ≤ 0.05.",
    "Lead: fgsea on the phenotype t ranking; NES, p, BH within collection, leading edge.",
    "Supplementary: limma on singscores; change per SD, t, p, BH within phenotype and model.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

summary |> filter(fdr_le_0.10 > 0)

sessionInfo()
