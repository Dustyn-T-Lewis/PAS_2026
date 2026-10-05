#!/usr/bin/env Rscript
# One table for every primary test: nominal p ≤ 0.05 against what chance gives, Π where it exists,
# FDR survivors, and limma's estimate of the share of true effects. Reads finished outputs and
# refits nothing.

source(here::here("R", "helpers.R"))

# propTrueNull reads the null share off the p-value histogram, which needs many tests.
MIN_TESTS_PI <- 100

out <- stage_paths("05_Summary", "01_Signal")
xl <- \(path, sheet) readxl::read_excel(here(path), sheet, guess_max = 1e5)
inputs <- c(
  proteins = here("02_Differential_Expression/02_Contrasts/c_data/02_contrasts.xlsx"),
  pathways = here("03_Pathway_Enrichment/02_Contrasts/c_data/02_contrasts.xlsx"),
  modules = here("04_Network/02_Contrasts/c_data/02_contrasts.xlsx"),
  protein_assoc = here("02_Differential_Expression/04_Associate/c_data/04_associate.xlsx"),
  pathway_assoc = here("03_Pathway_Enrichment/04_Associate/c_data/04_associate.xlsx"),
  module_assoc = here("04_Network/04_Associate/c_data/04_associate.xlsx")
)

families <- bind_rows(
  xl(inputs[["proteins"]], "results") |>
    transmute(
      level = "protein", question = "contrast", universe, term = contrast, p = P.Value,
      fdr = adj.P.Val, pi = pi_score
    ),
  xl(inputs[["pathways"]], "fgsea") |>
    transmute(
      level = paste("pathway", collection), question = "contrast", universe, term = contrast,
      p, fdr = padj
    ),
  xl(inputs[["modules"]], "results") |>
    transmute(
      level = "module", question = "contrast", universe = "full", term = contrast,
      p = P.Value, fdr = adj.P.Val
    ),
  xl(inputs[["protein_assoc"]], "association") |>
    transmute(
      level = "protein", question = model, universe = "full", term = phenotype,
      p = P.Value, fdr = adj.P.Val
    ),
  xl(inputs[["pathway_assoc"]], "phenotype_gsea") |>
    transmute(
      level = paste("pathway", collection), question = model, universe = "full",
      term = phenotype, p, fdr = padj
    ),
  xl(inputs[["module_assoc"]], "association") |>
    transmute(
      level = "module", question = model, universe = "full", term = phenotype,
      p = P.Value, fdr = adj.P.Val
    )
)

signal <- families |>
  reframe(
    tier_counts(p, fdr, if (all(is.na(pi))) NULL else pi),
    share_true_effects = if (n() >= MIN_TESTS_PI) 1 - limma::propTrueNull(p) else NA_real_,
    .by = c(level, question, universe, term)
  ) |>
  mutate(excess = p_le_0.05 - expected_by_chance, .after = expected_by_chance)

plot_index <- write_pdf(file.path(out$reports, "01_signal.pdf"), list(
  "nominal hits against chance, protein contrasts" = \() {
    x <- filter(signal, level == "protein", question == "contrast")
    par(mar = c(4, 12, 3, 1))
    barplot(rbind(x$expected_by_chance, x$p_le_0.05),
      beside = TRUE, horiz = TRUE, las = 1,
      names.arg = str_glue("{x$term} ({x$universe})"), cex.names = 0.7,
      col = c("grey80", "#B2182B"), xlab = "proteins at p ≤ 0.05",
      legend.text = c("chance", "observed"), args.legend = list(x = "bottomright", bty = "n"),
      main = "Observed against chance, with true-effect share from propTrueNull"
    )
  }
))

write_workbook(
  list(signal = signal, plot_index = plot_index),
  c(
    "Every family of primary tests: p ≤ 0.05, chance, excess, Π, FDR tiers, true-effect share.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)
signal |> filter(level == "protein", question == "contrast")
sessionInfo()
