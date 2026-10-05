#!/usr/bin/env Rscript
# F05 a. Pathways whose link to a phenotype differs by sex: phenotype fgsea, female fit and
# male-minus-female slope, the strongest non-redundant slope differences (03_Pathway_Enrichment/
# 04_Associate).

source(here::here("R", "panels.R"))

PHENOTYPE_SHORT <- c(
  conductance_sol_pyr = "Soleus pyruvate", capillary_density = "Capillary density",
  csa_mean = "CSA", conductance_sol_oct = "Soleus octanoyl"
)
MODEL_SHORT <- c("females only" = "F", "slope difference, male minus female" = "Δ")
term_label <- \(x, n = 40) str_trunc(enrichVolcano::clean_label(x, width = 200), n)
N_TERMS <- 14

inputs <- c(
  phenotype_gsea = here("03_Pathway_Enrichment", "04_Associate", "c_data", "04_associate.xlsx")
)
phenotype_gsea <- inputs[["phenotype_gsea"]] |>
  readxl::read_excel("phenotype_gsea", guess_max = 1e5) |>
  filter(phenotype %in% names(PHENOTYPE_SHORT), model %in% names(MODEL_SHORT))

top_terms <- phenotype_gsea |>
  filter(model == "slope difference, male minus female", dedup_status == "kept") |>
  slice_min(padj, n = N_TERMS, with_ties = FALSE) |>
  pull(term)
columns <- as.vector(t(outer(PHENOTYPE_SHORT, MODEL_SHORT, \(p, m) str_glue("{p} {m}"))))
data <- phenotype_gsea |>
  filter(term %in% top_terms) |>
  mutate(
    column = factor(str_glue("{PHENOTYPE_SHORT[phenotype]} {MODEL_SHORT[model]}"), columns),
    label = term_label(term), database = collection
  ) |>
  mutate(order = mean(NES[model == "slope difference, male minus female"]), .by = term) |>
  arrange(order) |>
  mutate(term = factor(term, unique(term)))
CONTRAST_COLOURS[columns] <- rep(c("#D55E00", "#7B3294"), times = length(PHENOTYPE_SHORT))
panel <- pathway_matrix(data,
  title = "Pathways whose link to phenotype differs by sex",
  subtitle = paste(
    "Proteins ranked by phenotype t, fgsea; columns in pairs, F (female fit) then Δ (male minus",
    "female slope):\nsoleus pyruvate, capillary density, CSA, soleus octanoyl;",
    "black ring: FDR ≤ 0.05"
  ),
  column_names = set_names(rep(c("F", "Δ"), times = length(PHENOTYPE_SHORT)), columns)
)
panel_data <- select(data, phenotype, model, collection, term, NES, p, padj)
panel_note <- "Phenotype fgsea for the strongest slope-difference terms (04_Associate)."
save_panel(panel, "F05", "main", "A_phenotype_links", 150, 110)
