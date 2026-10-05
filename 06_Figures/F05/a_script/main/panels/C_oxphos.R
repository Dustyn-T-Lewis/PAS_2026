#!/usr/bin/env Rscript
# F05 c. Hallmark oxidative phosphorylation in the protein ranking against soleus pyruvate
# conductance: fgsea running score, pooled and slope-difference rankings.

source(here::here("R", "panels.R"))

PHENOTYPE_SHORT <- c(
  conductance_sol_pyr = "Soleus pyruvate", capillary_density = "Capillary density",
  csa_mean = "CSA", conductance_sol_oct = "Soleus octanoyl"
)
MODEL_SHORT <- c("females only" = "F", "slope difference, male minus female" = "Δ")
term_label <- \(x, n = 40) str_trunc(enrichVolcano::clean_label(x, width = 200), n)
OXPHOS <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"

inputs <- c(
  phenotype_gsea = here("03_Pathway_Enrichment", "04_Associate", "c_data", "04_associate.xlsx"),
  association = here("02_Differential_Expression", "04_Associate", "c_data", "04_associate.xlsx"),
  gene_sets = here("03_Pathway_Enrichment", "00_Gene_Sets", "c_data", "gene_sets.rds")
)
read_sheet <- \(file, sheet) readxl::read_excel(inputs[[file]], sheet, guess_max = 1e5)

gene_sets <- readRDS(inputs[["gene_sets"]])
oxphos <- gene_sets$full$Hallmark[[OXPHOS]]
ranking_models <- c(
  "pooled slope" = "Pooled within-group slope",
  "slope difference, male minus female" = "Male minus female slope"
)
oxphos_stats <- read_sheet("phenotype_gsea", "phenotype_gsea") |>
  filter(term == OXPHOS, phenotype == "conductance_sol_pyr", model %in% names(ranking_models)) |>
  select(model, NES, padj)
association <- read_sheet("association", "association") |>
  filter(phenotype == "conductance_sol_pyr", model %in% names(ranking_models))
curves <- imap(ranking_models, \(model_label, model) {
  ranked <- association |>
    filter(model == !!model) |>
    transmute(protein, gene, contrast = model, logFC = log2fc_per_sd, t, P.Value, adj.P.Val) |>
    enrichVolcano::as_da(species = "Mus musculus")
  stats <- sort(set_names(ranked$rank, ranked$gene), decreasing = TRUE)
  stat <- filter(oxphos_stats, model == !!model)
  fgsea::plotEnrichment(oxphos, stats, ticksSize = 0.15)$data |>
    mutate(panel = str_glue("{model_label}\nNES {round(stat$NES, 2)}, FDR {signif(stat$padj, 2)}"))
}) |>
  list_rbind()
panel <- ggplot(curves, aes(rank, ES)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_line(aes(colour = panel), linewidth = 0.5) +
  scale_colour_manual(values = c("grey40", "#7B3294"), name = NULL) +
  labs(
    x = "Protein rank by association t (positive slope to negative)", y = "Enrichment score",
    title = "Oxidative phosphorylation against soleus pyruvate conductance",
    subtitle = "fgsea running score; Hallmark OXPHOS in the protein ranking"
  ) +
  theme_figure() +
  theme(legend.position = "inside", legend.position.inside = c(0.35, 0.25))

panel_data <- oxphos_stats
panel_note <- "OXPHOS against soleus pyruvate conductance, pooled and slope-difference rankings."
save_panel(panel, "F05", "main", "C_oxphos", 110, 90)
