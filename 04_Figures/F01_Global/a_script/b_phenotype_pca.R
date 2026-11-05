# F01 b: the phenotypes, scaled, in principal components, each trait's loading as a coloured
# arrow, and the PERMANOVA for sex and treatment.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))
source(here("04_Figures", "F01_Global", "a_script", "pca_panel.R"))

# One colour per phenotype for loading arrows, away from the group oranges and blues.
PHENOTYPE_COLOURS <- set_names(
  c("#1B9E77", "#7570B3", "#E6AB02", "#A6761D", "#666666", "#111111"), PHENOTYPE_LABELS
)

inputs <- c(pca = output_path("02_Statistics", "04_PCA", "pca.rds"))
pca <- readRDS(inputs[["pca"]])
scores <- pca$phenotype_scores
loadings <- pca$phenotype_loadings |>
  transmute(
    x = PC1 / max(abs(c(PC1, PC2))), y = PC2 / max(abs(c(PC1, PC2))),
    label = PHENOTYPE_LABELS[phenotype]
  )

plot <- pca_panel(
  scores, pca$phenotype_variance, "Phenotype PCA",
  str_glue("{nrow(loadings)} traits, scaled; {nrow(scores)} mice; PERMANOVA top right"),
  pca$phenotype_permanova, loadings, PHENOTYPE_COLOURS
)

save_figure(plot, "F01_Global", "b_phenotype_pca", 57, 64)
list(
  plot = plot,
  data = bind_rows(
    mutate(select(scores, sample_id, group, PC1, PC2), what = "score"),
    mutate(pca$phenotype_loadings, what = "loading", sample_id = phenotype)
  )
)
