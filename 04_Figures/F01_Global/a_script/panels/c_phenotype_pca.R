# F01 c: the six phenotypes, scaled, in principal components, with their loadings and PERMANOVA.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(pca = output_path("02_Statistics", "04_PCA", "pca.rds"))
pca <- readRDS(inputs[["pca"]])
scores <- pca$phenotype_scores
reach <- 0.8 * max(abs(c(scores$PC1, scores$PC2)))
loadings <- pca$phenotype_loadings |>
  transmute(
    x = PC1 / max(abs(c(PC1, PC2))) * reach, y = PC2 / max(abs(c(PC1, PC2))) * reach,
    label = PHENOTYPE_LABELS[phenotype]
  )

plot <- pca_panel(
  scores, pca$phenotype_variance, "Phenotype PCA", "Six traits, scaled",
  permanova_text(pca$phenotype_permanova), loadings
) +
  theme(legend.position = "none")

save_panel(plot, "F01_Global", "c_phenotype_pca", 54, 62)
list(
  plot = plot,
  data = bind_rows(
    mutate(select(scores, sample_id, group, PC1, PC2), what = "score"),
    mutate(pca$phenotype_loadings, what = "loading", sample_id = phenotype)
  ),
  inputs = inputs
)
