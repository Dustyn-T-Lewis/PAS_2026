# F01 d: the proteome of all twenty mice in principal components, PERMANOVA for sex and treatment,
# and the pellet's MitoCarta share fitted as a vector.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(pca = output_path("02_Statistics", "04_PCA", "pca.rds"))
pca <- readRDS(inputs[["pca"]])
scores <- pca$proteome_scores
reach <- 0.8 * max(abs(c(scores$PC1, scores$PC2)))
mito_arrow <- with(pca$mito_vector, tibble(
  x = PC1 * reach, y = PC2 * reach, label = sprintf("MitoCarta share, R² %.2f", r2)
))

plot <- pca_panel(
  scores, pca$proteome_variance, "Proteome PCA", str_glue("{nrow(scores)} mice"),
  permanova_text(pca$permanova), mito_arrow
)

save_panel(plot, "F01_Global", "d_proteome_pca", 60, 62)
list(
  plot = plot,
  data = select(scores, sample_id, group, sex, treatment, PC1, PC2, pct_mito),
  inputs = inputs
)
