# F01 c: the proteome of every mouse in principal components, the PERMANOVA for sex and
# treatment, and the direction in which the pellet's MitoCarta share rises, fitted as a vector.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(
  pca = output_path("02_Statistics", "04_PCA", "pca.rds"),
  fit = output_path("02_Statistics", "01_Fit", "fit.rds")
)
pca <- readRDS(inputs[["pca"]])
n_proteins <- nrow(readRDS(inputs[["fit"]])$E_norm$full)
scores <- pca$proteome_scores
reach <- 0.8 * max(abs(c(scores$PC1, scores$PC2)))
mito_label <- with(
  pca$mito_vector, sprintf("Pellet MitoCarta share rises (R² %.2f, p %.2g)", r2, p)
)
mito_arrow <- with(pca$mito_vector, tibble(x = PC1 * reach, y = PC2 * reach, label = mito_label))

plot <- pca_panel(
  scores, pca$proteome_variance, "Proteome PCA",
  str_glue(
    "{format(n_proteins, big.mark = ',')} proteins; {nrow(scores)} mice; PERMANOVA top right"
  ),
  pca$permanova, mito_arrow, set_names("#1B9E77", mito_label)
)

save_panel(plot, "F01_Global", "c_proteome_pca", 57, 64)
list(
  plot = plot,
  data = select(scores, sample_id, group, sex, treatment, PC1, PC2, pct_mito),
  inputs = inputs
)
