# F01 c: the proteome of every mouse in principal components, the PERMANOVA for sex and
# treatment, and the direction in which the pellet's MitoCarta share rises, fitted as a vector
# (envfit; its R² and p are in the README and in 04_PCA's workbook).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))
source(here("04_Figures", "F01_Global", "a_script", "pca_panel.R"))

inputs <- c(
  pca = output_path("02_Statistics", "04_PCA", "pca.rds"),
  fit = output_path("02_Statistics", "01_Fit", "fit.rds")
)
pca <- readRDS(inputs[["pca"]])
n_proteins <- nrow(readRDS(inputs[["fit"]])$E_norm$full)
scores <- pca$proteome_scores
mito_label <- "Pellet MitoCarta share (envfit)"
mito_arrow <- with(pca$mito_vector, tibble(x = PC1, y = PC2, label = mito_label))

plot <- pca_panel(
  scores, pca$proteome_variance, "Proteome PCA",
  str_glue(
    "{format(n_proteins, big.mark = ',')} proteins; {nrow(scores)} mice; PERMANOVA top right"
  ),
  pca$permanova, mito_arrow, set_names("grey30", mito_label)
) +
  # Panel b beside it already carries the group key.
  guides(fill = "none")

save_figure(plot, "F01_Global", "c_proteome_pca", 57, 64)
list(
  plot = plot,
  data = bind_rows(
    mutate(select(scores, sample_id, group, sex, treatment, PC1, PC2), what = "score"),
    mutate(pca$mito_vector, sample_id = mito_label, what = "envfit arrow")
  )
)
