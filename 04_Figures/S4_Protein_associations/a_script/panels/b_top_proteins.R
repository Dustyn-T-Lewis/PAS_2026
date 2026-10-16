# S4 b: for each phenotype and sex, the protein with the smallest within-treatment p, and its
# r-equivalent under the pooled and within-treatment models (01_Fit, full proteome). A protein that
# keeps its association within treatment tracks the trait among mice given the same treatment.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(fit = output_path("02_Statistics", "01_Fit", "fit.rds"))
associations <- readRDS(inputs[["fit"]])$protein_associations |>
  filter(universe == "full", model %in% names(MODEL_LABELS))
top <- associations |>
  filter(model == "treatment") |>
  slice_min(P.Value, n = 1, by = c(phenotype, sex), with_ties = FALSE) |>
  select(phenotype, sex, feature, gene)
shown <- associations |>
  semi_join(top, by = c("phenotype", "sex", "feature")) |>
  mutate(
    sex = factor(str_to_sentence(sex), c("Female", "Male")),
    row = factor(
      sprintf("%s: %s", PHENOTYPE_LABELS[phenotype], gene),
      rev(unique(sprintf("%s: %s", PHENOTYPE_LABELS[phenotype], gene)[order(phenotype)]))
    ),
    model = factor(MODEL_LABELS[model], MODEL_LABELS)
  )

plot <- ggplot(shown, aes(r_equiv, row)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.25) +
  geom_line(aes(group = row), colour = "grey70", linewidth = 0.4) +
  geom_point(aes(colour = model, shape = adj.P.Val <= ALPHA), size = 1.6) +
  facet_wrap(~sex, scales = "free_y") +
  scale_colour_manual(
    values = set_names(MODEL_COLOURS[names(MODEL_LABELS)], MODEL_LABELS), name = NULL
  ) +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c(`TRUE` = paste("FDR ≤", ALPHA), `FALSE` = paste("FDR >", ALPHA)), name = NULL
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(
    x = "r-equivalent of the limma t", y = NULL, title = "Strongest protein per phenotype",
    subtitle = "Chosen by within-treatment p; pooled and within-treatment estimates"
  ) +
  theme_figure() +
  theme(axis.text.y = element_text(size = 5.3))

save_panel(plot, "S4_Protein_associations", "b_top_proteins", 180, 70)
list(
  plot = plot,
  data = select(
    shown, phenotype, sex, model, feature, gene, estimate, r, r_equiv, P.Value, adj.P.Val
  ),
  inputs = inputs
)
