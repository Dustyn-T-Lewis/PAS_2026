# S6 b: for each phenotype and sex, the three proteins with the smallest within-treatment p, and
# their r-equivalents under the pooled and within-treatment models (01_Fit, full proteome). A
# protein that keeps its association within treatment tracks the trait among mice given the same
# treatment.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(fit = output_path("02_Statistics", "01_Fit", "fit.rds"))
TOP_N <- 3
# The association models in greys, darkest for within treatment, so a model never reads as a
# direction, a sex or a contrast.
MODEL_COLOURS <- c(pooled = "grey60", treatment = "grey15")
DRAWN_MODELS <- c("pooled", "treatment")
associations <- readRDS(inputs[["fit"]])$protein_associations |>
  filter(universe == "full", model %in% DRAWN_MODELS)
top <- associations |>
  filter(model == "treatment") |>
  slice_min(P.Value, n = TOP_N, by = c(phenotype, sex), with_ties = FALSE) |>
  mutate(rank = row_number(), .by = c(phenotype, sex)) |>
  select(phenotype, sex, feature, rank)
# One row per phenotype, sex and protein, smallest p on top; each sex's column keeps its own rows.
shown <- associations |>
  inner_join(top, by = c("phenotype", "sex", "feature")) |>
  mutate(
    trait = factor(PHENOTYPE_LABELS[phenotype], PHENOTYPE_LABELS),
    sex = factor(str_to_sentence(sex), c("Female", "Male")),
    model = factor(MODEL_LABELS[model], MODEL_LABELS[DRAWN_MODELS]),
    row = paste(phenotype, sex, feature),
    label = protein_label(gene, feature)
  ) |>
  arrange(trait, sex, desc(rank)) |>
  mutate(row = factor(row, unique(row)))

plot <- ggplot(shown, aes(r_equiv, row)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.25) +
  geom_line(aes(group = row), colour = "grey75", linewidth = 0.4) +
  geom_point(aes(colour = model, shape = adj.P.Val <= ALPHA),
    size = 1.5, stroke = 0.5, position = position_dodge(width = 0.4, orientation = "y")
  ) +
  ggh4x::facet_grid2(trait ~ sex, scales = "free_y", independent = "y", switch = "y") +
  scale_y_discrete(labels = set_names(shown$label, shown$row)) +
  scale_colour_manual(
    values = set_names(MODEL_COLOURS[DRAWN_MODELS], levels(shown$model)),
    name = NULL
  ) +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c(`TRUE` = paste("FDR ≤", ALPHA), `FALSE` = paste("FDR >", ALPHA)), name = NULL
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(
    x = "r-equivalent of the limma t", y = NULL,
    title = sprintf("Top %d proteins per phenotype", TOP_N),
    subtitle = "Chosen by within-treatment p; pooled and within-treatment estimates"
  ) +
  theme_figure() +
  theme(
    axis.text.y = element_text(size = 5.3), panel.grid.major.y = element_blank(),
    strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1),
    panel.spacing.y = unit(1.2, "mm"), panel.spacing.x = unit(4, "mm"),
    panel.border = element_rect(fill = NA, colour = "grey85", linewidth = 0.25)
  )

save_figure(plot, "F04_Phenotype_pathways", "S6b_top_proteins", 180, 100)
list(
  plot = plot,
  data = shown |>
    arrange(trait, sex, rank, model) |>
    select(phenotype, sex, rank, model, feature, gene, r_equiv, p = P.Value, fdr = adj.P.Val)
)
