# S4 a: proteins associated with each phenotype at p <= 0.05 within each sex, under each
# association model, against the count expected by chance (01_Fit, full proteome).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(fit = output_path("02_Statistics", "01_Fit", "fit.rds"))
ASSOCIATION_LABELS <- c(
  MODEL_LABELS["pooled"],
  pooled_mito = "Pooled + pellet MitoCarta share", MODEL_LABELS["treatment"]
)
counts <- readRDS(inputs[["fit"]])$association_summary |>
  filter(level == "protein", universe == "full") |>
  mutate(
    trait = factor(TRAIT_SHORT[phenotype], TRAIT_SHORT),
    sex = factor(str_to_sentence(sex), c("Female", "Male")),
    model = factor(ASSOCIATION_LABELS[model], ASSOCIATION_LABELS)
  )

plot <- ggplot(counts, aes(trait, p_le_0.05, fill = model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_hline(aes(yintercept = expected_by_chance), linetype = "dashed", linewidth = 0.25) +
  facet_wrap(~sex) +
  scale_fill_manual(
    values = set_names(MODEL_COLOURS[names(ASSOCIATION_LABELS)], ASSOCIATION_LABELS), name = NULL
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL, y = "proteins at p ≤ 0.05", title = "Protein–phenotype associations",
    subtitle = "Full proteome, within sex; dashed: expected by chance"
  ) +
  theme_figure() +
  theme(panel.grid.major.x = element_blank())

save_panel(plot, "S4_Protein_associations", "a_counts", 180, 60)
list(plot = plot, data = counts, inputs = inputs)
