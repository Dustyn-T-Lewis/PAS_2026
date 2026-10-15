# S4 a: proteins associated with each phenotype at p <= 0.05 within each sex, under each
# association model, against the count expected by chance (01_Fit, full proteome).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(fit = output_path("02_Statistics", "01_Fit", "01_fit.xlsx"))
ASSOCIATION_LABELS <- c(
  pooled = "Pooled", pooled_mito = "Pooled + pellet MitoCarta share",
  treatment = "Within treatment"
)
counts <- readxl::read_excel(inputs[["fit"]], "association_summary") |>
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
  scale_fill_manual(values = c("grey40", "grey70", "#7B3294"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL, y = "proteins at p ≤ 0.05", title = "Protein–phenotype associations",
    subtitle = "Full proteome, ten mice per sex; dashed: expected by chance"
  ) +
  theme_figure() +
  theme(panel.grid.major.x = element_blank(), strip.text = element_text(face = "bold"))

save_panel(plot, "S4_Protein_associations", "a_counts", 180, 60)
list(plot = plot, data = counts, inputs = inputs)
