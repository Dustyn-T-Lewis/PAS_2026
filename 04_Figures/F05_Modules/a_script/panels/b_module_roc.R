# F05 b: how well each module eigengene alone separates PAS from vehicle mice within a sex, for the
# modules at exact Wilcoxon p <= 0.05 in either sex: pROC's curve and AUC from 03_Classification.
# A module lower in PAS (AUC < 0.5) is drawn the other way up, as 1 - its curve, with AUC 1 - AUC,
# and labelled so.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(
  classification = output_path("02_Statistics", "03_Classification", "classification.rds")
)
cls <- readRDS(inputs[["classification"]])

shown <- cls$roc |>
  filter(level == "module", p <= ALPHA) |>
  arrange(comparison, p) |>
  mutate(
    cell = sprintf("%s, %s", str_to_sentence(comparison), feature),
    cell = factor(cell, cell),
    lower = auc < 0.5,
    text = sprintf(
      "AUC %.2f\n%s in PAS\np = %.3f", if_else(lower, 1 - auc, auc),
      if_else(lower, "lower", "higher"), p
    )
  )
curves <- cls$module_curves |>
  inner_join(select(shown, comparison, feature, cell, lower), by = c("comparison", "feature")) |>
  mutate(
    fpr = if_else(lower, specificity, 1 - specificity),
    tpr = if_else(lower, 1 - sensitivity, sensitivity)
  )
# The area under each curve, closed along the bottom edge.
areas <- curves |>
  group_by(cell, comparison) |>
  reframe(fpr = c(fpr, 1, 0), tpr = c(tpr, 0, 0))

plot <- ggplot(curves, aes(fpr, tpr)) +
  geom_polygon(aes(fill = I(SEX_COLOURS[str_to_sentence(comparison)])), areas, alpha = 0.3) +
  geom_abline(linetype = "dashed", colour = "grey60", linewidth = 0.25) +
  geom_path(linewidth = 0.4, colour = "grey15") +
  geom_text(
    aes(0.97, 0.05, label = text), shown,
    hjust = 1, vjust = 0, size = 1.7, lineheight = 0.9, inherit.aes = FALSE
  ) +
  facet_wrap(~cell, ncol = 3) +
  coord_equal() +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  labs(
    x = "1 − specificity", y = "sensitivity", title = "ROC of module eigengenes, PAS vs vehicle",
    subtitle = sprintf("Eigengene alone, within sex; Wilcoxon p ≤ %s", ALPHA)
  ) +
  theme_figure() +
  theme(strip.text = element_text(size = 5.5), panel.spacing.x = unit(2, "mm"))

save_panel(plot, "F05_Modules", "b_module_roc", 110, 78)
list(plot = plot, data = select(shown, -cell, -text), inputs = inputs)
