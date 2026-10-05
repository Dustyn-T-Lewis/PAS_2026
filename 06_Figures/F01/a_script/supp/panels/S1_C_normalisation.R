#!/usr/bin/env Rscript
# S1 c. log2 abundance per mouse before and after quantile normalisation.

source(here::here("R", "panels.R"))

inputs <- c(
  proteins = here("01_Preprocess", "02_Quantification", "c_data", "proteins.rds"),
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
de <- readRDS(inputs[["contrasts"]])
targets <- as_tibble(de$targets) |> mutate(group = factor(group, FIGURE_GROUPS))
before <- readRDS(inputs[["proteins"]])$proteins$E[rownames(de$E_norm$full), ]
normalisation <- list("Before (dpcQuant)" = before, "After (quantile)" = de$E_norm$full) |>
  imap(\(m, stage) {
    as_tibble(m, rownames = "protein") |>
      pivot_longer(-protein, names_to = "sample_id", values_to = "abundance") |>
      mutate(stage = stage)
  }) |>
  list_rbind() |>
  left_join(select(targets, sample_id, group), by = "sample_id") |>
  mutate(
    stage = factor(stage, c("Before (dpcQuant)", "After (quantile)")),
    sample_id = factor(sample_id, targets$sample_id[order(targets$group)])
  )
panel <- ggplot(normalisation, aes(sample_id, abundance, fill = group)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2, width = 0.7) +
  facet_wrap(~stage, ncol = 1) +
  scale_fill_manual(values = GROUP_COLOURS, name = NULL) +
  coord_cartesian(ylim = quantile(normalisation$abundance, c(0.01, 0.99), na.rm = TRUE)) +
  labs(
    x = NULL, y = "log2 abundance", title = "Normalisation",
    subtitle = "Full universe; the sex shift before normalisation is removed by it"
  ) +
  theme_figure() +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    strip.text = element_text(size = 6), legend.position = "bottom"
  )
panel_data <- summarise(normalisation,
  median = median(abundance, na.rm = TRUE), iqr = IQR(abundance, na.rm = TRUE),
  .by = c(stage, sample_id, group)
)
panel_note <- "Median and IQR of log2 abundance per mouse before and after quantile normalisation."
save_panel(panel, "F01", "supp", "S1_C_normalisation", 100, 70)
