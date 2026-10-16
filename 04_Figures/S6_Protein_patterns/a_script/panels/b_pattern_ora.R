# S6 b: what each response pattern of panel a holds. Over-representation (fgsea::fora) of each
# pattern's proteins against the detected proteome, full collections, from 02_Enrichment: the
# three sets with the smallest p per pattern, outlined where BH FDR <= 0.05.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"))
enr <- readRDS(inputs[["enrichment"]])

top <- enr$pattern_ora |>
  slice_min(p, n = 3, by = pattern, with_ties = FALSE) |>
  mutate(
    pattern = factor(pattern, FTEST_PATTERNS),
    label = sprintf("%s (%d/%d)", set_label(term, 32), overlap, size),
    row = reorder(paste(pattern, term), -p),
    significant = fdr <= ALPHA
  ) |>
  arrange(pattern, p)

plot <- ggplot(top, aes(-log10(p), row)) +
  geom_col(aes(fill = collection, colour = significant), width = 0.7, linewidth = 0.4) +
  geom_vline(xintercept = -log10(ALPHA), linetype = "dashed", linewidth = 0.25) +
  facet_wrap(~pattern, ncol = 1, scales = "free_y") +
  scale_y_discrete(labels = deframe(distinct(top, row, label))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = DATABASE_COLOURS, name = NULL) +
  scale_colour_manual(
    values = c(`TRUE` = "black", `FALSE` = "transparent"), breaks = "TRUE",
    labels = paste("FDR ≤", ALPHA), name = NULL
  ) +
  labs(
    x = "-log10 p", y = NULL, title = "What each pattern holds",
    subtitle = paste("ORA vs detected proteome; dashed: p =", ALPHA)
  ) +
  theme_figure() +
  theme(
    panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 5.2),
    strip.text = element_text(hjust = 1, size = 5.5), strip.clip = "off",
    legend.box = "vertical"
  ) +
  guides(
    fill = guide_legend(order = 1, nrow = 2),
    colour = guide_legend(order = 2, override.aes = list(fill = "white"))
  )

save_panel(plot, "S6_Protein_patterns", "b_pattern_ora", 64, 145)
list(plot = plot, data = select(top, -row), inputs = inputs)
