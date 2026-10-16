# F01 f: each sex's PAS signature (proteins at Π <= 0.05, up and down) placed in the other sex's
# ranking, drawn as limma's barcode (Wu & Smyth 2012): members as ticks, the up list above the
# axis and the down list below, with a tricube moving average of membership along the rank. Each
# lane carries the fry test across mice (02_Enrichment).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"))
enr <- readRDS(inputs[["enrichment"]])
ranks <- list(male = enr$ranks$PAS_vs_VEH_male$full, female = enr$ranks$PAS_vs_VEH_female$full)
SPAN <- 0.12
MIN_HEIGHT <- 0.34

lane_name <- \(source, ranked_in) {
  sprintf("%s signature\nin %s ranking", str_to_sentence(source), ranked_in)
}
lanes <- enr$signatures |>
  mutate(
    lane = factor(
      lane_name(source, ranked_in), lane_name(c("female", "male"), c("male", "female"))
    ),
    side = if_else(direction == "up", 1, -1),
    ranks = map2(members, ranked_in, \(m, target) {
      position <- rank(-ranks[[target]], ties.method = "first")
      unname(position[str_split_1(m, ";")])
    })
  )
n_total <- length(ranks$male)
ticks <- select(lanes, lane, direction, side, rank = ranks) |>
  unnest(rank)
worms <- lanes |>
  mutate(worm = map(ranks, \(r) {
    smooth <- limma::tricubeMovingAverage(as.integer(seq_len(n_total) %in% r), span = SPAN)
    tibble(rank = seq_len(n_total), height = MIN_HEIGHT + (1 - MIN_HEIGHT) * smooth / max(smooth))
  })) |>
  select(lane, direction, side, worm) |>
  unnest(worm) |>
  mutate(y = side * height, y0 = side * MIN_HEIGHT)
labels <- lanes |>
  mutate(
    text = sprintf(
      "%s (n = %d), %.0f%% same way, fry p %.2g", direction, n_ranked, 100 * same_way, p
    ),
    x = if_else(direction == "up", n_total * 0.995, n_total * 0.005),
    hjust = if_else(direction == "up", 1, 0), y = side * 1.38, vjust = if_else(side > 0, 1, 0)
  )

plot <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
  geom_segment(
    aes(rank, side * 0.06, xend = rank, yend = side * 0.3, colour = direction), ticks,
    linewidth = 0.2, alpha = 0.8, show.legend = FALSE
  ) +
  geom_ribbon(
    aes(rank, ymin = y0, ymax = y, fill = direction, colour = direction), worms,
    alpha = 0.25, linewidth = 0, key_glyph = "rect"
  ) +
  geom_line(aes(rank, y, colour = direction), worms, linewidth = 0.4, show.legend = FALSE) +
  geom_text(
    aes(x, y, label = text, hjust = hjust, vjust = vjust, colour = direction), labels,
    size = 1.9, fontface = "bold", show.legend = FALSE
  ) +
  facet_wrap(~lane, ncol = 1, strip.position = "left") +
  scale_colour_manual(
    values = c(up = DIRECTION_COLOURS[["Up"]], down = DIRECTION_COLOURS[["Down"]]),
    labels = c(up = "Up in its own sex", down = "Down in its own sex"),
    aesthetics = c("colour", "fill"), name = NULL,
    guide = guide_legend(nrow = 2, reverse = TRUE)
  ) +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = c(-1.42, 1.42), breaks = NULL, expand = expansion(0)) +
  labs(
    x = sprintf("rank in the other sex (n = %s)", format(n_total, big.mark = ",")),
    y = NULL, title = "Cross-sex signatures",
    subtitle = sprintf("Each sex's Π ≤ %s proteins, ranked in the other sex", ALPHA)
  ) +
  theme_figure() +
  theme(
    # Inside, not outside: patchwork then puts this axis in the same row as its neighbours'.
    strip.placement = "inside",
    strip.text.y.left = element_text(angle = 90, size = 5.5, face = "plain"),
    panel.grid.major = element_blank(), panel.spacing.y = unit(2, "mm")
  )

save_panel(plot, "F01_Global", "f_signature_barcode", 57, 64)
list(plot = plot, data = select(enr$signatures, -members), inputs = inputs)
