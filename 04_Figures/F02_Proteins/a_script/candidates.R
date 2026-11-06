# The candidate drawings F02's panels share, sourced by its panel scripts: the STRING network and
# every mouse for the strongest candidates.

GROUP_SHORT <- c(VEH_F = "F VEH", PAS_F = "F PAS", VEH_M = "M VEH", PAS_M = "M PAS")
# The every-mouse panels show the candidates passing at least this many routes; the other
# candidate thresholds come with candidates.rds from 02_Statistics/06_Candidates.
CANDIDATE_DOT_ROUTES <- 3

# One sex's candidates (those passing community_routes or more routes) as a STRING network:
# shaded hulls for communities, each named beside it by its ORA term; edges at the STRING score
# 06_Candidates used; node size routes passed, colour log2 fold change. Communities are kept
# together by adding invisible within-community edges to the layout only. Proteins without a
# STRING partner are left out of the drawing, counted under it and flagged by has_partner in the
# data, with the routes each candidate passed (route_*). cand is candidates.rds from
# 02_Statistics/06_Candidates. Returns the plot and the candidates drawn.
candidate_network <- function(cand, sex, title) {
  limits <- cand$thresholds
  shown <- filter(cand$candidates, sex == !!sex, routes >= limits$community_routes)
  nodes <- shown |>
    distinct(gene, .keep_all = TRUE) |>
    mutate(group = if_else(is.na(name), NA_character_, str_wrap(set_label(name), 24)))
  links <- cand$edges |>
    filter(sex == !!sex, a %in% nodes$gene, b %in% nodes$gene) |>
    distinct(a, b)
  linked <- nodes |> filter(!is.na(community) | gene %in% c(links$a, links$b))
  alone <- setdiff(nodes$gene, linked$gene)
  within <- linked |>
    filter(!is.na(community)) |>
    reframe(as_tibble(t(combn(gene, 2)), .name_repair = \(x) c("a", "b")), .by = community) |>
    select(a, b)
  layout_graph <- igraph::graph_from_data_frame(
    bind_rows(mutate(links, w = 1), mutate(within, w = 4)),
    directed = FALSE, vertices = select(linked, gene)
  )
  set.seed(SEED)
  xy <- igraph::layout_with_fr(layout_graph, weights = igraph::E(layout_graph)$w)
  positions <- tibble(gene = igraph::V(layout_graph)$name, x = xy[, 1], y = xy[, 2]) |>
    left_join(linked, by = "gene")
  segments <- links |>
    left_join(select(positions, a = gene, x, y), by = "a") |>
    left_join(select(positions, b = gene, xend = x, yend = y), by = "b")
  # Gene and community names repelled together, so neither covers a node or the other; each
  # community name starts above its hull's highest member.
  labels <- bind_rows(
    transmute(positions, x, y, label = gene, size = 1.9, colour = "black"),
    positions |>
      filter(!is.na(group)) |>
      summarise(x = mean(x), y = max(y), .by = group) |>
      transmute(x, y, label = group, size = 1.6, colour = "grey30")
  )
  limit <- max(abs(nodes$logFC))
  plot <- ggplot(positions, aes(x, y)) +
    ggforce::geom_mark_hull(
      aes(fill = group, filter = !is.na(group)),
      expand = unit(3.5, "mm"), radius = unit(3, "mm"), alpha = 0.18, colour = NA
    ) +
    geom_segment(aes(x, y, xend = xend, yend = yend), segments,
      colour = "grey60", linewidth = 0.3
    ) +
    geom_point(aes(colour = logFC, size = routes)) +
    ggrepel::geom_text_repel(aes(label = label), labels,
      size = labels$size, colour = labels$colour, lineheight = 0.9, seed = SEED,
      max.overlaps = Inf, box.padding = 0.15, point.padding = 0.1,
      min.segment.length = 0.3, segment.size = 0.2, bg.color = "white", bg.r = 0.15
    ) +
    scale_colour_gradient2(
      low = DIRECTION_COLOURS[["Down"]], high = DIRECTION_COLOURS[["Up"]],
      limits = c(-limit, limit), name = "log2 FC", guide = colour_bar(order = 1)
    ) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    # One size key for both sexes: the limits span every candidate passing community_routes.
    scale_size_continuous(
      range = c(1.6, 3.6), limits = c(limits$community_routes, max(cand$candidates$routes)),
      breaks = seq(limits$community_routes, max(cand$candidates$routes)), name = "routes",
      guide = guide_legend(order = 2)
    ) +
    # Room around the network for the community names.
    scale_x_continuous(expand = expansion(mult = 0.1)) +
    scale_y_continuous(expand = expansion(mult = 0.1)) +
    coord_cartesian(clip = "off") +
    labs(
      x = NULL, y = NULL, title = title,
      subtitle = sprintf(
        "%d proteins passing %d or more routes; edges: STRING score ≥ %d",
        nrow(shown), limits$community_routes, limits$string_score
      ),
      caption = if (length(alone)) {
        sprintf("%d proteins without a STRING partner (listed in the workbook)", length(alone))
      }
    ) +
    theme_figure() +
    theme(
      axis.text = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), legend.box = "vertical",
      plot.caption = element_text(size = 5.2, hjust = 0, colour = "grey30")
    )
  list(
    plot = plot,
    data = shown |>
      mutate(has_partner = gene %in% linked$gene) |>
      select(gene, logFC, routes, statistics:pathway, community, name, has_partner) |>
      rename_with(\(x) paste0("route_", x), statistics:pathway)
  )
}

# Every mouse for one sex's candidates passing CANDIDATE_DOT_ROUTES or more routes, most routes
# first and smallest p within: each protein's log2 abundance minus its own sex's vehicle mean, all
# four groups, open circles where the protein was not detected in that mouse. cands comes from
# 06_Candidates, fit from 01_Fit and quant from 02_Quantification. Returns the plot and its points:
# one row per protein and mouse, the value drawn and whether the protein was detected there.
candidate_dots <- function(cands, fit, quant, sex, title) {
  shown <- cands |>
    filter(sex == !!sex, routes >= CANDIDATE_DOT_ROUTES) |>
    arrange(desc(routes), P.Value)
  targets <- mutate(as_tibble(fit$targets), group = factor(group, FIGURE_GROUPS))
  values <- fit$E_norm$full[shown$feature, targets$sample_id, drop = FALSE]
  observations <- quant$proteins$other$n.observations
  detected <- observations[shown$feature, targets$sample_id, drop = FALSE] > 0
  genes <- set_names(
    sprintf("%s (%d)", protein_label(shown$gene, shown$feature), shown$routes), shown$feature
  )
  points <- as_tibble(centre_on_vehicle(values, targets), rownames = "feature") |>
    pivot_longer(-feature, names_to = "sample_id", values_to = "value") |>
    left_join(select(targets, sample_id, group), by = "sample_id") |>
    mutate(
      seen = detected[cbind(feature, sample_id)],
      protein = factor(genes[feature], unique(genes))
    )
  plot <- ggplot(points, aes(group, value, colour = group)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
    stat_summary(
      fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.25, colour = "grey30"
    ) +
    geom_point(aes(shape = seen),
      position = position_jitter(width = 0.12, seed = SEED), size = 0.9
    ) +
    facet_wrap(~protein, ncol = 7, scales = "free_y") +
    scale_colour_manual(values = GROUP_COLOURS, labels = GROUP_LABELS, name = NULL) +
    scale_shape_manual(
      values = c(`TRUE` = 16, `FALSE` = 1),
      labels = c(`TRUE` = "detected", `FALSE` = "not detected"),
      name = NULL
    ) +
    scale_x_discrete(labels = \(g) str_wrap(GROUP_SHORT[g], 3)) +
    labs(
      x = NULL, y = "log2 vs own-sex vehicle", title = title,
      subtitle = sprintf(
        "Proteins passing %d or more routes (routes in brackets); bar: group mean",
        CANDIDATE_DOT_ROUTES
      )
    ) +
    theme_figure() +
    theme(
      axis.text.x = element_text(size = 4.5, lineheight = 0.9),
      strip.text = element_text(size = 5.5), panel.spacing = unit(1.2, "mm")
    ) +
    guides(colour = guide_legend(nrow = 1, order = 1), shape = guide_legend(order = 2))
  list(
    plot = plot,
    data = points |>
      left_join(select(shown, feature, gene, routes), by = "feature") |>
      transmute(
        feature, gene, routes, sample_id, group,
        log2_vs_own_sex_vehicle = value,
        detected = seen
      )
  )
}
