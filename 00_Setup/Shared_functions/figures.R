# Figure style, after HRvLR's: colours sampled from the design panel so a contrast looks the same
# in the schematic and every plot, a minimal theme, and the plots more than one script draws.
# Sourced by helpers.R.

library(ggplot2)
library(patchwork)

SEX_COLOURS <- c(Female = "#D55E00", Male = "#0072B2")
SEX_SHAPES <- c(F = 16, M = 17)
FIGURE_GROUPS <- c("VEH_F", "PAS_F", "VEH_M", "PAS_M")
GROUP_COLOURS <- c(VEH_F = "#FDB863", PAS_F = "#D55E00", VEH_M = "#92C5DE", PAS_M = "#0072B2")
# Up and down in green and magenta (PiYG), so direction never borrows a group's orange or blue.
DIRECTION_COLOURS <- c(Down = "#C51B7D", Up = "#4D9221")
DIRECTION_LIGHT <- c(Down = "#F1B6DA", Up = "#B8E186")
PHENOTYPE_LABELS <- c(
  csa_mean = "CSA", capillary_density = "Capillary density",
  satellite_cells_per_fiber = "Satellite cells / fibre",
  conductance_sol_pyr = "Soleus pyruvate", conductance_pla_pyr = "Plantaris pyruvate",
  conductance_sol_oct = "Soleus octanoyl"
)
CONTRAST_LABELS <- c(
  PAS_vs_VEH_female = "PAS vs VEH, female", PAS_vs_VEH_male = "PAS vs VEH, male",
  treatment_by_sex = "Treatment × sex", sex_in_VEH = "Sex in VEH", sex_in_PAS = "Sex in PAS"
)
# Each contrast in three tones, as the design panel draws it: line for marks and box borders, fill
# for the box tint, ink for text on white.
CONTRAST_PALETTE <- tibble::tribble(
  ~contrast, ~line, ~fill, ~ink,
  "PAS_vs_VEH_female", "#D55E00", "#FBDCC8", "#A34400",
  "PAS_vs_VEH_male", "#0072B2", "#CCE3F0", "#005282",
  "treatment_by_sex", "#7B3294", "#E8D5EE", "#5A2370",
  "sex_in_VEH", "#737373", "#E5E5E5", "#525252",
  "sex_in_PAS", "#A6A6A6", "#F0F0F0", "#737373"
)
contrast_tone <- \(tone) set_names(CONTRAST_PALETTE[[tone]], CONTRAST_PALETTE$contrast)
DATABASE_COLOURS <- c("GO:BP" = "#999999", Reactome = "#56B4E9", MitoCarta = "#009E73")

# HRvLR's figure theme: minimal, a light major grid, bold title and a plain grey subtitle, with the
# panel letter in a left margin column level with the title.
theme_figure <- function() {
  theme_minimal(base_size = 7, base_family = "Helvetica") +
    theme(
      plot.title = element_text(face = "bold", size = 7.5, margin = margin(b = 1)),
      plot.subtitle = element_text(size = 6, colour = "grey30", margin = margin(b = 3)),
      plot.title.position = "plot",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2),
      legend.key.size = unit(3, "mm"),
      plot.tag = element_text(face = "bold", size = 9, vjust = 1, margin = margin(r = 2)),
      plot.tag.location = "margin", plot.tag.position = "left"
    )
}

# Facet strips as the design panel's boxes: the contrast's tint inside a border in its colour, the
# name in bold black; a level outside the palette gets no box.
contrast_strips <- function(levels, side = c("x", "y")) {
  side <- match.arg(side)
  boxes <- ggh4x::elem_list_rect(
    fill = unname(coalesce(contrast_tone("fill")[levels], "transparent")),
    colour = unname(coalesce(contrast_tone("line")[levels], "transparent")),
    linewidth = 0.6
  )
  text <- ggh4x::elem_list_text(colour = "black", face = "bold")
  if (side == "x") {
    ggh4x::strip_themed(background_x = boxes, text_x = text)
  } else {
    ggh4x::strip_themed(background_y = boxes, text_y = text)
  }
}

# A finished figure as PDF (cairo, which carries Π and ≤) and as a 600 dpi PNG.
save_figure <- function(figure, dir, name, width_mm, height_mm) {
  save <- \(ext, ...) {
    ggsave(file.path(dir, paste0(name, ext)), figure,
      width = width_mm, height = height_mm, units = "mm", ...
    )
  }
  save(".pdf", device = grDevices::cairo_pdf)
  save(".png", dpi = 600)
}

# A step's preview PDF: one landscape page per named function, each drawing with a package's own
# plot or base graphics; cairo, so Π and ≤ render. Returns the rows of the workbook's plot_index.
preview_pdf <- function(file, pages) {
  cairo_pdf(file, width = 7, height = 5, onefile = TRUE)
  on.exit(dev.off())
  defaults <- par(no.readonly = TRUE)
  walk(pages, \(draw) {
    draw()
    par(defaults)
  })
  tibble(file = basename(file), page = seq_along(pages), shows = names(pages))
}

# PAS against vehicle per phenotype and sex as a % change with its 95% CI, on a log ratio axis so
# rises and falls of the same size look the same size.
phenotype_effect_plot <- function(effects, title, subtitle) {
  effects |>
    mutate(
      sex = factor(if_else(sex == "F", "Female", "Male"), c("Male", "Female")),
      trait = factor(PHENOTYPE_LABELS[phenotype], rev(PHENOTYPE_LABELS)),
      across(c(pct_change, ci_low, ci_high), \(x) 1 + x / 100)
    ) |>
    ggplot(aes(pct_change, trait, colour = sex)) +
    geom_vline(xintercept = 1, colour = "grey50", linewidth = 0.3) +
    geom_pointrange(aes(xmin = ci_low, xmax = ci_high, shape = P.Value <= ALPHA),
      position = position_dodge(width = 0.6), size = 0.25, linewidth = 0.4, fill = "white"
    ) +
    scale_x_continuous(
      transform = "log2", breaks = c(0.25, 0.5, 1, 2, 4, 8),
      labels = \(x) sprintf("%+.0f%%", 100 * (x - 1))
    ) +
    scale_colour_manual(values = SEX_COLOURS, breaks = c("Female", "Male"), name = NULL) +
    scale_shape_manual(
      values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "p ≤ 0.05", `FALSE` = "p > 0.05"),
      name = NULL
    ) +
    labs(
      x = "PAS vs vehicle, % change (95% CI, log scale)", y = NULL, title = title,
      subtitle = subtitle
    ) +
    theme_figure()
}
