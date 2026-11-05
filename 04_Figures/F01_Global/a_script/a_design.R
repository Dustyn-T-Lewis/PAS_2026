# F01 a: the design. The image is assets/design.png (built from assets/design_spec/).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(
  design = here("04_Figures", "F01_Global", "assets", "design.png"),
  fit = output_path("02_Statistics", "01_Fit", "fit.rds")
)
targets <- readRDS(inputs[["fit"]])$targets
groups <- count(targets, group, name = "mice")
per_group <- unique(range(groups$mice))
WIDTH_MM <- 66
# The image resized to 600 dpi across the panel, so the PDF does not embed it at full size.
design <- magick::image_read(inputs[["design"]]) |>
  magick::image_resize(sprintf("%dx", ceiling(WIDTH_MM / 25.4 * 600)))

plot <- ggplot() +
  annotation_custom(grid::rasterGrob(design,
    y = grid::unit(1, "npc"), vjust = 1, interpolate = TRUE
  )) +
  labs(
    title = "Design",
    subtitle = str_glue("{nrow(targets)} mice, {str_flatten(per_group, ' to ')} per group")
  ) +
  theme_figure() +
  theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank())

save_figure(plot, "F01_Global", "a_design", WIDTH_MM, 64)
list(plot = plot, data = groups)
