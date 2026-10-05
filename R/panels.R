# Figures are built from panel scripts. Each panel script runs alone: it reads saved outputs, draws
# one panel, saves it to b_reports/<part>/panels/ and leaves `panel` (the ggplot), `panel_data` (its
# source table), `panel_note` (one line for the workbook) and `inputs` behind. A stitcher runs the
# panel scripts in order, assembles the figure and writes one workbook of panel data.

source(here::here("R", "helpers.R"))
source(here::here("R", "figures.R"))

figure_dirs <- function(figure, part) {
  root <- here("06_Figures", figure)
  dirs <- list(
    scripts = file.path(root, "a_script", part, "panels"),
    reports = file.path(root, "b_reports", part),
    panels = file.path(root, "b_reports", part, "panels"),
    data = file.path(root, "c_data")
  )
  walk(dirs[-1], dir.create, recursive = TRUE, showWarnings = FALSE)
  dirs
}

save_panel <- function(plot, figure, part, name, width_mm, height_mm) {
  ggsave(file.path(figure_dirs(figure, part)$panels, paste0(name, ".pdf")), plot,
    width = width_mm, height = height_mm, units = "mm", device = grDevices::cairo_pdf
  )
  invisible(plot)
}

run_panels <- function(figure, part, names) {
  scripts <- file.path(figure_dirs(figure, part)$scripts, paste0(names, ".R"))
  map(set_names(scripts, names), \(script) {
    env <- new.env(parent = globalenv())
    sys.source(script, envir = env)
    list(plot = env$panel, data = env$panel_data, note = env$panel_note, inputs = env$inputs)
  })
}

write_panel_data <- function(panels, figure, name) {
  inputs <- unlist(unname(map(panels, "inputs")))
  write_workbook(
    map(panels, "data"), unname(map_chr(panels, "note")),
    file.path(figure_dirs(figure, "main")$data, paste0(name, "_data.xlsx")),
    inputs[!duplicated(inputs)]
  )
}
