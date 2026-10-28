# Supplementary Figure 2, proteins per contrast: each comparison's Π <= 0.05 proteins across all
# mice beside the sets each block over-represents, and every protein's fold change and p. Each
# panel is its own script in panels/; this file runs them, stitches them under one key and writes
# one workbook.
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

PANELS <- c("a_pi_female", "b_pi_male", "c_pi_interaction", "d_volcanoes")
panels <- run_panels("S2_Proteins", PANELS)

heatmaps <- panels[1:3]
key <- pi_heatmap_key(unique(list_rbind(map(heatmaps, c("data", "ora")))$collection))
figure <- wrap_elements(full = wrap_plots(map(heatmaps, "plot"), nrow = 1)) /
  wrap_elements(full = key) / wrap_elements(full = panels$d_volcanoes$plot) +
  plot_layout(heights = c(126, 9, 62))
out <- stage_paths("04_Figures", "S2_Proteins")
save_figure(figure, out$reports, "S2_composite", width_mm = 180, height_mm = 197)

ORA_HOLDS <- paste(
  "The panel's ORA bars: per block, the three sets with the smallest p among those shown (one per",
  "overlap group at p ≤ 0.05), from 02_Enrichment's block_ora."
)
write_figure_workbook(
  panels,
  c(
    "Panel a: female PAS vs vehicle proteins at Π ≤ 0.05: log2 FC, p, FDR, Π.",
    ORA_HOLDS,
    "Panel b: male PAS vs vehicle proteins at Π ≤ 0.05: log2 FC, p, FDR, Π.",
    ORA_HOLDS,
    "Panel c: treatment × sex proteins at Π ≤ 0.05: log2 FC, p, FDR, Π.",
    ORA_HOLDS,
    "Panel d: every full-proteome protein per contrast: log2 FC, p, FDR, Π, colour class."
  ),
  file.path(out$data, "S2.xlsx")
)
