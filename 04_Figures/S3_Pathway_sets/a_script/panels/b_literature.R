# S3 b: the literature-informed and ROS sets, grouped by the hypothesis that named them, with the
# direction each hypothesis expected; fgsea per contrast, BH across the list (02_Enrichment).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"))
literature <- readRDS(inputs[["enrichment"]])$literature
cells <- literature |>
  transmute(
    contrast, NES,
    fdr = list_fdr, size, expected,
    row = factor(paste(hypothesis, set), rev(unique(paste(hypothesis, set)))),
    label = sprintf("%s [%s]", set_label(set, 38), expected),
    row_group = str_wrap(hypothesis, 12)
  )

plot <- set_dot_matrix(cells,
  title = "Literature-informed and ROS sets",
  subtitle = "Grouped by hypothesis; [expected direction]"
)

save_panel(plot, "S3_Pathway_sets", "b_literature", 100, 120)
list(plot = plot, data = select(cells, -row), inputs = inputs)
