# Step 02 of every comparison: which pathways its contrast shifted.
# Run by 03_Comparisons/<comparison>/a_script/run_comparison.R.

source(here::here("00_Setup", "Shared_functions", "comparisons.R"))

# Pathways: the comparison's contrast sliced from 02_Enrichment. fgsea on GO:BP, Reactome, Hallmark
# and KEGG in the full proteome and MitoCarta in the mitochondrial universe, BH within collection
# over every tested set, each beside the same set's singscore contrast; GO Slim compartments and
# the literature list alongside. GO Slim generic BP families (every slim ancestor a term reaches)
# summarise the GO:BP terms; they never choose what is tested.
report_pathways <- function(comparison) {
  out <- comparison_paths(comparison, "02_pathways")
  inputs <- c(enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"))
  enr <- readRDS(inputs[["enrichment"]])
  contrast <- COMPARISONS[[comparison]]$contrast
  collections <- deframe(distinct(enr$results, collection, universe))

  results <- filter(enr$results, contrast == !!contrast) |>
    select(-contrast) |>
    arrange(factor(collection, names(collections)), p)
  tier_table <- filter(enr$summary, contrast == !!contrast) |>
    select(-contrast)
  compartments <- filter(enr$compartments, contrast == !!contrast) |>
    select(-contrast)
  literature <- filter(enr$literature, contrast == !!contrast) |>
    select(-contrast)
  family_summary <- results |>
    filter(collection == "GO:BP", !is.na(families)) |>
    separate_longer_delim(families, "; ") |>
    rename(family = families) |>
    summarise(
      terms = n(), up_p = sum(p <= ALPHA & NES > 0), down_p = sum(p <= ALPHA & NES < 0),
      up_fdr = sum(padj <= ALPHA & NES > 0), down_fdr = sum(padj <= ALPHA & NES < 0),
      .by = family
    ) |>
    arrange(desc(up_p + down_p))

  gsea_table <- \(db) {
    u <- collections[[db]]
    \() {
      gsea_table_page(filter(results, collection == db), enr$sets[[u]], enr$ranks[[contrast]][[u]])
    }
  }
  plot_index <- preview_pdf(
    file.path(out$reports, "02_pathways.pdf"),
    set_names(map(names(collections), gsea_table), paste("strongest sets,", names(collections)))
  )

  write_workbook(
    list(
      summary = tier_table, fgsea = results, families = family_summary, compartments = compartments,
      literature = literature, plot_index = plot_index
    ),
    c(
      "Per universe and collection: p ≤ 0.05, FDR tiers.",
      paste(
        "Every tested set, by p: NES, p, BH FDR within collection, size, redundancy flag, camera",
        "p and FDR, singscore contrast (log2 FC, p, FDR, r with group), GO Slim families,",
        "leading edge."
      ),
      "Per GO Slim family: GO:BP terms tested and at p and FDR ≤ 0.05 by direction.",
      "GO Slim cellular compartments: fgsea, full proteome.",
      "Literature-informed and ROS sets, fgsea in the full proteome, BH across the list.",
      "Which PDF page shows what."
    ),
    out$workbook
  )
  tier_table
}
