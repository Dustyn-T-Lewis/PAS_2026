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
  inputs <- c(
    enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"),
    gene_sets = output_path("01_Processing", "04_Gene_Sets", "gene_sets.rds")
  )
  enr <- readRDS(inputs[["enrichment"]])
  gs <- readRDS(inputs[["gene_sets"]])
  contrast <- COMPARISONS[[comparison]]$contrast
  label <- COMPARISONS[[comparison]]$label
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
  families <- gs$goslim |>
    transmute(term, family = str_split(slim, "; ")) |>
    unnest(family) |>
    filter(family != "")
  family_summary <- results |>
    filter(collection == "GO:BP") |>
    inner_join(families, by = "term", relationship = "many-to-many") |>
    summarise(
      terms = n(), up_p = sum(p <= ALPHA & NES > 0), down_p = sum(p <= ALPHA & NES < 0),
      up_fdr = sum(padj <= ALPHA & NES > 0), down_fdr = sum(padj <= ALPHA & NES < 0),
      .by = family
    ) |>
    arrange(desc(up_p + down_p))

  # fgsea's own table of running scores: the five strongest sets each way per collection.
  gsea_table <- \(db) {
    \() {
      u <- collections[[db]]
      top <- results |>
        filter(collection == db) |>
        slice_min(p, n = 5, by = direction, with_ties = FALSE) |>
        arrange(desc(NES)) |>
        transmute(pathway = term, NES, pval = p, padj)
      short <- str_trunc(enrichVolcano::clean_label(top$pathway, width = 200), 36)
      print(fgsea::plotGseaTable(set_names(enr$sets[[u]][top$pathway], short),
        enr$ranks[[contrast]][[u]], mutate(top, pathway = short),
        gseaParam = 1
      ))
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
      "Per universe and collection: p ≤ 0.05, chance, FDR tiers.",
      paste(
        "Every tested set, by p: NES, p, BH FDR within collection, size, redundancy flag,",
        "singscore contrast (log2 FC, p, FDR, r with group), GO Slim families, leading edge."
      ),
      "Per GO Slim family: GO:BP terms tested and at p and FDR ≤ 0.05 by direction.",
      "GO Slim cellular compartments: fgsea, full proteome.",
      "Literature-informed and ROS sets, fgsea in the full proteome, BH across the list.",
      "Which PDF page shows what."
    ),
    out$workbook, inputs
  )
  tier_table
}
