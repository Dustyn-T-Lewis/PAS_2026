#!/usr/bin/env Rscript
# The gene-set collections every pathway step tests: Hallmark, Reactome, KEGG MEDICUS and full
# GO:BP through enrichVolcano for the full proteome, MitoCarta pathways for the mitochondrial
# universe, and the GO Slim cellular components for the compartment test. Each GO:BP term carries
# its GO Slim ancestors. Tests nothing.

source(here::here("R", "helpers.R"))

out <- stage_paths("03_Pathway_Enrichment", "00_Gene_Sets")
inputs <- c(mitocarta = here("00_Input", "Mouse.MitoCarta3.0.xls"))

full <- enrichVolcano::load_gene_sets(c("Hallmark", "Reactome", "KEGG", "GO:BP"),
  species = "Mus musculus"
)
mitocarta <- readxl::read_excel(inputs[["mitocarta"]], sheet = 4)
mito <- list(MitoCarta = set_names(
  str_split(mitocarta$Genes, ",\\s*"),
  paste0("MITOCARTA_", str_to_upper(str_replace_all(mitocarta$MitoPathway, "[^A-Za-z0-9]+", "_")))
))

# GO Slim from the obo file enrichVolcano pins, so both use one definition of the slim.
obo <- readLines(system.file("extdata", "goslim_generic.obo.gz", package = "enrichVolcano"))
terms <- split(obo, cumsum(obo == "[Term]"))
slim_ids <- \(namespace) {
  terms |>
    keep(\(x) any(x == paste("namespace:", namespace))) |>
    map_chr(\(x) str_remove(x[str_starts(x, "id: ")][1], "^id: "))
}
slim_bp <- slim_ids("biological_process")
slim_name <- AnnotationDbi::Term(slim_bp)

# Each compartment holds the genes annotated to it or any descendant, as enrichVolcano builds GO
# Slim BP sets.
cc <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db,
  keys = slim_ids("cellular_component"), columns = "SYMBOL", keytype = "GOALL"
)) |>
  distinct(go_id = GOALL, gene = SYMBOL) |>
  filter(!is.na(gene)) |>
  mutate(set = paste0("GOSLIM_", str_to_upper(str_replace_all(
    AnnotationDbi::Term(go_id), "[^A-Za-z0-9]+", "_"
  ))))
compartments <- list(Compartments = split(cc$gene, cc$set))

go_ids <- msigdbr::msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP") |>
  distinct(term = gs_name, go_id = gs_exact_source)
ancestors <- as.list(GO.db::GOBPANCESTOR)
goslim <- go_ids |>
  filter(term %in% names(full[["GO:BP"]])) |>
  mutate(slim = map_chr(go_id, \(id) {
    hit <- intersect(c(id, ancestors[[id]]), slim_bp)
    paste(sort(slim_name[hit]), collapse = "; ")
  }))

# Named before any result is read, and reported whatever they show.
ROS_SETS <- c(
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "GOBP_RESPONSE_TO_OXIDATIVE_STRESS",
  "GOBP_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS", "GOBP_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS",
  "GOBP_CELLULAR_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES",
  "REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES", "MITOCARTA_ROS_AND_GLUTATHIONE_METABOLISM"
)
# Named on 2026-10-05 from the PAS literature, after the contrasts had been seen, so they are
# literature-informed rather than a priori; every one is reported whatever it shows. Burke et al.
# 2026 (Nat Commun, PMID 42431880): PAS preserves muscle "by enhancing cellular energy status and
# translational capacity". Wang et al. 2019 (EMBO Rep, PMID 31318145): succinate shifts fibres
# towards oxidative metabolism through SUCNR1 and lowers glycolysis.
LITERATURE_SETS <- tribble(
  ~set, ~hypothesis, ~expected,
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "Energy status", "up",
  "MITOCARTA_TCA_CYCLE", "Energy status", "up",
  "MITOCARTA_FATTY_ACID_OXIDATION", "Energy status", "up",
  "GOBP_CYTOPLASMIC_TRANSLATION", "Translational capacity", "up",
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION", "Translational capacity", "up",
  "HALLMARK_MTORC1_SIGNALING", "Translational capacity", "up",
  "MITOCARTA_MITOCHONDRIAL_CENTRAL_DOGMA", "Translational capacity", "up",
  "HALLMARK_GLYCOLYSIS", "Oxidative fibre shift", "down",
  "REACTOME_STRIATED_MUSCLE_CONTRACTION", "Oxidative fibre shift", "up",
  "GOBP_PROTEASOME_MEDIATED_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS",
  "Atrophy defence", "down",
  "REACTOME_AUTOPHAGY", "Atrophy defence", "down",
  "MITOCARTA_LYSINE_METABOLISM", "Pipecolate handling", "unclear"
) |>
  bind_rows(tibble(set = ROS_SETS, hypothesis = "ROS (named in advance)", expected = "unclear"))
all_sets <- c(unlist(unname(full), recursive = FALSE), mito$MitoCarta)
ros <- tibble(
  set = ROS_SETS, present = ROS_SETS %in% names(all_sets),
  members = lengths(all_sets[ROS_SETS])
)

catalog <- c(full, mito, compartments) |>
  imap(\(sets, db) tibble(collection = db, set = names(sets), members = lengths(sets))) |>
  list_rbind() |>
  left_join(select(goslim, set = term, go_id, slim), by = "set")

saveRDS(
  list(
    full = full, mito = mito, compartments = compartments, goslim = goslim,
    slim_terms = goslim$term[goslim$go_id %in% slim_bp], ros = ROS_SETS,
    literature = LITERATURE_SETS,
    versions = attr(full, "versions")
  ),
  file.path(out$data, "gene_sets.rds")
)
write_workbook(
  list(
    catalog = catalog, ros_sets = ros,
    literature_sets = mutate(LITERATURE_SETS, present = set %in% names(all_sets)),
    versions = enframe(map_chr(attr(full, "versions"), as.character), "source", "version")
  ),
  c(
    "Every set: collection, member count before trimming to the data, GO ID and GO Slim ancestors.",
    "The ROS sets named in advance, whether each exists and its size.",
    "Literature-informed sets (chosen after results, all reported), hypothesis, expectation.",
    "Gene-set sources and versions."
  ),
  out$workbook, inputs
)

count(catalog, collection)
ros

sessionInfo()
