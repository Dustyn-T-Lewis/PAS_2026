#!/usr/bin/env Rscript
# Every protein grouping the arcs test, fixed before any result is read: Hallmark, Reactome, KEGG
# MEDICUS and GO:BP through enrichVolcano; MitoCarta pathways for the mitochondrial universe; GO
# Slim cellular components; CORUM complexes; and clusters of the STRING network among the detected
# proteins. None of them looks at the samples. Tests nothing.

source(here::here("R", "helpers.R"))

out <- stage_paths("00_Trunk", "04_Gene_Sets")
CACHE <- here("00_Trunk", "04_Gene_Sets", "cache")
CORUM_FILE <- file.path(CACHE, "corum_complete.txt")
MIN_DETECTED <- 5
inputs <- c(
  mitocarta = input_path("Mouse.MitoCarta3.0.xls"),
  proteins = here("00_Trunk", "02_Quantification", "c_data", "proteins.rds")
)
dir.create(CACHE, showWarnings = FALSE)

proteins <- readRDS(inputs[["proteins"]])$proteins
detected <- setdiff(unique(str_split_i(proteins$genes$Genes, ";", 1)), c("", NA))

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

# Each compartment holds the genes annotated to it or any descendant.
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

# CORUM 5 (Tsitsiridis et al. 2023, NAR): human and mouse complexes, human subunits mapped to
# mouse orthologs, kept when at least five subunits are detected. Identical detected member sets
# are kept once.
if (!file.exists(CORUM_FILE)) {
  download.file(paste0(
    "https://mips.helmholtz-muenchen.de/fastapi-corum/public/file/download_current_file",
    "?file_id=complete&file_format=txt"
  ), CORUM_FILE, quiet = TRUE)
}
inputs <- c(inputs, corum = CORUM_FILE)
corum <- readr::read_tsv(CORUM_FILE, show_col_types = FALSE) |>
  filter(organism %in% c("Human", "Mouse")) |>
  transmute(complex_id, complex_name, organism, gene = str_split(subunits_gene_name, ";")) |>
  unnest(gene) |>
  filter(gene != "")
orthologs <- babelgene::orthologs(unique(corum$gene[corum$organism == "Human"]),
  species = "Mus musculus", human = TRUE
) |>
  distinct(human_symbol, mouse = symbol)
complex_members <- corum |>
  left_join(orthologs, by = c(gene = "human_symbol"), relationship = "many-to-many") |>
  mutate(mouse = if_else(organism == "Mouse", gene, mouse)) |>
  filter(mouse %in% detected) |>
  distinct(complex_id, complex_name, mouse) |>
  summarise(members = list(sort(mouse)), .by = c(complex_id, complex_name)) |>
  filter(lengths(members) >= MIN_DETECTED) |>
  distinct(members, .keep_all = TRUE) |>
  mutate(set = paste0(
    "CORUM_", str_to_upper(str_replace_all(complex_name, "[^A-Za-z0-9]+", "_")), "_", complex_id
  ))
complexes <- list(CORUM = set_names(complex_members$members, complex_members$set))

# STRING v12 network among the detected proteins (score >= 400), clustered by walktrap; clusters
# of at least five proteins, named by their best over-represented Hallmark or Reactome term.
string_db <- STRINGdb::STRINGdb$new(
  version = "12.0", species = 10090, score_threshold = 400, network_type = "full",
  input_directory = CACHE
)
string_map <- string_db$map(data.frame(gene = detected), "gene", removeUnmappedRows = TRUE)
string_clusters <- string_db$get_clusters(unique(string_map$STRING_id), algorithm = "walktrap") |>
  map(\(ids) sort(unique(string_map$gene[string_map$STRING_id %in% ids]))) |>
  keep(\(g) length(g) >= MIN_DETECTED)
cluster_labels <- label_groups(string_clusters, full[c("Hallmark", "Reactome")], detected)
clusters <- list(STRING = set_names(
  string_clusters,
  sprintf("STRING_%02d_%s", seq_along(string_clusters), cluster_labels$label_id)
))

# Named before any result is read, and reported whatever they show.
ROS_SETS <- c(
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "GOBP_RESPONSE_TO_OXIDATIVE_STRESS",
  "GOBP_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS", "GOBP_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS",
  "GOBP_CELLULAR_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES",
  "REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES", "MITOCARTA_ROS_AND_GLUTATHIONE_METABOLISM"
)
# Named on 2026-10-05 from the PAS literature, after the v1 contrasts had been seen, so they are
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

collections <- c(full, mito, compartments, complexes, clusters)
catalog <- collections |>
  imap(\(sets, db) {
    tibble(
      collection = db, set = names(sets), members = lengths(sets),
      detected = map_int(sets, \(s) sum(s %in% detected))
    )
  }) |>
  list_rbind() |>
  left_join(select(goslim, set = term, go_id, slim), by = "set")
all_sets <- unlist(unname(collections), recursive = FALSE)
versions <- c(
  map_chr(attr(full, "versions"), as.character),
  CORUM = "5.3, released 2026-04-14 (md5 in input_manifest)",
  STRING = "12.0, score >= 400, walktrap"
)

saveRDS(
  list(
    full = full, mito = mito, compartments = compartments, complexes = complexes,
    clusters = clusters, goslim = goslim, slim_terms = goslim$term[goslim$go_id %in% slim_bp],
    ros = ROS_SETS, literature = LITERATURE_SETS, versions = versions
  ),
  file.path(out$data, "gene_sets.rds")
)
write_workbook(
  list(
    catalog = catalog,
    corum = transmute(complex_members, set, complex_name,
      detected = lengths(members), members = map_chr(members, \(g) paste(g, collapse = ";"))
    ),
    string_clusters = mutate(cluster_labels, set = names(clusters$STRING), .before = 1),
    literature_sets = mutate(LITERATURE_SETS, present = set %in% names(all_sets)),
    versions = enframe(versions, "source", "version")
  ),
  c(
    "Every set: collection, members, members detected here, GO ID and GO Slim ancestors.",
    "CORUM complexes with at least five detected subunits, and those subunits.",
    "STRING clusters: size, best over-represented term (FDR <= 0.05) and members.",
    "Literature-informed sets (chosen after the v1 results, all reported), hypothesis, expectation.",
    "Gene-set sources and versions."
  ),
  out$workbook, inputs
)

count(catalog, collection)

sessionInfo()
