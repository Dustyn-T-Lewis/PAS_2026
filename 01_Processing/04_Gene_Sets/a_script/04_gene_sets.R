# Every protein grouping the later steps test. Tested pathways: GO:BP, Reactome (disease and
# infection branch left out), Hallmark and KEGG LEGACY in the full proteome, MitoCarta in the
# mitochondrial universe, each trimmed to TEST_WINDOW detected genes, the same list for fgsea and
# singscore. Also GO Slim cellular components, CORUM complexes, STRING clusters, the literature and
# ROS sets, and a table of every detected protein with every tested group it belongs to. Hallmark
# also names clusters and modules and supplies four of the literature and ROS sets. None of this
# looks at the samples.

source(here::here("00_Setup", "Shared_functions", "helpers.R"))

out <- stage_paths("01_Processing", "04_Gene_Sets", reports = FALSE)
MIN_DETECTED <- 5
inputs <- c(
  mitocarta = input_path("Mouse.MitoCarta3.0.xls"),
  corum = input_path("corum_complete.txt"),
  STRING_FILES,
  proteins = here("01_Processing", "02_Quantification", "c_data", "proteins.rds")
)

quant <- readRDS(inputs[["proteins"]])
proteins <- quant$proteins
gene <- first_symbol(proteins$genes$Genes)
detected <- setdiff(unique(gene), c("", NA))
# enrichVolcano maps accessions to symbols itself, so the tested window counts the genes it ranks.
detected_in <- map(quant$universes, \(keep) {
  da <- enrichVolcano::as_da(
    tibble(
      protein = rownames(proteins)[keep], gene = gene[keep], logFC = 0, t = 0, P.Value = 1,
      adj.P.Val = 1, abundance = 0
    ),
    contrast = "none", species = "Mus musculus"
  )
  setdiff(unique(da$gene), c("", NA))
})

full <- enrichVolcano::load_gene_sets(c("Hallmark", "Reactome", "GO:BP"), species = "Mus musculus")
# KEGG LEGACY is human-only in MSigDB; msigdbr maps it to mouse orthologs, as it does every other
# collection here. enrichVolcano offers only KEGG MEDICUS, so this one is read directly.
kegg <- msigdbr::msigdbr(
  species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_LEGACY"
) |>
  distinct(gs_name, gene_symbol)
kegg <- split(kegg$gene_symbol, kegg$gs_name)
mitocarta <- readxl::read_excel(inputs[["mitocarta"]], sheet = 4)
mito <- list(MitoCarta = set_names(
  str_split(mitocarta$Genes, ",\\s*"),
  paste0("MITOCARTA_", set_id(mitocarta$MitoPathway))
))

# The tested lists: a set enters when its detected genes in the universe fall inside TEST_WINDOW,
# the window enrichVolcano applies to the ranked genes, so fgsea and singscore test the same sets.
in_window <- function(sets, universe) {
  n <- map_int(sets, \(s) length(intersect(s, detected_in[[universe]])))
  sets[n >= TEST_WINDOW[[universe]][1] & n <= TEST_WINDOW[[universe]][2]]
}
reactome <- full$Reactome[!str_detect(names(full$Reactome), REACTOME_DISEASE)]
tested <- list(
  full = list(
    `GO:BP` = in_window(full[["GO:BP"]], "full"), Reactome = in_window(reactome, "full"),
    Hallmark = in_window(full$Hallmark, "full"), KEGG = in_window(kegg, "full")
  ),
  mito = list(MitoCarta = in_window(mito$MitoCarta, "mito"))
)

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
  mutate(set = paste0("GOSLIM_", set_id(AnnotationDbi::Term(go_id))))
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

# CORUM 5.3 (Tsitsiridis et al. 2023, NAR; release 2026-04-14), every mammalian complex. Mouse
# subunits are used as named; rat, bovine, dog and pig subunits go to their human orthologs and
# every human symbol to its mouse ortholog (babelgene, HCOP); species babelgene lacks and the
# generic "Mammalia" entries carry human-style symbols and are read as human. A complex is kept
# when at least five subunits are detected here; identical detected member sets are kept once.
BABELGENE_SPECIES <- c(
  Rat = "Rattus norvegicus", Bovine = "Bos taurus", Dog = "Canis lupus familiaris",
  Pig = "Sus scrofa"
)
corum <- readr::read_tsv(inputs[["corum"]], show_col_types = FALSE) |>
  transmute(complex_id, complex_name, organism, gene = str_split(subunits_gene_name, ";")) |>
  unnest(gene) |>
  filter(!is.na(gene), gene != "")
to_human <- imap(BABELGENE_SPECIES, \(species, organism) {
  genes <- unique(corum$gene[corum$organism == organism])
  babelgene::orthologs(genes, species = species, human = FALSE) |>
    distinct(organism = organism, gene = symbol, human = human_symbol)
}) |>
  list_rbind()
corum <- corum |>
  left_join(to_human, by = c("organism", "gene"), relationship = "many-to-many") |>
  mutate(human = if_else(organism == "Mouse", NA_character_, coalesce(human, str_to_upper(gene))))
to_mouse <- babelgene::orthologs(
  unique(na.omit(corum$human)),
  species = "Mus musculus", human = TRUE
) |>
  distinct(human = human_symbol, mouse = symbol)
complex_members <- corum |>
  left_join(to_mouse, by = "human", relationship = "many-to-many") |>
  mutate(mouse = if_else(organism == "Mouse", gene, mouse)) |>
  filter(mouse %in% detected) |>
  distinct(complex_id, complex_name, organism, mouse) |>
  summarise(members = list(sort(unique(mouse))), .by = c(complex_id, complex_name, organism)) |>
  filter(lengths(members) >= MIN_DETECTED) |>
  arrange(complex_id) |>
  distinct(members, .keep_all = TRUE) |>
  mutate(set = paste0(
    "CORUM_",
    set_id(if_else(
      is.na(complex_name) | complex_name == "", paste("complex", complex_id), complex_name
    )),
    "_", complex_id
  ))
complexes <- list(CORUM = set_names(complex_members$members, complex_members$set))

# STRING v12 network among the detected proteins (score >= 400), clustered by walktrap; clusters
# of at least five proteins, named by their best over-represented Hallmark or Reactome term.
stringdb <- string_db()
string_map <- stringdb$map(data.frame(gene = detected), "gene", removeUnmappedRows = TRUE)
string_clusters <- stringdb$get_clusters(unique(string_map$STRING_id), algorithm = "walktrap") |>
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
# Named on 2026-10-05 from the PAS literature, after the first contrast results had been seen,
# so they are literature-informed rather than a priori; every one is reported whatever it shows.
# Burke et al. 2026 (Nat Commun, PMID 42431880): PAS preserves muscle "by enhancing cellular
# energy status and translational capacity". Wang et al. 2019 (EMBO Rep, PMID 31318145):
# succinate shifts fibres towards oxidative metabolism through SUCNR1 and lowers glycolysis.
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
all_sets <- unlist(unname(c(full, mito)), recursive = FALSE)
literature_sets <- all_sets[intersect(LITERATURE_SETS$set, names(all_sets))]

collections <- c(full, list(KEGG = kegg), mito, compartments, complexes, clusters)
tested_names <- unlist(map(tested, \(u) map(u, names)))
catalog <- collections |>
  imap(\(sets, db) {
    tibble(
      collection = db, set = names(sets), members = lengths(sets),
      detected = map_int(sets, \(s) sum(s %in% detected)),
      tested = set %in% c(tested_names, names(literature_sets)) |
        db %in% c("Compartments", "CORUM", "STRING")
    )
  }) |>
  list_rbind() |>
  left_join(select(goslim, set = term, go_id, slim), by = "set")

# Every detected protein and each tested group it belongs to; a protein inherits its first gene
# symbol's memberships. Long for filtering, wide for reading.
groups <- c(
  tested$full, tested$mito, compartments, complexes, clusters
)
protein_genes <- tibble(
  protein = rownames(proteins), gene = gene, is_mito = proteins$genes$is_mito
) |>
  filter(!is.na(gene), gene != "")
membership <- groups |>
  imap(\(sets, db) {
    tibble(
      collection = db, term = rep(names(sets), lengths(sets)),
      gene = unlist(sets, use.names = FALSE)
    )
  }) |>
  list_rbind() |>
  distinct() |>
  inner_join(protein_genes, by = "gene", relationship = "many-to-many") |>
  select(protein, gene, is_mito, collection, term) |>
  arrange(gene, protein, collection, term)
membership_wide <- membership |>
  summarise(
    n = n(), terms = paste(term, collapse = "; "),
    .by = c(protein, gene, is_mito, collection)
  ) |>
  pivot_wider(names_from = collection, values_from = c(n, terms), values_fill = list(n = 0L)) |>
  right_join(protein_genes, by = c("protein", "gene", "is_mito")) |>
  mutate(across(starts_with("n_"), \(x) coalesce(x, 0L))) |>
  select(
    protein, gene, is_mito,
    map(names(groups), \(db) paste0(c("n_", "terms_"), db)) |> unlist()
  ) |>
  arrange(gene, protein)
versions <- c(
  map_chr(attr(full, "versions"), as.character),
  CORUM = "5.3, released 2026-04-14, all mammalian complexes (md5 in input_manifest)",
  STRING = "12.0, score >= 400, walktrap"
)

saveRDS(
  list(
    full = full, mito = mito, tested = tested, compartments = compartments,
    complexes = complexes, clusters = clusters, membership = membership, goslim = goslim,
    slim_terms = goslim$term[goslim$go_id %in% slim_bp], ros = ROS_SETS,
    literature = LITERATURE_SETS, literature_sets = literature_sets, versions = versions
  ),
  file.path(out$data, "gene_sets.rds")
)
write_workbook(
  list(
    catalog = catalog,
    membership = membership_wide,
    membership_long = membership,
    corum = transmute(complex_members, set, complex_name, organism,
      detected = lengths(members), members = map_chr(members, \(g) paste(g, collapse = ";"))
    ),
    string_clusters = mutate(cluster_labels, set = names(clusters$STRING), .before = 1),
    literature_sets = mutate(LITERATURE_SETS, present = set %in% names(literature_sets)),
    versions = enframe(versions, "source", "version")
  ),
  c(
    "Every set: collection, members, detected, whether tested, GO ID and GO Slim ancestors.",
    "Every detected protein: count and names of the tested groups it belongs to, per collection.",
    "The same, one row per protein and group.",
    "CORUM complexes with at least five detected subunits, and those subunits.",
    "STRING clusters: size, best over-represented term (FDR <= 0.05) and members.",
    paste(
      "Literature-informed sets (chosen after the first contrast results, all reported),",
      "hypothesis, expectation."
    ),
    "Gene-set sources and versions."
  ),
  out$workbook, inputs
)

count(catalog, collection, tested)

sessionInfo()
