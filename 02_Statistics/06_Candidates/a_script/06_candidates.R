# Candidate proteins for follow-up: six routes, each answering a different worry (statistics,
# every mouse, both sexes, phenotype, module hub, pathway leading edge), and how many each protein
# passes. Discovery-stage candidates, not validated biomarkers. A protein counts in a sex only when
# detected in MIN_DETECTED_PER_GROUP or more mice of each group there, so its values are mostly
# measured. Candidates passing COMMUNITY_ROUTES or more are grouped by STRING; a community's name
# is its smallest-p ORA term, a label rather than a test.

source(here::here("00_Setup", "Shared_functions", "helpers.R"))

MIN_DETECTED_PER_GROUP <- 4
ROUTE_P <- 0.01
BOTH_SEXES_P <- 0.10
# Five against five, AUC moves in steps of 1/25: 0.92 is at most two of the 25 pairs out of order.
AUC_CUTS <- c(0.08, 0.92)
HUB_KME <- 0.8
ORA_MIN_SIZE <- 5
STRING_SCORE <- 700
COMMUNITY_ROUTES <- 2

out <- stage_paths("02_Statistics", "06_Candidates")
inputs <- c(
  fit = output_path("02_Statistics", "01_Fit", "fit.rds"),
  proteins = output_path("01_Processing", "02_Quantification", "proteins.rds"),
  design = output_path("01_Processing", "03_Design", "design.rds"),
  enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"),
  modules = output_path("01_Processing", "06_Modules", "modules.rds"),
  gene_sets = output_path("01_Processing", "04_Gene_Sets", "gene_sets.rds")
)
fit <- readRDS(inputs[["fit"]])
quant <- readRDS(inputs[["proteins"]])
d <- readRDS(inputs[["design"]])
enr <- readRDS(inputs[["enrichment"]])
mods <- readRDS(inputs[["modules"]])
gs <- readRDS(inputs[["gene_sets"]])
targets <- as_tibble(fit$targets)
expr <- fit$E_norm$full[, targets$sample_id]
detected <- quant$proteins$other$n.observations[rownames(expr), targets$sample_id] > 0
sv1 <- d$design[targets$sample_id, "sv1"]
SEX_CONTRAST <- map_chr(SEX_COMPARISONS, "contrast")

proteins <- filter(fit$protein_contrasts, universe == "full")
both <- proteins |>
  filter(contrast %in% SEX_CONTRAST) |>
  select(feature, contrast, logFC, P.Value) |>
  pivot_wider(names_from = contrast, values_from = c(logFC, P.Value)) |>
  filter(
    P.Value_PAS_vs_VEH_female <= BOTH_SEXES_P, P.Value_PAS_vs_VEH_male <= BOTH_SEXES_P,
    sign(logFC_PAS_vs_VEH_female) == sign(logFC_PAS_vs_VEH_male)
  ) |>
  pull(feature)

routes <- imap(SEX_CONTRAST, \(cn, sex) {
  code <- SEX_COMPARISONS[[sex]]$sex
  veh <- targets$group == paste0("VEH_", code)
  pas <- targets$group == paste0("PAS_", code)
  measured <- rowSums(detected[, veh]) >= MIN_DETECTED_PER_GROUP &
    rowSums(detected[, pas]) >= MIN_DETECTED_PER_GROUP
  in_sex <- veh | pas
  # The surrogate variable is removed within the sex before mice are compared one by one.
  adjusted <- limma::removeBatchEffect(expr[, in_sex], covariates = sv1[in_sex])
  case <- factor(pas[in_sex], c(FALSE, TRUE))
  # pROC's trapezoid sum lands a hair off the 1/25 steps, so it is rounded before the cuts.
  auc <- apply(adjusted, 1, \(v) {
    round(as.numeric(pROC::auc(case, v, levels = levels(case), direction = "<", quiet = TRUE)), 4)
  })
  associated <- fit$protein_associations |>
    filter(universe == "full", model == "treatment", sex == !!sex, P.Value <= ROUTE_P) |>
    pull(feature)
  responsive <- fit$score_contrasts |>
    filter(level == "module", contrast == cn, P.Value <= ALPHA) |>
    pull(feature)
  hubs <- filter(mods$modules, module %in% responsive, kME >= HUB_KME)$protein
  leading <- enr$results |>
    filter(contrast == cn, padj <= ALPHA, dedup_status == "kept") |>
    pull(leading_edge) |>
    str_split(";") |>
    unlist()
  proteins |>
    filter(contrast == cn, measured[feature]) |>
    transmute(
      sex, feature, gene, logFC, logFC_shrunk, P.Value, adj.P.Val, lfsr, pi_score,
      auc = auc[feature],
      statistics = P.Value <= ROUTE_P,
      every_mouse = auc >= AUC_CUTS[2] | auc <= AUC_CUTS[1],
      both_sexes = feature %in% both,
      phenotype = feature %in% associated,
      hub = feature %in% hubs,
      pathway = gene %in% leading & P.Value <= ALPHA
    ) |>
    mutate(
      routes = statistics + every_mouse + both_sexes + phenotype + hub + pathway,
      detected_veh = rowSums(detected[feature, veh, drop = FALSE]),
      detected_pas = rowSums(detected[feature, pas, drop = FALSE])
    )
}) |>
  list_rbind() |>
  filter(routes >= 1) |>
  arrange(sex, desc(routes), P.Value)

string <- string_db()
string_ids <- string$map(data.frame(gene = unique(routes$gene)), "gene", removeUnmappedRows = TRUE)
links <- string$get_interactions(string_ids$STRING_id) |>
  filter(combined_score >= STRING_SCORE) |>
  transmute(
    a = string_ids$gene[match(from, string_ids$STRING_id)],
    b = string_ids$gene[match(to, string_ids$STRING_id)]
  ) |>
  filter(a != b) |>
  distinct()
detected_genes <- unique(proteins$gene[!is.na(proteins$gene) & proteins$gene != ""])
full_sets <- unlist(unname(gs$full[c("Hallmark", "Reactome", "GO:BP")]), recursive = FALSE)
name_community <- function(members) {
  ora <- fgsea::fora(full_sets, members, detected_genes, minSize = ORA_MIN_SIZE, maxSize = 500) |>
    filter(overlap >= 2) |>
    arrange(pval)
  if (nrow(ora) == 0) NA_character_ else ora$pathway[1]
}
communities <- map(names(SEX_CONTRAST), \(sex) {
  genes <- unique(routes$gene[routes$sex == sex & routes$routes >= COMMUNITY_ROUTES])
  edges <- filter(links, a %in% genes, b %in% genes)
  graph <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = genes)
  membership <- igraph::membership(igraph::cluster_walktrap(graph))
  tibble(sex = sex, gene = names(membership), community = as.integer(membership)) |>
    mutate(size = n(), .by = community) |>
    mutate(community = if_else(size >= 2, community, NA_integer_))
}) |>
  list_rbind()
community_names <- communities |>
  filter(!is.na(community)) |>
  summarise(members = list(gene), .by = c(sex, community)) |>
  mutate(name = map_chr(members, name_community)) |>
  select(sex, community, name)
candidates <- routes |>
  left_join(select(communities, sex, gene, community), by = c("sex", "gene")) |>
  left_join(community_names, by = c("sex", "community"))
network_edges <- links |>
  inner_join(distinct(communities, sex, gene), by = c(a = "gene"), relationship = "many-to-many") |>
  semi_join(communities, by = c("sex", b = "gene"))

summary <- candidates |>
  summarise(
    any_route = n(), two_or_more = sum(routes >= 2), three_or_more = sum(routes >= 3),
    across(c(statistics, every_mouse, both_sexes, phenotype, hub, pathway), sum),
    .by = sex
  )

plot_index <- preview_pdf(file.path(out$reports, "06_candidates.pdf"), list(
  "routes passed per candidate" = \() {
    par(mfrow = c(1, 2))
    walk(names(SEX_CONTRAST), \(s) {
      barplot(table(candidates$routes[candidates$sex == s]),
        xlab = "routes passed", ylab = "proteins", main = s
      )
    })
  }
))

saveRDS(
  list(
    candidates = candidates, edges = network_edges,
    thresholds = list(p = ROUTE_P, string_score = STRING_SCORE, community_routes = COMMUNITY_ROUTES)
  ),
  file.path(out$data, "candidates.rds")
)
write_workbook(
  list(summary = summary, candidates = candidates, edges = network_edges, plot_index = plot_index),
  c(
    "Per sex: proteins passing any route, two or more, three or more, and each route's count.",
    paste(
      "Every protein passing at least one route in a sex: effect, p, FDR, lfsr, Π, AUC on",
      "sv1-adjusted values, each route, routes passed, detection per group, STRING community",
      "and its name."
    ),
    str_glue(
      "STRING v12 links (combined score >= {STRING_SCORE}) among candidates passing",
      " {COMMUNITY_ROUTES} or more routes."
    ),
    "Which PDF page shows what."
  ),
  out$workbook
)

summary
