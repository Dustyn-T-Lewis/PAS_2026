#!/usr/bin/env Rscript
# Whether the proteins PAS moved are linked in STRING, and which regulators the literature says
# would move them, using each database's own documented test: STRINGdb's PPI enrichment against
# the detected proteome, and INDRA CoGEx's discrete and signed regulator analyses. Service answers
# are cached, keyed on the lists they were asked about.

source(here::here("R", "helpers.R"))
library(httr2)

STRING_VERSION <- "12.0"
STRING_SCORE <- 400
FDR_LIST_CUT <- 0.10
MIN_EVIDENCE <- 2
FDR_CUT <- 0.05
N_TOP <- 10
# An interaction "up" means a larger PAS effect in females, not a rise in any tissue, so only the
# simple effects carry the direction signed_analysis scores.
SIGNED_CONTRASTS <- c("PAS_vs_VEH_female", "PAS_vs_VEH_male")

out <- stage_paths("04_Network", "06_Mechanism")
inputs <- c(
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)
CACHE <- file.path(out$data, "mechanism_cache.rds")
STRING_FILES <- here("04_Network", "stringdb_cache")

de <- readRDS(inputs[["contrasts"]])
detected <- tibble(protein = rownames(de$genes), symbol = str_split_i(de$genes$Genes, ";", 1))

female_fdr <- de$results |>
  filter(universe == "full", contrast == "PAS_vs_VEH_female", adj.P.Val < FDR_LIST_CUT) |>
  transmute(list = "female_FDR_0.10", contrast, direction = "both", protein)
pi_lists <- de$pi_lists |>
  filter(contrast %in% PRIMARY) |>
  transmute(list = str_glue("{contrast}_pi_{direction}"), contrast, direction, protein)
lists <- bind_rows(female_fdr, pi_lists) |> mutate(list = as.character(list))
stopifnot(nrow(female_fdr) > 0, all(lists$protein %in% detected$protein))
key <- rlang::hash(arrange(lists, list, protein))

# INDRA is built on human genes; babelgene's orthology table maps the mouse symbols.
orthologs <- babelgene::orthologs(setdiff(unique(detected$symbol), c("", NA)),
  species = "mouse", human = FALSE
) |>
  distinct(symbol, human_symbol)
human_of <- function(p) {
  unique(orthologs$human_symbol[orthologs$symbol %in% detected$symbol[detected$protein %in% p]])
}
human_background <- unique(orthologs$human_symbol)

indra <- function(endpoint, body) {
  request("https://discovery.indra.bio/api") |>
    req_url_path_append(endpoint) |>
    req_body_json(body) |>
    req_timeout(900) |>
    req_perform() |>
    resp_body_string() |>
    # CoGEx writes undefined p-values as a bare NaN, which JSON does not allow.
    str_replace_all(fixed(":NaN"), ":null") |>
    jsonlite::fromJSON()
}

if (!file.exists(CACHE) || !identical(readRDS(CACHE)$key, key)) {
  dir.create(STRING_FILES, showWarnings = FALSE)
  string_db <- STRINGdb::STRINGdb$new(
    version = STRING_VERSION, species = 10090, score_threshold = STRING_SCORE,
    network_type = "full", input_directory = STRING_FILES
  )
  string_map <- string_db$map(as.data.frame(detected), "protein", removeUnmappedRows = TRUE)
  string_db$set_background(unique(string_map$STRING_id))

  string <- lists |>
    summarise(proteins = list(protein), .by = list) |>
    mutate(result = map(proteins, \(p) {
      ids <- unique(string_map$STRING_id[string_map$protein %in% p])
      e <- string_db$get_ppi_enrichment(ids)
      tibble(n_mapped = length(ids), edges = e$edges, expected = e$lambda, p_value = e$enrichment)
    })) |>
    select(-proteins) |>
    tidyr::unnest(result)

  discrete <- indra("discrete_analysis", list(
    gene_list = human_of(female_fdr$protein), background_gene_list = human_background,
    method = "fdr_bh", alpha = FDR_CUT, keep_insignificant = TRUE,
    minimum_evidence_count = MIN_EVIDENCE, minimum_belief = 0, indra_path_analysis = TRUE
  ))

  signed <- map(SIGNED_CONTRASTS, \(cn) {
    genes <- \(dir) human_of(pi_lists$protein[pi_lists$contrast == cn & pi_lists$direction == dir])
    indra("signed_analysis", list(
      positive_genes = genes("up"), negative_genes = genes("down"), alpha = FDR_CUT,
      keep_insignificant = TRUE, minimum_evidence_count = MIN_EVIDENCE, minimum_belief = 0
    )) |>
      as_tibble() |>
      transmute(
        contrast = cn, regulator = name, curie, correct, incorrect, ambiguous,
        p_value = binom_pvalue
      )
  }) |>
    list_rbind()

  saveRDS(
    list(
      string = string, string_mapped = n_distinct(string_map$protein), discrete = discrete,
      signed = signed, fetched = Sys.Date(), key = key
    ),
    CACHE
  )
}

cache <- readRDS(CACHE)

string_results <- cache$string |>
  left_join(count(lists, list, contrast, direction, name = "n_proteins"), by = "list") |>
  mutate(fold = edges / expected, FDR = p.adjust(p_value, "BH")) |>
  relocate(list, contrast, direction, n_proteins)

regulator_sections <- c(
  "indra-upstream", "indra-upstream-kinases", "indra-upstream-tfs", "indra-downstream"
)
discrete_results <- map(regulator_sections, \(s) {
  as_tibble(cache$discrete[[s]]) |>
    transmute(analysis = s, regulator = name, curie, p_value = p, q_value = q)
}) |>
  list_rbind() |>
  arrange(analysis, p_value)

# signed_analysis reports raw binomial p-values; BH runs here, within each contrast.
signed_results <- cache$signed |>
  filter(!is.na(p_value)) |>
  mutate(q_value = p.adjust(p_value, "BH"), .by = contrast) |>
  arrange(contrast, p_value)

summary_counts <- bind_rows(
  summarise(discrete_results,
    tested = n(), q_below_cut = sum(q_value < FDR_CUT), min_q = min(q_value), .by = analysis
  ) |>
    mutate(input = "female_FDR_0.10", .before = 1),
  summarise(signed_results,
    tested = n(), q_below_cut = sum(q_value < FDR_CUT), min_q = min(q_value), .by = contrast
  ) |>
    transmute(
      input = str_glue("{contrast}_pi_up_down"), analysis = "signed", tested,
      q_below_cut, min_q
    )
)

coverage <- lists |>
  summarise(
    n_proteins = n(), n_human = length(human_of(protein)),
    .by = c(list, contrast, direction)
  )

settings <- tibble(
  setting = c(
    "STRING version", "STRING min score", "STRING background", "INDRA background",
    "INDRA min evidence", "fetched"
  ),
  value = c(
    STRING_VERSION, STRING_SCORE, str_glue("{cache$string_mapped} detected proteins"),
    str_glue("{length(human_background)} human orthologs of detected proteins"),
    MIN_EVIDENCE, format(cache$fetched)
  )
)

top_bars <- function(d, title) {
  d <- slice_min(d, p_value, n = N_TOP, with_ties = FALSE)
  barplot(rev(-log10(d$p_value)),
    names.arg = rev(str_trunc(d$regulator, 28)), horiz = TRUE, las = 1,
    cex.names = 0.6, xlab = "-log10 p", main = title, cex.main = 0.8
  )
  mtext(str_glue("smallest q = {signif(min(d$q_value), 2)}"), side = 3, line = 0, cex = 0.6)
}
# A list with no expected edges has no fold to draw; the table keeps it.
drawn <- filter(string_results, expected > 0)
plot_index <- write_pdf(file.path(out$reports, "06_mechanism.pdf"), list(
  "STRING PPI enrichment per list" = \() {
    par(mar = c(4, 14, 3, 1))
    y <- rev(seq_len(nrow(drawn)))
    plot(log2(drawn$fold), y,
      pch = if_else(drawn$FDR < FDR_CUT, 19, 1), yaxt = "n", ylab = "",
      xlab = "log2 observed / expected STRING edges (STRINGdb, detected background)",
      main = "STRING PPI enrichment (filled: FDR ≤ 0.05)", xlim = range(0, log2(drawn$fold))
    )
    axis(2, y, str_glue("{drawn$list} ({drawn$n_mapped})"), las = 1, cex.axis = 0.7)
    abline(v = 0, lty = 2)
  },
  "INDRA signed regulator analysis, female and male" = \() {
    par(mfrow = c(1, 2), mar = c(4, 10, 3, 1))
    walk(SIGNED_CONTRASTS, \(cn) {
      top_bars(filter(signed_results, contrast == cn), str_glue("signed: {cn}"))
    })
  },
  "INDRA regulator enrichment, female FDR list" = \() {
    par(mfrow = c(2, 2), mar = c(4, 10, 3, 1))
    walk(regulator_sections, \(s) top_bars(filter(discrete_results, analysis == s), s))
  }
))

write_workbook(
  list(
    string = string_results, indra_summary = summary_counts,
    indra_discrete = filter(discrete_results, p_value < FDR_CUT),
    indra_signed = filter(signed_results, p_value < FDR_CUT),
    coverage = coverage, settings = settings, plot_index = plot_index
  ),
  c(
    "STRING PPI enrichment per list: proteins, mapped, edges, expected, fold, p, BH FDR.",
    "INDRA regulators tested and reaching q ≤ 0.05, per analysis.",
    "INDRA discrete regulator enrichment on the female FDR list, p ≤ 0.05.",
    "INDRA signed analysis on the Π lists, p ≤ 0.05, BH within contrast.",
    "Proteins and human orthologs per list.",
    "Service versions, thresholds and fetch date.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

string_results
summary_counts

sessionInfo()
