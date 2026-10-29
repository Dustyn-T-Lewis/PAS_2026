# What went into each pellet. A mitochondria-enriched pellet also carries sarcoplasmic reticulum,
# myofibril, nerve, blood and vessel material in varying amounts, and a PAS effect on any of those
# could pass for a change in the tissue. Each mouse gets a score per marker group: the mean, over
# the group's proteins detected in at least MIN_MICE_DETECTED mice, of each protein's log2
# abundance centred and scaled across mice. Haemoglobin, filtered out of the proteome as a
# contaminant, is scored as its share of DIA-NN's normalised precursor signal, with the q-value
# filters 01_Filtering applies. Each score and the pellet MitoCarta share are tested for PAS within
# each sex (Welch t), with the difference and its 95% CI.

source(here::here("00_Setup", "Shared_functions", "helpers.R"))

MIN_MICE_DETECTED <- 15
MARKERS <- list(
  `Sarcoplasmic reticulum` = c(
    "Atp2a1", "Atp2a2", "Casq1", "Casq2", "Trdn", "Jph1", "Jph2", "Ryr1", "Srl", "Stac3",
    "Cacna1s", "Asph", "Sln", "Mypn"
  ),
  `Endoplasmic reticulum` = c("Canx", "Calr", "Hspa5", "Pdia3", "Pdia6", "P4hb", "Ssr1", "Sec61a1"),
  Myofibril = c(
    "Myh1", "Myh2", "Myh4", "Myh7", "Acta1", "Actn2", "Actn3", "Tnni1", "Tnni2", "Tnnt1",
    "Tnnt3", "Tnnc1", "Tnnc2", "Tpm1", "Tpm2", "Tpm3", "Myl1", "Myl2", "Mylpf", "Myom1", "Ttn",
    "Neb"
  ),
  Nerve = c("Mbp", "Mpz", "Prx", "Pmp22", "Nefl", "Nefm", "Nefh", "Prph", "S100b"),
  `Red cell` = c("Slc4a1", "Spta1", "Sptb", "Epb42", "Epb41", "Add2", "Gypa", "Rhag"),
  Vessel = c("Pecam1", "Cdh5", "Vwf", "Cav1", "Aoc3", "Cd34", "Eng", "Kdr", "Emcn")
)
HAEMOGLOBIN <- "Haemoglobin (% of signal)"
MITO_SHARE <- "MitoCarta share (% of signal)"

out <- stage_paths("02_Statistics", "05_Composition")
inputs <- c(
  proteins = output_path("01_Processing", "02_Quantification", "proteins.rds"),
  report = input_path("report.parquet")
)
quant <- readRDS(inputs[["proteins"]])
targets <- as_tibble(quant$proteins$targets)
expr <- quant$E_norm$full[, targets$sample_id]
detected <- quant$proteins$other$n.observations[rownames(expr), targets$sample_id] > 0
genes <- first_symbol(quant$proteins$genes[rownames(expr), "Genes"])

usable <- rowSums(detected) >= MIN_MICE_DETECTED
scaled <- t(scale(t(expr)))
marker_members <- imap(MARKERS, \(symbols, group) {
  tibble(group = group, protein = rownames(expr)[usable & genes %in% symbols]) |>
    mutate(gene = genes[protein])
}) |>
  list_rbind()
scores <- marker_members |>
  summarise(proteins = list(protein), n_proteins = n(), .by = group) |>
  mutate(score = map(proteins, \(p) colMeans(scaled[p, , drop = FALSE]))) |>
  select(group, n_proteins, score) |>
  unnest_longer(score, indices_to = "sample_id")

report <- nanoparquet::read_parquet(inputs[["report"]],
  col_select = c("Run", "Genes", "Precursor.Normalised", "Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value")
) |>
  filter(
    Q.Value <= Q_CUTOFF, Lib.Q.Value <= Q_CUTOFF, Lib.PG.Q.Value <= Q_CUTOFF, Run %in% targets$run
  )
haemoglobin <- report |>
  summarise(
    score = 100 * sum(Precursor.Normalised[str_detect(Genes, "^Hb[ab]")]) /
      sum(Precursor.Normalised),
    .by = Run
  ) |>
  transmute(
    group = HAEMOGLOBIN, n_proteins = NA_integer_,
    sample_id = targets$sample_id[match(Run, targets$run)], score
  )

composition <- bind_rows(
  scores,
  haemoglobin,
  tibble(
    group = MITO_SHARE, n_proteins = NA_integer_,
    sample_id = names(quant$pct_mito), score = unname(quant$pct_mito)
  )
) |>
  left_join(select(targets, sample_id, sex, treatment, group_label = group), by = "sample_id")

tests <- composition |>
  reframe(
    {
      w <- t.test(score[treatment == "PAS"], score[treatment == "VEH"])
      tibble(
        mean_pas = mean(score[treatment == "PAS"]), mean_veh = mean(score[treatment == "VEH"]),
        difference = unname(w$estimate[1] - w$estimate[2]),
        ci_low = w$conf.int[1], ci_high = w$conf.int[2], p = w$p.value
      )
    },
    .by = c(group, sex)
  ) |>
  mutate(sex = if_else(sex == "F", "female", "male"))

plot_index <- preview_pdf(file.path(out$reports, "05_composition.pdf"), list(
  "marker scores by group" = \() {
    groups <- unique(composition$group)
    par(mfrow = c(2, 4), mar = c(5, 4, 3, 1))
    walk(groups, \(g) {
      d <- filter(composition, group == g)
      boxplot(score ~ factor(group_label, FIGURE_GROUPS), d,
        las = 2, xlab = "", ylab = "score", main = g, cex.main = 0.9
      )
      stripchart(score ~ factor(group_label, FIGURE_GROUPS), d,
        vertical = TRUE, add = TRUE, pch = 16
      )
    })
  }
))

saveRDS(
  list(
    composition = composition, tests = tests,
    groups = list(markers = names(MARKERS), haemoglobin = HAEMOGLOBIN, mito_share = MITO_SHARE)
  ),
  file.path(out$data, "composition.rds")
)
write_workbook(
  list(
    tests = tests, scores = composition, members = marker_members, plot_index = plot_index
  ),
  c(
    "Per marker group and sex: PAS and vehicle means, difference with 95% CI, Welch p.",
    "Per mouse and marker group: score (mean scaled log2 abundance, or % of signal).",
    str_glue("Proteins scored per marker group (detected in at least {MIN_MICE_DETECTED} mice)."),
    "Which PDF page shows what."
  ),
  out$workbook
)

count(marker_members, group, name = "proteins")
tests
