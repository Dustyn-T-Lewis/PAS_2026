# Co-abundance modules by WGCNA's standard workflow on the proteins detected in every mouse, after
# removing sex and run depth (treatment kept), so the network follows co-variation within each sex.
# Values limpa infers for undetected proteins would group proteins by when they go missing, so only
# fully detected proteins enter; eigengenes come from the same adjusted matrix, so only within-sex
# contrasts are read on them.

source(here::here("00_Setup", "Shared_functions", "helpers.R"))
library(WGCNA)
disableWGCNAThreads()

out <- stage_paths("01_Processing", "06_Modules")
inputs <- c(
  proteins = output_path("01_Processing", "02_Quantification", "proteins.rds"),
  filtering = output_path("01_Processing", "01_Filtering", "precursors_filtered.rds"),
  gene_sets = output_path("01_Processing", "04_Gene_Sets", "gene_sets.rds")
)
RSQ_CUT <- 0.85
# WGCNA's FAQ power for a signed network of 20 to 30 samples when no power reaches the cut.
FALLBACK_POWER <- 16
N_HUBS <- 10
N_RESAMPLES <- 100
RESAMPLE_PER_GROUP <- 4

quant <- readRDS(inputs[["proteins"]])
gs <- readRDS(inputs[["gene_sets"]])
targets <- quant$proteins$targets
run_qc <- readRDS(inputs[["filtering"]])$run_qc
run_depth <- run_depth_log2(run_qc, targets$sample_id)
detected <- quant$proteins$other$n.observations[, targets$sample_id] > 0
complete <- rownames(detected)[rowSums(detected) == ncol(detected)]
expr <- quant$E_norm$full[complete, targets$sample_id]
adjusted <- limma::removeBatchEffect(
  expr,
  batch = targets$sex, covariates = run_depth, design = model.matrix(~treatment, targets)
)
datExpr <- t(adjusted)

# pickSoftThreshold prints its fit table whatever verbose says.
invisible(capture.output(
  sft <- pickSoftThreshold(datExpr,
    networkType = "signed", corFnc = bicor, corOptions = list(maxPOutliers = 0.05),
    RsquaredCut = RSQ_CUT, powerVector = c(1:10, seq(12, 30, 2)), verbose = 0
  )
))
power <- if (is.na(sft$powerEstimate)) FALLBACK_POWER else sft$powerEstimate
build <- \(x) {
  blockwiseModules(x,
    power = power, networkType = "signed", TOMType = "signed", corType = "bicor",
    maxPOutliers = 0.05, minModuleSize = 20, mergeCutHeight = 0.25, deepSplit = 2,
    numericLabels = FALSE, maxBlockSize = 5000, saveTOMs = FALSE, verbose = 0
  )
}
net <- build(datExpr)

module_eigengenes <- moduleEigengenes(datExpr, net$colors, excludeGrey = TRUE)$eigengenes
kme <- signedKME(datExpr, module_eigengenes, corFnc = "bicor", corOptions = "maxPOutliers = 0.05")
colnames(kme) <- str_remove(colnames(kme), "^kME")
eigengenes <- rename_with(module_eigengenes, \(x) str_remove(x, "^ME"))

gene <- first_symbol(quant$proteins$genes[complete, "Genes"])
modules <- tibble(
  protein = complete, gene = gene, is_mito = quant$proteins$genes[complete, "is_mito"],
  module = unname(net$colors)
) |>
  mutate(kME = as.matrix(kme)[cbind(match(protein, rownames(kme)), match(module, colnames(kme)))])
module_genes <- modules |>
  filter(module != "grey", gene != "") |>
  summarise(genes = list(unique(gene)), .by = module)
labels <- label_groups(
  set_names(module_genes$genes, module_genes$module), gs$full[c("Hallmark", "Reactome")],
  setdiff(unique(gene), c("", NA))
) |>
  mutate(module = module_genes$module, .before = 1)

# Stability: rebuild on subsamples of RESAMPLE_PER_GROUP mice per group and match each module to
# its best overlap in the rebuilt network.
set.seed(SEED)
stability <- map(seq_len(N_RESAMPLES), \(i) {
  keep <- split(seq_len(nrow(datExpr)), targets$group) |>
    map(\(idx) sample(idx, RESAMPLE_PER_GROUP)) |>
    unlist()
  colors <- build(datExpr[keep, ])$colors
  map(setdiff(unique(net$colors), "grey"), \(m) {
    members <- names(net$colors)[net$colors == m]
    best <- max(map_dbl(setdiff(unique(colors), "grey"), \(k) {
      other <- names(colors)[colors == k]
      length(intersect(members, other)) / length(union(members, other))
    }), 0)
    tibble(module = m, resample = i, jaccard = best)
  }) |>
    list_rbind()
}) |>
  list_rbind()
stability_summary <- stability |>
  summarise(
    median_jaccard = median(jaccard), q25 = quantile(jaccard, 0.25),
    q75 = quantile(jaccard, 0.75), .by = module
  )

# Pellet mitochondrial share, overall and within group, beside each eigengene: a module that only
# tracks how mitochondria-rich a pellet was would show it here.
within_group <- \(x) x - ave(x, targets$group)
sizes <- count(modules, module, name = "proteins") |>
  left_join(summarise(modules, mito_share = mean(is_mito), .by = module), by = "module") |>
  left_join(select(labels, module, label, label_fdr), by = "module") |>
  left_join(stability_summary, by = "module") |>
  left_join(
    tibble(
      module = colnames(eigengenes),
      r_pellet_mito = map_dbl(eigengenes, \(e) cor(e, quant$pct_mito)),
      r_pellet_mito_within_group = map_dbl(
        eigengenes, \(e) cor(within_group(e), within_group(quant$pct_mito))
      )
    ),
    by = "module"
  ) |>
  arrange(desc(proteins))
hubs <- modules |>
  filter(module != "grey") |>
  slice_max(kME, n = N_HUBS, by = module, with_ties = FALSE)
scale_free <- tibble(
  power = sft$fitIndices$Power, r2 = with(sft$fitIndices, -sign(slope) * SFT.R.sq),
  mean_connectivity = sft$fitIndices$mean.k.
)

plot_index <- preview_pdf(file.path(out$reports, "06_modules.pdf"), list(
  "soft threshold: scale-free fit and connectivity" = \() {
    par(mfrow = c(1, 2))
    plot(scale_free$power, scale_free$r2,
      type = "b", pch = 16, xlab = "Soft-threshold power", ylab = "Signed scale-free R²"
    )
    abline(h = RSQ_CUT, v = power, lty = 2)
    plot(scale_free$power, scale_free$mean_connectivity,
      type = "b", pch = 16, xlab = "Soft-threshold power", ylab = "Mean connectivity"
    )
    abline(v = power, lty = 2)
  },
  "module dendrogram" = \() {
    plotDendroAndColors(net$dendrograms[[1]], net$colors[net$blockGenes[[1]]], "Module",
      dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05
    )
  },
  "eigengene network" = \() {
    plotEigengeneNetworks(orderMEs(module_eigengenes), "",
      marDendro = c(0, 4, 1, 2), marHeatmap = c(3, 4, 1, 2), xLabelsAngle = 90
    )
  }
))

saveRDS(
  list(eigengenes = eigengenes, modules = modules, sizes = sizes),
  file.path(out$data, "modules.rds")
)
write_workbook(
  list(
    sizes = sizes, modules = modules, hubs = hubs, scale_free = scale_free,
    eigengenes = as_tibble(eigengenes, rownames = "sample_id"), stability = stability,
    plot_index = plot_index
  ),
  c(
    paste(
      "Per module: proteins, MitoCarta share, FDR-gated label, resampling stability (median",
      "and quartiles of best-match Jaccard), r with pellet mitochondrial share."
    ),
    "Every fully detected protein: module and signed membership (kME, bicor, adjusted data).",
    str_glue("The {N_HUBS} highest-kME proteins of each module."),
    "Scale-free fit and mean connectivity by soft-threshold power.",
    "Module eigengenes by mouse, from the sex- and run-depth-adjusted matrix.",
    str_glue("Best-match Jaccard of each module in each of the {N_RESAMPLES} subsample rebuilds."),
    "Which PDF page shows what."
  ),
  out$workbook
)

tibble(power = power, proteins = length(complete))
sizes
