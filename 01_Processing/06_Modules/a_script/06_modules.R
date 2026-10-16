# Co-abundance modules by WGCNA's standard workflow on the proteins detected in every mouse, built
# after removing sex and measured run depth (limma::removeBatchEffect, treatment kept), so the
# network follows co-variation within each sex instead of the sex difference or how deep a run
# went (Johnson et al. 2020, Nat Med; Parsana et al. 2019, Genome Biol). Values limpa infers for
# undetected proteins are tied within a protein and would group proteins by when they go missing,
# so only fully detected proteins enter. Eigengenes come from the same adjusted matrix the modules
# were built on, as WGCNA's authors advise; sex differences on them are therefore removed by
# design and only within-sex contrasts are read. Stability: the network is rebuilt on 100
# subsamples of 16 mice (four per group) and each module's best-match Jaccard is reported. Each
# module is named only by an over-represented term at FDR <= 0.05.

source(here::here("00_Setup", "Shared_functions", "helpers.R"))
library(WGCNA)
WGCNA::disableWGCNAThreads()

out <- stage_paths("01_Processing", "06_Modules")
inputs <- c(
  proteins = here("01_Processing", "02_Quantification", "c_data", "proteins.rds"),
  filtering = here("01_Processing", "01_Filtering", "c_data", "precursors_filtered.rds"),
  gene_sets = here("01_Processing", "04_Gene_Sets", "c_data", "gene_sets.rds")
)
RSQ_CUT <- 0.85
# WGCNA's FAQ power for a signed network of 20 to 30 samples when no power reaches the cut.
FALLBACK_POWER <- 16
N_HUBS <- 10
N_RESAMPLES <- 100

quant <- readRDS(inputs[["proteins"]])
gs <- readRDS(inputs[["gene_sets"]])
targets <- quant$proteins$targets
run_qc <- readRDS(inputs[["filtering"]])$run_qc
run_depth <- log2(run_qc$precursors[match(targets$sample_id, run_qc$sample_id)])
detected <- quant$proteins$other$n.observations[, targets$sample_id] > 0
complete <- rownames(detected)[rowSums(detected) == ncol(detected)]
expr <- quant$E_norm$full[complete, targets$sample_id]
adjusted <- limma::removeBatchEffect(
  expr,
  batch = targets$sex, covariates = run_depth, design = model.matrix(~treatment, targets)
)
datExpr <- t(adjusted)
stopifnot(goodSamplesGenes(datExpr, verbose = 0)$allOK)

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

eigengenes <- moduleEigengenes(datExpr, net$colors, excludeGrey = TRUE)$eigengenes
colnames(eigengenes) <- str_remove(colnames(eigengenes), "^ME")
kme <- signedKME(datExpr, moduleEigengenes(datExpr, net$colors, excludeGrey = TRUE)$eigengenes,
  corFnc = "bicor", corOptions = "maxPOutliers = 0.05"
)
colnames(kme) <- str_remove(colnames(kme), "^kME")

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

# Stability: rebuild on subsamples of 16 mice (four per group) and match each module to its best
# overlap in the rebuilt network.
set.seed(SEED)
stability <- map(seq_len(N_RESAMPLES), \(i) {
  keep <- unlist(map(split(seq_len(nrow(datExpr)), targets$group), \(idx) sample(idx, 4)))
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
    plotEigengeneNetworks(orderMEs(moduleEigengenes(datExpr, net$colors)$eigengenes), "",
      marDendro = c(0, 4, 1, 2), marHeatmap = c(3, 4, 1, 2), xLabelsAngle = 90
    )
  }
))

saveRDS(
  list(
    eigengenes = eigengenes, modules = modules, kme = kme, sizes = sizes, hubs = hubs,
    power = power, colors = net$colors, stability = stability, proteins_used = complete
  ),
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
    "Best-match Jaccard of each module in each of the 100 subsample rebuilds.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

sizes

sessionInfo()
