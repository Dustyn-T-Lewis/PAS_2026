# Co-abundance modules by WGCNA's standard workflow, built after removing sex and keeping treatment
# (limma::removeBatchEffect), so the network follows co-variation within each sex instead of the
# sex difference that dominates this proteome (Johnson et al. 2020, Nat Med; 2022, Nat Neurosci).
# Eigengenes are then taken from the unadjusted matrix, so 01_Fit can test sex and treatment on
# them. Each module is named only by an over-represented term at FDR <= 0.05.

source(here::here("00_Setup", "Shared_functions", "helpers.R"))
library(WGCNA)
WGCNA::disableWGCNAThreads()

out <- stage_paths("01_Processing", "06_Modules")
inputs <- c(
  proteins = here("01_Processing", "02_Quantification", "c_data", "proteins.rds"),
  gene_sets = here("01_Processing", "04_Gene_Sets", "c_data", "gene_sets.rds")
)
RSQ_CUT <- 0.85
# WGCNA's FAQ power for a signed network of 20 to 30 samples when no power reaches the cut.
FALLBACK_POWER <- 16
N_HUBS <- 10

quant <- readRDS(inputs[["proteins"]])
gs <- readRDS(inputs[["gene_sets"]])
targets <- quant$proteins$targets
expr <- quant$E_norm$full
adjusted <- limma::removeBatchEffect(
  expr,
  batch = targets$sex, design = model.matrix(~treatment, targets)
)
datExpr <- t(adjusted)
stopifnot(goodSamplesGenes(datExpr, verbose = 0)$allOK)

# pickSoftThreshold prints its fit table whatever verbose says.
invisible(capture.output(
  sft <- pickSoftThreshold(datExpr,
    networkType = "signed", corFnc = bicor, corOptions = list(maxPOutliers = 0.05),
    RsquaredCut = RSQ_CUT, powerVector = c(1:10, seq(12, 20, 2)), verbose = 0
  )
))
power <- if (is.na(sft$powerEstimate)) FALLBACK_POWER else sft$powerEstimate
net <- blockwiseModules(datExpr,
  power = power, networkType = "signed", TOMType = "signed", corType = "bicor",
  maxPOutliers = 0.05, minModuleSize = 20, mergeCutHeight = 0.25, deepSplit = 2,
  numericLabels = FALSE, maxBlockSize = 5000, saveTOMs = FALSE, verbose = 0
)

eigengenes <- moduleEigengenes(t(expr), net$colors, excludeGrey = TRUE)$eigengenes
colnames(eigengenes) <- str_remove(colnames(eigengenes), "^ME")
kme <- signedKME(datExpr, moduleEigengenes(datExpr, net$colors, excludeGrey = TRUE)$eigengenes,
  corFnc = "bicor", corOptions = "maxPOutliers = 0.05"
)
colnames(kme) <- str_remove(colnames(kme), "^kME")

gene <- first_symbol(quant$proteins$genes$Genes)
modules <- tibble(
  protein = rownames(expr), gene = gene, is_mito = quant$proteins$genes$is_mito,
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

# Pellet mitochondrial share, overall and within group, beside each eigengene: a module that only
# tracks how mitochondria-rich a pellet was would show it here.
within_group <- \(x) x - ave(x, targets$group)
sizes <- count(modules, module, name = "proteins") |>
  left_join(summarise(modules, mito_share = mean(is_mito), .by = module), by = "module") |>
  left_join(select(labels, module, label, label_fdr), by = "module") |>
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
    plotEigengeneNetworks(orderMEs(moduleEigengenes(t(expr), net$colors)$eigengenes), "",
      marDendro = c(0, 4, 1, 2), marHeatmap = c(3, 4, 1, 2), xLabelsAngle = 90
    )
  }
))

saveRDS(
  list(
    eigengenes = eigengenes, modules = modules, kme = kme, sizes = sizes, hubs = hubs,
    power = power, colors = net$colors
  ),
  file.path(out$data, "modules.rds")
)
write_workbook(
  list(
    sizes = sizes, modules = modules, hubs = hubs, scale_free = scale_free,
    eigengenes = as_tibble(eigengenes, rownames = "sample_id"), plot_index = plot_index
  ),
  c(
    "Per module: proteins, MitoCarta share, FDR-gated label, r with pellet mitochondrial share.",
    "Every protein: module and signed membership (kME, bicor, sex-adjusted data).",
    str_glue("The {N_HUBS} highest-kME proteins of each module."),
    "Scale-free fit and mean connectivity by soft-threshold power.",
    "Module eigengenes by mouse, from the unadjusted matrix.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

sizes

sessionInfo()
