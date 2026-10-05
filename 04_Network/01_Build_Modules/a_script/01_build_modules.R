#!/usr/bin/env Rscript
# Co-abundance modules by WGCNA's standard workflow: soft-threshold power, a signed network built
# blockwise, module eigengenes and signed module membership (kME). Settings are those of the first
# build, so the modules every later step reads stay the same.

source(here::here("R", "helpers.R"))
library(WGCNA)
WGCNA::disableWGCNAThreads()

out <- stage_paths("04_Network", "01_Build_Modules")
inputs <- c(
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds")
)

RSQ_CUT <- 0.85
# WGCNA's FAQ power for a signed network of 20 to 30 samples when no power reaches the cut.
FALLBACK_POWER <- 16
N_HUBS <- 10

de <- readRDS(inputs[["contrasts"]])
datExpr <- t(de$E_norm$full)
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
eigengenes <- orderMEs(moduleEigengenes(datExpr, net$colors, excludeGrey = TRUE)$eigengenes)
colnames(eigengenes) <- str_remove(colnames(eigengenes), "^ME")
kme <- signedKME(datExpr, eigengenes, corFnc = "bicor", corOptions = "maxPOutliers = 0.05")
colnames(kme) <- str_remove(colnames(kme), "^kME")

modules <- tibble(
  protein = colnames(datExpr), gene = str_split_i(de$genes$Genes, ";", 1),
  is_mito = de$genes$is_mito, module = unname(net$colors)
) |>
  # Grey has no eigengene, so its membership is NA.
  mutate(kME = as.matrix(kme)[cbind(match(protein, rownames(kme)), match(module, colnames(kme)))])
sizes <- count(modules, module, name = "proteins") |>
  left_join(summarise(modules, mito_share = mean(is_mito), .by = module), by = "module") |>
  arrange(desc(proteins))
hubs <- modules |>
  filter(module != "grey") |>
  slice_max(kME, n = N_HUBS, by = module, with_ties = FALSE)

# The first build's modules, checked so downstream results stay comparable.
stopifnot(
  setequal(sizes$module, c(
    "turquoise", "blue", "brown", "yellow", "grey", "green", "red",
    "black", "pink", "magenta"
  )),
  sizes$proteins[sizes$module == "turquoise"] == 520, sizes$proteins[sizes$module == "grey"] == 153
)

scale_free <- tibble(
  power = sft$fitIndices$Power, r2 = with(sft$fitIndices, -sign(slope) * SFT.R.sq),
  mean_connectivity = sft$fitIndices$mean.k.
)
plot_index <- write_pdf(file.path(out$reports, "01_build_modules.pdf"), list(
  "soft threshold: scale-free fit and connectivity by power" = \() {
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    plot(scale_free$power, scale_free$r2,
      type = "b", pch = 19, xlab = "power",
      ylab = "signed scale-free R²", main = str_glue("power used: {power}")
    )
    abline(h = RSQ_CUT, lty = 2, col = "#B2182B")
    plot(scale_free$power, scale_free$mean_connectivity,
      type = "b", pch = 19, xlab = "power",
      ylab = "mean connectivity", main = "mean connectivity"
    )
  },
  "protein dendrogram with module colours" = \() {
    plotDendroAndColors(net$dendrograms[[1]], net$colors[net$blockGenes[[1]]], "module",
      dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
      main = "Protein dendrogram and modules"
    )
  },
  "eigengene overview heatmap" = \() {
    title <- "Module eigengenes by mouse, rows clustered (display only)"
    ComplexHeatmap::draw(overview_heatmap(t(eigengenes), de$targets$group, title))
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
    eigengenes = as_tibble(eigengenes, rownames = "mouse"), plot_index = plot_index
  ),
  c(
    "Proteins per module and their MitoCarta share.",
    "Every protein: module and signed module membership (kME, bicor).",
    str_glue("The {N_HUBS} highest-kME proteins of each module."),
    "Scale-free fit and mean connectivity by soft-threshold power.",
    "Module eigengenes by mouse.",
    "Which PDF page shows what."
  ),
  out$workbook, inputs
)

sizes

sessionInfo()
