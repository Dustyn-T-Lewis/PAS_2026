#!/usr/bin/env Rscript
# Co-abundance modules by WGCNA's standard workflow, built after removing sex and keeping treatment
# (limma::removeBatchEffect), so the network follows co-variation within each sex instead of the
# sex difference that dominates this proteome (Johnson et al. 2020, Nat Med; 2022, Nat Neurosci).
# Eigengenes are then taken from the unadjusted matrix, so 07_Fit can test sex and treatment on
# them. Each module is named only by an over-represented term at FDR <= 0.05.

source(here::here("R", "helpers.R"))
library(WGCNA)
WGCNA::disableWGCNAThreads()

out <- stage_paths("00_Trunk", "06_Modules")
inputs <- c(
  proteins = here("00_Trunk", "02_Quantification", "c_data", "proteins.rds"),
  gene_sets = here("00_Trunk", "04_Gene_Sets", "c_data", "gene_sets.rds")
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

fit_plot <- scale_free |>
  pivot_longer(-power, names_to = "measure") |>
  mutate(measure = if_else(measure == "r2", "Signed scale-free R²", "Mean connectivity")) |>
  ggplot(aes(power, value)) +
  geom_line(linewidth = 0.3) +
  geom_point(size = 1.2) +
  geom_vline(xintercept = power, linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~measure, scales = "free_y") +
  labs(
    x = "Soft-threshold power", y = NULL, title = "Soft-threshold choice",
    subtitle = str_glue("Power used: {power}; R² cut {RSQ_CUT}")
  ) +
  theme_figure()
eigengene_plot <- as_tibble(eigengenes, rownames = "sample_id") |>
  pivot_longer(-sample_id, names_to = "module", values_to = "eigengene") |>
  left_join(select(targets, sample_id, group, sex), by = "sample_id") |>
  mutate(group = factor(group, FIGURE_GROUPS)) |>
  ggplot(aes(group, eigengene, colour = group, shape = sex)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.25, colour = "grey40") +
  geom_point(position = position_jitter(width = 0.12, seed = SEED), size = 1.4) +
  facet_wrap(~module) +
  scale_colour_manual(values = GROUP_COLOURS, guide = "none") +
  scale_shape_manual(values = SEX_SHAPES, guide = "none") +
  labs(
    x = NULL, y = "Eigengene (unadjusted data)", title = "Module eigengenes by mouse",
    subtitle = "Modules built with sex removed; eigengenes keep sex and treatment for testing"
  ) +
  theme_figure()
plot_index <- write_report(file.path(out$reports, "06_modules.pdf"), list(
  "soft threshold: scale-free fit and connectivity" = fit_plot,
  "eigengenes by group" = eigengene_plot,
  "eigengene heatmap" = group_heatmap(t(eigengenes), targets,
    title = "Module eigengenes", subtitle = "Rows z-scored; clustered for display only"
  )
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
