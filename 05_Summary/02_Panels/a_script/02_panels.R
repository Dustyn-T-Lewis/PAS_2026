#!/usr/bin/env Rscript
# Whether a combination of features separates the two groups of each comparison in mice a model has
# not seen. Two methods: mixOmics sPLS-DA/DIABLO (features kept per block fixed in advance) and
# nestedcv's elastic net (penalty tuned in an inner leave-one-out loop). Each runs on the features
# of one level alone, on those features with the phenotypes, and on the phenotypes alone. Scored by
# balanced error rate in held-out mice: leave-one-out for DIABLO (mixOmics perf takes no custom
# folds), leave-pair-out for the elastic net. The null is every one of the 126 ways to split the ten
# mice into two groups of five, refitted the same way, so p is exact. Results are cached.

source(here::here("R", "helpers.R"))

KEEP_X <- list(proteins = 20, sets = 10, modules = 3, phenotypes = 3)
DESIGN_WEIGHT <- 0.1
CORES <- max(1, parallel::detectCores() - 2)
LEVELS <- list(
  proteins = "proteins", sets = "sets", modules = "modules",
  combined = c("proteins", "sets", "modules")
)

out <- stage_paths("05_Summary", "02_Panels")
inputs <- c(
  contrasts = here("02_Differential_Expression", "02_Contrasts", "c_data", "contrasts.rds"),
  scores = here("03_Pathway_Enrichment", "01_Scores", "c_data", "set_scores.rds"),
  modules = here("04_Network", "01_Build_Modules", "c_data", "modules.rds"),
  phenotypes = here("00_Input", "phenotypes.xlsx")
)
CACHE <- file.path(out$data, "panels_cache.rds")

de <- readRDS(inputs[["contrasts"]])
sc <- readRDS(inputs[["scores"]])
mods <- readRDS(inputs[["modules"]])
samples <- rownames(de$targets)
ph <- read_phenotypes(samples)
blocks <- list(
  proteins = t(de$E_norm$full),
  sets = t(rbind(sc$scores$full, sc$scores$mito)),
  modules = as.matrix(mods$eigengenes),
  phenotypes = `rownames<-`(as.matrix(ph[PHENOTYPES]), samples)
)
stopifnot(all(map_lgl(blocks, \(b) identical(rownames(b), samples))))

configs <- bind_rows(
  expand_grid(level = names(LEVELS), input = c("features", "features + phenotypes")),
  tibble(level = "phenotypes", input = "phenotypes only")
) |>
  mutate(block_names = map2(level, input, \(l, i) {
    if (i == "phenotypes only") "phenotypes" else c(LEVELS[[l]], if (i != "features") "phenotypes")
  }))

diablo_ber <- function(X, y) {
  if (length(X) == 1) {
    fit <- mixOmics::splsda(X[[1]], y, ncomp = 1, keepX = KEEP_X[[names(X)]])
    pf <- mixOmics::perf(fit, validation = "loo", dist = "centroids.dist", progressBar = FALSE)
    return(list(ber = pf$error.rate$BER[1, 1], fit = fit))
  }
  design <- matrix(DESIGN_WEIGHT, length(X), length(X), dimnames = list(names(X), names(X)))
  diag(design) <- 0
  fit <- suppressMessages(mixOmics::block.splsda(X, y,
    ncomp = 1, keepX = KEEP_X[names(X)],
    design = design
  ))
  pf <- mixOmics::perf(fit, validation = "loo", dist = "centroids.dist", progressBar = FALSE)
  list(ber = pf$WeightedVote.error.rate$centroids.dist["Overall.BER", 1], fit = fit)
}
# Leave-pair-out (Airola et al. 2011): each outer fold holds one mouse from each group, so training
# stays balanced at four against four. Leave-one-out unbalances it and drives a penalised model's
# null error towards 1, not 0.5.
glmnet_ber <- function(X, y) {
  x <- do.call(cbind, unname(X))
  # Each mouse is predicted in five pairs; nestedcv labels predictions by row name.
  rownames(x) <- NULL
  pairs <- expand_grid(a = which(y == levels(y)[1]), b = which(y == levels(y)[2]))
  set.seed(SEED)
  fit <- muffle(
    nestedcv::nestcv.glmnet(y, x,
      family = "binomial", outer_folds = map2(pairs$a, pairs$b, c),
      n_inner_folds = length(y) - 2, alphaSet = 0.5, min_1se = 1, finalCV = FALSE,
      verbose = FALSE
    ),
    "convergence|one multinomial|observations"
  )
  wrong <- fit$output$predy != fit$output$testy
  list(
    ber = mean(tapply(wrong, fit$output$testy, mean)),
    empty = mean(map_int(fit$outer_result, \(r) length(r$coef) - 1L) == 0),
    fit = fit
  )
}

splits <- combn(10, 5)
splits <- splits[, splits[1, ] == 1]
grid <- expand_grid(comparison = names(COMPARISONS), configs)
# Each method's observed fit and its 126-split null; cached per method, keyed on what it saw.
fit_grid <- function(score) {
  pmap(list(grid$comparison, grid$block_names), \(cmp, bn) {
    groups <- COMPARISONS[[cmp]]
    keep <- ph$group %in% groups
    X <- map(blocks[bn], \(b) b[keep, , drop = FALSE])
    y <- factor(ph$group[keep], levels = groups)
    relabel <- \(s) factor(if_else(seq_along(y) %in% s, groups[1], groups[2]), levels = groups)
    null <- parallel::mclapply(seq_len(ncol(splits)), \(i) score(X, relabel(splits[, i]))$ber,
      mc.cores = CORES
    )
    list(observed = score(X, y), null = unlist(null))
  }, .progress = "panels")
}
keys <- list(
  diablo = rlang::hash(list(blocks, configs$block_names, KEEP_X, DESIGN_WEIGHT, "loo")),
  glmnet = rlang::hash(list(blocks, configs$block_names, "leave-pair-out", 0.5))
)
cache <- if (file.exists(CACHE)) readRDS(CACHE) else list()
for (m in names(keys)) {
  if (!identical(cache[[m]]$key, keys[[m]])) {
    fitter <- list(diablo = diablo_ber, glmnet = glmnet_ber)[[m]]
    cache[[m]] <- list(key = keys[[m]], runs = fit_grid(fitter))
    saveRDS(cache, CACHE)
  }
}
runs <- grid |>
  mutate(result = map2(cache$diablo$runs, cache$glmnet$runs, \(d, g) {
    list(diablo = d$observed, glmnet = g$observed, null = cbind(diablo = d$null, glmnet = g$null))
  }))

performance <- runs |>
  mutate(rows = map(result, \(r) {
    tibble(
      method = c("DIABLO / sPLS-DA", "elastic net (nestedcv)"),
      ber = c(r$diablo$ber, r$glmnet$ber),
      null_median = c(median(r$null[, "diablo"]), median(r$null[, "glmnet"])),
      p_value = c(
        mean(r$null[, "diablo"] <= r$diablo$ber), mean(r$null[, "glmnet"] <= r$glmnet$ber)
      ),
      empty_models = c(NA, r$glmnet$empty)
    )
  })) |>
  select(comparison, level, input, rows) |>
  unnest(rows)

gene <- set_names(str_split_i(de$genes$Genes, ";", 1), rownames(de$genes))
selected <- runs |>
  mutate(features = map(result, \(r) {
    sv <- mixOmics::selectVar(r$diablo$fit, comp = 1)
    sv <- if (!is.null(sv$value)) list(single = sv) else sv[names(r$diablo$fit$X)]
    diablo <- map(sv, \(v) {
      tibble(method = "DIABLO / sPLS-DA", feature = rownames(v$value), weight = v$value$value.var)
    }) |>
      list_rbind()
    folds <- r$glmnet$fit$outer_result
    glmnet <- map(folds, \(o) names(o$coef)[-1]) |>
      unlist() |>
      table() |>
      enframe("feature", "folds_selected") |>
      transmute(
        method = "elastic net (nestedcv)", feature = as.character(feature),
        weight = as.numeric(folds_selected) / length(folds)
      )
    bind_rows(diablo, glmnet)
  })) |>
  select(comparison, level, input, features) |>
  unnest(features) |>
  mutate(gene = gene[feature], .after = feature)

plot_index <- map(names(COMPARISONS), \(cmp) {
  r <- filter(runs, comparison == cmp)
  labels <- str_glue("{r$level}\n{r$input}")
  write_pdf(file.path(out$reports, str_glue("02_panels_{cmp}.pdf")), list(
    "balanced error rate against the 126-split null, both methods" = \() {
      par(mfrow = c(1, 2), mar = c(9, 4, 3, 1))
      iwalk(c(diablo = "DIABLO / sPLS-DA", glmnet = "elastic net (nestedcv)"), \(title, m) {
        stripchart(map(r$result, \(x) x$null[, m]),
          vertical = TRUE, method = "jitter", pch = 1,
          col = "grey65", cex = 0.4, group.names = labels, las = 2, cex.axis = 0.6,
          ylim = c(0, 1.1),
          ylab = "balanced error rate", main = str_glue("{cmp}: {title}"), cex.main = 0.85
        )
        points(seq_len(nrow(r)), map_dbl(r$result, \(x) x[[m]]$ber), pch = 19, col = "#B2182B")
      })
    }
  )) |>
    mutate(comparison = cmp, .before = 1)
}) |>
  list_rbind()

write_workbook(
  list(
    performance = performance, selected_features = selected, plot_index = plot_index,
    settings = tibble(
      setting = c(
        "keepX", "DIABLO design weight", "elastic net alpha", "lambda rule",
        "validation", "null"
      ),
      value = c(
        paste(names(KEEP_X), unlist(KEEP_X), sep = " = ", collapse = "; "), DESIGN_WEIGHT,
        0.5, "lambda.1se, inner cross-validation",
        "DIABLO: leave-one-out; elastic net: leave-pair-out (25 folds); balanced error rate",
        str_glue("all {ncol(splits)} five-and-five splits; p = share at or below observed")
      )
    )
  ),
  c(
    "Per comparison, level, input, method: BER, null median, exact p, empty elastic nets.",
    "Features used: DIABLO loadings (component 1); share of elastic-net folds selecting each.",
    "Which PDF page shows what.",
    "Model settings."
  ),
  out$workbook, inputs
)
performance |> filter(p_value <= ALPHA)
sessionInfo()
