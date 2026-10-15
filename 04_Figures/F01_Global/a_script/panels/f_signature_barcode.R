# F01 f: each sex's PAS signature (proteins at p <= 0.05, up and down) placed in the other sex's
# ranking, as limma's barcode; each lane carries the fry test across mice (02_Enrichment).
source(here::here("00_Setup", "Shared_functions", "helpers.R"))

inputs <- c(enrichment = output_path("02_Statistics", "02_Enrichment", "enrichment.rds"))
enr <- readRDS(inputs[["enrichment"]])
ranks <- list(male = enr$ranks$PAS_vs_VEH_male$full, female = enr$ranks$PAS_vs_VEH_female$full)
stopifnot(nrow(enr$signatures) == 4)

plot <- signature_barcode(enr$signatures, ranks) +
  labs(
    x = sprintf(
      "rank in the other sex (n = %s)", format(length(ranks$male), big.mark = ",")
    ),
    y = NULL, title = "Cross-sex signatures",
    subtitle = "Each sex's p ≤ 0.05 proteins, ranked in the other sex"
  )

save_panel(plot, "F01_Global", "f_signature_barcode", 60, 62)
list(plot = plot, data = select(enr$signatures, -members), inputs = inputs)
