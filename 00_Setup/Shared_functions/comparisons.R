# Helpers the four steps in 03_Comparisons/code/ share: where a comparison writes and which mice
# it looks at. Sourced by every step; it sources helpers.R.

source(here::here("00_Setup", "Shared_functions", "helpers.R"))

# A comparison writes every step into one folder: b_reports/ and c_data/, files named by step.
comparison_paths <- function(comparison, step) {
  root <- here("03_Comparisons", COMPARISONS[[comparison]]$dir)
  walk(file.path(root, c("b_reports", "c_data")), dir.create, showWarnings = FALSE)
  list(
    root = root, data = file.path(root, "c_data"), reports = file.path(root, "b_reports"),
    workbook = file.path(root, "c_data", paste0(step, ".xlsx"))
  )
}

# The mice a comparison looks at: one sex for the PAS comparisons, all twenty for the sex
# difference.
comparison_mice <- function(comparison, targets) {
  sex_code <- COMPARISONS[[comparison]]$sex
  if (is.na(sex_code)) targets$sample_id else targets$sample_id[targets$sex == sex_code]
}
