# 01_Filtering

Turns DIA-NN precursors into a filtered, contaminant-free precursor matrix for quantification.
No mouse is dropped and nothing is filtered on partial missingness, because `dpcQuant()` fits its
detection curve to those gaps.

- **Method:** `limpa::readDIANN()` with five q-value columns (`Q.Value`, `Lib.Q.Value`, `Lib.PG.Q.Value`, `Global.Q.Value`, `Global.PG.Q.Value`) at 0.01, the union of limpa's with- and without-MBR recipes; `removeNARows()` for all-missing precursors; removal of cRAP-only groups; `filterNonProteotypicPeptides()` then `filterCompoundProteins()`; removal of groups whose gene symbols all match an anchored keratin, immunoglobulin and haemoglobin pattern. `filterSingletonPeptides()` is not applied.
- **Inputs:** `00_Input/report.parquet`, `00_Input/metadata.csv`, `00_Input/report_stats.xlsx`
- **Outputs:** `c_data/01_filter.xlsx` (sheets: filter_log, contaminants_removed, crap_detected, contaminant_share, contaminant_per_run, run_qc), `c_data/precursors_filtered.rds` (read by `02_Quantification`), `b_reports/01_filter_figures.pdf` (raw signal against protein groups identified, per run)
- **Result:** 24,801 precursors after the q-value filter down to 20,355 in 1,923 protein groups, holding 74.2% of signal. Most of the loss is shared peptides (98.3% to 74.9%). 14 cRAP-only groups and 8 contaminant groups (0.003% of signal) were removed. Median plasma-protein share per group is 0.71 to 0.97%. Raw run signal spans 0.40 to 2.35 times the cohort median; no run is dropped.
- **Note:** The script stops if any cRAP group shares an accession with a mouse group, or if any of 12 named muscle proteins (`Ckm`, `Cycs`, `Mb`, `Sod1` and others) is lost, which catches an over-broad contaminant pattern. Albumin, fibrinogen and complement are kept, because the ischemic limb may hold extravasated plasma protein; `contaminant_share` reports that share by group. Run QC uses uncorrected `Precursor.Quantity`.
