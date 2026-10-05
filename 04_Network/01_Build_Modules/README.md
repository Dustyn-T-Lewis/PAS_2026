# 01_Build_Modules

Builds the co-abundance modules every later network step reads, using WGCNA's standard workflow on
the normalised full universe.

- **Method:** `WGCNA::pickSoftThreshold` (signed, `bicor`, maxPOutliers 0.05, scale-free R² cut 0.85; fallback power 16 if none reaches it), then `blockwiseModules` (signed TOM, minModuleSize 20, mergeCutHeight 0.25, deepSplit 2); eigengenes without grey; signed kME (`signedKME`, bicor); 10 hubs per module by kME.
- **Inputs:** `02_Differential_Expression/02_Contrasts/c_data/contrasts.rds`
- **Outputs:** `c_data/01_build_modules.xlsx` (sheets: sizes, modules, hubs, scale_free, eigengenes, plot_index), `c_data/modules.rds` (read by 02 to 05 and `05_Summary/02_Panels`), `b_reports/01_build_modules.pdf`
- **Result:** power 14 (first power with signed R² ≥ 0.85, R² 0.864). Nine modules plus grey: turquoise 520, blue 509, brown 170, yellow 160, green 112, red 77, black 56, pink 48, magenta 37, grey 153. Yellow is the most mitochondrial (64% MitoCarta), magenta the least (5%).
- **Note:** the script stops if the module set, turquoise size (520) or grey size (153) differ from the first build, so downstream results stay comparable.
