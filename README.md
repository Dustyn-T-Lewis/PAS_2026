# PAS_2026

Gastrocnemius proteomics of pipecolic acid plus succinate (PAS) against vehicle in mice with a ligated femoral artery. 20 mice, five per sex and treatment, one mitochondria-enriched pellet each, DIA mass spectrometry; 1,962 proteins, 595 in MitoCarta. Six phenotypes per mouse: CSA, capillary density and satellite cells by histology, and respiratory conductance in soleus and plantaris.

The project asks three questions:

1. Does PAS change females? PAS_F − VEH_F.
2. Does PAS change males? PAS_M − VEH_M.
3. Does the PAS response differ by sex? (PAS_F − VEH_F) − (PAS_M − VEH_M).

Differences between females and males at baseline are not a question here. They dominate the proteome (PERMANOVA R² 0.39), so they appear as context in `02_Statistics/04_PCA` and [Figure 1](04_Figures/F01_Global/README.md).

![The pipeline: inputs, processing, statistics run once on all 20 mice, the three comparisons and the figures](pipeline.png)

## Folders, in run order

| Folder | Holds |
|---|---|
| [`00_Setup`](00_Setup/README.md) | the study data and reference files (`Input/`) and the helpers several scripts source (`Shared_functions/`) |
| [`01_Processing`](01_Processing/README.md) | filtering, quantification, design, gene sets, set scores, modules |
| [`02_Statistics`](02_Statistics/README.md) | every model and test, run once on all 20 mice: the fit, pathway tests, classification, PCA, pellet composition, candidates, sensitivity |
| [`03_Comparisons`](03_Comparisons/README.md) | the three questions: four steps written once in `code/`, run for `Female`, `Male` and `Sex_Difference` |
| [`04_Figures`](04_Figures/README.md) | figures drawn from saved results only |

Data flows one way: processing feeds statistics, each comparison slices the statistics for its contrast, and the figures read the statistics. No comparison reads another. Every step reads only files earlier steps saved, so it reruns alone once those exist.

## Every step

```
01_Processing/NN_Step/, 02_Statistics/NN_Step/
├── a_script/   one R script
├── b_reports/  one preview PDF from package plots, where there is something to draw
├── c_data/     one workbook, plus any .rds a later step reads
└── README.md   what the step reads, writes and finds

04_Figures/F0N_Name/
├── a_script/   F0N.R (the composite) and one script per panel; supp/ holds S#.R and S#a_….R
├── b_reports/  the composite and each panel as PDF and PNG; supp/ the same for the supplements
├── c_data/     F0N.xlsx and S#.xlsx, one sheet per panel with the table it draws
└── README.md   the main figure and its supplements, with panel tables

03_Comparisons/
├── code/       01_proteins.R … 04_classification.R, each step written once
└── <comparison>/
    ├── a_script/  run_comparison.R: runs the steps for this contrast
    ├── b_reports/ one preview PDF per step
    ├── c_data/    one workbook per step, named by step
    └── README.md  the comparison's results
```

Every step workbook opens on `read_me`, which says what each sheet holds, and ends with `package_versions`, the packages attached at run time; a figure workbook holds only the panel sheets. Every step PDF is a short preview drawn with the package's own plots (limma, limpa, WGCNA, fgsea, pROC, vegan) or base graphics; the styled figures are in `04_Figures`. Which protein belongs to which pathway, complex or cluster is in `01_Processing/04_Gene_Sets/c_data/04_gene_sets.xlsx`, sheet `membership`.

## Conventions

- **Model:** one limpa/limma fit on all 20 mice: four group means and one surrogate variable from `sva`, which tracks run depth (r = 0.85) and does not differ by treatment. Proteins use limpa `dpcDE` with precision and sample weights, each contrast read as its own coefficient; set scores use limma with eBayes; the eight module eigengenes use ordinary least-squares t, since eBayes cannot estimate a prior from eight features.
- **Normalisation:** chosen per universe by replicate spread (NormalyzerDE measures) among methods that keep limpa's standard errors valid: none for the full proteome, median centring for MitoCarta.
- **Multiple testing:** BH within universe (full or mitochondrial) and contrast for proteins; within contrast and collection for fgsea; over every scored set in the universe for set scores. p ≤ 0.05 counts are reported as exploratory.
- **Π** = p^|log2 FC| (Xiao et al. 2014). It ranks on p and fold change together, controls no error rate, and is counted on its own. TREAT (fold change above 1.2) and ashr's local false sign rate sit beside it as error-controlled checks.
- **r** sits beside every test: r with group (PAS = 1) for a contrast; r across the sex's ten mice, and within vehicle and PAS, for an association.
- **Phenotype effect** means what PAS did to a phenotype, as % change against the same sex's vehicle. **Phenotype association** means a feature moving with a phenotype within one sex, reported pooled across PAS and vehicle mice and within treatment; it is co-movement, not a PAS effect.
- **Pathways:** fgsea through enrichVolcano on GO:BP, Reactome, Hallmark and KEGG LEGACY (full proteome) and MitoCarta (mitochondrial universe). The same sets are scored per mouse with singscore, and that mouse-level test sits beside every fgsea row. camera, which estimates each set's inter-protein correlation, sits beside every contrast fgsea row too. fgsea permutes genes, so where they disagree the camera and singscore results are the conservative ones.
- **Seed:** 20260910 for every stochastic step.

## Run

```sh
Rscript -e 'renv::restore()'
Rscript run_all.R
```

About 7.5 minutes on a laptop from a fresh clone. R 4.6.0. `renv.lock` pins every package; enrichVolcano installs from github.com/Dustyn-T-Lewis/enrichVolcano. CORUM and the STRING network are read from `00_Setup/Input/`, but it needs internet: STRINGdb checks string-db.org when it starts, and the STRING connectivity test in 02_Enrichment runs on its server. limpa's quantification is checkpointed in the repository at each detection slope; `01_Processing/02_Quantification/README.md` says how to refit it.

## Open with Dr. Ismaeel

Whether the measured limb is the ligated one; the muscle, marker and units of the histology; the capillary density definition; and the respiratory conductance units (`00_Setup/Input/phenotypes.xlsx`, sheet `notes`). The proteome comes from gastrocnemius, while respirometry comes from soleus and plantaris.

## Version 1

The exploratory first version is on branch `Dustyn/v1-exploratory-analyses`. Its protein contrasts are identical to these.

## Reference

Burke BI, Valentino TR, Ismaeel A, et al. Exercise-associated microbial metabolites prevent skeletal muscle atrophy in adult female mice. Nat Commun 2026. PMID 42431880.
