# PAS_2026

Gastrocnemius proteomics of pipecolic acid plus succinate (PAS) against vehicle in mice with a ligated femoral artery: 20 mice, five per sex and treatment, one mitochondria-enriched pellet each, DIA mass spectrometry, 1,962 proteins (595 in MitoCarta). Six phenotypes per mouse: CSA, capillary density and satellite cells by histology, and three respiratory conductances from soleus and plantaris. The project asks three questions: does PAS change females (PAS_F − VEH_F), does PAS change males (PAS_M − VEH_M), and does the PAS response differ by sex ((PAS_F − VEH_F) − (PAS_M − VEH_M))?

![The pipeline: inputs, processing, statistics run once on all 20 mice, the three comparisons and the figures](pipeline.png)

| Folder | Holds |
|---|---|
| [`00_Setup`](00_Setup/README.md) | input data, reference files, shared R helpers, schematic specs |
| [`01_Processing`](01_Processing/README.md) | precursors to protein values, the design, gene sets, set scores, modules |
| [`02_Statistics`](02_Statistics/README.md) | every model and test, run once on all 20 mice |
| [`03_Comparisons`](03_Comparisons/README.md) | the three questions, one folder each |
| [`04_Figures`](04_Figures/README.md) | Figures 1 to 5 and Supplementary Figures 1 to 6 |

## Conventions

- Each step folder holds `a_script/` (the script), `b_reports/` (a preview PDF or the figure files) and `c_data/` (one workbook plus any `.rds` a later step reads). Step workbooks open on `read_me` and end with `package_versions`.
- Data flows one way: processing, statistics, comparisons, figures. Each step reads only saved outputs of earlier steps.
- One limpa/limma fit on all 20 mice: four group means plus one `sva` surrogate variable, which tracks run depth (r = 0.85) and not treatment.
- Normalisation: none for the full proteome, median centring for MitoCarta, chosen by replicate spread.
- BH correction within universe and contrast for proteins, within contrast and collection for fgsea. Counts at p ≤ 0.05 are exploratory.
- Π = p^|log2 FC| (Xiao et al. 2014) ranks proteins and controls no error rate.
- Phenotype effect: what PAS did to a phenotype, as % change against the same sex's vehicle. Phenotype association: a feature moving with a phenotype within one sex. An association is co-movement, not a PAS effect.
- Baseline female against male differences are context, not a question. They dominate the proteome (PERMANOVA R² 0.39).
- Seed 20260910 for every stochastic step.

## Run

```sh
Rscript -e 'renv::restore()'
Rscript run_all.R
```

About 7 minutes on a laptop, 12 from a fresh clone (limpa refits its three cached quantifications). R 4.6.0; `renv.lock` pins every package, and enrichVolcano installs from github.com/Dustyn-T-Lewis/enrichVolcano. The run needs internet because STRINGdb contacts string-db.org.

## Version 1

The exploratory first version is on branch `Dustyn/v1-exploratory-analyses`. Its protein contrasts are identical to these.

## Reference

Burke BI, Valentino TR, Ismaeel A, et al. Exercise-associated microbial metabolites prevent skeletal muscle atrophy in adult female mice. Nat Commun 2026. PMID 42431880.
