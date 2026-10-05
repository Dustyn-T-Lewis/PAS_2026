# PAS_2026

Skeletal muscle proteomics of pipecolic acid plus succinate (PAS), a pair of exercise-associated
microbial metabolites, against vehicle in mice with a ligated femoral artery. Twenty mice, five per
sex and treatment, one gastrocnemius pellet each (mitochondria-enriched), DIA mass spectrometry,
plus six muscle phenotypes. The question: what does PAS do to the muscle proteome, and does it do
it differently in females and males?

1,842 proteins after filtering (568 in MitoCarta). Five contrasts: PAS vs vehicle in females and in
males, the treatment-by-sex interaction, and sex in each treatment.

## Results in one paragraph

Sex dominates the proteome (PERMANOVA R² 0.38, p = 0.001); PAS does not move it as a whole
(R² 0.04, p = 0.27). PAS moves 152 female and 133 male proteins at p ≤ 0.05 against 92 expected by
chance; 4 female proteins pass FDR ≤ 0.05, none in males. In females PAS raises muscle contraction,
myogenesis and myofibril pathways and lowers cell cycle; in males it lowers OXPHOS, EMT and actin
pathways. The two responses are independent (signature carry-over, fry p > 0.5), and contractile
proteins move up in PAS females and down in PAS males. PAS raises capillary density, CSA and
satellite cells in females; the protein signature marks treated females but does not track those
traits within groups.

## Layout

```
00_Input/                    data; runs nothing
01_Preprocess/               filtering, limpa quantification
02_Differential_Expression/  design, contrasts, classification, phenotype association
03_Pathway_Enrichment/       gene sets, scores, contrasts, classification, association, plots,
                             compartments, figure term selection
04_Network/                  WGCNA modules, contrasts, classification, association, STRING, INDRA
05_Summary/                  signal, panels, biclusters, sex patterns, global structure, magnitude,
                             signature carry-over, signature against phenotype
06_Figures/                  F01-F05 with supplements S1-S4
R/                           helpers.R (analysis), figures.R (style, shared panels), panels.R
Supplementary/               submission files, once drafted
run_all.R                    every step in order, each in its own R session
```

Every stage and sub-stage has a `README.md`, `a_script/` (one script), `b_reports/` (PDFs) and
`c_data/` (one workbook, plus `.rds` files a later step reads). Stages 01 to 05 are an exhaustive
exploration: every workbook holds nominal results with p, FDR and Π side by side and a contrast
or comparison column on every row, so any table filters in Excel. Each workbook opens with
`read_me` and ends with `input_manifest` (md5 of inputs) and `package_versions`.

Figures follow one pattern: `a_script/main/F0x.R` stitches the panel scripts in
`a_script/main/panels/`; supplements sit in `a_script/supp/` the same way; renders mirror this
under `b_reports/`; `c_data/` holds one workbook of panel data per figure.

## Conventions

- Significant: p ≤ 0.05 or FDR ≤ 0.05 (BH within contrast and collection); proteins also at
  FDR ≤ 0.10. Π = p^|log2 FC| (Xiao 2014) ranks proteins and controls no error rate.
- Sex differences are claimed only from the treatment-by-sex contrast. A hit in one sex alone is
  "detected in" that sex.
- Phenotypes come from soleus, plantaris and IHC; the proteome from gastrocnemius. Links between
  them are associations.

## Running

```sh
Rscript -e 'renv::restore()'   # once: installs the pinned packages
Rscript run_all.R
```

R 4.6.0. Main packages: limpa 1.4.0, limma 3.68.1, enrichVolcano 2.3.0
(github.com/Dustyn-T-Lewis/enrichVolcano), fgsea 1.38.0, WGCNA 1.74, mixOmics 6.36.0,
nestedcv 0.9.0, vegan 2.7.3, RRPP 2.2.0, singscore 1.32.0, STRINGdb 2.24.0, ggplot2 4.0.3,
patchwork 1.3.2. `renv.lock` pins the rest. `dpcQuant()` is checkpointed; STRING, INDRA and the
classification panels cache their results keyed on their inputs (first run of `02_Panels` takes
about 25 minutes on 12 cores).

## References

Burke BI, Valentino TR, Ismaeel A, et al. Exercise-associated microbial metabolites prevent
skeletal muscle atrophy in adult female mice. Nat Commun 2026. PMID 42431880.
