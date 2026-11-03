# F02_Proteins

- Figure 2: within each sex, which proteins are the strongest follow-up candidates, how they group by known interactions, and how every mouse looks for the strongest.
- Supplementary Figure 3: each comparison's Π ≤ 0.05 proteins across all 20 mice, cut into co-expression clusters, and the volcano plots.
- Supplementary Figure 4: how the protein results move when one processing or model choice changes.

Candidates and their six routes are defined in [`02_Statistics/06_Candidates`](../../02_Statistics/06_Candidates/README.md). They are discovery-stage candidates, not validated biomarkers.

## Figure 2

![Figure 2](b_reports/F02.png)

### a. Female candidate network

![Panel a](b_reports/a_female_network.png)

The 58 female candidates passing two or more routes; edges are STRING v12 links at combined score ≥ 700, colour is log2 fold change and size the number of routes. 35 have no STRING partner; the other 23 form eight communities, the largest named for TGF-β response regulation (Hsp90ab1, Cdc37, Cops3, Cops8, Sugt1, Fkbp1a) and skeletal muscle contraction (Atp2a1, Cacna1s, Srl, Casq1). A community's name may describe only some of its members.

### b. Male candidate network

![Panel b](b_reports/b_male_network.png)

The same for males: 47 protein groups on 46 nodes, since two Cav1 groups share one gene. 20 have no partner; the rest form seven communities, among them muscle contraction, EMT, nervous system development and heme metabolism. Many members are nerve, red-cell and extracellular matrix proteins, classes that fall in male PAS pellets ([05_Composition](../../02_Statistics/05_Composition/README.md)), so much of this network may record pellet content rather than a change in the fibres.

### c. Female candidates, every mouse

![Panel c](b_reports/c_female_candidates.png)

The 14 female proteins passing three or more routes, each mouse's log2 abundance minus its own-sex vehicle mean; an open circle marks a mouse where the protein was not detected. Nine fall and five rise (p 2.1e-05 to 0.032); only Dusp3 reaches FDR ≤ 0.05 (0.022). The routes chose these proteins from these same mice, so the female separation shows consistency, not confirmation.

### d. Male candidates, every mouse

![Panel d](b_reports/d_male_candidates.png)

The same for the 11 male protein groups passing three or more routes (Cav1 twice, labelled by accession). Nine fall and two rise (Mdh1, Comtd1); none reaches FDR ≤ 0.05 (smallest Cd36, 0.094). The list mixes myofibril, extracellular matrix and red-cell proteins, so the composition caveat from panel b applies.

## Supplementary Figure 3

![Supplementary Figure 3](b_reports/supp/S3.png)

The heatmaps show log2 abundance minus the own-sex vehicle mean (green above, magenta below, grey not detected), rows split by direction and cut into clusters as described in [`02_Enrichment`](../../02_Statistics/02_Enrichment/README.md). Bars beside each cluster are its three smallest-p ORA sets at p ≤ 0.05, on one −log10 p axis for panels a to c. The clusters come from proteins chosen for responding, so the ORA describes them and tests little.

### a. Female PAS against vehicle

![Panel a](b_reports/supp/S3a_pi_female.png)

The 107 female Π ≤ 0.05 proteins: 44 up in one cluster (C1, 43) and 63 down in C2 (51) and C3 (10), with three unassigned. No female set reaches FDR ≤ 0.05; C2's best is Reactome cell cycle (p 0.001, FDR 0.17). The female PAS columns separate by construction; the male columns are the informative part.

### b. Male PAS against vehicle

![Panel b](b_reports/supp/S3b_pi_male.png)

The 106 male Π ≤ 0.05 proteins: 27 up (C1, C2) and 79 down (C3, C4), with seven unassigned. C3 holds the muscle sets (Reactome striated muscle contraction FDR 5e-7). Two female mice, VEH_F5 and PAS_F2, show a solid green stripe here and in panel c from high fast-muscle sarcomere proteins; they score highest in their groups on the myofibril marker set.

### c. Treatment × sex

![Panel c](b_reports/supp/S3c_pi_interaction.png)

The 142 interaction Π ≤ 0.05 proteins: 97 with a larger female effect (C1 to C3) and 45 with a larger male effect (C4, C5). C3 holds Hallmark myogenesis (FDR 6e-9) and C2 Hallmark EMT (FDR 0.002). A positive difference can come from a female rise, a male fall or both, and the columns show which.

### d. Volcano plots

![Panel d](b_reports/supp/S3d_volcanoes.png)

Every protein per contrast, coloured where Π ≤ 0.05, with p = 0.05 dashed and FDR ≤ 0.05 ringed and named. At p ≤ 0.05: 177 female, 193 male, 162 interaction. Only Cops3, Dusp3 and Cops8 pass FDR ≤ 0.05, all female and down.

## Supplementary Figure 4

![Supplementary Figure 4](b_reports/supp/S4.png)

The variants are listed in [`02_Statistics/07_Sensitivity`](../../02_Statistics/07_Sensitivity/README.md).

### a. Protein results under each variant

![Panel a](b_reports/supp/S4a_counts.png)

Per variant and contrast: proteins at p ≤ 0.05, FDR ≤ 0.10 and FDR ≤ 0.05, and how many of the baseline's 50 smallest p stay in the top 50. Grey shade is the absolute log2 ratio of the count to baseline (one added to both), darker for further from baseline in either direction. Every variant finds 132 to 212 proteins at p ≤ 0.05; the female FDR ≤ 0.05 count runs from 1 to 13, and no variant gives a male protein at FDR ≤ 0.05.

### b. Baseline proteins under each variant

![Panel b](b_reports/supp/S4b_named.png)

The 22 baseline proteins at FDR ≤ 0.10 (18 female, 4 male), shaded by their FDR under each variant, with the log2 FC as text. No protein changes sign. Cops3, Dusp3 and Cops8 stay at FDR ≤ 0.10 under all seven variants; the four male proteins start near the line (FDR 0.094 to 0.096) and Slc2a4 holds under none.

## Run

`Rscript 04_Figures/F02_Proteins/a_script/F02.R`, `a_script/supp/S3.R` and `a_script/supp/S4.R`. Each panel script also runs alone.

## Workbook

`c_data/F02.xlsx`, `S3.xlsx` and `S4.xlsx` hold one sheet per panel with the table it draws. Each S3 heatmap has two sheets: `<panel>_heatmap` (proteins in drawn order with cluster) and `<panel>_ora` (the bars).
