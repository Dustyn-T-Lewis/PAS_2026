STYLE (shared by all seven schematics): A flat, clean scientific schematic in the look of a BioRender figure, landscape 16:9, on a pure white background. Shapes are rounded rectangles with soft pastel fills and thin outlines in a darker shade of the same hue; no gradients, shadows, 3D effects, textures or clip-art. All text is a clean sans-serif (Helvetica or Roboto style) in dark charcoal grey; titles bold, body text regular, set large enough to read when the image is shown at half size, with generous white space and boxes aligned on a tidy grid. Neutral boxes are very light grey with a mid-grey outline. Pipeline step cards are very pale blue-grey with a slate blue-grey outline. Female is orange (pale peach fill, burnt orange outline), male is blue (pale sky-blue fill, steel blue outline), and the treatment-by-sex interaction is purple (pale lilac fill, plum outline); these three colours are used for nothing else. Arrows are thin solid charcoal lines with small filled arrowheads. There are no icons, pictograms, animals, people, logos or decorations of any kind. Every text label is spelled exactly as written below, appears exactly once, and is fully legible; there is no other text anywhere in the image.

CONTENT: the whole PAS_2026 pipeline, read left to right.

Title at top left, bold: "PAS_2026 proteomics pipeline". Small dark grey rounded chip with white text in the top-right corner: "Rscript run_all.R".

Four columns of equal height, left to right. Each column is one large neutral rounded panel with its folder name in bold at the top; inside it a vertical stack of small rounded cards, one label per card. Columns 1 to 3 are each joined to the next by exactly one thick arrow from panel edge to panel edge.

Column 1 panel, header "00_Setup". Seven small neutral cards: "report.parquet", "report_stats.xlsx", "metadata.csv", "phenotypes.xlsx", "MitoCarta 3.0", "CORUM 5.3", "STRING v12". Below them, separated by a thin line, one card: "Shared_functions".

Column 2 panel, header "01_Processing". Six pale blue-grey step cards: "01 Filtering", "02 Quantification", "03 Design", "04 Gene sets", "05 Scores", "06 Modules". Under the cards, small text: "1,962 proteins, 595 in MitoCarta".

Column 3 panel, header "02_Statistics", with the small subtitle "run once on all 20 mice". Seven pale blue-grey step cards: "01 Fit", "02 Enrichment", "03 Classification", "04 PCA", "05 Composition", "06 Candidates", "07 Sensitivity".

Column 4 is split into two panels stacked one above the other, each reached by its own arrow from the column 3 panel (two arrows in total leave column 3).
- Upper panel, header "03_Comparisons", small subtitle "four steps written once in code/". Three cards side by side or stacked: an orange card with bold "Female" and below it "PAS_F − VEH_F"; a blue card with bold "Male" and below it "PAS_M − VEH_M"; a purple card with bold "Sex_Difference" and below it "(PAS_F − VEH_F) − (PAS_M − VEH_M)". Group names are written with plain underscores, never subscripts.
- Lower panel, header "04_Figures", small subtitle "drawn from saved results only; supplements filed under each figure". Six small neutral cards in two rows of three: "F01_Global", "F02_Proteins", "F03_Pathways", "F04_Phenotype_pathways", "F05_Complexes", "S1–S6".

A full-width very light grey bar along the bottom with the text: "every step folder:  a_script/   b_reports/   c_data/   README.md".
