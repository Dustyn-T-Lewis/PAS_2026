# 05 · Summary

Pulls the earlier stages together and asks the questions no single level answers: how much signal
each family of tests holds against chance, whether feature combinations separate groups in unseen
mice, which responses are shared or differ by sex, how far PAS moves the whole proteome, and whether
each sex's PAS signature carries over to the other sex or tracks phenotype. Steps 01 and 04 only read
finished workbooks; the rest fit their own models. Results are nominal, with statistics for
filtering; the working rule is p ≤ 0.05 and FDR ≤ 0.05.

| Sub-stage | Does | Writes |
|---|---|---|
| `01_Signal` | nominal hits against chance, Π, FDR tiers and true-effect share for every test family | workbook, PDF |
| `02_Panels` | DIABLO and elastic-net panels scored in held-out mice, exact permutation p | `panels_cache.rds`, workbook, 4 PDFs |
| `03_Biclusters` | plaid biclusters over 100 seeds and how often each recurs | workbook, PDF |
| `04_Sex_Patterns` | shared, detected-in-one-sex or sex-differential label per protein, pathway, module | workbook, PDF |
| `05_Global` | PCA, PERMANOVA, dispersion, envfit of the MitoCarta share | `global.rds`, workbook |
| `06_Magnitude` | each mouse's distance from its sex's vehicle centroid, tested by RRPP | `magnitude.rds`, workbook |
| `07_Signature` | each sex's PAS signature tested in the other sex with `limma::fry` | `signature.rds`, workbook |
| `08_Signature_Phenotype` | PAS effect on each phenotype; singscore signature score against phenotype | `signature_phenotype.rds`, workbook |

Every workbook opens with a `read_me` sheet (what each sheet holds) and ends with `input_manifest`
(files read, with md5) and `package_versions`. `05` to `08` write no PDF; their rds files feed
`06_Figures`.

Run order: `01` to `08` in number order, after stages 02 to 04 are complete (`01` reads every
contrast and associate workbook). `08` needs `07`'s `signature.rds`; the others are independent.
