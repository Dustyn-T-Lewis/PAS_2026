# 00_Setup

Files every script needs before it starts. Nothing here runs on its own.

| Folder | Holds |
|---|---|
| [`Input`](Input/README.md) | the DIA-NN report, sample sheet, phenotypes, MitoCarta 3.0, CORUM 5.3 and STRING v12 |
| [`Shared_functions`](Shared_functions/README.md) | R code sourced by more than one script |
| [`Schematics`](Schematics/) | the specs behind the seven schematic images |

`Schematics/` has one folder per image: `pipeline` (`pipeline.png`), `processing` (`01_Processing/schematic.png`), `statistics` (`02_Statistics/schematic.png`), `comparisons` (`03_Comparisons/schematic.png`), and `female`, `male` and `sex_difference` (each comparison folder's `schematic.png`). Each holds `SPEC.md` (every label drawn), `EVIDENCE.md` (the source of each element) and `DONOTDRAW.md`.
