# SwDGVirome: viral communities in Baltic Sea water and a deep granitic aquifer

## Overview

This repository contains the analysis notebooks for four samples from Baltic Sea water and the Äspö Hard Rock Laboratory (Sweden). The discovery catalogue contains 2,488 vOTUs; the primary catalogue contains 962 representatives. Comparisons among the four unreplicated samples are descriptive: sampling depth covaries with other environmental properties and is not isolated as a causal factor.

Predicted temperate fractions in Figures 2 and 3 are calculated as T/(T + V), using the summed TPM of predicted temperate and virulent vOTUs within the relevant sample or abundance group. Figure S2 instead reports count-based fractions for the discovery catalogue. Table S3 summarises quality-stratified results and taxonomic assignment coverage for the primary catalogue, with 962 rows of vOTU detail including sample-specific abundance groups. Fractions are displayed on a 0–1 scale.

## Installation

Requirements: R >= 4.2, Quarto CLI and its Pandoc installation, and the R dependencies declared in `DESCRIPTION`. The current workflow has been exercised with R 4.5. The notebooks use fonts including Arial and Times; local font availability can affect rendering.

```bash
git clone --branch dev_pan https://github.com/apan1027/SwDGVirome.git
cd SwDGVirome
```

Install dependencies from the repository root in R:

```r
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::pak()
```

## Inputs and execution

The repository tracks selected required inputs; most generated outputs and local archives are ignored.

- `analyses/data/00-raw/d00-resource/p0057v2.sqlite` and `imgvr_source.tsv`: input database and IMG/VR metadata for step 01.
- `analyses/data/00-raw/d11-figs8-public-reference/`: the two archived CSV inputs for Figure S8 and Table S2.
- `analyses/data/16-amg-targeted-validation/`: three archived inputs used by step 10. This historical directory is still required; there is no additional step 16 to run.

Run step 01 before the analyses that use its TSE objects. Step 09 additionally requires the output of step 04. Step 04 queries KEGG online, so successful annotation requires network access to KEGG. Check the query logs as well as render completion: failed queries can produce missing annotations or fallback classifications.

```bash
quarto render analyses/01-tse-construction.qmd
quarto render "analyses/02-fig2-viral diversity.qmd"
```

See [analyses/README.md](analyses/README.md) for the complete 01–11 notebook list and supplementary-table mapping. For GitHub-readable reports, render the relevant notebook with `--to gfm` and include its referenced preview images.

## Sample groups

| Code | Description | Sampling depth |
|------|-------------|----------------|
| BS | Baltic Sea | 0 m |
| SA | Shallow aquifer | 71 m |
| IA | Intermediate aquifer | 196 m |
| DA | Deep aquifer | 450 m |

## Outputs and reproducibility scope

Generated files are placed in `analyses/data/<step-name>/`. Selected Markdown reports, preview images and source tables are tracked for inspection. Figure 2b uses random rarefaction without a fixed seed, so rerun curves may vary.

Step 10 redraws genomic context and presents archived structure-comparison evidence; it does not rerun structure prediction or alignment. Step 11 redraws published coverage records and exports an archived reference-match table; it does not rerun mapping or sequence alignment. Step 06 generates the Table S3 quality and taxonomic-coverage CSV files and exports `Table_S3_quality_taxonomy.xlsx`, including the 962-row vOTU detail table. The selected CSVs are tracked; the generated Excel workbook remains a local output.

## Contact

Project maintainer: Cunli Pan (cunli.pan@tum.de).
