# Depth Stratifies Viral Survival Strategies in a Deep Granitic Aquifer

## Overview
This repository contains the analytical workflow for investigating depth-associated stratification of viral communities (0–450 m) in deep groundwater at the Äspö Hard Rock Laboratory (Sweden).

## Main results (summary)
- 2,488 viral operational taxonomic units (vOTUs) were identified across four depth zones.
- Viral lifestyle metrics suggest increasing prevalence of lysogeny with depth.
- Auxiliary metabolic genes (AMGs) show depth-stratified patterns consistent with environmental and host metabolic differences.

## Repository structure
- `analyses/`: Quarto analysis notebooks used to generate figures/tables.
- `R/`: Helper functions (e.g., `R/hello.R`).
- `analyses/data/00-raw/d00-resource/`: **Input data (not tracked)**; see “Input data” below.
- `renv.lock`, `renv/`: Reproducible R environment metadata.

## Workflow overview
Run notebooks in numeric order.

1. `01-tse-construction.qmd` — Build the `TreeSummarizedExperiment` object from the SQLite database (assays: counts/TPM/RPKM etc.; annotations stored in `rowData`, sample metadata in `colData`, auxiliary tables in `metadata`).
2. `02-fig2-viral-diversity.qmd` — Viral diversity analyses and Figure 2 panels (composition, rarefaction, Venn, lifestyle ratios, diversity indices, source/novelty, host phyla).
3. `03-fig3-core-virome.qmd` — Core virome definition and Figure 3 panels.
4. `04-fig4-abundant-rare-viral.qmd` — Abundant vs rare stratification and Figure 4 panels.
5. `05-fig5-heatmap.qmd` — AMG/KEGG pathway summaries and Figure 5 heatmaps.
6. `06-figs1-host-composition.qmd` through `11-figS6-amg-DA.qmd` — Supplementary figures (S1–S6).


## Sample Groups

| Code | Description | Depth Range |
|------|-------------|-------------|
| BS | Baltic Sea | Surface (0 m) |
| SA | Shallow Aquifer | 71 m |
| IA | Intermediate Aquifer | 196 m |
| DA | Deep Aquifer | 450 m |


## Dependencies
This project uses `renv` for reproducible R package management (see `renv.lock`).

## Installation
### Option 1: Use renv (recommended for exact reproducibility)
In R (from the project root):
```r
# Install renv
install.packages("renv")

# Restore exact package versions from renv.lock
renv::restore()
```
### Option 2: Manual installation (if renv fails）
```r
# Bioconductor packages
BiocManager::install(c("TreeSummarizedExperiment", "mia"))

# CRAN packages
install.packages(c("tidyverse", "vegan", "ggplot2", "pheatmap", 
                   "VennDiagram", "ComplexHeatmap", "openxlsx"))
```
**Troubleshooting:**
- If `renv::restore()` fails, try Option 2
- Check R version compatibility (see `renv.lock` for required R version)

### Key Packages
Core dependencies (complete list in `renv.lock`):
- **Data structures**: TreeSummarizedExperiment, mia, SummarizedExperiment
- **Database**: DBI, RSQLite
- **Tidyverse**: tidyverse, dplyr, tidyr, ggplot2
- **Visualization**: ComplexHeatmap, pheatmap, VernDiagram, patchwork
- **Analysis**: vegan, KEGGREST

### System Requirements
- R ≥ 4.0 (see `renv.lock` for exact version)
- Quarto CLI (for rendering `.qmd` files)
- 16+ GB RAM recommended

## Usage

### Quick Start
1. Restore packages: `renv::restore()` (see Installation)
2. Place input data files (see Input Data section)
3. Render analyses in order (see below)

### Step 1: Construct TSE object (required first)
````bash
quarto render analyses/01-tse-construction.qmd
````
This creates `data/01-tse-construction/tse.rds` used by all subsequent analyses.

### 2. Run main analyses (Fig. 2–5)
```bash
quarto render analyses/02-fig2-viral-diversity.qmd
quarto render analyses/03-fig3-core-virome.qmd
quarto render analyses/04-fig4-abundant-rare-viral.qmd
quarto render analyses/05-fig5-heatmap.qmd
```

### 3. Generate supplementary figures (Fig. S1–S6)
```bash
quarto render analyses/06-figs1-host-composition.qmd
quarto render analyses/07-figs2-venn-diagram.qmd
quarto render analyses/08-figs3-lifestyle.qmd
quarto render analyses/09-figs4-abundant-patterns.qmd
quarto render analyses/10-figs5-rare-patterns.qmd
quarto render analyses/11-figS6-amg-DA.qmd
```
**Note:** Scripts 02-11 all require `tse.rds` from Step 1.

## Input data requirements
Place the following files in:
analyses/data/00-raw/d00-resource/
p0057v2.sqlite – main SQLite database

imgvr_source.tsv – IMG/VR ecosystem source metadata

**Note:** The exact path is managed by `projthis` package. If your actual folder name differs, adjust the path above accordingly.


## Output Files

Each analysis generates outputs in `data/<script-name>/`:

**Example for `02-fig2-viral-diversity.qmd`:**
````
data/02-fig2-viral-diversity/
├── Fig2a.png                   # Family composition (300 DPI)
├── Fig2a_data.csv              # Source data
├── Fig2b_combined.png          # Rarefaction curves
├── ...
└── session_info.txt            # R session details
````

**Output formats:**
- Figures: PNG (300 DPI) and/or PDF
- Data tables: CSV or XLSX
- Session info: Complete package versions

**Note:** Exact output path printed in console when rendering.


## Common Issues

### Error: `file.exists(sqlite_path) is not TRUE`
**Cause:** Input SQLite database not found.

**Solution:**
1. Verify `p0057v2.sqlite` is in the resource directory
2. Check console output for expected path
3. Ensure file permissions allow reading

### Error: `cannot open file 'tse.rds'`
**Cause:** TSE object not built yet.

**Solution:** Run `01-tse-construction.qmd` first.

### Error: `Package 'XXX' not found`
**Cause:** Missing dependencies.

**Solution:**
````r
# Try renv restore
renv::restore()

# Or install manually
BiocManager::install("XXX")
````









