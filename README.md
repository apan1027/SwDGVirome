# Depth Stratifies Viral Survival Strategies in a Deep Granitic Aquifer

## Overview
This repository contains the analytical workflow for investigating depth-associated stratification of viral communities (0–450 m) in deep groundwater at the Äspö Hard Rock Laboratory (Sweden).

## Key Results
- **Viral Diversity**: 2,488 viral operational taxonomic units (vOTUs) identified across four depth zones.
- **Lifestyle Stratification**: Viral lifestyle metrics suggest increasing prevalence of lysogeny with depth.
- **Metabolic Potential**: Auxiliary metabolic genes (AMGs) show depth-stratified patterns consistent with environmental and host metabolic differences.

## Repository Structure
- `analyses/`: Quarto analysis notebooks used to generate figures/tables.
- `R/`: Helper functions and internal package code.
- `analyses/data/`: Processed data storage (not tracked by git).
- `analyses/data/00-raw/d00-resource/`: Input data directory (see "Input Data" below).
- `DESCRIPTION`: Project metadata and complete dependency list.

## Installation

This project is structured as an R package to ensure reproducibility and easy dependency management.

### Prerequisites
- **R** (>= 4.2)
- **Quarto CLI** (for rendering analysis reports)

### Installation

You can install the package and all its dependencies directly from GitHub using [pak](https://pak.r-lib.org/):

```r
# Install pak if not already available
if (!require("pak")) install.packages("pak")

# Install directly from GitHub
pak::pak("apan1027/SwDGVirome")
```

### Local Development / Analysis
If you want to run the analysis notebooks locally:

1. **Clone the repository:**
   ```bash
   git clone https://github.com/apan1027/SwDGVirome.git
   cd SwDGVirome
   ```

2. **Initialize dependencies:**
   Open the project in your IDE (RStudio/Positron) and run:
   ```r
   pak::pak()
   ```

## Workflow Overview

The analysis pipeline involves sequential execution of Quarto documents. For a detailed list of scripts and execution order, please refer to the [Analyses README](analyses/README.md).

**Key Steps:**
1.  **Data Construction**: Build the core data object from raw database files.
2.  **Visualisation**: Generate main and supplementary figures for the manuscript.

## Usage

### 1. Prepare Input Data
Ensure the following files are placed in the resource directory `analyses/data/00-raw/d00-resource/`:
*   `p0057v2.sqlite`: Main project database.
*   `imgvr_source.tsv`: IMG/VR ecosystem source metadata.

### 2. Run Analysis
You can run the full workflow or individual steps using Quarto.

**Example:**
```bash
# Step 1: Construct the TSE object (Required)
quarto render analyses/01-tse-construction.qmd

# Step 2: Run downstream analyses
quarto render analyses/02-fig2-viral-diversity.qmd
```
See [analyses/README.md](analyses/README.md) for the complete execution order.

## Sample Groups

| Code | Description | Depth Range |
|------|-------------|-------------|
| **BS** | Baltic Sea | Surface (0 m) |
| **SA** | Shallow Aquifer | 71 m |
| **IA** | Intermediate Aquifer | 196 m |
| **DA** | Deep Aquifer | 450 m |

## Output
Outputs are organized by script name in the `analyses/data/` directory. Each folder contains high-resolution figures, source data tables, and R session information.

## Contact
**Project Maintainer**: Cunli Pan (cunli.pan@tum.de)
**Repository**: [apan1027/SwDGVirome](https://github.com/apan1027/SwDGVirome)


