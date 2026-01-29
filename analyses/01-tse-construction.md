# 01-tse-construction
Cunli Pan, Jinlong Ru
2025-12-20

- [<span class="toc-section-number">1</span> Tasks](#tasks)
  - [<span class="toc-section-number">1.1</span> Task 1: Read SQLite
    Database](#task-1-read-sqlite-database)
  - [<span class="toc-section-number">1.2</span> Task 2: Build Assay
    Matrices](#task-2-build-assay-matrices)
  - [<span class="toc-section-number">1.3</span> Task 3: Build colData
    (Sample Metadata)](#task-3-build-coldata-sample-metadata)
  - [<span class="toc-section-number">1.4</span> Task 4: Build feature
    metadata (rowData)](#task-4-build-feature-metadata-rowdata)
    - [<span class="toc-section-number">1.4.1</span> 4.2 Ecosystem
      classification](#42-ecosystem-classification)
    - [<span class="toc-section-number">1.4.2</span> 4.3 Host
      predictions](#43-host-predictions)
    - [<span class="toc-section-number">1.4.3</span> 4.4 Contig
      linkage](#44-contig-linkage)
  - [<span class="toc-section-number">1.5</span> Task 5: Task 5:
    Construct and validate TSE
    object](#task-5-task-5-construct-and-validate-tse-object)
    - [<span class="toc-section-number">1.5.1</span> 5.2 Construct
      TSE](#52-construct-tse)
    - [<span class="toc-section-number">1.5.2</span> 5.3 Save and
      validate](#53-save-and-validate)

**Updated: 2026-01-29 15:40:39 CET.**

The purpose of this document is to construct the foundational
`TreeSummarizedExperiment` (TSE) object by integrating viral abundance
profiles, taxonomy, host predictions, and ecosystem classification
metadata from the project’s SQLite database.

<details class="code-fold">
<summary>Code</summary>

``` r
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(data.table)
  library(DBI)
  library(RSQLite)
  library(S4Vectors)
  library(IRanges)
  library(SummarizedExperiment)
  library(TreeSummarizedExperiment)
})

# Load package utility functions
devtools::load_all(here::here())
```

</details>

## Tasks

### Task 1: Read SQLite Database

<details class="code-fold">
<summary>Code</summary>

``` r
# 1.1 Connect to database
sqlite_path <- file.path(path_resource, "p0057v2.sqlite")
stopifnot(file.exists(sqlite_path))

con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
message("✅ Connected to SQLite: ", basename(sqlite_path))
```

</details>

    ✅ Connected to SQLite: p0057v2.sqlite

<details class="code-fold">
<summary>Code</summary>

``` r
# 1.2 List all tables
tbls <- DBI::dbListTables(con)

# 1.3 Read core tables
sample_meta_raw <- DBI::dbReadTable(con, "sample_metadata")
contig_anno_raw <- DBI::dbReadTable(con, "contig_annotation")

# 1.4 Read abundance tables
abundance_map <- c(
  counts           = "abundance_count",
  tpm              = "abundance_tpm",
  rpkm             = "abundance_rpkm",
  reads_per_base   = "abundance_reads_per_base",
  trimmed_mean     = "abundance_trimmed_mean",
  covered_fraction = "abundance_covered_fraction"
)

assay_tbls <- abundance_map[abundance_map %in% tbls]
if (!length(assay_tbls)) {
  DBI::dbDisconnect(con)
  stop("❌ No abundance tables found in SQLite")
}

assay_dfs <- lapply(assay_tbls, function(tbl_name) {
  DBI::dbReadTable(con, tbl_name)
})
names(assay_dfs) <- names(assay_tbls)

# 1.5 Read AMG tables
amg_dramv_raw <- DBI::dbReadTable(con, "amg_dramv")
amg_vibrant_raw <- DBI::dbReadTable(con, "amg_vibrant")

# 1.6 Read Host prediction tables
host_genus_raw <- DBI::dbReadTable(con, "host_genus")
host_genome_raw <- DBI::dbReadTable(con, "host_genome")

# 1.7 Read Protein annotation tables
anno_prot_dramv <- DBI::dbReadTable(con, "anno_prot_dramv")
anno_prot_eggnog <- DBI::dbReadTable(con, "anno_prot_eggnog")
anno_prot_phrog <- DBI::dbReadTable(con, "anno_prot_phrog")

# 1.8 Read ARG tables
anno_arg_abricate <- DBI::dbReadTable(con, "anno_arg_abricate")
anno_arg_deeparg <- DBI::dbReadTable(con, "anno_arg_deeparg")

# 1.9 Read IMG/VR ecosystem source (external file)
imgvr_path <- file.path(path_resource, "imgvr_source.tsv")
if (file.exists(imgvr_path)) {
  imgvr_raw <- read_tsv(imgvr_path, show_col_types = FALSE)
} else {
  warning("⚠️  imgvr_source.tsv not found")
  imgvr_raw <- NULL
}

# 1.10 Close connection
DBI::dbDisconnect(con)

# Summary
message("✅ Data loaded: ",
        nrow(sample_meta_raw), " samples, ",
        nrow(contig_anno_raw), " contigs, ",
        length(assay_tbls), " abundance tables, ",
        nrow(amg_dramv_raw) + nrow(amg_vibrant_raw), " AMG entries")
```

</details>

    ✅ Data loaded: 4 samples, 2493 contigs, 6 abundance tables, 275 AMG entries

### Task 2: Build Assay Matrices

Convert abundance data frames to matrices (vOTU × sample).

<details class="code-fold">
<summary>Code</summary>

``` r
# 2.1 Convert abundance data frames to matrices
# Using safe_to_matrix() from R/utils.R
assay_mats <- lapply(assay_dfs, safe_to_matrix)

# 2.2 Determine common samples and vOTUs
all_samples <- lapply(assay_mats, colnames)
common_samples <- Reduce(base::intersect, all_samples)

all_features <- lapply(assay_mats, rownames)
common_features <- Reduce(base::intersect, c(list(unique(contig_anno_raw$vOTU_id)), all_features))

stopifnot(length(common_samples) > 0, length(common_features) > 0)

# 2.3 Subset and align assays
assay_mats <- lapply(assay_mats, function(m) {
  m[common_features, common_samples, drop = FALSE]
})

# Verify all assays have same dimensions
dims <- sapply(assay_mats, dim)
stopifnot(length(unique(dims[1,])) == 1, length(unique(dims[2,])) == 1)

message("✅ Assays aligned: ", nrow(assay_mats[[1]]), " vOTUs x ",
        ncol(assay_mats[[1]]), " samples")
```

</details>

    ✅ Assays aligned: 2488 vOTUs x 4 samples

### Task 3: Build colData (Sample Metadata)

<details class="code-fold">
<summary>Code</summary>

``` r
# 3.1 Build colData with sample_group mapping
colData_df <- sample_meta_raw %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples)) %>%
  mutate(
    sample_group = case_when(
      sample_source == "Baltic Sea" ~ "BS",
      sample_source == "507B" ~ "SA",
      sample_source == "1327B" ~ "IA",
      sample_source == "TASF" ~ "DA",
      TRUE ~ as.character(sample_source)
    ),
    sample_group = factor(sample_group, levels = c("BS", "SA", "IA", "DA"))
  ) %>%
  column_to_rownames("sample_id")

# 3.2 Verify alignment
stopifnot(identical(rownames(colData_df), colnames(assay_mats[[1]])))

message("✅ Sample metadata built: ", nrow(colData_df), " samples (",
        paste(levels(colData_df$sample_group), collapse = ", "), ")")
```

</details>

    ✅ Sample metadata built: 4 samples (BS, SA, IA, DA)

### Task 4: Build feature metadata (rowData)

Complete vOTU-level metadata including taxonomy, quality metrics,
ecosystem classification, host predictions, and contig linkage. \####
4.1 Base taxonomy and quality metrics

<details class="code-fold">
<summary>Code</summary>

``` r
# 4.1.1 Select representative contig for each vOTU
# Strategy: Prioritize contigs with taxonomy info > taxonomy_priority > completeness > length
rep_contigs <- contig_anno_raw %>%
  mutate(
    completeness_num = suppressWarnings(as.numeric(checkv_completeness)),
    length_num = suppressWarnings(as.numeric(contig_length)),
    has_taxonomy = (!is.na(family) & family != "") | (!is.na(genus) & genus != "")
  ) %>%
  filter(vOTU_id %in% common_features) %>%
  group_by(vOTU_id) %>%
  arrange(
    desc(has_taxonomy),
    taxonomy_priority,
    desc(completeness_num),
    desc(length_num)
  ) %>%
  dplyr::slice_head(n = 1) %>%
  ungroup()

# 4.1.2 Calculate vOTU-level stats
votu_stats <- contig_anno_raw %>%
  filter(vOTU_id %in% common_features) %>%
  group_by(vOTU_id) %>%
  summarise(
    n_contigs = n_distinct(contig_id, na.rm = TRUE),
    total_length = sum(suppressWarnings(as.numeric(contig_length)), na.rm = TRUE),
    .groups = "drop"
  )

# 4.1.3 Extract taxonomy and quality from representative contigs
rowData_base <- rep_contigs %>%
  select(
    vOTU_id,
    # Taxonomy (8 columns)
    realm, kingdom, phylum, class, order, family, genus, species,
    # Quality metrics
    checkv_completeness, checkv_contamination, checkv_completeness_method,
    genomad_virus_score, genomad_n_hallmarks,
    # Lifestyle
    lifestyle, bacphlip_lifestyle,
    # vContact3
    vc_id, vc_genus, vc_novel_genus
  ) %>%
  left_join(votu_stats, by = "vOTU_id") %>%
  arrange(match(vOTU_id, rownames(assay_mats[[1]]))) %>%
  column_to_rownames("vOTU_id")

# 4.1.4 Add quality tier classification
rowData_base <- rowData_base %>%
  mutate(
    completeness_num = suppressWarnings(as.numeric(checkv_completeness)),
    quality_tier = case_when(
      completeness_num >= 90 ~ "Complete",
      completeness_num >= 50 ~ "High-quality",
      completeness_num >= 30 ~ "Medium-quality",
      TRUE ~ "Low-quality"
    )
  ) %>%
  select(-completeness_num)

stopifnot(identical(rownames(rowData_base), rownames(assay_mats[[1]])))

message("✅ 4.1 Base metadata: ", nrow(rowData_base), " vOTUs with taxonomy and quality")
```

</details>

    ✅ 4.1 Base metadata: 2488 vOTUs with taxonomy and quality

#### 4.2 Ecosystem classification

<details class="code-fold">
<summary>Code</summary>

``` r
 #| label: task4.2-ecosystem
#| message: true
#| warning: true

if (!is.null(imgvr_raw)) {

  # 4.2.1 Using extract_level5_position5() from R/utils.R

  # 4.2.2 Extract Level 5 from imgvr
  imgvr_eco <- imgvr_raw %>%
    select(contig_id, eco_class = `Ecosystem classification`) %>%
    mutate(level5 = sapply(eco_class, extract_level5_position5)) %>%
    filter(!is.na(level5)) %>%
    select(contig_id, eco_level5 = level5) %>%
    distinct()

  # 4.2.3 Get contig_id for each vOTU
  votu_to_contig <- rep_contigs %>%
    filter(vOTU_id %in% rownames(rowData_base)) %>%
    select(vOTU_id, contig_id)

  # 4.2.4 Merge ecosystem info
  votu_eco <- votu_to_contig %>%
    left_join(imgvr_eco, by = "contig_id") %>%
    select(vOTU_id, eco_level5)

  # 4.2.5 Add to rowData_base
  rowData_base <- rowData_base %>%
    rownames_to_column("vOTU_id") %>%
    left_join(votu_eco, by = "vOTU_id") %>%
    arrange(match(vOTU_id, rownames(assay_mats[[1]]))) %>%
    column_to_rownames("vOTU_id")

  n_with_eco <- sum(!is.na(rowData_base$eco_level5))
  message("✅ 4.2 Ecosystem: ", n_with_eco, "/", nrow(rowData_base),
          " vOTUs (", round(n_with_eco/nrow(rowData_base)*100, 1), "%)")

} else {
  warning("⚠️  Skipping ecosystem annotation (imgvr_source.tsv not loaded)")
}
```

</details>

    ✅ 4.2 Ecosystem: 1638/2488 vOTUs (65.8%)

<details class="code-fold">
<summary>Code</summary>

``` r
stopifnot(identical(rownames(rowData_base), rownames(assay_mats[[1]])))
```

</details>

#### 4.3 Host predictions

<details class="code-fold">
<summary>Code</summary>

``` r
# 4.3.1 Aggregate host genus predictions
host_genus_agg <- host_genus_raw %>%
  filter(vOTU_id %in% rownames(rowData_base)) %>%
  group_by(vOTU_id) %>%
  summarise(
    hostGenusset = paste(unique(Host.genus), collapse = ";"),
    hostnumbergenus = n_distinct(Host.genus),
    .groups = "drop"
  ) %>%
  mutate(
    hostSG = case_when(
      hostnumbergenus == 1 ~ "specialist",
      hostnumbergenus > 1 ~ "generalist",
      TRUE ~ NA_character_
    )
  )

# 4.3.2 Aggregate host species predictions
host_species_agg <- host_genome_raw %>%
  filter(vOTU_id %in% rownames(rowData_base)) %>%
  group_by(vOTU_id) %>%
  summarise(
    hostSpeciesset = paste(unique(Host.taxonomy), collapse = ";"),
    hostnumberspecies = n_distinct(Host.taxonomy),
    .groups = "drop"
  )

# 4.3.3 Merge host info into rowData
rowData_base <- rowData_base %>%
  rownames_to_column("vOTU_id") %>%
  left_join(host_genus_agg, by = "vOTU_id") %>%
  left_join(host_species_agg, by = "vOTU_id") %>%
  arrange(match(vOTU_id, rownames(assay_mats[[1]]))) %>%
  column_to_rownames("vOTU_id")

stopifnot(identical(rownames(rowData_base), rownames(assay_mats[[1]])))

n_with_host <- sum(!is.na(rowData_base$hostGenusset))
message("✅ 4.3 Host predictions: ", n_with_host, "/", nrow(rowData_base),
        " vOTUs (", round(n_with_host/nrow(rowData_base)*100, 1), "%)")
```

</details>

    ✅ 4.3 Host predictions: 909/2488 vOTUs (36.5%)

#### 4.4 Contig linkage

<details class="code-fold">
<summary>Code</summary>

``` r
# 4.4.1 Build vOTU → contig mapping
bridge_votu_contig <- contig_anno_raw %>%
  filter(!is.na(vOTU_id), !is.na(contig_id),
         vOTU_id != "", contig_id != "",
         vOTU_id %in% rownames(rowData_base)) %>%
  distinct(vOTU_id, contig_id) %>%
  arrange(vOTU_id, contig_id)

# 4.4.2 Create CharacterList
v2c_map <- split(bridge_votu_contig$contig_id, bridge_votu_contig$vOTU_id)
v2c_map <- lapply(v2c_map, unique)

lookup <- function(id) {
  x <- v2c_map[[id]]
  if (is.null(x)) character(0) else as.character(x)
}

votu_ids_ordered <- rownames(rowData_base)
contigs_aligned <- lapply(votu_ids_ordered, lookup)
contig_ids_list <- IRanges::CharacterList(contigs_aligned)
names(contig_ids_list) <- votu_ids_ordered

# 4.4.3 Add to rowData
rowData_base$contig_ids <- contig_ids_list

stopifnot(identical(rownames(rowData_base), rownames(assay_mats[[1]])))

message("✅ 4.4 Contig linkage: ", nrow(bridge_votu_contig), " edges, ",
        sum(lengths(contig_ids_list) > 1), " multi-contig vOTUs")
```

</details>

    ✅ 4.4 Contig linkage: 2493 edges, 5 multi-contig vOTUs

<details class="code-fold">
<summary>Code</summary>

``` r
# Final rowData summary
message("✅ Task 4 completed: rowData with ", ncol(rowData_base), " columns")
```

</details>

    ✅ Task 4 completed: rowData with 28 columns

### Task 5: Task 5: Construct and validate TSE object

Prepare metadata, construct TreeSummarizedExperiment, and validate
output. \#### 5.1 Prepare metadata list

<details class="code-fold">
<summary>Code</summary>

``` r
# Prepare complete metadata list
metadata_list <- list(
  # Construction info
  construction_info = list(
    created_date  = Sys.Date(),
    created_time  = Sys.time(),
    sqlite_path   = normalizePath(sqlite_path),
    n_samples_raw = nrow(sample_meta_raw),
    n_samples     = length(common_samples),
    n_votus_raw   = length(unique(contig_anno_raw$vOTU_id)),
    n_votus       = nrow(rowData_base),
    assays_loaded = names(assay_mats),
    strategy_representative_contig = "has_taxonomy > taxonomy_priority > completeness > length"
  ),

  # Protein annotations
  anno_prot_dramv = anno_prot_dramv,
  anno_prot_eggnog = anno_prot_eggnog,
  anno_prot_phrog = anno_prot_phrog,

  # AMG tables
  amg_dramv = amg_dramv_raw,
  amg_vibrant = amg_vibrant_raw,

  # ARG tables
  anno_arg_abricate = anno_arg_abricate,
  anno_arg_deeparg = anno_arg_deeparg,

  # Host prediction edge tables
  host_genus_edges = host_genus_raw,
  host_genome_edges = host_genome_raw,

  # Contig linkage info
  linkage_info = list(
    n_votu_contig_edges = nrow(bridge_votu_contig),
    multi_contig_votus = names(v2c_map)[lengths(v2c_map) > 1],
    n_multi_contig_votus = sum(lengths(v2c_map) > 1)
  )
)

# Add complete contig annotation
metadata_list$contig_annotation <- contig_anno_raw %>%
  filter(vOTU_id %in% rownames(rowData_base))

# Add complete imgvr source if available
if (!is.null(imgvr_raw)) {
  metadata_list$imgvr_source <- imgvr_raw
}

message("✅ 5.1 Metadata prepared: ", length(metadata_list), " components")
```

</details>

    ✅ 5.1 Metadata prepared: 13 components

#### 5.2 Construct TSE

<details class="code-fold">
<summary>Code</summary>

``` r
# 5.2.1 Final integrity checks
stopifnot(
  identical(rownames(assay_mats[[1]]), rownames(rowData_base)),
  identical(colnames(assay_mats[[1]]), rownames(colData_df))
)

message("✅ 5.2 All dimensions verified")
```

</details>

    ✅ 5.2 All dimensions verified

<details class="code-fold">
<summary>Code</summary>

``` r
# 5.2.2 Construct TSE
tse <- TreeSummarizedExperiment(
  assays  = S4Vectors::SimpleList(assay_mats),
  colData = S4Vectors::DataFrame(colData_df),
  rowData = S4Vectors::DataFrame(rowData_base),
  metadata = metadata_list
)

message("✅ TSE constructed successfully!")
```

</details>

    ✅ TSE constructed successfully!

<details class="code-fold">
<summary>Code</summary>

``` r
message("   Dimensions: ", nrow(tse), " vOTUs × ", ncol(tse), " samples")
```

</details>

       Dimensions: 2488 vOTUs × 4 samples

<details class="code-fold">
<summary>Code</summary>

``` r
message("   Assays: ", paste(names(assays(tse)), collapse = ", "))
```

</details>

       Assays: counts, tpm, rpkm, reads_per_base, trimmed_mean, covered_fraction

<details class="code-fold">
<summary>Code</summary>

``` r
message("   Metadata components: ", length(metadata(tse)))
```

</details>

       Metadata components: 13

#### 5.3 Save and validate

<details class="code-fold">
<summary>Code</summary>

``` r
# 5.3.1 Save TSE
output_path <- path_target("tse.rds")
saveRDS(tse, output_path)

# 5.3.2 Save session info
session_info_path <- path_target("session_info.txt")
writeLines(capture.output(sessionInfo()), session_info_path)

message("✅ Files saved:")
```

</details>

    ✅ Files saved:

<details class="code-fold">
<summary>Code</summary>

``` r
message("   - ", basename(output_path))
```

</details>

       - tse.rds

<details class="code-fold">
<summary>Code</summary>

``` r
message("   - ", basename(session_info_path))
```

</details>

       - session_info.txt

<details class="code-fold">
<summary>Code</summary>

``` r
# 5.3.3 Reload and validate
tse_reload <- readRDS(output_path)

# Quick validation checks
validation_passed <- all(
  identical(dim(tse_reload), dim(tse)),
  identical(assayNames(tse_reload), assayNames(tse)),
  "amg_dramv" %in% names(metadata(tse_reload)),
  "host_genome_edges" %in% names(metadata(tse_reload)),
  "contig_ids" %in% colnames(rowData(tse_reload)),
  "sample_group" %in% colnames(colData(tse_reload))
)

if (validation_passed) {
  message("✅ Validation passed: TSE ready for downstream analysis")
} else {
  stop("❌ Validation failed")
}
```

</details>

    ✅ Validation passed: TSE ready for downstream analysis

<details class="code-fold">
<summary>Code</summary>

``` r
# Display final summary
message("\n=== Final TSE Summary ===")
```

</details>


    === Final TSE Summary ===

<details class="code-fold">
<summary>Code</summary>

``` r
message("Dimensions: ", nrow(tse), " vOTUs × ", ncol(tse), " samples")
```

</details>

    Dimensions: 2488 vOTUs × 4 samples

<details class="code-fold">
<summary>Code</summary>

``` r
message("Sample groups: ", paste(levels(colData(tse)$sample_group), collapse = ", "))
```

</details>

    Sample groups: BS, SA, IA, DA

<details class="code-fold">
<summary>Code</summary>

``` r
message("Assays: ", paste(assayNames(tse), collapse = ", "))
```

</details>

    Assays: counts, tpm, rpkm, reads_per_base, trimmed_mean, covered_fraction

<details class="code-fold">
<summary>Code</summary>

``` r
message("rowData fields: ", ncol(rowData(tse)))
```

</details>

    rowData fields: 28

<details class="code-fold">
<summary>Code</summary>

``` r
message("Metadata tables: ", length(metadata(tse)))
```

</details>

    Metadata tables: 13
