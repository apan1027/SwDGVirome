06-figs2-lifestyle
================
today

**Updated: 2026-09-24 14:14:33 CET.**

``` r
suppressPackageStartupMessages({
  library(here)
  library(conflicted)
  library(tidyverse)
  library(data.table)
  library(patchwork)
  devtools::load_all()
})
```

## Tasks

This section is designated for your workflow implementation. **Please
write all your analysis and processing code here.**

Maintaining your code within this section ensures a clean project
structure and facilitates reproducibility. We recommend organizing your
workflow into modular tasks—such as data cleaning, analysis, and
visualization—using distinct code chunks.

<div class="callout-note">

## Example structure

This is an example structure to get you started. You can modify it as
per your needs.

</div>

## Tasks

### Task 1: Load the complete discovery catalogue

``` r
# Figure S2 is a sensitivity analysis of predicted lifestyle
# against contig length and CheckV quality.
#
# Therefore, unlike the main analyses, it uses the complete
# >=3 kb discovery catalogue so that the 3–5 kb group is retained.

tse_discovery <- readRDS(
  here(
    "data",
    "01-tse-construction",
    "tse_discovery_2488.rds"
  )
)

votu_metadata <- SummarizedExperiment::rowData(
  tse_discovery
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("vOTU_id")

stopifnot(nrow(votu_metadata) == 2488)

stopifnot(
  all(
    c(
      "representative_length_bp",
      "checkv_quality",
      "lifestyle"
    ) %in% colnames(votu_metadata)
  )
)

message(
  "Complete discovery catalogue loaded: ",
  nrow(votu_metadata),
  " vOTUs"
)
```

\###Task 2: Define length and CheckV quality groups

``` r
votu_metadata <- votu_metadata %>%
  dplyr::mutate(
    is_temperate = lifestyle == "temperate",

    length_group = cut(
      representative_length_bp,
      breaks = c(
        3000,
        5000,
        10000,
        20000,
        Inf
      ),
      right = FALSE,
      labels = c(
        "3–5 kb",
        "5–10 kb",
        "10–20 kb",
        "≥20 kb"
      )
    ),

    quality_group = factor(
      checkv_quality,
      levels = c(
        "Complete",
        "High-quality",
        "Medium-quality",
        "Low-quality",
        "Not-determined"
      )
    )
  )

stopifnot(!any(is.na(votu_metadata$length_group)))
stopifnot(!any(is.na(votu_metadata$quality_group)))

message(
  "Length and CheckV quality groups created"
)
```

### Task 3: Calculate predicted temperate fractions

``` r
# Fraction predicted temperate within each contig-length group
length_summary <- votu_metadata %>%
  dplyr::count(
    length_group,
    is_temperate,
    name = "n"
  ) %>%
  tidyr::complete(
    length_group,
    is_temperate,
    fill = list(n = 0)
  ) %>%
  dplyr::group_by(length_group) %>%
  dplyr::summarise(
    n_vOTU = sum(n),
    n_temperate = sum(n[is_temperate]),
    temperate_fraction = n_temperate / n_vOTU,
    .groups = "drop"
  )

# Fraction predicted temperate within each CheckV-quality group
quality_summary <- votu_metadata %>%
  dplyr::count(
    quality_group,
    is_temperate,
    name = "n"
  ) %>%
  tidyr::complete(
    quality_group,
    is_temperate,
    fill = list(n = 0)
  ) %>%
  dplyr::group_by(quality_group) %>%
  dplyr::summarise(
    n_vOTU = sum(n),
    n_temperate = sum(n[is_temperate]),
    temperate_fraction = n_temperate / n_vOTU,
    .groups = "drop"
  )

print(length_summary, n = Inf)
```

    ## # A tibble: 4 × 4
    ##   length_group n_vOTU n_temperate temperate_fraction
    ##   <fct>         <int>       <int>              <dbl>
    ## 1 3–5 kb         1526          12            0.00786
    ## 2 5–10 kb         689           9            0.0131
    ## 3 10–20 kb        207           7            0.0338
    ## 4 ≥20 kb           66           7            0.106

``` r
print(quality_summary, n = Inf)
```

    ## # A tibble: 5 × 4
    ##   quality_group  n_vOTU n_temperate temperate_fraction
    ##   <fct>           <int>       <int>              <dbl>
    ## 1 Complete            4           0             0
    ## 2 High-quality       12           2             0.167
    ## 3 Medium-quality     39           4             0.103
    ## 4 Low-quality      2311          25             0.0108
    ## 5 Not-determined    122           4             0.0328

\###Task 4: Validate the Figure S2 numbers

``` r
# Expected number of vOTUs in each length category
stopifnot(
  identical(
    length_summary$n_vOTU,
    c(
      1526L,
      689L,
      207L,
      66L
    )
  )
)

# Expected number predicted as temperate
stopifnot(
  identical(
    length_summary$n_temperate,
    c(
      12L,
      9L,
      7L,
      7L
    )
  )
)

# Expected number of vOTUs in each CheckV category
stopifnot(
  identical(
    quality_summary$n_vOTU,
    c(
      4L,
      12L,
      39L,
      2311L,
      122L
    )
  )
)

# Expected number predicted as temperate
stopifnot(
  identical(
    quality_summary$n_temperate,
    c(
      0L,
      2L,
      4L,
      25L,
      4L
    )
  )
)

stopifnot(
  sum(length_summary$n_vOTU) == 2488
)

stopifnot(
  sum(length_summary$n_temperate) == 35
)

stopifnot(
  sum(quality_summary$n_vOTU) == 2488
)

stopifnot(
  sum(quality_summary$n_temperate) == 35
)

message(
  "All Figure S2 numerical checks passed"
)
```

\###Task 5: Create Figure S2

``` r
panel_theme <-
  theme_classic(
    base_size = 18,
    base_family = "Times"
  ) +
  theme(
    axis.line = element_line(
      linewidth = 0.8,
      colour = "black"
    ),
    axis.ticks = element_line(
      linewidth = 0.7,
      colour = "black"
    ),
    axis.text = element_text(
      colour = "black"
    ),
    plot.margin = margin(
      10,
      10,
      10,
      10
    )
  )

# Figure S2a:
# predicted temperate fraction by representative contig length
p_figs2a <- ggplot(
  length_summary,
  aes(
    x = length_group,
    y = temperate_fraction
  )
) +
  geom_col(
    width = 0.68,
    fill = "#3274A1",
    colour = "black",
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label = paste0(
        n_temperate,
        "/",
        n_vOTU
      )
    ),
    vjust = -0.45,
    size = 5,
    family = "Times"
  ) +
  scale_y_continuous(
    limits = c(0, 0.20),
    breaks = seq(0, 0.20, 0.05),
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
    x = NULL,
    y = "Fraction of vOTUs predicted temperate"
  ) +
  panel_theme

# Figure S2b:
# predicted temperate fraction by CheckV quality category
p_figs2b <- ggplot(
  quality_summary,
  aes(
    x = quality_group,
    y = temperate_fraction
  )
) +
  geom_col(
    width = 0.68,
    fill = "#3274A1",
    colour = "black",
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label = paste0(
        n_temperate,
        "/",
        n_vOTU
      )
    ),
    vjust = -0.45,
    size = 5,
    family = "Times"
  ) +
  scale_y_continuous(
    limits = c(0, 0.185),
    breaks = seq(0, 0.15, 0.05),
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
    x = NULL,
    y = "Fraction of vOTUs predicted temperate"
  ) +
  panel_theme +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

# Combined preview
p_figs2_combined <-
  p_figs2a +
  p_figs2b +
  patchwork::plot_annotation(
    tag_levels = "a"
  )

print(p_figs2_combined)
```

![](06-figs2-lifestyle_files/figure-commonmark/task5-create-figure-1.png)<!-- -->

### Task 6: Save source data, PNG and TIFF files

### Task 7: Quality-stratified replication strategies in the primary catalogue

``` r
# Keep this analysis separate from the 2,488-vOTU Figure S2 objects.
local({
  primary <- readRDS(
    here::here("data", "01-tse-construction", "tse_primary_962.rds")
  )
  metadata <- as.data.frame(SummarizedExperiment::rowData(primary))
  abundance <- SummarizedExperiment::assay(primary, "tpm")
  sample_names <- c("BS", "SA", "IA", "DA")
  sample_index <- match(
    sample_names,
    as.character(SummarizedExperiment::colData(primary)$sample_group)
  )
  stopifnot(
    nrow(primary) == 962L,
    ncol(primary) == 4L,
    !anyNA(sample_index),
    !anyDuplicated(rownames(primary))
  )
  abundance <- abundance[, sample_index, drop = FALSE]
  colnames(abundance) <- sample_names
  quality <- as.character(metadata$checkv_quality)
  strategy <- as.character(metadata$lifestyle)
  quality_levels <- c(
    "Complete", "High-quality", "Medium-quality",
    "Low-quality", "Not-determined"
  )
  stopifnot(
    !anyNA(quality), all(quality %in% quality_levels),
    !anyNA(strategy), all(strategy %in% c("temperate", "virulent")),
    !anyNA(abundance), all(is.finite(abundance)), all(abundance >= 0)
  )
  mq_or_better <- quality %in% quality_levels[1:3]
  stopifnot(sum(mq_or_better) == 55L)

  select_group <- function(group) {
    if (group == "All primary") {
      rep(TRUE, nrow(primary))
    } else if (group == "Medium-quality or better") {
      mq_or_better
    } else {
      quality == group
    }
  }

  # Table S3a: catalogue counts, not TPM-weighted fractions.
  count_groups <- c(quality_levels, "All primary", "Medium-quality or better")
  quality_counts <- do.call(rbind, lapply(count_groups, function(group) {
    keep <- select_group(group)
    n_total <- sum(keep)
    n_temperate <- sum(keep & strategy == "temperate")
    data.frame(
      quality_group = group,
      n_vOTU = n_total,
      n_temperate = n_temperate,
      n_virulent = sum(keep & strategy == "virulent"),
      n_unassigned = 0L,
      temperate_count_fraction = n_temperate / n_total
    )
  }))

  # Table S3b: T/(T+V) within each sample and quality subset.
  tpm_groups <- c("All primary", "Medium-quality or better", quality_levels)
  results <- list()
  for (group in tpm_groups) {
    keep <- select_group(group)
    for (sample in sample_names) {
      a <- abundance[, sample]
      t <- sum(a[keep & strategy == "temperate"])
      v <- sum(a[keep & strategy == "virulent"])
      denominator <- t + v
      results[[length(results) + 1L]] <- data.frame(
        quality_group = group,
        sample = sample,
        n_catalogue = sum(keep),
        n_detected = sum(keep & a > 0),
        n_temperate_detected = sum(keep & a > 0 & strategy == "temperate"),
        n_virulent_detected = sum(keep & a > 0 & strategy == "virulent"),
        temperate_TPM = t,
        virulent_TPM = v,
        classified_TPM = denominator,
        temperate_fraction = if (denominator > 0) t / denominator else NA_real_,
        fraction_of_primary_TPM = denominator / sum(a)
      )
    }
  }
  sample_fractions <- do.call(rbind, results)
  all_rows <- sample_fractions[sample_fractions$quality_group == "All primary", ]
  mq_rows <- sample_fractions[
    sample_fractions$quality_group == "Medium-quality or better", ]
  comparison <- data.frame(
    sample = sample_names,
    primary_temperate_fraction = all_rows$temperate_fraction,
    MQ_or_better_temperate_fraction = mq_rows$temperate_fraction,
    MQ_or_better_n_detected = mq_rows$n_detected,
    MQ_or_better_n_temperate_detected = mq_rows$n_temperate_detected,
    MQ_or_better_fraction_of_primary_TPM = mq_rows$fraction_of_primary_TPM
  )

  # Validate the quality totals and reconcile the mutually exclusive groups.
  stopifnot(
    identical(quality_counts$n_vOTU[1:5], c(4L, 12L, 39L, 874L, 33L)),
    identical(quality_counts$n_temperate[1:5], c(0L, 2L, 4L, 16L, 1L)),
    sum(quality_counts$n_vOTU[1:5]) == 962L,
    sum(quality_counts$n_temperate[1:5]) == 23L,
    all(quality_counts$n_vOTU == quality_counts$n_temperate + quality_counts$n_virulent)
  )
  for (sample in sample_names) {
    rows <- sample_fractions$sample == sample &
      sample_fractions$quality_group %in% quality_levels
    stopifnot(abs(sum(sample_fractions$classified_TPM[rows]) -
      sum(abundance[, sample])) < 1e-7)
  }

  table_dir <- path_target("quality_sensitivity")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(quality_counts,
    file.path(table_dir, "Table_S3a_quality_counts.csv"), row.names = FALSE)
  utils::write.csv(sample_fractions,
    file.path(table_dir, "Table_S3b_quality_sample_fractions.csv"),
    row.names = FALSE, na = "")
  utils::write.csv(comparison,
    file.path(table_dir, "Table_S3_comparison.csv"), row.names = FALSE)
  source_data <- data.frame(
    vOTU_id = rownames(primary), checkv_quality = quality,
    final_strategy = strategy, MQ_or_better = as.integer(mq_or_better),
    abundance, check.names = FALSE
  )
  utils::write.csv(source_data,
    file.path(table_dir, "Table_S3_source_vOTU_TPM.csv"), row.names = FALSE)

  print(knitr::kable(quality_counts, digits = 4,
    caption = "Table S3a. Quality-stratified predicted replication-strategy counts."))
  print(knitr::kable(comparison, digits = 4,
    caption = "Table S3b. Sensitivity of TPM-weighted predicted temperate fractions."))
  message("Table S3 CSV files written to: ", normalizePath(table_dir))
  message("A missing fraction means T + V = 0 (not estimable), not a zero fraction.")
})
```

    ##
    ##
    ## Table: Table S3a. Quality-stratified predicted replication-strategy counts.
    ##
    ## |quality_group            | n_vOTU| n_temperate| n_virulent| n_unassigned| temperate_count_fraction|
    ## |:------------------------|------:|-----------:|----------:|------------:|------------------------:|
    ## |Complete                 |      4|           0|          4|            0|                   0.0000|
    ## |High-quality             |     12|           2|         10|            0|                   0.1667|
    ## |Medium-quality           |     39|           4|         35|            0|                   0.1026|
    ## |Low-quality              |    874|          16|        858|            0|                   0.0183|
    ## |Not-determined           |     33|           1|         32|            0|                   0.0303|
    ## |All primary              |    962|          23|        939|            0|                   0.0239|
    ## |Medium-quality or better |     55|           6|         49|            0|                   0.1091|
    ##
    ##
    ## Table: Table S3b. Sensitivity of TPM-weighted predicted temperate fractions.
    ##
    ## |sample | primary_temperate_fraction| MQ_or_better_temperate_fraction| MQ_or_better_n_detected| MQ_or_better_n_temperate_detected| MQ_or_better_fraction_of_primary_TPM|
    ## |:------|--------------------------:|-------------------------------:|-----------------------:|---------------------------------:|------------------------------------:|
    ## |BS     |                     0.0193|                          0.0812|                      39|                                 3|                               0.0700|
    ## |SA     |                     0.1616|                          0.4914|                      45|                                 5|                               0.1321|
    ## |IA     |                     0.0386|                          0.1954|                      18|                                 3|                               0.1041|
    ## |DA     |                     0.1892|                          0.5324|                      17|                                 1|                               0.3058|

### Task 8: Export Table S3 with taxonomic coverage

``` r
# Add after Task 7 in 06-figs2-lifestyle.qmd; keep Task 7 unchanged.
local({
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Please install openxlsx before running Task 8.")
  }
  table_dir <- path_target("quality_sensitivity")
  input_file <- path_source("01-tse-construction", "tse_primary_962.rds")
  primary <- readRDS(input_file)
  rd <- as.data.frame(SummarizedExperiment::rowData(primary))
  sample_names <- c("BS", "SA", "IA", "DA")
  ranks <- c("class", "order", "family", "genus")
  ids <- rownames(primary)
  sample_index <- match(sample_names,
    as.character(SummarizedExperiment::colData(primary)$sample_group))
  stopifnot(nrow(primary) == 962L, ncol(primary) == 4L,
    !anyDuplicated(ids), !anyNA(sample_index),
    all(c(ranks, "checkv_quality", "lifestyle",
      "representative_contig_id", "representative_length_bp") %in% names(rd)))
  tpm <- SummarizedExperiment::assay(primary, "tpm")[, sample_index, drop = FALSE]
  colnames(tpm) <- sample_names
  stopifnot(identical(rownames(tpm), ids), identical(rownames(rd), ids),
    !anyNA(tpm), all(is.finite(tpm)), all(tpm >= 0), all(colSums(tpm) > 0))
  # Read the quality results already generated by Task 7.
  read_task7 <- function(filename) {
    p <- file.path(table_dir, filename)
    if (!file.exists(p)) stop("Run Task 7 first. Missing file: ", p)
    utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  }
  quality_counts <- read_task7("Table_S3a_quality_counts.csv")
  sample_fractions <- read_task7("Table_S3b_quality_sample_fractions.csv")
  comparison <- read_task7("Table_S3_comparison.csv")
  source <- read_task7("Table_S3_source_vOTU_TPM.csv")
  stopifnot(nrow(source) == 962L, !anyDuplicated(source$vOTU_id),
    setequal(source$vOTU_id, ids))
  source <- source[match(ids, source$vOTU_id), , drop = FALSE]
  rownames(source) <- NULL
  quality <- as.character(rd$checkv_quality)
  strategy <- as.character(rd$lifestyle)
  mq <- quality %in% c("Complete", "High-quality", "Medium-quality")
  same_numbers <- function(x, y) {
    isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 1e-10))
  }
  stopifnot(identical(source$vOTU_id, ids),
    identical(source$checkv_quality, quality),
    identical(source$final_strategy, strategy),
    same_numbers(source$MQ_or_better, as.integer(mq)),
    same_numbers(as.matrix(source[sample_names]), tpm))
  quality_mask <- function(group) {
    if (group == "All primary") return(rep(TRUE, length(ids)))
    if (group == "Medium-quality or better") return(mq)
    quality == group
  }
   # Catch stale Task 7 summaries before combining them with current taxonomy.
  for (i in seq_len(nrow(quality_counts))) {
    keep <- quality_mask(quality_counts$quality_group[i])
    stopifnot(quality_counts$n_vOTU[i] == sum(keep),
      quality_counts$n_temperate[i] == sum(keep & strategy == "temperate"),
      quality_counts$n_virulent[i] == sum(keep & strategy == "virulent"))
  }
  for (i in seq_len(nrow(sample_fractions))) {
    keep <- quality_mask(sample_fractions$quality_group[i])
    a <- tpm[, sample_fractions$sample[i]]
    tt <- sum(a[keep & strategy == "temperate"])
    vv <- sum(a[keep & strategy == "virulent"])
    stopifnot(same_numbers(sample_fractions$temperate_TPM[i], tt),
      same_numbers(sample_fractions$virulent_TPM[i], vv),
      same_numbers(sample_fractions$temperate_fraction[i],
        if (tt + vv > 0) tt / (tt + vv) else NA_real_))
  }
  for (i in seq_len(nrow(comparison))) {
    s <- comparison$sample[i]
    full <- sample_fractions$quality_group == "All primary" & sample_fractions$sample == s
    high <- sample_fractions$quality_group == "Medium-quality or better" & sample_fractions$sample == s
    stopifnot(sum(full) == 1L, sum(high) == 1L,
      same_numbers(comparison$primary_temperate_fraction[i], sample_fractions$temperate_fraction[full]),
      same_numbers(comparison$MQ_or_better_temperate_fraction[i], sample_fractions$temperate_fraction[high]))
  }
  # Sample-specific groups use the total TPM of all 962 representatives.
  relative_tpm <- sweep(tpm, 2, colSums(tpm), "/")
  group <- matrix("Not detected", nrow(tpm), ncol(tpm), dimnames = dimnames(tpm))
  group[relative_tpm > 0 & relative_tpm < 0.001] <- "Rare"
  group[relative_tpm >= 0.001 & relative_tpm < 0.01] <- "Intermediate"
  group[relative_tpm >= 0.01] <- "Abundant"
  stopifnot(all(colSums(group == "Abundant") == c(4L, 37L, 31L, 20L)),
    all(colSums(group == "Rare") == c(552L, 708L, 118L, 165L)),
    all(colSums(tpm > 0) == c(864L, 786L, 156L, 211L)))
  # Preserve raw labels. A [Rank]_ placeholder is not an assignment at that rank.
  raw_taxonomy <- lapply(rd[ranks], as.character)
  assigned <- lapply(raw_taxonomy, function(x) {
    x <- trimws(x)
    missing <- is.na(x) | x == "" |
      grepl("^(unclassified|unknown|unassigned)([ _-]|$)|^(NA|N/A|-)$", x, ignore.case = TRUE)
    placeholder <- grepl("^\\[(Class|Order|Family|Genus)\\]_", x, ignore.case = TRUE)
    !missing & !placeholder
  })
  clean_taxonomy <- lapply(ranks, function(rank) {
    ifelse(assigned[[rank]], trimws(raw_taxonomy[[rank]]), "Unclassified")
  })
  names(clean_taxonomy) <- ranks
  detail <- data.frame(vOTU_id = ids, clean_taxonomy, stringsAsFactors = FALSE)
  for (s in sample_names) detail[[paste0(s, "_group")]] <- group[, s]
  detail$checkv_quality <- quality
  detail$final_strategy <- strategy
  detail$MQ_or_better <- as.integer(mq)
  detail$representative_contig_id <- as.character(rd$representative_contig_id)
  detail$representative_length_bp <- rd$representative_length_bp
  for (s in sample_names) {
    detail[[paste0(s, "_TPM")]] <- tpm[, s]
    detail[[paste0(s, "_relative_TPM")]] <- relative_tpm[, s]
  }
  for (rank in ranks) {
    detail[[paste0(rank, "_assigned")]] <- assigned[[rank]]
    detail[[paste0(rank, "_raw")]] <- raw_taxonomy[[rank]]
  }
  # All-primary overview + four samples x two subsets x four ranks = 36 rows.
  rank_summary <- function(keep, sample, subset) {
    n <- sum(keep)
    do.call(rbind, lapply(ranks, function(rank) {
      a <- sum(assigned[[rank]][keep])
      data.frame(sample = sample, subset = subset, rank = rank, n_vOTU = n,
        n_assigned = a, assigned_fraction = if (n > 0) a / n else NA_real_,
        n_unclassified = n - a,
        unclassified_fraction = if (n > 0) (n - a) / n else NA_real_)
    }))
  }
  summaries <- list(rank_summary(rep(TRUE, length(ids)), "All samples (unique)", "Primary catalogue"))
  for (s in sample_names) {
    for (g in c("Abundant", "Rare")) {
      summaries[[length(summaries) + 1L]] <- rank_summary(group[, s] == g, s, g)
    }
  }
  coverage <- do.call(rbind, summaries)
  rownames(coverage) <- NULL
  stopifnot(nrow(detail) == 962L, nrow(coverage) == 36L,
    all(coverage$n_assigned + coverage$n_unclassified == coverage$n_vOTU),
    all(abs(coverage$assigned_fraction + coverage$unclassified_fraction - 1) < 1e-12))
  notes <- data.frame(Item = c("Title", "Source", "Input MD5", "Scope", "Quality results",
    "Taxonomic coverage", "Assignment definition", "Raw labels", "Rank independence",
    "Abundance groups", "Relative TPM", "Sample-specific membership", "Fractions",
    "Unclassified", "Group column filtering", "Missing values"),
    Description = c(
    "Table S3. Quality assessment, replication-strategy sensitivity and taxonomic assignment coverage.",
    "analyses/data/01-tse-construction/tse_primary_962.rds; representative-contig annotations and TPM assay.",
    unname(tools::md5sum(input_file)),
    "Primary catalogue: 962 vOTUs with representative contigs >=5 kb. The 2,488-vOTU discovery set is not used in this table.",
    "Quality counts and TPM comparisons are the existing Task 7 outputs. Predicted temperate fraction is T/(T+V). Medium-quality or better includes Complete, High-quality and Medium-quality representatives.",
    "Taxonomic coverage is reported separately from quality stratification. For each rank, the denominator is all vOTUs in the indicated catalogue or sample-specific abundance subset; these are count fractions, not TPM shares or numbers of distinct taxa.",
    "Assigned means a nonempty, non-placeholder annotation in the archived rank field. No new taxonomic classification or ICTV validation was performed.",
    "The *_raw columns preserve source labels. Missing labels and [Class]_, [Order]_, [Family]_ or [Genus]_ placeholders are Unclassified at the corresponding rank; suffixes are not transferred to another rank.",
    "Ranks are counted independently. A genus annotation does not automatically fill a missing family or order annotation.",
    "Abundant: relative TPM >=0.01; Rare: >0 and <0.001; Intermediate: >=0.001 and <0.01; Not detected: TPM=0.",
    "Each sample's relative TPM is the vOTU TPM divided by total TPM across all 962 representatives. TPM denotes length- and library-size-normalized metagenomic abundance, not gene expression.",
    "Groups are assigned separately in BS, SA, IA and DA. One vOTU can belong to different groups in different samples; sample counts must not be summed as unique catalogue counts.",
    "All fraction and relative_TPM columns are numeric on a 0-1 scale. The displayed decimals do not change stored precision.",
    "Unclassified denotes no valid assignment at that rank in these archived outputs; it is not evidence of biological novelty or absence of a taxonomic group.",
    "Filter BS_group, SA_group, IA_group or DA_group in the vOTU data sheet, then inspect the cleaned rank columns and *_assigned indicators.",
    "A missing predicted temperate fraction means T+V=0, not zero temperate abundance. Blank *_raw cells represent missing source labels."
  ), stringsAsFactors = FALSE)

  utils::write.csv(coverage, file.path(table_dir, "Table_S3_taxonomic_coverage.csv"), row.names = FALSE, na = "")
  utils::write.csv(detail, file.path(table_dir, "Table_S3_vOTU_quality_taxonomy.csv"), row.names = FALSE, na = "")
  wb <- openxlsx::createWorkbook(creator = "SwDGVirome")
  openxlsx::modifyBaseFont(wb, fontName = "Arial", fontSize = 11)
  header <- openxlsx::createStyle(fgFill = "#E8EEF3", textDecoration = "bold",
    wrapText = TRUE, valign = "center")
  fraction_style <- openxlsx::createStyle(numFmt = "0.0000")
  integer_style <- openxlsx::createStyle(numFmt = "0")
  add_sheet <- function(name, data, widths = 22) {
    openxlsx::addWorksheet(wb, name, gridLines = FALSE)
    openxlsx::writeData(wb, name, data, headerStyle = header, withFilter = TRUE, keepNA = FALSE)
    openxlsx::setColWidths(wb, name, cols = seq_along(data), widths = widths)
    openxlsx::setRowHeights(wb, name, rows = 1, heights = 58)
    openxlsx::freezePane(wb, name, firstActiveRow = 2, firstActiveCol = 2)
    numeric_cols <- which(vapply(data, is.numeric, logical(1)))
    fraction_cols <- grep("fraction|relative_TPM|_TPM$", names(data))
    for (cols in list(setdiff(numeric_cols, fraction_cols), fraction_cols)) {
      if (length(cols)) openxlsx::addStyle(wb, name,
        if (all(cols %in% fraction_cols)) fraction_style else integer_style,
        rows = seq_len(nrow(data)) + 1, cols = cols, gridExpand = TRUE)
    }
  }
  add_sheet("Quality counts", quality_counts, c(30, rep(23, ncol(quality_counts) - 1)))
  add_sheet("TPM comparison", comparison, 26)
  full_start <- nrow(comparison) + 5L
  openxlsx::writeData(wb, "TPM comparison", "Full quality-stratified sample results", startRow = full_start - 1)
  openxlsx::writeData(wb, "TPM comparison", sample_fractions, startRow = full_start,
    headerStyle = header, keepNA = FALSE)
  openxlsx::setColWidths(wb, "TPM comparison", cols = seq_along(sample_fractions), widths = 26)
  openxlsx::setColWidths(wb, "TPM comparison", cols = c(1, 2, 3, 5, 6),
    widths = c(30, 32, 40, 42, 44))
  openxlsx::setRowHeights(wb, "TPM comparison", rows = full_start, heights = 58)
  openxlsx::addStyle(wb, "TPM comparison", fraction_style,
    rows = full_start + seq_len(nrow(sample_fractions)),
    cols = grep("TPM|fraction", names(sample_fractions)), gridExpand = TRUE)
  add_sheet("Taxonomic coverage", coverage, c(26, 25, 16, rep(23, 5)))
  excel_detail <- detail
  for (nm in grep("_raw$", names(excel_detail), value = TRUE)) {
    excel_detail[[nm]][!is.na(excel_detail[[nm]]) & excel_detail[[nm]] == ""] <- NA_character_
  }
  add_sheet("vOTU data", excel_detail, 24)
  # Keep small nonzero relative abundances visible without percentage formatting.
  openxlsx::addStyle(wb, "vOTU data", openxlsx::createStyle(numFmt = "0.000000"),
    rows = seq_len(nrow(detail)) + 1,
    cols = grep("_relative_TPM$", names(detail)), gridExpand = TRUE, stack = TRUE)
  openxlsx::setColWidths(wb, "vOTU data", cols = 1, widths = 16)
  openxlsx::setColWidths(wb, "vOTU data", cols = match("representative_contig_id", names(detail)), widths = 42)
  openxlsx::setColWidths(wb, "vOTU data", cols = grep("_raw$", names(detail)), widths = 48)
  add_sheet("Notes", notes, c(28, 115))
  openxlsx::addStyle(wb, "Notes", openxlsx::createStyle(wrapText = TRUE, valign = "top"),
    rows = seq_len(nrow(notes)) + 1, cols = 1:2, gridExpand = TRUE)
  openxlsx::setRowHeights(wb, "Notes", rows = seq_len(nrow(notes)) + 1, heights = 64)
  # No drawings are used; remove openxlsx's unused links to missing drawings.
  for (i in seq_along(wb$worksheets)) {
    wb$worksheets[[i]]$drawing <- character(0)
    wb$worksheets_rels[[i]] <- character(0)
  }
  xlsx_file <- file.path(table_dir, "Table_S3_quality_taxonomy.xlsx")
  openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)
  # Check saved values, including the full 962-row filtered-detail sheet.
  for (sheet in c("Taxonomic coverage", "vOTU data")) {
    expected <- if (sheet == "Taxonomic coverage") coverage else detail
    actual <- openxlsx::read.xlsx(xlsx_file, sheet = sheet, check.names = FALSE, na.strings = character())
    for (nm in names(expected)) {
      if (is.character(expected[[nm]])) {
        expected[[nm]][is.na(expected[[nm]]) | expected[[nm]] == ""] <- ""
        actual[[nm]][is.na(actual[[nm]]) | actual[[nm]] == ""] <- ""
      }
    }
    stopifnot(isTRUE(all.equal(actual, expected, check.attributes = FALSE, tolerance = 1e-10)))
  }
  print(knitr::kable(coverage, digits = 4,
    caption = "Table S3. Taxonomic assignment coverage (count fractions)."))
  message("Verified: 962 vOTUs; 36 coverage rows; five Excel sheets.")
  message("Table S3 Excel: ", normalizePath(xlsx_file))
})
```

    ##
    ##
    ## Table: Table S3. Taxonomic assignment coverage (count fractions).
    ##
    ## |sample               |subset            |rank   | n_vOTU| n_assigned| assigned_fraction| n_unclassified| unclassified_fraction|
    ## |:--------------------|:-----------------|:------|------:|----------:|-----------------:|--------------:|---------------------:|
    ## |All samples (unique) |Primary catalogue |class  |    962|        951|            0.9886|             11|                0.0114|
    ## |All samples (unique) |Primary catalogue |order  |    962|          8|            0.0083|            954|                0.9917|
    ## |All samples (unique) |Primary catalogue |family |    962|         55|            0.0572|            907|                0.9428|
    ## |All samples (unique) |Primary catalogue |genus  |    962|         21|            0.0218|            941|                0.9782|
    ## |BS                   |Abundant          |class  |      4|          4|            1.0000|              0|                0.0000|
    ## |BS                   |Abundant          |order  |      4|          0|            0.0000|              4|                1.0000|
    ## |BS                   |Abundant          |family |      4|          0|            0.0000|              4|                1.0000|
    ## |BS                   |Abundant          |genus  |      4|          0|            0.0000|              4|                1.0000|
    ## |BS                   |Rare              |class  |    552|        544|            0.9855|              8|                0.0145|
    ## |BS                   |Rare              |order  |    552|          1|            0.0018|            551|                0.9982|
    ## |BS                   |Rare              |family |    552|         29|            0.0525|            523|                0.9475|
    ## |BS                   |Rare              |genus  |    552|          9|            0.0163|            543|                0.9837|
    ## |SA                   |Abundant          |class  |     37|         36|            0.9730|              1|                0.0270|
    ## |SA                   |Abundant          |order  |     37|          0|            0.0000|             37|                1.0000|
    ## |SA                   |Abundant          |family |     37|          3|            0.0811|             34|                0.9189|
    ## |SA                   |Abundant          |genus  |     37|          1|            0.0270|             36|                0.9730|
    ## |SA                   |Rare              |class  |    708|        704|            0.9944|              4|                0.0056|
    ## |SA                   |Rare              |order  |    708|          4|            0.0056|            704|                0.9944|
    ## |SA                   |Rare              |family |    708|         42|            0.0593|            666|                0.9407|
    ## |SA                   |Rare              |genus  |    708|         15|            0.0212|            693|                0.9788|
    ## |IA                   |Abundant          |class  |     31|         31|            1.0000|              0|                0.0000|
    ## |IA                   |Abundant          |order  |     31|          1|            0.0323|             30|                0.9677|
    ## |IA                   |Abundant          |family |     31|          1|            0.0323|             30|                0.9677|
    ## |IA                   |Abundant          |genus  |     31|          0|            0.0000|             31|                1.0000|
    ## |IA                   |Rare              |class  |    118|        116|            0.9831|              2|                0.0169|
    ## |IA                   |Rare              |order  |    118|          0|            0.0000|            118|                1.0000|
    ## |IA                   |Rare              |family |    118|          3|            0.0254|            115|                0.9746|
    ## |IA                   |Rare              |genus  |    118|          2|            0.0169|            116|                0.9831|
    ## |DA                   |Abundant          |class  |     20|         19|            0.9500|              1|                0.0500|
    ## |DA                   |Abundant          |order  |     20|          1|            0.0500|             19|                0.9500|
    ## |DA                   |Abundant          |family |     20|          1|            0.0500|             19|                0.9500|
    ## |DA                   |Abundant          |genus  |     20|          2|            0.1000|             18|                0.9000|
    ## |DA                   |Rare              |class  |    165|        164|            0.9939|              1|                0.0061|
    ## |DA                   |Rare              |order  |    165|          1|            0.0061|            164|                0.9939|
    ## |DA                   |Rare              |family |    165|          2|            0.0121|            163|                0.9879|
    ## |DA                   |Rare              |genus  |    165|          2|            0.0121|            163|                0.9879|

## Files written

``` r
projthis::proj_dir_info(
  path_target(),
  tz = "CET"
) %>%
  knitr::kable()
```

| path                | type      | size | modification_time   |
|:--------------------|:----------|-----:|:--------------------|
| data                | directory |  224 | 2026-09-19 19:35:47 |
| figures             | directory |  288 | 2026-09-19 19:35:47 |
| quality_sensitivity | directory |  288 | 2026-09-24 13:39:53 |
| tiff_panels         | directory |  256 | 2026-09-19 19:35:47 |
