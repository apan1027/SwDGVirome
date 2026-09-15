# 07-figs3-abundant-patterns
Cunli Pan, Jinlong Ru
2025-12-20

- [<span class="toc-section-number">0.1</span> Task 1: Load
  TSE](#task-1-load-tse)
- [<span class="toc-section-number">0.2</span> Figure s3a: Abundant
  Viral Taxa (Family Level Bubble
  Plot)](#figure-s3a-abundant-viral-taxa-family-level-bubble-plot)
- [<span class="toc-section-number">0.3</span> Task 2: Figure s3b:
  Predicted Bacterial Host of Abundant Viruses
  (Heatmap)](#task-2-figure-s3b-predicted-bacterial-host-of-abundant-viruses-heatmap)
- [<span class="toc-section-number">0.4</span> Export Figure S3 panels
  as TIFF](#export-figure-s3-panels-as-tiff)

**Updated: 2026-09-15 19:45:08 CET.**

The purpose of this document is to characterize the “abundant
biosphere,” analyzing the distribution patterns, persistence, and host
associations of the most prevalent viral populations.

<details class="code-fold">
<summary>Code</summary>

``` r
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(data.table)
  library(TreeSummarizedExperiment)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(forcats)
})
```

</details>

    Warning: package 'S4Vectors' was built under R version 4.5.3

    Warning: package 'Biobase' was built under R version 4.5.3

<details class="code-fold">
<summary>Code</summary>

``` r
# Load package utility functions
devtools::load_all(here::here())
```

</details>

### Task 1: Load TSE

<details class="code-fold">
<summary>Code</summary>

``` r
tse <- readRDS(
  path_source(
    "01-tse-construction",
    "tse_primary_962.rds"
  )
)

tpm_mat <- SummarizedExperiment::assay(
  tse,
  "tpm"
)

primary_votu_ids <- rownames(tpm_mat)

stopifnot(
  nrow(tpm_mat) == 962,
  ncol(tpm_mat) == 4,
  length(primary_votu_ids) == 962
)

all_metadata <- S4Vectors::metadata(tse)

# metadata(tse) still contains records from the full discovery catalogue.
# Restrict every metadata table explicitly to the 962 primary vOTUs.
contig_anno <- all_metadata$contig_annotation %>%
  dplyr::filter(
    vOTU_id %in% primary_votu_ids
  )

imgvr <- all_metadata$imgvr_source %>%
  dplyr::semi_join(
    contig_anno %>%
      dplyr::select(contig_id) %>%
      dplyr::distinct(),
    by = "contig_id"
  )

host_genome <- all_metadata$host_genome_edges %>%
  dplyr::filter(
    vOTU_id %in% primary_votu_ids
  )

stopifnot(
  all(contig_anno$vOTU_id %in% primary_votu_ids),
  all(host_genome$vOTU_id %in% primary_votu_ids)
)

message(
  "Primary TSE loaded: ",
  nrow(tpm_mat),
  " vOTUs and ",
  ncol(tpm_mat),
  " samples"
)
```

</details>

    Primary TSE loaded: 962 vOTUs and 4 samples

<details class="code-fold">
<summary>Code</summary>

``` r
message(
  "Filtered metadata: ",
  nrow(contig_anno),
  " contig records, ",
  nrow(imgvr),
  " IMG/VR records and ",
  nrow(host_genome),
  " host-prediction records"
)
```

</details>

    Filtered metadata: 967 contig records, 967 IMG/VR records and 2761 host-prediction records

### Figure s3a: Abundant Viral Taxa (Family Level Bubble Plot)

<details class="code-fold">
<summary>Code</summary>

``` r
# Extract data from TSE
tpm_mat <- assays(tse)$tpm

family_df <- data.frame(
  vOTU_id = rownames(tse),
  family = as.character(rowData(tse)$family),
  stringsAsFactors = FALSE
)

sample_df <- data.frame(
  sample_id = colnames(tse),
  sample_group = as.character(colData(tse)$sample_group),
  stringsAsFactors = FALSE
)

# Convert to long format and calculate relative abundance
abundance_long <- as.data.frame(tpm_mat) %>%
  rownames_to_column("vOTU_id") %>%
  pivot_longer(cols = -vOTU_id, names_to = "sample_id", values_to = "tpm") %>%
  dplyr::left_join(sample_df, by = "sample_id") %>%
  dplyr::left_join(family_df, by = "vOTU_id")

# Clean family names
abundance_long <- abundance_long %>%
  mutate(
    family = str_remove(family, regex("^(\\[?Family\\]?_?|f__|family_)", ignore_case = TRUE)),
    family = str_trim(family),
    family = if_else(is.na(family) | family == "", "Unclassified virus", family)
  )

# Calculate relative abundance per sample_group
sample_totals <- abundance_long %>%
  group_by(sample_group) %>%
  summarise(total_TPM = sum(tpm, na.rm = TRUE), .groups = "drop")

abundance_rel <- abundance_long %>%
  dplyr::left_join(sample_totals, by = "sample_group") %>%
  mutate(rel_abundance = tpm / total_TPM)

# Identify abundant vOTUs (≥1%)
abundant_vOTUs <- abundance_rel %>%
  dplyr::filter(rel_abundance >= 0.01) %>%
  dplyr::select(sample_group, vOTU_id, family, rel_abundance) %>%
  arrange(sample_group, desc(rel_abundance))

# Prepare plot data
df_sum <- abundant_vOTUs %>%
  group_by(family, sample_group) %>%
  summarise(total_abund = sum(rel_abundance), .groups = "drop")

# Keep top 14 families
top_families <- df_sum %>%
  group_by(family) %>%
  summarise(mean_abund = mean(total_abund), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 14) %>%
  pull(family)

df_sum <- df_sum %>%
  dplyr::filter(family %in% top_families) %>%
  mutate(
    sample_group = factor(sample_group, levels = c("BS", "SA", "IA", "DA")),
    family = fct_relevel(family, "Unclassified virus", after = 0)
  )

# Create bubble plot
p_bubble <- ggplot(df_sum, aes(x = total_abund, y = fct_rev(factor(family)), color = sample_group)) +
  geom_point(aes(size = total_abund), alpha = 0.85) +
  scale_size_continuous(
    range = c(2, 10),
    name = "Relative Abundance",
    breaks = pretty(df_sum$total_abund, n = 4)
  ) +
  scale_x_continuous(
    breaks = seq(0, 1.0, by = 0.2),
    limits = c(0, 1.0),
    labels = number_format(accuracy = 0.1)
  ) +
  scale_color_manual(
    values = c("BS" = "#E69F00", "SA" = "#56B4E9", "IA" = "#009E73", "DA" = "#F0E442"),
    name = "Sample"
  ) +
  labs(x = "Relative abundance", y = NULL) +
  theme_minimal(base_size = 22) +
  theme(
    text = element_text(family = "Times"),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.5, color = "black"),
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "italic", size = 18, color = "black"),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 20, face = "plain", color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p_bubble)
```

</details>

![](07-figs3-abundant-patterns_files/figure-commonmark/FigS3a-bubble-1.png)

<details class="code-fold">
<summary>Code</summary>

``` r
# Save
ggsave(path_target("FigS3a_bubble.png"), p_bubble, width = 8.5, height = 6, dpi = 300)
write_csv(df_sum, path_target("FigS3a_bubble_data.csv"))
write_csv(abundant_vOTUs, path_target("FigS3a_abundant_full.csv"))

message("FigS3a completed: ", nrow(abundant_vOTUs), " abundant vOTU-sample pairs")
```

</details>

    FigS3a completed: 92 abundant vOTU-sample pairs

\#\| label: summary \#\| message: true

cat(“========================================================================”)
cat(“FIGURE s3 COMPLETED!”)
cat(“========================================================================”)
cat(“✅ Figure s3a: Abundant Viral Taxa (Bubble Plot)”) cat(” -
Criterion: ≥1% relative abundance per sample group“) cat(” - Display:
Top 14 families“) cat(” - Style: Exactly matching old code“) cat(”“)
cat(”✅ Figure s3b: Ecosystem Distribution (Vertical Bar Chart)“)
cat(” - Criterion: ≥1% relative abundance per sample“) cat(” -
Denominator: Total TPM of ALL viruses“) cat(” - Style: Exactly matching
old code“)
cat(”========================================================================“)

### Task 2: Figure s3b: Predicted Bacterial Host of Abundant Viruses (Heatmap)

<details class="code-fold">
<summary>Code</summary>

``` r
host_genome <- metadata(tse)$host_genome_edges

CONFIDENCE_THRESHOLD <- 80

host_filtered <- host_genome %>%
  dplyr::filter(Confidence.score >= CONFIDENCE_THRESHOLD) %>%
  mutate(
    Phylum = str_extract(Host.taxonomy, "(?<=p__)[^;]+"),
    Phylum = if_else(is.na(Phylum), "Unclassified", Phylum)
  )

if (!exists("abundant_votus_list")) {
  sample_df <- data.frame(
    sample_id = colnames(tse),
    sample_group = as.character(colData(tse)$sample_group),
    stringsAsFactors = FALSE
  )

  tpm_mat <- assays(tse)$tpm

  abundance_long <- as.data.frame(tpm_mat) %>%
    rownames_to_column("vOTU_id") %>%
    pivot_longer(cols = -vOTU_id, names_to = "sample_id", values_to = "TPM") %>%
    dplyr::left_join(sample_df, by = "sample_id")

  sample_totals <- abundance_long %>%
    group_by(sample_group) %>%
    summarise(total_TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

  abundance_rel <- abundance_long %>%
    dplyr::left_join(sample_totals, by = "sample_group") %>%
    mutate(rel_abundance = TPM / total_TPM)

  abundant_votus_list <- abundance_rel %>%
    dplyr::filter(rel_abundance >= 0.01) %>%
    pull(vOTU_id) %>%
    unique()
}

host_abundant <- host_filtered %>%
  dplyr::filter(vOTU_id %in% abundant_votus_list)

# Calculate weights for multi-Phylum vOTUs
host_count <- host_abundant %>%
  group_by(vOTU_id) %>%
  summarise(n_hosts = dplyr::n_distinct(Phylum), .groups = "drop")

host_weighted <- host_abundant %>%
  dplyr::select(vOTU_id, Phylum) %>%
  dplyr::distinct() %>%
  dplyr::left_join(host_count, by = "vOTU_id") %>%
  mutate(weight = 1 / n_hosts)

sample_df <- data.frame(
  sample_id = colnames(tse),
  sample_group = as.character(colData(tse)$sample_group),
  stringsAsFactors = FALSE
)

tpm_abundant <- as.data.frame(assays(tse)$tpm) %>%
  rownames_to_column("vOTU_id") %>%
  dplyr::filter(vOTU_id %in% abundant_votus_list) %>%
  pivot_longer(
    cols = -vOTU_id,
    names_to = "sample_id",
    values_to = "tpm"
  ) %>%
  dplyr::left_join(sample_df, by = "sample_id") %>%
  dplyr::filter(!is.na(sample_group)) %>%
  dplyr::semi_join(
    abundant_vOTUs %>%
      dplyr::select(sample_group, vOTU_id) %>%
      dplyr::distinct(),
    by = c("sample_group", "vOTU_id")
  )

host_abundance <- host_weighted %>%
  dplyr::left_join(tpm_abundant, by = "vOTU_id", relationship = "many-to-many")

heatmap_data <- host_abundance %>%
  dplyr::filter(tpm > 0) %>%
  group_by(sample_group, Phylum) %>%
  summarise(
    n_vOTUs_weighted = sum(weight),
    presence = if_else(sum(weight) > 0, 1, 0),
    .groups = "drop"
  )

all_samples <- c("BS", "SA", "IA", "DA")
all_phyla <- unique(host_weighted$Phylum)

heatmap_data_complete <- expand.grid(
  sample_group = all_samples,
  Phylum = all_phyla,
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(heatmap_data, by = c("sample_group", "Phylum")) %>%
  mutate(
    n_vOTUs_weighted = if_else(is.na(n_vOTUs_weighted), 0, n_vOTUs_weighted),
    presence = if_else(is.na(presence), 0, presence)
  )

phylum_order <- heatmap_data_complete %>%
  group_by(Phylum) %>%
  summarise(total_presence = sum(presence), .groups = "drop") %>%
  arrange(desc(total_presence)) %>%
  pull(Phylum)

heatmap_data_complete <- heatmap_data_complete %>%
  mutate(
    sample_group = factor(sample_group, levels = c("BS", "SA", "IA", "DA")),
    Phylum = factor(Phylum, levels = phylum_order)
  )

p_heatmap <- ggplot(heatmap_data_complete, aes(x = sample_group, y = Phylum)) +
  geom_tile(aes(fill = factor(presence)), color = "white", linewidth = 1) +
  scale_fill_manual(values = c("0" = "gray90", "1" = "#FC8D62"), guide = "none") +
  theme_minimal(base_size = 16) +
  theme(
    text = element_text(family = "Times"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18, color = "black"),
    axis.text.y = element_text(color = "black", size = 18, face = "italic"),
    axis.title = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p_heatmap)
```

</details>

![](07-figs3-abundant-patterns_files/figure-commonmark/FigS3c-host-heatmap-1.png)

<details class="code-fold">
<summary>Code</summary>

``` r
ggsave(path_target("FigS3b_host_heatmap.png"), p_heatmap, width = 8, height = 3.0, dpi = 300)
write_csv(heatmap_data_complete, path_target("FigS3b_host_heatmap_data.csv"))

summary_data <- heatmap_data_complete %>%
  group_by(sample_group) %>%
  summarise(
    n_phyla_present = sum(presence > 0),
    n_vOTUs_weighted = sum(n_vOTUs_weighted),
    .groups = "drop"
  )

write_csv(summary_data, path_target("FigS3c_summary.csv"))

message("FigS3b completed: ", dplyr::n_distinct(host_weighted$vOTU_id),
        " abundant vOTUs with host predictions (weighted)")
```

</details>

    FigS3b completed: 40 abundant vOTUs with host predictions (weighted)

### Export Figure S3 panels as TIFF
