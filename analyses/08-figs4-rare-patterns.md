# 08-figs4-rare-patterns
Cunli Pan, Jinlong Ru
2025-12-20

- [<span class="toc-section-number">0.1</span> Task 1: Load TSE and
  Extract Metadata](#task-1-load-tse-and-extract-metadata)
- [<span class="toc-section-number">0.2</span> Task 2: Figure S4a -
  Bubble Plot (Rare)](#task-2-figure-s4a---bubble-plot-rare)
- [<span class="toc-section-number">0.3</span> Task 4: Figure S4b - Host
  Heatmap](#task-4-figure-s4b---host-heatmap)
- [<span class="toc-section-number">0.4</span> Export Figure S4 panels
  as TIFF](#export-figure-s4-panels-as-tiff)

**Updated: 2026-09-07 17:29:01 CET.**

The purpose of this document is to investigate the “rare biosphere,”
examining the diversity, ecological drivers, and potential persistence
of low-abundance viral populations.

<details class="code-fold">
<summary>Code</summary>

``` r
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(TreeSummarizedExperiment)
  library(ggplot2)
  library(patchwork)
  library(scales)
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

### Task 1: Load TSE and Extract Metadata

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

### Task 2: Figure S4a - Bubble Plot (Rare)

<details class="code-fold">
<summary>Code</summary>

``` r
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

abundance_long <- as.data.frame(tpm_mat) %>%
  rownames_to_column("vOTU_id") %>%
  pivot_longer(cols = -vOTU_id, names_to = "sample_id", values_to = "tpm") %>%
  dplyr::left_join(sample_df, by = "sample_id") %>%
  dplyr::left_join(family_df, by = "vOTU_id")

abundance_long <- abundance_long %>%
  mutate(
    family = str_remove(family, regex("^(\\[?Family\\]?_?|f__|family_)", ignore_case = TRUE)),
    family = str_trim(family),
    family = if_else(is.na(family) | family == "", "Unclassified virus", family)
  )

sample_totals <- abundance_long %>%
  group_by(sample_group) %>%
  summarise(total_TPM = sum(tpm, na.rm = TRUE), .groups = "drop")

abundance_rel <- abundance_long %>%
  dplyr::left_join(sample_totals, by = "sample_group") %>%
  mutate(rel_abundance = tpm / total_TPM)

# Identify RARE vOTUs (< 0.1% = < 0.001)
rare_vOTUs <- abundance_rel %>%
  dplyr::filter(rel_abundance < 0.001 & rel_abundance > 0) %>%
  dplyr::select(sample_group, vOTU_id, family, rel_abundance) %>%
  arrange(sample_group, desc(rel_abundance))

df_sum <- rare_vOTUs %>%
  group_by(family, sample_group) %>%
  summarise(total_abund = sum(rel_abundance), .groups = "drop")

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
    family = forcats::fct_relevel(family, "Unclassified virus", after = 0)
  )

p_bubble <- ggplot(df_sum, aes(x = total_abund, y = forcats::fct_rev(factor(family)), color = sample_group)) +
  geom_point(aes(size = total_abund), alpha = 0.85) +
  scale_size_continuous(
    range = c(2, 10),
    name = "Relative Abundance",
    breaks = pretty(df_sum$total_abund, n = 4)
  ) +
  scale_x_continuous(
    breaks = seq(0, 0.8, by = 0.2),
    limits = c(0, max(df_sum$total_abund) * 1.1),
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

![](08-figs4-rare-patterns_files/figure-commonmark/figs4a-bubble-1.png)

<details class="code-fold">
<summary>Code</summary>

``` r
ggsave(path_target("FigS4a_bubble.png"), p_bubble, width = 8.5, height = 6, dpi = 300)
write_csv(df_sum, path_target("FigS4a_bubble_data.csv"))
write_csv(rare_vOTUs, path_target("FigS4a_rare_full.csv"))

message("FigS4a completed: ", nrow(rare_vOTUs), " rare vOTU-sample pairs")
```

</details>

    FigS4a completed: 1543 rare vOTU-sample pairs

### Task 4: Figure S4b - Host Heatmap

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

# Use rare vOTU list (from S4a/S4b)
if (!exists("rare_votus_list")) {
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

  rare_votus_list <- abundance_rel %>%
    dplyr::filter(rel_abundance < 0.001 & rel_abundance > 0) %>%
    pull(vOTU_id) %>%
    unique()
}

# Filter host predictions for rare vOTUs only
host_rare <- host_filtered %>%
  dplyr::filter(vOTU_id %in% rare_votus_list)

# ============================================================================
# NEW: Calculate weights for vOTUs with multiple Phylum predictions
# ============================================================================

# Step 1: Count how many different Phyla each vOTU has
host_count <- host_rare %>%
  dplyr::select(vOTU_id, Phylum) %>%
  dplyr::distinct() %>%
  group_by(vOTU_id) %>%
  summarise(n_hosts = dplyr::n_distinct(Phylum), .groups = "drop")

# Step 2: Calculate weights (1 / n_hosts)
host_weighted <- host_rare %>%
  dplyr::select(vOTU_id, Phylum) %>%
  dplyr::distinct() %>%
  dplyr::left_join(host_count, by = "vOTU_id") %>%
  mutate(weight = 1 / n_hosts)

# Diagnostic: Check multi-Phylum vOTUs
multi_phylum_votus <- host_count %>% dplyr::filter(n_hosts > 1)
if (nrow(multi_phylum_votus) > 0) {
  message(sprintf("  %d rare vOTUs have predictions for multiple Phyla (weighted)", nrow(multi_phylum_votus)))
} else {
  message("  All rare vOTUs have predictions for only one Phylum")
}
```

</details>

      17 rare vOTUs have predictions for multiple Phyla (weighted)

<details class="code-fold">
<summary>Code</summary>

``` r
# ============================================================================
# Prepare abundance data
# ============================================================================

sample_df <- data.frame(
  sample_id = colnames(tse),
  sample_group = as.character(colData(tse)$sample_group),
  stringsAsFactors = FALSE
)

tpm_rare <- as.data.frame(assays(tse)$tpm) %>%
  rownames_to_column("vOTU_id") %>%
  dplyr::filter(vOTU_id %in% rare_votus_list) %>%
  pivot_longer(
    cols = -vOTU_id,
    names_to = "sample_id",
    values_to = "tpm"
  ) %>%
  dplyr::left_join(sample_df, by = "sample_id") %>%
  dplyr::filter(!is.na(sample_group)) %>%
  dplyr::semi_join(
    rare_vOTUs %>%
      dplyr::select(sample_group, vOTU_id) %>%
      dplyr::distinct(),
    by = c("sample_group", "vOTU_id")
  )

# ============================================================================
# Create heatmap data (using WEIGHTED vOTU counts)
# ============================================================================

# Join with weighted host predictions
host_abundance <- host_weighted %>%  # ← Changed: use host_weighted
  dplyr::left_join(tpm_rare, by = "vOTU_id", relationship = "many-to-many")

# Calculate weighted vOTU counts per Phylum
heatmap_data <- host_abundance %>%
  dplyr::filter(tpm > 0) %>%
  group_by(sample_group, Phylum) %>%
  summarise(
    n_vOTUs_weighted = sum(weight),  # ← Changed: use weighted sum
    presence = if_else(sum(weight) > 0, 1, 0),  # ← Changed: based on weighted sum
    .groups = "drop"
  )

# Complete the grid (only rare vOTUs' Phyla)
all_samples <- c("BS", "SA", "IA", "DA")
all_phyla <- unique(host_weighted$Phylum)  # ← Changed: use host_weighted

heatmap_data_complete <- expand.grid(
  sample_group = all_samples,
  Phylum = all_phyla,
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(heatmap_data, by = c("sample_group", "Phylum")) %>%
  mutate(
    n_vOTUs_weighted = if_else(is.na(n_vOTUs_weighted), 0, n_vOTUs_weighted),  # ← Changed
    presence = if_else(is.na(presence), 0, presence)
  )

# Order phyla by total presence
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

# Create compact heatmap (horizontal layout)
p_heatmap <- ggplot(heatmap_data_complete, aes(x = sample_group, y = Phylum)) +
  geom_tile(aes(fill = factor(presence)), color = "white", linewidth = 1) +
  scale_fill_manual(
    values = c("0" = "gray90", "1" = "#66C2A5"),
    guide = "none"
  ) +
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

![](08-figs4-rare-patterns_files/figure-commonmark/figs4b-host-heatmap-1.png)

<details class="code-fold">
<summary>Code</summary>

``` r
ggsave(path_target("FigS4b_host_heatmap.png"), p_heatmap, width = 8, height = 10, dpi = 300)
write_csv(heatmap_data_complete, path_target("FigS4b_host_heatmap_data.csv"))

summary_data <- heatmap_data_complete %>%
  group_by(sample_group) %>%
  summarise(
    n_phyla_present = sum(presence > 0),
    n_vOTUs_weighted = sum(n_vOTUs_weighted),  # ← Changed: use weighted sum
    .groups = "drop"
  )

write_csv(summary_data, path_target("FigS4b_summary.csv"))

message("FigS4b completed: ", dplyr::n_distinct(host_weighted$vOTU_id),
        " rare vOTUs with host predictions (weighted)")
```

</details>

    FigS4b completed: 246 rare vOTUs with host predictions (weighted)

### Export Figure S4 panels as TIFF
