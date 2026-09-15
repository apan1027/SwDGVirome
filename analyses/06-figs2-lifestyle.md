# 06-figs2-lifestyle

2026-09-15

**Updated: 2026-09-15 15:11:59 CET.**

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

> [!NOTE]
>
> ### Example structure
>
> This is an example structure to get you started. You can modify it as
> per your needs.

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

\###Task 3: Calculate predicted temperate percentages

``` r
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
    pct_temperate =
      100 * n_temperate / n_vOTU,
    .groups = "drop"
  )

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
    pct_temperate =
      100 * n_temperate / n_vOTU,
    .groups = "drop"
  )

print(length_summary, n = Inf)
```

    # A tibble: 4 × 4
      length_group n_vOTU n_temperate pct_temperate
      <fct>         <int>       <int>         <dbl>
    1 3–5 kb         1526          12         0.786
    2 5–10 kb         689           9         1.31
    3 10–20 kb        207           7         3.38
    4 ≥20 kb           66           7        10.6

``` r
print(quality_summary, n = Inf)
```

    # A tibble: 5 × 4
      quality_group  n_vOTU n_temperate pct_temperate
      <fct>           <int>       <int>         <dbl>
    1 Complete            4           0          0
    2 High-quality       12           2         16.7
    3 Medium-quality     39           4         10.3
    4 Low-quality      2311          25          1.08
    5 Not-determined    122           4          3.28

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
# predicted temperate percentage by representative contig length
p_figs2a <- ggplot(
  length_summary,
  aes(
    x = length_group,
    y = pct_temperate
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
    limits = c(0, 20),
    breaks = seq(0, 20, 5),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
   x = NULL,
    y = "vOTUs predicted temperate (%)"
  ) +
  panel_theme

# Figure S2b:
# predicted temperate percentage by CheckV quality category
p_figs2b <- ggplot(
  quality_summary,
  aes(
    x = quality_group,
    y = pct_temperate
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
    limits = c(0, 18.5),
    breaks = seq(0, 15, 5),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
   x = NULL,
    y = "vOTUs predicted temperate (%)"
  ) +
  panel_theme +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

# Combined preview only.
# You can still use the two individual TIFF panels in Inkscape.
p_figs2_combined <-
  p_figs2a +
  p_figs2b +
  patchwork::plot_annotation(
    tag_levels = "a"
  )

print(p_figs2_combined)
```

![](06-figs2-lifestyle_files/figure-commonmark/task5-create-figure-1.png)

### Task 6: Save source data, PNG and TIFF files

## Files written

``` r
projthis::proj_dir_info(
  path_target(),
  tz = "CET"
) %>%
  knitr::kable()
```

| path        | type      | size | modification_time   |
|:------------|:----------|-----:|:--------------------|
| data        | directory |  128 | 2026-09-15 15:12:05 |
| figures     | directory |  160 | 2026-09-15 15:12:06 |
| tiff_panels | directory |  160 | 2026-09-15 15:12:06 |
