#' @title Swedish Water Virome Utility Functions
#' @description Utility functions for the Swedish Water Virome analysis pipeline.
#' @name utils
NULL

# =============================================================================
# Data Transformation Functions
# =============================================================================

#' Convert Data Frame to Matrix
#'
#' Converts a data frame to a numeric matrix with row names from an ID column.
#' Missing values are replaced with 0.
#'
#' @param df A data frame to convert.
#' @param id_col Character. Name of the column to use as row names.
#'   Default: "vOTU_id".
#'
#' @return A numeric matrix with row names set from `id_col`.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(vOTU_id = c("v1", "v2"), sample1 = c(10, 20))
#' mat <- safe_to_matrix(df)
#' }
#'
#' @export
safe_to_matrix <- function(df, id_col = "vOTU_id") {
  stopifnot(id_col %in% names(df))
  mat <- df |>
    tibble::as_tibble() |>
    tibble::column_to_rownames(id_col) |>
    as.matrix()
  storage.mode(mat) <- "double"
  mat[is.na(mat)] <- 0
  mat
}


#' Extract Level 5 from Ecosystem Classification String
#'
#' Parses a semicolon-separated ecosystem classification string and
#' returns the 5th level (position 5).
#'
#' @param eco_string Character. A semicolon-delimited ecosystem string
#'   (e.g., "Root;Level1;Level2;Level3;Level4;Level5").
#'
#' @return Character. The 5th level of the classification, or `NA_character_`
#'   if not available.
#'
#' @examples
#' extract_level5_position5("A;B;C;D;E;F")
#' # Returns "E"
#'
#' @export
extract_level5_position5 <- function(eco_string) {
  if (is.na(eco_string) || eco_string == "") {
    return(NA_character_)
  }
  parts <- strsplit(as.character(eco_string), ";")[[1]]
  parts <- trimws(parts)
  if (length(parts) >= 5 && parts[5] != "") {
    return(parts[5])
  } else {
    return(NA_character_)
  }
}


# =============================================================================
# Visualization Theme Functions
# =============================================================================

#' Custom ggplot2 Theme for Figure 2
#'
#' A consistent minimal theme for Figure 2 panels with Times font family.
#'
#' @return A ggplot2 theme object.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_fig2()
#' }
#'
#' @export
theme_fig2 <- function() {
  ggplot2::theme_minimal(base_size = 22) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Times"),
      axis.line = ggplot2::element_line(color = "black", linewidth = 1),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.8),
      axis.text = ggplot2::element_text(color = "black", size = 20, family = "Times"),
      axis.title = ggplot2::element_text(size = 23, face = "plain", family = "Times"),
      axis.title.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10, family = "Times"),
      legend.title = ggplot2::element_text(size = 12, face = "bold", family = "Times"),
      plot.margin = ggplot2::margin(15, 20, 10, 10)
    )
}


# =============================================================================
# Diversity Analysis Functions
# =============================================================================

#' Calculate Rarefaction Curves
#'
#' Computes rarefaction curves for diversity metrics (observed richness,
#' Chao1 richness, Shannon diversity) across multiple sequencing depths.
#'
#' @param count_mat A numeric matrix of counts (rows = features, columns = samples).
#' @param n_points Integer. Number of sampling depths to evaluate. Default: 20.
#' @param max_depth Integer. Maximum sampling depth. Default: 50000.
#'
#' @return A data frame with columns: `sample_id`, `sequencing_depth`,
#'   `observed_richness`, `chao1_richness`, `shannon_diversity`.
#'
#' @examples
#' \dontrun{
#' rarefaction_data <- calculate_rarefaction(count_mat, n_points = 20)
#' }
#'
#' @export
calculate_rarefaction <- function(count_mat, n_points = 20, max_depth = 50000) {
  results <- list()

  for (i in seq_len(ncol(count_mat))) {
    sample_id <- colnames(count_mat)[i]
    sample_counts <- count_mat[, i]
    sample_counts <- sample_counts[sample_counts > 0]

    if (length(sample_counts) == 0) next

    total_reads <- sum(sample_counts)
    if (total_reads < 100) next

    max_sampling_depth <- min(total_reads, max_depth)

    # Use logarithmic spacing (original version)
    if (max_sampling_depth <= 1000) {
      depths <- seq(100, max_sampling_depth, length.out = n_points)
    } else {
      depths <- unique(round(exp(seq(log(100), log(max_sampling_depth), length.out = n_points))))
      if (!max_sampling_depth %in% depths) {
        depths <- c(depths, max_sampling_depth)
      }
    }

    # Calculate diversity metrics at each depth
    chao1_values <- numeric(length(depths))
    shannon_values <- numeric(length(depths))
    observed_richness <- numeric(length(depths))

    for (j in 1:length(depths)) {
      depth <- depths[j]
      rarefied_sample <- vegan::rrarefy(sample_counts, depth)

      observed_richness[j] <- sum(rarefied_sample > 0)
      chao1_estimate <- vegan::estimateR(rarefied_sample)
      chao1_values[j] <- chao1_estimate[2]
      shannon_values[j] <- vegan::diversity(rarefied_sample, index = "shannon")
    }

    results[[sample_id]] <- data.frame(
      sample_id = sample_id,
      sequencing_depth = depths,
      observed_richness = observed_richness,
      chao1_richness = chao1_values,
      shannon_diversity = shannon_values,
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(results)
}


#' Get Viruses by Sample Group
#'
#' Returns virus IDs present (TPM > 0) in a specified sample group.
#'
#' @param tpm_mat A numeric matrix of TPM values (rows = vOTUs, columns = samples).
#' @param sample_groups A named character vector mapping sample IDs to group names.
#' @param group_name Character. The target group name to filter.
#'
#' @return Character vector of vOTU IDs present in the specified group.
#'
#' @examples
#' \dontrun{
#' bs_viruses <- get_viruses_by_group(tpm_mat, sample_groups, "BS")
#' }
#'
#' @export
get_viruses_by_group <- function(tpm_mat, sample_groups, group_name) {
  samples <- names(sample_groups)[sample_groups == group_name]
  group_tpm <- tpm_mat[, samples, drop = FALSE]
  rownames(group_tpm)[rowSums(group_tpm) > 0]
}


# =============================================================================
# KEGG Annotation Functions
# =============================================================================

#' Get KEGG Information for a Database ID
#'
#' Queries the KEGG database for information about a KO number or EC number.
#'
#' @param id Character. A KEGG KO number (e.g., "K00001") or EC number
#'   (e.g., "1.1.1.1").
#'
#' @return A tibble with columns: `dbid`, `db_type`, `Name`, `Definition`,
#'   `Pathway`, `Module`, `Class`.
#'
#' @examples
#' \dontrun{
#' get_kegg_info("K00001")
#' }
#'
#' @export
get_kegg_info <- function(id) {
  # KO number (Kxxxxx)
  if (stringr::str_detect(id, "^K\\d{5}")) {
    entry <- tryCatch(KEGGREST::keggGet(id)[[1]], error = function(e) NULL)
    if (!is.null(entry)) {
      return(tibble::tibble(
        dbid = id,
        db_type = "KEGG_KO",
        Name = entry$NAME[1] %||% NA_character_,
        Definition = entry$DEFINITION %||% NA_character_,
        Pathway = if (!is.null(entry$PATHWAY)) paste(names(entry$PATHWAY), collapse = "; ") else NA_character_,
        Module = if (!is.null(entry$MODULE)) paste(names(entry$MODULE), collapse = "; ") else NA_character_,
        Class = if (!is.null(entry$CLASS)) paste(entry$CLASS, collapse = "; ") else NA_character_
      ))
    }
  }

  # EC number (x.x.x.x)
  if (stringr::str_detect(id, "^\\d+\\.\\d+\\.\\d+\\.\\d+")) {
    entry <- tryCatch(KEGGREST::keggGet(paste0("ec:", id))[[1]], error = function(e) NULL)
    if (!is.null(entry)) {
      return(tibble::tibble(
        dbid = id,
        db_type = "EC_Number",
        Name = entry$NAME[1] %||% NA_character_,
        Definition = entry$DEFINITION %||% NA_character_,
        Pathway = if (!is.null(entry$PATHWAY)) paste(names(entry$PATHWAY), collapse = "; ") else NA_character_,
        Module = NA_character_,
        Class = NA_character_
      ))
    }
  }

  # Unknown ID
  tibble::tibble(
    dbid = id,
    db_type = "Other_DB",
    Name = NA_character_,
    Definition = NA_character_,
    Pathway = NA_character_,
    Module = NA_character_,
    Class = NA_character_
  )
}


#' Get KEGG Pathway Information
#'
#' Queries KEGG for pathway details and classifies into top-level categories.
#'
#' @param map_id Character. A KEGG map ID (e.g., "map00010").
#' @param index Integer. Current index for progress messages.
#' @param total Integer. Total number of pathways for progress messages.
#'
#' @return A tibble with columns: `Pathway_ID`, `Pathway_Name`, `Pathway_Class`,
#'   `Pathway_Top_Category`.
#'
#' @examples
#' \dontrun{
#' get_pathway_info("map00010", 1, 100)
#' }
#'
#' @export
get_pathway_info <- function(map_id, index, total) {
  message(sprintf("[%d/%d] Querying: %s", index, total, map_id))

  tryCatch({
    entry <- KEGGREST::keggGet(map_id)[[1]]
    name <- if (!is.null(entry$NAME)) paste(entry$NAME, collapse = "; ") else NA
    full_class <- if (!is.null(entry$CLASS)) paste(entry$CLASS, collapse = "; ") else NA

    # Seven major KEGG categories with hardcoded 1.0 Global maps
    top_class <- dplyr::case_when(
      # Global and overview maps (all belong to Metabolism)
      map_id %in% c(
        "map01100", "map01110", "map01120", "map01200", "map01210",
        "map01212", "map01230", "map01232", "map01250", "map01240",
        "map01220", "map01310", "map01320"
      ) ~ "Metabolism",
      # Class-based matching
      stringr::str_detect(full_class, "Metabolism") ~ "Metabolism",
      stringr::str_detect(full_class, "Genetic Information") ~ "Genetic Information Processing",
      stringr::str_detect(full_class, "Environmental Information") ~ "Environmental Information Processing",
      stringr::str_detect(full_class, "Cellular Processes") ~ "Cellular Processes",
      stringr::str_detect(full_class, "Organismal Systems") ~ "Organismal Systems",
      stringr::str_detect(full_class, "Human Diseases") ~ "Human Diseases",
      stringr::str_detect(full_class, "Drug Development") ~ "Drug Development",
      TRUE ~ "Other"
    )

    message(sprintf("  Success: %s -> %s", name, top_class))

    tibble::tibble(
      Pathway_ID = map_id,
      Pathway_Name = name,
      Pathway_Class = full_class,
      Pathway_Top_Category = top_class
    )

  }, error = function(e) {
    # Fallback: classify by pathway ID prefix
    prefix <- substr(map_id, 4, 5)
    top_class <- dplyr::case_when(
      prefix == "01" ~ "Metabolism",
      prefix == "02" ~ "Genetic Information Processing",
      prefix == "03" ~ "Environmental Information Processing",
      prefix == "04" ~ "Cellular Processes",
      prefix == "05" ~ "Organismal Systems",
      prefix == "06" ~ "Human Diseases",
      prefix == "07" ~ "Drug Development",
      TRUE ~ "Other"
    )

    message(sprintf("  Failed: %s -> Fallback: %s", e$message, top_class))

    tibble::tibble(
      Pathway_ID = map_id,
      Pathway_Name = paste("Unknown", map_id),
      Pathway_Class = NA_character_,
      Pathway_Top_Category = top_class
    )
  })
}


# =============================================================================
# Host Composition Functions
# =============================================================================

#' Calculate Host Relative Abundance
#'
#' Calculates the relative abundance of predicted bacterial hosts at a
#' specified taxonomic level, relative to total viral TPM.
#'
#' @param host_filtered A data frame with columns `vOTU_id` and the specified
#'   taxonomic level.
#' @param tpm_long A long-format data frame with columns `vOTU_id`, `sample_id`, `TPM`.
#' @param sample_total_tpm A data frame with columns `sample_id`, `total_TPM`.
#' @param sample_metadata A data frame with sample metadata including `sample_id`.
#' @param tax_level Character. Taxonomic level to summarize. Default: "Phylum".
#'
#' @return A data frame with relative abundance per sample and taxonomic group.
#'
#' @export
calculate_host_abundance <- function(host_filtered, tpm_long, sample_total_tpm,
                                     sample_metadata, tax_level = "Phylum") {

  host_count <- host_filtered |>
    dplyr::group_by(.data$vOTU_id) |>
    dplyr::summarise(n_hosts = dplyr::n_distinct(.data[[tax_level]]), .groups = "drop")

  host_weighted <- host_filtered |>
    dplyr::left_join(host_count, by = "vOTU_id") |>
    dplyr::mutate(weight = 1 / .data$n_hosts) |>
    dplyr::select("vOTU_id", dplyr::all_of(tax_level), "weight") |>
    dplyr::distinct()

  host_abundance <- tpm_long |>
    dplyr::inner_join(host_weighted, by = "vOTU_id", relationship = "many-to-many") |>
    dplyr::mutate(weighted_TPM = .data$TPM * .data$weight)

  tax_abundance <- host_abundance |>
    dplyr::group_by(.data$sample_id, .data[[tax_level]]) |>
    dplyr::summarise(tax_TPM = sum(.data$weighted_TPM), .groups = "drop")

  tax_rel_abundance <- tax_abundance |>
    dplyr::left_join(sample_total_tpm, by = "sample_id") |>
    dplyr::mutate(rel_abundance = .data$tax_TPM / .data$total_TPM) |>
    dplyr::select("sample_id", dplyr::all_of(tax_level), "rel_abundance") |>
    dplyr::left_join(sample_metadata, by = "sample_id")

  tax_rel_abundance
}


#' Calculate Normalized Host Abundance
#'
#' Calculates 100% normalized relative abundance of predicted bacterial hosts
#' at a specified taxonomic level.
#'
#' @inheritParams calculate_host_abundance
#'
#' @return A data frame with 100% normalized relative abundance per sample.
#'
#' @export
calculate_host_abundance_normalized <- function(host_filtered, tpm_long,
                                                sample_metadata, tax_level = "Phylum") {

  host_count <- host_filtered |>
    dplyr::group_by(.data$vOTU_id) |>
    dplyr::summarise(n_hosts = dplyr::n_distinct(.data[[tax_level]]), .groups = "drop")

  host_weighted <- host_filtered |>
    dplyr::left_join(host_count, by = "vOTU_id") |>
    dplyr::mutate(weight = 1 / .data$n_hosts) |>
    dplyr::select("vOTU_id", dplyr::all_of(tax_level), "weight") |>
    dplyr::distinct()

  host_abundance <- tpm_long |>
    dplyr::inner_join(host_weighted, by = "vOTU_id", relationship = "many-to-many") |>
    dplyr::mutate(weighted_TPM = .data$TPM * .data$weight)

  tax_abundance <- host_abundance |>
    dplyr::group_by(.data$sample_id, .data[[tax_level]]) |>
    dplyr::summarise(tax_TPM = sum(.data$weighted_TPM), .groups = "drop")

  tax_rel_abundance <- tax_abundance |>
    dplyr::group_by(.data$sample_id) |>
    dplyr::mutate(
      sample_assigned_total = sum(.data$tax_TPM),
      rel_abundance = .data$tax_TPM / .data$sample_assigned_total
    ) |>
    dplyr::ungroup() |>
    dplyr::select("sample_id", dplyr::all_of(tax_level), "rel_abundance") |>
    dplyr::left_join(sample_metadata, by = "sample_id")

  tax_rel_abundance
}


#' Create Stacked Bar Chart for Host Composition
#'
#' Generates a publication-ready stacked bar chart showing host composition
#' at the specified taxonomic level.
#'
#' @param tax_rel_abundance A data frame with columns `sample_group`,
#'   the specified taxonomic level, and `rel_abundance`.
#' @param tax_level Character. Taxonomic level for display. Default: "Phylum".
#' @param top_n Integer. Number of top taxa to display separately. Default: 10.
#' @param normalize Logical. If TRUE, scale y-axis to data range. Default: FALSE.
#'
#' @return A ggplot object.
#'
#' @export
create_stacked_bar <- function(tax_rel_abundance, tax_level = "Phylum",
                               top_n = 10, normalize = FALSE) {

  if (tax_level == "Class") top_n <- 15

  tax_totals <- tax_rel_abundance |>
    dplyr::group_by(.data[[tax_level]]) |>
    dplyr::summarise(total_abund = sum(.data$rel_abundance), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$total_abund))

  top_taxa <- tax_totals |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(.data[[tax_level]])

  plot_data <- tax_rel_abundance |>
    dplyr::mutate(
      Taxa_display = dplyr::if_else(.data[[tax_level]] %in% top_taxa,
                                    .data[[tax_level]], "Others")
    ) |>
    dplyr::group_by(.data$sample_group, .data$Taxa_display) |>
    dplyr::summarise(rel_abundance = sum(.data$rel_abundance), .groups = "drop") |>
    dplyr::mutate(
      Taxa_display = factor(.data$Taxa_display, levels = rev(c(top_taxa, "Others"))),
      sample_group = factor(.data$sample_group, levels = c("BS", "SA", "IA", "DA"))
    )

  n_colors <- length(top_taxa) + 1

  if (tax_level == "Phylum") {
    colors <- if (n_colors <= 12) {
      RColorBrewer::brewer.pal(max(3, n_colors), "Set3")
    } else {
      grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_colors)
    }
  } else {
    colors <- if (n_colors <= 12) {
      RColorBrewer::brewer.pal(max(3, n_colors), "Paired")
    } else {
      grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_colors)
    }
  }

  names(colors) <- c(top_taxa, "Others")
  colors["Others"] <- "#CCCCCC"

  ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data$sample_group,
    y = .data$rel_abundance,
    fill = .data$Taxa_display
  )) +
    ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.7, color = NA) +
    ggplot2::scale_fill_manual(values = colors, name = tax_level) +
    ggplot2::scale_y_continuous(
      labels = function(x) sprintf("%.1f", x),
      limits = if (!normalize) c(0, 1) else NULL,
      breaks = if (!normalize) seq(0, 1, 0.2) else ggplot2::waiver(),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(x = NULL, y = "Relative abundance of predicted bacterial hosts") +
    ggplot2::theme_minimal(base_size = 22) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Times"),
      axis.text = ggplot2::element_text(color = "black", size = 20),
      axis.title.y = ggplot2::element_text(size = 20, face = "plain"),
      axis.title.x = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black", linewidth = 1),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.8),
      legend.position = "right",
      legend.key.size = ggplot2::unit(0.4, "cm"),
      legend.text = ggplot2::element_text(size = 10, face = "italic"),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 10, 10, 10),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}
