# Usage:
# Rscript analyses/R/export_table_s3_amg.R ANALYSES_DIR BASE_TABLE_S3.xlsx OUTPUT.xlsx
# Regenerate the complete table, preserving the author's formatting and Notes.
# The base workbook is never overwritten. No Spark, API or network is required.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop(
  "Usage: Rscript export_table_s3_amg.R ANALYSES_DIR BASE_TABLE_S3.xlsx OUTPUT.xlsx")
analysis_dir <- normalizePath(args[1], mustWork = TRUE)
base_xlsx <- normalizePath(args[2], mustWork = TRUE)
dir.create(dirname(args[3]), recursive = TRUE, showWarnings = FALSE)
output_xlsx <- file.path(normalizePath(dirname(args[3])), basename(args[3]))
if (base_xlsx == output_xlsx) stop("Use a separate output path; the base workbook is preserved.")

# Execute the current QMD chunks instead of maintaining a second calculation.
# Only the Table S3 chunks run; Figure S2 and the other figures are not rerendered.
qmd <- readLines(file.path(analysis_dir, "06-figs2-lifestyle.qmd"), warn = FALSE)
table_chunk <- function(label) {
  label_line <- which(trimws(qmd) == paste0("#| label: ", label))
  if (length(label_line) != 1L) stop("Expected one Table S3 chunk: ", label)
  ends <- which(seq_along(qmd) > label_line & trimws(qmd) == "```")
  if (!length(ends)) stop("Unclosed chunk: ", label)
  parse(text = qmd[seq.int(label_line + 1L, min(ends) - 1L)])
}
suppressPackageStartupMessages(library(TreeSummarizedExperiment))
setwd(analysis_dir)
here::i_am("06-figs2-lifestyle.qmd")
run_env <- new.env(parent = globalenv())
run_env$path_target <- function(...) file.path(analysis_dir, "data/06-figs2-lifestyle", ...)
run_env$path_source <- function(...) file.path(analysis_dir, "data", ...)
for (label in c("table-s3-quality-sensitivity", "table-s3-taxonomic-coverage-export")) {
  eval(table_chunk(label), envir = run_env)
}
source(file.path(analysis_dir, "R", "amg_curation.R"), local = TRUE)
amg <- build_amg_evidence(analysis_dir)
append_amg_to_s3(base_xlsx, output_xlsx, amg,
  table_dir = file.path(analysis_dir, "data/06-figs2-lifestyle/quality_sensitivity"),
  replace_existing = TRUE)
message("Regenerated and verified author-formatted Table S3: ", output_xlsx)
