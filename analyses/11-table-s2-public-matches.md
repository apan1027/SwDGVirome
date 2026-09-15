# 11-table-s2-public-matches

2026-09-15

**Updated: 2026-09-15 15:12:16 CET.**

``` r
suppressPackageStartupMessages({
  library(here)
  library(conflicted)
  library(tidyverse)
  library(data.table)
  devtools::load_all()
})
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Please install openxlsx before rendering this QMD.")
}
```

## Tasks

### Task 1: Import existing sequence-match records

This step exports archived alignment results, without rerunning
alignments, read mapping, or cross-study species dereplication.

``` r
input_file <- file.path(path_data, "public_reference_alignments_primary962.csv")
if (!file.exists(input_file)) stop("Missing input: ", input_file)
matches <- utils::read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("vOTU_id", "reference_accession", "qlen", "slen", "ani",
              "af_query", "af_target", "representative_length_bp", "in_primary_962")
stopifnot(all(required %in% names(matches)))
stopifnot(
  !anyNA(matches[required]),
  nrow(matches) == 9L,
  length(unique(matches$vOTU_id)) == 9L,
  length(unique(matches$reference_accession)) == 6L,
  all(matches$in_primary_962),
  all(matches$qlen >= 5000),
  all(matches$qlen == matches$representative_length_bp),
  all(matches$slen > 0),
  all(matches$ani >= 95 & matches$ani <= 100),
  all(matches$af_query >= 85 & matches$af_query <= 100),
  all(matches$af_target > 0 & matches$af_target <= 100)
)
# Percentages are already on a 0-100 scale: do not multiply by 100.
table_s2 <- matches[c("vOTU_id", "reference_accession", "qlen", "slen",
                      "ani", "af_query", "af_target")]
names(table_s2) <- c("vOTU", "Reference accession", "Query length (bp)",
                    "Reference length (bp)", "ANI (%)", "Query aligned (%)",
                    "Reference aligned (%)")
```

### Task 2: Preview Table S2

``` r
table_title <- "Table S2. Sequence matches linking the primary catalogue to published Äspö contigs."
knitr::kable(table_s2, digits = 2, caption = table_title)
```

| vOTU | Reference accession | Query length (bp) | Reference length (bp) | ANI (%) | Query aligned (%) | Reference aligned (%) |
|:---|:---|---:|---:|---:|---:|---:|
| vOTU0002 | MT141479.1 | 9763 | 27425 | 98.49 | 94.07 | 33.03 |
| vOTU0005 | MT141479.1 | 5132 | 27425 | 98.62 | 87.30 | 16.32 |
| vOTU0154 | MT144629.1 | 13937 | 26937 | 99.78 | 91.48 | 47.34 |
| vOTU0180 | MT143973.1 | 12866 | 56417 | 100.00 | 100.00 | 22.81 |
| vOTU0189 | MT143973.1 | 12366 | 56417 | 100.00 | 100.00 | 21.92 |
| vOTU0220 | MT144091.1 | 11425 | 14852 | 100.00 | 99.76 | 76.74 |
| vOTU0443 | MT143973.1 | 7901 | 56417 | 100.00 | 100.00 | 14.00 |
| vOTU0774 | MT143893.1 | 5597 | 43525 | 99.62 | 99.84 | 12.84 |
| vOTU0825 | MT144653.1 | 5471 | 21334 | 100.00 | 100.00 | 25.64 |

Table S2. Sequence matches linking the primary catalogue to published
Äspö contigs.

### Task 3: Export Excel

``` r
notes <- data.frame(
  Item = c("Title", "Source records", "Reference publication", "ANI",
           "Aligned fractions", "Selection", "Interpretation", "Display"),
  Description = c(
    table_title,
    "Archived public_reference_alignments_primary962.csv; existing sequence-match records, not new alignments.",
    "Holmfeldt et al. (2021). The Fennoscandian Shield deep terrestrial virosphere suggests slow motion 'boom and burst' cycles. Communications Biology 4:307. https://doi.org/10.1038/s42003-021-01810-1",
    "Average nucleotide identity (%). Query is a representative in the present 962-vOTU primary catalogue; reference is a deposited Äspö contig.",
    "Query aligned (%) and Reference aligned (%) use the query and reference lengths as separate denominators. Values are on a 0-100 scale.",
    "Existing alignment records meeting ANI >=95% and query-aligned fraction >=85% were retained. No new alignment or cross-study species dereplication was performed for this comparison.",
    "Multiple current representatives can match different parts of the same longer reference. Nine matched representatives and six references do not establish either nine or six shared viral species. Published reference coverage across external libraries is shown in Fig. S8.",
    "Numerical values retain the precision of the source records. Percentage columns display two decimal places; a zero displays as 0. Lengths are integer base pairs."
  ), stringsAsFactors = FALSE
)
wb <- openxlsx::createWorkbook(creator = "SwDGVirome")
openxlsx::modifyBaseFont(wb, fontSize = 11, fontName = "Arial", fontColour = "#000000")
openxlsx::addWorksheet(wb, "Table S2", gridLines = FALSE)
openxlsx::addWorksheet(wb, "Notes", gridLines = FALSE)
header_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#E8EEF3",
                                      wrapText = TRUE, valign = "center")
openxlsx::writeData(wb, "Table S2", table_s2, headerStyle = header_style, withFilter = TRUE)
openxlsx::setColWidths(wb, "Table S2", cols = 1:7, widths = c(15, 23, 21, 23, 15, 23, 25))
openxlsx::setRowHeights(wb, "Table S2", rows = 1, heights = 32)
openxlsx::freezePane(wb, "Table S2", firstRow = TRUE)
openxlsx::addStyle(wb, "Table S2", openxlsx::createStyle(numFmt = "#,##0"),
                   rows = 2:10, cols = 3:4, gridExpand = TRUE)
openxlsx::addStyle(wb, "Table S2", openxlsx::createStyle(numFmt = "0.00;-0.00;0"),
                   rows = 2:10, cols = 5:7, gridExpand = TRUE)
openxlsx::writeData(wb, "Notes", notes, headerStyle = header_style)
openxlsx::setColWidths(wb, "Notes", cols = 1:2, widths = c(25, 110))
openxlsx::addStyle(wb, "Notes", openxlsx::createStyle(wrapText = TRUE, valign = "top"),
                   rows = 2:9, cols = 1:2, gridExpand = TRUE)
openxlsx::setRowHeights(wb, "Notes", rows = 2:9, heights = c(30, 35, 65, 45, 45, 55, 75, 45))
# These sheets contain no drawings. Remove openxlsx's unused default links
# so other spreadsheet readers do not look for nonexistent drawing files.
for (i in seq_along(wb$worksheets)) {
  wb$worksheets[[i]]$drawing <- character(0)
  wb$worksheets_rels[[i]] <- character(0)
}
xlsx_file <- path_target("Table_S2_public_reference_matches.xlsx")
openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)
# Read the saved Excel file back to verify every value and header.
roundtrip <- openxlsx::read.xlsx(xlsx_file, sheet = "Table S2", check.names = FALSE)
stopifnot(isTRUE(all.equal(roundtrip, table_s2, check.attributes = FALSE)))
message("Verified Excel export: 9 representatives, 6 references; ", xlsx_file)
```

    Verified Excel export: 9 representatives, 6 references; data/11-table-s2-public-matches/Table_S2_public_reference_matches.xlsx

## Files written

These files have been written to data/11-table-s2-public-matches:

``` r
projthis::proj_dir_info(path_target(), tz = "CET") %>% knitr::kable()
```

| path                                   | type |  size | modification_time   |
|:---------------------------------------|:-----|------:|:--------------------|
| Table_S2_public_reference_matches.xlsx | file | 8.85K | 2026-09-15 15:12:23 |
