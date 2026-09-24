# Build a traceable ORF-level evidence table without rerunning plots or web queries.
# Computational assessments are distinct from measured outputs and author approval.

amg_key <- function(x) paste(x$vOTU_id, x$protein_id, sep = "|")
amg_present <- function(x) !is.na(x) & trimws(as.character(x)) != "" &
  !tolower(trimws(as.character(x))) %in% c("na", "n/a", "none", "-", "[]")
amg_collapse <- function(x) paste(sort(unique(as.character(x[amg_present(x)]))), collapse = "; ")
amg_text <- function(x) {
  if (length(x) == 0L || is.na(x)) return("")
  if (is.numeric(x)) return(format(x, digits = 15, scientific = NA, trim = TRUE))
  as.character(x)
}
amg_record <- function(x, fields) {
  fields <- intersect(fields, names(x))
  parts <- vapply(fields, function(nm) {
    value <- amg_text(x[[nm]])
    if (value == "") "" else paste0(nm, "=", value)
  }, character(1))
  paste(parts[nzchar(parts)], collapse = "; ")
}

build_amg_evidence <- function(analysis_dir, output_dir = NULL) {
  analysis_dir <- normalizePath(analysis_dir, mustWork = TRUE)
  if (is.null(output_dir)) output_dir <- file.path(analysis_dir, "data/04-fig4-heatmap/amg_curation")
  raw_dir <- file.path(analysis_dir, "data/00-raw/d04-amg-curation")
  primary_file <- file.path(analysis_dir, "data/01-tse-construction/tse_primary_962.rds")
  read_input <- function(name) {
    p <- file.path(raw_dir, name)
    if (!file.exists(p)) stop("Missing AMG evidence input: ", p)
    utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
  }
  manifest <- read_input("source_manifest.csv")
  observed_md5 <- unname(tools::md5sum(file.path(raw_dir, manifest$file)))
  if (anyNA(observed_md5) || !identical(observed_md5, manifest$md5)) {
    stop("An archived AMG input is missing or differs from source_manifest.csv.")
  }
  if (!requireNamespace("TreeSummarizedExperiment", quietly = TRUE))
    stop("Install TreeSummarizedExperiment before reading the primary TSE.")
  tse <- readRDS(primary_file)
  rd <- as.data.frame(SummarizedExperiment::rowData(tse))
  md <- S4Vectors::metadata(tse)
  stopifnot(nrow(tse) == 962L, ncol(tse) == 4L, !anyDuplicated(rownames(tse)))
  candidates <- lapply(c("amg_dramv", "amg_vibrant"), function(nm) {
    d <- md[[nm]]
    d$source_record <- seq_len(nrow(d))
    d <- d[d$vOTU_id %in% rownames(tse), , drop = FALSE]
    d$annotation_tool <- if (nm == "amg_dramv") "DRAM-v" else "VIBRANT"
    d
  })
  names(candidates) <- c("amg_dramv", "amg_vibrant")
  keys <- unique(do.call(rbind, lapply(candidates, function(d) d[c("vOTU_id", "protein_id")])) )
  keys <- keys[order(keys$vOTU_id, keys$protein_id), , drop = FALSE]
  rownames(keys) <- NULL
  stopifnot(nrow(candidates$amg_dramv) == 72L, nrow(candidates$amg_vibrant) == 124L,
            nrow(keys) == 167L, length(unique(keys$vOTU_id)) == 124L,
            !anyDuplicated(amg_key(keys)), !anyDuplicated(keys$protein_id))
  frozen <- read_input("frozen_candidates.csv")
  frozen <- frozen[as.logical(frozen$in_primary_962), , drop = FALSE]
  stopifnot(!anyDuplicated(amg_key(frozen)), setequal(amg_key(frozen), amg_key(keys)))
  frozen <- frozen[match(amg_key(keys), amg_key(frozen)), , drop = FALSE]
  carrier <- rd[match(keys$vOTU_id, rownames(rd)), , drop = FALSE]
  stopifnot(identical(as.character(carrier$representative_contig_id), frozen$amg_contig_id),
            all(carrier$representative_length_bp == frozen$amg_contig_length))
  all_dram <- md$anno_prot_dramv
  all_dram$source_record <- seq_len(nrow(all_dram))
  focal <- all_dram[match(amg_key(keys), amg_key(all_dram)), , drop = FALSE]
  stopifnot(!anyNA(focal$protein_id), identical(amg_key(focal), amg_key(keys)),
            all(focal$contig_id == carrier$representative_contig_id))
  candidate_dram <- all_dram[amg_key(all_dram) %in% amg_key(keys), , drop = FALSE]
  stopifnot(nrow(candidate_dram) == 167L, !anyDuplicated(amg_key(candidate_dram)))
  start <- pmin(as.numeric(focal$start_position), as.numeric(focal$end_position))
  end <- pmax(as.numeric(focal$start_position), as.numeric(focal$end_position))
  stopifnot(!anyNA(start), !anyNA(end), all(start >= 1), all(end <= carrier$representative_length_bp))
  aa_with_stop <- (end - start + 1) / 3
  stopifnot(all(abs(aa_with_stop - round(aa_with_stop)) < 1e-9))
  distance <- pmin(start - 1, carrier$representative_length_bp - end)
  source_name <- "data/01-tse-construction/tse_primary_962.rds"
  evidence <- list()
  add_evidence <- function(i, type, result, source, record, subject = "", scope = "Direct extraction from archived input") {
    if (!nzchar(result)) return(invisible(NULL))
    evidence[[length(evidence) + 1L]] <<- data.frame(
      evidence_id = sprintf("E%05d", length(evidence) + 1L),
      vOTU_id = keys$vOTU_id[i], protein_id = keys$protein_id[i],
      evidence_type = type, subject_id = subject, result = result,
      source_file = source, source_record = as.character(record),
      verification_scope = scope, stringsAsFactors = FALSE)
  }
  # Values preserve their source names and scales. Missing values are omitted,
  # never replaced by zero. These records are not additional candidate genes.
  add_table <- function(tab, protein_col, fields, type, source, scope, subject_col = NULL) {
    if (!protein_col %in% names(tab)) stop("Missing identity field in ", source)
    for (j in seq_len(nrow(tab))) {
      i <- match(as.character(tab[[protein_col]][j]), keys$protein_id)
      if (is.na(i)) next
      if ("vOTU_id" %in% names(tab) && amg_present(tab$vOTU_id[j]) && tab$vOTU_id[j] != keys$vOTU_id[i])
        stop("Conflicting vOTU identity in ", source, " row ", j)
      record <- if ("source_record" %in% names(tab)) tab$source_record[j] else j
      row_scope <- if (is.function(scope)) scope(tab[j, , drop = FALSE], i) else scope
      add_evidence(i, type, amg_record(tab[j, , drop = FALSE], fields), source, record,
                   if (is.null(subject_col)) keys$protein_id[i] else amg_text(tab[[subject_col]][j]), row_scope)
    }
  }
  left_n <- right_n <- left_hits <- right_hits <- integer(nrow(keys))
  for (i in seq_len(nrow(keys))) {
    add_evidence(i, "Carrier quality and viral predictions", amg_record(carrier[i, , drop = FALSE],
      c("representative_contig_id", "representative_length_bp", "checkv_quality", "checkv_completeness",
        "checkv_completeness_method", "checkv_contamination", "checkv_provirus", "genomad_virus_score",
        "genomad_n_hallmarks", "virsorter2_score", "vibrant_type")), source_name,
      paste0("rowData: ", keys$vOTU_id[i]), carrier$representative_contig_id[i])
    d <- all_dram[all_dram$vOTU_id == keys$vOTU_id[i] & all_dram$contig_id == focal$contig_id[i], , drop = FALSE]
    d <- d[order(as.numeric(d$start_position), as.numeric(d$end_position)), , drop = FALSE]
    pos <- which(d$protein_id == keys$protein_id[i])
    stopifnot(length(pos) == 1L)
    l <- if (pos > 1) seq.int(max(1, pos - 5), pos - 1) else integer()
    r <- if (pos < nrow(d)) seq.int(pos + 1, min(nrow(d), pos + 5)) else integer()
    viral_match <- amg_present(d$viral_id) | amg_present(d$vogdb_id)
    left_n[i] <- length(l); right_n[i] <- length(r)
    left_hits[i] <- sum(viral_match[l]); right_hits[i] <- sum(viral_match[r])
    for (side in c("left", "right")) {
      ix <- if (side == "left") l else r
      txt <- if (!length(ix)) "No annotated ORF available on this side of the contig." else
        paste(vapply(ix, function(j) amg_record(d[j, , drop = FALSE],
          c("protein_id", "start_position", "end_position", "ko_id", "kegg_hit", "viral_id", "vogdb_id", "is_transposon")), character(1)), collapse = " | ")
      add_evidence(i, paste("Nearest five ORFs:", side), txt, source_name,
        paste0("metadata$anno_prot_dramv rows ", paste(d$source_record[ix], collapse = ",")),
        focal$contig_id[i], "Coordinate-ordered neighbours; viral database matches are not proof of viral origin or AMG activity")
    }
  }
  for (nm in names(candidates)) {
    d <- candidates[[nm]]
    for (j in seq_len(nrow(d))) {
      i <- match(amg_key(d[j, , drop = FALSE]), amg_key(keys))
      add_evidence(i, paste(d$annotation_tool[j], "candidate annotation"),
        amg_record(d[j, , drop = FALSE], c("dbid", "db_desc", "auxiliary_score", "amg_flags", "Pfam", "Pfam.name", "metabolism")),
        source_name, paste0("metadata$", nm, " row ", d$source_record[j]), keys$protein_id[i])
    }
  }
  add_table(candidate_dram, "protein_id", c("start_position", "end_position", "strandedness", "rank", "ko_id", "kegg_hit",
    "auxiliary_score", "amg_flags", "is_transposon", "pfam_hits", "viral_id", "viral_hit", "viral_identity", "viral_bitScore",
    "viral_eVal", "vogdb_id", "vogdb_hits", "vogdb_categories", "source_record"), "DRAM-v protein annotation",
    paste0(source_name, "::metadata$anno_prot_dramv"), "Direct extraction; viral_identity is a fraction (0-1)")
  egg <- md$anno_prot_eggnog
  add_table(egg, "protein_id", c("seed_ortholog", "evalue", "score", "COG_category", "Description", "Preferred_name", "EC", "KEGG_ko", "PFAMs"),
    "eggNOG protein annotation", paste0(source_name, "::metadata$anno_prot_eggnog"), "Direct extraction; automatic annotation")
  phrog <- md$anno_prot_phrog
  phrog$source_record <- seq_len(nrow(phrog))
  phrog <- phrog[phrog$protein_id %in% keys$protein_id, , drop = FALSE]
  phrog_i <- match(phrog$protein_id, keys$protein_id)
  # DRAM-v spans include a possible stop codon. A one-residue difference is allowed;
  # larger discrepancies invalidate use of that PHROG record as focal-ORF evidence.
  phrog$length_compatible <- is.finite(as.numeric(phrog$qlen)) &
    abs(as.numeric(phrog$qlen) - aa_with_stop[phrog_i]) <= 1
  phrog_bad <- unique(phrog$protein_id[!phrog$length_compatible])
  for (j in seq_len(nrow(phrog))) {
    i <- phrog_i[j]
    add_evidence(i, "PHROG annotation (identity check required)", amg_record(phrog[j, , drop = FALSE],
      c("target", "annot", "category", "evalue", "bits", "pident", "qlen", "tlen", "qcov", "tcov", "length_compatible")),
      source_name, paste0("metadata$anno_prot_phrog row ", phrog$source_record[j]), phrog$target[j],
      if (phrog$protein_id[j] == "FCH3CGLCCXYL4FAA_NODE_14_9") "Query length compatible, but target-coverage/hit record anomalous; not treated as a verified functional conflict pending provenance checks" else if (phrog$length_compatible[j]) "Reported query length compatible with ORF coordinates; identity is not confirmed by length alone; pident is percent; coverage is a fraction" else
        "Reported query length incompatible with ORF coordinates: excluded from focal-protein functional assessment")
  }
  # Archived targeted measurements. Old recommendation/verdict columns are not used.
  specs <- list(
    list("targeted_mec_control_assignments.csv", "protein_id", c("tool","database","identifier","assignment","category_or_flag","statistics"), "Mec source-specific annotation", "Original database-specific scores preserved; assignments are predictions, not measured function"),
    list("targeted_domain_scan_hmmer.csv", "target_protein", c("target_len","pfam_name","pfam_accession","model_len","full_evalue","full_score","domain_ievalue","domain_score","hmm_from","hmm_to","ali_from","ali_to","model_coverage","target_coverage"), "Targeted HMMER domain scan", "Archived measurement; full-sequence and domain scores remain separate"),
    list("targeted_tmalign_results.csv", "query", c("reference","comparison","tm_score_norm_query","tm_score_norm_reference","aligned_length","rmsd_angstrom","seq_identity_over_alignment","tool","access_utc"), "Targeted structural alignment", "Query identity and original TM-align output cross-checked for the 1SUR and 1RX2 comparisons"),
    list("targeted_structure_assessment.csv", "target", c("kind","method","n_residues","mean_plddt","median_plddt","min_plddt","max_plddt","frac_plddt_gt70","frac_plddt_gt90","source","access_utc"), "Structure model quality", "Candidate-sequence model; predicted structure does not establish activity"),
    list("targeted_annotation_hits.csv", "protein_id", c("database","hit_id","hit_description","evalue","bitscore","percent_identity","query_coverage","target_coverage","query_len","target_len","is_best_hit"), "Archived targeted annotation hit", "Archived extraction; recheck identity/length against current source records; not automatically adopted as a functional conclusion"),
    list("targeted_tree_placement.csv", "query_tip", c("gene_set","rooting","ancestor_level","n_other_tips","support_shalrt_ufboot","clade_composition","nearest_characterised_reference","nearest_characterised_distance","tree_file"), "Targeted phylogenetic placement", "Archived targeted tree analysis; distinct from historical Figure S6 trees"),
    list("additional_candidate_orfs.csv", "protein_id", c("protein_length_aa","observed_aa_length","start_position","end_position","strand","starts_with_M","internal_stop","ambiguous_residues","sequence_available","length_consistent_with_coords","sequence_source_file"), "Targeted sequence identity", "Archived ORF identity and sequence checks"),
    list("additional_axisB_family_assignment.csv", "protein_id", c("protein_length_aa","best_pfam_family","best_pfam_accession","pfam_evalue","pfam_bitscore","pfam_n_envelopes","pfam_hmm_spans","pfam_envelope_spans","pfam_domain_ievalues","pfam_domain_bitscores","pfam_model_coverage","pfam_query_coverage","domain_completeness","other_pfam_cutga","permissive_only_pfam","distinct_kegg_orthologues"), "Additional HMMER domain evidence", "Raw HMMER domain outputs cross-checked; pfam_bitscore is the full-sequence score; old functional verdicts excluded"),
    list("additional_axisC_function_summary.csv", "protein_id", c("candidate_length_aa","best_reference_accession","reference_length_aa","sequence_identity_pct","query_coverage_pct","reference_coverage_pct","key_positions_checked","key_positions_conserved","key_residues_present","key_residues_absent","tm_score_norm_reference","tm_score_norm_query","rmsd_angstrom","aligned_residues","structure_status","experimental_data_for_this_orf"), "Additional functional comparison", "Archived summary; raw rfbA/acpP TM-align output unavailable in this package; low-identity residue mappings require caution"),
    list("additional_axisD_candidate_placement.csv", "protein_id", c("gene","tree_file","tree_label","support_on_subtending_branch","support_parent_clade","n_tips_parent_clade","terminal_branch_length","meaningful_support"), "Historical phylogenetic placement", "Raw tree labels preserved; model, bootstrap scheme and full historical run provenance not recovered")
  )
  for (s in specs) {
    tab <- read_input(s[[1]])
    scope <- s[[5]]
    fields <- s[[3]]
    if (s[[1]] == "targeted_tree_placement.csv") {
      tab$original_query_tip <- tab$query_tip
      tab$query_tip <- sub("^QUERY_", "", tab$query_tip)
      fields <- c("original_query_tip", fields)
      stopifnot(all(tab$query_tip %in% keys$protein_id))
    }
    if (s[[1]] == "targeted_annotation_hits.csv") {
      scope <- function(row, i) {
        if (!grepl("PHROG", row$database, ignore.case = TRUE))
          return("Archived annotation extraction; predictions are not functional measurements")
        if (!is.finite(as.numeric(row$query_len)) || abs(as.numeric(row$query_len) - aa_with_stop[i]) > 1)
          return("PHROG query length incompatible with ORF coordinates: excluded from focal-protein functional assessment")
        if (row$protein_id == "FCH3CGLCCXYL4FAA_NODE_14_9")
          return("PHROG query length compatible, but target-coverage/hit information anomalous; not treated as a verified functional conflict pending provenance checks")
        "PHROG query length compatible; length alone does not confirm identity; coverage remains in source units"
      }
    }
    if (s[[1]] == "targeted_tmalign_results.csv") {
      scope <- function(row, i) {
        checked <- (row$query == "TARATASF_NODE_16_52" && row$reference == "1SUR") ||
          (row$query == "FCH3CGLCCXYL4FAA_NODE_14_9" && row$reference == "1RX2")
        if (checked) "Query identity and original TM-align output cross-checked for this comparison" else
          "Archived TM-align summary; original comparison log not independently cross-checked in this review"
      }
    }
    add_table(tab, s[[2]], fields, s[[4]],
      file.path("data/00-raw/d04-amg-curation", s[[1]]), scope)
  }
  # Per-residue evidence retains alignment uncertainty instead of treating a
  # projected residue mismatch as proof that the candidate lacks activity.
  sites <- read_input("targeted_catalytic_residue_assessment.csv")
  sites$protein_id <- ifelse(sites$gene_set == "cysH", "TARATASF_NODE_16_52",
    ifelse(sites$gene_set == "folA", "FCH3CGLCCXYL4FAA_NODE_14_9", NA_character_))
  add_table(sites, "protein_id", setdiff(names(sites), "protein_id"),
    "Targeted reference-residue comparison", "data/00-raw/d04-amg-curation/targeted_catalytic_residue_assessment.csv",
    "Archived alignment projection; conservation does not establish enzymatic activity")
  add_table(read_input("additional_axisC_catalytic_residues.csv"), "protein_id",
    c("ref_accession","ref_id","ref_pos","ref_aa","feature_type","description","evidence","candidate_pos","candidate_aa","status","same_residue_within_3aa","key_functional_position","pair_identity_pct","pair_query_cov_pct","pair_ref_cov_pct"),
    "Additional reference-residue comparison", "data/00-raw/d04-amg-curation/additional_axisC_catalytic_residues.csv",
    "Candidate residue letters checked against archived FASTA; distant alignment columns and functional equivalence remain uncertain")
  focus_sites <- read_input("targeted_cysH_catalytic_focus.csv")
  i <- match("TARATASF_NODE_16_52", keys$protein_id)
  for (j in seq_len(nrow(focus_sites))) add_evidence(i, "CysH reference catalytic-site comparison",
    amg_record(focus_sites[j, , drop = FALSE], names(focus_sites)),
    "data/00-raw/d04-amg-curation/targeted_cysH_catalytic_focus.csv", j, focus_sites$reference[j],
    "Actual candidate alignment; residue identity alone does not measure catalysis")
  proxy <- read_input("targeted_mec_proxy_pfam_scan.csv")
  i <- match("TARATASF_NODE_16_18", keys$protein_id)
  for (j in seq_len(nrow(proxy))) add_evidence(i, "Mec reference-proxy domain scan",
    amg_record(proxy[j, , drop = FALSE], names(proxy)),
    "data/00-raw/d04-amg-curation/targeted_mec_proxy_pfam_scan.csv", j, proxy$query[j],
    "REFERENCE PROXY ONLY: not a recovered candidate sequence or direct candidate domain scan")

  summary <- data.frame(vOTU_id = keys$vOTU_id, protein_id = keys$protein_id,
    annotation_label = frozen$gene_label, annotation_tools = "", original_KO_EC = "",
    reporting_status = "Unresolved candidate AMG",
    function_assessment = "Automated annotation recorded; specific biochemical function not established by an individual targeted assessment.",
    auxiliary_role = "Not established", review_scope = "Annotation and genomic-context screening",
    reason = "Available annotation and context records were checked; no claim of experimentally validated activity or auxiliary role is made.",
    contig_id = as.character(carrier$representative_contig_id),
    contig_length_bp = as.numeric(carrier$representative_length_bp),
    checkv_quality = as.character(carrier$checkv_quality),
    checkv_completeness_pct = as.numeric(carrier$checkv_completeness),
    genomad_virus_score = as.numeric(carrier$genomad_virus_score),
    genomad_n_hallmarks = as.numeric(carrier$genomad_n_hallmarks),
    virsorter2_score = as.numeric(carrier$virsorter2_score),
    dramv_auxiliary_score = as.numeric(focal$auxiliary_score), dramv_flags = as.character(focal$amg_flags),
    bp_to_nearest_contig_end = distance, left_ORFs_with_viral_database_match = left_hits,
    left_ORFs_examined = left_n, right_ORFs_with_viral_database_match = right_hits,
    right_ORFs_examined = right_n, context_limitations = "", annotation_cautions = "",
    stringsAsFactors = FALSE)
  for (i in seq_len(nrow(keys))) {
    matches <- lapply(candidates, function(d) d[amg_key(d) == amg_key(keys[i, , drop = FALSE]), , drop = FALSE])
    summary$annotation_tools[i] <- amg_collapse(unlist(lapply(matches, `[[`, "annotation_tool")))
    summary$original_KO_EC[i] <- amg_collapse(unlist(lapply(matches, `[[`, "dbid")))
    context <- character()
    if (distance[i] < 5000) context <- c(context, "Within 5 kb of contig boundary (context limitation, not an exclusion rule)")
    if (left_n[i] < 5 || right_n[i] < 5) context <- c(context, "Fewer than five neighbouring ORFs available on at least one side")
    if (left_hits[i] == 0 || right_hits[i] == 0) context <- c(context, "No recorded viral-database match on at least one side")
    summary$context_limitations[i] <- if (length(context)) paste(context, collapse = "; ") else "No boundary or match-coverage flag under these checks"
    caution <- character()
    if (isTRUE(as.logical(frozen$annotation_conflict[i]))) caution <- c(caution, "Competing-annotation flag in the archived automated screen; not itself a resolved conflict")
    if (keys$protein_id[i] %in% phrog_bad) caution <- c(caution, "PHROG query-length mismatch: incompatible records excluded from functional assessment")
    if (keys$protein_id[i] == "FCH3CGLCCXYL4FAA_NODE_14_9") caution <- c(caution,
      "Archived PHROG target-coverage/hit information anomalous; not treated as a verified functional conflict")
    ko <- unlist(strsplit(summary$original_KO_EC[i], "; ", fixed = TRUE))
    if (length(ko) > 1) caution <- c(caution, "Multiple candidate annotation identifiers retained")
    summary$annotation_cautions[i] <- if (length(caution)) paste(caution, collapse = "; ") else "No additional flag under the recorded checks"
  }
  decisions <- read_input("curation_decisions.csv")
  stopifnot(nrow(decisions) == 20L, !anyDuplicated(amg_key(decisions)),
            all(amg_key(decisions) %in% amg_key(keys)))
  ix <- match(amg_key(decisions), amg_key(summary))
  for (nm in c("function_assessment", "reporting_status", "reason")) {
    stopifnot(all(amg_present(decisions[[nm]])))
    summary[[nm]][ix] <- decisions[[nm]]
  }
  summary$review_scope[ix] <- "Targeted computational evidence review"
  add_table(decisions, "protein_id", c("function_assessment", "reporting_status", "reason", "evidence_sources", "assessment_date", "author_confirmation"),
    "Current evidence assessment (interpretation)", "data/00-raw/d04-amg-curation/curation_decisions.csv",
    "Interpretation of archived results; separate from raw measurements and author approval")
  evidence <- do.call(rbind, evidence)
  rownames(evidence) <- NULL
  summary$n_evidence_records <- vapply(keys$protein_id, function(id) sum(evidence$protein_id == id), integer(1))
  stopifnot(nrow(summary) == 167L, !anyDuplicated(amg_key(summary)),
            !anyDuplicated(evidence$evidence_id), all(evidence$protein_id %in% keys$protein_id),
            all(summary$auxiliary_role == "Not established"), all(summary$n_evidence_records > 0L),
            max(nchar(evidence$result)) < 32000L)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write <- function(x, name) utils::write.csv(x, file.path(output_dir, name), row.names = FALSE, na = "")
  write(summary, "Table_S3_AMG_curation.csv")
  write(evidence, "Table_S3_AMG_evidence.csv")
  write(phrog, "AMG_PHROG_identity_check.csv")
  input_paths <- c(primary_file, file.path(raw_dir, c(manifest$file, "source_manifest.csv", "curation_decisions.csv")))
  write(data.frame(input = sub(paste0(analysis_dir, "/"), "", input_paths, fixed = TRUE),
                   md5 = unname(tools::md5sum(input_paths))), "AMG_input_provenance.csv")
  message("AMG table: ", nrow(summary), " ORFs / ", length(unique(summary$vOTU_id)), " vOTUs; ",
          nrow(evidence), " evidence records; ", length(phrog_bad), " ORFs with incompatible PHROG query lengths.")
  list(curation = summary, evidence = evidence, phrog_identity = phrog, output_dir = output_dir)
}

add_amg_sheets <- function(wb, amg) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Install openxlsx.")
  if (any(c("AMG curation", "AMG evidence") %in% names(wb)))
    stop("AMG sheets already exist; start from the preserved pre-AMG workbook.")
  header <- openxlsx::createStyle(fgFill = "#E8EEF3", textDecoration = "bold", wrapText = TRUE,
    valign = "center", border = "Bottom", borderColour = "#AEBBC6")
  body <- openxlsx::createStyle(wrapText = TRUE, valign = "top", border = "Bottom", borderColour = "#E5E7EB")
  note <- openxlsx::createStyle(wrapText = TRUE, valign = "top", fontColour = "#404040", fontSize = 10)
  for (sheet in c("AMG curation", "AMG evidence")) {
    d <- if (sheet == "AMG curation") amg$curation else amg$evidence
    openxlsx::addWorksheet(wb, sheet, gridLines = FALSE)
    notes <- if (sheet == "AMG curation") c(
      "Table S3. Candidate-AMG evidence assessment within the primary viral catalogue.",
      "One row per candidate ORF: 167 ORFs on 124 of the 962 primary-catalogue vOTUs. Original annotation labels are retained; they do not establish enzymatic activity.",
      "All candidates have annotation/context screening; 20 have targeted computational evidence review. Unresolved means the available evidence does not justify a stronger claim. No auxiliary metabolic role is experimentally established.",
      "Viral-database matches in up to five neighbouring ORFs are distinct from geNomad hallmark counts. A location within 5 kb of a contig end is a context flag, not an automatic exclusion rule.",
      "Blank numeric cells mean unavailable, not zero. Completeness is in percent; prediction scores retain source scales. Filter AMG evidence by protein_id to inspect the counted records. Assessments are computational interpretations, not experimental validation."
    ) else c(
      "Table S3. Detailed evidence records for candidate AMG annotations.",
      "One row per source record or defined context summary; multiple rows can describe the same candidate. These are not additional genes. Filter protein_id to trace the corresponding curation entry.",
      "Result fields preserve source metric names: TM-score normalization is explicit; RMSD is in angstroms; pident is percent; viral_identity and query/reference coverage are fractions unless the source field explicitly ends in _pct or percent_identity.",
      "Source_record refers to a data-row number (excluding the CSV header) or a stated TSE metadata row/key. Missing fields are unavailable; they are not zero. Proxy/reference evidence is explicitly labelled and is not treated as a candidate-sequence measurement.",
      "Rows flagged as incompatible PHROG query lengths are retained for transparency and excluded from functional conclusions. Archived summaries without original alignment logs remain labelled as such. Current assessment records are interpretations, not new measurements."
    )
    for (r in seq_along(notes)) {
      openxlsx::writeData(wb, sheet, notes[r], startRow = r, colNames = FALSE)
      openxlsx::mergeCells(wb, sheet, cols = 1:min(6, ncol(d)), rows = r)
      openxlsx::addStyle(wb, sheet, note, rows = r, cols = 1:min(6,ncol(d)), gridExpand = TRUE)
    }
    openxlsx::setRowHeights(wb, sheet, rows = 1:5, heights = c(22,30,38,38,38))
    openxlsx::writeData(wb, sheet, d, startRow = 7, withFilter = TRUE, headerStyle = header, keepNA = FALSE)
    openxlsx::setRowHeights(wb, sheet, rows = 7, heights = 58)
    openxlsx::addStyle(wb, sheet, body, rows = 8:(nrow(d)+7), cols = seq_len(ncol(d)), gridExpand = TRUE)
    widths <- if (sheet == "AMG curation") c(15,39,29,18,26,34,67,23,35,80,38,20,22,23,22,22,22,22,22,25,25,20,25,20,72,72,22) else c(15,15,39,37,39,140,67,43,80)
    stopifnot(length(widths) == ncol(d))
    openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(d)), widths = widths)
    # Fit wrapped evidence text without imposing the tallest row on all records.
    line_counts <- sapply(seq_len(ncol(d)), function(j)
      pmax(1, ceiling(nchar(ifelse(is.na(d[[j]]), "", as.character(d[[j]]))) / (widths[j] * 0.8))))
    heights <- pmin(409, pmax(36, apply(line_counts, 1, max) * 14 + 10))
    openxlsx::setRowHeights(wb, sheet, rows = 8:(nrow(d)+7), heights = heights)
    openxlsx::freezePane(wb, sheet, firstActiveRow = 8, firstActiveCol = if (sheet == "AMG curation") 3 else 4)
    for (nm in names(d)[vapply(d, is.numeric, logical(1))]) {
      fmt <- if (grepl("score|completeness", nm)) "0.0000" else "0"
      openxlsx::addStyle(wb, sheet, openxlsx::createStyle(numFmt = fmt), rows = 8:(nrow(d)+7),
                         cols = match(nm, names(d)), stack = TRUE)
    }
  }
  # Keep the author's Notes page last, without changing its content or formatting.
  order <- names(wb)
  order <- c(setdiff(order, c("AMG curation", "AMG evidence", "Notes")), "AMG curation", "AMG evidence", intersect("Notes",order))
  openxlsx::worksheetOrder(wb) <- match(order, names(wb))
  invisible(wb)
}

append_amg_to_s3 <- function(base_xlsx, output_xlsx, amg,
                             table_dir = NULL, replace_existing = FALSE) {
  if (!file.exists(base_xlsx)) stop("Missing base Table S3: ", base_xlsx)
  if (normalizePath(base_xlsx) == normalizePath(output_xlsx, mustWork = FALSE))
    stop("Use a new output file to preserve the existing author-formatted Table S3.")
  wb <- openxlsx::loadWorkbook(base_xlsx)
  stopifnot(all(c("Quality counts", "TPM comparison", "Taxonomic coverage", "vOTU data") %in% names(wb)))
  # This export is for the existing data-only Table S3. openxlsx creates default
  # drawing/printer/VML links even when these parts do not exist. Remove only
  # those empty defaults; do not silently discard a future workbook's graphics.
  base_parts <- utils::unzip(base_xlsx, list = TRUE)$Name
  # Inspect the saved XML: loadWorkbook itself can synthesize printer links.
  worksheet_parts <- base_parts[grepl("^xl/worksheets/sheet[0-9]+\\.xml$", base_parts)]
  printer_link <- vapply(worksheet_parts, function(part) {
    xml <- paste(readLines(unz(base_xlsx, part), warn = FALSE), collapse = "")
    grepl('<pageSetup\\b[^>]*\\br:id\\s*=', xml, perl = TRUE)
  }, logical(1))
  if (any(grepl("^xl/(drawings|media|charts)/|^xl/comments", base_parts)) || any(printer_link))
    stop("Base Table S3 contains graphics, comments or printer parts; use a preservation-aware export.")
  blocks <- list()
  if (!is.null(table_dir)) {
    specs <- list(
      c("Quality counts", "Table_S3a_quality_counts.csv", "1"),
      c("TPM comparison", "Table_S3_comparison.csv", "1"),
      c("TPM comparison", "Table_S3b_quality_sample_fractions.csv", "9"),
      c("Taxonomic coverage", "Table_S3_taxonomic_coverage.csv", "1"),
      c("vOTU data", "Table_S3_vOTU_quality_taxonomy.csv", "1")
    )
    for (spec in specs) {
      p <- file.path(table_dir, spec[2])
      if (!file.exists(p)) stop("Missing regenerated Table S3 data: ", p)
      d <- utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE,
                           na.strings = "")
      start <- as.integer(spec[3])
      old_header <- openxlsx::read.xlsx(base_xlsx, sheet = spec[1],
        rows = start, cols = seq_along(d), colNames = FALSE)
      if (!identical(as.character(unlist(old_header, use.names = FALSE)), names(d)))
        stop("The author workbook's table layout changed in ", spec[1],
             "; align the column headers before refreshing.")
      # Overwrite values only. Preserve cell styles, row heights, filters,
      # author footnotes outside the data blocks, and the complete Notes sheet.
      openxlsx::deleteData(wb, spec[1],
        rows = seq.int(start, start + nrow(d)), cols = seq_along(d), gridExpand = TRUE)
      openxlsx::writeData(wb, spec[1], d, startRow = start,
                          colNames = TRUE, keepNA = FALSE)
      blocks[[length(blocks) + 1L]] <- list(sheet = spec[1], row = start, data = d)
    }
  }
  if (replace_existing) {
    for (sheet in intersect(c("AMG curation", "AMG evidence"), names(wb)))
      openxlsx::removeWorksheet(wb, sheet)
  }
  add_amg_sheets(wb, amg)
  prepare_s3_data_export(wb)
  dir.create(dirname(output_xlsx), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, output_xlsx, overwrite = TRUE)
  blocks <- c(blocks, list(list(sheet = "AMG curation", row = 7L, data = amg$curation),
                           list(sheet = "AMG evidence", row = 7L, data = amg$evidence)))
  for (block in blocks) {
    actual <- openxlsx::read.xlsx(output_xlsx, sheet = block$sheet,
      rows = seq.int(block$row, block$row + nrow(block$data)),
      cols = seq_along(block$data), check.names = FALSE, na.strings = character())
    expected <- block$data
    for (nm in names(expected)) {
      if (is.character(expected[[nm]])) {
        expected[[nm]][is.na(expected[[nm]])] <- ""
        actual[[nm]][is.na(actual[[nm]])] <- ""
      }
    }
    stopifnot(isTRUE(all.equal(actual, expected, check.attributes = FALSE,
                               tolerance = 1e-10)))
  }
  invisible(output_xlsx)
}

prepare_s3_data_export <- function(wb) {
  for (i in seq_along(wb$worksheets)) {
    ws <- wb$worksheets[[i]]
    if (length(ws$sheet_data$rows)) {
      last <- paste0(openxlsx::int2col(max(ws$sheet_data$cols)), max(ws$sheet_data$rows))
      ws$dimension <- paste0('<dimension ref="A1:', last, '"/>')
    }
    ws$drawing <- character(0)
    ws$legacyDrawing <- character(0)
    ws$pageSetup <- gsub(' r:id="[^"]*"', "", ws$pageSetup)
    rels <- wb$worksheets_rels[[i]]
    empty_default <- grepl('/relationships/(drawing|printerSettings|vmlDrawing)"', rels)
    wb$worksheets_rels[[i]] <- rels[!empty_default]
  }
  invisible(wb)
}
