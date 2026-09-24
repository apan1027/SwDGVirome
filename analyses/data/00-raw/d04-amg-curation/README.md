# Candidate-AMG evidence inputs for Table S3

The helper `analyses/R/amg_curation.R` combines the current primary TSE with these
archived inputs. It exports evidence; it does not rerun HMMER, ESMFold, TM-align,
gene calling or phylogenetic reconstruction. No external service is required.

## Inventory and units

- Filter `metadata(tse)$amg_dramv` and `$amg_vibrant` explicitly to the 962 primary
  representatives; TSE metadata does not automatically follow row subsetting.
- The 72 DRAM-v and 124 VIBRANT annotation rows represent 167 unique
  `(vOTU_id, protein_id)` pairs on 124 vOTUs. They are not 196 distinct genes.
- Original KO/EC labels are preserved. An annotation label, a KEGG pathway or a
  structural match does not by itself establish an auxiliary metabolic role.
- The Fig. 4 pathway aggregation uses carrier-vOTU/annotation units; it is not
  an independent gene-abundance measurement. Its deduplication can collapse
  multiple ORFs on the same carrier with the same annotation. This ORF inventory
  deliberately preserves all 167 pairs.

## What is recomputed

Current CheckV and viral-classifier fields are taken from primary TSE rowData.
DRAM-v coordinates identify each candidate and up to five neighbouring ORFs on
either side. A neighbouring ORF with a recorded RefSeq-viral or VOG match counts
as a **viral-database match**, not a confirmed structural hallmark. Counts are
separate from the original contig-level geNomad hallmark count. Boundary distance
is `min(start - 1, contig_length - end)` in bp, using 1-based inclusive
coordinates. The 5-kb flag describes limited context and is not an exclusion rule.

PHROG query lengths are compared with the coordinate span / 3, allowing one
residue for a terminal stop codon. Incompatible records are retained as evidence
of a source inconsistency but excluded from functional assessment. Length
compatibility alone does not prove identity. The focal folA PHROG coverage/hit
record has an additional unresolved provenance warning. An annotation-record
problem does not establish that the corresponding gene is absent or inactive.

## Evidence and interpretation

All 167 candidates receive annotation/context screening. Twenty have archived
targeted evidence, and their current, conservative interpretations are specified
in `curation_decisions.csv`. The other 147 retain unresolved specific function
and auxiliary role; they have not received an equivalent targeted assessment.
No experimentally demonstrated auxiliary role is claimed for any candidate.
Missing evidence is blank or explicitly unavailable, never a zero result.

The two focal sulfur-labelled ORFs are not retained as evidence for the claimed
sulfur-metabolism function, but remain in the complete candidate inventory.
Their assessment is not extrapolated to all cysH/mec-labelled ORFs. The mec domain
scan uses RefSeq proxy YP_009909783.1, not a recovered candidate sequence.
RNR/DHFR candidates are not automatically excluded merely because their possible
function supports viral replication. DNA methyltransferase family support is
distinguished from a DNMT3A-specific assignment.

`source_manifest.csv` records original repository-relative paths, file sizes,
MD5 and SHA-256 for the 23 unchanged archived CSVs. The helper checks their MD5s.
`curation_decisions.csv` is a new interpretation table, not a copied measurement
or author approval record. The output `AMG_input_provenance.csv` records checksums
for the actual TSE, archived inputs, manifest and decisions used in each run.

`AMG evidence` contains one source record or defined context summary per row.
Source-row numbers exclude CSV headers; TSE rows/keys are identified explicitly.
Each summary row gives its evidence-record count; filter the evidence sheet by protein_id to inspect those records. Old automated verdicts
are not adopted as current curation decisions.

## Verification limits

The original logs and candidate-model sequence identities were independently
checked for cysH versus PDB 1SUR and folA versus PDB 1RX2. Other comparisons are
labelled according to their actual verification level. Original rfbA/acpP
TM-align logs were not recovered, so those numbers remain archived summaries.
Reference residue projections, particularly at low identity, do not measure
enzyme activity. Historical S6 tree settings/run provenance remain incomplete;
archived support values do not resolve that limitation. Tool/database versions
and reference provenance are retained in the corresponding source CSVs.

Current assessments were assembled on 2026-09-24. Author confirmation is not
invented; the interpretation records explicitly state that it is not recorded.
