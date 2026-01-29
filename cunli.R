get_viruses_by_group <- function(tpm_mat, sample_groups, group_name) {
  samples <- names(sample_groups)[sample_groups == group_name]
  group_tpm <- tpm_mat[, samples, drop = FALSE]
  rownames(group_tpm)[rowSums(group_tpm) > 0]
}
