#' Auto-detect QC column names in obs
#'
#' Maps common QC column names from various pipelines to standardized names.
#' @param obs data.frame containing cell metadata
#' @param qc_cols optional named vector to override detection: c(pct_mito = "my_mito_col").
#'   Pass an empty list (or omit) to auto-detect.
#' @return named list with detected column names for each QC metric, or NULL if not found
#' @export
detect_qc_cols <- function(obs, qc_cols = NULL) {
  if (!is.data.frame(obs)) {
    stop("'obs' must be a data.frame", call. = FALSE)
  }

  # Manual override: only when non-empty and named
  if (!is.null(qc_cols) && length(qc_cols) > 0) {
    if (is.null(names(qc_cols))) {
      stop("'qc_cols' must be a named vector, e.g., c(pct_mito = 'my_col')", call. = FALSE)
    }
    result <- list()
    for (std_name in names(qc_cols)) {
      result[[std_name]] <- qc_cols[[std_name]]
    }
    return(result)
  }

  # Empty list or NULL → auto-detect
  # All entries in lowercase: matching is case-insensitive via tolower(colnames)
  dictionary <- list(
    pct_mito = c(
      "pct_mito", "percent.mito", "percent_mito", "mito_ratio",
      "pct_mitochondrial", "subsets_mito_percent", "percent.mt",
      "mt_ratio", "mito_percent", "subsets_mitochondrial_percent"
    ),
    n_genes = c(
      "n_genes_by_obs", "n_genes", "n_feature_by_obs", "n_features",
      "nfeature_rna", "n_genes_by_counts", "gene_count",
      "n_feature_rna", "n_gene"
    ),
    n_UMI = c(
      "n_counts", "total_counts", "n_umis", "n_umi", "total_count",
      "n_count", "ncount_rna", "n_umis_by_counts", "total_umi",
      "n_umis", "n_count_rna", "n_umi_by_counts"
    ),
    doublet_score = c(
      "doublet_score", "doublet_scores", "doublet_score_pred",
      "scdblfinder.score", "scrublet_score", "doubletfinder_score",
      "doublet_finder_score"
    )
  )

  colnames_lower <- tolower(colnames(obs))

  result <- list()
  for (std_name in names(dictionary)) {
    idx <- match(dictionary[[std_name]], colnames_lower)
    found <- idx[!is.na(idx)]
    if (length(found) > 0) {
      result[[std_name]] <- colnames(obs)[found[1]]
    }
    # else: key is omitted from result (entry not present)
  }

  return(result)
}
