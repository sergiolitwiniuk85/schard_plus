#' Auto-detect QC column names in obs
#'
#' Maps common QC column names from various pipelines to standardized names.
#' @param obs data.frame containing cell metadata
#' @param qc_cols optional named vector to override detection: c(pct_mito = "my_mito_col")
#' @return named list with detected column names for each QC metric, or NULL if not found
#' @export
detect_qc_cols <- function(obs, qc_cols = NULL) {
  if (!is.data.frame(obs)) {
    stop("'obs' must be a data.frame", call. = FALSE)
  }

  if (!is.null(qc_cols)) {
    if (is.null(names(qc_cols))) {
      stop("'qc_cols' must be a named vector, e.g., c(pct_mito = 'my_col')", call. = FALSE)
    }
    result <- list()
    for (std_name in names(qc_cols)) {
      result[[std_name]] <- qc_cols[[std_name]]
    }
    return(result)
  }

  dictionary <- list(
    pct_mito = c("pct_mito", "percent.mito", "percent_mito", "mito_ratio", "pct_mitochondrial"),
    n_genes = c("n_genes_by_obs", "n_genes", "n_feature_by_obs", "n_features"),
    n_UMI = c("n_counts", "total_counts", "n_umis", "n_umi", "total_count", "n_count"),
    doublet_score = c("doublet_score", "doublet_scores", "doublet_score_pred")
  )

  colnames_lower <- tolower(colnames(obs))

  result <- list()
  for (std_name in names(dictionary)) {
    idx <- match(dictionary[[std_name]], colnames_lower)
    found <- idx[!is.na(idx)]
    if (length(found) > 0) {
      result[[std_name]] <- colnames(obs)[found[1]]
    } else {
      result[[std_name]] <- NULL
    }
  }

  return(result)
}
