#' Compute standardized QC metrics from obs
#'
#' Extracts or computes QC metrics from cell metadata using auto-detected column names.
#' Returns a standardized data.frame.
#'
#' @param obs data.frame with cell metadata (from schard::h5ad2list)
#' @param qc_cols optional named vector to override QC column detection
#' @return data.frame with columns: cell_barcode, pct_mito, n_genes, n_UMI, doublet_score
#' @export
compute_qc_metrics <- function(obs, qc_cols = NULL) {
  if (!is.data.frame(obs) || nrow(obs) == 0) {
    stop("'obs' must be a non-empty data.frame", call. = FALSE)
  }

  detected <- detect_qc_cols(obs, qc_cols)

  result <- data.frame(
    cell_barcode = rownames(obs) %||% paste0("cell_", seq_len(nrow(obs))),
    pct_mito = if (!is.null(detected$pct_mito)) obs[[detected$pct_mito]] else NA_real_,
    n_genes = if (!is.null(detected$n_genes)) as.numeric(obs[[detected$n_genes]]) else NA_real_,
    n_UMI = if (!is.null(detected$n_UMI)) as.numeric(obs[[detected$n_UMI]]) else NA_real_,
    doublet_score = if (!is.null(detected$doublet_score)) as.numeric(obs[[detected$doublet_score]]) else NA_real_,
    stringsAsFactors = FALSE
  )

  missing <- names(detected)[vapply(detected, is.null, logical(1))]
  if (length(missing) > 0) {
    warning(
      "Could not detect QC columns for: ", paste(missing, collapse = ", "),
      ". Use 'qc_cols' to specify manually.", call. = FALSE
    )
  }

  cbind(result, obs)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Apply QC filter thresholds
#'
#' @param obs data.frame with cell metadata (must have detected QC columns)
#' @param thresholds list with threshold values, e.g.:
#'   list(pct_mito = 10, n_genes = c(500, 6000), n_UMI = c(1000, NA), doublet_score = 0.25)
#'   Each entry can be: a single number (upper bound), or a 2-element vector c(lower, upper).
#'   Use NA for one-sided bounds, e.g., c(500, NA) = lower bound only.
#' @param qc_cols optional named vector to override QC column detection
#' @return list with:
#'   - pass: logical vector (TRUE = passes all thresholds)
#'   - fail_reason: character vector with reason for failing ("" if pass)
#' @export
apply_filter_thresholds <- function(obs, thresholds = list(), qc_cols = NULL) {
  if (!is.data.frame(obs) || nrow(obs) == 0) {
    stop("'obs' must be a non-empty data.frame", call. = FALSE)
  }

  qc <- compute_qc_metrics(obs, qc_cols)
  n <- nrow(qc)
  pass <- rep(TRUE, n)
  fail_reason <- rep("", n)

  for (metric in names(thresholds)) {
    vals <- qc[[metric]]
    thresh <- thresholds[[metric]]

    if (is.null(vals) || all(is.na(vals))) next

    if (length(thresh) == 1) {
      failed <- !is.na(vals) & vals > thresh
      if (any(failed)) {
        pass[failed] <- FALSE
        fail_reason[failed] <- paste0(fail_reason[failed], sprintf("%s>%.1f; ", metric, thresh))
      }
    } else if (length(thresh) == 2) {
      if (!is.na(thresh[1])) {
        failed <- !is.na(vals) & vals < thresh[1]
        if (any(failed)) {
          pass[failed] <- FALSE
          fail_reason[failed] <- paste0(fail_reason[failed], sprintf("%s<%.1f; ", metric, thresh[1]))
        }
      }
      if (!is.na(thresh[2])) {
        failed <- !is.na(vals) & vals > thresh[2]
        if (any(failed)) {
          pass[failed] <- FALSE
          fail_reason[failed] <- paste0(fail_reason[failed], sprintf("%s>%.1f; ", metric, thresh[2]))
        }
      }
    }
  }

  list(pass = pass, fail_reason = fail_reason)
}

#' Compute per-replicate impact of filtering
#'
#' @param obs data.frame with cell metadata (including replicate column)
#' @param replicate_col character, column name identifying biological replicates
#' @param pass logical vector from apply_filter_thresholds
#' @param loss_threshold numeric, % cell loss that triggers a warning (default 25)
#' @param qc_cols optional named vector for QC column detection
#' @return data.frame with columns: replicate, total, passing, failing, pct_pass, pct_lost, warning
#' @export
impact_by_replicate <- function(obs, replicate_col, pass, loss_threshold = 25) {
  if (!replicate_col %in% colnames(obs)) {
    stop("Column '", replicate_col, "' not found in obs", call. = FALSE)
  }
  if (length(pass) != nrow(obs)) {
    stop("Length of 'pass' (", length(pass), ") must match nrow(obs) (", nrow(obs), ")", call. = FALSE)
  }

  replicates <- obs[[replicate_col]]
  total <- as.numeric(table(replicates))
  passing <- as.numeric(by(pass, replicates, sum))
  failing <- total - passing
  pct_pass <- passing / total
  pct_lost <- 1 - pct_pass
  warning <- pct_lost * 100 > loss_threshold

  result <- data.frame(
    replicate = names(table(replicates)),
    total = total,
    passing = passing,
    failing = failing,
    pct_pass = round(pct_pass * 100, 1),
    pct_lost = round(pct_lost * 100, 1),
    warning = warning,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  attr(result, "loss_threshold") <- loss_threshold
  result
}
