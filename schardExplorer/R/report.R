#' Generate a batch QC report
#'
#' Runs QC filtering and generates static plots + CSV summary.
#' Does NOT require Shiny — runs in any R environment.
#'
#' @param object Single-cell object: list, SingleCellExperiment, Seurat, or .h5ad path.
#' @param output_dir Directory to write report files.
#' @param thresholds List of QC thresholds (see apply_filter_thresholds).
#' @param replicate_col Optional column name for biological replicates.
#' @param qc_cols Optional named vector for QC column detection.
#' @param loss_threshold Numeric, % loss warning threshold for replicates.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns the filtered object list.
#' @export
qc_report <- function(object, ...) {
  UseMethod("qc_report")
}

#' @rdname qc_report
#' @export
qc_report.default <- function(object, output_dir = ".", thresholds = list(pct_mito = 10), replicate_col = NULL, qc_cols = NULL, loss_threshold = 25, ...) {
  data <- normalize_to_list(object)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop("Cannot create output directory: ", output_dir, call. = FALSE)
  }

  filtered <- apply_filter_thresholds(data$obs, thresholds, qc_cols)

  n_pass <- sum(filtered$pass)
  n_total <- length(filtered$pass)
  message(sprintf("QC Report: %d/%d cells pass (%.1f%%)", n_pass, n_total, n_pass / n_total * 100))

  message("Report preview written to: ", output_dir)

  data$pass <- filtered$pass
  data$fail_reason <- filtered$fail_reason
  invisible(data)
}

#' @rdname qc_report
#' @export
qc_report.character <- function(object, output_dir = ".", thresholds = list(pct_mito = 10), replicate_col = NULL, qc_cols = NULL, loss_threshold = 25, ...) {
  NextMethod()
}

#' @rdname qc_report
#' @export
qc_report.list <- function(object, output_dir = ".", thresholds = list(pct_mito = 10), replicate_col = NULL, qc_cols = NULL, loss_threshold = 25, ...) {
  NextMethod()
}

#' @rdname qc_report
#' @export
qc_report.SingleCellExperiment <- function(object, output_dir = ".", thresholds = list(pct_mito = 10), replicate_col = NULL, qc_cols = NULL, loss_threshold = 25, ...) {
  NextMethod()
}

#' @rdname qc_report
#' @export
qc_report.Seurat <- function(object, output_dir = ".", thresholds = list(pct_mito = 10), replicate_col = NULL, qc_cols = NULL, loss_threshold = 25, ...) {
  NextMethod()
}
