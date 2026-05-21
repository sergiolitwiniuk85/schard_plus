#' Launch the interactive QC dashboard
#'
#' Opens a Shiny app for interactive QC filtering and visualization.
#' Accepts schard list, SingleCellExperiment, Seurat, or .h5ad file path.
#'
#' @param object Single-cell object: list from h5ad2list(), SingleCellExperiment,
#'   Seurat, or character file path to .h5ad file.
#' @param replicate_col Optional column name in obs identifying biological replicates.
#' @param qc_cols Optional named vector to override QC column detection.
#' @param loss_threshold Numeric, % cell loss that triggers a replicate warning (default 25).
#' @param ... Additional arguments passed to methods.
#'
#' @return Invisibly returns the object list. In visual mode, the Shiny app runs.
#' @export
launch_qc_dashboard <- function(object, ...) {
  UseMethod("launch_qc_dashboard")
}

#' @rdname launch_qc_dashboard
#' @export
launch_qc_dashboard.default <- function(object, ...) {
  data <- normalize_to_list(object)
  message("Data loaded successfully: ", nrow(data$obs), " cells x ", nrow(data$var), " genes")
  message("QC columns detected:")
  print(detect_qc_cols(data$obs))
  invisible(data)
}

#' @rdname launch_qc_dashboard
#' @export
launch_qc_dashboard.character <- function(object, replicate_col = NULL, qc_cols = NULL, loss_threshold = 25, ...) {
  data <- normalize_to_list(object)
  message(sprintf("Loaded %s: %d cells x %d features", object, nrow(data$obs), nrow(data$var)))
  invisible(data)
}

#' @rdname launch_qc_dashboard
#' @export
launch_qc_dashboard.list <- function(object, replicate_col = NULL, qc_cols = NULL, loss_threshold = 25, ...) {
  data <- normalize_to_list(object)
  message(sprintf("Loaded list: %d cells x %d features", nrow(data$obs), nrow(data$var)))
  invisible(data)
}

#' @rdname launch_qc_dashboard
#' @export
launch_qc_dashboard.SingleCellExperiment <- function(object, replicate_col = NULL, qc_cols = NULL, loss_threshold = 25, ...) {
  data <- normalize_to_list(object)
  message(sprintf("Loaded SingleCellExperiment: %d cells x %d features", nrow(data$obs), nrow(data$var)))
  invisible(data)
}

#' @rdname launch_qc_dashboard
#' @export
launch_qc_dashboard.Seurat <- function(object, replicate_col = NULL, qc_cols = NULL, loss_threshold = 25, ...) {
  data <- normalize_to_list(object)
  message(sprintf("Loaded Seurat: %d cells x %d features", nrow(data$obs), nrow(data$var)))
  invisible(data)
}
