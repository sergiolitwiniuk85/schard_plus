#' Validate and normalize an object to a list
#'
#' Accepts list (from h5ad2list), SingleCellExperiment, Seurat, or character file path.
#' Returns a named list with components: X, obs, var, obsm, uns.
#' @param object An object to validate
#' @return A list with standardized components
#' @keywords internal
normalize_to_list <- function(object) {
  if (is.character(object)) {
    if (!file.exists(object)) {
      stop("File not found: ", object, call. = FALSE)
    }
    return(schard::h5ad2list(object))
  }
  if (inherits(object, "SingleCellExperiment")) {
    return(schard::h5ad2list(object))
  }
  if (inherits(object, "Seurat")) {
    return(schard::h5ad2list(object))
  }
  if (is.list(object) && !is.null(object$obs)) {
    return(object)
  }
  stop(
    "Unsupported object type. Supported: list (from h5ad2list), ",
    "SingleCellExperiment, Seurat, or character file path.",
    call. = FALSE
  )
}

#' Format a proportion as percentage string
#' @param x numeric between 0 and 1
#' @return character string like "45.2%"
#' @keywords internal
format_pct <- function(x) {
  paste0(round(x * 100, 1), "%")
}
