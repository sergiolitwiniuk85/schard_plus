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
    return(h5ad2list(object, load.obsm = TRUE))
  }
  if (inherits(object, "SingleCellExperiment")) {
    loadRequiredPackages("SingleCellExperiment")
    obsm <- SingleCellExperiment::reducedDims(object)
    names(obsm) <- sub("^", "X_", tolower(names(obsm)))
    return(list(
      X = SummarizedExperiment::assay(object),
      obs = as.data.frame(SummarizedExperiment::colData(object)),
      var = as.data.frame(SummarizedExperiment::rowData(object)),
      obsm = obsm,
      uns = list()
    ))
  }
  if (inherits(object, "Seurat")) {
    loadRequiredPackages("Seurat")
    assay <- Seurat::DefaultAssay(object)
    obsm <- list()
    for (dr in Seurat::Reductions(object)) {
      obsm[[paste0("X_", tolower(dr))]] <- Seurat::Embeddings(object, reduction = dr)
    }
    return(list(
      X = Seurat::GetAssayData(object, assay = assay, slot = "counts"),
      obs = object[[]],
      var = as.data.frame(object[[assay]][[]]),
      obsm = obsm,
      uns = list()
    ))
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
