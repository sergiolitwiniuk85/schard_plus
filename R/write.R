#' @import Matrix
#' @importFrom rhdf5 h5createFile h5createGroup h5write h5writeAttribute
#' @importFrom rhdf5 H5Fopen H5Fclose H5Gopen H5Gclose H5Dopen H5Dclose
#' @importFrom rhdf5 H5Lexists

.write_factor <- function(fid, path, x) {
  rhdf5::h5createGroup(fid, path)
  rhdf5::h5write(levels(x), fid, paste0(path, "/categories"))
  rhdf5::h5write(as.integer(x) - 1L, fid, paste0(path, "/codes"))
}

.write_csr <- function(fid, group_name, mat) {
  if (!inherits(mat, "dgCMatrix")) {
    mat <- as(as(mat, "Matrix"), "dgCMatrix")
  }
  csr <- as(Matrix::t(mat), "RsparseMatrix")

  rhdf5::h5createGroup(fid, group_name)
  rhdf5::h5write(as.numeric(csr@x), fid, paste0(group_name, "/data"))
  rhdf5::h5write(as.integer(csr@j), fid, paste0(group_name, "/indices"))
  rhdf5::h5write(as.integer(csr@p), fid, paste0(group_name, "/indptr"))

  gid <- rhdf5::H5Gopen(fid, group_name)
  rhdf5::h5writeAttribute("csr_matrix", gid, "encoding-type")
  rhdf5::h5writeAttribute("0.1.0", gid, "encoding-version")
  rhdf5::  h5writeAttribute(dim(csr), gid, "shape")
  rhdf5::H5Gclose(gid)
}

.write_dataframe <- function(fid, group_name, df) {
  rhdf5::h5createGroup(fid, group_name)

  index_col <- "_index"
  rnames <- rownames(df)
  if (is.null(rnames)) {
    rnames <- as.character(seq_len(nrow(df)))
  }
  rhdf5::h5write(rnames, fid, paste0(group_name, "/", index_col))

  col_order <- character(0)
  for (col in colnames(df)) {
    val <- df[[col]]
    path <- paste0(group_name, "/", col)
    col_order <- c(col_order, col)

    if (is.factor(val)) {
      .write_factor(fid, path, val)
    } else if (is.logical(val)) {
      rhdf5::h5write(as.integer(val), fid, path)
    } else if (is.character(val)) {
      rhdf5::h5write(val, fid, path)
    } else if (is.integer(val)) {
      rhdf5::h5write(val, fid, path)
    } else if (is.numeric(val)) {
      rhdf5::h5write(val, fid, path)
    }
  }

  gid <- rhdf5::H5Gopen(fid, group_name)
  rhdf5::h5writeAttribute("dataframe", gid, "encoding-type")
  rhdf5::h5writeAttribute("0.2.0", gid, "encoding-version")
  rhdf5::h5writeAttribute(index_col, gid, "_index")
  if (length(col_order) > 0) {
    rhdf5::h5writeAttribute(col_order, gid, "column-order")
  }
  rhdf5::H5Gclose(gid)
}

.write_obsm_list <- function(fid, obsm) {
  if (length(obsm) == 0) return(invisible())
  rhdf5::h5createGroup(fid, "obsm")
  for (nm in names(obsm)) {
    mat <- obsm[[nm]]
    if (is.null(mat)) next
    path <- paste0("obsm/", nm)
    rhdf5::h5write(mat, fid, path)
    did <- rhdf5::H5Dopen(fid, path)
    rhdf5::h5writeAttribute("array", did, "encoding-type")
    rhdf5::h5writeAttribute("0.2.0", did, "encoding-version")
    rhdf5::H5Dclose(did)
  }
}

.write_uns_list <- function(fid, uns, parent = "uns") {
  if (length(uns) == 0) return(invisible())
  if (!rhdf5::H5Lexists(fid, parent)) {
    rhdf5::h5createGroup(fid, parent)
  }
  for (nm in names(uns)) {
    val <- uns[[nm]]
    if (is.null(val)) next
    path <- paste0(parent, "/", nm)
    if (is.list(val) && !is.data.frame(val)) {
      .write_uns_list(fid, val, path)
    } else if (is.array(val) || is.matrix(val)) {
      rhdf5::h5write(val, fid, path)
    } else if (is.atomic(val) && length(val) > 0) {
      rhdf5::h5write(val, fid, path)
    }
  }
}

.write_h5ad_core <- function(file, X, obs, var, obsm, uns) {
  if (file.exists(file)) file.remove(file)

  rhdf5::h5createFile(file)
  rhdf5::h5closeAll()

  fid <- rhdf5::H5Fopen(file, "H5F_ACC_RDWR")
  tryCatch({
    .write_csr(fid, "X", X)
    .write_dataframe(fid, "obs", obs)
    .write_dataframe(fid, "var", var)
    .write_obsm_list(fid, obsm)
    .write_uns_list(fid, uns)
  }, finally = {
    rhdf5::H5Fclose(fid)
    rhdf5::h5closeAll()
  })
  invisible(file)
}

#' Write single-cell data to H5AD format
#'
#' Writes an R object to a scanpy-compatible H5AD (AnnData) file.
#' Methods exist for SingleCellExperiment, Seurat, and list objects.
#'
#' @param object An object to write. Supported classes:
#'   \code{SingleCellExperiment}, \code{Seurat}, or a \code{list} with
#'   components \code{X}, \code{obs}, \code{var}, \code{obsm}, and \code{uns}.
#' @param file Path to output \code{.h5ad} file.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns \code{file}.
#' @export
write_h5ad <- function(object, file, ...) {
  UseMethod("write_h5ad")
}

#' @rdname write_h5ad
write_h5ad.SingleCellExperiment <- function(object, file, ...) {
  loadRequiredPackages("SingleCellExperiment")

  mat <- Matrix::t(SummarizedExperiment::assay(object, 1))
  obs_df <- as.data.frame(SummarizedExperiment::colData(object))
  var_df <- as.data.frame(SummarizedExperiment::rowData(object))

  obsm <- as.list(SingleCellExperiment::reducedDims(object))
  obsm <- obsm[!vapply(obsm, is.null, logical(1))]

  uns <- S4Vectors::metadata(object)
  if (length(uns) == 0) uns <- list()

  .write_h5ad_core(file, mat, obs_df, var_df, obsm, uns)
}

#' @rdname write_h5ad
write_h5ad.Seurat <- function(object, file, ...) {
  loadRequiredPackages("Seurat")

  mat <- Matrix::t(Seurat::GetAssayData(object, slot = "counts"))

  obs_df <- object[[]]

  assay_obj <- Seurat::GetAssay(object)
  var_df <- as.data.frame(assay_obj[[]])
  if (nrow(var_df) > 0) {
    rownames(var_df) <- rownames(assay_obj)
  }

  obsm <- list()
  for (nm in Seurat::Reductions(object)) {
    emb <- tryCatch(
      Seurat::Embeddings(object, reduction = nm),
      error = function(e) NULL
    )
    if (!is.null(emb)) {
      obsm[[nm]] <- emb
    }
  }

  .write_h5ad_core(file, mat, obs_df, var_df, obsm, list())
}

#' @rdname write_h5ad
write_h5ad.list <- function(object, file, ...) {
  .write_h5ad_core(
    file,
    object[["X"]],
    if (is.null(object[["obs"]])) data.frame() else object[["obs"]],
    if (is.null(object[["var"]])) data.frame() else object[["var"]],
    if (is.null(object[["obsm"]])) list() else object[["obsm"]],
    if (is.null(object[["uns"]])) list() else object[["uns"]]
  )
}

#' @rdname write_h5ad
write_h5ad.default <- function(object, file, ...) {
  stop(
    "write_h5ad is not implemented for class ",
    paste(class(object), collapse = ", ")
  )
}
