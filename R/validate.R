#' Validate H5AD file integrity
#'
#' Checks that an H5AD file has consistent dimensions between the expression
#' matrix and its annotations (obs, var). Useful for diagnosing corrupt or
#' malformed files that fail during loading.
#'
#' Unlike \code{\link{h5ad2list}}, this function does NOT load the full
#' expression matrix (\code{X}). It only reads the shape attribute and the
#' dimensions of obs/var, making it fast even for large files.
#'
#' @param filename character, path to the H5AD file.
#' @param check.obsm logical, whether to check that all entries in obsm have
#'   the correct number of rows (matching obs). Default \code{FALSE} because
#'   obsm may contain arbitrary matrices.
#'
#' @return A list with:
#' \describe{
#'   \item{valid}{\code{TRUE} if all checks pass, \code{FALSE} otherwise.}
#'   \item{errors}{character vector of issues found (empty if valid).}
#'   \item{dims}{named list with the dimensions of \code{X}, \code{obs}, \code{var},
#'     and optionally \code{raw/X} and \code{raw/var}.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' h5ad_validate("adata.h5ad")
#' h5ad_validate("adata.h5ad", check.obsm = TRUE)
#' }
h5ad_validate <- function(filename, check.obsm = FALSE) {
  errors <- character(0)

  # --- file existence ---
  if (!file.exists(filename)) {
    return(list(
      valid = FALSE,
      errors = paste0("File not found: ", filename),
      dims = list()
    ))
  }

  # --- try to read structure ---
  h5struct <- tryCatch(
    rhdf5::h5ls(filename),
    error = function(e) {
      stop("Cannot read '", filename, "' as HDF5: ", conditionMessage(e), call. = FALSE)
    }
  )

  # --- check required groups ---
  required <- c("/obs", "/var", "/X")
  for (grp in required) {
    if (!any(h5struct$group == grp)) {
      errors <- c(errors, paste0("Missing required group: ", grp))
    }
  }
  if (length(errors) > 0) {
    return(list(valid = FALSE, errors = errors, dims = list()))
  }

  # --- helper: read dataset dim (NULL if not found or not a dataset) ---
  .get_dim <- function(group, name) {
    row <- h5struct$group == group & h5struct$name == name
    if (!any(row)) return(NULL)
    raw_dim <- h5struct$dim[row]
    if (nchar(raw_dim) == 0 || is.na(raw_dim)) return(NULL)
    as.numeric(strsplit(raw_dim, "\\s+")[[1]])
  }

  # --- read X shape from attributes (fast, no data loaded) ---
  x_shape <- tryCatch(
    as.integer(rhdf5::h5readAttributes(filename, "X")$shape),
    error = function(e) NULL
  )

  # Alternative: read dim from data subgroup (sparse encoding)
  if (is.null(x_shape)) {
    x_dim <- .get_dim("/X", "data")
    if (!is.null(x_dim) && length(x_dim) == 1) {
      # sparse: data is 1D, need to read from indptr instead
      # fallback: read attributes again with different approach
      x_attr <- tryCatch(rhdf5::h5readAttributes(filename, "X"), error = function(e) NULL)
      if (!is.null(x_attr)) {
        x_shape <- as.integer(x_attr$shape %||% x_attr$h5sparse_shape)
      }
    } else if (!is.null(x_dim) && length(x_dim) == 2) {
      x_shape <- x_dim
    }
  }

  if (is.null(x_shape)) {
    # last resort: read indptr length to guess nrows+1
    indptr_dim <- .get_dim("/X", "indptr")
    if (!is.null(indptr_dim) && length(indptr_dim) == 1) {
      # indptr length = n_rows + 1 (CSR) or n_cols + 1 (CSC)
      # We can't know orientation, so skip shape inference
    }
  }

  if (is.null(x_shape)) {
    errors <- c(errors, "Cannot determine X matrix shape. The file may be corrupted or in an unsupported format.")
  }

  # --- read obs nrows ---
  obs_dim <- .get_dim("/obs", "_index")
  if (is.null(obs_dim)) {
    # try any column in obs
    obs_cols <- h5struct$name[h5struct$group == "/obs"]
    obs_cols <- setdiff(obs_cols, "__categories")
    if (length(obs_cols) > 0) {
      obs_dim <- .get_dim("/obs", obs_cols[1])
    }
  }
  n_obs <- if (!is.null(obs_dim) && length(obs_dim) == 1) obs_dim[1] else NULL

  # --- read var nrows ---
  var_dim <- .get_dim("/var", "_index")
  if (is.null(var_dim)) {
    var_cols <- h5struct$name[h5struct$group == "/var"]
    var_cols <- setdiff(var_cols, "__categories")
    if (length(var_cols) > 0) {
      var_dim <- .get_dim("/var", var_cols[1])
    }
  }
  n_var <- if (!is.null(var_dim) && length(var_dim) == 1) var_dim[1] else NULL

  # --- dimension validation ---
  dims <- list(
    X = x_shape,
    obs = n_obs,
    var = n_var
  )

  if (!is.null(x_shape) && !is.null(n_obs) && length(x_shape) == 2) {
    if (x_shape[1] != n_obs) {
      errors <- c(errors, sprintf(
        "X has %d rows (observations), but obs has %d rows.",
        x_shape[1], n_obs
      ))
    }
  }

  if (!is.null(x_shape) && !is.null(n_var) && length(x_shape) == 2) {
    if (x_shape[2] != n_var) {
      errors <- c(errors, sprintf(
        "X has %d columns (features), but var has %d rows.",
        x_shape[2], n_var
      ))
    }
  }

  # --- raw/X and raw/var (optional) ---
  has_raw <- any(h5struct$group == "/raw/X")
  if (has_raw) {
    raw_x_shape <- tryCatch(
      as.integer(rhdf5::h5readAttributes(filename, "raw/X")$shape),
      error = function(e) NULL
    )

    raw_var_dim <- .get_dim("/raw/var", "_index")
    if (is.null(raw_var_dim)) {
      raw_var_cols <- h5struct$name[h5struct$group == "/raw/var"]
      raw_var_cols <- setdiff(raw_var_cols, "__categories")
      if (length(raw_var_cols) > 0) {
        raw_var_dim <- .get_dim("/raw/var", raw_var_cols[1])
      }
    }

    n_raw_var <- if (!is.null(raw_var_dim) && length(raw_var_dim) == 1) raw_var_dim[1] else NULL
    dims$raw_X <- raw_x_shape
    dims$raw_var <- n_raw_var

    if (!is.null(raw_x_shape) && !is.null(n_obs) && length(raw_x_shape) == 2) {
      if (raw_x_shape[1] != n_obs) {
        errors <- c(errors, sprintf(
          "raw/X has %d rows (observations), but obs has %d rows.",
          raw_x_shape[1], n_obs
        ))
      }
    }
    if (!is.null(raw_x_shape) && !is.null(n_raw_var) && length(raw_x_shape) == 2) {
      if (raw_x_shape[2] != n_raw_var) {
        errors <- c(errors, sprintf(
          "raw/X has %d columns (features), but raw/var has %d rows.",
          raw_x_shape[2], n_raw_var
        ))
      }
    }
  }

  # --- obsm validation (optional, requires loading shape) ---
  if (check.obsm && !is.null(n_obs)) {
    obsm_groups <- unique(h5struct$group[startsWith(h5struct$group, "/obsm/")])
    for (grp in obsm_groups) {
      obsm_shape <- tryCatch(
        as.integer(rhdf5::h5readAttributes(filename, grp)$shape),
        error = function(e) NULL
      )
      if (!is.null(obsm_shape) && length(obsm_shape) == 2) {
        entry_name <- sub("^/obsm/", "", grp)
        if (obsm_shape[1] != n_obs && obsm_shape[2] != n_obs) {
          errors <- c(errors, sprintf(
            "obsm/%s has shape (%d, %d); expected one dimension to match obs (%d).",
            entry_name, obsm_shape[1], obsm_shape[2], n_obs
          ))
        }
        dims[[paste0("obsm/", entry_name)]] <- obsm_shape
      }
    }
  }

  # also read X shape from raw group attributes as a fallback
  if (has_raw && is.null(dims$raw_X)) {
    raw_x_attr <- tryCatch(rhdf5::h5readAttributes(filename, "raw/X"), error = function(e) NULL)
    if (!is.null(raw_x_attr)) {
      dims$raw_X <- as.integer(raw_x_attr$shape %||% raw_x_attr$h5sparse_shape)
    }
  }

  list(
    valid = length(errors) == 0,
    errors = errors,
    dims = dims
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
