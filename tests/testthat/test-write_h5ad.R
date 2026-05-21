library(Matrix)

make_test_list <- function() {
  set.seed(1)
  x <- as(Matrix::rsparsematrix(5, 10, density = 0.4), "dgCMatrix")
  rownames(x) <- sprintf("gene_%03d", 1:5)
  colnames(x) <- sprintf("cell_%03d", 1:10)

  obs <- data.frame(
    batch = factor(sample(c("a", "b"), 10, replace = TRUE)),
    treatment = sample(c("ctrl", "drug"), 10, replace = TRUE),
    row.names = sprintf("cell_%03d", 1:10)
  )
  var <- data.frame(
    chr = sample(c("chr1", "chr2"), 5, replace = TRUE),
    row.names = sprintf("gene_%03d", 1:5)
  )
  obsm <- list(
    X_pca = matrix(rnorm(50), nrow = 10,
      dimnames = list(sprintf("cell_%03d", 1:10), paste0("PC", 1:5)))
  )

  list(X = x, obs = obs, var = var, obsm = obsm)
}

test_that("Basic round-trip with list", {
  lst <- make_test_list()
  tmp <- tempfile(fileext = ".h5ad")

  write_h5ad(lst, tmp)
  out <- h5ad2list(tmp, load.obsm = TRUE)

  expect_true(file.exists(tmp))
  expect_type(out, "list")
  expect_equal(dim(out$X), dim(lst$X))
  expect_equal(rownames(out$X), rownames(lst$X))
  expect_equal(colnames(out$X), colnames(lst$X))
  expect_equal(as.matrix(out$X), as.matrix(lst$X))
  expect_equal(nrow(out$obs), nrow(lst$obs))
  expect_equal(out$obs$batch, as.character(lst$obs$batch))
  expect_equal(out$obs$treatment, lst$obs$treatment)
  expect_equal(out$var$chr, lst$var$chr)
})

test_that("CSR encoding correctness", {
  set.seed(2)
  mat <- as(Matrix::rsparsematrix(5, 10, density = 0.4), "dgCMatrix")
  rownames(mat) <- sprintf("gene_%03d", 1:5)
  colnames(mat) <- sprintf("cell_%03d", 1:10)

  lst <- list(
    X = mat,
    obs = data.frame(row.names = colnames(mat)),
    var = data.frame(row.names = rownames(mat))
  )
  tmp <- tempfile(fileext = ".h5ad")
  write_h5ad(lst, tmp)

  x_attrs <- rhdf5::h5readAttributes(tmp, "/X")
  expect_equal(c(x_attrs[["encoding-type"]]), "csr_matrix")
  expect_equal(c(x_attrs[["encoding-version"]]), "0.1.0")
  expect_equal(c(x_attrs[["shape"]]), rev(dim(mat)))

  xdata <- rhdf5::h5read(tmp, "/X")
  expect_true("indptr" %in% names(xdata))
  expect_true("indices" %in% names(xdata))
  expect_true("data" %in% names(xdata))
  expect_length(xdata$indptr, ncol(mat) + 1)
})

test_that("Metadata encoding: factor levels preserved", {
  set.seed(3)
  mat <- as(Matrix::rsparsematrix(3, 5, density = 0.5), "dgCMatrix")
  obs <- data.frame(
    cell_type = factor(sample(c("T-cell", "B-cell", "NK"), 5, replace = TRUE),
      levels = c("T-cell", "B-cell", "NK")),
    condition = c("ctrl", "drug", "ctrl", "drug", "ctrl"),
    is_active = c(TRUE, FALSE, TRUE, TRUE, FALSE),
    value = c(1.1, 2.2, 3.3, 4.4, 5.5),
    row.names = paste0("cell_", 1:5)
  )
  var <- data.frame(row.names = paste0("gene_", 1:3))

  lst <- list(X = mat, obs = obs, var = var)
  tmp <- tempfile(fileext = ".h5ad")
  write_h5ad(lst, tmp)

  out <- h5ad2list(tmp)

  expect_equal(as.character(out$obs$cell_type), as.character(obs$cell_type))
  expect_equal(out$obs$condition, obs$condition)
  expect_equal(out$obs$value, obs$value)
  expect_equal(as.logical(out$obs$is_active), obs$is_active)
})

test_that("obsm round-trip", {
  set.seed(4)
  mat <- as(Matrix::rsparsematrix(3, 5, density = 0.5), "dgCMatrix")
  obs <- data.frame(row.names = paste0("cell_", 1:5))
  var <- data.frame(row.names = paste0("gene_", 1:3))

  pca <- matrix(rnorm(25), nrow = 5,
    dimnames = list(paste0("cell_", 1:5), paste0("PC", 1:5)))

  lst <- list(X = mat, obs = obs, var = var, obsm = list(X_pca = pca))
  tmp <- tempfile(fileext = ".h5ad")
  write_h5ad(lst, tmp)

  out <- h5ad2list(tmp, load.obsm = TRUE)

  expect_true("X_pca" %in% names(out$obsm))
  expect_equal(dim(out$obsm$X_pca), dim(pca))
  expect_equal(as.vector(out$obsm$X_pca), as.vector(pca))
})

test_that("write_h5ad errors on invalid object", {
  tmp <- tempfile(fileext = ".h5ad")
  expect_error(write_h5ad(iris, tmp), "write_h5ad is not implemented")
})

test_that("write_h5ad errors on invalid path", {
  set.seed(5)
  mat <- as(Matrix::rsparsematrix(3, 5, density = 0.5), "dgCMatrix")
  lst <- list(X = mat, obs = data.frame(row.names = paste0("cell_", 1:5)),
    var = data.frame(row.names = paste0("gene_", 1:3)))
  expect_error(write_h5ad(lst, "/nonexistent/dir/out.h5ad"))
})

test_that("obs with zero rows handled gracefully", {
  mat <- as(Matrix::rsparsematrix(0, 3, density = 0), "dgCMatrix")
  obs <- data.frame(row.names = character(0))
  var <- data.frame(row.names = paste0("gene_", 1:3))
  lst <- list(X = mat, obs = obs, var = var)
  tmp <- tempfile(fileext = ".h5ad")

  expect_error(write_h5ad(lst, tmp), NA)
})

test_that("obs/var with no additional columns handled gracefully", {
  mat <- as(Matrix::rsparsematrix(3, 5, density = 0.5), "dgCMatrix")
  obs <- data.frame(row.names = paste0("cell_", 1:5))
  var <- data.frame(row.names = paste0("gene_", 1:3))
  lst <- list(X = mat, obs = obs, var = var)
  tmp <- tempfile(fileext = ".h5ad")
  write_h5ad(lst, tmp)

  out <- h5ad2list(tmp)
  expect_equal(nrow(out$obs), 5)
  expect_equal(nrow(out$var), 3)
})

test_that("list with missing components uses defaults", {
  mat <- as(Matrix::rsparsematrix(3, 5, density = 0.5), "dgCMatrix")
  lst <- list(X = mat)
  tmp <- tempfile(fileext = ".h5ad")
  expect_error(write_h5ad(lst, tmp), NA)
})
