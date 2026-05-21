test_that("save_h5ad app builds correctly", {
  testthat::skip_if_not_installed("shinytest2")
  testthat::skip_on_cran()

  obj <- make_test_list()
  app <- run_dashboard_app(obj)
  expect_s3_class(app, "shiny.appobj")
})

# Helper: fixture with correct dims (R convention: X is n_var × n_obs)
make_filter_fixture <- function() {
  n_obs <- 50L
  n_var <- 100L
  list(
    X = Matrix::rsparsematrix(n_var, n_obs, 0.1),
    obs = data.frame(
      cell_barcode = paste0("cell_", seq_len(n_obs)),
      pct_mito = runif(n_obs, 0, 20),
      stringsAsFactors = FALSE
    ),
    var = data.frame(gene = paste0("gene_", seq_len(n_var)), stringsAsFactors = FALSE),
    obsm = list(X_umap = matrix(runif(n_obs * 2), n_obs, 2))
  )
}

test_that("subset logic works for filtered cells", {
  lst <- make_filter_fixture()
  keep <- c(rep(TRUE, 30), rep(FALSE, 20))

  filtered_list <- list(
    X = lst$X[, keep, drop = FALSE],
    obs = lst$obs[keep, , drop = FALSE],
    var = lst$var,
    obsm = lapply(lst$obsm, function(m) {
      if (is.matrix(m) || inherits(m, "Matrix")) m[keep, , drop = FALSE] else m
    })
  )

  expect_equal(ncol(filtered_list$X), 30)
  expect_equal(nrow(filtered_list$obs), 30)
  expect_equal(nrow(filtered_list$var), nrow(lst$var))
  expect_equal(nrow(filtered_list$obsm$X_umap), 30)
})

test_that("save_h5ad round-trips correctly", {
  lst <- make_filter_fixture()
  keep <- c(rep(TRUE, 30), rep(FALSE, 20))

  filtered_list <- list(
    X = lst$X[, keep, drop = FALSE],
    obs = lst$obs[keep, , drop = FALSE],
    var = lst$var,
    obsm = lapply(lst$obsm, function(m) {
      if (is.matrix(m) || inherits(m, "Matrix")) m[keep, , drop = FALSE] else m
    })
  )

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp))

  expect_error(schard::write_h5ad(filtered_list, tmp), NA)

  # Round-trip: read it back
  back <- schard::h5ad2list(tmp, load.obsm = TRUE)
  expect_equal(ncol(back$X), 30)
  expect_equal(nrow(back$obs), 30)
  expect_equal(nrow(back$var), nrow(lst$var))
  expect_true("X_umap" %in% names(back$obsm))
  expect_equal(nrow(back$obsm$X_umap), 30)
})

test_that("zero cells subset has zero cols/rows", {
  lst <- make_filter_fixture()
  keep <- rep(FALSE, 50)

  filtered_list <- list(
    X = lst$X[, keep, drop = FALSE],
    obs = lst$obs[keep, , drop = FALSE],
    var = lst$var,
    obsm = lst$obsm
  )

  expect_equal(ncol(filtered_list$X), 0)
  expect_equal(nrow(filtered_list$obs), 0)
})
