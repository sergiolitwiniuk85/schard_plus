test_that("qc_report writes all expected files", {
  obs <- make_test_obs()
  obj <- list(obs = obs, var = data.frame(gene = 1:5))
  obj$obsm <- list(X_umap = matrix(runif(100), 50, 2))

  tmp <- file.path(tempdir(), "qc_test_batch")
  result <- qc_report(obj, tmp, thresholds = list(pct_mito = 10))

  expect_true(file.exists(file.path(tmp, "filtered_cells.csv")))
  expect_true(file.exists(file.path(tmp, "qc_summary.html")))

  expect_true(is.list(result))
  expect_true("pass" %in% names(result))

  unlink(tmp, recursive = TRUE)
})

test_that("qc_report works with replicate tracking", {
  obs <- make_test_obs()
  obj <- list(obs = obs)
  tmp <- file.path(tempdir(), "qc_test_replicates")

  result <- qc_report(obj, tmp, thresholds = list(n_genes = c(2000, 5000)),
                       replicate_col = "donor_id")

  expect_true(file.exists(file.path(tmp, "impact_by_replicate.csv")))
  expect_true(file.exists(file.path(tmp, "impact_by_replicate.png")))

  unlink(tmp, recursive = TRUE)
})

test_that("qc_report errors on non-writable directory", {
  obs <- make_test_obs()
  obj <- list(obs = obs)

  expect_error(qc_report(obj, "/root_qc_test"))
})
