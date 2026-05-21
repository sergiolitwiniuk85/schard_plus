test_that("launch_qc_dashboard accepts character path", {
  fixture <- system.file("extdata", "test_matrix.h5ad", package = "schard")
  if (file.exists(fixture)) {
    result <- launch_qc_dashboard(fixture)
    expect_true(is.list(result))
    expect_true("obs" %in% names(result))
  } else {
    skip("schard test fixture not available")
  }
})

test_that("launch_qc_dashboard accepts list", {
  obs <- make_test_obs()
  data <- list(obs = obs, var = data.frame(gene = 1:5), X = NULL, obsm = NULL, uns = NULL)
  result <- launch_qc_dashboard(data)
  expect_true(is.list(result))
  expect_equal(nrow(result$obs), 50)
})

test_that("launch_qc_dashboard errors on invalid input", {
  expect_error(launch_qc_dashboard(42))
  expect_error(launch_qc_dashboard("nonexistent.h5ad"))
})

test_that("qc_report accepts list and returns filtered data", {
  obs <- make_test_obs()
  data <- list(obs = obs, var = data.frame(gene = 1:5))
  result <- qc_report(data, tempdir(), thresholds = list(n_genes = c(2000, 5000)))
  expect_true(is.list(result))
  expect_true("pass" %in% names(result))
  expect_true("fail_reason" %in% names(result))
})

test_that("qc_report creates output directory if missing", {
  obs <- make_test_obs()
  data <- list(obs = obs)
  tmp <- file.path(tempdir(), "qc_test_dir")
  result <- qc_report(data, tmp)
  expect_true(dir.exists(tmp))
  unlink(tmp, recursive = TRUE)
})
