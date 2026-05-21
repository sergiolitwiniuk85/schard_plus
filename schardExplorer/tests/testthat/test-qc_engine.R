make_obs <- function() {
  data.frame(
    row.names = paste0("cell_", 1:20),
    n_genes_by_obs = c(rep(3000, 10), rep(500, 10)),
    total_counts = c(rep(10000, 15), rep(500, 5)),
    pct_mito = c(rep(5, 15), rep(25, 5)),
    doublet_score = runif(20, 0, 0.5),
    cluster = sample(LETTERS[1:4], 20, replace = TRUE),
    donor_id = sample(paste0("donor_", 1:2), 20, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

test_that("compute_qc_metrics returns standardized columns", {
  obs <- make_obs()
  result <- compute_qc_metrics(obs)
  expect_true(is.data.frame(result))
  expect_true("cell_barcode" %in% colnames(result))
  expect_true("pct_mito" %in% colnames(result))
  expect_true("n_genes" %in% colnames(result))
  expect_true("n_UMI" %in% colnames(result))
  expect_true("doublet_score" %in% colnames(result))
})

test_that("compute_qc_metrics errors on empty obs", {
  expect_error(compute_qc_metrics(data.frame()))
  expect_error(compute_qc_metrics("not_a_df"))
})

test_that("apply_filter_thresholds correctly identifies pass/fail", {
  obs <- make_obs()
  result <- apply_filter_thresholds(obs, list(pct_mito = 10, n_genes = c(1000, NA)))
  expect_type(result, "list")
  expect_true("pass" %in% names(result))
  expect_true("fail_reason" %in% names(result))
  expect_length(result$pass, 20)
  expect_length(result$fail_reason, 20)
})

test_that("apply_filter_thresholds handles no thresholds", {
  obs <- make_obs()
  result <- apply_filter_thresholds(obs, list())
  expect_true(all(result$pass))
})

test_that("apply_filter_thresholds errors on invalid input", {
  expect_error(apply_filter_thresholds("bad"))
  expect_error(apply_filter_thresholds(data.frame()))
})

test_that("impact_by_replicate computes per-replicate stats", {
  obs <- make_obs()
  pass <- rep(c(TRUE, FALSE), each = 10)
  result <- impact_by_replicate(obs, "donor_id", pass)
  expect_true(is.data.frame(result))
  expect_true("replicate" %in% colnames(result))
  expect_true("total" %in% colnames(result))
  expect_true("passing" %in% colnames(result))
  expect_true("pct_lost" %in% colnames(result))
  expect_true("warning" %in% colnames(result))
})

test_that("impact_by_replicate errors on mismatched lengths", {
  obs <- make_obs()
  expect_error(impact_by_replicate(obs, "donor_id", c(TRUE, FALSE)))
})

test_that("impact_by_replicate errors on missing column", {
  obs <- make_obs()
  pass <- rep(TRUE, 20)
  expect_error(impact_by_replicate(obs, "nonexistent", pass))
})
