test_that("detect_qc_cols finds standard column names", {
  obs <- data.frame(
    pct_mito = 1:10,
    n_genes_by_obs = 1:10,
    total_counts = 1:10,
    doublet_score = 1:10
  )
  result <- detect_qc_cols(obs)
  expect_equal(result$pct_mito, "pct_mito")
  expect_equal(result$n_genes, "n_genes_by_obs")
  expect_equal(result$n_UMI, "total_counts")
  expect_equal(result$doublet_score, "doublet_score")
})

test_that("detect_qc_cols handles alternative names", {
  obs <- data.frame(
    percent.mito = 1:10,
    n_genes = 1:10,
    n_umis = 1:10,
    doublet_scores = 1:10
  )
  result <- detect_qc_cols(obs)
  expect_equal(result$pct_mito, "percent.mito")
  expect_equal(result$n_genes, "n_genes")
  expect_equal(result$n_UMI, "n_umis")
  expect_equal(result$doublet_score, "doublet_scores")
})

test_that("detect_qc_cols returns NULL for missing columns", {
  obs <- data.frame(cell_id = 1:10)
  result <- detect_qc_cols(obs)
  expect_null(result$pct_mito)
  expect_null(result$n_genes)
  expect_null(result$n_UMI)
  expect_null(result$doublet_score)
})

test_that("detect_qc_cols accepts user override", {
  obs <- data.frame(my_mito = 1:10, my_genes = 1:10)
  result <- detect_qc_cols(obs, qc_cols = c(pct_mito = "my_mito", n_genes = "my_genes"))
  expect_equal(result$pct_mito, "my_mito")
  expect_equal(result$n_genes, "my_genes")
})

test_that("detect_qc_cols errors on invalid input", {
  expect_error(detect_qc_cols("not_a_df"))
  expect_error(detect_qc_cols(data.frame(x = 1:3), qc_cols = c("unnamed")))
})
