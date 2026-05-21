# This test verifies the Shiny app starts correctly
# It does NOT test interactivity (that requires a browser)

test_that("run_dashboard_app returns a shiny.appobj", {
  testthat::skip_if_not_installed("shinytest2")
  testthat::skip_on_cran()

  obs <- make_test_obs()
  obj <- list(obs = obs, obsm = list(X_umap = matrix(runif(100), 50, 2)))

  app <- run_dashboard_app(obj)
  expect_s3_class(app, "shiny.appobj")
  expect_s3_class(app, "shiny.appobj")
})
