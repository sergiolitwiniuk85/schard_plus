test_that("h5ad2list loads fixture", {
  f <- load_fixture("test_matrix")
  skip_if(is.null(f) || !file.exists(f))

  lst <- h5ad2list(f, load.obsm = TRUE)

  expect_type(lst, "list")
  expect_true("X" %in% names(lst))
  expect_true("obs" %in% names(lst))
  expect_true("var" %in% names(lst))
  expect_true("obsm" %in% names(lst))

  expect_equal(dim(lst$X), c(5, 10))
  expect_equal(nrow(lst$obs), 10)
  expect_equal(nrow(lst$var), 5)

  expect_true("_index" %in% colnames(lst$obs))
  expect_true("batch" %in% colnames(lst$obs))

  expect_s4_class(lst$X, "dgCMatrix")
})

test_that("h5ad2data.frame loads obs", {
  f <- load_fixture("test_matrix")
  skip_if(is.null(f) || !file.exists(f))

  obs <- h5ad2data.frame(f, "obs")

  expect_s3_class(obs, "data.frame")
  expect_equal(nrow(obs), 10)
  expect_true("batch" %in% colnames(obs))
})

test_that("h5ad2Matrix loads X", {
  f <- load_fixture("test_matrix")
  skip_if(is.null(f) || !file.exists(f))

  mtx <- h5ad2Matrix(f, "X")

  expect_s4_class(mtx, "dgCMatrix")
  expect_equal(dim(mtx), c(5, 10))
})

test_that("h5ad2images handles missing images gracefully", {
  f <- load_fixture("test_matrix")
  skip_if(is.null(f) || !file.exists(f))

  expect_error(imgs <- schard:::h5ad2images(f), NA)
  imgs <- schard:::h5ad2images(f)
  expect_type(imgs, "list")
  expect_length(imgs, 0)
})
