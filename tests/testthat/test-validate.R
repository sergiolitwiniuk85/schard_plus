test_that("h5ad_validate passes on valid fixture", {
  f <- load_fixture("test_matrix")
  skip_if(is.null(f) || !file.exists(f))

  result <- h5ad_validate(f)

  expect_true(result$valid)
  expect_length(result$errors, 0)
  expect_equal(result$dims$X, c(10, 5))  # H5AD stores (n_obs, n_var)
  expect_equal(result$dims$obs, 10)
  expect_equal(result$dims$var, 5)
})

test_that("h5ad_validate fails on non-existent file", {
  result <- h5ad_validate("nonexistent.h5ad")

  expect_false(result$valid)
  expect_match(result$errors[1], "not found")
})

test_that("h5ad_validate fails on non-HDF5 file", {
  tmp <- tempfile(fileext = ".h5ad")
  writeLines("not an h5ad file", tmp)
  on.exit(unlink(tmp))

  expect_error(h5ad_validate(tmp), "Cannot read")
})

test_that("h5ad_validate reports dim mismatch", {
  # Create a minimal corrupt H5AD: fake HDF5 with matching obs/var
  # but X with wrong shape
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp))

  # Build a minimal HDF5 that looks like an H5AD
  # but with /X/data having mismatched dimensions vs /obs /var
  rhdf5::h5createFile(tmp)
  rhdf5::h5createGroup(tmp, "obs")
  rhdf5::h5createGroup(tmp, "var")
  rhdf5::h5createGroup(tmp, "X")

  # Write obs with 10 rows
  rhdf5::h5write(letters[1:10], tmp, "obs/_index")
  # Write var with 5 rows
  rhdf5::h5write(LETTERS[1:5], tmp, "var/_index")

  # Write X with WRONG shape: 8x5 instead of 10x5
  # Use a minimal sparse matrix
  rhdf5::h5write(c(1, 2, 3), tmp, "X/data")
  rhdf5::h5write(c(0, 2, 3, 3), tmp, "X/indices")
  rhdf5::h5write(c(0, 1, 2, 3), tmp, "X/indptr")
  # Write attributes via open HDF5 handles
  fid <- rhdf5::H5Fopen(tmp)
  gid <- rhdf5::H5Gopen(fid, "X")
  rhdf5::h5writeAttribute(c(8, 5), h5obj = gid, name = "shape")
  rhdf5::h5writeAttribute("csr_matrix", h5obj = gid, name = "encoding-type")
  rhdf5::H5Gclose(gid)
  rhdf5::H5Fclose(fid)

  result <- h5ad_validate(tmp)

  expect_false(result$valid)
  expect_match(result$errors[1], "X has 8 rows.*obs has 10")
})
