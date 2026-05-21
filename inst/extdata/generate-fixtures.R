local({
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install("rhdf5", ask = FALSE, update = FALSE)
  }
})

library(rhdf5)

set.seed(42)

n_obs <- 10
n_var <- 5

X_data <- as.numeric(sample(1:100, 20, replace = TRUE))
X_indices <- as.integer(sample(0:(n_var - 1), 20, replace = TRUE))
csr_counts <- tabulate(sample(1:n_obs, 20, replace = TRUE), n_obs)
X_indptr <- as.integer(c(0, cumsum(csr_counts)))

obs_names <- sprintf("cell_%03d", seq_len(n_obs))
var_names <- sprintf("gene_%03d", seq_len(n_var))
batch <- as.character(sample(c("a", "b"), n_obs, replace = TRUE))

obsm <- matrix(rnorm(n_obs * 5), nrow = n_obs)

out <- "inst/extdata/test_matrix.h5ad"
if (file.exists(out)) file.remove(out)

h5createFile(out)
h5closeAll()

h5createGroup(out, "obs")
h5write(obs_names, out, "obs/_index")
h5write(batch, out, "obs/batch")

h5createGroup(out, "var")
h5write(var_names, out, "var/_index")

h5createGroup(out, "X")
h5write(X_data, out, "X/data")
h5write(X_indices, out, "X/indices")
h5write(X_indptr, out, "X/indptr")

h5createGroup(out, "obsm")
h5write(obsm, out, "obsm/X_pca")

fid <- H5Fopen(out)

gid <- H5Gopen(fid, "obs")
h5writeAttribute("dataframe", gid, "encoding-type")
h5writeAttribute("0.2.0", gid, "encoding-version")
h5writeAttribute("_index", gid, "_index")
H5Gclose(gid)

gid <- H5Gopen(fid, "var")
h5writeAttribute("dataframe", gid, "encoding-type")
h5writeAttribute("0.2.0", gid, "encoding-version")
h5writeAttribute("_index", gid, "_index")
H5Gclose(gid)

gid <- H5Gopen(fid, "X")
h5writeAttribute("csr_matrix", gid, "encoding-type")
h5writeAttribute("0.1.0", gid, "encoding-version")
h5writeAttribute(c(n_var, n_obs), gid, "shape")
H5Gclose(gid)

did <- H5Dopen(fid, "obsm/X_pca")
h5writeAttribute("array", did, "encoding-type")
h5writeAttribute("0.2.0", did, "encoding-version")
H5Dclose(did)

H5Fclose(fid)
h5closeAll()

cat(sprintf("Created %s (%d obs x %d var)\n", out, n_obs, n_var))
