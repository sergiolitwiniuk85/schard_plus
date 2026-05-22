library(testthat)
library(schard)

load_fixture <- function(name) {
  system.file("extdata", paste0(name, ".h5ad"), package = "schard")
}

# Test helpers from schardExplorer
make_test_obs <- function() {
  data.frame(
    cell_barcode = paste0("cell_", 1:50),
    n_genes = rpois(50, 3000),
    n_counts = rpois(50, 10000),
    pct_mito = runif(50, 0, 20),
    doublet_score = runif(50, 0, 0.5),
    cluster = sample(LETTERS[1:5], 50, replace = TRUE),
    donor_id = sample(paste0("donor_", 1:3), 50, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

make_test_list <- function() {
  obs <- make_test_obs()
  list(
    obs = obs,
    var = data.frame(gene = paste0("gene_", 1:100)),
    X = Matrix::rsparsematrix(50, 100, 0.1),
    obsm = list(X_umap = matrix(runif(100), 50, 2))
  )
}
