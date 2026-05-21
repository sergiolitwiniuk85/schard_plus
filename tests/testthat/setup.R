library(testthat)
library(schard)

load_fixture <- function(name) {
  system.file("extdata", paste0(name, ".h5ad"), package = "schard")
}
