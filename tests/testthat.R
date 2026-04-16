library(testthat)
library(OpenMx)
if (requireNamespace("metaSEM")) {
  library(metaSEM)
}
library(metaDyn)
test_check("metaDyn")
