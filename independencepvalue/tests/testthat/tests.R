# Generated from create-independencepvalue.Rmd: do not edit by hand  
testthat::test_that("classical_p_val works", {
  set.seed(1)
  X <- matrix(rnorm(50), 10, 5)
  testthat::expect_equal(classical_p_val(S=cov(X), n=10, CP=c(rep(1, each = 3),rep(2, each = 2)), k=1, mc_iter=100), 0.83)
})

testthat::test_that("block_diag works", {
  set.seed(1)
  X <- matrix(rnorm(50), 10, 5)
  testthat::expect_equal(length(unique(block_diag(cor(X), c=0.5))), 3)
})

testthat::test_that("selective_p_val works", {
  set.seed(1)
  X <- matrix(rnorm(50), 10, 5)
  corX <- cor(cbind(X))
  block_diag_structure <- block_diag(corX, c= 0.5)
  testthat::expect_equal(round(selective_p_val(S=cov(X), n=10, CP=block_diag_structure, c=0.5, k=1, d0=5, mc_iter=100),2), 0.27)
})

