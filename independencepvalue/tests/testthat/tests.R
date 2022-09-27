# Generated from create-independencepvalue.Rmd: do not edit by hand  
testthat::test_that("test_stat_CCA() works", {
  set.seed(1)
  X <- matrix(rnorm(50), 10, 5)
  testthat::expect_equal(
    round(test_stat_CCA(S=cov(X), CP=rep(1:2, times=c(3, 2)), k=1)$statistic, 2),
    0.62)
  set.seed(1)
  X <- matrix(rnorm(50), 10, 5)
  # testing group 2 should give identical results:
  testthat::expect_equal(
    round(test_stat_CCA(S=cov(X), CP=rep(1:2, times=c(3, 2)), k=2)$statistic, 2),
    0.62)
})

testthat::test_that("MC_function_classical() works", {
  set.seed(123)
  p <- 5
  n <- 20
  nsim <- 1e4
  from_wisharts <- sapply(1:nsim,
                          function(i) MC_function_classical(p = p, rp = 1, n = n))
  from_beta <- 1 - rbeta(nsim, (p - 1) / 2, (n - p) / 2)
  probs <- seq(0.05, 0.95, length = 10)
  qw <- quantile(from_wisharts,probs = probs)
  qb <- quantile(from_beta,probs = probs)
  testthat::expect_true( all(abs(qw - qb) < 0.01) )
})

testthat::test_that("classical_p_val() works", {
  set.seed(1)
  X <- matrix(rnorm(50), 10, 5)
  testthat::expect_equal(
    classical_p_val(S=cov(X), n=10, CP=rep(1:2, times=c(3, 2)), k=1, mc_iter=100),
    0.83)
  set.seed(1)
  X <- matrix(rnorm(50), 10, 5)
  # testing group 2 should give identical results:
  testthat::expect_equal(
    classical_p_val(S=cov(X), n=10, CP=rep(1:2, times=c(3, 2)), k=2, mc_iter=100),
    0.83)
})

testthat::test_that("classical_p_val() works for $r=1$", {
  set.seed(1)
  X <- matrix(rnorm(50), 10, 5)
  testthat::expect_equal(
    round(classical_p_val(S=cov(X), n=10, CP=rep(1:2, times=c(1, 4)), k=1, mc_iter=100),2),
    0.31)
  set.seed(1)
  X <- matrix(rnorm(50), 10, 5)
  # testing group 2 should give identical results:
  testthat::expect_equal(
    round(classical_p_val(S=cov(X), n=10, CP=rep(1:2, times=c(1, 4)), k=2, mc_iter=100),2),
    0.31)
})

testthat::test_that("block_diag() works", {
  set.seed(1)
  X <- matrix(rnorm(50), 10, 5)
  testthat::expect_equal(length(unique(block_diag(cor(X), c=0.5))), 3)
})

testthat::test_that("selective_p_val() works", {
  set.seed(4)
  X <- matrix(rnorm(600), 30, 20)
  corX <- cor(cbind(X))
  block_diag_structure <- block_diag(corX, c=0.3)
  # Testing for r=1. 
  testthat::expect_equal(round(selective_p_val(S=cov(X), CP=block_diag_structure, k=6, n=30, c=0.3, d0=5, mc_iter=100),2), 0.81)
  # Testing for small r. 
  testthat::expect_equal(round(selective_p_val(S=cov(X), CP=block_diag_structure, k=2, n=30, c=0.3, d0=5, mc_iter=100),2), 0.15)
  # Testing for large r. 
  set.seed(1)
  testthat::expect_equal(round(selective_p_val(S=cov(X), CP=block_diag_structure, k=3, n=30, c=0.3, d0=5, mc_iter=100),2), 0.77)
})

testthat::test_that("selective_p_val_beta() and selective_p_val_MC() produce similar results", {
set.seed(4)
X <- matrix(rnorm(600), 30, 20)
corX <- cor(cbind(X))
block_diag_structure <- block_diag(corX, c=0.3)
set.seed(1)
p_beta <- selective_p_val(S=cov(X), CP=block_diag_structure, k=6, n=30, c=0.3, d0=5, mc_iter=100)
set.seed(1)
p_MC <- selective_p_val(S=cov(X), CP=block_diag_structure, k=6, n=30, c=0.3, d0=0, mc_iter=100)
testthat::expect_true(abs(p_beta - p_MC) < 0.1)
})

testthat::test_that("selective_p_val_numerical() and selective_p_val_MC() produce similar results", {
set.seed(4)
X <- matrix(rnorm(600), 30, 20)
corX <- cor(cbind(X))
block_diag_structure <- block_diag(corX, c=0.3)
set.seed(1)
p_numerical <- selective_p_val(S=cov(X), CP=block_diag_structure, k=2, n=30, c=0.3, d0=5, mc_iter=100)
set.seed(1)
p_MC <- selective_p_val(S=cov(X), CP=block_diag_structure, k=2, n=30, c=0.3, d0=0, mc_iter=1000)
testthat::expect_true(abs(p_numerical - p_MC) < 0.1)
})

