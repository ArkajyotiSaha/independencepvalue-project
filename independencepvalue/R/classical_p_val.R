# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Function to test the independence of two pre-specified groups of variables
#'
#' This tests the null hypothesis of independence between two pre-specified groups of Gaussian variables. This function approximates the p-value corresponding to the Wilk's lambda statistics with Monte Carlo simulation with \code{mc_iter} iterations.
#' 
#' @param S the covariance matrix of the data matrix \eqn{X}, where \eqn{X} = (\eqn{X_1} \eqn{X_2})
#' @param n the number of data points
#' @param CP the pre-specified grouping of the variables
#' @param k the group to be tested for independence with the remaining variables
#' @param mc_iter the number of Monte Carlo iterations to approximate the p-value
#' @return The p-value for the test of independence. 
#' @examples
#' # Simulates a 10 x 3 X_1 from N(0, I)
#' set.seed(1)
#' X_1 <- matrix(rnorm(30), 10, 3)
#'
##' # Simulates a 10 x 2 X_2 from N(0, I) independently of X_1
#' set.seed(2)
#' X_2 <- matrix(rnorm(20), 10, 2)
#'
#' # Compute the covariance matrix of X = (X_1 X_2).
#' covX <- cov(cbind(X_1, X_2))
#' # tests for a difference in means between X_1 and X_2
#' classical_p_val(S=covX, n=10, CP=c(rep(1, each = 3),rep(2, each = 2)), k=1, mc_iter=100)
#' @export 
classical_p_val <- function(S, n, CP, k, mc_iter= 1000){
  p <- nrow(S)
  ptemp <- sum(CP==k)
  if(2*ptemp >= p){
    p1 <- ptemp
    S11_x <- as.matrix(S[which(CP==k), which(CP==k)])#S11
    S22_x <- as.matrix(S[which(CP!=k), which(CP!=k)])#S22
    S12_x <- as.matrix(S[which(CP==k), which(CP!=k)])#S_12 
  }
  if(2*ptemp < p){
    p1 = p - ptemp
    S11_x <- as.matrix(S[which(CP!=k), which(CP!=k)])#S11
    S22_x <- as.matrix(S[which(CP==k), which(CP==k)])#S22
    S12_x <- as.matrix(S[which(CP!=k), which(CP==k)])#S_12 
  }
  #ensures that p2 <= p1.
  p2 <- p - p1
  if(p1 > p2){rp <- p2}#Define r(P)
  if(p1 <= p2){rp <- p1} 
  inv_S11_x <- solve(S11_x)#inv(S11)
  inv_S11_x_half <- amen::mhalf(inv_S11_x)#S_11^{-1/2} 
  inv_S22_x <- solve(S22_x)#inv(S22)
  inv_S22_x_half <- amen::mhalf(inv_S22_x)#S_22^{-1/2} 
  tilde_S12_x <- inv_S11_x_half %*% S12_x %*% inv_S22_x_half#S_12^W = covariance matrix of whitened X_1 and X_2
  svdecom <- svd(tilde_S12_x)#compact SVD 
  singular_values <- svdecom$d#Lambda
  test_stat <- prod(1-singular_values^2)#test statistic 
  if(rp == 1){
    classic_p_val <- 1 - stats::pbeta(singular_values^2, (p - 3)/2 + 1, (n - p - 2)/2 + 1) 
    # for r(p) = 1, we have simple evaluation of p_{LRT}.
  }
  if(rp > 1){
    sip <- future.apply::future_sapply(1:mc_iter, function(i) MC_function_classical(p, rp, n), future.seed = TRUE)
    #MC simulation to approximate the p-value. Simulates W and T. Returns vector of statistic. 
    classic_p_val <- mean(test_stat >= sip)#computes p-value. The test statistic is smaller if it is away from null. 
  }
  return(classic_p_val)
}
