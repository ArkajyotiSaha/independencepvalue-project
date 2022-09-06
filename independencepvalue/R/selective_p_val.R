# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Function to test the independence of a data-dependent group of variables with the remaining variables
#'
#' This tests the null hypothesis of independence between two groups of Gaussian variables, where the groups are obtained via thresholding with \code{block_diag}.
#' 
#' @param S the covariance matrix of the data matrix \eqn{X}, where \eqn{X} = (\eqn{X_1} \eqn{X_2})
#' @param n the number of data points
#' @param CP the grouping of the variables, an outcome of \code{block_diag}
#' @param c the threshold
#' @param k the group to be tested for independence with the remaining variables
#' @param d0 a natural number; ff the number of canonical correlations is greater than \code{d0}, Monte Carlo simulation will be used to approximate the p-value for computational convenience; default value is 5
#' @param mc_iter the number of Monte Carlo iterations used to approximate the p-value
#' @return The selective p-value for the test of independence.
#' @examples
#' # Simulates a 10 x 5 X from N(0, I)
#' set.seed(1)
#' X <- matrix(rnorm(50), 10, 5)
#'
#' # Compute the correlation matrix of X.
#' corX <- cor(X)
#' # Use 'block_diag' to obtain any block diagonal structure
#' block_diag_structure <- block_diag(corX, c= 0.5)
#' # test for independence of the variables in group 1 with the remaining variables
#' selective_p_val(S=cov(X), n=10, CP=block_diag_structure, c=0.5, k=1, d0=5, mc_iter=100)
#' @export 
selective_p_val <- function(S, n, CP, c, k, d0=5, mc_iter= 1000){
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
  inv_S11_x <- solve(S11_x)
  inv_S11_x_half <- amen::mhalf(inv_S11_x)
  inv_S22_x <- solve(S22_x)
  inv_S22_x_half <- amen::mhalf(inv_S22_x)
  tilde_S12_x <- inv_S11_x_half %*% S12_x %*% inv_S22_x_half
  svdecom <- svd(tilde_S12_x)
  singular_values <- svdecom$d
  test_stat <- prod(1-singular_values^2)
  A_x <- inv_S11_x_half %*% svdecom$u
  Gamma_x <- t(svdecom$v) %*% inv_S22_x_half
  L_x <- S11_x %*% A_x
  R_x <- Gamma_x %*% S22_x
  p1 <- max(nrow(L_x), ncol(R_x))
  p2 <- nrow(R_x)
  A <- matrix(0, (2 * p1 * p2 + 2 *p2), p2)
  for(i in 1:p1){
    for(j in 1:p2){
      A[2*((i-1) * p2 + j)-1, ] <- L_x[i,] * R_x[,j]/sqrt(S11_x[i,i]*S22_x[j, j])
      A[2*((i-1) * p2 + j), ] <- -L_x[i,] * R_x[,j]/sqrt(S11_x[i,i]*S22_x[j, j])
    }
  }
  if(p2 > 1){
    for(k in 1:(p2-1)){
      A[2*p1*p2 + k,k:(k+1)] <- c(-1, 1)
    }
  }
  A[2 * p1 * p2 + p2, p2] <- -1
  for(k in 1:(p2)){
    A[2*p1*p2 + p2 + k,k] <- 1
  }
  b <- c(rep(c, 2 * p1 * p2), rep(0, p2), rep(1, p2))
  du <- 0
  if(p2 <= d0){
    P <- rcdd::makeH(A, b, x = NULL)
    PV_d <- rcdd::scdd(P)
    if(p2 == 1){
      c1 <- max(PV_d$output[,3])^2
      I1 <- 0
      if(singular_values^2 < max(PV_d$output[,3])^2)
      {
        c1 <- singular_values^2
        I1_tot <- stats::pbeta(c( c1, max(PV_d$output[,3])^2), (p - 3)/2 + 1, (n - p - 2)/2 + 1)
        I1 <- I1_tot[2] - I1_tot[1]
      }
      I_denom_tot <- stats::pbeta(c(min(PV_d$output[,3])^2, max(PV_d$output[,3])^2), (p - 3)/2 + 1, (n - p - 2)/2 + 1)
      I_denom <- I_denom_tot[2] - I_denom_tot[1]
      du <- I1/I_denom
    }
    if(p2 > 1){
      V_d <- as.matrix(PV_d$output[ , - c(1, 2)])
      Pi <- volesti::Vpolytope(V = V_d)
      triang = try(geometry::delaunayn(Pi@V))
      prod_res <- 1
      for(i in 1:p2){
        prod_res <- prod_res * gamma((n - i)/2)/(gamma((n - p2 -  i)/2) * gamma((p1 - i + 1)/2) * gamma((p2 - i + 1)/2))
      }
      alpha <- log(pi^(p2/2) * 2^p *  prod_res)
      f_tot <- function(x) {
        dt <- dCCA(n, p, p2, alpha, x)
        return(c(dt, dt* (prod(1 - x^2) < test_stat)))
      }
      part_int <- function(i, Pi, triang, f_tot){
        if(stats::var(round(Pi@V[triang[i,], ][,1], 5)) == 0) Pi@V[triang[i,], ][,1][length(Pi@V[triang[i,], ][,1])] <- Pi@V[triang[i,], ][,1][length(Pi@V[triang[i,], ][,1])]*(1 - 10^(-5))
        return(SimplicialCubature::adaptIntegrateSimplex(f_tot, t(Pi@V[triang[i,], ]), fDim = 2)$integral)
      }
      par_I_tot_list <- future.apply::future_sapply(1:nrow(triang), part_int, Pi, triang, f_tot, future.seed=TRUE)
      du <- sum(par_I_tot_list[2,])/sum(par_I_tot_list[1,])
    }
  }
  if(du <= 0 || du >= 1 || p2 > d0){
    sip <- future.apply::future_sapply(1:mc_iter, function(i) MC_function_selective(p, p2, n, A, b), future.seed=TRUE)
    zp <- sum(as.numeric(sip[2,]))
    if(zp < 100){
      sip <- future.apply::future_sapply(1:(min((mc_iter * 100/zp), 100000)), function(i) MC_function_selective(p, p2, n, A, b), future.seed=TRUE)
    }
    du <- mean(test_stat >= sip[1,][sip[2,] == TRUE])
  }
  select_p_val <- du
  return(c(select_p_val))
}
