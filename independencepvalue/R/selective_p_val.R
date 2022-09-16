# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Function to test the independence of a data-dependent group of variables with the remaining variables
#'
#' Given a covariance matrix \eqn{S} of \eqn{p} Gaussian variables and a grouping obtained via thresholding at \eqn{c} with \code{block_diag}, this tests the null hypothesis of independence between two groups of Gaussian variables.
#' 
#' @param S a \eqn{p \times p} covariance matrix
#' @param CP a vector of length \eqn{p} with \eqn{i^{th}} element denoting the 
#' group \eqn{i^{th}} variable belongs to
#' @param k the group to be tested for independence with the remaining variables, i.e. \eqn{P = [i : CP[i]==k]}
#' @param n sample size
#' @param c a threshold
#' @param d0 a natural number; if the number of canonical correlations is greater than \code{d0}, Monte Carlo simulation will be used to approximate the p-value for computational convenience; default value is 5
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
selective_p_val <- function(S, CP, k, n, c, d0=5, mc_iter= 1000){
  test_hyp <- test_stat_CCA(S, CP, k)
  p <- nrow(S)
  p2 <- nrow(test_hyp$S22)
  p1 <- p - p2
  if(p2 == 1){
      diag_S <- diag(1/sqrt(diag(S)))
      R <- diag_S %*% S %*% diag_S
      g_u <- min(1, c^2 * test_hyp$statistic/max(abs(R[CP==k, CP!=k]))^2)
      I_denom_tot <- stats::pbeta(0, g_u, min(g_u, 1 - test_hyp$statistic), (p - 3)/2 + 1, (n - p - 2)/2 + 1)
      du <- (I_denom_tot[2] - I_denom_tot[1])/(I_denom_tot[2] - I_denom_tot[3])
  }
  else{
    L <- matrix(0, (2 * p1 * p2 + 2 *p2), p2)
    for(i in 1:p1){
      for(j in 1:p2){
        temp <- test_hyp$left_SV[i,] * test_hyp$right_SV[j,]/sqrt(test_hyp$S11[i,i]*test_hyp$S22[j, j])
        L[2*((i-1) * p2 + j)-1, ] <- temp
        L[2*((i-1) * p2 + j), ] <- - temp
      }
    }
    for(k in 1:(p2-1)){
      L[2*p1*p2 + k,k:(k+1)] <- c(-1, 1)
    }
    L[2 * p1 * p2 + p2, p2] <- -1
    for(k in 1:(p2)){
      L[2*p1*p2 + p2 + k,k] <- 1
    }
    g <- c(rep(c, 2 * p1 * p2), rep(0, p2), rep(1, p2))
    du <- 0
    if(p2 <= d0){
      P <- rcdd::makeH(L, g, x = NULL)
      PV_d <- rcdd::scdd(P)
      V_d <- as.matrix(PV_d$output[ , - c(1, 2)])
      Pi <- volesti::Vpolytope(V = V_d)
      triang = try(geometry::delaunayn(Pi@V), silent = TRUE)
      if(!inherits(triang,'try-error')){
        prod_res <- 1
        for(i in 1:p2){
          prod_res <- prod_res * gamma((n - i)/2)/(gamma((n - p2 -  i)/2) * gamma((p1 - i + 1)/2) * gamma((p2 - i + 1)/2))
        }
        alpha <- log(pi^(p2/2) * 2^p *  prod_res)
        f_tot <- function(x){
          dt <- dCCA(n, p, p2, alpha, x)
          return(c(dt, dt* (prod(1 - x^2) <= test_hyp$statistic)))
          }
        part_int <- function(i, Pi, triang, f_tot){
          if(stats::var(round(Pi@V[triang[i,], ][,1], 5)) == 0) {
            Pi@V[triang[i,], ][,1][length(Pi@V[triang[i,], ][,1])] <- Pi@V[triang[i,],][,1][length(Pi@V[triang[i,], ][,1])]*(1 - 10^(-5))
          }
          return(SimplicialCubature::adaptIntegrateSimplex(f_tot, t(Pi@V[triang[i,], ]), fDim = 2)$integral)
        }
        par_I_tot_list <- future.apply::future_sapply(1:nrow(triang), part_int, Pi, triang, f_tot, future.seed=TRUE)
        du <- sum(par_I_tot_list[2,])/sum(par_I_tot_list[1,])
        if(sum(par_I_tot_list[1,]) == 0){du <- 0}
      }
    }
    if(du <= 0 || du >= 1 || p2 > d0){
      sip <- future.apply::future_sapply(1:mc_iter, function(i) MC_function_selective(p, p2, n, L, g), future.seed=TRUE)
      zp <- sum(as.numeric(sip[2,]))
      if(zp < 100){
        sip <- future.apply::future_sapply(1:(min((mc_iter * 100/zp), 100000)), function(i) MC_function_selective(p, p2, n, L, g), future.seed=TRUE)
      }
      du <- mean(test_hyp$statistic >= sip[1,][sip[2,] == TRUE])
    }
  }
  select_p_val <- du
  return(c(select_p_val))
}
