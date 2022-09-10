# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Sample from the distribution of canonical correlations under the null 
#' 
#' Samples from the `rp`-dimensional joint distribution of the canonical 
#' correlations between two groups of variables of size `p1` and `p2`, where 
#' `rp` is the smaller of `p1` and `p2` and `p=p1+p2`, and `n` is the sample size.
#' 
#' @param p p_1+p_2
#' @param rp min(p_1, p_2)
#' @param n sample size
#' @return A vector of length `rp` sampled from the joint distribution described
#' above.
#' @keywords internal
sample_canoncor_under_null <- function(p, rp, n) {
  tilde_W_X <- stats::rWishart(1, p - rp, diag(rp))#simulate W
  tilde_T_X <- stats::rWishart(1, n - (p - rp) - 1, diag(rp))#simulate T
  tilde_F_X <- tilde_W_X[,,1] %*% solve(tilde_T_X[,,1])# inv(W)T
  return(eigen(tilde_F_X)$values) # eigenvalues of (inv(W)T)
}
