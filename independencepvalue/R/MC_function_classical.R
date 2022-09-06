# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Function for MC simulation for classical inference
#'
#' Given \eqn{p = p_1 + p2}, and \code{rp} = \eqn{min(p_1, p_2)}, simulates from Wilks' lambda distribution with appropriate parameters.
#' 
#' @param p p_1+p_2
#' @param rp min(p_1, p_2)
#' @param n sample size
#' @return A sample from Wilks' lambda distribution
#' @keywords internal
MC_function_classical <- function(p, rp, n){
  tilde_W_X <- stats::rWishart(1, p - rp, diag(rp))#simulate W
  tilde_T_X <- stats::rWishart(1, n - (p - rp) - 1, diag(rp))#simulate T
  tilde_F_X <- tilde_W_X[,,1] %*% solve(tilde_T_X[,,1])# inv(W)T
  eigen_tilde_F_X <- eigen(tilde_F_X) # eigen decomposition of (inv(W)T)
  F_X_eigenvalues <- eigen_tilde_F_X$values # eigenvalues of (inv(W)T)
  statistic <- 1/prod(1+F_X_eigenvalues)
  # prod(1 - lambda_i^2) = prod(1 - Psi_i/(1 + Psi_i)) = prod(1/(1 + Psi_i))
  return(statistic)
}
