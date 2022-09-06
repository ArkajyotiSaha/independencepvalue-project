# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Function for MC simulation for selective inference
#'
#' For a given group of variables \eqn{P}, simulates from the joint distribution of the canonical correlations and check if \eqn{P} is recovered by applying `block_diag` on the perturbed covariance matrix. This condition can be characterized as the vector of canonical correlations \eqn{\Lambda} belonging in a convex polytope \eqn{A\Lambda \leq b}, where the matrix \eqn{A} and the vector \eqn{b} are functions of the threshold \eqn{c} and the sample covariance matrix.
#' 
#' @param p p_1+p_2
#' @param rp min(p_1, p_2)
#' @param n sample size
#' @param A matrix with \code{rp} rows
#' @param b numeric vector
#' @return
#' \item{statistic}{Test statistic corresponding to simulated \eqn{\Lambda}.}
#' \item{status}{Logical vector if the simulated \eqn{\Lambda} satisfies \eqn{A\Lambda \leq b}.}
#' @keywords internal
MC_function_selective <- function(p, rp, n, A, b){
  if(nrow(A)!=length(b)){stop("error: number of rows of matrix A must be equal to the numebr of rows of vector b")}
  tilde_W_X <- stats::rWishart(1, p - rp, diag(rp))#simulate W
  tilde_T_X <- stats::rWishart(1, n - (p - rp) - 1, diag(rp))#simulate T
    
  tilde_F_X <- tilde_W_X[,,1] %*% solve(tilde_T_X[,,1])# inv(W)T
  eigen_tilde_F_X <- eigen(tilde_F_X) # eigen decomposition of (inv(W)T)
  F_X_eigenvalues <- eigen_tilde_F_X$values # eigenvalues of (inv(W)T)
  statistic <- 1/prod(1+F_X_eigenvalues)
  # prod(1 - lambda_i^2) = prod(1 - Psi_i/(1 + Psi_i)) = prod(1/(1 + Psi_i))
  status <- all(A %*% sqrt(F_X_eigenvalues/(1+F_X_eigenvalues)) <= b)
  return(list(statistic=statistic, status=status))
}
