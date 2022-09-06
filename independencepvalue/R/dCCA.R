# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Function to evaluate  the joint density of the canonical correlations
#'
#' This function evaluates the joint density of the canonical correlations at a specific value \eqn{Lambda}.
#' 
#' @param n sample size
#' @param p total number of variables
#' @param rp minimum of the size of the two groups of variables. 
#' @param a a scaling constant
#' @param Lambda the \code{rp} length vector where the joint density is to be evaluatedl \eqn{1 > Lambda[1]\ge Lambda[2] \ge ... \ge Lambda[rp] > 0} 
#' @return Scaled joined density of canonial correlations, evaluated at \eqn{Lambda}.
#' @keywords internal
dCCA <- function(n, p, rp, a, Lambda){
  c1 <- (p - 2*rp) * sum(log(Lambda))
  c2 <- ((n - p - 2)/2) * sum(log(1-Lambda^2))
  dip <- matrix(1, length(Lambda), length(Lambda))
  dip[lower.tri(dip)] <- stats::dist(Lambda^2)
  c3 <- sum(log(matrixStats::colProds(dip)))
  return(exp(a+c1+c2+c3))
}
