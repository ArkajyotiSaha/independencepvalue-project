# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Function to obtain block diagonal structure through thresholding
#'
#' Given the correlation matrix matrix of a data matrix \eqn{X}, this function discovrs the block-diagonal structure by thresholding the absolute values of the entries of the correlation matrix at c. We create an adjacency matrix with the elements being 1 iff the corresponding member of the correlation matrix has an absolute value greater than or equal to \eqn{c}. This is equivalent to performing a single linkage hierarcical clustering on the variables, with the distance matrix given by \eqn{1 - |R|} and truncating the tree at \eqn{1-c}.
#' 
#' @param R the \eqn{p \time p} correlation matrix of the data
#' @param c the threshold
#' @return A \eqn{p} length integer vector whose \eqn{i^{th}} element denotes the group \eqn{i^{th}} variable belongs to.
#' @examples
#' # Simulates a 10 x 5 X from N(0, I)
#' set.seed(1)
#' X <- matrix(rnorm(50), 10, 5)
#' # Compute the coreelation matrix of X.
#' corX <- cor(X)
#' # Compute the block diagonal structure at c=0.5
#' block_diag(R=corX, c=0.5)
#' @export 
block_diag <- function(R, c){
  dis_R <- 1 - abs(R)
  test <- stats::as.dist(dis_R, diag = TRUE)
  clust_result <- stats::hclust(test, method = "single")
  return(stats::cutree(clust_result, h=(1 - c)))
}
