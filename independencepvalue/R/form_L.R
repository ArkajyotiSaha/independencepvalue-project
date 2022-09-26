# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Form L matrix described in Theorem 2
#' 
#' @param test_hyp output of `test_stat_CCA()`
#' 
#' @keywords internal
form_L <- function(test_hyp) {
  p1 <- nrow(test_hyp$S11)
  p2 <- nrow(test_hyp$S22)
  
  L <- matrix(0, (2 * p1 * p2 + 2 * p2), p2)
  for (i in 1:p1) {
    for (j in 1:p2) {
      temp <-
        test_hyp$left_SV[i, ] * test_hyp$right_SV[j, ] / sqrt(test_hyp$S11[i, i] *
                                                                test_hyp$S22[j, j])
      L[2 * ((i - 1) * p2 + j) - 1,] <- temp
      L[2 * ((i - 1) * p2 + j),] <- -temp
    }
  }
  for (k in 1:(p2 - 1)) {
    L[2 * p1 * p2 + k, k:(k + 1)] <- c(-1, 1)
  }
  L[2 * p1 * p2 + p2, p2] <- -1
  for (k in 1:(p2)) {
    L[2 * p1 * p2 + p2 + k, k] <- 1
  }
  return(L)
}
