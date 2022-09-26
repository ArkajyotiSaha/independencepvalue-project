# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Compute selective p-value using numerical integration
#' 
#' @param S a \eqn{p \times p} covariance matrix
#' @param CP a vector of length \eqn{p} with \eqn{i^{th}} element denoting the 
#' group \eqn{i^{th}} variable belongs to
#' @param k the group to be tested for independence with the remaining variables, i.e. \eqn{P = [i : CP[i]==k]}
#' @param n sample size
#' @param c a threshold
#' @param test_hyp output of `test_stat_CCA()`
#' 
#' @keywords internal
selective_p_val_beta <- function(S, CP, k, n, c, test_hyp) {
  p <- nrow(S)
  diag_S <- diag(1 / sqrt(diag(S)))
  R <- diag_S %*% S %*% diag_S
  g_u <-
    min(1, c ^ 2 * test_hyp$statistic / max(abs(R[CP == k, CP != k])) ^ 2)
  message("ARKA: Wouldn't it be simpler to write (p - 1)/2 and (n - p)/2 in the next line?")
  I_denom_tot <- stats::pbeta(0, g_u, min(g_u, 1 - test_hyp$statistic), (p - 3) / 2 + 1, (n - p - 2) / 2 + 1)
  (I_denom_tot[2] - I_denom_tot[1]) / (I_denom_tot[2] - I_denom_tot[3])
}
