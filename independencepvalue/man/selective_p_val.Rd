% Generated by roxygen2: do not edit by hand
% Please edit documentation in create-independencepvalue.Rmd.
\name{selective_p_val}
\alias{selective_p_val}
\title{Test independence of a data-dependent group of variables with the rest}
\usage{
selective_p_val(
  S,
  CP,
  k,
  n,
  c,
  d0 = 5,
  tol = 1e-05,
  maxeval = 1e+05,
  mc_iter = 1000
)
}
\arguments{
\item{S}{a \eqn{p \times p} covariance matrix}

\item{CP}{a vector of length \eqn{p} with \eqn{i^{th}} element denoting the
group \eqn{i^{th}} variable belongs to}

\item{k}{the group to be tested for independence with the remaining variables, i.e. \eqn{P = [i : CP[i]==k]}}

\item{n}{sample size}

\item{c}{a threshold}

\item{d0}{a natural number; if the number of canonical correlations is greater than \code{d0}, Monte Carlo simulation will be used to approximate the p-value for computational convenience; default value is 5}

\item{tol}{the relative tolerance used to approximate the p-value using \code{selective_p_val_integrate()}; default value is 1e-05}

\item{maxeval}{the maximum number of function evaluations used to approximate the p-value using \code{selective_p_val_integrate()}; we recommend using a high value of this to obtain an approximation with high accuracy; default value is 10,000}

\item{mc_iter}{the number of Monte Carlo iterations used to approximate the p-value; we recommend using a high value of this to obtain an approximation with high accuracy; default value is 1,000}
}
\value{
The selective p-value for the test of independence.
}
\description{
Given a covariance matrix \code{S} of \code{p} Gaussian variables and a grouping obtained
via thresholding absolute correlations at \code{1-c} using \code{block_diag()}, this
function tests the null hypothesis of independence between two groups of
Gaussian variables.
}
\examples{
# Simulates a 10 x 5 X from N(0, I)
set.seed(1)
X <- matrix(rnorm(50), 10, 5)

# Compute the correlation matrix of X.
corX <- cor(X)
# Use 'block_diag' to obtain any block diagonal structure
block_diag_structure <- block_diag(corX, c= 0.5)
# test for independence of the variables in group 1 with the remaining variables
selective_p_val(S=cov(X), n=10, CP=block_diag_structure, c=0.5, 
k=1, d0=5, tol = 1e-05, maxeval = 10000, mc_iter=100)
}