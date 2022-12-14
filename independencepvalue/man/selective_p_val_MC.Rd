% Generated by roxygen2: do not edit by hand
% Please edit documentation in create-independencepvalue.Rmd.
\name{selective_p_val_MC}
\alias{selective_p_val_MC}
\title{Compute selective p-value using Monte Carlo}
\usage{
selective_p_val_MC(n, L, g, test_hyp, mc_iter)
}
\arguments{
\item{n}{sample size}

\item{L}{matrix used to define conditioning set}

\item{g}{vector used to define conditioning set}

\item{test_hyp}{output of \code{test_stat_CCA()}}

\item{mc_iter}{number of Monte Carlo iterations}
}
\description{
Compute selective p-value using Monte Carlo
}
\keyword{internal}