% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vendi.R
\name{score_K}
\alias{score_K}
\title{score_K}
\usage{
score_K(K, q = 1, p = NULL)
}
\arguments{
\item{K}{A \eqn{n x n} similarity matrix (must be Positive Semi-definite)}

\item{q}{Vendi score order, must be non-negative}

\item{p}{A vector of probabilities (\code{Default=NULL}) to compute a
probability-weighted Vendi Score. \code{p=NULL} returns an un-weighted score}
}
\value{
The Vendi score of order \eqn{q} for the provided data
}
\description{
Computes the Vendi score when provided with a pre-computed
 similarity matrix.
}
\examples{
K <- matrix(data=c(1,1,0,1,1,0,0,0,1), nrow=3, ncol=3)
VS <- score_K(K=K, q=1, p=NULL)
}
