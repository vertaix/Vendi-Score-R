% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vendi.R
\name{score}
\alias{score}
\title{score}
\usage{
score(X, k, q = 1, p = NULL)
}
\arguments{
\item{X}{A \eqn{N} length datastructure that supports indexing so that
index i corresponds to the \code{i}th sample where \eqn{N} is the number of samples.
We support matrices, dataframes, vectors, lists, and higher-dimensional arrays.}

\item{k}{A positive semi-definite pairwise similarity function that takes two samples
and computes a similarity between \eqn{0} and \eqn{1}}

\item{q}{Vendi score order, must be non-negative}

\item{p}{A vector of probabilities (\code{Default=NULL}) to compute a
probability-weighted Vendi Score. \code{p=NULL} returns an un-weighted scores}
}
\value{
The Vendi score of order \eqn{q} for the provided data
}
\description{
Computes the Vendi score of a dataset given a pair-wise similarity
 function and an order \eqn{q}
}
\examples{
samples <- matrix(data=c(1,1,0,1,1,0,0,0,1), nrow=3, ncol=3)
rbf_kernel <- function(x, y, gamma = 1) exp(-gamma * sum((x - y)^2))
VS <- score(X=samples, k=rbf_kernel, q=1, p=NULL)
}
