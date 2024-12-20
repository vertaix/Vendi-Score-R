% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vendi.R
\name{score_cosine}
\alias{score_cosine}
\title{score_cosine}
\usage{
score_cosine(samples, q = 1, p = NULL)
}
\arguments{
\item{samples}{A \eqn{N x D} normalized matrix of data ($L2$ norm is $1$ for each sample)
where \eqn{N} is the number of samples and \eqn{D} is the dimensions}

\item{q}{Vendi score order, must be non-negative}

\item{p}{A vector of probabilities (\code{Default=NULL}) to compute a
probability-weighted Vendi Score. $\code{p=NULL}$ returns an un-weighted scores}
}
\value{
The Vendi score of order \eqn{q} for the data
}
\description{
Uses the cosine similarity trick to compute the Vendi Score efficiently
 for normalized numerical data
}
\examples{
samples <- matrix(sample(20, size=20, replace=TRUE, prob=NULL), nrow=10, ncol=2)
norm_samples <- t(apply(samples, 1, function(row) row / sqrt(sum(row^2))))
VS <- score_cosine(samples=norm_samples, q=1, p=NULL)
}
