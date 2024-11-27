#' @title weight_K
#'
#' @description Computes a weighted similarity matrix.
#'
#' @param K A \eqn{n \times n} normalized similarity matrix (must be Positive Semi-definite)
#' @param p, an optional vector of length n must be non-negative and sum to 1. \code{Default=NULL}
#'
#' @return A weighted similarity matrix \code{K}
#'
#' @export
#'
weight_K <- function(K, p=NULL){
  if(is.null(p)){
    Ksim <- K/nrow(K)
  } else{ #compute the weighted vendi-score similarity matrix here
    if(sum(p)!=1){
      stop('Invalid Weights, must sum to 1')
    }
    Ksim <- K * outer(p**0.5, p**0.5)
  }
  return(Ksim)
}

#' @title entropy_q
#'
#' @description Returns the entropy of the eigenvalues of a given similarity matrix.
#'  This is equivalent to the log of the Vendi Score
#'
#' @param K A \eqn{n x n} normalized similarity matrix (must be Positive Semi-definite)
#' @param q Vendi score order, must be non-negative. \code{Default=1.}
#'
#' @return A single value describing the eigenvalues entropy of \code{K}
#'
#' @export

entropy_q <- function(K, q=1){
  ev <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  #we filter out exceedingly small values below the allowed precision
  #  Helps with numerical stability
  ev <- ev[ev>1e-15]

  if(q==1){
    return(sum(-1*(ev*log(ev))))
  }
  else if(q=='inf'){
    return(-log(max(ev)))
  }
  else{
    return(log(sum(ev**q))/(1-q))
  }
}

#' @title score_cosine
#'
#' @description Uses the cosine similarity trick to compute the Vendi Score efficiently
#'  for normalized numerical data
#' @param samples A \eqn{N x D} normalized matrix of data ($L2$ norm is $1$ for each sample)
#'  where \eqn{N} is the number of samples and \eqn{D} is the dimensions
#' @param q Vendi score order, must be non-negative
#' @param p A vector of probabilities (\code{Default=NULL}) to compute a
#'  probability-weighted Vendi Score. $\code{p=NULL}$ returns an un-weighted scores
#' @return The Vendi score of order \eqn{q} for the data
#' @examples
#' samples <- matrix(sample(20, size=20, replace=TRUE, prob=NULL), nrow=10, ncol=2)
#' norm_samples <- t(apply(samples, 1, function(row) row / sqrt(sum(row^2))))
#' VS <- score_cosine(samples=norm_samples, q=1, p=NULL)
#' @export

score_cosine <- function(samples, q=1, p=NULL){
  if(!is.matrix(samples)){
    stop('Data must be in matrix form')
  }
  if(nrow(samples)>ncol(samples)){
    if(is.null(p)){
      norm_samples <- samples/ (nrow(samples)**0.5)
      Ksim <- t(norm_samples) %*% norm_samples
    } else{
      norm_samples <- samples * p**0.5
      Ksim <- t(norm_samples) %*% norm_samples
    }
  } else{
    Ksim <- samples %*% t(samples)
    Ksim <- weight_K(Ksim, p=p)
  }
  return(exp(entropy_q(Ksim, q=q)))
}


.extract_sample <- function(X, i) {
  #Used to get the ith sample from X
  if (is.matrix(X) || is.data.frame(X)) {
    return(X[i, ])  # Extract row i for matrix/data.frame
  } else if (is.vector(X) || is.list(X)) {
    return(X[[i]])  # Extract ith element for list/vector
  } else if (is.array(X)) {
    # Generalize extraction for higher-dimensional arrays
    dims <- dim(X)
    slice <- as.list(rep(TRUE, length(dims)))  # Create a list of "TRUE" to index all other dimensions
    slice[[1]] <- i  # Select the ith slice along the first dimension
    return(do.call(`[`, c(list(X), slice)))  # Use do.call to subset the array
  }
}

#' @title score
#'
#' @description Computes the Vendi score of a dataset given a pair-wise similarity
#'  function and an order \eqn{q}
#' @param X A \eqn{N} length datastructure that supports indexing so that
#'  index i corresponds to the \code{i}th sample where \eqn{N} is the number of samples.
#'  We support matrices, dataframes, vectors, lists, and higher-dimensional arrays.
#' @param k A positive semi-definite pairwise similarity function that takes two samples
#'  and computes a similarity between \eqn{0} and \eqn{1}
#' @param q Vendi score order, must be non-negative
#' @param p A vector of probabilities (\code{Default=NULL}) to compute a
#'  probability-weighted Vendi Score. \code{p=NULL} returns an un-weighted scores
#' @return The Vendi score of order \eqn{q} for the provided data
#' @examples
#' samples <- matrix(data=c(1,1,0,1,1,0,0,0,1), nrow=3, ncol=3)
#' rbf_kernel <- function(x, y, gamma = 1) exp(-gamma * sum((x - y)^2))
#' VS <- score(X=samples, k=rbf_kernel, q=1, p=NULL)
#' @export

score <- function(X, k, q=1, p=NULL){
  if (is.matrix(X) || is.data.frame(X)) {
    n <- nrow(X)  # Number of rows (samples)
  } else if (is.vector(X) || is.list(X)) {
    n <- length(X)  # Number of samples
  } else if (is.array(X)) {
    n <- dim(X)[1]
  } else {
    stop("Unsupported data type for X")
  }

  Ksim <- matrix(nrow = n, ncol = n)
  for(i in 1:nrow(Ksim)){
    for(j in i:ncol(Ksim)){
      sample_i <- .extract_sample(X,i)
      sample_j <- .extract_sample(X,j)
      Ksim[i,j] = k(sample_i, sample_j)
      Ksim[j,i] = Ksim[i,j]
    }
  }
  Ksim <- weight_K(Ksim, p=p)
  return(exp(entropy_q(Ksim, q=q)))
}

#' @title score_K
#'
#' @description Computes the Vendi score when provided with a pre-computed
#'  similarity matrix.
#' @param K A \eqn{n x n} similarity matrix (must be Positive Semi-definite)
#' @param q Vendi score order, must be non-negative
#' @param p A vector of probabilities (\code{Default=NULL}) to compute a
#'  probability-weighted Vendi Score. \code{p=NULL} returns an un-weighted score
#' @return The Vendi score of order \eqn{q} for the provided data
#' @examples
#' K <- matrix(data=c(1,1,0,1,1,0,0,0,1), nrow=3, ncol=3)
#' VS <- score_K(K=K, q=1, p=NULL)
#' @export

score_K <- function(K, q=1, p=NULL){
  if (!is.matrix(K)){
    stop("Unsupported data type for K")
  }
  Ksim <- weight_K(K, p=p)
  return(exp(entropy_q(Ksim, q=q)))
}
