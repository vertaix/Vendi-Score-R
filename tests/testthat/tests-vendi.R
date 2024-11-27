library(VendiScore)
library(testthat)

K <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow=3)
test_that('Entropy_Q Returns Correct Output', {
  # For 3 elements with 0 similarity, vendi score should always return number of unique elements
  qs <- c(0.1, 0.5, 1., 2., -1)
  for(q in qs){
    if(q<0){
      q='inf'
    }
    result <- entropy_q(K/nrow(K), q=q)
    expect_equal(3, exp(result))
  }
})

samples <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow=3)
test_that('score_cosine Returns Correct Output', {
  #Similar to the previous test
  qs <- c(0.1, 0.5, 1., 2., -1)
  for(q in qs){
    if(q<0){
      q='inf'
    }
    result <- score_cosine(samples, q=q)
    expect_equal(3, result)
  }
})

samples <- matrix(c(1, 0, 0, 1, 0, 1, 1, 0), nrow=4)
k <- function(x, y){
  #For two-dimesional inputs, k(x,y) is 0.5*number of matching dimensions
  sim = 0
  if(x[1]==y[1]){
    sim=sim+0.5
  }
  if(x[2]==y[2]){
    sim=sim+0.5
  }
  return(sim)
}
test_that('score Returns Correct Output', {
  #We have a simple pair-wise metric above
  ## There are exactly two 'species'
  qs <- c(0.1, 0.5, 1., 2., -1)
  for(q in qs){
    if(q<0){
      q='inf'
    }
    result <- score(samples, k=k, q=q)
    expect_equal(2, result)
  }
})

samples1 <- matrix(c(1, 0, 0, 1, 0, 1, 1, 0, 0, 0), nrow=5, byrow=TRUE)
samples2 <- matrix(c(1, 0, 0, 1, 0, 1, 1, 0, 1, 0), nrow=5, byrow=TRUE)
test_that('score Returns Correct Output', {
  #A more complicated test now:
  ## In samples1, we add another sample that is distinct from all others with some similarity
  ## In samples2, we add a sample that is duplicate of an existing sample
  ## Vendi score should be always be less than 3 (given previous example had score 2)
  ## Samples1 should be more diverse than samples2 since it has one less exact-duplicate
  qs <- c(0.1, 0.5, 1., 2., -1)
  for(q in qs){
    if(q<0){
      q='inf'
    }
    result1 <- score(samples1, k=k, q=q)
    result2 <- score(samples2, k=k, q=q)
    expect_gt(3, result1)
    expect_gt(3, result2)
    expect_gt(result1, result2)
  }
})

a_3d <- array(c(0, 0, 1, 0, 0, 0,
                0, 1, 0, 0, 0, 0,
                0, 0, 0, 1, 0, 0),
              dim = c(3, 3, 2))
k_3d <- function(x,y){
  #Simple function that counts number of matching dimensions for 2d array objects
  sim <- 0
  for(i in 1:nrow(x)){
    for(j in 1:ncol(x)){
      if(x[i,j]==y[i,j]){
        sim <- sim + 1/6
      }
    }
  }
  return(sim)
}
test_that('Handling Higher-Dimensional Arrays', {
  #a_3d is a 3D array.
  ## score can handle such data structures
  ## similarity function yields that all samples are distinct
  expect_equal(1.98, score(a_3d, k_3d, q=1.), tolerance=1e-2)
})

weights <- c(0.1, 0.4, 0.5)
K <- matrix(c(1,0,0,0,1,0.5,0,0.5, 1), nrow=3, byrow=TRUE)
true_Kw <- matrix(c(0.1, 0, 0, 0, 0.4, 0.5*0.4**0.5 * 0.5**0.5, 0, 0.5*0.4**0.5 * 0.5**0.5, 0.5), nrow=3, byrow=TRUE)
test_that('weight_K Returns Correct Output', {
  #K is a matrix and weights are the corresponding weights we want to use in weighting
  pred_weightK <- weight_K(K, weights)
  for(i in 1:3){
    for(j in 1:3){
        expect_equal(pred_weightK[i,j], true_Kw[i,j])
    }
  }
})

K <- matrix(c(1,0,0,0,1,0,0,0, 1), nrow=3, byrow=TRUE)
test_that('score_K Returns Correct Output', {
  #K is a matrix and weights are the corresponding weights we want to use in weighting
  vendi_score <- score_K(K, q=1, p=NULL)
  expect_equal(3, 3)
})



