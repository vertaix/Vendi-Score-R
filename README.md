# R Implementation of the Vendi Score: A Diversity Evaluation Metric for Machine Learning and Science

This repository contains the R implementation of the Vendi Score (VS), a metric for evaluating diversity in machine learning and the natural sciences. The Vendi Scores are a family of diversity metrics that are flexible, interpretable, and unsupervised. Defined as the entropy associated with the eigenvalues of a sample similarity matrix $\mathbf{K}$, the Vendi Score only requires a pair-wise similarity measure. The Vendi Score is defined in this [paper](https://arxiv.org/abs/2210.02410). The order $q$ of the Vendi Score weight the importance of rare and common elements in the diversity computation, as described in this [paper](https://arxiv.org/abs/2310.12952).

Check out our Python implementation of the Vendi Scores [here](https://github.com/vertaix/Vendi-Score/tree/main)!

## Installation

The Vendi Scores in R require no additional dependencies and can be directly installed with devtools.

``` r
devtools::install_github("vertaix/Vendi-Score-R")
```

## Computing the Vendi Scores

The Vendi Scores have 3 inputs: your data, a pair-wise similarity metric $k$, and an order $q$.

The data can be any data frame, matrix, vector, list or higher-dimensional array for which index $i$ corresponds to the $i$ th sample. The pair-wise similarity function $k$ should be symmetric and $k(x,x)=1$. The order $q$ can be any non-negative value.

``` r
library(VendiScore)
library(datasets)
data(iris)
iris_mat <- data.matrix(iris)
iris_mat <- iris_mat[,colnames(iris_mat)!='Species']

rbf_kernel <- function(x, y, gamma = 0.1) exp(-gamma * sum((x - y)^2))

score(iris_mat, rbf_kernel, q=1.)
# 3.169735
```

A score of about $3$ was expected since we have data from $3$ species, but we did not need to have class labels to measure the diversity of the dataset.

We can also use the cosine kernel trick to speed up Vendi Score computation for larger datasets in numerical form that can use a cosine similarity kernel. Data must be normalized in this case.

``` r
norm_samples <- t(apply(iris_mat, 1, function(row) row / sqrt(sum(row^2))))
VS <- score_cosine(samples=norm_samples, q=1)
# 1.20783
```

In cases where already have pre-computed a similarity matrix:
``` r
K <- matrix(data=c(1,1,0,1,1,0,0,0,1), nrow=3, ncol=3)
VS <- score_K(K, q=1.)
# 1.88988
```
We provide documentation for all functions in the package.

Check out our vignette for a demonstration of the advantages of the Vendi Score over metrics like average similarity.

## Citation

``` bibtex
@article{friedman2022vendi,
  title={The Vendi Score: A Diversity Evaluation Metric for Machine Learning},
  author={Friedman, Dan and Dieng, Adji Bousso},
  journal={arXiv preprint arXiv:2210.02410},
  year={2022}
}
```

``` bibtex
@article{pasarkar2023cousins,
      title={Cousins Of The Vendi Score: A Family Of Similarity-Based Diversity Metrics For Science And Machine Learning}, 
      author={Pasarkar, Amey P and Dieng, Adji Bousso},
      journal={arXiv preprint arXiv:2310.12952},
      year={2023},
}
```
