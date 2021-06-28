
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sparsediscrim

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Codecov test
coverage](https://codecov.io/gh/topepo/sparsediscrim/branch/main/graph/badge.svg)](https://codecov.io/gh/topepo/sparsediscrim?branch=main)
[![R-CMD-check](https://github.com/topepo/sparsediscrim/workflows/R-CMD-check/badge.svg)](https://github.com/topepo/sparsediscrim/actions)
<!-- badges: end -->

The R package `sparsediscrim` provides a collection of sparse and
regularized discriminant analysis classifiers that are especially useful
for when applied to small-sample, high-dimensional data sets.

The package was archived in 2018 and was re-released in 2021. The
package code was forked from [John Ramey’s
repo](https://github.com/ramhiser/sparsediscrim) and subsequently
modified.

## Installation

You can install the stable version on
[CRAN](https://cran.r-project.org/package=sparsediscrim):

``` r
install.packages('sparsediscrim', dependencies = TRUE)
```

If you prefer to download the latest version, instead type:

``` r
library(devtools)
install_github('topepo/sparsediscrim')
```

## Usage

The formula and non-formula interfaces can be used:

``` r
library(sparsediscrim)

data(parabolic, package = "modeldata")

qda_mod <- qda_shrink_mean(class ~ ., data = parabolic)
# or
qda_mod <- qda_shrink_mean(x = parabolic[, 1:2], y = parabolic$class)

qda_mod
#> Shrinkage-Mean-Based Diagonal QDA
#> 
#> Sample Size: 500 
#> Number of Features: 2 
#> 
#> Classes and Prior Probabilities:
#>   Class1 (48.8%), Class2 (51.2%)

# Prediction uses the `type` argument: 

parabolic_grid <-
   expand.grid(X1 = seq(-5, 5, length = 100),
               X2 = seq(-5, 5, length = 100))


parabolic_grid$qda <- predict(qda_mod, parabolic_grid, type = "prob")$Class1

library(ggplot2)
ggplot(parabolic, aes(x = X1, y = X2)) +
   geom_point(aes(col = class), alpha = .5) +
   geom_contour(data = parabolic_grid, aes(z = qda), col = "black", breaks = .5) +
   theme_bw() +
   theme(legend.position = "top") +
   coord_equal()
```

<img src="man/figures/README-usage-1.png" width="100%" />

## Classifiers

The `sparsediscrim` package features the following classifier (the R
function is included within parentheses):

-   [High-Dimensional Regularized Discriminant
    Analysis](https://arxiv.org/abs/1602.01182) (`rda_high_dim()`) from
    Ramey et al. (2015)

The `sparsediscrim` package also includes a variety of additional
classifiers intended for small-sample, high-dimensional data sets. These
include:

| Classifier                                                    | Author                                                                                         | R Function                                                                                           |
|---------------------------------------------------------------|------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------|
| Diagonal Linear Discriminant Analysis                         | [Dudoit et al. (2002)](https://www.tandfonline.com/doi/abs/10.1198/016214502753479248)         | [`lda_diag()`](https://topepo.github.io/sparsediscrim/reference/lda_diag.html)                       |
| Diagonal Quadratic Discriminant Analysis                      | [Dudoit et al. (2002)](https://www.tandfonline.com/doi/abs/10.1198/016214502753479248)         | [`qda_diag()`](https://topepo.github.io/sparsediscrim/reference/qda_diag.html)                       |
| Shrinkage-based Diagonal Linear Discriminant Analysis         | [Pang et al. (2009)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1541-0420.2009.01200.x) | [`lda_shrink_cov()`](https://topepo.github.io/sparsediscrim/reference/lda_shrink_cov.html)           |
| Shrinkage-based Diagonal Quadratic Discriminant Analysis      | [Pang et al. (2009)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1541-0420.2009.01200.x) | [`qda_shrink_cov()`](https://topepo.github.io/sparsediscrim/reference/qda_shrink_cov.html)           |
| Shrinkage-mean-based Diagonal Linear Discriminant Analysis    | [Tong et al. (2012)](https://academic.oup.com/bioinformatics/article/28/4/531/211887)          | [`lda_shrink_mean()`](https://topepo.github.io/sparsediscrim/reference/lda_shrink_mean.html)         |
| Shrinkage-mean-based Diagonal Quadratic Discriminant Analysis | [Tong et al. (2012)](https://academic.oup.com/bioinformatics/article/28/4/531/211887)          | [`qda_shrink_mean()`](https://topepo.github.io/sparsediscrim/reference/qda_shrink_mean.html)         |
| Minimum Distance Empirical Bayesian Estimator (MDEB)          | [Srivistava and Kubokawa (2007)](http://www.utstat.utoronto.ca/~srivasta/exp1.pdf)             | [`lda_emp_bayes()`](https://topepo.github.io/sparsediscrim/reference/lda_emp_bayes.html)             |
| Minimum Distance Rule using Modified Empirical Bayes (MDMEB)  | [Srivistava and Kubokawa (2007)](http://www.utstat.utoronto.ca/~srivasta/exp1.pdf)             | [`lda_emp_bayes_eigen()`](https://topepo.github.io/sparsediscrim/reference/lda_emp_bayes_eigen.html) |
| Minimum Distance Rule using Moore-Penrose Inverse (MDMP)      | [Srivistava and Kubokawa (2007)](http://www.utstat.utoronto.ca/~srivasta/exp1.pdf)             | [`lda_eigen()`](https://topepo.github.io/sparsediscrim/reference/lda_eigen.html)                     |

We also include modifications to Linear Discriminant Analysis (LDA) with
regularized covariance-matrix estimators:

-   Moore-Penrose Pseudo-Inverse
    ([`lda_pseudo()`](https://topepo.github.io/sparsediscrim/reference/lda_pseudo.html))
-   Schafer-Strimmer estimator
    ([`lda_schafer()`](https://topepo.github.io/sparsediscrim/reference/lda_schafer.html))
-   Thomaz-Kitani-Gillies estimator
    ([`lda_thomaz()`](https://topepo.github.io/sparsediscrim/reference/lda_thomaz.html))
