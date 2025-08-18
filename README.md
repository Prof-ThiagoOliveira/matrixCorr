
<!-- README.md is generated from README.Rmd. Please edit that file -->

# matrixCorr

<!-- badges: start -->

[![R-CMD-check.yaml](https://github.com/Prof-ThiagoOliveira/matrixCorr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Prof-ThiagoOliveira/matrixCorr/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/Prof-ThiagoOliveira/kendall_tau_rank_corr/graph/badge.svg)](https://app.codecov.io/gh/Prof-ThiagoOliveira/kendall_tau_rank_corr)
<!-- badges: end -->

`matrixCorr` provides fast and scalable computation of correlation
matrices using optimized C++ routines. It supports rank-based and
parametric methods, including Kendall’s tau (tau-a and tau-b), with
planned support for Spearman and Pearson correlations.

## Features

- High-performance C++ backend using `Rcpp`
- Support for tied data via Kendall’s tau-b
- Easy-to-use `kendall_tau()` function for matrices and data frames
- ggplot2-based heatmap visualization
- Designed for large-scale, high-dimensional data

## Installation

``` r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("Prof-ThiagoOliveira/matrixCorr")
```

## Example

### **Kendall’s Tau Example**

``` r
library(matrixCorr)

# Simulated data
set.seed(42)
mat <- cbind(A = rnorm(100), B = rnorm(100), C = rnorm(100))

# Compute Kendall's tau correlation matrix
ktau <- kendall_tau(mat)

# Print matrix
print(ktau)

# Visualize with ggplot2
plot(ktau)
```

------------------------------------------------------------------------

### \*Spearman’s Rho Example\*\*

``` r
library(matrixCorr)

# Simulated data with some ties
set.seed(123)
mat <- cbind(
  A = sample(1:10, 100, replace = TRUE),
  B = sample(1:10, 100, replace = TRUE),
  C = rnorm(100)
)

# Compute Spearman's rho correlation matrix
spearman <- spearman_rho(mat)

# Print matrix
print(spearman)

# Visualize with ggplot2
plot(spearman)
```

------------------------------------------------------------------------

### **Pearson Correlation Example**

``` r
library(matrixCorr)

# Simulated continuous data
set.seed(999)
mat <- cbind(
  A = rnorm(100),
  B = 0.5 * rnorm(100) + 0.5,
  C = runif(100)
)

# Compute Pearson correlation matrix
pcorr <- pearson_corr(mat)

# Print matrix
print(pcorr)

# Visualize with ggplot2
plot(pcorr)
```

## License

MIT [Thiago de Paula Oliveira](https://orcid.org/0000-0002-4555-2584)
