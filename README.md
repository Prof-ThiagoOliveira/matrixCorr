
<!-- README.md is generated from README.Rmd. Please edit that file -->

# matrixCorr <img src="man/figures/logo.png" align="right" height="160" />

<!-- badges: start -->

[![CRAN
version](https://www.r-pkg.org/badges/version/matrixCorr)](https://CRAN.R-project.org/package=matrixCorr)
![CRAN
downloads](https://cranlogs.r-pkg.org/badges/grand-total/matrixCorr)
[![R-CMD-check.yaml](https://github.com/Prof-ThiagoOliveira/matrixCorr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Prof-ThiagoOliveira/matrixCorr/actions/workflows/R-CMD-check.yaml)
[![test-coverage.yaml](https://github.com/Prof-ThiagoOliveira/matrixCorr/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/Prof-ThiagoOliveira/matrixCorr/actions/workflows/test-coverage.yaml)
<!-- badges: end -->

`matrixCorr` computes correlation and related association matrices from
small to high-dimensional data using simple, consistent functions and
sensible defaults. It includes shrinkage and robust options for noisy or
**p \>= n** settings, plus convenient print/plot/summary methods.
Performance-critical paths are implemented in C++ with BLAS/OpenMP and
memory-aware symmetric updates. The API accepts base matrices and data
frames and returns standard R objects via a consistent S3 interface.

Contributions from other researchers who want to add new correlation
methods are very welcome. A central goal of `matrixCorr` is to keep
efficient correlation and agreement estimation in one package with a
common interface and consistent outputs, so methods can be extended,
compared, and used without repeated translation across packages.

Supported measures include Pearson, Spearman, Kendall, distance
correlation, partial correlation, robust biweight mid-correlation,
percentage bend, Winsorized, skipped correlation, and latent
categorical/ordinal correlations (tetrachoric, polychoric, polyserial,
and biserial), plus repeated-measures correlation (`rmcorr()`);
agreement tools cover Bland-Altman (two-method and repeated-measures),
Lin’s concordance correlation coefficient (including repeated-measures
LMM/REML extensions), and intraclass correlation for both wide and
repeated-measures designs.

## Features

- High-performance C++ backend using `Rcpp`
- General correlations such as `pearson_corr()`, `spearman_rho()`,
  `kendall_tau()`
- Robust correlation metrics (`bicor()`, `pbcor()`, `wincor()`,
  `skipped_corr()`)
- Distance correlation (`dcor()`)
- Partial correlation (`pcorr()`)
- Latent categorical/ordinal correlations (`tetrachoric()`,
  `polychoric()`, `polyserial()`, `biserial()`)
- Repeated-measures correlation (`rmcorr()`)
- Shrinkage for $p >> n$ (`shrinkage_corr()`)
- Agreement metrics
  - Bland-Altman (two-method `ba()` and repeated-measures `ba_rm()`),
  - Lin’s concordance correlation coefficient (pairwise `ccc()`,
    repeated-measures LMM/REML `ccc_rm_reml()` and non-parametric
    `ccc_rm_ustat()`),
  - Intraclass correlation (wide-data `icc()` with pairwise and overall
    scope, repeated-measures REML `icc_rm_reml()`)
- Interactive Shiny viewers for matrix-style outputs with a dedicated
  repeated-measures correlation viewer (`view_rmcorr_shiny()`)

## Installation

``` r
# Install from CRAN
install.packages("matrixCorr")

# Development version from GitHub
# install.packages("remotes")
remotes::install_github("Prof-ThiagoOliveira/matrixCorr")
```

## Thread settings (`n_threads`)

Most computational functions expose an `n_threads` argument and default to:

``` r
getOption("matrixCorr.threads", 1L)
```

To define a package-wide default once per session:

``` r
options(matrixCorr.threads = parallel::detectCores(logical = FALSE))
```

## Quick start

### Wide-data correlation workflow

``` r
library(matrixCorr)

set.seed(1)
X <- as.data.frame(matrix(rnorm(300 * 6), ncol = 6))
names(X) <- paste0("V", 1:6)

R_pear <- pearson_corr(X, ci = TRUE)
R_bicor <- bicor(X)

print(R_pear, digits = 2)
#> Pearson correlation matrix
#>   method      : pearson
#>   dimensions  : 6 x 6
#>   ci          : yes
#> 
#>       V1    V2    V3    V4    V5    V6
#> V1  1.00  0.02  0.04 -0.02 -0.07  0.01
#> V2  0.02  1.00  0.04  0.03 -0.05  0.13
#> V3  0.04  0.04  1.00 -0.06  0.08 -0.14
#> V4 -0.02  0.03 -0.06  1.00  0.07  0.03
#> V5 -0.07 -0.05  0.08  0.07  1.00  0.04
#> V6  0.01  0.13 -0.14  0.03  0.04  1.00
summary(R_pear)
#> Pearson correlation summary
#>   method      : pearson
#>   dimensions  : 6 x 6
#>   pairs       : 15
#>   n_complete  : 300
#>   estimate    : -0.1410 to 0.1272
#>   most_negative: V3-V6 (-0.1410)
#>   most_positive: V2-V6 (0.1272)
#>   ci          : 95%
#>   ci_method   : fisher_z
#>   ci_width    : 0.222 to 0.226
#>   cross_zero  : 13 pair(s)
#> 
#> Strongest pairs by |estimate|
#> 
#>  item1 item2 estimate n_complete lwr    upr   
#>  V3    V6    -0.1410  300        -0.250 -0.028
#>  V2    V6     0.1272  300         0.014  0.237
#>  V3    V5     0.0776  300        -0.036  0.189
#>  V4    V5     0.0724  300        -0.041  0.184
#>  V1    V5    -0.0650  300        -0.177  0.049
#> ... 10 more rows not shown (omitted)
#> Use as.data.frame()/tidy()/as.matrix() to inspect the full result.
plot(R_bicor)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" alt="" width="100%" />

The same matrix-style workflow extends to Spearman, Kendall, distance
correlation, partial correlation, shrinkage correlation, latent
correlation, and the robust estimators `pbcor()`, `wincor()`, and
`skipped_corr()`.

### Agreement and repeated-measures workflow

``` r
set.seed(6)
S <- 24
Tm <- 4
id <- factor(rep(seq_len(S), each = 2 * Tm))
method <- factor(rep(rep(c("A", "B"), each = Tm), times = S))
time <- rep(rep(seq_len(Tm), times = 2), times = S)

u <- rnorm(S, 0, 0.9)[as.integer(id)]
um <- rnorm(S * 2, 0, 0.25)
um <- um[(as.integer(id) - 1L) * 2L + as.integer(method)]
y <- u + um + (method == "B") * 0.2 + rnorm(length(id), 0, 0.35)

dat_rm <- data.frame(y, id, method, time)

fit_ccc_rm <- ccc_rm_reml(
  dat_rm,
  response = "y",
  subject = "id",
  method = "method",
  time = "time"
)

summary(fit_ccc_rm)
#> 
#> Repeated-measures concordance (REML)
#> 
#> Concordance estimates
#> 
#>  item1 item2 estimate n_subjects n_obs SB     se_ccc
#>  A     B     0.8996   24         96    0.0554 0.0184
#> 
#> Variance components
#> 
#>  sigma2_subject sigma2_subject_method sigma2_subject_time sigma2_error
#>  0.7941         0                     0                   0.1329      
#> 
#> AR(1) diagnostics
#> 
#>  ar1_rho ar1_rho_lag1 ar1_rho_mom ar1_pairs ar1_pval use_ar1 ar1_recommend
#>  0.0179  0.0179       0.0179      144       0.8298   FALSE   FALSE
```

Agreement and reliability methods use the same general inspection
pattern, but they target different quantities. The package includes
Bland-Altman analysis, concordance correlation, and intraclass
correlation for both wide and repeated-measures designs.

## Vignettes

The package documentation is organised as a set of workflow vignettes.
The README is intentionally brief; the vignettes are the main user
guide.

Start here:

- `vignette("v01-matrixCorr-introduction", package = "matrixCorr")`
  introduces the package structure, common object behaviour, and shared
  display conventions.
- `vignette("v02-wide-correlation-workflows", package = "matrixCorr")`
  covers Pearson, Spearman, Kendall, distance correlation, and the
  general wide-data matrix workflow.
- `vignette("v03-robust-and-highdim-correlation", package = "matrixCorr")`
  covers robust estimators, shrinkage, and high-dimensional settings.
- `vignette("v04-latent-and-mixed-scale-correlation", package = "matrixCorr")`
  covers tetrachoric, polychoric, polyserial, and biserial correlation.
- `vignette("v05-agreement-and-icc-wide", package = "matrixCorr")`
  covers Bland-Altman analysis, concordance, and intraclass correlation
  for wide data.
- `vignette("v06-repeated-measures-workflows", package = "matrixCorr")`
  covers repeated-measures correlation, repeated agreement, and repeated
  reliability workflows.

If you want a compact overview of the available estimators, start with
the introduction vignette and then move to the workflow family that
matches your data layout and scientific question.

## Contributing

Issues and pull requests are welcome. Please see `CONTRIBUTING.md` for
guidelines and `cran-comments.md`/`DESCRIPTION` for package metadata.

## License

MIT [Thiago de Paula Oliveira](https://orcid.org/0000-0002-4555-2584)

See inst/LICENSE for the full MIT license text.
