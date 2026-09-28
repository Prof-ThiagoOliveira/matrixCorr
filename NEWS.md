# matrixCorr 0.12.3

## New features

* Distance-correlation API correction: `dcor()` now returns the conventional
  sample distance correlation \(R_n\). Previous versions returned a
  U-centred bias-corrected squared distance-correlation statistic. The 
  latter is now exposed as `bcdcor()`,
  with signed finite-sample values preserved. Use `dcor(..., squared = TRUE)`
  for conventional \(R_n^2\).
* Added `robust_ccc()` for pairwise MCD-based robust concordance correlation.
  The estimator substitutes joint reweighted minimum covariance determinant
  location and scatter estimates into Lin's concordance coefficient and
  exposes `alpha`, `mcd_nsamp`, and `seed` controls.
* `robust_ccc()` supports the package's standard matrix, sparse, and edge-list
  outputs; `error`, `complete`, and `pairwise` missing-data modes; and the
  shared `print()`, `summary()`, `plot()`, `estimate()`, `coef()`, `tidy()`,
  `ci()`, and `confint()` interfaces.
* Optional paired percentile-bootstrap confidence intervals are available via
  `ci = TRUE`, with bootstrap success counts retained in result diagnostics.
