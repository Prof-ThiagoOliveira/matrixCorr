## R CMD check results

0 errors | 0 warnings | 0 notes

## Comments

This is a new release of `matrixCorr`.
* Adds optional inference for distance correlation.
* Adds optional large-sample confidence intervals for biweight
  mid-correlation.
* Adds optional large-sample p-values and confidence intervals for
  biserial correlation.
* Adds optional p-values and percentile-bootstrap confidence intervals
  for percentage bend correlation and pairwise Winsorized correlation.
* Adds intraclass correlation support for wide data via `icc()`,
  including pairwise and overall ANOVA-based coefficients.
* Adds repeated-measures intraclass correlation via `icc_rm_reml()`
  using the package's existing REML variance-components backend.
* Adds a vignette set covering the main package workflows, including
  wide-data correlation, robust and high-dimensional methods, latent
  correlation, agreement and intraclass correlation, and
  repeated-measures workflows.



## Test environments

Local checks and tests were run on Windows with R 4.5.x.
