## R CMD check results

0 errors | 0 warnings | 0 notes

## Comments

This is a new release of `matrixCorr`.
* updated cpp generic scripts to remove unnecessary overhead
* fixed S3 summary/print formatting for correlation summaries so the
  `Strongest pairs by |estimate|` section now consistently respects the
  user-requested `digits` (and `ci_digits` for CI bounds)
* added regression tests covering the `summary()`/`print()` digits behavior



## Test environments

Local checks and tests were run on Windows with R 4.5.x.
