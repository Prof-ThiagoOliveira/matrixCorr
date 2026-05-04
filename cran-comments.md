## R CMD check results

0 errors | 0 warnings | 0 notes

## Comments

This is a new release of `matrixCorr`.
* updated cpp generic scripts to remove unnecessary overhead
* fixed S3 summary/print formatting for correlation summaries so the
  `Strongest pairs by |estimate|` section now consistently respects the
  user-requested `digits` (and `ci_digits` for CI bounds)
* added regression tests covering the `summary()`/`print()` digits behavior
* addressed an MKL-specific segfault reported on the CRAN MKL checks for
  `ccc_rm_reml()`/`icc_rm_reml()`. The REML C++ backend now enforces the
  documented single-thread OpenMP default when `n_threads` is omitted, and the
  BLAS thread guard is also applied to Intel MKL when runtime controls are
  available. This avoids nested package OpenMP plus threaded MKL BLAS/LAPACK
  execution in the repeated-measures REML path.



## Test environments

Local checks and tests were run on Windows with R 4.5.x. A source build and
targeted installed-package checks for the failing repeated-measures REML
examples completed successfully locally. A full Windows check compiled and
installed the package and completed examples successfully; the local run was
stopped while running the longer testthat suite.
