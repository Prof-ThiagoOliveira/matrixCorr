## R CMD check results

0 errors | 0 warnings | 0 notes

## Comments

This is a new release of `matrixCorr`.
* updated cpp generic scripts to remove unnecessary overhead
* fixed S3 summary/print formatting for correlation summaries so the
  `Strongest pairs by |estimate|` section now consistently respects the
  user-requested `digits` (and `ci_digits` for CI bounds)
* added regression tests covering the `summary()`/`print()` digits behavior
* added `cohen_kappa()`, an exported implementation of pairwise unweighted
  Cohen's kappa for nominal ratings, supporting both two-vector and matrix
  workflows
* added `weighted_kappa()`, an exported implementation of pairwise weighted
  Cohen's kappa for ordered categorical ratings, supporting two-vector and
  matrix workflows with unweighted, linear, quadratic, and custom symmetric
  agreement weights
* added `multirater_kappa()`, an exported implementation of panel-level
  multi-rater nominal kappa for three or more raters, including Fleiss'
  fixed-marginal kappa and Randolph's free-marginal kappa, with optional
  exact Fleiss estimation and asymptotic or jackknife inference
* added `krippendorff_alpha()`, an exported implementation of
  Krippendorff's alpha for panel-level reliability/agreement with nominal,
  ordinal, interval, and ratio disagreement, supporting ratings and
  counts input, missing ratings, customary bootstrap CI, and analytical
  jackknife inference
* added `prob_agree()`, an exported implementation of the Stevens and
  Anderson-Cook (2017) probability of agreement for binomial reliability
  curves.

## Bug

* A dedicated GitHub Actions workflow has been added to test the MKL-related
  segmentation fault reported by CRAN. The workflow is available here:
  https://github.com/Prof-ThiagoOliveira/matrixCorr/blob/main/.github/workflows/mkl-check.yaml
  This workflow runs matrixCorr under R-devel on Ubuntu with Intel oneMKL
  configured as the BLAS/LAPACK backend. It runs the repeated-measures
  CCC/ICC REML reproducer and the repeated-measures CCC and ICC test files.
  The package passes these Ubuntu + MKL checks in both a single-threaded
  configuration (OMP_NUM_THREADS = 1, MKL_NUM_THREADS = 1) and a
  multi-threaded configuration (OMP_NUM_THREADS = 2, MKL_NUM_THREADS = 2),
  without any segmentation fault.
* addressed an MKL-specific segfault reported on the CRAN MKL checks for
  `ccc_rm_reml()`/`icc_rm_reml()`. The REML C++ backend now enforces the
  documented single-thread OpenMP default when `n_threads` is omitted, and the
  BLAS thread guard is also applied to Intel MKL when runtime controls are
  available. This avoids nested package OpenMP plus threaded MKL BLAS/LAPACK
  execution in the repeated-measures REML path.
* added defensive C++ checks for a second repeated-measures REML segfault path
  where optional method/time inputs could be represented as 0-length vectors in
  low-level calls. Subject-level grouping now keeps method/time arrays aligned
  with row indices using `-1` sentinels when those optional inputs are absent,
  and `ccc_vc_cpp()` now validates empty designs and optional matrix dimensions
  before creating Armadillo views. A regression test covers the absent
  method/time-vector path.

## Test environments

Local checks and tests were run on Windows with R 4.5.x. A source build and
targeted installed-package checks for the failing repeated-measures REML
examples completed successfully locally. A full Windows check compiled and
installed the package and completed examples successfully; the local run was
stopped while running the longer testthat suite.
