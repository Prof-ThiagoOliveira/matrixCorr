## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* Added a runtime guard so the BLAS thread locker is skipped on MKL builds (previous CRAN MKL checks hit a segfault inside the guard before any examples ran). Set `MATRIXCORR_DISABLE_BLAS_GUARD=1` to opt out explicitly.
