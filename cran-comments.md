## Resubmission after archive

This is a resubmission of matrixCorr, which was archived on 2025-10-19
due to check issues.

Changes since 0.5.1:
* Fixed a segfault in the `ccc_lmm_reml()` example that appeared on the
  MKL BLAS additional-issues checks.
* Simplified the `configure` script to avoid bash-specific `local`
  declarations and use only POSIX `/bin/sh` features.
