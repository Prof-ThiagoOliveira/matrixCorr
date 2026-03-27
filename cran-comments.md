## R CMD check results

0 errors | 0 warnings | 0 notes

## Comments

This is a new release of `matrixCorr`.

Main changes in this submission:

- Shortened several public function names for a more concise API:
  `ba()`, `ba_rm()`, `ccc_rm_reml()`, `ccc_rm_ustat()`, `pcorr()`,
  `bicor()`, and `dcor()` are now the primary exported names.
- Added `pcorr(method = "glasso")` for sparse graphical-lasso partial
  correlation estimation, alongside the existing sample, ridge, and OAS
  options.
- Added robust correlation methods `pbcor()`, `wincor()`, and
  `skipped_corr()`.
- Added latent correlation methods such as `tetrachoric()`, 
`polychoric()`, `polyserial()`, and `biserial()`.
- Added standard `summary()` methods across the matrix-style 
correlation objects.
- Added repeated-measures correlation via `rmcorr()` for within-subject
  association, with pairwise and matrix outputs.
- Updated `summary(ba_rm())` console printing to show repeated-measures
  Bland-Altman results in grouped blocks, keeping all reported outputs
  while improving readability in the console.
- Added a dedicated repeated-measures correlation Shiny viewer via
  `view_rmcorr_shiny()`.

## Test environments

Local checks and tests were run on Windows with R 4.5.x.
