# Package-level documentation and namespace setup

#' matrixCorr
#'
#' Correlation, agreement, and intraclass correlation estimators with a
#' consistent S3 interface.
#'
#' `print()` methods provide compact previews and `summary()` methods provide
#' richer but still bounded digests. Truncation affects console display only;
#' full results remain available through direct extraction and coercion helpers
#' such as `as.matrix()`, `as.data.frame()`, and `tidy()` where implemented.
#'
#' Package-wide display defaults can be controlled with options:
#' \itemize{
#'   \item `matrixCorr.threads` (default `1L`; used as the default for
#'   `n_threads` in thread-aware estimators)
#'   \item `matrixCorr.print_max_rows` (default `20L`)
#'   \item `matrixCorr.print_topn` (default `5L`)
#'   \item `matrixCorr.print_max_vars` (default `NULL`, derived from console width)
#'   \item `matrixCorr.print_show_ci` (default `"yes"`)
#'   \item `matrixCorr.summary_max_rows` (default `12L`)
#'   \item `matrixCorr.summary_topn` (default `5L`)
#'   \item `matrixCorr.summary_max_vars` (default `10L`)
#'   \item `matrixCorr.summary_show_ci` (default `"yes"`)
#' }
#'
#' Display methods validate user-facing arguments with `cli` conditions.
#' Confidence-interval visibility is explicit: `show_ci` accepts only `"yes"`
#' or `"no"`.
#'
#' A common setup is:
#' \preformatted{
#' options(matrixCorr.threads = parallel::detectCores(logical = FALSE))
#' }
#'
#' @name matrixCorr-package
#' @aliases matrixCorr RcppExports
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib matrixCorr, .registration = TRUE
## usethis namespace: end

# Silence NOTES about NSE vars in ggplot2 etc.
utils::globalVariables(c(
  "Var1", "Var2", "Tau", "Rho", "Pearson", "CCC",
  "label", "PCor", "r", "dCor", "bicor", "ci_label",
  "Tetrachoric", "Polychoric", "Polyserial", "Biserial",
  "diffs", "hi_l", "hi_u", "lab", "lo_l", "lo_u", "loaL",
  "loaU", "md", "md_l", "md_u", "means", "j", "k"
))
