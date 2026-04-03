# Internal exact percentage-bend helpers. These stay in R until the compiled
# backend reproduces WRS2::pbcor() exactly across wide matrices.
.mc_pbcor_scores <- function(z, beta = 0.2) {
  med <- stats::median(z)
  omega <- sort(abs(z - med))[floor((1 - beta) * length(z))]
  if (!is.finite(omega) || omega <= 0) return(rep(NA_real_, length(z)))

  psi <- (z - med) / omega
  i1 <- sum(psi < -1)
  i2 <- sum(psi > 1)
  sx <- z
  sx[psi < -1 | psi > 1] <- 0
  denom <- length(z) - i1 - i2
  if (denom <= 0) return(rep(NA_real_, length(z)))

  theta <- (sum(sx) + omega * (i2 - i1)) / denom
  a <- (z - theta) / omega
  a[a <= -1] <- -1
  a[a >= 1] <- 1
  a
}

.mc_pbcor_pair_exact <- function(x, y, beta = 0.2) {
  ax <- .mc_pbcor_scores(x, beta = beta)
  ay <- .mc_pbcor_scores(y, beta = beta)
  if (anyNA(ax) || anyNA(ay)) return(NA_real_)

  den_x <- sum(ax^2)
  den_y <- sum(ay^2)
  if (!is.finite(den_x) || !is.finite(den_y) || den_x <= 0 || den_y <= 0) {
    return(NA_real_)
  }

  sum(ax * ay) / sqrt(den_x * den_y)
}

.mc_pbcor_matrix_exact <- function(X, beta = 0.2) {
  n <- nrow(X)
  p <- ncol(X)
  if (p < 2L) stop("At least two numeric columns are required.", call. = FALSE)
  scores <- matrix(0, nrow = n, ncol = p)
  valid <- rep(FALSE, p)

  for (j in seq_len(p)) {
    a <- .mc_pbcor_scores(X[, j], beta = beta)
    if (!anyNA(a)) {
      ss <- sum(a^2)
      if (is.finite(ss) && ss > 0) {
        scores[, j] <- a
        valid[j] <- TRUE
      }
    }
  }

  ss <- colSums(scores^2)
  out <- as.matrix(crossprod(scores))
  denom <- sqrt(outer(ss, ss))
  out <- as.matrix(out / denom)
  out <- pmin(1, pmax(-1, out))
  out <- matrix(out, nrow = p, ncol = p, dimnames = list(colnames(X), colnames(X)))

  if (any(!valid)) {
    bad <- which(!valid)
    out[bad, ] <- NA_real_
    out[, bad] <- NA_real_
  }
  diag(out) <- 1
  out
}

.mc_pbcor_matrix_pairwise_exact <- function(X, beta = 0.2, min_n = 5L) {
  p <- ncol(X)
  if (p < 2L) stop("At least two numeric columns are required.", call. = FALSE)
  out <- matrix(NA_real_, nrow = p, ncol = p)
  dimnames(out) <- list(colnames(X), colnames(X))

  finite_n <- colSums(is.finite(X))
  diag(out) <- 1

  for (j in seq_len(p - 1L)) {
    for (k in seq.int(j + 1L, p)) {
      ok <- is.finite(X[, j]) & is.finite(X[, k])
      if (sum(ok) < min_n) next
      out[j, k] <- .mc_pbcor_pair_exact(X[ok, j], X[ok, k], beta = beta)
      out[k, j] <- out[j, k]
    }
  }

  if (any(finite_n < 2L)) {
    bad <- which(finite_n < 2L)
    out[bad, ] <- NA_real_
    out[, bad] <- NA_real_
    diag(out) <- 1
  }
  out
}

#' Percentage bend correlation
#'
#' @description
#' Computes all pairwise percentage bend correlations for the numeric columns of
#' a matrix or data frame. Percentage bend correlation limits the influence of
#' extreme marginal observations by bending standardised deviations into the
#' interval \eqn{[-1, 1]}, yielding a Pearson-like measure that is robust to
#' outliers and heavy tails.
#'
#' @param data A numeric matrix or data frame containing numeric columns.
#' @param beta Bending constant in \code{[0, 0.5)} that sets the cutoff used to
#'   bend standardised deviations toward the interval \eqn{[-1, 1]}. Larger
#'   values cause more observations to be bent and increase resistance to
#'   marginal outliers. Default \code{0.2}. See Details.
#' @param na_method One of \code{"error"} (default) or \code{"pairwise"}.
#'   With \code{"pairwise"}, each correlation is computed on the overlapping
#'   complete rows for the column pair.
#' @param n_threads Integer \eqn{\geq 1}. Kept for API consistency with the
#'   other robust correlation wrappers. It is currently validated but not used
#'   by the exact percentage-bend implementation.
#' @param x An object of class \code{pbcor}.
#' @param digits Integer; number of digits to print.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param ... Additional arguments passed to the underlying print or plot helper.
#' @param title Character; plot title.
#' @param low_color,high_color,mid_color Colors used in the heatmap.
#' @param value_text_size Numeric text size for overlaid cell values.
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
#'
#' @return A symmetric correlation matrix with class \code{pbcor} and
#'   attributes \code{method = "percentage_bend_correlation"},
#'   \code{description}, and \code{package = "matrixCorr"}.
#'
#' @details
#' For a column \eqn{x = (x_i)_{i=1}^n}, let \eqn{m = \mathrm{med}(x)} and let
#' \eqn{\omega_\beta(x)} be the \eqn{\lfloor (1-\beta)n \rfloor}-th order
#' statistic of \eqn{|x_i - m|}. The constant \code{beta} determines the cutoff
#' \eqn{\omega_\beta(x)} used to standardise deviations from the median. As
#' \code{beta} increases, the selected cutoff becomes smaller, so a larger
#' fraction of observations is truncated to the bounds \eqn{-1} and \eqn{1}.
#' This makes the correlation more resistant to marginal outliers. The one-step
#' percentage-bend location is
#' \deqn{
#' \hat\theta_{pb}(x) =
#' \frac{\sum_{i: |\psi_i| \le 1} x_i + \omega_\beta(x)(i_2 - i_1)}
#'      {n - i_1 - i_2},
#' \qquad
#' \psi_i = \frac{x_i - m}{\omega_\beta(x)},
#' }
#' where \eqn{i_1 = \sum \mathbf{1}(\psi_i < -1)} and
#' \eqn{i_2 = \sum \mathbf{1}(\psi_i > 1)}.
#'
#' The bent scores are
#' \deqn{
#' a_i = \max\{-1, \min(1, (x_i - \hat\theta_{pb}(x))/\omega_\beta(x))\},
#' }
#' and likewise \eqn{b_i} for a second variable \eqn{y}. The percentage bend
#' correlation is
#' \deqn{
#' r_{pb}(x,y) =
#' \frac{\sum_i a_i b_i}
#'      {\sqrt{\sum_i a_i^2}\sqrt{\sum_i b_i^2}}.
#' }
#'
#' When \code{na_method = "error"}, bent scores are computed once per column and
#' the matrix is formed from their cross-products. When
#' \code{na_method = "pairwise"}, each pair is recomputed on its complete-case
#' overlap, which can break positive semidefiniteness as with pairwise Pearson
#' correlation.
#'
#' @references
#' Wilcox, R. R. (1994). The percentage bend correlation coefficient.
#' Psychometrika, 59(4), 601-616. \doi{10.1007/BF02294395}
#'
#' @seealso [wincor()], [skipped_corr()], [bicor()]
#'
#' @examples
#' set.seed(10)
#' X <- matrix(rnorm(150 * 4), ncol = 4)
#' X[sample(length(X), 8)] <- X[sample(length(X), 8)] + 10
#'
#' R <- pbcor(X)
#' print(R, digits = 2)
#' summary(R)
#' plot(R)
#'
#' # Interactive viewing (requires shiny)
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   view_corr_shiny(R)
#' }
#'
#' @author Thiago de Paula Oliveira
#' @export
pbcor <- function(data,
                  beta = 0.2,
                  na_method = c("error", "pairwise"),
                  n_threads = getOption("matrixCorr.threads", 1L)) {
  na_method <- match.arg(na_method)
  check_scalar_numeric(beta,
                       arg = "beta",
                       lower = 0,
                       upper = 0.5,
                       closed_lower = TRUE,
                       closed_upper = FALSE)
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")

  numeric_data <- if (na_method == "error") {
    validate_corr_input(data)
  } else {
    validate_corr_input(data, check_na = FALSE)
  }
  colnames_data <- colnames(numeric_data)

  res <- if (na_method == "error") {
    .mc_pbcor_matrix_exact(numeric_data, beta = beta)
  } else {
    .mc_pbcor_matrix_pairwise_exact(numeric_data, beta = beta, min_n = 5L)
  }

  colnames(res) <- rownames(res) <- colnames_data
  .mc_structure_corr_matrix(
    res,
    class_name = "pbcor",
    method = "percentage_bend_correlation",
    description = paste0(
      "Percentage bend correlation; beta = ", beta,
      "; NA mode = ", na_method, "."
    )
  )
}

#' @rdname pbcor
#' @method print pbcor
#' @export
print.pbcor <- function(x, digits = 4, n = NULL, topn = NULL,
                        max_vars = NULL, width = NULL,
                        show_ci = NULL, ...) {
  .mc_print_corr_matrix(
    x,
    header = "Percentage bend correlation matrix",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
}

#' @rdname pbcor
#' @method plot pbcor
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
#' @export
plot.pbcor <- function(x,
                       title = "Percentage bend correlation heatmap",
                       low_color = "indianred1",
                       high_color = "steelblue1",
                       mid_color = "white",
                       value_text_size = 4,
                       show_value = TRUE,
                       ...) {
  .mc_plot_corr_matrix(
    x,
    class_name = "pbcor",
    fill_name = "pbcor",
    title = title,
    low_color = low_color,
    high_color = high_color,
    mid_color = mid_color,
    value_text_size = value_text_size,
    show_value = show_value,
    ...
  )
}

#' @rdname pbcor
#' @method summary pbcor
#' @param object An object of class \code{pbcor}.
#' @export
summary.pbcor <- function(object, n = NULL, topn = NULL,
                          max_vars = NULL, width = NULL,
                          show_ci = NULL, ...) {
  .mc_summary_corr_matrix(object, topn = topn)
}
