#' @title Pairwise Winsorized correlation
#'
#' @description
#' Computes all pairwise Winsorized correlation coefficients for the numeric
#' columns of a matrix or data frame using a high-performance 'C++' backend.
#'
#' This function Winsorizes each margin at proportion \code{tr} and then
#' computes ordinary Pearson correlation on the Winsorized values. It is a
#' simple robust alternative to Pearson correlation when the main concern is
#' unusually large or small observations in the marginal distributions.
#'
#' @param data A numeric matrix or a data frame with at least two numeric
#' columns. All non-numeric columns will be excluded.
#' @param tr Winsorization proportion in \code{[0, 0.5)}. For a sample of size
#'   \eqn{n}, let \eqn{g = \lfloor tr \cdot n \rfloor}; the \eqn{g} smallest
#'   observations are set to the \eqn{(g+1)}-st order statistic and the
#'   \eqn{g} largest observations are set to the \eqn{(n-g)}-th order
#'   statistic. Default \code{0.2}.
#' @param na_method One of \code{"error"} (default) or \code{"pairwise"}.
#' @param n_threads Integer \eqn{\geq 1}. Number of OpenMP threads. Defaults to
#'   \code{getOption("matrixCorr.threads", 1L)}.
#' @param x An object of class \code{wincor}.
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
#' @return A symmetric correlation matrix with class \code{wincor} and
#'   attributes \code{method = "winsorized_correlation"}, \code{description},
#'   and \code{package = "matrixCorr"}.
#'
#' @details
#' Let \eqn{X \in \mathbb{R}^{n \times p}} be a numeric matrix with rows as
#' observations and columns as variables. For a column
#' \eqn{x = (x_i)_{i=1}^n}, write the order statistics as
#' \eqn{x_{(1)} \le \cdots \le x_{(n)}} and let
#' \eqn{g = \lfloor tr \cdot n \rfloor}. The Winsorized values can be written as
#' \deqn{
#' x_i^{(w)} \;=\; \max\!\bigl\{x_{(g+1)},\, \min(x_i, x_{(n-g)})\bigr\}.
#' }
#' For two columns \eqn{x} and \eqn{y}, the Winsorized correlation is the
#' ordinary Pearson correlation computed from \eqn{x^{(w)}} and \eqn{y^{(w)}}:
#' \deqn{
#' r_w(x,y) \;=\;
#' \frac{\sum_{i=1}^n (x_i^{(w)}-\bar x^{(w)})(y_i^{(w)}-\bar y^{(w)})}
#'      {\sqrt{\sum_{i=1}^n (x_i^{(w)}-\bar x^{(w)})^2}\;
#'       \sqrt{\sum_{i=1}^n (y_i^{(w)}-\bar y^{(w)})^2}}.
#' }
#'
#' In matrix form, let \eqn{X^{(w)}} contain the Winsorized columns and define
#' the centred, unit-norm columns
#' \deqn{
#' z_{\cdot j} =
#' \frac{x_{\cdot j}^{(w)} - \bar x_j^{(w)} \mathbf{1}}
#'      {\sqrt{\sum_{i=1}^n (x_{ij}^{(w)}-\bar x_j^{(w)})^2}},
#' \qquad j=1,\ldots,p.
#' }
#' If \eqn{Z = [z_{\cdot 1}, \ldots, z_{\cdot p}]}, then the Winsorized
#' correlation matrix is
#' \deqn{
#' R_w \;=\; Z^\top Z.
#' }
#'
#' Winsorization acts on each margin separately, so it guards against marginal
#' outliers and heavy tails but does not target unusual points in the joint
#' cloud. This implementation Winsorizes each column in 'C++', centres and
#' normalises it, and forms the complete-data matrix from cross-products. With
#' \code{na_method = "pairwise"}, each pair is recomputed on its overlap of
#' non-missing rows. As with Pearson correlation, the complete-data path yields
#' a symmetric positive semidefinite matrix, whereas pairwise deletion can
#' break positive semidefiniteness.
#'
#' \strong{Computational complexity.} In the complete-data path, Winsorizing the
#' columns requires sorting within each column, and forming the cross-product
#' matrix costs \eqn{O(n p^2)} with \eqn{O(p^2)} output storage.
#'
#' @references
#' Wilcox, R. R. (1993). Some results on a Winsorized correlation coefficient.
#' British Journal of Mathematical and Statistical Psychology, 46(2), 339-349.
#' \doi{10.1111/j.2044-8317.1993.tb01020.x}
#'
#' Wilcox, R. R. (2012). Introduction to Robust Estimation and Hypothesis
#' Testing (3rd ed.). Academic Press.
#'
#' @seealso [pbcor()], [skipped_corr()], [bicor()]
#'
#' @examples
#' set.seed(11)
#' X <- matrix(rnorm(180 * 4), ncol = 4)
#' X[sample(length(X), 6)] <- X[sample(length(X), 6)] - 12
#'
#' R <- wincor(X, tr = 0.2)
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
wincor <- function(data,
                   tr = 0.2,
                   na_method = c("error", "pairwise"),
                   n_threads = getOption("matrixCorr.threads", 1L)) {
  na_method <- match.arg(na_method)
  check_scalar_numeric(tr,
                       arg = "tr",
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
    wincor_matrix_cpp(numeric_data, tr = tr, n_threads = n_threads)
  } else {
    wincor_matrix_pairwise_cpp(numeric_data, tr = tr, min_n = 5L, n_threads = n_threads)
  }

  colnames(res) <- rownames(res) <- colnames_data
  .mc_structure_corr_matrix(
    res,
    class_name = "wincor",
    method = "winsorized_correlation",
    description = paste0(
      "Winsorized correlation; tr = ", tr,
      "; NA mode = ", na_method, "."
    )
  )
}

#' @rdname wincor
#' @method print wincor
#' @export
print.wincor <- function(x, digits = 4, n = NULL, topn = NULL,
                         max_vars = NULL, width = NULL,
                         show_ci = NULL, ...) {
  .mc_print_corr_matrix(
    x,
    header = "Winsorized correlation matrix",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
}

#' @rdname wincor
#' @method plot wincor
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
#' @export
plot.wincor <- function(x,
                        title = "Winsorized correlation heatmap",
                        low_color = "indianred1",
                        high_color = "steelblue1",
                        mid_color = "white",
                        value_text_size = 4,
                        show_value = TRUE,
                        ...) {
  .mc_plot_corr_matrix(
    x,
    class_name = "wincor",
    fill_name = "wincor",
    title = title,
    low_color = low_color,
    high_color = high_color,
    mid_color = mid_color,
    value_text_size = value_text_size,
    show_value = show_value,
    ...
  )
}

#' @rdname wincor
#' @method summary wincor
#' @param object An object of class \code{wincor}.
#' @export
summary.wincor <- function(object, n = NULL, topn = NULL,
                           max_vars = NULL, width = NULL,
                           show_ci = NULL, ...) {
  .mc_summary_corr_matrix(object, topn = topn)
}
