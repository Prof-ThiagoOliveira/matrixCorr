#' @title Pairwise Kendall's Tau Rank Correlation
#'
#' @description
#' Computes all pairwise Kendall's tau rank correlation coefficients for the
#' numeric columns of a matrix or data frame using a high-performance
#' 'C++'.
#'
#' This function uses a fast and scalable algorithm implemented in 'C++'
#' to compute both Kendall's tau-a (when no ties are present) and tau-b
#' (when ties are detected), making it suitable for large datasets. It follows
#' the Knight (1966) \eqn{O(n \log n)} scheme, where a single sort on one
#' variable, in-block sorting of the paired variable within tie groups, and a
#' global merge-sort–based inversion count with closed-form tie corrections.
#'
#' @param data A numeric matrix or a data frame with at least two numeric
#' columns. All non-numeric columns will be excluded. Each column must have
#' at least two non-missing values and contain no NAs.
#'
#' @return A symmetric numeric matrix where the \code{(i, j)}-th element is
#' the Kendall's tau correlation between the \code{i}-th and \code{j}-th
#' numeric columns of the input.
#'
#' @details
#' Kendall's tau is a rank-based measure of association between two variables.
#' For a dataset with \eqn{n} observations on variables \eqn{X} and \eqn{Y},
#' let \eqn{n_0 = n(n - 1)/2} be the number of unordered pairs, \eqn{C} the
#' number of concordant pairs, and \eqn{D} the number of discordant pairs.
#' Let \eqn{T_x = \sum_g t_g (t_g - 1)/2} and \eqn{T_y = \sum_h u_h (u_h - 1)/2}
#' be the numbers of tied pairs within \eqn{X} and within \eqn{Y}, respectively,
#' where \eqn{t_g} and \eqn{u_h} are tie-group sizes in \eqn{X} and \eqn{Y}.
#'
#' The tie-robust Kendall's tau-b is defined as:
#' \deqn{ \tau_b = \frac{C - D}{\sqrt{(n_0 - T_x)\,(n_0 - T_y)}}. }
#'
#' When there are no ties, the function computes the faster \emph{tau-a}
#' version:
#' \deqn{ \tau_a = \frac{C - D}{n(n-1)/2}. }
#'
#' Note that if \eqn{T_{xy}} denotes the number of pairs tied on both
#' variables, then pairs partition as \eqn{C + D + T_x + T_y + T_{xy} = n_0}.
#'
#' The function automatically selects \emph{tau-a} or \emph{tau-b} depending on
#' the presence of ties. In degenerate cases where a variable is constant
#' (so \eqn{n_0 = T_x} or \eqn{n_0 = T_y}), the tau-b denominator is zero and
#' the correlation is undefined (returned as \code{NA}).
#'
#' Performance is maximised by computing correlations in 'C++' directly from
#' the matrix columns using the Knight \eqn{O(n \log n)} procedure; columns are
#' discretised once and pairwise computations are parallelised where available.
#'
#' @note Missing values are not allowed. Columns with fewer than two
#' observations are excluded.
#'
#' @references
#' Kendall, M. G. (1938). A New Measure of Rank Correlation. \emph{Biometrika},
#' 30(1/2), 81–93.
#'
#' Knight, W. R. (1966). A Computer Method for Calculating Kendall’s Tau with
#' Ungrouped Data. \emph{Journal of the American Statistical Association},
#' 61(314), 436–439.
#'
#' @examples
#' # Basic usage with a matrix
#' mat <- cbind(a = rnorm(100), b = rnorm(100), c = rnorm(100))
#' kt <- kendall_tau(mat)
#' print(kt)
#' plot(kt)
#'
#' # With a large data frame
#' df <- data.frame(x = rnorm(1e4), y = rnorm(1e4), z = rnorm(1e4))
#' kendall_tau(df)
#'
#' # Including ties
#' tied_df <- data.frame(
#'   v1 = rep(1:5, each = 20),
#'   v2 = rep(5:1, each = 20),
#'   v3 = rnorm(100)
#' )
#' kt <- kendall_tau(tied_df)
#' print(kt)
#' plot(kt)
#'
#' @seealso \code{\link{print.kendall_matrix}},
#' \code{\link{print.kendall_matrix}}
#' @author Thiago de Paula Oliveira
#' @export
kendall_tau <- function(data) {
  numeric_data <- validate_corr_input(data)
  colnames_data <- colnames(numeric_data)
  result <- kendall_matrix_cpp(numeric_data)
  colnames(result) <- rownames(result) <- colnames_data
  attr(result, "method") <- "kendall"
  attr(result, "description") <-
    "Pairwise Kendall's tau (auto tau-a/tau-b) correlation matrix"
  attr(result, "package") <- "matrixCorr"
  class(result) <- c("kendall_matrix", "matrix")
  return(result)
}


#' @rdname kendall_tau
#' @method print kendall_matrix
#' @title Print Method for \code{kendall_matrix} Objects
#'
#' @description Prints a summary of the Kendall's tau correlation matrix,
#' including description and method metadata.
#'
#' @param x An object of class \code{kendall_matrix}.
#' @param digits Integer; number of decimal places to print
#' @param max_rows Optional integer; maximum number of rows to display.
#'  If \code{NULL}, all rows are shown.
#' @param max_cols Optional integer; maximum number of columns to display.
#' If \code{NULL}, all columns are shown.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @return Invisibly returns the \code{kendall_matrix} object.
#' @author Thiago de Paula Oliveira
#' @export
print.kendall_matrix <- function(x, digits = 4, max_rows = NULL,
                                 max_cols = NULL, ...) {
  method <- attr(x, "method")
  header <- if (!is.null(method)) {
    sprintf("Kendall correlation (%s):", method)
  } else {
    "Kendall correlation:"
  }
  cat(header, "\n")

  m <- as.matrix(x)
  attributes(m) <- attributes(m)[c("dim", "dimnames")]

  # Truncation for large matrices
  if (!is.null(max_rows) || !is.null(max_cols)) {
    nr <- nrow(m); nc <- ncol(m)
    r  <- if (is.null(max_rows)) nr else min(nr, max_rows)
    c  <- if (is.null(max_cols)) nc else min(nc, max_cols)
    mm <- round(m[seq_len(r), seq_len(c), drop = FALSE], digits)
    print(mm, ...)
    if (nr > r || nc > c) {
      cat(sprintf("... omitted: %d rows, %d cols\n", nr - r, nc - c))
    }
  } else {
    print(round(m, digits), ...)
  }

  invisible(x)
}

#' @rdname kendall_tau
#' @method plot kendall_matrix
#' @title Plot Method for \code{kendall_matrix} Objects
#'
#' @description Generates a ggplot2-based heatmap of the Kendall's tau
#' correlation matrix.
#'
#' @param x An object of class \code{kendall_matrix}.
#' @param title Plot title. Default is \code{"Kendall's Tau Correlation
#' Heatmap"}.
#' @param low_color Color for the minimum tau value. Default is
#' \code{"indianred1"}.
#' @param high_color Color for the maximum tau value. Default is
#' \code{"steelblue1"}.
#' @param mid_color Color for zero correlation. Default is \code{"white"}.
#' @param value_text_size Font size for displaying correlation values. Default
#' is \code{4}.
#' @param ... Additional arguments passed to \code{ggplot2::theme()} or other
#' \code{ggplot2} layers.
#'
#' @return A \code{ggplot} object representing the heatmap.
#' @import ggplot2
#' @author Thiago de Paula Oliveira
#' @export
plot.kendall_matrix <- function(x, title = "Kendall's Tau correlation heatmap",
                                low_color = "indianred1", high_color = "steelblue1",
                                mid_color = "white",
                                value_text_size = 4, ...) {

  if (!inherits(x, "kendall_matrix")) {
    stop("x must be of class 'kendall_matrix'.")
  }

  mat <- as.matrix(x)
  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("Var1", "Var2", "Tau")

  # Reverse Y-axis so the diagonal appears from top-left to bottom-right
  df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))

  p <- ggplot2::ggplot(df, ggplot2::aes(Var2, Var1, fill = Tau)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Tau)),
                       size = value_text_size, color = "black") +
    ggplot2::scale_fill_gradient2(
      low = low_color, high = high_color, mid = mid_color,
      midpoint = 0, limit = c(-1, 1), name = "Tau"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      ...
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = NULL, y = NULL)

  return(p)
}
