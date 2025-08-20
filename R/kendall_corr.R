#' @title Fast Pairwise Kendall's Tau Rank Correlation
#'
#' @description
#' Computes all pairwise Kendall's tau rank correlation coefficients for the
#' numeric columns of a matrix or data frame using a high-performance
#' \code{C++}.
#'
#' This function uses a fast and scalable algorithm implemented in \code{C++}
#' to compute both Kendall's tau-a (when no ties are present) and tau-b
#' (when ties are detected), making it suitable for large datasets. Internally,
#' it calls a highly optimized function that uses a combination of merge-sort-
#' based inversion counting and a Fenwick Tree (binary indexed tree)
#' for efficient tie handling.
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
#' For a dataset with \eqn{n} observations of two variables \eqn{X} and
#' \eqn{Y}, Kendall's tau coefficient is defined as:
#'
#' \deqn{ \tau = \frac{C - D}{\sqrt{(C + D + T_x)(C + D + T_y)}} }
#'
#' where:
#' \itemize{
#'   \item \eqn{C} is the number of concordant pairs defined by
#'   \eqn{(x_i - x_j)(y_i - y_j) > 0}
#'   \item \eqn{D} is the number of discordant pairs defined by
#'   \eqn{(x_i - x_j)(y_i - y_j) < 0}
#'   \item \eqn{T_x}, \eqn{T_y} are the number of tied pairs in \eqn{X} and
#'   \eqn{Y}, respectively
#' }
#'
#' When there are no ties, the function computes the faster \emph{tau-a}
#' version:
#' \deqn{ \tau_a = \frac{C - D}{n(n-1)/2} }
#'
#' The function automatically selects \emph{tau-a} or \emph{tau-b} depending
#' on the presence of ties. Performance is maximized by computing correlations
#' in \code{C++} directly from the matrix columns.
#'
#' @note Missing values are not allowed. Columns with fewer than two
#' observations are excluded.
#'
#' @references
#' Kendall, M. G. (1938). A New Measure of Rank Correlation. Biometrika,
#' 30(1/2), 81–93.
#'
#' @examples
#' # Basic usage with a matrix
#' mat <- cbind(a = rnorm(100), b = rnorm(100), c = rnorm(100))
#' kt <- kendall_tau(mat)
#' print(kt)
#' plot(kt)
#'
#' \dontrun{
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
#' }
#'
#' @seealso \code{\link{print.kendall_matrix}},
#' \code{\link{print.kendall_matrix}}
#' @author Thiago de Paula Oliveira \email{toliveira@abacusbio.com}
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
#' @param digits Integer; number of decimal places to print in the concordance
#' @param max_rows Optional integer; maximum number of rows to display.
#'  If \code{NULL}, all rows are shown.
#' @param max_cols Optional integer; maximum number of columns to display.
#' If \code{NULL}, all columns are shown.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @return Invisibly returns the \code{kendall_matrix} object.
#' @author Thiago de Paula Oliveira
#' @export
print.kendall_matrix <- function(x, digits = 3, max_rows = NULL,
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

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
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
