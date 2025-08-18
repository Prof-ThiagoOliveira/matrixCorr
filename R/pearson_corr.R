#' @title Fast pairwise Pearson correlation
#'
#' @description
#' Computes all pairwise Pearson correlation coefficients for the numeric
#' columns of a matrix or data frame using a high-performance \code{C++}
#' backend.
#'
#' This function uses a direct Pearson formula implementation in C++ to achieve
#' fast and scalable correlation computations, especially for large datasets.
#'
#' @param data A numeric matrix or a data frame with at least two numeric
#' columns. All non-numeric columns will be excluded. Each column must have
#' at least two non-missing values and contain no NAs.
#'
#' @return A symmetric numeric matrix where the \code{(i, j)}-th element is
#' the Pearson correlation between the \code{i}-th and \code{j}-th
#' numeric columns of the input.
#'
#' @details
#' Pearson's correlation measures the linear relationship between two variables.
#' It is defined as:
#'
#' \deqn{ r = \frac{\sum (x_i - \bar{x})(y_i - \bar{y})}{\sqrt{\sum (x_i -
#' \bar{x})^2 \sum (y_i - \bar{y})^2}} }
#'
#' This implementation avoids calculating means explicitly and instead
#' uses a numerically stable online form.
#'
#' @note Missing values are not allowed. Columns with fewer than two
#' observations are excluded.
#'
#' @references
#' Pearson, K. (1895). "Notes on regression and inheritance in the case of
#' two parents". Proceedings of the Royal Society of London, 58, 240–242.
#'
#' @examples
#' mat <- cbind(a = rnorm(100), b = rnorm(100), c = rnorm(100))
#' pr <- pearson_corr(mat)
#' print(pr)
#' plot(pr)
#'
#' \dontrun{
#' df <- data.frame(x = rnorm(1e4), y = rnorm(1e4), z = rnorm(1e4))
#' pearson_corr(df)
#' }
#'
#' @useDynLib matrixCorr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @seealso \code{\link{print.pearson_corr}}, \code{\link{plot.pearson_corr}}
#' @author Thiago de Paula Oliveira \email{toliveira@abacusbio.com}
#' @export
pearson_corr <- function(data) {
  numeric_data <- validate_corr_input(data)
  colnames_data <- colnames(numeric_data)
  result <- pearson_matrix_cpp(numeric_data)
  colnames(result) <- rownames(result) <- colnames_data
  attr(result, "method") <- "pearson"
  attr(result, "description") <- "Pairwise Pearson correlation matrix"
  attr(result, "package") <- "matrixCorr"
  class(result) <- c("pearson_corr", "matrix")
  return(result)
}


#' @rdname pearson_corr
#' @method print pearson_corr
#' @title Print Method for \code{pearson_corr} Objects
#'
#' @description Prints a summary of the Pearson correlation matrix,
#' including description and method metadata.
#'
#' @param x An object of class \code{pearson_corr}.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @return Invisibly returns the \code{pearson_corr} object.
#' @author Thiago de Paula Oliveira
#' @export
print.pearson_corr <- function(x, ...) {
  cat("pearson_corr object\n")
  desc <- attr(x, "description")
  if (!is.null(desc)) cat("Description:", desc, "\n")
  method <- attr(x, "method")
  if (!is.null(method)) cat("Method:", method, "\n")
  cat("\nCorrelation matrix:\n")
  NextMethod("print", x, ...)
  invisible(x)
}

#' @rdname pearson_corr
#' @method plot pearson_corr
#' @title Plot Method for \code{pearson_corr} Objects
#'
#' @description Generates a ggplot2-based heatmap of the Pearson correlation
#' matrix.
#'
#' @param x An object of class \code{pearson_corr}.
#' @param title Plot title. Default is \code{"Pearson correlation heatmap"}.
#' @param low_color Color for the minimum correlation. Default is
#' \code{"indianred1"}.
#' @param high_color Color for the maximum correlation. Default is
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
plot.pearson_corr <-
  function(x, title = "Pearson correlation heatmap",
           low_color = "indianred1", high_color = "steelblue1",
           mid_color = "white", value_text_size = 4, ...) {

  if (!inherits(x, "pearson_corr")) {
    stop("x must be of class 'pearson_corr'.")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }

  mat <- as.matrix(x)
  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("Var1", "Var2", "Pearson")

  df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))

  p <- ggplot2::ggplot(df, ggplot2::aes(Var2, Var1, fill = Pearson)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Pearson)),
                       size = value_text_size, color = "black") +
    ggplot2::scale_fill_gradient2(
      low = low_color, high = high_color, mid = mid_color,
      midpoint = 0, limit = c(-1, 1), name = "Pearson"
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

