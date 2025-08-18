#' @title Fast pairwise Spearman's rank correlation
#'
#' @description
#' Computes all pairwise Spearman's rank correlation coefficients for the
#' numeric columns of a matrix or data frame using a high-performance
#' \code{C++} backend.
#'
#' This function ranks the data and computes Pearson correlation on ranks,
#' which is equivalent to Spearman’s rho. It supports large datasets and
#' is optimized in C++ for performance.
#'
#' @param data A numeric matrix or a data frame with at least two numeric
#' columns. All non-numeric columns will be excluded. Each column must have
#' at least two non-missing values and contain no NAs.
#'
#' @return A symmetric numeric matrix where the \code{(i, j)}-th element is
#' the Spearman correlation between the \code{i}-th and \code{j}-th
#' numeric columns of the input.
#'
#' @details
#' Spearman's correlation assesses how well the relationship between two
#' variables can be described using a monotonic function. It is computed by
#' ranking the data and then applying the Pearson correlation formula to the
#' ranks.
#'
#' @note Missing values are not allowed. Columns with fewer than two
#' observations are excluded.
#'
#' @references
#' Spearman, C. (1904). The proof and measurement of association between
#' two things. International Journal of Epidemiology, 39(5), 1137-1150.
#'
#' @examples
#' mat <- cbind(a = rnorm(100), b = rnorm(100), c = rnorm(100))
#' sp <- spearman_rho(mat)
#' print(sp)
#' plot(sp)
#'
#' \dontrun{
#' df <- data.frame(x = rnorm(1e4), y = rnorm(1e4), z = rnorm(1e4))
#' spearman_rho(df)
#'
#' tied_df <- data.frame(
#'   v1 = rep(1:5, each = 20),
#'   v2 = rep(5:1, each = 20),
#'   v3 = rnorm(100)
#' )
#' sp <- spearman_rho(tied_df)
#' print(sp)
#' plot(sp)
#' }
#'
#' @useDynLib matrixCorr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @seealso \code{\link{print.spearman_rho}}, \code{\link{plot.spearman_rho}}
#' @author Thiago de Paula Oliveira \email{toliveira@abacusbio.com}
#' @export
spearman_rho <- function(data) {
  numeric_data <- validate_corr_input(data)
  colnames_data <- colnames(numeric_data)
  result <- spearman_matrix_cpp(numeric_data)
  colnames(result) <- rownames(result) <- colnames_data
  attr(result, "method") <- "spearman"
  attr(result, "description") <- "Pairwise Spearman's rank correlation matrix"
  attr(result, "package") <- "matrixCorr"
  class(result) <- c("spearman_rho", "matrix")
  return(result)
}

#' @rdname spearman_rho
#' @method print spearman_rho
#' @title Print Method for \code{spearman_rho} Objects
#'
#' @description Prints a summary of the Spearman's correlation matrix,
#' including description and method metadata.
#'
#' @param x An object of class \code{spearman_rho}.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @return Invisibly returns the \code{spearman_rho} object.
#' @author Thiago de Paula Oliveira
#' @export
print.spearman_rho <- function(x, ...) {
  cat("spearman_rho object\n")
  desc <- attr(x, "description")
  if (!is.null(desc)) cat("Description:", desc, "\n")
  method <- attr(x, "method")
  if (!is.null(method)) cat("Method:", method, "\n")
  cat("\nCorrelation matrix:\n")
  NextMethod("print", x, ...)
  invisible(x)
}


#' @rdname spearman_rho
#' @method plot spearman_rho
#' @title Plot Method for \code{spearman_rho} Objects
#'
#' @description Generates a ggplot2-based heatmap of the Spearman's rank
#' correlation matrix.
#'
#' @param x An object of class \code{spearman_rho}.
#' @param title Plot title. Default is \code{"Spearman's rank correlation
#' heatmap"}.
#' @param low_color Color for the minimum rho value. Default is
#'  \code{"indianred1"}.
#' @param high_color Color for the maximum rho value. Default is
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
plot.spearman_rho <-
  function(x, title = "Spearman's rank correlation heatmap",
           low_color = "indianred1", high_color = "steelblue1",
           mid_color = "white", value_text_size = 4, ...) {

  if (!inherits(x, "spearman_rho")) {
    stop("x must be of class 'spearman_rho'.")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }

  mat <- as.matrix(x)
  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("Var1", "Var2", "Rho")

  df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))

  p <- ggplot2::ggplot(df, ggplot2::aes(Var2, Var1, fill = Rho)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Rho)),
                       size = value_text_size, color = "black") +
    ggplot2::scale_fill_gradient2(
      low = low_color, high = high_color, mid = mid_color,
      midpoint = 0, limit = c(-1, 1), name = "Rho"
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
