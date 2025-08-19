#' @title Pairwise Lin's concordance correlation coefficient
#'
#' @description
#' Computes all pairwise Lin's Concordance Correlation Coefficients (CCC)
#' from the numeric columns of a matrix or data frame. CCC measures both
#' precision (Pearson correlation) and accuracy (closeness to the 45-degree line).
#' This function is backed by a high-performance \code{C++} implementation.
#'
#' @details
#' Lin's CCC is defined as:
#' \deqn{
#' \rho_c = \frac{2 \cdot \mathrm{cov}(X, Y)}{\sigma_X^2 + \sigma_Y^2 +
#' (\mu_X - \mu_Y)^2}
#' }{
#' rho_c = 2 * cov(X, Y) / [var(X) + var(Y) + (mean(X) - mean(Y))^2]
#' }
#'
#' This formula combines the Pearson correlation coefficient
#' \eqn{r = \mathrm{cov}(X, Y) / (\sigma_X \sigma_Y)}
#' with a bias correction factor:
#' \deqn{
#' C_b = \frac{2 \sigma_X \sigma_Y}{\sigma_X^2 + \sigma_Y^2 + (\mu_X - \mu_Y)^2}
#' }
#'
#' Confidence intervals are not provided in the matrix version to retain speed,
#' but can be computed separately for individual variable pairs using the
#' scalar form of Lin's CCC.
#'
#' Missing values are not allowed; input must be clean numeric data.
#'
#' @param data A numeric matrix or data frame with at least two numeric columns.
#' Non-numeric columns will be ignored.
#' @param ci Logical; if TRUE, return lower and upper confidence bounds
#' @param conf_level Confidence level for CI, default = 0.95
#' @param verbose Logical; if TRUE, prints how many threads are used
#'
#' @return A symmetric numeric matrix with class \code{"ccc"} and attributes:
#' \itemize{
#'   \item \code{method}: The method used ("Lin's concordance")
#'   \item \code{description}: Description string
#' }
#'  If \code{ci = FALSE}, returns matrix of class \code{"ccc"}.
#'         If \code{ci = TRUE}, returns a list with elements: \code{est},
#'         \code{lwr.ci}, \code{upr.ci}.
#'
#' @seealso \code{\link{print.ccc}}, \code{\link{plot.ccc}}
#'
#' @examples
#' mat <- cbind(x = rnorm(100), y = rnorm(100), z = rnorm(100))
#' result <- ccc(mat)
#' print(result)
#' plot(result)
#'
#' @importFrom stats var cov cor
#' @importFrom graphics plot
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2
#' @importFrom ggplot2 theme_minimal element_text coord_fixed labs theme
#' @author Thiago de Paula Oliveira
#' @references
#' Lin L (1989). A concordance correlation coefficient to evaluate
#' reproducibility. Biometrics 45: 255-268.
#' @references
#' Lin L (2000). A note on the concordance correlation coefficient.
#' Biometrics 56: 324-325.
#' @references
#' Bland J, Altman D (1986). Statistical methods for assessing agreement
#' between two methods of clinical measurement. The Lancet 327: 307-310.
#' @export
ccc <- function(data, ci = FALSE, conf_level = 0.95, verbose = FALSE) {
  numeric_data <- validate_corr_input(data)
  mat <- as.matrix(numeric_data)
  colnames_data <- colnames(numeric_data)

  if (verbose) {
    cat("Using", openmp_threads(), "OpenMP threads\n")
  }

  if (ci) {
    out <- ccc_with_ci_cpp(mat, conf_level)
    out$est <- `dimnames<-`(out$est, list(colnames_data, colnames_data))
    out$lwr.ci <- `dimnames<-`(out$lwr.ci, list(colnames_data, colnames_data))
    out$upr.ci <- `dimnames<-`(out$upr.ci, list(colnames_data, colnames_data))
  } else {
    est <- ccc_cpp(mat)
    out <- `dimnames<-`(est, list(colnames_data, colnames_data))
  }

  attr(out, "method") <- "Lin's concordance"
  attr(out, "description") <- if (ci) {
    "Pairwise Lin's concordance with confidence intervals"
  } else {
    "Pairwise Lin's concordance correlation matrix"
  }
  attr(out, "package") <- "matrixCorr"
  class(out) <- "ccc"

  return(out)
}

#' @rdname ccc
#' @method print ccc
#' @param digits Integer; number of decimal places to print in the concordance
#' matrix (default is 3).
#' @export
print.ccc <- function(x, digits = 3, ...) {
  cat("ccc object\n")
  desc <- attr(x, "description")
  if (!is.null(desc)) cat("Description:", desc, "\n")
  method <- attr(x, "method")
  if (!is.null(method)) cat("Method:", method, "\n\n")

  if (is.list(x) && all(c("est", "lwr.ci", "upr.ci") %in% names(x))) {
    cat("Estimates:\n")
    print(round(x$est, digits = digits))
    cat("\n95% Confidence Intervals:\nLower:\n")
    print(round(x$lwr.ci, digits = 2))
    cat("\nUpper:\n")
    print(round(x$upr.ci, digits = 2))

  } else if (inherits(x, "matrix")) {
    cat("Concordance matrix:\n")
    print(round(x, digits))
  } else {
    stop("Invalid object format for class 'ccc'.")
  }

  invisible(x)
}


#' @rdname ccc
#' @method plot ccc
#' @param x An object of class \code{"ccc"} (either a matrix or a list with
#' confidence intervals).
#' @param digits Number of digits to display (for printing).
#' @param title Title for the plot.
#' @param low_color Color for low CCC values.
#' @param high_color Color for high CCC values.
#' @param mid_color Color for mid CCC values (typically near 0).
#' @param value_text_size Text size for numbers in the heatmap.
#' @param ... Additional arguments passed to underlying functions
#' (like \code{theme} or \code{print}).
#' @importFrom ggplot2 .data
#' @export
plot.ccc <- function(x,
                     title = "Lin's Concordance Correlation Heatmap",
                     low_color = "indianred1",
                     high_color = "steelblue1",
                     mid_color = "white",
                     value_text_size = 4, ...) {

  if (!inherits(x, "ccc")) {
    stop("x must be of class 'ccc'.")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }

  mat <- if (is.list(x) && !is.null(x$est)) x$est else as.matrix(x)

  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("Var1", "Var2", "CCC")
  df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))

  p <- ggplot2::ggplot(df, ggplot2::aes(Var2, Var1, fill = CCC)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$CCC)),
                       size = value_text_size, color = "black") +
    ggplot2::scale_fill_gradient2(
      low = low_color, high = high_color, mid = mid_color,
      midpoint = 0, limit = c(-1, 1), name = "CCC"
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
