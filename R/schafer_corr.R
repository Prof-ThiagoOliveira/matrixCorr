#' @title Schäfer-Strimmer shrinkage correlation
#'
#' @description
#' Computes a shrinkage correlation matrix using the Schäfer-Strimmer approach,
#' with an analytic, data-driven shrinkage intensity \eqn{\hat\lambda}. The
#' method shrinks off-diagonal correlations towards zero while keeping the
#' diagonal at one, stabilising estimates when \eqn{p} is large relative to
#' \eqn{n}.
#'
#' This function uses a high-performance \code{C++} backend that forms
#' \eqn{X^\top X} via BLAS \code{SYRK}, applies centring via a rank-1 update,
#' converts to Pearson correlation, estimates \eqn{\hat\lambda}, and shrinks
#' the off-diagonals:
#' \eqn{R_{\mathrm{shr}} = (1-\hat\lambda)R + \hat\lambda I}.
#'
#' @param data A numeric matrix or a data frame with at least two numeric
#' columns. All non-numeric columns will be excluded. Columns must be numeric
#' and contain no \code{NA}s.
#'
#' @return A symmetric numeric matrix of class \code{schafer_corr} where entry
#' \code{(i, j)} is the shrunk correlation between the \code{i}-th and
#' \code{j}-th numeric columns. Attributes:
#' \itemize{
#'   \item \code{method} = \code{"schafer_shrinkage"}
#'   \item \code{description} = \code{"Schäfer-Strimmer shrinkage correlation
#'   matrix"}
#'   \item \code{package} = \code{"matrixCorr"}
#' }
#' Columns with zero variance are set to \code{NA} across row/column (including
#' the diagonal), matching \code{pearson_corr()} behaviour.
#'
#' @details
#' Let \eqn{R} be the sample Pearson correlation matrix. The Schäfer-Strimmer
#' shrinkage estimator targets the identity in correlation space and uses
#' \eqn{\hat\lambda = \frac{\sum_{i<j}\widehat{\mathrm{Var}}(r_{ij})}
#' {\sum_{i<j} r_{ij}^2}} (clamped to \eqn{[0,1]}), where
#' \eqn{\widehat{\mathrm{Var}}(r_{ij}) \approx \frac{(1-r_{ij}^2)^2}{n-1}}.
#' The returned estimator is \eqn{R_{\mathrm{shr}} = (1-\hat\lambda)R +
#' \hat\lambda I}.
#'
#' @note
#' No missing values are permitted. Columns with fewer than two observations
#' or zero variance are flagged as \code{NA} (row/column).
#'
#' @references
#' Schäfer, J. & Strimmer, K. (2005). A shrinkage approach to large-scale
#' covariance matrix estimation and implications for functional genomics.
#' \emph{Statistical Applications in Genetics and Molecular Biology}, 4(1).
#'
#' @examples
#' set.seed(1)
#' n <- 50; p <- 200
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("V", seq_len(p))
#' Rshr <- schafer_corr(X)
#' print(Rshr, digits = 2, max_rows = 6, max_cols = 6)
#' \donttest{
#'   plot(Rshr, title = "Schäfer-Strimmer shrinkage correlation (subset)")
#' }
#'
#' @seealso \code{\link{print.schafer_corr}}, \code{\link{plot.schafer_corr}},
#'   \code{\link{pearson_corr}}
#' @author Thiago de Paula Oliveira \email{toliveira@@abacusbio.com}
#' @export
schafer_corr <- function(data) {
  numeric_data <- validate_corr_input(data)
  colnames_data <- colnames(numeric_data)

  # call the C++ backend
  result <- sss_cor_cpp(numeric_data)

  # dimnames and metadata
  colnames(result) <- rownames(result) <- colnames_data
  attr(result, "method") <- "schafer_shrinkage"
  attr(result, "description") <- "Schafer-Strimmer shrinkage correlation matrix"
  attr(result, "package") <- "matrixCorr"

  class(result) <- c("schafer_corr", "matrix")
  result
}

#' @rdname schafer_corr
#' @method print schafer_corr
#' @title Print Method for \code{schafer_corr} Objects
#'
#' @description Prints a summary of the shrinkage correlation matrix with
#' optional truncation for large objects.
#'
#' @param x An object of class \code{schafer_corr}.
#' @param digits Integer; number of decimal places to print.
#' @param max_rows Optional integer; maximum number of rows to display.
#' If \code{NULL}, all rows are shown.
#' @param max_cols Optional integer; maximum number of columns to display.
#' If \code{NULL}, all columns are shown.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.schafer_corr <- function(x, digits = 3, max_rows = NULL,
                               max_cols = NULL, ...) {
  cat("Schafer-Strimmer shrinkage correlation matrix:\n")
  m <- as.matrix(x)
  attributes(m) <- attributes(m)[c("dim", "dimnames")]

  # Truncate display for large matrices
  if (!is.null(max_rows) || !is.null(max_cols)) {
    nr <- nrow(m); nc <- ncol(m)
    r  <- if (is.null(max_rows)) nr else min(nr, max_rows)
    c  <- if (is.null(max_cols)) nc else min(nc, max_cols)
    m2 <- round(m[seq_len(r), seq_len(c), drop = FALSE], digits)
    print(m2, ...)
    if (nr > r || nc > c) {
      cat(sprintf("... omitted: %d rows, %d cols\n", nr - r, nc - c))
    }
  } else {
    print(round(m, digits), ...)
  }

  invisible(x)
}

#' @rdname schafer_corr
#' @method plot schafer_corr
#' @title Plot Method for \code{schafer_corr} Objects
#'
#' @description Generates a \pkg{ggplot2}-based heatmap of the shrinkage
#' correlation matrix.
#'
#' @param x An object of class \code{schafer_corr}.
#' @param title Plot title. Default is
#' \code{"Schafer-Strimmer shrinkage correlation heatmap"}.
#' @param low_color Colour for the minimum correlation. Default is
#' \code{"indianred1"}.
#' @param high_color Colour for the maximum correlation. Default is
#' \code{"steelblue1"}.
#' @param mid_color Colour for zero correlation. Default is \code{"white"}.
#' @param value_text_size Font size for correlation values. Default is \code{4}.
#' @param ... Additional arguments passed to \code{ggplot2::theme()} or other
#' \code{ggplot2} layers.
#'
#' @return A \code{ggplot} object.
#' @import ggplot2
#' @export
plot.schafer_corr <-
  function(x,
           title = "Schafer-Strimmer shrinkage correlation heatmap",
           low_color = "indianred1",
           high_color = "steelblue1",
           mid_color = "white",
           value_text_size = 4, ...) {

    if (!inherits(x, "schafer_corr")) {
      stop("x must be of class 'schafer_corr'.")
    }

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required for plotting.")
    }

    mat <- as.matrix(x)
    df <- as.data.frame(as.table(mat))
    colnames(df) <- c("Var1", "Var2", "Shrinkage")

    df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))

    p <- ggplot2::ggplot(df, ggplot2::aes(Var2, Var1, fill = Shrinkage)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Shrinkage)),
                         size = value_text_size, color = "black") +
      ggplot2::scale_fill_gradient2(
        low = low_color, high = high_color, mid = mid_color,
        midpoint = 0, limit = c(-1, 1), name = "Shrinkage r"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank(),
        ...
      ) +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = title, x = NULL, y = NULL)

    p
  }
