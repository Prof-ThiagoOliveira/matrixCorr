#' @title Pairwise Pearson correlation
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
#' \strong{Statistical formulation.} Let \eqn{X \in \mathbb{R}^{n \times p}} be a
#' numeric matrix with rows as observations and columns as variables, and let
#' \eqn{\mu = \tfrac{1}{n}\mathbf{1}^\top X} be the vector of column means.
#' Define the centred cross-product matrix
#' \deqn{S \;=\; (X - \mathbf{1}\mu)^\top (X - \mathbf{1}\mu)
#'       \;=\; X^\top X \;-\; n\,\mu\,\mu^\top.}
#' The (unbiased) sample covariance is then
#' \deqn{\widehat{\Sigma} \;=\; \tfrac{1}{n-1}\,S,}
#' and the vector of sample standard deviations is
#' \deqn{s_i \;=\; \sqrt{\widehat{\Sigma}_{ii}}, \qquad i=1,\ldots,p.}
#' The Pearson correlation matrix \eqn{R} is obtained by standardising
#' \eqn{\widehat{\Sigma}}, given by:
#' \deqn{R \;=\; D^{-1/2}\,\widehat{\Sigma}\,D^{-1/2}, \qquad
#'       D \;=\; \mathrm{diag}(\widehat{\Sigma}_{11},\ldots,\widehat{\Sigma}_{pp}),}
#' equivalently, entrywise
#' \deqn{R_{ij} \;=\; \frac{\widehat{\Sigma}_{ij}}{s_i\,s_j}, \quad i \neq j,
#'       \qquad R_{ii} \;=\; 1.}
#'
#' If \eqn{s_i = 0} (zero variance),
#' the \eqn{i}-th row and column are set to \code{NA}. Tiny negative values on
#' the covariance diagonal caused by floating-point rounding are reduced to
#' zero before taking square roots. No missing values are permitted in \eqn{X}.
#'
#' The implementation forms \eqn{X^\top X} via a
#' symmetric rank-\eqn{k} update (BLAS \code{SYRK}) on the upper triangle, then
#' applies the rank-1 correction \eqn{-\,n\,\mu\,\mu^\top} to obtain \eqn{S}
#' without explicitly materialising \eqn{X - \mathbf{1}\mu}. After scaling by
#' \eqn{1/(n-1)}, triangular normalisation by \eqn{D^{-1/2}} yields \eqn{R},
#' which is then symmetrised. The dominant cost is \eqn{O(n p^2)} flops with
#' \eqn{O(p^2)} memory.
#'
#' This implementation avoids calculating means explicitly and instead
#' uses a numerically stable form.
#'
#' @note Missing values are not allowed. Columns with fewer than two
#' observations are excluded.
#'
#' @references
#' Pearson, K. (1895). "Notes on regression and inheritance in the case of
#' two parents". Proceedings of the Royal Society of London, 58, 240–242.
#'
#' @examples
#' ## MVN with AR(1) correlation
#' set.seed(123)
#' p <- 6; n <- 300; rho <- 0.5
#' # true correlation
#' Sigma <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
#' L <- chol(Sigma)
#' # MVN(n, 0, Sigma)
#' X <- matrix(rnorm(n * p), n, p) %*% L
#' colnames(X) <- paste0("V", seq_len(p))
#'
#' pr <- pearson_corr(X)
#' print(pr, digits = 2)
#' plot(pr)
#'
#' ## Compare the sample estimate to the truth
#' Rhat <- cor(X)
#' # estimated
#' round(Rhat[1:4, 1:4], 2)
#' # true
#' round(Sigma[1:4, 1:4], 2)
#' off <- upper.tri(Sigma, diag = FALSE)
#' # MAE on off-diagonals
#' mean(abs(Rhat[off] - Sigma[off]))
#'
#' ## Larger n reduces sampling error
#' n2 <- 2000
#' X2 <- matrix(rnorm(n2 * p), n2, p) %*% L
#' Rhat2 <- cor(X2)
#' off <- upper.tri(Sigma, diag = FALSE)
#' ## mean absolute error (MAE) of the off-diagonal correlations
#' mean(abs(Rhat2[off] - Sigma[off]))
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
#' @param digits Integer; number of decimal places to print in the concordance
#' @param max_rows Optional integer; maximum number of rows to display.
#'  If \code{NULL}, all rows are shown.
#' @param max_cols Optional integer; maximum number of columns to display.
#' If \code{NULL}, all columns are shown.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @return Invisibly returns the \code{pearson_corr} object.
#' @author Thiago de Paula Oliveira
#' @export
print.pearson_corr <- function(x, digits = 3, max_rows = NULL,
                               max_cols = NULL, ...) {
  cat("Pearson correlation matrix:\n")
  m <- as.matrix(x)
  attributes(m) <- attributes(m)[c("dim", "dimnames")]

  # Truncation for large matrices
  if (!is.null(max_rows) || !is.null(max_cols)) {
    nr <- nrow(m); nc <- ncol(m)
    r  <- if (is.null(max_rows)) nr else min(nr, max_rows)
    c  <- if (is.null(max_cols)) nc else min(nc, max_cols)
    m  <- m[seq_len(r), seq_len(c), drop = FALSE]
    m  <- round(m, digits)
    print(m, ...)
    if (nr > r || nc > c) {
      cat(sprintf("... omitted: %d rows, %d cols\n", nr - r, nc - c))
    }
  } else {
    print(round(m, digits), ...)
  }

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

