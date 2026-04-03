#' @title Pairwise Distance Correlation (dCor)
#'
#' @description
#' Computes pairwise distance correlations for the numeric columns of a matrix
#' or data frame using a high-performance 'C++' backend. Distance correlation
#' detects general dependence, including non-linear relationships.
#'
#' @param data A numeric matrix or a data frame with at least two numeric
#' columns. All non-numeric columns are dropped. Columns must be numeric.
#' @param check_na Logical (default \code{TRUE}). When \code{TRUE}, inputs must
#' be free of \code{NA}/\code{NaN}/\code{Inf}. Set to \code{FALSE} only if you
#' have already handled missingness upstream.
#'
#' @return A symmetric numeric matrix where the \code{(i, j)} entry is the
#' unbiased distance correlation between the \code{i}-th and \code{j}-th
#' numeric columns. The object has class \code{dcor} with attributes
#' \code{method = "distance_correlation"}, \code{description}, and
#' \code{package = "matrixCorr"}.
#'
#' @details
#' Let \eqn{x \in \mathbb{R}^n} and \eqn{D^{(x)}} be the pairwise distance matrix
#' with zero diagonal: \eqn{D^{(x)}_{ii} = 0}, \eqn{D^{(x)}_{ij} = |x_i - x_j|} for
#' \eqn{i \neq j}. Define row sums \eqn{r^{(x)}_i = \sum_{k \neq i} D^{(x)}_{ik}} and
#' grand sum \eqn{S^{(x)} = \sum_{i \neq k} D^{(x)}_{ik}}. The U-centred matrix is
#' \deqn{A^{(x)}_{ij} =
#'       \begin{cases}
#'         D^{(x)}_{ij} - \dfrac{r^{(x)}_i + r^{(x)}_j}{n - 2}
#'         + \dfrac{S^{(x)}}{(n - 1)(n - 2)}, & i \neq j,\\[6pt]
#'         0, & i = j~.
#'       \end{cases}}
#' For two variables \eqn{x,y}, the unbiased distance covariance and variances are
#' \deqn{\widehat{\mathrm{dCov}}^2_u(x,y) = \frac{2}{n(n-3)} \sum_{i<j} A^{(x)}_{ij} A^{(y)}_{ij}
#' \;=\; \frac{1}{n(n-3)} \sum_{i \neq j} A^{(x)}_{ij} A^{(y)}_{ij},}
#' with \eqn{\widehat{\mathrm{dVar}}^2_u(x)} defined analogously from \eqn{A^{(x)}}.
#' The unbiased distance correlation is
#' \deqn{\widehat{\mathrm{dCor}}_u(x,y) =
#'       \frac{\widehat{\mathrm{dCov}}_u(x,y)}
#'            {\sqrt{\widehat{\mathrm{dVar}}_u(x)\,\widehat{\mathrm{dVar}}_u(y)}} \in [0,1].}
#'
#' \strong{Computation.} All heavy lifting (distance matrices, U-centering,
#' and unbiased scaling) is implemented in C++ (`ustat_dcor_matrix_cpp`), so
#' the R wrapper only validates/coerces the input. OpenMP parallelises the
#' upper-triangular loops. The implementation includes a Huo-Szekely style
#' univariate \eqn{O(n \log n)} dispatch for pairwise terms. We also have an exact
#' unbiased \eqn{O(n^2)} fallback retained for robustness in small-sample or
#' non-finite-path cases; no external dependencies are used.
#'
#' @note Requires \eqn{n \ge 4}. Columns with (near) zero unbiased distance
#' variance yield \code{NA} in their row/column. Typical per-pair cost uses
#' the \eqn{O(n \log n)} fast path, with \eqn{O(n^2)} fallback when needed.
#'
#' @references
#' Székely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007).
#' Measuring and testing dependence by correlation of distances.
#' \emph{Annals of Statistics}, 35(6), 2769–2794.
#'
#' Székely, G. J., & Rizzo, M. L. (2013).
#' The distance correlation t-test of independence.
#' \emph{Journal of Multivariate Analysis}, 117, 193-213.
#'
#' @examples
#' ##Independent variables -> dCor ~ 0
#' set.seed(1)
#' X <- cbind(a = rnorm(200), b = rnorm(200))
#' D <- dcor(X)
#' print(D, digits = 3)
#' summary(D)
#'
#' ## Non-linear dependence: Pearson ~ 0, but unbiased dCor > 0
#' set.seed(42)
#' n <- 200
#' x <- rnorm(n)
#' y <- x^2 + rnorm(n, sd = 0.2)
#' XY <- cbind(x = x, y = y)
#' D2 <- dcor(XY)
#' # Compare Pearson vs unbiased distance correlation
#' round(c(pearson = cor(XY)[1, 2], dcor = D2["x", "y"]), 3)
#' summary(D2)
#' plot(D2, title = "Unbiased distance correlation (non-linear example)")
#'
#' ## Small AR(1) multivariate normal example
#' set.seed(7)
#' p <- 5; n <- 150; rho <- 0.6
#' Sigma <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
#' X3 <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
#' colnames(X3) <- paste0("V", seq_len(p))
#' D3 <- dcor(X3)
#' print(D3[1:3, 1:3], digits = 2)
#'
#' # Interactive viewing (requires shiny)
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   view_corr_shiny(D)
#' }
#'
#' @author Thiago de paula Oliveira
#'
#' @export
dcor <- function(data, check_na = TRUE) {
  numeric_data <- validate_corr_input(data, check_na = check_na)
  colnames_data <- colnames(numeric_data)

  dcor_matrix <- ustat_dcor_matrix_cpp(numeric_data)
  colnames(dcor_matrix) <- rownames(dcor_matrix) <- colnames_data

  .mc_structure_corr_matrix(
    dcor_matrix,
    class_name = "dcor",
    method = "distance_correlation",
    description = "Pairwise distance correlation matrix (unbiased)"
  )
}

#' @rdname dcor
#' @method print dcor
#' @title Print Method for \code{dcor} Objects
#'
#' @param x An object of class \code{dcor}.
#' @param digits Integer; number of decimal places to print.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.dcor <- function(x, digits = 4, n = NULL, topn = NULL,
                       max_vars = NULL, width = NULL,
                       show_ci = NULL, ...) {
  .mc_print_corr_matrix(
    x,
    header = "Distance correlation (dCor) matrix",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
}

#' @rdname dcor
#' @method plot dcor
#' @title Plot Method for \code{dcor} Objects
#'
#' @param x An object of class \code{dcor}.
#' @param title Plot title. Default is \code{"Distance correlation heatmap"}.
#' @param low_color Colour for zero correlation. Default is \code{"white"}.
#' @param high_color Colour for strong correlation. Default is \code{"steelblue1"}.
#' @param value_text_size Font size for displaying values. Default is \code{4}.
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
#' @param ... Additional arguments passed to \code{ggplot2::theme()} or other
#' \code{ggplot2} layers.
#'
#' @return A \code{ggplot} object representing the heatmap.
#' @import ggplot2
#' @export
plot.dcor <-
  function(x, title = "Distance correlation heatmap",
           low_color = "white", high_color = "steelblue1",
           value_text_size = 4, show_value = TRUE, ...) {

    check_inherits(x, "dcor")
    check_bool(show_value, arg = "show_value")
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg ggplot2} is required for plotting.")
    }

    mat <- as.matrix(x)
    df <- as.data.frame(as.table(mat))
    colnames(df) <- c("Var1", "Var2", "dCor")

    df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))

    p <- ggplot2::ggplot(df, ggplot2::aes(Var2, Var1, fill = dCor)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient(
        low = low_color, high = high_color,
        limits = c(0, 1), name = "dCor"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank(),
        ...
      ) +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = title, x = NULL, y = NULL)

    if (isTRUE(show_value) && !is.null(value_text_size) && is.finite(value_text_size)) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2f", dCor)),
        size = value_text_size,
        color = "black"
      )
    }

    p
  }

#' @rdname dcor
#' @method summary dcor
#' @param object An object of class \code{dcor}.
#' @export
summary.dcor <- function(object, n = NULL, topn = NULL,
                         max_vars = NULL, width = NULL,
                         show_ci = NULL, ...) {
  .mc_summary_corr_matrix(object, topn = topn)
}
