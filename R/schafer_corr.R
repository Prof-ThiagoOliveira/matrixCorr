#' @title Schafer-Strimmer shrinkage correlation
#'
#' @description
#' Computes a Schafer-Strimmer shrinkage correlation matrix for numeric data
#' using a high-performance 'C++' backend. This stabilises Pearson correlation
#' estimates by shrinking off-diagonal entries towards zero.
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
#'   \item \code{description} = \code{"Schafer-Strimmer shrinkage correlation
#'   matrix"}
#'   \item \code{package} = \code{"matrixCorr"}
#' }
#' Columns with zero variance are set to \code{NA} across row/column (including
#' the diagonal), matching \code{pearson_corr()} behaviour.
#'
#' @details
#' Let \eqn{R} be the sample Pearson correlation matrix. The Schafer-Strimmer
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
#' Schafer, J. & Strimmer, K. (2005). A shrinkage approach to large-scale
#' covariance matrix estimation and implications for functional genomics.
#' \emph{Statistical Applications in Genetics and Molecular Biology}, 4(1).
#'
#' @examples
#' ## Multivariate normal with AR(1) dependence (Toeplitz correlation)
#' set.seed(1)
#' n <- 80; p <- 40; rho <- 0.6
#' d <- abs(outer(seq_len(p), seq_len(p), "-"))
#' Sigma <- rho^d
#'
#' X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
#' colnames(X) <- paste0("V", seq_len(p))
#'
#' Rshr <- schafer_corr(X)
#' print(Rshr, digits = 2, n = 6, max_vars = 6)
#' summary(Rshr)
#' plot(Rshr)
#'
#' ## Shrinkage typically moves the sample correlation closer to the truth
#' Rraw <- stats::cor(X)
#' off  <- upper.tri(Sigma, diag = FALSE)
#' mae_raw <- mean(abs(Rraw[off] - Sigma[off]))
#' mae_shr <- mean(abs(Rshr[off] - Sigma[off]))
#' print(c(MAE_raw = mae_raw, MAE_shrunk = mae_shr))
#' plot(Rshr, title = "Schafer-Strimmer shrinkage correlation")
#'
#' # Interactive viewing (requires shiny)
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   view_corr_shiny(Rshr)
#' }
#'
#' @seealso \code{\link{print.schafer_corr}}, \code{\link{plot.schafer_corr}},
#'   \code{\link{pearson_corr}}
#' @author Thiago de Paula Oliveira
#' @export
schafer_corr <- function(data) {
  numeric_data <- validate_corr_input(data)
  colnames_data <- colnames(numeric_data)

  # call the C++ backend
  result <- sss_cor_cpp(numeric_data)

  # dimnames and metadata
  colnames(result) <- rownames(result) <- colnames_data
  .mc_structure_corr_matrix(
    result,
    class_name = "schafer_corr",
    method = "schafer_shrinkage",
    description = "Schafer-Strimmer shrinkage correlation matrix"
  )
}

#' @rdname schafer_corr
#' @method print schafer_corr
#' @title Print Method for \code{schafer_corr} Objects
#'
#' @param x An object of class \code{schafer_corr}.
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
print.schafer_corr <- function(x, digits = 4, n = NULL, topn = NULL,
                               max_vars = NULL, width = NULL,
                               show_ci = NULL, ...) {
  .mc_print_corr_matrix(
    x,
    header = "Schafer-Strimmer shrinkage correlation matrix",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
}

#' @rdname schafer_corr
#' @method plot schafer_corr
#' @title Plot Method for \code{schafer_corr} Objects
#'
#' @param x An object of class \code{schafer_corr}.
#' @param title Plot title.
#' @param cluster Logical; if TRUE, reorder rows/cols by hierarchical clustering
#'        on distance \eqn{1 - r}.
#' @param hclust_method Linkage method for \code{hclust}; default \code{"complete"}.
#' @param triangle One of \code{"full"}, \code{"upper"}, \code{"lower"}.
#' Default to \code{upper}.
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles (subject to \code{value_text_limit}).
#' @param show_values Deprecated compatibility alias for \code{show_value}. If
#'   supplied, it overrides \code{show_value}.
#' @param value_text_limit Integer threshold controlling when values are drawn.
#' @param value_text_size Font size for values if shown.
#' @param palette Character; \code{"diverging"} (default) or \code{"viridis"}.
#' @param ... Additional arguments passed to \code{ggplot2::theme()}.
#'
#' @return A \code{ggplot} object.
#' @import ggplot2
#' @export
plot.schafer_corr <- function(
    x,
    title = "Schafer-Strimmer shrinkage correlation",
    cluster = TRUE,
    hclust_method = "complete",
    triangle = c("upper", "lower", "full"),
    show_value = TRUE,
    show_values = NULL,
    value_text_limit = 60,
    value_text_size = 3,
    palette = c("diverging", "viridis"),
    ...
) {
  check_inherits(x, "schafer_corr")
  if (!is.null(show_values)) {
    show_value <- show_values
  }
  check_bool(show_value, arg = "show_value")
  triangle <- match.arg(triangle)
  palette  <- match.arg(palette)

  mat <- as.matrix(x)
  p   <- ncol(mat)

  # clustering to reveal blocks
  if (cluster && p > 1) {
    d <- stats::as.dist(1 - pmax(pmin(mat, 1), -1))
    hc <- stats::hclust(d, method = hclust_method)
    ord <- hc$order
    mat <- mat[ord, ord, drop = FALSE]
  }

  if (triangle == "upper") {
    mat[lower.tri(mat)] <- NA_real_
  } else if (triangle == "lower") {
    mat[upper.tri(mat)] <- NA_real_
  }
  rn <- rownames(mat); if (is.null(rn)) rn <- as.character(seq_len(p))
  cn <- colnames(mat); if (is.null(cn)) cn <- as.character(seq_len(p))
  df <- data.frame(
    Var1 = factor(rep(rn, each = p), levels = rev(rn)),
    Var2 = factor(rep(cn, times = p), levels = cn),
    r    = as.vector(mat),
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$r), , drop = FALSE]   # drop masked NAs

  # Palette
  if (palette == "diverging") {
    fill_scale <- ggplot2::scale_fill_gradient2(
      low = "indianred1", mid = "white", high = "steelblue1",
      midpoint = 0, limits = c(-1, 1), name = "r"
    )
  } else {
    if (!requireNamespace("viridisLite", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg viridisLite} is required for {.arg palette} = \"viridis\".")
    }
    fill_scale <- ggplot2::scale_fill_gradientn(
      colours = viridisLite::viridis(256, option = "B"),
      limits = c(-1, 1), name = "r"
    )
  }

  p_ <- ggplot2::ggplot(df, ggplot2::aes(x = Var2, y = Var1, fill = r)) +
    ggplot2::geom_raster() +
    fill_scale +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      ...
    )

  # Draw numbers only for small matrices
  if (show_value && p <= value_text_limit) {
    p_ <- p_ + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", r)),
      size = value_text_size, colour = "black"
    )
  } else if (p > 60) {
    # Hide tick labels for very large matrices
    p_ <- p_ + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )
  }

  p_
}

#' @rdname schafer_corr
#' @method summary schafer_corr
#' @param object An object of class \code{schafer_corr}.
#' @export
summary.schafer_corr <- function(object, n = NULL, topn = NULL,
                                 max_vars = NULL, width = NULL,
                                 show_ci = NULL, ...) {
  .mc_summary_corr_matrix(object, topn = topn)
}
