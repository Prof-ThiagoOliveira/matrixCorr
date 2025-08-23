#' @title Pairwise Lin's concordance correlation coefficient
#'
#' @description
#' Computes all pairwise Lin's Concordance Correlation Coefficients (CCC)
#' from the numeric columns of a matrix or data frame. CCC measures both
#' precision (Pearson correlation) and accuracy (closeness to the 45-degree line).
#' This function is backed by a high-performance 'C++' implementation.
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
#' # Example with multivariate normal data
#' Sigma <- matrix(c(1, 0.5, 0.3,
#'                   0.5, 1, 0.4,
#'                   0.3, 0.4, 1), nrow = 3)
#' mu <- c(0, 0, 0)
#' set.seed(123)
#' mat_mvn <- MASS::mvrnorm(n = 100, mu = mu, Sigma = Sigma)
#' result_mvn <- ccc(mat_mvn)
#' print(result_mvn)
#' plot(result_mvn)
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

  if (verbose) cat("Using", openmp_threads(), "OpenMP threads\n")

  if (ci) {
    ccc_lin <- ccc_with_ci_cpp(mat, conf_level)
    ccc_lin$est    <- `dimnames<-`(ccc_lin$est,
                                   list(colnames_data, colnames_data))
    ccc_lin$lwr.ci <- `dimnames<-`(ccc_lin$lwr.ci,
                                   list(colnames_data, colnames_data))
    ccc_lin$upr.ci <- `dimnames<-`(ccc_lin$upr.ci,
                                   list(colnames_data, colnames_data))

    attr(ccc_lin, "method") <- "Lin's concordance"
    attr(ccc_lin, "description") <-
      "Pairwise Lin's concordance with confidence intervals"
    attr(ccc_lin, "package") <- "matrixCorr"
    class(ccc_lin) <- c("ccc", "ccc_ci")   # list with CIs
  } else {
    est <- ccc_cpp(mat)
    ccc_lin <- `dimnames<-`(est, list(colnames_data, colnames_data))

    attr(ccc_lin, "method") <- "Lin's concordance"
    attr(ccc_lin, "description") <- "Pairwise Lin's concordance correlation matrix"
    attr(ccc_lin, "package") <- "matrixCorr"
    class(ccc_lin) <- c("ccc", "matrix")   # matrix printing still available
  }

  ccc_lin
}


#' @rdname ccc
#' @method print ccc
#' @param digits Integer; number of decimal places to print in the concordance
#' matrix (default is 4).
#' @export
print.ccc <- function(x, digits = 4, ...) {
  # helper to strip non-essential attributes
  strip_attrs <- function(m) {
    m <- as.matrix(m)
    attributes(m) <- attributes(m)[c("dim", "dimnames")]
    m
  }

  if (inherits(x, "ccc_ci") || (is.list(x) && all(c("est", "lwr.ci", "upr.ci") %in% names(x)))) {
    cat("Concordance matrix (estimates):\n")
    est <- strip_attrs(x$est)
    print(round(est, digits), ...)

    cat("\nConfidence intervals:\nLower:\n")
    lwr <- strip_attrs(x$lwr.ci)
    print(round(lwr, 2), ...)

    cat("\nUpper:\n")
    upr <- strip_attrs(x$upr.ci)
    print(round(upr, 2), ...)
  } else if (is.matrix(x)) {
    cat("Concordance matrix:\n")
    M <- strip_attrs(x)
    print(round(M, digits), ...)
  } else {
    stop("Invalid object format for class 'ccc'.")
  }

  invisible(x)
}

#' @rdname ccc
#' @method plot ccc
#' @param x An object of class \code{"ccc"} (either a matrix or a list with CIs).
#' @param title Title for the plot.
#' @param low_color Color for low CCC values.
#' @param high_color Color for high CCC values.
#' @param mid_color Color for mid CCC values.
#' @param value_text_size Text size for CCC values in the heatmap.
#' @param ci_text_size Text size for confidence intervals.
#' @param ... Passed to \code{ggplot2::theme()}.
#' @export
plot.ccc <- function(x,
                     title = "Lin's Concordance Correlation Heatmap",
                     low_color = "indianred1",
                     high_color = "steelblue1",
                     mid_color = "white",
                     value_text_size = 4,
                     ci_text_size = 3,
                     ...) {

  if (!inherits(x, "ccc"))
    stop("x must be of class 'ccc'.")

  # --- Build long data with proper alignment by (Var1, Var2) ---
  est_mat <- if (is.list(x) && !is.null(x$est)) x$est else unclass(x)
  df_est  <- as.data.frame(as.table(est_mat))
  names(df_est) <- c("Var1", "Var2", "CCC")

  if (is.list(x) && !is.null(x$lwr.ci) && !is.null(x$upr.ci)) {
    df_lwr <- as.data.frame(as.table(x$lwr.ci)); names(df_lwr)[3] <- "lwr"
    df_upr <- as.data.frame(as.table(x$upr.ci)); names(df_upr)[3] <- "upr"
    df <- Reduce(function(a, b) merge(a, b, by = c("Var1","Var2"), all = TRUE),
                 list(df_est, df_lwr, df_upr))

    # Blank CI on the diagonal (show 1.00 but no CI there)
    diag_idx <- df$Var1 == df$Var2
    df$lwr[diag_idx] <- NA_real_
    df$upr[diag_idx] <- NA_real_
    df$ci_label <- ifelse(is.na(df$lwr) | is.na(df$upr),
                          NA_character_,
                          sprintf("(%.2f, %.2f)", df$lwr, df$upr))
  } else {
    df <- df_est
    df$ci_label <- NA_character_
  }

  # Reverse Y axis for heatmap look, keep X in natural order
  lev_row <- unique(df_est$Var1)
  lev_col <- unique(df_est$Var2)
  df$Var1 <- factor(df$Var1, levels = rev(lev_row))
  df$Var2 <- factor(df$Var2, levels = lev_col)

  # Main estimate label
  df$label <- sprintf("%.2f", df$CCC)

  # --- Plot ---
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Var2, y = Var1, fill = CCC)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = label), size = value_text_size) +
    ggplot2::scale_fill_gradient2(
      low = low_color, high = high_color, mid = mid_color,
      midpoint = 0, limit = c(-1, 1), name = "CCC"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      ...
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL)

  if (any(!is.na(df$ci_label))) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = ci_label, y = as.numeric(Var1) - 0.25),
      size = ci_text_size, color = "gray30", na.rm = TRUE
    )
  }

  p
}
