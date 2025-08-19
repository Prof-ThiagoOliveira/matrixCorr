#' @title Partial correlation matrix (sample / ridge / OAS)
#'
#' @description
#' Computes the Gaussian partial correlation matrix from a numeric data frame
#' or matrix, using one of three covariance estimators: unbiased sample
#' covariance, ridge-regularised covariance, or OAS shrinkage to a scaled
#' identity (recommended when \eqn{p \gg n}).
#'
#' @param data A numeric matrix or data frame with at least two numeric columns.
#'   Non-numeric columns are ignored.
#' @param method Character; one of \code{"oas"}, \code{"ridge"}, \code{"sample"}.
#'   Default \code{"oas"}.
#' @param lambda Numeric \eqn{\ge 0}; ridge penalty added to the covariance
#'   diagonal when \code{method = "ridge"}. Ignored otherwise. Default \code{1e-3}.
#' @param return_cov_precision Logical; if \code{TRUE}, also return the
#'   covariance (\code{cov}) and precision (\code{precision}) matrices used to
#'   form the partial correlations. Default to \code{FALSE}
#'
#' @return An object of class \code{"partial_corr"} (a list) with elements:
#'   \itemize{
#'     \item \code{pcor}: \eqn{p \times p} partial correlation matrix.
#'     \item \code{cov} (if requested): covariance matrix used.
#'     \item \code{precision} (if requested): precision matrix \eqn{\Theta}.
#'     \item \code{method}: the estimator used (\code{"oas"}, \code{"ridge"},
#'     or \code{"sample"}).
#'     \item \code{lambda}: ridge penalty (or \code{NA_real_}).
#'     \item \code{rho}: OAS shrinkage weight in \eqn{[0,1]} (or \code{NA_real_}).
#'     \item \code{jitter}: diagonal jitter added (if any) to ensure positive
#'     definiteness.
#'   }
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(2000), nrow = 100, ncol = 20)
#' pc <- partial_correlation(X, method = "oas")
#' print(pc, digits = 2)
#' plot(pc)
#'
#' @references
#' Chen, Y., Wiesel, A., & Hero, A. O. III (2011).
#' Robust Shrinkage Estimation of High-dimensional Covariance Matrices.
#' IEEE Transactions on Signal Processing.
#'
#' @export
partial_correlation <- function(data, method = c("oas","ridge","sample"),
                                lambda = 1e-3, return_cov_precision = FALSE) {
  method <- match.arg(method)

  numeric_data <-
    if (is.matrix(data) && is.double(data) && all(is.finite(data))) {
      data
    } else {
      # drops non-numeric
      validate_corr_input(data)
    }

  res <- partial_correlation_cpp(numeric_data, method, lambda, return_cov_precision)

  # set dimnames (cheap; attributes only)
  dn <- list(colnames(numeric_data), colnames(numeric_data))
  if (!is.null(res$pcor)) {
    dimnames(res$pcor) <- dn
    if (!is.null(res$cov))       dimnames(res$cov)       <- dn
    if (!is.null(res$precision)) dimnames(res$precision) <- dn
  } else {
    pcor <- res[[1]]; dimnames(pcor) <- dn; res <- list(pcor = pcor)
  }

  res$method <- method
  res$lambda <- if (identical(method, "ridge")) lambda else NA_real_
  res$rho    <- if (identical(method, "oas"))   res$rho %||% NA_real_ else NA_real_
  res$jitter <- res$jitter %||% NA_real_
  class(res) <- c("partial_corr", "list")
  res
}


# small helper for older R versions without %||%
`%||%` <- function(a, b) if (!is.null(a)) a else b


#' @title Print method for \code{partial_corr}
#' @description Prints only the partial correlation matrix (no attribute spam),
#'   with an optional one-line header stating the estimator used.
#' @param x An object of class \code{partial_corr}.
#' @param digits Integer; number of decimal places for display (default 3).
#' @param show_method Logical; print a one-line header with \code{method}
#'   (and \code{lambda}/\code{rho} if available). Default \code{TRUE}.
#' @param max_rows,max_cols Optional integer limits for display; if provided,
#'   the printed matrix is truncated with a note about omitted rows/cols.
#' @param ... Passed to \code{print}.
#' @return Invisibly returns \code{x}.
#' @export
print.partial_correlation <-
  function(x,
           digits = 3,
           show_method = TRUE,
           max_rows = NULL,
           max_cols = NULL,
           ...) {
  stopifnot(inherits(x, "partial_corr"))
  if (!is.matrix(x$pcor)) stop("`x$pcor` must be a matrix.")

  if (isTRUE(show_method)) {
    meth <- x$method %||% NA_character_
    hdr <- switch(
      tolower(as.character(meth)),
      "oas"   = {
        rho <- if (!is.null(x$rho) && is.finite(x$rho)) sprintf(", OAS rho=%.3f", x$rho) else ""
        paste0("Partial correlation (OAS", rho, ")")
      },
      "ridge" = {
        lam <- if (!is.null(x$lambda) && is.finite(x$lambda)) sprintf(", lambda=%.3g", x$lambda) else ""
        paste0("Partial correlation (ridge", lam, ")")
      },
      "sample" = "Partial correlation (sample covariance)",
      "Partial correlation"
    )
    cat(hdr, "\n")
  } else {
    cat("Partial correlation matrix:\n")
  }

  # Prepare matrix
  M <- x$pcor
  attributes(M) <- attributes(M)[c("dim", "dimnames")]

  if (!is.null(max_rows) || !is.null(max_cols)) {
    nr <- nrow(M); nc <- ncol(M)
    r  <- if (is.null(max_rows)) nr else min(nr, max_rows)
    c  <- if (is.null(max_cols)) nc else min(nc, max_cols)
    mm <- round(M[seq_len(r), seq_len(c), drop = FALSE], digits)
    print(mm, ...)
    if (nr > r || nc > c) {
      cat(sprintf("... omitted: %d rows, %d cols\n", nr - r, nc - c))
    }
  } else {
    print(round(M, digits), ...)
  }

  invisible(x)
}


#' @rdname partial_correlation
#' @method plot partial_corr
#' @title Plot Method for \code{partial_corr} Objects
#'
#' @description
#' Produces a \pkg{ggplot2}-based heatmap of the partial correlation matrix
#' stored in \code{x$pcor}. Optionally masks the diagonal and/or reorders
#' variables via hierarchical clustering of \eqn{1 - |pcor|}.
#'
#' @param x An object of class \code{partial_corr}.
#' @param title Plot title. By default, constructed from the estimator in
#' \code{x$method}.
#' @param low_color Colour for low (negative) values. Default
#' \code{"indianred1"}.
#' @param high_color Colour for high (positive) values. Default
#' \code{"steelblue1"}.
#' @param mid_color Colour for zero. Default \code{"white"}.
#' @param value_text_size Font size for cell labels. Default \code{4}.
#' @param mask_diag Logical; if \code{TRUE}, the diagonal is masked
#' (set to \code{NA}) and not labelled. Default \code{TRUE}.
#' @param reorder Logical; if \code{TRUE}, variables are reordered by
#' hierarchical clustering of \eqn{1 - |pcor|}. Default \code{FALSE}.
#' @param ... Additional arguments passed to \code{ggplot2::theme()} or other
#'   \pkg{ggplot2} layers.
#'
#' @return A \code{ggplot} object.
#' @import ggplot2
#' @importFrom stats as.dist hclust
#' @export
plot.partial_corr <- function(
    x,
    title = NULL,
    low_color  = "indianred1",
    high_color = "steelblue1",
    mid_color  = "white",
    value_text_size = 4,
    mask_diag = TRUE,
    reorder   = FALSE,
    ...
) {
  if (!inherits(x, "partial_corr")) stop("x must be of class 'partial_corr'.")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plotting.")

  M <- x$pcor
  if (!is.matrix(M)) stop("`x$pcor` must be a matrix.")

  # Ensure dimnames for labelling
  if (is.null(colnames(M))) colnames(M) <- paste0("V", seq_len(ncol(M)))
  if (is.null(rownames(M))) rownames(M) <- colnames(M)

  # Optional reordering by hierarchical clustering of 1 - |pcor|
  if (isTRUE(reorder) && nrow(M) >= 2L) {
    D  <- stats::as.dist(1 - pmin(1, abs(M)))
    hc <- stats::hclust(D, method = "average")
    ord <- hc$order
    M <- M[ord, ord, drop = FALSE]
  }

  # Default title constructed from x$method (+ tuning, if present)
  if (is.null(title)) {
    method <- tolower(as.character(x$method %||% ""))
    extra <- switch(
      method,
      "oas"   = if (is.finite(x$rho %||% NA_real_)) sprintf(" (OAS, rho=%.3f)", x$rho) else " (OAS)",
      "ridge" = if (is.finite(x$lambda %||% NA_real_)) sprintf(" (ridge, lambda=%.3g)", x$lambda) else " (ridge)",
      "sample" = " (sample)",
      ""
    )
    title <- paste0("Partial correlation heatmap", extra)
  }

  # Build plotting frame
  df <- as.data.frame(as.table(M))
  names(df) <- c("Var1", "Var2", "PCor")

  # Mask diagonal if requested
  if (isTRUE(mask_diag)) {
    df$PCor[df$Var1 == df$Var2] <- NA_real_
  }

  # Reverse y-axis order for a tidy heatmap
  df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))

  ggplot2::ggplot(df, ggplot2::aes(x = Var2, y = Var1, fill = PCor)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.na(PCor), "", sprintf("%.2f", PCor))),
      size = value_text_size, colour = "black"
    ) +
    ggplot2::scale_fill_gradient2(
      low = low_color, high = high_color, mid = mid_color,
      midpoint = 0, limits = c(-1, 1), name = "Partial r",
      na.value = "grey95"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid  = ggplot2::element_blank(),
      ...
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = NULL, y = NULL)
}
