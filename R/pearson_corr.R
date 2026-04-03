#' @title Pairwise Pearson correlation
#'
#' @description
#' Computes pairwise Pearson correlations for the numeric columns of a matrix
#' or data frame using a high-performance 'C++' backend. Optional Fisher-z
#' confidence intervals are available.
#'
#' @param data A numeric matrix or a data frame with at least two numeric
#' columns. All non-numeric columns will be excluded. Each column must have
#' at least two non-missing values.
#' @param check_na Logical (default \code{TRUE}). If \code{TRUE}, inputs must be
#' free of \code{NA}/\code{NaN}/\code{Inf}. Set to \code{FALSE} only when the
#' caller already handled missingness.
#' @param ci Logical (default \code{FALSE}). If \code{TRUE}, attach pairwise
#' Fisher-\eqn{z} confidence intervals for the off-diagonal Pearson
#' correlations.
#' @param conf_level Confidence level used when \code{ci = TRUE}. Default is
#' \code{0.95}.
#'
#' @return A symmetric numeric matrix where the \code{(i, j)}-th element is
#' the Pearson correlation between the \code{i}-th and \code{j}-th
#' numeric columns of the input. When \code{ci = TRUE}, the object also
#' carries a \code{ci} attribute with elements \code{est}, \code{lwr.ci},
#' \code{upr.ci}, and \code{conf.level}. When pairwise-complete evaluation is
#' used, pairwise sample sizes are stored in \code{attr(x, "diagnostics")$n_complete}.
#'
#' @details
#' Let \eqn{X \in \mathbb{R}^{n \times p}} be a
#' numeric matrix with rows as observations and columns as variables, and let
#' \eqn{\mathbf{1} \in \mathbb{R}^n} denote the all-ones vector. Define the column
#' means \eqn{\mu = (1/n)\,\mathbf{1}^\top X} and the centred cross-product
#' matrix
#' \deqn{ S \;=\; (X - \mathbf{1}\mu)^\top (X - \mathbf{1}\mu)
#'       \;=\; X^\top \!\Big(I_n - \tfrac{1}{n}\mathbf{1}\mathbf{1}^\top\Big) X
#'       \;=\; X^\top X \;-\; n\,\mu\,\mu^\top. }
#' The (unbiased) sample covariance is
#' \deqn{ \widehat{\Sigma} \;=\; \tfrac{1}{n-1}\,S, }
#' and the sample standard deviations are \eqn{s_i = \sqrt{\widehat{\Sigma}_{ii}}}.
#' The Pearson correlation matrix is obtained by standardising \eqn{\widehat{\Sigma}}, and it is given by
#' \deqn{ R \;=\; D^{-1/2}\,\widehat{\Sigma}\,D^{-1/2}, \qquad
#'       D \;=\; \mathrm{diag}(\widehat{\Sigma}_{11},\ldots,\widehat{\Sigma}_{pp}), }
#' equivalently, entrywise \eqn{R_{ij} = \widehat{\Sigma}_{ij}/(s_i s_j)} for
#' \eqn{i \neq j} and \eqn{R_{ii} = 1}. With \eqn{1/(n-1)} scaling,
#' \eqn{\widehat{\Sigma}} is unbiased for the covariance; the induced
#' correlations are biased in finite samples.
#'
#' The implementation forms \eqn{X^\top X} via a BLAS
#' symmetric rank-\eqn{k} update (SYRK) on the upper triangle, then applies the
#' rank-1 correction \eqn{-\,n\,\mu\,\mu^\top} to obtain \eqn{S} without
#' explicitly materialising \eqn{X - \mathbf{1}\mu}. After scaling by
#' \eqn{1/(n-1)}, triangular normalisation by \eqn{D^{-1/2}} yields \eqn{R},
#' which is then symmetrised to remove round-off asymmetry. Tiny negative values
#' on the covariance diagonal due to floating-point rounding are truncated to
#' zero before taking square roots.
#'
#' If a variable has zero variance (\eqn{s_i = 0}), the
#' corresponding row and column of \eqn{R} are set to \code{NA}. When
#' \code{check_na = FALSE}, each \eqn{(i,j)} correlation is recomputed on the
#' pairwise complete-case overlap of columns \eqn{i} and \eqn{j}.
#'
#' When \code{ci = TRUE}, Fisher-\eqn{z} confidence intervals are computed from
#' the observed pairwise Pearson correlation \eqn{r_{ij}} and the pairwise
#' complete-case sample size \eqn{n_{ij}}:
#' \deqn{
#' z_{ij} = \operatorname{atanh}(r_{ij}), \qquad
#' \operatorname{SE}(z_{ij}) = \frac{1}{\sqrt{n_{ij} - 3}}.
#' }
#' With \eqn{z_{1-\alpha/2} = \Phi^{-1}(1 - \alpha/2)}, the confidence limits are
#' \deqn{
#' \tanh\!\bigl(z_{ij} - z_{1-\alpha/2}\operatorname{SE}(z_{ij})\bigr)
#' \;\;\text{and}\;\;
#' \tanh\!\bigl(z_{ij} + z_{1-\alpha/2}\operatorname{SE}(z_{ij})\bigr).
#' }
#' Confidence intervals are reported only when \eqn{n_{ij} > 3}.
#'
#' \strong{Computational complexity.} The dominant cost is \eqn{O(n p^2)} flops
#' with \eqn{O(p^2)} memory.
#'
#' @note Missing values are not allowed when \code{check_na = TRUE}. Columns
#' with fewer than two observations are excluded.
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
#' summary(pr)
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
#' # Interactive viewing (requires shiny)
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   view_corr_shiny(pr)
#' }
#'
#' @useDynLib matrixCorr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @seealso \code{\link{print.pearson_corr}}, \code{\link{plot.pearson_corr}}
#' @author Thiago de Paula Oliveira
#' @export
pearson_corr <- function(data, check_na = TRUE, ci = FALSE, conf_level = 0.95) {
  check_bool(ci, arg = "ci")
  if (isTRUE(ci)) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  }

  numeric_data <- validate_corr_input(data, check_na = check_na)
  colnames_data <- colnames(numeric_data)
  diagnostics <- NULL
  ci_attr <- NULL

  if (isTRUE(check_na) && !isTRUE(ci)) {
    result <- pearson_matrix_cpp(numeric_data)
  } else {
    pairwise <- pearson_matrix_pairwise_cpp(
      numeric_data,
      return_ci = ci,
      conf_level = conf_level
    )
    result <- pairwise$est
    diagnostics <- list(n_complete = pairwise$n_complete)
    dimnames(diagnostics$n_complete) <- list(colnames_data, colnames_data)
    if (isTRUE(ci)) {
      ci_attr <- list(
        est = unclass(result),
        lwr.ci = unclass(pairwise$lwr),
        upr.ci = unclass(pairwise$upr),
        conf.level = pairwise$conf_level
      )
      dimnames(ci_attr$est) <- list(colnames_data, colnames_data)
      dimnames(ci_attr$lwr.ci) <- list(colnames_data, colnames_data)
      dimnames(ci_attr$upr.ci) <- list(colnames_data, colnames_data)
    }
  }

  colnames(result) <- rownames(result) <- colnames_data
  out <- .mc_structure_corr_matrix(
    result,
    class_name = "pearson_corr",
    method = "pearson",
    description = "Pairwise Pearson correlation matrix",
    diagnostics = diagnostics
  )
  if (!is.null(ci_attr)) {
    attr(out, "ci") <- ci_attr
    attr(out, "conf.level") <- conf_level
  }
  out
}

.mc_pearson_ci_attr <- function(x) {
  attr(x, "ci", exact = TRUE)
}

.mc_pearson_pairwise_summary <- function(object,
                                         digits = 4,
                                         ci_digits = 3,
                                         show_ci = NULL) {
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )
  check_inherits(object, "pearson_corr")

  est <- as.matrix(object)
  rn <- rownames(est); cn <- colnames(est)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(est)))
  if (is.null(cn)) cn <- as.character(seq_len(ncol(est)))

  ci <- .mc_pearson_ci_attr(object)
  diag_attr <- attr(object, "diagnostics", exact = TRUE)
  include_ci <- identical(show_ci, "yes") && !is.null(ci)

  rows <- vector("list", nrow(est) * (ncol(est) - 1L) / 2L)
  k <- 0L
  for (i in seq_len(nrow(est) - 1L)) {
    for (j in (i + 1L):ncol(est)) {
      k <- k + 1L
      rec <- list(
        var1 = rn[i],
        var2 = cn[j],
        estimate = round(est[i, j], digits)
      )
      if (is.list(diag_attr) && is.matrix(diag_attr$n_complete)) {
        rec$n_complete <- as.integer(diag_attr$n_complete[i, j])
      }
      if (include_ci) {
        rec$lwr <- if (!is.null(ci$lwr.ci) && is.finite(ci$lwr.ci[i, j])) round(ci$lwr.ci[i, j], ci_digits) else NA_real_
        rec$upr <- if (!is.null(ci$upr.ci) && is.finite(ci$upr.ci[i, j])) round(ci$upr.ci[i, j], ci_digits) else NA_real_
      }
      rows[[k]] <- rec
    }
  }

  df <- do.call(rbind.data.frame, rows)
  rownames(df) <- NULL
  if ("estimate" %in% names(df)) df$estimate <- as.numeric(df$estimate)
  if ("lwr" %in% names(df)) df$lwr <- as.numeric(df$lwr)
  if ("upr" %in% names(df)) df$upr <- as.numeric(df$upr)
  if ("n_complete" %in% names(df)) df$n_complete <- as.integer(df$n_complete)

  out <- structure(df, class = c("summary.pearson_corr", "data.frame"))
  attr(out, "overview") <- .mc_summary_corr_matrix(object)
  attr(out, "has_ci") <- include_ci
  attr(out, "conf.level") <- if (is.null(ci)) NA_real_ else ci$conf.level
  attr(out, "digits") <- digits
  attr(out, "ci_digits") <- ci_digits
  out
}


#' @rdname pearson_corr
#' @method print pearson_corr
#' @title Print Method for \code{pearson_corr} Objects
#'
#' @param x An object of class \code{pearson_corr}.
#' @param digits Integer; number of decimal places to print in the concordance
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param ci_digits Integer; digits for Pearson confidence limits.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @return Invisibly returns the \code{pearson_corr} object.
#' @author Thiago de Paula Oliveira
#' @export
print.pearson_corr <- function(x, digits = 4, n = NULL,
                               topn = NULL,
                               max_vars = NULL,
                               width = NULL,
                               ci_digits = 3,
                               show_ci = NULL,
                               ...) {
  .mc_print_corr_matrix(
    x,
    header = "Pearson correlation matrix",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
  invisible(x)
}

#' @rdname pearson_corr
#' @method plot pearson_corr
#' @title Plot Method for \code{pearson_corr} Objects
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
#' @param ci_text_size Text size for confidence intervals in the heatmap.
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
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
           mid_color = "white", value_text_size = 4,
           ci_text_size = 3, show_value = TRUE, ...) {
  check_bool(show_value, arg = "show_value")
  ci <- .mc_pearson_ci_attr(x)
  if (is.null(ci) || is.null(ci$lwr.ci) || is.null(ci$upr.ci)) {
    return(.mc_plot_corr_matrix(
      x, class_name = "pearson_corr", fill_name = "Pearson",
      title = title, low_color = low_color, high_color = high_color,
      mid_color = mid_color, value_text_size = value_text_size,
      show_value = show_value, ...
    ))
  }

  est_mat <- as.matrix(x)
  df_est <- as.data.frame(as.table(est_mat))
  names(df_est) <- c("Var1", "Var2", "pearson")

  df_lwr <- as.data.frame(as.table(ci$lwr.ci))
  names(df_lwr)[3] <- "lwr"
  df_upr <- as.data.frame(as.table(ci$upr.ci))
  names(df_upr)[3] <- "upr"
  df <- Reduce(
    function(a, b) merge(a, b, by = c("Var1", "Var2"), all = TRUE),
    list(df_est, df_lwr, df_upr)
  )

  diag_idx <- df$Var1 == df$Var2
  df$lwr[diag_idx] <- NA_real_
  df$upr[diag_idx] <- NA_real_
  df$ci_label <- ifelse(
    is.na(df$lwr) | is.na(df$upr),
    NA_character_,
    sprintf("[%.3f, %.3f]", df$lwr, df$upr)
  )

  lev_row <- unique(df_est$Var1)
  lev_col <- unique(df_est$Var2)
  df$Var1 <- factor(df$Var1, levels = rev(lev_row))
  df$Var2 <- factor(df$Var2, levels = lev_col)
  df$label <- sprintf("%.2f", df$pearson)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Var2, y = Var1, fill = .data$pearson)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(
      low = low_color,
      high = high_color,
      mid = mid_color,
      midpoint = 0,
      limits = c(-1, 1),
      name = "Pearson"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      ...
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = NULL, y = NULL)

  if (isTRUE(show_value)) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = label), size = value_text_size, color = "black")
  }

  if (isTRUE(show_value) && any(!is.na(df$ci_label))) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = ci_label, y = as.numeric(Var1) - 0.25),
      size = ci_text_size,
      color = "gray30",
      na.rm = TRUE
    )
  }

  p
}

#' @rdname pearson_corr
#' @method summary pearson_corr
#' @param object An object of class \code{pearson_corr}.
#' @param ci_digits Integer; digits for Pearson confidence limits in the
#'   pairwise summary.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @export
summary.pearson_corr <- function(object,
                                 n = NULL,
                                 topn = NULL,
                                 max_vars = NULL,
                                 width = NULL,
                                 ci_digits = 3,
                                 show_ci = NULL,
                                 ...) {
  check_inherits(object, "pearson_corr")
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )
  if (is.null(.mc_pearson_ci_attr(object))) {
    return(.mc_summary_corr_matrix(object, topn = topn))
  }
  .mc_pearson_pairwise_summary(
    object,
    ci_digits = ci_digits,
    show_ci = show_ci
  )
}

#' @rdname pearson_corr
#' @method print summary.pearson_corr
#' @param x An object of class \code{summary.pearson_corr}.
#' @export
print.summary.pearson_corr <- function(x, digits = NULL, n = NULL,
                                       topn = NULL, max_vars = NULL,
                                       width = NULL, show_ci = NULL, ...) {
  .mc_print_pairwise_summary_digest(
    x,
    title = "Pearson correlation summary",
    digits = .mc_coalesce(digits, 4),
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ci_method = "fisher_z",
    ...
  )
  invisible(x)
}

