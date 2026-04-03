#' @title Pairwise (or Two-Vector) Kendall's Tau Rank Correlation
#'
#' @description
#' Computes pairwise Kendall's tau correlations for numeric data using a
#' high-performance 'C++' backend. Optional confidence intervals are available
#' for matrix and data-frame input.
#'
#' @param data
#' For matrix/data frame mode, a numeric matrix or a data frame with at least
#' two numeric columns. All non-numeric columns are excluded. For two-vector
#' mode, a numeric vector \code{x}.
#' @param y Optional numeric vector \code{y} of the same length as \code{data}
#' when \code{data} is a vector. If supplied, the function computes the
#' Kendall correlation \emph{between \code{data} and \code{y}} using a
#' low-overhead scalar path and returns a single number.
#' @param check_na Logical (default \code{TRUE}). If \code{TRUE}, inputs must
#' be free of missing/undefined values. Use \code{FALSE} only when missingness
#' has already been handled upstream.
#' @param ci Logical (default \code{FALSE}). If \code{TRUE}, attach pairwise
#' confidence intervals for the off-diagonal Kendall correlations in
#' matrix/data-frame mode.
#' @param conf_level Confidence level used when \code{ci = TRUE}. Default is
#' \code{0.95}.
#' @param ci_method Confidence-interval engine used when \code{ci = TRUE}.
#' Supported Kendall methods are \code{"fieller"} (default),
#' \code{"brown_benedetti"}, and \code{"if_el"}.
#'
#' @return
#' \itemize{
#'   \item If \code{y} is \code{NULL} and \code{data} is a matrix/data frame: a
#'   symmetric numeric matrix where entry \code{(i, j)} is the Kendall's tau
#'   correlation between the \code{i}-th and \code{j}-th numeric columns. When
#'   \code{ci = TRUE}, the object also carries a \code{ci} attribute with
#'   elements \code{est}, \code{lwr.ci}, \code{upr.ci}, \code{conf.level}, and
#'   \code{ci.method}. Pairwise complete-case sample sizes are stored in
#'   \code{attr(x, "diagnostics")$n_complete}.
#'   \item If \code{y} is provided (two-vector mode): a single numeric scalar,
#'   the Kendall's tau correlation between \code{data} and \code{y}.
#' }
#'
#' @details
#' Kendall's tau is a rank-based measure of association between two variables.
#' For a dataset with \eqn{n} observations on variables \eqn{X} and \eqn{Y},
#' let \eqn{n_0 = n(n - 1)/2} be the number of unordered pairs, \eqn{C} the
#' number of concordant pairs, and \eqn{D} the number of discordant pairs.
#' Let \eqn{T_x = \sum_g t_g (t_g - 1)/2} and \eqn{T_y = \sum_h u_h (u_h - 1)/2}
#' be the numbers of tied pairs within \eqn{X} and within \eqn{Y}, respectively,
#' where \eqn{t_g} and \eqn{u_h} are tie-group sizes in \eqn{X} and \eqn{Y}.
#'
#' The tie-robust Kendall's tau-b is:
#' \deqn{ \tau_b = \frac{C - D}{\sqrt{(n_0 - T_x)\,(n_0 - T_y)}}. }
#' When there are no ties (\eqn{T_x = T_y = 0}), this reduces to tau-a:
#' \deqn{ \tau_a = \frac{C - D}{n(n-1)/2}. }
#'
#' The function automatically handles ties. In degenerate cases where a
#' variable is constant (\eqn{n_0 = T_x} or \eqn{n_0 = T_y}), the tau-b
#' denominator is zero and the correlation is undefined (returned as \code{NA}
#' off the diagonal).
#'
#' When \code{check_na = FALSE}, each \eqn{(i,j)} estimate is recomputed on the
#' pairwise complete-case overlap of columns \eqn{i} and \eqn{j}. Confidence
#' intervals use the observed pairwise-complete Kendall estimate and the same
#' pairwise complete-case overlap.
#'
#' With \code{ci_method = "fieller"}, the interval is built on the Fisher-style
#' transformed scale \eqn{z = \operatorname{atanh}(\hat\tau)} using Fieller's
#' asymptotic standard error
#' \deqn{ \operatorname{SE}(z) = \sqrt{\frac{0.437}{n - 4}}, }
#' where \eqn{n} is the pairwise complete-case sample size. The interval is then
#' mapped back with \code{tanh()} and clipped to \eqn{[-1, 1]} for numerical
#' safety. This is the default Kendall CI and is intended to be the fast,
#' production-oriented choice.
#'
#' With \code{ci_method = "brown_benedetti"}, the interval is computed on the
#' Kendall tau scale using the Brown-Benedetti large-sample variance for
#' Kendall's tau-b. This path is tie-aware, remains on the original Kendall
#' scale, and is intended as a conventional asymptotic alternative when a
#' direct tau-scale interval is preferred.
#'
#' With \code{ci_method = "if_el"}, the interval is computed in 'C++' using an
#' influence-function empirical-likelihood construction built from the
#' linearised Kendall estimating equation. The lower and upper limits are found
#' by solving the empirical-likelihood ratio equation against the
#' \eqn{\chi^2_1}-cutoff implied by \code{conf_level}. This method is slower
#' than \code{"fieller"} and is intended for specialised inference.
#'
#' \strong{Performance:}
#' \itemize{
#'   \item In the \strong{two-vector mode} (\code{y} supplied), the C++ backend uses a
#'   raw-double path with minimal overhead.
#'   \item In the \strong{matrix/data-frame mode}, the no-missing estimate-only path
#'   uses the Knight (1966) \eqn{O(n \log n)} algorithm. Pairwise-complete
#'   inference paths recompute each pair on its complete-case overlap; the
#'   \code{"brown_benedetti"} interval adds tie-aware large-sample variance
#'   calculations and \code{"if_el"} adds extra per-pair likelihood solving.
#' }
#'
#' @note Missing values are not allowed when \code{check_na = TRUE}. Columns
#' with fewer than two observations are excluded. Confidence intervals are not
#' available in the two-vector interface.
#'
#' @references
#' Kendall, M. G. (1938). A New Measure of Rank Correlation. \emph{Biometrika},
#' 30(1/2), 81-93.
#'
#' Knight, W. R. (1966). A Computer Method for Calculating Kendall's Tau with
#' Ungrouped Data. \emph{Journal of the American Statistical Association},
#' 61(314), 436-439.
#'
#' Fieller, E. C., Hartley, H. O., & Pearson, E. S. (1957). Tests for rank
#' correlation coefficients. I. \emph{Biometrika}, 44(3/4), 470-481.
#'
#' Brown, M. B., & Benedetti, J. K. (1977). Sampling behavior of tests for
#' correlation in two-way contingency tables. \emph{Journal of the American
#' Statistical Association}, 72(358), 309-315.
#'
#' Huang, Z., & Qin, G. (2023). Influence function-based confidence intervals
#' for the Kendall rank correlation coefficient. \emph{Computational
#' Statistics}, 38(2), 1041-1055.
#'
#' Croux, C., & Dehon, C. (2010). Influence functions of the Spearman and
#' Kendall correlation measures. \emph{Statistical Methods & Applications},
#' 19, 497-515.
#'
#' @examples
#' # Basic usage with a matrix
#' mat <- cbind(a = rnorm(100), b = rnorm(100), c = rnorm(100))
#' kt <- kendall_tau(mat)
#' print(kt)
#' summary(kt)
#' plot(kt)
#'
#' # Confidence intervals
#' kt_ci <- kendall_tau(mat[, 1:3], ci = TRUE)
#' print(kt_ci, show_ci = "yes")
#' summary(kt_ci)
#'
#' # Two-vector mode (scalar path)
#' x <- rnorm(1000); y <- 0.5 * x + rnorm(1000)
#' kendall_tau(x, y)
#'
#' # Including ties
#' tied_df <- data.frame(
#'   v1 = rep(1:5, each = 20),
#'   v2 = rep(5:1, each = 20),
#'   v3 = rnorm(100)
#' )
#' kt_tied <- kendall_tau(tied_df, ci = TRUE, ci_method = "fieller")
#' print(kt_tied, show_ci = "yes")
#'
#' # Interactive viewing (requires shiny)
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   view_corr_shiny(kt)
#' }
#'
#' @seealso \code{\link{print.kendall_matrix}}, \code{\link{plot.kendall_matrix}}
#' @author Thiago de Paula Oliveira
#' @export
kendall_tau <- function(data, y = NULL, check_na = TRUE, ci = FALSE,
                        conf_level = 0.95,
                        ci_method = c("fieller", "if_el", "brown_benedetti")) {
  check_bool(ci, arg = "ci")
  ci_method <- match.arg(ci_method)
  if (isTRUE(ci)) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  }

  if (!is.null(y)) {
    if (isTRUE(ci)) {
      abort_bad_arg("ci",
        message = "Confidence intervals are not available when {.arg y} is supplied."
      )
    }
    if (!is.numeric(data) || !is.numeric(y)) {
      abort_bad_arg("data",
        message = "and {.arg y} must be numeric vectors for two-vector mode."
      )
    }
    check_same_length(data, y, arg_x = "data", arg_y = "y")
    if (check_na && (any(!is.finite(data)) || any(!is.finite(y)))) {
      abort_bad_arg("data",
        message = "and {.arg y} must be free of NA/NaN/Inf when {.arg check_na} = TRUE.",
        .hint   = "Set `check_na = FALSE` only if missingness has been handled upstream."
      )
    }

    tau <- kendall_tau2_cpp(as.numeric(data), as.numeric(y))
    return(as.numeric(tau))
  }

  numeric_data <- validate_corr_input(data, check_na = check_na)
  colnames_data <- colnames(numeric_data)
  diagnostics <- NULL
  ci_attr <- NULL

  if (isTRUE(check_na) && !isTRUE(ci)) {
    result <- kendall_matrix_cpp(numeric_data)
  } else {
    pairwise <- kendall_matrix_pairwise_cpp(
      numeric_data,
      return_ci = ci,
      conf_level = conf_level,
      ci_method = ci_method
    )
    result <- pairwise$est
    diagnostics <- list(n_complete = pairwise$n_complete)
    dimnames(diagnostics$n_complete) <- list(colnames_data, colnames_data)
    if (isTRUE(ci)) {
      ci_attr <- list(
        est = unclass(result),
        lwr.ci = unclass(pairwise$lwr),
        upr.ci = unclass(pairwise$upr),
        conf.level = pairwise$conf_level,
        ci.method = pairwise$ci_method
      )
      dimnames(ci_attr$est) <- list(colnames_data, colnames_data)
      dimnames(ci_attr$lwr.ci) <- list(colnames_data, colnames_data)
      dimnames(ci_attr$upr.ci) <- list(colnames_data, colnames_data)
    }
  }

  colnames(result) <- rownames(result) <- colnames_data
  out <- .mc_structure_corr_matrix(
    result,
    class_name = "kendall_matrix",
    method = "kendall",
    description = "Pairwise Kendall's tau (auto tau-a/tau-b) correlation matrix",
    diagnostics = diagnostics
  )
  if (!is.null(ci_attr)) {
    attr(out, "ci") <- ci_attr
    attr(out, "conf.level") <- conf_level
    attr(out, "ci.method") <- ci_method
  }
  out
}

.mc_kendall_ci_attr <- function(x) {
  attr(x, "ci", exact = TRUE)
}

.mc_kendall_pairwise_summary <- function(object,
                                         digits = 4,
                                         ci_digits = 3,
                                         show_ci = NULL) {
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )
  check_inherits(object, "kendall_matrix")

  est <- as.matrix(object)
  rn <- rownames(est); cn <- colnames(est)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(est)))
  if (is.null(cn)) cn <- as.character(seq_len(ncol(est)))

  ci <- .mc_kendall_ci_attr(object)
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

  out <- structure(df, class = c("summary.kendall_matrix", "data.frame"))
  attr(out, "overview") <- .mc_summary_corr_matrix(object)
  attr(out, "has_ci") <- include_ci
  attr(out, "conf.level") <- if (is.null(ci)) NA_real_ else ci$conf.level
  attr(out, "ci.method") <- if (is.null(ci)) NA_character_ else ci$ci.method
  attr(out, "digits") <- digits
  attr(out, "ci_digits") <- ci_digits
  out
}

#' @rdname kendall_tau
#' @method print kendall_matrix
#' @title Print Method for \code{kendall_matrix} Objects
#'
#' @param x An object of class \code{kendall_matrix}.
#' @param digits Integer; number of decimal places to print.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param ci_digits Integer; digits for Kendall confidence limits.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @return Invisibly returns the \code{kendall_matrix} object.
#' @export
print.kendall_matrix <- function(x, digits = 4, n = NULL,
                                 topn = NULL, max_vars = NULL,
                                 width = NULL, ci_digits = 3,
                                 show_ci = NULL, ...) {
  .mc_print_corr_matrix(
    x,
    header = "Kendall correlation matrix",
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

#' @rdname kendall_tau
#' @method plot kendall_matrix
#' @title Plot Method for \code{kendall_matrix} Objects
#'
#' @param x An object of class \code{kendall_matrix}.
#' @param title Plot title. Default is \code{"Kendall's Tau correlation
#' heatmap"}.
#' @param low_color Color for the minimum tau value. Default is
#' \code{"indianred1"}.
#' @param high_color Color for the maximum tau value. Default is
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
#' @export
plot.kendall_matrix <- function(x, title = "Kendall's Tau correlation heatmap",
                                low_color = "indianred1", high_color = "steelblue1",
                                mid_color = "white", value_text_size = 4,
                                ci_text_size = 3, show_value = TRUE, ...) {
  check_bool(show_value, arg = "show_value")
  ci <- .mc_kendall_ci_attr(x)
  if (is.null(ci) || is.null(ci$lwr.ci) || is.null(ci$upr.ci)) {
    return(.mc_plot_corr_matrix(
      x, class_name = "kendall_matrix", fill_name = "Tau",
      title = title, low_color = low_color, high_color = high_color,
      mid_color = mid_color, value_text_size = value_text_size,
      show_value = show_value, ...
    ))
  }

  est_mat <- as.matrix(x)
  df_est <- as.data.frame(as.table(est_mat))
  names(df_est) <- c("Var1", "Var2", "tau")

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
  df$label <- sprintf("%.2f", df$tau)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Var2, y = Var1, fill = .data$tau)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(
      low = low_color,
      high = high_color,
      mid = mid_color,
      midpoint = 0,
      limits = c(-1, 1),
      name = "Tau"
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

#' @rdname kendall_tau
#' @method summary kendall_matrix
#' @param object An object of class \code{kendall_matrix}.
#' @param ci_digits Integer; digits for Kendall confidence limits in the
#'   pairwise summary.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @export
summary.kendall_matrix <- function(object,
                                   n = NULL,
                                   topn = NULL,
                                   max_vars = NULL,
                                   width = NULL,
                                   ci_digits = 3,
                                   show_ci = NULL,
                                   ...) {
  check_inherits(object, "kendall_matrix")
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )
  if (is.null(.mc_kendall_ci_attr(object))) {
    return(.mc_summary_corr_matrix(object, topn = topn))
  }
  .mc_kendall_pairwise_summary(
    object,
    ci_digits = ci_digits,
    show_ci = show_ci
  )
}

#' @rdname kendall_tau
#' @method print summary.kendall_matrix
#' @param x An object of class \code{summary.kendall_matrix}.
#' @export
print.summary.kendall_matrix <- function(x, digits = NULL, n = NULL,
                                         topn = NULL, max_vars = NULL,
                                         width = NULL, show_ci = NULL, ...) {
  ci_method <- attr(x, "ci.method", exact = TRUE)
  .mc_print_pairwise_summary_digest(
    x,
    title = "Kendall correlation summary",
    digits = .mc_coalesce(digits, 4),
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ci_method = ci_method,
    ...
  )
  invisible(x)
}
