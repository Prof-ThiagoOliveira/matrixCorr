#' @title Pairwise Winsorized correlation
#'
#' @description
#' Computes all pairwise Winsorized correlation coefficients for the numeric
#' columns of a matrix or data frame using a high-performance 'C++' backend.
#'
#' This function Winsorizes each margin at proportion \code{tr} and then
#' computes ordinary Pearson correlation on the Winsorized values. It is a
#' simple robust alternative to Pearson correlation when the main concern is
#' unusually large or small observations in the marginal distributions.
#'
#' @param data A numeric matrix or a data frame with at least two numeric
#' columns. All non-numeric columns will be excluded.
#' @param tr Winsorization proportion in \code{[0, 0.5)}. For a sample of size
#'   \eqn{n}, let \eqn{g = \lfloor tr \cdot n \rfloor}; the \eqn{g} smallest
#'   observations are set to the \eqn{(g+1)}-st order statistic and the
#'   \eqn{g} largest observations are set to the \eqn{(n-g)}-th order
#'   statistic. Default \code{0.2}.
#' @param na_method One of \code{"error"} (default) or \code{"pairwise"}.
#' @param n_threads Integer \eqn{\geq 1}. Number of OpenMP threads. Defaults to
#'   \code{getOption("matrixCorr.threads", 1L)}.
#' @param ci Logical (default \code{FALSE}). If \code{TRUE}, attach percentile
#'   bootstrap confidence intervals for each pairwise estimate.
#' @param p_value Logical (default \code{FALSE}). If \code{TRUE}, attach the
#'   method-specific large-sample test statistic and two-sided p-value for each
#'   pairwise estimate.
#' @param conf_level Confidence level used when \code{ci = TRUE}. Default
#'   \code{0.95}.
#' @param n_boot Integer \eqn{\geq 1}. Number of bootstrap resamples used when
#'   \code{ci = TRUE}. Default \code{500}.
#' @param seed Optional positive integer used to seed the bootstrap resampling
#'   when \code{ci = TRUE}. If \code{NULL}, the current random-number stream is
#'   used.
#' @param x An object of class \code{wincor}.
#' @param digits Integer; number of digits to print.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param ci_digits Integer; digits used for confidence limits in pairwise
#'   summaries.
#' @param p_digits Integer; digits used for p-values in pairwise summaries.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param ... Additional arguments passed to the underlying print or plot helper.
#' @param title Character; plot title.
#' @param low_color,high_color,mid_color Colors used in the heatmap.
#' @param value_text_size Numeric text size for overlaid cell values.
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
#'
#' @return A symmetric correlation matrix with class \code{wincor} and
#'   attributes \code{method = "winsorized_correlation"}, \code{description},
#'   and \code{package = "matrixCorr"}. When \code{ci = TRUE}, the returned
#'   object also carries a \code{ci} attribute with elements \code{est},
#'   \code{lwr.ci}, \code{upr.ci}, \code{conf.level}, and \code{ci.method},
#'   plus \code{attr(x, "conf.level")}. When \code{p_value = TRUE}, it also
#'   carries an \code{inference} attribute with elements \code{estimate},
#'   \code{statistic}, \code{parameter}, \code{p_value}, \code{n_obs}, and
#'   \code{alternative}. When either inferential option is requested, the
#'   object also carries \code{diagnostics$n_complete}.
#'
#' @details
#' Let \eqn{X \in \mathbb{R}^{n \times p}} be a numeric matrix with rows as
#' observations and columns as variables. For a column
#' \eqn{x = (x_i)_{i=1}^n}, write the order statistics as
#' \eqn{x_{(1)} \le \cdots \le x_{(n)}} and let
#' \eqn{g = \lfloor tr \cdot n \rfloor}. The Winsorized values can be written as
#' \deqn{
#' x_i^{(w)} \;=\; \max\!\bigl\{x_{(g+1)},\, \min(x_i, x_{(n-g)})\bigr\}.
#' }
#' For two columns \eqn{x} and \eqn{y}, the Winsorized correlation is the
#' ordinary Pearson correlation computed from \eqn{x^{(w)}} and \eqn{y^{(w)}}:
#' \deqn{
#' r_w(x,y) \;=\;
#' \frac{\sum_{i=1}^n (x_i^{(w)}-\bar x^{(w)})(y_i^{(w)}-\bar y^{(w)})}
#'      {\sqrt{\sum_{i=1}^n (x_i^{(w)}-\bar x^{(w)})^2}\;
#'       \sqrt{\sum_{i=1}^n (y_i^{(w)}-\bar y^{(w)})^2}}.
#' }
#'
#' In matrix form, let \eqn{X^{(w)}} contain the Winsorized columns and define
#' the centred, unit-norm columns
#' \deqn{
#' z_{\cdot j} =
#' \frac{x_{\cdot j}^{(w)} - \bar x_j^{(w)} \mathbf{1}}
#'      {\sqrt{\sum_{i=1}^n (x_{ij}^{(w)}-\bar x_j^{(w)})^2}},
#' \qquad j=1,\ldots,p.
#' }
#' If \eqn{Z = [z_{\cdot 1}, \ldots, z_{\cdot p}]}, then the Winsorized
#' correlation matrix is
#' \deqn{
#' R_w \;=\; Z^\top Z.
#' }
#'
#' Winsorization acts on each margin separately, so it guards against marginal
#' outliers and heavy tails but does not target unusual points in the joint
#' cloud. This implementation Winsorizes each column in 'C++', centres and
#' normalises it, and forms the complete-data matrix from cross-products. With
#' \code{na_method = "pairwise"}, each pair is recomputed on its overlap of
#' non-missing rows. As with Pearson correlation, the complete-data path yields
#' a symmetric positive semidefinite matrix, whereas pairwise deletion can
#' break positive semidefiniteness. If the Winsorized variance of a column is
#' zero, correlations involving that column are returned as \code{NA}.
#'
#' When \code{p_value = TRUE}, inference follows the method-specific test based
#' on
#' \deqn{
#' T_{ij} = r_{w,ij}\sqrt{\frac{n_{ij} - 2}{1 - r_{w,ij}^2}},
#' }
#' evaluated against a \eqn{t}-distribution with
#' \eqn{n_{ij} - 2g_{ij} - 2} degrees of freedom, where
#' \eqn{g_{ij} = \lfloor tr \cdot n_{ij} \rfloor} and \eqn{n_{ij}} is the
#' pairwise complete-case sample size for the corresponding column pair. The
#' p-value is reported only when the pair is not identical and the resulting
#' degrees of freedom are positive. When \code{ci = TRUE}, the interval is a
#' percentile bootstrap interval based on \eqn{n_{\mathrm{boot}}} resamples
#' drawn from the pairwise complete cases. If
#' \eqn{\tilde r_{w,(1)} \le \cdots \le \tilde r_{w,(B)}} denotes the sorted
#' bootstrap sample of finite estimates with \eqn{B} retained resamples, the
#' reported limits are
#' \deqn{
#' \tilde r_{w,(\ell)} \quad \text{and} \quad \tilde r_{w,(u)},
#' }
#' where \eqn{\ell = \lfloor (\alpha/2) B + 0.5 \rfloor} and
#' \eqn{u = \lfloor (1-\alpha/2) B + 0.5 \rfloor} for
#' \eqn{\alpha = 1 - \mathrm{conf\_level}}. Resamples that yield undefined
#' estimates are discarded before the percentile limits are formed.
#'
#' \strong{Computational complexity.} In the complete-data path, Winsorizing the
#' columns requires sorting within each column, and forming the cross-product
#' matrix costs \eqn{O(n p^2)} with \eqn{O(p^2)} output storage. When
#' \code{ci = TRUE}, the bootstrap cost is incurred separately for each column
#' pair.
#'
#' @references
#' Wilcox, R. R. (1993). Some results on a Winsorized correlation coefficient.
#' British Journal of Mathematical and Statistical Psychology, 46(2), 339-349.
#' \doi{10.1111/j.2044-8317.1993.tb01020.x}
#'
#' Wilcox, R. R. (2012). Introduction to Robust Estimation and Hypothesis
#' Testing (3rd ed.). Academic Press.
#'
#' @seealso [pbcor()], [skipped_corr()], [bicor()]
#'
#' @examples
#' set.seed(11)
#' X <- matrix(rnorm(180 * 4), ncol = 4)
#' X[sample(length(X), 6)] <- X[sample(length(X), 6)] - 12
#'
#' R <- wincor(X, tr = 0.2)
#' print(R, digits = 2)
#' summary(R)
#' plot(R)
#'
#' # Interactive viewing (requires shiny)
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   view_corr_shiny(R)
#' }
#'
#' @author Thiago de Paula Oliveira
#' @export
wincor <- function(data,
                   na_method = c("error", "pairwise"),
                   ci = FALSE,
                   p_value = FALSE,
                   conf_level = 0.95,
                   n_threads = getOption("matrixCorr.threads", 1L),
                   tr = 0.2,
                   n_boot = 500L,
                   seed = NULL) {
  na_method <- match.arg(na_method)
  check_scalar_numeric(tr,
                       arg = "tr",
                       lower = 0,
                       upper = 0.5,
                       closed_lower = TRUE,
                       closed_upper = FALSE)
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")
  check_bool(ci, arg = "ci")
  check_bool(p_value, arg = "p_value")
  check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  n_boot <- check_scalar_int_pos(n_boot, arg = "n_boot")
  if (!is.null(seed)) {
    seed <- check_scalar_int_pos(seed, arg = "seed")
  }

  numeric_data <- if (na_method == "error") {
    validate_corr_input(data)
  } else {
    validate_corr_input(data, check_na = FALSE)
  }
  colnames_data <- colnames(numeric_data)

  res <- if (na_method == "error") {
    wincor_matrix_cpp(numeric_data, tr = tr, n_threads = n_threads)
  } else {
    wincor_matrix_pairwise_cpp(numeric_data, tr = tr, min_n = 5L, n_threads = n_threads)
  }

  colnames(res) <- rownames(res) <- colnames_data
  payload <- NULL
  if (isTRUE(ci) || isTRUE(p_value)) {
    payload <- .mc_wincor_pairwise_payload(
      numeric_data,
      est = res,
      tr = tr,
      ci = ci,
      p_value = p_value,
      conf_level = conf_level,
      n_boot = n_boot,
      seed = seed
    )
  }

  out <- .mc_structure_corr_matrix(
    res,
    class_name = "wincor",
    method = "winsorized_correlation",
    description = paste0(
      "Winsorized correlation; tr = ", tr,
      "; NA mode = ", na_method, "."
    ),
    diagnostics = if (is.null(payload)) NULL else payload$diagnostics
  )
  if (!is.null(payload$ci)) {
    attr(out, "ci") <- payload$ci
    attr(out, "conf.level") <- conf_level
    attr(out, "n_boot") <- n_boot
  }
  if (!is.null(payload$inference)) {
    attr(out, "inference") <- payload$inference
  }
  out
}

.mc_winsorize_vec <- function(z, tr = 0.2) {
  n <- length(z)
  if (!n) {
    return(z)
  }
  g <- floor(tr * n)
  if (g <= 0L) {
    return(z)
  }

  zs <- sort(z)
  low <- zs[[g + 1L]]
  high <- zs[[n - g]]
  pmax(low, pmin(z, high))
}

.mc_wincor_pair_exact <- function(x, y, tr = 0.2) {
  xw <- .mc_winsorize_vec(x, tr = tr)
  yw <- .mc_winsorize_vec(y, tr = tr)

  xc <- xw - mean(xw)
  yc <- yw - mean(yw)
  den_x <- sum(xc^2)
  den_y <- sum(yc^2)
  if (!is.finite(den_x) || !is.finite(den_y) || den_x <= 0 || den_y <= 0) {
    return(NA_real_)
  }

  sum(xc * yc) / sqrt(den_x * den_y)
}

.mc_wincor_pair_inference <- function(x,
                                      y,
                                      tr = 0.2,
                                      ci = FALSE,
                                      p_value = FALSE,
                                      conf_level = 0.95,
                                      n_boot = 500L,
                                      seed = NULL) {
  est <- .mc_wincor_pair_exact(x, y, tr = tr)
  n_complete <- length(x)
  g <- floor(tr * n_complete)
  parameter <- n_complete - 2 * g - 2
  statistic <- NA_real_
  p_val <- NA_real_
  lwr <- NA_real_
  upr <- NA_real_

  if (isTRUE(p_value) && is.finite(est) && parameter > 0L && any(x != y)) {
    if (abs(est) >= 1) {
      statistic <- sign(est) * Inf
      p_val <- 0
    } else {
      statistic <- est * sqrt((n_complete - 2) / (1 - est^2))
      p_val <- 2 * stats::pt(abs(statistic), df = parameter, lower.tail = FALSE)
    }
  }

  if (isTRUE(ci) && is.finite(est) && n_complete > 1L) {
    boot_vals <- .mc_eval_with_seed(
      seed,
      replicate(
        n_boot,
        {
          idx <- sample.int(n_complete, size = n_complete, replace = TRUE)
          .mc_wincor_pair_exact(x[idx], y[idx], tr = tr)
        }
      )
    )
    ci_vals <- .mc_percentile_boot_ci(boot_vals, conf_level = conf_level)
    lwr <- ci_vals[[1L]]
    upr <- ci_vals[[2L]]
  }

  list(
    estimate = est,
    statistic = statistic,
    p_value = p_val,
    parameter = parameter,
    lwr = lwr,
    upr = upr,
    n_complete = n_complete
  )
}

.mc_wincor_ci_attr <- function(x) {
  attr(x, "ci", exact = TRUE)
}

.mc_wincor_inference_attr <- function(x) {
  attr(x, "inference", exact = TRUE)
}

.mc_wincor_pairwise_payload <- function(X,
                                        est,
                                        tr = 0.2,
                                        ci = FALSE,
                                        p_value = FALSE,
                                        conf_level = 0.95,
                                        n_boot = 500L,
                                        seed = NULL) {
  est <- as.matrix(est)
  p <- ncol(est)
  dn <- dimnames(est)
  n_complete <- .mc_corr_n_complete(X)
  dimnames(n_complete) <- dn

  ci_lwr <- ci_upr <- NULL
  if (isTRUE(ci)) {
    ci_lwr <- matrix(NA_real_, p, p, dimnames = dn)
    ci_upr <- matrix(NA_real_, p, p, dimnames = dn)
    diag(ci_lwr) <- 1
    diag(ci_upr) <- 1
  }

  stat_mat <- p_mat <- df_mat <- NULL
  if (isTRUE(p_value)) {
    stat_mat <- matrix(NA_real_, p, p, dimnames = dn)
    p_mat <- matrix(NA_real_, p, p, dimnames = dn)
    df_mat <- matrix(NA_real_, p, p, dimnames = dn)
    diag(stat_mat) <- Inf
    diag(p_mat) <- 0
  }

  pair_id <- 0L
  for (i in seq_len(p - 1L)) {
    for (j in seq.int(i + 1L, p)) {
      pair_id <- pair_id + 1L
      ok <- is.finite(X[, i]) & is.finite(X[, j])
      n_ij <- sum(ok)
      n_complete[i, j] <- n_complete[j, i] <- n_ij
      if (n_ij < 5L) {
        next
      }

      info <- .mc_wincor_pair_inference(
        X[ok, i],
        X[ok, j],
        tr = tr,
        ci = ci,
        p_value = p_value,
        conf_level = conf_level,
        n_boot = n_boot,
        seed = .mc_seed_offset(seed, pair_id - 1L)
      )

      if (isTRUE(ci)) {
        ci_lwr[i, j] <- ci_lwr[j, i] <- info$lwr
        ci_upr[i, j] <- ci_upr[j, i] <- info$upr
      }
      if (isTRUE(p_value)) {
        stat_mat[i, j] <- stat_mat[j, i] <- info$statistic
        p_mat[i, j] <- p_mat[j, i] <- info$p_value
        df_mat[i, j] <- df_mat[j, i] <- info$parameter
      }
    }
  }

  diagnostics <- list(n_complete = n_complete)
  ci_attr <- NULL
  if (isTRUE(ci)) {
    ci_attr <- list(
      est = unclass(est),
      lwr.ci = ci_lwr,
      upr.ci = ci_upr,
      conf.level = conf_level,
      ci.method = "percentile_bootstrap"
    )
  }

  inference_attr <- NULL
  if (isTRUE(p_value)) {
    inference_attr <- list(
      method = "winsorized_t_test",
      estimate = unclass(est),
      statistic = stat_mat,
      p_value = p_mat,
      parameter = df_mat,
      n_obs = n_complete,
      alternative = "two.sided"
    )
  }

  list(
    diagnostics = diagnostics,
    ci = ci_attr,
    inference = inference_attr
  )
}

.mc_wincor_pairwise_summary <- function(object,
                                        digits = 4,
                                        ci_digits = 3,
                                        p_digits = 4,
                                        show_ci = NULL) {
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )
  check_inherits(object, "wincor")

  est <- as.matrix(object)
  rn <- rownames(est); cn <- colnames(est)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(est)))
  if (is.null(cn)) cn <- as.character(seq_len(ncol(est)))

  ci <- .mc_wincor_ci_attr(object)
  inf <- .mc_wincor_inference_attr(object)
  diag_attr <- attr(object, "diagnostics", exact = TRUE)
  include_ci <- identical(show_ci, "yes") && !is.null(ci)
  include_p <- !is.null(inf) && !is.null(inf$p_value)

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
        rec$lwr <- if (is.finite(ci$lwr.ci[i, j])) round(ci$lwr.ci[i, j], ci_digits) else NA_real_
        rec$upr <- if (is.finite(ci$upr.ci[i, j])) round(ci$upr.ci[i, j], ci_digits) else NA_real_
      }
      if (include_p) {
        rec$statistic <- if (is.finite(inf$statistic[i, j])) round(inf$statistic[i, j], digits) else NA_real_
        rec$p_value <- if (is.finite(inf$p_value[i, j])) round(inf$p_value[i, j], p_digits) else NA_real_
      }
      rows[[k]] <- rec
    }
  }

  df <- do.call(rbind.data.frame, rows)
  rownames(df) <- NULL
  num_cols <- intersect(c("estimate", "lwr", "upr", "statistic", "p_value"), names(df))
  int_cols <- intersect(c("n_complete"), names(df))
  for (nm in num_cols) df[[nm]] <- as.numeric(df[[nm]])
  for (nm in int_cols) df[[nm]] <- as.integer(df[[nm]])

  out <- .mc_finalize_summary_df(df, class_name = "summary.wincor")
  attr(out, "overview") <- .mc_summary_corr_matrix(object)
  attr(out, "has_ci") <- include_ci
  attr(out, "has_p") <- include_p
  attr(out, "conf.level") <- if (is.null(ci)) NA_real_ else ci$conf.level
  attr(out, "digits") <- digits
  attr(out, "ci_digits") <- ci_digits
  attr(out, "p_digits") <- p_digits
  attr(out, "inference_method") <- if (is.null(inf)) NA_character_ else inf$method
  attr(out, "n_boot") <- if (is.null(ci)) NA_integer_ else attr(object, "n_boot", exact = TRUE) %||% NA_integer_
  out
}

#' @rdname wincor
#' @method print wincor
#' @export
print.wincor <- function(x, digits = 4, n = NULL, topn = NULL,
                         max_vars = NULL, width = NULL,
                         show_ci = NULL, ...) {
  .mc_print_corr_matrix(
    x,
    header = "Winsorized correlation matrix",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
}

#' @rdname wincor
#' @method plot wincor
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
#' @export
plot.wincor <- function(x,
                        title = "Winsorized correlation heatmap",
                        low_color = "indianred1",
                        high_color = "steelblue1",
                        mid_color = "white",
                        value_text_size = 4,
                        show_value = TRUE,
                        ...) {
  .mc_plot_corr_matrix(
    x,
    class_name = "wincor",
    fill_name = "wincor",
    title = title,
    low_color = low_color,
    high_color = high_color,
    mid_color = mid_color,
    value_text_size = value_text_size,
    show_value = show_value,
    ...
  )
}

#' @rdname wincor
#' @method summary wincor
#' @param object An object of class \code{wincor}.
#' @export
summary.wincor <- function(object, n = NULL, topn = NULL,
                           max_vars = NULL, width = NULL,
                           ci_digits = 3,
                           p_digits = 4,
                           show_ci = NULL, ...) {
  check_inherits(object, "wincor")
  if (is.null(.mc_wincor_ci_attr(object)) && is.null(.mc_wincor_inference_attr(object))) {
    return(.mc_summary_corr_matrix(object, topn = topn))
  }
  .mc_wincor_pairwise_summary(
    object,
    ci_digits = ci_digits,
    p_digits = p_digits,
    show_ci = show_ci
  )
}

#' @rdname wincor
#' @method print summary.wincor
#' @param x An object of class \code{summary.wincor}.
#' @export
print.summary.wincor <- function(x, digits = NULL, n = NULL,
                                 topn = NULL, max_vars = NULL,
                                 width = NULL, show_ci = NULL, ...) {
  extra_items <- c(
    inference = attr(x, "inference_method", exact = TRUE),
    n_boot = if (isTRUE(attr(x, "has_ci", exact = TRUE))) {
      as.character(attr(x, "n_boot", exact = TRUE))
    } else {
      NA_character_
    }
  )
  .mc_print_pairwise_summary_digest(
    x,
    title = "Winsorized correlation summary",
    digits = .mc_coalesce(digits, 4),
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ci_method = "percentile_bootstrap",
    extra_items = extra_items,
    ...
  )
  invisible(x)
}
