# Internal exact percentage-bend helpers. These stay in R until the compiled
# backend reproduces the established wide-matrix reference behaviour exactly.
.mc_pbcor_scores <- function(z, beta = 0.2) {
  med <- stats::median(z)
  omega <- sort(abs(z - med))[floor((1 - beta) * length(z))]
  if (!is.finite(omega) || omega <= 0) return(rep(NA_real_, length(z)))

  psi <- (z - med) / omega
  i1 <- sum(psi < -1)
  i2 <- sum(psi > 1)
  sx <- z
  sx[psi < -1 | psi > 1] <- 0
  denom <- length(z) - i1 - i2
  if (denom <= 0) return(rep(NA_real_, length(z)))

  theta <- (sum(sx) + omega * (i2 - i1)) / denom
  a <- (z - theta) / omega
  a[a <= -1] <- -1
  a[a >= 1] <- 1
  a
}

.mc_pbcor_pair_exact <- function(x, y, beta = 0.2) {
  ax <- .mc_pbcor_scores(x, beta = beta)
  ay <- .mc_pbcor_scores(y, beta = beta)
  if (anyNA(ax) || anyNA(ay)) return(NA_real_)

  den_x <- sum(ax^2)
  den_y <- sum(ay^2)
  if (!is.finite(den_x) || !is.finite(den_y) || den_x <= 0 || den_y <= 0) {
    return(NA_real_)
  }

  sum(ax * ay) / sqrt(den_x * den_y)
}

.mc_pbcor_matrix_exact <- function(X, beta = 0.2) {
  n <- nrow(X)
  p <- ncol(X)
  if (p < 2L) stop("At least two numeric columns are required.", call. = FALSE)
  scores <- matrix(0, nrow = n, ncol = p)
  valid <- rep(FALSE, p)

  for (j in seq_len(p)) {
    a <- .mc_pbcor_scores(X[, j], beta = beta)
    if (!anyNA(a)) {
      ss <- sum(a^2)
      if (is.finite(ss) && ss > 0) {
        scores[, j] <- a
        valid[j] <- TRUE
      }
    }
  }

  ss <- colSums(scores^2)
  out <- as.matrix(crossprod(scores))
  denom <- sqrt(outer(ss, ss))
  out <- as.matrix(out / denom)
  out <- pmin(1, pmax(-1, out))
  out <- matrix(out, nrow = p, ncol = p, dimnames = list(colnames(X), colnames(X)))

  if (any(!valid)) {
    bad <- which(!valid)
    out[bad, ] <- NA_real_
    out[, bad] <- NA_real_
  }
  diag(out) <- 1
  out
}

.mc_pbcor_matrix_pairwise_exact <- function(X, beta = 0.2, min_n = 5L) {
  p <- ncol(X)
  if (p < 2L) stop("At least two numeric columns are required.", call. = FALSE)
  out <- matrix(NA_real_, nrow = p, ncol = p)
  dimnames(out) <- list(colnames(X), colnames(X))

  finite_n <- colSums(is.finite(X))
  diag(out) <- 1

  for (j in seq_len(p - 1L)) {
    for (k in seq.int(j + 1L, p)) {
      ok <- is.finite(X[, j]) & is.finite(X[, k])
      if (sum(ok) < min_n) next
      out[j, k] <- .mc_pbcor_pair_exact(X[ok, j], X[ok, k], beta = beta)
      out[k, j] <- out[j, k]
    }
  }

  if (any(finite_n < 2L)) {
    bad <- which(finite_n < 2L)
    out[bad, ] <- NA_real_
    out[, bad] <- NA_real_
    diag(out) <- 1
  }
  out
}

.mc_pbcor_pair_inference <- function(x,
                                     y,
                                     beta = 0.2,
                                     ci = FALSE,
                                     p_value = FALSE,
                                     conf_level = 0.95,
                                     n_boot = 500L,
                                     seed = NULL) {
  est <- .mc_pbcor_pair_exact(x, y, beta = beta)
  n_complete <- length(x)
  statistic <- NA_real_
  p_val <- NA_real_
  parameter <- NA_real_
  lwr <- NA_real_
  upr <- NA_real_

  if (isTRUE(p_value) && is.finite(est) && n_complete > 2L) {
    parameter <- n_complete - 2
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
          .mc_pbcor_pair_exact(x[idx], y[idx], beta = beta)
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

.mc_pbcor_ci_attr <- function(x) {
  attr(x, "ci", exact = TRUE)
}

.mc_pbcor_inference_attr <- function(x) {
  attr(x, "inference", exact = TRUE)
}

.mc_pbcor_pairwise_payload <- function(X,
                                       est,
                                       beta = 0.2,
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

      info <- .mc_pbcor_pair_inference(
        X[ok, i],
        X[ok, j],
        beta = beta,
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
      method = "percentage_bend_t_test",
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

.mc_pbcor_pairwise_summary <- function(object,
                                       digits = 4,
                                       ci_digits = 3,
                                       p_digits = 4,
                                       show_ci = NULL) {
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )
  check_inherits(object, "pbcor")

  est <- as.matrix(object)
  rn <- rownames(est); cn <- colnames(est)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(est)))
  if (is.null(cn)) cn <- as.character(seq_len(ncol(est)))

  ci <- .mc_pbcor_ci_attr(object)
  inf <- .mc_pbcor_inference_attr(object)
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

  out <- .mc_finalize_summary_df(df, class_name = "summary.pbcor")
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

#' Percentage bend correlation
#'
#' @description
#' Computes all pairwise percentage bend correlations for the numeric columns of
#' a matrix or data frame. Percentage bend correlation limits the influence of
#' extreme marginal observations by bending standardised deviations into the
#' interval \eqn{[-1, 1]}, yielding a Pearson-like measure that is robust to
#' outliers and heavy tails.
#'
#' @param data A numeric matrix or data frame containing numeric columns.
#' @param beta Bending constant in \code{[0, 0.5)} that sets the cutoff used to
#'   bend standardised deviations toward the interval \eqn{[-1, 1]}. Larger
#'   values cause more observations to be bent and increase resistance to
#'   marginal outliers. Default \code{0.2}. See Details.
#' @param na_method One of \code{"error"} (default) or \code{"pairwise"}.
#'   With \code{"pairwise"}, each correlation is computed on the overlapping
#'   complete rows for the column pair.
#' @param n_threads Integer \eqn{\geq 1}. Number of OpenMP threads used for the
#'   point-estimate matrix computation. Defaults to
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
#' @param x An object of class \code{pbcor}.
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
#' @return A symmetric correlation matrix with class \code{pbcor} and
#'   attributes \code{method = "percentage_bend_correlation"},
#'   \code{description}, and \code{package = "matrixCorr"}. When
#'   \code{ci = TRUE}, the returned object also carries a \code{ci} attribute
#'   with elements \code{est}, \code{lwr.ci}, \code{upr.ci},
#'   \code{conf.level}, and \code{ci.method}, plus
#'   \code{attr(x, "conf.level")}. When \code{p_value = TRUE}, it also carries
#'   an \code{inference} attribute with elements \code{estimate},
#'   \code{statistic}, \code{parameter}, \code{p_value}, \code{n_obs}, and
#'   \code{alternative}. When either inferential option is requested, the
#'   object also carries \code{diagnostics$n_complete}.
#'
#' @details
#' Let \eqn{X \in \mathbb{R}^{n \times p}} be a numeric matrix with rows as
#' observations and columns as variables. For a column
#' \eqn{x = (x_i)_{i=1}^n}, let \eqn{m_x = \mathrm{med}(x)} and define
#' \eqn{\omega_\beta(x)} as the \eqn{\lfloor (1-\beta)n \rfloor}-th order
#' statistic of \eqn{|x_i - m_x|}. Larger values of \code{beta} reduce
#' \eqn{\omega_\beta(x)}, so more observations are bent to the bounds
#' \eqn{-1} and \eqn{1}.
#'
#' The one-step percentage-bend location is
#' \deqn{
#' \hat\theta_{pb}(x) =
#' \frac{\sum_{i: |\psi_i| \le 1} x_i + \omega_\beta(x)(i_2 - i_1)}
#'      {n - i_1 - i_2},
#' \qquad
#' \psi_i = \frac{x_i - m_x}{\omega_\beta(x)},
#' }
#' where \eqn{i_1 = \sum_{i=1}^n \mathbf{1}(\psi_i < -1)} and
#' \eqn{i_2 = \sum_{i=1}^n \mathbf{1}(\psi_i > 1)}. The bent scores are
#' \deqn{
#' a_i = \max\!\left\{-1,\; \min\!\left(1,\frac{x_i - \hat\theta_{pb}(x)}
#' {\omega_\beta(x)}\right)\right\},
#' }
#' and likewise \eqn{b_i} for a second column \eqn{y}. The percentage bend
#' correlation for the pair \eqn{(x,y)} is
#' \deqn{
#' r_{pb}(x,y) =
#' \frac{\sum_{i=1}^n a_i b_i}
#'      {\sqrt{\sum_{i=1}^n a_i^2}\sqrt{\sum_{i=1}^n b_i^2}}.
#' }
#'
#' In the complete-data path, the bent score vectors are computed once per
#' column and collected into a matrix \eqn{A = [a_{\cdot 1}, \ldots, a_{\cdot p}]},
#' after which the correlation matrix is formed from their cross-products:
#' \deqn{
#' R_{pb} = D_A^{-1/2} A^\top A D_A^{-1/2},
#' \qquad
#' D_A = \mathrm{diag}(a_{\cdot 1}^\top a_{\cdot 1}, \ldots,
#' a_{\cdot p}^\top a_{\cdot p}).
#' }
#' If a column yields an undefined bent-score denominator, the corresponding row
#' and column are returned as \code{NA}. With \code{na_method = "pairwise"},
#' each pair is recomputed on its complete-case overlap. As with pairwise
#' Pearson correlation, this pairwise path can break positive semidefiniteness.
#'
#' When \code{p_value = TRUE}, the method-specific test statistic for a pairwise
#' estimate \eqn{r_{pb}} based on \eqn{n_{ij}} complete observations is
#' \deqn{
#' T_{ij} = r_{pb,ij}\sqrt{\frac{n_{ij} - 2}{1 - r_{pb,ij}^2}},
#' }
#' and the reported p-value is the two-sided Student-\eqn{t} tail probability
#' with \eqn{n_{ij}-2} degrees of freedom. When \code{ci = TRUE}, the interval
#' is a percentile bootstrap interval based on \eqn{n_{\mathrm{boot}}}
#' resamples drawn from the pairwise complete cases. If
#' \eqn{\tilde r_{pb,(1)} \le \cdots \le \tilde r_{pb,(B)}} denotes the sorted
#' bootstrap sample of finite estimates with \eqn{B} retained resamples, the
#' reported limits are
#' \deqn{
#' \tilde r_{pb,(\ell)} \quad \text{and} \quad \tilde r_{pb,(u)},
#' }
#' where \eqn{\ell = \lfloor (\alpha/2) B + 0.5 \rfloor} and
#' \eqn{u = \lfloor (1-\alpha/2) B + 0.5 \rfloor} for
#' \eqn{\alpha = 1 - \mathrm{conf\_level}}. Resamples that yield undefined
#' estimates are discarded before the percentile limits are formed.
#'
#' \strong{Computational complexity.} In the complete-data path, forming the
#' bent scores requires sorting within each column and the cross-product step
#' costs \eqn{O(n p^2)} with \eqn{O(p^2)} output storage. When
#' \code{ci = TRUE}, the bootstrap cost is incurred separately for each column
#' pair.
#'
#' @references
#' Wilcox, R. R. (1994). The percentage bend correlation coefficient.
#' Psychometrika, 59(4), 601-616. \doi{10.1007/BF02294395}
#'
#' @seealso [wincor()], [skipped_corr()], [bicor()]
#'
#' @examples
#' set.seed(10)
#' X <- matrix(rnorm(150 * 4), ncol = 4)
#' X[sample(length(X), 8)] <- X[sample(length(X), 8)] + 10
#'
#' R <- pbcor(X)
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
pbcor <- function(data,
                  na_method = c("error", "pairwise"),
                  ci = FALSE,
                  p_value = FALSE,
                  conf_level = 0.95,
                  n_threads = getOption("matrixCorr.threads", 1L),
                  beta = 0.2,
                  n_boot = 500L,
                  seed = NULL) {
  na_method <- match.arg(na_method)
  check_scalar_numeric(beta,
                       arg = "beta",
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
  dn <- .mc_square_dimnames(colnames_data)

  res <- if (na_method == "error") {
    pbcor_matrix_cpp(numeric_data, beta = beta, n_threads = n_threads)
  } else {
    pbcor_matrix_pairwise_cpp(numeric_data, beta = beta, min_n = 5L, n_threads = n_threads)
  }
  res <- .mc_set_matrix_dimnames(res, colnames_data)

  payload <- NULL
  if (isTRUE(ci) || isTRUE(p_value)) {
    payload <- .mc_pbcor_pairwise_payload(
      numeric_data,
      est = res,
      beta = beta,
      ci = ci,
      p_value = p_value,
      conf_level = conf_level,
      n_boot = n_boot,
      seed = seed
    )
  }

  .mc_structure_corr_matrix(
    res,
    class_name = "pbcor",
    method = "percentage_bend_correlation",
    description = paste0(
      "Percentage bend correlation; beta = ", beta,
      "; NA mode = ", na_method, "."
    ),
    diagnostics = if (is.null(payload)) NULL else payload$diagnostics,
    extra_attrs = c(
      if (!is.null(payload$ci)) {
        list(
          ci = payload$ci,
          conf.level = conf_level,
          n_boot = n_boot
        )
      },
      if (!is.null(payload$inference)) {
        list(inference = payload$inference)
      }
    )
  )
}

#' @rdname pbcor
#' @method print pbcor
#' @export
print.pbcor <- function(x, digits = 4, n = NULL, topn = NULL,
                        max_vars = NULL, width = NULL,
                        show_ci = NULL, ...) {
  .mc_print_corr_matrix(
    x,
    header = "Percentage bend correlation matrix",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
}

#' @rdname pbcor
#' @method plot pbcor
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
#' @export
plot.pbcor <- function(x,
                       title = "Percentage bend correlation heatmap",
                       low_color = "indianred1",
                       high_color = "steelblue1",
                       mid_color = "white",
                       value_text_size = 4,
                       show_value = TRUE,
                       ...) {
  .mc_plot_corr_matrix(
    x,
    class_name = "pbcor",
    fill_name = "pbcor",
    title = title,
    low_color = low_color,
    high_color = high_color,
    mid_color = mid_color,
    value_text_size = value_text_size,
    show_value = show_value,
    ...
  )
}

#' @rdname pbcor
#' @method summary pbcor
#' @param object An object of class \code{pbcor}.
#' @export
summary.pbcor <- function(object, n = NULL, topn = NULL,
                          max_vars = NULL, width = NULL,
                          ci_digits = 3,
                          p_digits = 4,
                          show_ci = NULL, ...) {
  check_inherits(object, "pbcor")
  if (is.null(.mc_pbcor_ci_attr(object)) && is.null(.mc_pbcor_inference_attr(object))) {
    return(.mc_summary_corr_matrix(object, topn = topn))
  }
  .mc_pbcor_pairwise_summary(
    object,
    ci_digits = ci_digits,
    p_digits = p_digits,
    show_ci = show_ci
  )
}

#' @rdname pbcor
#' @method print summary.pbcor
#' @param x An object of class \code{summary.pbcor}.
#' @export
print.summary.pbcor <- function(x, digits = NULL, n = NULL,
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
    title = "Percentage bend correlation summary",
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
