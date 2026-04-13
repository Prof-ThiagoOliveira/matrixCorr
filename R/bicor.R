#' Biweight mid-correlation (bicor)
#'
#' @description
#' Computes pairwise biweight mid-correlations for numeric data. Bicor is a
#' robust, Pearson-like correlation that down-weights outliers and heavy-tailed
#' observations. Optional large-sample confidence intervals are available as a
#' derived feature.
#'
#'
#' @param data A numeric matrix or a data frame containing numeric columns.
#'   Factors, logicals and common time classes are dropped in the data-frame
#'   path. Missing values are not allowed unless \code{na_method = "pairwise"}.
#' @param c_const Positive numeric. Tukey biweight tuning constant applied to the
#'   \emph{raw} MAD; default \code{9} (Langfelder & Horvath’s convention).
#' @param max_p_outliers Numeric in \code{(0, 1]}. Optional cap on the maximum
#'   proportion of outliers \emph{on each side}; if \code{< 1}, side-specific
#'   rescaling maps those quantiles to \code{|u|=1}. Use \code{1} to disable.
#' @param pearson_fallback Character scalar indicating the fallback policy.
#'   One of:
#'   \itemize{
#'     \item \code{"hybrid"} (default): if a column has MAD = 0, that column
#'       uses Pearson standardisation, yielding a hybrid correlation.
#'     \item \code{"none"}: return \code{NA} if a column has MAD = 0 or becomes
#'       degenerate after weighting.
#'     \item \code{"all"}: force ordinary Pearson for all columns.
#'   }
#' @param na_method One of \code{"error"} (default, fastest) or \code{"pairwise"}.
#'   With \code{"pairwise"}, each \eqn{(j,k)} correlation is computed on the
#'   intersection of non-missing rows for the pair.
#' @param mad_consistent Logical; if \code{TRUE}, use the normal-consistent MAD
#'   (\code{MAD_raw * 1.4826}) in the bicor weights. Default \code{FALSE} to
#'   match Langfelder & Horvath (2012).
#' @param w Optional non-negative numeric vector of length \code{nrow(data)}
#'   giving \emph{row weights}. When supplied, weighted medians/MADs are used
#'   and Tukey weights are multiplied by \code{w} before normalisation.
#' @param sparse_threshold Optional numeric \eqn{\geq 0}. If supplied, sets
#'   entries with \code{|r| < sparse_threshold} to 0 and returns a sparse
#'   \code{"ddiMatrix"} (requires \pkg{Matrix}).
#' @param ci Logical (default \code{FALSE}). If \code{TRUE}, attach approximate
#'   Fisher-z confidence intervals for the off-diagonal biweight
#'   mid-correlations. This confidence interval is provided as an additional
#'   large-sample approximation.
#' @param conf_level Confidence level used when \code{ci = TRUE}. Default is
#'   \code{0.95}.
#' @param n_threads Integer \eqn{\geq 1}. Number of OpenMP threads. Defaults to
#'   \code{getOption("matrixCorr.threads", 1L)}.
#'
#' @return A symmetric correlation matrix with class \code{bicor}
#'   (or a \code{dgCMatrix} if \code{sparse_threshold} is used), with attributes:
#'   \code{method = "biweight_mid_correlation"}, \code{description},
#'   and \code{package = "matrixCorr"}. Downstream code should be prepared to
#'   handle either a dense numeric matrix or a sparse \code{dgCMatrix}. When
#'   \code{ci = TRUE}, the object also carries a \code{ci} attribute with
#'   elements \code{est}, \code{lwr.ci}, \code{upr.ci}, \code{conf.level}, and
#'   \code{ci.method}, together with an \code{inference} attribute containing
#'   the standard large-sample summary matrices \code{estimate},
#'   \code{statistic}, \code{p_value}, \code{Z}, and \code{n_obs}. Pairwise
#'   complete-case counts are stored in \code{attr(x, "diagnostics")$n_complete}.
#' Internally, all medians/MADs, Tukey weights, optional pairwise-NA handling,
#' and OpenMP loops are implemented in the C++ helpers
#' (\code{bicor_*_cpp()}), so the R wrapper mostly validates arguments and
#' dispatches to the appropriate backend.
#'
#' @details
#'
#' For a column \eqn{x = (x_a)_{a=1}^m}, let \eqn{\mathrm{med}(x)} be the median and
#' \eqn{\mathrm{MAD}(x) = \mathrm{med}(|x - \mathrm{med}(x)|)} the (raw) median
#' absolute deviation. If \code{mad_consistent = TRUE}, the consistent scale
#' \eqn{\mathrm{MAD}^\star(x) = 1.4826\,\mathrm{MAD}(x)} is used. With tuning constant
#' \eqn{c>0}, define
#' \deqn{u_a = \frac{x_a - \mathrm{med}(x)}{c\,\mathrm{MAD}^{(\star)}(x)}.}
#' The Tukey biweight gives per-observation weights
#' \deqn{w_a = (1 - u_a^2)^2\,\mathbf{1}\{|u_a| < 1\}.}
#' Robust standardisation of a column is
#' \deqn{\tilde x_a =
#' \frac{(x_a - \mathrm{med}(x))\,w_a}{
#'       \sqrt{\sum_{b=1}^m \big[(x_b - \mathrm{med}(x))\,w_b\big]^2}}.}
#' For two columns \eqn{x,y}, the biweight mid-correlation is
#' \deqn{\mathrm{bicor}(x,y) = \sum_{a=1}^m \tilde x_a\,\tilde y_a \in [-1,1].}
#'
#' \strong{Capping the maximum proportion of outliers (\code{max_p_outliers}).}
#' If \code{max_p_outliers < 1}, let \eqn{q_L = Q_x(\text{max\_p\_outliers})} and
#' \eqn{q_U = Q_x(1-\text{max\_p\_outliers})} be the lower/upper quantiles of \eqn{x}.
#' If the corresponding \eqn{|u|} at either quantile exceeds 1, \eqn{u} is rescaled
#' \emph{separately} on the negative and positive sides so that those quantiles land at
#' \eqn{|u|=1}. This guarantees that all observations between the two quantiles receive
#' positive weight. Note the bound applies per side, so up to \eqn{2\,\text{max\_p\_outliers}}
#' of observations can be treated as outliers overall.
#'
#' \strong{Fallback when for zero MAD / degeneracy (\code{pearson_fallback}).}
#' If a column has \eqn{\mathrm{MAD}=0} or the robust denominator becomes zero,
#' the following rules apply:
#' \itemize{
#'   \item \code{"none"} when correlations involving that column are \code{NA} (diagonal
#'         remains 1).
#'   \item \code{"hybrid"} when only the affected column switches to Pearson standardisation
#'         \eqn{\bar x_a = (x_a - \overline{x}) / \sqrt{\sum_b (x_b - \overline{x})^2}},
#'         yielding the hybrid correlation
#'         \deqn{\mathrm{bicor}_{\mathrm{hyb}}(x,y) = \sum_a \bar x_a\,\tilde y_a,}
#'         with the other column still robust-standardised.
#'   \item \code{"all"} when all columns use ordinary Pearson standardisation; the result
#'         equals \code{stats::cor(…, method="pearson")} when the NA policy matches.
#' }
#'
#' \strong{Handling missing values (\code{na_method}).}
#' \itemize{
#'   \item \code{"error"} (default): inputs must be finite; this yields a symmetric,
#'         positive semidefinite (PSD) matrix since \eqn{R = \tilde X^\top \tilde X}.
#'   \item \code{"pairwise"}: each \eqn{R_{jk}} is computed on the intersection of
#'         rows where both columns are finite. Pairs with fewer than 5 overlapping
#'         rows return \code{NA} (guarding against instability). Pairwise deletion can
#'         break PSD, as in the Pearson case.
#' }
#'
#' \strong{Row weights (\code{w}).}
#' When \code{w} is supplied (non-negative, length \eqn{m}), the weighted median
#' \eqn{\mathrm{med}_w(x)} and weighted MAD
#' \eqn{\mathrm{MAD}_w(x) = \mathrm{med}_w(|x - \mathrm{med}_w(x)|)} are used to form
#' \eqn{u}. The Tukey weights are then multiplied by the observation weights prior
#' to normalisation:
#' \deqn{\tilde x_a =
#' \frac{(x_a - \mathrm{med}_w(x))\,w_a\,w^{(\mathrm{obs})}_a}{
#'       \sqrt{\sum_b \big[(x_b - \mathrm{med}_w(x))\,w_b\,w^{(\mathrm{obs})}_b\big]^2}},}
#' where \eqn{w^{(\mathrm{obs})}_a \ge 0} are the user-supplied row weights and \eqn{w_a}
#' are the Tukey biweights built from the weighted median/MAD. Weighted pairwise
#' behaves analogously on each column pair's overlap.
#'
#' \strong{MAD choice (\code{mad_consistent}).}
#' Setting \code{mad_consistent = TRUE} multiplies the raw MAD by 1.4826 inside
#' \eqn{u}. Equivalently, it uses an effective tuning constant
#' \eqn{c^\star = c \times 1.4826}. The default \code{FALSE} reproduces the convention
#' in Langfelder & Horvath (2012).
#'
#' \strong{Optional sparsification (\code{sparse_threshold}).}
#' If provided, entries with \eqn{|r| < \text{sparse\_threshold}} are set to 0 and the
#' result is returned as a \code{"ddiMatrix"} (diagonal is forced to 1). This is a
#' post-processing step that does not alter the per-pair estimates.
#'
#' \strong{Computation and threads.}
#' Columns are robust-standardised in parallel and the matrix is formed as
#' \eqn{R = \tilde X^\top \tilde X}. \code{n_threads} selects the number of OpenMP
#' threads; by default it uses \code{getOption("matrixCorr.threads", 1L)}.
#'
#' \strong{Large-sample inference.}
#' For a pairwise estimate \eqn{r} computed from \eqn{n_{jk}} observed rows, the
#' standard large-sample summaries use
#' \deqn{
#' Z_{jk} = \frac{1}{2}\log\!\left(\frac{1 + r}{1 - r}\right)\sqrt{n_{jk} - 2}
#' }
#' and
#' \deqn{
#' T_{jk} = \sqrt{n_{jk} - 2}\;\frac{|r|}{\sqrt{1-r^2}}.
#' }
#' The reported p-value is the two-sided Student-\eqn{t} tail probability with
#' \eqn{n_{jk}-2} degrees of freedom. When \code{ci = TRUE}, the package also
#' reports an approximate Fisher-z confidence interval obtained from
#' \deqn{
#' z_{jk} = \operatorname{atanh}(r), \qquad
#' \operatorname{SE}(z_{jk}) = \frac{1}{\sqrt{n_{jk}-3}},
#' }
#' followed by back-transformation with \code{tanh()}. Confidence intervals are currently available only for dense,
#' unweighted outputs.
#'
#' \strong{Basic properties.}
#' \eqn{\mathrm{bicor}(a x + b,\; c y + d) = \mathrm{sign}(ac)\,\mathrm{bicor}(x,y)}.
#' With no missing data (and with per-column hybrid/robust standardisation), the
#' output is symmetric and PSD. As with Pearson, affine equivariance does not hold
#' for the associated biweight midcovariance.
#'
#' @references
#' Langfelder, P. & Horvath, S. (2012).
#' Fast R Functions for Robust Correlations and Hierarchical Clustering.
#' Journal of Statistical Software, 46(11), 1–17. \doi{10.18637/jss.v046.i11}
#'
#' @importFrom Matrix Matrix
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(2000 * 40), 2000, 40)
#' R <- bicor(X, c_const = 9, max_p_outliers = 1,
#'                        pearson_fallback = "hybrid")
#' print(attr(R, "method"))
#' summary(R)
#' R_ci <- bicor(X[, 1:5], ci = TRUE)
#' summary(R_ci)
#'
#' # Interactive viewing (requires shiny)
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   view_corr_shiny(R)
#' }
#'
#' @author Thiago de Paula Oliveira
#' @export
bicor <- function(
    data,
    na_method        = c("error", "pairwise"),
    ci               = FALSE,
    conf_level       = 0.95,
    n_threads        = getOption("matrixCorr.threads", 1L),
    c_const          = 9,
    max_p_outliers   = 1,
    pearson_fallback = c("hybrid", "none", "all"),
    mad_consistent   = FALSE,
    w                = NULL,
    sparse_threshold = NULL
) {
  if (missing(na_method) &&
      isFALSE(ci) &&
      missing(n_threads) &&
      identical(c_const, 9) &&
      identical(max_p_outliers, 1) &&
      missing(pearson_fallback) &&
      isFALSE(mad_consistent) &&
      is.null(w) &&
      is.null(sparse_threshold)) {
    numeric_data <- validate_corr_input(data)
    colnames_data <- colnames(numeric_data)
    prev_threads <- get_omp_threads()
    on.exit(set_omp_threads(as.integer(prev_threads)), add = TRUE)
    res <- bicor_matrix_cpp(
      numeric_data,
      c_const = 9,
      maxPOutliers = 1,
      pearson_fallback = 1L,
      n_threads = getOption("matrixCorr.threads", 1L)
    )
    if (!is.null(colnames_data)) {
      dimnames(res) <- .mc_square_dimnames(colnames_data)
    }
    return(.mc_structure_corr_matrix(
      res,
      class_name = "bicor",
      method = "biweight_mid_correlation",
      description = paste0(
        "Median/MAD-based biweight mid-correlation (bicor); max_p_outliers = 1",
        ", MAD = raw; fallback = hybrid; NA mode = error."
      )
    ))
  }

  pf <- match.arg(pearson_fallback)
  pf_int <- switch(pf, "none" = 0L, "hybrid" = 1L, "all" = 2L)
  na_method <- match.arg(na_method)

  # --- checks
  check_scalar_nonneg(c_const, arg = "c_const", strict = TRUE)
  check_scalar_numeric(max_p_outliers,
                       arg          = "max_p_outliers",
                       lower        = 0,
                       upper        = 1,
                       closed_lower = FALSE,
                       closed_upper = TRUE)
  check_bool(mad_consistent, arg = "mad_consistent")
  check_bool(ci, arg = "ci")
  if (isTRUE(ci)) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  }
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")

  if (!is.null(sparse_threshold)) {
    check_scalar_nonneg(sparse_threshold, arg = "sparse_threshold", strict = FALSE)
  }

  # --- validate/coerce input (allow NA only in pairwise mode)
  numeric_data <- if (na_method == "error") {
    validate_corr_input(data)
  } else {
    validate_corr_input(data, check_na = FALSE)
  }
  colnames_data <- colnames(numeric_data)
  w <- check_weights(w, n = nrow(numeric_data), arg = "w")
  if (isTRUE(ci) && !is.null(w)) {
    abort_bad_arg(
      "ci",
      message = "Confidence intervals are unavailable when {.arg w} is supplied."
    )
  }
  if (isTRUE(ci) && !is.null(sparse_threshold)) {
    abort_bad_arg(
      "ci",
      message = "Confidence intervals are unavailable when {.arg sparse_threshold} is supplied."
    )
  }

  # --- MAD consistency via effective c
  c_eff <- if (isTRUE(mad_consistent)) c_const * 1.4826 else c_const
  dn <- .mc_square_dimnames(colnames_data)

  # --- choose backend
  prev_threads <- get_omp_threads()
  on.exit(set_omp_threads(as.integer(prev_threads)), add = TRUE)
  if (is.null(w) && na_method == "error") {
    res <- bicor_matrix_cpp(
      numeric_data,
      c_const          = c_eff,
      maxPOutliers     = max_p_outliers,
      pearson_fallback = pf_int,
      n_threads        = n_threads
    )
  } else if (is.null(w) && na_method == "pairwise") {
    res <- bicor_matrix_pairwise_cpp(
      numeric_data,
      c_const          = c_eff,
      maxPOutliers     = max_p_outliers,
      pearson_fallback = pf_int,
      min_n            = 5L,
      n_threads        = n_threads
    )
  } else if (!is.null(w) && na_method == "error") {
    res <- bicor_matrix_weighted_cpp(
      numeric_data, w,
      c_const          = c_eff,
      maxPOutliers     = max_p_outliers,
      pearson_fallback = pf_int,
      n_threads        = n_threads
    )
  } else { # weighted + pairwise
    res <- bicor_matrix_weighted_pairwise_cpp(
      numeric_data, w,
      c_const          = c_eff,
      maxPOutliers     = max_p_outliers,
      pearson_fallback = pf_int,
      min_n            = 5L,
      n_threads        = n_threads
    )
  }
  res <- .mc_set_matrix_dimnames(res, colnames_data)

  # --- names & metadata
  desc <- paste0(
    "Median/MAD-based biweight mid-correlation (bicor); max_p_outliers = ", max_p_outliers,
    ", MAD = ", if (mad_consistent) "normal-consistent (1.4826 * raw)" else "raw",
    "; fallback = ", pf, "; NA mode = ", na_method, "."
  )
  diagnostics <- NULL
  ci_attr <- NULL
  inference_attr <- NULL
  if (isTRUE(ci)) {
    n_complete <- .mc_bicor_n_complete(numeric_data)
    diagnostics <- list(n_complete = .mc_set_matrix_dimnames(n_complete, colnames_data))
    inference_attr <- .mc_bicor_large_sample_inference(
      est = res,
      n_complete = diagnostics$n_complete
    )
    ci_attr <- .mc_bicor_fisher_ci(
      est = res,
      n_complete = diagnostics$n_complete,
      conf_level = conf_level
    )
  }

  # --- default: dense matrix with S3 class (original behaviour)
  if (is.null(sparse_threshold)) {
    return(.mc_structure_corr_matrix(
      res,
      class_name = "bicor",
      method = "biweight_mid_correlation",
      description = desc,
      diagnostics = diagnostics,
      extra_attrs = c(
        if (!is.null(ci_attr)) {
          list(
            ci = ci_attr,
            conf.level = conf_level
          )
        },
        if (!is.null(inference_attr)) {
          list(inference = inference_attr)
        }
      )
    ))
  }

  # --- optional sparse thresholding (S4 Matrix, no S3 class mutation)
  # zero-out small entries safely, keeping NAs intact
  idx <- !is.na(res) & (abs(res) < sparse_threshold)
  if (any(idx)) res[idx] <- 0
  diag(res) <- 1

  # mark degenerate columns (only diagonal finite) with NA diagonal for sparse reporting
  degenerate_idx <- which(rowSums(is.finite(res)) <= 1)
  if (length(degenerate_idx)) {
    res[cbind(degenerate_idx, degenerate_idx)] <- NA_real_
  }

  sparse_obj <- Matrix::forceSymmetric(Matrix::Matrix(res, sparse = TRUE), uplo = "U")
  sparse_info <- list(class = class(sparse_obj), threshold = sparse_threshold)
  res_dense <- as.matrix(sparse_obj)
  res_dense <- .mc_structure_corr_matrix(
    res_dense,
    class_name = "bicor",
    method = "biweight_mid_correlation",
    description = paste0(desc, " Sparse threshold = ", sparse_threshold, "."),
    extra_attrs = list(sparse_info = sparse_info)
  )
  res_dense
}

.mc_bicor_n_complete <- function(x) {
  x <- as.matrix(x)
  fin <- is.finite(x)
  n_complete <- crossprod(unclass(fin))
  storage.mode(n_complete) <- "integer"
  n_complete
}

.mc_bicor_large_sample_inference <- function(est, n_complete) {
  est <- as.matrix(est)
  n_complete <- as.matrix(n_complete)
  p <- ncol(est)
  dn <- dimnames(est)

  statistic <- matrix(NA_real_, p, p, dimnames = dn)
  p_value <- matrix(NA_real_, p, p, dimnames = dn)
  Z <- matrix(NA_real_, p, p, dimnames = dn)
  parameter <- matrix(NA_real_, p, p, dimnames = dn)

  valid <- is.finite(est) & is.finite(n_complete) & n_complete > 2
  exact <- valid & abs(est) >= 1
  approx <- valid & !exact

  parameter[valid] <- n_complete[valid] - 2
  if (any(exact)) {
    statistic[exact] <- Inf
    p_value[exact] <- 0
    Z[exact] <- sign(est[exact]) * Inf
  }
  if (any(approx)) {
    eps <- sqrt(.Machine$double.eps)
    est_clip <- pmax(pmin(est[approx], 1 - eps), -1 + eps)
    df <- n_complete[approx] - 2
    statistic[approx] <- sqrt(df) * abs(est_clip) / sqrt(1 - est_clip^2)
    p_value[approx] <- 2 * stats::pt(statistic[approx], df = df, lower.tail = FALSE)
    Z[approx] <- atanh(est_clip) * sqrt(df)
  }

  diag(parameter) <- NA_real_
  diag(p_value) <- 0
  diag(statistic) <- Inf
  diag(Z) <- Inf

  list(
    method = "large_sample_bicor",
    estimate = unclass(est),
    statistic = statistic,
    p_value = p_value,
    Z = Z,
    parameter = parameter,
    n_obs = n_complete,
    alternative = "two.sided"
  )
}

.mc_bicor_fisher_ci <- function(est, n_complete, conf_level = 0.95) {
  est <- as.matrix(est)
  n_complete <- as.matrix(n_complete)
  p <- ncol(est)
  dn <- dimnames(est)
  lwr <- matrix(NA_real_, p, p, dimnames = dn)
  upr <- matrix(NA_real_, p, p, dimnames = dn)
  diag(lwr) <- 1
  diag(upr) <- 1

  valid <- is.finite(est) & is.finite(n_complete) & n_complete > 3
  exact <- valid & abs(est) >= 1
  approx <- valid & !exact
  if (any(exact)) {
    lwr[exact] <- est[exact]
    upr[exact] <- est[exact]
  }
  if (any(approx)) {
    eps <- sqrt(.Machine$double.eps)
    est_clip <- pmax(pmin(est[approx], 1 - eps), -1 + eps)
    z <- atanh(est_clip)
    se <- 1 / sqrt(n_complete[approx] - 3)
    crit <- stats::qnorm(0.5 * (1 + conf_level))
    lwr[approx] <- tanh(z - crit * se)
    upr[approx] <- tanh(z + crit * se)
  }

  list(
    est = unclass(est),
    lwr.ci = lwr,
    upr.ci = upr,
    conf.level = conf_level,
    ci.method = "fisher_z_bicor"
  )
}

.mc_bicor_ci_attr <- function(x) {
  attr(x, "ci", exact = TRUE)
}

.mc_bicor_inference_attr <- function(x) {
  attr(x, "inference", exact = TRUE)
}

.mc_bicor_pairwise_summary <- function(object,
                                       digits = 4,
                                       ci_digits = 3,
                                       p_digits = 4,
                                       show_ci = NULL) {
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )
  check_inherits(object, "bicor")

  est <- as.matrix(object)
  rn <- rownames(est); cn <- colnames(est)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(est)))
  if (is.null(cn)) cn <- as.character(seq_len(ncol(est)))

  ci <- .mc_bicor_ci_attr(object)
  inf <- .mc_bicor_inference_attr(object)
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
        rec$fisher_z <- if (is.finite(inf$Z[i, j])) round(inf$Z[i, j], digits) else NA_real_
        rec$p_value <- if (is.finite(inf$p_value[i, j])) round(inf$p_value[i, j], p_digits) else NA_real_
      }
      rows[[k]] <- rec
    }
  }

  df <- do.call(rbind.data.frame, rows)
  rownames(df) <- NULL
  num_cols <- intersect(c("estimate", "lwr", "upr", "statistic", "fisher_z", "p_value"), names(df))
  int_cols <- intersect(c("n_complete"), names(df))
  for (nm in num_cols) df[[nm]] <- as.numeric(df[[nm]])
  for (nm in int_cols) df[[nm]] <- as.integer(df[[nm]])

  out <- .mc_finalize_summary_df(df, class_name = "summary.bicor")
  attr(out, "overview") <- .mc_summary_corr_matrix(object)
  attr(out, "has_ci") <- include_ci
  attr(out, "has_p") <- include_p
  attr(out, "conf.level") <- if (is.null(ci)) NA_real_ else ci$conf.level
  attr(out, "digits") <- digits
  attr(out, "ci_digits") <- ci_digits
  attr(out, "p_digits") <- p_digits
  attr(out, "inference_method") <- if (is.null(inf)) NA_character_ else inf$method
  out
}

#' @rdname bicor
#' @export
diag.bicor <- function(x, ...) {
  base::diag(as.matrix(x), ...)
}



#' @rdname bicor
#' @method print bicor
#' @title Print Method for \code{bicor} Objects
#'
#' @param x An object of class \code{bicor}.
#' @param digits Integer; number of decimal places used for the matrix.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param ci_digits Integer; digits for bicor confidence limits.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param na_print Character; how to display missing values.
#' @param ... Additional arguments passed to \code{print()}.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.bicor <- function(x,
                        digits   = 4,
                        n = NULL,
                        topn = NULL,
                        max_vars = NULL,
                        width    = NULL,
                        ci_digits = 3,
                        show_ci = NULL,
                        na_print = "NA",
                        ...) {
  check_inherits(x, "bicor")
  .mc_print_corr_matrix(
    x,
    header = "Biweight mid-correlation matrix (bicor)",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    na.print = na_print,
    ...
  )
}

#' @rdname bicor
#' @method plot bicor
#' @title Plot Method for \code{bicor} Objects
#'
#' @param x An object of class \code{bicor}.
#' @param title Plot title. Default is \code{"Biweight mid-correlation heatmap"}.
#' @param reorder Character; one of \code{"none"} (default) or \code{"hclust"}.
#'   If \code{"hclust"}, variables are reordered by complete-linkage clustering
#'   on the distance \eqn{d = 1 - r}, after replacing \code{NA} by 0 for
#'   clustering purposes only.
#' @param triangle One of \code{"full"} (default), \code{"lower"}, or \code{"upper"}
#'   to display the full matrix or a single triangle.
#' @param low_color,mid_color,high_color Colours for the gradient in
#'   \code{scale_fill_gradient2}. Defaults are \code{"indianred1"},
#'   \code{"white"}, \code{"steelblue1"}.
#' @param value_text_size Numeric; font size for cell labels. Set to \code{NULL}
#'   to suppress labels (recommended for large matrices).
#' @param ci_text_size Text size for confidence-interval labels in the heatmap.
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
#' @param na_fill Fill colour for \code{NA} cells. Default \code{"grey90"}.
#' @param ... Additional arguments passed to \code{ggplot2::theme()} or other layers.
#'
#' @return A \code{ggplot} object.
#' @import ggplot2
#' @importFrom stats as.dist hclust
#' @export
plot.bicor <- function(
    x,
    title = "Biweight mid-correlation heatmap",
    reorder  = c("none", "hclust"),
    triangle = c("full", "lower", "upper"),
    low_color = "indianred1",
    mid_color = "white",
    high_color = "steelblue1",
    value_text_size = 3,
    ci_text_size = 2.5,
    show_value = TRUE,
    na_fill = "grey90",
    ...
) {
  check_inherits(x, "bicor")
  check_bool(show_value, arg = "show_value")

  reorder  <- match.arg(reorder)
  triangle <- match.arg(triangle)

  mat <- as.matrix(x)

  # Optional reordering via robust-inspired clustering on 1 - r
  if (reorder == "hclust" && ncol(mat) > 1) {
    R_fill <- mat
    R_fill[is.na(R_fill)] <- 0   # for clustering only
    R_fill[R_fill > 1] <- 1; R_fill[R_fill < -1] <- -1
    d <- stats::as.dist(1 - R_fill)
    ord <- stats::hclust(d, method = "complete")$order
    mat <- mat[ord, ord, drop = FALSE]
  }

  # Prepare long format with indices for triangle filtering
  rn <- rownames(mat); cn <- colnames(mat)
  if (is.null(rn)) rn <- seq_len(nrow(mat)) else rn <- as.character(rn)
  if (is.null(cn)) cn <- seq_len(ncol(mat)) else cn <- as.character(cn)
  dimnames(mat) <- list(rn, cn)

  idx_r <- rep(rn, each = length(cn))
  idx_c <- rep(cn, times = length(rn))

  tm <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(tm) <- c("Var1", "Var2", "bicor")
  tm$Var1 <- as.character(tm$Var1)
  tm$Var2 <- as.character(tm$Var2)
  miss_r <- is.na(tm$Var1) | tm$Var1 == ""
  miss_c <- is.na(tm$Var2) | tm$Var2 == ""
  if (any(miss_r)) tm$Var1[miss_r] <- idx_r[miss_r]
  if (any(miss_c)) tm$Var2[miss_c] <- idx_c[miss_c]

  df <- tm
  df$Var1 <- factor(as.character(df$Var1), levels = rev(rn))
  df$Var2 <- factor(as.character(df$Var2), levels = cn)
  ci <- .mc_bicor_ci_attr(x)
  if (!is.null(ci) && is.matrix(ci$lwr.ci) && is.matrix(ci$upr.ci)) {
    lwr_mat <- ci$lwr.ci
    upr_mat <- ci$upr.ci
    if (reorder == "hclust" && ncol(mat) > 1) {
      lwr_mat <- lwr_mat[ord, ord, drop = FALSE]
      upr_mat <- upr_mat[ord, ord, drop = FALSE]
    }
    dimnames(lwr_mat) <- list(rn, cn)
    dimnames(upr_mat) <- list(rn, cn)
    df_lwr <- as.data.frame(as.table(lwr_mat), stringsAsFactors = FALSE)
    df_upr <- as.data.frame(as.table(upr_mat), stringsAsFactors = FALSE)
    names(df_lwr) <- c("Var1", "Var2", "lwr")
    names(df_upr) <- c("Var1", "Var2", "upr")
    df <- Reduce(
      function(a, b) merge(a, b, by = c("Var1", "Var2"), all.x = TRUE),
      list(df, df_lwr, df_upr)
    )
    df$Var1 <- factor(as.character(df$Var1), levels = rev(rn))
    df$Var2 <- factor(as.character(df$Var2), levels = cn)
  } else {
    df$lwr <- NA_real_
    df$upr <- NA_real_
  }

  # Triangle filtering
  if (triangle != "full") {
    # map factors back to integer indices
    i <- as.integer(df$Var1)         # note: Var1 is reversed
    j <- as.integer(df$Var2)
    # Because Var1 is reversed, convert to row indices in original order
    n <- length(levels(df$Var1))
    i_true <- n - i + 1
    keep <- switch(triangle,
                   lower = i_true >= j,
                   upper = i_true <= j)
    df <- df[keep, , drop = FALSE]
  }
  df$ci_label <- ifelse(
    is.na(df$lwr) | is.na(df$upr) | (as.character(df$Var1) == as.character(df$Var2)),
    NA_character_,
    sprintf("[%.3f, %.3f]", df$lwr, df$upr)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(Var2, Var1, fill = bicor)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(
      low = low_color, mid = mid_color, high = high_color,
      midpoint = 0, limits = c(-1, 1), na.value = na_fill, name = "bicor"
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
      ggplot2::aes(label = ifelse(is.na(bicor), "NA", sprintf("%.2f", bicor))),
      size = value_text_size, color = "black"
    )
  }
  if (isTRUE(show_value) && any(!is.na(df$ci_label))) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = ci_label, y = as.numeric(Var1) - 0.22),
      size = ci_text_size,
      color = "gray30",
      na.rm = TRUE
    )
  }

  return(p)
}

#' @rdname bicor
#' @method summary bicor
#' @param object An object of class \code{bicor}.
#' @param ci_digits Integer; digits for bicor confidence limits in the pairwise
#'   summary.
#' @param p_digits Integer; digits for bicor p-values in the pairwise summary.
#' @export
summary.bicor <- function(object, n = NULL, topn = NULL,
                          max_vars = NULL, width = NULL,
                          ci_digits = 3,
                          p_digits = 4,
                          show_ci = NULL, ...) {
  check_inherits(object, "bicor")
  if (is.null(.mc_bicor_ci_attr(object))) {
    return(.mc_summary_corr_matrix(object, topn = topn))
  }
  .mc_bicor_pairwise_summary(
    object,
    ci_digits = ci_digits,
    p_digits = p_digits,
    show_ci = show_ci
  )
}

#' @rdname bicor
#' @method print summary.bicor
#' @param x An object of class \code{summary.bicor}.
#' @export
print.summary.bicor <- function(x, digits = NULL, n = NULL,
                                topn = NULL, max_vars = NULL,
                                width = NULL, show_ci = NULL, ...) {
  extra_items <- c(inference = attr(x, "inference_method", exact = TRUE))
  .mc_print_pairwise_summary_digest(
    x,
    title = "Biweight mid-correlation summary",
    digits = .mc_coalesce(digits, 4),
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ci_method = "fisher_z_bicor",
    extra_items = extra_items,
    ...
  )
  invisible(x)
}

