#' Pairwise Chatterjee Rank Correlation
#'
#' @description
#' Computes the directed Chatterjee rank correlation coefficient for numeric
#' vectors or for all directed column pairs of a numeric matrix/data frame.
#' The matrix orientation is \code{result[i, j] = xi(V_i, V_j)}, where
#' \code{V_i} is the predictor/sorting variable and \code{V_j} is the
#' response/ranked variable.
#'
#' @param data A numeric matrix or a data frame with at least two numeric
#' columns, or a numeric predictor vector when \code{y} is supplied. All
#' non-numeric columns will be excluded in matrix/data-frame mode. Each
#' retained column must have at least two non-missing values.
#' @param y Optional numeric response vector for two-vector mode. When supplied,
#'   \code{xi_corr(data, y)} estimates \eqn{\xi(data, y)}, where \code{data}
#'   is the sorting/predictor variable and \code{y} is the response/ranked
#'   variable.
#' @param na_method Character scalar controlling missing-data handling.
#'   \code{"error"} rejects missing, \code{NaN}, and infinite values.
#'   \code{"pairwise"} recomputes each directed pair on its own pairwise
#'   complete-case overlap. This is permissive, but different directed entries
#'   may be based on different rows. \code{"complete"} performs listwise
#'   deletion once across the retained numeric columns and then computes all
#'   directed entries on the same complete sample.
#' @param ci Logical (default \code{FALSE}). If \code{TRUE}, attach
#'   m-out-of-n bootstrap confidence intervals for the directed Chatterjee
#'   estimates. In two-vector mode the scalar estimate carries a \code{ci}
#'   attribute.
#' @param conf_level Confidence level used when \code{ci = TRUE}. Default is
#'   \code{0.95}.
#' @param ci_method Confidence interval method:
#'   \itemize{
#'   \item \code{"auto"} (default): use \code{"dette_kroll"} when the
#'   pair-specific complete-case sample size is less than or equal to
#'   \code{large_sample_cutoff}, and \code{"n_choose_m"} for larger samples.
#'   \item \code{"dette_kroll"}: m-out-of-n variance estimation followed by a
#'   normal interval, following Dette and Kroll.
#'   \item \code{"n_choose_m"}: non-parametric m-out-of-n basic interval based
#'   on the centred and scaled bootstrap statistic, following the large-sample
#'   n-choose-m approach used by Dalitz, Arning and Goebbels.
#'   }
#' @param bootstrap_reps Number of m-out-of-n bootstrap replicates used when
#'   \code{ci = TRUE}. Larger values reduce Monte Carlo noise but increase
#'   runtime.
#' @param m Optional subsample size for the m-out-of-n bootstrap. If
#'   \code{NULL}, method-specific defaults are used: \code{floor(sqrt(n))} for
#'   \code{"dette_kroll"} and \code{round(2 * sqrt(n))} for
#'   \code{"n_choose_m"}, each bounded to \code{[2, n - 1]}.
#' @param large_sample_cutoff For \code{ci_method = "auto"}, use
#'   \code{"n_choose_m"} when pair-specific \code{n_complete} is greater than
#'   this cutoff; otherwise use \code{"dette_kroll"}.
#' @param bias_correction Finite-sample normalisation:
#'   \itemize{
#'   \item \code{"none"} (default): return Chatterjee's original finite-sample
#'   statistic.
#'   \item \code{"upper_bound"}: divide by the raw finite-sample upper-bound
#'   estimate \eqn{\xi_n(Y,Y)} and return
#'   \eqn{\max(-1, \xi_n(X,Y) / \xi_n(Y,Y))}. This can reduce the finite-sample
#'   downward bias of the raw statistic, but it is not the default.
#'   }
#' @param tie_method How ties in the sorting variable \eqn{X} are broken:
#'   \itemize{
#'   \item \code{"random"} (default): randomly permute observations within tied
#'   \eqn{X} groups, matching Chatterjee's random tie-breaking definition.
#'   \item \code{"first"}: preserve stable input order within tied \eqn{X}
#'   groups, useful for deterministic diagnostics and reproducible comparisons.
#'   }
#' @param seed Optional positive integer seed for reproducible tie breaking and
#'   bootstrap resampling.
#' @param n_threads Integer \eqn{\geq 1}. Number of OpenMP threads. Defaults to
#'   \code{getOption("matrixCorr.threads", 1L)}. Random tie breaking and
#'   confidence-interval paths are evaluated serially to preserve R RNG
#'   semantics.
#' @param output Output representation for the computed directed estimates.
#'   \itemize{
#'   \item \code{"matrix"} (default): full dense directed matrix; best when you
#'   need matrix algebra, heatmaps, or full compatibility with existing code.
#'   \item \code{"sparse"}: sparse matrix from \pkg{Matrix} containing retained
#'   directed entries; useful when thresholding leaves relatively few values.
#'   \item \code{"edge_list"}: long-form data frame with columns
#'   \code{row}, \code{col}, and \code{value}; convenient for filtering,
#'   joins, and graph/network workflows. Because \eqn{\xi(X,Y)} is directed,
#'   both directions are represented when retained.
#'   }
#' @param threshold Non-negative absolute-value filter for non-matrix outputs:
#'   keep entries with \code{abs(value) >= threshold}. Use
#'   \code{threshold > 0} when you want only stronger directed associations
#'   with \code{output = "sparse"} or \code{"edge_list"}. Must be \code{0}
#'   when \code{output = "matrix"}.
#' @param diag Logical; whether to include diagonal entries in
#'   \code{"sparse"} and \code{"edge_list"} outputs.
#' @param ... Compatibility arguments. The deprecated \code{check_na} alias is
#'   accepted.
#'
#' @return
#' A directed numeric matrix where the \code{(i, j)}-th element is
#' Chatterjee's \eqn{\xi} from the \code{i}-th numeric column to the
#' \code{j}-th numeric column. The dense matrix inherits from
#' \code{c("corr_matrix", "chatterjee_xi", "corr_result", "matrix")}. When
#' \code{ci = TRUE}, the object also carries a \code{ci} attribute with
#' elements \code{est}, \code{lwr.ci}, \code{upr.ci}, \code{conf.level},
#' \code{ci.method}, \code{se}, \code{m}, and \code{bootstrap_reps}. When
#' pairwise-complete evaluation is used, pairwise sample sizes are stored in
#' \code{attr(x, "diagnostics")$n_complete}. In two-vector mode, a numeric
#' scalar is returned; when \code{ci = TRUE}, it carries confidence-interval
#' attributes.
#'
#' @details
#' Let \eqn{(X_k,Y_k), k=1,\ldots,n}, be complete paired observations.
#' Chatterjee's rank correlation is directed: \eqn{\xi(X,Y)} measures how well
#' \eqn{Y} is functionally determined by \eqn{X}. It is not a symmetric
#' association measure, so in matrix mode
#' \deqn{\mathrm{result}_{ij} = \xi(X_{\cdot i}, X_{\cdot j})}
#' need not equal \eqn{\mathrm{result}_{ji}}. The row variable is always the
#' sorting/predictor variable and the column variable is the response/ranked
#' variable.
#'
#' For a given directed pair, observations are sorted by \eqn{X}. Let
#' \eqn{r_i} be the number of response values \eqn{Y_j \le Y_i}, and let
#' \eqn{l_i} be the number of response values \eqn{Y_j \ge Y_i}, evaluated for
#' the \eqn{i}-th observation in sorted-by-\eqn{X} order. The tied-response
#' finite-sample estimator is
#' \deqn{
#' \xi_n = 1 -
#' \frac{n \sum_{i=1}^{n-1} |r_{i+1} - r_i|}
#'      {2 \sum_{i=1}^n l_i (n - l_i)}.
#' }
#' If all response values are equal, the denominator is zero and the estimate
#' is \code{NA}. When there are no ties in \eqn{Y}, this reduces to the familiar
#' rank-difference expression
#' \deqn{
#' \xi_n = 1 -
#' \frac{3 \sum_{i=1}^{n-1}|r_{i+1}-r_i|}{n^2 - 1}.
#' }
#' The raw finite-sample statistic is not forced to one on the diagonal; with
#' no ties and \code{X = Y}, it equals \eqn{(n - 2)/(n + 1)}. Use
#' \code{bias_correction = "upper_bound"} only when this finite-sample
#' upper-bound normalisation is desired.
#'
#' Ties in the sorting variable \eqn{X} are handled before computing adjacent
#' rank differences. Chatterjee's definition uses random tie breaking, provided
#' here by \code{tie_method = "random"}. For reproducible diagnostics or when
#' the input order should define tied-\eqn{X} ordering, use
#' \code{tie_method = "first"}. Ties in the response variable \eqn{Y} are
#' handled by the general \eqn{r_i,l_i} formula above and do not require random
#' tie breaking.
#'
#' \strong{Confidence intervals.} The ordinary n-out-of-n bootstrap is not
#' used for Chatterjee's coefficient. Both available intervals use
#' m-out-of-n subsampling without replacement, consistent with the bootstrap
#' family considered for Chatterjee's rank correlation.
#' With \code{ci_method = "dette_kroll"}, subsamples of size \eqn{m} are drawn
#' without replacement and the limiting standard deviation is estimated as
#' \deqn{\widehat{\sigma} = \sqrt{m}\,\mathrm{sd}(\xi_m^*).}
#' The reported interval is the normal interval
#' \deqn{
#' \xi_n \pm z_{1-\alpha/2}\widehat{\sigma}/\sqrt{n},
#' \quad \alpha = 1-\texttt{conf\_level}.
#' }
#' If \code{m = NULL}, this method uses \eqn{m=\lfloor\sqrt{n}\rfloor}, bounded
#' to \eqn{[2,n-1]}.
#'
#' With \code{ci_method = "n_choose_m"}, subsamples are also drawn without
#' replacement, but the interval is obtained by inverting the empirical
#' quantiles of the centred and scaled statistic
#' \deqn{T^* = \sqrt{m}(\xi_m^* - \xi_n).}
#' If \eqn{q_{\alpha/2}} and \eqn{q_{1-\alpha/2}} are the bootstrap quantiles,
#' the basic m-out-of-n limits are
#' \deqn{
#' \xi_n - q_{1-\alpha/2}/\sqrt{n}
#' \quad\text{and}\quad
#' \xi_n - q_{\alpha/2}/\sqrt{n}.
#' }
#' If finite Monte Carlo quantiles fall entirely on one side of zero, the
#' inversion is anchored at zero so that the reported interval contains the
#' observed estimate; the bounds are not clipped to the population parameter
#' range. If \code{m = NULL}, this method uses
#' \eqn{m=\mathrm{round}(2\sqrt{n})}, bounded to \eqn{[2,n-1]}, following the
#' implementation rule used by Dalitz, Arning and Goebbels.
#'
#' \strong{Which CI method should be used?} The Dette-Kroll interval is the
#' more conservative default for small to moderate complete-case sample sizes
#' in this implementation, because it uses the m-out-of-n bootstrap only to
#' estimate a normal-approximation standard error. The non-parametric
#' n-choose-m interval is useful for larger samples when a direct
#' m-out-of-n bootstrap interval is preferred, but Dalitz, Arning and Goebbels
#' report that this interval may need fairly large \eqn{n} to approach nominal
#' coverage in some settings. Therefore \code{ci_method = "auto"} uses
#' \code{"dette_kroll"} when pair-specific \code{n_complete <=}
#' \code{large_sample_cutoff} and \code{"n_choose_m"} when
#' \code{n_complete > large_sample_cutoff}. Increase
#' \code{bootstrap_reps} for final analyses to reduce Monte Carlo error, and
#' consider setting \code{m} explicitly when a study protocol requires a
#' fixed subsample size.
#'
#' \strong{Computation.} For complete finite matrices without confidence
#' intervals, the C++ backend precomputes each column's sorting order and
#' response ranks and reuses them across directed pairs. This costs
#' \eqn{O(p\,n\log n + p^2 n)} time and \eqn{O(pn + p^2)} memory. Pairwise
#' missing-data evaluation and bootstrap confidence intervals recompute on the
#' relevant complete-case samples and are correspondingly more expensive.
#'
#' @references
#' Chatterjee, S. (2021). A New Coefficient of Correlation. \emph{Journal of
#' the American Statistical Association}, 116, 2009-2022.
#'
#' Dette, H. and Kroll, M. (2025). A simple bootstrap for Chatterjee's rank
#' correlation. \emph{Biometrika}, 112(1), asae045.
#' \doi{10.1093/biomet/asae045}.
#'
#' Lin, Z. and Han, F. (2025). Limit theorems of Chatterjee's rank correlation.
#' arXiv:2204.08031v4. \doi{10.48550/arXiv.2204.08031}.
#'
#' Dalitz, C., Arning, J. and Goebbels, S. (2024). A Simple Bias Reduction for
#' Chatterjee's Correlation.
#'
#' @examples
#' ## Example 1: independence versus functional dependence
#' ## Chatterjee's xi targets whether Y is determined by X.
#' set.seed(1)
#' n <- 300
#' x <- runif(n, -1, 1)
#' y_independent <- rnorm(n)
#' y_function <- sin(2 * pi * x)
#' c(
#'   independent = xi_corr(x, y_independent, tie_method = "first"),
#'   functional = xi_corr(x, y_function, tie_method = "first")
#' )
#'
#' ## Example 2: non-monotone functional dependence is directed
#' ## x determines x^2, but x^2 does not determine the sign of x.
#' ## The reverse raw finite-sample estimate can therefore be much smaller,
#' ## and may be negative when sorted response ranks oscillate strongly.
#' x <- seq(-1, 1, length.out = 300)
#' y <- x^2
#' c(
#'   xi_x_to_y = xi_corr(x, y, tie_method = "first"),
#'   xi_y_to_x = xi_corr(y, x, tie_method = "first")
#' )
#'
#' ## Example 3: the raw finite-sample diagonal is not forced to one
#' ## The optional upper-bound normalisation rescales this finite-sample ceiling.
#' z <- 1:20
#' c(
#'   raw = xi_corr(z, z, tie_method = "first"),
#'   upper_bound = xi_corr(
#'     z, z,
#'     tie_method = "first",
#'     bias_correction = "upper_bound"
#'   )
#' )
#'
#' ## Example 4: directed matrix workflow
#' X <- cbind(
#'   x = x,
#'   square = x^2,
#'   sine = sin(2 * pi * x),
#'   noise = rnorm(length(x))
#' )
#' xi <- xi_corr(X, tie_method = "first")
#' print(xi, digits = 3)
#' summary(xi)
#' estimate(xi)
#' tidy(xi)
#'
#' ## Example 5: Dette-Kroll bootstrap interval for a moderate sample
#' \donttest{
#' xi_ci <- xi_corr(
#'   X[, c("x", "sine")],
#'   ci = TRUE,
#'   ci_method = "dette_kroll",
#'   bootstrap_reps = 49,
#'   seed = 1,
#'   tie_method = "first"
#' )
#' summary(xi_ci)
#' ci(xi_ci)
#' confint(xi_ci)
#' plot(xi_ci)
#' }
#'
#' ## Example 6: n-choose-m interval for larger samples
#' \donttest{
#' set.seed(2)
#' x_large <- runif(1200, -1, 1)
#' y_large <- sin(2 * pi * x_large) + rnorm(1200, sd = 0.2)
#' xi_corr(
#'   x_large,
#'   y_large,
#'   ci = TRUE,
#'   ci_method = "n_choose_m",
#'   bootstrap_reps = 99,
#'   seed = 2,
#'   tie_method = "first"
#' )
#' }
#'
#' @useDynLib matrixCorr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @author Thiago de Paula Oliveira
#' @export
xi_corr <- function(data,
                    y = NULL,
                    na_method = c("error", "pairwise", "complete"),
                    ci = FALSE,
                    conf_level = 0.95,
                    ci_method = c("auto", "dette_kroll", "n_choose_m"),
                    bootstrap_reps = 999L,
                    m = NULL,
                    large_sample_cutoff = 1000L,
                    bias_correction = c("none", "upper_bound"),
                    tie_method = c("random", "first"),
                    seed = NULL,
                    n_threads = getOption("matrixCorr.threads", 1L),
                    output = c("matrix", "sparse", "edge_list"),
                    threshold = 0,
                    diag = TRUE,
                    ...) {
  output_cfg <- .mc_validate_thresholded_output_request(
    output = output,
    threshold = threshold,
    diag = diag
  )
  n_threads_missing <- missing(n_threads)

  if (...length() == 0L && missing(na_method)) {
    na_cfg <- list(na_method = "error", check_na = TRUE)
  } else {
    legacy_args <- .mc_extract_legacy_aliases(list(...), allowed = "check_na")
    na_cfg <- resolve_na_args(
      na_method = na_method,
      check_na = legacy_args$check_na %||% NULL,
      na_method_missing = missing(na_method),
      allowed = c("error", "pairwise", "complete")
    )
  }

  check_bool(ci, arg = "ci")
  if (isTRUE(ci)) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
    bootstrap_reps <- check_scalar_int_pos(bootstrap_reps, arg = "bootstrap_reps")
    large_sample_cutoff <- check_scalar_int_pos(large_sample_cutoff, arg = "large_sample_cutoff")
    if (!is.null(m)) {
      m <- check_scalar_int_pos(m, arg = "m")
      if (m < 2L) {
        abort_bad_arg("m", message = "must be at least 2.")
      }
    }
  }
  ci_method <- match.arg(ci_method)
  bias_correction <- match.arg(bias_correction)
  tie_method <- match.arg(tie_method)
  if (!is.null(seed)) {
    seed <- check_scalar_int_pos(seed, arg = "seed")
  }

  prev_threads <- .mc_prepare_omp_threads(
    n_threads,
    n_threads_missing = n_threads_missing
  )
  if (!is.null(prev_threads)) {
    on.exit(.mc_exit_omp_threads(prev_threads), add = TRUE)
  }

  eval_seeded <- function(expr) .mc_eval_with_seed(seed, expr)

  if (!is.null(y)) {
    if (!identical(output_cfg$output, "matrix")) {
      abort_bad_arg(
        "output",
        message = "non-matrix output modes are unavailable in two-vector mode."
      )
    }
    if (!is.numeric(data) || !is.numeric(y)) {
      abort_bad_arg(
        "data",
        message = "and {.arg y} must be numeric vectors for two-vector mode."
      )
    }
    check_same_length(data, y, arg_x = "data", arg_y = "y")
    data <- as.numeric(data)
    y <- as.numeric(y)
    if (na_cfg$check_na && (any(!is.finite(data)) || any(!is.finite(y)))) {
      abort_bad_arg(
        "data",
        message = "and {.arg y} must be free of NA/NaN/Inf when {.arg na_method} is {.val error}.",
        .hint = "Set `na_method = \"pairwise\"` or `na_method = \"complete\"` to use complete paired observations."
      )
    }
    if (!identical(na_cfg$na_method, "error")) {
      keep <- is.finite(data) & is.finite(y)
      if (sum(keep) < 2L) {
        abort_bad_arg(
          "data",
          message = "must retain at least 2 complete paired observations; retained {sum(keep)}."
        )
      }
      data <- data[keep]
      y <- y[keep]
    }

    if (!isTRUE(ci)) {
      return(eval_seeded(chatterjee_xi_vec_cpp(
        data,
        y,
        tie_method = tie_method,
        bias_correction = bias_correction
      )))
    }

    pair <- eval_seeded(chatterjee_xi_matrix_pairwise_cpp(
      cbind(data, y),
      return_ci = TRUE,
      conf_level = conf_level,
      ci_method = ci_method,
      bootstrap_reps = bootstrap_reps,
      m = if (is.null(m)) NULL else as.integer(m),
      large_sample_cutoff = large_sample_cutoff,
      tie_method = tie_method,
      bias_correction = bias_correction,
      n_threads = as.integer(n_threads)
    ))
    out <- as.numeric(pair$est[1L, 2L])
    attr(out, "ci") <- list(
      lwr = as.numeric(pair$lwr[1L, 2L]),
      upr = as.numeric(pair$upr[1L, 2L]),
      conf.level = pair$conf_level,
      ci.method = as.character(pair$ci_method[1L, 2L]),
      se = as.numeric(pair$se[1L, 2L]),
      m = as.integer(pair$m[1L, 2L]),
      bootstrap_reps = as.integer(pair$bootstrap_reps),
      n_complete = as.integer(pair$n_complete[1L, 2L])
    )
    attr(out, "conf.level") <- pair$conf_level
    attr(out, "ci.method") <- as.character(pair$ci_method[1L, 2L])
    class(out) <- c("chatterjee_xi_scalar", "chatterjee_xi", class(out))
    return(out)
  }

  numeric_data <- validate_corr_input(data, check_na = na_cfg$check_na)
  diagnostics_extra <- NULL
  if (identical(na_cfg$na_method, "complete")) {
    cc <- .mc_complete_case_matrix(numeric_data, min_n = 2L, arg = "data")
    numeric_data <- cc$data
    diagnostics_extra <- cc$diagnostics
  }

  colnames_data <- colnames(numeric_data)
  dn <- if (is.null(colnames_data)) NULL else .mc_square_dimnames(colnames_data)
  diagnostics <- NULL
  ci_attr <- NULL
  extra_attrs <- list()

  if (!identical(na_cfg$na_method, "pairwise") && !isTRUE(ci)) {
    result <- eval_seeded(chatterjee_xi_matrix_cpp(
      numeric_data,
      tie_method = tie_method,
      bias_correction = bias_correction,
      n_threads = as.integer(n_threads)
    ))
  } else {
    pairwise <- eval_seeded(chatterjee_xi_matrix_pairwise_cpp(
      numeric_data,
      return_ci = ci,
      conf_level = conf_level,
      ci_method = ci_method,
      bootstrap_reps = bootstrap_reps,
      m = if (is.null(m)) NULL else as.integer(m),
      large_sample_cutoff = large_sample_cutoff,
      tie_method = tie_method,
      bias_correction = bias_correction,
      n_threads = as.integer(n_threads)
    ))
    result <- pairwise$est
    diagnostics <- list(
      n_complete = .mc_set_matrix_dimnames(pairwise$n_complete, colnames_data)
    )
    if (isTRUE(ci)) {
      ci_method_payload <- pairwise$ci_method
      dimnames(ci_method_payload) <- dn
      se_payload <- .mc_set_matrix_dimnames(pairwise$se, colnames_data)
      m_payload <- .mc_set_matrix_dimnames(pairwise$m, colnames_data)
      diagnostics$ci_method <- ci_method_payload
      diagnostics$m <- m_payload
      diagnostics$bootstrap_reps <- as.integer(pairwise$bootstrap_reps)
      diagnostics$se <- se_payload
      ci_attr <- list(
        est = .mc_set_matrix_dimnames(unclass(result), colnames_data),
        lwr.ci = .mc_set_matrix_dimnames(unclass(pairwise$lwr), colnames_data),
        upr.ci = .mc_set_matrix_dimnames(unclass(pairwise$upr), colnames_data),
        conf.level = pairwise$conf_level,
        ci.method = ci_method_payload,
        se = se_payload,
        m = m_payload,
        bootstrap_reps = as.integer(pairwise$bootstrap_reps)
      )
      extra_attrs$ci.method <- if (identical(ci_method, "auto")) "auto" else ci_method
    }
  }
  diagnostics <- .mc_merge_diagnostics(diagnostics, diagnostics_extra)

  out <- .mc_structure_corr_matrix(
    result,
    class_name = "chatterjee_xi",
    method = "chatterjee_xi",
    description = "Pairwise Chatterjee rank correlation coefficient matrix (directed)",
    symmetric = FALSE,
    diagnostics = diagnostics,
    dimnames = dn,
    extra_attrs = c(
      list(
        bias_correction = bias_correction,
        tie_method = tie_method
      ),
      if (!is.null(ci_attr)) {
        c(
          list(
            ci = ci_attr,
            conf.level = conf_level
          ),
          extra_attrs
        )
      } else {
        list()
      }
    )
  )
  .mc_finalize_corr_output_fast(
    out,
    output = output_cfg$output,
    threshold = output_cfg$threshold,
    diag = output_cfg$diag
  )
}

#' @rdname xi_corr
#' @method print chatterjee_xi
#' @param x An object of class \code{chatterjee_xi}.
#' @param digits Integer; number of decimal places to print.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; \code{NULL}
#'   derives this from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param ci_digits Integer; digits for confidence limits.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @export
print.chatterjee_xi <- function(x,
                                digits = 4,
                                n = NULL,
                                topn = NULL,
                                max_vars = NULL,
                                width = NULL,
                                ci_digits = 3,
                                show_ci = NULL,
                                ...) {
  .mc_print_corr_matrix(
    x,
    header = "Chatterjee xi directed-dependence matrix",
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

#' @rdname xi_corr
#' @method print chatterjee_xi_scalar
#' @export
print.chatterjee_xi_scalar <- function(x,
                                       digits = 4,
                                       ci_digits = 4,
                                       ...) {
  ci <- attr(x, "ci", exact = TRUE)
  cat("Chatterjee xi directed-dependence estimate\n")
  cat("  method      : chatterjee_xi\n")
  cat("  estimate    : ", formatC(unname(as.numeric(x)[1L]), digits = digits, format = "f"), "\n", sep = "")
  if (is.list(ci)) {
    cat("  ci          : yes\n")
    cat("  ci_method   : ", ci$ci.method %||% NA_character_, "\n", sep = "")
    cat("  conf_level  : ", ci$conf.level %||% attr(x, "conf.level", exact = TRUE) %||% NA_real_, "\n", sep = "")
    cat("  n_complete  : ", ci$n_complete %||% NA_integer_, "\n", sep = "")
    cat("  m           : ", ci$m %||% NA_integer_, "\n", sep = "")
    cat("  boot_reps   : ", ci$bootstrap_reps %||% NA_integer_, "\n", sep = "")
    cat(
      "  interval    : [",
      formatC(ci$lwr %||% NA_real_, digits = ci_digits, format = "f"),
      ", ",
      formatC(ci$upr %||% NA_real_, digits = ci_digits, format = "f"),
      "]\n",
      sep = ""
    )
  } else {
    cat("  ci          : no\n")
  }
  invisible(x)
}

#' @rdname xi_corr
#' @method summary chatterjee_xi_scalar
#' @param object An object of class \code{chatterjee_xi_scalar}.
#' @export
summary.chatterjee_xi_scalar <- function(object,
                                         digits = 4,
                                         ci_digits = 4,
                                         ...) {
  ci <- attr(object, "ci", exact = TRUE)
  out <- data.frame(
    estimate = unname(as.numeric(object)[1L]),
    lwr = if (is.list(ci)) ci$lwr %||% NA_real_ else NA_real_,
    upr = if (is.list(ci)) ci$upr %||% NA_real_ else NA_real_,
    conf_level = if (is.list(ci)) ci$conf.level %||% attr(object, "conf.level", exact = TRUE) %||% NA_real_ else NA_real_,
    ci_method = if (is.list(ci)) ci$ci.method %||% NA_character_ else NA_character_,
    se = if (is.list(ci)) ci$se %||% NA_real_ else NA_real_,
    m = if (is.list(ci)) ci$m %||% NA_integer_ else NA_integer_,
    bootstrap_reps = if (is.list(ci)) ci$bootstrap_reps %||% NA_integer_ else NA_integer_,
    n_complete = if (is.list(ci)) ci$n_complete %||% NA_integer_ else NA_integer_,
    stringsAsFactors = FALSE
  )
  class(out) <- c("summary.chatterjee_xi_scalar", "summary.matrixCorr", "data.frame")
  attr(out, "method") <- "chatterjee_xi"
  attr(out, "summary_title") <- "Chatterjee xi directed-dependence estimate"
  attr(out, "digits") <- digits
  attr(out, "ci_digits") <- ci_digits
  out
}

#' @rdname xi_corr
#' @method print summary.chatterjee_xi_scalar
#' @export
print.summary.chatterjee_xi_scalar <- function(x,
                                               digits = NULL,
                                               ci_digits = NULL,
                                               ...) {
  digits <- digits %||% attr(x, "digits", exact = TRUE) %||% 4L
  ci_digits <- ci_digits %||% attr(x, "ci_digits", exact = TRUE) %||% digits
  cat(attr(x, "summary_title", exact = TRUE) %||% "Chatterjee xi directed-dependence estimate", "\n")
  cat("  method      : chatterjee_xi\n")
  cat("  estimate    : ", formatC(x$estimate[[1L]], digits = digits, format = "f"), "\n", sep = "")
  if (is.finite(x$lwr[[1L]]) || is.finite(x$upr[[1L]])) {
    cat("  ci          : yes\n")
    cat("  ci_method   : ", x$ci_method[[1L]], "\n", sep = "")
    cat("  conf_level  : ", x$conf_level[[1L]], "\n", sep = "")
    cat("  n_complete  : ", x$n_complete[[1L]], "\n", sep = "")
    cat("  m           : ", x$m[[1L]], "\n", sep = "")
    cat("  boot_reps   : ", x$bootstrap_reps[[1L]], "\n", sep = "")
    cat(
      "  interval    : [",
      formatC(x$lwr[[1L]], digits = ci_digits, format = "f"),
      ", ",
      formatC(x$upr[[1L]], digits = ci_digits, format = "f"),
      "]\n",
      sep = ""
    )
    if ("se" %in% names(x) && is.finite(x$se[[1L]])) {
      cat("  se          : ", formatC(x$se[[1L]], digits = ci_digits, format = "f"), "\n", sep = "")
    }
  } else {
    cat("  ci          : no\n")
  }
  invisible(x)
}
