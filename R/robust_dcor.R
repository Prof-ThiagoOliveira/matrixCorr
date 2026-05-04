#' @title Robust Distance Correlation
#'
#' @description
#' Computes robust distance correlations by applying the biloop transformation
#' to each numeric variable and then computing unbiased distance correlation on
#' the transformed variables.
#'
#' @param data A numeric matrix or a data frame with at least two numeric
#' columns. All non-numeric columns are dropped. Columns must be numeric.
#' @param na_method Character scalar controlling missing-data handling.
#'   \code{"error"} rejects missing, \code{NaN}, and infinite values.
#'   \code{"pairwise"} recomputes each association on its own pairwise
#'   complete-case overlap.
#' @param p_value Logical (default \code{FALSE}). If \code{TRUE}, attach
#' permutation p-values for pairwise independence.
#' @param inference Character scalar. Use \code{"none"} (default) for estimates
#' only or \code{"permutation"} with \code{p_value = TRUE}.
#' @param n_perm Positive integer; number of permutations used when
#' \code{p_value = TRUE}. Default \code{999}.
#' @param seed Optional positive integer seed for permutation inference.
#' @param c_const Positive numeric tuning constant for the biloop
#' transformation. Default \code{4}.
#' @param n_threads Integer \eqn{\geq 1}. Number of OpenMP threads. Defaults to
#' \code{getOption("matrixCorr.threads", 1L)}.
#' @param output Output representation for the computed estimates:
#' \code{"matrix"}, \code{"sparse"}, or \code{"edge_list"}.
#' @param threshold Non-negative absolute-value filter for non-matrix outputs.
#' Must be \code{0} when \code{output = "matrix"}.
#' @param diag Logical; whether to include diagonal entries in
#' \code{"sparse"} and \code{"edge_list"} outputs.
#' @param ... Deprecated compatibility aliases. Currently only
#' \code{check_na} is supported.
#' @param x,object An object of class \code{robust_dcor}.
#' @param digits Integer; number of decimal places to print.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns.
#' @param width Optional display width.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param title Plot title.
#' @param low_color,high_color,mid_color Colors used in the heatmap.
#' @param value_text_size,ci_text_size Text sizes for plot labels.
#' @param show_value Logical; overlay numeric values if \code{TRUE}.
#'
#' @return A symmetric correlation result. Dense output inherits from
#' \code{robust_dcor}, \code{corr_matrix}, and \code{matrix}. Sparse and
#' edge-list outputs use the package-standard \code{corr_sparse} and
#' \code{corr_edge_list} representations.
#'
#' @section Relationship to [dcor()]:
#' `robust_dcor()` is not a replacement for `dcor()`. It estimates distance
#' correlation after a bounded robust transformation. Large differences between
#' `dcor()` and `robust_dcor()` can indicate that the classical dependence
#' signal is driven by tail observations or outliers.
#'
#' @details
#' For each variable \eqn{x = (x_1,\ldots,x_n)}, the wrapper first robustly
#' standardises values using the median
#' \deqn{m_x = \mathrm{median}(x)}
#' and raw median absolute deviation
#' \deqn{s_x = \mathrm{median}\{|x_i - m_x|\}.}
#' If this scale is zero or non-finite, the implementation falls back to
#' \eqn{\mathrm{IQR}(x) / 1.34898} and then the ordinary sample standard
#' deviation. Columns with no positive finite fallback scale are treated as
#' degenerate.
#'
#' The standardised value \eqn{z_i = (x_i - m_x) / s_x} is mapped to two
#' bounded coordinates by the biloop transformation
#' \deqn{\theta_i = 2\pi\tanh(z_i / c),}
#' \deqn{u_i =
#' \begin{cases}
#' c\{1-\cos(\theta_i)\}, & z_i > 0,\\
#' -c\{1-\cos(\theta_i)\}, & z_i < 0,\\
#' 0, & z_i = 0,
#' \end{cases}}
#' and
#' \deqn{v_i = \sin(\theta_i).}
#' Pairwise distances are Euclidean distances in this transformed
#' two-dimensional space,
#' \deqn{D^{(j)}_{ab} =
#' \left[(u_{aj}-u_{bj})^2 + (v_{aj}-v_{bj})^2\right]^{1/2},\quad a \ne b,}
#' with zero diagonal.
#'
#' The transformed distance matrix is U-centred as
#' \deqn{A^{(j)}_{ab} =
#' D^{(j)}_{ab} - \frac{r^{(j)}_a + r^{(j)}_b}{n - 2}
#' + \frac{S^{(j)}}{(n - 1)(n - 2)},\quad a \ne b,}
#' where \eqn{r^{(j)}_a = \sum_{b \ne a}D^{(j)}_{ab}} and
#' \eqn{S^{(j)} = \sum_{a \ne b}D^{(j)}_{ab}}. The diagonal of
#' \eqn{A^{(j)}} is zero.
#'
#' For variables \eqn{j} and \eqn{k}, the unbiased robust distance covariance
#' is
#' \deqn{\widehat{\mathrm{dCov}}^2_u(j,k) =
#' \frac{2}{n(n-3)}\sum_{a < b} A^{(j)}_{ab} A^{(k)}_{ab}.}
#' The corresponding robust distance correlation is the usual non-negative
#' distance-correlation ratio based on this covariance and the two transformed
#' distance variances. Small negative numerical artifacts are clipped to zero.
#'
#' @note
#' `robust_dcor()` is more robust to extreme observations than classical
#' `dcor()`, but it may downweight genuine tail dependence. Classical
#' `dcor()` may be preferable when tail dependence is the scientific target.
#' Comparing both methods is recommended.
#'
#' @references
#' Leyder, J., Raymaekers, J., & Rousseeuw, P. J. (2025). Robust distance
#' correlation through bounded transformations.
#'
#' @seealso [dcor()], [wincor()], [pbcor()], [skipped_corr()]
#'
#' @examples
#' ## Non-linear dependence: both estimators detect association.
#' set.seed(1)
#' n <- 200
#' x <- rnorm(n)
#' y <- x^2 + rnorm(n, sd = 0.2)
#'
#' X <- cbind(x = x, y = y)
#'
#' classical <- dcor(X)
#' robust <- robust_dcor(X)
#' round(c(
#'   dcor = classical["x", "y"],
#'   robust_dcor = robust["x", "y"]
#' ), 3)
#'
#' ## One diagonal outlier can inflate classical dCor more than robust dCor.
#' set.seed(45)
#' x <- rnorm(20)
#' y <- rnorm(20)
#' x[1] <- 10
#' y[1] <- 10
#' X_out <- cbind(x = x, y = y)
#'
#' classical <- dcor(X_out)
#' robust <- robust_dcor(X_out)
#' round(c(
#'   dcor = classical["x", "y"],
#'   robust_dcor = robust["x", "y"]
#' ), 3)
#'
#' print(classical)
#' print(robust)
#' summary(robust)
#' plot(robust)
#'
#' ## Several variables.
#' set.seed(7)
#' z <- rnorm(120)
#' X_multi <- cbind(
#'   linear = z + rnorm(120, sd = 0.3),
#'   nonlinear = z^2 + rnorm(120, sd = 0.3),
#'   noise = rnorm(120),
#'   outlier = rnorm(120)
#' )
#' X_multi[1, "outlier"] <- 12
#' X_multi[1, "noise"] <- 12
#'
#' robust_multi <- robust_dcor(X_multi)
#'
#' print(robust_multi)
#' summary(robust_multi)
#' plot(robust_multi)
#'
#' @export
robust_dcor <- function(data,
                        na_method = c("error", "pairwise"),
                        p_value = FALSE,
                        inference = c("none", "permutation"),
                        n_perm = 999L,
                        seed = NULL,
                        c_const = 4,
                        n_threads = getOption("matrixCorr.threads", 1L),
                        output = c("matrix", "sparse", "edge_list"),
                        threshold = 0,
                        diag = TRUE,
                        ...) {
  output_cfg <- .mc_validate_output_args(
    output = output,
    threshold = threshold,
    diag = diag
  )

  legacy_args <- .mc_extract_legacy_aliases(list(...), allowed = "check_na")
  na_cfg <- resolve_na_args(
    na_method = na_method,
    check_na = legacy_args$check_na %||% NULL,
    na_method_missing = missing(na_method)
  )

  check_bool(p_value, arg = "p_value")
  inference <- match.arg(inference)
  check_scalar_numeric(c_const, arg = "c_const", lower = 0, closed_lower = FALSE)
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")

  if (isTRUE(p_value) && !identical(inference, "permutation")) {
    abort_bad_arg(
      "inference",
      message = "must be \"permutation\" when p_value = TRUE for robust_dcor().",
      .hint = "The classical distance-correlation t-test is not used for biloop-transformed robust dCor."
    )
  }

  if (isTRUE(p_value)) {
    n_perm <- check_scalar_int_pos(n_perm, arg = "n_perm")
    if (!is.null(seed)) seed <- check_scalar_int_pos(seed, arg = "seed")
  }

  numeric_data <- validate_corr_input(data, check_na = na_cfg$check_na)
  colnames_data <- colnames(numeric_data)
  dn <- .mc_square_dimnames(colnames_data)

  prev_threads <- .mc_prepare_omp_threads(
    n_threads,
    n_threads_missing = missing(n_threads)
  )
  if (!is.null(prev_threads)) {
    on.exit(.mc_exit_omp_threads(prev_threads), add = TRUE)
  }

  diagnostics <- NULL
  inference_attr <- NULL
  if (isTRUE(na_cfg$check_na) && !isTRUE(p_value)) {
    est <- robust_dcor_matrix_cpp(
      numeric_data,
      c_const = c_const,
      n_threads = n_threads
    )
  } else {
    res <- robust_dcor_matrix_pairwise_cpp(
      numeric_data,
      c_const = c_const,
      return_inference = isTRUE(p_value),
      n_perm = n_perm,
      seed = if (is.null(seed)) NULL else as.integer(seed),
      n_threads = n_threads
    )
    est <- res$est
    diagnostics <- list(
      n_complete = .mc_set_matrix_dimnames(unclass(res$n_complete), colnames_data),
      scale_method = "median_raw_mad_with_fallback"
    )

    if (isTRUE(p_value)) {
      inference_attr <- list(
        method = "permutation",
        estimate = .mc_set_matrix_dimnames(unclass(res$estimate), colnames_data),
        statistic = NULL,
        parameter = matrix(
          n_perm,
          nrow = ncol(numeric_data),
          ncol = ncol(numeric_data),
          dimnames = dn
        ),
        p_value = .mc_set_matrix_dimnames(unclass(res$p_value), colnames_data),
        alternative = "greater",
        n_perm = n_perm
      )
    }
  }

  est <- .mc_set_matrix_dimnames(est, colnames_data)

  out <- .mc_structure_corr_matrix(
    est,
    class_name = "robust_dcor",
    method = "robust_distance_correlation",
    description = paste0(
      "Biloop-transformed robust distance correlation; c_const = ",
      c_const,
      "; NA mode = ",
      na_cfg$na_method,
      "."
    ),
    symmetric = TRUE,
    diagnostics = diagnostics,
    dimnames = dn,
    extra_attrs = c(
      list(
        transform = "biloop",
        c_const = c_const,
        standardise = "median_raw_mad",
        classical_ref = "dcor"
      ),
      if (!is.null(inference_attr)) list(inference = inference_attr)
    )
  )

  .mc_finalize_corr_output_fast(
    out,
    output = output_cfg$output,
    threshold = output_cfg$threshold,
    diag = output_cfg$diag
  )
}

#' @rdname robust_dcor
#' @method print robust_dcor
#' @export
print.robust_dcor <- function(x,
                              digits = 4,
                              n = NULL,
                              topn = NULL,
                              max_vars = NULL,
                              width = NULL,
                              show_ci = NULL,
                              ...) {
  .mc_print_corr_matrix(
    x,
    header = "Robust distance correlation matrix",
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

#' @rdname robust_dcor
#' @method summary robust_dcor
#' @export
summary.robust_dcor <- function(object, topn = NULL, show_ci = NULL, ...) {
  NextMethod("summary")
}

#' @rdname robust_dcor
#' @method plot robust_dcor
#' @export
plot.robust_dcor <- function(x,
                             title = NULL,
                             low_color = "indianred1",
                             high_color = "steelblue1",
                             mid_color = "white",
                             value_text_size = 4,
                             ci_text_size = 3,
                             show_value = TRUE,
                             ...) {
  .mc_plot_corr_result(
    x,
    title = title %||% "Robust distance correlation heatmap",
    low_color = low_color,
    high_color = high_color,
    mid_color = mid_color,
    value_text_size = value_text_size,
    ci_text_size = ci_text_size,
    show_value = show_value,
    ...
  )
}
