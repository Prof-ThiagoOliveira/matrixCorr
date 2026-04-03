#' Deprecated Compatibility Wrappers
#'
#' Temporary wrappers for functions renamed in `matrixCorr` 1.0.0. These
#' wrappers preserve the pre-1.0 entry points while warning that they will be
#' removed in 2.0.0.
#'
#' @param group1,group2 Numeric vectors of equal length.
#' @param two Positive scalar; the multiple of the standard deviation used to
#'   define the limits of agreement.
#' @param mode Integer; 1 uses `group1 - group2`, 2 uses `group2 - group1`.
#' @param conf_level Confidence level.
#' @param verbose Logical; print brief progress or diagnostic output.
#' @param data A `data.frame`, matrix, or repeated-measures dataset accepted by
#'   the corresponding replacement function.
#' @param response Numeric response vector or column name, depending on the
#'   target method.
#' @param subject Subject identifier or subject column name.
#' @param method Method label or method column name.
#' @param time Replicate/time index or time column name.
#' @param include_slope Logical; whether to estimate proportional bias.
#' @param use_ar1 Logical; whether to use AR(1) within-subject correlation.
#' @param ar1_rho AR(1) parameter.
#' @param max_iter,tol EM control parameters.
#' @param c_const Positive numeric Tukey biweight tuning constant.
#' @param max_p_outliers Numeric in `(0, 1]`; optional cap on the maximum
#'   proportion of outliers on each side.
#' @param pearson_fallback Character fallback policy used by `bicor()`.
#' @param na_method Missing-data policy used by `bicor()`.
#' @param mad_consistent Logical; if `TRUE`, uses the consistency-corrected MAD.
#' @param w Optional vector of case weights.
#' @param sparse_threshold Optional threshold controlling sparse output.
#' @param n_threads Integer number of OpenMP threads.
#' @param rind Character; column identifying subjects for `ccc_rm_reml()`.
#' @param interaction Logical; forwarded to `ccc_rm_reml()`.
#' @param Dmat Optional distance matrix forwarded to `ccc_rm_reml()`.
#' @param Dmat_type Character selector controlling how `Dmat` is constructed.
#' @param Dmat_weights Optional weights used when `Dmat_type` requires them.
#' @param Dmat_rescale Logical; whether to rescale `Dmat`.
#' @param ci_mode Character selector for the confidence-interval scale used by
#'   `ccc_rm_reml()`.
#' @param digits Display precision forwarded to `ccc_rm_reml()`.
#' @param use_message Logical; whether the deprecated wrapper emits a lifecycle
#'   message.
#' @param ar Character selector for the within-subject residual correlation
#'   model.
#' @param ar_rho Numeric AR(1) parameter.
#' @param slope Character selector for the proportional-bias slope structure.
#' @param slope_var Optional covariance matrix for custom slopes.
#' @param slope_Z Optional design matrix for custom slopes.
#' @param drop_zero_cols Logical; whether zero-variance design columns are
#'   dropped.
#' @param vc_select Character selector controlling variance-component
#'   selection.
#' @param vc_alpha Significance level used in variance-component selection.
#' @param vc_test_order Character vector controlling the variance-component
#'   test order.
#' @param include_subj_method Optional logical override for the
#'   subject-by-method component.
#' @param include_subj_time Optional logical override for the subject-by-time
#'   component.
#' @param sb_zero_tol Numerical tolerance used when stabilising the scale-bias
#'   term.
#' @param delta Numeric power exponent for U-statistics distances.
#' @param lambda Numeric regularisation strength used by `pcorr()`.
#' @param return_cov_precision Logical; if `TRUE`, also return covariance and
#'   precision matrices.
#' @param ci Logical; if `TRUE`, request confidence intervals when supported by
#'   the replacement function.
#' @param check_na Logical validation flag used by `dcor()`.
#'
#' @details
#' Renamed functions:
#' \itemize{
#'   \item `bland_altman()` -> `ba()`
#'   \item `bland_altman_repeated()` -> `ba_rm()`
#'   \item `biweight_mid_corr()` -> `bicor()`
#'   \item `distance_corr()` -> `dcor()`
#'   \item `partial_correlation()` -> `pcorr()`
#'   \item `ccc_lmm_reml()` -> `ccc_rm_reml()`
#'   \item `ccc_pairwise_u_stat()` -> `ccc_rm_ustat()`
#' }
#'
#' The deprecated wrappers will be removed in `matrixCorr` 2.0.0.
#'
#' @name deprecated-matrixCorr
NULL

.mc_deprecate <- function(old,
                          new,
                          details = NULL,
                          remove_in = "2.0.0") {
  msg <- paste0(
    "`", old, "()` is deprecated in matrixCorr 1.0.0; use `", new,
    "()` instead. It will be removed in ", remove_in, ".",
    if (!is.null(details)) paste0(" ", details) else ""
  )

  .Deprecated(new = new, package = "matrixCorr", msg = msg)
  invisible(NULL)
}

#' @rdname deprecated-matrixCorr
#' @export
bland_altman <- function(group1,
                         group2,
                         two = 1.96,
                         mode = 1L,
                         conf_level = 0.95,
                         verbose = FALSE) {
  .mc_deprecate(
    old = "bland_altman",
    new = "ba",
    details = "Argument `two` was renamed to `loa_multiplier`."
  )

  ba(
    group1 = group1,
    group2 = group2,
    loa_multiplier = two,
    mode = mode,
    conf_level = conf_level,
    verbose = verbose
  )
}

#' @rdname deprecated-matrixCorr
#' @export
bland_altman_repeated <- function(data = NULL, response, subject, method, time,
                                  two = 1.96, conf_level = 0.95,
                                  include_slope = FALSE,
                                  use_ar1 = FALSE, ar1_rho = NA_real_,
                                  max_iter = 200L, tol = 1e-6,
                                  verbose = FALSE) {
  .mc_deprecate(
    old = "bland_altman_repeated",
    new = "ba_rm",
    details = "Argument `two` was renamed to `loa_multiplier`."
  )

  ba_rm(
    data = data,
    response = response,
    subject = subject,
    method = method,
    time = time,
    loa_multiplier = two,
    conf_level = conf_level,
    include_slope = include_slope,
    use_ar1 = use_ar1,
    ar1_rho = ar1_rho,
    max_iter = max_iter,
    tol = tol,
    verbose = verbose
  )
}

#' @rdname deprecated-matrixCorr
#' @export
biweight_mid_corr <- function(
    data,
    c_const = 9,
    max_p_outliers = 1,
    pearson_fallback = c("hybrid", "none", "all"),
    na_method = c("error", "pairwise"),
    mad_consistent = FALSE,
    w = NULL,
    sparse_threshold = NULL,
    n_threads = getOption("matrixCorr.threads", 1L)
) {
  .mc_deprecate(old = "biweight_mid_corr", new = "bicor")

  bicor(
    data = data,
    c_const = c_const,
    max_p_outliers = max_p_outliers,
    pearson_fallback = pearson_fallback,
    na_method = na_method,
    mad_consistent = mad_consistent,
    w = w,
    sparse_threshold = sparse_threshold,
    n_threads = n_threads
  )
}

#' @rdname deprecated-matrixCorr
#' @export
distance_corr <- function(data, check_na = TRUE) {
  .mc_deprecate(old = "distance_corr", new = "dcor")
  dcor(data = data, check_na = check_na)
}

#' @rdname deprecated-matrixCorr
#' @export
partial_correlation <- function(data,
                                method = c("oas", "ridge", "sample"),
                                lambda = 1e-3,
                                return_cov_precision = FALSE,
                                ci = FALSE,
                                conf_level = 0.95) {
  .mc_deprecate(
    old = "partial_correlation",
    new = "pcorr",
    details = "This wrapper preserves the pre-1.0 default `method = \"oas\"`."
  )

  pcorr(
    data = data,
    method = match.arg(method),
    lambda = lambda,
    return_cov_precision = return_cov_precision,
    ci = ci,
    conf_level = conf_level
  )
}

#' @rdname deprecated-matrixCorr
#' @export
ccc_lmm_reml <- function(data, response, rind,
                         method = NULL, time = NULL, interaction = FALSE,
                         max_iter = 100, tol = 1e-6,
                         Dmat = NULL,
                         Dmat_type = c("time-avg", "typical-visit", "weighted-avg", "weighted-sq"),
                         Dmat_weights = NULL,
                         Dmat_rescale = TRUE,
                         ci = FALSE, conf_level = 0.95,
                         ci_mode = c("auto", "raw", "logit"),
                         verbose = FALSE, digits = 4, use_message = TRUE,
                         ar = c("none", "ar1"),
                         ar_rho = NA_real_,
                         slope = c("none", "subject", "method", "custom"),
                         slope_var = NULL,
                         slope_Z = NULL,
                         drop_zero_cols = TRUE,
                         vc_select = c("auto", "none"),
                         vc_alpha = 0.05,
                         vc_test_order = c("subj_time", "subj_method"),
                         include_subj_method = NULL,
                         include_subj_time = NULL,
                         sb_zero_tol = 1e-10) {
  .mc_deprecate(old = "ccc_lmm_reml", new = "ccc_rm_reml")

  ccc_rm_reml(
    data = data,
    response = response,
    rind = rind,
    method = method,
    time = time,
    interaction = interaction,
    max_iter = max_iter,
    tol = tol,
    Dmat = Dmat,
    Dmat_type = Dmat_type,
    Dmat_weights = Dmat_weights,
    Dmat_rescale = Dmat_rescale,
    ci = ci,
    conf_level = conf_level,
    ci_mode = ci_mode,
    verbose = verbose,
    digits = digits,
    use_message = use_message,
    ar = ar,
    ar_rho = ar_rho,
    slope = slope,
    slope_var = slope_var,
    slope_Z = slope_Z,
    drop_zero_cols = drop_zero_cols,
    vc_select = vc_select,
    vc_alpha = vc_alpha,
    vc_test_order = vc_test_order,
    include_subj_method = include_subj_method,
    include_subj_time = include_subj_time,
    sb_zero_tol = sb_zero_tol
  )
}

#' @rdname deprecated-matrixCorr
#' @export
ccc_pairwise_u_stat <- function(data,
                                response,
                                method,
                                subject,
                                time = NULL,
                                Dmat = NULL,
                                delta = 1,
                                ci = FALSE,
                                conf_level = 0.95,
                                n_threads = getOption("matrixCorr.threads", 1L),
                                verbose = FALSE) {
  .mc_deprecate(old = "ccc_pairwise_u_stat", new = "ccc_rm_ustat")

  ccc_rm_ustat(
    data = data,
    response = response,
    method = method,
    subject = subject,
    time = time,
    Dmat = Dmat,
    delta = delta,
    ci = ci,
    conf_level = conf_level,
    n_threads = n_threads,
    verbose = verbose
  )
}
