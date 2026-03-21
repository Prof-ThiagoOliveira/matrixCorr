#' @title Short Aliases for Agreement and Repeated-Measures Methods
#' @description Convenience aliases for selected agreement and repeated-measures
#'   functions in \pkg{matrixCorr}. These aliases do not change computation or
#'   outputs; they simply forward to the corresponding full function names.
#'
#' @param group1,group2 Numeric vectors of equal length.
#' @param two Positive scalar; the multiple of the standard deviation used to
#'   define the limits of agreement.
#' @param mode Integer; 1 uses \code{group1 - group2}, 2 uses \code{group2 - group1}.
#' @param conf_level Confidence level.
#' @param verbose Logical; print brief progress or diagnostic output.
#' @param data A \code{data.frame} containing the repeated-measures dataset.
#' @param response Numeric response vector or column name, depending on the target method.
#' @param subject Subject identifier or subject column name.
#' @param method Method label or method column name.
#' @param time Replicate/time index or time column name.
#' @param include_slope Logical; whether to estimate proportional bias.
#' @param use_ar1 Logical; whether to use AR(1) within-subject correlation.
#' @param ar1_rho AR(1) parameter.
#' @param max_iter,tol EM control parameters.
#' @param rind Character; column identifying subjects for \code{ccc_rm_reml()}.
#' @param interaction,max_iter,Dmat,Dmat_type,Dmat_weights,Dmat_rescale,ci,
#'   ci_mode,digits,use_message,ar,ar_rho,slope,slope_var,slope_Z,
#'   drop_zero_cols,vc_select,vc_alpha,vc_test_order,include_subj_method,
#'   include_subj_time,sb_zero_tol Additional arguments forwarded to
#'   \code{ccc_lmm_reml()}.
#' @param Dmat Optional numeric weight matrix.
#' @param delta Numeric power exponent for U-statistics distances.
#' @param n_threads Integer number of OpenMP threads.
#'
#' @return The same object returned by the corresponding full function.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(40)
#' y <- x + rnorm(40, sd = 0.5)
#' ba(x, y)
#'
#' @export
ba <- function(group1,
               group2,
               two = 1.96,
               mode = 1L,
               conf_level = 0.95,
               verbose = FALSE) {
  bland_altman(
    group1 = group1,
    group2 = group2,
    two = two,
    mode = mode,
    conf_level = conf_level,
    verbose = verbose
  )
}

#' @rdname agreement_aliases
#' @export
ba_rm <- function(data = NULL, response, subject, method, time,
                  two = 1.96, conf_level = 0.95,
                  include_slope = FALSE,
                  use_ar1 = FALSE, ar1_rho = NA_real_,
                  max_iter = 200L, tol = 1e-6,
                  verbose = FALSE) {
  bland_altman_repeated(
    data = data,
    response = response,
    subject = subject,
    method = method,
    time = time,
    two = two,
    conf_level = conf_level,
    include_slope = include_slope,
    use_ar1 = use_ar1,
    ar1_rho = ar1_rho,
    max_iter = max_iter,
    tol = tol,
    verbose = verbose
  )
}

#' @rdname agreement_aliases
#' @export
ccc_rm_reml <- function(data, response, rind,
                        method = NULL, time = NULL, interaction = FALSE,
                        max_iter = 100, tol = 1e-6,
                        Dmat = NULL,
                        Dmat_type = c("time-avg","typical-visit","weighted-avg","weighted-sq"),
                        Dmat_weights = NULL,
                        Dmat_rescale = TRUE,
                        ci = FALSE, conf_level = 0.95,
                        ci_mode = c("auto","raw","logit"),
                        verbose = FALSE, digits = 4, use_message = TRUE,
                        ar = c("none", "ar1"),
                        ar_rho = NA_real_,
                        slope = c("none", "subject", "method", "custom"),
                        slope_var = NULL,
                        slope_Z = NULL,
                        drop_zero_cols = TRUE,
                        vc_select = c("auto","none"),
                        vc_alpha = 0.05,
                        vc_test_order = c("subj_time","subj_method"),
                        include_subj_method = NULL,
                        include_subj_time = NULL,
                        sb_zero_tol = 1e-10) {
  ccc_lmm_reml(
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

#' @rdname agreement_aliases
#' @export
ccc_rm_ustat <- function(data,
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
  ccc_pairwise_u_stat(
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
