#' @title Repeated-Measures Coefficient of Individual Agreement
#'
#' @description
#' Computes pairwise repeated-measures coefficients of individual agreement
#' (CIA) from long-format matched repeated-measures data using the categorical
#' condition ANOVA formulation of Haber, Gao, and Barnhart (2010).
#'
#' @details
#' `cia_rm()` is for matched repeated measurements under conditions such as
#' raters, visits, laboratories, treatments, or time points. It is not a
#' technical-replicate estimator. If the same subject has true technical
#' replicates within each subject-method cell, use [cia()] instead.
#'
#' Let `Y_ijk` denote the measurement on subject `i`, with method `j`, under
#' condition `k`, and consider the model
#' \deqn{
#' Y_{ijk} = \mu + \alpha_i + \beta_j + \gamma_k +
#' (\alpha\beta)_{ij} + (\alpha\gamma)_{ik} + (\beta\gamma)_{jk} + e_{ijk}.
#' }
#' Here subject is random, method and time/condition are fixed,
#' subject-method and subject-time are random interactions, method-time is a
#' fixed interaction, and `e_ijk` is the residual.
#'
#' For each method pair and condition `k`, the repeated-measures CIA is
#' \deqn{
#' \psi_k = \frac{2\sigma_e^2}{d_k^2 + 2\sigma_{\alpha\beta}^2 + 2\sigma_e^2},
#' }
#' where `d_k` is the condition-specific mean difference between the two
#' methods, \eqn{\sigma_{\alpha\beta}^2} is the subject-method variance
#' component, and \eqn{\sigma_e^2} is the residual variance component.
#'
#' The function fits this estimator separately to every method pair. The
#' homogeneity of agreement across conditions is tested by the method-time
#' interaction:
#' \deqn{
#' F = \frac{MS_{\mathrm{method}\times\mathrm{time}}}{MS_{\mathrm{error}}}.
#' }
#' The returned `homogeneity_F` attribute stores this test statistic for each
#' method pair, and `homogeneity_p` stores the corresponding upper-tail
#' `F`-test p-value using degrees of freedom `df_method_time` and `df_error`.
#' Larger `homogeneity_F` values and smaller `homogeneity_p` values indicate
#' stronger evidence that agreement changes across conditions, meaning that the
#' condition-specific CIA curves should be interpreted directly rather than
#' relying on a single pooled/common coefficient. Conversely, a small
#' `homogeneity_F` and a non-small `homogeneity_p` indicate that the data do
#' not show strong evidence of method-by-condition heterogeneity, so a common
#' estimate may be a reasonable summary if it is scientifically meaningful.
#'
#' When `homogeneous = TRUE`, the function also reports a common or pooled CIA
#' for each method pair from the reduced model that pools the method-time sum of
#' squares into the residual term. That common estimate is meaningful only when
#' agreement is reasonably homogeneous across conditions.
#'
#' Confidence intervals can be computed by a delta-method normal approximation
#' or by a subject-level percentile bootstrap. The delta-method interval is the
#' default. For each method pair and condition `k`, let
#' \eqn{\hat\psi_k} be the ANOVA estimator, let \eqn{g_k} be its numerical
#' gradient with respect to the subject-level moment vector, and let
#' \eqn{\Sigma} be the empirical covariance matrix of that moment vector divided
#' by the number of subjects. The delta-method standard error is
#' \deqn{
#' \widehat{\mathrm{se}}(\hat\psi_k) =
#' \sqrt{\max\left(g_k^\top \Sigma g_k, 0\right)}.
#' }
#' The raw normal interval is
#' \deqn{
#' \hat\psi_k \pm z_{1-\alpha/2}\,\widehat{\mathrm{se}}(\hat\psi_k).
#' }
#' Under the bootstrap option, the interval is given by the empirical
#' percentile limits from subject-resampled estimates. The argument
#' `estimator` controls whether the reported result is the literal
#' method-of-moments ratio or a bounded variance-component variant. Under
#' `estimator = "vc_constrained"`, the estimated subject-method variance
#' component is constrained to be non-negative before converting it back to CIA,
#' and the reported interval limits are also truncated to the CIA parameter
#' space,
#' \deqn{
#' \mathrm{lwr} = \max(0, \mathrm{lwr}_{\mathrm{raw}}), \qquad
#' \mathrm{upr} = \min(1, \mathrm{upr}_{\mathrm{raw}}).
#' }
#' This keeps the reported interval inside the CIA parameter space `[0, 1]`.
#' Use `estimator = "mom_unconstrained"` to inspect the literal raw
#' method-of-moments estimator and the corresponding unbounded interval on the
#' estimator scale.
#'
#' The current implementation requires exactly one observation in every
#' subject-method-time cell. It therefore targets the balanced repeated-measures
#' ANOVA setting from the cited paper and returns an error otherwise.
#'
#' @param data A data frame in long format.
#' @param response Character scalar naming the numeric response column.
#' @param subject Character scalar naming the subject identifier column.
#' @param method Character scalar naming the method column. Required.
#' @param time Character scalar naming the repeated condition column. Required.
#' @param ci Logical; if `TRUE`, attach confidence intervals.
#' @param conf_level Confidence level for intervals when `ci = TRUE`.
#' @param n_threads Positive integer thread hint passed to the C++ backend.
#' @param estimator One of `"vc_constrained"` or `"mom_unconstrained"`.
#'   `"vc_constrained"` is the default bounded variance-component estimator. It
#'   constrains the estimated subject-method variance component to be
#'   non-negative before converting it back to CIA. `"mom_unconstrained"` is
#'   the raw method-of-moments CIA ratio and may exceed 1.
#' @param inference One of `"delta"`, `"bootstrap"`, or `"none"`. The delta
#'   method is the default when `ci = TRUE`.
#' @param verbose Logical; if `TRUE`, emit progress messages for CI work.
#' @param digits Integer print precision carried on the returned object.
#' @param use_message Logical; if `TRUE`, use package messaging helpers when
#'   `verbose = TRUE`.
#' @param homogeneous Logical; if `TRUE`, also compute the pooled/common CIA
#'   from the reduced model without method-time interaction.
#' @param B Number of subject-bootstrap resamples when `ci = TRUE`.
#' @param seed Optional positive integer seed for reproducible bootstrap
#'   resampling.
#' @param ... Unused.
#'
#' @return
#' If `ci = FALSE`, the result is a list of class `c("cia_rm", "cia")`.
#' If `ci = TRUE`, the result is a list of class
#' `c("cia_rm", "cia_ci", "cia")`.
#'
#' Main components:
#' \itemize{
#'   \item `est`: condition-specific pairwise CIA estimates.
#'   \item `common`: pairwise overall summary CIA estimates from the reduced
#'   homogeneous model.
#'   \item `condition`: the ordered condition/time labels used in the output.
#'   \item `se`: standard errors for `est` when `ci = TRUE`.
#'   \item `lwr.ci` and `upr.ci`: lower and upper confidence limits for `est`
#'   when `ci = TRUE`.
#'   \item `common.se`: standard errors for `common` when `ci = TRUE`.
#'   \item `common.lwr.ci` and `common.upr.ci`: lower and upper confidence
#'   limits for `common` when `ci = TRUE`.
#' }
#'
#' When there are three or more methods, additional overall multi-method
#' components are returned:
#' \itemize{
#'   \item `overall`: condition-specific overall CIA across all methods.
#'   \item `overall.common`: the overall summary CIA across all methods when
#'   `homogeneous = TRUE`.
#'   \item `overall.se`: standard errors for `overall` when `ci = TRUE`.
#'   \item `overall.lwr.ci` and `overall.upr.ci`: lower and upper confidence
#'   limits for `overall` when `ci = TRUE`.
#'   \item `overall.common.se`: standard errors for `overall.common` when
#'   `ci = TRUE`.
#'   \item `overall.common.lwr.ci` and `overall.common.upr.ci`: lower and upper
#'   confidence limits for `overall.common` when `ci = TRUE`.
#' }
#'
#' The object carries pairwise ANOVA diagnostics as attributes, including
#' `sigma2_error`, `sigma2_subject_method`, `repeatability`,
#' `homogeneity_F`, `homogeneity_p`, `df_method_time`, `df_error`,
#' `n_obs`, `n_subjects`, `n_methods`, `n_times`, `time_levels`, and
#' `method_levels`. Here `homogeneity_F` is the method-time interaction test
#' statistic for each method pair, and `homogeneity_p` is the corresponding
#' p-value for the null hypothesis of homogeneous agreement across conditions.
#'
#' @seealso [ccc_rm_reml()], [icc_rm_reml()], and [cia()].
#'
#' @references
#' Haber M, Gao J, Barnhart HX. (2010). Evaluation of agreement between
#' measurement methods from data with matched repeated measurements via the
#' coefficient of individual agreement. \emph{Journal of Data Science}, 8,
#' 457-469.
#'
#' Barnhart HX, Kosinski AS, Haber M. (2007). Assessing individual agreement.
#' \emph{Journal of Biopharmaceutical Statistics}, 17(4), 697-719.
#'
#' @examples
#' # Example 1
#' set.seed(1)
#' dat_rater <- expand.grid(
#'   id = factor(sprintf("s%02d", 1:8)),
#'   method = factor(c("Device_A", "Device_B")),
#'   time = factor(c("Rater_1", "Rater_2", "Rater_3")),
#'   KEEP.OUT.ATTRS = FALSE
#' )
#' subj_eff <- rnorm(nlevels(dat_rater$id), sd = 1)[dat_rater$id]
#' method_shift <- c(Device_A = 0, Device_B = 0.25)[dat_rater$method]
#' rater_shift <- c(Rater_1 = 0, Rater_2 = 0.15, Rater_3 = -0.10)[dat_rater$time]
#' interaction_shift <- ifelse(
#'   dat_rater$method == "Device_B" & dat_rater$time == "Rater_3",
#'   0.12, 0
#' )
#' dat_rater$y <- subj_eff + method_shift + rater_shift +
#'   interaction_shift + rnorm(nrow(dat_rater), sd = 0.25)
#'
#' fit_rater <- cia_rm(
#'   dat_rater,
#'   response = "y",
#'   subject = "id",
#'   method = "method",
#'   time = "time",
#'   homogeneous = TRUE
#' )
#' print(fit_rater)
#' summary(fit_rater)
#' plot(fit_rater)
#'
#' # Example 2
#' set.seed(2)
#' dat_time <- expand.grid(
#'   id = factor(sprintf("s%02d", 1:10)),
#'   method = factor(c("Assay_A", "Assay_B", "Assay_C", "Assay_D")),
#'   time = factor(
#'     c("baseline", "week2", "month1", "month2", "month3"),
#'     levels = c("baseline", "week2", "month1", "month2", "month3")
#'   ),
#'   KEEP.OUT.ATTRS = FALSE
#' )
#' subj_eff <- rnorm(nlevels(dat_time$id), sd = 0.9)[dat_time$id]
#' method_shift <- c(Assay_A = 0, Assay_B = 0.20, Assay_C = -0.10, Assay_D = 0.35)[dat_time$method]
#' time_shift <- c(
#'   baseline = 0, week2 = 0.10, month1 = 0.22, month2 = 0.32, month3 = 0.42
#' )[dat_time$time]
#' interaction_shift <- ifelse(dat_time$method == "Assay_B" & dat_time$time == "month3", 0.10, 0) +
#'   ifelse(dat_time$method == "Assay_D" & dat_time$time == "month2", -0.08, 0)
#' dat_time$y <- subj_eff + method_shift + time_shift +
#'   interaction_shift + rnorm(nrow(dat_time), sd = 0.30)
#'
#' fit_time <- cia_rm(
#'   dat_time,
#'   response = "y",
#'   subject = "id",
#'   method = "method",
#'   time = "time",
#'   ci = TRUE
#' )
#' print(fit_time)
#' plot(fit_time)
#'
#' # Example 3
#' set.seed(3)
#' dat_days <- expand.grid(
#'   id = factor(sprintf("s%02d", 1:12)),
#'   method = factor(c("Sensor_A", "Sensor_B", "Sensor_C")),
#'   time = 1:15,
#'   KEEP.OUT.ATTRS = FALSE
#' )
#' subj_eff <- rnorm(nlevels(dat_days$id), sd = 0.8)[dat_days$id]
#' method_shift <- c(Sensor_A = 0, Sensor_B = 0.15, Sensor_C = -0.08)[dat_days$method]
#' day_trend <- 0.05 * dat_days$time
#' interaction_shift <- ifelse(dat_days$method == "Sensor_B", 0.01 * dat_days$time, 0) +
#'   ifelse(dat_days$method == "Sensor_C", -0.005 * dat_days$time, 0)
#' dat_days$y <- subj_eff + method_shift + day_trend +
#'   interaction_shift + rnorm(nrow(dat_days), sd = 0.22)
#'
#' fit_days <- cia_rm(
#'   dat_days,
#'   response = "y",
#'   subject = "id",
#'   method = "method",
#'   time = "time",
#'   ci = TRUE
#' )
#' plot(fit_days, facet_by_pair = TRUE)
#'
#' @author Thiago de Paula Oliveira
#' @export
cia_rm <- function(data,
                   response,
                   subject,
                   method = NULL,
                   time = NULL,
                   ci = FALSE,
                   conf_level = 0.95,
                   n_threads = getOption("matrixCorr.threads", 1L),
                   estimator = c("vc_constrained", "mom_unconstrained"),
                   inference = c("delta", "bootstrap", "none"),
                   verbose = FALSE,
                   digits = 4,
                   use_message = TRUE,
                   homogeneous = FALSE,
                   B = 1000L,
                   seed = NULL,
                   ...) {
  estimator <- match.arg(estimator)
  inference <- match.arg(inference)
  check_bool(ci, arg = "ci")
  check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")
  check_bool(verbose, arg = "verbose")
  digits <- check_scalar_int_pos(digits, arg = "digits")
  check_bool(use_message, arg = "use_message")
  check_bool(homogeneous, arg = "homogeneous")
  B <- check_scalar_int_pos(B, arg = "B")
  if (B < 2L) {
    abort_bad_arg("B", message = "must be >= 2.")
  }
  if (!is.null(seed)) {
    seed <- check_scalar_int_pos(seed, arg = "seed")
  }
  if (is.null(method)) {
    abort_bad_arg("method", message = "must be supplied.")
  }
  if (is.null(time)) {
    abort_bad_arg("time", message = "must be supplied.")
  }
  if (isTRUE(ci) && identical(inference, "none")) {
    abort_bad_arg("inference", message = "must not be {.val none} when {.arg ci} is TRUE.")
  }
  if (!isTRUE(ci)) {
    inference <- "none"
  }

  prep <- .mc_cia_rm_prepare(
    data = data,
    response = response,
    subject = subject,
    method = method,
    time = time
  )

  raw <- cia_rm_anova_cpp(
    y = prep$y,
    subject = prep$subject_code,
    method = prep$method_code,
    time = prep$time_code,
    n_subjects = prep$n_subjects,
    n_methods = prep$n_methods,
    n_times = prep$n_times,
    homogeneous = homogeneous,
    constrain_vc = identical(estimator, "vc_constrained"),
    n_threads = n_threads
  )

  out <- .mc_cia_rm_build_object(
    prep = prep,
    raw = raw,
    ci = ci,
    conf_level = conf_level,
    homogeneous = homogeneous,
    digits = digits
  )
  attr(out, "estimator") <- estimator

  if (!isTRUE(ci)) {
    return(out)
  }

  if (isTRUE(verbose)) {
    .mc_cia_rm_note(
      sprintf(
        "Computing %s confidence intervals for repeated-measures CIA.",
        if (identical(inference, "delta")) "delta-method" else "bootstrap"
      ),
      use_message = use_message
    )
  }

  if (identical(inference, "delta")) {
    delta <- .mc_cia_rm_delta(
      prep = prep,
      raw = raw,
      conf_level = conf_level,
      homogeneous = homogeneous,
      estimator = estimator
    )
    out$se <- delta$se
    out$lwr.ci <- delta$lwr
    out$upr.ci <- delta$upr
    out$common.se <- delta$common_se
    out$common.lwr.ci <- delta$common_lwr
    out$common.upr.ci <- delta$common_upr
    if (prep$n_methods >= 3L) {
      out$overall.se <- delta$overall_se
      out$overall.lwr.ci <- delta$overall_lwr
      out$overall.upr.ci <- delta$overall_upr
      out$overall.common.se <- delta$overall_common_se
      out$overall.common.lwr.ci <- delta$overall_common_lwr
      out$overall.common.upr.ci <- delta$overall_common_upr
    }
    attr(out, "ci.method") <- "delta_normal"
    attr(out, "estimator") <- estimator
  } else if (identical(inference, "bootstrap")) {
    boot <- .mc_cia_rm_bootstrap(
      prep = prep,
      homogeneous = homogeneous,
      B = B,
      conf_level = conf_level,
      seed = seed,
      n_threads = n_threads,
      estimator = estimator
    )
    out$lwr.ci <- boot$lwr
    out$upr.ci <- boot$upr
    out$common.lwr.ci <- boot$common_lwr
    out$common.upr.ci <- boot$common_upr
    if (prep$n_methods >= 3L) {
      out$overall.lwr.ci <- boot$overall_lwr
      out$overall.upr.ci <- boot$overall_upr
      out$overall.common.lwr.ci <- boot$overall_common_lwr
      out$overall.common.upr.ci <- boot$overall_common_upr
    }
    attr(out, "ci.method") <- "subject_bootstrap_percentile"
    attr(out, "estimator") <- estimator
    attr(out, "bootstrap") <- boot[c("n_successful", "B")]
  }

  out
}

.mc_cia_rm_prepare <- function(data, response, subject, method, time) {
  df <- as.data.frame(data)
  cols <- c(response, subject, method, time)
  names(cols) <- c("response", "subject", "method", "time")
  for (nm in names(cols)) {
    check_scalar_character(cols[[nm]], arg = nm)
  }
  check_required_cols(df, unname(cols), df_arg = "data")

  df <- df[, unname(cols), drop = FALSE]
  names(df) <- names(cols)
  time_raw <- df$time

  df$response <- suppressWarnings(as.numeric(df$response))
  if (any(!is.finite(df$response))) {
    abort_bad_arg("response", message = "must contain only finite numeric values.")
  }
  if (anyNA(df$subject) || anyNA(df$method) || anyNA(df$time)) {
    abort_bad_arg("data", message = "must not contain missing values in subject, method, or time.")
  }

  df$subject <- droplevels(factor(df$subject))
  df$method <- droplevels(factor(df$method))
  df$time <- droplevels(factor(df$time))
  time_meta <- .mc_cia_rm_time_meta(time_raw = time_raw, time_factor = df$time)

  if (nlevels(df$method) < 2L) {
    abort_bad_arg("method", message = "must retain at least two method levels.")
  }
  if (nlevels(df$time) < 2L) {
    abort_bad_arg("time", message = "must retain at least two time/condition levels.")
  }
  if (nlevels(df$subject) < 2L) {
    abort_bad_arg("subject", message = "must retain at least two subjects.")
  }

  counts <- stats::xtabs(~ subject + method + time, data = df)
  if (any(counts != 1L)) {
    abort_bad_arg(
      "data",
      message = paste(
        "cia_rm() currently implements the balanced categorical repeated-measures ANOVA estimator;",
        "each subject-method-time cell must contain exactly one observation."
      )
    )
  }

  ord <- order(df$subject, df$method, df$time)
  df <- df[ord, , drop = FALSE]
  rownames(df) <- NULL

  subject_levels <- levels(df$subject)
  method_levels <- levels(df$method)
  time_levels <- levels(df$time)
  subject_code <- as.integer(df$subject)
  method_code <- as.integer(df$method)
  time_code <- as.integer(df$time)
  subject_rows <- split(seq_len(nrow(df)), subject_code)

  list(
    data = df,
    y = df$response,
    subject_code = subject_code,
    method_code = method_code,
    time_code = time_code,
    subject_levels = subject_levels,
    method_levels = method_levels,
    time_levels = time_levels,
    n_subjects = length(subject_levels),
    n_methods = length(method_levels),
    n_times = length(time_levels),
    n_obs = nrow(df),
    subject_rows = subject_rows,
    rows_per_subject = length(method_levels) * length(time_levels),
    time_var = cols[["time"]],
    time_scale = time_meta$scale,
    time_plot_values = time_meta$values
  )
}

.mc_cia_rm_time_meta <- function(time_raw, time_factor) {
  time_levels <- levels(time_factor)
  level_index <- match(time_levels, as.character(time_factor))
  raw_levels <- time_raw[level_index]
  scale <- "categorical"
  values <- as.character(time_levels)

  if (inherits(raw_levels, "Date")) {
    scale <- "continuous"
    values <- as.Date(raw_levels)
  } else if (inherits(raw_levels, c("POSIXct", "POSIXlt"))) {
    scale <- "continuous"
    values <- as.POSIXct(raw_levels)
  } else if (is.numeric(time_raw) || is.integer(time_raw)) {
    scale <- "continuous"
    values <- as.numeric(raw_levels)
  } else {
    num_vals <- suppressWarnings(as.numeric(as.character(time_levels)))
    if (all(is.finite(num_vals))) {
      scale <- "continuous"
      values <- num_vals
    } else {
      date_vals <- tryCatch(
        suppressWarnings(as.Date(as.character(time_levels))),
        error = function(...) rep(as.Date(NA), length(time_levels))
      )
      if (all(!is.na(date_vals))) {
        scale <- "continuous"
        values <- date_vals
      }
    }
  }

  list(scale = scale, values = values)
}

.mc_cia_rm_array <- function(x, method_levels, time_levels) {
  arr <- array(
    as.numeric(x),
    dim = c(length(method_levels), length(method_levels), length(time_levels)),
    dimnames = list(method_levels, method_levels, time_levels)
  )
  arr
}

.mc_cia_rm_matrix <- function(x, method_levels) {
  mat <- as.matrix(x)
  dimnames(mat) <- .mc_square_dimnames(method_levels)
  mat
}

.mc_cia_rm_build_object <- function(prep, raw, ci, conf_level, homogeneous, digits) {
  est <- .mc_cia_rm_array(raw$est, prep$method_levels, prep$time_levels)
  common <- .mc_cia_rm_matrix(raw$common, prep$method_levels)
  common_sigma2_error <- .mc_cia_rm_matrix(raw$common_sigma2_error, prep$method_levels)
  common_sigma2_subject_method <- .mc_cia_rm_matrix(raw$common_sigma2_subject_method, prep$method_levels)
  common_repeatability <- .mc_cia_rm_matrix(raw$common_repeatability, prep$method_levels)
  has_overall <- prep$n_methods >= 3L
  overall <- if (has_overall) stats::setNames(as.numeric(raw$overall), prep$time_levels) else NULL
  overall_common <- if (has_overall) {
    if (isTRUE(homogeneous)) as.numeric(raw$overall_common) else NA_real_
  } else {
    NULL
  }

  if (!isTRUE(homogeneous)) {
    common[row(common) != col(common)] <- NA_real_
    common_sigma2_error[row(common_sigma2_error) != col(common_sigma2_error)] <- NA_real_
    common_sigma2_subject_method[row(common_sigma2_subject_method) != col(common_sigma2_subject_method)] <- NA_real_
    common_repeatability[row(common_repeatability) != col(common_repeatability)] <- NA_real_
  }

  out <- if (isTRUE(ci)) {
    structure(
      list(
        est = est,
        se = array(NA_real_, dim = dim(est), dimnames = dimnames(est)),
        lwr.ci = array(NA_real_, dim = dim(est), dimnames = dimnames(est)),
        upr.ci = array(NA_real_, dim = dim(est), dimnames = dimnames(est)),
        common = common,
        common.se = matrix(NA_real_, nrow = prep$n_methods, ncol = prep$n_methods, dimnames = .mc_square_dimnames(prep$method_levels)),
        common.lwr.ci = matrix(NA_real_, nrow = prep$n_methods, ncol = prep$n_methods, dimnames = .mc_square_dimnames(prep$method_levels)),
        common.upr.ci = matrix(NA_real_, nrow = prep$n_methods, ncol = prep$n_methods, dimnames = .mc_square_dimnames(prep$method_levels)),
        condition = prep$time_levels
      ),
      class = c("cia_rm", "cia_ci", "cia")
    )
  } else {
    structure(
      list(
        est = est,
        common = common,
        condition = prep$time_levels
      ),
      class = c("cia_rm", "cia")
    )
  }

  if (has_overall) {
    out$overall <- overall
    out$overall.common <- overall_common
    if (isTRUE(ci)) {
      out$overall.se <- stats::setNames(rep(NA_real_, prep$n_times), prep$time_levels)
      out$overall.lwr.ci <- stats::setNames(rep(NA_real_, prep$n_times), prep$time_levels)
      out$overall.upr.ci <- stats::setNames(rep(NA_real_, prep$n_times), prep$time_levels)
      out$overall.common.se <- NA_real_
      out$overall.common.lwr.ci <- NA_real_
      out$overall.common.upr.ci <- NA_real_
    }
  }

  attr(out, "method") <- "Repeated-measures coefficient of individual agreement"
  attr(out, "description") <- "Pairwise repeated-measures CIA from categorical-condition ANOVA"
  attr(out, "package") <- "matrixCorr"
  if (isTRUE(ci)) attr(out, "conf.level") <- conf_level
  attr(out, "sigma2_error") <- .mc_cia_rm_matrix(raw$sigma2_error, prep$method_levels)
  attr(out, "sigma2_subject_method") <- .mc_cia_rm_matrix(raw$sigma2_subject_method, prep$method_levels)
  attr(out, "repeatability") <- .mc_cia_rm_matrix(raw$repeatability, prep$method_levels)
  attr(out, "method_time_diff") <- .mc_cia_rm_array(raw$method_time_diff, prep$method_levels, prep$time_levels)
  attr(out, "ms_subject_method") <- .mc_cia_rm_matrix(raw$ms_subject_method, prep$method_levels)
  attr(out, "ms_method_time") <- .mc_cia_rm_matrix(raw$ms_method_time, prep$method_levels)
  attr(out, "ms_error") <- .mc_cia_rm_matrix(raw$ms_error, prep$method_levels)
  attr(out, "ss_subject_method") <- .mc_cia_rm_matrix(raw$ss_subject_method, prep$method_levels)
  attr(out, "ss_method_time") <- .mc_cia_rm_matrix(raw$ss_method_time, prep$method_levels)
  attr(out, "ss_error") <- .mc_cia_rm_matrix(raw$ss_error, prep$method_levels)
  attr(out, "df_subject_method") <- .mc_cia_rm_matrix(raw$df_subject_method, prep$method_levels)
  attr(out, "df_method_time") <- .mc_cia_rm_matrix(raw$df_method_time, prep$method_levels)
  attr(out, "df_error") <- .mc_cia_rm_matrix(raw$df_error, prep$method_levels)
  attr(out, "homogeneity_F") <- .mc_cia_rm_matrix(raw$homogeneity_F, prep$method_levels)
  attr(out, "homogeneity_p") <- .mc_cia_rm_matrix(raw$homogeneity_p, prep$method_levels)
  attr(out, "common_sigma2_error") <- common_sigma2_error
  attr(out, "common_sigma2_subject_method") <- common_sigma2_subject_method
  attr(out, "common_repeatability") <- common_repeatability
  attr(out, "overall_sigma2_error") <- if (has_overall) as.numeric(raw$overall_sigma2_error) else NULL
  attr(out, "overall_sigma2_subject_method") <- if (has_overall) as.numeric(raw$overall_sigma2_subject_method) else NULL
  attr(out, "overall_repeatability") <- if (has_overall) as.numeric(raw$overall_repeatability) else NULL
  attr(out, "overall_homogeneity_F") <- if (has_overall) as.numeric(raw$overall_homogeneity_F) else NULL
  attr(out, "overall_homogeneity_p") <- if (has_overall) as.numeric(raw$overall_homogeneity_p) else NULL
  attr(out, "overall_df_subject_method") <- if (has_overall) as.numeric(raw$overall_df_subject_method) else NULL
  attr(out, "overall_df_method_time") <- if (has_overall) as.numeric(raw$overall_df_method_time) else NULL
  attr(out, "overall_df_error") <- if (has_overall) as.numeric(raw$overall_df_error) else NULL
  attr(out, "overall_common_sigma2_error") <- if (has_overall) {
    if (isTRUE(homogeneous)) as.numeric(raw$overall_common_sigma2_error) else NA_real_
  } else {
    NULL
  }
  attr(out, "overall_common_sigma2_subject_method") <- if (has_overall) {
    if (isTRUE(homogeneous)) as.numeric(raw$overall_common_sigma2_subject_method) else NA_real_
  } else {
    NULL
  }
  attr(out, "overall_common_repeatability") <- if (has_overall) {
    if (isTRUE(homogeneous)) as.numeric(raw$overall_common_repeatability) else NA_real_
  } else {
    NULL
  }
  attr(out, "n_obs") <- prep$n_obs
  attr(out, "n_subjects") <- prep$n_subjects
  attr(out, "n_methods") <- prep$n_methods
  attr(out, "n_times") <- prep$n_times
  attr(out, "time_levels") <- prep$time_levels
  attr(out, "time_var") <- prep$time_var
  attr(out, "time_scale") <- prep$time_scale
  attr(out, "time_plot_values") <- prep$time_plot_values
  attr(out, "method_levels") <- prep$method_levels
  attr(out, "homogeneous") <- homogeneous
  attr(out, "digits") <- digits
  out
}

.mc_cia_rm_note <- function(text, use_message = TRUE) {
  if (isTRUE(use_message)) {
    cli::cli_inform(text)
  } else {
    base::message(text)
  }
}

.mc_cia_rm_cube <- function(prep) {
  aperm(
    array(prep$y, dim = c(prep$n_times, prep$n_methods, prep$n_subjects)),
    c(3L, 2L, 1L)
  )
}

.mc_cia_rm_pair_diff_matrix <- function(cube, method_a, method_b) {
  cube[, method_a, , drop = FALSE][, 1L, ] - cube[, method_b, , drop = FALSE][, 1L, ]
}

.mc_cia_rm_subject_cell_matrix <- function(cube) {
  n <- dim(cube)[1L]
  p <- dim(cube)[2L] * dim(cube)[3L]
  out <- matrix(NA_real_, nrow = n, ncol = p)
  for (i in seq_len(n)) {
    out[i, ] <- as.vector(cube[i, , ])
  }
  out
}

.mc_cia_rm_mean_pair_sqdiff <- function(x) {
  x <- as.numeric(x)
  J <- length(x)
  if (J < 2L) {
    return(NA_real_)
  }
  acc <- 0
  count <- 0L
  for (i in seq_len(J - 1L)) {
    for (j in (i + 1L):J) {
      acc <- acc + (x[i] - x[j])^2
      count <- count + 1L
    }
  }
  acc / count
}

.mc_cia_rm_overall_payload_from_moments <- function(m,
                                                    J,
                                                    K,
                                                    n,
                                                    homogeneous = TRUE,
                                                    constrain_vc = FALSE) {
  P <- J * K
  mu <- m[seq_len(P)]
  M <- matrix(m[P + seq_len(P * P)], nrow = P, ncol = P)
  mu_mat <- matrix(mu, nrow = J, ncol = K)
  Sx <- (n / (n - 1)) * (M - tcrossprod(mu))

  Pj <- diag(J) - matrix(1 / J, nrow = J, ncol = J)
  Pk <- diag(K) - matrix(1 / K, nrow = K, ncol = K)
  T_method <- matrix(0, nrow = J, ncol = P)
  for (k in seq_len(K)) {
    idx <- ((k - 1L) * J + 1L):(k * J)
    T_method[, idx] <- diag(J) / K
  }
  S_method <- T_method %*% Sx %*% t(T_method)
  H <- kronecker(Pk, Pj)

  ss_subject_method <- K * (n - 1) * sum(diag(Pj %*% S_method %*% Pj))
  ss_error <- (n - 1) * sum(diag(H %*% Sx %*% H))

  method_mean <- rowMeans(mu_mat)
  time_mean <- colMeans(mu_mat)
  grand_mean <- mean(mu_mat)
  interaction <- mu_mat - method_mean - rep(time_mean, each = J) + grand_mean
  ss_method_time <- n * sum(interaction^2)

  df_subject_method <- (n - 1) * (J - 1)
  df_method_time <- (J - 1) * (K - 1)
  df_error <- (n - 1) * (J - 1) * (K - 1)

  ms_subject_method <- ss_subject_method / df_subject_method
  ms_method_time <- ss_method_time / df_method_time
  ms_error <- ss_error / df_error

  sigma2_error <- ms_error
  sigma2_subject_method <- (ms_subject_method - ms_error) / K
  if (isTRUE(constrain_vc)) {
    sigma2_error <- max(sigma2_error, 0)
    sigma2_subject_method <- max(sigma2_subject_method, 0)
  }

  overall <- vapply(
    seq_len(K),
    function(k) {
      mean_sq_diff <- .mc_cia_rm_mean_pair_sqdiff(mu_mat[, k])
      denom <- mean_sq_diff + 2 * sigma2_subject_method + 2 * sigma2_error
      if (is.finite(denom) && denom > 0) {
        2 * sigma2_error / denom
      } else {
        NA_real_
      }
    },
    numeric(1)
  )

  homogeneity_F <- NA_real_
  homogeneity_p <- NA_real_
  if (is.finite(ms_method_time) && is.finite(ms_error)) {
    if (ms_error > 0) {
      homogeneity_F <- ms_method_time / ms_error
      homogeneity_p <- stats::pf(homogeneity_F, df1 = df_method_time, df2 = df_error, lower.tail = FALSE)
    } else if (ms_method_time > 0) {
      homogeneity_F <- Inf
      homogeneity_p <- 0
    }
  }

  overall_common <- NA_real_
  common_sigma2_error <- NA_real_
  common_sigma2_subject_method <- NA_real_
  common_repeatability <- NA_real_
  if (isTRUE(homogeneous)) {
    ss_error_reduced <- ss_method_time + ss_error
    df_error_reduced <- df_method_time + df_error
    ms_error_reduced <- ss_error_reduced / df_error_reduced
    common_sigma2_error <- ms_error_reduced
    common_sigma2_subject_method <- (ms_subject_method - ms_error_reduced) / K
    if (isTRUE(constrain_vc)) {
      common_sigma2_error <- max(common_sigma2_error, 0)
      common_sigma2_subject_method <- max(common_sigma2_subject_method, 0)
    }
    mean_sq_diff_common <- .mc_cia_rm_mean_pair_sqdiff(method_mean)
    denom_common <- mean_sq_diff_common + 2 * common_sigma2_subject_method + 2 * common_sigma2_error
    if (is.finite(denom_common) && denom_common > 0) {
      overall_common <- 2 * common_sigma2_error / denom_common
    }
    if (is.finite(common_sigma2_error) && common_sigma2_error >= 0) {
      common_repeatability <- 1.96 * sqrt(2 * common_sigma2_error)
    }
  }

  repeatability <- if (is.finite(sigma2_error) && sigma2_error >= 0) {
    1.96 * sqrt(2 * sigma2_error)
  } else {
    NA_real_
  }

  list(
    overall = overall,
    overall_common = overall_common,
    sigma2_error = sigma2_error,
    sigma2_subject_method = sigma2_subject_method,
    repeatability = repeatability,
    ms_subject_method = ms_subject_method,
    ms_method_time = ms_method_time,
    ms_error = ms_error,
    ss_subject_method = ss_subject_method,
    ss_method_time = ss_method_time,
    ss_error = ss_error,
    df_subject_method = df_subject_method,
    df_method_time = df_method_time,
    df_error = df_error,
    homogeneity_F = homogeneity_F,
    homogeneity_p = homogeneity_p,
    common_sigma2_error = common_sigma2_error,
    common_sigma2_subject_method = common_sigma2_subject_method,
    common_repeatability = common_repeatability
  )
}

.mc_cia_rm_overall_moment_stats <- function(X, J, K, constrain_vc = FALSE, homogeneous = TRUE) {
  n <- nrow(X)
  P <- ncol(X)
  mu <- colMeans(X)
  M <- crossprod(X) / n
  g <- cbind(
    X,
    t(vapply(seq_len(n), function(i) as.vector(tcrossprod(X[i, ])), numeric(P * P)))
  )
  Sigma_mean <- if (n >= 2L) stats::cov(g) / n else matrix(NA_real_, ncol(g), ncol(g))
  m <- c(mu, as.vector(M))
  payload <- .mc_cia_rm_overall_payload_from_moments(
    m = m,
    J = J,
    K = K,
    n = n,
    homogeneous = homogeneous,
    constrain_vc = constrain_vc
  )
  c(
    list(
      n = n,
      J = J,
      K = K,
      mu = mu,
      M = M,
      Sigma_mean = Sigma_mean
    ),
    payload
  )
}

.mc_cia_rm_pair_moment_stats <- function(D, constrain_vc = FALSE) {
  n <- nrow(D)
  K <- ncol(D)
  d <- colMeans(D)
  M <- Reduce(`+`, lapply(seq_len(n), function(i) tcrossprod(D[i, ]))) / n
  C <- (n / (n - 1)) * (M - tcrossprod(d))
  off_mean <- mean(C[row(C) != col(C)])
  diag_mean <- mean(diag(C))
  sigma2_subject_method <- off_mean / 2
  sigma2_error <- (diag_mean - off_mean) / 2
  if (isTRUE(constrain_vc)) {
    sigma2_subject_method <- max(sigma2_subject_method, 0)
  }
  denom <- d^2 + 2 * sigma2_subject_method + 2 * sigma2_error
  est <- ifelse(is.finite(denom) & denom > 0, 2 * sigma2_error / denom, NA_real_)
  d_bar <- mean(d)
  P <- diag(K) - matrix(1 / K, nrow = K, ncol = K)
  Q <- M - tcrossprod(d)
  ss_time <- n * as.numeric(crossprod(d - d_bar))
  ss_error_reduced <- n * sum(P * Q)
  ms_error_reduced <- (ss_time + ss_error_reduced) / (n * (K - 1))
  mean_a2 <- mean(M)
  ms_subject_method <- K * (n / (n - 1)) * (mean_a2 - d_bar^2)
  common_sigma2_error <- ms_error_reduced / 2
  common_sigma2_subject_method <- (ms_subject_method - ms_error_reduced) / (2 * K)
  if (isTRUE(constrain_vc)) {
    common_sigma2_subject_method <- max(common_sigma2_subject_method, 0)
  }
  common_denom <- d_bar^2 + 2 * common_sigma2_subject_method + 2 * common_sigma2_error
  common_est <- if (is.finite(common_denom) && common_denom > 0) {
    2 * common_sigma2_error / common_denom
  } else {
    NA_real_
  }

  g <- cbind(
    D,
    t(vapply(seq_len(n), function(i) as.vector(tcrossprod(D[i, ])), numeric(K * K)))
  )
  Sigma_mean <- if (n >= 2L) stats::cov(g) / n else matrix(NA_real_, ncol(g), ncol(g))

  list(
    n = n,
    K = K,
    d = d,
    M = M,
    Sigma_mean = Sigma_mean,
    sigma2_error = sigma2_error,
    sigma2_subject_method = sigma2_subject_method,
    est = est,
    common_sigma2_error = common_sigma2_error,
    common_sigma2_subject_method = common_sigma2_subject_method,
    common = common_est
  )
}

.mc_cia_rm_pair_est_from_moments <- function(m, K, n, constrain_vc = FALSE) {
  d <- m[seq_len(K)]
  M <- matrix(m[K + seq_len(K * K)], nrow = K, ncol = K)
  C <- (n / (n - 1)) * (M - tcrossprod(d))
  off_mean <- mean(C[row(C) != col(C)])
  diag_mean <- mean(diag(C))
  sigma2_subject_method <- off_mean / 2
  sigma2_error <- (diag_mean - off_mean) / 2
  if (isTRUE(constrain_vc)) {
    sigma2_subject_method <- max(sigma2_subject_method, 0)
  }
  denom <- d^2 + 2 * sigma2_subject_method + 2 * sigma2_error
  ifelse(is.finite(denom) & denom > 0, 2 * sigma2_error / denom, NA_real_)
}

.mc_cia_rm_common_est_from_moments <- function(m, K, n, constrain_vc = FALSE) {
  d <- m[seq_len(K)]
  M <- matrix(m[K + seq_len(K * K)], nrow = K, ncol = K)
  d_bar <- mean(d)
  P <- diag(K) - matrix(1 / K, nrow = K, ncol = K)
  Q <- M - tcrossprod(d)
  ss_time <- n * as.numeric(crossprod(d - d_bar))
  ss_error_reduced <- n * sum(P * Q)
  ms_error_reduced <- (ss_time + ss_error_reduced) / (n * (K - 1))
  mean_a2 <- mean(M)
  ms_subject_method <- K * (n / (n - 1)) * (mean_a2 - d_bar^2)
  sigma2_error <- ms_error_reduced / 2
  sigma2_subject_method <- (ms_subject_method - ms_error_reduced) / (2 * K)
  if (isTRUE(constrain_vc)) {
    sigma2_subject_method <- max(sigma2_subject_method, 0)
  }
  denom <- d_bar^2 + 2 * sigma2_subject_method + 2 * sigma2_error
  if (is.finite(denom) && denom > 0) {
    2 * sigma2_error / denom
  } else {
    NA_real_
  }
}

.mc_cia_rm_num_grad <- function(fun, m) {
  p <- length(m)
  g <- numeric(p)
  for (idx in seq_len(p)) {
    step <- 1e-6 * max(1, abs(m[[idx]]))
    plus <- m
    minus <- m
    plus[[idx]] <- plus[[idx]] + step
    minus[[idx]] <- minus[[idx]] - step
    g[[idx]] <- (fun(plus) - fun(minus)) / (2 * step)
  }
  g
}

.mc_cia_rm_apply_bounds <- function(lwr, upr, estimator) {
  if (identical(estimator, "vc_constrained")) {
    lwr <- pmax(0, lwr)
    upr <- pmin(1, upr)
  }
  list(lwr = lwr, upr = upr)
}

.mc_cia_rm_delta <- function(prep, raw, conf_level, homogeneous, estimator) {
  delta_raw <- cia_rm_delta_cpp(
    y = prep$y,
    subject = prep$subject_code,
    method = prep$method_code,
    time = prep$time_code,
    n_subjects = prep$n_subjects,
    n_methods = prep$n_methods,
    n_times = prep$n_times,
    homogeneous = homogeneous,
    constrain_vc = identical(estimator, "vc_constrained"),
    conf_level = conf_level,
    n_threads = 1L
  )

  out <- list(
    se = .mc_cia_rm_array(delta_raw$se, prep$method_levels, prep$time_levels),
    lwr = .mc_cia_rm_array(delta_raw$lwr, prep$method_levels, prep$time_levels),
    upr = .mc_cia_rm_array(delta_raw$upr, prep$method_levels, prep$time_levels),
    common_se = .mc_cia_rm_matrix(delta_raw$common_se, prep$method_levels),
    common_lwr = .mc_cia_rm_matrix(delta_raw$common_lwr, prep$method_levels),
    common_upr = .mc_cia_rm_matrix(delta_raw$common_upr, prep$method_levels),
    overall_se = NULL,
    overall_lwr = NULL,
    overall_upr = NULL,
    overall_common_se = NULL,
    overall_common_lwr = NULL,
    overall_common_upr = NULL
  )

  if (prep$n_methods >= 3L) {
    out$overall_se <- stats::setNames(as.numeric(delta_raw$overall_se), prep$time_levels)
    out$overall_lwr <- stats::setNames(as.numeric(delta_raw$overall_lwr), prep$time_levels)
    out$overall_upr <- stats::setNames(as.numeric(delta_raw$overall_upr), prep$time_levels)
    out$overall_common_se <- as.numeric(delta_raw$overall_common_se)
    out$overall_common_lwr <- as.numeric(delta_raw$overall_common_lwr)
    out$overall_common_upr <- as.numeric(delta_raw$overall_common_upr)
  }

  out
}

.mc_cia_rm_bootstrap <- function(prep, homogeneous, B, conf_level, seed = NULL, n_threads = 1L, estimator = "vc_constrained") {
  .mc_eval_with_seed(seed, {
    alpha <- (1 - conf_level) / 2
    est_arr <- array(NA_real_, dim = c(prep$n_methods, prep$n_methods, prep$n_times, B))
    common_arr <- array(NA_real_, dim = c(prep$n_methods, prep$n_methods, B))
    overall_arr <- if (prep$n_methods >= 3L) {
      matrix(NA_real_, nrow = prep$n_times, ncol = B, dimnames = list(prep$time_levels, NULL))
    } else {
      NULL
    }
    overall_common_arr <- if (prep$n_methods >= 3L) rep(NA_real_, B) else NULL

    for (b in seq_len(B)) {
      sampled <- sample.int(prep$n_subjects, size = prep$n_subjects, replace = TRUE)
      boot_rows <- integer(prep$n_obs)
      pos <- 1L
      for (i in seq_along(sampled)) {
        idx <- prep$subject_rows[[sampled[[i]]]]
        rng <- pos:(pos + length(idx) - 1L)
        boot_rows[rng] <- idx
        pos <- pos + length(idx)
      }

      y <- prep$y[boot_rows]
      method_code <- prep$method_code[boot_rows]
      time_code <- prep$time_code[boot_rows]
      subject_code <- rep(seq_len(prep$n_subjects), each = prep$rows_per_subject)

      boot_raw <- tryCatch(
        cia_rm_anova_cpp(
          y = y,
          subject = subject_code,
          method = method_code,
          time = time_code,
          n_subjects = prep$n_subjects,
          n_methods = prep$n_methods,
          n_times = prep$n_times,
          homogeneous = homogeneous,
          constrain_vc = identical(estimator, "vc_constrained"),
          n_threads = n_threads
        ),
        error = function(...) NULL
      )
      if (is.null(boot_raw)) {
        next
      }
      est_arr[, , , b] <- .mc_cia_rm_array(boot_raw$est, prep$method_levels, prep$time_levels)
      common_arr[, , b] <- .mc_cia_rm_matrix(boot_raw$common, prep$method_levels)
      if (prep$n_methods >= 3L) {
        overall_arr[, b] <- as.numeric(boot_raw$overall)
        overall_common_arr[[b]] <- as.numeric(boot_raw$overall_common)
      }
    }

    lwr <- array(NA_real_, dim = c(prep$n_methods, prep$n_methods, prep$n_times), dimnames = list(prep$method_levels, prep$method_levels, prep$time_levels))
    upr <- lwr
    common_lwr <- matrix(NA_real_, prep$n_methods, prep$n_methods, dimnames = .mc_square_dimnames(prep$method_levels))
    common_upr <- common_lwr
    overall_lwr <- if (prep$n_methods >= 3L) stats::setNames(rep(NA_real_, prep$n_times), prep$time_levels) else NULL
    overall_upr <- if (prep$n_methods >= 3L) stats::setNames(rep(NA_real_, prep$n_times), prep$time_levels) else NULL
    overall_common_lwr <- if (prep$n_methods >= 3L) NA_real_ else NULL
    overall_common_upr <- if (prep$n_methods >= 3L) NA_real_ else NULL
    for (k in seq_len(prep$n_times)) {
      diag(lwr[, , k]) <- 1
      diag(upr[, , k]) <- 1
    }
    diag(common_lwr) <- 1
    diag(common_upr) <- 1

    for (i in seq_len(prep$n_methods)) {
      for (j in seq_len(prep$n_methods)) {
        if (i == j) next
        for (k in seq_len(prep$n_times)) {
          vals <- est_arr[i, j, k, ]
          ok <- is.finite(vals)
          if (sum(ok) >= 2L) {
            qs <- stats::quantile(vals[ok], probs = c(alpha, 1 - alpha), names = FALSE, type = 6, na.rm = TRUE)
            bounded <- .mc_cia_rm_apply_bounds(as.numeric(qs[[1L]]), as.numeric(qs[[2L]]), estimator)
            lwr[i, j, k] <- bounded$lwr
            upr[i, j, k] <- bounded$upr
          }
        }
        vals_common <- common_arr[i, j, ]
        ok_common <- is.finite(vals_common)
        if (sum(ok_common) >= 2L) {
          qs <- stats::quantile(vals_common[ok_common], probs = c(alpha, 1 - alpha), names = FALSE, type = 6, na.rm = TRUE)
          bounded_common <- .mc_cia_rm_apply_bounds(as.numeric(qs[[1L]]), as.numeric(qs[[2L]]), estimator)
          common_lwr[i, j] <- bounded_common$lwr
          common_upr[i, j] <- bounded_common$upr
        }
      }
    }

    if (prep$n_methods >= 3L) {
      for (k in seq_len(prep$n_times)) {
        vals <- overall_arr[k, ]
        ok <- is.finite(vals)
        if (sum(ok) >= 2L) {
          qs <- stats::quantile(vals[ok], probs = c(alpha, 1 - alpha), names = FALSE, type = 6, na.rm = TRUE)
          bounded <- .mc_cia_rm_apply_bounds(as.numeric(qs[[1L]]), as.numeric(qs[[2L]]), estimator)
          overall_lwr[[k]] <- bounded$lwr
          overall_upr[[k]] <- bounded$upr
        }
      }
      ok_common <- is.finite(overall_common_arr)
      if (sum(ok_common) >= 2L) {
        qs <- stats::quantile(overall_common_arr[ok_common], probs = c(alpha, 1 - alpha), names = FALSE, type = 6, na.rm = TRUE)
        bounded_common <- .mc_cia_rm_apply_bounds(as.numeric(qs[[1L]]), as.numeric(qs[[2L]]), estimator)
        overall_common_lwr <- bounded_common$lwr
        overall_common_upr <- bounded_common$upr
      }
    }

    list(
      lwr = lwr,
      upr = upr,
      common_lwr = common_lwr,
      common_upr = common_upr,
      overall_lwr = overall_lwr,
      overall_upr = overall_upr,
      overall_common_lwr = overall_common_lwr,
      overall_common_upr = overall_common_upr,
      n_successful = sum(vapply(seq_len(B), function(b) any(is.finite(est_arr[, , , b])), logical(1))),
      B = B
    )
  })
}

.mc_cia_rm_flatten <- function(object) {
  est <- object$est
  se <- object$se
  method_levels <- attr(object, "method_levels", exact = TRUE)
  time_levels <- attr(object, "time_levels", exact = TRUE)
  sigma2_error <- attr(object, "sigma2_error", exact = TRUE)
  sigma2_subject_method <- attr(object, "sigma2_subject_method", exact = TRUE)
  repeatability <- attr(object, "repeatability", exact = TRUE)
  homogeneity_F <- attr(object, "homogeneity_F", exact = TRUE)
  homogeneity_p <- attr(object, "homogeneity_p", exact = TRUE)
  common <- object$common
  common_se <- object$common.se
  common_sigma2_error <- attr(object, "common_sigma2_error", exact = TRUE)
  common_sigma2_subject_method <- attr(object, "common_sigma2_subject_method", exact = TRUE)
  common_repeatability <- attr(object, "common_repeatability", exact = TRUE)
  has_ci <- inherits(object, "cia_ci")

  rows <- vector("list", 0L)
  for (i in seq_len(length(method_levels) - 1L)) {
    for (j in (i + 1L):length(method_levels)) {
      for (k in seq_along(time_levels)) {
        row <- data.frame(
          method_1 = method_levels[[i]],
          method_2 = method_levels[[j]],
          time = time_levels[[k]],
          cia = as.numeric(est[i, j, k]),
          sigma2_error = as.numeric(sigma2_error[i, j]),
          sigma2_subject_method = as.numeric(sigma2_subject_method[i, j]),
          repeatability = as.numeric(repeatability[i, j]),
          homogeneity_F = as.numeric(homogeneity_F[i, j]),
          homogeneity_p = as.numeric(homogeneity_p[i, j]),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        if (!is.null(se)) {
          row$se <- as.numeric(se[i, j, k])
        }
        if (has_ci) {
          row$lwr.ci <- as.numeric(object$lwr.ci[i, j, k])
          row$upr.ci <- as.numeric(object$upr.ci[i, j, k])
        }
        rows[[length(rows) + 1L]] <- row
      }

      if (is.finite(common[i, j])) {
        row_common <- data.frame(
          method_1 = method_levels[[i]],
          method_2 = method_levels[[j]],
          time = "common",
          cia = as.numeric(common[i, j]),
          sigma2_error = as.numeric(common_sigma2_error[i, j]),
          sigma2_subject_method = as.numeric(common_sigma2_subject_method[i, j]),
          repeatability = as.numeric(common_repeatability[i, j]),
          homogeneity_F = as.numeric(homogeneity_F[i, j]),
          homogeneity_p = as.numeric(homogeneity_p[i, j]),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        if (!is.null(common_se)) {
          row_common$se <- as.numeric(common_se[i, j])
        }
        if (has_ci) {
          row_common$lwr.ci <- as.numeric(object$common.lwr.ci[i, j])
          row_common$upr.ci <- as.numeric(object$common.upr.ci[i, j])
        }
        rows[[length(rows) + 1L]] <- row_common
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.mc_cia_rm_overall_df <- function(object) {
  overall <- object$overall
  if (is.null(overall)) {
    return(NULL)
  }
  has_ci <- inherits(object, "cia_ci")
  out <- data.frame(
    time = names(overall) %||% attr(object, "time_levels", exact = TRUE),
    cia = as.numeric(overall),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!is.null(object$overall.se)) {
    out$se <- as.numeric(object$overall.se)
  }
  if (has_ci && !is.null(object$overall.lwr.ci) && !is.null(object$overall.upr.ci)) {
    out$lwr.ci <- as.numeric(object$overall.lwr.ci)
    out$upr.ci <- as.numeric(object$overall.upr.ci)
  }
  out
}

.mc_cia_rm_overall_common_df <- function(object) {
  overall_common <- object$overall.common
  if (is.null(overall_common) || !is.finite(overall_common)) {
    return(NULL)
  }
  has_ci <- inherits(object, "cia_ci")
  out <- data.frame(
    time = "common",
    cia = as.numeric(overall_common),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!is.null(object$overall.common.se)) {
    out$se <- as.numeric(object$overall.common.se)
  }
  if (has_ci && !is.null(object$overall.common.lwr.ci) && !is.null(object$overall.common.upr.ci)) {
    out$lwr.ci <- as.numeric(object$overall.common.lwr.ci)
    out$upr.ci <- as.numeric(object$overall.common.upr.ci)
  }
  out
}

#' @rdname cia_rm
#' @method summary cia_rm
#' @param object A `cia_rm` object.
#' @param ci_digits Integer; number of digits for confidence interval bounds.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns.
#' @param width Optional display width.
#' @param show_ci One of `"yes"` or `"no"`.
#' @param ... Unused.
#' @export
summary.cia_rm <- function(object,
                           digits = 4,
                           ci_digits = 3,
                           n = NULL,
                           topn = NULL,
                           max_vars = NULL,
                           width = NULL,
                           show_ci = NULL,
                           ...) {
  out <- .mc_cia_rm_flatten(object)
  num_cols <- c("cia", "sigma2_error", "sigma2_subject_method", "repeatability", "homogeneity_F", "homogeneity_p")
  if ("se" %in% names(out)) {
    num_cols <- c("se", num_cols)
  }
  for (nm in intersect(num_cols, names(out))) {
    out[[nm]] <- round(out[[nm]], digits)
  }
  for (nm in intersect(c("lwr.ci", "upr.ci"), names(out))) {
    out[[nm]] <- round(out[[nm]], ci_digits)
  }
  out <- .mc_finalize_summary_df(out, class_name = "summary.cia_rm", repeated = TRUE)
  attr(out, "conf.level") <- attr(object, "conf.level", exact = TRUE)
  attr(out, "has_ci") <- inherits(object, "cia_ci")
  attr(out, "digits") <- digits
  attr(out, "ci_digits") <- ci_digits
  attr(out, "method") <- attr(object, "method", exact = TRUE)
  attr(out, "description") <- attr(object, "description", exact = TRUE)
  attr(out, "n_subjects_total") <- attr(object, "n_subjects", exact = TRUE)
  attr(out, "n_methods_total") <- attr(object, "n_methods", exact = TRUE)
  attr(out, "n_times_total") <- attr(object, "n_times", exact = TRUE)
  attr(out, "homogeneous") <- attr(object, "homogeneous", exact = TRUE)
  overall_out <- .mc_cia_rm_overall_df(object)
  if (!is.null(overall_out)) {
    for (nm in intersect(c("cia", "se"), names(overall_out))) {
      overall_out[[nm]] <- round(overall_out[[nm]], digits)
    }
    for (nm in intersect(c("lwr.ci", "upr.ci"), names(overall_out))) {
      overall_out[[nm]] <- round(overall_out[[nm]], ci_digits)
    }
  }
  overall_common_out <- .mc_cia_rm_overall_common_df(object)
  if (!is.null(overall_common_out)) {
    for (nm in intersect(c("cia", "se"), names(overall_common_out))) {
      overall_common_out[[nm]] <- round(overall_common_out[[nm]], digits)
    }
    for (nm in intersect(c("lwr.ci", "upr.ci"), names(overall_common_out))) {
      overall_common_out[[nm]] <- round(overall_common_out[[nm]], ci_digits)
    }
  }
  attr(out, "overall_section") <- overall_out
  attr(out, "overall_common_section") <- overall_common_out
  out
}

#' @rdname cia_rm
#' @method print cia_rm
#' @param x A `cia_rm` object.
#' @export
print.cia_rm <- function(x,
                         digits = 4,
                         n = NULL,
                         topn = NULL,
                         max_vars = NULL,
                         width = NULL,
                         show_ci = NULL,
                         ...) {
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("print_show_ci", "yes")
  )

  legacy <- .mc_extract_legacy_display_args(
    list(...),
    n = n,
    topn = topn,
    max_vars = max_vars
  )
  cfg <- .mc_resolve_display_args(
    context = "print",
    n = legacy$n,
    topn = legacy$topn,
    max_vars = legacy$max_vars,
    width = width,
    show_ci = show_ci
  )

  has_ci <- inherits(x, "cia_ci")
  cl <- suppressWarnings(as.numeric(attr(x, "conf.level", exact = TRUE)))
  if (length(cl) != 1L || !is.finite(cl)) cl <- NA_real_

  digest <- c(
    n_subjects = .mc_count_fmt(attr(x, "n_subjects", exact = TRUE)),
    n_methods = .mc_count_fmt(attr(x, "n_methods", exact = TRUE)),
    n_times = .mc_count_fmt(attr(x, "n_times", exact = TRUE))
  )
  if (isTRUE(attr(x, "homogeneous", exact = TRUE))) {
    digest <- c(digest, common = "yes")
  }
  if (has_ci && identical(cfg$show_ci, "yes")) {
    digest <- c(digest, ci = sprintf("%g%%", 100 * cl))
  }
  .mc_print_named_digest(
    digest,
    header = "Repeated-measures coefficient of individual agreement"
  )
  cat("\n")

  preview <- .mc_cia_rm_flatten(x)
  keep <- c("method_1", "method_2", "time", "cia", "lwr.ci", "upr.ci")
  preview <- preview[, intersect(keep, names(preview)), drop = FALSE]
  if ("cia" %in% names(preview)) {
    preview$cia <- round(preview$cia, digits)
  }
  for (nm in intersect(c("lwr.ci", "upr.ci"), names(preview))) {
    preview[[nm]] <- round(preview[[nm]], digits)
  }
  if (!has_ci || identical(cfg$show_ci, "no")) {
    preview <- preview[, setdiff(names(preview), c("lwr.ci", "upr.ci")), drop = FALSE]
  }

  condition_rows <- preview$time != "common"
  common_rows <- preview$time == "common"
  overall_preview <- .mc_cia_rm_overall_df(x)
  overall_common_preview <- .mc_cia_rm_overall_common_df(x)
  if (!is.null(overall_preview)) {
    overall_preview$cia <- round(overall_preview$cia, digits)
    for (nm in intersect(c("lwr.ci", "upr.ci"), names(overall_preview))) {
      overall_preview[[nm]] <- round(overall_preview[[nm]], digits)
    }
    if (!has_ci || identical(cfg$show_ci, "no")) {
      overall_preview <- overall_preview[, setdiff(names(overall_preview), c("lwr.ci", "upr.ci")), drop = FALSE]
    }
  }
  if (!is.null(overall_common_preview)) {
    overall_common_preview$cia <- round(overall_common_preview$cia, digits)
    for (nm in intersect(c("lwr.ci", "upr.ci"), names(overall_common_preview))) {
      overall_common_preview[[nm]] <- round(overall_common_preview[[nm]], digits)
    }
    if (!has_ci || identical(cfg$show_ci, "no")) {
      overall_common_preview <- overall_common_preview[, setdiff(names(overall_common_preview), c("lwr.ci", "upr.ci")), drop = FALSE]
    }
  }

  if (any(condition_rows)) {
    cat("Condition-specific CIA estimates\n\n")
    do.call(
      .mc_print_preview_table,
      c(
        list(
          df = preview[condition_rows, , drop = FALSE],
          n = cfg$n,
          topn = cfg$topn,
          max_vars = cfg$max_vars,
          width = cfg$width,
          context = "print",
          full_hint = TRUE,
          summary_hint = TRUE
        ),
        legacy$dots
      )
    )
  }
  if (any(common_rows)) {
    cat("\nCommon CIA estimates\n\n")
    do.call(
      .mc_print_preview_table,
      c(
        list(
          df = preview[common_rows, , drop = FALSE],
          n = cfg$n,
          topn = cfg$topn,
          max_vars = cfg$max_vars,
          width = cfg$width,
          context = "print",
          full_hint = TRUE,
          summary_hint = TRUE
        ),
        legacy$dots
      )
    )
  }
  if (!is.null(overall_preview) && nrow(overall_preview)) {
    cat("\nOverall multi-method CIA\n\n")
    do.call(
      .mc_print_preview_table,
      c(
        list(
          df = overall_preview,
          n = cfg$n,
          topn = cfg$topn,
          max_vars = cfg$max_vars,
          width = cfg$width,
          context = "print",
          full_hint = TRUE,
          summary_hint = TRUE
        ),
        legacy$dots
      )
    )
  }
  if (!is.null(overall_common_preview) && nrow(overall_common_preview)) {
    cat("\nCommon overall CIA\n\n")
    do.call(
      .mc_print_preview_table,
      c(
        list(
          df = overall_common_preview,
          n = cfg$n,
          topn = cfg$topn,
          max_vars = cfg$max_vars,
          width = cfg$width,
          context = "print",
          full_hint = TRUE,
          summary_hint = TRUE
        ),
        legacy$dots
      )
    )
  }
  invisible(x)
}

#' @rdname cia_rm
#' @method print summary.cia_rm
#' @param x A `summary.cia_rm` object.
#' @export
print.summary.cia_rm <- function(x,
                                 digits = NULL,
                                 n = NULL,
                                 topn = NULL,
                                 max_vars = NULL,
                                 width = NULL,
                                 show_ci = NULL,
                                 ...) {
  show_ci <- .mc_resolve_show_ci(show_ci, context = "summary")
  has_ci <- isTRUE(attr(x, "has_ci"))
  cl <- suppressWarnings(as.numeric(attr(x, "conf.level", exact = TRUE)))
  if (length(cl) != 1L || !is.finite(cl)) cl <- NA_real_

  legacy <- .mc_extract_legacy_display_args(
    list(...),
    n = n,
    topn = topn,
    max_vars = max_vars
  )
  cfg <- .mc_resolve_display_args(
    context = "summary",
    n = legacy$n,
    topn = legacy$topn,
    max_vars = legacy$max_vars,
    width = width,
    show_ci = show_ci
  )

  cat(
    "\n",
    .mc_header_with_ci(
      "Repeated-measures coefficient of individual agreement",
      cl,
      if (has_ci) cfg$show_ci else "no"
    ),
    "\n\n",
    sep = ""
  )

  estimates <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  keep_est <- c("method_1", "method_2", "time", "cia", "se", "lwr.ci", "upr.ci")
  estimates <- estimates[, intersect(keep_est, names(estimates)), drop = FALSE]
  if (!has_ci || identical(cfg$show_ci, "no")) {
    estimates <- estimates[, setdiff(names(estimates), c("lwr.ci", "upr.ci")), drop = FALSE]
  }

  diagnostics <- unique(
    as.data.frame(
      x[, c(
        "method_1", "method_2",
        "sigma2_error", "sigma2_subject_method",
        "repeatability", "homogeneity_F", "homogeneity_p"
      ), drop = FALSE],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )
  overall_estimates <- attr(x, "overall_section", exact = TRUE)
  overall_common_estimates <- attr(x, "overall_common_section", exact = TRUE)

  cat("CIA estimates\n\n")
  do.call(
    .mc_print_preview_table,
    c(
      list(
        df = estimates,
        n = cfg$n,
        topn = cfg$topn,
        max_vars = cfg$max_vars,
        width = cfg$width,
        context = "summary",
        full_hint = FALSE,
        summary_hint = FALSE
      ),
      legacy$dots
    )
  )

  if (!is.null(overall_estimates) && nrow(overall_estimates)) {
    cat("\nOverall multi-method CIA\n\n")
    do.call(
      .mc_print_preview_table,
      c(
        list(
          df = overall_estimates,
          n = cfg$n,
          topn = cfg$topn,
          max_vars = cfg$max_vars,
          width = cfg$width,
          context = "summary",
          full_hint = FALSE,
          summary_hint = FALSE
        ),
        legacy$dots
      )
    )
  }

  if (!is.null(overall_common_estimates) && nrow(overall_common_estimates)) {
    cat("\nCommon overall CIA\n\n")
    do.call(
      .mc_print_preview_table,
      c(
        list(
          df = overall_common_estimates,
          n = cfg$n,
          topn = cfg$topn,
          max_vars = cfg$max_vars,
          width = cfg$width,
          context = "summary",
          full_hint = FALSE,
          summary_hint = FALSE
        ),
        legacy$dots
      )
    )
  }

  if (nrow(diagnostics)) {
    cat("\nPair diagnostics\n\n")
    do.call(
      .mc_print_preview_table,
      c(
        list(
          df = diagnostics,
          n = cfg$n,
          topn = cfg$topn,
          max_vars = cfg$max_vars,
          width = cfg$width,
          context = "summary",
          full_hint = FALSE,
          summary_hint = FALSE
        ),
        legacy$dots
      )
    )
  }
  invisible(x)
}

#' @rdname cia_rm
#' @method plot cia_rm
#' @param title Optional plot title.
#' @param facet_by_pair Logical; if `TRUE`, draw one panel per method pair.
#' @param facet_scales Passed to `ggplot2::facet_wrap(scales = ...)` when
#'   `facet_by_pair = TRUE`. One of `"fixed"` or `"free_y"`.
#' @param show_ci Logical; if `TRUE` and confidence intervals are present, show
#'   them on the plot. By default this is inferred from the object class.
#' @param show_common Logical; if `TRUE` and common estimates are available,
#'   draw dashed horizontal reference lines for the pooled/common CIA.
#' @param point_size Numeric point size passed to `ggplot2::geom_point()`.
#' @param line_size Numeric line width passed to `ggplot2::geom_line()`.
#' @param ci_alpha Numeric alpha used for confidence ribbons in continuous
#'   plots.
#' @param ci_linewidth Numeric line width used for confidence-interval
#'   error bars in categorical plots.
#' @param ... Additional arguments passed to `ggplot2::theme()`.
#' @export
plot.cia_rm <- function(x,
                        title = NULL,
                        facet_by_pair = FALSE,
                        facet_scales = c("fixed", "free_y"),
                        show_ci = NULL,
                        show_common = TRUE,
                        point_size = 2.2,
                        line_size = 0.7,
                        ci_alpha = 0.16,
                        ci_linewidth = 0.45,
                        ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} is required for plotting.")
  }
  if (!inherits(x, "cia_rm")) {
    cli::cli_abort("{.arg x} must inherit from {.cls cia_rm}.")
  }
  check_bool(facet_by_pair, arg = "facet_by_pair")
  facet_scales <- match.arg(facet_scales)
  check_bool(show_common, arg = "show_common")
  if (!is.null(show_ci)) {
    check_bool(show_ci, arg = "show_ci")
  }

  if (is.null(show_ci)) {
    show_ci <- inherits(x, "cia_ci")
  }

  df <- .mc_cia_rm_flatten(x)
  df <- df[df$time != "common", , drop = FALSE]
  if (!nrow(df)) {
    cli::cli_abort("No condition-specific CIA rows are available for plotting.")
  }

  df$pair <- paste(df$method_1, df$method_2, sep = " vs ")
  time_levels <- attr(x, "time_levels", exact = TRUE)
  time_scale <- attr(x, "time_scale", exact = TRUE)
  time_values <- attr(x, "time_plot_values", exact = TRUE)
  time_var <- attr(x, "time_var", exact = TRUE)
  if (!length(time_levels) || !length(time_values) || length(time_levels) != length(time_values)) {
    time_levels <- unique(df$time)
    time_values <- time_levels
    time_scale <- "categorical"
  }

  if (identical(time_scale, "continuous")) {
    df$time_value <- time_values[match(df$time, time_levels)]
    ord <- order(df$pair, df$time_value)
    df <- df[ord, , drop = FALSE]
  } else {
    df$time_value <- factor(df$time, levels = time_levels)
  }

  common_df <- df[FALSE, c("pair", "cia"), drop = FALSE]
  if (isTRUE(show_common) && is.matrix(x$common)) {
    common_rows <- .mc_cia_rm_flatten(x)
    common_rows <- common_rows[common_rows$time == "common", , drop = FALSE]
    if (nrow(common_rows)) {
      common_rows$pair <- paste(common_rows$method_1, common_rows$method_2, sep = " vs ")
      common_df <- common_rows[, c("pair", "cia"), drop = FALSE]
    }
  }

  has_ci_cols <- all(c("lwr.ci", "upr.ci") %in% names(df))
  show_ci <- isTRUE(show_ci) && has_ci_cols
  if (is.null(title)) {
    title <- if (identical(time_scale, "continuous")) {
      "Repeated-measures CIA across condition levels"
    } else {
      "Repeated-measures CIA by condition"
    }
  }

  subtitle <- if (identical(time_scale, "continuous")) {
    "Condition scale detected as continuous/discrete"
  } else {
    "Condition scale detected as categorical"
  }

  if (identical(time_scale, "continuous")) {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data$time_value,
        y = .data$cia,
        colour = .data$pair,
        group = .data$pair
      )
    )
    if (show_ci) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$lwr.ci, ymax = .data$upr.ci, fill = .data$pair),
        alpha = ci_alpha,
        colour = NA,
        inherit.aes = TRUE
      )
    }
    if (nrow(common_df)) {
      p <- p + ggplot2::geom_hline(
        data = common_df,
        ggplot2::aes(yintercept = .data$cia, colour = .data$pair),
        linetype = "dashed",
        linewidth = 0.45,
        inherit.aes = FALSE
      )
    }
    p <- p +
      ggplot2::geom_line(linewidth = line_size) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        ...
      ) +
      ggplot2::labs(
        title = title,
        subtitle = subtitle,
        x = time_var,
        y = "CIA",
        colour = "Method pair",
        fill = "Method pair"
      )

    if (isTRUE(facet_by_pair)) {
      p <- p +
        ggplot2::facet_wrap(ggplot2::vars(.data$pair), scales = facet_scales) +
        ggplot2::guides(colour = "none", fill = "none")
    }
    return(p)
  }

  dodge <- if (isTRUE(facet_by_pair)) {
    ggplot2::position_identity()
  } else {
    ggplot2::position_dodge(width = 0.45)
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$time_value,
      y = .data$cia,
      colour = .data$pair
    )
  )
  if (show_ci) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$lwr.ci, ymax = .data$upr.ci),
      width = 0.15,
      linewidth = ci_linewidth,
      position = dodge
    )
  }
  if (nrow(common_df)) {
    p <- p + ggplot2::geom_hline(
      data = common_df,
      ggplot2::aes(yintercept = .data$cia, colour = .data$pair),
      linetype = "dashed",
      linewidth = 0.45,
      inherit.aes = FALSE
    )
  }
  p <- p +
    ggplot2::geom_point(size = point_size, position = dodge) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
      ...
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = time_var,
      y = "CIA",
      colour = "Method pair"
    )

  if (isTRUE(facet_by_pair)) {
    p <- p +
      ggplot2::facet_wrap(ggplot2::vars(.data$pair), scales = facet_scales) +
      ggplot2::guides(colour = "none")
  }
  p
}

