#' @title Bland-Altman statistics with confidence intervals
#'
#' @description
#' Computes Bland-Altman mean difference and limits of agreement (LoA)
#' between two numeric measurement vectors, including t-based confidence
#' intervals for the mean difference and for each LoA using 'C++' backend.
#' If `group2` is omitted and `group1` is a numeric matrix or data frame with
#' at least two numeric columns, `ba()` computes all pairwise contrasts across
#' methods and returns a pairwise Bland-Altman matrix object.
#'
#' Note: Lin's concordance correlation coefficient (CCC) is a complementary,
#' single-number summary of agreement (precision + accuracy). It is useful for
#' quick screening or reporting an overall CI, but may miss systematic or
#' magnitude-dependent bias; consider reporting CCC alongside Bland-Altman.
#'
#' @details
#' Given paired measurements \eqn{(x_i, y_i)}, Bland-Altman analysis uses
#' \eqn{d_i = x_i - y_i} (or \eqn{y_i - x_i} if \code{mode = 2}) and
#' \eqn{m_i = (x_i + y_i)/2}. The mean difference \eqn{\bar d} estimates bias.
#' The limits of agreement (LoA) are \eqn{\bar d \pm z \cdot s_d}, where
#' \eqn{s_d} is the sample standard deviation of \eqn{d_i} and \eqn{z}
#' (argument \code{loa_multiplier}) is typically 1.96 for nominal 95% LoA.
#'
#' When `group2` is omitted and `group1` is a wide numeric matrix or data
#' frame, the same two-method calculation is applied to every unordered column
#' pair. The returned pairwise matrices follow the requested `mode`: with
#' `mode = 1`, upper-triangle entries represent row minus column; with
#' `mode = 2`, upper-triangle entries represent column minus row.
#'
#' Confidence intervals use Student's \eqn{t} distribution with \eqn{n-1}
#' degrees of freedom, with
#' \itemize{
#'   \item Mean-difference CI given by \eqn{\bar d \pm t_{n-1,\,1-\alpha/2}\,
#'   s_d/\sqrt{n}}; and
#'   \item LoA CI given by \eqn{(\bar d \pm z\, s_d) \;\pm\;
#'   t_{n-1,\,1-\alpha/2}\, s_d\,\sqrt{3/n}}.
#' }
#'
#' Assumptions include approximately normal differences and roughly constant
#' variability across the measurement range; if differences increase with
#' magnitude, consider a transformation before analysis. Missing values are
#' removed pairwise (rows with an \code{NA} in either input are dropped before
#' calling the C++ backend).
#'
#' Probability of agreement, available through \code{\link{prob_agree}}, is a
#' tolerance-based companion to Bland-Altman analysis. Bland-Altman reports bias
#' and limits of agreement for paired differences. `prob_agree()` instead uses
#' the sampling distribution of estimated differences to quantify the
#' probability that two estimated quantities or curves agree within a
#' user-specified practical tolerance.
#'
#' @param group1 Numeric vector of paired measurements, or a numeric matrix/data
#'   frame with at least two numeric columns when `group2` is omitted.
#' @param group2 Optional numeric vector of paired measurements. If omitted,
#'   `group1` must be a numeric matrix or data frame and all unordered numeric
#'   column pairs are analysed.
#' @param loa_multiplier Positive scalar; the multiple of the standard deviation used to
#'   define the LoA (default 1.96 for nominal 95\% agreement). The confidence
#'   intervals always use \eqn{t_{n-1,\,1-\alpha/2}} regardless of this choice.
#' @param mode Integer; 1 uses \code{group1 - group2}, 2 uses \code{group2 - group1}.
#' @param conf_level Confidence level for CIs (default 0.95).
#' @param n_threads Integer \eqn{\geq 1}. Number of OpenMP threads. Defaults to
#'   \code{getOption("matrixCorr.threads", 1L)}.
#' @param verbose Logical; if TRUE, prints how many OpenMP threads are used.
#'
#' @return
#' If `group1` and `group2` are both supplied, an object of class
#' \code{"ba"} (list) with elements:
#' \itemize{
#'   \item \code{means}, \code{diffs}: numeric vectors
#'   \item \code{groups}: data.frame used after NA removal
#'   \item \code{n_obs}: integer, number of complete pairs used.
#'   \item \code{based.on}: compatibility alias for \code{n_obs}.
#'   \item \code{lower.limit}, \code{mean.diffs}, \code{upper.limit}
#'   \item \code{lines}: named numeric vector (lower, mean, upper)
#'   \item \code{CI.lines}: named numeric vector for CIs of those lines
#'   \item \code{loa_multiplier}, \code{critical.diff}
#' }
#'
#' If `group2` is omitted and `group1` is a wide numeric matrix or data frame,
#' an object of class \code{"ba_matrix"} with pairwise components:
#' \itemize{
#'   \item \code{bias}: pairwise matrix of Bland-Altman mean differences.
#'   \item \code{sd_loa}: pairwise matrix of SDs of differences.
#'   \item \code{loa_lower}, \code{loa_upper}: pairwise matrices of LoA
#'   endpoints.
#'   \item \code{width}: pairwise matrix of LoA widths.
#'   \item \code{n}: integer pairwise matrix of complete-case counts.
#'   \item \code{mean_ci_low}, \code{mean_ci_high}: pairwise matrices of CI
#'   bounds for the mean difference.
#'   \item \code{loa_lower_ci_low}, \code{loa_lower_ci_high}: pairwise matrices
#'   of CI bounds for the lower LoA.
#'   \item \code{loa_upper_ci_low}, \code{loa_upper_ci_high}: pairwise matrices
#'   of CI bounds for the upper LoA.
#'   \item \code{methods}: character vector of analysed method names.
#'   \item \code{loa_multiplier}, \code{mode}: calculation settings reused for
#'   every pair.
#' }
#'
#' @seealso \code{\link{print.ba}}, \code{\link{plot.ba}},
#'  \code{\link{ccc}}, \code{\link{prob_agree}}, \code{\link{ccc_rm_ustat}},
#'  \code{\link{ccc_rm_reml}}
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100, 100, 10)
#' y <- x + rnorm(100, 0, 8)
#' fit_ba <- ba(x, y)
#' print(fit_ba)
#' plot(fit_ba)
#'
#' # Pairwise Bland-Altman across 3 methods
#' set.seed(7)
#' wide3 <- data.frame(
#'   ref = rnorm(80, 100, 8),
#'   m2  = rnorm(80, 101, 8),
#'   m3  = rnorm(80,  99, 9)
#' )
#' fit_ba3 <- ba(wide3)
#' print(fit_ba3)
#' summary(fit_ba3)
#' plot(fit_ba3)
#'
#' @references
#' Bland JM, Altman DG (1986). Statistical methods for assessing agreement
#' between two methods of clinical measurement. *The Lancet*, 307-310.
#' @references
#' Bland JM, Altman DG (1999). Measuring agreement in method comparison studies.
#' *Statistical Methods in Medical Research*, 8(2), 135-160.
#'
#' @author Thiago de Paula Oliveira
#'
#' @export
ba <- function(group1,
               group2,
               loa_multiplier = 1.96,
               mode = 1L,
               conf_level = 0.95,
               n_threads = getOption("matrixCorr.threads", 1L),
               verbose = FALSE) {
  if (missing(group2)) {
    return(.mc_ba_fit_matrix(
      data = group1,
      loa_multiplier = loa_multiplier,
      mode = mode,
      conf_level = conf_level,
      n_threads = n_threads,
      verbose = verbose
    ))
  }

  .mc_ba_fit_two(
    group1 = group1,
    group2 = group2,
    loa_multiplier = loa_multiplier,
    mode = mode,
    conf_level = conf_level,
    n_threads = n_threads,
    verbose = verbose
  )
}

.mc_ba_fit_two <- function(group1,
                           group2,
                           loa_multiplier,
                           mode,
                           conf_level,
                           n_threads,
                           verbose) {
  # -- validate ---------------------------------------------------------------
  if (!is.numeric(group1) || !is.numeric(group2)) {
    abort_bad_arg("group1",
      message = "and {.arg group2} must be numeric vectors."
    )
  }
  check_same_length(group1, group2, arg_x = "group1", arg_y = "group2")
  if (!is.numeric(loa_multiplier) || length(loa_multiplier) != 1L || !is.finite(loa_multiplier) || loa_multiplier <= 0) {
    abort_bad_arg("loa_multiplier",
      message = "`{arg}` must be a positive scalar."
    )
  }
  if (!rlang::is_scalar_integerish(mode) || !(as.integer(mode) %in% c(1L, 2L))) {
    abort_bad_arg("mode",
      message = "must be 1 or 2."
    )
  }
  mode <- as.integer(mode)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    abort_bad_arg("conf_level",
      message = "`conf_level` must be in (0, 1)."
    )
  }
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")
  check_bool(verbose, arg = "verbose")

  called.with <- length(group1)
  complete_pairs <- !is.na(group1) & !is.na(group2)
  n_pairs <- sum(complete_pairs)
  if (n_pairs < 2L) {
    abort_bad_arg("group1",
      message = "must provide at least two complete pairs with {.arg group2} after removing missing values (found {n_pairs}).",
      n_pairs = n_pairs
    )
  }
  inform_if_verbose(
    "Using {n_threads} OpenMP thread{?s}.",
    .verbose = verbose
  )

  # -- compute in C++ ---------------------------------------------------------
  prev_threads <- get_omp_threads()
  on.exit(set_omp_threads(as.integer(prev_threads)), add = TRUE)
  ba_out <- bland_altman_cpp(
    group1 = group1,
    group2 = group2,
    loa_multiplier = loa_multiplier,
    mode = mode,
    conf_level = conf_level,
    n_threads = n_threads,
    include_legacy_payload = FALSE
  )
  ba_out <- structure(
    list(
      means = as.numeric(ba_out$means),
      diffs = as.numeric(ba_out$diffs),
      n_obs = as.integer(ba_out$based.on),
      lower.limit = as.numeric(ba_out$lower.limit),
      mean.diffs = as.numeric(ba_out$mean.diffs),
      upper.limit = as.numeric(ba_out$upper.limit),
      loa_multiplier = as.numeric(ba_out$loa_multiplier),
      critical.diff = as.numeric(ba_out$critical.diff),
      ci = .mc_ba_normalize_ci_lines(ba_out$CI.lines),
      mode = mode
    ),
    class = "ba"
  )
  attr(ba_out, "method") <- "Bland-Altman"
  attr(ba_out, "description") <- "Mean difference and limits of agreement with CIs"
  attr(ba_out, "package") <- "matrixCorr"
  attr(ba_out, "conf.level") <- conf_level
  attr(ba_out, "called.with") <- length(group1)
  ba_out
}

.mc_ba_fit_matrix <- function(data,
                              loa_multiplier,
                              mode,
                              conf_level,
                              n_threads,
                              verbose) {
  if (is.data.frame(data)) {
    is_num <- vapply(data, is.numeric, logical(1))
    if (!any(is_num)) {
      abort_bad_arg("group1",
        message = "must contain at least two numeric columns when {.arg group2} is omitted."
      )
    }
    data <- data[is_num]
  }

  if (!is.matrix(data) && !is.data.frame(data)) {
    abort_bad_arg("group1",
      message = "must be a numeric matrix or data frame when {.arg group2} is omitted."
    )
  }
  if (is.matrix(data) && !is.numeric(data)) {
    abort_bad_arg("group1",
      message = "must be a numeric matrix when {.arg group2} is omitted."
    )
  }

  data_wide <- as.data.frame(data, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(data_wide) < 2L) {
    abort_bad_arg("group1",
      message = "must contain at least two numeric columns when {.arg group2} is omitted."
    )
  }
  methods <- names(data_wide)
  if (is.null(methods) || anyNA(methods) || any(methods == "")) {
    methods <- paste0("method", seq_len(ncol(data_wide)))
    names(data_wide) <- methods
  }

  if (!is.numeric(loa_multiplier) || length(loa_multiplier) != 1L || !is.finite(loa_multiplier) || loa_multiplier <= 0) {
    abort_bad_arg("loa_multiplier",
      message = "`{arg}` must be a positive scalar."
    )
  }
  if (!rlang::is_scalar_integerish(mode) || !(as.integer(mode) %in% c(1L, 2L))) {
    abort_bad_arg("mode",
      message = "must be 1 or 2."
    )
  }
  mode <- as.integer(mode)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    abort_bad_arg("conf_level",
      message = "`conf_level` must be in (0, 1)."
    )
  }
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")
  check_bool(verbose, arg = "verbose")

  if (isTRUE(verbose)) {
    inform_if_verbose(
      "Pairwise Bland-Altman: {ncol(data_wide)} methods -> {choose(ncol(data_wide), 2L)} pairs.",
      .verbose = TRUE
    )
  }

  mm <- length(methods)
  bias <- .make_named_matrix_ba(methods)
  sd_loa <- .make_named_matrix_ba(methods)
  loa_lower <- .make_named_matrix_ba(methods)
  loa_upper <- .make_named_matrix_ba(methods)
  width <- .make_named_matrix_ba(methods)
  n_mat <- .make_named_matrix_ba(methods, NA_integer_, "integer")
  mean_ci_low <- .make_named_matrix_ba(methods)
  mean_ci_high <- .make_named_matrix_ba(methods)
  lo_ci_low <- .make_named_matrix_ba(methods)
  lo_ci_high <- .make_named_matrix_ba(methods)
  hi_ci_low <- .make_named_matrix_ba(methods)
  hi_ci_high <- .make_named_matrix_ba(methods)

  for (i in 1:(mm - 1L)) for (j in (i + 1L):mm) {
    keep <- !is.na(data_wide[[i]]) & !is.na(data_wide[[j]])
    n_pair <- sum(keep)
    n_mat[i, j] <- n_pair
    n_mat[j, i] <- n_pair
    if (n_pair < 2L) {
      next
    }

    fit <- .mc_ba_fit_two(
      group1 = data_wide[[i]],
      group2 = data_wide[[j]],
      loa_multiplier = loa_multiplier,
      mode = mode,
      conf_level = conf_level,
      n_threads = n_threads,
      verbose = FALSE
    )

    md <- as.numeric(fit$mean.diffs)
    sdv <- as.numeric(fit$critical.diff) / as.numeric(fit$loa_multiplier)
    lo <- as.numeric(fit$lower.limit)
    hi <- as.numeric(fit$upper.limit)

    bias[i, j] <- md
    bias[j, i] <- -md
    sd_loa[i, j] <- sdv
    sd_loa[j, i] <- sdv
    loa_lower[i, j] <- lo
    loa_lower[j, i] <- -hi
    loa_upper[i, j] <- hi
    loa_upper[j, i] <- -lo
    width[i, j] <- hi - lo
    width[j, i] <- hi - lo

    mean_ci_low[i, j] <- as.numeric(fit$CI.lines[["mean.diff.ci.lower"]])
    mean_ci_high[i, j] <- as.numeric(fit$CI.lines[["mean.diff.ci.upper"]])
    lo_ci_low[i, j] <- as.numeric(fit$CI.lines[["lower.limit.ci.lower"]])
    lo_ci_high[i, j] <- as.numeric(fit$CI.lines[["lower.limit.ci.upper"]])
    hi_ci_low[i, j] <- as.numeric(fit$CI.lines[["upper.limit.ci.lower"]])
    hi_ci_high[i, j] <- as.numeric(fit$CI.lines[["upper.limit.ci.upper"]])

    mean_ci_low[j, i] <- -mean_ci_high[i, j]
    mean_ci_high[j, i] <- -mean_ci_low[i, j]
    lo_ci_low[j, i] <- -hi_ci_high[i, j]
    lo_ci_high[j, i] <- -hi_ci_low[i, j]
    hi_ci_low[j, i] <- -lo_ci_high[i, j]
    hi_ci_high[j, i] <- -lo_ci_low[i, j]
  }

  out <- structure(
    list(
      bias = bias,
      sd_loa = sd_loa,
      loa_lower = loa_lower,
      loa_upper = loa_upper,
      width = width,
      n = n_mat,
      mean_ci_low = mean_ci_low,
      mean_ci_high = mean_ci_high,
      loa_lower_ci_low = lo_ci_low,
      loa_lower_ci_high = lo_ci_high,
      loa_upper_ci_low = hi_ci_low,
      loa_upper_ci_high = hi_ci_high,
      methods = methods,
      loa_multiplier = loa_multiplier,
      mode = mode,
      data_wide = data_wide
    ),
    class = c("ba_matrix", "list")
  )
  attr(out, "method") <- "Bland-Altman"
  attr(out, "description") <- "Pairwise mean differences and limits of agreement with CIs"
  attr(out, "package") <- "matrixCorr"
  attr(out, "conf.level") <- conf_level
  out
}

.mc_ba_normalize_ci_lines <- function(x) {
  out <- as.numeric(x)
  nms <- names(x)
  if (is.null(nms) && length(out) == 6L) {
    nms <- c(
      "lower.limit.ci.lower",
      "lower.limit.ci.upper",
      "mean.diff.ci.lower",
      "mean.diff.ci.upper",
      "upper.limit.ci.lower",
      "upper.limit.ci.upper"
    )
  }
  stats::setNames(out, nms)
}

.mc_ba_groups <- function(x) {
  means <- as.numeric(unclass(x)$means)
  diffs <- as.numeric(unclass(x)$diffs)
  mode <- as.integer(unclass(x)$mode %||% 1L)
  if (identical(mode, 2L)) {
    group1 <- means - diffs / 2
    group2 <- means + diffs / 2
  } else {
    group1 <- means + diffs / 2
    group2 <- means - diffs / 2
  }
  data.frame(group1 = group1, group2 = group2, check.names = FALSE)
}

.mc_ba_lines <- function(x) {
  c(
    lower.limit = as.numeric(unclass(x)$lower.limit),
    mean.diffs = as.numeric(unclass(x)$mean.diffs),
    upper.limit = as.numeric(unclass(x)$upper.limit)
  )
}

.mc_ba_ci_lines <- function(x) {
  ci <- .mc_ba_normalize_ci_lines(unclass(x)$ci %||% numeric())
  if (length(ci)) {
    return(ci)
  }
  numeric()
}

#' @export
names.ba <- function(x) {
  unique(c(
    names(unclass(x)),
    "based.on",
    "groups",
    "lines",
    "CI.lines"
  ))
}

#' @export
`$.ba` <- function(x, name) {
  core <- unclass(x)
  switch(name,
    based.on = core$n_obs,
    groups = .mc_ba_groups(x),
    lines = .mc_ba_lines(x),
    CI.lines = .mc_ba_ci_lines(x),
    core[[name]]
  )
}

#' @export
`[[.ba` <- function(x, i, ...) {
  if (is.numeric(i)) {
    nms <- names(x)
    if (length(i) != 1L || is.na(i) || i < 1L || i > length(nms)) {
      abort_bad_arg("i", message = "must be a valid component index.")
    }
    i <- nms[[i]]
  }
  `$.ba`(x, i)
}

#' @rdname ba
#' @method print ba
#' @param x A \code{"ba"} object.
#' @param digits Number of digits for estimates (default 3).
#' @param ci_digits Number of digits for CI bounds (default 3).
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param ... Unused.
#' @export
print.ba <- function(x, digits = 3, ci_digits = 3,
                     n = NULL, topn = NULL,
                     max_vars = NULL, width = NULL,
                     show_ci = NULL, ...) {
  check_inherits(x, "ba")
  show_ci <- .mc_resolve_show_ci(show_ci, context = "print")

  n   <- as.integer(x$based.on)
  loa_multiplier <- as.numeric(x$loa_multiplier)
  cl  <- suppressWarnings(as.numeric(attr(x, "conf.level")))
  if (!is.finite(cl)) cl <- NA_real_

  # core numbers (as scalars)
  bias   <- as.numeric(x$mean.diffs)
  loa_lo <- as.numeric(x$lower.limit)
  loa_hi <- as.numeric(x$upper.limit)
  sd_d   <- as.numeric(x$critical.diff) / loa_multiplier
  loa_width <- loa_hi - loa_lo

  # CIs (robust extraction by name)
  cil <- function(nm) as.numeric(x$CI.lines[[nm]])
  bias_l <- cil("mean.diff.ci.lower"); bias_u <- cil("mean.diff.ci.upper")
  lo_l   <- cil("lower.limit.ci.lower"); lo_u <- cil("lower.limit.ci.upper")
  hi_l   <- cil("upper.limit.ci.lower"); hi_u <- cil("upper.limit.ci.upper")

  .mc_print_named_digest(
    c(
      based_on = n,
      loa_rule = sprintf("mean +/- %.3g * SD", loa_multiplier),
      if (identical(show_ci, "yes") && is.finite(cl)) c(ci = sprintf("%g%%", 100 * cl)),
      sd_diff = formatC(sd_d, format = "f", digits = digits),
      width = formatC(loa_width, format = "f", digits = digits)
    ),
    header = "Bland-Altman preview:"
  )
  cat("\n")

  # nicely aligned three-row table
  df <- data.frame(
    quantity = c("Mean difference", "Lower LoA", "Upper LoA"),
    estimate = c(bias, loa_lo, loa_hi),
    lwr      = c(bias_l, lo_l, hi_l),
    upr      = c(bias_u, lo_u, hi_u),
    check.names = FALSE,
    row.names   = NULL
  )
  df$estimate <- formatC(df$estimate, format = "f", digits = digits)
  if (identical(show_ci, "yes")) {
    df$lwr <- formatC(df$lwr, format = "f", digits = ci_digits)
    df$upr <- formatC(df$upr, format = "f", digits = ci_digits)
  } else {
    df$lwr <- NULL
    df$upr <- NULL
  }

  .mc_print_preview_table(
    df,
    n = .mc_coalesce(n, .mc_display_option("print_max_rows", 20L)),
    topn = .mc_coalesce(topn, .mc_display_option("print_topn", 5L)),
    max_vars = .mc_coalesce(max_vars, .mc_display_option("print_max_vars", NULL)),
    width = .mc_coalesce(width, getOption("width", 80L)),
    context = "print",
    full_hint = TRUE,
    summary_hint = TRUE,
    ...
  )
  invisible(x)
}

#' @rdname ba
#' @method summary ba
#' @param object A \code{"ba"} object.
#' @export
summary.ba <- function(object,
                       digits = 3,
                       ci_digits = 3,
                       ...) {
  check_inherits(object, "ba")

  cil <- function(nm) as.numeric(object$CI.lines[[nm]])

  out <- data.frame(
    n_obs          = as.integer(object$based.on),
    bias           = round(as.numeric(object$mean.diffs), digits),
    sd_loa         = round(as.numeric(object$critical.diff) / as.numeric(object$loa_multiplier), digits),
    loa_low        = round(as.numeric(object$lower.limit), digits),
    loa_up         = round(as.numeric(object$upper.limit), digits),
    width          = round(as.numeric(object$upper.limit - object$lower.limit), digits),
    loa_multiplier = round(as.numeric(object$loa_multiplier), digits),
    bias_lwr       = round(cil("mean.diff.ci.lower"), ci_digits),
    bias_upr       = round(cil("mean.diff.ci.upper"), ci_digits),
    lo_lwr         = round(cil("lower.limit.ci.lower"), ci_digits),
    lo_upr         = round(cil("lower.limit.ci.upper"), ci_digits),
    up_lwr         = round(cil("upper.limit.ci.lower"), ci_digits),
    up_upr         = round(cil("upper.limit.ci.upper"), ci_digits),
    check.names = FALSE
  )

  out <- .mc_finalize_summary_df(out, class_name = "summary.ba")
  attr(out, "conf.level") <- suppressWarnings(as.numeric(attr(object, "conf.level")))
  attr(out, "digits") <- digits
  attr(out, "ci_digits") <- ci_digits
  out
}

#' @rdname ba
#' @method print summary.ba
#' @param ... Passed to \code{\link[base]{print.data.frame}}.
#' @export
print.summary.ba <- function(x, digits = NULL, n = NULL,
                             topn = NULL, max_vars = NULL,
                             width = NULL, show_ci = NULL, ...) {
  cl <- suppressWarnings(as.numeric(attr(x, "conf.level")))
  show_ci <- .mc_resolve_show_ci(show_ci, context = "summary")
  header <- .mc_header_with_ci("Bland-Altman (two methods)", cl, show_ci)
  .mc_print_sectioned_table(
    x,
    sections = list(
      list(
        title = "Agreement estimates",
        cols = c("n_obs", "bias", "sd_loa", "loa_low", "loa_up", "width", "loa_multiplier")
      ),
      list(
        title = "Confidence intervals",
        cols = c("bias_lwr", "bias_upr", "lo_lwr", "lo_upr", "up_lwr", "up_upr")
      )
    ),
    header = header,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    print_overview = FALSE,
    ...
  )
  invisible(x)
}

#' @rdname ba
#' @method plot ba
#' @param x A \code{"ba"} object.
#' @param title Plot title.
#' @param subtitle Optional subtitle. If NULL, shows n and LoA summary.
#' @param point_alpha Point transparency.
#' @param point_size Point size.
#' @param line_size Line width for mean/LoA.
#' @param shade_ci Logical; if TRUE, draw shaded CI bands instead of 6 dashed
#' lines.
#' @param shade_alpha Transparency of CI bands.
#' @param smoother One of "none", "loess", "lm" to visualize proportional bias.
#' @param symmetrize_y Logical; if TRUE, y-axis centered at mean difference
#' with symmetric limits.
#' @param show_value Logical; included for a consistent plotting interface.
#'   Bland-Altman plots do not overlay numeric cell values, so this argument
#'   currently has no effect.
#' @param ... Passed to \code{ggplot2::theme()} (ggplot path) or \code{plot()}.
#' @importFrom graphics abline lines par rect
#' @export
plot.ba <- function(x,
                    title = "Bland-Altman Plot",
                    subtitle = NULL,
                    point_alpha = 0.7,
                    point_size  = 2.2,
                    line_size   = 0.8,
                    shade_ci    = TRUE,
                    shade_alpha = 0.08,
                    smoother    = c("none", "loess", "lm"),
                    symmetrize_y = TRUE,
                    show_value = TRUE,
                    ...) {
  check_inherits(x, "ba")
  check_bool(show_value, arg = "show_value")
  smoother <- match.arg(smoother)

  means <- as.numeric(x$means)
  diffs <- as.numeric(x$diffs)

  # scalars
  md    <- as.numeric(x$mean.diffs)
  loaL  <- as.numeric(x$lower.limit)
  loaU  <- as.numeric(x$upper.limit)
  loa_multiplier <- as.numeric(x$loa_multiplier)
  n     <- as.integer(x$based.on)
  sd_d  <- as.numeric(x$critical.diff) / loa_multiplier
  cl    <- suppressWarnings(as.numeric(attr(x, "conf.level")))
  ci    <- function(nm) as.numeric(x$CI.lines[[nm]])

  if (is.null(subtitle)) {
    subtitle <- if (is.finite(cl)) {
      sprintf("n = %d; mean diff = %.2f; LoA = [%.2f, %.2f]; %g%% CI shown",
              n, md, loaL, loaU, 100*cl)
    } else {
      sprintf("n = %d; mean diff = %.2f; LoA = [%.2f, %.2f]",
              n, md, loaL, loaU)
    }
  }

  # y limits (symmetric around md if requested)
  y_rng <- range(c(diffs, loaL, loaU), na.rm = TRUE)
  if (isTRUE(symmetrize_y)) {
    half <- max(abs(y_rng - md))
    y_rng <- c(md - half, md + half)
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    ## ---------- Base fallback ----------
    plot(means, diffs,
         xlab = "Mean of methods", ylab = "Difference between methods",
         main = title, sub = subtitle,
         pch = 16, cex = point_size / 2.2,
         col = grDevices::gray(0.2, alpha = point_alpha),
         ylim = y_rng, ...)
    # shaded CI bands (rectangles)
    if (shade_ci) {
      usr <- par("usr")
      rect(xleft = usr[1], xright = usr[2],
           ybottom = ci("mean.diff.ci.lower"), ytop = ci("mean.diff.ci.upper"),
           border = NA, col = grDevices::gray(0.2, alpha = shade_alpha))
      rect(xleft = usr[1], xright = usr[2],
           ybottom = ci("lower.limit.ci.lower"), ytop = ci("lower.limit.ci.upper"),
           border = NA, col = grDevices::gray(0.2, alpha = shade_alpha))
      rect(xleft = usr[1], xright = usr[2],
           ybottom = ci("upper.limit.ci.lower"), ytop = ci("upper.limit.ci.upper"),
           border = NA, col = grDevices::gray(0.2, alpha = shade_alpha))
    } else {
      abline(h = c(ci("mean.diff.ci.lower"), ci("mean.diff.ci.upper")), lty = 2)
      abline(h = c(ci("lower.limit.ci.lower"), ci("lower.limit.ci.upper")), lty = 2)
      abline(h = c(ci("upper.limit.ci.lower"), ci("upper.limit.ci.upper")), lty = 2)
    }
    # reference & main lines
    abline(h = 0,   col = "grey70", lty = 3)
    abline(h = md,  lwd = line_size * 1.4)
    abline(h = loaL, lwd = line_size * 1.2)
    abline(h = loaU, lwd = line_size * 1.2)

    # optional smoother
    if (smoother != "none") {
      fit <- if (smoother == "lm") stats::lm(diffs ~ means) else stats::loess(diffs ~ means)
      xs <- seq(min(means), max(means), length.out = 200)
      ys <- stats::predict(fit, newdata = data.frame(means = xs))
      lines(xs, ys, lty = 1, col = "grey40")
    }
    return(invisible(NULL))
  }

  ## ---------- ggplot path ----------
  df <- data.frame(means = means, diffs = diffs)
  xr <- range(means, na.rm = TRUE)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = means, y = diffs))

  if (shade_ci) {
    p <- p +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                        ymin = ci("mean.diff.ci.lower"), ymax = ci("mean.diff.ci.upper"),
                        alpha = shade_alpha) +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                        ymin = ci("lower.limit.ci.lower"), ymax = ci("lower.limit.ci.upper"),
                        alpha = shade_alpha) +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                        ymin = ci("upper.limit.ci.lower"), ymax = ci("upper.limit.ci.upper"),
                        alpha = shade_alpha)
  } else {
    p <- p +
      ggplot2::geom_hline(yintercept = c(ci("mean.diff.ci.lower"),
                                         ci("mean.diff.ci.upper"),
                                         ci("lower.limit.ci.lower"),
                                         ci("lower.limit.ci.upper"),
                                         ci("upper.limit.ci.lower"),
                                         ci("upper.limit.ci.upper")),
                          linetype = "dashed")
  }

  p <- p +
    ggplot2::geom_point(alpha = point_alpha, size = point_size) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dotted", color = "grey40") +
    ggplot2::geom_hline(yintercept = md,   linewidth = line_size) +
    ggplot2::geom_hline(yintercept = loaL, linewidth = line_size) +
    ggplot2::geom_hline(yintercept = loaU, linewidth = line_size)

  if (smoother == "lm") {
    p <- p + ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.7)
  } else if (smoother == "loess") {
    p <- p + ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 0.7, span = 0.9)
  }

  p <- p +
    ggplot2::coord_cartesian(ylim = y_rng, expand = TRUE) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), ...) +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = "Mean of methods", y = "Difference between methods")

  p
}

#' @rdname ba
#' @method print ba_matrix
#' @param x A \code{"ba_matrix"} object.
#' @param style Show the pairwise result as \code{"pairs"} or \code{"matrices"}.
#' @export
print.ba_matrix <- function(x,
                            digits = 3,
                            ci_digits = 3,
                            n = NULL,
                            topn = NULL,
                            max_vars = NULL,
                            width = NULL,
                            show_ci = NULL,
                            style = c("pairs", "matrices"),
                            ...) {
  check_inherits(x, "ba_matrix")
  style <- match.arg(style)
  show_ci <- .mc_resolve_show_ci(show_ci, context = "print")
  cl <- suppressWarnings(as.numeric(attr(x, "conf.level")))

  if (style == "matrices") {
    .mc_print_named_digest(
      c(
        methods = length(x$methods),
        mode = if (identical(as.integer(x$mode), 1L)) "row - column" else "column - row",
        if (identical(show_ci, "yes") && is.finite(cl)) c(ci = sprintf("%g%%", 100 * cl)),
        components = "bias, sd_loa, loa_low, loa_up, width, n"
      ),
      header = "Bland-Altman matrix preview:"
    )
    cat("\nPrimary component: bias\n\n")
    .mc_print_preview_matrix(
      x$bias,
      digits = digits,
      n = .mc_coalesce(n, .mc_display_option("print_max_rows", 20L)),
      topn = .mc_coalesce(topn, .mc_display_option("print_topn", 5L)),
      max_vars = .mc_coalesce(max_vars, .mc_display_option("print_max_vars", NULL)),
      width = .mc_coalesce(width, getOption("width", 80L)),
      context = "print",
      full_hint = TRUE,
      summary_hint = TRUE,
      ...
    )
    return(invisible(x))
  }

  sm <- summary(x, digits = digits, ci_digits = ci_digits)
  cols <- c("method1", "method2", "bias", "sd_loa", "loa_low", "loa_up", "width", "n_obs")
  if (identical(show_ci, "yes")) {
    cols <- c(cols, "bias_lwr", "bias_upr", "lo_lwr", "lo_upr", "up_lwr", "up_upr")
  }
  df <- as.data.frame(sm)[, cols[cols %in% names(sm)], drop = FALSE]
  header <- .mc_header_with_ci(
    sprintf("Bland-Altman (%s)", if (identical(as.integer(x$mode), 1L)) "row - column" else "column - row"),
    cl,
    show_ci
  )
  .mc_print_summary_table(
    df,
    header = header,
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    print_overview = FALSE,
    ...
  )
  invisible(x)
}

#' @rdname ba
#' @method summary ba_matrix
#' @export
summary.ba_matrix <- function(object,
                              digits = 3,
                              ci_digits = 3,
                              ...) {
  check_inherits(object, "ba_matrix")
  methods <- object$methods
  m <- length(methods)

  rows <- vector("list", m * (m - 1L) / 2L)
  k <- 0L
  for (i in 1:(m - 1L)) for (j in (i + 1L):m) {
    k <- k + 1L
    rows[[k]] <- list(
      method1 = methods[i],
      method2 = methods[j],
      bias = round(object$bias[i, j], digits),
      sd_loa = round(object$sd_loa[i, j], digits),
      loa_low = round(object$loa_lower[i, j], digits),
      loa_up = round(object$loa_upper[i, j], digits),
      width = round(object$width[i, j], digits),
      n_obs = suppressWarnings(as.integer(object$n[i, j])),
      loa_multiplier = round(as.numeric(object$loa_multiplier), digits),
      bias_lwr = round(object$mean_ci_low[i, j], ci_digits),
      bias_upr = round(object$mean_ci_high[i, j], ci_digits),
      lo_lwr = round(object$loa_lower_ci_low[i, j], ci_digits),
      lo_upr = round(object$loa_lower_ci_high[i, j], ci_digits),
      up_lwr = round(object$loa_upper_ci_low[i, j], ci_digits),
      up_upr = round(object$loa_upper_ci_high[i, j], ci_digits)
    )
  }

  out <- .mc_finalize_summary_df(
    do.call(rbind.data.frame, rows),
    class_name = "summary.ba_matrix"
  )
  attr(out, "conf.level") <- suppressWarnings(as.numeric(attr(object, "conf.level")))
  attr(out, "mode") <- as.integer(object$mode)
  out
}

#' @rdname ba
#' @method print summary.ba_matrix
#' @export
print.summary.ba_matrix <- function(x, digits = NULL, n = NULL,
                                    topn = NULL, max_vars = NULL,
                                    width = NULL, show_ci = NULL, ...) {
  cl <- suppressWarnings(as.numeric(attr(x, "conf.level")))
  mode <- as.integer(attr(x, "mode") %||% 1L)
  show_ci <- .mc_resolve_show_ci(show_ci, context = "summary")
  header <- .mc_header_with_ci(
    sprintf("Bland-Altman (%s)", if (identical(mode, 1L)) "pairwise row - column" else "pairwise column - row"),
    cl,
    show_ci
  )
  .mc_print_sectioned_table(
    x,
    sections = list(
      list(
        title = "Agreement estimates",
        cols = c("item1", "item2", "n_obs", "bias", "sd_loa", "loa_low", "loa_up", "width", "loa_multiplier")
      ),
      list(
        title = "Confidence intervals",
        cols = c("bias_lwr", "bias_upr", "lo_lwr", "lo_upr", "up_lwr", "up_upr")
      )
    ),
    header = header,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    print_overview = FALSE,
    ...
  )
  invisible(x)
}

#' @rdname ba
#' @method plot ba_matrix
#' @param pairs Optional character vector of pair labels to display.
#' @param against Optional single method name; if supplied, only contrasts
#'   involving that method are plotted.
#' @param facet_scales Either \code{"free_y"} (default) or \code{"fixed"}.
#' @export
plot.ba_matrix <- function(x,
                           pairs = NULL,
                           against = NULL,
                           facet_scales = c("free_y", "fixed"),
                           title = "Bland-Altman (pairwise)",
                           point_alpha = 0.6,
                           point_size = 1.8,
                           line_size = 0.7,
                           shade_ci = TRUE,
                           shade_alpha = 0.08,
                           smoother = c("none", "loess", "lm"),
                           show_value = TRUE,
                           ...) {
  check_inherits(x, "ba_matrix")
  check_bool(show_value, arg = "show_value")
  facet_scales <- match.arg(facet_scales)
  smoother <- match.arg(smoother)
  methods <- x$methods
  mode <- as.integer(x$mode %||% 1L)
  m <- length(methods)
  lab_pair <- function(i, j) {
    if (identical(mode, 2L)) {
      paste(methods[j], "\u2212", methods[i])
    } else {
      paste(methods[i], "\u2212", methods[j])
    }
  }

  idx_upper <- which(upper.tri(matrix(NA_real_, m, m)), arr.ind = TRUE)
  all_pairs <- data.frame(
    i = idx_upper[, 1],
    j = idx_upper[, 2],
    lab = lab_pair(idx_upper[, 1], idx_upper[, 2]),
    stringsAsFactors = FALSE
  )

  if (!is.null(against)) {
    if (!against %in% methods) {
      abort_bad_arg("against",
        message = "must be one of: {methods*?{, }{ and }}.",
        methods = methods
      )
    }
    hit <- match(against, methods)
    all_pairs <- subset(all_pairs, i == hit | j == hit)
  } else if (!is.null(pairs)) {
    all_pairs <- subset(all_pairs, lab %in% pairs)
    if (!nrow(all_pairs)) {
      abort_bad_arg("pairs",
        message = "None of the requested pairs matched available contrasts."
      )
    }
  }
  pairs_order <- all_pairs$lab

  bands <- data.frame(
    pair = lab_pair(idx_upper[, 1], idx_upper[, 2]),
    md = as.vector(x$bias[idx_upper]),
    loaL = as.vector(x$loa_lower[idx_upper]),
    loaU = as.vector(x$loa_upper[idx_upper]),
    md_l = as.vector(x$mean_ci_low[idx_upper]),
    md_u = as.vector(x$mean_ci_high[idx_upper]),
    lo_l = as.vector(x$loa_lower_ci_low[idx_upper]),
    lo_u = as.vector(x$loa_lower_ci_high[idx_upper]),
    hi_l = as.vector(x$loa_upper_ci_low[idx_upper]),
    hi_u = as.vector(x$loa_upper_ci_high[idx_upper]),
    stringsAsFactors = FALSE
  )
  bands <- bands[bands$pair %in% pairs_order, , drop = FALSE]
  bands$pair <- factor(bands$pair, levels = pairs_order)

  pts <- do.call(rbind, lapply(seq_len(nrow(all_pairs)), function(k) {
    i <- all_pairs$i[k]
    j <- all_pairs$j[k]
    a <- x$data_wide[[i]]
    b <- x$data_wide[[j]]
    keep <- !is.na(a) & !is.na(b)
    if (!any(keep)) {
      return(NULL)
    }
    means <- (a[keep] + b[keep]) / 2
    diffs <- if (identical(mode, 2L)) b[keep] - a[keep] else a[keep] - b[keep]
    data.frame(pair = all_pairs$lab[k], means = means, diffs = diffs, stringsAsFactors = FALSE)
  }))
  if (!is.null(pts) && nrow(pts)) {
    pts$pair <- factor(pts$pair, levels = pairs_order)
  } else {
    pts <- NULL
  }

  p <- ggplot2::ggplot() + ggplot2::facet_wrap(~ pair, scales = facet_scales)
  p <- p +
    ggplot2::geom_blank(data = bands, ggplot2::aes(x = 0, y = md)) +
    ggplot2::geom_blank(data = bands, ggplot2::aes(x = 0, y = loaL)) +
    ggplot2::geom_blank(data = bands, ggplot2::aes(x = 0, y = loaU))

  if (shade_ci) {
    p <- p +
      ggplot2::geom_rect(data = bands, ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = md_l, ymax = md_u),
                         inherit.aes = FALSE, alpha = shade_alpha) +
      ggplot2::geom_rect(data = bands, ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = lo_l, ymax = lo_u),
                         inherit.aes = FALSE, alpha = shade_alpha) +
      ggplot2::geom_rect(data = bands, ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = hi_l, ymax = hi_u),
                         inherit.aes = FALSE, alpha = shade_alpha)
  } else {
    p <- p +
      ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = md_l), linetype = "dashed") +
      ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = md_u), linetype = "dashed") +
      ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = lo_l), linetype = "dashed") +
      ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = lo_u), linetype = "dashed") +
      ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = hi_l), linetype = "dashed") +
      ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = hi_u), linetype = "dashed")
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dotted", colour = "grey40") +
    ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = md), linewidth = line_size) +
    ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = loaL), linewidth = line_size) +
    ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = loaU), linewidth = line_size)

  if (!is.null(pts)) {
    p <- p + ggplot2::geom_point(data = pts, ggplot2::aes(x = means, y = diffs),
                                 alpha = point_alpha, size = point_size)
    if (smoother == "lm") {
      p <- p + ggplot2::geom_smooth(data = pts, ggplot2::aes(x = means, y = diffs),
                                    method = "lm", se = FALSE, linewidth = 0.7)
    } else if (smoother == "loess") {
      p <- p + ggplot2::geom_smooth(data = pts, ggplot2::aes(x = means, y = diffs),
                                    method = "loess", se = FALSE, linewidth = 0.7, span = 0.9)
    }
  }

  p +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), ...) +
    ggplot2::labs(
      title = title,
      x = if (!is.null(pts)) "Mean of methods" else NULL,
      y = if (identical(mode, 2L)) "Difference (column - row)" else "Difference (row - column)"
    )
}
