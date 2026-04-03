#' @title Bland-Altman statistics with confidence intervals
#'
#' @description
#' Computes Bland-Altman mean difference and limits of agreement (LoA)
#' between two numeric measurement vectors, including t-based confidence
#' intervals for the mean difference and for each LoA using 'C++' backend.
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
#' @param group1,group2 Numeric vectors of equal length.
#' @param loa_multiplier Positive scalar; the multiple of the standard deviation used to
#'   define the LoA (default 1.96 for nominal 95\% agreement). The confidence
#'   intervals always use \eqn{t_{n-1,\,1-\alpha/2}} regardless of this choice.
#' @param mode Integer; 1 uses \code{group1 - group2}, 2 uses \code{group2 - group1}.
#' @param conf_level Confidence level for CIs (default 0.95).
#' @param verbose Logical; if TRUE, prints how many OpenMP threads are used.
#'
#' @return An object of class \code{"ba"} (list) with elements:
#' \itemize{
#'   \item \code{means}, \code{diffs}: numeric vectors
#'   \item \code{groups}: data.frame used after NA removal
#'   \item \code{based.on}: integer, number of pairs used
#'   \item \code{lower.limit}, \code{mean.diffs}, \code{upper.limit}
#'   \item \code{lines}: named numeric vector (lower, mean, upper)
#'   \item \code{CI.lines}: named numeric vector for CIs of those lines
#'   \item \code{loa_multiplier}, \code{critical.diff}
#' }
#'
#' @seealso \code{\link{print.ba}}, \code{\link{plot.ba}},
#'  \code{\link{ccc}},\code{\link{ccc_rm_ustat}},
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
               verbose = FALSE) {
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
  if (isTRUE(verbose)) cat("Using", ba_openmp_threads(), "OpenMP threads\n")

  # -- compute in C++ ---------------------------------------------------------
  ba_out <- bland_altman_cpp(group1, group2, loa_multiplier, mode, conf_level)

  ba_out <- structure(ba_out, class = "ba")
  attr(ba_out, "method")      <- "Bland-Altman"
  attr(ba_out, "description") <- "Mean difference and limits of agreement with CIs"
  attr(ba_out, "package")     <- "matrixCorr"
  attr(ba_out, "conf.level")  <- conf_level
  attr(ba_out, "called.with") <- length(group1)
  ba_out
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
    n              = as.integer(object$based.on),
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

  out <- structure(out, class = c("summary.ba", "data.frame"))
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
        cols = c("n", "bias", "sd_loa", "loa_low", "loa_up", "width", "loa_multiplier")
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
