#' @title Bland-Altman statistics with confidence intervals
#'
#' @description
#' Computes Bland-Altman mean difference and limits of agreement (LoA)
#' between two numeric measurement vectors, including t-based confidence
#' intervals for the mean difference and for each LoA using 'C++' backend.
#' Two usage modes:
#' \itemize{
#'   \item \strong{Two-vector mode (default):} If \code{group1} and \code{group2}
#'         are numeric vectors, computes standard Bland-Altman statistics with CIs
#'         and returns a \code{"ba"} object.
#'   \item \strong{Multi-method mode:} If \code{group2 = NULL} and \code{group1}
#'         is a numeric \code{data.frame}/\code{data.table}/matrix with one column
#'         per method, computes the full pairwise Bland-Altman matrix across all
#'         columns and returns a \code{"ba_matrix"} object.
#' }
#'
#' Note: Lin's concordance correlation coefficient (CCC) is a complementary,
#' single-number summary of agreement (precision + accuracy) for continuous
#' variables. It is useful for
#' quick screening or reporting an overall CI, but may miss systematic or
#' magnitude-dependent bias; consider reporting CCC alongside Bland-Altman.
#'
#' @details
#' Given paired measurements \eqn{(x_i, y_i)}, Bland-Altman analysis uses
#' \eqn{d_i = x_i - y_i} (or \eqn{y_i - x_i} if \code{mode = 2}) and
#' \eqn{m_i = (x_i + y_i)/2}. The mean difference \eqn{\bar d} estimates bias.
#' The limits of agreement (LoA) are \eqn{\bar d \pm z \cdot s_d}, where
#' \eqn{s_d} is the sample standard deviation of \eqn{d_i} and \eqn{z}
#' (argument \code{two}) is typically 1.96 for nominal 95% LoA.
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
#' removed pairwise (any row with an \code{NA} in either input is dropped).
#'
#' @param group1 Numeric vector; or a numeric data.frame/data.table/matrix
#'   (columns = methods) when \code{group2 = NULL}.
#' @param group2 Optional numeric vector. If provided, the function performs a
#'   two-vector Bland-Altman comparison against \code{group1}.
#' @param two Positive scalar; the multiple of SD for the LoA (default 1.96).
#' @param mode Integer.
#'   When both \code{group1} and \code{group2} are numeric vectors,
#'   \code{mode = 1} computes \code{group1 - group2} and
#'   \code{mode = 2} computes \code{group2 - group1}.
#'   When \code{group2 = NULL} and \code{group1} is multi-column
#'   (matrix/data.frame; multi-method mode), \code{mode} is ignored and
#'   pairwise results are reported in a fixed \emph{row - column}
#'   orientation in the returned matrices.

#' @param conf_level Confidence level for CIs (default 0.95).
#' @param verbose Logical; if TRUE, prints how many OpenMP threads are used.
#'
#' @return Either a \code{"ba"} object (two-vector mode) or a
#' \code{"ba_matrix"} object containing:
#' \itemize{
#'   \item \code{means}, \code{diffs}: numeric vectors
#'   \item \code{groups}: data.frame used after NA removal
#'   \item \code{based.on}: integer, number of pairs used
#'   \item \code{lower.limit}, \code{mean.diffs}, \code{upper.limit}
#'   \item \code{lines}: named numeric vector (lower, mean, upper)
#'   \item \code{CI.lines}: named numeric vector for CIs of those lines
#'   \item \code{two}, \code{critical.diff}
#' }
#'
#' @seealso \code{\link{print.ba}}, \code{\link{plot.ba}},
#'  \code{\link{ccc}},\code{\link{ccc_pairwise_u_stat}},
#'  \code{\link{ccc_lmm_reml}}
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100, 100, 10)
#' y <- x + rnorm(100, 0, 8)
#' ba <- bland_altman(x, y)
#' print(ba)
#' plot(ba)
#'
#' # ---- Multi-method example (matrix/data.frame input) -----------------------
#' set.seed(42)
#' n <- 120
#' true <- rnorm(n, 50, 7)
#' A <- true + rnorm(n, 0, 3)           # baseline method
#' B <- true + 2 + rnorm(n, 0, 3)       # constant +2 bias
#' C <- 1.10 * true + rnorm(n, 0, 4)    # proportional bias, larger SD
#' D <- true + rnorm(n, 0, 6)           # higher random error
#' X <- data.frame(MethodA = A, MethodB = B, MethodC = C, MethodD = D)
#'
#' # group2 = NULL and multi-column group1 -> returns a "ba_matrix" object
#' baM <- bland_altman(X)
#'
#' # Summaries (orientation: row - column)
#' print(baM, show = c("bias", "width"))
#' round(baM$bias, 2)    # mean differences (row - column)
#' round(baM$width, 2)   # LoA widths (upper - lower)
#'
#' # Plot - faceted BA scatter by pair
#' plot(baM, data = X, smoother = "lm")
#'
#' # Focus on a reference method versus all others
#' plot(baM, data = X, against = "MethodA", facet_scales = "free_y")
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
bland_altman <- function(group1,
                         group2 = NULL,
                         two = 1.96,
                         mode = 1L,
                         conf_level = 0.95,
                         verbose = FALSE) {

  # --- Multi-method mode dispatch -------------------------------------------
  # group2 is NULL and group1 is a multi-column object -> compute BA matrix
  if (is.null(group2)) {
    if (is.data.frame(group1) || is.matrix(group1)) {
      return(
        bland_altman_matrix(
          group1, two = two, conf_level = conf_level, verbose = verbose
        )
      )
    } else {
      stop(paste("When 'group2' is NULL, 'group1' must be a",
                 "data.frame/data.table or matrix with at least",
                 "two numeric columns."))
    }
  }

  # --- Two-vector mode (original behaviour) ---------------------------------
  # validate
  if (!is.numeric(group1) || !is.numeric(group2))
    stop(paste("In two-vector mode, 'group1' and 'group2' must be numeric",
               "vectors. For multi-method use, supply 'group1' as a",
               "multi-column object and set group2 = NULL."))
  if (length(group1) != length(group2))
    stop("'group1' and 'group2' must have the same length.")
  if (!is.numeric(two) || length(two) != 1L || two <= 0)
    stop("'two' must be a positive scalar.")
  mode <- as.integer(mode)
  if (!mode %in% c(1L, 2L))
    stop("'mode' must be 1 or 2.")
  if (!(is.numeric(conf_level) && length(conf_level) == 1L && conf_level > 0 &&
        conf_level < 1))
    stop("'conf_level' must be in (0, 1).")

  if (isTRUE(verbose)) cat("Using", ba_openmp_threads(), "OpenMP threads\n")

  # compute (C++; pairwise NA removal is handled inside bland_altman_cpp)
  ba_out <- bland_altman_cpp(group1, group2, two, mode, conf_level)

  attr(ba_out, "method")      <- "Bland-Altman"
  attr(ba_out, "description") <- "Mean difference and limits of agreement with CIs"
  attr(ba_out, "package")     <- "matrixCorr"
  attr(ba_out, "conf.level")  <- conf_level
  attr(ba_out, "called.with") <- length(group1)

  class(ba_out) <- c("ba", "bland_altman")
  ba_out
}

#-------------------------------------------------------------------------------
# Pairwise Bland-Altman matrix across multiple methods
#-------------------------------------------------------------------------------
#' @title Pairwise Bland-Altman matrix across multiple methods
#' @description Internal helper used by \code{bland_altman()} when \code{group2 = NULL}
#' and \code{group1} is a multi-column object. Computes all pairwise Bland-Altman
#' results using the two-vector \code{bland_altman()} path.
#' @param X Numeric data.frame/data.table or matrix; columns are methods.
#' @param two Positive scalar; LoA multiple of SD (default 1.96).
#' @param conf_level Confidence level for CIs (default 0.95).
#' @param verbose Logical; if TRUE, prints a single line about threads.
#' @return An object of class \code{"ba_matrix"}: list of matrices with row/col
#'   names = method names. Upper triangle corresponds to (row - column).
#' @keywords internal
bland_altman_matrix <- function(X,
                                two = 1.96,
                                conf_level = 0.95,
                                verbose = FALSE) {
  X <- as.data.frame(X, check.names = FALSE)
  if (!all(vapply(X, is.numeric, logical(1))))
    stop("All columns of 'group1' must be numeric for multi-method Bland-Altman.")
  m <- ncol(X)
  if (m < 2L) stop("Need at least two methods (two columns) when group2 = NULL.")
  methods <- colnames(X); if (is.null(methods)) methods <- paste0("M", seq_len(m))
  if (isTRUE(verbose)) cat("Multi-method mode: using", ba_openmp_threads(),
                           "OpenMP threads\n")

  mk <- function(val = NA_real_) {
    out <- matrix(val, m, m)
    dimnames(out) <- list(methods, methods)
    out
    }
  bias         <- mk()
  sd_diff      <- mk()
  loa_lower    <- mk()
  loa_upper    <- mk()
  width        <- mk()
  n_mat        <- mk(NA_integer_)
  mean_ci_low  <- mk()
  mean_ci_high <- mk()
  lo_ci_low    <- mk()
  lo_ci_high   <- mk()
  hi_ci_low    <- mk()
  hi_ci_high   <- mk()

  for (j in 1:(m-1)) for (k in (j+1):m) {
    xj <- X[[j]]; xk <- X[[k]]
    ok <- is.finite(xj) & is.finite(xk)     # pairwise NA removal to mirror C++
    xj <- xj[ok]; xk <- xk[ok]
    if (length(xj) < 3L) next

    # Use the two-vector path (C++ backend)
    ba <- bland_altman(xj, xk, two = two, mode = 1L, conf_level = conf_level,
                       verbose = FALSE)

    md   <- as.numeric(ba$mean.diffs)
    loaL <- as.numeric(ba$lower.limit)
    loaU <- as.numeric(ba$upper.limit)
    n    <- as.integer(ba$based.on)
    sdd  <- as.numeric(ba$critical.diff) / two

    bias[j,k]      <- md;     bias[k,j]      <- -md
    sd_diff[j,k]   <- sdd;    sd_diff[k,j]   <- sdd
    loa_lower[j,k] <- loaL;   loa_lower[k,j] <- -loaU
    loa_upper[j,k] <- loaU;   loa_upper[k,j] <- -loaL
    width[j,k]     <- loaU - loaL; width[k,j] <- width[j,k]
    n_mat[j,k]     <- n;      n_mat[k,j]     <- n

    cil <- function(nm) as.numeric(ba$CI.lines[[nm]])
    mean_ci_low[j,k]  <- cil("mean.diff.ci.lower")
    mean_ci_low[k,j]  <- -cil("mean.diff.ci.upper")
    mean_ci_high[j,k] <- cil("mean.diff.ci.upper")
    mean_ci_high[k,j] <- -cil("mean.diff.ci.lower")
    lo_ci_low[j,k]    <- cil("lower.limit.ci.lower")
    lo_ci_low[k,j]    <- -cil("upper.limit.ci.upper")
    lo_ci_high[j,k]   <- cil("lower.limit.ci.upper")
    lo_ci_high[k,j]   <- -cil("upper.limit.ci.lower")
    hi_ci_low[j,k]    <- cil("upper.limit.ci.lower")
    hi_ci_low[k,j]    <- -cil("lower.limit.ci.upper")
    hi_ci_high[j,k]   <- cil("upper.limit.ci.upper")
    hi_ci_high[k,j]   <- -cil("lower.limit.ci.lower")
  }

  out <- list(
    bias = bias,
    sd_diff = sd_diff,
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
    two = two,
    conf_level = conf_level,
    methods = methods,
    na_strategy = "pairwise"
  )
  class(out) <- c("ba_matrix","list")
  out
}

#' @rdname bland_altman
#' @method print ba
#' @param x A \code{"ba"} object.
#' @param digits Number of digits for estimates (default 3).
#' @param ci_digits Number of digits for CI bounds (default 3).
#' @param ... Unused.
#' @export
print.ba <- function(x, digits = 3, ci_digits = 3, ...) {
  if (!inherits(x, "ba")) stop("Object is not of class 'ba'.")

  n   <- as.integer(x$based.on)
  two <- as.numeric(x$two)
  cl  <- suppressWarnings(as.numeric(attr(x, "conf.level")))
  if (!is.finite(cl)) cl <- NA_real_

  # core numbers (as scalars)
  bias   <- as.numeric(x$mean.diffs)
  loa_lo <- as.numeric(x$lower.limit)
  loa_hi <- as.numeric(x$upper.limit)
  sd_d   <- as.numeric(x$critical.diff) / two
  width  <- loa_hi - loa_lo

  # CIs (robust extraction by name)
  cil <- function(nm) as.numeric(x$CI.lines[[nm]])
  bias_l <- cil("mean.diff.ci.lower"); bias_u <- cil("mean.diff.ci.upper")
  lo_l   <- cil("lower.limit.ci.lower"); lo_u <- cil("lower.limit.ci.upper")
  hi_l   <- cil("upper.limit.ci.lower"); hi_u <- cil("upper.limit.ci.upper")

  # header
  if (is.finite(cl)) {
    cat(sprintf("Bland-Altman (n = %d) - LoA = mean +/- %.3g * SD, %g%% CI\n\n",
                n, two, 100*cl))
  } else {
    cat(sprintf("Bland-Altman (n = %d) - LoA = mean +/- %.3g * SD\n\n", n, two))
  }

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
  df$lwr      <- formatC(df$lwr,      format = "f", digits = ci_digits)
  df$upr      <- formatC(df$upr,      format = "f", digits = ci_digits)

  print(df, row.names = FALSE, right = FALSE)
  cat(sprintf("\nSD(differences): %s   LoA width: %s\n",
              formatC(sd_d, format = "f", digits = digits),
              formatC(width, format = "f", digits = digits)))
  invisible(x)
}

#' @rdname bland_altman
#' @method print ba_matrix
#' @param x A \code{"ba_matrix"} object.
#' @param digits Number of digits for estimates (default 3).
#' @param show Character vector indicating which matrices to print for a
#'   \code{"ba_matrix"} object. Allowed values are \code{"bias"},
#'   \code{"width"}, and \code{"sd_diff"}. \code{"bias"} prints the matrix of
#'   mean differences (row - column); \code{"width"} prints the matrix of limits
#'   of agreement widths (upper - lower); \code{"sd_diff"} prints the matrix of
#'   standard deviations of the pairwise differences. May include one or more
#'   of these values in any order. Defaults to \code{c("bias","width","sd_diff")}.
#' @param ... Unused.
#' @export
print.ba_matrix <- function(x, digits = 3, show = c("bias","width","sd_diff"), ...) {
  stopifnot(inherits(x, "ba_matrix"))
  show <- match.arg(show, choices = c("bias","width","sd_diff"), several.ok = TRUE)
  cat(sprintf("Pairwise Bland-Altman matrices (two = %.3g, conf.level = %.2f)\n",
              x$two, x$conf_level))
  if ("bias" %in% show) {
    cat("\nMean difference (row - column):\n")
    print(round(x$bias, digits))
  }
  if ("width" %in% show) {
    cat("\nLoA width (upper - lower):\n")
    print(round(x$width, digits))
  }
  if ("sd_diff" %in% show) {
    cat("\nSD of differences:\n")
    print(round(x$sd_diff, digits))
  }
  invisible(x)
}

#' @rdname bland_altman
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
                    ...) {
  if (!inherits(x, "ba")) stop("x must be of class 'ba'.")
  smoother <- match.arg(smoother)

  means <- as.numeric(x$means)
  diffs <- as.numeric(x$diffs)

  # scalars
  md    <- as.numeric(x$mean.diffs)
  loaL  <- as.numeric(x$lower.limit)
  loaU  <- as.numeric(x$upper.limit)
  two   <- as.numeric(x$two)
  n     <- as.integer(x$based.on)
  sd_d  <- as.numeric(x$critical.diff) / two
  cl    <- suppressWarnings(as.numeric(attr(x, "conf.level")))
  ci    <- function(nm) as.numeric(x$CI.lines[[nm]])

  if (is.null(subtitle)) {
    subtitle <- if (is.finite(cl)) {
      sprintf("n = %d  *  mean diff = %.2f  *  LoA = [%.2f, %.2f]  *  %g%% CI shown",
              n, md, loaL, loaU, 100*cl)
    } else {
      sprintf("n = %d  *  mean diff = %.2f  *  LoA = [%.2f, %.2f]",
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
           ybottom = ci("lower.limit.ci.lower"),
           ytop = ci("lower.limit.ci.upper"),
           border = NA, col = grDevices::gray(0.2, alpha = shade_alpha))
      rect(xleft = usr[1], xright = usr[2],
           ybottom = ci("upper.limit.ci.lower"),
           ytop = ci("upper.limit.ci.upper"),
           border = NA, col = grDevices::gray(0.2, alpha = shade_alpha))
    } else {
      abline(h = c(ci("mean.diff.ci.lower"), ci("mean.diff.ci.upper")), lty = 2)
      abline(h = c(ci("lower.limit.ci.lower"), ci("lower.limit.ci.upper")),
             lty = 2)
      abline(h = c(ci("upper.limit.ci.lower"), ci("upper.limit.ci.upper")),
             lty = 2)
    }
    # reference & main lines
    abline(h = 0,   col = "grey70", lty = 3)
    abline(h = md,  lwd = line_size * 1.4)
    abline(h = loaL, lwd = line_size * 1.2)
    abline(h = loaU, lwd = line_size * 1.2)

    # optional smoother
    if (smoother != "none") {
      fit <- if (smoother == "lm") {
        stats::lm(diffs ~ means)
      } else stats::loess(diffs ~ means)
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
                        ymin = ci("mean.diff.ci.lower"),
                        ymax = ci("mean.diff.ci.upper"),
                        alpha = shade_alpha) +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                        ymin = ci("lower.limit.ci.lower"),
                        ymax = ci("lower.limit.ci.upper"),
                        alpha = shade_alpha) +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                        ymin = ci("upper.limit.ci.lower"),
                        ymax = ci("upper.limit.ci.upper"),
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
    ggplot2::geom_hline(yintercept = 0, size = 0.4, linetype = "dotted",
                        color = "grey40") +
    ggplot2::geom_hline(yintercept = md,   size = line_size) +
    ggplot2::geom_hline(yintercept = loaL, size = line_size) +
    ggplot2::geom_hline(yintercept = loaU, size = line_size)

  if (smoother == "lm") {
    p <- p + ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.7)
  } else if (smoother == "loess") {
    p <- p + ggplot2::geom_smooth(method = "loess", se = FALSE,
                                  linewidth = 0.7, span = 0.9)
  }

  p <- p +
    ggplot2::coord_cartesian(ylim = y_rng, expand = TRUE) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), ...) +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = "Mean of methods", y = "Difference between methods")

  p
}

#' @rdname bland_altman
#' @method plot ba_matrix
#' @param x A \code{"ba_matrix"} object returned by \code{bland_altman()}
#'   when called with \code{group2 = NULL}.
#' @param data Optional wide numeric \code{data.frame}/matrix containing the
#'   original method columns used to build \code{x}. Required if \code{x} does
#'   not carry \code{attr(x, "data")}. The column names must match
#'   \code{x$methods}. This is used to compute point-level means and
#'   differences for each pairwise panel; the summaries stored in \code{x}
#'   (bias/LoA) are insufficient to draw the scatter without the raw data.
#' @param pairs Optional character vector of pair labels to display. Labels
#'   must follow the internal convention \emph{row - column}, e.g.
#'   \code{"MethodA - MethodB"}. If \code{NULL} (default), all
#'   upper-triangular pairs are shown.
#' @param against Optional single method name. If supplied, that method is
#'   plotted against every other method (one facet per pair). Takes precedence
#'   over \code{pairs} when both are provided.
#' @param facet_scales Facet scaling behaviour. Use \code{"free_y"} (default)
#'   to allow each panel its own y-axis (useful when LoA differ markedly across
#'   pairs), or \code{"fixed"} for a common y-axis across panels.
#' @param title Plot title.
#' @param ... Additional theme tweaks passed to \code{ggplot2::theme()}
#'   (e.g., \code{panel.grid = ggplot2::element_blank()}). Not forwarded to
#'   geoms.
#' @importFrom stats complete.cases
#' @export
plot.ba_matrix <- function(x,
                           data = NULL,
                           pairs = NULL,
                           against = NULL,
                           facet_scales = c("free_y","fixed"),
                           title = "Bland-Altman (pairwise, faceted)",
                           point_alpha = 0.6,
                           point_size  = 1.8,
                           line_size   = 0.7,
                           shade_ci    = TRUE,
                           shade_alpha = 0.08,
                           smoother    = c("none","loess","lm"),
                           ...) {
  if (!inherits(x, "ba_matrix")) stop("x must be of class 'ba_matrix'.")
  facet_scales <- match.arg(facet_scales)
  smoother <- match.arg(smoother)

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for faceted plotting. Please install it.")

  # Obtain wide raw data
  if (is.null(data)) data <- attr(x, "data")
  if (is.null(data))
    stop("Faceted BA plot needs the wide raw data. Supply it via `data =` or store it in attr(x, 'data').")
  data <- as.data.frame(data, check.names = FALSE)

  methods <- x$methods
  if (!all(methods %in% colnames(data)))
    stop("`data` does not contain all methods recorded in `x$methods`.")

  # Build the list of pairs to plot (row - column orientation)
  idx_upper <- which(upper.tri(matrix(NA_real_, length(methods), length(methods))), arr.ind = TRUE)
  all_pairs <- data.frame(
    j = idx_upper[,1],
    k = idx_upper[,2],
    lab = paste(methods[idx_upper[,1]], "\u2212", methods[idx_upper[,2]]), # "A - B"
    stringsAsFactors = FALSE
  )

  if (!is.null(against)) {
    if (!against %in% methods) stop("`against` must be one of: ", paste(methods, collapse = ", "))
    j_idx <- match(against, methods)
    keep  <- all_pairs$j == j_idx | all_pairs$k == j_idx
    all_pairs <- all_pairs[keep, , drop = FALSE]
    # ensure label orientation uses row - column;
    # if against is column, we'll still show "row - column"
  } else if (!is.null(pairs)) {
    all_pairs <- subset(all_pairs, lab %in% pairs)
    if (!nrow(all_pairs)) stop("None of the requested `pairs` matched the available pairs.")
  }

  # Construct per-panel data
  build_panel <- function(j, k, lab) {
    g1 <- data[[methods[j]]]; g2 <- data[[methods[k]]]
    ok <- is.finite(g1) & is.finite(g2)
    g1 <- g1[ok]; g2 <- g2[ok]
    if (length(g1) < 2L) return(NULL)

    means <- (g1 + g2) / 2
    diffs <- (g1 - g2)              # row - column orientation

    # panel stats from x (no recomputation)
    md   <- x$bias[j, k]
    loaL <- x$loa_lower[j, k]
    loaU <- x$loa_upper[j, k]

    # CIs (if available)
    md_l <- x$mean_ci_low[j, k];  md_u <- x$mean_ci_high[j, k]
    lo_l <- x$loa_lower_ci_low[j, k]; lo_u <- x$loa_lower_ci_high[j, k]
    hi_l <- x$loa_upper_ci_low[j, k]; hi_u <- x$loa_upper_ci_high[j, k]

    pts <- data.frame(
      means = means, diffs = diffs, pair = lab, stringsAsFactors = FALSE
    )

    bands <- data.frame(
      pair = lab,
      md = md, loaL = loaL, loaU = loaU,
      md_l = md_l, md_u = md_u,
      lo_l = lo_l, lo_u = lo_u,
      hi_l = hi_l, hi_u = hi_u,
      stringsAsFactors = FALSE
    )
    list(points = pts, bands = bands)
  }

  res <- lapply(seq_len(nrow(all_pairs)), function(i) {
    build_panel(all_pairs$j[i], all_pairs$k[i], all_pairs$lab[i])
  })
  res <- Filter(Negate(is.null), res)
  if (!length(res)) stop("No panel had at least two complete pairs after NA removal.")

  pts   <- do.call(rbind, lapply(res, `[[`, "points"))
  bands <- do.call(rbind, lapply(res, `[[`, "bands"))

  p <- ggplot2::ggplot(pts, ggplot2::aes(x = means, y = diffs)) +
    ggplot2::geom_point(alpha = point_alpha, size = point_size) +
    ggplot2::facet_wrap(~ pair, scales = facet_scales) +
    ggplot2::geom_hline(yintercept = 0, size = 0.35, linetype = "dotted",
                        colour = "grey40")

  # CI shading per panel  ----  add inherit.aes = FALSE
  if (shade_ci) {
    p <- p +
      ggplot2::geom_rect(
        data = bands,
        ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = md_l, ymax = md_u),
        inherit.aes = FALSE, alpha = shade_alpha
      ) +
      ggplot2::geom_rect(
        data = bands,
        ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = lo_l, ymax = lo_u),
        inherit.aes = FALSE, alpha = shade_alpha
      ) +
      ggplot2::geom_rect(
        data = bands,
        ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = hi_l, ymax = hi_u),
        inherit.aes = FALSE, alpha = shade_alpha
      )
  } else {
    p <- p +
      ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = md_l),
                          linetype = "dashed") +
      ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = md_u),
                          linetype = "dashed") +
      ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = lo_l),
                          linetype = "dashed") +
      ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = lo_u),
                          linetype = "dashed") +
      ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = hi_l),
                          linetype = "dashed") +
      ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = hi_u),
                          linetype = "dashed")
  }

  # --- Mean/LoA lines per panel (remove inherit.aes, use size =) ---
  p <- p +
    ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = md),
                        size = line_size) +
    ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = loaL),
                        size = line_size) +
    ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = loaU),
                        size = line_size)

  # Optional smoother per panel  ----  use size= for compatibility
  if (smoother == "lm") {
    p <- p + ggplot2::geom_smooth(method = "lm", se = FALSE, size = 0.7)
  } else if (smoother == "loess") {
    p <- p + ggplot2::geom_smooth(method = "loess", se = FALSE, size = 0.7,
                                  span = 0.9)
  }

  p +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), ...) +
    ggplot2::labs(
      title = title,
      x = "Mean of methods", y = "Difference between methods (row - column)"
    )
}
