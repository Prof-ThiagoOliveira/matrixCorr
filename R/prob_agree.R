#' Probability of Agreement for Reliability Curves
#'
#' @description
#' `prob_agree()` implements the probability of agreement following Stevens and
#' Anderson-Cook (2017). It fits binomial reliability curves by population or
#' generation and estimates, for each pair, the probability that the fitted
#' reliabilities differ by no more than a user-specified practically negligible
#' amount.
#'
#' This is not a correlation coefficient and should not be interpreted as
#' linear association. It is a tolerance-based agreement measure. It is
#' computed from the sampling distribution of the estimated difference between
#' two reliability curves.
#'
#' @details
#' For each pair of fitted reliability curves \eqn{\hat\pi_1(a)} and
#' \eqn{\hat\pi_2(a)}, Stevens and Anderson-Cook define
#' \deqn{\theta(a) = P\{|\tilde\pi_1(a) - \tilde\pi_2(a)| \le \delta(a) \mid a\}.}
#' The models are fit with a binomial GLM using either a logit or probit link,
#' with \eqn{g(\pi_i) = \beta_0 + \beta_1 a_i}. The C++ backend evaluates the
#' two-population large-sample normal approximation described in the paper and
#' Supplementary Material A; when more than two groups are supplied, `prob_agree()`
#' applies that calculation to all pairwise group comparisons.
#'
#' Missing rows in `response`, `predictor`, or `group` are removed before model
#' fitting. The response must be binary, coded as `0`/`1`, `FALSE`/`TRUE`, or a
#' two-level factor.
#'
#' @param data A data frame containing the response, predictor, and group
#'   variables.
#' @param response Character scalar naming the binary pass/fail response column.
#' @param predictor Character scalar naming the numeric age, time, or operating
#'   condition column.
#' @param group Character scalar naming the population/generation column. All
#'   pairwise comparisons among the observed non-missing levels are evaluated.
#' @param delta Positive scalar or vector tolerance for symmetric limits
#'   `-delta` and `delta`. A length-two vector with a negative lower bound and a
#'   positive upper bound, for example `c(-0.03, 0.05)`, is treated as an
#'   asymmetric tolerance interval. No default is provided because the tolerance
#'   is part of the scientific estimand.
#' @param limits Numeric vector `c(lower, upper)` for asymmetric tolerance
#'   limits. Mutually exclusive with `delta`.
#' @param link Link used for the reliability curve: `"logit"` or `"probit"`.
#' @param newdata Optional data frame containing predictor values where the
#'   probability of agreement is evaluated. If `NULL`, a regular grid over the
#'   observed predictor range is used.
#' @param grid_size Integer number of evaluation points when `newdata = NULL`.
#' @param ci Logical; if `TRUE`, attach pointwise confidence limits.
#' @param conf_level Confidence level for intervals.
#' @param max_iter Maximum number of IRLS iterations for each fitted curve.
#' @param tol Convergence tolerance for IRLS coefficient updates.
#' @param verbose Logical; if `TRUE`, emits a short `cli` message.
#'
#' @return A data frame with class
#' `c("prob_agree_curve", "prob_agree", "data.frame")` and columns `group1`,
#' `group2`, the predictor, `prob_agree`, and optional `lwr.ci`, `upr.ci`.
#'
#' @references
#' Stevens, N. T. and Anderson-Cook, C. M. (2017). Comparing the Reliability of
#' Related Populations With the Probability of Agreement. *Technometrics*,
#' 59(3), 371-380. doi:10.1080/00401706.2016.1214180.
#'
#' @examples
#' # Stevens and Anderson-Cook's probability of agreement evaluates whether
#' # two related populations have reliability curves that are similar enough
#' # to be treated as practically homogeneous. At each age, the agreement
#' # hypothesis is that the difference between the two fitted reliabilities
#' # lies inside the user-specified practical tolerance interval.
#' set.seed(1)
#' n <- 160
#' dat <- data.frame(
#'   age = c(runif(n, 0, 60), runif(n, 0, 45)),
#'   generation = rep(c("Gen1", "Gen2"), each = n)
#' )
#' eta <- ifelse(
#'   dat$generation == "Gen1",
#'   4.3 - 0.045 * dat$age,
#'   4.0 - 0.040 * dat$age
#' )
#' dat$pass <- rbinom(nrow(dat), size = 1, prob = plogis(eta))
#'
#' fit_pa <- prob_agree(
#'   dat,
#'   response = "pass",
#'   predictor = "age",
#'   group = "generation",
#'   delta = 0.05,
#'   link = "logit",
#'   ci = TRUE
#' )
#' print(fit_pa)
#' summary(fit_pa)
#' plot(fit_pa)
#'
#' # Four generations are compared as all pairwise two-generation contrasts.
#' set.seed(2)
#' n4 <- 120
#' dat4 <- data.frame(
#'   age = rep(runif(n4, 0, 55), times = 4),
#'   generation = rep(paste0("Gen", 1:4), each = n4)
#' )
#' shifts <- c(4.2, 4.0, 3.8, 3.6)
#' slopes <- c(-0.040, -0.042, -0.044, -0.046)
#' gen_id <- match(dat4$generation, paste0("Gen", 1:4))
#' eta4 <- shifts[gen_id] + slopes[gen_id] * dat4$age
#' dat4$pass <- rbinom(nrow(dat4), size = 1, prob = plogis(eta4))
#'
#' fit4 <- prob_agree(
#'   dat4,
#'   response = "pass",
#'   predictor = "age",
#'   group = "generation",
#'   limits = c(-0.03, 0.05),
#'   link = "logit",
#'   ci = FALSE
#' )
#' print(fit4)
#' plot(fit4)
#'
#' @export
prob_agree <- function(data,
                       response,
                       predictor,
                       group,
                       delta = NULL,
                       limits = NULL,
                       link = c("logit", "probit"),
                       newdata = NULL,
                       grid_size = 100L,
                       ci = TRUE,
                       conf_level = 0.95,
                       max_iter = 50L,
                       tol = 1e-8,
                       verbose = FALSE) {
  link <- match.arg(link)
  check_bool(ci, arg = "ci")
  check_bool(verbose, arg = "verbose")
  check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  grid_size <- check_scalar_int_pos(grid_size, arg = "grid_size")
  max_iter <- check_scalar_int_pos(max_iter, arg = "max_iter")
  check_scalar_numeric(tol, arg = "tol", lower = 0, finite = TRUE)

  if (!is.data.frame(data)) {
    abort_bad_arg("data", message = "must be a data frame.")
  }
  for (arg in c("response", "predictor", "group")) {
    val <- get(arg, envir = environment())
    if (!is.character(val) || length(val) != 1L || is.na(val) || !nzchar(val)) {
      abort_bad_arg(arg, message = "must be a single column name.")
    }
    if (!val %in% names(data)) {
      abort_bad_arg(arg, message = "must name a column in {.arg data}.")
    }
  }

  y <- .pa_response01(data[[response]])
  x <- data[[predictor]]
  if (!is.numeric(x)) {
    abort_bad_arg("predictor", message = "must name a numeric column.")
  }
  g <- factor(data[[group]])
  keep <- stats::complete.cases(y, x, g)
  if (sum(keep) < 6L) {
    abort_bad_arg("data", message = "must contain at least six complete rows across the groups.")
  }
  y <- y[keep]
  x <- as.numeric(x[keep])
  g <- droplevels(g[keep])
  if (nlevels(g) < 2L) {
    abort_bad_arg("group", message = "must retain at least two levels after removing missing rows.")
  }

  eval_x <- .pa_eval_predictor(data[[predictor]], predictor, newdata, grid_size)
  tol_obj <- .pa_tolerance(delta, limits, length(eval_x))
  link_code <- switch(link, logit = 1L, probit = 2L)
  levels_g <- levels(g)
  pair_mat <- utils::combn(levels_g, 2L)
  pair_outputs <- vector("list", ncol(pair_mat))
  pair_models <- vector("list", ncol(pair_mat))

  for (j in seq_len(ncol(pair_mat))) {
    pair <- pair_mat[, j]
    idx <- g %in% pair
    pair_g <- droplevels(g[idx])
    raw <- prob_agree_fit_cpp(
      response = y[idx],
      predictor = x[idx],
      group = as.integer(pair_g),
      eval_predictor = eval_x,
      link_code = link_code,
      lower = tol_obj$lower,
      upper = tol_obj$upper,
      ci = ci,
      ci_method = if (ci) 1L else 0L,
      conf_level = conf_level,
      max_iter = max_iter,
      tol = tol
    )
    pair_outputs[[j]] <- .pa_build_output(
      raw,
      ci = ci,
      predictor = predictor,
      group1 = pair[[1L]],
      group2 = pair[[2L]]
    )
    pair_models[[j]] <- list(
      group1 = pair[[1L]],
      group2 = pair[[2L]],
      coefficients = stats::setNames(
        list(
          stats::setNames(as.numeric(raw$beta1), c("(Intercept)", predictor)),
          stats::setNames(as.numeric(raw$beta2), c("(Intercept)", predictor))
        ),
        pair
      ),
      vcov = stats::setNames(list(raw$vcov1, raw$vcov2), pair),
      n_obs = stats::setNames(as.integer(c(raw$n1, raw$n2)), pair),
      converged = stats::setNames(as.logical(raw$converged), pair),
      iterations = stats::setNames(as.integer(raw$iterations), pair)
    )
  }
  out <- do.call(rbind, pair_outputs)
  row.names(out) <- NULL
  class(out) <- c("prob_agree_curve", "prob_agree", "data.frame")

  attr(out, "method") <- "probability_of_agreement"
  attr(out, "description") <- "Stevens-Anderson-Cook probability of agreement"
  attr(out, "package") <- "matrixCorr"
  attr(out, "estimator") <- "large_sample_normal"
  attr(out, "tolerance") <- tol_obj$tolerance
  attr(out, "link") <- link
  attr(out, "response") <- response
  attr(out, "predictor") <- predictor
  attr(out, "group") <- group
  attr(out, "group.levels") <- levels_g
  attr(out, "pairs") <- data.frame(
    group1 = pair_mat[1L, ],
    group2 = pair_mat[2L, ],
    pair = paste(pair_mat[1L, ], pair_mat[2L, ], sep = " vs "),
    stringsAsFactors = FALSE
  )
  attr(out, "pair_models") <- stats::setNames(pair_models, attr(out, "pairs")$pair)
  attr(out, "ci.method") <- if (ci) "delta_method" else "none"
  attr(out, "conf.level") <- conf_level
  attr(out, "reference") <- "Stevens and Anderson-Cook (2017), Technometrics, doi:10.1080/00401706.2016.1214180"

  inform_if_verbose(
    "Computed probability of agreement for {ncol(pair_mat)} pairwise comparison{?s} at {length(eval_x)} evaluation point{?s}.",
    .verbose = verbose
  )
  out
}

.pa_response01 <- function(x) {
  if (is.logical(x)) return(as.numeric(x))
  if (is.factor(x)) {
    if (nlevels(droplevels(x)) != 2L) abort_bad_arg("response", message = "must be binary.")
    return(as.numeric(droplevels(x)) - 1)
  }
  if (is.numeric(x) || is.integer(x)) {
    out <- as.numeric(x)
    vals <- stats::na.omit(unique(out))
    if (!all(vals %in% c(0, 1))) abort_bad_arg("response", message = "must contain only 0/1 values.")
    return(out)
  }
  abort_bad_arg("response", message = "must be binary numeric, logical, or a two-level factor.")
}

.pa_eval_predictor <- function(x, predictor, newdata, grid_size) {
  if (is.null(newdata)) {
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng)) || rng[[1L]] == rng[[2L]]) {
      abort_bad_arg("predictor", message = "must have a finite non-zero range when {.arg newdata} is NULL.")
    }
    return(seq(rng[[1L]], rng[[2L]], length.out = grid_size))
  }
  out <- if (is.data.frame(newdata)) {
    if (!predictor %in% names(newdata)) abort_bad_arg("newdata", message = "must contain the predictor column {.val {predictor}}.")
    newdata[[predictor]]
  } else {
    newdata
  }
  if (!is.numeric(out) || any(!is.finite(out))) abort_bad_arg("newdata", message = "must provide finite numeric predictor values.")
  as.numeric(out)
}

.pa_tolerance <- function(delta, limits, n) {
  if (!is.null(delta) && !is.null(limits)) cli::cli_abort("{.arg delta} and {.arg limits} are mutually exclusive.")
  if (is.null(delta) && is.null(limits)) cli::cli_abort("A practical tolerance must be supplied via {.arg delta} or {.arg limits}.")
  if (!is.null(delta)) {
    delta <- as.numeric(delta)
    delta_input <- delta
    if (length(delta) == 2L && any(delta < 0)) {
      if (any(!is.finite(delta)) || delta[[1L]] > delta[[2L]]) {
        abort_bad_arg("delta", message = "must be a finite c(lower, upper) vector with lower <= upper when used as asymmetric tolerance.")
      }
      return(list(
        lower = rep(delta[[1L]], n),
        upper = rep(delta[[2L]], n),
        tolerance = list(delta = NULL, limits = delta)
      ))
    }
    if (!(length(delta) %in% c(1L, n)) || any(!is.finite(delta)) || any(delta < 0)) {
      abort_bad_arg("delta", message = "must be a non-negative finite scalar or vector, or c(lower, upper) for asymmetric tolerance.")
    }
    if (length(delta) == 1L) delta <- rep(delta, n)
    return(list(
      lower = -delta,
      upper = delta,
      tolerance = list(delta = delta_input, limits = NULL)
    ))
  }
  limits <- as.numeric(limits)
  if (length(limits) != 2L || any(!is.finite(limits)) || limits[[1L]] > limits[[2L]]) {
    abort_bad_arg("limits", message = "must be a finite numeric vector c(lower, upper) with lower <= upper.")
  }
  list(
    lower = rep(limits[[1L]], n),
    upper = rep(limits[[2L]], n),
    tolerance = list(delta = NULL, limits = limits)
  )
}

.pa_build_output <- function(raw, ci, predictor, group1, group2) {
  out <- data.frame(
    group1 = group1,
    group2 = group2,
    predictor_value = as.numeric(raw$predictor),
    prob_agree = as.numeric(raw$prob_agree),
    check.names = FALSE
  )
  names(out)[[3L]] <- predictor
  if (isTRUE(ci)) {
    out$lwr.ci <- as.numeric(raw$lwr.ci)
    out$upr.ci <- as.numeric(raw$upr.ci)
  }
  structure(out, class = c("prob_agree_curve", "prob_agree", "data.frame"))
}

#' @method print prob_agree
#' @export
print.prob_agree <- function(x, digits = 4, ...) {
  rng <- range(x$prob_agree, na.rm = TRUE)
  tol <- attr(x, "tolerance", exact = TRUE)
  tol_txt <- if (!is.null(tol$delta)) {
    if (length(tol$delta) == 1L) sprintf("symmetric delta = %s", format(tol$delta, digits = digits)) else "symmetric vector delta"
  } else {
    sprintf("limits = [%s, %s]", format(tol$limits[[1L]], digits = digits), format(tol$limits[[2L]], digits = digits))
  }
  groups <- attr(x, "group.levels", exact = TRUE)
  pairs <- attr(x, "pairs", exact = TRUE)
  n_eval <- length(unique(x[[attr(x, "predictor", exact = TRUE)]]))
  cat("Probability-of-agreement result\n")
  cat("  method      : Stevens-Anderson-Cook large-sample normal\n")
  cat("  groups      :", paste(groups, collapse = ", "), "\n")
  cat("  comparisons :", nrow(pairs), "\n")
  cat("  link        :", attr(x, "link", exact = TRUE), "\n")
  cat("  tolerance   :", tol_txt, "\n")
  cat("  evaluations :", n_eval, "per comparison\n")
  cat("  range       :", paste(format(round(rng, digits), nsmall = digits), collapse = " to "), "\n")
  cat("  ci          :", if (all(c("lwr.ci", "upr.ci") %in% names(x))) "yes" else "no", "\n\n")
  print(.pa_preview_df(x, digits = digits), row.names = FALSE, ...)
  invisible(x)
}

#' @method summary prob_agree
#' @export
summary.prob_agree <- function(object, digits = 4, ...) {
  out <- as.data.frame(object)
  num_cols <- vapply(out, is.numeric, logical(1))
  out[num_cols] <- lapply(out[num_cols], round, digits = digits)
  structure(out, class = c("summary.prob_agree", "data.frame"))
}

#' @method print summary.prob_agree
#' @export
print.summary.prob_agree <- function(x, ...) {
  cat("Probability-of-agreement summary\n\n")
  print.data.frame(.pa_preview_df(x), row.names = FALSE, ...)
  n_omitted <- nrow(x) - nrow(.pa_preview_df(x))
  if (n_omitted > 0L) {
    cat("\n", n_omitted, " additional evaluation point", if (n_omitted == 1L) "" else "s", " not shown\n", sep = "")
  }
  invisible(x)
}

.pa_preview_df <- function(x, n = 10L, digits = 4) {
  out <- as.data.frame(x)
  if (nrow(out) > n) out <- utils::head(out, n)
  num_cols <- vapply(out, is.numeric, logical(1))
  out[num_cols] <- lapply(out[num_cols], round, digits = digits)
  out
}

#' Plot probability-of-agreement results
#'
#' @param x An object returned by `prob_agree()`.
#' @param threshold Optional probability threshold shown in curve plots.
#' @param title Optional plot title.
#' @param style Plot style: `"auto"` and `"curve"` draw pairwise probability
#'   curves; `"facet"` draws one curve panel per comparison; `"heatmap"` shows
#'   probability by predictor value and comparison.
#' @param show_ci Logical; for curve plots, controls whether confidence ribbons
#'   are shown. By default ribbons are shown for one comparison and suppressed
#'   for multiple comparisons.
#' @param ... Additional arguments passed to [ggplot2::theme()].
#'
#' @return A `ggplot` object.
#' @method plot prob_agree
#' @export
plot.prob_agree <- function(x,
                            threshold = NULL,
                            title = "Probability of agreement",
                            style = c("auto", "curve", "facet", "heatmap"),
                            show_ci = NULL,
                            ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) cli::cli_abort("Package {.pkg ggplot2} is required for plotting.")
  style <- match.arg(style)
  predictor <- attr(x, "predictor", exact = TRUE) %||% names(x)[[2L]]
  df <- as.data.frame(x)
  df <- df[order(df[[predictor]]), , drop = FALSE]
  df$pair <- paste(df$group1, df$group2, sep = " vs ")
  n_pairs <- length(unique(df$pair))
  if (identical(style, "auto")) style <- "curve"

  if (identical(style, "heatmap")) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[predictor]], y = .data$pair, fill = .data$prob_agree)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(
        limits = c(0, 1),
        option = "C",
        name = "Pr(agreement)"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), ...) +
      ggplot2::labs(title = title, x = predictor, y = "Comparison")
    if (!is.null(threshold)) {
      check_scalar_numeric(threshold, arg = "threshold", lower = 0, upper = 1)
    }
    return(p)
  }

  if (identical(style, "facet") && is.null(show_ci)) {
    show_ci <- all(c("lwr.ci", "upr.ci") %in% names(df))
  }
  if (is.null(show_ci)) show_ci <- n_pairs == 1L
  check_bool(show_ci, arg = "show_ci")
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[predictor]], y = .data$prob_agree, colour = .data$pair))
  has_ci_ribbon <- isTRUE(show_ci) && all(c("lwr.ci", "upr.ci") %in% names(df))
  if (has_ci_ribbon) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lwr.ci, ymax = .data$upr.ci, fill = .data$pair),
      alpha = 0.15,
      colour = NA
    )
  }
  p <- p +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), ...) +
    ggplot2::labs(title = title, x = predictor, y = "Probability of agreement", colour = "Comparison")
  if (identical(style, "facet")) {
    p <- p +
      ggplot2::facet_wrap(ggplot2::vars(.data$pair)) +
      ggplot2::guides(colour = "none", fill = "none")
  }
  if (has_ci_ribbon) {
    p <- p + ggplot2::labs(fill = "Comparison")
  }
  if (n_pairs == 1L) {
    p <- p + ggplot2::guides(colour = "none", fill = "none")
  }
  if (!is.null(threshold)) {
    check_scalar_numeric(threshold, arg = "threshold", lower = 0, upper = 1)
    p <- p + ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", colour = "grey35")
  }
  p
}
