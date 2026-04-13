#' @title Repeated-Measures Correlation (rmcorr)
#'
#' @description
#' Computes repeated-measures correlation for two or more continuous responses
#' observed repeatedly within subjects. Supply a \code{data.frame} plus column
#' names, or pass the response matrix/data frame and subject vector directly.
#' The repeated observations are indexed only by \code{subject}, thus no explicit
#' time variable is modeled, and the method targets a common within-subject
#' linear association after removing subject-specific means.
#'
#' @param data Optional \code{data.frame}/matrix containing the repeated-measures
#'   dataset.
#' @param response Either:
#'   \itemize{
#'     \item a character vector of at least two column names in \code{data}, or
#'     \item a numeric matrix/data frame with at least two columns.
#'   }
#'   If exactly two responses are supplied, the function returns a pairwise
#'   repeated-measures correlation object of class \code{"rmcorr"}. If three or
#'   more responses are supplied, a symmetric matrix of class
#'   \code{"rmcorr_matrix"} is returned.
#' @param subject Subject identifier (factor/character/integer/numeric) or a
#'   single character string naming the subject column in \code{data}.
#' @param conf_level Confidence level used for Wald confidence intervals on the
#'   repeated-measures correlation (default \code{0.95}).
#' @param na_method Character scalar controlling missing-data handling.
#'   \code{"error"} rejects missing values in the selected response columns or
#'   \code{subject}. \code{"pairwise"} uses pairwise complete cases and drops
#'   subjects with fewer than two complete pairs for the specific contrast.
#' @param n_threads Integer \eqn{\ge 1}. Number of OpenMP threads used by the
#'   C++ backend in matrix mode. Defaults to
#'   \code{getOption("matrixCorr.threads", 1L)}.
#' @param keep_data Logical (default \code{FALSE}). If \code{TRUE}, returned
#'   objects retain a compact copy of the numeric response matrix plus encoded
#'   subject identifiers. For pairwise fits this enables \code{plot()} to
#'   reconstruct subject-level fitted lines lazily; for matrix outputs it also
#'   allows \code{view_rmcorr_shiny()} to rebuild pairwise scatterplots from the
#'   returned object alone.
#' @param verbose Logical; if \code{TRUE}, prints a short note about the thread
#'   count used in matrix mode.
#' @param ... Deprecated compatibility aliases. The legacy \code{check_na}
#'   argument is still accepted here temporarily.
#'
#' @return
#' Either a \code{"rmcorr"} object (exactly two responses) or a
#' \code{"rmcorr_matrix"} object (pairwise results when \eqn{\ge}3 responses).
#'
#' \strong{If \code{"rmcorr"} (exactly two responses)}, outputs include:
#' \itemize{
#'   \item \code{estimate}; repeated-measures correlation estimate.
#'   \item \code{p_value}; two-sided p-value for the common within-subject slope.
#'   \item \code{lwr}, \code{upr}; confidence interval limits for
#'         \code{estimate}.
#'   \item \code{slope}; common within-subject slope.
#'   \item \code{df}; residual degrees of freedom \eqn{N - S - 1}.
#'   \item \code{n_obs}; number of complete observations retained after
#'         dropping subjects with fewer than two repeated pairs.
#'   \item \code{n_subjects}; number of contributing subjects.
#'   \item \code{responses}; names of the fitted response variables.
#'   \item compatibility aliases \code{r}, \code{conf_int}, and \code{based.on}
#'         are reconstructed on access without duplicate storage.
#'   \item when \code{keep_data = TRUE}, compact source data are retained so
#'         \code{plot()} can lazily reconstruct \code{data_long},
#'         \code{intercepts}, and fitted lines; these are otherwise not stored.
#' }
#'
#' \strong{If \code{"rmcorr_matrix"} (\eqn{\ge}3 responses)}, outputs are:
#' \itemize{
#'   \item a symmetric numeric matrix of pairwise repeated-measures correlations.
#'   \item attributes \code{method}, \code{description}, and
#'         \code{package = "matrixCorr"}.
#'   \item \code{diagnostics}; a list with square matrices for \code{slope},
#'         \code{p_value}, \code{df}, \code{n_complete}, \code{n_subjects},
#'         \code{conf_low}, and \code{conf_high}, plus scalar
#'         \code{conf_level}.
#' }
#'
#' @details
#' Repeated-measures correlation estimates the common \emph{within-subject}
#' linear association between two variables measured repeatedly on the same
#' subjects. It differs from agreement methods such as Lin's CCC or
#' Bland-Altman analysis because those target concordance or interchangeability,
#' whereas repeated-measures correlation targets the strength of the
#' subject-centred association.
#'
#' For subject \eqn{i = 1,\ldots,S} and repeated observations
#' \eqn{j = 1,\ldots,n_i}, let \eqn{x_{ij}} and \eqn{y_{ij}} denote the two
#' responses. Define subject-specific means
#' \deqn{\bar x_i = \frac{1}{n_i}\sum_{j=1}^{n_i} x_{ij}, \qquad
#'       \bar y_i = \frac{1}{n_i}\sum_{j=1}^{n_i} y_{ij}.}
#' The repeated-measures correlation uses within-subject centred values
#' \deqn{x^\ast_{ij} = x_{ij} - \bar x_i, \qquad
#'       y^\ast_{ij} = y_{ij} - \bar y_i}
#' and computes
#' \deqn{r_{\mathrm{rm}} =
#'       \frac{\sum_i \sum_j x^\ast_{ij} y^\ast_{ij}}
#'            {\sqrt{\left(\sum_i \sum_j {x^\ast_{ij}}^2\right)
#'                   \left(\sum_i \sum_j {y^\ast_{ij}}^2\right)}}.}
#'
#' Equivalently, this is the correlation implied by an ANCOVA model with a
#' common slope and subject-specific intercepts:
#' \deqn{y_{ij} = \alpha_i + \beta x_{ij} + \varepsilon_{ij}.}
#' The returned slope is
#' \deqn{\hat\beta =
#'       \frac{\sum_i \sum_j x^\ast_{ij} y^\ast_{ij}}
#'            {\sum_i \sum_j {x^\ast_{ij}}^2},}
#' and the subject-specific fitted intercepts are
#' \eqn{\hat\alpha_i = \bar y_i - \hat\beta \bar x_i}. Residual degrees of
#' freedom are \eqn{N - S - 1}, where \eqn{N = \sum_i n_i} after filtering to
#' complete observations and retaining only subjects with at least two repeated
#' pairs.
#'
#' Confidence intervals are computed with a Fisher \eqn{z}-transformation of
#' \eqn{r_{\mathrm{rm}}} and then back-transformed to the correlation scale. In
#' matrix mode, the same estimator is applied to every pair of selected
#' response columns.
#'
#' @references
#' Bakdash, J. Z., & Marusich, L. R. (2017). Repeated Measures Correlation.
#' \emph{Frontiers in Psychology}, 8, 456.
#' \doi{10.3389/fpsyg.2017.00456}
#'
#' @examples
#' set.seed(2026)
#' n_subjects <- 20
#' n_rep <- 4
#' subject <- rep(seq_len(n_subjects), each = n_rep)
#' subj_eff_x <- rnorm(n_subjects, sd = 1.5)
#' subj_eff_y <- rnorm(n_subjects, sd = 2.0)
#' within_signal <- rnorm(n_subjects * n_rep)
#'
#' dat <- data.frame(
#'   subject = subject,
#'   x = subj_eff_x[subject] + within_signal + rnorm(n_subjects * n_rep, sd = 0.2),
#'   y = subj_eff_y[subject] + 0.8 * within_signal + rnorm(n_subjects * n_rep, sd = 0.3),
#'   z = subj_eff_y[subject] - 0.4 * within_signal + rnorm(n_subjects * n_rep, sd = 0.4)
#' )
#'
#' fit_xy <- rmcorr(dat, response = c("x", "y"), subject = "subject", keep_data = TRUE)
#' print(fit_xy)
#' summary(fit_xy)
#' plot(fit_xy)
#'
#' fit_mat <- rmcorr(dat, response = c("x", "y", "z"), subject = "subject")
#' print(fit_mat, digits = 3)
#' summary(fit_mat)
#' plot(fit_mat)
#'
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   fit_mat_view <- rmcorr(
#'     dat,
#'     response = c("x", "y", "z"),
#'     subject = "subject",
#'     keep_data = TRUE
#'   )
#'   view_rmcorr_shiny(fit_mat_view)
#' }
#' @author Thiago de Paula Oliveira
#' @export
rmcorr <- function(data = NULL, response, subject,
                   na_method = c("error", "pairwise"), conf_level = 0.95,
                   n_threads = getOption("matrixCorr.threads", 1L),
                   keep_data = FALSE,
                   verbose = FALSE,
                   ...) {
  # legacy positional signature support: rmcorr(response, subject, ...)
  if (!missing(data) && missing(subject)) {
    data_nrow <- tryCatch(NROW(data), error = function(...) NA_integer_)
    response_is_subject_vector <-
      !missing(response) &&
      isTRUE(length(response) == data_nrow) &&
      !(
        is.character(response) &&
        length(response) >= 2L &&
        inherits(data, "data.frame") &&
        all(response %in% names(data))
      )

    if (!inherits(data, "data.frame") || response_is_subject_vector) {
      subject <- response
      response <- data
      data <- NULL
    }
  }

  if (...length() == 0L &&
      missing(na_method) &&
      isFALSE(keep_data) &&
      isFALSE(verbose)) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
    n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")

    resolved <- .mc_rmcorr_resolve_inputs(data = data, response = response, subject = subject)
    response_mat <- resolved$response
    subject_vec <- resolved$subject
    response_names <- colnames(response_mat)

    if (ncol(response_mat) < 2L) {
      abort_bad_arg(
        "response",
        message = "must supply at least two numeric response columns."
      )
    }
    if (nrow(response_mat) != length(subject_vec)) {
      abort_bad_arg(
        "response",
        message = "and {.arg subject} must have the same number of observations."
      )
    }

    .mc_check_latent_missing(
      list(response = response_mat, subject = subject_vec),
      check_na = TRUE,
      arg = "response"
    )

    subj_info <- .mc_rmcorr_encode_subject(subject_vec, arg = "subject")

    if (ncol(response_mat) == 2L) {
      stats <- rmcorr_pair_cpp(
        x = response_mat[, 1L],
        y = response_mat[, 2L],
        subject = subj_info$code,
        conf_level = conf_level
      )
      return(.mc_rmcorr_build_pair_object(
        stats = stats,
        response_mat = response_mat,
        subject = subject_vec,
        response_names = response_names,
        subject_name = resolved$subject_name,
        na_method = "error",
        conf_level = conf_level,
        keep_data = FALSE
      ))
    }

    prev_threads <- get_omp_threads()
    on.exit(set_omp_threads(as.integer(prev_threads)), add = TRUE)

    out <- rmcorr_matrix_cpp(
      x = response_mat,
      y = response_mat,
      subject = subj_info$code,
      symmetric = TRUE,
      conf_level = conf_level,
      n_threads = n_threads
    )

    est <- out$estimate
    dimnames(est) <- list(response_names, response_names)
    diagnostics <- list(
      slope = .mc_rmcorr_set_dimnames(out$slope, response_names),
      p_value = .mc_rmcorr_set_dimnames(out$p_value, response_names),
      df = .mc_rmcorr_set_dimnames(out$df, response_names),
      n_complete = .mc_rmcorr_set_dimnames(out$n_complete, response_names),
      n_subjects = .mc_rmcorr_set_dimnames(out$n_subjects, response_names),
      conf_low = .mc_rmcorr_set_dimnames(out$conf_low, response_names),
      conf_high = .mc_rmcorr_set_dimnames(out$conf_high, response_names),
      conf_level = conf_level
    )
    return(.mc_structure_corr_matrix(
      est,
      class_name = "rmcorr_matrix",
      method = "rmcorr",
      description = "Repeated-measures correlation matrix",
      diagnostics = diagnostics
    ))
  }

  legacy_args <- .mc_extract_legacy_aliases(list(...), allowed = "check_na")
  na_cfg <- resolve_na_args(
    na_method = na_method,
    check_na = legacy_args$check_na %||% NULL,
    na_method_missing = missing(na_method)
  )

  check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")
  check_bool(keep_data, arg = "keep_data")
  check_bool(verbose, arg = "verbose")

  resolved <- .mc_rmcorr_resolve_inputs(data = data, response = response, subject = subject)
  response_mat <- resolved$response
  subject_vec <- resolved$subject
  response_names <- colnames(response_mat)

  if (ncol(response_mat) < 2L) {
    abort_bad_arg(
      "response",
      message = "must supply at least two numeric response columns."
    )
  }

  if (nrow(response_mat) != length(subject_vec)) {
    abort_bad_arg(
      "response",
      message = "and {.arg subject} must have the same number of observations."
    )
  }

  if (na_cfg$check_na) {
    .mc_check_latent_missing(
      list(response = response_mat, subject = subject_vec),
      check_na = TRUE,
      arg = "response"
    )
  }

  subj_info <- .mc_rmcorr_encode_subject(subject_vec, arg = "subject")

  if (ncol(response_mat) == 2L) {
    stats <- rmcorr_pair_cpp(
      x = response_mat[, 1L],
      y = response_mat[, 2L],
      subject = subj_info$code,
      conf_level = conf_level
    )
    return(.mc_rmcorr_build_pair_object(
      stats = stats,
      response_mat = response_mat,
      subject = subject_vec,
      response_names = response_names,
      subject_name = resolved$subject_name,
      na_method = na_cfg$na_method,
      conf_level = conf_level,
      keep_data = keep_data
    ))
  }

  if (verbose) {
    cat("Using", n_threads, "OpenMP threads\n")
  }

  prev_threads <- get_omp_threads()
  on.exit(set_omp_threads(as.integer(prev_threads)), add = TRUE)

  out <- rmcorr_matrix_cpp(
    x = response_mat,
    y = response_mat,
    subject = subj_info$code,
    symmetric = TRUE,
    conf_level = conf_level,
    n_threads = n_threads
  )

  est <- out$estimate
  dimnames(est) <- list(response_names, response_names)

  diagnostics <- list(
    slope = .mc_rmcorr_set_dimnames(out$slope, response_names),
    p_value = .mc_rmcorr_set_dimnames(out$p_value, response_names),
    df = .mc_rmcorr_set_dimnames(out$df, response_names),
    n_complete = .mc_rmcorr_set_dimnames(out$n_complete, response_names),
    n_subjects = .mc_rmcorr_set_dimnames(out$n_subjects, response_names),
    conf_low = .mc_rmcorr_set_dimnames(out$conf_low, response_names),
    conf_high = .mc_rmcorr_set_dimnames(out$conf_high, response_names),
    conf_level = conf_level
  )

  out_obj <- .mc_structure_corr_matrix(
    est,
    class_name = "rmcorr_matrix",
    method = "rmcorr",
    description = "Repeated-measures correlation matrix",
    diagnostics = diagnostics
  )
  if (keep_data) {
    attr(out_obj, "source_data") <- .mc_rmcorr_compact_source_data(
      response_mat = response_mat,
      subject = subject_vec,
      subject_name = resolved$subject_name,
      conf_level = conf_level
    )
  }
  out_obj
}

.mc_rmcorr_resolve_inputs <- function(data, response, subject) {
  pull_data_col <- function(df, x, arg) {
    if (is.character(x) && length(x) >= 1L) {
      missing_cols <- setdiff(x, names(df))
      if (length(missing_cols) > 0L) {
        abort_bad_arg(
          "data",
          message = "is missing required column(s): {paste(missing_cols, collapse = ', ')}.",
          missing_cols = missing_cols
        )
      }
      if (identical(arg, "response")) {
        return(df[x])
      }
      if (length(x) != 1L) {
        abort_bad_arg(
          arg,
          message = "must be a vector or a single column name in {.arg data}."
        )
      }
      return(df[[x]])
    }
    x
  }

  if (!is.null(data)) {
    if (!inherits(data, "data.frame")) {
      data <- as.data.frame(data, stringsAsFactors = FALSE)
    }
    response <- pull_data_col(data, response, "response")
    subject_name <-
      if (is.character(subject) && length(subject) == 1L) subject else "subject"
    subject <- pull_data_col(data, subject, "subject")
  } else {
    subject_name <- "subject"
  }

  response_mat <- .mc_extract_continuous_matrix(response, arg = "response", min_cols = 2L)

  if (!is.atomic(subject) && !is.factor(subject)) {
    abort_bad_arg(
      "subject",
      message = "must be a vector or a single column name in {.arg data}."
    )
  }

  list(
    response = response_mat,
    subject = subject,
    subject_name = subject_name
  )
}

.mc_rmcorr_encode_subject <- function(subject, arg = "subject") {
  if (is.numeric(subject)) {
    bad <- is.na(subject) | is.nan(subject) | is.infinite(subject)
    if (any(bad)) subject[bad] <- NA_real_
  }
  f <- factor(subject)
  if (nlevels(f) < 2L) {
    abort_bad_arg(
      arg,
      message = "must contain at least two distinct subjects."
    )
  }
  list(code = as.integer(f), levels = levels(f))
}

.mc_rmcorr_pair_data <- function(response_mat, subject) {
  dat <- data.frame(
    .response1 = as.numeric(response_mat[, 1L]),
    .response2 = as.numeric(response_mat[, 2L]),
    .subject = subject,
    check.names = FALSE
  )
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  if (!nrow(dat)) return(dat)
  counts <- table(dat$.subject)
  keep <- names(counts[counts >= 2L])
  dat <- dat[dat$.subject %in% keep, , drop = FALSE]
  dat$.subject <- droplevels(factor(dat$.subject))
  rownames(dat) <- NULL
  dat
}

.mc_rmcorr_build_pair_object <- function(stats, response_mat, subject,
                                         response_names, subject_name,
                                         na_method,
                                         conf_level, keep_data = FALSE) {
  dat <- .mc_rmcorr_pair_data(response_mat, subject)
  slope <- as.numeric(stats$slope)
  valid <- isTRUE(stats$valid) && is.finite(stats$estimate)

  subj_counts <- if (nrow(dat)) as.integer(table(dat$.subject)) else integer()
  out <- structure(
    list(
      estimate = as.numeric(stats$estimate),
      p_value = as.numeric(stats$p_value),
      lwr = as.numeric(stats$conf_low),
      upr = as.numeric(stats$conf_high),
      conf_level = conf_level,
      slope = slope,
      df = as.numeric(stats$df),
      n_obs = as.integer(stats$n_complete),
      n_subjects = as.integer(stats$n_subjects),
      na_method = na_method,
      responses = response_names,
      subject_name = subject_name,
      obs_per_subject_min = if (length(subj_counts)) min(subj_counts) else NA_integer_,
      obs_per_subject_max = if (length(subj_counts)) max(subj_counts) else NA_integer_,
      valid = valid
    ),
    class = "rmcorr",
    method = "rmcorr",
    description = "Repeated-measures correlation",
    package = "matrixCorr"
  )
  if (isTRUE(keep_data)) {
    attr(out, "source_data") <- .mc_rmcorr_compact_source_data(
      response_mat = response_mat,
      subject = subject,
      subject_name = subject_name,
      conf_level = conf_level
    )
  }
  out
}

.mc_rmcorr_set_dimnames <- function(x, nm) {
  x <- as.matrix(x)
  dimnames(x) <- list(nm, nm)
  x
}

.mc_rmcorr_pick_diag_vals <- function(diag_attr, object, field) {
  x <- diag_attr[[field]]
  m <- as.matrix(object)
  if (!is.matrix(x) || !identical(dim(x), dim(m))) return(numeric())
  selector <- if (isTRUE(nrow(m) == ncol(m)) && nrow(m) > 1L) {
    upper.tri(m, diag = FALSE)
  } else {
    matrix(TRUE, nrow(m), ncol(m))
  }
  vals <- x[selector]
  vals[is.finite(vals)]
}

.mc_rmcorr_compact_source_data <- function(response_mat, subject, subject_name,
                                           conf_level) {
  f <- factor(subject)
  list(
    response = unname(response_mat),
    response_names = colnames(response_mat),
    subject_code = as.integer(f),
    subject_levels = levels(f),
    subject_name = subject_name,
    conf_level = conf_level
  )
}

.mc_rmcorr_pair_plot_payload <- function(x) {
  source_data <- attr(x, "source_data", exact = TRUE)
  if (!is.list(source_data) || is.null(source_data$response) || is.null(source_data$subject_code)) {
    return(NULL)
  }

  subject <- factor(
    source_data$subject_code,
    levels = seq_along(source_data$subject_levels),
    labels = source_data$subject_levels
  )
  dat <- .mc_rmcorr_pair_data(source_data$response, subject)
  if (!nrow(dat)) {
    return(list(data_long = dat, intercepts = numeric(), fitted = dat[0, c(".subject"), drop = FALSE]))
  }

  slope <- as.numeric(x$estimate) * 0 + as.numeric(x$slope)
  intercepts <- tapply(
    dat$.response2 - slope * dat$.response1,
    dat$.subject,
    mean
  )
  fitted <- do.call(
    rbind,
    lapply(names(intercepts), function(s) {
      idx <- dat$.subject == s
      x_min <- min(dat$.response1[idx])
      x_max <- max(dat$.response1[idx])
      b0 <- as.numeric(intercepts[[s]])
      data.frame(
        .subject = s,
        x = c(x_min, x_max),
        y = c(b0 + slope * x_min, b0 + slope * x_max),
        check.names = FALSE
      )
    })
  )
  list(data_long = dat, intercepts = intercepts, fitted = fitted)
}

.mc_rmcorr_component_names <- function(x) {
  out <- c(
    "estimate", "p_value", "lwr", "upr", "conf_level", "slope", "df",
    "n_obs", "n_subjects", "na_method", "responses", "subject_name",
    "obs_per_subject_min", "obs_per_subject_max", "valid"
  )
  out <- c(out, "r", "conf_int", "based.on")
  if (!is.null(attr(x, "source_data", exact = TRUE))) {
    out <- c(out, "source_data", "data_long", "intercepts", "fitted")
  }
  unique(out)
}

#' @export
names.rmcorr <- function(x) {
  .mc_rmcorr_component_names(x)
}

#' @export
`$.rmcorr` <- function(x, name) {
  core <- unclass(x)
  switch(name,
    r = core$estimate,
    conf_int = c(lower = core$lwr, upper = core$upr),
    based.on = core$n_obs,
    source_data = attr(x, "source_data", exact = TRUE),
    data_long = {
      payload <- .mc_rmcorr_pair_plot_payload(x)
      if (is.null(payload)) NULL else payload$data_long
    },
    intercepts = {
      payload <- .mc_rmcorr_pair_plot_payload(x)
      if (is.null(payload)) NULL else payload$intercepts
    },
    fitted = {
      payload <- .mc_rmcorr_pair_plot_payload(x)
      if (is.null(payload)) NULL else payload$fitted
    },
    core[[name]]
  )
}

#' @export
`[[.rmcorr` <- function(x, i, ...) {
  if (is.numeric(i)) {
    nms <- names(x)
    if (length(i) != 1L || is.na(i) || i < 1L || i > length(nms)) {
      abort_bad_arg("i", message = "must be a valid component index.")
    }
    i <- nms[[i]]
  }
  `$.rmcorr`(x, i)
}

#' @title Methods for Pairwise rmcorr Objects
#' @description Print, summarize, and plot methods for pairwise
#'   repeated-measures correlation objects of class \code{"rmcorr"} and
#'   \code{"summary.rmcorr"}.
#' @param x An object of class \code{"rmcorr"} or \code{"summary.rmcorr"}.
#' @param object An object of class \code{"rmcorr"}.
#' @param digits Number of significant digits to print.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param title Optional plot title for \code{plot.rmcorr()}. Defaults to a
#'   title containing the estimated repeated-measures correlation.
#' @param point_alpha Alpha transparency for scatterplot points.
#' @param line_width Line width for subject-specific fitted lines.
#' @param show_legend Logical; if \code{TRUE}, shows the subject legend in the
#'   pairwise scatterplot.
#' @param show_value Logical; included for a consistent plotting interface.
#'   Pairwise repeated-measures plots do not overlay numeric cell values, so this
#'   argument currently has no effect.
#' @param ... Additional arguments passed to downstream methods. For
#'   \code{plot.rmcorr()}, these are passed to \code{ggplot2::theme()}.
#' @rdname rmcorr_methods
#' @method print rmcorr
#' @export
print.rmcorr <- function(x, digits = 4, n = NULL, topn = NULL,
                         max_vars = NULL, width = NULL,
                         show_ci = NULL, ...) {
  check_inherits(x, "rmcorr")
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("print_show_ci", "yes")
  )
  digest <- c(
    method = attr(x, "method"),
    responses = sprintf("%s vs %s", x$responses[[1L]], x$responses[[2L]]),
    n_obs = x$n_obs,
    n_subjects = x$n_subjects,
    estimate = format(signif(x$estimate, digits = digits)),
    if (identical(show_ci, "yes")) {
      c(ci = sprintf(
        "[%s, %s] (%g%%)",
        format(signif(x$lwr, digits = digits)),
        format(signif(x$upr, digits = digits)),
        100 * x$conf_level
      ))
    }
  )
  .mc_print_named_digest(digest, header = "Repeated-measures correlation preview:")
  invisible(x)
}

#' @rdname rmcorr_methods
#' @method summary rmcorr
#' @export
summary.rmcorr <- function(object, n = NULL, topn = NULL,
                           max_vars = NULL, width = NULL,
                           show_ci = NULL, ...) {
  check_inherits(object, "rmcorr")
  out <- list(
    class = "rmcorr",
    method = attr(object, "method"),
    description = attr(object, "description"),
    responses = object$responses,
    n_obs = object$n_obs,
    n_subjects = object$n_subjects,
    df = object$df,
    estimate = object$estimate,
    slope = object$slope,
    p_value = object$p_value,
    lwr = object$lwr,
    upr = object$upr,
    conf_level = object$conf_level,
    obs_per_subject_min = object$obs_per_subject_min %||% NA_integer_,
    obs_per_subject_max = object$obs_per_subject_max %||% NA_integer_,
    n_intercepts = object$n_subjects %||% NA_integer_,
    ci_method = "fisher_z",
    ci_width = if (all(is.finite(c(object$lwr, object$upr)))) object$upr - object$lwr else NA_real_,
    valid = isTRUE(object$valid)
  )
  .mc_finalize_summary_list(out, class_name = "summary.rmcorr")
}

#' @rdname rmcorr_methods
#' @method print summary.rmcorr
#' @export
print.summary.rmcorr <- function(x, digits = 4, n = NULL, topn = NULL,
                                 max_vars = NULL, width = NULL,
                                 show_ci = NULL, ...) {
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )
  digest <- c(
    method = x$method,
    responses = sprintf("%s vs %s", x$responses[[1L]], x$responses[[2L]]),
    n_obs = x$n_obs,
    n_subjects = x$n_subjects,
    obs_per_subject = .mc_format_scalar_or_range(
      x$obs_per_subject_min,
      x$obs_per_subject_max,
      digits = digits,
      integer = TRUE
    ),
    df = format(signif(x$df, digits = digits)),
    estimate = format(signif(x$estimate, digits = digits)),
    slope = format(signif(x$slope, digits = digits)),
    p_value = format(signif(x$p_value, digits = digits)),
    if (identical(show_ci, "yes")) {
      c(
        interval = sprintf(
          "[%s, %s]",
          format(signif(x$lwr, digits = digits)),
          format(signif(x$upr, digits = digits))
        ),
        ci_level = sprintf("%g%%", 100 * x$conf_level),
        ci_method = x$ci_method,
        ci_width = format(signif(x$ci_width, digits = digits))
      )
    }
  )
  if (is.finite(x$n_intercepts) && x$n_intercepts > 0L) {
    digest <- c(digest, fitted_subjects = x$n_intercepts)
  }
  if (!isTRUE(x$valid)) {
    digest <- c(digest, valid = "no")
  }
  digest <- digest[!is.na(digest) & nzchar(digest)]
  .mc_print_named_digest(digest, header = "Repeated-measures correlation summary:")
  invisible(x)
}

#' @rdname rmcorr_methods
#' @method plot rmcorr
#' @export
plot.rmcorr <- function(x,
                        title = NULL,
                        point_alpha = 0.8,
                        line_width = 0.8,
                        show_legend = FALSE,
                        show_value = TRUE,
                        ...) {
  check_inherits(x, "rmcorr")
  check_bool(show_value, arg = "show_value")
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} is required for plotting.")
  }
  payload <- .mc_rmcorr_pair_plot_payload(x)
  if (is.null(payload)) {
    cli::cli_abort(
      c(
        "The repeated-measures correlation fit does not retain plotting data.",
        "i" = "Refit with {.code keep_data = TRUE} to enable {.fn plot} for pairwise repeated-measures correlations."
      )
    )
  }
  if (!nrow(payload$data_long) || !isTRUE(x$valid)) {
    cli::cli_abort("The repeated-measures correlation fit is not valid for plotting.")
  }

  title <- title %||% sprintf(
    "Repeated-measures correlation: r = %.2f",
    x$r
  )

  p <- ggplot2::ggplot(
    payload$data_long,
    ggplot2::aes(
      x = .data$.response1,
      y = .data$.response2,
      colour = .data$.subject
    )
  ) +
    ggplot2::geom_point(alpha = point_alpha) +
    ggplot2::geom_line(
      data = payload$fitted,
      ggplot2::aes(x = .data$x, y = .data$y, colour = .data$.subject, group = .data$.subject),
      linewidth = line_width,
      inherit.aes = FALSE
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = if (show_legend) "right" else "none",
      ...
    ) +
    ggplot2::labs(
      title = title,
      x = x$responses[[1L]],
      y = x$responses[[2L]],
      colour = x$subject_name
    )

  p
}

#' @title Methods for rmcorr Matrix Objects
#' @description Print, summarize, and plot methods for repeated-measures
#'   correlation matrix objects of class \code{"rmcorr_matrix"} and
#'   \code{"summary.rmcorr_matrix"}.
#' @param x An object of class \code{"rmcorr_matrix"} or
#'   \code{"summary.rmcorr_matrix"}.
#' @param object An object of class \code{"rmcorr_matrix"}.
#' @param digits Number of significant digits to print.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param title Plot title for \code{plot.rmcorr_matrix()}.
#' @param low_color,high_color,mid_color Colours used for negative, positive,
#'   and midpoint values in the heatmap.
#' @param value_text_size Size of the overlaid numeric value labels in the
#'   heatmap.
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on heatmap tiles when the plot type supports them.
#' @param ... Additional arguments passed to downstream methods.
#' @rdname rmcorr_matrix_methods
#' @method print rmcorr_matrix
#' @export
print.rmcorr_matrix <- function(x, digits = 4, n = NULL,
                                topn = NULL, max_vars = NULL,
                                width = NULL, show_ci = NULL, ...) {
  .mc_print_corr_matrix(
    x,
    header = "Repeated-measures correlation matrix",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
}

#' @rdname rmcorr_matrix_methods
#' @method plot rmcorr_matrix
#' @export
plot.rmcorr_matrix <- function(x, title = "Repeated-measures correlation heatmap",
                               low_color = "indianred1",
                               high_color = "steelblue1",
                               mid_color = "white",
                               value_text_size = 4,
                               show_value = TRUE, ...) {
  .mc_plot_corr_matrix(
    x,
    class_name = "rmcorr_matrix",
    fill_name = "rmcorr",
    title = title,
    low_color = low_color,
    high_color = high_color,
    mid_color = mid_color,
    value_text_size = value_text_size,
    show_value = show_value,
    ...
  )
}

#' @rdname rmcorr_matrix_methods
#' @method summary rmcorr_matrix
#' @export
summary.rmcorr_matrix <- function(object, n = NULL, topn = NULL,
                                  max_vars = NULL, width = NULL,
                                  show_ci = NULL, ...) {
  check_inherits(object, "rmcorr_matrix")
  diag_attr <- attr(object, "diagnostics")
  m <- as.matrix(object)
  symmetric <- isTRUE(nrow(m) == ncol(m)) && isTRUE(isSymmetric(unclass(m)))
  vals <- if (symmetric && nrow(m) > 1L) {
    m[upper.tri(m, diag = FALSE)]
  } else {
    as.vector(m)
  }
  vals <- vals[is.finite(vals)]
  n_pairs <- if (symmetric) {
    stats::setNames(nrow(m) * (nrow(m) - 1L) / 2L, NULL)
  } else {
    stats::setNames(nrow(m) * ncol(m), NULL)
  }

  out <- list(
    class = class(object)[1L],
    method = attr(object, "method"),
    description = attr(object, "description"),
    n_variables = if (symmetric) nrow(m) else NA_integer_,
    n_pairs = as.integer(n_pairs),
    n_rows = nrow(m),
    n_cols = ncol(m),
    symmetric = symmetric,
    n_missing = sum(is.na(m)),
    estimate_min = if (length(vals)) min(vals) else NA_real_,
    estimate_max = if (length(vals)) max(vals) else NA_real_,
    header = "Repeated-measures correlation summary"
  )
  out$df_min <- {
    vals <- .mc_rmcorr_pick_diag_vals(diag_attr, object, "df")
    if (length(vals)) min(vals) else NA_real_
  }
  out$df_max <- {
    vals <- .mc_rmcorr_pick_diag_vals(diag_attr, object, "df")
    if (length(vals)) max(vals) else NA_real_
  }
  out$n_complete_min <- {
    vals <- .mc_rmcorr_pick_diag_vals(diag_attr, object, "n_complete")
    if (length(vals)) min(vals) else NA_real_
  }
  out$n_complete_max <- {
    vals <- .mc_rmcorr_pick_diag_vals(diag_attr, object, "n_complete")
    if (length(vals)) max(vals) else NA_real_
  }
  out$n_subjects_min <- {
    vals <- .mc_rmcorr_pick_diag_vals(diag_attr, object, "n_subjects")
    if (length(vals)) min(vals) else NA_real_
  }
  out$n_subjects_max <- {
    vals <- .mc_rmcorr_pick_diag_vals(diag_attr, object, "n_subjects")
    if (length(vals)) max(vals) else NA_real_
  }
  out$p_value_min <- {
    vals <- .mc_rmcorr_pick_diag_vals(diag_attr, object, "p_value")
    if (length(vals)) min(vals) else NA_real_
  }
  out$p_value_max <- {
    vals <- .mc_rmcorr_pick_diag_vals(diag_attr, object, "p_value")
    if (length(vals)) max(vals) else NA_real_
  }
  out$conf_level <- diag_attr$conf_level %||% NA_real_
  out$top_results <- .mc_summary_top_pairs(object, digits = 4, topn = .mc_coalesce(topn, .mc_display_option("summary_topn", 5L)))
  .mc_finalize_summary_list(out, class_name = "summary.rmcorr_matrix")
}

#' @rdname rmcorr_matrix_methods
#' @method print summary.rmcorr_matrix
#' @export
print.summary.rmcorr_matrix <- function(x, digits = 4, n = NULL,
                                        topn = NULL, max_vars = NULL,
                                        width = NULL, show_ci = NULL, ...) {
  cfg <- .mc_resolve_display_args(
    context = "summary",
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci
  )
  digest <- c(
    method = x$method,
    dimensions = sprintf("%d x %d", x$n_rows, x$n_cols),
    pairs = .mc_count_fmt(x$n_pairs),
    estimate = .mc_format_scalar_or_range(x$estimate_min, x$estimate_max, digits = digits)
  )
  if (!is.na(x$n_variables)) digest <- c(digest, n_variables = x$n_variables)
  if (is.finite(x$conf_level)) digest <- c(digest, conf_level = format(signif(x$conf_level, digits = digits)))
  if (is.finite(x$n_complete_min) && is.finite(x$n_complete_max)) {
    digest <- c(digest, n_complete = .mc_format_scalar_or_range(x$n_complete_min, x$n_complete_max, digits = digits))
  }
  if (is.finite(x$df_min) && is.finite(x$df_max)) {
    digest <- c(digest, df = .mc_format_scalar_or_range(x$df_min, x$df_max, digits = digits))
  }
  if (is.finite(x$n_subjects_min) && is.finite(x$n_subjects_max)) {
    digest <- c(digest, n_subjects = .mc_format_scalar_or_range(x$n_subjects_min, x$n_subjects_max, digits = digits))
  }
  if (is.finite(x$p_value_min) && is.finite(x$p_value_max)) {
    digest <- c(digest, p_value = .mc_format_scalar_or_range(x$p_value_min, x$p_value_max, digits = digits))
  }
  if (isTRUE(x$n_missing > 0L)) {
    digest <- c(digest, missing = .mc_count_fmt(x$n_missing))
  }
  .mc_print_named_digest(digest, header = x$header %||% "Repeated-measures correlation summary")
  if (is.data.frame(x$top_results) && nrow(x$top_results)) {
    .mc_print_ranked_pairs_preview(
      x$top_results,
      header = "Strongest pairs by |estimate|",
      topn = cfg$topn,
      max_vars = cfg$max_vars,
      width = cfg$width,
      show_ci = cfg$show_ci,
      ...
    )
  }
  invisible(x)
}

print.summary_rmcorr <- print.summary.rmcorr
print.summary_rmcorr_matrix <- print.summary.rmcorr_matrix
