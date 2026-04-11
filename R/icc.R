#' @title Intraclass Correlation for Wide Data
#'
#' @description
#' Computes intraclass correlation coefficients for the numeric columns of a
#' matrix or data frame using the classical ANOVA mean-square formulas. The
#' output can be either a pairwise matrix across columns or an overall
#' all-column coefficient table.
#'
#' @details
#' Each column is treated as a measurement channel, method, or rater and each
#' row is treated as a subject.
#'
#' The function supports two distinct analysis targets.
#'
#' - `scope = "pairwise"` answers: "how reliable is each specific column pair?"
#'   Each estimate is based on exactly two columns and the output is a symmetric
#'   matrix.
#' - `scope = "overall"` answers: "how reliable is the full set of columns when
#'   analysed jointly?" The output is the standard six-form overall ANOVA table
#'   (`ICC1`, `ICC2`, `ICC3`, `ICC1k`, `ICC2k`, `ICC3k`).
#'
#' These two scopes do not target the same quantity. The overall coefficients
#' are computed from the full multi-column ANOVA decomposition and are not
#' obtained by averaging or otherwise aggregating the pairwise matrix.
#'
#' The three main choice arguments determine the classical ICC form.
#'
#' - `model` controls the rater structure:
#'   - `"oneway"` uses the one-way random-effects formulation.
#'   - `"twoway_random"` uses the two-way random-effects formulation.
#'   - `"twoway_mixed"` uses the two-way mixed-effects formulation.
#' - `type` controls whether systematic column mean differences are penalized:
#'   - `"consistency"` targets consistency across columns.
#'   - `"agreement"` targets absolute agreement across columns.
#' - `unit` controls whether reliability refers to one measurement or to the
#'   average of multiple measurements:
#'   - `"single"` returns the single-measure coefficient.
#'   - `"average"` returns the average-measure coefficient.
#'
#' The supported mappings are:
#'
#' - `model = "oneway", type = "consistency", unit = "single"` gives `ICC1`.
#' - `model = "oneway", type = "consistency", unit = "average"` gives `ICC1k`.
#' - `model = "twoway_random", type = "agreement", unit = "single"` gives `ICC2`.
#' - `model = "twoway_random", type = "agreement", unit = "average"` gives `ICC2k`.
#' - `model = "twoway_random", type = "consistency", unit = "single"` gives `ICC3`.
#' - `model = "twoway_random", type = "consistency", unit = "average"` gives `ICC3k`.
#' - `model = "twoway_mixed", type = "agreement", unit = "single"` gives the
#'   mixed-effects absolute-agreement analogue with the same classical point
#'   formula as `ICC2`.
#' - `model = "twoway_mixed", type = "agreement", unit = "average"` gives the
#'   corresponding average-measure analogue.
#' - `model = "twoway_mixed", type = "consistency", unit = "single"` gives `ICC3`.
#' - `model = "twoway_mixed", type = "consistency", unit = "average"` gives `ICC3k`.
#'
#' The combination `model = "oneway", type = "agreement"` is not defined here
#' and returns an error.
#'
#' For `scope = "pairwise"`, the point estimates are computed in `C++` directly
#' from the two-column ANOVA mean squares for each complete pair. For
#' `unit = "average"`, the implementation uses `k = 2` because each estimate is
#' based on exactly two columns.
#'
#' For `scope = "overall"`, the point estimates are computed jointly from the
#' full wide matrix using the classical ANOVA decomposition over all columns.
#' Here the average-measure coefficients use `k = ncol(data)` after any row
#' filtering required by `check_na`.
#'
#' Missing-data handling depends on `scope`:
#'
#' - with `check_na = TRUE`, missing values are rejected before estimation;
#' - with `check_na = FALSE` and `scope = "pairwise"`, each pair uses its own
#'   complete-case overlap;
#' - with `check_na = FALSE` and `scope = "overall"`, rows are restricted to
#'   complete cases across all columns because the overall ANOVA requires a
#'   common wide matrix.
#'
#' When `ci = TRUE`, confidence intervals are obtained from the classical
#' F-based ANOVA formulas corresponding to the selected coefficient. For
#' `scope = "pairwise"`, non-estimable off-diagonal pairs return `NA`. For
#' `scope = "overall"`, the coefficient table includes interval columns for all
#' six standard rows.
#'
#' @param data A numeric matrix or data frame with at least two numeric columns.
#' @param model Character scalar selecting the classical ICC model.
#'   `"oneway"` uses the one-way random-effects formulation,
#'   `"twoway_random"` uses the two-way random-effects formulation, and
#'   `"twoway_mixed"` uses the two-way mixed-effects formulation.
#' @param type Character scalar selecting the reliability target.
#'   `"consistency"` does not penalize systematic column mean differences in
#'   the same way as absolute agreement, whereas `"agreement"` targets absolute
#'   agreement and therefore incorporates those differences.
#' @param unit Character scalar selecting whether the coefficient refers to a
#'   single measurement (`"single"`) or to the mean of multiple measurements
#'   (`"average"`). For `scope = "pairwise"`, the average-measure coefficient
#'   always uses `k = 2`. For `scope = "overall"`, it uses the full number of
#'   analysed columns.
#' @param scope Character scalar selecting the analysis target.
#'   `"pairwise"` returns the symmetric matrix of two-column ICCs.
#'   `"overall"` returns the six-row overall ANOVA table for the full set of
#'   columns. Choose `"pairwise"` when the question is about specific column
#'   pairs; choose `"overall"` when the question is about the full set of
#'   columns analysed jointly.
#' @param ci Logical; if `TRUE`, return confidence intervals.
#' @param conf_level Confidence level for the interval output. Ignored when
#'   `ci = FALSE`.
#' @param check_na Logical. If `TRUE`, missing values are rejected before
#'   estimation. If `FALSE`, `scope = "pairwise"` uses pair-specific complete
#'   cases, while `scope = "overall"` uses complete rows across all analysed
#'   columns.
#' @param n_threads Integer number of OpenMP threads.
#' @param verbose Logical; if `TRUE`, report how many threads are requested.
#'
#' @return
#' For `scope = "pairwise"`, if `ci = FALSE`, a symmetric matrix of class
#' `icc`. If `ci = TRUE`, a list with elements `est`, `lwr.ci`, and `upr.ci`
#' and class `c("icc", "icc_ci")`.
#'
#' For `scope = "overall"`, a list of class `c("icc_overall", "icc")` with a
#' coefficient table, ANOVA table, and mean-square metadata. The coefficient
#' table always includes the standard six overall coefficients. Confidence
#' interval columns are attached when `ci = TRUE`.
#'
#' All outputs carry attributes describing the selected model, type, unit, and
#' method.
#'
#' @seealso [ccc()], [ccc_rm_reml()], [ba()], [ba_rm()], [rmcorr()]
#'
#' @references
#' Shrout PE, Fleiss JL (1979). Intraclass correlations: uses in assessing
#' rater reliability. Psychological Bulletin, 86(2), 420-428.
#'
#' McGraw KO, Wong SP (1996). Forming inferences about some intraclass
#' correlation coefficients. Psychological Methods, 1(1), 30-46.
#'
#' @examples
#' set.seed(123)
#' n <- 40
#' subj <- rnorm(n, sd = 1)
#' dat <- data.frame(
#'   m1 = subj + rnorm(n, sd = 0.3),
#'   m2 = 0.2 + subj + rnorm(n, sd = 0.4),
#'   m3 = -0.1 + subj + rnorm(n, sd = 0.5)
#' )
#'
#' fit_icc <- icc(dat,
#'   model = "twoway_random",
#'   type = "agreement",
#'   unit = "single",
#'   scope = "pairwise"
#' )
#' print(fit_icc)
#' summary(fit_icc)
#'
#' fit_icc_overall <- icc(dat, scope = "overall", ci = TRUE)
#' print(fit_icc_overall)
#' summary(fit_icc_overall)
#'
#' @export
icc <- function(data,
                model = c("oneway", "twoway_random", "twoway_mixed"),
                type = c("consistency", "agreement"),
                unit = c("single", "average"),
                scope = c("pairwise", "overall"),
                ci = FALSE,
                conf_level = 0.95,
                check_na = TRUE,
                n_threads = getOption("matrixCorr.threads", 1L),
                verbose = FALSE) {
  model <- match.arg(model)
  type <- match.arg(type)
  unit <- match.arg(unit)
  scope <- match.arg(scope)

  check_bool(ci, arg = "ci")
  check_bool(check_na, arg = "check_na")
  check_bool(verbose, arg = "verbose")
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")
  if (isTRUE(ci)) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  }

  if (identical(model, "oneway") && identical(type, "agreement")) {
    abort_bad_arg("type",
      message = "cannot be {.val agreement} when {.arg model} is {.val oneway}.",
      .hint = "Use type = \"consistency\" for the one-way formulation."
    )
  }

  numeric_data <- validate_corr_input(data, check_na = check_na)
  mat <- as.matrix(numeric_data)
  colnames_data <- colnames(numeric_data)

  if (verbose) cat("Using", n_threads, "OpenMP threads\n")

  form_code <- .mc_icc_form_code(model, type)
  average_unit <- identical(unit, "average")
  selected_coefficient <- .mc_icc_selected_coefficient(model, type, unit)

  if (identical(scope, "overall")) {
    if (check_na) {
      overall_data <- validate_corr_input(data, check_na = TRUE)
    } else {
      overall_data <- validate_corr_input(data, check_na = FALSE)
      keep <- stats::complete.cases(overall_data)
      overall_data <- overall_data[keep, , drop = FALSE]
      if (nrow(overall_data) < 2L) {
        abort_bad_arg("data",
          message = "must retain at least two complete rows across all columns when {.arg scope} is {.val overall}."
        )
      }
    }

    mat <- as.matrix(overall_data)
    if (verbose) cat("Using native ANOVA backend\n")

    fit <- icc_overall_cpp(
      mat,
      return_ci = ci,
      conf_level = conf_level
    )

    coefficients <- fit$coefficients
    coefficients$selected <- coefficients$coefficient == selected_coefficient
    anova <- fit$anova

    out <- structure(
      list(
        coefficients = coefficients,
        anova = anova,
        mean_squares = fit$mean_squares
      ),
      class = c("icc_overall", "icc")
    )

    attr(out, "method") <- "Overall intraclass correlation table"
    attr(out, "description") <- "Overall intraclass correlation coefficients from the classical ANOVA decomposition"
    attr(out, "package") <- "matrixCorr"
    attr(out, "scope") <- scope
    attr(out, "model") <- model
    attr(out, "type") <- type
    attr(out, "unit") <- unit
    attr(out, "selected_coefficient") <- selected_coefficient
    attr(out, "selected_row") <- coefficients[coefficients$selected, , drop = FALSE]
    attr(out, "conf.level") <- conf_level
    attr(out, "ci.method") <- if (isTRUE(ci)) "anova_f" else NULL
    attr(out, "diagnostics") <- list(
      n_complete = nrow(mat),
      n_subjects = fit$n_subjects,
      n_raters = fit$n_raters,
      dropped_rows = if (check_na) 0L else nrow(validate_corr_input(data, check_na = FALSE)) - nrow(mat)
    )
    return(out)
  }

  fit <- icc_matrix_cpp(
    mat,
    form_code = form_code,
    average_unit = average_unit,
    pairwise_complete = !check_na,
    return_ci = ci,
    conf_level = conf_level,
    n_threads = n_threads
  )

  fit$est <- `dimnames<-`(fit$est, list(colnames_data, colnames_data))
  fit$n_complete <- `dimnames<-`(fit$n_complete, list(colnames_data, colnames_data))
  if (isTRUE(ci)) {
    fit$lwr.ci <- `dimnames<-`(fit$lwr.ci, list(colnames_data, colnames_data))
    fit$upr.ci <- `dimnames<-`(fit$upr.ci, list(colnames_data, colnames_data))
  }

  method_label <- .mc_icc_method_label(model = model, type = type, unit = unit)
  description <- paste0("Pairwise intraclass correlation matrix (", method_label, ")")

  if (isTRUE(ci)) {
    out <- structure(
      list(est = fit$est, lwr.ci = fit$lwr.ci, upr.ci = fit$upr.ci),
      class = c("icc", "icc_ci")
    )
  } else {
    out <- structure(fit$est, class = c("icc", "matrix"))
  }

  attr(out, "method") <- method_label
  attr(out, "description") <- description
  attr(out, "package") <- "matrixCorr"
  attr(out, "scope") <- scope
  attr(out, "model") <- model
  attr(out, "type") <- type
  attr(out, "unit") <- unit
  attr(out, "k") <- 2L
  attr(out, "conf.level") <- conf_level
  attr(out, "ci.method") <- if (isTRUE(ci)) "anova_f" else NULL
  attr(out, "diagnostics") <- list(n_complete = fit$n_complete)
  out
}

.mc_icc_form_code <- function(model, type) {
  if (identical(model, "oneway")) {
    return(0L)
  }
  if (identical(type, "agreement")) {
    return(1L)
  }
  2L
}

.mc_icc_method_label <- function(model, type, unit) {
  base <- switch(model,
    oneway = "one-way",
    twoway_random = "two-way random",
    twoway_mixed = "two-way mixed"
  )
  paste("Intraclass correlation", sprintf("(%s, %s, %s)", base, type, unit))
}

.mc_icc_selected_coefficient <- function(model, type, unit) {
  base <- if (identical(model, "oneway")) {
    "ICC1"
  } else if (identical(type, "agreement")) {
    "ICC2"
  } else {
    "ICC3"
  }
  if (identical(unit, "average")) paste0(base, "k") else base
}

#' @rdname icc
#' @method print icc
#' @param x An intraclass-correlation object returned by `icc()`, or a summary
#'   object returned by `summary()` for that fit.
#' @param digits Integer; number of digits to print.
#' @param ci_digits Integer; number of digits for CI bounds.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives
#'   this from console width.
#' @param width Optional display width; defaults to `getOption("width")`.
#' @param show_ci One of `"yes"` or `"no"`.
#' @param ... Passed to the underlying print helper.
#' @export
print.icc <- function(x,
                      digits = 4,
                      ci_digits = 4,
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

  is_ci_obj <- inherits(x, "icc_ci") ||
    (is.list(x) && all(c("est", "lwr.ci", "upr.ci") %in% names(x)))

  est <- if (is_ci_obj) as.matrix(x$est) else as.matrix(x)
  .mc_print_corr_matrix(
    x,
    header = "Intraclass correlation matrix",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    mat = est,
    ...
  )
  invisible(x)
}

#' @rdname icc
#' @method summary icc
#' @param object An intraclass-correlation object returned by `icc()`.
#' @export
summary.icc <- function(object,
                        digits = 4,
                        ci_digits = 2,
                        n = NULL,
                        topn = NULL,
                        max_vars = NULL,
                        width = NULL,
                        show_ci = NULL,
                        ...) {
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )

  is_ci_obj <- inherits(object, "icc_ci") ||
    (is.list(object) && all(c("est", "lwr.ci", "upr.ci") %in% names(object)))

  if (is_ci_obj) {
    est <- as.matrix(object$est)
    lwr <- as.matrix(object$lwr.ci)
    upr <- as.matrix(object$upr.ci)
    conf_level <- suppressWarnings(as.numeric(attr(object, "conf.level")))
  } else if (is.matrix(object)) {
    est <- as.matrix(object)
    lwr <- matrix(NA_real_, nrow(est), ncol(est), dimnames = dimnames(est))
    upr <- lwr
    conf_level <- suppressWarnings(as.numeric(attr(object, "conf.level")))
  } else {
    abort_bad_arg("object",
      message = "must be a matrix or a list with elements `est`, `lwr.ci`, and `upr.ci`."
    )
  }
  if (length(conf_level) != 1L || !is.finite(conf_level)) conf_level <- NA_real_

  rn <- rownames(est); cn <- colnames(est)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(est)))
  if (is.null(cn)) cn <- as.character(seq_len(ncol(est)))

  has_any_ci <- any(is.finite(lwr) | is.finite(upr))
  include_ci <- identical(show_ci, "yes") && has_any_ci

  rows <- vector("list", nrow(est) * (ncol(est) - 1L) / 2L)
  k <- 0L
  for (i in seq_len(nrow(est) - 1L)) {
    for (j in (i + 1L):ncol(est)) {
      k <- k + 1L
      rec <- list(
        method1 = rn[i],
        method2 = cn[j],
        estimate = round(est[i, j], digits)
      )
      if (include_ci) {
        rec$lwr <- if (is.na(lwr[i, j])) NA_real_ else round(lwr[i, j], ci_digits)
        rec$upr <- if (is.na(upr[i, j])) NA_real_ else round(upr[i, j], ci_digits)
      }
      rows[[k]] <- rec
    }
  }

  df <- do.call(rbind.data.frame, rows)
  rownames(df) <- NULL
  df$estimate <- as.numeric(df$estimate)
  if (include_ci) {
    df$lwr <- as.numeric(df$lwr)
    df$upr <- as.numeric(df$upr)
  }

  df <- structure(df, class = c("summary.icc", "data.frame"))
  attr(df, "overview") <- .mc_summary_corr_matrix(est, topn = topn)
  attr(df, "conf.level") <- if (is.finite(conf_level)) conf_level else NA_real_
  attr(df, "has_ci") <- isTRUE(include_ci)
  attr(df, "digits") <- digits
  attr(df, "ci_digits") <- ci_digits
  attr(df, "ci_method") <- attr(object, "ci.method", exact = TRUE) %||% NA_character_
  df
}

#' @rdname icc
#' @method print summary.icc
#' @export
print.summary.icc <- function(x, digits = NULL, n = NULL,
                              topn = NULL, max_vars = NULL,
                              width = NULL, show_ci = NULL, ...) {
  .mc_print_pairwise_summary_digest(
    x,
    title = "Intraclass correlation summary",
    digits = .mc_coalesce(digits, 4),
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ci_method = .mc_coalesce(attr(x, "ci_method"), "anova_f"),
    ...
  )
  invisible(x)
}

#' @rdname icc
#' @method print icc_overall
#' @param x An intraclass-correlation object returned by `icc()`, or a summary
#'   object returned by `summary()` for that fit.
#' @export
print.icc_overall <- function(x,
                              digits = 4,
                              ci_digits = 4,
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

  coeff <- summary(x,
    digits = digits,
    ci_digits = ci_digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci
  )

  diag_info <- attr(x, "diagnostics", exact = TRUE) %||% list()
  digest <- c(method = attr(x, "method", exact = TRUE) %||% "Overall intraclass correlation")
  if (is.finite(diag_info$n_subjects %||% NA_real_)) {
    digest <- c(digest, subjects = .mc_count_fmt(diag_info$n_subjects))
  }
  if (is.finite(diag_info$n_raters %||% NA_real_)) {
    digest <- c(digest, raters = .mc_count_fmt(diag_info$n_raters))
  }
  if (!is.null(attr(x, "selected_coefficient", exact = TRUE))) {
    digest <- c(digest, selected = attr(x, "selected_coefficient", exact = TRUE))
  }
  .mc_print_named_digest(digest, header = "Overall intraclass correlation")

  cfg <- .mc_resolve_display_args(
    context = "print",
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci
  )
  cat("\nCoefficient table\n\n")
  .mc_print_preview_table(
    if (identical(show_ci, "no")) coeff[, setdiff(names(coeff), c("lwr", "upr")), drop = FALSE] else coeff,
    n = cfg$n,
    topn = cfg$topn,
    max_vars = cfg$max_vars,
    width = cfg$width,
    context = "print",
    full_hint = TRUE,
    summary_hint = FALSE,
    ...
  )
  invisible(x)
}

#' @rdname icc
#' @method summary icc_overall
#' @param object An intraclass-correlation object returned by `icc()`.
#' @export
summary.icc_overall <- function(object,
                                digits = 4,
                                ci_digits = 2,
                                n = NULL,
                                topn = NULL,
                                max_vars = NULL,
                                width = NULL,
                                show_ci = NULL,
                                ...) {
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )

  coeff <- object$coefficients
  num_cols <- intersect(c("estimate", "statistic", "df1", "df2", "p_value"), names(coeff))
  for (nm in num_cols) coeff[[nm]] <- round(as.numeric(coeff[[nm]]), digits)
  if (all(c("lwr", "upr") %in% names(coeff))) {
    coeff$lwr <- round(as.numeric(coeff$lwr), ci_digits)
    coeff$upr <- round(as.numeric(coeff$upr), ci_digits)
  }
  coeff$selected <- as.logical(coeff$selected)

  anova <- object$anova
  if (!is.null(anova)) {
    anova$df <- round(as.numeric(anova$df), digits)
    anova$sum_sq <- round(as.numeric(anova$sum_sq), digits)
    anova$mean_sq <- round(as.numeric(anova$mean_sq), digits)
    anova$statistic <- round(as.numeric(anova$statistic), digits)
    anova$p_value <- round(as.numeric(anova$p_value), digits)
  }

  out <- structure(coeff, class = c("summary.icc_overall", "data.frame"))
  attr(out, "anova") <- anova
  attr(out, "mean_squares") <- object$mean_squares
  attr(out, "conf.level") <- attr(object, "conf.level", exact = TRUE)
  attr(out, "has_ci") <- all(c("lwr", "upr") %in% names(coeff))
  attr(out, "digits") <- digits
  attr(out, "ci_digits") <- ci_digits
  attr(out, "selected_coefficient") <- attr(object, "selected_coefficient", exact = TRUE)
  attr(out, "diagnostics") <- attr(object, "diagnostics", exact = TRUE)
  out
}

#' @rdname icc
#' @method print summary.icc_overall
#' @param x An intraclass-correlation object returned by `icc()`, or a summary
#'   object returned by `summary()` for that fit.
#' @export
print.summary.icc_overall <- function(x,
                                      digits = NULL,
                                      n = NULL,
                                      topn = NULL,
                                      max_vars = NULL,
                                      width = NULL,
                                      show_ci = NULL,
                                      ...) {
  show_ci <- .mc_resolve_show_ci(show_ci, context = "summary")
  digits <- .mc_coalesce(digits, attr(x, "digits", exact = TRUE) %||% 4)
  cfg <- .mc_resolve_display_args(
    context = "summary",
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci
  )
  diag_info <- attr(x, "diagnostics", exact = TRUE) %||% list()
  ms <- attr(x, "mean_squares", exact = TRUE)

  digest <- character()
  if (!is.null(attr(x, "selected_coefficient", exact = TRUE))) {
    digest <- c(digest, selected = attr(x, "selected_coefficient", exact = TRUE))
  }
  if (is.finite(diag_info$n_subjects %||% NA_real_)) {
    digest <- c(digest, subjects = .mc_count_fmt(diag_info$n_subjects))
  }
  if (is.finite(diag_info$n_raters %||% NA_real_)) {
    digest <- c(digest, raters = .mc_count_fmt(diag_info$n_raters))
  }
  if (is.finite(diag_info$n_complete %||% NA_real_)) {
    digest <- c(digest, n_complete = .mc_count_fmt(diag_info$n_complete))
  }
  if (is.numeric(ms) && length(ms)) {
    digest <- c(
      digest,
      ms_subject = formatC(unname(ms["ms_subject"]), format = "f", digits = digits),
      ms_rater = formatC(unname(ms["ms_rater"]), format = "f", digits = digits),
      ms_error = formatC(unname(ms["ms_error"]), format = "f", digits = digits)
    )
  }
  digest <- digest[!is.na(digest) & nzchar(digest)]
  .mc_print_named_digest(digest, header = "Overall intraclass correlation summary")

  cat("\nICC coefficients\n\n")
  .mc_print_preview_table(
    if (identical(show_ci, "no")) x[, setdiff(names(x), c("lwr", "upr")), drop = FALSE] else x,
    n = cfg$n,
    topn = cfg$topn,
    max_vars = cfg$max_vars,
    width = cfg$width,
    context = "summary",
    full_hint = FALSE,
    summary_hint = FALSE,
    ...
  )

  anova <- attr(x, "anova", exact = TRUE)
  if (is.data.frame(anova) && nrow(anova)) {
    cat("\nANOVA decomposition\n\n")
    .mc_print_preview_table(
      anova,
      n = cfg$n,
      topn = cfg$topn,
      max_vars = cfg$max_vars,
      width = cfg$width,
      context = "summary",
      full_hint = FALSE,
      summary_hint = FALSE,
      ...
    )
  }
  invisible(x)
}

#' @title Repeated-Measures Intraclass Correlation via REML
#'
#' @description
#' Computes pairwise repeated-measures intraclass correlation coefficients from
#' long-format data using the same REML and Woodbury-identity backend used by
#' the repeated-measures agreement models.
#'
#' @details
#' The repeated-measures model is fit separately for each method pair using the
#' same REML and Woodbury-identity backend used by the repeated-measures
#' concordance estimator.
#'
#' \strong{Kernel \eqn{D_m} and the fixed-bias term \eqn{S_B}.}
#' Let \eqn{d = L^\top \hat\beta} stack the within-time, pairwise method
#' differences, grouped by time. The symmetric positive semidefinite kernel
#' \eqn{D_m \succeq 0} selects which functional of the bias profile is targeted
#' by \eqn{S_B}. Internally, the code rescales any supplied or constructed
#' \eqn{D_m} to satisfy \eqn{1^\top D_m 1 = n_t} for stability and comparability.
#'
#' \itemize{
#'   \item \code{Dmat_type = "time-avg"} targets the square of the
#'   time-averaged bias.
#'   \item \code{Dmat_type = "typical-visit"} targets the average of squared
#'   per-time biases.
#'   \item \code{Dmat_type = "weighted-avg"} targets the square of a weighted
#'   time average.
#'   \item \code{Dmat_type = "weighted-sq"} targets the weighted average of
#'   squared per-time biases.
#' }
#'
#' As in the repeated-measures concordance implementation, \eqn{S_B} is the
#' fitted fixed-effect dispersion term induced by \eqn{D_m}. It enters the
#' denominator only for `type = "agreement"`.
#'
#' \strong{Time-averaging and shrinkage factors.}
#' The fitted variance components reported in the summary are
#' \eqn{\sigma_S^2}, \eqn{\sigma_{S\times M}^2}, \eqn{\sigma_{S\times T}^2},
#' and \eqn{\sigma_E^2}, stored respectively as `sigma2_subject`,
#' `sigma2_subject_method`, `sigma2_subject_time`, and `sigma2_error`.
#'
#' The repeated-measures ICC always uses \eqn{\sigma_S^2} in the numerator only.
#' The denominator uses the same time-averaging logic already used by the
#' repeated-measures concordance backend, through two pair-specific factors
#' \eqn{\bar{\kappa}_g} and \eqn{\bar{\kappa}_e}.
#'
#' If `time` is absent, the implementation sets
#' \deqn{ \bar{\kappa}_g = 0, \qquad \bar{\kappa}_e = 1. }
#'
#' If `time` is present and the target is a single visit
#' (\code{Dmat_type \%in\% c("typical-visit", "weighted-sq")}), the implementation
#' leaves the time-varying terms unshrunk:
#' \deqn{ \bar{\kappa}_g = 1, \qquad \bar{\kappa}_e = 1. }
#'
#' If `time` is present and the target is time-averaged
#' (\code{Dmat_type \%in\% c("time-avg", "weighted-avg")}), then for each observed
#' subject-method unit with \eqn{T} distinct observed visits:
#' \itemize{
#'   \item with equal weights,
#'   \deqn{ \kappa_g = \frac{1}{T}, \qquad
#'          \kappa_e = \frac{T + 2\sum_{h=1}^{T-1}(T-h)\rho^h}{T^2}, }
#'   with \eqn{\kappa_e=\kappa_g} when residuals are iid;
#'   \item with normalized visit weights \eqn{w_1,\ldots,w_T},
#'   \deqn{ \kappa_g = \sum_{t=1}^{T} w_t^2, \qquad
#'          \kappa_e = \sum_{t=1}^{T}\sum_{s=1}^{T} w_t w_s \rho^{|t-s|}, }
#'   with \eqn{\kappa_e=\kappa_g} when residuals are iid.
#' }
#'
#' With unbalanced \eqn{T}, the implementation averages the per-unit
#' \eqn{\kappa} values across the observations contributing to the pair and then
#' clamps both \eqn{\bar{\kappa}_g} and \eqn{\bar{\kappa}_e} to
#' \eqn{[10^{-12},\,1]} for numerical stability.
#'
#' \strong{Repeated-measures ICC.}
#' Let \eqn{\sigma_{S\times M,\mathrm{eff}}^2} denote the effective
#' subject-by-method variance term, equal to \eqn{\sigma_{S\times M}^2} when
#' the subject-by-method random effect is included and \eqn{0} otherwise. Let
#' \eqn{\sigma_{S\times T,\mathrm{eff}}^2} denote the effective subject-by-time
#' variance term, equal to \eqn{\sigma_{S\times T}^2} when the subject-by-time
#' random effect is included and \eqn{0} otherwise.
#'
#' For `type = "consistency"`, the reported ICC is
#' \deqn{ \mathrm{ICC}_{\mathrm{consistency}} \;=\;
#'       \frac{\sigma_S^2}{
#'       \sigma_S^2 + \sigma_{S\times M,\mathrm{eff}}^2 +
#'       \bar{\kappa}_g\,\sigma_{S\times T,\mathrm{eff}}^2 +
#'       \bar{\kappa}_e\,\sigma_E^2}. }
#'
#' For `type = "agreement"`, the denominator additionally includes \eqn{S_B}:
#' \deqn{ \mathrm{ICC}_{\mathrm{agreement}} \;=\;
#'       \frac{\sigma_S^2}{
#'       \sigma_S^2 + \sigma_{S\times M,\mathrm{eff}}^2 +
#'       \bar{\kappa}_g\,\sigma_{S\times T,\mathrm{eff}}^2 +
#'       \bar{\kappa}_e\,\sigma_E^2 + S_B}. }
#'
#' This differs from repeated-measures concordance because ICC uses only
#' \eqn{\sigma_S^2} in the numerator, whereas the concordance numerator also
#' includes the time-averaged subject-time term. Extra random-effect variances
#' \eqn{\{\sigma_{Z,j}^2\}} from `slope` / `slope_Z` are estimated by the shared
#' backend but are not included in the ICC denominator.
#'
#' \strong{CIs / SEs (delta method for ICC).}
#' Confidence intervals are built from the same REML fit by a large-sample delta
#' method. If `ci_mode = "raw"`, a Wald interval is formed on the ICC scale,
#' \deqn{ \widehat{\mathrm{ICC}} \;\pm\; z_{1-\alpha/2}\,
#'       \widehat{\mathrm{se}}\{\widehat{\mathrm{ICC}}\}, }
#' and truncated to \eqn{[0,1]}. If `ci_mode = "logit"`, the same Wald
#' construction is applied on the logit scale and then back-transformed. If
#' `ci_mode = "auto"`, the backend selects between the raw-scale and logit-scale
#' interval per estimate.
#'
#' \strong{Choosing \eqn{\rho} for AR(1).}
#' When \code{ar="ar1"} and \code{ar_rho = NA}, \eqn{\rho} is estimated by
#' profiling the REML log-likelihood at \eqn{(\hat\beta,\hat G,\hat\sigma_E^2)}.
#' With very few visits per subject, \eqn{\rho} can be weakly identified; consider
#' sensitivity checks over a plausible range.
#'
#' @section Notes on stability and performance:
#' All per-subject solves are \eqn{r\times r} with
#' \eqn{r = 1 + n_m + n_t + q_Z}, so cost scales with the number of subjects and
#' the fixed-effects dimension rather than the total number of observations.
#' Solvers use symmetric positive-definite paths with a small diagonal ridge and
#' pseudo-inverse fallback, which helps for very small or unbalanced subsets and
#' near-boundary estimates. For \code{AR(1)}, observations are ordered by time
#' within subject; \code{NA} time codes break the run, and gaps between factor
#' levels are treated as regular steps.
#'
#' Heteroscedastic slopes across \eqn{Z} columns are supported. Each \eqn{Z}
#' column has its own variance component \eqn{\sigma_{Z,j}^2}, but
#' cross-covariances among \eqn{Z} columns are set to zero.
#'
#' @section Threading and BLAS guards:
#' The C++ backend uses OpenMP loops while also forcing vendor BLAS libraries to
#' run single-threaded so that overall CPU usage stays predictable. On OpenBLAS
#' and Apple's Accelerate this is handled automatically. On Intel MKL builds the
#' guard is disabled by default, but you can also opt out manually by setting
#' \code{MATRIXCORR_DISABLE_BLAS_GUARD=1} in the environment before loading the
#' package.
#'
#' @param data A data frame.
#' @param response Character. Response variable name.
#' @param subject Character. Subject ID variable name.
#' @param method Character or \code{NULL}. Optional column name of method factor
#'   (added to fixed effects).
#' @param time Character or \code{NULL}. Optional column name of time factor
#'   (added to fixed effects).
#' @param type Character scalar; one of \code{c("consistency","agreement")}.
#'   For \code{"consistency"}, systematic fixed method differences are not
#'   penalized through \eqn{S_B}. For \code{"agreement"}, the fixed-effect
#'   dispersion term \eqn{S_B} is added to the denominator.
#' @param interaction Logical. Include \code{method:time} interaction?
#'   (default \code{FALSE}).
#' @param max_iter Integer. Maximum iterations for variance-component updates
#'   (default \code{100}).
#' @param tol Numeric. Convergence tolerance on parameter change
#'   (default \code{1e-6}).
#'
#' @param Dmat Optional \eqn{n_t \times n_t} numeric matrix to weight/aggregate
#'   time-specific fixed biases in the \eqn{S_B} quadratic form. If supplied, it
#'   is used (after optional mass rescaling; see \code{Dmat_rescale}) whenever at
#'   least two \emph{present} time levels exist; otherwise it is ignored. \strong{If
#'   \code{Dmat} is \code{NULL}}, a canonical kernel \eqn{D_m} is \emph{constructed}
#'   from \code{Dmat_type} and \code{Dmat_weights} (see below). \code{Dmat} should be
#'   symmetric positive semidefinite; small asymmetries are symmetrized internally.
#'
#' @param Dmat_type Character, one of \code{c("time-avg","typical-visit",
#'   "weighted-avg","weighted-sq")}. Only used when \code{Dmat = NULL}.
#'   It selects the aggregation target for time-specific fixed biases in
#'   \eqn{S_B}. Options are:
#'
#'   \itemize{
#'     \item \code{"time-avg"}: square of the time-averaged bias, \eqn{D_m=(1/n_t)\,11^\top}.
#'     \item \code{"typical-visit"}: average of squared per-time biases, \eqn{D_m=I_{n_t}}.
#'     \item \code{"weighted-avg"}: square of a weighted average, \eqn{D_m=n_t\,w\,w^\top} with \eqn{\sum w=1}.
#'     \item \code{"weighted-sq"}: weighted average of squared biases, \eqn{D_m=n_t\,\mathrm{diag}(w)} with \eqn{\sum w=1}.
#'   }
#'   Pick \code{"time-avg"} for ICC targeting the time-averaged measurement; pick
#'   \code{"typical-visit"} for ICC targeting a randomly sampled visit (typical
#'   occasion). Default \code{"time-avg"}.
#'
#' @param Dmat_weights Optional numeric weights \eqn{w} used when
#'   \code{Dmat_type \%in\% c("weighted-avg","weighted-sq")}. Must be nonnegative and
#'   finite. If \code{names(w)} are provided, they should match the \emph{full} time
#'   levels in \code{data}; they are aligned to the \emph{present} time subset per fit.
#'   If unnamed, the length must equal the number of present time levels. In all cases
#'   \eqn{w} is internally normalized to sum to 1.
#'
#' @param Dmat_rescale Logical. When \code{TRUE} (default), the supplied/built
#'   \eqn{D_m} is rescaled to satisfy the simple mass rule
#'   \eqn{1^\top D_m 1 = n_t}. This keeps the \eqn{S_B} denominator invariant and
#'   harmonizes with the time-averaging factors used for the variance terms.
#'
#' @param ci Logical. If \code{TRUE}, return a CI container for the repeated ICC.
#' @param conf_level Numeric in \eqn{(0,1)}. Confidence level when
#'   \code{ci = TRUE} (default \code{0.95}).
#' @param ci_mode Character scalar; one of \code{c("auto","raw","logit")}.
#'   Controls how confidence intervals are computed when \code{ci = TRUE}.
#'   If \code{"raw"}, a Wald CI is formed on the ICC scale and truncated to
#'   \code{[0,1]}. If \code{"logit"}, a Wald CI is computed on the
#'   \eqn{\mathrm{logit}(\mathrm{ICC})} scale and back-transformed to the
#'   original scale. If \code{"auto"} (default), the backend chooses per estimate
#'   between the raw-scale and logit-scale interval.
#' @param verbose Logical. If \code{TRUE}, prints a structured summary of the
#'   fitted variance components and \eqn{S_B} for each fit. Default \code{FALSE}.
#' @param digits Integer \eqn{(\ge 0)}. Number of decimal places to use in the
#'   printed summary when \code{verbose = TRUE}. Default \code{4}.
#' @param use_message Logical. When \code{verbose = TRUE}, choose the printing
#'   mechanism, where \code{TRUE} uses \code{message()} (respects \code{sink()},
#'   easily suppressible via \code{suppressMessages()}), whereas \code{FALSE}
#'   uses \code{cat()} to \code{stdout}. Default \code{TRUE}.
#'
#' @param ar Character. Residual correlation structure: \code{"none"} (iid) or
#'   \code{"ar1"} for subject-level AR(1) correlation within contiguous time
#'   runs. Default \code{c("none")}.
#' @param ar_rho Numeric in \eqn{(-0.999,\,0.999)} or \code{NA}.
#'   If \code{ar = "ar1"} and \code{ar_rho} is finite, it is treated as fixed.
#'   If \code{ar = "ar1"} and \code{ar_rho = NA}, \eqn{\rho} is estimated by
#'   profiling a 1-D objective (REML when available; an approximation otherwise).
#'   Default \code{NA_real_}.
#' @param slope Character. Optional extra random-effect design \eqn{Z}.
#'   With \code{"subject"} a single random slope is added (one column in \eqn{Z});
#'   with \code{"method"} one column per method level is added; with
#'   \code{"custom"} you provide \code{slope_Z} directly. Default
#'   \code{c("none","subject","method","custom")}.
#' @param slope_var For \code{slope \%in\% c("subject","method")}, a character
#'   string giving the name of a column in \code{data} used as the slope regressor
#'   (e.g., centered time). It is looked up inside \code{data}; do not
#'   pass the vector itself. NAs are treated as zeros in \eqn{Z}.
#' @param slope_Z For \code{slope = "custom"}, a numeric matrix with \eqn{n}
#'   rows (same order as \code{data}) providing the full extra random-effect
#'   design \eqn{Z}. \strong{Each column of \code{slope_Z} has its own variance
#'   component} \eqn{\sigma_{Z,j}^2}; columns are treated as \emph{uncorrelated}
#'   (diagonal block in \eqn{G}). Ignored otherwise.
#' @param drop_zero_cols Logical. When \code{slope = "method"}, drop all-zero
#'   columns of \eqn{Z} after subsetting (useful in pairwise fits). Default
#'   \code{TRUE}.
#'
#' @param vc_select Character scalar; one of \code{c("auto","none")}.
#'   Controls how the subject by method \eqn{\sigma^2_{A\times M}} ("subj_method") and
#'   subject by time \eqn{\sigma^2_{A\times T}} ("subj_time") variance components are
#'   included. If \code{"auto"} (default), the function performs boundary-aware
#'   REML likelihood-ratio tests (LRTs; null on the boundary at zero with a
#'   half-\eqn{\chi^2_1} reference) to decide whether to retain each component,
#'   in the order given by \code{vc_test_order}. If \code{"none"}, no testing
#'   is done and inclusion is taken from \code{include_subj_method}/\code{include_subj_time}
#'   (or, if \code{NULL}, from the mere presence of the corresponding factor in
#'   the design). In pairwise fits, the decision is made independently for each
#'   method pair.
#'
#' @param vc_alpha Numeric scalar in \eqn{(0,1)}; default \code{0.05}.
#'   Per-component significance level for the boundary-aware REML LRTs used when
#'   \code{vc_select = "auto"}. The tests are one-sided for variance components
#'   on the boundary and are not multiplicity-adjusted.
#'
#' @param vc_test_order Character vector (length 2) with a permutation of
#'   \code{c("subj_time","subj_method")}; default \code{c("subj_time","subj_method")}. Specifies the order
#'   in which the two variance components are tested when \code{vc_select = "auto"}.
#'   The component tested first may be dropped before testing the second. If a
#'   factor is absent in the design (e.g., no time factor so "subj_time" is undefined),
#'   the corresponding test is skipped.
#'
#' @param include_subj_method,include_subj_time Logical scalars or \code{NULL}.
#'   When \code{vc_select = "none"}, these control whether the
#'   \eqn{\sigma^2_{A\times M}} ("subj_method") and \eqn{\sigma^2_{A\times T}} ("subj_time")
#'   random effects are included (\code{TRUE}) or excluded (\code{FALSE}) in the model.
#'   If \code{NULL} (default), inclusion defaults to the presence of the
#'   corresponding factor in the data (i.e., at least two method/time levels).
#'   When \code{vc_select = "auto"}, these arguments are ignored (automatic
#'   selection is used instead).
#'
#' @param sb_zero_tol Non-negative numeric scalar; default \code{1e-10}.
#'   Numerical threshold for the fixed-effect dispersion term \eqn{S_B}.
#'   After computing \eqn{\widehat{S_B}} and its delta-method variance, if
#'   \eqn{\widehat{S_B} \le} \code{sb_zero_tol} or non-finite, the procedure
#'   treats \eqn{S_B} as fixed at zero in the delta step.
#' @param x For `print()`, an object returned by `icc_rm_reml()` or
#'   `summary.icc_rm_reml()`.
#' @param object For `summary()`, an object returned by `icc_rm_reml()`.
#' @param ci_digits Integer; number of digits for confidence interval bounds in
#'   printed method summaries.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives
#'   this from console width.
#' @param width Optional display width; defaults to `getOption("width")`.
#' @param show_ci One of `"yes"` or `"no"`.
#' @param ... Passed to the underlying display helpers.
#'
#' @return
#' A repeated-measures pairwise ICC object. Without confidence intervals the
#' result is a symmetric matrix of class `c("icc_rm_reml", "icc", "matrix")`.
#' With confidence intervals it is a list with `est`, `lwr.ci`, and `upr.ci`
#' and class `c("icc_rm_reml", "icc_ci", "icc")`. Both carry the fitted
#' variance-component matrices as attributes.
#'
#' @seealso `icc()`, [ccc()], [ccc_rm_reml()], [ba()], [ba_rm()], [rmcorr()]
#'
#' @references
#' Shrout PE, Fleiss JL (1979). Intraclass correlations: uses in assessing
#' rater reliability. Psychological Bulletin, 86(2), 420-428.
#'
#' McGraw KO, Wong SP (1996). Forming inferences about some intraclass
#' correlation coefficients. Psychological Methods, 1(1), 30-46.
#'
#' @examples
#' set.seed(321)
#' n_id <- 20
#' n_time <- 3
#' id <- factor(rep(seq_len(n_id), each = 2 * n_time))
#' method <- factor(rep(rep(c("A", "B"), each = n_time), times = n_id))
#' time <- factor(rep(seq_len(n_time), times = 2 * n_id))
#'
#' subj <- rnorm(n_id, sd = 1)[as.integer(id)]
#' subj_method <- rnorm(n_id * 2, sd = 0.25)
#' sm <- subj_method[(as.integer(id) - 1L) * 2L + as.integer(method)]
#' y <- subj + sm + 0.3 * (method == "B") + rnorm(length(id), sd = 0.35)
#'
#' dat_rm <- data.frame(y = y, id = id, method = method, time = time)
#'
#' fit_icc_rm <- icc_rm_reml(
#'   dat_rm,
#'   response = "y",
#'   subject = "id",
#'   method = "method",
#'   time = "time",
#'   type = "consistency"
#' )
#' print(fit_icc_rm)
#' summary(fit_icc_rm)
#'
#' @export
icc_rm_reml <- function(data, response, subject,
                        method = NULL, time = NULL,
                        type = c("consistency", "agreement"),
                        interaction = FALSE,
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
  type <- match.arg(type)
  ar <- match.arg(ar)
  slope <- match.arg(slope)
  Dmat_type <- match.arg(Dmat_type)
  vc_select <- match.arg(vc_select)
  vc_test_order <- match.arg(vc_test_order, several.ok = TRUE)
  ci_mode <- match.arg(ci_mode)
  ci_mode_int <- switch(ci_mode, raw = 0L, logit = 1L, auto = 2L)

  check_bool(interaction, arg = "interaction")
  max_iter <- check_scalar_int_pos(max_iter, arg = "max_iter")
  check_scalar_nonneg(tol, arg = "tol", strict = TRUE)
  check_bool(ci, arg = "ci")
  check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  check_bool(verbose, arg = "verbose")
  check_bool(use_message, arg = "use_message")
  check_bool(Dmat_rescale, arg = "Dmat_rescale")
  check_bool(drop_zero_cols, arg = "drop_zero_cols")
  check_prob_scalar(vc_alpha, arg = "vc_alpha", open_ends = TRUE)

  if (identical(ar, "ar1")) {
    if (length(ar_rho) != 1L) {
      abort_bad_arg("ar_rho",
        message = "must be length 1 (or NA to estimate)."
      )
    }
    if (!is.na(ar_rho)) {
      check_scalar_numeric(ar_rho,
        arg = "ar_rho",
        lower = -0.999,
        upper = 0.999,
        closed_lower = FALSE,
        closed_upper = FALSE
      )
    }
  } else {
    ar_rho <- NA_real_
  }

  df <- as.data.frame(data)
  req_cols <- c(response, subject, method, time, slope_var)
  req_cols <- req_cols[!vapply(req_cols, is.null, logical(1))]
  check_required_cols(df, req_cols, df_arg = "data")

  df[[response]] <- as.numeric(df[[response]])
  if (anyNA(df[[response]])) {
    abort_bad_arg("response",
      message = "must reference a numeric column in {.arg data}."
    )
  }
  df[[subject]] <- factor(df[[subject]])
  if (!is.null(method)) df[[method]] <- factor(df[[method]])
  if (!is.null(time)) df[[time]] <- factor(df[[time]])
  all_time_lvls <- if (!is.null(time)) levels(df[[time]]) else character(0)

  terms_rhs <- "1"
  if (!is.null(method)) terms_rhs <- c(terms_rhs, method)
  if (!is.null(time)) terms_rhs <- c(terms_rhs, time)
  if (!is.null(method) && !is.null(time) && interaction) {
    terms_rhs <- c(terms_rhs, sprintf("%s:%s", method, time))
  }
  fml <- as.formula(paste("~", paste(terms_rhs, collapse = " + ")))

  extra_label <- switch(slope,
    subject = "random slope (subject)",
    method = "random slope (by method)",
    custom = "custom random effect",
    NULL
  )

  if (is.null(method) || nlevels(df[[method]]) < 2L) {
    cli::cli_abort(
      c(
        "At least two method levels are required to compute pairwise intraclass correlation.",
        "i" = "Supply {.arg method} with >= 2 levels; the overall coefficient is not computed."
      )
    )
  }

  metric_mode <- if (identical(type, "agreement")) 2L else 1L

  out <- ccc_lmm_reml_pairwise(
    df = df,
    fml = fml,
    response = response,
    rind = subject,
    method = method,
    time = time,
    slope = slope,
    slope_var = slope_var,
    slope_Z = slope_Z,
    drop_zero_cols = drop_zero_cols,
    Dmat = Dmat,
    ar = ar,
    ar_rho = ar_rho,
    max_iter = max_iter,
    tol = tol,
    conf_level = conf_level,
    verbose = verbose,
    digits = digits,
    use_message = use_message,
    extra_label = extra_label,
    ci = ci,
    ci_mode_int = ci_mode_int,
    all_time_lvls = all_time_lvls,
    Dmat_type = Dmat_type,
    Dmat_weights = Dmat_weights,
    Dmat_rescale = Dmat_rescale,
    vc_select = vc_select,
    vc_alpha = vc_alpha,
    vc_test_order = vc_test_order,
    include_subj_method = include_subj_method,
    include_subj_time = include_subj_time,
    sb_zero_tol = sb_zero_tol,
    metric_mode = metric_mode,
    estimate_key = "metric",
    se_key = "se_metric",
    out_class = "icc_rm_reml",
    out_method = paste("Repeated-measures intraclass correlation", sprintf("(%s)", type)),
    out_description = "Pairwise repeated-measures intraclass correlation from random-effects REML",
    summary_title = "Repeated-measures intraclass correlation matrix",
    vc_engine_name = "icc_rm_reml",
    vc_se_label = "SE(ICC)"
  )
  attr(out, "type") <- type
  attr(out, "model") <- "reml_variance_components"
  out
}

#' @rdname icc_rm_reml
#' @method print icc_rm_reml
#' @param x For `print()`, an object returned by `icc_rm_reml()` or
#'   `summary.icc_rm_reml()`.
#' @param ci_digits Integer; number of digits for confidence interval bounds in
#'   printed method summaries.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives
#'   this from console width.
#' @param width Optional display width; defaults to `getOption("width")`.
#' @param show_ci One of `"yes"` or `"no"`.
#' @param ... Passed to the underlying display helpers.
#' @export
print.icc_rm_reml <- function(x,
                              digits = 4,
                              ci_digits = 4,
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
  est <- if (is.list(x) && !is.null(x$est)) as.matrix(x$est) else as.matrix(x)
  .mc_print_corr_matrix(
    x,
    header = "Repeated-measures intraclass correlation matrix",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    mat = est,
    ...
  )
  invisible(x)
}

#' @rdname icc_rm_reml
#' @method summary icc_rm_reml
#' @param object For `summary()`, an object returned by `icc_rm_reml()`.
#' @param ci_digits Integer; number of digits for confidence interval bounds in
#'   printed method summaries.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives
#'   this from console width.
#' @param width Optional display width; defaults to `getOption("width")`.
#' @param show_ci One of `"yes"` or `"no"`.
#' @param ... Passed to the underlying display helpers.
#' @export
summary.icc_rm_reml <- function(object,
                                digits = 4,
                                ci_digits = 2,
                                n = NULL,
                                topn = NULL,
                                max_vars = NULL,
                                width = NULL,
                                show_ci = NULL,
                                ...) {
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )

  base_summary <- summary.icc(
    object,
    digits = digits,
    ci_digits = ci_digits,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci
  )

  est_mat <- if (is.list(object) && !is.null(object$est)) as.matrix(object$est) else as.matrix(object)
  n_pairs <- nrow(base_summary)

  extract_pairs_num <- function(val) {
    out <- numeric(n_pairs)
    k <- 0L
    for (i in seq_len(nrow(est_mat) - 1L)) {
      for (j in (i + 1L):ncol(est_mat)) {
        k <- k + 1L
        out[k] <- if (is.null(val)) {
          NA_real_
        } else if (is.matrix(val)) {
          suppressWarnings(as.numeric(val[i, j]))
        } else {
          vv <- suppressWarnings(as.numeric(val))
          if (length(vv) == 1L) vv else NA_real_
        }
      }
    }
    out
  }

  extract_pairs_extra_list <- function(val) {
    if (is.matrix(val) && typeof(val) == "list") {
      out <- vector("list", n_pairs)
      k <- 0L
      for (i in seq_len(nrow(est_mat) - 1L)) {
        for (j in (i + 1L):ncol(est_mat)) {
          k <- k + 1L
          out[[k]] <- val[[i, j]]
        }
      }
      return(out)
    }
    if (is.null(val)) return(vector("list", n_pairs))
    vv <- suppressWarnings(as.numeric(val))
    replicate(n_pairs, vv, simplify = FALSE)
  }

  out <- base_summary
  out$sigma2_subject <- round(extract_pairs_num(attr(object, "sigma2_subject")), digits)
  out$sigma2_subject_method <- round(extract_pairs_num(attr(object, "sigma2_subject_method")), digits)
  out$sigma2_subject_time <- round(extract_pairs_num(attr(object, "sigma2_subject_time")), digits)
  out$sigma2_error <- round(extract_pairs_num(attr(object, "sigma2_error")), digits)
  out$SB <- round(extract_pairs_num(attr(object, "SB")), digits)
  out$se_icc <- round(extract_pairs_num(attr(object, "se_metric")), digits)
  out$residual_model <- rep(attr(object, "residual_model", exact = TRUE) %||% "iid", n_pairs)

  extra_list <- extract_pairs_extra_list(attr(object, "sigma2_extra"))
  max_k <- 0L
  for (v in extra_list) if (!is.null(v)) max_k <- max(max_k, length(v))
  if (max_k > 0L) {
    for (kk in seq_len(max_k)) {
      out[[paste0("sigma2_extra", kk)]] <- vapply(extra_list, function(v) {
        if (is.null(v) || length(v) < kk) return(NA_real_)
        round(as.numeric(v[[kk]]), digits)
      }, numeric(1))
    }
  }

  extract_pairs_logi <- function(val) {
    out <- logical(n_pairs)
    k <- 0L
    for (i in seq_len(nrow(est_mat) - 1L)) {
      for (j in (i + 1L):ncol(est_mat)) {
        k <- k + 1L
        out[k] <- if (is.null(val)) FALSE else isTRUE(if (is.matrix(val)) val[i, j] else val)
      }
    }
    out
  }

  if (!is.null(attr(object, "ar_rho"))) {
    out$ar1_rho <- round(extract_pairs_num(attr(object, "ar_rho")), digits)
  }
  if (!is.null(attr(object, "ar1_rho_lag1"))) {
    vals <- round(extract_pairs_num(attr(object, "ar1_rho_lag1")), digits)
    out$ar1_rho_lag1 <- vals
    out$ar1_rho_mom <- vals
  }
  if (!is.null(attr(object, "ar1_pairs"))) {
    out$ar1_pairs <- as.integer(extract_pairs_num(attr(object, "ar1_pairs")))
  }
  if (!is.null(attr(object, "ar1_pval"))) {
    out$ar1_pval <- round(extract_pairs_num(attr(object, "ar1_pval")), digits)
  }
  if (!is.null(attr(object, "use_ar1"))) {
    out$use_ar1 <- extract_pairs_logi(attr(object, "use_ar1"))
    out$ar1_recommend <- out$use_ar1
  }

  attr(out, "conf.level") <- attr(base_summary, "conf.level")
  attr(out, "has_ci") <- attr(base_summary, "has_ci")
  attr(out, "digits") <- digits
  attr(out, "ci_digits") <- ci_digits
  structure(out, class = c("summary.icc_rm_reml", "data.frame"))
}

#' @rdname icc_rm_reml
#' @method print summary.icc_rm_reml
#' @param x For `print()`, an object returned by `icc_rm_reml()` or
#'   `summary.icc_rm_reml()`.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives
#'   this from console width.
#' @param width Optional display width; defaults to `getOption("width")`.
#' @param show_ci One of `"yes"` or `"no"`.
#' @param ... Passed to the underlying display helpers.
#' @export
print.summary.icc_rm_reml <- function(x, digits = NULL, n = NULL,
                                      topn = NULL, max_vars = NULL,
                                      width = NULL, show_ci = NULL, ...) {
  show_ci <- .mc_resolve_show_ci(show_ci, context = "summary")
  has_ci <- isTRUE(attr(x, "has_ci")) || all(c("lwr", "upr") %in% names(x))
  cl <- suppressWarnings(as.numeric(attr(x, "conf.level")))
  if (!is.finite(cl)) cl <- NA_real_

  .mc_print_sectioned_table(
    x,
    sections = list(
      list(
        title = "ICC estimates",
        cols = c("method1", "method2", "estimate", "lwr", "upr", "se_icc", "residual_model")
      ),
      list(
        title = "Variance components",
        cols = c("sigma2_subject", "sigma2_subject_method", "sigma2_subject_time",
                 "sigma2_error", "SB", grep("^sigma2_extra", names(x), value = TRUE))
      ),
      list(
        title = "AR(1) diagnostics",
        cols = c("ar1_rho", "ar1_rho_lag1", "ar1_rho_mom", "ar1_pairs",
                 "ar1_pval", "use_ar1", "ar1_recommend")
      )
    ),
    header = .mc_header_with_ci("Repeated-measures intraclass correlation (REML)", cl, if (has_ci) show_ci else "no"),
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
  invisible(x)
}
