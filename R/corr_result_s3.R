#' S3 Summary for Dense Correlation Results
#'
#' Representation-first summary for dense correlation outputs.
#'
#' @param object A dense correlation result (`corr_matrix`).
#' @param topn Optional number of head/tail rows when preview is truncated.
#' @param show_ci One of `"yes"` or `"no"`.
#' @param ... Unused.
#'
#' @details
#' Arithmetic, comparisons, logical operators, `Math` functions such as
#' `abs()` and `round()`, and indexed replacement operate on the estimates
#' and return the usual base R result. Stored inference applies to the
#' original estimates. Use `estimate(object)` for numeric calculations.
#'
#' Pearson summaries without a stored inference payload report
#' `fisher_z = atanh(estimate)`, with non-finite transforms stored as `NA`.
#' For each pair whose finite complete-case count is greater than 3, they also
#' report the Fisher normal approximation for independent bivariate normal
#' samples, `statistic = fisher_z / (1 / sqrt(n_complete - 3))`, and a
#' two-sided normal p-value. That approximation is not an exact *t* test.
#' Pairs with `n_complete <= 3` leave the statistic and p-value missing, and
#' one such pair does not suppress inference for other eligible pairs.
#' A confidence-interval request does not by itself supply p-values.
#'
#' For every other estimator, the summary retrieves statistic, degrees of
#' freedom, and p-values only from inference the estimator already stored. It
#' does not calculate a test from the coefficient. The Fisher transform is not
#' a universal inference column. Biweight midcorrelation keeps its stored
#' Student-*t*-like statistic and exposes the stored standardized transform as
#' `z_statistic`.
#'
#' When the object already stores the matrices, the summary appends
#' `p_value_adjusted`, `reject`, `skipped_n`, and `skipped_prop` after the
#' existing columns. Those values are copied. They are not recomputed.
#'
#' @return A standardized summary data frame with class
#'   `c("summary.corr_result", "data.frame")` (plus compatibility classes).
#' @method summary corr_matrix
#' @export
summary.corr_matrix <- function(object, topn = NULL, show_ci = NULL, ...) {
  .mc_summary_corr_result(
    object,
    output_class = "summary.corr_matrix",
    topn = topn,
    show_ci = show_ci
  )
}

#' @export
Ops.corr_matrix <- function(e1, e2) {
  op <- .Primitive(.Generic)
  if (inherits(e1, "corr_matrix")) e1 <- .mc_plain_numeric_matrix(e1)
  if (missing(e2)) return(op(e1))
  if (inherits(e2, "corr_matrix")) e2 <- .mc_plain_numeric_matrix(e2)
  op(e1, e2)
}

#' @export
Math.corr_matrix <- function(x, ...) {
  op <- get(.Generic, envir = baseenv(), inherits = FALSE)
  op(.mc_plain_numeric_matrix(x), ...)
}

#' @export
`[<-.corr_matrix` <- function(x, ..., value) {
  `[<-`(.mc_plain_numeric_matrix(x), ..., value = value)
}

#' @export
`[[<-.corr_matrix` <- function(x, ..., value) {
  `[[<-`(.mc_plain_numeric_matrix(x), ..., value = value)
}

#' S3 Summary for Sparse Correlation Results
#'
#' Representation-first summary for sparse correlation outputs.
#'
#' @param object A sparse correlation result.
#' @param topn Optional number of head/tail rows when preview is truncated.
#' @param show_ci One of `"yes"` or `"no"`.
#' @param ... Unused.
#'
#' @details
#' Use `base::summary(object)` for the statistical summary regardless of
#' whether Matrix is attached. `Matrix::summary(object)` reports the stored
#' sparse entries as row indices, column indices and values.
#'
#' See [summary.corr_matrix()] for which inference columns are retrieved from
#' stored results, which Pearson-only Fisher normal approximation is
#' calculated when no inference payload was stored, and which stored
#' adjusted-p, rejection, and skipped-correlation diagnostic columns are
#' appended when those matrices are already present.
#'
#' @return A standardized summary data frame with class
#'   `c("summary.corr_result", "data.frame")` (plus compatibility classes).
#' @method summary corr_sparse
#' @export
summary.corr_sparse <- function(object, topn = NULL, show_ci = NULL, ...) {
  .mc_summary_corr_result(
    object,
    output_class = "summary.corr_sparse",
    topn = topn,
    show_ci = show_ci
  )
}

#' S3 Summary for Edge-List Correlation Results
#'
#' Representation-first summary for edge-list outputs.
#'
#' @param object An edge-list correlation result.
#' @param topn Optional number of head/tail rows when preview is truncated.
#' @param show_ci One of `"yes"` or `"no"`.
#' @param ... Unused.
#'
#' @details
#' See [summary.corr_matrix()] for which inference columns are retrieved from
#' stored results, which Pearson-only Fisher normal approximation is
#' calculated when no inference payload was stored, and which stored
#' adjusted-p, rejection, and skipped-correlation diagnostic columns are
#' appended when those matrices are already present.
#'
#' @return A standardized summary data frame with class
#'   `c("summary.corr_result", "data.frame")` (plus compatibility classes).
#' @method summary corr_edge_list
#' @export
summary.corr_edge_list <- function(object, topn = NULL, show_ci = NULL, ...) {
  .mc_summary_corr_result(
    object,
    output_class = "summary.corr_edge_list",
    topn = topn,
    show_ci = show_ci
  )
}

#' Edge-list data-frame view
#' @param x A `corr_edge_list` object.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param ... Unused.
#' @return A data frame with columns `row`, `col`, `value`.
#' @method as.data.frame corr_edge_list
#' @export
as.data.frame.corr_edge_list <- function(x, row.names = NULL, optional = FALSE, ...) {
  .mc_corr_as_edge_df(x)
}

#' Matrix sparse fallback for correlation sparse output
#' @keywords internal
#' @exportS3Method summary sparseMatrix
#' @noRd
summary.sparseMatrix <- function(object, topn = NULL, show_ci = NULL, ...) {
  if (identical(attr(object, "corr_output_class", exact = TRUE), "corr_sparse") ||
      identical(attr(object, "output", exact = TRUE), "sparse")) {
    return(summary.corr_sparse(object, topn = topn, show_ci = show_ci, ...))
  }
  NextMethod()
}

#' @keywords internal
#' @exportS3Method summary dsCMatrix
#' @noRd
summary.dsCMatrix <- function(object, topn = NULL, show_ci = NULL, ...) {
  if (identical(attr(object, "corr_output_class", exact = TRUE), "corr_sparse") ||
      identical(attr(object, "output", exact = TRUE), "sparse")) {
    return(summary.corr_sparse(object, topn = topn, show_ci = show_ci, ...))
  }
  NextMethod()
}

#' @keywords internal
#' @exportS3Method summary dgCMatrix
#' @noRd
summary.dgCMatrix <- function(object, topn = NULL, show_ci = NULL, ...) {
  if (identical(attr(object, "corr_output_class", exact = TRUE), "corr_sparse") ||
      identical(attr(object, "output", exact = TRUE), "sparse")) {
    return(summary.corr_sparse(object, topn = topn, show_ci = show_ci, ...))
  }
  NextMethod()
}

#' @title Print Edge-List Correlation Results
#' @param x An edge-list correlation result.
#' @param digits Number of digits for numeric values.
#' @param n Optional preview row threshold.
#' @param topn Optional number of head/tail rows when preview is truncated.
#' @param max_vars Optional maximum number of visible columns in preview.
#' @param width Optional output width.
#' @param show_ci One of `"yes"` or `"no"`.
#' @param ... Unused.
#' @method print corr_edge_list
#' @export
print.corr_edge_list <- function(x,
                                 digits = 4,
                                 n = NULL,
                                 topn = NULL,
                                 max_vars = NULL,
                                 width = NULL,
                                 show_ci = NULL,
                                 ...) {
  .mc_print_corr_records(
    x,
    output_label = "edge list",
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

#' @title Print Standardized Correlation Summaries
#' @param x A `summary.corr_result` object.
#' @param digits Number of digits for numeric values.
#' @param n Optional preview row threshold.
#' @param topn Optional number of head/tail rows when preview is truncated.
#' @param max_vars Optional maximum number of visible columns in preview.
#' @param width Optional output width.
#' @param show_ci One of `"yes"` or `"no"`.
#' @param ... Unused.
#' @return Invisibly returns `x`.
#' @method print summary.corr_result
#' @export
print.summary.corr_result <- function(x,
                                      digits = 4,
                                      n = NULL,
                                      topn = NULL,
                                      max_vars = NULL,
                                      width = NULL,
                                      show_ci = NULL,
                                      ...) {
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

  overview <- attr(x, "overview", exact = TRUE)
  header <- attr(x, "summary_title", exact = TRUE) %||%
    sprintf("%s summary", attr(x, "method", exact = TRUE) %||%
      (overview$method %||% "correlation"))
  digest <- c(
    output = attr(x, "output", exact = TRUE) %||% (overview$output %||% NA_character_),
    dimensions = if (is.list(overview)) sprintf("%d x %d", overview$n_rows, overview$n_cols) else NA_character_,
    retained_pairs = if (is.list(overview)) .mc_count_fmt(overview$n_pairs_retained) else NA_character_,
    threshold = if (is.list(overview)) formatC(overview$threshold, format = "f", digits = digits) else NA_character_,
    diag = if (is.list(overview)) if (isTRUE(overview$diag)) "included" else "excluded" else NA_character_,
    estimate = if (is.list(overview)) {
      .mc_format_scalar_or_range(overview$estimate_min, overview$estimate_max, digits = digits)
    } else {
      NA_character_
    }
  )
  if (is.list(overview) && is.finite(overview$threshold_sets) && overview$threshold_sets > 0L) {
    digest <- c(digest, thresholds = sprintf("%s set(s)", .mc_count_fmt(overview$threshold_sets)))
  }
  if (is.list(overview)) {
    skipped_n_txt <- .mc_format_scalar_or_range(overview$skipped_n_min, overview$skipped_n_max, digits = digits, integer = TRUE)
    if (!is.null(skipped_n_txt)) {
      digest <- c(digest, skipped_n = skipped_n_txt)
    }
    skipped_prop_txt <- .mc_format_scalar_or_range(overview$skipped_prop_min, overview$skipped_prop_max, digits = digits, integer = FALSE)
    if (!is.null(skipped_prop_txt)) {
      digest <- c(digest, skipped_prop = skipped_prop_txt)
    }
    if (!is.null(overview$inference_method) &&
        length(overview$inference_method) == 1L &&
        !is.na(overview$inference_method) &&
        nzchar(overview$inference_method)) {
      digest <- c(digest, inference = overview$inference_method)
    }
  }
  if (identical(cfg$show_ci, "yes")) {
    summary_df <- as.data.frame(x, stringsAsFactors = FALSE)
    digest <- c(
      digest,
      .mc_ci_digest(
        summary_df,
        conf_level = suppressWarnings(as.numeric(attr(x, "conf.level", exact = TRUE))),
        ci_method = attr(x, "ci_method", exact = TRUE) %||%
          (attr(x, "overview", exact = TRUE)$ci_method %||% NULL),
        digits = 3
      )
    )
    constant_m <- .mc_single_display_value(summary_df$m)
    if (!is.null(constant_m)) {
      digest <- c(digest, m = constant_m)
    }
    constant_boot <- .mc_single_display_value(summary_df$bootstrap_reps)
    if (!is.null(constant_boot)) {
      digest <- c(digest, boot_reps = constant_boot)
    }
    se_vals <- suppressWarnings(as.numeric(summary_df$se))
    se_vals <- se_vals[is.finite(se_vals)]
    if (length(se_vals)) {
      se_txt <- .mc_format_scalar_or_range(
        min(se_vals),
        max(se_vals),
        digits = 3,
        integer = FALSE
      )
      if (!is.null(se_txt)) {
        digest <- c(digest, se = se_txt)
      }
    }
  }
  digest <- digest[!is.na(digest) & nzchar(digest)]
  .mc_print_named_digest(digest, header = header)
  cat("\n")

  preview <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  ci_digits <- .mc_coalesce(attr(x, "ci_digits", exact = TRUE), digits)
  preview <- .mc_compact_corr_summary_table(
    preview,
    show_ci = cfg$show_ci,
    conf_level = suppressWarnings(as.numeric(attr(x, "conf.level", exact = TRUE))),
    digits = digits,
    ci_digits = ci_digits
  )
  preview <- .mc_format_summary_numeric_table(
    preview,
    digits = digits,
    ci_digits = ci_digits
  )
  .mc_print_preview_table(
    preview,
    n = cfg$n,
    topn = cfg$topn,
    max_vars = cfg$max_vars,
    width = cfg$width,
    context = "summary",
    full_hint = TRUE,
    summary_hint = FALSE
  )
  invisible(x)
}

.mc_single_display_value <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  vals <- unique(x[!is.na(x)])
  if (length(vals) != 1L) {
    return(NULL)
  }
  if (is.numeric(vals)) {
    vals <- vals[is.finite(vals)]
    if (length(vals) != 1L) {
      return(NULL)
    }
    return(format(vals[[1L]], trim = TRUE, scientific = FALSE))
  }
  vals <- as.character(vals)
  if (!nzchar(vals[[1L]])) {
    return(NULL)
  }
  vals[[1L]]
}

.mc_compact_corr_summary_table <- function(df,
                                           show_ci = "yes",
                                           conf_level = NA_real_,
                                           digits = 4,
                                           ci_digits = digits) {
  if (!is.data.frame(df) || !nrow(df)) {
    return(df)
  }

  out <- df
  if ("n_complete" %in% names(out)) {
    names(out)[match("n_complete", names(out))] <- "n"
  }

  if (identical(show_ci, "yes") && all(c("lwr", "upr") %in% names(out))) {
    ci_label <- if (is.finite(conf_level)) {
      sprintf("%g%% CI", 100 * conf_level)
    } else {
      "CI"
    }
    lwr <- suppressWarnings(as.numeric(out$lwr))
    upr <- suppressWarnings(as.numeric(out$upr))
    out[[ci_label]] <- ifelse(
      is.finite(lwr) & is.finite(upr),
      sprintf(
        "[%s, %s]",
        formatC(lwr, format = "f", digits = ci_digits),
        formatC(upr, format = "f", digits = ci_digits)
      ),
      NA_character_
    )
  }

  always_drop <- c(
    "lwr", "upr", "fisher_z", "z_statistic", "statistic", "df", "se",
    "bootstrap_reps", "skipped_n", "skipped_prop"
  )
  if (identical(show_ci, "no")) {
    always_drop <- c(always_drop, grep("CI$", names(out), value = TRUE))
  }

  if ("ci_method" %in% names(out) && !is.null(.mc_single_display_value(out$ci_method))) {
    always_drop <- c(always_drop, "ci_method")
  }
  if ("m" %in% names(out) && !is.null(.mc_single_display_value(out$m))) {
    always_drop <- c(always_drop, "m")
  }
  out <- out[, setdiff(names(out), always_drop), drop = FALSE]

  preferred <- c(
    "item1", "item2", "estimate", "n",
    grep("CI$", names(out), value = TRUE),
    "p_value", "p_value_adjusted", "ci_method", "m"
  )
  preferred <- unique(preferred[preferred %in% names(out)])
  out[, c(preferred, setdiff(names(out), preferred)), drop = FALSE]
}

#' S3 Plot for Dense Correlation Results
#'
#' @param x A dense correlation result (`corr_matrix`).
#' @param title Optional plot title.
#' @param low_color Fill color for -1.
#' @param high_color Fill color for +1.
#' @param mid_color Fill color for 0.
#' @param value_text_size Text size for optional overlaid values.
#' @param ci_text_size Text size for optional confidence-interval labels.
#' @param show_value Logical; overlay values if `TRUE`.
#' @param ... Additional theme arguments.
#'
#' @return A `ggplot` heatmap.
#' @method plot corr_matrix
#' @export
plot.corr_matrix <- function(x,
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
    title = title %||% "Correlation heatmap",
    low_color = low_color,
    high_color = high_color,
    mid_color = mid_color,
    value_text_size = value_text_size,
    ci_text_size = ci_text_size,
    show_value = show_value,
    ...
  )
}

#' S3 Plot for Sparse Correlation Results
#'
#' @inheritParams plot.corr_matrix
#' @param x A sparse correlation result.
#' @method plot corr_sparse
#' @export
plot.corr_sparse <- function(x,
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
    title = title %||% "Retained sparse correlation heatmap",
    low_color = low_color,
    high_color = high_color,
    mid_color = mid_color,
    value_text_size = value_text_size,
    ci_text_size = ci_text_size,
    show_value = show_value,
    ...
  )
}

#' S3 Plot for Edge-List Correlation Results
#'
#' @inheritParams plot.corr_matrix
#' @param x An edge-list correlation result.
#' @method plot corr_edge_list
#' @export
plot.corr_edge_list <- function(x,
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
    title = title %||% "Edge-list correlation heatmap",
    low_color = low_color,
    high_color = high_color,
    mid_color = mid_color,
    value_text_size = value_text_size,
    ci_text_size = ci_text_size,
    show_value = show_value,
    ...
  )
}

#' Matrix sparse fallback for correlation sparse output
#' @keywords internal
#' @exportS3Method plot sparseMatrix
#' @noRd
plot.sparseMatrix <- function(x,
                              title = NULL,
                              low_color = "indianred1",
                              high_color = "steelblue1",
                              mid_color = "white",
                              value_text_size = 4,
                              ci_text_size = 3,
                              show_value = TRUE,
                              ...) {
  if (identical(attr(x, "corr_output_class", exact = TRUE), "corr_sparse") ||
      identical(attr(x, "output", exact = TRUE), "sparse")) {
    return(plot.corr_sparse(
      x,
      title = title,
      low_color = low_color,
      high_color = high_color,
      mid_color = mid_color,
      value_text_size = value_text_size,
      ci_text_size = ci_text_size,
      show_value = show_value,
      ...
    ))
  }
  graphics::plot(as.matrix(x), ...)
}

#' @keywords internal
#' @exportS3Method plot dsCMatrix
#' @noRd
plot.dsCMatrix <- function(x,
                           title = NULL,
                           low_color = "indianred1",
                           high_color = "steelblue1",
                           mid_color = "white",
                           value_text_size = 4,
                           ci_text_size = 3,
                           show_value = TRUE,
                           ...) {
  if (identical(attr(x, "corr_output_class", exact = TRUE), "corr_sparse") ||
      identical(attr(x, "output", exact = TRUE), "sparse")) {
    return(plot.corr_sparse(
      x,
      title = title,
      low_color = low_color,
      high_color = high_color,
      mid_color = mid_color,
      value_text_size = value_text_size,
      ci_text_size = ci_text_size,
      show_value = show_value,
      ...
    ))
  }
  NextMethod()
}

#' @keywords internal
#' @exportS3Method plot dgCMatrix
#' @noRd
plot.dgCMatrix <- function(x,
                           title = NULL,
                           low_color = "indianred1",
                           high_color = "steelblue1",
                           mid_color = "white",
                           value_text_size = 4,
                           ci_text_size = 3,
                           show_value = TRUE,
                           ...) {
  if (identical(attr(x, "corr_output_class", exact = TRUE), "corr_sparse") ||
      identical(attr(x, "output", exact = TRUE), "sparse")) {
    return(plot.corr_sparse(
      x,
      title = title,
      low_color = low_color,
      high_color = high_color,
      mid_color = mid_color,
      value_text_size = value_text_size,
      ci_text_size = ci_text_size,
      show_value = show_value,
      ...
    ))
  }
  NextMethod()
}

#' @keywords internal
#' @noRd
.mc_summary_corr_result <- function(object,
                                    output_class,
                                    topn = NULL,
                                    show_ci = NULL) {
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )
  all_pairs <- .mc_corr_pairs(object, diag = TRUE)
  raw_pairs <- .mc_corr_pairs(object, diag = FALSE, indexed_pairs = all_pairs)
  estimator_class <- attr(object, "corr_estimator_class", exact = TRUE) %||%
    .mc_corr_estimator_class(object)
  estimator_class <- if (is.character(estimator_class) &&
    length(estimator_class) >= 1L &&
    !is.na(estimator_class[[1L]]) &&
    nzchar(estimator_class[[1L]])) {
    estimator_class[[1L]]
  } else {
    "corr_estimator"
  }
  pairs <- data.frame(
    item1 = as.character(raw_pairs$row %||% character()),
    item2 = as.character(raw_pairs$col %||% character()),
    estimate = as.numeric(raw_pairs$value %||% numeric()),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  dm <- .mc_corr_dim(object)
  ii <- as.integer(raw_pairs$i)
  jj <- as.integer(raw_pairs$j)
  idx_ok <- is.finite(ii) & is.finite(jj)
  idx <- list(i = ii, j = jj, ok = idx_ok)

  diag_attr <- attr(object, "diagnostics", exact = TRUE)
  if (is.list(diag_attr) && is.matrix(diag_attr$n_complete)) {
    pairs$n_complete <- as.integer(round(.mc_extract_pair_matrix(
      diag_attr$n_complete, idx, dm, symmetric = FALSE
    )))
  }
  if (identical(estimator_class, "chatterjee_xi") && is.list(diag_attr)) {
    if (is.matrix(diag_attr$ci_method)) {
      method_vals <- rep(NA_character_, nrow(pairs))
      method_vals[idx_ok] <- as.character(diag_attr$ci_method[cbind(ii[idx_ok], jj[idx_ok])])
      pairs$ci_method <- method_vals
    }
    if (is.matrix(diag_attr$m)) {
      pairs$m <- as.integer(round(.mc_extract_pair_matrix(
        diag_attr$m, idx, dm, symmetric = FALSE
      )))
    }
    if (is.matrix(diag_attr$se)) {
      pairs$se <- as.numeric(.mc_extract_pair_matrix(
        diag_attr$se, idx, dm, symmetric = FALSE
      ))
    }
    if (!is.null(diag_attr$bootstrap_reps)) {
      pairs$bootstrap_reps <- rep(as.integer(diag_attr$bootstrap_reps[[1L]]), nrow(pairs))
    }
  }

  ci_attr <- attr(object, "ci", exact = TRUE)
  if (is.list(ci_attr)) {
    if (is.matrix(ci_attr$lwr.ci)) {
      pairs$lwr <- as.numeric(.mc_extract_pair_matrix(ci_attr$lwr.ci, idx, dm, symmetric = FALSE))
    }
    if (is.matrix(ci_attr$upr.ci)) {
      pairs$upr <- as.numeric(.mc_extract_pair_matrix(ci_attr$upr.ci, idx, dm, symmetric = FALSE))
    }
  }

  inf <- attr(object, "inference", exact = TRUE)
  if (is.list(inf)) {
    if (is.matrix(inf$bcdcor)) {
      pairs$bcdcor <- as.numeric(.mc_extract_pair_matrix(inf$bcdcor, idx, dm, symmetric = FALSE))
    }
    if (is.matrix(inf$statistic)) {
      pairs$statistic <- as.numeric(.mc_extract_pair_matrix(inf$statistic, idx, dm, symmetric = FALSE))
    }
    param_mat <- inf$parameter %||% inf$df
    if (is.matrix(param_mat)) {
      pairs$df <- as.numeric(.mc_extract_pair_matrix(param_mat, idx, dm, symmetric = FALSE))
    }
    if (is.matrix(inf$p_value)) {
      pairs$p_value <- as.numeric(.mc_extract_pair_matrix(inf$p_value, idx, dm, symmetric = FALSE))
    }
    if (identical(estimator_class, "bicor") && is.matrix(inf$Z)) {
      pairs$z_statistic <- as.numeric(.mc_extract_pair_matrix(inf$Z, idx, dm, symmetric = FALSE))
    }
  } else if (identical(estimator_class, "pearson_corr")) {
    z <- suppressWarnings(atanh(pairs$estimate))
    z[!is.finite(z)] <- NA_real_
    pairs$fisher_z <- z
    if ("n_complete" %in% names(pairs) && nrow(pairs)) {
      n <- as.numeric(pairs$n_complete)
      eligible <- is.finite(n) & n > 3
      if (any(eligible)) {
        pairs$statistic <- rep(NA_real_, nrow(pairs))
        pairs$p_value <- rep(NA_real_, nrow(pairs))
        se <- 1 / sqrt(n[eligible] - 3)
        pairs$statistic[eligible] <- z[eligible] / se
        pairs$p_value[eligible] <- 2 * stats::pnorm(
          abs(pairs$statistic[eligible]),
          lower.tail = FALSE
        )
      }
    }
  }

  if (identical(estimator_class, "pearson_corr") && !"fisher_z" %in% names(pairs)) {
    z <- suppressWarnings(atanh(pairs$estimate))
    z[!is.finite(z)] <- NA_real_
    pairs$fisher_z <- z
  }

  if (is.list(inf) && is.matrix(inf$p_value_adjusted)) {
    pairs$p_value_adjusted <- as.numeric(.mc_extract_pair_matrix(
      inf$p_value_adjusted, idx, dm, symmetric = FALSE
    ))
  }
  if (is.list(inf) && is.matrix(inf$reject)) {
    pairs$reject <- .mc_extract_pair_matrix(
      inf$reject, idx, dm, logical = TRUE, symmetric = FALSE
    )
  }
  if (is.list(diag_attr) && is.matrix(diag_attr$skipped_n)) {
    pairs$skipped_n <- as.integer(round(.mc_extract_pair_matrix(
      diag_attr$skipped_n, idx, dm, symmetric = FALSE
    )))
  }
  if (is.list(diag_attr) && is.matrix(diag_attr$skipped_prop)) {
    pairs$skipped_prop <- as.numeric(.mc_extract_pair_matrix(
      diag_attr$skipped_prop, idx, dm, symmetric = FALSE
    ))
  }

  if (nrow(pairs)) {
    pairs <- pairs[order(abs(pairs$estimate), decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  }
  rownames(pairs) <- NULL

  overview <- .mc_corr_overview(object, indexed_pairs = all_pairs)
  summary_classes <- c(
    "summary.corr_result",
    output_class,
    paste0("summary.", estimator_class),
    if (estimator_class %in% c("tetrachoric_corr", "polychoric_corr", "polyserial_corr", "biserial_corr")) {
      "summary.latent_corr"
    } else {
      character()
    },
    "summary.matrixCorr",
    "data.frame"
  )
  attr(pairs, "class") <- NULL
  class(pairs) <- unique(summary_classes)

  topn_use <- .mc_coalesce(topn, .mc_display_option("summary_topn", 5L))
  topn_use <- .mc_validate_optional_count(topn_use, arg = "topn", allow_null = FALSE)
  top_tbl <- if (nrow(pairs)) {
    pairs[utils::head(seq_len(nrow(pairs)), topn_use), , drop = FALSE]
  } else {
    pairs
  }
  attr(pairs, "overview") <- overview
  attr(pairs, "method") <- overview$method
  attr(pairs, "output") <- overview$output
  attr(pairs, "threshold") <- overview$threshold
  attr(pairs, "diag") <- overview$diag
  attr(pairs, "has_ci") <- any(c("lwr", "upr") %in% names(pairs))
  attr(pairs, "has_p") <- "p_value" %in% names(pairs)
  attr(pairs, "conf.level") <- overview$conf.level
  ci_method_attr <- attr(object, "ci", exact = TRUE)$ci.method %||%
    attr(object, "ci.method", exact = TRUE) %||%
    NA_character_
  if (is.list(ci_method_attr)) {
    ci_method_attr <- NA_character_
  } else if (length(ci_method_attr) != 1L) {
    method_values <- unique(as.character(ci_method_attr))
    method_values <- method_values[!is.na(method_values) & nzchar(method_values)]
    ci_method_attr <- if (length(method_values) == 1L) {
      method_values
    } else if (length(method_values) > 1L) {
      "mixed"
    } else {
      NA_character_
    }
  }
  ci_method_attr <- as.character(ci_method_attr)[[1L]]
  attr(pairs, "ci_method") <- ci_method_attr
  if ((is.null(attr(pairs, "ci_method", exact = TRUE)) ||
       is.na(attr(pairs, "ci_method", exact = TRUE)) ||
       !nzchar(attr(pairs, "ci_method", exact = TRUE))) &&
      isTRUE(attr(pairs, "has_ci"))) {
    ci_map <- c(
      pearson_corr = "fisher_z",
      spearman_rho = "jackknife_euclidean_likelihood",
      kendall_matrix = "fieller"
    )
    attr(pairs, "ci_method") <- unname(ci_map[estimator_class] %||% NA_character_)
  }
  attr(pairs, "top_results") <- top_tbl
  attr(pairs, "show_ci") <- show_ci
  attr(pairs, "summary_title") <- .mc_corr_summary_title(
    estimator_class,
    has_ci = isTRUE(attr(pairs, "has_ci")),
    has_p = isTRUE(attr(pairs, "has_p"))
  )
  pairs
}

#' @keywords internal
#' @noRd
.mc_corr_summary_title <- function(estimator_class, has_ci = FALSE, has_p = FALSE) {
  estimator_class <- if (is.character(estimator_class) &&
    length(estimator_class) >= 1L &&
    !is.na(estimator_class[[1L]]) &&
    nzchar(estimator_class[[1L]])) {
    estimator_class[[1L]]
  } else {
    "corr_estimator"
  }
  if (estimator_class %in% c("tetrachoric_corr", "polychoric_corr", "polyserial_corr", "biserial_corr")) {
    if (!isTRUE(has_ci) && !isTRUE(has_p)) {
      return("Latent correlation summary")
    }
  }
  map <- c(
    pearson_corr = "Pearson correlation summary",
    chatterjee_xi = "Chatterjee xi directed-dependence summary",
    spearman_rho = if (isTRUE(has_ci) || isTRUE(has_p)) "Spearman correlation summary" else "Correlation summary",
    kendall_matrix = "Kendall correlation summary",
    cohen_kappa = "Cohen's kappa agreement summary",
    weighted_kappa = "Weighted Cohen's kappa agreement summary",
    bicor = "Biweight mid-correlation summary",
    dcor = "Distance correlation summary",
    bcdcor = "Bias-corrected squared distance correlation summary",
    lancaster_corr = "Lancaster correlation summary",
    robust_dcor = "Robust distance correlation summary",
    pbcor = "Percentage bend correlation summary",
    wincor = "Winsorized correlation summary",
    skipped_corr = "Skipped correlation summary",
    shrinkage_corr = "Correlation summary",
    schafer_corr = "Correlation summary",
    tetrachoric_corr = "Tetrachoric correlation summary",
    polychoric_corr = "Polychoric correlation summary",
    polyserial_corr = "Polyserial correlation summary",
    biserial_corr = "Biserial correlation summary",
    partial_corr_matrix = "Partial correlation summary"
  )
  label <- unname(map[estimator_class])
  if (length(label) == 0L || is.na(label) || !nzchar(label)) {
    "Correlation summary"
  } else {
    label
  }
}

#' Summary Accessor for Correlation Summaries
#' @param x A `summary.corr_result` object.
#' @param name A column name or summary metadata key.
#' @return A summary column (if present) or summary metadata entry.
#' @method $ summary.corr_result
#' @export
`$.summary.corr_result` <- function(x, name) {
  if (name %in% names(x)) {
    return(NextMethod("$"))
  }
  if (identical(name, "top_results")) {
    return(attr(x, "top_results", exact = TRUE))
  }
  overview <- attr(x, "overview", exact = TRUE)
  if (is.list(overview) && name %in% names(overview)) {
    return(overview[[name]])
  }
  attr(x, name, exact = TRUE)
}

#' Summary Accessor for Correlation Summaries
#' @param x A `summary.corr_result` object.
#' @param i A column name or summary metadata key.
#' @param ... Unused.
#' @return A summary column (if present) or summary metadata entry.
#' @method [[ summary.corr_result
#' @export
`[[.summary.corr_result` <- function(x, i, ...) {
  if (is.character(i) && length(i) == 1L && !i %in% names(x)) {
    return(`$.summary.corr_result`(x, i))
  }
  NextMethod("[[")
}

#' @keywords internal
#' @noRd
.mc_method_label <- function(x) {
  method <- attr(x, "method", exact = TRUE) %||%
    attr(x, "corr_estimator_class", exact = TRUE) %||%
    "correlation"
  if (!is.character(method) || length(method) != 1L || is.na(method) || !nzchar(method)) {
    return("Correlation")
  }
  paste0(toupper(substr(method, 1L, 1L)), substr(method, 2L, nchar(method)))
}

#' @keywords internal
#' @noRd
.mc_print_corr_records <- function(x,
                                   output_label,
                                   digits = 4,
                                   n = NULL,
                                   topn = NULL,
                                   max_vars = NULL,
                                   width = NULL,
                                   show_ci = NULL,
                                   ...) {
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

  df <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  header <- sprintf("%s correlation %s", .mc_method_label(x), output_label)
  digest <- c(
    method = attr(x, "method", exact = TRUE) %||% NA_character_,
    dimensions = sprintf("%d x %d", nrow(df), ncol(df))
  )
  digest <- digest[!is.na(digest) & nzchar(digest)]
  .mc_print_named_digest(digest, header = header)
  cat("\n")

  if ("value" %in% names(df)) {
    df$value <- formatC(as.numeric(df$value), format = "f", digits = digits)
  }
  .mc_print_preview_table(
    df,
    n = cfg$n,
    topn = cfg$topn,
    max_vars = cfg$max_vars,
    width = cfg$width,
    context = "print",
    full_hint = TRUE,
    summary_hint = TRUE
  )
}

#' @keywords internal
#' @noRd
.mc_plot_grid_ci <- function(ci_attr, dm, ii, jj) {
  n <- length(ii)
  lwr <- rep(NA_real_, n)
  upr <- rep(NA_real_, n)
  has_ci <- is.list(ci_attr) &&
    is.matrix(ci_attr$lwr.ci) &&
    is.matrix(ci_attr$upr.ci) &&
    identical(dim(ci_attr$lwr.ci), dm) &&
    identical(dim(ci_attr$upr.ci), dm)
  if (isTRUE(has_ci) && n) {
    ok <- is.finite(ii) & is.finite(jj) &
      ii >= 1L & ii <= dm[[1L]] &
      jj >= 1L & jj <= dm[[2L]]
    if (any(ok)) {
      lwr[ok] <- ci_attr$lwr.ci[cbind(ii[ok], jj[ok])]
      upr[ok] <- ci_attr$upr.ci[cbind(ii[ok], jj[ok])]
    }
  }
  list(lwr = as.numeric(lwr), upr = as.numeric(upr))
}

.mc_plot_grid_mirror <- function(edf, off, row, col, i = NULL, j = NULL) {
  if (!any(off)) {
    return(edf)
  }
  mirror <- edf[off, , drop = FALSE]
  if (!is.null(i)) {
    mirror$i <- i
    mirror$j <- j
  }
  mirror$row <- row
  mirror$col <- col
  rbind(edf, mirror)
}

.mc_corr_plot_grid <- function(x) {
  dm <- .mc_corr_dim(x)
  dn <- .mc_corr_dimnames(x)
  ambiguous <- .mc_axis_names_positional(dn[[1L]], dn[[2L]])
  if (ambiguous) {
    edf <- .mc_resolve_indexed_pairs(x)
  } else {
    edf <- .mc_corr_as_edge_df(x)
  }
  if (!nrow(edf)) {
    return(data.frame(
      row = character(),
      col = character(),
      value = numeric(),
      lwr = numeric(),
      upr = numeric(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  if (ambiguous) {
    ii <- edf$i
    jj <- edf$j
  } else {
    edf$row <- as.character(edf$row)
    edf$col <- as.character(edf$col)
    edf$value <- as.numeric(edf$value)
    ii <- if (!is.null(dn[[1L]])) match(edf$row, dn[[1L]]) else suppressWarnings(as.integer(edf$row))
    jj <- if (!is.null(dn[[2L]])) match(edf$col, dn[[2L]]) else suppressWarnings(as.integer(edf$col))
  }
  bounds <- .mc_plot_grid_ci(attr(x, "ci", exact = TRUE), dm, ii, jj)
  edf$lwr <- bounds$lwr
  edf$upr <- bounds$upr

  if (isTRUE(attr(x, "corr_symmetric", exact = TRUE))) {
    if (ambiguous) {
      off <- edf$i != edf$j
      edf <- .mc_plot_grid_mirror(
        edf,
        off,
        row = if (is.null(dn[[1L]])) as.character(edf$j[off]) else dn[[1L]][edf$j[off]],
        col = if (is.null(dn[[2L]])) as.character(edf$i[off]) else dn[[2L]][edf$i[off]],
        i = edf$j[off],
        j = edf$i[off]
      )
    } else {
      off <- edf$row != edf$col
      edf <- .mc_plot_grid_mirror(
        edf,
        off,
        row = edf$col[off],
        col = edf$row[off]
      )
    }
  }

  if (ambiguous) {
    edf$row <- factor(edf$i, levels = rev(seq_len(dm[[1L]])))
    edf$col <- factor(edf$j, levels = seq_len(dm[[2L]]))
    out <- edf[, c("row", "col", "value", "lwr", "upr"), drop = FALSE]
    attr(out, "mc_axis_labels") <- list(
      row = rev(dn[[1L]] %||% as.character(seq_len(dm[[1L]]))),
      col = dn[[2L]] %||% as.character(seq_len(dm[[2L]]))
    )
    return(out)
  }

  if (!is.null(dn[[1L]])) {
    edf$row <- factor(edf$row, levels = rev(dn[[1L]]))
  } else {
    edf$row <- factor(edf$row, levels = rev(as.character(seq_len(dm[[1L]]))))
  }
  if (!is.null(dn[[2L]])) {
    edf$col <- factor(edf$col, levels = dn[[2L]])
  } else {
    edf$col <- factor(edf$col, levels = as.character(seq_len(dm[[2L]])))
  }
  edf
}

#' @keywords internal
#' @noRd
.mc_plot_corr_result <- function(x,
                                 title,
                                 low_color,
                                 high_color,
                                 mid_color,
                                 value_text_size,
                                 ci_text_size,
                                 show_value,
                                 ...) {
  dots <- list(...)
  if ("type" %in% names(dots)) {
    abort_bad_arg(
      "type",
      message = "is not used for matrix-style correlation heatmaps.",
      .hint = "For estimator-specific plot modes such as {.arg type} = {.val estimate}, use a non-matrix result object."
    )
  }
  check_bool(show_value, arg = "show_value")
  df <- .mc_corr_plot_grid(x)
  axis_labels <- attr(df, "mc_axis_labels", exact = TRUE)
  if (!"ci_label" %in% names(df)) {
    has_ci_cols <- all(c("lwr", "upr") %in% names(df))
    if (isTRUE(has_ci_cols)) {
      same_cell <- .mc_same_axis_square(x) &
        (as.character(df$row) == as.character(df$col))
      df$ci_label <- ifelse(
        is.finite(df$lwr) & is.finite(df$upr) & !same_cell,
        sprintf("[%.3f, %.3f]", df$lwr, df$upr),
        NA_character_
      )
    } else {
      df$ci_label <- NA_character_
    }
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$col, .data$row, fill = .data$value)) +
    ggplot2::geom_tile(color = "white")

  estimator_class <- .mc_corr_estimator_class(x)
  if (identical(estimator_class, "dcor")) {
    squared <- isTRUE(attr(x, "squared", exact = TRUE))
    fill_name <- if (squared) "dCor^2" else "dCor"
    if (identical(title, "Correlation heatmap")) {
      title <- if (squared) "Squared distance correlation heatmap" else "Distance correlation heatmap"
    }
    p <- p + ggplot2::scale_fill_gradient(
      low = mid_color,
      high = high_color,
      limits = c(0, 1),
      name = fill_name
    )
  } else {
    fill_name <- if (identical(estimator_class, "bcdcor")) {
      if (identical(title, "Correlation heatmap")) {
        title <- "Bias-corrected squared distance correlation heatmap"
      }
      "bc-dCor^2"
    } else {
      "Correlation"
    }
    p <- p + ggplot2::scale_fill_gradient2(
      low = low_color,
      high = high_color,
      mid = mid_color,
      midpoint = 0,
      limits = if (identical(estimator_class, "bcdcor")) NULL else c(-1, 1),
      name = fill_name
    )
  }

  p <- p +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      ...
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = NULL, y = NULL)
  p <- .mc_add_positional_axis_scales(p, axis_labels)

  if (isTRUE(show_value) && is.finite(value_text_size) && !is.null(value_text_size)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", .data$value)),
      size = value_text_size,
      color = "black"
    )
  }
  if (isTRUE(show_value) &&
      is.finite(ci_text_size) &&
      !is.null(ci_text_size) &&
      any(!is.na(df$ci_label))) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$ci_label, y = as.numeric(.data$row) - 0.25),
      size = ci_text_size,
      color = "gray30",
      na.rm = TRUE
    )
  }
  p
}

.mc_summary_corr_matrix <- function(object,
                                    header = "Correlation summary",
                                    topn = NULL) {
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
    correct = attr(object, "correct"),
    n_variables = if (symmetric) nrow(m) else NA_integer_,
    n_pairs = as.integer(n_pairs),
    n_rows = nrow(m),
    n_cols = ncol(m),
    symmetric = symmetric,
    n_missing = sum(is.na(m)),
    estimate_min = if (length(vals)) min(vals) else NA_real_,
    estimate_max = if (length(vals)) max(vals) else NA_real_,
    threshold_sets = {
      thr <- attr(object, "thresholds")
      if (is.list(thr)) length(thr) else NA_integer_
    },
    has_ci = .mc_has_ci_payload(object),
    has_p = .mc_has_p_payload(object)
  )
  extremes <- .mc_pair_extremes(object, digits = 4)
  out$most_negative <- extremes$most_negative
  out$most_positive <- extremes$most_positive

  diag_attr <- attr(object, "diagnostics")
  if (is.list(diag_attr)) {
    selector <- if (symmetric && nrow(m) > 1L) upper.tri(m, diag = FALSE) else matrix(TRUE, nrow(m), ncol(m))
    pair_pick <- function(x, transform = identity, default = NA_integer_) {
      if (!is.matrix(x) || !identical(dim(x), dim(m))) return(default)
      vals <- transform(x[selector])
      vals <- vals[!is.na(vals)]
      vals <- vals[is.finite(vals)]
      if (!length(vals)) return(default)
      vals
    }

    zero_vals <- pair_pick(diag_attr$zero_cells)
    ncomp_vals <- pair_pick(diag_attr$n_complete)
    skipped_n_vals <- pair_pick(diag_attr$skipped_n, default = integer())
    skipped_prop_vals <- pair_pick(diag_attr$skipped_prop, default = numeric())
    boundary_vals <- pair_pick(diag_attr$boundary, as.logical, NA)
    near_vals <- pair_pick(diag_attr$near_boundary, as.logical, NA)
    corrected_vals <- pair_pick(diag_attr$corrected, as.logical, NA)
    conv_vals <- pair_pick(diag_attr$converged, as.logical, NA)

    out$zero_cell_pairs <- if (is.atomic(zero_vals) && length(zero_vals)) sum(zero_vals > 0, na.rm = TRUE) else NA_integer_
    out$boundary_pairs <- if (is.atomic(boundary_vals) && length(boundary_vals)) sum(boundary_vals, na.rm = TRUE) else NA_integer_
    out$near_boundary_pairs <- if (is.atomic(near_vals) && length(near_vals)) sum(near_vals, na.rm = TRUE) else NA_integer_
    out$corrected_pairs <- if (is.atomic(corrected_vals) && length(corrected_vals)) sum(corrected_vals, na.rm = TRUE) else NA_integer_
    out$converged_pairs <- if (is.atomic(conv_vals) && length(conv_vals)) sum(conv_vals, na.rm = TRUE) else NA_integer_
    out$n_complete_min <- if (is.atomic(ncomp_vals) && length(ncomp_vals)) min(ncomp_vals, na.rm = TRUE) else NA_integer_
    out$n_complete_max <- if (is.atomic(ncomp_vals) && length(ncomp_vals)) max(ncomp_vals, na.rm = TRUE) else NA_integer_
    out$skipped_n_min <- if (is.atomic(skipped_n_vals) && length(skipped_n_vals)) min(skipped_n_vals, na.rm = TRUE) else NA_integer_
    out$skipped_n_max <- if (is.atomic(skipped_n_vals) && length(skipped_n_vals)) max(skipped_n_vals, na.rm = TRUE) else NA_integer_
    out$skipped_prop_min <- if (is.atomic(skipped_prop_vals) && length(skipped_prop_vals)) min(skipped_prop_vals, na.rm = TRUE) else NA_real_
    out$skipped_prop_max <- if (is.atomic(skipped_prop_vals) && length(skipped_prop_vals)) max(skipped_prop_vals, na.rm = TRUE) else NA_real_
      out$optimizer_tol <- if (is.null(diag_attr$optimizer_tol)) NA_real_ else diag_attr$optimizer_tol
    } else {
    out$zero_cell_pairs <- NA_integer_
    out$boundary_pairs <- NA_integer_
    out$near_boundary_pairs <- NA_integer_
    out$corrected_pairs <- NA_integer_
    out$converged_pairs <- NA_integer_
    out$n_complete_min <- NA_integer_
    out$n_complete_max <- NA_integer_
    out$skipped_n_min <- NA_integer_
    out$skipped_n_max <- NA_integer_
    out$skipped_prop_min <- NA_real_
    out$skipped_prop_max <- NA_real_
      out$optimizer_tol <- NA_real_
    }

  ci_attr <- attr(object, "ci", exact = TRUE)
  out$ci_conf_level <- NA_real_
  out$ci_method <- NA_character_
  out$ci_width_min <- NA_real_
  out$ci_width_max <- NA_real_
  out$ci_cross_zero <- NA_integer_
  if (is.list(ci_attr) &&
      is.matrix(ci_attr$lwr.ci) &&
      is.matrix(ci_attr$upr.ci) &&
      identical(dim(ci_attr$lwr.ci), dim(m)) &&
      identical(dim(ci_attr$upr.ci), dim(m))) {
    selector <- if (symmetric && nrow(m) > 1L) upper.tri(m, diag = FALSE) else matrix(TRUE, nrow(m), ncol(m))
    lwr_vals <- ci_attr$lwr.ci[selector]
    upr_vals <- ci_attr$upr.ci[selector]
    keep_ci <- is.finite(lwr_vals) & is.finite(upr_vals)
    out$ci_conf_level <- suppressWarnings(as.numeric(ci_attr$conf.level %||% attr(object, "conf.level", exact = TRUE)))
    out$ci_method <- as.character(ci_attr$ci.method %||% attr(object, "ci.method", exact = TRUE) %||% NA_character_)
    if (any(keep_ci)) {
      widths <- upr_vals[keep_ci] - lwr_vals[keep_ci]
      out$ci_width_min <- min(widths)
      out$ci_width_max <- max(widths)
      out$ci_cross_zero <- sum(lwr_vals[keep_ci] <= 0 & upr_vals[keep_ci] >= 0, na.rm = TRUE)
    }
  }

  out$header <- header
  out$top_results <- .mc_summary_top_pairs(
    object,
    digits = 4,
    topn = .mc_coalesce(topn, .mc_display_option("summary_topn", 5L))
  )
  out <- .mc_finalize_summary_list(
    out,
    class_name = if (identical(header, "Latent correlation summary")) {
      "summary.latent_corr"
    } else {
      "summary.matrixCorr"
    }
  )
  if (identical(class(out)[1L], "summary.matrixCorr")) {
    class(out) <- c("summary.matrixCorr", "summary.corr_matrix")
  } else {
    class(out) <- c("summary.latent_corr", "summary.matrixCorr", "summary.corr_matrix")
  }
  out
}

#' Correlation overview abstraction
#' @keywords internal
#' @noRd
.mc_corr_overview <- function(x, indexed_pairs = NULL) {
  dm <- .mc_corr_dim(x)
  edges <- if (is.null(indexed_pairs)) .mc_corr_pairs(x, diag = TRUE) else indexed_pairs
  method <- attr(x, "method", exact = TRUE) %||% .mc_corr_estimator_class(x)
  output <- attr(x, "output", exact = TRUE) %||% "matrix"
  thr <- attr(x, "threshold", exact = TRUE) %||% 0
  use_diag <- attr(x, "diag", exact = TRUE) %||% TRUE
  symmetric <- isTRUE(attr(x, "corr_symmetric", exact = TRUE))
  pair_count <- if (!isTRUE(symmetric)) {
    nrow(edges)
  } else if (!nrow(edges)) {
    0L
  } else {
    keep <- edges$i != edges$j
    if (!anyNA(keep)) sum(keep) else nrow(edges[keep, , drop = FALSE])
  }
  has_ci <- !is.null(attr(x, "ci", exact = TRUE))
  inf_attr <- attr(x, "inference", exact = TRUE)
  has_p <- is.list(inf_attr) && is.matrix(inf_attr$p_value)
  conf.level <- attr(x, "conf.level", exact = TRUE) %||%
    (attr(x, "ci", exact = TRUE)$conf.level %||% NA_real_)
  ci_method <- attr(x, "ci", exact = TRUE)$ci.method %||%
    attr(x, "ci.method", exact = TRUE) %||%
    NA_character_
  inference_method <- if (is.list(inf_attr) && is.character(inf_attr$method) &&
                          length(inf_attr$method) == 1L) {
    inf_attr$method
  } else {
    NA_character_
  }

  matrix_output <- identical(output, "matrix") || inherits(x, "corr_matrix")
  use_matrix_overview <- isTRUE(matrix_output)
  selector <- NULL
  if (isTRUE(use_matrix_overview)) {
    m <- .mc_corr_as_dense_matrix(x)
    selector <- if (symmetric && dm[[1L]] == dm[[2L]]) {
      upper.tri(m, diag = FALSE)
    } else {
      matrix(TRUE, nrow = dm[[1L]], ncol = dm[[2L]])
    }
    vals <- if (is.matrix(m) && identical(dim(m), dm)) {
      as.numeric(m[selector])
    } else {
      numeric()
    }
    n_missing <- if (is.matrix(m)) sum(is.na(m)) else NA_integer_
  } else {
    vals <- as.numeric(edges$value)
    if (symmetric && dm[[1L]] == dm[[2L]] && nrow(edges)) {
      vals <- vals[edges$i != edges$j]
    }
    n_missing <- sum(is.na(edges$value))
  }
  vals <- vals[is.finite(vals)]
  # Source-wide diagnostic range: strict upper triangle when the result is
  # marked symmetric and square, otherwise the full rectangle. The range does
  # not follow the retained edge or threshold filter.
  pick_range <- function(mat) {
    if (!is.matrix(mat) || !identical(dim(mat), dm)) {
      return(c(NA_real_, NA_real_))
    }
    idx <- selector
    if (is.null(idx)) {
      idx <- if (symmetric && dm[[1L]] == dm[[2L]]) {
        upper.tri(mat, diag = FALSE)
      } else {
        matrix(TRUE, nrow = dm[[1L]], ncol = dm[[2L]])
      }
    }
    vv <- as.numeric(mat[idx])
    vv <- vv[is.finite(vv)]
    if (!length(vv)) {
      return(c(NA_real_, NA_real_))
    }
    c(min(vv), max(vv))
  }
  diag_attr <- attr(x, "diagnostics", exact = TRUE)
  skipped_n <- if (is.list(diag_attr)) pick_range(diag_attr$skipped_n) else c(NA_real_, NA_real_)
  skipped_prop <- if (is.list(diag_attr)) pick_range(diag_attr$skipped_prop) else c(NA_real_, NA_real_)
  threshold_sets <- {
    threshold_sets_attr <- attr(x, "thresholds", exact = TRUE)
    if (is.list(threshold_sets_attr)) length(threshold_sets_attr) else NA_integer_
  }
  thr_num <- if (is.list(thr)) numeric() else suppressWarnings(as.numeric(thr))
  if (!length(thr_num) || !is.finite(thr_num[[1L]])) {
    thr_num <- 0
  } else {
    thr_num <- thr_num[[1L]]
  }
  conf_num <- suppressWarnings(as.numeric(conf.level))
  if (!length(conf_num) || !is.finite(conf_num[[1L]])) {
    conf_num <- NA_real_
  } else {
    conf_num <- conf_num[[1L]]
  }
  if (length(ci_method) != 1L || is.list(ci_method)) {
    ci_method <- NA_character_
  }
  if (length(inference_method) != 1L || is.list(inference_method)) {
    inference_method <- NA_character_
  }

  list(
    class = .mc_corr_estimator_class(x),
    method = method,
    output = output,
    threshold = thr_num,
    diag = isTRUE(use_diag),
    n_variables = if (symmetric) as.integer(dm[[1L]]) else NA_integer_,
    n_rows = as.integer(dm[[1L]]),
    n_cols = as.integer(dm[[2L]]),
    n_pairs = pair_count,
    n_pairs_retained = pair_count,
    symmetric = symmetric,
    n_missing = n_missing,
    threshold_sets = threshold_sets,
    estimate_min = if (length(vals)) min(vals) else NA_real_,
    estimate_max = if (length(vals)) max(vals) else NA_real_,
    has_ci = has_ci,
    has_p = has_p,
    conf.level = conf_num,
    ci_method = ci_method,
    inference_method = inference_method,
    skipped_n_min = skipped_n[[1L]],
    skipped_n_max = skipped_n[[2L]],
    skipped_prop_min = skipped_prop[[1L]],
    skipped_prop_max = skipped_prop[[2L]]
  )
}

.mc_round_finite <- function(x, digits) {
  x <- as.numeric(x)
  if (is.null(digits)) {
    return(x)
  }
  out <- rep(NA_real_, length(x))
  keep <- is.finite(x)
  out[keep] <- round(x[keep], digits)
  out
}

.mc_pairwise_matrix_summary <- function(object,
                                        class_name,
                                        digits,
                                        ci_digits = NULL,
                                        p_digits = NULL,
                                        show_ci = NULL,
                                        ci_attr = .mc_ci_attr(object),
                                        inference_attr = .mc_inference_attr(object),
                                        diagnostics_attr = attr(object, "diagnostics", exact = TRUE),
                                        include_n_complete = TRUE,
                                        include_ci = TRUE,
                                        include_p = TRUE,
                                        extra_columns = NULL,
                                        extra_attrs = list(),
                                        overview = NULL) {
  est <- as.matrix(object)
  nr <- nrow(est)
  nc <- ncol(est)
  rn <- rownames(est)
  cn <- colnames(est)
  if (is.null(rn)) rn <- as.character(seq_len(nr))
  if (is.null(cn)) cn <- as.character(seq_len(nc))

  n_pairs <- nr * (nc - 1L) / 2L
  var1 <- character(n_pairs)
  var2 <- character(n_pairs)
  estimate <- numeric(n_pairs)

  has_n_complete <- isTRUE(include_n_complete) &&
    is.list(diagnostics_attr) &&
    is.matrix(diagnostics_attr$n_complete)
  n_complete <- if (has_n_complete) integer(n_pairs) else NULL

  show_ci <- if (is.null(show_ci)) "yes" else show_ci
  has_ci <- isTRUE(include_ci) &&
    identical(show_ci, "yes") &&
    is.list(ci_attr)
  lwr <- if (has_ci) numeric(n_pairs) else NULL
  upr <- if (has_ci) numeric(n_pairs) else NULL

  extra_columns <- extra_columns %||% list()
  extra_values <- vector("list", length(extra_columns))
  names(extra_values) <- names(extra_columns)
  if (length(extra_columns)) {
    for (nm in names(extra_columns)) {
      spec <- extra_columns[[nm]]
      type <- spec$type %||% "numeric"
      extra_values[[nm]] <- switch(
        type,
        integer = integer(n_pairs),
        logical = logical(n_pairs),
        character = character(n_pairs),
        numeric(n_pairs)
      )
    }
  }

  k <- 0L
  for (i in seq_len(nr - 1L)) {
    for (j in (i + 1L):nc) {
      k <- k + 1L
      var1[[k]] <- rn[[i]]
      var2[[k]] <- cn[[j]]
      estimate[[k]] <- round(est[i, j], digits)
      if (!is.null(n_complete)) {
        n_complete[[k]] <- as.integer(diagnostics_attr$n_complete[i, j])
      }
      if (has_ci) {
        lwr[[k]] <- if (!is.null(ci_attr$lwr.ci) && is.finite(ci_attr$lwr.ci[i, j])) {
          round(ci_attr$lwr.ci[i, j], ci_digits)
        } else {
          NA_real_
        }
        upr[[k]] <- if (!is.null(ci_attr$upr.ci) && is.finite(ci_attr$upr.ci[i, j])) {
          round(ci_attr$upr.ci[i, j], ci_digits)
        } else {
          NA_real_
        }
      }
      if (length(extra_columns)) {
        for (nm in names(extra_columns)) {
          spec <- extra_columns[[nm]]
          mat <- spec$matrix
          value <- if (is.matrix(mat)) {
            mat[i, j]
          } else if (is.function(spec$value)) {
            spec$value(i, j)
          } else {
            NA
          }
          type <- spec$type %||% "numeric"
          extra_values[[nm]][[k]] <- switch(
            type,
            integer = as.integer(value),
            logical = isTRUE(value),
            character = as.character(value),
            .mc_round_finite(value, if ("digits" %in% names(spec)) spec$digits else digits)
          )
        }
      }
    }
  }

  df <- data.frame(
    var1 = var1,
    var2 = var2,
    estimate = as.numeric(estimate),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!is.null(n_complete)) df$n_complete <- as.integer(n_complete)
  if (!is.null(lwr)) df$lwr <- as.numeric(lwr)
  if (!is.null(upr)) df$upr <- as.numeric(upr)
  if (length(extra_values)) {
    for (nm in names(extra_values)) {
      df[[nm]] <- extra_values[[nm]]
    }
  }
  rownames(df) <- NULL

  out <- .mc_finalize_summary_df(df, class_name = class_name)
  attr(out, "overview") <- overview %||% .mc_summary_corr_matrix(object)
  attr(out, "has_ci") <- has_ci
  if (!is.null(include_p)) {
    attr(out, "has_p") <- isTRUE(include_p)
  }
  attr(out, "conf.level") <- if (is.null(ci_attr)) NA_real_ else ci_attr$conf.level
  attr(out, "digits") <- digits
  if (!is.null(ci_digits)) attr(out, "ci_digits") <- ci_digits
  if (!is.null(p_digits)) attr(out, "p_digits") <- p_digits
  if (length(extra_attrs)) {
    for (nm in names(extra_attrs)) {
      attr(out, nm) <- extra_attrs[[nm]]
    }
  }
  out
}
