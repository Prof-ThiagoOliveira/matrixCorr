#' S3 Summary for Dense Correlation Results
#'
#' Representation-first summary for dense correlation outputs.
#'
#' @param object A dense correlation result (`corr_matrix`).
#' @param topn Optional number of head/tail rows when preview is truncated.
#' @param show_ci One of `"yes"` or `"no"`.
#' @param ... Unused.
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

#' S3 Summary for Sparse Correlation Results
#'
#' Representation-first summary for sparse correlation outputs.
#'
#' @param object A sparse correlation result.
#' @param topn Optional number of head/tail rows when preview is truncated.
#' @param show_ci One of `"yes"` or `"no"`.
#' @param ... Unused.
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

#' S3 Summary for Packed-Upper Correlation Results
#'
#' Representation-first summary for packed upper-triangle outputs.
#'
#' @param object A packed upper-triangle correlation result.
#' @param topn Optional number of head/tail rows when preview is truncated.
#' @param show_ci One of `"yes"` or `"no"`.
#' @param ... Unused.
#'
#' @return A standardized summary data frame with class
#'   `c("summary.corr_result", "data.frame")` (plus compatibility classes).
#' @method summary corr_packed_upper
#' @export
summary.corr_packed_upper <- function(object, topn = NULL, show_ci = NULL, ...) {
  .mc_summary_corr_result(
    object,
    output_class = "summary.corr_packed_upper",
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

#' @title Print Packed-Upper Correlation Results
#' @param x A packed-upper correlation result.
#' @param digits Number of digits for numeric values.
#' @param n Optional preview row threshold.
#' @param topn Optional number of head/tail rows when preview is truncated.
#' @param max_vars Optional maximum number of visible columns in preview.
#' @param width Optional output width.
#' @param show_ci One of `"yes"` or `"no"`.
#' @param ... Unused.
#' @return Invisibly returns `x`.
#' @method print corr_packed_upper
#' @export
print.corr_packed_upper <- function(x,
                                    digits = 4,
                                    n = NULL,
                                    topn = NULL,
                                    max_vars = NULL,
                                    width = NULL,
                                    show_ci = NULL,
                                    ...) {
  .mc_print_corr_records(
    x,
    output_label = "packed upper",
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

#' @title Print Edge-List Correlation Results
#' @param x An edge-list correlation result.
#' @inheritParams print.corr_packed_upper
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
    digest <- c(
      digest,
      .mc_ci_digest(
        as.data.frame(x, stringsAsFactors = FALSE),
        conf_level = suppressWarnings(as.numeric(attr(x, "conf.level", exact = TRUE))),
        ci_method = attr(x, "ci_method", exact = TRUE) %||%
          (attr(x, "overview", exact = TRUE)$ci_method %||% NULL),
        digits = 3
      )
    )
  }
  digest <- digest[!is.na(digest) & nzchar(digest)]
  .mc_print_named_digest(digest, header = header)
  cat("\n")

  preview <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  if (identical(cfg$show_ci, "no")) {
    preview <- preview[, setdiff(names(preview), c("lwr", "upr")), drop = FALSE]
  }
  if (nrow(preview)) {
    if ("estimate" %in% names(preview)) {
      preview$estimate <- formatC(as.numeric(preview$estimate), format = "f", digits = digits)
    }
    for (nm in intersect(c("lwr", "upr", "statistic", "df", "p_value", "fisher_z"), names(preview))) {
      preview[[nm]] <- formatC(as.numeric(preview[[nm]]), format = "f", digits = digits)
    }
  }
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
  top_tbl <- attr(x, "top_results", exact = TRUE)
  if (is.data.frame(top_tbl) && nrow(top_tbl)) {
    if (identical(cfg$show_ci, "no")) {
      top_tbl <- top_tbl[, setdiff(names(top_tbl), c("lwr", "upr")), drop = FALSE]
    }
    cat("\nStrongest pairs by |estimate|\n\n")
    .mc_print_preview_table(
      top_tbl,
      n = max(2L, nrow(top_tbl) + 1L),
      topn = min(cfg$topn, max(1L, nrow(top_tbl))),
      max_vars = cfg$max_vars,
      width = cfg$width,
      context = "summary",
      full_hint = TRUE,
      summary_hint = FALSE
    )
  }
  invisible(x)
}

#' S3 Plot for Dense Correlation Results
#'
#' @param x A dense correlation result (`corr_matrix`).
#' @param title Optional plot title.
#' @param low_color Fill color for -1.
#' @param high_color Fill color for +1.
#' @param mid_color Fill color for 0.
#' @param value_text_size Text size for optional overlaid values.
#' @param ci_text_size Reserved for compatibility.
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
    show_value = show_value,
    ...
  )
}

#' S3 Plot for Packed-Upper Correlation Results
#'
#' @inheritParams plot.corr_matrix
#' @param x A packed-upper correlation result.
#' @method plot corr_packed_upper
#' @export
plot.corr_packed_upper <- function(x,
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
    title = title %||% "Packed upper-triangle correlation heatmap",
    low_color = low_color,
    high_color = high_color,
    mid_color = mid_color,
    value_text_size = value_text_size,
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
  raw_pairs <- .mc_corr_pairs(object, diag = FALSE)
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
  dn <- .mc_corr_dimnames(object)
  rn <- dn[[1L]]
  cn <- dn[[2L]]
  ii <- if (is.null(rn)) suppressWarnings(as.integer(pairs$item1)) else match(pairs$item1, rn)
  jj <- if (is.null(cn)) suppressWarnings(as.integer(pairs$item2)) else match(pairs$item2, cn)
  idx_ok <- is.finite(ii) & is.finite(jj)

  add_from_matrix <- function(mat) {
    if (!is.matrix(mat) || length(dm) != 2L || !identical(dim(mat), dm) || !nrow(pairs)) {
      return(rep(NA_real_, nrow(pairs)))
    }
    out <- rep(NA_real_, nrow(pairs))
    out[idx_ok] <- mat[cbind(ii[idx_ok], jj[idx_ok])]
    out
  }

  diag_attr <- attr(object, "diagnostics", exact = TRUE)
  if (is.list(diag_attr) && is.matrix(diag_attr$n_complete)) {
    pairs$n_complete <- as.integer(round(add_from_matrix(diag_attr$n_complete)))
  }

  ci_attr <- attr(object, "ci", exact = TRUE)
  if (is.list(ci_attr)) {
    if (is.matrix(ci_attr$lwr.ci)) pairs$lwr <- as.numeric(add_from_matrix(ci_attr$lwr.ci))
    if (is.matrix(ci_attr$upr.ci)) pairs$upr <- as.numeric(add_from_matrix(ci_attr$upr.ci))
  }

  inf <- attr(object, "inference", exact = TRUE)
  if (is.list(inf)) {
    if (is.matrix(inf$statistic)) pairs$statistic <- as.numeric(add_from_matrix(inf$statistic))
    param_mat <- inf$parameter %||% inf$df
    if (is.matrix(param_mat)) pairs$df <- as.numeric(add_from_matrix(param_mat))
    if (is.matrix(inf$p_value)) pairs$p_value <- as.numeric(add_from_matrix(inf$p_value))
  } else {
    needs_fallback <- estimator_class %in% c(
      "bicor", "dcor", "pearson_corr", "spearman_rho",
      "kendall_matrix", "pbcor", "wincor", "skipped_corr",
      "partial_corr_matrix"
    )
    if (isTRUE(needs_fallback)) {
      # Backward-compatible fallback for CI-only summaries that historically
      # exposed Fisher-z and p-value columns.
      z <- suppressWarnings(atanh(pairs$estimate))
      z[!is.finite(z)] <- NA_real_
      pairs$fisher_z <- z
      if ("n_complete" %in% names(pairs) && all(is.finite(pairs$n_complete))) {
        se <- 1 / sqrt(pmax(1, as.numeric(pairs$n_complete) - 3))
        pairs$statistic <- z / se
        pairs$p_value <- 2 * stats::pnorm(abs(pairs$statistic), lower.tail = FALSE)
      }
    }
  }

  if (!"fisher_z" %in% names(pairs)) {
    z <- suppressWarnings(atanh(pairs$estimate))
    z[!is.finite(z)] <- NA_real_
    pairs$fisher_z <- z
  }

  if (nrow(pairs)) {
    pairs <- pairs[order(abs(pairs$estimate), decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  }
  rownames(pairs) <- NULL

  overview <- .mc_corr_overview(object)
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
  class(pairs) <- unique(summary_classes)

  topn_use <- .mc_coalesce(topn, .mc_display_option("summary_topn", 5L))
  topn_use <- .mc_validate_optional_count(topn_use, arg = "topn", allow_null = FALSE)
  top_tbl <- if (nrow(pairs)) {
    pairs[utils::head(order(abs(pairs$estimate), decreasing = TRUE, na.last = TRUE), topn_use), , drop = FALSE]
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
  attr(pairs, "ci_method") <- attr(object, "ci", exact = TRUE)$ci.method %||%
    attr(object, "ci.method", exact = TRUE) %||%
    NA_character_
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
    spearman_rho = if (isTRUE(has_ci) || isTRUE(has_p)) "Spearman correlation summary" else "Correlation summary",
    kendall_matrix = "Kendall correlation summary",
    bicor = "Biweight mid-correlation summary",
    dcor = "Distance correlation summary",
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
.mc_corr_plot_grid <- function(x) {
  dm <- .mc_corr_dim(x)
  dn <- .mc_corr_dimnames(x)
  edf <- .mc_corr_as_edge_df(x)
  if (!nrow(edf)) {
    return(data.frame(
      row = character(),
      col = character(),
      value = numeric(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  symm <- isTRUE(attr(x, "corr_symmetric", exact = TRUE))
  if (symm) {
    off <- edf$row != edf$col
    if (any(off)) {
      edf <- rbind(
        edf,
        data.frame(
          row = edf$col[off],
          col = edf$row[off],
          value = edf$value[off],
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      )
    }
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
                                 show_value,
                                 ...) {
  check_bool(show_value, arg = "show_value")
  df <- .mc_corr_plot_grid(x)
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$col, .data$row, fill = .data$value)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(
      low = low_color,
      high = high_color,
      mid = mid_color,
      midpoint = 0,
      limits = c(-1, 1),
      name = "Correlation"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = NULL, y = NULL)

  if (isTRUE(show_value) && is.finite(value_text_size) && !is.null(value_text_size)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", .data$value)),
      size = value_text_size,
      color = "black"
    )
  }
  p
}

