#' Internal display helpers
#'
#' Shared console preview and summary formatting for matrixCorr S3 methods.
#'
#' @keywords internal
#' @noRd
NULL

.mc_coalesce <- function(x, y) if (is.null(x)) y else x

.mc_display_option <- function(name, default) {
  getOption(paste0("matrixCorr.", name), default = default)
}

.mc_validate_yes_no <- function(x,
                                arg = as.character(substitute(x)),
                                default = NULL) {
  if (is.null(x)) {
    x <- default
  }
  if (!is.character(x) || length(x) != 1L || is.na(x) || !x %in% c("yes", "no")) {
    abort_bad_arg(arg, message = "must be one of \"yes\" or \"no\".")
  }
  x
}

.mc_resolve_show_ci <- function(show_ci,
                                context = c("print", "summary")) {
  context <- match.arg(context)
  .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option(paste0(context, "_show_ci"), "yes")
  )
}

.mc_extract_legacy_display_args <- function(dots,
                                            n = NULL,
                                            topn = NULL,
                                            max_vars = NULL) {
  dots <- if (is.null(dots)) list() else dots
  if (is.null(n) && !is.null(dots$max_rows)) {
    n <- dots$max_rows
    if (is.null(topn) && is.numeric(n) && is.finite(n)) {
      topn <- max(1L, floor(as.integer(n) / 2L))
    }
  }
  if (is.null(max_vars) && !is.null(dots$max_cols)) {
    max_vars <- dots$max_cols
  }
  dots$max_rows <- NULL
  dots$max_cols <- NULL
  list(dots = dots, n = n, topn = topn, max_vars = max_vars)
}

.mc_validate_optional_count <- function(x,
                                        arg = as.character(substitute(x)),
                                        allow_null = TRUE) {
  if (is.null(x)) {
    if (allow_null) return(NULL)
    abort_bad_arg(arg, message = "must not be NULL.")
  }
  if (!rlang::is_scalar_integerish(x) || is.na(x) || x <= 0) {
    abort_bad_arg(arg, message = "must be a positive integer.")
  }
  as.integer(x)
}

.mc_resolve_display_args <- function(context = c("print", "summary"),
                                     n = NULL,
                                     topn = NULL,
                                     max_vars = NULL,
                                     width = NULL,
                                     show_ci = NULL) {
  context <- match.arg(context)
  topn_supplied <- !is.null(topn)
  defaults <- if (identical(context, "print")) {
    list(
      n = .mc_display_option("print_max_rows", 20L),
      topn = .mc_display_option("print_topn", 5L),
      max_vars = .mc_display_option("print_max_vars", NULL),
      show_ci = .mc_display_option("print_show_ci", "yes")
    )
  } else {
    list(
      n = .mc_display_option("summary_max_rows", 12L),
      topn = .mc_display_option("summary_topn", 5L),
      max_vars = .mc_display_option("summary_max_vars", 10L),
      show_ci = .mc_display_option("summary_show_ci", "yes")
    )
  }

  n <- .mc_validate_optional_count(.mc_coalesce(n, defaults$n), arg = "n", allow_null = FALSE)
  topn <- .mc_validate_optional_count(.mc_coalesce(topn, defaults$topn), arg = "topn", allow_null = FALSE)
  max_vars <- .mc_validate_optional_count(.mc_coalesce(max_vars, defaults$max_vars), arg = "max_vars")
  width <- .mc_validate_optional_count(.mc_coalesce(width, getOption("width", 80L)), arg = "width", allow_null = FALSE)
  show_ci <- .mc_validate_yes_no(show_ci, arg = "show_ci", default = defaults$show_ci)

  if (n <= 1L) {
    abort_bad_arg("n", message = "must be greater than 1.")
  }

  if (topn * 2L > n) {
    if (isTRUE(topn_supplied)) {
      cli::cli_warn(
        c(
          "Display truncation parameters are inconsistent.",
          "i" = "{.arg topn} is larger than half of {.arg n}; reducing the head/tail preview size."
        ),
        class = c("matrixCorr_warning", "matrixCorr_display_warning")
      )
    }
    topn <- max(1L, floor(n / 2L))
  }

  list(
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci
  )
}

.mc_count_fmt <- function(x) {
  format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

.mc_pick_head_tail <- function(total, threshold, topn) {
  if (total <= threshold) {
    return(list(index = seq_len(total), omitted = 0L, truncated = FALSE))
  }
  keep <- unique(c(seq_len(topn), seq.int(total - topn + 1L, total)))
  list(index = keep, omitted = total - length(keep), truncated = TRUE)
}

.mc_width_to_max_vars <- function(labels,
                                  width,
                                  min_visible = 2L) {
  labels <- as.character(labels)
  if (!length(labels)) return(0L)
  name_width <- max(nchar(labels, type = "width"), na.rm = TRUE)
  approx_cell <- max(8L, name_width + 1L)
  available <- max(min_visible, floor((width - 12L) / approx_cell))
  min(length(labels), available)
}

.mc_pick_leading_trailing <- function(total, max_visible) {
  if (total <= max_visible) {
    return(list(index = seq_len(total), omitted = 0L, truncated = FALSE))
  }
  left <- ceiling(max_visible / 2)
  right <- floor(max_visible / 2)
  keep <- unique(c(seq_len(left), seq.int(total - right + 1L, total)))
  list(index = keep, omitted = total - length(keep), truncated = TRUE)
}

.mc_preview_columns <- function(labels, width, max_vars = NULL) {
  if (!length(labels)) return(list(index = integer(), omitted = 0L, truncated = FALSE))
  visible <- if (is.null(max_vars)) {
    .mc_width_to_max_vars(labels = labels, width = width)
  } else {
    min(length(labels), max_vars)
  }
  visible <- max(2L, visible)
  .mc_pick_leading_trailing(length(labels), visible)
}

.mc_format_matrix_values <- function(m, digits = 4) {
  if (is.numeric(m)) {
    return(matrix(
      ifelse(is.na(m), NA_character_, formatC(m, format = "f", digits = digits)),
      nrow = nrow(m),
      ncol = ncol(m),
      dimnames = dimnames(m)
    ))
  }
  matrix(as.character(m), nrow = nrow(m), ncol = ncol(m), dimnames = dimnames(m))
}

.mc_build_matrix_preview <- function(m,
                                     digits = 4,
                                     n = 20L,
                                     topn = 5L,
                                     max_vars = NULL,
                                     width = getOption("width", 80L)) {
  row_info <- .mc_pick_head_tail(nrow(m), threshold = n, topn = topn)
  col_labels <- .mc_coalesce(colnames(m), as.character(seq_len(ncol(m))))
  col_info <- .mc_preview_columns(col_labels, width = width, max_vars = max_vars)

  fmt <- .mc_format_matrix_values(m, digits = digits)
  preview <- fmt[row_info$index, col_info$index, drop = FALSE]

  if (col_info$truncated) {
    left <- ceiling(length(col_info$index) / 2)
    preview <- cbind(
      preview[, seq_len(left), drop = FALSE],
      "..." = rep("...", nrow(preview)),
      preview[, seq.int(left + 1L, ncol(preview)), drop = FALSE]
    )
  }

  if (row_info$truncated) {
    top_keep <- topn
    preview <- rbind(
      preview[seq_len(top_keep), , drop = FALSE],
      matrix("...", nrow = 1L, ncol = ncol(preview),
             dimnames = list("...", colnames(preview))),
      preview[seq.int(nrow(preview) - topn + 1L, nrow(preview)), , drop = FALSE]
    )
  }

  list(
    data = as.data.frame(preview, stringsAsFactors = FALSE, check.names = FALSE),
    row_info = row_info,
    col_info = col_info
  )
}

.mc_build_table_preview <- function(df,
                                    n = 20L,
                                    topn = 5L,
                                    max_vars = NULL,
                                    width = getOption("width", 80L)) {
  row_info <- .mc_pick_head_tail(nrow(df), threshold = n, topn = topn)
  col_info <- .mc_preview_columns(names(df), width = width, max_vars = max_vars)
  preview <- df[row_info$index, col_info$index, drop = FALSE]

  if (col_info$truncated) {
    left <- ceiling(length(col_info$index) / 2)
    preview <- data.frame(
      preview[, seq_len(left), drop = FALSE],
      "..." = rep("...", nrow(preview)),
      preview[, seq.int(left + 1L, ncol(preview)), drop = FALSE],
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  if (row_info$truncated) {
    ellipsis_row <- as.list(rep("...", ncol(preview)))
    names(ellipsis_row) <- names(preview)
    preview <- rbind(
      preview[seq_len(topn), , drop = FALSE],
      as.data.frame(ellipsis_row, stringsAsFactors = FALSE, check.names = FALSE),
      preview[seq.int(nrow(preview) - topn + 1L, nrow(preview)), , drop = FALSE]
    )
    rownames(preview) <- NULL
  }

  list(
    data = preview,
    row_info = row_info,
    col_info = col_info
  )
}

.mc_display_footer <- function(row_omitted = 0L,
                               col_omitted = 0L,
                               context = c("print", "summary"),
                               full_hint = TRUE,
                               summary_hint = FALSE) {
  context <- match.arg(context)
  lines <- character()
  truncated <- (row_omitted > 0L) || (col_omitted > 0L)
  if (row_omitted > 0L) {
    lines <- c(lines, sprintf("... %s more rows not shown (omitted)", .mc_count_fmt(row_omitted)))
  }
  if (col_omitted > 0L) {
    lines <- c(lines, sprintf("... %s more variables not shown (omitted)", .mc_count_fmt(col_omitted)))
  }
  if (summary_hint && identical(context, "print") && truncated) {
    lines <- c(lines, "Use summary() for a richer digest.")
  }
  if (full_hint && truncated) {
    lines <- c(lines, "Use as.data.frame()/tidy()/as.matrix() to inspect the full result.")
  }
  lines
}

.mc_ci_column_index <- function(df) {
  if (!ncol(df)) return(integer())
  nm <- names(df)
  which(grepl("(^lwr$|^upr$|(^|[._])lwr$|(^|[._])upr$|^ci$|^ci_method$|^conf[._]level$)", nm))
}

.mc_header_with_ci <- function(base, conf_level = NA_real_, show_ci = "yes") {
  if (identical(show_ci, "yes") && is.finite(conf_level)) {
    sprintf("%s (%g%% CI)", base, 100 * conf_level)
  } else {
    base
  }
}

.mc_print_preview_table <- function(df,
                                    n = 20L,
                                    topn = 5L,
                                    max_vars = NULL,
                                    width = getOption("width", 80L),
                                    context = c("print", "summary"),
                                    full_hint = TRUE,
                                    summary_hint = FALSE,
                                    ...) {
  context <- match.arg(context)
  preview <- .mc_build_table_preview(df, n = n, topn = topn, max_vars = max_vars, width = width)
  print.data.frame(preview$data, row.names = FALSE, right = FALSE, ...)
  foot <- .mc_display_footer(
    row_omitted = preview$row_info$omitted,
    col_omitted = preview$col_info$omitted,
    context = context,
    full_hint = full_hint,
    summary_hint = summary_hint
  )
  if (length(foot)) cat(paste0(foot, collapse = "\n"), "\n", sep = "")
  invisible(preview)
}

.mc_print_preview_matrix <- function(m,
                                     digits = 4,
                                     n = 20L,
                                     topn = 5L,
                                     max_vars = NULL,
                                     width = getOption("width", 80L),
                                     context = c("print", "summary"),
                                     full_hint = TRUE,
                                     summary_hint = FALSE,
                                     ...) {
  context <- match.arg(context)
  preview <- .mc_build_matrix_preview(
    m,
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width
  )
  print.data.frame(preview$data, right = TRUE, ...)
  foot <- .mc_display_footer(
    row_omitted = preview$row_info$omitted,
    col_omitted = preview$col_info$omitted,
    context = context,
    full_hint = full_hint,
    summary_hint = summary_hint
  )
  if (length(foot)) cat(paste0(foot, collapse = "\n"), "\n", sep = "")
  invisible(preview)
}

.mc_line_item <- function(label, value) {
  cat(sprintf("  %-12s: %s\n", label, value))
}

.mc_print_named_digest <- function(items, header = NULL) {
  if (!is.null(header)) cat(header, "\n", sep = "")
  for (nm in names(items)) {
    .mc_line_item(nm, items[[nm]])
  }
  invisible(items)
}

.mc_format_scalar_or_range <- function(x_min,
                                       x_max = x_min,
                                       digits = 4,
                                       integer = FALSE) {
  if (!is.finite(x_min) || !is.finite(x_max)) {
    return(NULL)
  }
  fmt_one <- function(x) {
    if (isTRUE(integer)) {
      .mc_count_fmt(x)
    } else {
      formatC(x, format = "f", digits = digits)
    }
  }
  if (isTRUE(all.equal(x_min, x_max))) {
    return(fmt_one(x_min))
  }
  sprintf("%s to %s", fmt_one(x_min), fmt_one(x_max))
}

.mc_pair_label <- function(var1, var2) {
  sprintf("%s-%s", var1, var2)
}

.mc_top_pair_by <- function(df, which_idx, digits = 4) {
  if (!nrow(df) || is.na(which_idx) || which_idx < 1L || which_idx > nrow(df)) {
    return(NULL)
  }
  row <- df[which_idx, , drop = FALSE]
  if (!all(c("item1", "item2", "estimate") %in% names(row))) {
    return(NULL)
  }
  sprintf(
    "%s (%s)",
    .mc_pair_label(row$item1[[1L]], row$item2[[1L]]),
    formatC(row$estimate[[1L]], format = "f", digits = digits)
  )
}

.mc_pair_extremes <- function(object, digits = 4) {
  tbl <- .mc_format_pair_table(as.matrix(object), digits = digits, sort_abs = FALSE)
  if (!nrow(tbl)) {
    return(list(
      most_negative = NULL,
      most_positive = NULL
    ))
  }
  est <- tbl$estimate
  list(
    most_negative = .mc_top_pair_by(tbl, which.min(est), digits = digits),
    most_positive = .mc_top_pair_by(tbl, which.max(est), digits = digits)
  )
}

.mc_ci_digest <- function(df,
                          conf_level = NA_real_,
                          ci_method = NULL,
                          digits = 3) {
  if (!all(c("lwr", "upr") %in% names(df))) {
    return(setNames(character(), character()))
  }
  keep <- is.finite(df$lwr) & is.finite(df$upr)
  if (!any(keep)) {
    return(setNames(character(), character()))
  }
  widths <- df$upr[keep] - df$lwr[keep]
  out <- c(ci = if (is.finite(conf_level)) sprintf("%g%%", 100 * conf_level) else "available")
  if (is.character(ci_method) && length(ci_method) == 1L && !is.na(ci_method) && nzchar(ci_method)) {
    out <- c(out, ci_method = ci_method)
  }
  width_txt <- .mc_format_scalar_or_range(min(widths), max(widths), digits = digits, integer = FALSE)
  if (!is.null(width_txt)) {
    out <- c(out, ci_width = width_txt)
  }
  cross_zero <- sum(df$lwr[keep] <= 0 & df$upr[keep] >= 0, na.rm = TRUE)
  if (cross_zero > 0L) {
    out <- c(out, cross_zero = sprintf("%s pair(s)", .mc_count_fmt(cross_zero)))
  }
  out
}

.mc_pairwise_top_results <- function(df,
                                     topn = 5L,
                                     show_ci = "yes") {
  if (!nrow(df) || !"estimate" %in% names(df)) {
    return(df[0, , drop = FALSE])
  }
  topn <- min(.mc_validate_optional_count(topn, arg = "topn", allow_null = FALSE), nrow(df))
  ord <- order(abs(df$estimate), decreasing = TRUE, na.last = NA)
  out <- df[utils::head(ord, topn), , drop = FALSE]
  keep <- intersect(
    c(
      "item1", "item2", "estimate", "lwr", "upr", "statistic", "df",
      "p_value", "n_complete", "n_subjects", "n_obs",
      "p_value_adjusted", "reject", "skipped_n", "skipped_prop",
      "lwr", "upr", "p_value", "p_value_adjusted", "reject"
    ),
    names(out)
  )
  if (identical(show_ci, "no")) {
    keep <- setdiff(keep, c("lwr", "upr"))
  }
  out <- out[, keep, drop = FALSE]
  rownames(out) <- NULL
  out
}

.mc_standardize_summary_pairs <- function(df) {
  if (!is.data.frame(df)) return(df)
  nms <- names(df)
  if ("var1" %in% nms) nms[nms == "var1"] <- "item1"
  if ("var2" %in% nms) nms[nms == "var2"] <- "item2"
  if ("method1" %in% nms) nms[nms == "method1"] <- "item1"
  if ("method2" %in% nms) nms[nms == "method2"] <- "item2"
  names(df) <- nms
  df
}

.mc_standardize_summary_counts <- function(df,
                                           repeated = FALSE,
                                           rename_n = FALSE) {
  if (!is.data.frame(df)) return(df)
  if ("based.on" %in% names(df) && !"n_obs" %in% names(df)) {
    names(df)[names(df) == "based.on"] <- "n_obs"
  }
  if (isTRUE(rename_n) && "n" %in% names(df) && !"n_obs" %in% names(df)) {
    names(df)[names(df) == "n"] <- "n_obs"
  }
  if (isTRUE(repeated) && "n_complete" %in% names(df) && !"n_obs" %in% names(df)) {
    names(df)[names(df) == "n_complete"] <- "n_obs"
  }
  df
}

.mc_reorder_summary_columns <- function(df,
                                        repeated = FALSE) {
  if (!is.data.frame(df)) return(df)
  core <- c("item1", "item2", "estimate", "lwr", "upr", "statistic", "df", "p_value")
  counts <- if (isTRUE(repeated)) c("n_subjects", "n_obs") else "n_complete"
  keep <- c(core, counts)
  ordered <- c(intersect(keep, names(df)), setdiff(names(df), keep))
  df[, ordered, drop = FALSE]
}

.mc_finalize_summary_df <- function(df,
                                    class_name,
                                    repeated = FALSE,
                                    rename_n = FALSE) {
  extra_attrs <- attributes(df)
  extra_attrs[c("names", "row.names", "class")] <- NULL
  df <- .mc_standardize_summary_pairs(df)
  df <- .mc_standardize_summary_counts(df, repeated = repeated, rename_n = rename_n)
  df <- .mc_reorder_summary_columns(df, repeated = repeated)
  attributes(df) <- c(
    attributes(df),
    extra_attrs,
    list(class = c(class_name, "summary.matrixCorr", "data.frame"))
  )
  df
}

.mc_finalize_summary_list <- function(x,
                                      class_name) {
  class(x) <- c(class_name, "summary.matrixCorr")
  x
}

.mc_print_ranked_pairs_preview <- function(df,
                                           header = "Strongest pairs by |estimate|",
                                           topn = 5L,
                                           max_vars = NULL,
                                           width = getOption("width", 80L),
                                           show_ci = "yes",
                                           ...) {
  top <- .mc_pairwise_top_results(df, topn = topn, show_ci = show_ci)
  if (!nrow(top)) {
    return(invisible(top))
  }
  cat("\n", header, "\n\n", sep = "")
  .mc_print_preview_table(
    top,
    n = max(2L, nrow(top) + 1L),
    topn = max(1L, min(topn, nrow(top))),
    max_vars = max_vars,
    width = width,
    context = "summary",
    full_hint = FALSE,
    summary_hint = FALSE,
    ...
  )
  foot <- .mc_display_footer(
    row_omitted = max(0L, nrow(df) - nrow(top)),
    col_omitted = 0L,
    context = "summary",
    full_hint = nrow(df) > nrow(top),
    summary_hint = FALSE
  )
  if (length(foot)) cat(paste0(foot, collapse = "\n"), "\n", sep = "")
  invisible(top)
}

.mc_corr_summary_digest_items <- function(x, digits = 4, show_ci = "yes") {
  items <- c(
    method = .mc_coalesce(x$method, NA_character_),
    dimensions = sprintf("%d x %d", x$n_rows, x$n_cols),
    pairs = .mc_count_fmt(x$n_pairs)
  )
  items <- items[!is.na(items) & nzchar(items)]

  if (isTRUE(x$has_ci) && identical(show_ci, "yes")) {
    items <- c(items, ci = "yes")
  }
  if (isTRUE(x$has_p)) {
    items <- c(items, p_values = "yes")
  }
  if (isTRUE(x$n_missing > 0L)) {
    items <- c(items, missing = .mc_count_fmt(x$n_missing))
  }
  n_complete_txt <- .mc_format_scalar_or_range(x$n_complete_min, x$n_complete_max, digits = digits, integer = TRUE)
  if (!is.null(n_complete_txt)) {
    items <- c(items, n_complete = n_complete_txt)
  }
  skipped_n_txt <- .mc_format_scalar_or_range(x$skipped_n_min, x$skipped_n_max, digits = digits, integer = TRUE)
  if (!is.null(skipped_n_txt) && !(identical(skipped_n_txt, "0"))) {
    items <- c(items, skipped_n = skipped_n_txt)
  }
  skipped_prop_txt <- .mc_format_scalar_or_range(x$skipped_prop_min, x$skipped_prop_max, digits = digits, integer = FALSE)
  if (!is.null(skipped_prop_txt) && !(identical(skipped_prop_txt, formatC(0, format = "f", digits = digits)))) {
    items <- c(items, skipped_prop = skipped_prop_txt)
  }
  if (!is.null(x$correct) &&
      length(x$correct) == 1L &&
      !is.na(x$correct) &&
      is.finite(x$correct) &&
      !isTRUE(all.equal(x$correct, 0))) {
    items <- c(items, correct = formatC(x$correct, format = "f", digits = digits))
  }
  if (is.finite(x$threshold_sets) && x$threshold_sets > 0L) {
    items <- c(items, thresholds = sprintf("%s set(s)", .mc_count_fmt(x$threshold_sets)))
  }
  if (is.finite(x$zero_cell_pairs) && x$zero_cell_pairs > 0L) {
    items <- c(items, zero_cells = sprintf("%s pair(s)", .mc_count_fmt(x$zero_cell_pairs)))
  }
  if (is.finite(x$corrected_pairs) && x$corrected_pairs > 0L) {
    items <- c(items, corrected = sprintf("%s pair(s)", .mc_count_fmt(x$corrected_pairs)))
  }
  if (is.finite(x$boundary_pairs) && x$boundary_pairs > 0L) {
    items <- c(items, boundary = sprintf("%s pair(s)", .mc_count_fmt(x$boundary_pairs)))
  }
  if (is.finite(x$near_boundary_pairs) && x$near_boundary_pairs > 0L) {
    items <- c(items, near_bound = sprintf("%s pair(s)", .mc_count_fmt(x$near_boundary_pairs)))
  }
  if (is.finite(x$converged_pairs) && x$converged_pairs > 0L) {
    items <- c(items, converged = sprintf("%s pair(s)", .mc_count_fmt(x$converged_pairs)))
  }
  if (is.finite(x$optimizer_tol)) {
    items <- c(items, opt_tol = format(signif(x$optimizer_tol, digits = digits)))
  }
  est_range <- .mc_format_scalar_or_range(x$estimate_min, x$estimate_max, digits = digits, integer = FALSE)
  if (!is.null(est_range)) {
    items <- c(items, estimate = est_range)
  }
  if (!is.null(x$most_negative) && nzchar(x$most_negative)) {
    items <- c(items, most_negative = x$most_negative)
  }
  if (!is.null(x$most_positive) && nzchar(x$most_positive)) {
    items <- c(items, most_positive = x$most_positive)
  }
  items
}

.mc_print_pairwise_summary_digest <- function(x,
                                              title,
                                              digits = 4,
                                              n = NULL,
                                              topn = NULL,
                                              max_vars = NULL,
                                              width = NULL,
                                              show_ci = NULL,
                                              ci_method = NULL,
                                              extra_items = NULL,
                                              ...) {
  show_ci <- .mc_resolve_show_ci(show_ci, context = "summary")
  cfg <- .mc_resolve_display_args(
    context = "summary",
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci
  )
  digits <- .mc_coalesce(attr(x, "digits"), digits)
  overview <- attr(x, "overview", exact = TRUE)

  digest <- if (inherits(overview, "summary.matrixCorr")) {
    .mc_corr_summary_digest_items(overview, digits = digits, show_ci = cfg$show_ci)
  } else {
    setNames(character(), character())
  }
  if (identical(cfg$show_ci, "yes") && isTRUE(attr(x, "has_ci"))) {
    digest <- c(
      digest,
      .mc_ci_digest(
        as.data.frame(x, stringsAsFactors = FALSE),
        conf_level = suppressWarnings(as.numeric(attr(x, "conf.level"))),
        ci_method = ci_method,
        digits = .mc_coalesce(attr(x, "ci_digits"), 3)
      )
    )
  }
  if (length(extra_items)) {
    digest <- c(digest, extra_items)
  }
  if (length(digest)) {
    digest <- digest[!duplicated(names(digest), fromLast = TRUE)]
  }
  digest <- digest[!is.na(digest) & nzchar(digest) & !(names(digest) == "multiplicity" & digest == "none")]

  .mc_print_named_digest(digest, header = title)
  .mc_print_ranked_pairs_preview(
    .mc_reorder_summary_columns(
      .mc_standardize_summary_counts(
        .mc_standardize_summary_pairs(as.data.frame(x, stringsAsFactors = FALSE)),
        repeated = any(c("n_subjects", "n_obs") %in% names(x))
      ),
      repeated = any(c("n_subjects", "n_obs") %in% names(x))
    ),
    header = "Strongest pairs by |estimate|",
    topn = cfg$topn,
    max_vars = cfg$max_vars,
    width = cfg$width,
    show_ci = cfg$show_ci,
    ...
  )
  invisible(x)
}

.mc_print_overview_if_present <- function(x, ...) {
  overview <- attr(x, "overview", exact = TRUE)
  show_overview <- attr(x, "show_overview", exact = TRUE)
  if (isTRUE(show_overview) || is.null(show_overview)) {
    if (!is.null(overview)) {
      class(overview) <- unique(c(class(overview), "summary.matrixCorr"))
      print.summary.matrixCorr(overview, ...)
    }
  }
  invisible(NULL)
}

.mc_print_summary_table <- function(x,
                                    header,
                                    digits = NULL,
                                    n = NULL,
                                    topn = NULL,
                                    max_vars = NULL,
                                    width = NULL,
                                    show_ci = NULL,
                                    print_overview = TRUE,
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

  if (isTRUE(print_overview)) {
    do.call(
      .mc_print_overview_if_present,
      c(
        list(
          x = x,
          digits = if (is.null(digits)) 4 else digits,
          n = cfg$n,
          topn = cfg$topn,
          max_vars = cfg$max_vars,
          width = cfg$width,
          show_ci = cfg$show_ci
        ),
        legacy$dots
      )
    )
  }

  cat("\n", header, "\n\n", sep = "")
  table_x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (identical(cfg$show_ci, "no")) {
    ci_cols <- .mc_ci_column_index(table_x)
    if (length(ci_cols)) {
      table_x <- table_x[, -ci_cols, drop = FALSE]
    }
  }
  do.call(
    .mc_print_preview_table,
    c(
      list(
        df = table_x,
        n = cfg$n,
        topn = cfg$topn,
        max_vars = cfg$max_vars,
        width = cfg$width,
        context = "summary",
        full_hint = TRUE,
        summary_hint = FALSE
      ),
      legacy$dots
    )
  )
  invisible(x)
}

.mc_print_sectioned_table <- function(x,
                                      sections,
                                      header,
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
  table_x <- as.data.frame(x, stringsAsFactors = FALSE)
  ci_cols <- .mc_ci_column_index(table_x)

  if (!is.null(header)) {
    cat("\n", header, "\n\n", sep = "")
  }

  printed <- 0L
  truncated_any <- FALSE
  for (section in sections) {
    cols <- intersect(section$cols, names(table_x))
    if (identical(cfg$show_ci, "no") && length(ci_cols)) {
      cols <- setdiff(cols, names(table_x)[ci_cols])
    }
    if (!length(cols)) next
    if (printed > 0L) cat("\n")
    cat(section$title, "\n\n", sep = "")
    preview <- do.call(
      .mc_print_preview_table,
      c(
        list(
          df = table_x[, cols, drop = FALSE],
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
    truncated_any <- truncated_any ||
      isTRUE(preview$row_info$omitted > 0L) ||
      isTRUE(preview$col_info$omitted > 0L)
    printed <- printed + 1L
  }
  if (printed > 0L && truncated_any) {
    cat("Use as.data.frame()/tidy()/as.matrix() to inspect the full result.\n")
  }
  invisible(x)
}

.mc_format_pair_table <- function(m,
                                  digits = 4,
                                  sort_abs = TRUE) {
  m <- as.matrix(m)
  rn <- .mc_coalesce(rownames(m), as.character(seq_len(nrow(m))))
  cn <- .mc_coalesce(colnames(m), as.character(seq_len(ncol(m))))
  symmetric <- isTRUE(nrow(m) == ncol(m)) && isTRUE(isSymmetric(unclass(m)))
  selector <- if (symmetric && nrow(m) > 1L) {
    upper.tri(m, diag = FALSE)
  } else {
    matrix(TRUE, nrow(m), ncol(m))
  }
  idx <- which(selector, arr.ind = TRUE)
  out <- data.frame(
    item1 = rn[idx[, 1L]],
    item2 = cn[idx[, 2L]],
    estimate = as.numeric(m[idx]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (sort_abs && nrow(out)) {
    out <- out[order(abs(out$estimate), decreasing = TRUE), , drop = FALSE]
    rownames(out) <- NULL
  }
  out$estimate <- round(out$estimate, digits)
  out
}

.mc_summary_top_pairs <- function(object,
                                  digits = 4,
                                  topn = 5L) {
  topn <- .mc_validate_optional_count(topn, arg = "topn", allow_null = FALSE)
  tbl <- .mc_format_pair_table(as.matrix(object), digits = digits, sort_abs = TRUE)
  utils::head(tbl, topn)
}

.mc_has_ci_payload <- function(x) {
  ci <- attr(x, "ci", exact = TRUE)
  if (is.list(ci)) return(TRUE)
  is.list(x) && all(c("est", "lwr.ci", "upr.ci") %in% names(x))
}

.mc_has_p_payload <- function(x) {
  inf <- attr(x, "inference", exact = TRUE)
  if (is.list(inf) && !is.null(inf$p_value)) return(TRUE)
  diag_attr <- attr(x, "diagnostics", exact = TRUE)
  is.list(diag_attr) && any(c("p_value", "p_value_adjusted") %in% names(diag_attr))
}
