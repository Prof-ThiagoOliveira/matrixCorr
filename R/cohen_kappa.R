#' @title Pairwise Cohen's Kappa for Nominal Ratings
#'
#' @description
#' Computes unweighted Cohen's kappa for either a pair of nominal rating
#' vectors or all pairwise combinations of nominal columns in a matrix or data
#' frame.
#'
#' @param data In matrix mode, a matrix or data frame whose rows are
#'   observational units and whose columns are raters or classifiers. Supported
#'   column types are factor, ordered factor, character, logical, integer, and
#'   numeric, all treated as nominal categories here. If the ratings are truly
#'   ordinal and disagreements should be weighted by distance, use
#'   [weighted_kappa()] instead. In two-vector mode, the first nominal rating
#'   vector.
#' @param y Optional second nominal rating vector. When supplied, the function
#'   returns a single Cohen's kappa estimate for \code{data} and \code{y}.
#' @param na_method Character scalar controlling missing-data handling.
#'   \code{"error"} rejects missing values. \code{"complete"} applies listwise
#'   deletion across retained columns before computing the matrix. In matrix
#'   mode, \code{"pairwise"} uses pair-specific complete observations.
#' @param ci Logical; if \code{TRUE}, attach large-sample delta-method
#'   confidence intervals.
#' @param p_value Logical; if \code{TRUE}, attach large-sample delta-method
#'   Wald test statistics and p-values.
#' @param conf_level Confidence level used when \code{ci = TRUE}. Default is
#'   \code{0.95}.
#' @param n_threads Integer \eqn{\ge 1}. Number of OpenMP threads used by the
#'   matrix C++ backend. Defaults to
#'   \code{getOption("matrixCorr.threads", 1L)}.
#' @param output Output representation for matrix mode.
#'   \itemize{
#'   \item \code{"matrix"} (default): full dense matrix.
#'   \item \code{"sparse"}: sparse matrix from \pkg{Matrix} containing only
#'   retained entries after thresholding.
#'   \item \code{"edge_list"}: long-form data frame representation.
#'   }
#' @param threshold Non-negative absolute-value filter for non-matrix outputs:
#'   retain entries with \code{abs(value) >= threshold}. Must be \code{0} when
#'   \code{output = "matrix"}.
#' @param diag Logical; whether to include diagonal entries in sparse and
#'   edge-list outputs.
#' @param ... Reserved for future extensions. Unsupported extra arguments are
#'   rejected.
#'
#' @details
#' Cohen's kappa is an agreement coefficient for two raters assigning the same
#' units to mutually exclusive nominal categories. For contingency-table cell
#' proportions \eqn{p_{ij}},
#' \deqn{
#' \kappa = \frac{p_o - p_e}{1 - p_e},
#' }
#' where \eqn{p_o = \sum_i p_{ii}} is the observed agreement and
#' \eqn{p_e = \sum_i p_{i+} p_{+i}} is the chance agreement implied by the
#' marginal category proportions.
#'
#' This implementation is strictly the original unweighted nominal-scale Cohen's
#' kappa. If the categories are ordinal and near disagreements should count as
#' less severe than distant disagreements, use [weighted_kappa()] instead.
#'
#' In matrix mode, columns are treated as raters/classifiers and rows as shared
#' observational units. All pairwise Cohen's kappas between columns are
#' computed. Category labels are encoded in R before dispatch to C++, using a
#' common label mapping across all columns so that identical labels in different
#' columns correspond to the same category code.
#'
#' Missing values are encoded as \code{NA_integer_} before entering the C++
#' backend. With \code{na_method = "error"}, missing values are rejected before
#' computation. With \code{na_method = "complete"}, listwise deletion is
#' applied across retained columns. With \code{na_method = "pairwise"}, each
#' pair uses its own complete observations. Pairwise complete counts are stored
#' in \code{attr(x, "diagnostics")$n_complete}.
#'
#' \strong{Confidence intervals and standard errors.}
#' The implementation uses the exact large-sample formula. Let \eqn{p_{ab}}
#' be the empirical cell proportions,
#' \eqn{q_a = \sum_b p_{ab}} the row margins, \eqn{r_b = \sum_a p_{ab}} the
#' column margins, and \eqn{D = 1 - p_e}. For each cell \eqn{(a,b)}, define
#' the influence contribution
#' \deqn{ g_{ab}
#'       \;=\;
#'       \frac{\mathbf{1}(a=b)\,D + (p_o - 1)\,(r_a + q_b)}{D^2}. }
#' The variance estimator used by the code is
#' \deqn{ \widehat{\mathrm{Var}}(\hat\kappa)
#'       \;=\;
#'       \frac{1}{n}\left(
#'       \sum_{a,b} p_{ab} g_{ab}^2 -
#'       \Big(\sum_{a,b} p_{ab} g_{ab}\Big)^2
#'       \right), }
#' with \eqn{n} the number of complete paired ratings and
#' \eqn{\widehat{\mathrm{se}}(\hat\kappa)=\sqrt{\widehat{\mathrm{Var}}(\hat\kappa)}}.
#' The confidence interval is the Wald interval
#' \deqn{ \hat\kappa \pm z_{1-\alpha/2}\,\widehat{\mathrm{se}}(\hat\kappa), }
#' truncated to \eqn{[-1, 1]} in the returned result.
#'
#' @return
#' If \code{y} is supplied, a scalar S3 object of class \code{"cohen_kappa"}
#' backed by a numeric value, with attributes \code{diagnostics}, and
#' optionally \code{ci}, \code{inference}, and \code{conf.level}. Otherwise a
#' symmetric matrix-style result with estimator class \code{cohen_kappa}.
#'
#' @references
#' Cohen, J. (1960). A coefficient of agreement for nominal scales.
#' \emph{Educational and Psychological Measurement}, 20(1), 37-46.
#' \doi{10.1177/001316446002000104}
#'
#' @seealso
#' [weighted_kappa()] for two-rater ordered-category agreement with
#' distance-sensitive disagreement weights; [multirater_kappa()] for
#' panel-level nominal agreement among three or more raters.
#'
#' @examples
#' x <- factor(c("A", "A", "B", "B", "A", "B"))
#' y <- factor(c("A", "B", "B", "B", "A", "A"))
#' cohen_kappa(x, y)
#'
#' raters <- data.frame(
#'   r1 = factor(c("low", "low", "high", "high", "mid")),
#'   r2 = factor(c("low", "mid", "high", "high", "mid")),
#'   r3 = c("low", "low", "high", "mid", "mid")
#' )
#' ck <- cohen_kappa(raters)
#' print(ck)
#' summary(ck)
#' plot(ck)
#'
#' @author Thiago de Paula Oliveira
#' @export
cohen_kappa <- function(data,
                        y = NULL,
                        na_method = c("error", "pairwise", "complete"),
                        ci = FALSE,
                        p_value = FALSE,
                        conf_level = 0.95,
                        n_threads = getOption("matrixCorr.threads", 1L),
                        output = c("matrix", "sparse", "edge_list"),
                        threshold = 0,
                        diag = TRUE,
                        ...) {
  .mc_extract_legacy_aliases(list(...), allowed = character())
  output_cfg <- .mc_validate_thresholded_output_request(
    output = output,
    threshold = threshold,
    diag = diag
  )
  na_method <- validate_na_method(
    na_method,
    arg = "na_method",
    allowed = c("error", "pairwise", "complete")
  )
  check_bool(ci, arg = "ci")
  check_bool(p_value, arg = "p_value")
  if (isTRUE(ci)) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  }

  if (!is.null(y)) {
    if (!identical(output_cfg$output, "matrix")) {
      abort_bad_arg(
        "output",
        message = "must be {.val matrix} when {.arg y} is supplied."
      )
    }
    x_labels <- .mc_nominal_vector_labels(data, arg = "data")
    y_labels <- .mc_nominal_vector_labels(y, arg = "y")
    check_same_length(x_labels, y_labels, arg_x = "data", arg_y = "y")

    if (identical(na_method, "error") && (anyNA(x_labels) || anyNA(y_labels))) {
      abort_bad_arg(
        "data",
        message = "and {.arg y} must not contain missing values when {.arg na_method} is {.val error}.",
        .hint = "Use {.arg na_method} = {.val pairwise} or {.val complete} to drop incomplete pairs."
      )
    }

    if (!identical(na_method, "error")) {
      keep <- stats::complete.cases(x_labels, y_labels)
      x_labels <- x_labels[keep]
      y_labels <- y_labels[keep]
    }

    enc <- .mc_encode_nominal_pair_labels(x_labels, y_labels)
    fit <- cohen_kappa_pair_cpp(
      enc$x,
      enc$y,
      n_levels = enc$n_levels,
      return_inference = isTRUE(ci) || isTRUE(p_value),
      conf_level = conf_level
    )
    return(.mc_cohen_kappa_scalar(
      fit,
      ci = ci,
      p_value = p_value,
      conf_level = conf_level
    ))
  }

  enc <- .mc_extract_nominal_columns(data, arg = "data", min_cols = 2L)
  X <- enc$matrix
  col_names <- enc$names
  dn <- .mc_square_dimnames(col_names)
  n_levels <- enc$n_levels

  listwise_meta <- NULL
  if (identical(na_method, "error") && anyNA(X)) {
    abort_bad_arg(
      "data",
      message = "contains missing values when {.arg na_method} is {.val error}.",
      .hint = "Use {.arg na_method} = {.val pairwise} or {.val complete} to handle incomplete units."
    )
  }
  if (identical(na_method, "complete")) {
    keep <- stats::complete.cases(X)
    listwise_meta <- list(
      na_method = "complete",
      n_original = nrow(X),
      n_complete_rows = sum(keep),
      n_removed = nrow(X) - sum(keep),
      complete_rows = keep,
      common_sample = TRUE
    )
    X <- X[keep, , drop = FALSE]
  }

  direct_threshold <- !isTRUE(ci) &&
    !isTRUE(p_value) &&
    !identical(na_method, "pairwise") &&
    isTRUE(output_cfg$thresholded)

  prev_threads <- .mc_prepare_omp_threads(
    n_threads,
    n_threads_missing = missing(n_threads)
  )
  if (!is.null(prev_threads)) {
    on.exit(.mc_exit_omp_threads(prev_threads), add = TRUE)
  }

  if (direct_threshold) {
    trip <- cohen_kappa_threshold_triplets_cpp(
      X,
      n_levels = n_levels,
      threshold = output_cfg$threshold,
      diag = output_cfg$diag,
      n_threads = n_threads
    )
    n_complete_mat <- matrix(
      as.integer(nrow(X)),
      nrow = ncol(X),
      ncol = ncol(X),
      dimnames = dn
    )
    diag(n_complete_mat) <- colSums(!is.na(X))
    return(.mc_finalize_triplets_output(
      triplets = trip,
      output = output_cfg$output,
      estimator_class = "cohen_kappa",
      method = "cohen_kappa",
      description = "Pairwise Cohen's kappa agreement matrix",
      threshold = output_cfg$threshold,
      diag = output_cfg$diag,
      source_dim = as.integer(c(ncol(X), ncol(X))),
      source_dimnames = dn,
      diagnostics = .mc_merge_diagnostics(
        list(n_complete = n_complete_mat),
        listwise_meta
      ),
      symmetric = TRUE
    ))
  }

  fit <- cohen_kappa_matrix_cpp(
    X,
    n_levels = n_levels,
    pairwise_complete = identical(na_method, "pairwise"),
    return_inference = isTRUE(ci) || isTRUE(p_value),
    conf_level = conf_level,
    n_threads = n_threads
  )
  fit$est <- .mc_set_matrix_dimnames(fit$est, col_names)
  fit$n_complete <- .mc_set_matrix_dimnames(fit$n_complete, col_names)
  fit$po <- .mc_set_matrix_dimnames(fit$po, col_names)
  fit$pe <- .mc_set_matrix_dimnames(fit$pe, col_names)
  if (!is.null(fit$se)) {
    fit$se <- .mc_set_matrix_dimnames(fit$se, col_names)
    fit$lwr <- .mc_set_matrix_dimnames(fit$lwr, col_names)
    fit$upr <- .mc_set_matrix_dimnames(fit$upr, col_names)
    fit$statistic <- .mc_set_matrix_dimnames(fit$statistic, col_names)
    fit$p_value <- .mc_set_matrix_dimnames(fit$p_value, col_names)
  }

  diagnostics <- .mc_merge_diagnostics(
    list(
      n_complete = fit$n_complete,
      observed_agreement = fit$po,
      expected_agreement = fit$pe
    ),
    if (identical(na_method, "pairwise")) {
      list(
        na_method = "pairwise",
        common_sample = FALSE
      )
    } else {
      listwise_meta
    }
  )
  ci_attr <- if (isTRUE(ci)) {
    list(
      est = fit$est,
      lwr.ci = fit$lwr,
      upr.ci = fit$upr,
      conf.level = conf_level,
      ci.method = "delta"
    )
  } else {
    NULL
  }
  inference_attr <- if (isTRUE(p_value)) {
    list(
      method = "wald_z_delta_kappa",
      estimate = fit$est,
      se = fit$se,
      statistic = fit$statistic,
      p_value = fit$p_value,
      n_obs = fit$n_complete
    )
  } else {
    NULL
  }

  out <- .mc_new_corr_matrix(
    mat = fit$est,
    estimator_class = "cohen_kappa",
    method = "cohen_kappa",
    description = "Pairwise Cohen's kappa agreement matrix",
    diagnostics = diagnostics,
    ci = ci_attr,
    conf.level = if (isTRUE(ci)) conf_level else NULL,
    symmetric = TRUE,
    extra_attrs = list(inference = inference_attr)
  )

  .mc_finalize_corr_output_fast(
    out,
    output = output_cfg$output,
    threshold = output_cfg$threshold,
    diag = output_cfg$diag
  )
}

.mc_is_plain_nominal_atomic <- function(x) {
  cls <- class(x)
  identical(cls, "character") ||
    identical(cls, "logical") ||
    identical(cls, "integer") ||
    identical(cls, "numeric")
}

.mc_nominal_vector_labels <- function(x,
                                      arg = as.character(substitute(x))) {
  if (is.list(x) && !is.data.frame(x)) {
    abort_bad_arg(arg, message = "must not be a list.")
  }
  if (inherits(x, "Date") || inherits(x, "POSIXt") || inherits(x, "difftime")) {
    abort_bad_arg(
      arg,
      message = "must not be a Date/POSIXct/POSIXlt/difftime object."
    )
  }
  if (is.factor(x)) {
    out <- as.character(x)
    out[is.na(x)] <- NA_character_
    return(out)
  }
  if (.mc_is_plain_nominal_atomic(x)) {
    if (is.numeric(x) && any(!is.na(x) & !is.finite(x))) {
      abort_bad_arg(
        arg,
        message = "must not contain NaN or infinite values."
      )
    }
    out <- as.character(x)
    out[is.na(x)] <- NA_character_
    return(out)
  }
  abort_bad_arg(
    arg,
    message = paste(
      "must be a factor, ordered factor, character, logical, integer,",
      "or numeric nominal vector."
    )
  )
}

.mc_extract_nominal_columns <- function(data,
                                        arg = as.character(substitute(data)),
                                        min_cols = 2L) {
  cols <- NULL
  col_names <- NULL
  if (is.data.frame(data)) {
    cols <- unclass(data)
    col_names <- names(data)
  } else if (is.matrix(data)) {
    cols <- lapply(seq_len(ncol(data)), function(j) data[, j])
    col_names <- colnames(data)
  } else {
    abort_bad_arg(
      arg,
      message = "must be a matrix or data frame in matrix mode."
    )
  }

  if (length(cols) < min_cols) {
    abort_bad_arg(
      arg,
      message = "must contain at least {min_cols} nominal column{?s}.",
      min_cols = min_cols
    )
  }

  labels <- vector("list", length(cols))
  for (j in seq_along(cols)) {
    col_arg <- if (!is.null(col_names) && nzchar(col_names[[j]] %||% "")) {
      sprintf("%s$%s", arg, col_names[[j]])
    } else {
      sprintf("%s[[%d]]", arg, j)
    }
    labels[[j]] <- .mc_nominal_vector_labels(cols[[j]], arg = col_arg)
  }
  names(labels) <- col_names
  .mc_encode_nominal_label_list(labels)
}

.mc_encode_nominal_pair_labels <- function(x_labels, y_labels) {
  levels_all <- unique(c(
    x_labels[!is.na(x_labels)],
    y_labels[!is.na(y_labels)]
  ))
  list(
    x = .mc_match_nominal_levels(x_labels, levels_all),
    y = .mc_match_nominal_levels(y_labels, levels_all),
    levels = levels_all,
    n_levels = length(levels_all)
  )
}

.mc_encode_nominal_label_list <- function(labels) {
  levels_all <- unique(unlist(
    lapply(labels, function(z) z[!is.na(z)]),
    use.names = FALSE
  ))
  X <- matrix(
    NA_integer_,
    nrow = length(labels[[1L]]),
    ncol = length(labels),
    dimnames = list(NULL, names(labels))
  )
  for (j in seq_along(labels)) {
    X[, j] <- .mc_match_nominal_levels(labels[[j]], levels_all)
  }
  list(
    matrix = X,
    levels = levels_all,
    n_levels = length(levels_all),
    names = names(labels)
  )
}

.mc_match_nominal_levels <- function(labels, levels_all) {
  out <- match(labels, levels_all)
  out[is.na(labels)] <- NA_integer_
  as.integer(out)
}

.mc_cohen_kappa_scalar <- function(fit,
                                   ci = FALSE,
                                   p_value = FALSE,
                                   conf_level = 0.95) {
  est <- as.numeric(fit$est)
  out <- structure(est, class = c("cohen_kappa", "numeric"))
  attr(out, "method") <- "cohen_kappa"
  attr(out, "description") <- "Pairwise Cohen's kappa agreement"
  attr(out, "package") <- "matrixCorr"
  attr(out, "estimate") <- est
  attr(out, "diagnostics") <- list(
    n_complete = as.integer(fit$n_complete),
    observed_agreement = as.numeric(fit$po),
    expected_agreement = as.numeric(fit$pe)
  )
  if (isTRUE(ci)) {
    attr(out, "ci") <- list(
      est = est,
      lwr.ci = as.numeric(fit$lwr),
      upr.ci = as.numeric(fit$upr),
      conf.level = conf_level,
      ci.method = "delta"
    )
    attr(out, "conf.level") <- conf_level
    attr(out, "conf_level") <- conf_level
  }
  if (isTRUE(p_value)) {
    attr(out, "inference") <- list(
      method = "wald_z_delta_kappa",
      estimate = est,
      se = as.numeric(fit$se),
      statistic = as.numeric(fit$statistic),
      p_value = as.numeric(fit$p_value),
      n_obs = as.integer(fit$n_complete)
    )
  }
  out
}

#' @rdname cohen_kappa
#' @method summary cohen_kappa
#' @param object A scalar or matrix-style \code{cohen_kappa} object.
#' @param digits Integer; number of decimal places for displayed values.
#' @param ci_digits Integer; number of decimal places for confidence limits.
#' @param p_digits Integer; number of decimal places for p-values.
#' @param ... Unused.
#' @export
summary.cohen_kappa <- function(object,
                                digits = 4,
                                ci_digits = 3,
                                p_digits = 4,
                                ...) {
  if (inherits(object, "corr_result")) {
    return(summary.corr_matrix(object, ...))
  }
  check_inherits(object, "cohen_kappa")

  diag_attr <- attr(object, "diagnostics", exact = TRUE)
  ci_attr <- attr(object, "ci", exact = TRUE)
  inf_attr <- attr(object, "inference", exact = TRUE)

  df <- data.frame(
    estimate = as.numeric(object),
    n_complete = as.integer(diag_attr$n_complete %||% NA_integer_),
    observed_agreement = as.numeric(diag_attr$observed_agreement %||% NA_real_),
    expected_agreement = as.numeric(diag_attr$expected_agreement %||% NA_real_),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (is.list(ci_attr)) {
    df$lwr <- as.numeric(ci_attr$lwr.ci %||% NA_real_)
    df$upr <- as.numeric(ci_attr$upr.ci %||% NA_real_)
  }
  if (is.list(inf_attr)) {
    df$se <- as.numeric(inf_attr$se %||% NA_real_)
    df$statistic <- as.numeric(inf_attr$statistic %||% NA_real_)
    df$p_value <- as.numeric(inf_attr$p_value %||% NA_real_)
  }

  out <- .mc_finalize_summary_df(df, class_name = "summary.cohen_kappa")
  attr(out, "digits") <- digits
  attr(out, "ci_digits") <- ci_digits
  attr(out, "p_digits") <- p_digits
  attr(out, "has_ci") <- is.list(ci_attr)
  attr(out, "has_p") <- is.list(inf_attr)
  attr(out, "conf.level") <- attr(object, "conf.level", exact = TRUE) %||%
    attr(object, "conf_level", exact = TRUE) %||%
    NA_real_
  attr(out, "ci_method") <- if (is.list(ci_attr)) ci_attr$ci.method %||% NA_character_ else NA_character_
  attr(out, "inference_method") <- if (is.list(inf_attr)) inf_attr$method %||% NA_character_ else NA_character_
  out
}

#' @rdname cohen_kappa
#' @method print summary.cohen_kappa
#' @param x A \code{summary.cohen_kappa} object.
#' @param n Optional preview row threshold.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns.
#' @param width Optional display width.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @export
print.summary.cohen_kappa <- function(x,
                                      digits = NULL,
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
  digits <- .mc_coalesce(digits, .mc_coalesce(attr(x, "digits", exact = TRUE), 4L))
  ci_digits <- .mc_coalesce(attr(x, "ci_digits", exact = TRUE), digits)
  p_digits <- .mc_coalesce(attr(x, "p_digits", exact = TRUE), digits)

  digest <- c(
    method = "cohen_kappa",
    estimate = formatC(as.numeric(x$estimate[[1L]]), format = "f", digits = digits),
    n_complete = .mc_count_fmt(x$n_complete[[1L]]),
    observed_agreement = formatC(as.numeric(x$observed_agreement[[1L]]), format = "f", digits = digits),
    expected_agreement = formatC(as.numeric(x$expected_agreement[[1L]]), format = "f", digits = digits)
  )
  if (isTRUE(attr(x, "has_ci", exact = TRUE))) {
    digest <- c(
      digest,
      ci = sprintf("%g%%", 100 * suppressWarnings(as.numeric(attr(x, "conf.level", exact = TRUE)))),
      ci_method = attr(x, "ci_method", exact = TRUE) %||% NA_character_
    )
  }
  if (isTRUE(attr(x, "has_p", exact = TRUE))) {
    digest <- c(
      digest,
      inference = attr(x, "inference_method", exact = TRUE) %||% NA_character_
    )
  }
  digest <- digest[!is.na(digest) & nzchar(digest)]
  .mc_print_named_digest(digest, header = "Cohen's kappa agreement summary")
  cat("\n")

  preview <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  if (identical(show_ci, "no")) {
    preview <- preview[, setdiff(names(preview), c("lwr", "upr")), drop = FALSE]
  }
  if ("estimate" %in% names(preview)) preview$estimate <- round(preview$estimate, digits)
  if ("observed_agreement" %in% names(preview)) preview$observed_agreement <- round(preview$observed_agreement, digits)
  if ("expected_agreement" %in% names(preview)) preview$expected_agreement <- round(preview$expected_agreement, digits)
  if ("lwr" %in% names(preview)) preview$lwr <- round(preview$lwr, ci_digits)
  if ("upr" %in% names(preview)) preview$upr <- round(preview$upr, ci_digits)
  if ("se" %in% names(preview)) preview$se <- round(preview$se, digits)
  if ("statistic" %in% names(preview)) preview$statistic <- round(preview$statistic, digits)
  if ("p_value" %in% names(preview)) preview$p_value <- round(preview$p_value, p_digits)

  cfg <- .mc_resolve_display_args(
    context = "summary",
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci
  )
  .mc_print_preview_table(
    preview,
    n = cfg$n,
    topn = cfg$topn,
    max_vars = cfg$max_vars,
    width = cfg$width,
    context = "summary",
    full_hint = FALSE,
    summary_hint = FALSE
  )
  invisible(x)
}

#' @rdname cohen_kappa
#' @method print cohen_kappa
#' @param x A scalar or matrix-style \code{cohen_kappa} object.
#' @param digits Integer; number of decimal places for displayed values.
#' @param n Optional preview row threshold.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns.
#' @param width Optional display width.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param ... Additional arguments passed to downstream print helpers.
#' @export
print.cohen_kappa <- function(x, digits = 4, n = NULL, topn = NULL,
                              max_vars = NULL, width = NULL,
                              show_ci = NULL, ...) {
  if (!inherits(x, "corr_result")) {
    return(invisible(print.summary.cohen_kappa(
      summary.cohen_kappa(x),
      digits = digits,
      n = n,
      topn = topn,
      max_vars = max_vars,
      width = width,
      show_ci = show_ci,
      ...
    )))
  }
  .mc_print_corr_matrix(
    x,
    header = "Cohen's kappa agreement matrix",
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

#' @rdname cohen_kappa
#' @method plot cohen_kappa
#' @param x A scalar or matrix-style \code{cohen_kappa} object.
#' @param title Optional plot title.
#' @param low_color Fill/color used for negative agreement.
#' @param high_color Fill/color used for positive agreement.
#' @param mid_color Unused placeholder for API consistency with matrix heatmaps.
#' @param value_text_size Text size for the estimate label.
#' @param ci_text_size Text size for the CI label.
#' @param show_value Logical; whether to print the estimate and CI labels.
#' @param ... Additional theme arguments.
#' @export
plot.cohen_kappa <- function(x,
                             title = NULL,
                             low_color = "indianred1",
                             high_color = "steelblue1",
                             mid_color = "white",
                             value_text_size = 4,
                             ci_text_size = 3,
                             show_value = TRUE,
                             ...) {
  if (inherits(x, "corr_result")) {
    return(.mc_plot_corr_result(
      x,
      title = title %||% "Cohen's kappa agreement heatmap",
      low_color = low_color,
      high_color = high_color,
      mid_color = mid_color,
      value_text_size = value_text_size,
      ci_text_size = ci_text_size,
      show_value = show_value,
      ...
    ))
  }
  check_inherits(x, "cohen_kappa")
  check_bool(show_value, arg = "show_value")
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} is required for plotting.")
  }

  ci_attr <- attr(x, "ci", exact = TRUE)
  est <- as.numeric(x)
  df <- data.frame(
    label = "Cohen's kappa",
    estimate = est,
    lwr = if (is.list(ci_attr)) as.numeric(ci_attr$lwr.ci %||% NA_real_) else NA_real_,
    upr = if (is.list(ci_attr)) as.numeric(ci_attr$upr.ci %||% NA_real_) else NA_real_,
    stringsAsFactors = FALSE
  )
  df$label_y <- 1
  df$estimate_label <- sprintf("%.2f", df$estimate)
  df$ci_label <- ifelse(
    is.finite(df$lwr) & is.finite(df$upr),
    sprintf("[%.2f, %.2f]", df$lwr, df$upr),
    NA_character_
  )
  point_color <- if (is.finite(est) && est < 0) low_color else high_color

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$estimate, y = .data$label_y)) +
    ggplot2::geom_vline(xintercept = 0, color = "gray85", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 1, color = "gray90", linewidth = 0.4, linetype = "dashed")
  if (is.finite(df$lwr[[1L]]) && is.finite(df$upr[[1L]])) {
    p <- p + ggplot2::geom_segment(
      ggplot2::aes(x = .data$lwr, xend = .data$upr, y = .data$label_y, yend = .data$label_y),
      linewidth = 1,
      color = point_color
    )
  }
  p <- p +
    ggplot2::geom_point(size = 3, color = point_color) +
    ggplot2::scale_x_continuous(limits = c(-1, 1)) +
    ggplot2::scale_y_continuous(breaks = 1, labels = "Cohen's kappa") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      ...
    ) +
    ggplot2::labs(
      title = title %||% "Cohen's kappa agreement",
      x = "Kappa"
    )

  if (isTRUE(show_value)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$estimate_label),
      nudge_y = 0.12,
      size = value_text_size,
      color = "black"
    )
    if (any(!is.na(df$ci_label))) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data$ci_label),
        nudge_y = -0.12,
        size = ci_text_size,
        color = "gray30",
        na.rm = TRUE
      )
    }
  }
  p
}
