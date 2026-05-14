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
#'   numeric, all treated as nominal categories. In two-vector mode, the first
#'   nominal rating vector.
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
#' kappa. Weighted kappa is intentionally not included here.
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
#' When \code{ci = TRUE} or \code{p_value = TRUE}, inferential quantities are
#' based on a large-sample multinomial delta-method approximation for the
#' unweighted kappa estimator.
#'
#' @return
#' If \code{y} is supplied, a numeric scalar with attributes
#' \code{diagnostics}, and optionally \code{ci}, \code{inference}, and
#' \code{conf.level}. Otherwise a symmetric correlation-style result with
#' estimator class \code{cohen_kappa}. Matrix outputs are created with
#' \code{.mc_new_corr_matrix()}; sparse and edge-list outputs dispatch through
#' the shared matrixCorr output-routing helpers.
#'
#' @references
#' Cohen, J. (1960). A coefficient of agreement for nominal scales.
#' \emph{Educational and Psychological Measurement}, 20(1), 37-46.
#' \doi{10.1177/001316446002000104}
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
  out <- est
  attr(out, "method") <- "cohen_kappa"
  attr(out, "description") <- "Pairwise Cohen's kappa agreement"
  attr(out, "package") <- "matrixCorr"
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

#' @method print cohen_kappa
#' @export
print.cohen_kappa <- function(x, digits = 4, n = NULL, topn = NULL,
                              max_vars = NULL, width = NULL,
                              show_ci = NULL, ...) {
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
