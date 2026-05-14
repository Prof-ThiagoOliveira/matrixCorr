#' @title Pairwise Weighted Cohen's Kappa for Ordered Ratings
#'
#' @description
#' Computes weighted Cohen's kappa for either a pair of ordinal rating vectors
#' or all pairwise combinations of ordinal columns in a matrix or data frame.
#'
#' @param data In matrix mode, a matrix or data frame whose rows are
#'   observational units and whose columns are raters or classifiers. Each
#'   column contains ordered category ratings. In two-vector mode, the first
#'   ordinal rating vector.
#' @param y Optional second ordinal rating vector. When supplied, the function
#'   returns a single weighted kappa estimate for \code{data} and \code{y}.
#' @param weights Weight specification. The default \code{"quadratic"} uses
#'   Fleiss-Cohen-style agreement weights. \code{"linear"} uses equal-spacing
#'   agreement weights. \code{"unweighted"} reduces to ordinary Cohen's kappa.
#'   Custom numeric weight matrices must be symmetric agreement weights with
#'   diagonal 1 and entries in \code{[0, 1]}.
#' @param levels Optional ordered category levels. Weighted kappa depends on
#'   category order and will not silently alphabetise arbitrary labels. If
#'   omitted, order is inferred only when all involved ratings are ordered
#'   factors with identical levels or when all involved ratings are numeric or
#'   integer and can be ordered by their sorted unique values.
#' @param na_method Character scalar controlling missing-data handling.
#'   \code{"error"} rejects missing values. \code{"complete"} applies listwise
#'   deletion across retained columns before computing the matrix. In matrix
#'   mode, \code{"pairwise"} uses pair-specific complete observations.
#' @param ci Logical; if \code{TRUE}, attach large-sample delta-method
#'   confidence intervals.
#' @param p_value Logical; if \code{TRUE}, attach large-sample delta-method
#'   Wald test statistics and p-values for \eqn{H_0: \kappa_w = 0}.
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
#' This implementation uses agreement/similarity weights internally:
#' \itemize{
#'   \item \code{"unweighted"}: \eqn{w_{ij} = 1(i = j)}
#'   \item \code{"linear"}: \eqn{w_{ij} = 1 - |i - j| / (K - 1)}
#'   \item \code{"quadratic"}: \eqn{w_{ij} = 1 - (i - j)^2 / (K - 1)^2}
#'   }
#' The \code{"quadratic"} scheme is the default and corresponds to
#' Fleiss-Cohen-style agreement weights. The \code{"linear"} scheme
#' corresponds to equal-spacing agreement weights.
#'
#' For matrix mode, rows are shared observational units and columns are raters
#' or classifiers. A common category map is resolved in R and the main
#' computation is performed in C++. Missing-data handling follows the usual
#' \pkg{matrixCorr} \code{na_method} conventions. Pairwise complete counts are
#' stored in \code{attr(x, "diagnostics")$n_complete}.
#'
#' Confidence intervals and p-values use a large-sample multinomial
#' delta-method approximation. Use [cohen_kappa()] for unordered nominal
#' categories where all disagreements are equally serious.
#'
#' @return
#' If \code{y} is supplied and neither \code{ci} nor \code{p_value} is
#' requested, a single numeric estimate. If \code{y} is supplied and inference
#' is requested, a small named list is returned. Otherwise a symmetric
#' matrix-style result with estimator class \code{weighted_kappa}.
#'
#' @references
#' Cohen, J. (1968). Weighted kappa: Nominal scale agreement with provision for
#' scaled disagreement or partial credit. \emph{Psychological Bulletin},
#' 70(4), 213-220. \doi{10.1037/h0026256}
#'
#' @seealso
#' [cohen_kappa()] for unweighted two-rater nominal agreement;
#' [multirater_kappa()] for panel-level nominal agreement among three or more
#' raters.
#'
#' @examples
#' raters <- data.frame(
#'   r1 = ordered(c("low", "low", "mid", "high", "high"),
#'                levels = c("low", "mid", "high")),
#'   r2 = ordered(c("low", "mid", "mid", "high", "high"),
#'                levels = c("low", "mid", "high")),
#'   r3 = ordered(c("low", "low", "high", "high", "mid"),
#'                levels = c("low", "mid", "high"))
#' )
#'
#' wk <- weighted_kappa(raters)
#' print(wk)
#' summary(wk)
#' plot(wk)
#'
#' x <- raters$r1
#' y <- raters$r2
#' weighted_kappa(x, y, weights = "linear")
#'
#' @author Thiago de Paula Oliveira
#' @export
weighted_kappa <- function(data,
                           y = NULL,
                           weights = c("quadratic", "linear", "unweighted"),
                           levels = NULL,
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

  enc <- .mc_encode_weighted_kappa_data(
    data = data,
    y = y,
    levels = levels,
    na_method = na_method
  )
  W <- .mc_resolve_weighted_kappa_weights(weights, enc$n_levels)

  if (!is.null(y)) {
    if (!identical(output_cfg$output, "matrix")) {
      abort_bad_arg(
        "output",
        message = "must be {.val matrix} when {.arg y} is supplied."
      )
    }

    x_codes <- enc$x
    y_codes <- enc$y
    if (identical(na_method, "error") && (anyNA(x_codes) || anyNA(y_codes))) {
      abort_bad_arg(
        "data",
        message = "and {.arg y} must not contain missing values when {.arg na_method} is {.val error}.",
        .hint = "Use {.arg na_method} = {.val pairwise} or {.val complete} to drop incomplete pairs."
      )
    }
    if (!identical(na_method, "error")) {
      keep <- stats::complete.cases(x_codes, y_codes)
      x_codes <- x_codes[keep]
      y_codes <- y_codes[keep]
    }

    fit <- weighted_kappa_pair_cpp(
      x_codes,
      y_codes,
      W,
      return_inference = isTRUE(ci) || isTRUE(p_value),
      conf_level = conf_level
    )
    if (!isTRUE(ci) && !isTRUE(p_value)) {
      return(as.numeric(fit$estimate))
    }

    out <- list(
      estimate = as.numeric(fit$estimate),
      n_complete = as.integer(fit$n_complete),
      observed_agreement = as.numeric(fit$observed_agreement),
      expected_agreement = as.numeric(fit$expected_agreement),
      levels = enc$levels,
      weights = W,
      weight_type = attr(W, "weight_type", exact = TRUE)
    )
    if (isTRUE(ci)) {
      out$lwr <- as.numeric(fit$lwr)
      out$upr <- as.numeric(fit$upr)
      out$conf_level <- conf_level
      out$ci_method <- "delta"
    }
    if (isTRUE(p_value)) {
      out$se <- as.numeric(fit$se)
      out$statistic <- as.numeric(fit$statistic)
      out$p_value <- as.numeric(fit$p_value)
    }
    return(out)
  }

  X <- enc$X
  colnames_data <- enc$diagnostics$colnames_data
  dn <- .mc_square_dimnames(colnames_data)

  if (identical(na_method, "error") && anyNA(X)) {
    abort_bad_arg(
      "data",
      message = "contains missing values when {.arg na_method} is {.val error}.",
      .hint = "Use {.arg na_method} = {.val pairwise} or {.val complete} to handle incomplete units."
    )
  }

  listwise_meta <- NULL
  if (identical(na_method, "complete")) {
    cc <- .mc_weighted_kappa_complete_cases(X, min_n = 2L)
    X <- cc$X
    listwise_meta <- cc$diagnostics
  }

  direct_threshold <- .mc_supports_direct_threshold_path(
    method = "weighted_kappa",
    na_method = na_method,
    ci = ci,
    output = output_cfg$output,
    threshold = output_cfg$threshold,
    symmetric = TRUE,
    exact = TRUE,
    has_inference = p_value
  )

  prev_threads <- .mc_prepare_omp_threads(
    n_threads,
    n_threads_missing = missing(n_threads)
  )
  if (!is.null(prev_threads)) {
    on.exit(.mc_exit_omp_threads(prev_threads), add = TRUE)
  }

  if (direct_threshold) {
    trip <- weighted_kappa_threshold_triplets_cpp(
      X,
      W,
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
      estimator_class = "weighted_kappa",
      method = "weighted_kappa",
      description = "Pairwise weighted Cohen's kappa agreement matrix",
      threshold = output_cfg$threshold,
      diag = output_cfg$diag,
      source_dim = as.integer(c(ncol(X), ncol(X))),
      source_dimnames = dn,
      diagnostics = list(
        n_complete = n_complete_mat,
        levels = enc$levels,
        weights = .mc_weighted_kappa_weights_for_diagnostics(W, enc$levels),
        weight_type = attr(W, "weight_type", exact = TRUE)
      ),
      symmetric = TRUE
    ))
  }

  fit <- weighted_kappa_matrix_cpp(
    X,
    W,
    pairwise_complete = identical(na_method, "pairwise"),
    return_inference = isTRUE(ci) || isTRUE(p_value),
    conf_level = conf_level,
    n_threads = n_threads
  )
  out <- .mc_structure_weighted_kappa_matrix(
    raw = fit,
    colnames_data = colnames_data,
    weights = W,
    levels = enc$levels,
    ci = ci,
    p_value = p_value,
    conf_level = conf_level
  )
  attr(out, "diagnostics") <- .mc_merge_diagnostics(
    attr(out, "diagnostics", exact = TRUE),
    if (identical(na_method, "pairwise")) {
      list(
        na_method = "pairwise",
        common_sample = FALSE
      )
    } else {
      listwise_meta
    }
  )

  .mc_finalize_corr_output_fast(
    out,
    output = output_cfg$output,
    threshold = output_cfg$threshold,
    diag = output_cfg$diag
  )
}

.mc_weighted_kappa_prepare_vector <- function(x,
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
    key <- as.character(x)
    key[is.na(x)] <- NA_character_
    return(list(
      key = key,
      is_ordered = is.ordered(x),
      ordered_levels = if (is.ordered(x)) levels(x) else NULL,
      numeric_values = NULL,
      numeric_infer = FALSE
    ))
  }

  if (is.numeric(x) || is.integer(x)) {
    if (any(!is.na(x) & !is.finite(x))) {
      abort_bad_arg(arg, message = "must not contain NaN or infinite values.")
    }
    key <- as.character(x)
    key[is.na(x)] <- NA_character_
    return(list(
      key = key,
      is_ordered = FALSE,
      ordered_levels = NULL,
      numeric_values = x,
      numeric_infer = TRUE
    ))
  }

  if (is.character(x) || is.logical(x)) {
    key <- as.character(x)
    key[is.na(x)] <- NA_character_
    return(list(
      key = key,
      is_ordered = FALSE,
      ordered_levels = NULL,
      numeric_values = NULL,
      numeric_infer = FALSE
    ))
  }

  abort_bad_arg(
    arg,
    message = paste(
      "must be an ordered factor, factor, character, logical, integer,",
      "or numeric vector."
    )
  )
}

.mc_weighted_kappa_normalize_levels <- function(levels) {
  if (is.null(levels)) {
    return(NULL)
  }
  if (is.list(levels) || is.matrix(levels)) {
    abort_bad_arg("levels", message = "must be an atomic vector of ordered categories.")
  }
  lev <- if (is.factor(levels)) as.character(levels) else levels
  if (!is.atomic(lev)) {
    abort_bad_arg("levels", message = "must be an atomic vector of ordered categories.")
  }
  if (!length(lev)) {
    abort_bad_arg("levels", message = "must contain at least one category.")
  }
  if (anyNA(lev)) {
    abort_bad_arg("levels", message = "must not contain missing values.")
  }
  lev_key <- as.character(lev)
  if (anyDuplicated(lev_key)) {
    abort_bad_arg("levels", message = "must not contain duplicate categories.")
  }
  lev
}

.mc_weighted_kappa_resolve_levels <- function(prepared, levels = NULL) {
  if (!is.null(levels)) {
    return(.mc_weighted_kappa_normalize_levels(levels))
  }

  if (all(vapply(prepared, `[[`, logical(1), "is_ordered"))) {
    base_levels <- prepared[[1L]]$ordered_levels
    same_levels <- vapply(
      prepared,
      function(z) identical(z$ordered_levels, base_levels),
      logical(1)
    )
    if (all(same_levels)) {
      return(base_levels)
    }
  }

  if (all(vapply(prepared, `[[`, logical(1), "numeric_infer"))) {
    vals <- unlist(lapply(prepared, function(z) z$numeric_values[!is.na(z$numeric_values)]), use.names = FALSE)
    if (!length(vals)) {
      abort_bad_arg(
        "levels",
        message = "must be supplied when no non-missing ratings are available to infer the category order."
      )
    }
    return(sort(unique(vals)))
  }

  abort_bad_arg(
    "levels",
    message = paste(
      "must be supplied when the category order cannot be inferred.",
      "Automatic order inference is limited to identical ordered-factor levels",
      "or numeric/integer ratings."
    )
  )
}

.mc_weighted_kappa_weights_for_diagnostics <- function(weights, levels) {
  out <- unclass(weights)
  dimnames(out) <- list(as.character(levels), as.character(levels))
  out
}

#' @keywords internal
#' @noRd
.mc_encode_weighted_kappa_data <- function(data,
                                           y = NULL,
                                           levels = NULL,
                                           na_method) {
  validate_na_method(
    na_method,
    arg = "na_method",
    allowed = c("error", "pairwise", "complete")
  )

  if (!is.null(y)) {
    check_same_length(data, y, arg_x = "data", arg_y = "y")
    prepared <- list(
      data = .mc_weighted_kappa_prepare_vector(data, arg = "data"),
      y = .mc_weighted_kappa_prepare_vector(y, arg = "y")
    )
    resolved_levels <- .mc_weighted_kappa_resolve_levels(prepared, levels = levels)
    level_key <- as.character(resolved_levels)

    x_codes <- match(prepared$data$key, level_key)
    y_codes <- match(prepared$y$key, level_key)
    x_codes[is.na(prepared$data$key)] <- NA_integer_
    y_codes[is.na(prepared$y$key)] <- NA_integer_

    bad_x <- !is.na(prepared$data$key) & is.na(x_codes)
    bad_y <- !is.na(prepared$y$key) & is.na(y_codes)
    if (any(bad_x)) {
      abort_bad_arg("data", message = "contains values not present in {.arg levels}.")
    }
    if (any(bad_y)) {
      abort_bad_arg("y", message = "contains values not present in {.arg levels}.")
    }

    return(list(
      x = as.integer(x_codes),
      y = as.integer(y_codes),
      levels = resolved_levels,
      n_levels = length(resolved_levels),
      diagnostics = list(mode = "pair")
    ))
  }

  cols <- NULL
  colnames_data <- NULL
  if (is.data.frame(data)) {
    cols <- unclass(data)
    colnames_data <- names(data)
  } else if (is.matrix(data)) {
    cols <- lapply(seq_len(ncol(data)), function(j) data[, j])
    colnames_data <- colnames(data)
  } else {
    abort_bad_arg(
      "data",
      message = "must be a matrix or data frame in matrix mode."
    )
  }

  if (length(cols) < 2L) {
    abort_bad_arg(
      "data",
      message = "must contain at least 2 rating columns."
    )
  }

  prepared <- vector("list", length(cols))
  for (j in seq_along(cols)) {
    col_arg <- if (!is.null(colnames_data) && nzchar(colnames_data[[j]] %||% "")) {
      sprintf("data$%s", colnames_data[[j]])
    } else {
      sprintf("data[[%d]]", j)
    }
    prepared[[j]] <- .mc_weighted_kappa_prepare_vector(cols[[j]], arg = col_arg)
  }

  resolved_levels <- .mc_weighted_kappa_resolve_levels(prepared, levels = levels)
  level_key <- as.character(resolved_levels)
  X <- matrix(
    NA_integer_,
    nrow = length(prepared[[1L]]$key),
    ncol = length(prepared),
    dimnames = list(NULL, colnames_data)
  )
  for (j in seq_along(prepared)) {
    codes <- match(prepared[[j]]$key, level_key)
    codes[is.na(prepared[[j]]$key)] <- NA_integer_
    bad <- !is.na(prepared[[j]]$key) & is.na(codes)
    if (any(bad)) {
      abort_bad_arg(
        "data",
        message = "contains values not present in {.arg levels}."
      )
    }
    X[, j] <- as.integer(codes)
  }

  list(
    X = X,
    levels = resolved_levels,
    n_levels = length(resolved_levels),
    diagnostics = list(
      mode = "matrix",
      colnames_data = colnames_data
    )
  )
}

#' @keywords internal
#' @noRd
.mc_resolve_weighted_kappa_weights <- function(weights, n_levels) {
  check_scalar_numeric(
    n_levels,
    arg = "n_levels",
    lower = 1,
    closed_lower = TRUE
  )

  if (is.character(weights)) {
    if (length(weights) != 1L) {
      if (identical(weights, c("quadratic", "linear", "unweighted"))) {
        weights <- "quadratic"
      } else {
        abort_bad_arg("weights", message = "must be a single character weight scheme or a numeric matrix.")
      }
    }
    if (is.na(weights)) {
      abort_bad_arg("weights", message = "must be a single character weight scheme or a numeric matrix.")
    }
    alias_map <- c(
      quadratic = "quadratic",
      fleiss_cohen = "quadratic",
      "fleiss-cohen" = "quadratic",
      squared = "quadratic",
      linear = "linear",
      equal_spacing = "linear",
      "equal-spacing" = "linear",
      equal = "linear",
      unweighted = "unweighted"
    )
    key <- tolower(weights[[1L]])
    resolved <- unname(alias_map[key])
    if (!length(resolved) || is.na(resolved)) {
      abort_bad_arg(
        "weights",
        message = paste(
          "must be one of {.val quadratic}, {.val linear},",
          "{.val unweighted}, or a numeric square matrix."
        )
      )
    }

    K <- as.integer(n_levels)
    if (identical(resolved, "unweighted")) {
      W <- diag(1, nrow = K, ncol = K)
    } else {
      idx <- outer(seq_len(K), seq_len(K), "-")
      if (K <= 1L) {
        W <- matrix(1, nrow = K, ncol = K)
      } else if (identical(resolved, "linear")) {
        W <- 1 - abs(idx) / (K - 1)
      } else {
        W <- 1 - (idx^2) / ((K - 1)^2)
      }
    }
    storage.mode(W) <- "double"
    attr(W, "weight_type") <- resolved
    return(W)
  }

  if (!is.numeric(weights) || !is.matrix(weights)) {
    abort_bad_arg(
      "weights",
      message = "must be a single character weight scheme or a numeric square matrix."
    )
  }
  check_matrix_dims(weights, arg = "weights")
  if (nrow(weights) != ncol(weights)) {
    abort_bad_arg("weights", message = "must be a square matrix.")
  }
  if (nrow(weights) != as.integer(n_levels)) {
    abort_bad_arg(
      "weights",
      message = "must have dimension {n_levels} x {n_levels}.",
      n_levels = as.integer(n_levels)
    )
  }
  if (any(!is.finite(weights))) {
    abort_bad_arg("weights", message = "must contain only finite values.")
  }
  tol <- .Machine$double.eps^0.5
  if (any(abs(diag(weights) - 1) > tol)) {
    abort_bad_arg("weights", message = "must have diagonal entries equal to 1.")
  }
  if (any(weights < -tol | weights > 1 + tol)) {
    abort_bad_arg("weights", message = "must have all values in the interval [0, 1].")
  }
  if (!isTRUE(isSymmetric(weights, tol = tol))) {
    abort_bad_arg(
      "weights",
      message = "must be symmetric for this matrix-oriented implementation.",
      .hint = "Asymmetric custom agreement weights are not yet supported."
    )
  }

  W <- unclass(weights)
  storage.mode(W) <- "double"
  attr(W, "weight_type") <- "custom"
  W
}

#' @keywords internal
#' @noRd
.mc_weighted_kappa_complete_cases <- function(X, min_n = 2L) {
  keep <- stats::complete.cases(X)
  list(
    X = X[keep, , drop = FALSE],
    keep = keep,
    diagnostics = list(
      na_method = "complete",
      n_original = nrow(X),
      n_complete_rows = sum(keep),
      n_removed = nrow(X) - sum(keep),
      complete_rows = keep,
      common_sample = TRUE,
      min_n = as.integer(min_n)
    )
  )
}

#' @keywords internal
#' @noRd
.mc_structure_weighted_kappa_matrix <- function(raw,
                                                colnames_data,
                                                weights,
                                                levels,
                                                ci,
                                                p_value,
                                                conf_level) {
  est <- .mc_set_matrix_dimnames(raw$est, colnames_data)
  n_complete <- .mc_set_matrix_dimnames(raw$n_complete, colnames_data)
  observed_agreement <- .mc_set_matrix_dimnames(raw$observed_agreement, colnames_data)
  expected_agreement <- .mc_set_matrix_dimnames(raw$expected_agreement, colnames_data)

  ci_attr <- NULL
  inference_attr <- NULL
  if (!is.null(raw$se)) {
    se <- .mc_set_matrix_dimnames(raw$se, colnames_data)
    lwr <- .mc_set_matrix_dimnames(raw$lwr, colnames_data)
    upr <- .mc_set_matrix_dimnames(raw$upr, colnames_data)
    statistic <- .mc_set_matrix_dimnames(raw$statistic, colnames_data)
    p_mat <- .mc_set_matrix_dimnames(raw$p_value, colnames_data)
    if (isTRUE(ci)) {
      ci_attr <- list(
        est = est,
        lwr.ci = lwr,
        upr.ci = upr,
        conf.level = conf_level,
        ci.method = "delta"
      )
    }
    if (isTRUE(p_value)) {
      inference_attr <- list(
        estimate = est,
        se = se,
        statistic = statistic,
        p_value = p_mat,
        n_obs = n_complete
      )
    }
  }

  .mc_new_corr_matrix(
    mat = est,
    estimator_class = "weighted_kappa",
    method = "weighted_kappa",
    description = "Pairwise weighted Cohen's kappa agreement matrix",
    diagnostics = list(
      n_complete = n_complete,
      observed_agreement = observed_agreement,
      expected_agreement = expected_agreement,
      levels = levels,
      weights = .mc_weighted_kappa_weights_for_diagnostics(weights, levels),
      weight_type = attr(weights, "weight_type", exact = TRUE)
    ),
    ci = ci_attr,
    conf.level = if (isTRUE(ci)) conf_level else NULL,
    symmetric = TRUE,
    extra_attrs = list(inference = inference_attr)
  )
}

#' @rdname weighted_kappa
#' @method print weighted_kappa
#' @param x A matrix-style \code{weighted_kappa} object.
#' @param digits Integer; number of decimal places for displayed values.
#' @param n Optional preview row threshold.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns.
#' @param width Optional display width.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @export
print.weighted_kappa <- function(x,
                                 digits = 4,
                                 n = NULL,
                                 topn = NULL,
                                 max_vars = NULL,
                                 width = NULL,
                                 show_ci = NULL,
                                 ...) {
  .mc_print_corr_matrix(
    x,
    header = "Weighted Cohen's kappa agreement matrix",
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
