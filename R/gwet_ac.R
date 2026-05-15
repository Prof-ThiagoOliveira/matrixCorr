#' @title Gwet's AC1 and AC2 Agreement Coefficients
#'
#' @description
#' Estimates Gwet's chance-corrected agreement coefficient for either
#' unweighted nominal agreement (AC1) or weighted agreement (AC2) in two-rater
#' pairwise form or panel-level multi-rater form.
#'
#' @param data Input data. If \code{y} is supplied, \code{data} is the first
#'   rater vector. Otherwise:
#'   for \code{input = "pairwise"}, rows are units and columns are raters;
#'   for \code{input = "ratings"}, rows are items and columns are raters;
#'   for \code{input = "counts"}, rows are items and columns are categories
#'   containing rater counts.
#' @param y Optional second rater vector for scalar two-rater agreement.
#' @param weights Weight specification. \code{"unweighted"} yields AC1.
#'   Weighted schemes yield AC2. Supported built-in schemes are
#'   \code{"linear"}, \code{"quadratic"}, \code{"ordinal"},
#'   \code{"radical"}, \code{"ratio"}, \code{"circular"}, and
#'   \code{"bipolar"}. A custom numeric agreement-weight matrix may also be
#'   supplied.
#' @param levels Optional explicit category labels. This is mainly useful when
#'   unused categories should still be retained in the calculation.
#' @param input One of \code{"pairwise"}, \code{"ratings"}, or \code{"counts"}.
#' @param na_method Missing-data rule. For scalar two-rater and pairwise matrix
#'   modes, \code{"error"} rejects missing values, \code{"complete"} applies
#'   listwise deletion in matrix mode, and \code{"pairwise"} uses complete
#'   pairs within each pair of raters. For \code{input = "ratings"},
#'   \code{"available"} retains items with at least \code{min_raters}
#'   observed ratings. For \code{input = "counts"}, missing values are not
#'   allowed.
#' @param min_raters Minimum number of observed ratings required for an item to
#'   be retained in panel mode. Must be at least \code{2}.
#' @param by_category Logical; if \code{TRUE}, attach one-vs-all category-wise
#'   AC1 estimates for panel-level fits. This is available only for the
#'   unweighted coefficient.
#' @param ci Logical; if \code{TRUE}, attach confidence intervals.
#' @param p_value Logical; if \code{TRUE}, attach test statistics and p-values.
#' @param conf_level Confidence level for intervals. Default is \code{0.95}.
#' @param se_method Standard-error method. \code{"asymptotic"} uses the
#'   analytic large-sample formula. \code{"jackknife"} is retained as an
#'   alternative for panel-level fits. \code{"none"} disables inferential
#'   output.
#' @param n_threads Integer \eqn{\ge 1}. Number of OpenMP threads used by the
#'   C++ backend.
#' @param output One of \code{"matrix"}, \code{"sparse"}, or
#'   \code{"edge_list"} for pairwise matrix mode.
#' @param threshold Numeric threshold used only for thresholded pairwise matrix
#'   output.
#' @param diag Logical; include diagonal entries in thresholded pairwise matrix
#'   output.
#' @param verbose Logical; if \code{TRUE}, emit short progress messages.
#' @param ... Reserved for future extensions. Unsupported extra arguments are
#'   rejected.
#'
#' @details
#' For two raters, let \eqn{w_{ab}} denote the agreement weight for cell
#' \eqn{(a,b)}, \eqn{p_{ab}} the cell proportion, and
#' \eqn{\pi_k = (p_{k+} + p_{+k})/2}. Then
#' \deqn{ P_o = \sum_{a,b} w_{ab} p_{ab}, }
#' \deqn{ P_e =
#'        \frac{\sum_{a,b} w_{ab}}{q(q-1)}
#'        \sum_{k=1}^q \pi_k(1-\pi_k), }
#' and
#' \deqn{ AC = \frac{P_o - P_e}{1 - P_e}. }
#' With identity weights \eqn{w_{ab} = 1(a=b)} this is Gwet's AC1, the
#' nominal-agreement coefficient. With non-identity agreement weights this is
#' Gwet's AC2, which extends the same chance-correction idea to ordered or
#' partially creditable disagreement by letting near disagreements receive
#' larger agreement weights than distant disagreements.
#'
#' For panel-level counts, if \eqn{r_{ik}} is the number of raters assigning
#' item \eqn{i} to category \eqn{k}, \eqn{r_i = \sum_k r_{ik}}, and
#' \eqn{w_{kl}} is the agreement-weight matrix, then the item-level agreement is
#' \deqn{ A_i =
#'        \frac{1}{r_i(r_i-1)}
#'        \sum_k r_{ik}\left(\sum_l w_{kl}r_{il}-1\right). }
#' The observed agreement is the mean of \eqn{A_i} over retained items. The
#' chance term uses the average item-wise category proportions
#' \eqn{\pi_k = m^{-1}\sum_i r_{ik}/r_i}, where \eqn{m} is the number of
#' retained items.
#'
#' \strong{Confidence intervals and inference.}
#'
#' The primary inferential path is the analytic large-sample method. For two
#' raters, define
#' \deqn{ s_{ab}
#'        =
#'        w_{ab}
#'        -
#'        \frac{2(1-AC)\sum_{u,v} w_{uv}}{q(q-1)}
#'        \left(1 - \frac{\pi_a + \pi_b}{2}\right). }
#' The backend variance estimator is
#' \deqn{ \widehat{\mathrm{Var}}(\widehat{AC})
#'        =
#'        \frac{1}{n(1-P_e)^2}
#'        \left[
#'        \sum_{a,b} p_{ab} s_{ab}^2
#'        -
#'        \left\{P_o - 2(1-AC)P_e\right\}^2
#'        \right], }
#' where \eqn{n} is the number of complete paired ratings.
#'
#' For panel-level counts, write
#' \deqn{ AC_i = \frac{A_i - P_e}{1 - P_e}, }
#' and define the item-specific chance term
#' \deqn{ P_{e,i}
#'        =
#'        \frac{\sum_{k,l} w_{kl}}{q(q-1)}
#'        \sum_{k=1}^q \frac{r_{ik}}{r_i}(1-\pi_k). }
#' The asymptotic linearised contribution used by the backend is
#' \deqn{ \xi_i
#'        =
#'        AC_i
#'        -
#'        \frac{2(1-AC)}{1-P_e}(P_{e,i} - P_e), }
#' giving
#' \deqn{ \widehat{\mathrm{Var}}(\widehat{AC})
#'        =
#'        \frac{1}{m(m-1)}
#'        \sum_{i=1}^m (\xi_i - AC)^2, }
#' where \eqn{m} is the number of retained items.
#'
#' The reported standard error is
#' \eqn{\widehat{\mathrm{se}}(\widehat{AC}) =
#' \sqrt{\widehat{\mathrm{Var}}(\widehat{AC})}}, and the confidence interval is
#' the t interval
#' \deqn{ \widehat{AC} \pm t_{1-\alpha/2,\nu}\,
#'        \widehat{\mathrm{se}}(\widehat{AC}), }
#' truncated to \eqn{[-1,1]}. Here \eqn{\nu = n-1} for the two-rater path and
#' \eqn{\nu = m-1} for panel-level asymptotic inference. The reported test
#' statistic is the corresponding t ratio for \eqn{H_0: AC = 0}. For
#' panel-level fits, \code{se_method = "jackknife"} remains available as a
#' second option, using the leave-one-item-out variance
#' \deqn{ \widehat{\mathrm{Var}}_{\mathrm{JK}}(\widehat{AC})
#'        =
#'        \frac{m-1}{m}
#'        \sum_{i=1}^m
#'        \left(\widehat{AC}_{(-i)} - \overline{AC}_{(-\cdot)}\right)^2. }
#'
#' @return
#' If \code{y} is supplied, a scalar numeric object of class
#' \code{c("gwet_ac", "numeric")}. If \code{input = "pairwise"}, a
#' \code{corr_result}. If \code{input} is \code{"ratings"} or \code{"counts"},
#' a one-row data frame with class
#' \code{c("gwet_ac", "agreement_result", "data.frame")}.
#'
#' @references
#' Gwet, K. L. (2008). Computing inter-rater reliability and its variance in
#' the presence of high agreement. \emph{British Journal of Mathematical and
#' Statistical Psychology}, 61, 29-48. \doi{10.1348/000711006X126600}
#'
#' @seealso
#' [cohen_kappa()] for two-rater nominal kappa;
#' [multirater_kappa()] for panel-level nominal kappa;
#' [weighted_kappa()] for weighted Cohen-type agreement.
#'
#' @examples
#' x <- c("A", "A", "B", "B", "A", "C")
#' y <- c("A", "B", "B", "B", "A", "C")
#' gwet_ac(x, y)
#' gwet_ac(x, y, weights = "quadratic", levels = c("A", "B", "C"))
#'
#' raters <- data.frame(
#'   r1 = c("A", "A", "B", "C", "A", "B"),
#'   r2 = c("A", "B", "B", "C", "A", "B"),
#'   r3 = c("A", "A", "B", "B", "A", "C"),
#'   stringsAsFactors = FALSE
#' )
#'
#' fit_pw <- gwet_ac(raters)
#' fit_panel <- gwet_ac(raters, input = "ratings")
#' print(fit_pw)
#' print(fit_panel)
#'
#' @author Thiago de Paula Oliveira
#' @export
gwet_ac <- function(data,
                    y = NULL,
                    weights = c("unweighted", "linear", "quadratic", "ordinal", "radical", "ratio", "circular", "bipolar"),
                    levels = NULL,
                    input = c("pairwise", "ratings", "counts"),
                    na_method = c("error", "pairwise", "complete", "available"),
                    min_raters = 2L,
                    by_category = FALSE,
                    ci = FALSE,
                    p_value = FALSE,
                    conf_level = 0.95,
                    se_method = c("asymptotic", "jackknife", "none"),
                    n_threads = getOption("matrixCorr.threads", 1L),
                    output = c("matrix", "sparse", "edge_list"),
                    threshold = 0,
                    diag = TRUE,
                    verbose = FALSE,
                    ...) {
  .mc_extract_legacy_aliases(list(...), allowed = character())
  output_cfg <- .mc_validate_thresholded_output_request(
    output = output,
    threshold = threshold,
    diag = diag
  )
  input <- match_arg(input, c("pairwise", "ratings", "counts"), arg_name = "input")
  na_method <- validate_na_method(
    na_method,
    arg = "na_method",
    allowed = c("error", "pairwise", "complete", "available")
  )
  min_raters <- check_scalar_int_pos(min_raters, arg = "min_raters")
  if (min_raters < 2L) {
    abort_bad_arg("min_raters", message = "must be >= 2.")
  }
  check_bool(by_category, arg = "by_category")
  check_bool(ci, arg = "ci")
  check_bool(p_value, arg = "p_value")
  check_bool(verbose, arg = "verbose")
  check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")

  se_method <- match_arg(se_method, c("asymptotic", "jackknife", "none"), arg_name = "se_method")
  if ((isTRUE(ci) || isTRUE(p_value)) && identical(se_method, "none")) {
    abort_bad_arg("se_method", message = "must not be 'none' when ci or p_value is TRUE.")
  }

  if (!is.null(y)) {
    if (!identical(output_cfg$output, "matrix")) {
      abort_bad_arg("output", message = "must be {.val matrix} when {.arg y} is supplied.")
    }
    if (identical(na_method, "available")) {
      abort_bad_arg("na_method", message = "na_method = 'available' is not defined when {.arg y} is supplied.")
    }
    if (identical(se_method, "jackknife") && (isTRUE(ci) || isTRUE(p_value))) {
      abort_bad_arg("se_method", message = "{.val jackknife} is not defined for two-rater AC1/AC2.")
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

    enc <- .mc_encode_gwet_pair_labels(x_labels, y_labels, levels = levels)
    W <- .mc_resolve_gwet_weights(weights, enc$levels)
    fit <- gwet_ac_pair_cpp(
      enc$x,
      enc$y,
      weights = W,
      return_inference = isTRUE(ci) || isTRUE(p_value),
      conf_level = conf_level
    )
    return(.mc_gwet_ac_scalar(
      fit,
      weights = W,
      levels = enc$levels,
      ci = ci,
      p_value = p_value,
      conf_level = conf_level
    ))
  }

  if (identical(input, "pairwise")) {
    if (identical(na_method, "available")) {
      abort_bad_arg("na_method", message = "na_method = 'available' is not defined for pairwise matrix agreement.")
    }
    if (identical(se_method, "jackknife") && (isTRUE(ci) || isTRUE(p_value))) {
      abort_bad_arg("se_method", message = "{.val jackknife} is not defined for pairwise matrix AC1/AC2.")
    }
    return(.mc_gwet_ac_pairwise_matrix(
      data = data,
      weights = weights,
      levels = levels,
      na_method = na_method,
      ci = ci,
      p_value = p_value,
      conf_level = conf_level,
      n_threads = n_threads,
      output_cfg = output_cfg
    ))
  }

  if (!identical(output_cfg$output, "matrix") || output_cfg$thresholded) {
    abort_bad_arg(
      "output",
      message = "must be {.val matrix} with {.arg threshold} = 0 for panel-level agreement."
    )
  }
  if (identical(input, "ratings") && identical(na_method, "pairwise")) {
    abort_bad_arg(
      "na_method",
      message = "na_method = 'pairwise' is not defined for panel-level agreement; use 'available' to allow item-specific missing ratings."
    )
  }

  prep <- if (identical(input, "ratings")) {
    .mc_prepare_multirater_kappa_ratings(
      data = data,
      levels = levels,
      na_method = na_method,
      min_raters = min_raters
    )
  } else {
    .mc_prepare_multirater_kappa_counts(
      data = data,
      levels = levels,
      min_raters = min_raters
    )
  }
  W <- .mc_resolve_gwet_weights(weights, prep$categories)
  if (isTRUE(by_category) && !identical(attr(W, "weight_type", exact = TRUE), "unweighted")) {
    abort_bad_arg("by_category", message = "must be FALSE for weighted AC2 fits.")
  }

  inference_requested <- isTRUE(ci) || isTRUE(p_value)
  inference_code <- if (!inference_requested) {
    0L
  } else if (identical(se_method, "jackknife")) {
    2L
  } else {
    1L
  }

  prev_threads <- .mc_prepare_omp_threads(
    n_threads,
    n_threads_missing = missing(n_threads)
  )
  if (!is.null(prev_threads)) {
    on.exit(.mc_exit_omp_threads(prev_threads), add = TRUE)
  }

  raw <- if (identical(input, "ratings")) {
    gwet_ac_ratings_cpp(
      ratings = prep$X,
      n_levels = prep$n_levels,
      weights = W,
      na_code = switch(prep$na_method, error = 1L, complete = 2L, available = 3L),
      min_raters = min_raters,
      by_category = by_category,
      inference_code = inference_code,
      conf_level = conf_level,
      n_threads = n_threads
    )
  } else {
    gwet_ac_counts_cpp(
      counts = prep$counts,
      weights = W,
      by_category = by_category,
      inference_code = inference_code,
      conf_level = conf_level,
      n_threads = n_threads
    )
  }

  out <- .mc_build_gwet_ac_output(
    raw = raw,
    input = input,
    categories = prep$categories,
    weights = W,
    ci = ci,
    p_value = p_value,
    conf_level = conf_level,
    na_method = if (identical(input, "ratings")) prep$na_method else "available",
    min_raters = min_raters,
    se_method = if (isTRUE(inference_requested)) {
      if (!is.null(raw$se_method)) as.character(raw$se_method) else se_method
    } else {
      "none"
    },
    by_category = by_category,
    prep_diagnostics = prep$diagnostics
  )

  inform_if_verbose(
    "Computed {attr(out, 'coefficient', exact = TRUE)} for {attr(out, 'n_items', exact = TRUE)} item{?s} across {length(prep$categories)} categor{?y/ies}.",
    .verbose = verbose
  )
  out
}

.mc_encode_gwet_pair_labels <- function(x_labels, y_labels, levels = NULL) {
  if (is.null(levels)) {
    return(.mc_encode_nominal_pair_labels(x_labels, y_labels))
  }
  if (!is.atomic(levels) || length(levels) < 2L) {
    abort_bad_arg("levels", message = "must be an atomic vector with at least two categories.")
  }
  if (anyNA(levels)) {
    abort_bad_arg("levels", message = "must not contain missing values.")
  }
  level_keys <- as.character(levels)
  if (anyDuplicated(level_keys)) {
    abort_bad_arg("levels", message = "must not contain duplicate categories.")
  }
  observed <- c(x_labels[!is.na(x_labels)], y_labels[!is.na(y_labels)])
  if (length(observed) && any(!observed %in% level_keys)) {
    abort_bad_arg("levels", message = "must contain every non-missing observed category.")
  }
  list(
    x = .mc_match_nominal_levels(x_labels, level_keys),
    y = .mc_match_nominal_levels(y_labels, level_keys),
    levels = levels,
    n_levels = length(level_keys)
  )
}

.mc_extract_gwet_pairwise_columns <- function(data, levels = NULL) {
  if (is.null(levels)) {
    return(.mc_extract_nominal_columns(data, arg = "data", min_cols = 2L))
  }

  if (!is.data.frame(data) && !is.matrix(data)) {
    abort_bad_arg("data", message = "must be a matrix or data frame in matrix mode.")
  }
  if (!is.atomic(levels) || length(levels) < 2L) {
    abort_bad_arg("levels", message = "must be an atomic vector with at least two categories.")
  }
  if (anyNA(levels)) {
    abort_bad_arg("levels", message = "must not contain missing values.")
  }
  level_keys <- as.character(levels)
  if (anyDuplicated(level_keys)) {
    abort_bad_arg("levels", message = "must not contain duplicate categories.")
  }

  cols <- if (is.data.frame(data)) unclass(data) else lapply(seq_len(ncol(data)), function(j) data[, j])
  col_names <- if (is.data.frame(data)) names(data) else colnames(data)
  if (length(cols) < 2L) {
    abort_bad_arg("data", message = "must contain at least two nominal columns.")
  }
  labels <- vector("list", length(cols))
  for (j in seq_along(cols)) {
    col_arg <- if (!is.null(col_names) && nzchar(col_names[[j]] %||% "")) {
      sprintf("data$%s", col_names[[j]])
    } else {
      sprintf("data[[%d]]", j)
    }
    labels[[j]] <- .mc_nominal_vector_labels(cols[[j]], arg = col_arg)
  }
  observed <- unlist(lapply(labels, function(z) z[!is.na(z)]), use.names = FALSE)
  if (length(observed) && any(!observed %in% level_keys)) {
    abort_bad_arg("levels", message = "must contain every non-missing observed category.")
  }

  X <- matrix(
    NA_integer_,
    nrow = length(labels[[1L]]),
    ncol = length(labels),
    dimnames = list(NULL, col_names)
  )
  for (j in seq_along(labels)) {
    X[, j] <- .mc_match_nominal_levels(labels[[j]], level_keys)
  }
  list(
    matrix = X,
    levels = levels,
    n_levels = length(level_keys),
    names = col_names
  )
}

.mc_gwet_weight_support <- function(categories) {
  if (is.numeric(categories) || is.integer(categories)) {
    return(as.numeric(categories))
  }
  seq_along(categories)
}

.mc_resolve_gwet_weights <- function(weights, categories) {
  n_levels <- length(categories)
  check_scalar_numeric(n_levels, arg = "categories", lower = 1, closed_lower = TRUE)

  if (is.character(weights)) {
    if (length(weights) != 1L) {
      weights <- weights[[1L]]
    }
    key <- tolower(weights[[1L]])
    alias_map <- c(
      unweighted = "unweighted",
      identity = "unweighted",
      linear = "linear",
      quadratic = "quadratic",
      ordinal = "ordinal",
      radical = "radical",
      ratio = "ratio",
      circular = "circular",
      bipolar = "bipolar"
    )
    resolved <- unname(alias_map[key])
    if (!length(resolved) || is.na(resolved)) {
      abort_bad_arg(
        "weights",
        message = paste(
          "must be one of {.val unweighted}, {.val linear}, {.val quadratic},",
          "{.val ordinal}, {.val radical}, {.val ratio}, {.val circular},",
          "{.val bipolar}, or a numeric square matrix."
        )
      )
    }

    K <- as.integer(n_levels)
    if (identical(resolved, "unweighted")) {
      W <- diag(1, nrow = K, ncol = K)
    } else {
      support <- .mc_gwet_weight_support(categories)
      if (!is.numeric(support) || any(!is.finite(support))) {
        support <- seq_len(K)
      }
      if (length(support) != K) {
        support <- seq_len(K)
      }
      xmin <- min(support)
      xmax <- max(support)
      span <- xmax - xmin
      idx_i <- matrix(rep(support, each = K), nrow = K)
      idx_j <- matrix(rep(support, times = K), nrow = K)
      if (identical(resolved, "linear")) {
        W <- if (K <= 1L || span == 0) matrix(1, K, K) else 1 - abs(idx_i - idx_j) / abs(span)
      } else if (identical(resolved, "quadratic")) {
        W <- if (K <= 1L || span == 0) matrix(1, K, K) else 1 - ((idx_i - idx_j)^2) / (span^2)
      } else if (identical(resolved, "radical")) {
        W <- if (K <= 1L || span == 0) matrix(1, K, K) else 1 - sqrt(abs(idx_i - idx_j)) / sqrt(abs(span))
      } else if (identical(resolved, "ratio")) {
        denom <- ((xmax - xmin) / (xmax + xmin))^2
        if (!is.finite(denom) || denom <= 0) {
          support <- seq_len(K)
          idx_i <- matrix(rep(support, each = K), nrow = K)
          idx_j <- matrix(rep(support, times = K), nrow = K)
          denom <- ((max(support) - min(support)) / (max(support) + min(support)))^2
        }
        W <- 1 - (((idx_i - idx_j) / (idx_i + idx_j))^2) / denom
        W[!is.finite(W)] <- 0
      } else if (identical(resolved, "circular")) {
        U <- span + 1
        base <- sin(pi * (idx_i - idx_j) / U)^2
        W <- 1 - base / max(base)
        W[!is.finite(W)] <- 1
      } else if (identical(resolved, "bipolar")) {
        W <- matrix(0, K, K)
        for (i in seq_len(K)) {
          for (j in seq_len(K)) {
            if (i == j) {
              W[i, j] <- 1
            } else {
              denom <- (((support[i] + support[j]) - 2 * xmin) * (2 * xmax - (support[i] + support[j])))
              W[i, j] <- if (!is.finite(denom) || denom == 0) 0 else
                (support[i] - support[j])^2 / denom
            }
          }
        }
        W <- 1 - W / max(W)
        W[!is.finite(W)] <- 1
      } else {
        W <- matrix(0, K, K)
        for (i in seq_len(K)) {
          for (j in seq_len(K)) {
            nkl <- abs(i - j) + 1
            W[i, j] <- nkl * (nkl - 1) / 2
          }
        }
        W <- 1 - W / max(W)
      }
    }
    storage.mode(W) <- "double"
    diag(W) <- 1
    W[W < 0] <- 0
    W[W > 1] <- 1
    attr(W, "weight_type") <- resolved
    attr(W, "coefficient") <- if (identical(resolved, "unweighted")) "AC1" else "AC2"
    return(W)
  }

  if (!is.numeric(weights) || !is.matrix(weights)) {
    abort_bad_arg("weights", message = "must be a supported character scheme or a numeric square matrix.")
  }
  check_matrix_dims(weights, arg = "weights")
  if (nrow(weights) != ncol(weights)) {
    abort_bad_arg("weights", message = "must be a square matrix.")
  }
  if (nrow(weights) != n_levels) {
    abort_bad_arg("weights", message = "must have dimension {n_levels} x {n_levels}.", n_levels = n_levels)
  }
  if (any(!is.finite(weights))) {
    abort_bad_arg("weights", message = "must contain only finite values.")
  }
  if (any(weights < 0 | weights > 1)) {
    abort_bad_arg("weights", message = "must contain values in [0, 1].")
  }
  if (any(abs(diag(weights) - 1) > sqrt(.Machine$double.eps))) {
    abort_bad_arg("weights", message = "must have diagonal entries equal to 1.")
  }
  if (isTRUE(max(abs(weights - t(weights))) > sqrt(.Machine$double.eps))) {
    abort_bad_arg("weights", message = "must be symmetric.")
  }
  W <- unclass(weights)
  storage.mode(W) <- "double"
  attr(W, "weight_type") <- "custom"
  attr(W, "coefficient") <- if (isTRUE(all.equal(W, diag(n_levels), tolerance = sqrt(.Machine$double.eps)))) "AC1" else "AC2"
  W
}

.mc_gwet_ac_pairwise_matrix <- function(data,
                                        weights,
                                        levels,
                                        na_method,
                                        ci,
                                        p_value,
                                        conf_level,
                                        n_threads,
                                        output_cfg) {
  enc <- .mc_extract_gwet_pairwise_columns(data, levels = levels)
  X <- enc$matrix
  col_names <- enc$names
  dn <- .mc_square_dimnames(col_names)
  W <- .mc_resolve_gwet_weights(weights, enc$levels)

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

  if (direct_threshold) {
    trip <- gwet_ac_threshold_triplets_cpp(
      X,
      weights = W,
      drop_unused_levels = is.null(levels),
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
      estimator_class = "gwet_ac",
      method = "gwet_ac",
      description = sprintf("Pairwise Gwet %s agreement matrix", attr(W, "coefficient", exact = TRUE)),
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

  fit <- gwet_ac_matrix_cpp(
    X,
    weights = W,
    drop_unused_levels = is.null(levels),
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
      expected_agreement = fit$pe,
      weights = W,
      weight_type = attr(W, "weight_type", exact = TRUE)
    ),
    if (identical(na_method, "pairwise")) {
      list(na_method = "pairwise", common_sample = FALSE)
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
      ci.method = "asymptotic"
    )
  } else {
    NULL
  }
  inference_attr <- if (isTRUE(p_value)) {
    list(
      method = "t_asymptotic_gwet",
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
    estimator_class = "gwet_ac",
    method = "gwet_ac",
    description = sprintf("Pairwise Gwet %s agreement matrix", attr(W, "coefficient", exact = TRUE)),
    diagnostics = diagnostics,
    ci = ci_attr,
    conf.level = if (isTRUE(ci)) conf_level else NULL,
    symmetric = TRUE,
    extra_attrs = list(
      inference = inference_attr,
      weights = W,
      weight_type = attr(W, "weight_type", exact = TRUE),
      coefficient = attr(W, "coefficient", exact = TRUE),
      levels = enc$levels
    )
  )

  .mc_finalize_corr_output_fast(
    out,
    output = output_cfg$output,
    threshold = output_cfg$threshold,
    diag = output_cfg$diag
  )
}

.mc_gwet_ac_scalar <- function(fit,
                               weights,
                               levels,
                               ci = FALSE,
                               p_value = FALSE,
                               conf_level = 0.95) {
  est <- as.numeric(fit$est)
  out <- structure(est, class = c("gwet_ac", "numeric"))
  attr(out, "method") <- "gwet_ac"
  attr(out, "description") <- sprintf("Pairwise Gwet %s agreement", attr(weights, "coefficient", exact = TRUE))
  attr(out, "package") <- "matrixCorr"
  attr(out, "estimate") <- est
  attr(out, "levels") <- levels
  attr(out, "weights") <- weights
  attr(out, "weight_type") <- attr(weights, "weight_type", exact = TRUE)
  attr(out, "coefficient") <- attr(weights, "coefficient", exact = TRUE)
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
      ci.method = "asymptotic"
    )
    attr(out, "conf.level") <- conf_level
    attr(out, "conf_level") <- conf_level
  }
  if (isTRUE(p_value)) {
    attr(out, "inference") <- list(
      method = "t_asymptotic_gwet",
      estimate = est,
      se = as.numeric(fit$se),
      statistic = as.numeric(fit$statistic),
      p_value = as.numeric(fit$p_value),
      n_obs = as.integer(fit$n_complete)
    )
  }
  out
}

.mc_build_gwet_ac_output <- function(raw,
                                     input,
                                     categories,
                                     weights,
                                     ci,
                                     p_value,
                                     conf_level,
                                     na_method,
                                     min_raters,
                                     se_method,
                                     by_category,
                                     prep_diagnostics) {
  coefficient <- attr(weights, "coefficient", exact = TRUE)
  inference_requested <- isTRUE(ci) || isTRUE(p_value)
  out <- data.frame(
    method = "gwet_ac",
    ac = as.numeric(raw$estimate),
    observed_agreement = as.numeric(raw$observed_agreement),
    expected_agreement = as.numeric(raw$expected_agreement),
    n_items = as.integer(raw$n_items),
    n_categories = as.integer(raw$n_categories),
    n_raters_min = as.integer(raw$n_raters_min),
    n_raters_max = as.integer(raw$n_raters_max),
    balanced = as.logical(raw$balanced),
    n_ratings_total = as.integer(raw$n_ratings_total),
    se = if (!is.null(raw$se)) as.numeric(raw$se) else NA_real_,
    statistic = if (isTRUE(p_value) && !is.null(raw$statistic)) as.numeric(raw$statistic) else NA_real_,
    p_value = if (isTRUE(p_value) && !is.null(raw$p_value)) as.numeric(raw$p_value) else NA_real_,
    lwr.ci = if (isTRUE(ci) && !is.null(raw$lwr)) as.numeric(raw$lwr) else NA_real_,
    upr.ci = if (isTRUE(ci) && !is.null(raw$upr)) as.numeric(raw$upr) else NA_real_,
    conf.level = if (isTRUE(ci)) conf_level else NA_real_,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  by_category_df <- NULL
  if (isTRUE(by_category)) {
    by_category_df <- data.frame(
      category = as.character(categories),
      ac = as.numeric(raw$by_category_estimate),
      observed_agreement = as.numeric(raw$by_category_observed_agreement),
      expected_agreement = as.numeric(raw$by_category_expected_agreement),
      se = if (!is.null(raw$by_category_se)) as.numeric(raw$by_category_se) else NA_real_,
      statistic = if (isTRUE(p_value) && !is.null(raw$by_category_statistic)) as.numeric(raw$by_category_statistic) else NA_real_,
      p_value = if (isTRUE(p_value) && !is.null(raw$by_category_p_value)) as.numeric(raw$by_category_p_value) else NA_real_,
      lwr.ci = if (isTRUE(ci) && !is.null(raw$by_category_lwr)) as.numeric(raw$by_category_lwr) else NA_real_,
      upr.ci = if (isTRUE(ci) && !is.null(raw$by_category_upr)) as.numeric(raw$by_category_upr) else NA_real_,
      n_items = as.integer(raw$n_items),
      n_ratings_total = as.integer(raw$n_ratings_total),
      category_count = as.integer(raw$category_counts),
      category_proportion = as.numeric(raw$category_proportions),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  class(out) <- c("gwet_ac", "agreement_result", "data.frame")
  attr(out, "method") <- "gwet_ac"
  attr(out, "description") <- sprintf("Gwet %s chance-corrected agreement", coefficient)
  attr(out, "package") <- "matrixCorr"
  attr(out, "input") <- input
  attr(out, "categories") <- categories
  attr(out, "n_categories") <- length(categories)
  attr(out, "n_items") <- as.integer(raw$n_items)
  attr(out, "n_raters") <- list(
    min = as.integer(raw$n_raters_min),
    max = as.integer(raw$n_raters_max),
    balanced = as.logical(raw$balanced)
  )
  attr(out, "na_method") <- na_method
  attr(out, "min_raters") <- min_raters
  attr(out, "weights") <- weights
  attr(out, "weight_type") <- attr(weights, "weight_type", exact = TRUE)
  attr(out, "coefficient") <- coefficient
  attr(out, "diagnostics") <- c(
    list(
      category_counts = as.integer(raw$category_counts),
      category_proportions = as.numeric(raw$category_proportions),
      item_category_counts = {
        mat <- raw$item_category_counts
        if (!is.null(mat)) {
          mat <- as.matrix(mat)
          colnames(mat) <- as.character(categories)
          rownames(mat) <- paste0("item_", seq_len(nrow(mat)))
        }
        mat
      },
      item_agreement = as.numeric(raw$item_agreement),
      n_raters_by_item = as.integer(raw$n_raters_by_item)
    ),
    prep_diagnostics
  )
  attr(out, "by_category") <- by_category_df
  attr(out, "se.method") <- if (isTRUE(inference_requested)) se_method else "none"
  attr(out, "conf.level") <- if (isTRUE(ci)) conf_level else NULL
  attr(out, "reference") <- "Gwet (2008), British Journal of Mathematical and Statistical Psychology, 61, 29-48, doi:10.1348/000711006X126600"
  out
}

.mc_as_multirater_plot_input <- function(x) {
  tmp <- x
  tmp$kappa <- tmp$ac
  tmp$ac <- NULL
  by_cat <- attr(tmp, "by_category", exact = TRUE)
  if (is.data.frame(by_cat) && "ac" %in% names(by_cat)) {
    by_cat$kappa <- by_cat$ac
    by_cat$ac <- NULL
    attr(tmp, "by_category") <- by_cat
  }
  class(tmp) <- c("multirater_kappa", "agreement_result", "data.frame")
  attr(tmp, "kappa_method") <- attr(x, "coefficient", exact = TRUE) %||% "gwet_ac"
  attr(tmp, "exact") <- FALSE
  tmp
}

#' @rdname gwet_ac
#' @method summary gwet_ac
#' @param object A \code{gwet_ac} object.
#' @param digits Integer; number of decimal places for displayed values.
#' @param ci_digits Integer; number of decimal places for confidence limits.
#' @param p_digits Integer; number of decimal places for p-values.
#' @param ... Unused.
#' @export
summary.gwet_ac <- function(object,
                            digits = 4,
                            ci_digits = 3,
                            p_digits = 4,
                            ...) {
  if (inherits(object, "corr_result")) {
    return(summary.corr_matrix(object, ...))
  }
  if (inherits(object, "agreement_result") && is.data.frame(object)) {
    out <- as.data.frame(object, stringsAsFactors = FALSE, check.names = FALSE)
    preferred_order <- c(
      "method", "ac", "observed_agreement", "expected_agreement", "n_items",
      "n_categories", "n_raters_min", "n_raters_max", "balanced",
      "n_ratings_total", "se", "statistic", "p_value", "lwr.ci", "upr.ci",
      "conf.level"
    )
    out <- out[, c(intersect(preferred_order, names(out)), setdiff(names(out), preferred_order)), drop = FALSE]
    drop_if_all_na <- c("se", "statistic", "p_value", "lwr.ci", "upr.ci", "conf.level")
    keep_cols <- vapply(names(out), function(nm) {
      if (!nm %in% drop_if_all_na) {
        return(TRUE)
      }
      !all(is.na(out[[nm]]))
    }, logical(1))
    out <- out[, keep_cols, drop = FALSE]
    num_cols <- vapply(out, is.numeric, logical(1))
    out[num_cols] <- lapply(out[num_cols], round, digits = digits)
    out <- .mc_finalize_summary_df(out, class_name = "summary.gwet_ac")
    attr(out, "result_type") <- "multirater"
    attr(out, "digits") <- digits
    attr(out, "by_category") <- attr(object, "by_category", exact = TRUE)
    attr(out, "se.method") <- attr(object, "se.method", exact = TRUE)
    attr(out, "coefficient") <- attr(object, "coefficient", exact = TRUE)
    return(out)
  }

  check_inherits(object, "gwet_ac")
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
  out <- .mc_finalize_summary_df(df, class_name = "summary.gwet_ac")
  attr(out, "result_type") <- "scalar"
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
  attr(out, "coefficient") <- attr(object, "coefficient", exact = TRUE)
  out
}

#' @rdname gwet_ac
#' @method print summary.gwet_ac
#' @param x A \code{summary.gwet_ac} object.
#' @param n Optional preview row threshold.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns.
#' @param width Optional display width.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @export
print.summary.gwet_ac <- function(x,
                                  digits = NULL,
                                  n = NULL,
                                  topn = NULL,
                                  max_vars = NULL,
                                  width = NULL,
                                  show_ci = NULL,
                                  ...) {
  coefficient <- attr(x, "coefficient", exact = TRUE) %||% "AC"
  if (identical(attr(x, "result_type", exact = TRUE), "multirater")) {
    digits <- .mc_coalesce(digits, attr(x, "digits", exact = TRUE) %||% 4L)
    digest <- c(
      coefficient = coefficient,
      ac = formatC(x$ac[[1L]], format = "f", digits = digits),
      items = .mc_count_fmt(x$n_items[[1L]]),
      categories = .mc_count_fmt(x$n_categories[[1L]])
    )
    se_method <- attr(x, "se.method", exact = TRUE) %||% "none"
    if (!identical(se_method, "none")) {
      digest <- c(digest, se_method = se_method)
    }
    .mc_print_named_digest(digest, header = "Gwet agreement summary")
    cat("\n")
    .mc_print_summary_table(
      x,
      header = "Estimate table",
      digits = digits,
      n = n,
      topn = topn,
      max_vars = max_vars,
      width = width,
      show_ci = show_ci,
      print_overview = FALSE,
      ...
    )
    return(invisible(x))
  }

  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )
  digits <- .mc_coalesce(digits, .mc_coalesce(attr(x, "digits", exact = TRUE), 4L))
  ci_digits <- .mc_coalesce(attr(x, "ci_digits", exact = TRUE), digits)
  p_digits <- .mc_coalesce(attr(x, "p_digits", exact = TRUE), digits)

  digest <- c(
    coefficient = coefficient,
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
  .mc_print_named_digest(digest, header = "Gwet agreement summary")
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

#' @rdname gwet_ac
#' @method print gwet_ac
#' @param x A \code{gwet_ac} object.
#' @param show_by_category Logical; whether to print attached category-wise
#'   panel results when available.
#' @export
print.gwet_ac <- function(x,
                          digits = 4,
                          n = NULL,
                          topn = NULL,
                          max_vars = NULL,
                          width = NULL,
                          show_ci = NULL,
                          show_by_category = FALSE,
                          ...) {
  coefficient <- attr(x, "coefficient", exact = TRUE) %||% "AC"
  if (inherits(x, "corr_result")) {
    .mc_print_corr_matrix(
      x,
      header = sprintf("Gwet %s agreement matrix", coefficient),
      digits = digits,
      n = n,
      topn = topn,
      max_vars = max_vars,
      width = width,
      show_ci = show_ci,
      ...
    )
    return(invisible(x))
  }

  if (inherits(x, "agreement_result") && is.data.frame(x)) {
    check_bool(show_by_category, arg = "show_by_category")
    cat("Gwet agreement\n")
    cat("  coefficient        :", coefficient, "\n")
    cat("  estimate           :", formatC(x$ac[[1L]], format = "f", digits = digits), "\n")
    cat("  observed agreement :", formatC(x$observed_agreement[[1L]], format = "f", digits = digits), "\n")
    cat("  expected agreement :", formatC(x$expected_agreement[[1L]], format = "f", digits = digits), "\n")
    cat("  items              :", .mc_count_fmt(x$n_items[[1L]]), "\n")
    cat("  categories         :", .mc_count_fmt(x$n_categories[[1L]]), "\n")
    nr <- attr(x, "n_raters", exact = TRUE)
    cat(
      "  raters per item    :",
      paste0(.mc_count_fmt(nr$min), " to ", .mc_count_fmt(nr$max)),
      if (isTRUE(nr$balanced)) "(balanced)" else "(unbalanced)",
      "\n"
    )
    if (is.finite(x$lwr.ci[[1L]]) && is.finite(x$upr.ci[[1L]])) {
      cat(
        "  CI                 :",
        sprintf(
          "%g%% [%s, %s]",
          100 * suppressWarnings(as.numeric(attr(x, "conf.level", exact = TRUE))),
          formatC(x$lwr.ci[[1L]], format = "f", digits = digits),
          formatC(x$upr.ci[[1L]], format = "f", digits = digits)
        ),
        "\n"
      )
    }
    if (is.finite(x$p_value[[1L]])) {
      cat("  p-value            :", formatC(x$p_value[[1L]], format = "f", digits = digits), "\n")
    }
    if (isTRUE(show_by_category)) {
      by_cat <- attr(x, "by_category", exact = TRUE)
      if (is.data.frame(by_cat) && nrow(by_cat)) {
        cat("\nBy-category AC1\n\n")
        print.data.frame(by_cat, row.names = FALSE, ...)
      }
    }
    return(invisible(x))
  }

  invisible(print.summary.gwet_ac(
    summary.gwet_ac(x),
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  ))
}

#' @rdname gwet_ac
#' @method plot gwet_ac
#' @param title Optional plot title.
#' @param low_color Fill/color used for negative agreement.
#' @param high_color Fill/color used for positive agreement.
#' @param mid_color Unused placeholder for API consistency with matrix heatmaps.
#' @param value_text_size Text size for estimate labels in matrix plots.
#' @param ci_text_size Text size for CI labels in matrix plots.
#' @param show_value Logical; whether to print estimate and CI labels.
#' @param type Plot type for panel-level fits: \code{"agreement_map"},
#'   \code{"estimate"}, \code{"item_agreement"}, \code{"category_proportion"},
#'   or \code{"by_category"}.
#' @param bins Integer number of bins retained for compatibility with the
#'   panel-level item-agreement plot.
#' @export
plot.gwet_ac <- function(x,
                         title = NULL,
                         low_color = "indianred1",
                         high_color = "steelblue1",
                         mid_color = "white",
                         value_text_size = 4,
                         ci_text_size = 3,
                         show_value = TRUE,
                         type = c("agreement_map", "estimate", "item_agreement", "category_proportion", "by_category"),
                         bins = 30L,
                         ...) {
  coefficient <- attr(x, "coefficient", exact = TRUE) %||% "AC"
  if (inherits(x, "corr_result")) {
    return(.mc_plot_corr_result(
      x,
      title = title %||% sprintf("Gwet %s agreement heatmap", coefficient),
      low_color = low_color,
      high_color = high_color,
      mid_color = mid_color,
      value_text_size = value_text_size,
      ci_text_size = ci_text_size,
      show_value = show_value,
      ...
    ))
  }

  if (inherits(x, "agreement_result") && is.data.frame(x)) {
    return(plot.multirater_kappa(
      .mc_as_multirater_plot_input(x),
      type = match.arg(type),
      bins = bins,
      ...
    ))
  }

  check_inherits(x, "gwet_ac")
  check_bool(show_value, arg = "show_value")
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} is required for plotting.")
  }

  ci_attr <- attr(x, "ci", exact = TRUE)
  est <- as.numeric(x)
  df <- data.frame(
    label = coefficient,
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
    ggplot2::scale_y_continuous(breaks = 1, labels = coefficient) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      ...
    ) +
    ggplot2::labs(
      title = title %||% sprintf("Gwet %s agreement", coefficient),
      x = coefficient
    )

  if (isTRUE(show_value)) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$estimate_label),
        nudge_y = 0.18,
        size = 4
      )
    if (is.finite(df$lwr[[1L]]) && is.finite(df$upr[[1L]])) {
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes(label = .data$ci_label),
          nudge_y = -0.18,
          size = 3
        )
    }
  }
  p
}
