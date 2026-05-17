#' @title Multi-Rater Kappa for Nominal Ratings
#'
#' @description
#' Estimates panel-level chance-corrected agreement among multiple raters
#' assigning items to nominal categories.
#'
#' @param data Input ratings or counts data.
#'   For \code{input = "ratings"}, rows are items/subjects and columns are
#'   raters/classifiers. For \code{input = "counts"}, rows are items/subjects,
#'   columns are categories, and values are counts of raters assigning that
#'   category to that item.
#' @param levels Optional explicit category labels. For nominal multi-rater
#'   kappa, category order is not used in the estimator itself, but explicit
#'   levels are useful when unobserved categories should be retained in the
#'   output. If factor columns are mixed with non-factor columns, \code{levels}
#'   must be supplied to avoid ambiguous implicit ordering.
#' @param input One of \code{"ratings"} or \code{"counts"}.
#' @param method One of \code{"fleiss"} (fixed-marginal multi-rater kappa) or
#'   \code{"randolph"} (free-marginal multi-rater kappa).
#' @param exact Logical; if \code{TRUE}, use the exact Fleiss fixed-marginal
#'   variant. This is available only for \code{input = "ratings"} with
#'   complete item-by-rater data and \code{method = "fleiss"}.
#' @param na_method Missing-data rule for \code{input = "ratings"}.
#'   \code{"error"} rejects missing ratings. \code{"complete"} removes any item
#'   with one or more missing ratings. \code{"available"} retains items with at
#'   least \code{min_raters} observed ratings and allows the number of raters to
#'   vary by item. For \code{input = "counts"}, missing counts are not allowed.
#' @param min_raters Minimum number of observed ratings required for an item to
#'   be retained when \code{na_method = "available"} or when
#'   \code{input = "counts"} rows are filtered by total raters. Must be at
#'   least \code{2}.
#' @param by_category Logical; if \code{TRUE}, attach category-wise kappas based
#'   on target-vs-other binary collapses of the category count matrix.
#' @param ci Logical; if \code{TRUE}, attach confidence intervals using the
#'   selected \code{se_method}.
#' @param p_value Logical; if \code{TRUE}, attach z statistics and two-sided
#'   p-values using the selected \code{se_method}.
#' @param conf_level Confidence level for intervals. Default is \code{0.95}.
#' @param se_method Standard-error method. \code{"asymptotic"} uses the
#'   closed-form large-sample variance for the standard non-exact Fleiss kappa
#'   when a common number of raters is available. \code{"jackknife"} uses
#'   leave-one-item-out jackknife inference. \code{"none"} disables
#'   inferential quantities and may not be combined with \code{ci = TRUE} or
#'   \code{p_value = TRUE}. When inferential quantities are requested and
#'   \code{se_method} is left at its default, \code{multirater_kappa()}
#'   automatically falls back to \code{"jackknife"} if asymptotic inference is
#'   unavailable for the chosen estimator/data combination. \code{se_method}
#'   and \code{conf_level} are validated only when \code{ci} or \code{p_value}
#'   requests inferential output.
#' @param n_threads Integer \eqn{\ge 1}. Number of OpenMP threads used by the
#'   C++ backend.
#' @param verbose Logical; if \code{TRUE}, emit short progress messages.
#' @param ... Reserved for future extensions. Unsupported extra arguments are
#'   rejected.
#'
#' @details
#' `multirater_kappa()` returns one panel-level agreement coefficient, not a
#' pairwise matrix. Fleiss' kappa is the default fixed-marginal estimator:
#' \deqn{\kappa = \frac{\bar P - P_e}{1 - P_e},}
#' where \eqn{\bar P} is the mean item-level agreement and
#' \eqn{P_e = \sum_j p_j^2} uses the pooled marginal category proportions.
#'
#' Randolph's free-marginal alternative replaces the expected agreement with
#' \eqn{P_e = 1 / K}. This function is for nominal categories and does not use
#' category ordering. Use [weighted_kappa()] for two-rater ordered-category
#' agreement and [cohen_kappa()] for two-rater nominal agreement.
#'
#' When \code{exact = TRUE}, the exact Fleiss fixed-marginal estimate requires
#' the original item-by-rater rating matrix. In \pkg{matrixCorr}, closed-form
#' asymptotic inference is only available for the standard non-exact Fleiss
#' estimator with a common number of raters per item. In other settings, use
#' \code{se_method = "jackknife"} when inferential quantities are needed.
#'
#' With \code{input = "ratings"} and \code{na_method = "available"}, the
#' implementation supports unbalanced item-specific numbers of raters as a
#' generalisation beyond the strict equal-rater Fleiss setting. The returned
#' diagnostics record the observed per-item rater counts.
#'
#' \strong{Confidence intervals and standard errors.}
#' Two inference paths are implemented.
#'
#' If \code{se_method = "jackknife"}, the method computes leave-one-item-out
#' estimates \eqn{\hat\kappa_{(-1)},\ldots,\hat\kappa_{(-m)}} over the
#' \eqn{m} retained items and forms
#' \deqn{ \bar\kappa_{(.)} = \frac{1}{m}\sum_{i=1}^m \hat\kappa_{(-i)}, }
#' \deqn{ \widehat{\mathrm{Var}}_{\mathrm{JK}}(\hat\kappa)
#'       \;=\;
#'       \frac{m-1}{m}\sum_{i=1}^m
#'       \left(\hat\kappa_{(-i)} - \bar\kappa_{(.)}\right)^2. }
#' The standard error is
#' \eqn{\widehat{\mathrm{se}}(\hat\kappa)=
#' \sqrt{\widehat{\mathrm{Var}}_{\mathrm{JK}}(\hat\kappa)}} and the CI is
#' \deqn{ \hat\kappa \pm z_{1-\alpha/2}\,\widehat{\mathrm{se}}(\hat\kappa), }
#' truncated to \eqn{[-1, 1]}.
#'
#' If \code{se_method = "asymptotic"}, the closed-form variance is available
#' only for the standard non-exact Fleiss estimator
#' with a common number of raters per item. Let \eqn{m} be the number of items,
#' \eqn{n} the common number of raters, \eqn{p_j} the pooled category
#' proportions, and define
#' \deqn{ S = \sum_j p_j(1-p_j), \qquad
#'        C = \sum_j p_j(1-p_j)(1-2p_j). }
#' Then the variance estimator is
#' \deqn{ \widehat{\mathrm{Var}}(\hat\kappa)
#'       \;=\;
#'       \frac{2\,(S^2 - C)}{S^2\,m\,n\,(n-1)}. }
#' The standard error and CI are then obtained from the same Wald form
#' above and truncated to \eqn{[-1, 1]}. When asymptotic inference is
#' unavailable and inferential output is requested with the default setting,
#' \code{multirater_kappa()} automatically falls back to the jackknife path.
#'
#' @return
#' A one-row data frame with class
#' \code{c("multirater_kappa", "agreement_result", "data.frame")}. The object
#' carries estimator metadata and diagnostics in attributes, and may also
#' attach category-wise binary-collapsed results in
#' \code{attr(x, "by_category")}.
#'
#' @references
#' Fleiss, J. L. (1971). Measuring nominal scale agreement among many raters.
#' \emph{Psychological Bulletin}, 76, 378-382. \doi{10.1037/h0031619}
#'
#' Randolph, J. J. (2005). Free-Marginal Multirater Kappa: An Alternative to
#' Fleiss' Fixed-Marginal Multirater Kappa.
#'
#' @seealso
#' [cohen_kappa()] for unweighted two-rater nominal agreement;
#' [weighted_kappa()] for two-rater ordered-category agreement.
#'
#' @examples
#' raters <- data.frame(
#'   r1 = c("A", "A", "B", "C", "A", "B"),
#'   r2 = c("A", "B", "B", "C", "A", "B"),
#'   r3 = c("A", "A", "B", "B", "A", "C"),
#'   stringsAsFactors = FALSE
#' )
#'
#' fit <- multirater_kappa(raters)
#' print(fit)
#' summary(fit)
#' # The default plot is an item-by-category agreement map.
#' # Rows are items, ordered from stronger to weaker item-level agreement.
#' # Columns are categories. Each tile shows how many raters assigned that
#' # category to that item, and darker fill means a larger share of raters
#' # chose that category for that item.
#' plot(fit)
#'
#' @author Thiago de Paula Oliveira
#' @export
multirater_kappa <- function(data,
                             levels = NULL,
                             input = c("ratings", "counts"),
                             method = c("fleiss", "randolph"),
                             exact = FALSE,
                             na_method = c("error", "complete", "available"),
                             min_raters = 2L,
                             by_category = FALSE,
                             ci = FALSE,
                             p_value = FALSE,
                             conf_level = 0.95,
                             se_method = c("asymptotic", "jackknife", "none"),
                             n_threads = getOption("matrixCorr.threads", 1L),
                             verbose = FALSE,
                             ...) {
  .mc_extract_legacy_aliases(list(...), allowed = character())

  input <- match_arg(input, c("ratings", "counts"), arg_name = "input")
  method <- match_arg(method, c("fleiss", "randolph"), arg_name = "method")
  check_bool(exact, arg = "exact")
  na_method <- if (missing(na_method)) {
    "error"
  } else {
    validate_na_method(
      na_method,
      arg = "na_method",
      allowed = c("error", "complete", "available", "pairwise")
    )
  }
  if (identical(na_method, "pairwise")) {
    abort_bad_arg(
      "na_method",
      message = "na_method = 'pairwise' is not defined for panel-level multi-rater kappa; use 'available' to allow item-specific missing ratings."
    )
  }
  if (isTRUE(exact) && !identical(method, "fleiss")) {
    abort_bad_arg(
      "exact",
      message = "must be FALSE unless {.arg method} is {.val fleiss}."
    )
  }
  if (isTRUE(exact) && !identical(input, "ratings")) {
    abort_bad_arg(
      "exact",
      message = "must be FALSE when {.arg input} is {.val counts}.",
      .hint = "The exact Fleiss variant requires the original item-by-rater ratings because aggregated counts do not identify rater-specific marginals."
    )
  }
  if (isTRUE(exact) && identical(na_method, "available")) {
    abort_bad_arg(
      "exact",
      message = "must be FALSE when {.arg na_method} is {.val available}.",
      .hint = "The exact Fleiss variant is defined here only for complete item-by-rater panels."
    )
  }

  min_raters <- check_scalar_int_pos(min_raters, arg = "min_raters")
  if (min_raters < 2L) {
    abort_bad_arg("min_raters", message = "must be >= 2.")
  }
  check_bool(by_category, arg = "by_category")
  check_bool(ci, arg = "ci")
  check_bool(p_value, arg = "p_value")
  check_bool(verbose, arg = "verbose")
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")

  inference_requested <- isTRUE(ci) || isTRUE(p_value)
  se_method_missing <- missing(se_method)
  if (isTRUE(inference_requested)) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
    se_method <- match_arg(se_method, c("asymptotic", "jackknife", "none"), arg_name = "se_method")
    if (identical(se_method, "none")) {
      abort_bad_arg(
        "se_method",
        message = "must not be 'none' when ci or p_value is TRUE."
      )
    }
  } else {
    se_method <- "none"
    conf_level <- 0.95
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

  method_code <- switch(method, fleiss = 1L, randolph = 2L)
  inference_code <- if (isTRUE(inference_requested)) {
    switch(se_method, none = 0L, asymptotic = 1L, jackknife = 2L)
  } else {
    0L
  }

  prev_threads <- .mc_prepare_omp_threads(
    n_threads,
    n_threads_missing = missing(n_threads)
  )
  if (!is.null(prev_threads)) {
    on.exit(.mc_exit_omp_threads(prev_threads), add = TRUE)
  }

  run_backend <- function(inference_code) {
    if (identical(input, "ratings")) {
      multirater_kappa_ratings_cpp(
        ratings = prep$X,
        n_levels = prep$n_levels,
        method_code = method_code,
        na_code = switch(prep$na_method, error = 1L, complete = 2L, available = 3L),
        min_raters = min_raters,
        by_category = by_category,
        inference_code = inference_code,
        conf_level = conf_level,
        n_threads = n_threads,
        exact = exact
      )
    } else {
      multirater_kappa_counts_cpp(
        counts = prep$counts,
        method_code = method_code,
        by_category = by_category,
        inference_code = inference_code,
        conf_level = conf_level,
        n_threads = n_threads
      )
    }
  }

  raw <- run_backend(inference_code)
  effective_se_method <- if (isTRUE(inference_requested)) se_method else "none"

  if (isTRUE(inference_requested) &&
      identical(se_method, "asymptotic") &&
      !isTRUE(raw$asymptotic_available)) {
    if (isTRUE(se_method_missing)) {
      raw <- run_backend(2L)
      effective_se_method <- "jackknife"
      inform_if_verbose(
        "Asymptotic inference is unavailable for this estimator/data combination; using {.val jackknife} instead.",
        .verbose = verbose
      )
    } else {
      abort_bad_arg(
        "se_method",
        message = "{.val asymptotic} is not available for this estimator/data combination.",
        .hint = "Use {.arg se_method} = {.val jackknife} for exact Fleiss, Randolph, or unbalanced multi-rater panels."
      )
    }
  }

  out <- .mc_build_multirater_kappa_output(
    raw = raw,
    method = method,
    input = input,
    categories = prep$categories,
    ci = ci,
    p_value = p_value,
    conf_level = conf_level,
    na_method = if (identical(input, "ratings")) prep$na_method else na_method,
    min_raters = min_raters,
    se_method = effective_se_method,
    by_category = by_category,
    prep_diagnostics = prep$diagnostics
  )

  inform_if_verbose(
    "Computed {.val {method}}{if (isTRUE(exact)) ' exact' else ''} multi-rater kappa for {attr(out, 'n_items', exact = TRUE)} item{?s} across {length(prep$categories)} categor{?y/ies}.",
    .verbose = verbose
  )
  out
}

.mc_prepare_multirater_kappa_ratings <- function(data,
                                                 levels,
                                                 na_method,
                                                 min_raters) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    abort_bad_arg("data", message = "must be a matrix or data frame when {.arg input} is {.val ratings}.")
  }

  df <- if (is.data.frame(data)) data else as.data.frame(data, stringsAsFactors = FALSE)
  if (ncol(df) < 2L) {
    abort_bad_arg("data", message = "must contain at least two rater columns.")
  }

  cols <- lapply(seq_along(df), function(j) .mc_multirater_rating_column(df[[j]], arg = names(df)[[j]] %||% "data"))
  names(cols) <- names(df)
  categories <- .mc_resolve_multirater_levels(cols, levels = levels)
  if (length(categories) < 2L) {
    abort_bad_arg("levels", message = "must define at least two categories for multi-rater kappa.")
  }

  level_keys <- as.character(categories)
  X <- matrix(
    NA_integer_,
    nrow = nrow(df),
    ncol = ncol(df),
    dimnames = list(NULL, names(df))
  )

  for (j in seq_along(cols)) {
    vals <- cols[[j]]
    keys <- as.character(vals)
    codes <- match(keys, level_keys)
    codes[is.na(vals)] <- NA_integer_
    if (any(!is.na(vals) & is.na(codes))) {
      abort_bad_arg(
        "data",
        message = "contains values not present in {.arg levels}."
      )
    }
    X[, j] <- as.integer(codes)
  }

  keep <- rep(TRUE, nrow(X))
  if (identical(na_method, "error") && anyNA(X)) {
    abort_bad_arg(
      "data",
      message = "contains missing ratings when {.arg na_method} is {.val error}."
    )
  }
  if (identical(na_method, "complete")) {
    keep <- stats::complete.cases(X)
    X <- X[keep, , drop = FALSE]
  } else if (identical(na_method, "available")) {
    keep <- rowSums(!is.na(X)) >= min_raters
    X <- X[keep, , drop = FALSE]
  }

  if (nrow(X) < 2L) {
    abort_bad_arg(
      "data",
      message = "must retain at least two items after applying missing-data handling."
    )
  }

  list(
    X = X,
    categories = categories,
    n_levels = length(categories),
    retained_rows = keep,
    dropped_rows = !keep,
    na_method = na_method,
    diagnostics = list(
      n_original = nrow(df),
      n_retained = nrow(X),
      retained_rows = keep,
      dropped_rows = !keep,
      n_available_by_item = rowSums(!is.na(if (is.data.frame(data)) as.matrix(data) else data))
    )
  )
}

.mc_prepare_multirater_kappa_counts <- function(data,
                                                levels,
                                                min_raters) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    abort_bad_arg("data", message = "must be a matrix or data frame when {.arg input} is {.val counts}.")
  }

  X <- if (is.matrix(data)) data else as.matrix(data)
  if (!is.numeric(X)) {
    abort_bad_arg("data", message = "must contain numeric count columns when {.arg input} is {.val counts}.")
  }
  if (ncol(X) < 2L) {
    abort_bad_arg("data", message = "must contain at least two category columns.")
  }
  if (anyNA(X)) {
    abort_bad_arg("data", message = "must not contain missing counts.")
  }
  if (any(!is.finite(X))) {
    abort_bad_arg("data", message = "must contain only finite counts.")
  }
  if (any(X < 0)) {
    abort_bad_arg("data", message = "must contain non-negative counts.")
  }
  if (any(abs(X - round(X)) > sqrt(.Machine$double.eps))) {
    abort_bad_arg("data", message = "must contain integer-like count values.")
  }

  categories <- if (!is.null(levels)) {
    if (length(levels) != ncol(X)) {
      abort_bad_arg("levels", message = "must have length equal to ncol(data) for {.arg input} = {.val counts}.")
    }
    levels
  } else if (!is.null(colnames(X))) {
    colnames(X)
  } else {
    paste0("category_", seq_len(ncol(X)))
  }

  if (length(categories) < 2L) {
    abort_bad_arg("levels", message = "must define at least two categories for multi-rater kappa.")
  }

  row_sums <- rowSums(X)
  keep <- row_sums >= min_raters
  counts <- matrix(
    as.integer(round(X[keep, , drop = FALSE])),
    ncol = ncol(X),
    dimnames = list(NULL, as.character(categories))
  )

  if (nrow(counts) < 2L) {
    abort_bad_arg(
      "data",
      message = "must retain at least two items with at least {.val {min_raters}} raters."
    )
  }

  list(
    counts = counts,
    categories = categories,
    n_levels = length(categories),
    retained_rows = keep,
    dropped_rows = !keep,
    na_method = "available",
    diagnostics = list(
      n_original = nrow(X),
      n_retained = nrow(counts),
      retained_rows = keep,
      dropped_rows = !keep,
      row_sums = row_sums
    )
  )
}

.mc_resolve_multirater_levels <- function(data, levels) {
  cols <- if (is.list(data) && !is.data.frame(data)) data else as.list(data)

  if (!is.null(levels)) {
    if (!is.atomic(levels) || length(levels) < 1L) {
      abort_bad_arg("levels", message = "must be an atomic vector of category labels.")
    }
    if (anyNA(levels)) {
      abort_bad_arg("levels", message = "must not contain missing values.")
    }
    keys <- as.character(levels)
    if (anyDuplicated(keys)) {
      abort_bad_arg("levels", message = "must not contain duplicate categories.")
    }
    observed <- unlist(lapply(cols, function(z) as.character(z[!is.na(z)])), use.names = FALSE)
    if (length(observed) && any(!observed %in% keys)) {
      abort_bad_arg("levels", message = "must contain every non-missing observed category.")
    }
    return(levels)
  }

  all_factor <- all(vapply(cols, is.factor, logical(1)))
  any_factor <- any(vapply(cols, is.factor, logical(1)))
  all_numericish <- all(vapply(cols, function(z) is.numeric(z) || is.integer(z) || is.logical(z), logical(1)))

  if (all_factor) {
    return(unique(unlist(lapply(cols, levels), use.names = FALSE)))
  }
  if (any_factor) {
    abort_bad_arg(
      "levels",
      message = "must be supplied when factor and non-factor rating columns are mixed."
    )
  }
  if (all_numericish) {
    observed <- unlist(lapply(cols, function(z) z[!is.na(z)]), use.names = FALSE)
    return(sort(unique(observed)))
  }
  unique(unlist(lapply(cols, function(z) as.character(z[!is.na(z)])), use.names = FALSE))
}

.mc_multirater_rating_column <- function(x, arg) {
  if (is.list(x) && !is.atomic(x)) {
    abort_bad_arg(arg, message = "must not be a list column.")
  }
  if (is.matrix(x) || is.data.frame(x)) {
    abort_bad_arg(arg, message = "must be a vector-like rating column.")
  }
  if (inherits(x, "Date") || inherits(x, "POSIXct") || inherits(x, "POSIXt") || inherits(x, "difftime")) {
    abort_bad_arg(arg, message = "must not be a date-time or time-difference column.")
  }
  if (is.factor(x)) {
    return(x)
  }
  if (is.character(x) || is.logical(x) || is.numeric(x) || is.integer(x)) {
    return(x)
  }
  abort_bad_arg(
    arg,
    message = "must be factor, ordered factor, character, logical, integer, or numeric."
  )
}

.mc_build_multirater_kappa_output <- function(raw,
                                              method,
                                              input,
                                              categories,
                                              ci,
                                              p_value,
                                              conf_level,
                                              na_method,
                                              min_raters,
                                              se_method,
                                              by_category,
                                              prep_diagnostics) {
  inference_requested <- isTRUE(ci) || isTRUE(p_value)
  se <- if (!is.null(raw$se)) as.numeric(raw$se) else NA_real_
  statistic <- if (isTRUE(p_value) && !is.null(raw$statistic)) as.numeric(raw$statistic) else NA_real_
  p_val <- if (isTRUE(p_value) && !is.null(raw$p_value)) as.numeric(raw$p_value) else NA_real_
  lwr <- if (isTRUE(ci) && !is.null(raw$lwr)) as.numeric(raw$lwr) else NA_real_
  upr <- if (isTRUE(ci) && !is.null(raw$upr)) as.numeric(raw$upr) else NA_real_

  out <- data.frame(
    method = method,
    kappa = as.numeric(raw$estimate),
    observed_agreement = as.numeric(raw$observed_agreement),
    expected_agreement = as.numeric(raw$expected_agreement),
    n_items = as.integer(raw$n_items),
    n_categories = as.integer(raw$n_categories),
    n_raters_min = as.integer(raw$n_raters_min),
    n_raters_max = as.integer(raw$n_raters_max),
    balanced = as.logical(raw$balanced),
    n_ratings_total = as.integer(raw$n_ratings_total),
    se = se,
    statistic = statistic,
    p_value = p_val,
    lwr.ci = lwr,
    upr.ci = upr,
    conf.level = if (isTRUE(ci)) conf_level else NA_real_,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  by_category_df <- NULL
  if (isTRUE(by_category)) {
    by_category_df <- data.frame(
      category = as.character(categories),
      kappa = as.numeric(raw$by_category_estimate),
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

  class(out) <- c("multirater_kappa", "agreement_result", "data.frame")
  attr(out, "method") <- "multirater_kappa"
  attr(out, "description") <- "Multi-rater chance-corrected nominal agreement"
  attr(out, "package") <- "matrixCorr"
  attr(out, "input") <- input
  attr(out, "kappa_method") <- method
  attr(out, "exact") <- isTRUE(raw$exact)
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
  attr(out, "reference") <- switch(
    method,
    fleiss = "Fleiss (1971), Psychological Bulletin, 76, 378-382, doi:10.1037/h0031619",
    randolph = "Randolph (2005), Free-Marginal Multirater Kappa"
  )
  out
}

#' @method print multirater_kappa
#' @rdname multirater_kappa
#' @param x A \code{multirater_kappa} object.
#' @param digits Integer; number of decimal places for displayed values.
#' @param show_by_category Logical; whether to print the attached category-wise
#'   kappa table when available.
#' @param ... Additional arguments passed to [print.data.frame()].
#' @export
print.multirater_kappa <- function(x,
                                   digits = 4,
                                   show_by_category = FALSE,
                                   ...) {
  check_inherits(x, "multirater_kappa")
  check_bool(show_by_category, arg = "show_by_category")

  nr <- attr(x, "n_raters", exact = TRUE)
  method_label <- attr(x, "kappa_method", exact = TRUE)
  if (isTRUE(attr(x, "exact", exact = TRUE))) {
    method_label <- paste0(method_label, " (exact)")
  }
  cat("Multi-rater kappa agreement\n")
  cat("  method             :", method_label, "\n")
  cat("  kappa              :", formatC(x$kappa[[1L]], format = "f", digits = digits), "\n")
  cat("  observed agreement :", formatC(x$observed_agreement[[1L]], format = "f", digits = digits), "\n")
  cat("  expected agreement :", formatC(x$expected_agreement[[1L]], format = "f", digits = digits), "\n")
  cat("  items              :", .mc_count_fmt(x$n_items[[1L]]), "\n")
  cat("  categories         :", .mc_count_fmt(x$n_categories[[1L]]), "\n")
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
      cat("\nBy-category kappas\n\n")
      print.data.frame(by_cat, row.names = FALSE, ...)
    }
  }
  invisible(x)
}

#' @method summary multirater_kappa
#' @rdname multirater_kappa
#' @param object A \code{multirater_kappa} object.
#' @param digits Integer; number of decimal places for displayed values.
#' @param ... Unused.
#' @export
summary.multirater_kappa <- function(object,
                                     digits = 4,
                                     ...) {
  check_inherits(object, "multirater_kappa")

  out <- as.data.frame(object, stringsAsFactors = FALSE, check.names = FALSE)
  preferred_order <- c(
    "method",
    "kappa",
    "observed_agreement",
    "expected_agreement",
    "n_items",
    "n_categories",
    "n_raters_min",
    "n_raters_max",
    "balanced",
    "n_ratings_total",
    "se",
    "statistic",
    "p_value",
    "lwr.ci",
    "upr.ci",
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
  out <- .mc_finalize_summary_df(out, class_name = "summary.multirater_kappa")
  out <- out[, c(intersect(preferred_order, names(out)), setdiff(names(out), preferred_order)), drop = FALSE]
  class(out) <- unique(c("summary.multirater_kappa", "summary.agreement_result", class(out)))

  diag_attr <- attr(object, "diagnostics", exact = TRUE)
  item_agreement <- as.numeric(diag_attr$item_agreement %||% numeric())
  n_raters_by_item <- as.integer(diag_attr$n_raters_by_item %||% integer())
  attr(out, "diagnostics_summary") <- list(
    category_counts = as.integer(diag_attr$category_counts %||% integer()),
    category_proportions = as.numeric(diag_attr$category_proportions %||% numeric()),
    item_agreement_summary = if (length(item_agreement)) {
      c(
        min = min(item_agreement, na.rm = TRUE),
        mean = mean(item_agreement, na.rm = TRUE),
        median = stats::median(item_agreement, na.rm = TRUE),
        max = max(item_agreement, na.rm = TRUE)
      )
    } else {
      c(min = NA_real_, mean = NA_real_, median = NA_real_, max = NA_real_)
    },
    n_raters_summary = if (length(n_raters_by_item)) {
      c(
        min = min(n_raters_by_item, na.rm = TRUE),
        mean = mean(n_raters_by_item, na.rm = TRUE),
        median = stats::median(n_raters_by_item, na.rm = TRUE),
        max = max(n_raters_by_item, na.rm = TRUE)
      )
    } else {
      c(min = NA_real_, mean = NA_real_, median = NA_real_, max = NA_real_)
    }
  )
  attr(out, "by_category") <- attr(object, "by_category", exact = TRUE)
  attr(out, "kappa_method") <- attr(object, "kappa_method", exact = TRUE)
  attr(out, "exact") <- attr(object, "exact", exact = TRUE)
  attr(out, "conf.level") <- attr(object, "conf.level", exact = TRUE)
  attr(out, "se.method") <- attr(object, "se.method", exact = TRUE)
  attr(out, "digits") <- digits
  out
}

#' @method print summary.multirater_kappa
#' @rdname multirater_kappa
#' @param n Optional preview row threshold.
#' @param topn Optional number of leading/trailing rows when truncated.
#' @param max_vars Optional maximum number of visible columns.
#' @param width Optional display width.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @export
print.summary.multirater_kappa <- function(x,
                                           digits = NULL,
                                           n = NULL,
                                           topn = NULL,
                                           max_vars = NULL,
                                           width = NULL,
                                           show_ci = NULL,
                                           ...) {
  digits <- .mc_coalesce(digits, attr(x, "digits", exact = TRUE) %||% 4L)
  digest <- c(
    method = if (isTRUE(attr(x, "exact", exact = TRUE))) {
      paste0(attr(x, "kappa_method", exact = TRUE) %||% NA_character_, " (exact)")
    } else {
      attr(x, "kappa_method", exact = TRUE) %||% NA_character_
    },
    kappa = formatC(x$kappa[[1L]], format = "f", digits = digits),
    items = .mc_count_fmt(x$n_items[[1L]]),
    categories = .mc_count_fmt(x$n_categories[[1L]])
  )
  se_method <- attr(x, "se.method", exact = TRUE) %||% "none"
  if (!identical(se_method, "none")) {
    digest <- c(digest, se_method = se_method)
  }
  digest <- digest[!is.na(digest) & nzchar(digest)]
  .mc_print_named_digest(digest, header = "Multi-rater kappa agreement summary")
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
  invisible(x)
}

#' @method plot multirater_kappa
#' @rdname multirater_kappa
#' @param type Plot type: \code{"agreement_map"}, \code{"estimate"}, \code{"item_agreement"},
#'   \code{"category_proportion"}, or \code{"by_category"}.
#' @param bins Integer number of bins retained for compatibility; currently
#'   unused by the default item-agreement profile plot.
#' @export
plot.multirater_kappa <- function(x,
                                  type = c("agreement_map", "estimate", "item_agreement", "category_proportion", "by_category"),
                                  bins = 30L,
                                  ...) {
  check_inherits(x, "multirater_kappa")
  type <- match.arg(type)
  bins <- check_scalar_int_pos(bins, arg = "bins")
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} is required for plotting.")
  }

  diag_attr <- attr(x, "diagnostics", exact = TRUE)
  if (identical(type, "agreement_map")) {
    counts <- diag_attr$item_category_counts
    if (is.null(counts) || !is.matrix(counts) || !nrow(counts) || !ncol(counts)) {
      abort_bad_arg("type", message = "cannot be {.val agreement_map} because retained item-category counts are unavailable.")
    }
    item_agreement <- as.numeric(diag_attr$item_agreement %||% rep(NA_real_, nrow(counts)))
    n_raters_by_item <- as.numeric(diag_attr$n_raters_by_item %||% rowSums(counts))
    row_id <- seq_len(nrow(counts))
    dominant <- apply(counts, 1L, function(z) {
      idx <- which.max(z)
      if (!length(idx) || all(z == 0)) "" else colnames(counts)[[idx[[1L]]]]
    })
    ord <- order(item_agreement, decreasing = TRUE, na.last = TRUE)
    counts <- counts[ord, , drop = FALSE]
    item_agreement <- item_agreement[ord]
    n_raters_by_item <- n_raters_by_item[ord]
    dominant <- dominant[ord]
    row_id <- row_id[ord]

    df <- as.data.frame(as.table(counts), stringsAsFactors = FALSE)
    names(df) <- c("item", "category", "count")
    df$item_index <- match(df$item, rownames(counts))
    df$n_raters <- n_raters_by_item[df$item_index]
    df$share <- ifelse(df$n_raters > 0, df$count / df$n_raters, NA_real_)
    df$item_agreement <- item_agreement[df$item_index]
    df$dominant <- dominant[df$item_index]
    df$item_label <- sprintf(
      "item %s | a=%s | top=%s",
      row_id[df$item_index],
      formatC(df$item_agreement, format = "f", digits = 2),
      df$dominant
    )
    item_levels <- unique(df$item_label)
    df$item_label <- factor(df$item_label, levels = rev(item_levels))

    return(
      ggplot2::ggplot(df, ggplot2::aes(x = .data$category, y = .data$item_label, fill = .data$share)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.5) +
        ggplot2::geom_text(
          ggplot2::aes(label = ifelse(.data$count > 0, .data$count, "")),
          size = 3,
          color = "black"
        ) +
        ggplot2::scale_fill_gradient(
          low = "#F3E9DC",
          high = "#274C77",
          limits = c(0, 1),
          name = "Share of raters"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          ...
        ) +
        ggplot2::labs(
          title = "Item-by-category agreement map",
          subtitle = "Items ordered by item-level agreement; tile labels are rater counts",
          x = "Category"
        )
    )
  }

  if (identical(type, "estimate")) {
    df <- data.frame(
      label = "Panel kappa",
      x = x$kappa[[1L]],
      lwr = x$lwr.ci[[1L]],
      upr = x$upr.ci[[1L]],
      stringsAsFactors = FALSE
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$label)) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.5, color = "grey80") +
      ggplot2::geom_vline(xintercept = 1, linewidth = 0.4, color = "grey90", linetype = "dashed") +
      ggplot2::geom_vline(xintercept = -1, linewidth = 0.4, color = "grey90", linetype = "dashed")
    if (is.finite(df$lwr[[1L]]) && is.finite(df$upr[[1L]])) {
      p <- p + ggplot2::geom_segment(
        ggplot2::aes(x = .data$lwr, xend = .data$upr, y = .data$label, yend = .data$label),
        linewidth = 1.2,
        color = "#5B7C99"
      )
    }
    return(
      p +
        ggplot2::geom_point(size = 3.2, color = "#1F4E79") +
        ggplot2::coord_cartesian(xlim = c(-1, 1)) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          ...
        ) +
        ggplot2::labs(
          title = "Multi-rater kappa estimate",
          subtitle = sprintf(
            "%s method; %s items across %s categories",
            if (isTRUE(attr(x, "exact", exact = TRUE))) {
              paste0(attr(x, "kappa_method", exact = TRUE), " (exact)")
            } else {
              attr(x, "kappa_method", exact = TRUE)
            },
            .mc_count_fmt(x$n_items[[1L]]),
            .mc_count_fmt(x$n_categories[[1L]])
          ),
          x = "Kappa"
        )
    )
  }

  if (identical(type, "item_agreement")) {
    df <- data.frame(
      item = seq_along(diag_attr$item_agreement %||% numeric()),
      item_agreement = sort(as.numeric(diag_attr$item_agreement %||% numeric())),
      stringsAsFactors = FALSE
    )
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = .data$item, y = .data$item_agreement)) +
        ggplot2::geom_hline(
          yintercept = x$observed_agreement[[1L]],
          linewidth = 0.6,
          linetype = "dashed",
          color = "#B04A5A"
        ) +
        ggplot2::geom_linerange(
          ggplot2::aes(ymin = 0, ymax = .data$item_agreement),
          linewidth = 0.8,
          color = "#7AA6C2"
        ) +
        ggplot2::geom_point(size = 2.3, color = "#1F4E79") +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), ...) +
        ggplot2::labs(
          title = "Item-level agreement profile",
          subtitle = sprintf(
            "Dashed line shows mean observed agreement = %s",
            formatC(x$observed_agreement[[1L]], format = "f", digits = 3)
          ),
          x = "Items ordered by agreement",
          y = "Agreement"
        )
    )
  }

  if (identical(type, "category_proportion")) {
    df <- data.frame(
      category = factor(
        as.character(attr(x, "categories", exact = TRUE)),
        levels = rev(as.character(attr(x, "categories", exact = TRUE)))
      ),
      proportion = as.numeric(diag_attr$category_proportions %||% numeric()),
      stringsAsFactors = FALSE
    )
    return(
      ggplot2::ggplot(df, ggplot2::aes(y = .data$category, x = .data$proportion)) +
        ggplot2::geom_segment(
          ggplot2::aes(x = 0, xend = .data$proportion, yend = .data$category),
          linewidth = 1,
          color = "#C5A46D"
        ) +
        ggplot2::geom_point(size = 3, color = "#8C5A2B") +
        ggplot2::coord_cartesian(xlim = c(0, 1)) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), ...) +
        ggplot2::labs(
          title = "Marginal category proportions",
          x = "Proportion",
          y = "Category"
        )
    )
  }

  by_cat <- attr(x, "by_category", exact = TRUE)
  if (!is.data.frame(by_cat) || !nrow(by_cat)) {
    abort_bad_arg(
      "type",
      message = "cannot be {.val by_category} unless {.arg by_category} was requested in the original fit."
    )
  }
  by_cat$category <- factor(by_cat$category, levels = by_cat$category[order(by_cat$kappa, decreasing = TRUE)])
  p <- ggplot2::ggplot(by_cat, ggplot2::aes(y = .data$category, x = .data$kappa)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.5, color = "grey80") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), ...) +
    ggplot2::labs(
      title = "Category-wise binary-collapsed kappas",
      x = "Kappa",
      y = "Category"
    )
  if (all(c("lwr.ci", "upr.ci") %in% names(by_cat)) &&
      any(is.finite(by_cat$lwr.ci) & is.finite(by_cat$upr.ci))) {
    p <- p +
      ggplot2::geom_segment(
        data = by_cat[is.finite(by_cat$lwr.ci) & is.finite(by_cat$upr.ci), , drop = FALSE],
        ggplot2::aes(x = .data$lwr.ci, xend = .data$upr.ci, yend = .data$category),
        linewidth = 1,
        color = "#7E9F45"
      )
  }
  p +
    ggplot2::geom_point(size = 3, color = "#4C6A2F") +
    ggplot2::coord_cartesian(xlim = c(-1, 1))
}
