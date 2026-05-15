#' @title Krippendorff's Alpha Reliability Coefficient
#'
#' @description
#' Estimates Krippendorff's alpha as a panel-level reliability/agreement
#' coefficient among two or more coders, raters, observers, judges, or
#' instruments. Missing ratings are supported through \code{na_method}, and
#' nominal, ordinal, interval, and ratio disagreement functions are available.
#'
#' @param data Input data.
#'   For \code{input = "ratings"}, rows are units/items and columns are
#'   coders/raters. For \code{input = "counts"}, rows are units/items and
#'   columns are categories containing per-item category counts.
#' @param levels Optional category labels in analysis order. For ordinal
#'   character data, explicit \code{levels} are required unless the data are
#'   ordered factors. For interval and ratio alpha, \code{levels} must be
#'   coercible to finite numeric values.
#' @param input One of \code{"ratings"} or \code{"counts"}.
#' @param level Measurement level: \code{"nominal"}, \code{"ordinal"},
#'   \code{"interval"}, or \code{"ratio"}.
#' @param method Estimator. \code{"customary"} uses the coincidence-matrix
#'   estimator of Krippendorff/Hayes. \code{"analytical"} uses the improved
#'   estimator and jackknife inference direction described by Hughes.
#' @param na_method Missing-data rule. \code{"error"} rejects missing ratings,
#'   \code{"complete"} removes any item with one or more missing ratings, and
#'   \code{"available"} retains items with at least \code{min_raters} observed
#'   ratings. For \code{input = "counts"}, missing counts are not allowed and
#'   rows with total counts below \code{min_raters} are dropped.
#' @param min_raters Minimum number of observed ratings required for a retained
#'   item. Must be at least \code{2}.
#' @param ci Logical; if \code{TRUE}, attach confidence intervals.
#' @param p_value Logical; if \code{TRUE}, attach test statistics and p-values.
#'   This is available only for \code{method = "analytical"} with jackknife
#'   inference.
#' @param conf_level Confidence level for intervals. Default is \code{0.95}.
#' @param se_method Standard-error method. \code{"auto"} resolves to
#'   \code{"bootstrap"} for customary alpha when \code{ci = TRUE}, and to
#'   \code{"jackknife"} for analytical alpha when \code{ci = TRUE} or
#'   \code{p_value = TRUE}. \code{"bootstrap"} is allowed only for customary
#'   alpha. \code{"jackknife"} is the primary inference path for analytical
#'   alpha. \code{"none"} disables uncertainty output.
#' @param n_boot Number of bootstrap resamples used for customary alpha when
#'   \code{ci = TRUE}.
#' @param seed Optional positive integer seed used for reproducible bootstrap
#'   resampling.
#' @param n_threads Integer \eqn{\ge 1}. Number of OpenMP threads used by the
#'   C++ backend.
#' @param return_matrices Logical; if \code{TRUE}, attach the retained count
#'   matrix, observed coincidence matrix, expected coincidence matrix, and
#'   disagreement matrix in \code{attr(x, "matrices")}.
#' @param verbose Logical; if \code{TRUE}, emit short progress messages.
#' @param ... Reserved for future extensions. Unsupported extra arguments are
#'   rejected.
#'
#' @details
#' Krippendorff's alpha is a panel-level reliability coefficient, not a
#' pairwise correlation matrix. The customary coincidence-matrix estimator is
#' \deqn{ \alpha = 1 - \frac{D_o}{D_e}, }
#' where \eqn{D_o} is the observed disagreement and \eqn{D_e} is the expected
#' disagreement under chance assignment.
#'
#' For retained item \eqn{u} with \eqn{m_u} observed ratings and category
#' counts \eqn{n_{uc}}, the observed coincidence matrix updates by
#' \deqn{ O_{ck} \leftarrow O_{ck} +
#'        \frac{n_{uc}(n_{uk} - 1(c=k))}{m_u - 1}. }
#' Let \eqn{n_c = \sum_k O_{ck}} and \eqn{n = \sum_c n_c}. Then
#' \deqn{ D_o = \frac{1}{n}\sum_{c,k} O_{ck}\,\delta^2_{ck}, }
#' \deqn{ D_e = \frac{1}{n(n-1)}\sum_{c,k} n_c n_k\,\delta^2_{ck}, }
#' where \eqn{\delta^2_{ck}} is the level-specific disagreement matrix.
#'
#' The disagreement functions implemented here are:
#' \itemize{
#'   \item nominal: \eqn{\delta^2_{ck} = 1(c \ne k)}
#'   \item ordinal: Krippendorff's cumulative-margin disagreement based on the
#'   pooled category margins
#'   \item interval: \eqn{\delta^2_{ck} = (v_c - v_k)^2}
#'   \item ratio:
#'   \eqn{\delta^2_{ck} = (v_c - v_k)^2 / (v_c + v_k)^2}, with
#'   \eqn{\delta^2_{ck}=0} when both values are zero
#' }
#'
#' \strong{Analytical estimator and inference.}
#' Hughes (2024) recommends a different point estimator based on within-item
#' and pooled pairwise disagreement. Let \eqn{a} be the number of retained
#' items, \eqn{m_i} the number of ratings in item \eqn{i}, and
#' \eqn{N = \sum_i m_i}. Define
#' \deqn{ \bar n = \frac{N - \sum_i m_i^2 / N}{a - 1}. }
#' The within-item mean square is
#' \deqn{ MSE =
#'        \frac{1}{N}\sum_{i=1}^a
#'        \frac{\sum_{j<k}\delta^2(x_{ij}, x_{ik})}{m_i - 1}, }
#' the pooled total disagreement is
#' \deqn{ SST = \frac{1}{N}\sum_{r<s}\delta^2(z_r, z_s), }
#' where \eqn{z_1,\ldots,z_N} are the pooled observed ratings, and
#' \deqn{ SSA = SST - (N-a)\,MSE, \qquad
#'        MSA = \frac{SSA}{a-1}. }
#' With \eqn{\theta = MSA/MSE} and \eqn{\eta = \log(\theta)}, the analytical
#' alpha estimate is
#' \deqn{ \alpha =
#'        \frac{\theta - 1}{\theta + \bar n - 1}. }
#'
#' When analytical inference is requested, the implementation uses the
#' delete-one-item jackknife on \eqn{\eta}. If \eqn{\eta_{(-i)}} is the
#' leave-one-item-out transform and \eqn{a} is the number of retained items,
#' the jackknife pseudo-values are
#' \deqn{ p_i = a\,\hat\eta - (a-1)\,\hat\eta_{(-i)}, }
#' with standard error
#' \deqn{ \widehat{\mathrm{se}}(\hat\eta)
#'        = \sqrt{\mathrm{Var}(p_1,\ldots,p_a) / a}. }
#' The confidence interval is formed on the \eqn{\eta}-scale with a
#' \eqn{t_{a-1}} critical value and back-transformed to alpha using the full
#' data \eqn{\bar n}. The analytical p-value tests \eqn{H_0: \alpha = 0},
#' equivalently \eqn{H_0: \eta = 0}, via the t statistic
#' \deqn{ t = \hat\eta / \widehat{\mathrm{se}}(\hat\eta). }
#'
#' For customary alpha, percentile bootstrap intervals are available by
#' resampling retained items with replacement and recomputing the customary
#' estimate.
#'
#' @return
#' A one-row data frame with class
#' \code{c("krippendorff_alpha", "agreement_result", "data.frame")}. The
#' object stores method metadata, diagnostics, and optional matrices in
#' attributes.
#'
#' @references
#' Krippendorff, K. (2011/2013). Computing Krippendorff's alpha-reliability.
#'
#' Hayes, A. F. and Krippendorff, K. (2007). Answering the call for a standard
#' reliability measure for coding data. \emph{Communication Methods and
#' Measures}, 1(1), 77-89.
#'
#' Hughes, J. (2024). Toward improved inference for Krippendorff's Alpha
#' agreement coefficient. \emph{Journal of Statistical Planning and
#' Inference}, 233, 106170.
#'
#' @seealso
#' [multirater_kappa()] for nominal multi-rater kappa;
#' [gwet_ac()] for AC1/AC2 agreement coefficients.
#'
#' @examples
#' raters <- data.frame(
#'   r1 = c("A", "A", "B", "B", "C", "A"),
#'   r2 = c("A", "B", "B", "B", "C", "A"),
#'   r3 = c("A", "A", "B", "C", "C", "A"),
#'   stringsAsFactors = FALSE
#' )
#'
#' fit <- krippendorff_alpha(raters, level = "nominal", na_method = "available")
#' print(fit)
#' summary(fit)
#' plot(fit)
#'
#' @author Thiago de Paula Oliveira
#' @export
krippendorff_alpha <- function(data,
                               levels = NULL,
                               input = c("ratings", "counts"),
                               level = c("nominal", "ordinal", "interval", "ratio"),
                               method = c("customary", "analytical"),
                               na_method = c("error", "complete", "available"),
                               min_raters = 2L,
                               ci = FALSE,
                               p_value = FALSE,
                               conf_level = 0.95,
                               se_method = c("auto", "bootstrap", "jackknife", "none"),
                               n_boot = 1000L,
                               seed = NULL,
                               n_threads = getOption("matrixCorr.threads", 1L),
                               return_matrices = FALSE,
                               verbose = FALSE,
                               ...) {
  .mc_extract_legacy_aliases(list(...), allowed = character())

  input <- match_arg(input, c("ratings", "counts"), arg_name = "input")
  level <- match_arg(level, c("nominal", "ordinal", "interval", "ratio"), arg_name = "level")
  method <- match_arg(method, c("customary", "analytical"), arg_name = "method")
  na_method <- validate_na_method(
    na_method,
    arg = "na_method",
    allowed = c("error", "complete", "available")
  )
  min_raters <- check_scalar_int_pos(min_raters, arg = "min_raters")
  if (min_raters < 2L) {
    abort_bad_arg("min_raters", message = "must be >= 2.")
  }
  check_bool(ci, arg = "ci")
  check_bool(p_value, arg = "p_value")
  check_bool(return_matrices, arg = "return_matrices")
  check_bool(verbose, arg = "verbose")
  check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")
  if (!is.null(seed)) {
    seed <- check_scalar_int_pos(seed, arg = "seed")
  }

  se_method <- match_arg(
    se_method,
    c("auto", "bootstrap", "jackknife", "none"),
    arg_name = "se_method"
  )
  resolved_se_method <- .mc_krippendorff_resolve_se_method(
    method = method,
    se_method = se_method,
    ci = ci,
    p_value = p_value
  )

  if (identical(resolved_se_method, "bootstrap")) {
    n_boot <- check_scalar_int_pos(n_boot, arg = "n_boot")
  }

  prep <- if (identical(input, "ratings")) {
    .mc_prepare_agreement_ratings(
      data = data,
      levels = levels,
      level = level,
      na_method = na_method,
      min_raters = min_raters
    )
  } else {
    .mc_prepare_agreement_counts(
      data = data,
      levels = levels,
      level = level,
      min_raters = min_raters
    )
  }

  delta2 <- .mc_krippendorff_delta2(
    categories = prep$categories,
    level = level,
    coincidence_margins = colSums(prep$counts)
  )

  prev_threads <- .mc_prepare_omp_threads(
    n_threads,
    n_threads_missing = missing(n_threads)
  )
  if (!is.null(prev_threads)) {
    on.exit(.mc_exit_omp_threads(prev_threads), add = TRUE)
  }

  raw <- krippendorff_alpha_core_cpp(
    counts = prep$counts,
    delta2 = delta2,
    method_code = if (identical(method, "analytical")) 2L else 1L,
    return_matrices = return_matrices,
    n_threads = n_threads
  )

  estimate <- .mc_clamp_unit(as.numeric(raw$estimate))
  se_out <- NA_real_
  statistic_out <- NA_real_
  p_value_out <- NA_real_
  lwr <- NA_real_
  upr <- NA_real_
  ci_method <- NULL
  boot_vals <- NULL
  boot_kept <- NULL
  jackknife <- NULL

  if (identical(resolved_se_method, "bootstrap")) {
    boot_vals <- .mc_krippendorff_bootstrap_customary(
      counts = prep$counts,
      delta2 = delta2,
      n_boot = n_boot,
      seed = seed,
      n_threads = n_threads
    )
    boot_kept <- sum(is.finite(boot_vals))
    if (boot_kept >= 2L) {
      se_out <- stats::sd(boot_vals[is.finite(boot_vals)])
    }
    if (isTRUE(ci)) {
      ci_vals <- .mc_percentile_boot_ci(boot_vals, conf_level = conf_level)
      lwr <- .mc_clamp_unit(ci_vals[[1L]])
      upr <- .mc_clamp_unit(ci_vals[[2L]])
      ci_method <- "percentile_bootstrap"
    }
  } else if (identical(resolved_se_method, "jackknife")) {
    jackknife <- .mc_krippendorff_analytical_jackknife(
      counts = prep$counts,
      delta2 = delta2,
      full_raw = raw,
      conf_level = conf_level,
      ci = ci,
      p_value = p_value,
      n_threads = n_threads
    )
    se_out <- jackknife$se
    statistic_out <- jackknife$statistic
    p_value_out <- jackknife$p_value
    lwr <- jackknife$lwr
    upr <- jackknife$upr
    ci_method <- if (isTRUE(ci)) "jackknife_t_eta" else NULL
  }

  out <- data.frame(
    method = "krippendorff_alpha",
    alpha = estimate,
    level = level,
    estimator = method,
    observed_disagreement = as.numeric(raw$observed_disagreement),
    expected_disagreement = as.numeric(raw$expected_disagreement),
    n_items = as.integer(raw$n_items),
    n_categories = as.integer(raw$n_categories),
    n_raters_min = as.integer(raw$n_raters_min),
    n_raters_max = as.integer(raw$n_raters_max),
    balanced = as.logical(raw$balanced),
    n_ratings_total = as.integer(raw$n_ratings_total),
    se = se_out,
    statistic = statistic_out,
    p_value = p_value_out,
    lwr.ci = if (isTRUE(ci)) lwr else NA_real_,
    upr.ci = if (isTRUE(ci)) upr else NA_real_,
    conf.level = if (isTRUE(ci)) conf_level else NA_real_,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  class(out) <- c("krippendorff_alpha", "agreement_result", "data.frame")

  matrices_attr <- NULL
  if (isTRUE(return_matrices)) {
    coincidence <- as.matrix(raw$coincidence)
    expected_coincidence <- as.matrix(raw$expected_coincidence)
    colnames(coincidence) <- rownames(coincidence) <- as.character(prep$categories)
    colnames(expected_coincidence) <- rownames(expected_coincidence) <- as.character(prep$categories)
    colnames(delta2) <- rownames(delta2) <- as.character(prep$categories)
    counts_mat <- prep$counts
    colnames(counts_mat) <- as.character(prep$categories)
    matrices_attr <- list(
      counts = counts_mat,
      coincidence = coincidence,
      expected_coincidence = expected_coincidence,
      delta2 = delta2
    )
  }

  diagnostics <- c(
    list(
      category_counts = as.integer(raw$category_counts),
      category_proportions = as.numeric(raw$category_proportions),
      item_category_counts = prep$counts,
      n_raters_by_item = as.integer(raw$n_raters_by_item),
      item_disagreement = as.numeric(raw$item_disagreement),
      customary_alpha = as.numeric(raw$customary_alpha),
      analytical_alpha = as.numeric(raw$analytical_alpha),
      mse = as.numeric(raw$mse),
      sst = as.numeric(raw$sst),
      ssa = as.numeric(raw$ssa),
      msa = as.numeric(raw$msa),
      theta = as.numeric(raw$theta),
      eta = as.numeric(raw$eta),
      n_bar = as.numeric(raw$n_bar),
      ci_method = ci_method,
      boot_finite = boot_kept,
      eta_jackknife = if (!is.null(jackknife)) jackknife$eta_minus else NULL,
      eta_pseudovalues = if (!is.null(jackknife)) jackknife$pseudo else NULL
    ),
    prep$diagnostics
  )

  attr(out, "method") <- "krippendorff_alpha"
  attr(out, "description") <- "Krippendorff's alpha reliability coefficient"
  attr(out, "package") <- "matrixCorr"
  attr(out, "input") <- input
  attr(out, "level") <- level
  attr(out, "estimator") <- method
  attr(out, "categories") <- prep$categories
  attr(out, "n_categories") <- length(prep$categories)
  attr(out, "n_items") <- as.integer(raw$n_items)
  attr(out, "n_raters") <- list(
    min = as.integer(raw$n_raters_min),
    max = as.integer(raw$n_raters_max),
    balanced = as.logical(raw$balanced)
  )
  attr(out, "na_method") <- if (identical(input, "counts")) na_method else prep$na_method
  attr(out, "min_raters") <- min_raters
  attr(out, "se.method") <- resolved_se_method
  attr(out, "conf.level") <- if (isTRUE(ci)) conf_level else NULL
  attr(out, "n_boot") <- if (identical(resolved_se_method, "bootstrap")) n_boot else NULL
  attr(out, "diagnostics") <- diagnostics
  attr(out, "reference") <- if (identical(method, "analytical")) {
    paste(
      "Hughes (2024), Journal of Statistical Planning and Inference, 233, 106170;",
      "Krippendorff (2011/2013), Computing Krippendorff's alpha-reliability"
    )
  } else {
    "Krippendorff (2011/2013), Computing Krippendorff's alpha-reliability"
  }
  if (isTRUE(return_matrices)) {
    attr(out, "matrices") <- matrices_attr
  }

  inform_if_verbose(
    "Computed Krippendorff's alpha ({.val {method}}, {.val {level}}) for {attr(out, 'n_items', exact = TRUE)} item{?s} across {length(prep$categories)} categor{?y/ies}.",
    .verbose = verbose
  )
  out
}

.mc_agreement_rating_column <- function(x,
                                        arg = as.character(substitute(x))) {
  out <- .mc_weighted_kappa_prepare_vector(x, arg = arg)
  out$is_factor <- is.factor(x)
  out$factor_levels <- if (is.factor(x)) levels(x) else NULL
  out
}

.mc_resolve_agreement_levels <- function(prepared,
                                         levels = NULL,
                                         level = c("nominal", "ordinal", "interval", "ratio")) {
  level <- match.arg(level)
  observed <- unlist(lapply(prepared, function(z) z$key[!is.na(z$key)]), use.names = FALSE)

  if (!is.null(levels)) {
    lev <- .mc_weighted_kappa_normalize_levels(levels)
    lev_key <- as.character(lev)
    if (length(observed) && any(!observed %in% lev_key)) {
      abort_bad_arg("levels", message = "must contain every non-missing observed category.")
    }
    if (level %in% c("interval", "ratio")) {
      num <- suppressWarnings(as.numeric(lev_key))
      if (anyNA(num) || any(!is.finite(num))) {
        abort_bad_arg("levels", message = "must be coercible to finite numeric values for interval/ratio alpha.")
      }
      if (identical(level, "ratio") && any(num < 0)) {
        abort_bad_arg("levels", message = "must be non-negative for ratio alpha.")
      }
    }
    return(lev)
  }

  if (!length(observed)) {
    abort_bad_arg(
      "levels",
      message = "must be supplied when no non-missing ratings are available to infer categories."
    )
  }

  if (identical(level, "nominal")) {
    if (all(vapply(prepared, function(z) isTRUE(z$numeric_infer), logical(1)))) {
      vals <- unlist(lapply(prepared, function(z) z$numeric_values[!is.na(z$numeric_values)]), use.names = FALSE)
      return(sort(unique(vals)))
    }
    from_levels <- unlist(lapply(prepared, function(z) z$factor_levels), use.names = FALSE)
    if (length(from_levels)) {
      out <- unique(c(from_levels, observed))
      return(out)
    }
    return(unique(observed))
  }

  if (identical(level, "ordinal")) {
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
    numeric_ok <- suppressWarnings(as.numeric(observed))
    if (all(is.finite(numeric_ok))) {
      return(sort(unique(numeric_ok)))
    }
    abort_bad_arg(
      "levels",
      message = paste(
        "must be supplied for ordinal character/unordered-factor data.",
        "Automatic order inference is limited to identical ordered-factor levels or numeric values."
      )
    )
  }

  numeric_ok <- suppressWarnings(as.numeric(observed))
  if (anyNA(numeric_ok) || any(!is.finite(numeric_ok))) {
    abort_bad_arg(
      "levels",
      message = "must be supplied as numeric/coercible values for interval/ratio alpha."
    )
  }
  out <- sort(unique(numeric_ok))
  if (identical(level, "ratio") && any(out < 0)) {
    abort_bad_arg("data", message = "must not contain negative categories for ratio alpha.")
  }
  out
}

.mc_prepare_agreement_ratings <- function(data,
                                          levels = NULL,
                                          level = c("nominal", "ordinal", "interval", "ratio"),
                                          na_method = c("error", "complete", "available"),
                                          min_raters = 2L) {
  level <- match.arg(level)
  na_method <- match.arg(na_method)
  if (!is.data.frame(data) && !is.matrix(data)) {
    abort_bad_arg("data", message = "must be a matrix or data frame when {.arg input} is {.val ratings}.")
  }

  df <- if (is.data.frame(data)) data else as.data.frame(data, stringsAsFactors = FALSE)
  if (ncol(df) < 2L) {
    abort_bad_arg("data", message = "must contain at least two rater columns.")
  }

  cols <- vector("list", ncol(df))
  for (j in seq_along(cols)) {
    col_arg <- names(df)[[j]]
    if (is.null(col_arg) || !nzchar(col_arg)) {
      col_arg <- sprintf("data[[%d]]", j)
    } else {
      col_arg <- sprintf("data$%s", col_arg)
    }
    cols[[j]] <- .mc_agreement_rating_column(df[[j]], arg = col_arg)
  }
  categories <- .mc_resolve_agreement_levels(cols, levels = levels, level = level)
  category_keys <- as.character(categories)
  if (length(categories) < 2L) {
    abort_bad_arg("levels", message = "must define at least two categories.")
  }

  X <- matrix(
    NA_integer_,
    nrow = nrow(df),
    ncol = ncol(df),
    dimnames = list(NULL, names(df))
  )
  for (j in seq_along(cols)) {
    codes <- match(cols[[j]]$key, category_keys)
    codes[is.na(cols[[j]]$key)] <- NA_integer_
    bad <- !is.na(cols[[j]]$key) & is.na(codes)
    if (any(bad)) {
      abort_bad_arg("data", message = "contains values not present in {.arg levels}.")
    }
    X[, j] <- as.integer(codes)
  }

  n_available <- rowSums(!is.na(X))
  keep <- rep(TRUE, nrow(X))
  if (identical(na_method, "error") && anyNA(X)) {
    abort_bad_arg(
      "data",
      message = "contains missing ratings when {.arg na_method} is {.val error}."
    )
  }
  if (identical(na_method, "complete")) {
    keep <- stats::complete.cases(X)
  } else if (identical(na_method, "available")) {
    keep <- n_available >= min_raters
  }
  X_retained <- X[keep, , drop = FALSE]
  retained_available <- n_available[keep]
  if (length(retained_available)) {
    keep_raters <- retained_available >= min_raters
    X_retained <- X_retained[keep_raters, , drop = FALSE]
    idx_keep <- which(keep)
    keep[idx_keep[!keep_raters]] <- FALSE
  }

  if (nrow(X_retained) < 2L) {
    abort_bad_arg(
      "data",
      message = "must retain at least two items with at least {.val {min_raters}} observed ratings."
    )
  }

  counts <- matrix(
    0L,
    nrow = nrow(X_retained),
    ncol = length(categories),
    dimnames = list(paste0("item_", seq_len(nrow(X_retained))), as.character(categories))
  )
  for (i in seq_len(nrow(X_retained))) {
    xi <- X_retained[i, ]
    xi <- xi[!is.na(xi)]
    if (length(xi)) {
      counts[i, ] <- tabulate(xi, nbins = ncol(counts))
    }
  }

  list(
    counts = counts,
    categories = categories,
    n_levels = length(categories),
    retained_rows = keep,
    dropped_rows = !keep,
    na_method = na_method,
    diagnostics = list(
      n_original = nrow(df),
      n_retained = nrow(counts),
      retained_rows = keep,
      dropped_rows = !keep,
      n_available_by_item = n_available
    )
  )
}

.mc_prepare_agreement_counts <- function(data,
                                         levels = NULL,
                                         level = c("nominal", "ordinal", "interval", "ratio"),
                                         min_raters = 2L) {
  level <- match.arg(level)
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

  x_names <- colnames(X)
  if (!is.null(x_names) && anyDuplicated(x_names)) {
    abort_bad_arg("data", message = "must not contain duplicate category column names.")
  }

  categories <- NULL
  counts_aligned <- NULL
  if (!is.null(levels)) {
    categories <- .mc_resolve_agreement_levels(
      prepared = list(list(
        key = if (!is.null(x_names)) x_names else character(),
        is_ordered = FALSE,
        ordered_levels = NULL,
        numeric_values = NULL,
        numeric_infer = FALSE,
        is_factor = FALSE,
        factor_levels = NULL
      )),
      levels = levels,
      level = level
    )
    cat_keys <- as.character(categories)
    if (is.null(x_names)) {
      if (length(categories) != ncol(X)) {
        abort_bad_arg(
          "levels",
          message = "must have length equal to ncol(data) when count columns are unnamed."
        )
      }
      counts_aligned <- matrix(
        as.integer(round(X)),
        nrow = nrow(X),
        ncol = ncol(X),
        dimnames = list(NULL, cat_keys)
      )
    } else {
      idx <- match(x_names, cat_keys)
      if (anyNA(idx)) {
        abort_bad_arg("levels", message = "must contain every observed count column name.")
      }
      if (anyDuplicated(idx)) {
        abort_bad_arg("data", message = "maps multiple columns to the same category in {.arg levels}.")
      }
      counts_aligned <- matrix(
        0L,
        nrow = nrow(X),
        ncol = length(categories),
        dimnames = list(NULL, cat_keys)
      )
      counts_aligned[, idx] <- as.integer(round(X))
    }
  } else {
    if (is.null(x_names)) {
      if (level %in% c("interval", "ratio")) {
        abort_bad_arg(
          "levels",
          message = "must be supplied or inferred from count-column names for interval/ratio alpha."
        )
      }
      categories <- seq_len(ncol(X))
    } else {
      categories <- if (level %in% c("interval", "ratio")) {
        num <- suppressWarnings(as.numeric(x_names))
        if (anyNA(num) || any(!is.finite(num))) {
          abort_bad_arg(
            "levels",
            message = "count-column names must be coercible to finite numeric values for interval/ratio alpha."
          )
        }
        if (identical(level, "ratio") && any(num < 0)) {
          abort_bad_arg("data", message = "must not contain negative categories for ratio alpha.")
        }
        num
      } else if (identical(level, "ordinal")) {
        x_names
      } else {
        x_names
      }
    }
    categories <- .mc_resolve_agreement_levels(
      prepared = list(list(
        key = as.character(categories),
        is_ordered = FALSE,
        ordered_levels = NULL,
        numeric_values = if (is.numeric(categories)) categories else NULL,
        numeric_infer = is.numeric(categories),
        is_factor = FALSE,
        factor_levels = NULL
      )),
      levels = categories,
      level = level
    )
    counts_aligned <- matrix(
      as.integer(round(X)),
      nrow = nrow(X),
      ncol = ncol(X),
      dimnames = list(NULL, as.character(categories))
    )
  }

  row_sums <- rowSums(counts_aligned)
  keep <- row_sums >= min_raters
  counts_retained <- counts_aligned[keep, , drop = FALSE]
  rownames(counts_retained) <- paste0("item_", seq_len(nrow(counts_retained)))
  if (nrow(counts_retained) < 2L) {
    abort_bad_arg(
      "data",
      message = "must retain at least two items with at least {.val {min_raters}} ratings."
    )
  }

  list(
    counts = counts_retained,
    categories = categories,
    n_levels = length(categories),
    retained_rows = keep,
    dropped_rows = !keep,
    na_method = "available",
    diagnostics = list(
      n_original = nrow(counts_aligned),
      n_retained = nrow(counts_retained),
      retained_rows = keep,
      dropped_rows = !keep,
      row_sums = row_sums
    )
  )
}

.mc_krippendorff_delta2 <- function(categories,
                                    level = c("nominal", "ordinal", "interval", "ratio"),
                                    coincidence_margins = NULL) {
  level <- match.arg(level)
  K <- length(categories)
  out <- matrix(0, nrow = K, ncol = K, dimnames = list(as.character(categories), as.character(categories)))
  if (K <= 1L) {
    return(out)
  }

  if (identical(level, "nominal")) {
    out[row(out) != col(out)] <- 1
    return(out)
  }

  if (identical(level, "ordinal")) {
    if (is.null(coincidence_margins) || length(coincidence_margins) != K) {
      abort_bad_arg(
        "coincidence_margins",
        message = "must be supplied with one pooled category margin per category for ordinal alpha."
      )
    }
    ng <- as.numeric(coincidence_margins)
    if (any(!is.finite(ng)) || any(ng < 0)) {
      abort_bad_arg("coincidence_margins", message = "must contain finite non-negative margins.")
    }
    cs <- cumsum(ng)
    for (c in seq_len(K)) {
      if (c >= K) {
        next
      }
      for (k in seq.int(c + 1L, K)) {
        mass_ck <- cs[[k]] - if (c > 1L) cs[[c - 1L]] else 0
        val <- (mass_ck - (ng[[c]] + ng[[k]]) / 2)^2
        out[c, k] <- val
        out[k, c] <- val
      }
    }
    diag(out) <- 0
    return(out)
  }

  vals <- suppressWarnings(as.numeric(as.character(categories)))
  if (anyNA(vals) || any(!is.finite(vals))) {
    abort_bad_arg("levels", message = "must be coercible to finite numeric values for interval/ratio alpha.")
  }
  if (identical(level, "ratio") && any(vals < 0)) {
    abort_bad_arg("levels", message = "must be non-negative for ratio alpha.")
  }

  for (c in seq_len(K)) {
    for (k in seq_len(K)) {
      if (c == k) {
        next
      }
      if (identical(level, "interval")) {
        out[c, k] <- (vals[[c]] - vals[[k]])^2
      } else {
        denom <- (vals[[c]] + vals[[k]])^2
        out[c, k] <- if (denom <= 0) 0 else ((vals[[c]] - vals[[k]])^2) / denom
      }
    }
  }
  out
}

.mc_clamp_unit <- function(x, tol = 1e-12) {
  if (!is.finite(x)) {
    return(NA_real_)
  }
  if (x > 1 && x <= 1 + tol) {
    return(1)
  }
  if (x < -1 && x >= -1 - tol) {
    return(-1)
  }
  x
}

.mc_krippendorff_resolve_se_method <- function(method,
                                               se_method,
                                               ci = FALSE,
                                               p_value = FALSE) {
  inference_requested <- isTRUE(ci) || isTRUE(p_value)
  if (!inference_requested) {
    return("none")
  }

  if (identical(method, "customary")) {
    if (isTRUE(p_value)) {
      abort_bad_arg(
        "p_value",
        message = "must be FALSE for the customary estimator.",
        .hint = "P-values are available only for {.arg method} = {.val analytical} with jackknife inference."
      )
    }
    if (identical(se_method, "auto")) {
      return("bootstrap")
    }
    if (!identical(se_method, "bootstrap")) {
      abort_bad_arg(
        "se_method",
        message = "must be {.val bootstrap}, {.val auto}, or {.val none} for customary alpha."
      )
    }
    return("bootstrap")
  }

  if (identical(se_method, "auto")) {
    return("jackknife")
  }
  if (identical(se_method, "bootstrap")) {
    abort_bad_arg(
      "se_method",
      message = "{.val bootstrap} is not available for the analytical estimator.",
      .hint = "Use {.arg se_method} = {.val jackknife}."
    )
  }
  if (identical(se_method, "none")) {
    abort_bad_arg(
      "se_method",
      message = "must not be {.val none} when ci or p_value is TRUE."
    )
  }
  "jackknife"
}

.mc_krippendorff_bootstrap_customary <- function(counts,
                                                 delta2,
                                                 n_boot,
                                                 seed = NULL,
                                                 n_threads = 1L) {
  n_items <- nrow(counts)
  .mc_eval_with_seed(seed, {
    replicate(n_boot, {
      idx <- sample.int(n_items, size = n_items, replace = TRUE)
      raw_b <- krippendorff_alpha_core_cpp(
        counts = counts[idx, , drop = FALSE],
        delta2 = delta2,
        method_code = 1L,
        return_matrices = FALSE,
        n_threads = n_threads
      )
      .mc_clamp_unit(as.numeric(raw_b$estimate))
    })
  })
}

.mc_krippendorff_alpha_from_eta <- function(eta, n_bar) {
  theta <- exp(eta)
  .mc_clamp_unit((theta - 1) / (theta + n_bar - 1))
}

.mc_krippendorff_analytical_jackknife <- function(counts,
                                                  delta2,
                                                  full_raw,
                                                  conf_level = 0.95,
                                                  ci = FALSE,
                                                  p_value = FALSE,
                                                  n_threads = 1L) {
  out <- list(
    se = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    lwr = NA_real_,
    upr = NA_real_,
    eta_minus = rep(NA_real_, nrow(counts)),
    pseudo = rep(NA_real_, nrow(counts))
  )

  a <- as.integer(full_raw$n_items)
  eta_hat <- as.numeric(full_raw$eta)
  n_bar <- as.numeric(full_raw$n_bar)
  if (a < 3L || !is.finite(eta_hat) || !is.finite(n_bar)) {
    return(out)
  }

  eta_minus <- vapply(seq_len(a), function(i) {
    raw_i <- krippendorff_alpha_core_cpp(
      counts = counts[-i, , drop = FALSE],
      delta2 = delta2,
      method_code = 2L,
      return_matrices = FALSE,
      n_threads = n_threads
    )
    as.numeric(raw_i$eta)
  }, numeric(1))

  if (any(!is.finite(eta_minus))) {
    out$eta_minus <- eta_minus
    return(out)
  }

  pseudo <- a * eta_hat - (a - 1) * eta_minus
  se_eta <- sqrt(stats::var(pseudo) / a)
  out$eta_minus <- eta_minus
  out$pseudo <- pseudo
  out$se <- if (is.finite(se_eta) && se_eta > 0) se_eta else NA_real_
  if (!is.finite(out$se) || out$se <= 0) {
    return(out)
  }

  df <- a - 1
  tcrit <- stats::qt(1 - (1 - conf_level) / 2, df = df)
  if (isTRUE(ci)) {
    lo_eta <- eta_hat - tcrit * out$se
    hi_eta <- eta_hat + tcrit * out$se
    out$lwr <- .mc_krippendorff_alpha_from_eta(lo_eta, n_bar)
    out$upr <- .mc_krippendorff_alpha_from_eta(hi_eta, n_bar)
  }
  if (isTRUE(p_value)) {
    out$statistic <- eta_hat / out$se
    out$p_value <- 2 * stats::pt(abs(out$statistic), df = df, lower.tail = FALSE)
  }
  out
}

#' @method print krippendorff_alpha
#' @rdname krippendorff_alpha
#' @param x A \code{krippendorff_alpha} object.
#' @param digits Integer; number of decimal places for displayed numeric values.
#' @export
print.krippendorff_alpha <- function(x,
                                     digits = 4,
                                     ...) {
  check_inherits(x, "krippendorff_alpha")
  nr <- attr(x, "n_raters", exact = TRUE)
  has_ci <- is.finite(x$lwr.ci[[1L]]) && is.finite(x$upr.ci[[1L]])
  digest <- c(
    estimator = attr(x, "estimator", exact = TRUE),
    level = attr(x, "level", exact = TRUE),
    alpha = formatC(x$alpha[[1L]], format = "f", digits = digits),
    items = .mc_count_fmt(x$n_items[[1L]]),
    categories = .mc_count_fmt(x$n_categories[[1L]]),
    raters_per_item = paste0(
      .mc_count_fmt(nr$min), " to ", .mc_count_fmt(nr$max),
      if (isTRUE(nr$balanced)) " (balanced)" else " (unbalanced)"
    ),
    ci = if (has_ci) "yes" else "no",
    se_method = .mc_coalesce(attr(x, "se.method", exact = TRUE), "none")
  )
  .mc_print_named_digest(digest, header = "Krippendorff's alpha reliability coefficient")
  if (has_ci) {
    cat(
      sprintf(
        "  %-12s: %g%% [%s, %s]\n",
        "interval",
        100 * suppressWarnings(as.numeric(attr(x, "conf.level", exact = TRUE))),
        formatC(x$lwr.ci[[1L]], format = "f", digits = digits),
        formatC(x$upr.ci[[1L]], format = "f", digits = digits)
      )
    )
  }
  if (is.finite(x$p_value[[1L]])) {
    cat(sprintf("  %-12s: %s\n", "p_value", formatC(x$p_value[[1L]], format = "f", digits = digits)))
  }
  invisible(x)
}

#' @method summary krippendorff_alpha
#' @rdname krippendorff_alpha
#' @param object A \code{krippendorff_alpha} object.
#' @param digits Integer; number of decimal places for rounded numeric columns.
#' @param ... Unused.
#' @export
summary.krippendorff_alpha <- function(object,
                                       digits = 4,
                                       ...) {
  check_inherits(object, "krippendorff_alpha")
  out <- as.data.frame(object, stringsAsFactors = FALSE, check.names = FALSE)
  preferred_order <- c(
    "method", "alpha", "lwr.ci", "upr.ci", "se", "statistic", "p_value",
    "observed_disagreement", "expected_disagreement", "n_items",
    "n_categories", "n_ratings_total", "n_raters_min", "n_raters_max",
    "balanced", "level", "estimator"
  )
  out$se_method <- .mc_coalesce(attr(object, "se.method", exact = TRUE), "none")
  out <- out[, c(intersect(preferred_order, names(out)), "se_method", setdiff(names(out), c(preferred_order, "se_method"))), drop = FALSE]
  drop_if_all_na <- c("lwr.ci", "upr.ci", "se", "statistic", "p_value", "conf.level")
  keep_cols <- vapply(names(out), function(nm) {
    if (!nm %in% drop_if_all_na) {
      return(TRUE)
    }
    !all(is.na(out[[nm]]))
  }, logical(1))
  out <- out[, keep_cols, drop = FALSE]
  num_cols <- vapply(out, is.numeric, logical(1))
  out[num_cols] <- lapply(out[num_cols], round, digits = digits)
  out <- .mc_finalize_summary_df(out, class_name = "summary.krippendorff_alpha")
  class(out) <- unique(c("summary.krippendorff_alpha", "summary.agreement_result", class(out)))
  attr(out, "digits") <- digits
  attr(out, "conf.level") <- attr(object, "conf.level", exact = TRUE)
  attr(out, "se.method") <- attr(object, "se.method", exact = TRUE)
  attr(out, "estimator") <- attr(object, "estimator", exact = TRUE)
  attr(out, "level") <- attr(object, "level", exact = TRUE)
  out
}

#' @method print summary.krippendorff_alpha
#' @rdname krippendorff_alpha
#' @param n Optional preview row threshold.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns.
#' @param width Optional display width.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @export
print.summary.krippendorff_alpha <- function(x,
                                             digits = NULL,
                                             n = NULL,
                                             topn = NULL,
                                             max_vars = NULL,
                                             width = NULL,
                                             show_ci = NULL,
                                             ...) {
  digits <- .mc_coalesce(digits, .mc_coalesce(attr(x, "digits", exact = TRUE), 4L))
  show_ci <- .mc_resolve_show_ci(show_ci, context = "summary")
  digest <- c(
    estimator = attr(x, "estimator", exact = TRUE),
    level = attr(x, "level", exact = TRUE),
    alpha = formatC(x$alpha[[1L]], format = "f", digits = digits),
    items = .mc_count_fmt(x$n_items[[1L]]),
    categories = .mc_count_fmt(x$n_categories[[1L]])
  )
  se_method <- .mc_coalesce(attr(x, "se.method", exact = TRUE), "none")
  if (!identical(se_method, "none")) {
    digest <- c(digest, se_method = se_method)
  }
  .mc_print_named_digest(
    digest,
    header = .mc_header_with_ci(
      "Krippendorff's alpha summary",
      suppressWarnings(as.numeric(.mc_coalesce(attr(x, "conf.level", exact = TRUE), NA_real_))),
      show_ci = show_ci
    )
  )
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

#' @method plot krippendorff_alpha
#' @rdname krippendorff_alpha
#' @param type Plot type: \code{"estimate"}, \code{"item_disagreement"},
#'   \code{"category_proportion"}, or \code{"coincidence"}.
#' @export
plot.krippendorff_alpha <- function(x,
                                    type = c("estimate", "item_disagreement", "category_proportion", "coincidence"),
                                    ...) {
  check_inherits(x, "krippendorff_alpha")
  type <- match.arg(type)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} is required for plotting.")
  }

  diag_attr <- attr(x, "diagnostics", exact = TRUE)
  if (identical(type, "estimate")) {
    df <- data.frame(
      label = "Krippendorff's alpha",
      x = x$alpha[[1L]],
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
          title = "Krippendorff's alpha estimate",
          subtitle = sprintf(
            "%s estimator; %s level",
            attr(x, "estimator", exact = TRUE),
            attr(x, "level", exact = TRUE)
          ),
          x = "Alpha"
        )
    )
  }

  if (identical(type, "item_disagreement")) {
        item_dis <- as.numeric(.mc_coalesce(diag_attr$item_disagreement, numeric()))
    if (!length(item_dis) || !any(is.finite(item_dis))) {
      abort_bad_arg("type", message = "cannot be {.val item_disagreement} because item-level disagreement diagnostics are unavailable.")
    }
    df <- data.frame(
      item = seq_along(item_dis),
      item_disagreement = sort(item_dis),
      stringsAsFactors = FALSE
    )
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = .data$item, y = .data$item_disagreement)) +
        ggplot2::geom_hline(
          yintercept = mean(item_dis, na.rm = TRUE),
          linewidth = 0.6,
          linetype = "dashed",
          color = "#B04A5A"
        ) +
        ggplot2::geom_linerange(
          ggplot2::aes(ymin = 0, ymax = .data$item_disagreement),
          linewidth = 0.8,
          color = "#7AA6C2"
        ) +
        ggplot2::geom_point(size = 2.3, color = "#1F4E79") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), ...) +
        ggplot2::labs(
          title = "Item-level disagreement profile",
          subtitle = "Dashed line shows the mean item-level disagreement",
          x = "Items ordered by disagreement",
          y = "Disagreement"
        )
    )
  }

  if (identical(type, "category_proportion")) {
    props <- as.numeric(.mc_coalesce(diag_attr$category_proportions, numeric()))
    cats <- as.character(attr(x, "categories", exact = TRUE))
    df <- data.frame(
      category = factor(cats, levels = rev(cats)),
      proportion = props,
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

  mats <- attr(x, "matrices", exact = TRUE)
  if (is.null(mats) || !is.matrix(mats$coincidence)) {
    abort_bad_arg(
      "type",
      message = "cannot be {.val coincidence} unless {.arg return_matrices} was requested in the original fit."
    )
  }
  df <- as.data.frame(as.table(mats$coincidence), stringsAsFactors = FALSE)
  names(df) <- c("category1", "category2", "value")
  ggplot2::ggplot(df, ggplot2::aes(x = .data$category2, y = .data$category1, fill = .data$value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(label = formatC(.data$value, format = "f", digits = 1)),
      size = 3
    ) +
    ggplot2::scale_fill_gradient(low = "#F4EBD0", high = "#2E5EAA", name = "Coincidence") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      ...
    ) +
    ggplot2::labs(title = "Observed coincidence matrix")
}
