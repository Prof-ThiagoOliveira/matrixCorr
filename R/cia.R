#' @title Coefficient of Individual Agreement
#'
#' @description
#' `cia()` estimates the coefficient of individual agreement. CIA assesses
#' individual-level interchangeability by comparing between-method disagreement
#' with within-method replicate disagreement. Unlike CCC, CIA is not intended
#' to be driven by between-subject heterogeneity.
#'
#' The estimator requires replicated readings within method. A data set with
#' one observation per subject per method is insufficient for CIA because the
#' within-method disagreement term cannot be estimated.
#'
#' @param data A data frame in long format.
#' @param response Character scalar naming the numeric measurement column.
#' @param subject Character scalar naming the subject/unit identifier column.
#' @param method Character scalar naming the method/device/rater column.
#' @param replicate Character scalar naming the replicate identifier within each
#'   subject-method cell.
#' @param reference Optional character scalar naming the reference method. When
#'   `NULL` (default), the no-reference CIA is estimated.
#' @param scope One of `"pairwise"` or `"overall"`. `"pairwise"` returns the
#'   two-method
#'   CIA for each method pair, or the reference-scaled CIA for each
#'   method-versus-reference comparison when `reference` is supplied.
#'   `"overall"` returns one CIA coefficient across all retained methods.
#' @param estimator One of `"mom_unconstrained"` or `"vc_constrained"`.
#'   `"mom_unconstrained"` is the literal raw method-of-moments CIA ratio and
#'   may exceed 1. `"vc_constrained"` is a bounded package-level
#'   variance-component variant that constrains the estimated inter-method
#'   variance component to be non-negative before converting it back to CIA.
#' @param ci Logical; if `TRUE`, attach confidence intervals using the selected
#'   `inference` method.
#' @param conf_level Confidence level used when `ci = TRUE`. Default is `0.95`.
#' @param inference One of `"delta"`, `"bootstrap"`, or `"none"`. When
#'   confidence intervals are requested, pairwise CIA uses a delta-method normal
#'   interval by default and also supports subject-bootstrap percentile
#'   intervals. For overall CIA, delta-method normal intervals are available for
#'   the unconstrained moment estimator, and subject-bootstrap percentile
#'   intervals are also available.
#' @param B Number of subject-bootstrap resamples when `ci = TRUE` and
#'   `inference = "bootstrap"`.
#' @param seed Optional positive integer seed for reproducible bootstrap
#'   resampling.
#' @param n_threads Integer >= 1. Number of OpenMP threads passed to the C++
#'   backend.
#' @param verbose Logical; if `TRUE`, emit short progress messages.
#'
#' @details
#' Let `Y_ijk` denote replicate `k` on subject `i` measured by method `j`.
#'
#' Without a reference method, CIA is defined by
#' \deqn{
#' \mathrm{CIA}^N =
#' \frac{\sum_j E\{(Y_{ijk} - Y_{ijk'})^2\} / 2}
#'      {\sum_{j < j'} E\{(Y_{ij} - Y_{ij'})^2\} / (J - 1)}.
#' }
#'
#' With a reference method `J`, CIA is defined by
#' \deqn{
#' \mathrm{CIA}^R =
#' \frac{E\{(Y_{iJk} - Y_{iJk'})^2\}}
#'      {\sum_{j=1}^{J-1} E\{(Y_{ij} - Y_{iJ})^2\} / (J - 1)}.
#' }
#'
#' This implementation uses method-of-moments estimators built from
#' subject-level disagreement functions.
#'
#' For pairwise CIA, consider two methods `X` and `Y`. For each eligible
#' subject `i`, define
#' \deqn{
#' G_{1i} = \mathrm{mean}\{(X_{ik} - X_{ik'})^2\},
#' }
#' \deqn{
#' G_{2i} = \mathrm{mean}\{(Y_{il} - Y_{il'})^2\},
#' }
#' \deqn{
#' G_{3i} = \mathrm{mean}\{(X_{ik} - Y_{il})^2\},
#' }
#' where the first two means are over all distinct within-method replicate
#' pairs and the third mean is over all cross-method replicate combinations.
#' Let \eqn{\bar G_1}, \eqn{\bar G_2}, and \eqn{\bar G_3} denote the averages
#' of these subject-level quantities across eligible subjects. Then the
#' no-reference pairwise estimator is
#' \deqn{
#' \widehat{\mathrm{CIA}}_{XY}^{\mathrm{raw}} =
#' \frac{(\bar G_1 + \bar G_2)/2}{\bar G_3},
#' }
#' and the reference pairwise estimator is
#' \deqn{
#' \widehat{\mathrm{CIA}}_{XR}^{\mathrm{raw}} = \frac{\bar G_1}{\bar G_3},
#' }
#' where `X` is the reference method.
#'
#' For overall CIA, the current implementation follows the balanced replicated
#' formulas in the cited papers. This requires `n` common subjects, `J`
#' retained methods, and the same replicate count `K >= 2` in every retained
#' subject-method cell. If the data are not balanced in this sense,
#' `scope = "overall"` returns an informative error rather than using an
#' approximation.
#'
#' Write \eqn{\bar Y_{ij}} for the mean of the `K` replicates on subject `i`
#' and method `j`. Define the within-method mean square
#' \deqn{
#' A_{ij} = \frac{\sum_{k=1}^{K} (Y_{ijk} - \bar Y_{ij})^2}{K - 1}.
#' }
#'
#' For the no-reference overall estimator, define the subject-level within term
#' \deqn{
#' A_{i\cdot} = \frac{1}{J} \sum_{j=1}^{J} A_{ij},
#' }
#' the subject-wide mean
#' \deqn{
#' \bar Y_{i\cdot\cdot} = \frac{1}{J} \sum_{j=1}^{J} \bar Y_{ij},
#' }
#' and the subject-level denominator
#' \deqn{
#' B_{Ni} =
#' \frac{\sum_{j=1}^{J} (\bar Y_{ij} - \bar Y_{i\cdot\cdot})^2}{J - 1}
#' + \left(1 - \frac{1}{K}\right) A_{i\cdot}.
#' }
#' With \eqn{\bar A = n^{-1}\sum_i A_{i\cdot}} and
#' \eqn{\bar B = n^{-1}\sum_i B_{Ni}}, the overall no-reference estimator is
#' \deqn{
#' \widehat{\mathrm{CIA}}_{N}^{\mathrm{raw}} = \frac{\bar A}{\bar B}.
#' }
#'
#' For the reference overall estimator, let method `R` be the retained
#' reference method. Define
#' \deqn{
#' A_{iR} = A_{i,R}
#' }
#' and
#' \deqn{
#' B_{Ri} =
#' \frac{\sum_{j \ne R} (\bar Y_{ij} - \bar Y_{iR})^2}{J - 1}
#' + \left(1 - \frac{1}{K}\right)
#'   \frac{\sum_{j \ne R} A_{ij}}{J - 1}
#' + \left(1 - \frac{1}{K}\right) A_{iR}.
#' }
#' With \eqn{\bar A_R = n^{-1}\sum_i A_{iR}} and
#' \eqn{\bar B_R = n^{-1}\sum_i B_{Ri}}, the overall reference estimator is
#' \deqn{
#' \widehat{\mathrm{CIA}}_{R}^{\mathrm{raw}} = \frac{2 \bar A_R}{\bar B_R}.
#' }
#'
#' Two estimators are available. `estimator = "mom_unconstrained"` reports the
#' raw ratio estimator directly. `estimator = "vc_constrained"` applies the
#' non-negative variance-component boundary from the cited papers. In the
#' no-reference setting this uses
#' \eqn{\hat\tau^2_{\mathrm{raw}} = \bar B - \bar A}, and in the reference
#' setting it uses \eqn{\hat\tau^2_{\mathrm{raw}} = \bar B_R/2 - \bar A_R}.
#' The constrained estimator then sets
#' \eqn{\hat\tau^2 = \max(\hat\tau^2_{\mathrm{raw}}, 0)} and reports
#' \eqn{\widehat{\mathrm{CIA}} = \widehat{\sigma}_W^2 /
#' (\widehat{\sigma}_W^2 + \hat\tau^2)} on the corresponding scale. This
#' applies the boundary on the implied inter-method variance component rather
#' than clamping CIA directly.
#'
#' When confidence intervals are requested, pairwise CIA uses a large-sample
#' delta-method normal interval by default and also supports subject-bootstrap
#' percentile intervals. For overall CIA, delta-method normal intervals are
#' available for the unconstrained moment estimator, and subject-bootstrap
#' percentile intervals are also available.
#'
#' High CIA indicates stronger individual agreement. The FDA individual
#' bioequivalence boundary `IEC <= 2.4948` corresponds to `CIA >= 0.445`, and
#' `CIA >= 0.8` is sometimes used as a stronger practical rule. Such thresholds
#' are context-dependent and are not hard-coded by this function.
#'
#' Missing rows in the required columns are removed before estimation and the
#' counts are recorded in `attr(x, "diagnostics")`.
#'
#' @return
#' For `scope = "overall"`, a one-row data frame with class
#' `c("cia_overall", "cia", "data.frame")`. For `scope = "pairwise"` and
#' `ci = FALSE`, a dense matrix-style object using the package's standard
#' correlation-result infrastructure. For `scope = "pairwise"` and `ci = TRUE`,
#' a list with elements `est`, `lwr.ci`, and `upr.ci`, classed as
#' `c("cia", "cia_ci")`.
#'
#' @references
#' Barnhart HX, Kosinski AS, Haber M. (2007). Assessing individual agreement.
#' \emph{Journal of Biopharmaceutical Statistics}, 17(4), 697-719.
#' \doi{10.1080/10543400701329489}
#'
#' Barnhart HX, Haber M, Lokhnygina Y, Kosinski AS. (2007). Comparison of
#' concordance correlation coefficient and coefficient of individual agreement
#' in assessing agreement. \emph{Journal of Biopharmaceutical Statistics},
#' 17(4), 721-738.
#'
#' Pan Y, Gao J, Haber M, Barnhart HX. (2010). Estimation of coefficients of
#' individual agreement (CIA's) for quantitative and binary data using SAS and
#' R. \emph{Computer Methods and Programs in Biomedicine}.
#'
#' @seealso \code{\link{ccc}}, \code{\link{ba}}, and \code{\link{prob_agree}}.
#' \code{cia()} answers the question "Are methods interchangeable for
#' individual subjects once within-method replicate disagreement is taken into
#' account?". In contrast, \code{ccc()} answers "How well do two paired
#' measurements agree overall, combining correlation with location/scale
#' agreement?".
#'
#' @examples
#' set.seed(1)
#' subjects <- sprintf("s%02d", 1:30)
#' methods <- c("A", "B", "C")
#' replicates <- sprintf("r%02d", 1:20)
#' dat <- expand.grid(
#'   subject = subjects,
#'   method = methods,
#'   replicate = replicates,
#'   KEEP.OUT.ATTRS = FALSE,
#'   stringsAsFactors = FALSE
#' )
#' subject_effect <- stats::rnorm(length(subjects), sd = 3)
#' method_shift <- c(A = 0, B = 0.15, C = -0.10)
#' method_sd <- c(A = 0.45, B = 0.45, C = 0.30)
#' dat$value <- subject_effect[match(dat$subject, subjects)] +
#'   method_shift[dat$method] +
#'   stats::rnorm(nrow(dat), sd = method_sd[dat$method])
#'
#' fit_overall <- cia(
#'   dat,
#'   response = "value",
#'   subject = "subject",
#'   method = "method",
#'   replicate = "replicate",
#'   scope = "overall"
#' )
#' print(fit_overall)
#' summary(fit_overall)
#' plot(fit_overall)
#'
#' fit_pairwise <- cia(
#'   dat,
#'   response = "value",
#'   subject = "subject",
#'   method = "method",
#'   replicate = "replicate",
#'   scope = "pairwise"
#' )
#' print(fit_pairwise)
#' summary(fit_pairwise)
#' plot(fit_pairwise)
#'
#' @author Thiago de Paula Oliveira
#' @export
cia <- function(data,
                response,
                subject,
                method,
                replicate,
                reference = NULL,
                scope = c("pairwise", "overall"),
                estimator = c("mom_unconstrained", "vc_constrained"),
                ci = FALSE,
                conf_level = 0.95,
                inference = c("delta", "bootstrap", "none"),
                B = 1000L,
                seed = NULL,
                n_threads = getOption("matrixCorr.threads", 1L),
                verbose = FALSE) {
  scope <- match_arg(scope, c("pairwise", "overall"), arg_name = "scope")
  estimator <- match_arg(
    estimator,
    c("mom_unconstrained", "vc_constrained"),
    arg_name = "estimator"
  )
  inference <- if (!missing(inference)) {
    match_arg(inference, c("delta", "bootstrap", "none"), arg_name = "inference")
  } else {
    NULL
  }
  check_bool(ci, arg = "ci")
  check_bool(verbose, arg = "verbose")
  check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")
  B <- check_scalar_int_pos(B, arg = "B")
  if (B < 2L) {
    abort_bad_arg("B", message = "must be >= 2.")
  }
  if (!is.null(seed)) {
    seed <- check_scalar_int_pos(seed, arg = "seed")
  }
  inference <- .mc_cia_resolve_inference(
    ci = ci,
    scope = scope,
    estimator = estimator,
    inference = inference,
    inference_missing = missing(inference)
  )

  prep <- .mc_cia_prepare_long(
    data = data,
    response = response,
    subject = subject,
    method = method,
    replicate = replicate,
    reference = reference,
    verbose = verbose
  )

  raw <- cia_moments_cpp(
    y = prep$y,
    subject = prep$subject_code,
    method = prep$method_code,
    replicate = prep$replicate_code,
    n_methods = prep$n_methods,
    reference_method = prep$reference_index,
    has_reference = prep$has_reference,
    pairwise = identical(scope, "pairwise"),
    n_threads = n_threads
  )

  diagnostics <- .mc_cia_build_diagnostics(prep, raw)

  if (identical(scope, "overall")) {
    overall <- .mc_cia_overall_moments(
      prep,
      estimator = estimator,
      n_threads = n_threads
    )
    out <- .mc_cia_build_overall_from_moments(
      prep = prep,
      diagnostics = diagnostics,
      moments = overall,
      estimator = estimator
    )
    if (isTRUE(ci) && identical(inference, "delta")) {
      if (!identical(estimator, "mom_unconstrained")) {
        abort_bad_arg(
          "inference",
          message = paste(
            "Delta-method confidence intervals are implemented for the unconstrained moment estimator only;",
            "use {.val bootstrap} for the constrained estimator."
          )
        )
      }
      ci_res <- .mc_cia_overall_delta_ci(overall, conf_level = conf_level)
      out$lwr.ci <- ci_res$lwr
      out$upr.ci <- ci_res$upr
      out$se <- ci_res$se
      if (!is.null(ci_res$diagnostic)) {
        attr(out, "ci.diagnostic") <- ci_res$diagnostic
      }
      attr(out, "ci.method") <- "delta_normal"
    } else if (isTRUE(ci) && identical(inference, "bootstrap")) {
      boot <- .mc_cia_bootstrap_overall(
        prep = prep,
        estimator = estimator,
        B = B,
        conf_level = conf_level,
        seed = seed,
        n_threads = n_threads
      )
      out$lwr.ci <- boot$lwr
      out$upr.ci <- boot$upr
      attr(out, "bootstrap") <- boot[c("n_successful", "B")]
      attr(out, "ci.method") <- "subject_bootstrap_percentile"
    }
    attr(out, "conf.level") <- if (isTRUE(ci)) conf_level else NULL
    return(out)
  }

  .mc_cia_assert_estimable(
    prep,
    diagnostics,
    scope = scope,
    n_threads = n_threads
  )

  if (!isTRUE(ci)) {
    return(.mc_cia_build_pairwise_matrix(
      prep,
      diagnostics,
      estimator = estimator,
      n_threads = n_threads
    ))
  }

  boot <- NULL
  if (identical(inference, "bootstrap")) {
    boot <- .mc_cia_bootstrap_pairwise(
      prep = prep,
      estimator = estimator,
      B = B,
      conf_level = conf_level,
      seed = seed,
      n_threads = n_threads
    )
  }
  .mc_cia_build_pairwise_ci(
    prep,
    diagnostics,
    estimator = estimator,
    inference = inference,
    boot = boot,
    conf_level = conf_level,
    n_threads = n_threads
  )
}

.mc_cia_prepare_long <- function(data,
                                 response,
                                 subject,
                                 method,
                                 replicate,
                                 reference,
                                 verbose = FALSE) {
  if (!is.data.frame(data)) {
    abort_bad_arg("data", message = "must be a data frame.")
  }
  cols <- c(response = response, subject = subject, method = method, replicate = replicate)
  for (nm in names(cols)) {
    val <- cols[[nm]]
    if (!is.character(val) || length(val) != 1L || is.na(val) || !nzchar(val)) {
      abort_bad_arg(nm, message = "must be a single column name.")
    }
  }
  check_required_cols(data, unname(cols), df_arg = "data")

  dat <- data[, unname(cols), drop = FALSE]
  names(dat) <- names(cols)
  keep <- stats::complete.cases(dat)
  n_removed <- sum(!keep)
  if (n_removed > 0L) {
    inform_if_verbose(
      "Removed {n_removed} incomplete row{?s} before CIA estimation.",
      .verbose = verbose,
      n_removed = n_removed
    )
    dat <- dat[keep, , drop = FALSE]
  }
  if (!nrow(dat)) {
    abort_bad_arg("data", message = "does not contain any complete rows in the required CIA columns.")
  }

  y <- suppressWarnings(as.numeric(dat$response))
  if (any(!is.finite(y))) {
    abort_bad_arg("response", message = "must name a numeric column with finite values.")
  }

  subj_chr <- as.character(dat$subject)
  meth_chr <- as.character(dat$method)
  rep_chr <- as.character(dat$replicate)

  dup_key <- paste(subj_chr, meth_chr, rep_chr, sep = "\r")
  if (anyDuplicated(dup_key)) {
    abort_bad_arg(
      "data",
      message = "contains duplicated subject-method-replicate combinations.",
      .hint = "Each replicate identifier must be unique within each subject-method cell."
    )
  }

  method_levels <- unique(meth_chr)
  if (length(method_levels) < 2L) {
    abort_bad_arg("method", message = "must retain at least two methods after removing incomplete rows.")
  }
  subject_levels <- unique(subj_chr)
  if (length(subject_levels) < 1L) {
    abort_bad_arg("subject", message = "must retain at least one subject after removing incomplete rows.")
  }

  if (is.null(reference)) {
    reference_index <- 0L
    reference_name <- NA_character_
    has_reference <- FALSE
    coefficient <- "CIA_N"
  } else {
    if (!is.character(reference) || length(reference) != 1L || is.na(reference) || !nzchar(reference)) {
      abort_bad_arg("reference", message = "must be NULL or a single observed method name.")
    }
    reference_index <- match(reference, method_levels, nomatch = 0L)
    if (reference_index <= 0L) {
      abort_bad_arg(
        "reference",
        message = "must match one of the retained method labels.",
        .hint = "Check the spelling of the reference method after missing-row removal."
      )
    }
    reference_name <- reference
    has_reference <- TRUE
    coefficient <- "CIA_R"
  }

  subject_code <- as.integer(factor(subj_chr, levels = subject_levels))
  method_code <- as.integer(factor(meth_chr, levels = method_levels))
  replicate_levels <- unique(rep_chr)
  replicate_code <- as.integer(factor(rep_chr, levels = replicate_levels))
  subject_rows <- split(seq_len(nrow(dat)), subject_code)
  by_subject_method <- lapply(subject_rows, function(idx) split(y[idx], method_code[idx]))

  list(
    data = dat,
    y = y,
    subject_code = subject_code,
    method_code = method_code,
    replicate_code = replicate_code,
    subject_rows = subject_rows,
    by_subject_method = by_subject_method,
    method_levels = method_levels,
    subject_levels = subject_levels,
    replicate_levels = replicate_levels,
    n_methods = length(method_levels),
    n_subjects = length(subject_levels),
    has_reference = has_reference,
    reference_index = as.integer(reference_index),
    reference = reference_name,
    coefficient = coefficient,
    na_rows_removed = as.integer(n_removed),
    n_original_rows = as.integer(nrow(data)),
    n_complete_rows = as.integer(nrow(dat))
  )
}

.mc_cia_named_vector <- function(x, nm) {
  x <- as.numeric(x)
  names(x) <- nm
  x
}

.mc_cia_named_matrix <- function(x, nm) {
  x <- as.matrix(x)
  dimnames(x) <- .mc_square_dimnames(nm)
  x
}

.mc_cia_named_rect_matrix <- function(x, row_nm, col_nm) {
  x <- as.matrix(x)
  dimnames(x) <- list(row_nm, col_nm)
  x
}

.mc_cia_build_diagnostics <- function(prep, raw) {
  within_msd <- .mc_cia_named_vector(raw$within_msd, prep$method_levels)
  sigma_within <- .mc_cia_named_vector(raw$sigma_within, prep$method_levels)
  n_replicate_pairs <- setNames(as.integer(raw$n_replicate_pairs), prep$method_levels)
  between_msd <- .mc_cia_named_matrix(raw$between_msd, prep$method_levels)
  n_between_pairs <- matrix(
    as.integer(raw$n_between_pairs),
    nrow = prep$n_methods,
    ncol = prep$n_methods,
    dimnames = .mc_square_dimnames(prep$method_levels)
  )
  list(
    n_original_rows = prep$n_original_rows,
    n_complete_rows = prep$n_complete_rows,
    na_rows_removed = prep$na_rows_removed,
    n_subjects = as.integer(raw$n_subjects),
    n_obs = as.integer(raw$n_obs),
    method_levels = prep$method_levels,
    reference = prep$reference,
    coefficient = prep$coefficient,
    within_msd = within_msd,
    sigma_within = sigma_within,
    between_msd = between_msd,
    n_replicate_pairs = n_replicate_pairs,
    n_between_pairs = n_between_pairs,
    n_complete = n_between_pairs
  )
}

.mc_cia_resolve_inference <- function(ci, scope, estimator, inference, inference_missing = FALSE) {
  if (!isTRUE(ci)) {
    return("none")
  }
  if (isTRUE(inference_missing) || is.null(inference)) {
    if (identical(scope, "pairwise")) {
      return("delta")
    }
    return(if (identical(estimator, "mom_unconstrained")) "delta" else "bootstrap")
  }
  if (identical(inference, "none")) {
    abort_bad_arg("inference", message = "must not be {.val none} when {.arg ci} is TRUE.")
  }
  if (identical(scope, "overall") &&
      identical(inference, "delta") &&
      !identical(estimator, "mom_unconstrained")) {
    abort_bad_arg(
      "inference",
      message = paste(
        "Delta-method confidence intervals are implemented for the unconstrained moment estimator only;",
        "use {.val bootstrap} for the constrained estimator."
      )
    )
  }
  inference
}

.mc_cia_subject_within_msd <- function(x) {
  n <- length(x)
  if (n < 2L) {
    return(NA_real_)
  }
  d2 <- outer(x, x, FUN = "-")^2
  mean(d2[upper.tri(d2)])
}

.mc_cia_subject_between_msd <- function(x, y) {
  if (!length(x) || !length(y)) {
    return(NA_real_)
  }
  mean(outer(x, y, FUN = "-")^2)
}

.mc_cia_subject_values <- function(prep, subject_idx, method_idx) {
  vals <- prep$by_subject_method[[subject_idx]][[as.character(method_idx)]]
  if (is.null(vals)) numeric(0) else as.numeric(vals)
}

.mc_cia_overall_arrays <- function(prep, n_threads = 1L) {
  raw <- tryCatch(
    cia_overall_balanced_cpp(
      y = prep$y,
      subject = prep$subject_code,
      method = prep$method_code,
      replicate = prep$replicate_code,
      n_methods = prep$n_methods,
      reference_method = prep$reference_index,
      has_reference = prep$has_reference,
      n_threads = n_threads
    ),
    error = function(e) {
      abort_bad_arg("data", message = conditionMessage(e))
    }
  )

  Ybar <- .mc_cia_named_rect_matrix(
    raw$Ybar,
    row_nm = as.character(prep$subject_levels),
    col_nm = prep$method_levels
  )
  A <- .mc_cia_named_rect_matrix(
    raw$A,
    row_nm = as.character(prep$subject_levels),
    col_nm = prep$method_levels
  )

  list(
    Ybar = Ybar,
    A = A,
    n = as.integer(raw$n_subjects),
    J = as.integer(raw$n_methods),
    K = as.integer(raw$K),
    ref_idx = as.integer(raw$ref_idx),
    A_i = as.numeric(raw$A_i),
    B_i = as.numeric(raw$B_i),
    A_bar = as.numeric(raw$A_bar),
    B_bar = as.numeric(raw$B_bar),
    estimate_raw = as.numeric(raw$estimate_raw)
  )
}

.mc_cia_overall_payload <- function(A_bar, B_bar, reference = FALSE, estimator) {
  if (isTRUE(reference)) {
    numerator_term <- 2 * A_bar
    cia_raw <- suppressWarnings(numerator_term / B_bar)
    if (isTRUE(is.nan(cia_raw))) {
      cia_raw <- NA_real_
    }
    tau2_raw <- B_bar / 2 - A_bar
    boundary <- is.finite(tau2_raw) && tau2_raw < 0
    tau2 <- if (is.finite(tau2_raw)) max(tau2_raw, 0) else NA_real_
    cia <- if (identical(estimator, "mom_unconstrained")) {
      cia_raw
    } else {
      denom <- A_bar + tau2
      if (!is.finite(denom) || denom == 0) NA_real_ else A_bar / denom
    }
  } else {
    numerator_term <- A_bar
    cia_raw <- suppressWarnings(A_bar / B_bar)
    if (isTRUE(is.nan(cia_raw))) {
      cia_raw <- NA_real_
    }
    tau2_raw <- B_bar - A_bar
    boundary <- is.finite(tau2_raw) && tau2_raw < 0
    tau2 <- if (is.finite(tau2_raw)) max(tau2_raw, 0) else NA_real_
    cia <- if (identical(estimator, "mom_unconstrained")) {
      cia_raw
    } else {
      denom <- A_bar + tau2
      if (!is.finite(denom) || denom == 0) NA_real_ else A_bar / denom
    }
  }

  list(
    cia = as.numeric(cia),
    cia_raw = as.numeric(cia_raw),
    tau2 = as.numeric(tau2),
    tau2_raw = as.numeric(tau2_raw),
    boundary = isTRUE(boundary),
    numerator_term = as.numeric(numerator_term),
    denominator_term = as.numeric(B_bar),
    A_bar = as.numeric(A_bar)
  )
}

.mc_cia_overall_moments <- function(prep, estimator, n_threads = 1L) {
  arr <- .mc_cia_overall_arrays(prep, n_threads = n_threads)
  payload <- .mc_cia_overall_payload(
    A_bar = arr$A_bar,
    B_bar = arr$B_bar,
    reference = isTRUE(prep$has_reference),
    estimator = estimator
  )

  c(
    list(
      Ybar = arr$Ybar,
      A = arr$A,
      A_i = arr$A_i,
      B_i = arr$B_i,
      n = arr$n,
      J = arr$J,
      K = arr$K,
      ref_idx = arr$ref_idx
    ),
    payload
  )
}

.mc_cia_overall_delta_ci <- function(moments, conf_level = 0.95) {
  A_i <- as.numeric(moments$A_i)
  B_i <- as.numeric(moments$B_i)
  n <- as.integer(moments$n)
  if (n < 2L) {
    return(list(lwr = NA_real_, upr = NA_real_, se = NA_real_, var = NA_real_, diagnostic = "At least two subjects are required for a delta-method CI."))
  }

  A_bar <- mean(A_i)
  B_bar <- mean(B_i)
  if (!is.finite(B_bar) || B_bar <= 0) {
    return(list(lwr = NA_real_, upr = NA_real_, se = NA_real_, var = NA_real_, diagnostic = "The overall denominator term is non-finite or non-positive, so the delta-method CI is not available."))
  }

  var_A <- stats::var(A_i) / n
  var_B <- stats::var(B_i) / n
  cov_A_B <- stats::cov(A_i, B_i) / n
  var_ratio <-
    (1 / (B_bar^2)) * var_A +
    ((A_bar^2) / (B_bar^4)) * var_B -
    ((2 * A_bar) / (B_bar^3)) * cov_A_B
  ratio_multiplier <- if (moments$ref_idx > 0L) 2 else 1
  var_psi <- (ratio_multiplier^2) * var_ratio
  if (!is.finite(var_psi)) {
    return(list(lwr = NA_real_, upr = NA_real_, se = NA_real_, var = NA_real_, diagnostic = "The delta-method variance was non-finite."))
  }
  tol <- sqrt(.Machine$double.eps) * max(1, abs(var_psi))
  if (var_psi < 0 && abs(var_psi) <= tol) {
    var_psi <- 0
  }
  if (var_psi < 0) {
    return(list(lwr = NA_real_, upr = NA_real_, se = NA_real_, var = var_psi, diagnostic = "The delta-method variance was negative beyond floating-point tolerance."))
  }

  psi <- as.numeric(moments$cia_raw)
  se <- sqrt(var_psi)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  list(
    lwr = psi - z * se,
    upr = psi + z * se,
    se = se,
    var = var_psi,
    diagnostic = NULL
  )
}

.mc_cia_component_payload <- function(sigma2_within, between_term, estimator) {
  cia_raw <- suppressWarnings(sigma2_within / between_term)
  if (isTRUE(is.nan(cia_raw))) {
    cia_raw <- NA_real_
  }
  tau2_raw <- between_term - sigma2_within
  boundary <- is.finite(tau2_raw) && tau2_raw < 0
  tau2 <- if (is.finite(tau2_raw)) max(tau2_raw, 0) else NA_real_
  cia <- if (identical(estimator, "mom_unconstrained")) {
    cia_raw
  } else {
    denom <- sigma2_within + tau2
    if (!is.finite(denom) || denom == 0) NA_real_ else sigma2_within / denom
  }
  list(
    cia = as.numeric(cia),
    cia_raw = as.numeric(cia_raw),
    tau2_raw = as.numeric(tau2_raw),
    tau2 = as.numeric(tau2),
    boundary = isTRUE(boundary),
    sigma2_within = as.numeric(sigma2_within),
    between_term = as.numeric(between_term)
  )
}

.mc_cia_pairwise_delta_ci <- function(A,
                                      B,
                                      N,
                                      S11,
                                      S33,
                                      S13,
                                      S22 = NA_real_,
                                      S12 = NA_real_,
                                      S23 = NA_real_,
                                      reference_pair = FALSE,
                                      conf_level = 0.95,
                                      center) {
  if (N < 2L || !is.finite(A) || !is.finite(B) || B == 0) {
    return(list(lwr = NA_real_, upr = NA_real_, se = NA_real_, var = NA_real_))
  }

  if (isTRUE(reference_pair)) {
    var_A <- S11 / N
    var_B <- S33 / N
    cov_A_B <- S13 / N
  } else {
    var_A <- (S11 + S22 + 2 * S12) / (4 * N)
    var_B <- S33 / N
    cov_A_B <- (S13 + S23) / (2 * N)
  }

  var_cia <-
    (1 / (B^2)) * var_A +
    ((A^2) / (B^4)) * var_B -
    ((2 * A) / (B^3)) * cov_A_B
  se <- sqrt(max(var_cia, 0))
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  list(
    lwr = center - z * se,
    upr = center + z * se,
    se = se,
    var = var_cia
  )
}

.mc_cia_pairwise_payload <- function(prep,
                                     diagnostics = NULL,
                                     estimator,
                                     conf_level = NULL,
                                     return_delta = FALSE,
                                     n_threads = 1L) {
  raw <- cia_pairwise_stats_cpp(
    y = prep$y,
    subject = prep$subject_code,
    method = prep$method_code,
    replicate = prep$replicate_code,
    n_methods = prep$n_methods,
    reference_method = prep$reference_index,
    has_reference = prep$has_reference,
    n_threads = n_threads
  )
  nm <- prep$method_levels
  est <- matrix(NA_real_, prep$n_methods, prep$n_methods, dimnames = .mc_square_dimnames(nm))
  cia_raw <- est
  tau2_raw <- est
  tau2 <- est
  sigma2_within <- est
  lwr <- est
  upr <- est
  n_eligible <- matrix(0L, prep$n_methods, prep$n_methods, dimnames = .mc_square_dimnames(nm))
  boundary <- matrix(FALSE, prep$n_methods, prep$n_methods, dimnames = .mc_square_dimnames(nm))

  diag(est) <- 1
  diag(cia_raw) <- 1
  diag(tau2_raw) <- 0
  diag(tau2) <- 0
  diag(sigma2_within) <- 0
  diag(lwr) <- 1
  diag(upr) <- 1
  diag(n_eligible) <- prep$n_subjects

  if (isTRUE(prep$has_reference)) {
    ref_idx <- prep$reference_index
    for (j in seq_len(prep$n_methods)) {
      if (j == ref_idx) {
        next
      }
      N <- as.integer(raw$n_eligible[j, ref_idx])
      n_eligible[j, ref_idx] <- N
      n_eligible[ref_idx, j] <- N
      if (N <= 0L) {
        next
      }
      A <- as.numeric(raw$A_bar[j, ref_idx])
      B <- as.numeric(raw$B_bar[j, ref_idx])
      payload <- .mc_cia_component_payload(
        sigma2_within = A,
        between_term = B,
        estimator = estimator
      )
      est[j, ref_idx] <- payload$cia
      est[ref_idx, j] <- payload$cia
      cia_raw[j, ref_idx] <- payload$cia_raw
      cia_raw[ref_idx, j] <- payload$cia_raw
      tau2_raw[j, ref_idx] <- payload$tau2_raw
      tau2_raw[ref_idx, j] <- payload$tau2_raw
      tau2[j, ref_idx] <- payload$tau2
      tau2[ref_idx, j] <- payload$tau2
      sigma2_within[j, ref_idx] <- payload$sigma2_within
      sigma2_within[ref_idx, j] <- payload$sigma2_within
      boundary[j, ref_idx] <- payload$boundary
      boundary[ref_idx, j] <- payload$boundary
      if (isTRUE(return_delta)) {
        ci <- .mc_cia_pairwise_delta_ci(
          A = A,
          B = B,
          N = N,
          S11 = as.numeric(raw$S11[j, ref_idx]),
          S33 = as.numeric(raw$S33[j, ref_idx]),
          S13 = as.numeric(raw$S13[j, ref_idx]),
          reference_pair = TRUE,
          conf_level = conf_level,
          center = payload$cia
        )
        lwr[j, ref_idx] <- ci$lwr
        lwr[ref_idx, j] <- ci$lwr
        upr[j, ref_idx] <- ci$upr
        upr[ref_idx, j] <- ci$upr
      }
    }
    return(list(
      est = est,
      cia_raw = cia_raw,
      tau2_raw = tau2_raw,
      tau2 = tau2,
      sigma2_within = sigma2_within,
      boundary = boundary,
      n_eligible = n_eligible,
      lwr = lwr,
      upr = upr
    ))
  }

  for (j in seq_len(prep$n_methods - 1L)) {
    for (k in (j + 1L):prep$n_methods) {
      N <- as.integer(raw$n_eligible[j, k])
      n_eligible[j, k] <- N
      n_eligible[k, j] <- N
      if (N <= 0L) {
        next
      }
      A <- as.numeric(raw$A_bar[j, k])
      B <- as.numeric(raw$B_bar[j, k])
      payload <- .mc_cia_component_payload(
        sigma2_within = A,
        between_term = B,
        estimator = estimator
      )
      est[j, k] <- payload$cia
      est[k, j] <- payload$cia
      cia_raw[j, k] <- payload$cia_raw
      cia_raw[k, j] <- payload$cia_raw
      tau2_raw[j, k] <- payload$tau2_raw
      tau2_raw[k, j] <- payload$tau2_raw
      tau2[j, k] <- payload$tau2
      tau2[k, j] <- payload$tau2
      sigma2_within[j, k] <- payload$sigma2_within
      sigma2_within[k, j] <- payload$sigma2_within
      boundary[j, k] <- payload$boundary
      boundary[k, j] <- payload$boundary
      if (isTRUE(return_delta)) {
        ci <- .mc_cia_pairwise_delta_ci(
          A = A,
          B = B,
          N = N,
          S11 = as.numeric(raw$S11[j, k]),
          S22 = as.numeric(raw$S22[j, k]),
          S33 = as.numeric(raw$S33[j, k]),
          S12 = as.numeric(raw$S12[j, k]),
          S13 = as.numeric(raw$S13[j, k]),
          S23 = as.numeric(raw$S23[j, k]),
          reference_pair = FALSE,
          conf_level = conf_level,
          center = payload$cia
        )
        lwr[j, k] <- ci$lwr
        lwr[k, j] <- ci$lwr
        upr[j, k] <- ci$upr
        upr[k, j] <- ci$upr
      }
    }
  }
  list(
    est = est,
    cia_raw = cia_raw,
    tau2_raw = tau2_raw,
    tau2 = tau2,
    sigma2_within = sigma2_within,
    boundary = boundary,
    n_eligible = n_eligible,
    lwr = lwr,
    upr = upr
  )
}

.mc_cia_assert_estimable <- function(prep,
                                     diagnostics,
                                     scope = c("overall", "pairwise"),
                                     n_threads = 1L) {
  scope <- match.arg(scope)
  if (!any(diagnostics$n_replicate_pairs > 0L)) {
    abort_bad_arg(
      "data",
      message = "does not contain any subject-method cell with at least two replicated readings.",
      .hint = "CIA requires replicated readings within method."
    )
  }

  if (!identical(scope, "overall")) {
    est <- .mc_cia_pairwise_payload(
      prep,
      diagnostics = diagnostics,
      estimator = "vc_constrained",
      n_threads = n_threads
    )$est
    if (!any(is.finite(est[upper.tri(est, diag = FALSE)]))) {
      abort_bad_arg(
        "data",
        message = "does not contain enough replicated and cross-method information to estimate any pairwise CIA value."
      )
    }
  }
  invisible(NULL)
}

.mc_cia_build_overall_from_moments <- function(prep, diagnostics, moments, estimator) {
  out <- data.frame(
    coefficient = prep$coefficient,
    estimator = estimator,
    cia = moments$cia,
    cia_raw = moments$cia_raw,
    tau2 = moments$tau2,
    tau2_raw = moments$tau2_raw,
    boundary = moments$boundary,
    reference = if (isTRUE(prep$has_reference)) prep$reference else NA_character_,
    n_methods = prep$n_methods,
    n_subjects = moments$n,
    n_obs = diagnostics$n_obs,
    K = moments$K,
    numerator_term = moments$numerator_term,
    denominator_term = moments$denominator_term,
    within_msd = moments$numerator_term,
    between_msd = moments$denominator_term,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  class(out) <- c("cia_overall", "cia", "data.frame")
  attr(out, "method") <- "Coefficient of individual agreement"
  attr(out, "description") <- if (isTRUE(prep$has_reference)) {
    "Overall reference-scaled coefficient of individual agreement"
  } else {
    "Overall coefficient of individual agreement"
  }
  attr(out, "package") <- "matrixCorr"
  attr(out, "scope") <- "overall"
  attr(out, "reference") <- prep$reference
  attr(out, "coefficient") <- prep$coefficient
  attr(out, "estimator") <- estimator
  attr(out, "boundary") <- moments$boundary
  attr(out, "diagnostics") <- c(diagnostics, list(K = moments$K))
  attr(out, "overall_moments") <- moments
  out
}

.mc_cia_warn_small_pairwise_n <- function(payload) {
  off_diag <- upper.tri(payload$n_eligible)
  small <- payload$n_eligible[off_diag] > 0L & payload$n_eligible[off_diag] < 10L
  if (any(small, na.rm = TRUE)) {
    cli::cli_warn(
      "Some pairwise CIA estimates use fewer than 10 eligible subjects.",
      class = "matrixCorr_cia_small_n"
    )
  }
}

.mc_cia_pairwise_diagnostics <- function(diagnostics, payload) {
  diagnostics$n_complete <- payload$n_eligible
  diagnostics$n_eligible_subjects <- payload$n_eligible
  diagnostics
}

.mc_cia_build_pairwise_matrix <- function(prep, diagnostics, estimator, n_threads = 1L) {
  payload <- .mc_cia_pairwise_payload(
    prep,
    diagnostics = diagnostics,
    estimator = estimator,
    n_threads = n_threads
  )
  .mc_cia_warn_small_pairwise_n(payload)
  pair_diag <- .mc_cia_pairwise_diagnostics(diagnostics, payload)
  out <- .mc_new_corr_matrix(
    mat = payload$est,
    estimator_class = "cia",
    method = "Coefficient of individual agreement",
    description = if (isTRUE(prep$has_reference)) {
      "Reference-scaled pairwise coefficient of individual agreement"
    } else {
      "Pairwise coefficient of individual agreement"
    },
    diagnostics = pair_diag,
    symmetric = TRUE,
    extra_attrs = list(
      scope = "pairwise",
      reference = prep$reference,
      coefficient = prep$coefficient,
      estimator = estimator,
      cia_raw = payload$cia_raw,
      tau2_raw = payload$tau2_raw,
      tau2 = payload$tau2,
      sigma2_within = payload$sigma2_within,
      boundary = payload$boundary
    )
  )
  class(out) <- unique(c("cia", class(out)))
  out
}

.mc_cia_build_pairwise_ci <- function(prep,
                                      diagnostics,
                                      estimator,
                                      inference,
                                      boot = NULL,
                                      conf_level,
                                      n_threads = 1L) {
  payload <- .mc_cia_pairwise_payload(
    prep,
    diagnostics = diagnostics,
    estimator = estimator,
    conf_level = conf_level,
    return_delta = identical(inference, "delta"),
    n_threads = n_threads
  )
  .mc_cia_warn_small_pairwise_n(payload)
  out <- list(
    est = payload$est,
    lwr.ci = if (identical(inference, "delta")) payload$lwr else boot$lwr,
    upr.ci = if (identical(inference, "delta")) payload$upr else boot$upr
  )
  class(out) <- c("cia", "cia_ci")
  attr(out, "method") <- "Coefficient of individual agreement"
  attr(out, "description") <- if (isTRUE(prep$has_reference)) {
    "Reference-scaled pairwise coefficient of individual agreement with confidence intervals"
  } else {
    "Pairwise coefficient of individual agreement with confidence intervals"
  }
  attr(out, "package") <- "matrixCorr"
  attr(out, "scope") <- "pairwise"
  attr(out, "reference") <- prep$reference
  attr(out, "coefficient") <- prep$coefficient
  attr(out, "estimator") <- estimator
  attr(out, "cia_raw") <- payload$cia_raw
  attr(out, "tau2_raw") <- payload$tau2_raw
  attr(out, "tau2") <- payload$tau2
  attr(out, "sigma2_within") <- payload$sigma2_within
  attr(out, "boundary") <- payload$boundary
  attr(out, "diagnostics") <- .mc_cia_pairwise_diagnostics(diagnostics, payload)
  attr(out, "conf.level") <- conf_level
  attr(out, "ci.method") <- if (identical(inference, "delta")) {
    "delta_normal"
  } else {
    "subject_bootstrap_percentile"
  }
  if (!is.null(boot)) {
    attr(out, "bootstrap") <- boot[c("n_successful", "B")]
  }
  out
}

.mc_cia_bootstrap_vectors <- function(prep, sampled_subjects) {
  total <- sum(lengths(prep$subject_rows[sampled_subjects]))
  y <- numeric(total)
  subject <- integer(total)
  method <- integer(total)
  replicate <- integer(total)
  pos <- 1L
  for (draw in seq_along(sampled_subjects)) {
    idx <- prep$subject_rows[[sampled_subjects[[draw]]]]
    n_i <- length(idx)
    rng <- pos:(pos + n_i - 1L)
    y[rng] <- prep$y[idx]
    method[rng] <- prep$method_code[idx]
    replicate[rng] <- prep$replicate_code[idx]
    subject[rng] <- draw
    pos <- pos + n_i
  }
  list(y = y, subject = subject, method = method, replicate = replicate)
}

.mc_cia_bootstrap_prep <- function(prep, boot_vec) {
  subject_rows <- split(seq_along(boot_vec$y), boot_vec$subject)
  by_subject_method <- lapply(subject_rows, function(idx) split(boot_vec$y[idx], boot_vec$method[idx]))
  list(
    y = boot_vec$y,
    subject_code = boot_vec$subject,
    method_code = boot_vec$method,
    replicate_code = boot_vec$replicate,
    subject_rows = subject_rows,
    by_subject_method = by_subject_method,
    method_levels = prep$method_levels,
    subject_levels = seq_along(subject_rows),
    replicate_levels = prep$replicate_levels,
    n_methods = prep$n_methods,
    n_subjects = length(subject_rows),
    has_reference = prep$has_reference,
    reference_index = prep$reference_index,
    reference = prep$reference,
    coefficient = prep$coefficient,
    n_original_rows = length(boot_vec$y),
    n_complete_rows = length(boot_vec$y),
    na_rows_removed = 0L
  )
}

.mc_cia_bootstrap_overall <- function(prep, estimator, B, conf_level, seed = NULL, n_threads = 1L) {
  .mc_eval_with_seed(seed, {
    alpha <- (1 - conf_level) / 2
    est <- rep(NA_real_, B)
    for (b in seq_len(B)) {
      sampled <- sample.int(prep$n_subjects, size = prep$n_subjects, replace = TRUE)
      boot_vec <- .mc_cia_bootstrap_vectors(prep, sampled)
      boot_prep <- .mc_cia_bootstrap_prep(prep, boot_vec)
      est[[b]] <- tryCatch(
        .mc_cia_overall_moments(
          boot_prep,
          estimator = estimator,
          n_threads = n_threads
        )$cia,
        error = function(...) NA_real_
      )
    }
    keep <- is.finite(est)
    ci <- if (sum(keep) >= 2L) {
      stats::quantile(est[keep], probs = c(alpha, 1 - alpha), names = FALSE, type = 6, na.rm = TRUE)
    } else {
      c(NA_real_, NA_real_)
    }
    list(
      lwr = as.numeric(ci[[1L]]),
      upr = as.numeric(ci[[2L]]),
      n_successful = sum(keep),
      B = B
    )
  })
}

.mc_cia_bootstrap_pairwise <- function(prep, estimator, B, conf_level, seed = NULL, n_threads = 1L) {
  .mc_eval_with_seed(seed, {
    alpha <- (1 - conf_level) / 2
    arr <- array(NA_real_, dim = c(prep$n_methods, prep$n_methods, B))
    for (b in seq_len(B)) {
      sampled <- sample.int(prep$n_subjects, size = prep$n_subjects, replace = TRUE)
      boot_vec <- .mc_cia_bootstrap_vectors(prep, sampled)
      boot_prep <- .mc_cia_bootstrap_prep(prep, boot_vec)
      arr[, , b] <- .mc_cia_pairwise_payload(
        boot_prep,
        estimator = estimator,
        n_threads = n_threads
      )$est
    }

    lwr <- matrix(NA_real_, prep$n_methods, prep$n_methods, dimnames = .mc_square_dimnames(prep$method_levels))
    upr <- lwr
    keep_any <- 0L
    for (i in seq_len(prep$n_methods)) {
      for (j in seq_len(prep$n_methods)) {
        vals <- arr[i, j, ]
        ok <- is.finite(vals)
        if (sum(ok) >= 2L) {
          qs <- stats::quantile(vals[ok], probs = c(alpha, 1 - alpha), names = FALSE, type = 6, na.rm = TRUE)
          lwr[i, j] <- as.numeric(qs[[1L]])
          upr[i, j] <- as.numeric(qs[[2L]])
          keep_any <- keep_any + 1L
        }
      }
    }
    diag(lwr) <- 1
    diag(upr) <- 1
    list(
      lwr = lwr,
      upr = upr,
      n_successful = sum(vapply(seq_len(B), function(b) any(is.finite(arr[, , b])), logical(1))),
      B = B
    )
  })
}

.mc_cia_as_corr_result <- function(x) {
  if (inherits(x, "corr_result")) {
    return(x)
  }
  if (inherits(x, "cia_ci")) {
    ci_attr <- list(
      est = as.matrix(x$est),
      lwr.ci = as.matrix(x$lwr.ci),
      upr.ci = as.matrix(x$upr.ci),
      conf.level = attr(x, "conf.level", exact = TRUE),
      ci.method = attr(x, "ci.method", exact = TRUE) %||% "subject_bootstrap_percentile"
    )
    return(.mc_new_corr_matrix(
      mat = as.matrix(x$est),
      estimator_class = "cia",
      method = attr(x, "method", exact = TRUE) %||% "Coefficient of individual agreement",
      description = attr(x, "description", exact = TRUE) %||% "Pairwise coefficient of individual agreement with confidence intervals",
      diagnostics = attr(x, "diagnostics", exact = TRUE),
      ci = ci_attr,
      conf.level = attr(x, "conf.level", exact = TRUE),
      symmetric = TRUE,
      extra_attrs = list(
        scope = attr(x, "scope", exact = TRUE),
        reference = attr(x, "reference", exact = TRUE),
        coefficient = attr(x, "coefficient", exact = TRUE),
        estimator = attr(x, "estimator", exact = TRUE),
        cia_raw = attr(x, "cia_raw", exact = TRUE),
        tau2_raw = attr(x, "tau2_raw", exact = TRUE),
        tau2 = attr(x, "tau2", exact = TRUE),
        sigma2_within = attr(x, "sigma2_within", exact = TRUE),
        boundary = attr(x, "boundary", exact = TRUE)
      )
    ))
  }
  abort_bad_arg("x", message = "must be a CIA object.")
}

.mc_cia_has_boundary <- function(x) {
  boundary <- attr(x, "boundary", exact = TRUE)
  if (is.null(boundary) && is.data.frame(x) && "boundary" %in% names(x)) {
    boundary <- x$boundary
  }
  isTRUE(any(boundary, na.rm = TRUE))
}

.mc_cia_boundary_note <- function() {
  paste(
    "Note: the raw method-of-moments estimator implied negative inter-method variability",
    "for at least one comparison, so the reported CIA was constrained to the",
    "valid variance-component parameter space."
  )
}

.mc_cia_summary_add_pairwise_fields <- function(df, corr_obj) {
  raw_mat <- attr(corr_obj, "cia_raw", exact = TRUE)
  tau2_raw_mat <- attr(corr_obj, "tau2_raw", exact = TRUE)
  tau2_mat <- attr(corr_obj, "tau2", exact = TRUE)
  boundary_mat <- attr(corr_obj, "boundary", exact = TRUE)
  if (is.null(raw_mat) || !all(c("item1", "item2") %in% names(df))) {
    return(df)
  }
  idx <- cbind(
    match(df$item1, rownames(raw_mat)),
    match(df$item2, colnames(raw_mat))
  )
  df$cia_raw <- raw_mat[idx]
  df$tau2_raw <- tau2_raw_mat[idx]
  df$tau2 <- tau2_mat[idx]
  df$boundary <- as.logical(boundary_mat[idx])
  df
}

#' @rdname cia
#' @method summary cia
#' @param object A `cia` object.
#' @param digits Integer; number of decimal places for estimates.
#' @param ci_digits Integer; number of decimal places for confidence limits.
#' @param ... Unused.
#' @export
summary.cia <- function(object,
                        digits = 4,
                        ci_digits = 3,
                        ...) {
  if (inherits(object, "cia_overall") && is.data.frame(object)) {
    out <- as.data.frame(object, stringsAsFactors = FALSE, check.names = FALSE)
    out <- .mc_finalize_summary_df(out, class_name = "summary.cia")
    attr(out, "result_type") <- "overall"
    attr(out, "digits") <- digits
    attr(out, "ci_digits") <- ci_digits
    attr(out, "conf.level") <- attr(object, "conf.level", exact = TRUE) %||% NA_real_
    attr(out, "coefficient") <- attr(object, "coefficient", exact = TRUE) %||% NA_character_
    attr(out, "reference") <- attr(object, "reference", exact = TRUE) %||% NA_character_
    attr(out, "estimator") <- attr(object, "estimator", exact = TRUE) %||% NA_character_
    attr(out, "boundary") <- attr(object, "boundary", exact = TRUE) %||% FALSE
    return(out)
  }
  corr_obj <- .mc_cia_as_corr_result(object)
  out <- summary.corr_matrix(corr_obj, ...)
  if ("fisher_z" %in% names(out)) {
    out$fisher_z <- NULL
  }
  out <- .mc_cia_summary_add_pairwise_fields(out, corr_obj)
  top_tbl <- attr(out, "top_results", exact = TRUE)
  if (is.data.frame(top_tbl) && "fisher_z" %in% names(top_tbl)) {
    top_tbl$fisher_z <- NULL
  }
  if (is.data.frame(top_tbl)) {
    top_tbl <- .mc_cia_summary_add_pairwise_fields(top_tbl, corr_obj)
    attr(out, "top_results") <- top_tbl
  }
  class(out) <- unique(c("summary.cia", class(out)))
  attr(out, "estimator") <- attr(corr_obj, "estimator", exact = TRUE) %||% NA_character_
  attr(out, "boundary") <- attr(corr_obj, "boundary", exact = TRUE)
  attr(out, "summary_title") <- "Coefficient of individual agreement summary"
  out
}

#' @rdname cia
#' @method print cia
#' @param x A `cia` object.
#' @param digits Integer; number of decimal places for estimates.
#' @param n Optional preview row threshold.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns.
#' @param width Optional display width.
#' @param show_ci One of `"yes"` or `"no"`.
#' @param ... Additional arguments passed to downstream print helpers.
#' @export
print.cia <- function(x,
                      digits = 4,
                      n = NULL,
                      topn = NULL,
                      max_vars = NULL,
                      width = NULL,
                      show_ci = NULL,
                      ...) {
  if (inherits(x, "cia_overall") && is.data.frame(x)) {
    return(invisible(print.summary.cia(
      summary.cia(x),
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
    .mc_cia_as_corr_result(x),
    header = if (isTRUE(nzchar(attr(.mc_cia_as_corr_result(x), "reference", exact = TRUE) %||% ""))) {
      "Coefficient of individual agreement matrix"
    } else {
      "Coefficient of individual agreement matrix"
    },
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
  if (identical(attr(.mc_cia_as_corr_result(x), "estimator", exact = TRUE), "vc_constrained") &&
      .mc_cia_has_boundary(.mc_cia_as_corr_result(x))) {
    cat("\n", .mc_cia_boundary_note(), "\n", sep = "")
  }
  invisible(x)
}

#' @rdname cia
#' @method print summary.cia
#' @param x A `summary.cia` object.
#' @export
print.summary.cia <- function(x,
                              digits = NULL,
                              n = NULL,
                              topn = NULL,
                              max_vars = NULL,
                              width = NULL,
                              show_ci = NULL,
                              ...) {
  if (!identical(attr(x, "result_type", exact = TRUE), "overall")) {
    print.summary.corr_result(
      x,
      digits = digits,
      n = n,
      topn = topn,
      max_vars = max_vars,
      width = width,
      show_ci = show_ci,
      ...
    )
    if (identical(attr(x, "estimator", exact = TRUE), "vc_constrained") &&
        .mc_cia_has_boundary(x)) {
      cat("\n", .mc_cia_boundary_note(), "\n", sep = "")
    }
    return(invisible(x))
  }

  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )
  digits <- .mc_coalesce(digits, .mc_coalesce(attr(x, "digits", exact = TRUE), 4L))
  ci_digits <- .mc_coalesce(attr(x, "ci_digits", exact = TRUE), digits)

  digest <- c(
    coefficient = x$coefficient[[1L]],
    estimator = x$estimator[[1L]],
    cia = formatC(as.numeric(x$cia[[1L]]), format = "f", digits = digits),
    n_methods = .mc_count_fmt(x$n_methods[[1L]]),
    n_subjects = .mc_count_fmt(x$n_subjects[[1L]]),
    n_obs = .mc_count_fmt(x$n_obs[[1L]])
  )
  if (!is.na(x$reference[[1L]]) && nzchar(x$reference[[1L]])) {
    digest <- c(digest, reference = x$reference[[1L]])
  }
  if (isTRUE(identical(show_ci, "yes")) &&
      all(c("lwr.ci", "upr.ci") %in% names(x)) &&
      is.finite(x$lwr.ci[[1L]]) &&
      is.finite(x$upr.ci[[1L]])) {
    digest <- c(
      digest,
      ci = sprintf("%g%%", 100 * suppressWarnings(as.numeric(attr(x, "conf.level", exact = TRUE))))
    )
  }
  .mc_print_named_digest(digest, header = "Coefficient of individual agreement")
  cat("\n")

  preview <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  num_cols <- intersect(
    c(
      "cia",
      "cia_raw",
      "tau2",
      "tau2_raw",
      "K",
      "numerator_term",
      "denominator_term",
      "within_msd",
      "between_msd",
      "se"
    ),
    names(preview)
  )
  for (nm in num_cols) preview[[nm]] <- round(preview[[nm]], digits)
  if ("lwr.ci" %in% names(preview)) preview$lwr.ci <- round(preview$lwr.ci, ci_digits)
  if ("upr.ci" %in% names(preview)) preview$upr.ci <- round(preview$upr.ci, ci_digits)
  has_ci_cols <- all(c("lwr.ci", "upr.ci") %in% names(preview))
  if (identical(show_ci, "no") ||
      (isTRUE(has_ci_cols) &&
       !any(is.finite(preview$lwr.ci) | is.finite(preview$upr.ci)))) {
    preview <- preview[, setdiff(names(preview), c("lwr.ci", "upr.ci")), drop = FALSE]
  }

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
  if (identical(attr(x, "estimator", exact = TRUE), "vc_constrained") &&
      .mc_cia_has_boundary(x)) {
    cat("\n", .mc_cia_boundary_note(), "\n", sep = "")
  }
  invisible(x)
}

#' @rdname cia
#' @method plot cia
#' @param x A `cia` object.
#' @param title Optional plot title.
#' @param low_color Color used for lower agreement values.
#' @param high_color Color used for higher agreement values.
#' @param mid_color Midpoint color for pairwise heatmaps.
#' @param value_text_size Text size for overlaid estimate labels.
#' @param ci_text_size Text size for CI labels.
#' @param show_value Logical; whether to overlay numeric values.
#' @export
plot.cia <- function(x,
                     title = NULL,
                     low_color = "indianred1",
                     high_color = "steelblue1",
                     mid_color = "white",
                     value_text_size = 4,
                     ci_text_size = 3,
                     show_value = TRUE,
                     ...) {
  if (inherits(x, "cia_overall") && is.data.frame(x)) {
    check_bool(show_value, arg = "show_value")
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg ggplot2} is required for plotting.")
    }
    df <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
    est <- as.numeric(df$cia[[1L]])
    lwr <- if ("lwr.ci" %in% names(df)) as.numeric(df$lwr.ci[[1L]]) else NA_real_
    upr <- if ("upr.ci" %in% names(df)) as.numeric(df$upr.ci[[1L]]) else NA_real_
    x_min <- min(c(0, est, lwr), na.rm = TRUE)
    x_max <- max(c(1, est, upr), na.rm = TRUE)
    if (!is.finite(x_min) || !is.finite(x_max)) {
      x_min <- 0
      x_max <- max(1, est + 0.1)
    }
    if (identical(x_min, x_max)) {
      x_min <- x_min - 0.1
      x_max <- x_max + 0.1
    }
    pad <- max(0.04, 0.06 * (x_max - x_min))
    point_color <- if (is.finite(est) && est < 0.5) low_color else high_color
    x_limits <- c(x_min - pad, x_max + pad)
    band_df <- data.frame(
      xmin = c(0, 0.445, 0.8, 1),
      xmax = c(0.445, 0.8, 1, x_limits[[2L]]),
      ymin = 0.7,
      ymax = 1.3,
      fill = c("#f6d6d6", "#f4ead0", "#dbeedc", "#d9e7f7")
    )
    band_df <- band_df[band_df$xmin < x_limits[[2L]], , drop = FALSE]
    band_df$xmin <- pmax(band_df$xmin, x_limits[[1L]])
    band_df$xmax <- pmin(band_df$xmax, x_limits[[2L]])

    guide_df <- data.frame(
      x = c(0.445, 0.8),
      label = c("0.445", "0.800"),
      text = c("IEC-linked reference", "strong-agreement reference")
    )
    guide_df <- guide_df[guide_df$x >= x_limits[[1L]] & guide_df$x <= x_limits[[2L]], , drop = FALSE]

    label_df <- data.frame(
      x = c(0.2225, (0.445 + 0.8) / 2, 0.9),
      y = c(1.18, 1.18, 1.18),
      label = c("lower", "intermediate", "strong")
    )
    label_df <- label_df[label_df$x >= x_limits[[1L]] & label_df$x <= x_limits[[2L]], , drop = FALSE]

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$cia, y = 1)) +
      ggplot2::geom_rect(
        data = band_df,
        ggplot2::aes(
          xmin = .data$xmin,
          xmax = .data$xmax,
          ymin = .data$ymin,
          ymax = .data$ymax
        ),
        inherit.aes = FALSE,
        fill = band_df$fill,
        color = NA,
        alpha = 0.85
      ) +
      ggplot2::geom_segment(
        inherit.aes = FALSE,
        data = data.frame(x = x_limits[[1L]], xend = x_limits[[2L]], y = 1, yend = 1),
        ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
        linewidth = 6,
        color = "white",
        lineend = "round"
      ) +
      ggplot2::geom_point(size = 3.6, color = point_color)
    if (is.finite(lwr) && is.finite(upr)) {
      p <- p + ggplot2::geom_segment(
        ggplot2::aes(x = .data$lwr.ci, xend = .data$upr.ci, y = 1, yend = 1),
        linewidth = 1.4,
        color = point_color
      )
    }
    if (nrow(guide_df)) {
      p <- p + ggplot2::geom_vline(
        data = guide_df,
        ggplot2::aes(xintercept = .data$x),
        inherit.aes = FALSE,
        color = "gray45",
        linetype = "22",
        linewidth = 0.5
      ) +
        ggplot2::geom_text(
          data = guide_df,
          ggplot2::aes(x = .data$x, y = 0.76, label = .data$label),
          inherit.aes = FALSE,
          size = 3,
          color = "gray35",
          vjust = 1
        )
    }
    if (nrow(label_df)) {
      p <- p + ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        inherit.aes = FALSE,
        size = 3.2,
        color = "gray25",
        fontface = "bold"
      )
    }
    if (isTRUE(show_value)) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2f", .data$cia)),
        nudge_y = 0.16,
        size = value_text_size
      )
    }
    subtitle <- if (!is.na(df$reference[[1L]]) && nzchar(df$reference[[1L]])) {
      sprintf("Overall reference-scaled CIA using %s as the reference method", df$reference[[1L]])
    } else {
      "Overall no-reference CIA from within-method replicate disagreement"
    }
    p <- p +
      ggplot2::scale_x_continuous(limits = x_limits) +
      ggplot2::scale_y_continuous(NULL, breaks = NULL) +
      ggplot2::labs(
        title = title %||% "Coefficient of individual agreement",
        subtitle = subtitle,
        x = "CIA",
        y = NULL,
        caption = "Guide lines at 0.445 and 0.8 are contextual references, not fixed pass/fail thresholds."
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        plot.caption = ggplot2::element_text(hjust = 0)
      ) +
      ggplot2::coord_cartesian(clip = "off")
    return(p)
  }

  .mc_plot_corr_result(
    .mc_cia_as_corr_result(x),
    title = title %||% "Coefficient of individual agreement heatmap",
    low_color = low_color,
    high_color = high_color,
    mid_color = mid_color,
    value_text_size = value_text_size,
    ci_text_size = ci_text_size,
    show_value = show_value,
    ...
  )
}

#' @rdname cia
#' @method print cia_ci
#' @export
print.cia_ci <- function(x, ...) {
  print.cia(x, ...)
}
