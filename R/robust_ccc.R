#' MCD based robust concordance correlation coefficient
#'
#' @description
#' Computes pairwise robust concordance correlation coefficients from the
#' numeric columns of a matrix or data frame. The estimator replaces the
#' classical means, variances, and covariance in Lin's concordance coefficient
#' with joint minimum covariance determinant estimates.
#'
#' Classical concordance can be strongly distorted by outliers because it
#' depends directly on means, variances, and covariance. The minimum covariance
#' determinant is a high breakdown estimator of multivariate location and
#' scatter. It is designed primarily for data with an approximately
#' elliptically symmetric and unimodal core distribution. It is not universally
#' robust or distribution free.
#'
#' @param data A numeric matrix or data frame containing at least two numeric
#'   columns. Non numeric columns in a data frame are ignored.
#' @param alpha Numeric scalar in \eqn{[0.5, 1]}. This controls the fraction of
#'   observations used by the minimum covariance determinant estimator.
#'   Approximately \code{alpha * n} observations define the high quality
#'   subset. Lower values provide greater resistance to contamination. The
#'   approximate trimmed fraction is \code{1 - alpha}. Default \code{0.75}.
#' @param ci Logical. If \code{TRUE}, attach paired nonparametric percentile
#'   bootstrap confidence intervals. Default \code{FALSE}.
#' @param conf_level Confidence level used when \code{ci = TRUE}. Default
#'   \code{0.95}.
#' @param n_boot Positive integer giving the number of bootstrap resamples used
#'   when \code{ci = TRUE}. Default \code{500}.
#' @param seed Optional positive integer controlling both the stochastic
#'   FastMCD point estimates and bootstrap resampling. Stable offsets are used
#'   for different variable pairs. If \code{NULL}, the current R random number
#'   stream is used.
#' @param na_method Character scalar controlling missing data handling.
#'   \code{"error"} rejects missing, \code{NaN}, and infinite values.
#'   \code{"complete"} removes rows that are incomplete in any retained
#'   numeric column before all estimates are calculated. \code{"pairwise"}
#'   selects finite rows separately for each variable pair.
#' @param output Output representation. One of \code{"matrix"},
#'   \code{"sparse"}, or \code{"edge_list"}.
#' @param threshold Nonnegative absolute value filter for sparse and edge list
#'   outputs. It must be zero for matrix output.
#' @param diag Logical. Whether sparse and edge list outputs include diagonal
#'   entries.
#' @param mcd_nsamp Positive integer giving the number of subsets used for
#'   initial estimates in the FastMCD search. Passed to the \code{nsamp}
#'   argument of [robustbase::covMcd()]. Default \code{500L}. This is a search
#'   count, not the subset size controlled by \code{alpha} or the number of
#'   bootstrap resamples controlled by \code{n_boot}. Only numeric counts up to
#'   \code{.Machine$integer.max} are supported, not the character modes of
#'   \code{covMcd()}.
#'
#' @details
#' For each paired sample, let \eqn{\hat\mu_{X,MCD}} and
#' \eqn{\hat\mu_{Y,MCD}} be the components of the joint MCD location estimate.
#' Let \eqn{\hat\sigma^2_{X,MCD}}, \eqn{\hat\sigma^2_{Y,MCD}}, and
#' \eqn{\hat\sigma_{XY,MCD}} be entries of the joint MCD scatter estimate. The
#' robust concordance estimate is
#' \deqn{
#' \hat\rho_R =
#' \frac{2\hat\sigma_{XY,MCD}}
#' {\hat\sigma^2_{X,MCD} + \hat\sigma^2_{Y,MCD} +
#' (\hat\mu_{X,MCD} - \hat\mu_{Y,MCD})^2}.
#' }
#'
#' The joint location and scatter are fitted with \code{robustbase::covMcd()}
#' using \code{raw.only = FALSE} and \code{nsamp = mcd_nsamp}. Consequently, this
#' function uses the final reweighted MCD estimates rather than the raw MCD
#' subset estimates. No population moment rescaling is applied.
#'
#' The default \code{mcd_nsamp = 500L} retains the published software
#' convention. Larger counts increase the search effort and generally increase
#' computation time. Changing the count can change the numerical MCD solution,
#' and hence the coefficient and confidence interval, even with the same
#' \code{seed}. The same count is used for every pair and every bootstrap
#' refit. With a nonnull \code{seed}, repeated calls are reproducible for fixed
#' data, arguments, and seed.
#'
#' Reweighting is part of this estimator and is not a user option. Even when
#' \code{alpha = 1}, the final \code{covMcd()} location and scatter can still
#' be reweighted and must not be replaced by ordinary sample moments.
#'
#' At least three finite paired observations are required for basic
#' computability. An MCD fit can still be undefined for a very small or
#' geometrically degenerate sample. If either robust marginal scatter is
#' effectively zero, the corresponding coefficient is undefined and
#' \code{NA} is returned. A diagonal entry is one only when its robust scatter
#' is estimable and nondegenerate.
#'
#' Confidence intervals are a matrixCorr inference extension. They are not an
#' interval proposed by Bulut, Zobu, and Sağlam. Each bootstrap replicate
#' samples paired rows with replacement and recalculates the complete MCD
#' estimator. Failed or nonfinite replicate estimates are omitted before the
#' percentile limits are calculated. The \code{diagnostics} attribute contains
#' \code{n_boot_success} when intervals are requested. Interval calculation is
#' expensive because every variable pair requires \code{n_boot} additional MCD
#' fits. Its leading cost is
#' \eqn{O(p^2 B \mathop{\mathrm{cost}}(\text{bivariate MCD}))}.
#'
#' This estimator retains the form and interpretation of Lin's concordance
#' coefficient. It differs from the robust concordance variants of King and
#' Chinchilli, which use alternative distance or loss functions. It also differs
#' from the Bayesian robust concordance estimator of Feng, Baumgartner, and
#' Svetnik, which uses heavy tailed multivariate modelling.
#'
#' @return A standard matrixCorr correlation result. Dense output has class
#'   \code{c("corr_matrix", "robust_ccc", "corr_result", "matrix")}.
#'   Confidence intervals are stored in the standard \code{ci} attribute with
#'   \code{ci.method = "percentile_bootstrap"}. The result records
#'   \code{alpha = alpha}, \code{mcd_nsamp = mcd_nsamp}, and
#'   \code{mcd_estimate = "reweighted"}.
#'
#' @author Thiago de Paula Oliveira
#'
#' @references
#' Bulut H, Zobu M, Sağlam V (2026). A High-Breakdown MCD-Based Robust
#' Concordance Correlation Coefficient. Mathematics 14(1), 196.
#' \doi{10.3390/math14010196}
#'
#' Hubert M, Debruyne M, Rousseeuw PJ (2018). Minimum covariance determinant
#' and extensions. WIREs Computational Statistics 10, e1421.
#' \doi{10.1002/wics.1421}
#'
#' Rousseeuw PJ, Van Driessen K (1999). A Fast Algorithm for the Minimum
#' Covariance Determinant Estimator. Technometrics 41, 212 to 223.
#' \doi{10.1080/00401706.1999.10485670}
#'
#' King TS, Chinchilli VM (2001). Robust estimators of the concordance
#' correlation coefficient. Journal of Biopharmaceutical Statistics 11, 83 to
#' 105. \doi{10.1081/BIP-100107651}
#'
#' Feng D, Baumgartner R, Svetnik V (2015). A Robust Bayesian Estimate of the
#' Concordance Correlation Coefficient. Journal of Biopharmaceutical Statistics
#' 25, 490 to 507. \doi{10.1080/10543406.2014.920342}
#'
#' Crawford SB and colleagues (2007). Computer programs for the concordance
#' correlation coefficient. Computer Methods and Programs in Biomedicine 88,
#' 62 to 74. \doi{10.1016/j.cmpb.2007.07.003}
#'
#' @seealso [robustbase::covMcd()], [ccc()], [pbcor()], [wincor()]
#'
#' @examples
#' # Three measurement methods with an approximately normal core
#' set.seed(42)
#' x <- rnorm(40)
#' measurements <- cbind(
#'   reference = x,
#'   method_a = x + rnorm(40, sd = 0.2),
#'   method_b = 0.4 + 1.2 * x + rnorm(40, sd = 0.3)
#' )
#' fit <- robust_ccc(measurements, seed = 42)
#' print(fit)
#' summary(fit)
#' estimate(fit)
#' coef(fit)
#' tidy(fit)
#' plot(fit)
#'
#' # Compare classical and robust CCC before and after discordant contamination
#' contaminated <- measurements
#' contaminated[1, "reference"] <- 9
#' contaminated[1, "method_a"] <- -8
#' fit_contaminated <- robust_ccc(contaminated, seed = 42)
#' rbind(
#'   clean = c(
#'     classical = estimate(ccc(measurements))[1, 2],
#'     robust = estimate(fit)[1, 2]
#'   ),
#'   contaminated = c(
#'     classical = estimate(ccc(contaminated))[1, 2],
#'     robust = estimate(fit_contaminated)[1, 2]
#'   )
#' )
#'
#' # Alpha controls the retained subset fraction, not the search count
#' fit_half <- robust_ccc(contaminated, alpha = 0.5, seed = 42)
#' estimate(fit_half)
#'
#' # Increase the initial subset search count while retaining reweighted MCD
#' fit_more <- robust_ccc(contaminated, mcd_nsamp = 1000L, seed = 42)
#' attr(fit_more, "mcd_nsamp")
#' attr(fit_more, "mcd_estimate")
#'
#' # Global complete rows versus the finite overlap for each pair
#' incomplete <- measurements
#' incomplete[1, "reference"] <- NA
#' incomplete[2, "method_a"] <- Inf
#' fit_complete <- robust_ccc(incomplete, na_method = "complete", seed = 42)
#' fit_pairwise <- robust_ccc(incomplete, na_method = "pairwise", seed = 42)
#' estimate(fit_complete)
#' estimate(fit_pairwise)
#'
#' # Thresholded outputs retain entries with absolute CCC at least 0.9
#' # Every MCD pair is still fitted before output conversion
#' fit_sparse <- robust_ccc(
#'   measurements, seed = 42, output = "sparse", threshold = 0.9, diag = FALSE
#' )
#' fit_edges <- robust_ccc(
#'   measurements, seed = 42, output = "edge_list", threshold = 0.9, diag = FALSE
#' )
#' fit_sparse
#' tidy(fit_edges)
#'
#' # Paired percentile bootstrap intervals require many additional MCD fits
#' \donttest{
#' # A small count keeps this illustration short; the default is 500 resamples
#' fit_ci <- robust_ccc(
#'   measurements, ci = TRUE, conf_level = 0.95, n_boot = 100L, seed = 42
#' )
#' summary(fit_ci)
#' tidy(fit_ci)
#' ci(fit_ci)
#' confint(fit_ci)
#' plot(fit_ci)
#' attr(fit_ci, "diagnostics")$n_complete
#' # Off diagonal entries count finite bootstrap estimates
#' # Valid diagonal intervals are fixed at [1, 1] without refitting
#' attr(fit_ci, "diagnostics")$n_boot_success
#' }
#'
#' # Interactive viewing requires both suggested packages
#' if (interactive() &&
#'     requireNamespace("shiny", quietly = TRUE) &&
#'     requireNamespace("shinyWidgets", quietly = TRUE)) {
#'   view_corr_shiny(fit)
#' }
#'
#' @export
robust_ccc <- function(data,
                       alpha = 0.75,
                       ci = FALSE,
                       conf_level = 0.95,
                       n_boot = 500L,
                       seed = NULL,
                       na_method = c("error", "complete", "pairwise"),
                       output = c("matrix", "sparse", "edge_list"),
                       threshold = 0,
                       diag = TRUE,
                       mcd_nsamp = 500L) {
  check_scalar_numeric(
    alpha,
    arg = "alpha",
    lower = 0.5,
    upper = 1,
    closed_lower = TRUE,
    closed_upper = TRUE
  )
  check_bool(ci, arg = "ci")
  check_scalar_numeric(mcd_nsamp, arg = "mcd_nsamp", lower = 1,
                       upper = .Machine$integer.max)
  mcd_nsamp <- check_scalar_int_pos(mcd_nsamp, arg = "mcd_nsamp")
  if (!is.null(seed)) {
    seed <- check_scalar_int_pos(seed, arg = "seed")
  }
  na_method <- match.arg(na_method)

  desc <- paste0(
    "MCD robust concordance; alpha = ", alpha,
    "; reweighted MCD; NA mode = ", na_method, "."
  )

  .mc_inference_corr_wrapper(
    data = data,
    na_method = na_method,
    ci = ci,
    p_value = FALSE,
    conf_level = conf_level,
    n_threads = 1L,
    output = output,
    threshold = threshold,
    diag = diag,
    estimator_class = "robust_ccc",
    method = "mcd_robust_concordance",
    description = desc,
    kernel_matrix = function(x, n_threads) {
      .mc_robust_ccc_matrix(
        x, alpha = alpha, seed = seed, pairwise = FALSE,
        mcd_nsamp = mcd_nsamp
      )
    },
    kernel_pairwise = function(x, n_threads) {
      .mc_robust_ccc_matrix(
        x, alpha = alpha, seed = seed, pairwise = TRUE,
        mcd_nsamp = mcd_nsamp
      )
    },
    payload_builder = function(x, est, ci, p_value, conf_level, n_boot, seed) {
      .mc_robust_ccc_payload(
        x,
        est = est,
        alpha = alpha,
        conf_level = conf_level,
        n_boot = n_boot,
        seed = seed,
        mcd_nsamp = mcd_nsamp
      )
    },
    min_n = 3L,
    symmetric = TRUE,
    extra_attrs = list(
      alpha = alpha,
      mcd_estimate = "reweighted",
      mcd_nsamp = mcd_nsamp
    ),
    n_boot = n_boot,
    seed = seed
  )
}

.mc_robust_ccc_pair <- function(x, y, alpha = 0.75, seed = NULL,
                                mcd_nsamp = 500L) {
  if (length(x) != length(y) || length(x) < 3L ||
      any(!is.finite(x)) || any(!is.finite(y))) {
    return(NA_real_)
  }

  fit <- tryCatch(
    suppressWarnings(
      .mc_eval_with_seed(
        seed,
        robustbase::covMcd(
          cbind(x, y),
          alpha = alpha,
          raw.only = FALSE,
          nsamp = mcd_nsamp
        )
      )
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(NA_real_)
  }

  mu <- as.numeric(fit$center)
  scatter <- as.matrix(fit$cov)
  if (length(mu) != 2L || !identical(dim(scatter), c(2L, 2L)) ||
      any(!is.finite(mu)) || any(!is.finite(scatter))) {
    return(NA_real_)
  }

  marginal <- c(scatter[1L, 1L], scatter[2L, 2L])
  scatter_scale <- max(abs(scatter))
  scatter_tol <- 100 * .Machine$double.eps *
    max(scatter_scale, .Machine$double.xmin)
  if (any(marginal <= scatter_tol)) {
    return(NA_real_)
  }

  location_gap_sq <- (mu[[1L]] - mu[[2L]])^2
  denominator <- sum(marginal) + location_gap_sq
  denominator_tol <- 100 * .Machine$double.eps * max(
    abs(c(marginal, location_gap_sq)),
    .Machine$double.xmin
  )
  if (!is.finite(denominator) || denominator <= denominator_tol) {
    return(NA_real_)
  }

  value <- 2 * scatter[1L, 2L] / denominator
  if (!is.finite(value)) {
    return(NA_real_)
  }

  range_tol <- sqrt(.Machine$double.eps)
  if (value > 1) {
    if (value <= 1 + range_tol) return(1)
    return(NA_real_)
  }
  if (value < -1) {
    if (value >= -1 - range_tol) return(-1)
    return(NA_real_)
  }
  as.numeric(value)
}

.mc_robust_ccc_matrix <- function(X,
                                  alpha = 0.75,
                                  seed = NULL,
                                  pairwise = FALSE,
                                  min_n = 3L,
                                  mcd_nsamp = 500L) {
  X <- as.matrix(X)
  p <- ncol(X)
  out <- matrix(NA_real_, nrow = p, ncol = p, dimnames = list(colnames(X), colnames(X)))

  pair_id <- 0L
  if (p > 1L) {
    for (i in seq_len(p - 1L)) {
      for (j in seq.int(i + 1L, p)) {
        ok <- if (isTRUE(pairwise)) {
          is.finite(X[, i]) & is.finite(X[, j])
        } else {
          rep(TRUE, nrow(X))
        }
        if (sum(ok) >= min_n) {
          value <- .mc_robust_ccc_pair(
            X[ok, i],
            X[ok, j],
            alpha = alpha,
            seed = .mc_seed_offset(seed, pair_id),
            mcd_nsamp = mcd_nsamp
          )
          out[i, j] <- out[j, i] <- value
        }
        pair_id <- pair_id + 1L
      }
    }
  }

  diagonal_offset <- choose(p, 2L)
  for (i in seq_len(p)) {
    ok <- if (isTRUE(pairwise)) is.finite(X[, i]) else rep(TRUE, nrow(X))
    if (sum(ok) >= min_n) {
      out[i, i] <- .mc_robust_ccc_pair(
        X[ok, i],
        X[ok, i],
        alpha = alpha,
        seed = .mc_seed_offset(seed, diagonal_offset + i - 1L),
        mcd_nsamp = mcd_nsamp
      )
    }
  }
  out
}

.mc_robust_ccc_payload <- function(X,
                                   est,
                                   alpha = 0.75,
                                   conf_level = 0.95,
                                   n_boot = 500L,
                                   seed = NULL,
                                   mcd_nsamp = 500L) {
  est <- as.matrix(est)
  p <- ncol(est)
  dn <- dimnames(est)
  n_complete <- .mc_corr_n_complete(X)
  dimnames(n_complete) <- dn

  ci_lwr <- matrix(NA_real_, p, p, dimnames = dn)
  ci_upr <- matrix(NA_real_, p, p, dimnames = dn)
  n_boot_success <- matrix(0L, p, p, dimnames = dn)

  valid_diagonal <- is.finite(diag(est)) & diag(est) == 1
  diag(ci_lwr)[valid_diagonal] <- 1
  diag(ci_upr)[valid_diagonal] <- 1
  diag(n_boot_success)[valid_diagonal] <- as.integer(n_boot)

  pair_id <- 0L
  if (p > 1L) {
    for (i in seq_len(p - 1L)) {
      for (j in seq.int(i + 1L, p)) {
        ok <- is.finite(X[, i]) & is.finite(X[, j])
        n_ij <- sum(ok)
        n_complete[i, j] <- n_complete[j, i] <- n_ij

        if (n_ij >= 3L && is.finite(est[i, j])) {
          x <- X[ok, i]
          y <- X[ok, j]
          boot_values <- .mc_eval_with_seed(
            .mc_seed_offset(seed, pair_id),
            vapply(
              seq_len(n_boot),
              function(b) {
                idx <- sample.int(n_ij, size = n_ij, replace = TRUE)
                .mc_robust_ccc_pair(
                  x[idx], y[idx], alpha = alpha, seed = NULL,
                  mcd_nsamp = mcd_nsamp
                )
              },
              numeric(1L)
            )
          )
          success <- sum(is.finite(boot_values))
          limits <- .mc_percentile_boot_ci(boot_values, conf_level = conf_level)
          ci_lwr[i, j] <- ci_lwr[j, i] <- limits[[1L]]
          ci_upr[i, j] <- ci_upr[j, i] <- limits[[2L]]
          n_boot_success[i, j] <- n_boot_success[j, i] <- success
        }
        pair_id <- pair_id + 1L
      }
    }
  }

  list(
    diagnostics = list(
      n_complete = n_complete,
      n_boot_success = n_boot_success
    ),
    ci = list(
      est = unclass(est),
      lwr.ci = ci_lwr,
      upr.ci = ci_upr,
      conf.level = conf_level,
      ci.method = "percentile_bootstrap"
    ),
    inference = NULL
  )
}
