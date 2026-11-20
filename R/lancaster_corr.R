#' Lancaster correlation
#'
#' Estimate pairwise dependence through linear association and association
#' between squared magnitudes. The coefficient is the larger absolute
#' correlation and takes values in `[0, 1]`. For a bivariate normal population
#' it equals the absolute Pearson correlation. Zero can occur under dependence.
#'
#' @param data Numeric matrix or data frame with observations in rows.
#' @param method With `"rank"`, transform each margin to normal scores using
#'   `qnorm(rank(x) / (n + 1))` and correlate the scores and their squares.
#'   With `"moment"`, correlate the observations and their squared standardized
#'   values. Both transformations use the observations available for that pair.
#' @param na_method Missing-data policy: error, pairwise deletion, or deletion
#'   of incomplete rows across all columns.
#' @param ci Compute pointwise covariance-bootstrap confidence intervals.
#' @param p_value Compute Monte Carlo permutation tests of pairwise independence.
#' @param conf_level Confidence level for each interval.
#' @param ci_method Covariance-bootstrap interval from equation (18)
#'   (`"conservative"`, default) or equation (17) (`"regular"`).
#' @param n_boot Number of paired bootstrap draws for each interval.
#' @param n_perm Number of random permutations per independence test.
#' @param seed Optional positive integer seed. When omitted, resampling draws
#'   a seed from R's RNG, so `set.seed()` also controls reproducibility.
#' @param n_threads Maximum number of C++ worker threads.
#' @param output Output representation: matrix, sparse matrix, or edge list.
#' @param threshold Minimum absolute estimate retained in sparse or edge output.
#' @param diag Include the diagonal in sparse and edge-list output.
#' @details
#' **Estimation.** For a pair of transformed variables `u` and `v`, the
#' estimate is `max(abs(cor(u, v)), abs(cor(u^2, v^2)))`. The second
#' component can detect association in magnitude even when the first is zero.
#'
#' With `method = "rank"`, each observation with rank \eqn{r} receives the
#' normal score \code{qnorm(r / (n + 1))}, where \eqn{n} is the number of
#' observations used for that pair. These are the plotting positions of
#' Holzmann and Klar. Midpoint positions \code{(r - 0.5) / n} can give
#' different finite-sample estimates. The rank theory assumes continuous
#' margins. The implementation assigns average ranks to ties, including
#' repeated observations in bootstrap samples.
#'
#' With `method = "moment"`, each variable is centered and scaled before
#' squaring. The first component is therefore the ordinary Pearson correlation.
#' The population coefficient requires finite fourth moments and positive
#' variances for both components. The moment-based interval additionally
#' assumes finite eighth moments. For either estimator, a pair with fewer
#' than three usable observations or a constant input to either correlation
#' receives an `NA` estimate.
#'
#' **Confidence intervals.** Inference assumes independent, identically
#' distributed rows. Each bootstrap draw resamples paired observations and
#' recomputes both signed components. The regular interval uses the bootstrap
#' standard error of the component with the larger observed absolute value
#' (the first component in a tie). Its derivation assumes unequal absolute
#' population components. The default conservative interval uses their joint
#' bootstrap covariance to adjust the lower endpoint and retains the regular
#' upper endpoint. Both intervals are pointwise and clipped to `[0, 1]`.
#' Finite-sample coverage can fall below nominal, especially near independence.
#' Off-diagonal limits are `NA` if fewer than two bootstrap draws yield both
#' components. A defined diagonal estimate has both limits equal to one.
#'
#' **Independence test.** For each off-diagonal pair, the test permutes one
#' variable across rows and recomputes the Lancaster estimate. If `b` of
#' `n_perm` permutations yield an estimate at least as large as the observed
#' value (allowing for numerical roundoff), the p-value is
#' `(b + 1) / (n_perm + 1)`. Diagonal p-values are `NA`.
#'
#' @return A matrixCorr result supporting `print()`, `summary()`, `plot()`,
#'   `estimate()`, `coef()`, `tidy()`, `ci()`, and `confint()`. The diagnostics
#'   attribute contains pair counts and the signed `linear` and `quadratic`
#'   components. For sparse results use `base::summary()` for the statistical
#'   summary.
#' @references Holzmann H, Klar B (2025). Lancaster correlation: A new
#'   dependence measure linked to maximum correlation. Scandinavian Journal
#'   of Statistics, 52, 145–169. \doi{10.1111/sjos.12733}.
#'   Phipson B, Smyth GK (2010). Permutation P-values should never be zero.
#'   \doi{10.2202/1544-6115.1585}.
#' @examples
#' # Two sensors record deviations from a target. For each reading of A,
#' # sensor B records the same magnitude once in each direction.
#' a <- rep(c(-3, -2, -1, 1, 2, 3), each = 2)
#' readings <- cbind(A = a, B = a * rep(c(-1, 1), 6))
#' fit <- lancaster_corr(readings, method = "moment")
#' d <- attr(fit, "diagnostics")
#' round(c(pearson = estimate(pearson_corr(readings))["A", "B"],
#'         squared = d$quadratic["A", "B"],
#'         lancaster = estimate(fit)["A", "B"]), 3)
#' # Pearson is 0. The squared deviations agree exactly, giving Lancaster 1.
#'
#' # Which tree measurement is more closely associated with timber volume?
#' # Compare trunk diameter (called Girth in the trees data) and height.
#' timber <- datasets::trees
#' rank_fit <- lancaster_corr(timber, method = "rank", ci = TRUE,
#'                            p_value = TRUE, seed = 42)
#' pairs <- tidy(rank_fit)
#' pairs[pairs$item2 == "Volume", c("item1", "item2", "estimate",
#'                                  "lwr", "upr", "p_value")]
#' # Each interval describes one pair. Each p-value tests independence
#' # within that pair. Compare the estimates to assess relative strength.
#'
#' # A log scale preserves the ordering of these positive volumes.
#' # The rank-based estimates therefore agree on both measurement scales.
#' log_timber <- transform(timber, Volume = log(Volume))
#' all.equal(estimate(rank_fit), estimate(lancaster_corr(log_timber)))
#' @export
lancaster_corr <- function(data, method = c("rank", "moment"),
                           na_method = c("error", "pairwise", "complete"),
                           ci = FALSE, p_value = FALSE, conf_level = 0.95,
                           n_boot = 999L, n_perm = 999L, seed = NULL,
                           n_threads = getOption("matrixCorr.threads", 1L),
                           output = c("matrix", "sparse", "edge_list"),
                           threshold = 0, diag = TRUE,
                           ci_method = c("conservative", "regular")) {
  method <- match.arg(method)
  na_method <- match.arg(na_method)
  ci_method <- match.arg(ci_method)
  check_bool(ci, arg = "ci")
  check_bool(p_value, arg = "p_value")
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")
  if (ci) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
    n_boot <- check_scalar_int_pos(n_boot, arg = "n_boot")
    if (n_boot < 2L) abort_bad_arg("n_boot", message = "must be at least two.")
  }
  if (p_value) n_perm <- check_scalar_int_pos(n_perm, arg = "n_perm")
  cfg <- .mc_prepare_corr_output(output, threshold, diag)
  input <- .mc_prepare_corr_input(data, .mc_resolve_corr_na(na_method))
  if (ci || p_value) {
    seed <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else
      check_scalar_int_pos(seed, arg = "seed")
  }

  result <- lancaster_matrix_cpp(
    input$data, method == "rank", ci, p_value, if (ci) conf_level else 0.95,
    if (ci) n_boot else 0L, if (p_value) n_perm else 0L,
    if (ci || p_value) seed else 0L, n_threads, ci_method == "conservative"
  )
  result <- lapply(result, function(m) .mc_set_matrix_dimnames(m, input$colnames))
  diagnostics <- list(n_complete = result$n_complete, linear = result$first,
                      quadratic = result$second)
  attrs <- list()
  if (ci) {
    diagnostics$n_boot_success <- result$n_boot_success
    attrs$ci <- list(est = result$est, lwr.ci = result$lwr, upr.ci = result$upr,
                     conf.level = conf_level, ci.method = paste0("covariance_bootstrap_", ci_method))
    attrs$conf.level <- conf_level
    attrs$n_boot <- n_boot
  }
  if (p_value) {
    attrs$inference <- list(method = "permutation", statistic = result$est,
                            p_value = result$p_value, n_perm = n_perm,
                            alternative = "greater")
  }
  .mc_finalize_corr_result(
    result$est, "lancaster_corr", paste0("lancaster_", method),
    paste("Lancaster correlation,", method, "estimator"), cfg,
    diagnostics = .mc_merge_diagnostics(input$diagnostics, diagnostics),
    dimnames = input$dimnames, symmetric = TRUE, extra_attrs = attrs
  )
}

#' @method print lancaster_corr
#' @export
print.lancaster_corr <- function(x, ...) {
  .mc_print_corr_matrix(x, header = "Lancaster correlation matrix", ...)
}
