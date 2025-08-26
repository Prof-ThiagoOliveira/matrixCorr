#' @title Repeated-Measures Lin's Concordance Correlation Coefficient (CCC)
#'
#' @description
#' Computes all pairwise Lin's Concordance Correlation Coefficients (CCC)
#' across multiple methods (L \eqn{\geq} 2) for repeated-measures data.
#' Each subject must be measured by all methods across the same set of time
#' points or replicates.
#'
#' CCC measures both accuracy (how close measurements are to the line of
#' equality) and precision (Pearson correlation). Confidence intervals are
#' optionally computed using a U-statistics-based estimator with Fisher's Z
#' transformation
#'
#' @param data A data frame containing the repeated-measures dataset.
#' @param ry Character. Name of the numeric outcome column.
#' @param rmet Character. Name of the method column (factor with L
#' \eqn{\geq} 2 levels).
#' @param rtime Character or NULL. Name of the time/repetition column. If NULL,
#' one time point is assumed.
#' @param Dmat Optional numeric weight matrix (T \eqn{\times} T) for
#' timepoints. Defaults to identity.
#' @param delta Numeric. Power exponent used in the distance computations
#' between method trajectories
#' across time points. This controls the contribution of differences between
#' measurements:
#' \itemize{
#'   \item \code{delta = 1} (default) uses **absolute differences**.
#'   \item \code{delta = 2} uses **squared differences**, more sensitive to
#'   larger deviations.
#'   \item \code{delta = 0} reduces to a **binary distance** (presence/absence
#'   of disagreement), analogous to a repeated-measures version of the kappa
#'   statistic.
#' }
#' The choice of \code{delta} should reflect the penalty you want to assign to
#' measurement disagreement.
#
#' @param ci Logical. If TRUE, returns confidence intervals (default FALSE).
#' @param conf_level Confidence level for CI (default 0.95).
#' @param n_threads Integer (\eqn{\geq} 1). Number of OpenMP threads to use for computation.
#'   Defaults to \code{getOption("matrixCorr.threads", 1L)}.
#' @param verbose Logical. If TRUE, prints diagnostic output (default FALSE).
#'
#' @return
#' If \code{ci = FALSE}, a symmetric matrix of class \code{"ccc"} (estimates only).
#' If \code{ci = TRUE}, a list of class \code{"ccc"}, \code{"ccc_ci"} with elements:
#' \itemize{
#'   \item \code{est}: CCC estimate matrix
#'   \item \code{lwr.ci}: Lower bound matrix
#'   \item \code{upr.ci}: Upper bound matrix
#' }
#'
#' @details
#' This function computes pairwise Lin's Concordance Correlation Coefficient
#' (CCC) between methods in a repeated-measures design using a
#' U-statistics-based nonparametric estimator proposed by
#' Carrasco et al. (2013). It is computationally efficient and robust,
#' particularly for large-scale or balanced longitudinal designs.
#'
#' Lin's CCC is defined as
#' \deqn{
#' \rho_c = \frac{2 \cdot \mathrm{cov}(X, Y)}{\sigma_X^2 + \sigma_Y^2 +
#' (\mu_X - \mu_Y)^2}
#' }{
#' CCC = 2 * cov(X, Y) / [var(X) + var(Y) + (mean(X) - mean(Y))^2]
#' }
#' where:
#' \itemize{
#'   \item \eqn{X} and \eqn{Y} are paired measurements from two methods.
#'   \item \eqn{\mu_X}, \eqn{\mu_Y} are means, and \eqn{\sigma_X^2},
#'   \eqn{\sigma_Y^2} are variances.
#' }
#'
#' ## U-statistics Estimation
#' For repeated measures across \eqn{T} time points and \eqn{n} subjects we
#' assume
#' \itemize{
#'   \item all \eqn{n(n-1)} pairs of subjects are considered to compute a
#'   U-statistic estimator for within-method and cross-method distances.
#'   \item if `delta > 0`, pairwise distances are raised to a power before
#'   applying a time-weighted kernel matrix \eqn{D}.
#'   \item if `delta = 0`, the method reduces to a version similar to a
#'   repeated-measures kappa.
#' }
#'
#' ## Confidence Intervals
#' Confidence intervals are constructed using a **Fisher Z-transformation**
#' of the CCC. Specifically,
#' \itemize{
#'   \item The CCC is transformed using
#'   \eqn{Z = 0.5 \log((1 + \rho_c) / (1 - \rho_c))}.
#'   \item Standard errors are computed from the asymptotic variance of the
#'   U-statistic.
#'   \item Normal-based intervals are computed on the Z-scale and then
#'   back-transformed to the CCC scale.
#' }
#'
#' ## Assumptions
#' \itemize{
#'   \item The design must be **balanced**, where all subjects must have
#'   complete observations for all methods and time points.
#'   \item The method is **nonparametric** and does not require assumptions
#'   of normality or linear mixed effects.
#'   \item Weights (`Dmat`) allow differential importance of time points.
#' }
#'
#' For **unbalanced** or **complex hierarchical** data (e.g.,
#' missing timepoints, covariate adjustments), consider using
#' \code{\link{ccc_lmm_reml}}, which uses a variance components approach
#' via linear mixed models.
#'
#' @references
#' Lin L (1989). A concordance correlation coefficient to evaluate reproducibility.
#' \emph{Biometrics}, 45: 255-268.
#'
#' Lin L (2000). A note on the concordance correlation coefficient.
#' \emph{Biometrics}, 56: 324-325.
#'
#' Carrasco JL, Jover L (2003). Estimating the concordance correlation coefficient:
#' a new approach. \emph{Computational Statistics & Data Analysis}, 47(4): 519-539.
#'
#' @seealso \code{\link{ccc}}, \code{\link{ccc_lmm_reml}},
#' \code{\link{plot.ccc}}, \code{\link{print.ccc}}
#'
#' @examples
#' set.seed(123)
#' df <- expand.grid(subject = 1:10,
#'                   time = 1:2,
#'                   method = c("A", "B", "C"))
#' df$y <- rnorm(nrow(df), mean = match(df$method, c("A", "B", "C")), sd = 1)
#'
#' # CCC matrix (no CIs)
#' ccc1 <- ccc_pairwise_u_stat(df, ry = "y", rmet = "method", rtime = "time")
#' print(ccc1)
#' summary(ccc1)
#' plot(ccc1)
#'
#' # With confidence intervals
#' ccc2 <- ccc_pairwise_u_stat(df, ry = "y", rmet = "method", rtime = "time", ci = TRUE)
#' print(ccc2)
#' summary(ccc2)
#' plot(ccc2)
#'
#' #------------------------------------------------------------------------
#' # Choosing delta based on distance sensitivity
#' #------------------------------------------------------------------------
#' # Absolute distance (L1 norm) - robust
#' ccc_pairwise_u_stat(df, ry = "y", rmet = "method", rtime = "time", delta = 1)
#'
#' # Squared distance (L2 norm) - amplifies large deviations
#' ccc_pairwise_u_stat(df, ry = "y", rmet = "method", rtime = "time", delta = 2)
#'
#' # Presence/absence of disagreement (like kappa)
#' ccc_pairwise_u_stat(df, ry = "y", rmet = "method", rtime = "time", delta = 0)
#'
#' @author Thiago de Paula Oliveira
#' @export
ccc_pairwise_u_stat <- function(data,
                        ry,
                        rmet,
                        rtime = NULL,
                        Dmat = NULL,
                        delta = 1,
                        ci = FALSE,
                        conf_level = 0.95,
                        n_threads = getOption("matrixCorr.threads", 1L),
                        verbose = FALSE) {
  df <- as.data.frame(data)
  df[[rmet]] <- factor(df[[rmet]])
  method_levels <- levels(df[[rmet]])
  L <- length(method_levels)

  if (L < 2) stop("Need at least two methods (levels in rmet)")

  if (is.null(rtime)) {
    df$time_i <- 0
    ntime <- 1
  } else {
    df[[rtime]] <- factor(df[[rtime]])
    df$time_i <- as.integer(df[[rtime]]) - 1
    ntime <- length(levels(df[[rtime]]))
  }

  if (is.null(Dmat)) {
    Dmat <- diag(ntime)
  } else {
    Dmat <- as.matrix(Dmat)
    if (!all(dim(Dmat) == c(ntime, ntime))) stop("Dmat dimension mismatch")
  }

  cccr_est <- matrix(1, L, L, dimnames = list(method_levels, method_levels))
  cccr_lwr <- matrix(NA_real_, L, L, dimnames = list(method_levels, method_levels))
  cccr_upr <- matrix(NA_real_, L, L, dimnames = list(method_levels, method_levels))

  n_threads <- max(1L, as.integer(n_threads))
  set_omp_threads(n_threads)
  if (verbose) cat("Using", get_omp_threads(), "OpenMP threads\n")

  # Loop over all method pairs
  for (i in 1:(L - 1)) {
    for (j in (i + 1):L) {
      m1 <- method_levels[i]
      m2 <- method_levels[j]

      df_sub <- df[df[[rmet]] %in% c(m1, m2), , drop = FALSE]
      df_sub$met_i <- as.integer(factor(df_sub[[rmet]], levels = c(m1, m2))) - 1

      ns <- nrow(df_sub) / (2 * ntime)
      if (ns != floor(ns))
        stop(sprintf("Data inconsistent for method pair %s vs %s", m1, m2))

      res <- cccUst_rcpp(df_sub[[ry]],
                         df_sub$met_i,
                         df_sub$time_i,
                         0, 1,
                         ntime, ns,
                         Dmat, delta,
                         conf_level)

      cccr_est[i, j] <- cccr_est[j, i] <- res$CCC
      cccr_lwr[i, j] <- cccr_lwr[j, i] <- res$LL_CI
      cccr_upr[i, j] <- cccr_upr[j, i] <- res$UL_CI
    }
  }

  if (ci) {
    result <- list(est = cccr_est, lwr.ci = cccr_lwr, upr.ci = cccr_upr)
    attr(result, "method") <- "Repeated-measures Lin's concordance"
    attr(result, "description") <- "Repeated-measures CCC (pairwise) with confidence intervals"
    attr(result, "package") <- "matrixCorr"
    attr(result, "conf.level")  <- conf_level
    class(result) <- c("ccc", "ccc_ci")
  } else {
    result <- cccr_est
    attr(result, "method") <- "Repeated-measures Lin's concordance"
    attr(result, "description") <- "Repeated-measures CCC (pairwise matrix)"
    attr(result, "package") <- "matrixCorr"
    class(result) <- c("ccc", "matrix")
  }

  return(result)
}

#' @title Concordance Correlation via REML (Linear Mixed-Effects Model)
#'
#' @description
#' Compute Lin's Concordance Correlation Coefficient (CCC) from a linear
#' mixed-effects model fitted by REML. The fixed-effects part can include
#' \code{method} and/or \code{time} factors (and optionally their interaction),
#' while a subject-specific random intercept captures between-subject variation.
#' The implementation avoids any \eqn{n \times n}
#' inversions by working with small per-subject systems via the Woodbury
#' identity.
#'
#' \strong{Assumption:} time levels are treated as \emph{regular, equally spaced}
#' visits indexed by their order within subject. The AR(1) residual model is discrete-time
#' on visit index (not calendar time). \code{NA} time codes break the serial run. Gaps in
#' the factor levels are \emph{ignored} (adjacent observed visits are treated as lag-1).
#'
#' @param data A data frame.
#' @param ry Character. Response variable name.
#' @param rind Character. Subject ID variable name (random intercept).
#' @param rmet Character or \code{NULL}. Optional column name of method factor
#'   (added to fixed effects).
#' @param rtime Character or \code{NULL}. Optional column name of time factor
#'   (added to fixed effects).
#' @param interaction Logical. Include \code{method:time} interaction?
#'   (default \code{FALSE}).
#' @param max_iter Integer. Maximum iterations for variance-component updates
#'   (default \code{100}).
#' @param tol Numeric. Convergence tolerance on parameter change
#'   (default \code{1e-6}).
#'
#' @param Dmat Optional \eqn{n_t \times n_t} numeric matrix to weight/aggregate
#'   time-specific fixed biases in the \eqn{S_B} quadratic form. If supplied, it
#'   is used (after optional mass rescaling; see \code{Dmat_rescale}) whenever at
#'   least two \emph{present} time levels exist; otherwise it is ignored. \strong{If
#'   \code{Dmat} is \code{NULL}}, a canonical kernel \eqn{D_m} is \emph{constructed}
#'   from \code{Dmat_type} and \code{Dmat_weights} (see below). \code{Dmat} should be
#'   symmetric positive semidefinite; small asymmetries are symmetrized internally.
#'
#' @param Dmat_type Character, one of \code{c("time-avg","typical-visit",
#'   "weighted-avg","weighted-sq")}. Only used when \code{Dmat = NULL}.
#'   It selects the aggregation target for time-specific fixed biases in
#'   \eqn{S_B}. Options are:
#'
#'   \itemize{
#'     \item \code{"time-avg"}: square of the time-averaged bias, \eqn{D_m=(1/n_t)\,11^\top}.
#'     \item \code{"typical-visit"}: average of squared per-time biases, \eqn{D_m=I_{n_t}}.
#'     \item \code{"weighted-avg"}: square of a weighted average, \eqn{D_m=n_t\,w\,w^\top} with \eqn{\sum w=1}.
#'     \item \code{"weighted-sq"}: weighted average of squared biases, \eqn{D_m=n_t\,\mathrm{diag}(w)} with \eqn{\sum w=1}.
#'   }
#'   Pick \code{"time-avg"} for CCC targeting the time-averaged measurement; pick
#'   \code{"typical-visit"} for CCC targeting a randomly sampled visit (typical occasion).
#'   Default \code{"time-avg"}.
#'
#' @param Dmat_weights Optional numeric weights \eqn{w} used when
#'   \code{Dmat_type \%in\% c("weighted-avg","weighted-sq")}. Must be nonnegative and
#'   finite. If \code{names(w)} are provided, they should match the \emph{full} time
#'   levels in \code{data}; they are aligned to the \emph{present} time subset per fit.
#'   If unnamed, the length must equal the number of present time levels. In all cases
#'   \eqn{w} is internally normalized to sum to 1.
#'
#' @param Dmat_rescale Logical. When \code{TRUE} (default), the supplied/built
#'   \eqn{D_m} is rescaled to satisfy the simple mass rule
#'   \eqn{1^\top D_m 1 = n_t}. This keeps the \eqn{S_B} denominator invariant and
#'   harmonizes with the \eqn{\kappa}-shrinkage used for variance terms.
#'
#' @param ci Logical. If \code{TRUE}, return a CI container; limits are computed
#'   by a large-sample delta method for CCC (see \strong{CIs} note below).
#' @param conf_level Numeric in \eqn{(0,1)}. Confidence level when
#'   \code{ci = TRUE} (default \code{0.95}).
#' @param verbose Logical. If \code{TRUE}, prints a structured summary of the
#'   fitted variance components and \eqn{S_B} for each fit (overall or
#'   pairwise). Default \code{FALSE}.
#' @param digits Integer \eqn{(\ge 0)}. Number of decimal places to use in the
#'   printed summary when \code{verbose = TRUE}. Default \code{4}.
#' @param use_message Logical. When \code{verbose = TRUE}, choose the printing
#'   mechanism, where \code{TRUE} uses \code{message()} (respects \code{sink()},
#'   easily suppressible via \code{suppressMessages()}), whereas \code{FALSE}
#'   uses \code{cat()} to \code{stdout}. Default \code{TRUE}.
#'
#' @param ar Character. Residual correlation structure: \code{"none"} (iid) or
#'   \code{"ar1"} for subject-level AR(1) correlation within contiguous time
#'   runs. Default \code{c("none","ar1")}.
#' @param ar_rho Numeric of length 1 in \eqn{(-0.999,\,0.999)} or \code{NA}.
#'   When \code{ar = "ar1"} and \code{ar_rho} is finite, it is treated as fixed.
#'   When \code{ar = "ar1"} and \code{ar_rho = NA}, \emph{\eqn{\rho} is estimated}
#'   by profiling a 1-D objective (preferring REML if available; a proxy objective
#'   is used as a fallback). Default \code{NA_real_}.
#'
#' @param slope Character. Optional extra random-effect design \eqn{Z}.
#'   With \code{"subject"} a single random slope is added (one column in \eqn{Z});
#'   with \code{"method"} one column per method level is added; with
#'   \code{"custom"} you provide \code{slope_Z} directly. Default
#'   \code{c("none","subject","method","custom")}.
#' @param slope_var For \code{slope \%in\% c("subject","method")}, a character
#'   string giving the name of a column in \code{data} used as the slope regressor
#'   (e.g., centered time). It is looked up inside \code{data}; do not
#'   pass the vector itself. NAs are treated as zeros in \eqn{Z}.
#' @param slope_Z For \code{slope = "custom"}, a numeric matrix with \eqn{n}
#'   rows (same order as \code{data}) providing the full extra random-effect
#'   design \eqn{Z}. \strong{Note:} all columns of \code{slope_Z} share a
#'   single pooled variance component \eqn{\sigma_Z^2}; per-column variances
#'   are not yet supported (under development). Ignored otherwise.
#' @param drop_zero_cols Logical. When \code{slope = "method"}, drop all-zero
#'   columns of \eqn{Z} after subsetting (useful in pairwise fits). Default
#'   \code{TRUE}.
#'
#' @details
#' For measurement \eqn{y_{ij}} on subject \eqn{i} under fixed
#' levels (method, time), we fit
#' \deqn{ y = X\beta + Zu + \varepsilon,\qquad
#'        u \sim N(0,\,G),\ \varepsilon \sim N(0,\,R). }
#' Notation: \eqn{m} subjects, \eqn{n=\sum_i n_i} total rows;
#' \eqn{nm} method levels; \eqn{nt} time levels; \eqn{q_Z} extra
#' random-slope columns (if any); \eqn{r=1+nm+nt} (or \eqn{1+nm+nt+q_Z} with slopes).
#' Here \eqn{Z} is the subject-structured random-effects design and \eqn{G} is
#' block-diagonal at the subject level with the following \emph{per-subject}
#' parameterisation. Specifically,
#' \itemize{
#'   \item one random intercept with variance \eqn{\sigma_A^2};
#'   \item optionally, \emph{method} deviations (one column per method level)
#'         with a common variance \eqn{\sigma_{A\times M}^2} and zero
#'         covariances across levels (i.e., multiple of an identity);
#'   \item optionally, \emph{time} deviations (one column per time level)
#'         with a common variance \eqn{\sigma_{A\times T}^2} and zero
#'         covariances across levels;
#'   \item optionally, an \emph{extra} random effect aligned with \eqn{Z}
#'         (random slope), where variance \eqn{\sigma_Z^2} times an identity on the
#'         \eqn{Z}-columns (see \strong{Random-slope \eqn{Z}}).
#' }
#' The fixed-effects design is \code{~ 1 + rmet + rtime} and, if
#' \code{interaction=TRUE}, \code{+ rmet:rtime}.
#'
#' \strong{Residual correlation \eqn{R} (regular, equally spaced time).}
#' Write \eqn{R_i=\sigma_E^2\,C_i(\rho)}. With \code{ar="none"}, \eqn{C_i=I}.
#' With \code{ar="ar1"}, within-subject residuals follow a \emph{discrete} AR(1)
#' process along the visit index after sorting by increasing time level. Ties
#' retain input order, and any \code{NA} time code breaks the series so each
#' contiguous block of non-\code{NA} times forms a run. The correlation
#' between \emph{adjacent observed visits} in a run is \eqn{\rho}; we do not use
#' calendar-time gaps. Internally we work with the \emph{precision} of the AR(1)
#' correlation: for a run of length \eqn{L\ge 2}, the tridiagonal inverse has
#' \deqn{ (C^{-1})_{11}=(C^{-1})_{LL}=\frac{1}{1-\rho^2},\quad
#'        (C^{-1})_{tt}=\frac{1+\rho^2}{1-\rho^2}\ (2\le t\le L-1),\quad
#'        (C^{-1})_{t,t+1}=(C^{-1})_{t+1,t}=\frac{-\rho}{1-\rho^2}. }
#' The working inverse is \eqn{\,R_i^{-1}=\sigma_E^{-2}\,C_i(\rho)^{-1}}.
#'
#' \strong{Per-subject Woodbury system.} For subject \eqn{i} with \eqn{n_i}
#' rows, define the per-subject random-effects design \eqn{U_i} (columns:
#' intercept, method indicators, time indicators; dimension
#' \eqn{\,r=1+nm+nt\,}). The core never forms
#' \eqn{V_i = R_i + U_i G U_i^\top} explicitly. Instead,
#' \deqn{ M_i \;=\; G^{-1} \;+\; U_i^\top R_i^{-1} U_i, }
#' and accumulates GLS blocks via rank-\eqn{r} corrections using
#' \eqn{\,V_i^{-1} = R_i^{-1}-R_i^{-1}U_i M_i^{-1}U_i^\top R_i^{-1}\,}:
#' \deqn{ X^\top V^{-1} X \;=\; \sum_i \Big[
#'        X_i^\top R_i^{-1}X_i \;-\; (X_i^\top R_i^{-1}U_i)\,
#'        M_i^{-1}\,(U_i^\top R_i^{-1}X_i) \Big], }
#' \deqn{ X^\top V^{-1} y \;=\; \sum_i \Big[
#'        X_i^\top R_i^{-1}y_i \;-\; (X_i^\top R_i^{-1}U_i)\,M_i^{-1}\,
#'        (U_i^\top R_i^{-1}y_i) \Big]. }
#' Because \eqn{G^{-1}} is diagonal with positive entries, each \eqn{M_i} is
#' symmetric positive definite; solves/inversions use symmetric-PD routines with
#' a tiny diagonal "jitter" and a pseudo-inverse fallback when needed.
#'
#' \strong{Random-slope \eqn{Z}.}
#' Besides \eqn{U_i}, the function can include an extra design \eqn{Z_i} and a
#' corresponding variance \eqn{\sigma_Z^2}:
#' \itemize{
#'   \item \code{slope="subject"}: \eqn{Z} has one column (the regressor in
#'         \code{slope_var}); \eqn{Z_{i}} is the subject-\eqn{i} block.
#'   \item \code{slope="method"}: \eqn{Z} has one column per method level;
#'         row \eqn{t} uses the slope regressor if its method equals level \eqn{\ell},
#'         otherwise 0; all-zero columns can be dropped via
#'         \code{drop_zero_cols=TRUE} after subsetting.
#'   \item \code{slope="custom"}: \eqn{Z} is provided fully via \code{slope_Z}.
#' }
#' Computations simply augment \eqn{\tilde U_i=[U_i\ Z_i]} and the corresponding
#' inverse-variance block. The EM updates then include
#' \deqn{ \sigma_Z^{2\,(new)} \;=\; \frac{1}{m\,q_Z}
#'       \sum_i \sum_{j=1}^{q_Z}\!\Big( b_{i,\text{extra},j}^2 +
#'       (M_i^{-1})_{\text{extra},jj} \Big)
#'       \quad (\text{if } q_Z>0). }
#' \emph{Interpretation:} \eqn{\sigma_Z^2} represents additional within-subject
#' variability explained by the slope regressor(s) and is \emph{not} part of the CCC
#' denominator (agreement across methods/time).
#'
#' \strong{EM-style variance-component updates.} With current \eqn{\hat\beta},
#' form residuals \eqn{r_i = y_i - X_i\hat\beta}. The BLUPs and conditional
#' covariances are
#' \deqn{ b_i \;=\; M_i^{-1}\,(U_i^\top R_i^{-1} r_i), \qquad
#'       \mathrm{Var}(b_i\mid y) \;=\; M_i^{-1}. }
#' Let \eqn{e_i=r_i-U_i b_i}. Expected squares then yield closed-form updates:
#' \deqn{ \sigma_A^{2\,(new)} \;=\; \frac{1}{m}\sum_i \Big( b_{i,0}^2 +
#' (M_i^{-1})_{00} \Big), }
#' \deqn{ \sigma_{A\times M}^{2\,(new)} \;=\; \frac{1}{m\,nm}
#'       \sum_i \sum_{\ell=1}^{nm}\!\Big( b_{i,\ell}^2 +
#'       (M_i^{-1})_{\ell\ell} \Big)
#'       \quad (\text{if } nm>0), }
#' \deqn{ \sigma_{A\times T}^{2\,(new)} \;=\; \frac{1}{m\,nt}
#'       \sum_i \sum_{t=1}^{nt}\!\Big( b_{i,t}^2 + (M_i^{-1})_{tt} \Big)
#'       \quad (\text{if } nt>0), }
#' \deqn{ \sigma_E^{2\,(new)} \;=\; \frac{1}{n} \sum_i
#'       \Big( e_i^\top C_i(\rho)^{-1} e_i \;+\;
#'       \mathrm{tr}\!\big(M_i^{-1}U_i^\top C_i(\rho)^{-1} U_i\big) \Big), }
#' where \eqn{C_i(\rho)^{-1}} is the AR(1) precision built on the
#' time-ordered runs described above (and equals \eqn{I} for i.i.d. residuals).
#' Iterate until the \eqn{\ell_1} change across components is \eqn{<}
#' \code{tol} or \code{max_iter} is reached.
#'
#' \strong{Fixed-effect dispersion \eqn{S_B}: choosing the time-kernel \eqn{D_m}.}
#'
#' Let \eqn{d = L^\top \hat\beta} stack the within-time, pairwise method differences,
#' grouped by time as \eqn{d=(d_1^\top,\ldots,d_{n_t}^\top)^\top} with
#' \eqn{d_t \in \mathbb{R}^{P}} and \eqn{P = n_m(n_m-1)}. The symmetric
#' positive semidefinite kernel \eqn{D_m \succeq 0} selects which functional of the
#' bias profile \eqn{t \mapsto d_t} is targeted by \eqn{S_B}. Internally, the code
#' rescales any supplied/built \eqn{D_m} to satisfy \eqn{1^\top D_m 1 = n_t} for
#' stability and comparability.
#'
#' \itemize{
#'   \item \code{Dmat_type = "time-avg"} (square of the time-averaged bias).
#'     Let \deqn{ w \;=\; \frac{1}{n_t}\,\mathbf{1}_{n_t}, \qquad
#'                D_m \;\propto\; I_P \otimes (w\,w^\top), }
#'     so that
#'     \deqn{ d^\top D_m d \;\propto\; \sum_{p=1}^{P}
#'            \left( \frac{1}{n_t}\sum_{t=1}^{n_t} d_{t,p} \right)^{\!2}. }
#'     Methods have equal \eqn{\textit{time-averaged}}
#'     means within subject, i.e. \eqn{\sum_{t=1}^{n_t} d_{t,p}/n_t = 0} for all
#'     \eqn{p}. Appropriate when decisions depend on an average over time and
#'     opposite-signed biases are allowed to cancel.
#'
#'   \item \code{Dmat_type = "typical-visit"} (average of squared per-time biases).
#'     With equal visit probability, take
#'     \deqn{ D_m \;\propto\; I_P \otimes
#'            \mathrm{diag}\!\Big(\tfrac{1}{n_t},\ldots,\tfrac{1}{n_t}\Big), }
#'     yielding
#'     \deqn{ d^\top D_m d \;\propto\; \frac{1}{n_t}
#'            \sum_{t=1}^{n_t}\sum_{p=1}^{P} d_{t,p}^{\,2}. }
#'     Methods agree on a \eqn{\textit{typical}}
#'     occasion drawn uniformly from the visit set. Use when each visit matters
#'     on its own; alternating signs \eqn{d_{t,p}} do not cancel.
#'
#'   \item \code{Dmat_type = "weighted-avg"} (square of a weighted time average).
#'     For user weights \eqn{a=(a_1,\ldots,a_{n_t})^\top} with \eqn{a_t \ge 0}, set
#'     \deqn{ w \;=\; \frac{a}{\sum_{t=1}^{n_t} a_t}, \qquad
#'            D_m \;\propto\; I_P \otimes (w\,w^\top), }
#'     so that
#'     \deqn{ d^\top D_m d \;\propto\; \sum_{p=1}^{P}
#'            \left( \sum_{t=1}^{n_t} w_t\, d_{t,p} \right)^{\!2}. }
#'     Methods have equal \eqn{\textit{weighted
#'     time-averaged}} means, i.e. \eqn{\sum_{t=1}^{n_t} w_t\, d_{t,p} = 0} for all
#'     \eqn{p}. Use when some visits (e.g., baseline/harvest) are a priori more
#'     influential; opposite-signed biases may cancel according to \eqn{w}.
#'
#'   \item \code{Dmat_type = "weighted-sq"} (weighted average of squared per-time biases).
#'     With the same weights \eqn{w}, take
#'     \deqn{ D_m \;\propto\; I_P \otimes \mathrm{diag}(w_1,\ldots,w_{n_t}), }
#'     giving
#'     \deqn{ d^\top D_m d \;\propto\; \sum_{t=1}^{n_t} w_t
#'            \sum_{p=1}^{P} d_{t,p}^{\,2}. }
#'     Methods agree at visits sampled with
#'     probabilities \eqn{\{w_t\}}, counting each visit’s discrepancy on its own.
#'     Use when per-visit agreement is required but some visits should be
#'     emphasised more than others.
#' }
#'
#' \strong{Time-averaging for CCC (regular visits).}
#' The reported CCC targets agreement of the \emph{time-averaged} measurements
#' per method within subject by default (\code{Dmat_type="time-avg"}). Averaging over \eqn{T}
#' non-\code{NA} visits shrinks time-varying components by
#' \deqn{ \kappa_T^{(g)} \;=\; 1/T, \qquad
#'       \kappa_T^{(e)} \;=\; \{T + 2\sum_{k=1}^{T-1}(T-k)\rho^k\}/T^2, }
#' with \eqn{\kappa_T^{(e)}=1/T} when residuals are i.i.d. With unbalanced \eqn{T}, the
#' implementation averages the per-(subject,method) \eqn{\kappa} values across the
#' pairs contributing to CCC and then clamps \eqn{\kappa_T^{(e)}} to
#' \eqn{[10^{-12},\,1]} for numerical stability. Choosing
#' \code{Dmat_type="typical-visit"} makes \eqn{S_B} match the interpretation of a
#' randomly sampled occasion instead.
#'
#' \strong{Concordance correlation coefficient.} The CCC used is
#' \deqn{ \mathrm{CCC} \;=\;
#'       \frac{\sigma_A^2 + \kappa_T^{(g)}\,\sigma_{A\times T}^2}
#'            {\sigma_A^2 + \sigma_{A\times M}^2 +
#'             \kappa_T^{(g)}\,\sigma_{A\times T}^2 + S_B +
#'             \kappa_T^{(e)}\,\sigma_E^2}. }
#' Special cases: with no method factor, \eqn{S_B=\sigma_{A\times M}^2=0}; with
#' a single time level, \eqn{\sigma_{A\times T}^2=0} (no \eqn{\kappa}-shrinkage).
#' When \eqn{T=1} or \eqn{\rho=0}, both \eqn{\kappa}-factors equal 1. The extra
#' random-effect variance \eqn{\sigma_Z^2} (if used) is not included.
#'
#' \strong{CIs / SEs (delta method for CCC).}
#' Let
#' \deqn{ \theta \;=\; \big(\sigma_A^2,\ \sigma_{A\times M}^2,\
#' \sigma_{A\times T}^2,\ \sigma_E^2,\ S_B\big)^\top, }
#' and write \eqn{\mathrm{CCC}(\theta)=N/D} with
#' \eqn{N=\sigma_A^2+\kappa_T^{(g)}\sigma_{A\times T}^2} and
#' \eqn{D=\sigma_A^2+\sigma_{A\times M}^2+\kappa_T^{(g)}\sigma_{A\times T}^2+S_B+\kappa_T^{(e)}\sigma_E^2}.
#' The gradient components are
#' \deqn{ \frac{\partial\,\mathrm{CCC}}{\partial \sigma_A^2}
#'       \;=\; \frac{\sigma_{A\times M}^2 + S_B + \kappa_T^{(e)}\sigma_E^2}{D^2}, }
#' \deqn{ \frac{\partial\,\mathrm{CCC}}{\partial \sigma_{A\times M}^2}
#'       \;=\; -\,\frac{N}{D^2}, \qquad
#'        \frac{\partial\,\mathrm{CCC}}{\partial \sigma_{A\times T}^2}
#'       \;=\; \frac{\kappa_T^{(g)}\big(\sigma_{A\times M}^2 + S_B +
#'                                     \kappa_T^{(e)}\sigma_E^2\big)}{D^2}, }
#' \deqn{ \frac{\partial\,\mathrm{CCC}}{\partial \sigma_E^2}
#'       \;=\; -\,\frac{\kappa_T^{(e)}\,N}{D^2}, \qquad
#'        \frac{\partial\,\mathrm{CCC}}{\partial S_B}
#'       \;=\; -\,\frac{N}{D^2}. }
#'
#' \emph{Estimating \eqn{\mathrm{Var}(\hat\theta)}.}
#' The EM updates write each variance component as an average of per-subject
#' quantities. For subject \eqn{i},
#' \deqn{ t_{A,i} \;=\; b_{i,0}^2 + (M_i^{-1})_{00},\qquad
#'        t_{M,i} \;=\; \frac{1}{nm}\sum_{\ell=1}^{nm}
#'                        \Big(b_{i,\ell}^2 + (M_i^{-1})_{\ell\ell}\Big), }
#' \deqn{ t_{T,i} \;=\; \frac{1}{nt}\sum_{j=1}^{nt}
#'                        \Big(b_{i,j}^2 + (M_i^{-1})_{jj}\Big),\qquad
#'        s_i \;=\; \frac{e_i^\top C_i(\rho)^{-1} e_i +
#'        \mathrm{tr}\!\big(M_i^{-1}U_i^\top C_i(\rho)^{-1} U_i\big)}{n_i}, }
#' where \eqn{b_i = M_i^{-1}(U_i^\top R_i^{-1} r_i)} and
#' \eqn{M_i = G^{-1} + U_i^\top R_i^{-1} U_i}.
#' With \eqn{m} subjects, form the empirical covariance of the stacked
#' subject vectors and scale by \eqn{m} to approximate the covariance of the
#' means:
#' \deqn{ \widehat{\mathrm{Cov}}\!\left(
#'       \begin{bmatrix} t_{A,\cdot} \\ t_{M,\cdot} \\ t_{T,\cdot} \end{bmatrix}
#'       \right)
#'       \;\approx\; \frac{1}{m}\,
#'        \mathrm{Cov}_i\!\left(
#'       \begin{bmatrix} t_{A,i} \\ t_{M,i} \\ t_{T,i} \end{bmatrix}\right). }
#' (Drop rows/columns as needed when \code{nm==0} or \code{nt==0}.)
#'
#' The residual variance estimator is a weighted mean
#' \eqn{\hat\sigma_E^2=\sum_i w_i s_i} with \eqn{w_i=n_i/n}. Its variance is
#' approximated by the variance of a weighted mean of independent terms,
#' \deqn{ \widehat{\mathrm{Var}}(\hat\sigma_E^2)
#'       \;\approx\; \Big(\sum_i w_i^2\Big)\,\widehat{\mathrm{Var}}(s_i), }
#' where \eqn{\widehat{\mathrm{Var}}(s_i)} is the sample variance across
#' subjects. The method-dispersion term uses the quadratic-form delta already
#' computed for \eqn{S_B}:
#' \deqn{ \widehat{\mathrm{Var}}(S_B)
#'       \;=\; \frac{2\,\mathrm{tr}\!\big((A_{\!fix}\,\mathrm{Var}(\hat\beta))^2\big)
#'              \;+\; 4\,\hat\beta^\top A_{\!fix}\,\mathrm{Var}(\hat\beta)\,
#'              A_{\!fix}\,\hat\beta}
#'                    {\big[nm\,(nm-1)\,\max(nt,1)\big]^2}, }
#' with \eqn{A_{\!fix}=L\,D_m\,L^\top}.
#'
#' \emph{Putting it together.} Assemble
#' \eqn{\widehat{\mathrm{Var}}(\hat\theta)} by combining the
#' \eqn{(\sigma_A^2,\sigma_{A\times M}^2,\sigma_{A\times T}^2)} covariance
#' block from the subject-level empirical covariance, add the
#' \eqn{\widehat{\mathrm{Var}}(\hat\sigma_E^2)} and
#' \eqn{\widehat{\mathrm{Var}}(S_B)} terms on the diagonal,
#' and ignore cross-covariances across these blocks (a standard large-sample
#' simplification). Then
#' \deqn{ \widehat{\mathrm{se}}\{\widehat{\mathrm{CCC}}\}
#'       \;=\; \sqrt{\,\nabla \mathrm{CCC}(\hat\theta)^\top\,
#'                     \widehat{\mathrm{Var}}(\hat\theta)\,
#'                     \nabla \mathrm{CCC}(\hat\theta)\,}. }
#'
#' A two-sided \eqn{(1-\alpha)} normal CI is
#' \deqn{ \widehat{\mathrm{CCC}} \;\pm\; z_{1-\alpha/2}\,
#'       \widehat{\mathrm{se}}\{\widehat{\mathrm{CCC}}\}, }
#' truncated to \eqn{[0,1]} in the output for convenience. When \eqn{S_B} is
#' truncated at 0 or samples are very small/imbalanced, the normal CI can be
#' mildly anti-conservative near the boundary; a logit transform for CCC or a
#' subject-level (cluster) bootstrap can be used for sensitivity analysis.
#'
#' \strong{Choosing \eqn{\rho} for AR(1).}
#' When \code{ar="ar1"} and \code{ar_rho = NA}, \eqn{\rho} can be estimated by
#' profiling the REML log-likelihood computed at the fitted
#' \eqn{(\hat\beta,\hat G,\hat\sigma_E^2)}; when \code{ar_rho} is finite, that
#' fixed value is used. Very small numbers of visits per subject can make
#' \eqn{\rho} weakly identified; sensitivity checks over a plausible range (e.g., 0.3–0.8)
#' are recommended.
#'
#' @section Notes on stability and performance:
#' All per-subject solves are \eqn{\,r\times r} with \eqn{r=1+nm+nt+q_Z}, so cost
#' scales with the number of subjects and the fixed-effects dimension rather
#' than the total number of observations. Solvers use symmetric-PD paths with
#' a small diagonal ridge and pseudo-inverse fallback, which helps for
#' tiny/unbalanced subsets and near-boundary estimates. Very small samples or
#' extreme imbalance can still make \eqn{S_B} numerically delicate; negative
#' estimates are truncated to 0 by construction. For AR(1), observations are
#' ordered by time within subject before building the run-wise tridiagonal precision;
#' \code{NA} time codes break the run, and gaps between factor levels are treated as
#' regular steps (we do not use elapsed time).
#'
#' \emph{Heteroscedastic slopes are not supported yet.}
#' When \code{slope_Z = "custom"} the current implementation assumes a single
#' \eqn{\sigma_Z^2} for \emph{all} custom-slope columns. Distinct variances
#' (e.g., one for the global slope and another for method-specific
#' slopes) are not estimated. Column rescaling changes the implied prior on
#' \eqn{b_{i,\text{extra}}} but does not introduce separate variance components.
#'
#' @seealso \code{build_L_Dm_Z_cpp}
#' for constructing \eqn{L}/\eqn{D_m}/\eqn{Z}; \code{\link{ccc_pairwise_u_stat}}
#' for a U-statistic alternative; and \pkg{cccrm} for a reference approach via
#' \pkg{nlme}.
#'
#' @references
#' Lin L (1989). A concordance correlation coefficient to evaluate reproducibility.
#' \emph{Biometrics}, 45: 255-268.
#'
#' Lin L (2000). A note on the concordance correlation coefficient.
#' \emph{Biometrics}, 56: 324-325.
#'
#' Carrasco, J. L. et al. (2013). Estimation of the concordance
#' correlation coefficient for repeated measures using SAS and R.
#' \emph{Computer Methods and Programs in Biomedicine}, 109(3), 293–304.
#' \doi{10.1016/j.cmpb.2012.09.002}
#'
#' King et al. (2007). A Class of Repeated Measures Concordance
#' Correlation Coefficients.
#' \emph{Journal of Biopharmaceutical Statistics}, 17(4).
#' \doi{10.1080/10543400701329455}
#'
#' @examples
#' #--------------------------------------------------------------------
#' ## Two methods (no time)
#' #--------------------------------------------------------------------
#' set.seed(1)
#' n_subj <- 30
#' meth   <- factor(rep(c("A","B"), each = n_subj))
#' id     <- factor(rep(seq_len(n_subj), times = 2))
#' sigA <- 1.0; sigE <- 0.5
#' u  <- rnorm(n_subj, 0, sqrt(sigA))
#' y  <- c(u + rnorm(n_subj, 0, sqrt(sigE)),
#'          u + 0.2 + rnorm(n_subj, 0, sqrt(sigE)))
#' dat <- data.frame(y, id, method = meth)
#' ccc_rm1 <- ccc_lmm_reml(dat, ry = "y", rind = "id", rmet = "method")
#' print(ccc_rm1)
#' summary(ccc_rm1)
#'
#' # 95% CI container
#' ccc_rm2 <- ccc_lmm_reml(dat, ry = "y", rind = "id", rmet = "method", ci = TRUE)
#' ccc_rm2
#'
#' #--------------------------------------------------------------------
#' ## Two methods x time (balanced 2x2), with and without interaction
#' #--------------------------------------------------------------------
#' dat$time <- factor(rep(rep(c("t1","t2"), each = n_subj/2), times = 2))
#' ccc_lmm_reml(dat, "y", "id", rmet = "method", rtime = "time",
#'              interaction = FALSE)
#' ccc_lmm_reml(dat, "y", "id", rmet = "method", rtime = "time",
#'              interaction = TRUE, verbose = TRUE)
#'
#' #--------------------------------------------------------------------
#' ## Random slope by subject; AR(1) with rho estimated
#' #--------------------------------------------------------------------
#' dat$t_num <- as.integer(dat$time)
#' dat$t_c   <- ave(dat$t_num, dat$id, FUN = function(v) v - mean(v))
#' ccc_lmm_reml(dat, "y", "id", rmet = "method", rtime = "time",
#'              slope = "subject", slope_var = "t_c",
#'              ar = "ar1", ar_rho = NA_real_, verbose = TRUE)
#'
#' #--------------------------------------------------------------------
#' ## D_m choices: time-averaged (default) vs typical visit
#' #--------------------------------------------------------------------
#' ccc_timeavg   <- ccc_lmm_reml(dat, "y", "id", rmet = "method", rtime = "time",
#'                               Dmat_type = "time-avg")
#' ccc_typical   <- ccc_lmm_reml(dat, "y", "id", rmet = "method", rtime = "time",
#'                               Dmat_type = "typical-visit")
#'
#' # Weighted time-averaging (named weights aligned to time levels)
#' w <- c(t1 = 0.25, t2 = 0.75)
#' ccc_wavg <- ccc_lmm_reml(dat, "y", "id", rmet = "method", rtime = "time",
#'                          Dmat_type = "weighted-avg", Dmat_weights = w)
#'
#' # ------------------------------------------------------------------
#' # AR(1) residual correlation with fixed rho
#' # ------------------------------------------------------------------
#'   set.seed(10)
#'   n_subj <- 40
#'   n_time <- 6                      # ≥ 3 time points recommended for AR(1)
#'   id  <- factor(rep(seq_len(n_subj), each = n_time))
#'   tim <- factor(rep(seq_len(n_time),  times = n_subj))
#'   beta0 <- 0; beta_t <- 0.2
#'   rho_true <- 0.6; sigE <- 0.7
#'   y <- numeric(length(id))
#'   for (i in seq_len(n_subj)) {
#'     idx <- which(id == levels(id)[i])
#'     e <- stats::arima.sim(list(ar = rho_true), n = n_time, sd = sigE)
#'     y[idx] <- beta0 + beta_t*(seq_len(n_time) - mean(seq_len(n_time))) + e
#'   }
#'   dat_ar <- data.frame(y = y, id = id, time = tim)
#'   ccc_lmm_reml(dat_ar, ry = "y", rind = "id", rtime = "time",
#'                ar = "ar1", ar_rho = 0.6, verbose = TRUE)
#'
#' # ------------------------------------------------------------------
#' # Random slope by SUBJECT
#' # ------------------------------------------------------------------
#' set.seed(2)
#' n_subj <- 60; n_time <- 4
#' id  <- factor(rep(seq_len(n_subj), each = 2 * n_time))
#' tim <- factor(rep(rep(seq_len(n_time), times = 2), times = n_subj))
#' method <- factor(rep(rep(c("A","B"), each = n_time), times = n_subj))
#' subj <- as.integer(id)
#' slope_i <- rnorm(n_subj, 0, 0.15)
#' slope_vec <- slope_i[subj]
#' base <- rnorm(n_subj, 0, 1.0)[subj]
#' tnum <- as.integer(tim)
#' y <- base + 0.3*(method=="B") + slope_vec*(tnum - mean(seq_len(n_time))) +
#'      rnorm(length(id), 0, 0.5)
#' dat_s <- data.frame(y, id, method, time = tim)
#' dat_s$t_num <- as.integer(dat_s$time)
#' dat_s$t_c   <- ave(dat_s$t_num, dat_s$id, FUN = function(v) v - mean(v))
#' ccc_lmm_reml(dat_s, "y", "id", rmet = "method", rtime = "time",
#'              slope = "subject", slope_var = "t_c", verbose = TRUE)
#'
#' # ------------------------------------------------------------------
#' # Random slope by METHOD
#' # ------------------------------------------------------------------
#' set.seed(3)
#' n_subj <- 60; n_time <- 4
#' id  <- factor(rep(seq_len(n_subj), each = 2 * n_time))
#' tim <- factor(rep(rep(seq_len(n_time), times = 2), times = n_subj))
#' method <- factor(rep(rep(c("A","B"), each = n_time), times = n_subj))
#' slope_m <- ifelse(method=="B", 0.25, 0.00)
#' base <- rnorm(n_subj, 0, 1.0)[as.integer(id)]
#' tnum <- as.integer(tim)
#' y <- base + 0.3*(method=="B") + slope_m*(tnum - mean(seq_len(n_time))) +
#'      rnorm(length(id), 0, 0.5)
#' dat_m <- data.frame(y, id, method, time = tim)
#' dat_m$t_num <- as.integer(dat_m$time)
#' dat_m$t_c   <- ave(dat_m$t_num, dat_m$id, FUN = function(v) v - mean(v))
#' ccc_lmm_reml(dat_m, "y", "id", rmet = "method", rtime = "time",
#'              slope = "method", slope_var = "t_c", verbose = TRUE)
#'
#' # ------------------------------------------------------------------
#' # Random slopes for SUBJECT and METHOD (custom Z)
#' # ------------------------------------------------------------------
#' set.seed(4)
#' n_subj <- 50; n_time <- 4
#' id  <- factor(rep(seq_len(n_subj), each = 2 * n_time))
#' tim <- factor(rep(rep(seq_len(n_time), times = 2), times = n_subj))
#' method <- factor(rep(rep(c("A","B"), each = n_time), times = n_subj))
#' subj <- as.integer(id)
#' slope_subj <- rnorm(n_subj, 0, 0.12)[subj]
#' slope_B    <- ifelse(method=="B", 0.18, 0.00)
#' tnum <- as.integer(tim)
#' base <- rnorm(n_subj, 0, 1.0)[subj]
#' y <- base + 0.3*(method=="B") +
#'      (slope_subj + slope_B) * (tnum - mean(seq_len(n_time))) +
#'      rnorm(length(id), 0, 0.5)
#' dat_both <- data.frame(y, id, method, time = tim)
#' dat_both$t_num <- as.integer(dat_both$time)
#' dat_both$t_c   <- ave(dat_both$t_num, dat_both$id, FUN = function(v) v - mean(v))
#' MM <- model.matrix(~ 0 + method, data = dat_both)  # one col per method
#' # All custom-slope columns share the same variance component and are
#' # modelled as uncorrelated. Take care when using custom approach.
#' Z_custom <- cbind(
#'   subj_slope   = dat_both$t_c,
#'   MM * dat_both$t_c
#' )
#' ccc_lmm_reml(dat_both, "y", "id", rmet = "method", rtime = "time",
#'              slope = "custom", slope_Z = Z_custom, verbose = TRUE)
#'
#' @author Thiago de Paula Oliveira
#' @importFrom stats as.formula model.matrix setNames qnorm optimize
#' @export
ccc_lmm_reml <- function(data, ry, rind,
                         rmet = NULL, rtime = NULL, interaction = FALSE,
                         max_iter = 100, tol = 1e-6,
                         Dmat = NULL,
                         Dmat_type = c("time-avg","typical-visit","weighted-avg","weighted-sq"),
                         Dmat_weights = NULL,
                         Dmat_rescale = TRUE,
                         ci = FALSE, conf_level = 0.95,
                         verbose = FALSE, digits = 4, use_message = TRUE,
                         ar = c("none", "ar1"),
                         ar_rho = NA_real_,
                         slope = c("none", "subject", "method", "custom"),
                         slope_var = NULL,
                         slope_Z = NULL,
                         drop_zero_cols = TRUE) {

  ar         <- match.arg(ar)
  slope      <- match.arg(slope)
  Dmat_type  <- match.arg(Dmat_type)

  if (identical(ar, "ar1")) {
    if (length(ar_rho) != 1L) stop("ar_rho must be length 1 (or NA to estimate).")
    if (!is.na(ar_rho) && abs(ar_rho) >= 0.999) stop("ar_rho must be in (-0.999, 0.999).")
  }

  df <- as.data.frame(data)
  df[[ry]]   <- as.numeric(df[[ry]])
  df[[rind]] <- factor(df[[rind]])
  if (!is.null(rmet))  df[[rmet]]  <- factor(df[[rmet]])
  if (!is.null(rtime)) df[[rtime]] <- factor(df[[rtime]])
  all_time_lvls <- if (!is.null(rtime)) levels(df[[rtime]]) else character(0)

  rhs <- "1"
  if (!is.null(rmet))  rhs <- paste(rhs, "+", rmet)
  if (!is.null(rtime)) rhs <- paste(rhs, "+", rtime)
  if (!is.null(rmet) && !is.null(rtime) && interaction) rhs <- paste(rhs, "+", paste0(rmet, ":", rtime))
  fml <- as.formula(paste("~", rhs))

  extra_label <- switch(slope,
                        "subject" = "random slope (subject)",
                        "method"  = "random slope (by method)",
                        "custom"  = "custom random effect",
                        NULL)

  if (is.null(rmet) || nlevels(df[[rmet]]) < 2L) {
    return(ccc_lmm_reml_overall(df, fml, ry, rind, rmet, rtime,
                                slope, slope_var, slope_Z, drop_zero_cols,
                                Dmat, ar, ar_rho, max_iter, tol,
                                conf_level, verbose, digits, use_message,
                                extra_label, ci,
                                Dmat_type = Dmat_type,
                                Dmat_weights = Dmat_weights,
                                Dmat_rescale = Dmat_rescale))
  } else {
    return(ccc_lmm_reml_pairwise(df, fml, ry, rind, rmet, rtime,
                                 slope, slope_var, slope_Z, drop_zero_cols,
                                 Dmat, ar, ar_rho, max_iter, tol,
                                 conf_level, verbose, digits, use_message,
                                 extra_label, ci, all_time_lvls,
                                 Dmat_type = Dmat_type,
                                 Dmat_weights = Dmat_weights,
                                 Dmat_rescale = Dmat_rescale))
  }
}

#' @title num_or_na
#' @description Helper to safely coerce a value to numeric or return NA if invalid.
#' @keywords internal
num_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) != 1 || !is.finite(x)) NA_real_ else x
}

#' @title compute_ci_from_se
#' @description Compute confidence intervals from point estimate and standard error.
#' @keywords internal
compute_ci_from_se <- function(ccc, se, level) {
  if (!is.finite(ccc) || !is.finite(se)) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - level)/2)
  c(max(0, min(1, ccc - z * se)), max(0, min(1, ccc + z * se)))
}

#' @title .vc_message
#' @description Display variance component estimation details to the console.
#' @keywords internal
.vc_message <- function(ans, label, nm, nt, conf_level,
                        digits = 4, use_message = TRUE,
                        extra_label = NULL,
                        ar = c("none", "ar1"),
                        ar_rho = NA_real_) {
  ar <- match.arg(ar)

  fmt <- function(x) {
    if (is.null(x)) return("NA")
    x <- suppressWarnings(as.numeric(x))
    if (length(x) == 0 || all(!is.finite(x))) return("NA")
    if (length(x) == 1) {
      return(formatC(x, format = "f", digits = digits))
    }
    # vector: pretty-print as a comma-separated list
    paste0("[", paste(formatC(x, format = "f", digits = digits), collapse = ", "), "]")
  }

  colw <- 38L
  v <- function(s, x) sprintf("  %-*s : %s", colw, s, fmt(x))

  out <- c(
    sprintf("---- matrixCorr::ccc_lmm_reml - variance-components (%s) ----", label),
    sprintf("Design: methods nm = %d, times nt = %d", nm, nt)
  )
  if (identical(ar, "ar1")) {
    info <- if (is.na(ar_rho)) "AR(1) (rho estimated)" else sprintf("AR(1) with rho = %s", fmt(ar_rho))
    out <- c(out, paste("Residual correlation:", info))
  } else {
    out <- c(out, "Residual correlation: independent (iid)")
  }

  out <- c(out,
           "Estimates:",
           v("sigma_A^2 (subject)",             ans[["sigma2_subject"]]),
           v("sigma_A_M^2 (subject x method)",  ans[["sigma2_subject_method"]]),
           v("sigma_A_T^2 (subject x time)",    ans[["sigma2_subject_time"]]),
           v("sigma_E^2 (error)",               ans[["sigma2_error"]]))

  if (!is.null(ans[["sigma2_extra"]])) {
    lab <- if (is.null(extra_label)) "extra random effect" else extra_label
    out <- c(out, v(sprintf("sigma_Z^2 (%s)", lab), ans[["sigma2_extra"]]))
  }

  out <- c(out,
           v("S_B (fixed-effect dispersion)",   ans[["SB"]]),
           v("SE(CCC)",                         ans[["se_ccc"]]),
           "--------------------------------------------------------------------------")

  if (use_message) message(paste(out, collapse = "\n")) else cat(paste(out, collapse = "\n"), "\n")
}

#' @title build_LDZ
#' @description Internal helper to construct L, Dm, and Z matrices for random effects.
#' @keywords internal
build_LDZ <- function(colnames_X, method_levels, time_levels, Dsub, df_sub,
                      rmet_name, rtime_name, slope, interaction, slope_var,
                      drop_zero_cols) {
  slope_mode_cpp <- switch(slope, none = "none", subject = "subject", method = "method", custom = "none")
  if (!identical(slope, "custom")) {
    build_L_Dm_Z_cpp(
      colnames_X      = colnames_X,
      rmet_name       = if (is.null(rmet_name)) NULL else rmet_name,
      rtime_name      = if (is.null(rtime_name)) NULL else rtime_name,
      method_levels   = if (is.null(rmet_name)) character(0) else method_levels,
      time_levels     = if (is.null(rtime_name)) character(0) else time_levels,
      has_interaction = interaction,
      Dmat_global     = Dsub,
      slope_mode      = slope_mode_cpp,
      slope_var       = if (!is.null(slope_var)) df_sub[[slope_var]] else NULL,
      method_codes    = if (!is.null(rmet_name)) as.integer(df_sub[[rmet_name]]) else NULL,
      drop_zero_cols  = drop_zero_cols
    )
  } else {
    build_L_Dm_cpp(
      colnames_X      = colnames_X,
      rmet_name       = if (is.null(rmet_name)) NULL else rmet_name,
      rtime_name      = if (is.null(rtime_name)) NULL else rtime_name,
      method_levels   = if (is.null(rmet_name)) character(0) else method_levels,
      time_levels     = if (is.null(rtime_name)) character(0) else time_levels,
      has_interaction = interaction,
      Dmat_global     = Dsub
    )
  }
}

#' @title run_cpp
#' @description Wrapper for calling 'C++' backend for CCC estimation.
#' @keywords internal
run_cpp <- function(Xr, yr, subject, method_int, time_int, Laux, Z,
                    use_ar1, ar1_rho, max_iter, tol, conf_level) {
  ccc_vc_cpp(
    Xr = unname(Xr),
    yr = yr,
    subject = subject,
    method  = method_int,
    time    = time_int,
    nm = Laux$nm, nt = Laux$nt,
    max_iter = max_iter, tol = tol,
    conf_level = conf_level,
    Lr   = if (is.null(Laux$L)) NULL else unname(Laux$L),
    auxDr = if (is.null(Laux$Dm)) NULL else unname(Laux$Dm),
    Zr   = if (is.null(Z)) NULL else unname(Z),
    use_ar1 = use_ar1,
    ar1_rho = as.numeric(ar1_rho)
  )
}

#' @title estimate_rho
#' @description Estimate AR(1) correlation parameter rho by optimizing REML log-likelihood.
#' @keywords internal
estimate_rho <- function(Xr, yr, subject, method_int, time_int, Laux, Z,
                         rho_lo = -0.95, rho_hi = 0.95,
                         max_iter = 100, tol = 1e-6, conf_level = 0.95) {
  obj <- function(r) {
    fit <- run_cpp(Xr, yr, subject, method_int, time_int, Laux, Z,
                   use_ar1 = TRUE, ar1_rho = r,
                   max_iter = max_iter, tol = tol, conf_level = conf_level)
    ll <- suppressWarnings(as.numeric(fit[["reml_loglik"]]))
    if (!is.finite(ll)) return(Inf)  # never prefer invalid values
    -ll  # minimize negative REML loglik
  }
  oo <- optimize(obj, interval = c(rho_lo, rho_hi))
  list(rho = unname(oo$minimum), used_reml = TRUE)
}


#' @title ccc_lmm_reml_overall
#' @description Internal function to handle overall CCC estimation when `rmet` is NULL or has < 2 levels.
#' @keywords internal
ccc_lmm_reml_overall <- function(df, fml, ry, rind, rmet, rtime,
                                 slope, slope_var, slope_Z, drop_zero_cols,
                                 Dmat, ar, ar_rho, max_iter, tol,
                                 conf_level, verbose, digits, use_message,
                                 extra_label, ci,
                                 Dmat_type = c("time-avg","typical-visit","weighted-avg","weighted-sq"),
                                 Dmat_weights = NULL,
                                 Dmat_rescale = TRUE) {

  Dmat_type <- match.arg(Dmat_type)

  X <- model.matrix(fml, data = df)

  ## Determine present time levels in this (overall) fit
  lev_time_sub <- if (!is.null(rtime)) levels(df[[rtime]]) else character(0)

  ## Build/subset the time kernel Dsub only when there are ≥ 2 time levels
  if (!is.null(rtime) && length(lev_time_sub) >= 2L) {
    if (!is.null(Dmat)) {
      Dfull <- as.matrix(Dmat)
      ## In overall fits lev_time_sub are the global time levels; require conformity
      if (nrow(Dfull) != length(lev_time_sub) || ncol(Dfull) != length(lev_time_sub))
        stop("Dmat has incompatible dimension for present time levels (overall fit).")
      Dsub <- Dfull
      if (isTRUE(Dmat_rescale))
        Dsub <- .Dmat_normalise_mass(Dsub, length(lev_time_sub))
      ## soft symmetrisation for safety
      Dsub <- 0.5 * (Dsub + t(Dsub))
    } else {
      ## Construct from type/weights; align (named) weights if provided
      w_sub <- .align_weights_to_levels(Dmat_weights, lev_time_sub, lev_time_sub)
      Dsub  <- .Dmat_build_kernel(length(lev_time_sub),
                                  type    = Dmat_type,
                                  w       = w_sub,
                                  rescale = Dmat_rescale)
    }
  } else {
    Dsub <- NULL
  }

  ## infer "has_interaction" from model matrix columns
  has_interaction <- any(grepl(":", colnames(X), fixed = TRUE))

  Laux <- build_LDZ(
    colnames_X    = colnames(X),
    method_levels = if (is.null(rmet)) character(0) else levels(df[[rmet]]),
    time_levels   = if (is.null(rtime)) character(0) else levels(df[[rtime]]),
    Dsub          = Dsub,
    df_sub        = df,
    rmet_name     = rmet,
    rtime_name    = rtime,
    slope         = slope,
    interaction   = has_interaction,
    slope_var     = slope_var,
    drop_zero_cols= drop_zero_cols
  )

  method_int <- if (!is.null(rmet)  && nlevels(df[[rmet]])  >= 2L) as.integer(df[[rmet]])  else integer(0)
  time_int   <- if (!is.null(rtime) && nlevels(df[[rtime]]) >= 2L) as.integer(df[[rtime]]) else integer(0)
  Z_overall  <- if (!identical(slope, "custom")) Laux$Z else {
    if (is.null(slope_Z)) NULL else as.matrix(slope_Z)
  }

  ## estimate rho if requested
  rho_used <- if (identical(ar, "ar1") && is.na(ar_rho)) {
    er <- estimate_rho(X, df[[ry]], as.integer(df[[rind]]),
                       method_int, time_int, Laux, Z_overall,
                       max_iter = max_iter, tol = tol, conf_level = conf_level)
    if (!er$used_reml && isTRUE(verbose)) {
      message("Note: ccc_vc_cpp did not return 'reml_loglik'; rho estimated via a proxy objective.")
    }
    er$rho
  } else ar_rho

  ans <- tryCatch(
    run_cpp(X, df[[ry]], as.integer(df[[rind]]),
            method_int, time_int, Laux, Z_overall,
            use_ar1 = identical(ar, "ar1"),
            ar1_rho = if (identical(ar, "ar1")) rho_used else 0,
            max_iter = max_iter, tol = tol, conf_level = conf_level),
    error = function(e) {
      stop("ccc_vc_cpp failed on this dataset (near-singular tiny data): ",
           conditionMessage(e), call. = FALSE)
    }
  )

  if (isTRUE(verbose)) {
    .vc_message(ans, label = "Overall", nm = Laux$nm, nt = Laux$nt,
                conf_level = conf_level, digits = digits, use_message = use_message,
                extra_label = extra_label, ar = ar,
                ar_rho = if (identical(ar, "ar1")) rho_used else NA_real_)
  }

  lab <- "Overall"
  est_mat <- matrix(unname(ans$ccc), 1, 1, dimnames = list(lab, lab))

  if (isTRUE(ci)) {
    lwr_cpp <- num_or_na(ans[["lwr"]]); upr_cpp <- num_or_na(ans[["upr"]])
    if (is.na(lwr_cpp) || is.na(upr_cpp)) {
      se_cpp <- num_or_na(ans[["se_ccc"]])
      ci2 <- compute_ci_from_se(num_or_na(ans[["ccc"]]), se_cpp, conf_level)
      lwr_cpp <- ci2[1]; upr_cpp <- ci2[2]
    }
    out <- list(
      est    = est_mat,
      lwr.ci = matrix(lwr_cpp, 1, 1, dimnames = dimnames(est_mat)),
      upr.ci = matrix(upr_cpp, 1, 1, dimnames = dimnames(est_mat))
    )
    attr(out, "method")      <- "Variance Components REML"
    attr(out, "description") <- "Lin's CCC from random-effects LMM"
    attr(out, "package")     <- "matrixCorr"
    attr(out, "conf.level")  <- conf_level
    if (identical(ar, "ar1")) attr(out, "ar_rho") <- as.numeric(rho_used)

    ## variance-components attributes for overall (scalars)
    attr(out, "sigma2_subject")        <- num_or_na(ans[["sigma2_subject"]])
    attr(out, "sigma2_subject_method") <- num_or_na(ans[["sigma2_subject_method"]])
    attr(out, "sigma2_subject_time")   <- num_or_na(ans[["sigma2_subject_time"]])
    attr(out, "sigma2_error")          <- num_or_na(ans[["sigma2_error"]])
    attr(out, "sigma2_extra")          <- num_or_na(ans[["sigma2_extra"]])
    attr(out, "SB")                    <- num_or_na(ans[["SB"]])
    attr(out, "se_ccc")                <- num_or_na(ans[["se_ccc"]])

    class(out) <- c("ccc_lmm_reml", "matrixCorr_ccc_ci", "matrixCorr_ccc", "ccc")
    return(out)
  } else {
    out <- est_mat
    attr(out, "method")      <- "Variance Components REML"
    attr(out, "description") <- "Lin's CCC from random-effects LMM"
    attr(out, "package")     <- "matrixCorr"
    if (identical(ar, "ar1")) attr(out, "ar_rho") <- as.numeric(rho_used)

    ## variance-components attributes for overall (scalars)
    attr(out, "sigma2_subject")        <- num_or_na(ans[["sigma2_subject"]])
    attr(out, "sigma2_subject_method") <- num_or_na(ans[["sigma2_subject_method"]])
    attr(out, "sigma2_subject_time")   <- num_or_na(ans[["sigma2_subject_time"]])
    attr(out, "sigma2_error")          <- num_or_na(ans[["sigma2_error"]])
    attr(out, "sigma2_extra")          <- num_or_na(ans[["sigma2_extra"]])
    attr(out, "SB")                    <- num_or_na(ans[["SB"]])
    attr(out, "se_ccc")                <- num_or_na(ans[["se_ccc"]])

    class(out) <- c("ccc_lmm_reml", "matrixCorr_ccc", "ccc", "matrix")
    return(out)
  }
}

#' @title ccc_lmm_reml_pairwise
#' @description Internal function to handle pairwise CCC estimation for each method pair.
#' @keywords internal
ccc_lmm_reml_pairwise <- function(df, fml, ry, rind, rmet, rtime,
                                  slope, slope_var, slope_Z, drop_zero_cols,
                                  Dmat, ar, ar_rho, max_iter, tol,
                                  conf_level, verbose, digits, use_message,
                                  extra_label, ci, all_time_lvls,
                                  Dmat_type = c("time-avg","typical-visit","weighted-avg","weighted-sq"),
                                  Dmat_weights = NULL,
                                  Dmat_rescale = TRUE) {

  Dmat_type <- match.arg(Dmat_type)

  df[[rmet]] <- droplevels(df[[rmet]])
  method_levels <- levels(df[[rmet]])
  Lm <- length(method_levels)

  est_mat <- matrix(1,  Lm, Lm, dimnames = list(method_levels, method_levels))
  if (isTRUE(ci)) {
    lwr_mat <- matrix(NA_real_, Lm, Lm, dimnames = list(method_levels, method_levels))
    upr_mat <- matrix(NA_real_, Lm, Lm, dimnames = list(method_levels, method_levels))
  }

  # store rho per pair if estimated
  rho_mat <- matrix(NA_real_, Lm, Lm, dimnames = list(method_levels, method_levels))

  # variance component containers (per pair)
  vc_subject        <- matrix(NA_real_, Lm, Lm, dimnames = list(method_levels, method_levels))
  vc_subject_method <- matrix(NA_real_, Lm, Lm, dimnames = list(method_levels, method_levels))
  vc_subject_time   <- matrix(NA_real_, Lm, Lm, dimnames = list(method_levels, method_levels))
  vc_error          <- matrix(NA_real_, Lm, Lm, dimnames = list(method_levels, method_levels))
  vc_extra          <- matrix(NA_real_, Lm, Lm, dimnames = list(method_levels, method_levels))
  vc_SB             <- matrix(NA_real_, Lm, Lm, dimnames = list(method_levels, method_levels))
  vc_se_ccc         <- matrix(NA_real_, Lm, Lm, dimnames = list(method_levels, method_levels))

  for (i in 1:(Lm - 1L)) {
    for (j in (i + 1L):Lm) {
      m1 <- method_levels[i]; m2 <- method_levels[j]

      idx <- which(df[[rmet]] %in% c(m1, m2))
      subj_int   <- as.integer(df[[rind]][idx])
      y_sub      <- df[[ry]][idx]
      met_fac    <- droplevels(df[[rmet]][idx])        # exactly 2 levels
      time_fac   <- if (!is.null(rtime)) droplevels(df[[rtime]][idx]) else NULL

      Xp <- model.matrix(fml, data = df[idx, , drop = FALSE])

      # Present time levels in this pair
      lev_time_sub <- if (!is.null(time_fac)) levels(time_fac) else character(0)

      # Build/subset Dmat for this pair (only if ≥ 2 time levels)
      if (!is.null(rtime) && length(lev_time_sub) >= 2L) {
        if (!is.null(Dmat)) {
          Dfull <- as.matrix(Dmat)
          if (!is.null(all_time_lvls) && nrow(Dfull) == length(all_time_lvls) && ncol(Dfull) == length(all_time_lvls)) {
            pos  <- match(lev_time_sub, all_time_lvls)
            Dsub <- Dfull[pos, pos, drop = FALSE]
          } else if (nrow(Dfull) == length(lev_time_sub) && ncol(Dfull) == length(lev_time_sub)) {
            # already in the pair's time order
            Dsub <- Dfull
          } else {
            stop("Dmat has incompatible dimension for present time levels in a pairwise fit.")
          }
          if (isTRUE(Dmat_rescale))
            Dsub <- .Dmat_normalise_mass(Dsub, length(lev_time_sub))
          # soft symmetrisation for safety
          Dsub <- 0.5 * (Dsub + t(Dsub))
        } else {
          # Construct from type/weights; align (named) weights from global to present levels
          w_sub <- .align_weights_to_levels(Dmat_weights, lev_time_sub, all_time_lvls)
          Dsub  <- .Dmat_build_kernel(length(lev_time_sub),
                                      type    = Dmat_type,
                                      w       = w_sub,
                                      rescale = Dmat_rescale)
        }
      } else {
        Dsub <- NULL
      }

      # infer "has_interaction" from model matrix columns for this pair
      has_interaction <- any(grepl(":", colnames(Xp), fixed = TRUE))

      df_sub <- df[idx, , drop = FALSE]
      Laux <- build_LDZ(
        colnames_X    = colnames(Xp),
        method_levels = levels(met_fac),
        time_levels   = lev_time_sub,
        Dsub          = Dsub,
        df_sub        = df_sub,
        rmet_name     = rmet,
        rtime_name    = rtime,
        slope         = slope,
        interaction   = has_interaction,
        slope_var     = slope_var,
        drop_zero_cols= drop_zero_cols
      )

      method_int <- if (nlevels(met_fac)  >= 2L) as.integer(met_fac)  else integer(0)
      time_int   <- if (!is.null(time_fac) && nlevels(time_fac) >= 2L) as.integer(time_fac) else integer(0)
      Zp <- if (!identical(slope, "custom")) Laux$Z else {
        if (is.null(slope_Z)) NULL else as.matrix(slope_Z)[idx, , drop = FALSE]
      }

      # estimate rho for this pair if requested
      rho_used <- if (identical(ar, "ar1") && is.na(ar_rho)) {
        er <- estimate_rho(Xp, y_sub, subj_int, method_int, time_int, Laux, Zp,
                           max_iter = max_iter, tol = tol, conf_level = conf_level)
        if (!er$used_reml && isTRUE(verbose)) {
          message(sprintf("Note: REML log-lik not available for pair (%s,%s); rho estimated via a proxy objective.", m1, m2))
        }
        er$rho
      } else ar_rho
      rho_mat[i, j] <- rho_mat[j, i] <- as.numeric(rho_used)

      ans <- tryCatch(
        run_cpp(Xp, y_sub, subj_int, method_int, time_int, Laux, Zp,
                use_ar1 = identical(ar, "ar1"),
                ar1_rho = if (identical(ar, "ar1")) rho_used else 0,
                max_iter = max_iter, tol = tol, conf_level = conf_level),
        error = function(e) {
          warning(sprintf("ccc_vc_cpp failed for pair (%s, %s): %s",
                          m1, m2, conditionMessage(e)))
          NULL
        }
      )

      if (isTRUE(verbose) && !is.null(ans)) {
        .vc_message(ans, label = sprintf("Pair: %s vs %s", m1, m2),
                    nm = Laux$nm, nt = Laux$nt,
                    conf_level = conf_level, digits = digits,
                    use_message = use_message,
                    extra_label = extra_label, ar = ar,
                    ar_rho = if (identical(ar, "ar1")) rho_used else NA_real_)
      }

      val <- if (is.null(ans)) NA_real_ else unname(ans$ccc)
      est_mat[i, j] <- est_mat[j, i] <- val

      if (!is.null(ans)) {
        vc_subject[i, j]        <- vc_subject[j, i]        <- num_or_na(ans[["sigma2_subject"]])
        vc_subject_method[i, j] <- vc_subject_method[j, i] <- num_or_na(ans[["sigma2_subject_method"]])
        vc_subject_time[i, j]   <- vc_subject_time[j, i]   <- num_or_na(ans[["sigma2_subject_time"]])
        vc_error[i, j]          <- vc_error[j, i]          <- num_or_na(ans[["sigma2_error"]])
        extra <- ans[["sigma2_extra"]]
        vc_extra[i, j] <- vc_extra[j, i] <-
          if (length(extra) == 1L) as.numeric(extra) else mean(as.numeric(extra))
        vc_SB[i, j]             <- vc_SB[j, i]             <- num_or_na(ans[["SB"]])
        vc_se_ccc[i, j]         <- vc_se_ccc[j, i]         <- num_or_na(ans[["se_ccc"]])
      }

      if (isTRUE(ci)) {
        lwr_cpp <- num_or_na(if (!is.null(ans)) ans[["lwr"]] else NA_real_)
        upr_cpp <- num_or_na(if (!is.null(ans)) ans[["upr"]] else NA_real_)
        if (is.na(lwr_cpp) || is.na(upr_cpp)) {
          se_cpp <- num_or_na(if (!is.null(ans)) ans[["se_ccc"]] else NA_real_)
          ci2 <- compute_ci_from_se(num_or_na(val), se_cpp, conf_level)
          lwr_cpp <- ci2[1]; upr_cpp <- ci2[2]
        }
        lwr_mat[i, j] <- lwr_mat[j, i] <- lwr_cpp
        upr_mat[i, j] <- upr_mat[j, i] <- upr_cpp
      }
    }
  }

  diag(est_mat) <- 1
  if (isTRUE(ci)) {
    diag(lwr_mat) <- NA_real_
    diag(upr_mat) <- NA_real_
    out <- list(est = est_mat, lwr.ci = lwr_mat, upr.ci = upr_mat)
    attr(out, "method")      <- "Variance Components REML - pairwise"
    attr(out, "description") <- "Lin's CCC per method pair from random-effects LMM"
    attr(out, "package")     <- "matrixCorr"
    attr(out, "conf.level")  <- conf_level
    if (identical(ar, "ar1")) attr(out, "ar_rho") <- rho_mat

    # attach variance-component matrices
    attr(out, "sigma2_subject")        <- vc_subject
    attr(out, "sigma2_subject_method") <- vc_subject_method
    attr(out, "sigma2_subject_time")   <- vc_subject_time
    attr(out, "sigma2_error")          <- vc_error
    attr(out, "sigma2_extra")          <- vc_extra
    attr(out, "SB")                    <- vc_SB
    attr(out, "se_ccc")                <- vc_se_ccc

    class(out) <- c("ccc_lmm_reml", "matrixCorr_ccc_ci", "matrixCorr_ccc", "ccc")
    return(out)
  } else {
    out <- est_mat
    attr(out, "method")      <- "Variance Components REML - pairwise"
    attr(out, "description") <- "Lin's CCC per method pair from random-effects LMM"
    attr(out, "package")     <- "matrixCorr"
    if (identical(ar, "ar1")) attr(out, "ar_rho") <- rho_mat

    # attach variance-component matrices
    attr(out, "sigma2_subject")        <- vc_subject
    attr(out, "sigma2_subject_method") <- vc_subject_method
    attr(out, "sigma2_subject_time")   <- vc_subject_time
    attr(out, "sigma2_error")          <- vc_error
    attr(out, "sigma2_extra")          <- vc_extra
    attr(out, "SB")                    <- vc_SB
    attr(out, "se_ccc")                <- vc_se_ccc

    class(out) <- c("ccc_lmm_reml", "matrixCorr_ccc", "ccc", "matrix")
    return(out)
  }
}

#' @keywords internal
.Dmat_normalise_mass <- function(D, target_mass) {
  if (is.null(D)) return(NULL)
  one <- rep(1, nrow(D))
  mass <- as.numeric(t(one) %*% D %*% one)
  # leave as is; 'C++' will guard S_B
  if (!is.finite(mass) || mass <= 0) return(D)
  D * (target_mass / mass)
}

#' @keywords internal
.Dmat_build_kernel <- function(nt, type = c("time-avg","typical-visit","weighted-avg","weighted-sq"),
                               w = NULL, rescale = TRUE) {
  type <- match.arg(type)
  if (nt < 2) return(NULL)
  if (type %in% c("weighted-avg","weighted-sq")) {
    if (is.null(w)) w <- rep(1/nt, nt)
    w <- as.numeric(w)
    if (length(w) != nt) stop("Dmat_weights length must equal the number of present time levels.")
    if (!all(is.finite(w)) || any(w < 0)) stop("Dmat_weights must be non-negative and finite.")
    sw <- sum(w); if (sw <= 0) stop("Dmat_weights sums to zero.")
    w <- w / sw
  }
  D <- switch(type,
              "time-avg"      = (1/nt) * matrix(1, nt, nt),
              "typical-visit" = diag(nt),
              "weighted-avg"  = nt * (w %o% w),
              "weighted-sq"   = nt * diag(w)
  )
  if (rescale) D <- .Dmat_normalise_mass(D, nt)
  # symmetrise softly for safety
  0.5 * (D + t(D))
}

#' Align (optional named) weights to a subset of time levels
#' @keywords internal
.align_weights_to_levels <- function(w, present_lvls, all_lvls) {
  if (is.null(w)) return(NULL)
  if (!is.null(names(w))) {
    idx <- match(present_lvls, all_lvls)
    if (anyNA(idx)) stop("Present time levels not found in all_time_lvls.")
    w_all <- rep(NA_real_, length(all_lvls))
    w_all[seq_along(all_lvls)] <- w[all_lvls]
    w_sub <- w_all[idx]
    if (anyNA(w_sub)) stop("Missing weights for some present time levels.")
    as.numeric(w_sub)
  } else {
    if (length(w) != length(present_lvls))
      stop("Unnamed Dmat_weights must have length equal to the number of present time levels.")
    as.numeric(w)
  }
}

#' Print method for matrixCorr CCC objects
#'
#' @param x A `matrixCorr_ccc` or `matrixCorr_ccc_ci` object.
#' @param digits Number of digits for CCC estimates.
#' @param ci_digits Number of digits for CI bounds.
#' @param show_ci One of `"auto"`, `"yes"`, `"no"`.
#' @param ... Passed to underlying printers.
#' @export
#' @method print matrixCorr_ccc
print.matrixCorr_ccc <- function(x,
                                 digits = 4,
                                 ci_digits = 4,
                                 show_ci = c("auto", "yes", "no"),
                                 ...) {
  show_ci <- match.arg(show_ci)
  is_ci_obj <- inherits(x, "matrixCorr_ccc_ci") ||
    (is.list(x) && all(c("est", "lwr.ci", "upr.ci") %in% names(x)))

  if (is_ci_obj) {
    est <- as.matrix(x$est)
    lwr <- as.matrix(x$lwr.ci)
    upr <- as.matrix(x$lwr.ci); upr[,] <- x$upr.ci
  } else if (is.matrix(x)) {
    est <- as.matrix(x)
    lwr <- matrix(NA_real_, nrow(est), ncol(est), dimnames = dimnames(est))
    upr <- lwr
  } else {
    stop("Invalid object format for class 'ccc'.")
  }

  rn <- rownames(est); cn <- colnames(est)
  if (is.null(rn)) rn <- paste0("m", seq_len(nrow(est)))
  if (is.null(cn)) cn <- rn

  has_any_ci <- any(is.finite(lwr) | is.finite(upr))
  include_ci <- switch(show_ci, auto = has_any_ci, yes = TRUE, no = FALSE)

  cl <- suppressWarnings(as.numeric(attr(x, "conf.level")))
  if (include_ci && is.finite(cl)) {
    cat(sprintf("Concordance pairs (Lin's CCC, %g%% CI)\n\n", 100 * cl))
  } else {
    cat("Concordance pairs (Lin's CCC)\n\n")
  }

  if (nrow(est) == 1L && ncol(est) == 1L) {
    df <- data.frame(
      method1  = rn[1],
      method2  = cn[1],
      estimate = formatC(est[1,1], format = "f", digits = digits),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    if (include_ci) {
      df$lwr <- ifelse(is.na(lwr[1,1]), NA,
                       formatC(lwr[1,1], format = "f", digits = ci_digits))
      df$upr <- ifelse(is.na(upr[1,1]), NA,
                       formatC(upr[1,1], format = "f", digits = ci_digits))
    }
    print(df, row.names = FALSE, right = FALSE, ...)
    return(invisible(x))
  }

  rows <- vector("list", nrow(est) * (ncol(est) - 1L) / 2L); k <- 0L
  for (i in seq_len(nrow(est) - 1L)) {
    for (j in (i + 1L):ncol(est)) {
      k <- k + 1L
      row <- list(
        method1  = rn[i],
        method2  = cn[j],
        estimate = formatC(est[i, j], format = "f", digits = digits)
      )
      if (include_ci) {
        row$lwr <- ifelse(is.na(lwr[i, j]), NA,
                          formatC(lwr[i, j], format = "f", digits = ci_digits))
        row$upr <- ifelse(is.na(upr[i, j]), NA,
                          formatC(upr[i, j], format = "f", digits = ci_digits))
      }
      rows[[k]] <- row
    }
  }

  df <- do.call(rbind.data.frame, rows); rownames(df) <- NULL
  print(df, row.names = FALSE, right = FALSE, ...)
  invisible(x)
}

#' Print method for matrixCorr CCC objects with CIs
#'
#' @inheritParams print.matrixCorr_ccc
#' @export
#' @method print matrixCorr_ccc_ci
print.matrixCorr_ccc_ci <- function(x, ...) {
  print.matrixCorr_ccc(x, ...)
}

#' S3 print for legacy class `ccc_ci`
#'
#' For compatibility with objects that still carry class `"ccc_ci"`.
#' @inheritParams print.matrixCorr_ccc
#' @export
#' @method print ccc_ci
print.ccc_ci <- function(x, ...) {
  print.matrixCorr_ccc(x, ...)
}

#' @title Summary Method for `ccc_lmm_reml` Objects
#'
#' @description Produces a detailed summary of a `"ccc_lmm_reml"` object, including
#' Lin's CCC estimates and associated variance component estimates per method pair.
#'
#' @param object An object of class `"ccc_lmm_reml"`, as returned by [ccc_lmm_reml()].
#' @param digits Integer; number of decimal places to round CCC estimates and components.
#' @param ci_digits Integer; decimal places for confidence interval bounds.
#' @param show_ci Character string indicating whether to show confidence intervals:
#'   `"auto"` (default) shows only if non-NA CIs exist, `"yes"` always shows CIs,
#'   `"no"` never shows them.
#' @param ... Additional arguments (ignored).
#'
#' @return A data frame of class `"summary.ccc_lmm_reml"` with columns:
#'   \code{method1}, \code{method2}, \code{estimate}, and optionally \code{lwr}, \code{upr},
#'   as well as variance component estimates: \code{sigma2_subject}, \code{sigma2_subject_method},
#'   \code{sigma2_subject_time}, \code{sigma2_error}, \code{sigma2_extra}, \code{SB}, \code{se_ccc}.
#'
#' @export
#' @method summary ccc_lmm_reml
summary.ccc_lmm_reml <- function(object,
                                 digits = 4,
                                 ci_digits = 2,
                                 show_ci = c("auto", "yes", "no"),
                                 ...) {
  show_ci <- match.arg(show_ci)

  # Base CCC summary (handles CI and formatting choices)
  base_summary <- summary.ccc(object, digits = digits,
                              ci_digits = ci_digits, show_ci = show_ci)

  # Pull the estimate matrix to know the size/order of pairs
  est_mat <- if (is.list(object) && !is.null(object$est)) {
    as.matrix(object$est)
  } else {
    as.matrix(object)
  }

  rn <- rownames(est_mat); cn <- colnames(est_mat)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(est_mat)))
  if (is.null(cn)) cn <- as.character(seq_len(ncol(est_mat)))

  # Helper: extract VC values in the same order as summary.ccc builds rows
  extract_pairs_vec <- function(val) {
    # overall 1x1
    if (nrow(est_mat) == 1L && ncol(est_mat) == 1L) {
      if (is.null(val)) return(NA_real_)
      if (is.matrix(val)) {
        return(suppressWarnings(as.numeric(val[1, 1])))
      } else {
        return(suppressWarnings(as.numeric(val[1])))
      }
    }

    # pairwise >= 2 methods
    out <- numeric(nrow(est_mat) * (ncol(est_mat) - 1L) / 2L)
    k <- 0L

    # Named access if available; otherwise positional
    use_named <- is.matrix(val) && !is.null(rownames(val)) && !is.null(colnames(val))

    for (i in seq_len(nrow(est_mat) - 1L)) {
      for (j in (i + 1L):ncol(est_mat)) {
        k <- k + 1L
        if (is.null(val)) {
          out[k] <- NA_real_
        } else if (is.matrix(val)) {
          if (use_named) {
            out[k] <- suppressWarnings(as.numeric(val[rn[i], cn[j]]))
          } else {
            out[k] <- suppressWarnings(as.numeric(val[i, j]))
          }
        } else {
          # scalar or weird length; recycle if length 1, else NA
          vv <- suppressWarnings(as.numeric(val))
          out[k] <- if (length(vv) == 1L) vv else NA_real_
        }
      }
    }
    out
  }

  vc_names <- c("sigma2_subject",
                "sigma2_subject_method",
                "sigma2_subject_time",
                "sigma2_error",
                "sigma2_extra",
                "SB",
                "se_ccc")

  # Build VC columns aligned to the pair order
  vc_cols <- lapply(vc_names, function(nm) extract_pairs_vec(attr(object, nm)))
  names(vc_cols) <- vc_names

  # Bind and round
  out <- cbind(base_summary, lapply(vc_cols, function(x) round(x, digits)))
  for (nm in vc_names) out[[nm]] <- as.numeric(out[[nm]])

  class(out) <- c("summary.ccc_lmm_reml", "data.frame")
  attr(out, "conf.level") <- attr(base_summary, "conf.level")
  attr(out, "has_ci")     <- attr(base_summary, "has_ci")
  attr(out, "digits")     <- digits
  attr(out, "ci_digits")  <- ci_digits
  out
}
