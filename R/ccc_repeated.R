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
#'   (default \code{TRUE}).
#' @param max_iter Integer. Maximum iterations for variance-component updates
#'   (default \code{100}).
#' @param tol Numeric. Convergence tolerance on parameter change
#'   (default \code{1e-6}).
#' @param Dmat Optional \eqn{n_t \times n_t} numeric matrix to weight/contrast
#'   time levels when computing the fixed-effect dispersion term \eqn{S_B}.
#'   Defaults to the identity.
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
#'   When \code{ar = "ar1"} and \code{ar_rho = NA}, a one-dimensional profile
#'   optimization is used to estimate \eqn{\rho} (see \strong{AR(1)} in
#'   \strong{Details}). Default \code{NA_real_}.
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
#' Here \eqn{Z} is the subject-structured random-effects design and \eqn{G} is
#' block-diagonal at the subject level with the following \emph{per-subject}
#' parameterization:
#' \itemize{
#'   \item one random intercept with variance \eqn{\sigma_A^2};
#'   \item optionally, \emph{method} deviations (one column per method level)
#'         with a common variance \eqn{\sigma_{A\times M}^2} and zero
#'         covariances across levels (i.e., multiple of an identity);
#'   \item optionally, \emph{time} deviations (one column per time level)
#'         with a common variance \eqn{\sigma_{A\times T}^2} and zero
#'         covariances across levels;
#'   \item optionally, an \emph{extra} random effect aligned with \eqn{Z}
#'         (random slope): variance \eqn{\sigma_Z^2} times an identity on the
#'         \eqn{Z}-columns (see \strong{Random-slope \eqn{Z}}).
#' }
#' The fixed-effects design is \code{~ 1 + rmet + rtime} and, if
#' \code{interaction=TRUE}, \code{+ rmet:rtime}.
#'
#' \strong{Residual correlation \eqn{R} (regular, equally spaced time).}
#' With \code{ar="none"}, \eqn{R=\sigma_E^2 I}. With \code{ar="ar1"}, within-subject
#' residuals follow a \emph{discrete} AR(1) process along the visit index after sorting
#' by increasing time level. Ties retain input order, and any \code{NA} time code breaks
#' the series so each contiguous block of non-\code{NA} times forms a run. The correlation
#' between \emph{adjacent observed visits} in a run is \eqn{\rho}; we do not use calendar-time
#' gaps. Internally we work with the \emph{precision} of the AR(1) correlation:
#' for a run of length \eqn{L\ge 2}, the tridiagonal inverse has
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
#' Solves/inversions use symmetric-PD routines with a tiny diagonal "jitter" and
#' a pseudo-inverse fallback when needed.
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
#' Computations simply augment \eqn{U_i} to \eqn{\tilde U_i=[U_i\ Z_i]} and
#' \eqn{G^{-1}} to \eqn{\tilde G^{-1}} by appending a diagonal block
#' \eqn{\sigma_Z^{-2} I_{q_Z}}. The EM updates then include
#' \deqn{ \sigma_Z^{2\,(new)} \;=\; \frac{1}{m\,q_Z}
#'       \sum_i \sum_{j=1}^{q_Z}\!\Big( b_{i,\text{extra},j}^2 +
#'       (M_i^{-1})_{\text{extra},jj} \Big)
#'       \quad (\text{if } q_Z>0). }
#' \emph{Interpretation:} \eqn{\sigma_Z^2} represents additional within-subject
#' variability explained by the slope regressor(s). By design it \emph{does not}
#' enter the CCC formula below, where CCC targets agreement across methods/time, not
#' variability along a subject- or method-specific slope.
#'
#' \strong{EM-style variance-component updates.} With current \eqn{\hat\beta},
#' residuals \eqn{r_i = y_i - X_i\hat\beta} are formed. The BLUPs and
#' conditional covariances are
#' \deqn{ b_i \;=\; M_i^{-1}\,(U_i^\top R_i^{-1} r_i), \qquad
#'       \mathrm{Var}(b_i\mid y) \;=\; M_i^{-1}. }
#' Expected squares yield closed-form updates:
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
#'       \Big( r_i^\top C_i(\rho)^{-1} r_i +
#'       \mathrm{tr}\!\big(M_i^{-1}\,U_i^\top C_i(\rho)^{-1} U_i\big) \Big), }
#' where \eqn{C_i(\rho)^{-1}} is the AR(1) precision built on the
#' time-ordered runs described above (and equals \eqn{I} for iid residuals).
#' Iterate until the \eqn{\ell_1} change across components is \eqn{<}
#' \code{tol} or \code{max_iter} is reached.
#'
#' \strong{Fixed-effect dispersion \eqn{\mathbf{S_B}}.} Method dispersion is
#' computed from \eqn{\hat\beta} and \eqn{\mathrm{Var}(\hat\beta)} with a
#' contrast matrix \eqn{L} (columns encode pairwise method differences within
#' each time level) and an optional time-weighting matrix \eqn{\mathrm{D_m}}:
#' \deqn{ S_B \;=\;
#'  \frac{\big(L^\top \hat\beta\big)^\top\,\mathrm{D_m}\,
#'  \big(L^\top \hat\beta\big)
#'        \;-\; \mathrm{tr}\!\Big(\big(L\,\mathrm{D_m}\,L^\top\big)\,
#'        \mathrm{Var}(\hat\beta)\Big)}
#'       {\,nm\,(nm-1)\,\max(nt,1)\,}, }
#' truncated at 0. The helper \code{build_L_Dm_cpp} constructs \eqn{L} so it
#' aligns exactly with the columns of \eqn{X=\mathrm{model.matrix}(\cdot)}.
#' For exactly two methods (\eqn{nm=2}), a fast path builds \eqn{L} directly
#' from the design's column names, where, with interaction, the per-time
#' difference at time \eqn{j} is
#' \eqn{\beta_{\text{met2}}+\beta_{\text{met2:time}_j}} (baseline time uses
#' \eqn{\beta_{\text{met2}}}); while without interaction, the same
#'  \eqn{\beta_{\text{met2}}} is used for all times.
#'
#' \strong{Time-averaging for CCC (regular visits).}
#' The reported CCC targets agreement of the \emph{time-averaged} measurements
#' per method within subject. Averaging over \eqn{T} non-\code{NA}
#' \emph{visits} (equally spaced, discrete AR) leaves subject
#' (\eqn{\sigma_A^2}) and subject\eqn{\times}method (\eqn{\sigma_{A\times M}^2})
#' components unchanged, while shrinking time-varying components by
#' \deqn{ \kappa_T^{(g)} \;=\; 1/T \quad (\text{subject}\times\text{time}), }
#' \deqn{ \kappa_T^{(e)} \;=\; \frac{T + 2\sum_{k=1}^{T-1} (T-k)\,\rho^k}{T^2}
#'        \quad \text{(residual AR(1) on visit index; equals }1/T\text{ if iid).} }
#' With unbalanced \eqn{T}, the implementation uses the average of the
#' per-(subject,method) \eqn{\kappa} values across the pairs contributing
#' to CCC. (Numerical guard: after averaging,
#' \eqn{\kappa_T^{(e)}} is clamped to \eqn{[10^{-12},\,1]}.)


#'
#' \strong{Concordance correlation coefficient.} The CCC used is defined by
#' \deqn{ \mathrm{CCC} \;=\;
#'       \frac{\sigma_A^2 + \kappa_T^{(g)}\,\sigma_{A\times T}^2}
#'            {\sigma_A^2 + \sigma_{A\times M}^2 +
#'             \kappa_T^{(g)}\,\sigma_{A\times T}^2 + S_B +
#'             \kappa_T^{(e)}\,\sigma_E^2}. }
#' There are special cases when there is no method factor (or a single level),
#' then \eqn{S_B=0} and \eqn{\sigma_{A\times M}^2=0}; if there is no
#' time factor (or a single level), then \eqn{\sigma_{A\times T}^2=0} and
#' \eqn{\kappa_T^{(g)}} is unused. When \eqn{T=1} (one time point) or
#' \eqn{\rho=0}, both \eqn{\kappa}-factors equal 1 and the formula reduces to the
#' familiar single-occasion expression.
#' The extra random-effect variance \eqn{\sigma_Z^2} (if used) is \emph{not}
#' included, where CCC targets agreement across methods/time, not variability along
#' the user-specified slope.
#'
#' \strong{CIs / SEs (delta method for CCC).}
#' Let
#' \deqn{ \theta \;=\; \big(\sigma_A^2,\ \sigma_{A\times M}^2,\
#' \sigma_{A\times T}^2,\ \sigma_E^2,\ S_B\big)^\top, }
#' and write the concordance as
#' \deqn{ \mathrm{CCC}(\theta) \;=\; \frac{N}{D}
#'       \;=\; \frac{\sigma_A^2 + \kappa_T^{(g)}\,\sigma_{A\times T}^2}
#'                    {\sigma_A^2 + \sigma_{A\times M}^2 +
#'                    \kappa_T^{(g)}\,\sigma_{A\times T}^2 + S_B +
#'                    \kappa_T^{(e)}\,\sigma_E^2}. }
#'
#' A first-order (large-sample) standard error follows from the delta method:
#' \deqn{ \mathrm{Var}\{\widehat{\mathrm{CCC}}\}
#'       \;\approx\; \nabla \mathrm{CCC}(\hat\theta)^\top\,
#'                   \mathrm{Var}(\hat\theta)\,
#'                   \nabla \mathrm{CCC}(\hat\theta), }
#' with gradient components (using \eqn{N} and \eqn{D} as above)
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
#'        s_i \;=\; \frac{r_i^\top C_i(\rho)^{-1} r_i +
#'        \mathrm{tr}\!\big(M_i^{-1}U_i^\top C_i(\rho)^{-1} U_i\big)}{n_i}, }
#' where \eqn{b_i = M_i^{-1}(U_i^\top R_i^{-1} r_i)} and
#' \eqn{M_i = G^{-1} + U_i^\top R_i^{-1} U_i}.
#' With \eqn{m} subjects, we form the empirical covariance of the stacked
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
#' with \eqn{A_{\!fix}=L\,\mathrm{D_m}\,L^\top}.
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
#' This function \emph{requires a fixed} \eqn{\rho\in(-0.999,0.999)} when \code{ar="ar1"}; it does
#' not estimate \eqn{\rho} internally. Higher-level wrappers may profile \eqn{\rho} by refitting
#' across a grid and selecting by (R)EML or another criterion, but that is outside this core.
#'
#' @return
#' \itemize{
#'   \item If \code{rmet} is \code{NULL} or has a single level, an object of
#'     class \code{c("ccc","ccc_ci")} (when \code{ci=TRUE}) or
#'     \code{c("ccc","matrix")} with a \eqn{1\times 1} matrix containing the
#'     overall CCC estimate.
#'   \item If \code{rmet} has \eqn{L\geq 2} levels, a symmetric \eqn{L\times L}
#'     matrix with pairwise CCC estimates between methods (diagonal set to 1).
#'     When \code{ci=TRUE}, \code{lwr.ci} and \code{upr.ci} matrices are
#'     included.
#' }
#' In all cases, attributes \code{"method"}, \code{"description"},
#' \code{"package"}, and (if \code{ci=TRUE}) \code{"conf.level"} are set.
#' When \code{ar="ar1"}, an additional attribute \code{"ar_rho"} is attached:
#' a scalar (overall) or an \eqn{L\times L} matrix (pairwise) with the
#' \eqn{\rho} values used/estimated.
#'
#' @section Notes on stability and performance:
#' All per-subject solves are \eqn{\,r\times r} with \eqn{r=1+nm+nt+q_Z}, so cost
#' scales with the number of subjects and the fixed-effects dimension rather
#' than the total number of observations. Solvers use symmetric-PD paths with
#' a small diagonal ridge and pseudo-inverse fallback, which helps for
#' tiny/unbalanced subsets and near-boundary estimates. Very small samples or
#' extreme imbalance can still make \eqn{S_B} numerically delicate; negative
#' estimates are truncated to 0 by construction.
#' For AR(1), observations are first ordered by time within subject before building the
#' run-wise tridiagonal precision; \code{NA} time codes break the run, and gaps between
#' factor levels are treated as regular steps (we do not use elapsed time).
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
#' Carrasco JL, Jover L (2003). Estimating the concordance correlation coefficient:
#' a new approach. \emph{Computational Statistics & Data Analysis}, 47(4): 519-539.
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
#' ## Random slope by subject: create a centered numeric time within subject
#' #--------------------------------------------------------------------
#' dat$t_num <- as.integer(dat$time)
#' dat$t_c   <- ave(dat$t_num, dat$id, FUN = function(v) v - mean(v))
#' ccc_lmm_reml(dat, "y", "id", rmet = "method", rtime = "time",
#'              slope = "subject", slope_var = "t_c",
#'              ar = "ar1", ar_rho = NA_real_, verbose = TRUE)
#'
#' #--------------------------------------------------------------------
#' ## Three methods - pairwise CCCs
#' #--------------------------------------------------------------------
#' set.seed(2)
#' n_subj <- 40
#' id2     <- factor(rep(seq_len(n_subj), times = 3))
#' method2 <- factor(rep(c("A","B","C"), each = n_subj))
#' sigA <- 1.2; sigE <- 0.6
#' u  <- rnorm(n_subj, 0, sqrt(sigA))
#' mu <- c(A = 0.00, B = 0.15, C = -0.10)
#' e  <- rnorm(3 * n_subj, 0, sqrt(sigE))
#' y2 <- u[as.integer(id2)] + unname(mu[method2]) + e
#' dat3 <- data.frame(y = y2, id = id2, method = method2)
#' ccc_lmm_reml(dat3, "y", "id", rmet = "method", verbose = TRUE)
#'
#' # ------------------------------------------------------------------
#' # AR(1) residual correlation (fixed rho, threads forced to 1)
#' # When needed: repeated measures over time with serially correlated
#' # residuals within subject (e.g., values drift smoothly across visits).
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
#'     # small linear trend so AR(1) isn't swallowed by fixed effects
#'     y[idx] <- beta0 + beta_t*(seq_len(n_time) - mean(seq_len(n_time))) + e
#'   }
#'   dat_ar <- data.frame(y = y, id = id, time = tim)
#'   # Fit with AR(1) and rho fixed (nonzero). Estimation of rho is not
#'   # implemented yet; use a plausible value (e.g., 0.4–0.8) for sensitivity.
#'   ccc_lmm_reml(dat_ar, ry = "y", rind = "id", rtime = "time",
#'                ar = "ar1", ar_rho = 0.6, verbose = TRUE)
#'
#' # ------------------------------------------------------------------
#' # Random slope by SUBJECT
#' # When needed: each subject shows a systematic linear change over time
#' # (e.g., individual-specific trends), regardless of method.
#' # ------------------------------------------------------------------
#' set.seed(2)
#' n_subj <- 60; n_time <- 4
#' id  <- factor(rep(seq_len(n_subj), each = 2 * n_time))
#' tim <- factor(rep(rep(seq_len(n_time), times = 2), times = n_subj))
#' method <- factor(rep(rep(c("A","B"), each = n_time), times = n_subj))
#' # subject-specific slopes around 0
#' subj <- as.integer(id)
#' slope_i <- rnorm(n_subj, 0, 0.15)
#' slope_vec <- slope_i[subj]
#' base <- rnorm(n_subj, 0, 1.0)[subj]
#' tnum <- as.integer(tim)
#' y <- base + 0.3*(method=="B") + slope_vec*(tnum - mean(seq_len(n_time))) +
#'      rnorm(length(id), 0, 0.5)
#' dat_s <- data.frame(y, id, method, time = tim)
#' # center time within subject (recommended for random slopes)
#' dat_s$t_num <- as.integer(dat_s$time)
#' dat_s$t_c   <- ave(dat_s$t_num, dat_s$id, FUN = function(v) v - mean(v))
#' ccc_lmm_reml(dat_s, "y", "id", rmet = "method", rtime = "time",
#'              slope = "subject", slope_var = "t_c", verbose = TRUE)
#'
#' # ------------------------------------------------------------------
#' # Random slope by METHOD
#' # When needed: methods drift differently across time (e.g., biases that
#' # change with time are method-specific).
#' # ------------------------------------------------------------------
#' set.seed(3)
#' n_subj <- 60; n_time <- 4
#' id  <- factor(rep(seq_len(n_subj), each = 2 * n_time))
#' tim <- factor(rep(rep(seq_len(n_time), times = 2), times = n_subj))
#' method <- factor(rep(rep(c("A","B"), each = n_time), times = n_subj))
#' # method-specific slopes: A ~ 0, B ~ positive drift
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
#' # Random slopes for SUBJECT *and* METHOD (custom Z)
#' # When needed: subjects have their own time trends, AND methods carry
#' # additional method-specific trends. Supply a custom Z with both parts.
#' # ------------------------------------------------------------------
#' set.seed(4)
#' n_subj <- 50; n_time <- 4
#' id  <- factor(rep(seq_len(n_subj), each = 2 * n_time))
#' tim <- factor(rep(rep(seq_len(n_time), times = 2), times = n_subj))
#' method <- factor(rep(rep(c("A","B"), each = n_time), times = n_subj))
#' subj <- as.integer(id)
#' # subject slopes + extra slope on method B
#' slope_subj <- rnorm(n_subj, 0, 0.12)[subj]
#' slope_B    <- ifelse(method=="B", 0.18, 0.00)
#' tnum <- as.integer(tim)
#' base <- rnorm(n_subj, 0, 1.0)[subj]
#' y <- base + 0.3*(method=="B") +
#'      (slope_subj + slope_B) * (tnum - mean(seq_len(n_time))) +
#'      rnorm(length(id), 0, 0.5)
#' dat_both <- data.frame(y, id, method, time = tim)
#' # build Z = [subject_slope | method_A_slope | method_B_slope]
#' dat_both$t_num <- as.integer(dat_both$time)
#' dat_both$t_c   <- ave(dat_both$t_num, dat_both$id, FUN = function(v) v - mean(v))
#' MM <- model.matrix(~ 0 + method, data = dat_both)  # one col per method
#' Z_custom <- cbind(
#'   subj_slope   = dat_both$t_c,            # subject slope
#'   MM * dat_both$t_c                       # method-specific slopes
#' )
#' ccc_lmm_reml(dat_both, "y", "id", rmet = "method", rtime = "time",
#'              slope = "custom", slope_Z = Z_custom, verbose = TRUE)
#'
#' @author Thiago de Paula Oliveira
#' @importFrom stats as.formula model.matrix setNames qnorm optimize
#' @export
ccc_lmm_reml <- function(data, ry, rind,
                         rmet = NULL, rtime = NULL, interaction = TRUE,
                         max_iter = 100, tol = 1e-6,
                         Dmat = NULL, ci = FALSE, conf_level = 0.95,
                         verbose = FALSE, digits = 4, use_message = TRUE,
                         ar = c("none", "ar1"),
                         ar_rho = NA_real_,
                         slope = c("none", "subject", "method", "custom"),
                         slope_var = NULL,
                         slope_Z = NULL,
                         drop_zero_cols = TRUE) {

  ar    <- match.arg(ar)
  slope <- match.arg(slope)

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
                                extra_label, ci))
  } else {
    return(ccc_lmm_reml_pairwise(df, fml, ry, rind, rmet, rtime,
                                 slope, slope_var, slope_Z, drop_zero_cols,
                                 Dmat, ar, ar_rho, max_iter, tol,
                                 conf_level, verbose, digits, use_message,
                                 extra_label, ci, all_time_lvls))
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
  fmt <- function(x) if (is.na(x)) "NA" else formatC(x, format = "f", digits = digits)
  colw <- 38L
  v <- function(s, x) sprintf("  %-*s : %s", colw, s, fmt(as.numeric(x)))
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
#' @description Wrapper for calling C++ backend for CCC estimation.
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
  used_reml <- TRUE
  obj <- function(r) {
    fit <- run_cpp(Xr, yr, subject, method_int, time_int, Laux, Z,
                   use_ar1 = TRUE, ar1_rho = r,
                   max_iter = max_iter, tol = tol, conf_level = conf_level)
    ll <- fit[["reml_loglik"]]
    if (is.null(ll) || !is.finite(ll)) {
      used_reml <<- FALSE
      se_ccc <- suppressWarnings(as.numeric(fit[["se_ccc"]]))
      SB     <- suppressWarnings(as.numeric(fit[["SB"]]))
      if (!is.finite(se_ccc)) se_ccc <- 0
      if (!is.finite(SB))     SB     <- 0
      return(se_ccc + SB)
    }
    -as.numeric(ll)
  }
  oo <- optimize(obj, interval = c(rho_lo, rho_hi))
  list(rho = unname(oo$minimum), used_reml = used_reml)
}

#' @title ccc_lmm_reml_overall
#' @description Internal function to handle overall CCC estimation when `rmet` is NULL or has < 2 levels.
#' @keywords internal
ccc_lmm_reml_overall <- function(df, fml, ry, rind, rmet, rtime,
                                 slope, slope_var, slope_Z, drop_zero_cols,
                                 Dmat, ar, ar_rho, max_iter, tol,
                                 conf_level, verbose, digits, use_message,
                                 extra_label, ci) {

  X <- model.matrix(fml, data = df)

  # For overall case, pass Dmat only if time has >= 2 levels
  Dsub <- if (!is.null(Dmat) && !is.null(rtime) && nlevels(df[[rtime]]) >= 2L) {
    as.matrix(Dmat)
  } else NULL

  # infer "has_interaction" from model matrix columns
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

  # estimate rho if requested
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

    # variance-components attributes for overall (scalars)
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

    # variance-components attributes for overall (scalars)
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
                                  extra_label, ci, all_time_lvls) {

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

      # subset Dmat to present time levels in this pair (if any)
      lev_time_sub <- if (!is.null(time_fac)) levels(time_fac) else character(0)
      Dsub <- if (!is.null(Dmat) && length(lev_time_sub) >= 2L) {
        pos <- match(lev_time_sub, all_time_lvls)
        as.matrix(Dmat[pos, pos, drop = FALSE])
      } else NULL

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
        vc_extra[i, j]          <- vc_extra[j, i]          <- num_or_na(ans[["sigma2_extra"]])
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
