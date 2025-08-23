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
#' @param data A data frame.
#' @param ry Character. Response variable name.
#' @param rind Character. Subject ID variable name (random intercept).
#' @param rmet Character or \code{NULL}. Optional column name of method factor
#' (added to fixed effects).
#' @param rtime Character or \code{NULL}. Optional column name of time factor
#' (added to fixed effects).
#' @param interaction Logical. Include \code{method:time} interaction?
#' (default \code{TRUE}).
#' @param max_iter Integer. Maximum iterations for variance-component updates
#' (default \code{100}).
#' @param tol Numeric. Convergence tolerance on parameter change
#' (default \code{1e-6}).
#' @param Dmat Optional \eqn{n_t \times n_t} numeric matrix to weight/contrast
#' time levels when computing the fixed-effect dispersion term \eqn{S_B}.
#' Defaults to the identity.
#' @param ci Logical. If \code{TRUE}, return a CI container with NA limits
#' (see \strong{CIs} note below).
#' @param conf_level Numeric in \eqn{(0,1)}. Confidence level when
#' \code{ci=TRUE} (default \code{0.95}).
#' @param verbose Logical. If \code{TRUE}, prints a structured summary of the
#'   fitted variance components and \eqn{S_B} for each fit (overall or
#'   pairwise). Default \code{FALSE}.
#' @param digits Integer (\(\geq 0\)). Number of decimal places to use in the
#'   printed summary when \code{verbose=TRUE}. Default \code{4}.
#' @param use_message Logical. When \code{verbose=TRUE}, choose the printing
#'   mechanism, where \code{TRUE} uses \code{message()} (respects \code{sink()},
#'   easily suppressible via \code{suppressMessages()}), whereas \code{FALSE}
#'   uses \code{cat()} to \code{stdout}. Default \code{TRUE}.
#'
#' @details
#' For measurement \eqn{y_{ij}} on subject \eqn{i} under fixed
#' levels (method, time), we fit
#' \deqn{ y = X\beta + Zu + \varepsilon,\qquad
#'        u \sim N(0,\,G),\ \varepsilon \sim N(0,\,\sigma_E^2 I). }
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
#'         covariances across levels.
#' }
#' The fixed-effects design is \code{~ 1 + rmet + rtime} and, if
#' \code{interaction=TRUE}, \code{+ rmet:rtime}.
#'
#' \strong{Per-subject Woodbury system.} For subject \eqn{i} with \eqn{n_i}
#' rows, define the per-subject random-effects design \eqn{U_i} (columns:
#' intercept, method indicators, time indicators; dimension
#' \eqn{\,r=1+nm+nt\,}). The core never forms
#' \eqn{V_i = \sigma_E^2 I_{n_i} + U_i G U_i^\top} explicitly. Instead,
#' it builds
#' \deqn{ M_i \;=\; G^{-1} \;+\; \frac{1}{\sigma_E^2}\,U_i^\top U_i, }
#' and accumulates GLS blocks via rank-\eqn{r} corrections using
#' \eqn{\,V_i^{-1} = \tfrac{1}{\sigma_E^2}\big(I_{n_i} - U_i M_i^{-1}
#' U_i^\top / \sigma_E^2\big)}:
#' \deqn{ X^\top V^{-1} X \;=\; \sum_i \frac{1}{\sigma_E^2}\Big[
#'        X_i^\top X_i \;-\; (X_i^\top U_i)\,
#'        M_i^{-1}\,(U_i^\top X_i)/\sigma_E^2 \Big], }
#' \deqn{ X^\top V^{-1} y \;=\; \sum_i \frac{1}{\sigma_E^2}\Big[
#'        X_i^\top y_i \;-\; (X_i^\top U_i)\,M_i^{-1}\,
#'        (U_i^\top y_i)/\sigma_E^2 \Big]. }
#' Solves/inversions use symmetric-PD routines with a tiny diagonal "jitter" and
#' a pseudo-inverse fallback when needed.
#'
#' \strong{EM-style variance-component updates.} With current \eqn{\hat\beta},
#' residuals \eqn{r_i = y_i - X_i\hat\beta} are formed. The BLUPs and
#' conditional covariances are
#' \deqn{ b_i \;=\; M_i^{-1}\,(U_i^\top r_i \!/\! \sigma_E^2), \qquad
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
#'       \Big( \| r_i - U_i b_i \|^2 + \mathrm{tr}\!\big(M_i^{-1} U_i^\top
#'       U_i\big) \Big). }
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
#' truncated at 0. The helper \code{\link{build_L_Dm_cpp}} constructs \eqn{L} so it
#' aligns exactly with the columns of \eqn{X=\mathrm{model.matrix}(\cdot)}
#' passed to 'C++'.
#' For exactly two methods (\eqn{nm=2}), a fast path builds \eqn{L} directly
#' from the design's column names, where, with interaction, the per-time
#' difference at time \eqn{j} is
#' \eqn{\beta_{\text{met2}}+\beta_{\text{met2:time}_j}} (baseline time uses
#' \eqn{\beta_{\text{met2}}}); while without interaction, the same
#'  \eqn{\beta_{\text{met2}}} is used for all times.
#'
#' \strong{Concordance correlation coefficient.} The CCC used is defined by
#' \deqn{ \mathrm{CCC} \;=\;
#'       \frac{\sigma_A^2 + \sigma_{A\times T}^2}
#'            {\sigma_A^2 + \sigma_{A\times M}^2 +
#'             \sigma_{A\times T}^2 + S_B + \sigma_E^2}. }
#' There are special cases when there is no method factor (or a single level),
#' then \eqn{S_B=0} and \eqn{\sigma_{A\times M}^2=0}; if there is no
#' time factor (or a single level), then \eqn{\sigma_{A\times T}^2=0}.
#'
#' \strong{CIs / SEs (delta method for CCC).}
#' Let
#' \deqn{ \theta \;=\; \big(\sigma_A^2,\ \sigma_{A\times M}^2,\
#' \sigma_{A\times T}^2,\ \sigma_E^2,\ S_B\big)^\top, }
#' and write the concordance as
#' \deqn{ \mathrm{CCC}(\theta) \;=\; \frac{N}{D}
#'       \;=\; \frac{\sigma_A^2 + \sigma_{A\times T}^2}
#'                    {\sigma_A^2 + \sigma_{A\times M}^2 +
#'                    \sigma_{A\times T}^2 + S_B + \sigma_E^2}. }
#'
#' A first-order (large-sample) standard error follows from the delta method:
#' \deqn{ \mathrm{Var}\{\widehat{\mathrm{CCC}}\}
#'       \;\approx\; \nabla \mathrm{CCC}(\hat\theta)^\top\,
#'                   \mathrm{Var}(\hat\theta)\,
#'                   \nabla \mathrm{CCC}(\hat\theta), }
#' with gradient components (using \eqn{N} and \eqn{D} as above)
#' \deqn{ \frac{\partial\,\mathrm{CCC}}{\partial \sigma_A^2}
#'       \;=\; \frac{D - N}{D^2}
#'       \;=\; \frac{\sigma_{A\times M}^2 + S_B + \sigma_E^2}{D^2}, }
#' \deqn{ \frac{\partial\,\mathrm{CCC}}{\partial \sigma_{A\times M}^2}
#'       \;=\; -\,\frac{N}{D^2}, \qquad
#'        \frac{\partial\,\mathrm{CCC}}{\partial \sigma_{A\times T}^2}
#'       \;=\; \frac{D - N}{D^2}, }
#' \deqn{ \frac{\partial\,\mathrm{CCC}}{\partial \sigma_E^2}
#'       \;=\; -\,\frac{N}{D^2}, \qquad
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
#'        s_i \;=\; \frac{\|r_i - U_i b_i\|^2 +
#'        \mathrm{tr}(M_i^{-1}U_i^\top U_i)}{n_i}, }
#' where \eqn{b_i = M_i^{-1}(U_i^\top r_i/\sigma_E^2)} and
#' \eqn{M_i = G^{-1} + U_i^\top U_i/\sigma_E^2}.
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
#'
#' @section Notes on stability and performance:
#' All per-subject solves are \eqn{\,r\times r} with \eqn{r=1+nm+nt}, so cost
#' scales with the number of subjects and the fixed-effects dimension rather
#' than the total number of observations. Solvers use symmetric-PD paths with
#' a small diagonal ridge and pseudo-inverse fallback, which helps for
#' tiny/unbalanced subsets and near-boundary estimates. Very small samples or
#' extreme imbalance can still make \eqn{S_B} numerically delicate; negative
#' estimates are truncated to 0 by construction.
#'
#' @seealso \code{\link{build_L_Dm_cpp}} for constructing \eqn{L} and
#' \eqn{\mathrm{D_m}}; \code{\link{ccc_pairwise_u_stat}} for a U-statistic
#' alternative; and \pkg{cccrm} for a reference approach via \pkg{nlme}.
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
#' plot(ccc_rm1)
#'
#' # 95% CI container (limits currently NA by design)
#' ccc_rm2 <- ccc_lmm_reml(dat, ry = "y", rind = "id", rmet = "method", ci = TRUE)
#' print(ccc_rm2)
#' summary(ccc_rm2)
#' plot(ccc_rm2)
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
#' ## Three methods - pairwise CCCs
#' #--------------------------------------------------------------------
#' set.seed(2)
#' n_subj <- 40
#' id     <- factor(rep(seq_len(n_subj), times = 3))
#' method <- factor(rep(c("A","B","C"), each = n_subj))
#' sigA <- 1.2; sigE <- 0.6
#' u  <- rnorm(n_subj, 0, sqrt(sigA))
#' mu <- c(A = 0.00, B = 0.15, C = -0.10)
#' e  <- rnorm(3 * n_subj, 0, sqrt(sigE))
#' y  <- u[as.integer(id)] + unname(mu[method]) + e
#' dat3 <- data.frame(y, id, method)
#' ccc_lmm_reml(dat3, "y", "id", rmet = "method")
#'
#' # To get variance-components estimate per method combination, turn
#' # verbose to TRUE
#' ccc_lmm_reml(dat3, "y", "id", rmet = "method", verbose = TRUE)
#'
#' @author Thiago de Paula Oliveira
#' @importFrom stats as.formula model.matrix setNames qnorm
#' @export
ccc_lmm_reml <- function(data, ry, rind,
                         rmet = NULL, rtime = NULL, interaction = TRUE,
                         max_iter = 100, tol = 1e-6,
                         Dmat = NULL, ci = FALSE, conf_level = 0.95,
                         verbose = FALSE, digits = 4, use_message = TRUE) {

  ## ---- small helpers --------------------------------------------------------
  num_or_na <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (length(x) != 1 || !is.finite(x)) NA_real_ else x
  }
  compute_ci_from_se <- function(ccc, se, level) {
    if (!is.finite(ccc) || !is.finite(se)) return(c(NA_real_, NA_real_))
    z <- qnorm(1 - (1 - level)/2)
    c(max(0, min(1, ccc - z*se)), max(0, min(1, ccc + z*se)))
  }
  .vc_message <- function(ans, label, nm, nt, conf_level,
                          digits = 4, use_message = TRUE) {
    fmt <- function(x) if (is.na(x)) "NA" else formatC(x, format = "f", digits = digits)
    out <- c(
      sprintf("---- matrixCorr::ccc_lmm_reml - variance-components (%s) ----", label),
      sprintf("Design: methods nm = %d, times nt = %d", nm, nt),
      "Estimates:",
      sprintf("  sigma_A^2 (subject)            : %s", fmt(num_or_na(ans[["sigma2_subject"]]))),
      sprintf("  sigma_A_M^2 (subject x method) : %s", fmt(num_or_na(ans[["sigma2_subject_method"]]))),
      sprintf("  sigma_A_T^2 (subject x time)   : %s", fmt(num_or_na(ans[["sigma2_subject_time"]]))),
      sprintf("  sigma_E^2 (error)              : %s", fmt(num_or_na(ans[["sigma2_error"]]))),
      sprintf("  S_B (fixed-effect dispersion)  : %s", fmt(num_or_na(ans[["SB"]]))),
      sprintf("  SE(CCC)                        : %s", fmt(num_or_na(ans[["se_ccc"]]))),
      "--------------------------------------------------------------------------"
    )
    if (use_message) message(paste(out, collapse = "\n")) else cat(paste(out, collapse = "\n"), "\n")
  }

  ## ---- coerce once ----------------------------------------------------------
  df <- as.data.frame(data)
  df[[ry]]   <- as.numeric(df[[ry]])
  df[[rind]] <- factor(df[[rind]])
  if (!is.null(rmet))  df[[rmet]]  <- factor(df[[rmet]])
  if (!is.null(rtime)) df[[rtime]] <- factor(df[[rtime]])

  # keep the overall time levels for Dmat subsetting in pairwise path
  all_time_lvls <- if (!is.null(rtime)) levels(df[[rtime]]) else character(0)

  ## ---- fixed-effects formula ------------------------------------------------
  rhs <- "1"
  if (!is.null(rmet))  rhs <- paste(rhs, "+", rmet)
  if (!is.null(rtime)) rhs <- paste(rhs, "+", rtime)
  if (!is.null(rmet) && !is.null(rtime) && isTRUE(interaction))
    rhs <- paste(rhs, "+", paste0(rmet, ":", rtime))
  fml <- as.formula(paste("~", rhs))

  ## ---- Case 1: no (or single) method -> overall CCC (1x1) ------------------
  if (is.null(rmet) || nlevels(df[[rmet]]) < 2L) {
    X <- model.matrix(fml, data = df)

    # For overall case, pass Dmat if time has >=2 levels (else NULL)
    Dsub <- if (!is.null(Dmat) && !is.null(rtime) && nlevels(df[[rtime]]) >= 2L) {
      as.matrix(Dmat)
    } else NULL

    Laux <- build_L_Dm_cpp(
      colnames_X      = colnames(X),
      rmet_name       = if (is.null(rmet)) NULL else rmet,
      rtime_name      = if (is.null(rtime)) NULL else rtime,
      method_levels   = if (is.null(rmet)) character(0) else levels(df[[rmet]]),
      time_levels     = if (is.null(rtime)) character(0) else levels(df[[rtime]]),
      has_interaction = interaction,
      Dmat_global     = Dsub
    )

    # integer codes for method/time (0-length if absent or single level)
    method_int <- if (!is.null(rmet)  && nlevels(df[[rmet]])  >= 2L) as.integer(df[[rmet]])  else integer(0)
    time_int   <- if (!is.null(rtime) && nlevels(df[[rtime]]) >= 2L) as.integer(df[[rtime]]) else integer(0)

    ans <- tryCatch(
      ccc_vc_cpp(
        Xr = unname(X),
        yr = df[[ry]],
        subject = as.integer(df[[rind]]),
        method  = method_int,
        time    = time_int,
        nm = Laux$nm, nt = Laux$nt,
        max_iter = max_iter, tol = tol,
        conf_level = conf_level,
        Lr   = if (is.null(Laux$L))    NULL else unname(Laux$L),
        auxDr = if (is.null(Laux$Dm)) NULL else unname(Laux$Dm)
      ),
      error = function(e) {
        stop("ccc_vc_cpp failed on this dataset (near-singular tiny data): ",
             conditionMessage(e), call. = FALSE)
      }
    )

    if (isTRUE(verbose)) {
      .vc_message(ans, label = "Overall", nm = Laux$nm, nt = Laux$nt,
                  conf_level = conf_level, digits = digits, use_message = use_message)
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
      attr(out, "method")     <- "Variance Components REML"
      attr(out, "description")<- "Lin's CCC from random-effects LMM"
      attr(out, "package")    <- "matrixCorr"
      attr(out, "conf.level") <- conf_level
      class(out) <- c("ccc", "ccc_ci")
      return(out)
    } else {
      est <- est_mat
      attr(est, "method")      <- "Variance Components REML"
      attr(est, "description") <- "Lin's CCC from random-effects LMM"
      attr(est, "package")     <- "matrixCorr"
      class(est) <- c("ccc", "matrix")
      return(est)
    }
  }

  ## ---- Case 2: pairwise by method (L >= 2) ---------------------------------
  df[[rmet]] <- droplevels(df[[rmet]])
  method_levels <- levels(df[[rmet]])
  Lm <- length(method_levels)

  est_mat <- matrix(1,  Lm, Lm, dimnames = list(method_levels, method_levels))
  if (isTRUE(ci)) {
    lwr_mat <- matrix(NA_real_, Lm, Lm, dimnames = list(method_levels, method_levels))
    upr_mat <- matrix(NA_real_, Lm, Lm, dimnames = list(method_levels, method_levels))
  }

  for (i in 1:(Lm - 1L)) {
    for (j in (i + 1L):Lm) {
      m1 <- method_levels[i]; m2 <- method_levels[j]

      idx <- which(df[[rmet]] %in% c(m1, m2))
      # cheap “sub” views
      subj_int   <- as.integer(df[[rind]][idx])
      y_sub      <- df[[ry]][idx]
      met_fac    <- droplevels(df[[rmet]][idx])        # ensures exactly 2 levels
      time_fac   <- if (!is.null(rtime)) droplevels(df[[rtime]][idx]) else NULL

      Xp <- model.matrix(fml, data = df[idx, , drop = FALSE])

      # subset Dmat to present time levels in this pair (if any)
      lev_time_sub <- if (!is.null(time_fac)) levels(time_fac) else character(0)
      Dsub <- if (!is.null(Dmat) && length(lev_time_sub) >= 2L) {
        pos <- match(lev_time_sub, all_time_lvls)
        as.matrix(Dmat[pos, pos, drop = FALSE])
      } else NULL

      Laux <- build_L_Dm_cpp(
        colnames_X      = colnames(Xp),
        rmet_name       = rmet,
        rtime_name      = if (is.null(rtime)) NULL else rtime,
        method_levels   = levels(met_fac),
        time_levels     = lev_time_sub,
        has_interaction = interaction,
        Dmat_global     = Dsub
      )

      method_int <- if (nlevels(met_fac)  >= 2L) as.integer(met_fac)  else integer(0)
      time_int   <- if (!is.null(time_fac) && nlevels(time_fac) >= 2L) as.integer(time_fac) else integer(0)

      ans <- tryCatch(
        ccc_vc_cpp(
          Xr = unname(Xp),
          yr = y_sub,
          subject = subj_int,
          method  = method_int,
          time    = time_int,
          nm = Laux$nm, nt = Laux$nt,
          max_iter = max_iter, tol = tol,
          conf_level = conf_level,
          Lr   = if (is.null(Laux$L))    NULL else unname(Laux$L),
          auxDr = if (is.null(Laux$Dm)) NULL else unname(Laux$Dm)
        ),
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
                    use_message = use_message)
      }

      val <- if (is.null(ans)) NA_real_ else unname(ans$ccc)
      est_mat[i, j] <- est_mat[j, i] <- val

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
    class(out) <- c("ccc", "ccc_ci")
    return(out)
  } else {
    est <- est_mat
    attr(est, "method")      <- "Variance Components REML - pairwise"
    attr(est, "description") <- "Lin's CCC per method pair from random-effects LMM"
    attr(est, "package")     <- "matrixCorr"
    class(est) <- c("ccc", "matrix")
    return(est)
  }
}
