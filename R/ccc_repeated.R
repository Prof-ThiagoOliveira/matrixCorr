#' @title Repeated-Measures Lin's Concordance Correlation Coefficient (CCC)
#'
#' @description
#' Computes all pairwise Lin's Concordance Correlation Coefficients (CCC)
#' across multiple methods (L ≥ 2) for repeated-measures data. Each subject
#' must be measured by all methods across the same set of time points or
#' replicates.
#'
#' CCC measures both accuracy (how close measurements are to the line of
#' equality) and precision (Pearson correlation). Confidence intervals are
#' optionally computed using a U-statistics-based estimator with Fisher's Z
#' transformation
#'
#' @param data A data frame containing the repeated-measures dataset.
#' @param ry Character. Name of the numeric outcome column.
#' @param rmet Character. Name of the method column (factor with L ≥ 2 levels).
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
#' @param n_threads Integer (≥ 1). Number of OpenMP threads to use for computation.
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
#' This function computes pairwise Lin’s Concordance Correlation Coefficient
#' (CCC) between methods in a repeated-measures design using a
#' U-statistics-based nonparametric estimator proposed by
#' Carrasco et al. (2013). It is computationally efficient and robust,
#' particularly for large-scale or balanced longitudinal designs.
#'
#' Lin’s CCC is defined as
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
#' \emph{Biometrics}, 45: 255–268.
#'
#' Lin L (2000). A note on the concordance correlation coefficient.
#' \emph{Biometrics}, 56: 324–325.
#'
#' Carrasco JL, Jover L (2003). Estimating the concordance correlation coefficient:
#' a new approach. \emph{Computational Statistics & Data Analysis}, 47(4): 519–539.
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
#' plot(ccc1)
#'
#' # With confidence intervals
#' ccc2 <- ccc_pairwise_u_stat(df, ry = "y", rmet = "method", rtime = "time", ci = TRUE)
#' print(ccc2)
#' plot(ccc2)
#'
#' #------------------------------------------------------------------------
#' # Choosing delta based on distance sensitivity
#' #------------------------------------------------------------------------
#' # Absolute distance (L1 norm) – robust
#' ccc_pairwise_u_stat(df, ry = "y", rmet = "method", rtime = "time", delta = 1)
#'
#' # Squared distance (L2 norm) – amplifies large deviations
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
#' Compute Lin’s Concordance Correlation Coefficient (CCC) from a linear
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
#' @param max_iter Integer. Maximum iterations for variance–component updates
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
#' \strong{EM-style variance–component updates.} With current \eqn{\hat\beta},
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
#' each time level) and an optional time-weighting matrix \eqn{\mathrm{Dm}}:
#' \deqn{ S_B \;=\;
#'  \frac{\big(L^\top \hat\beta\big)^\top\,\mathrm{Dm}\,
#'  \big(L^\top \hat\beta\big)
#'        \;-\; \mathrm{tr}\!\Big(\big(L\,\mathrm{Dm}\,L^\top\big)\,
#'        \mathrm{Var}(\hat\beta)\Big)}
#'       {\,nm\,(nm-1)\,\max(nt,1)\,}, }
#' truncated at 0. The helper \code{\link{build_L_Dm}} constructs \eqn{L} so it
#' aligns exactly with the columns of \eqn{X=\mathrm{model.matrix}(\cdot)}
#' passed to 'C++'.
#' For exactly two methods (\eqn{nm=2}), a fast path builds \eqn{L} directly
#' from the design’s column names, where, with interaction, the per-time
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
#' \strong{CIs / SEs.}
#'
#' @return
#' \itemize{
#'   \item If \code{rmet} is \code{NULL} or has a single level, an object of
#'     class \code{c("ccc","ccc_ci")} (when \code{ci=TRUE}) or
#'     \code{c("ccc","matrix")} with a \eqn{1\times 1} matrix containing the
#'     overall CCC estimate.
#'   \item If \code{rmet} has \eqn{L\ge 2} levels, a symmetric \eqn{L\times L}
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
#' @seealso \code{\link{build_L_Dm}} for constructing \eqn{L} and
#' \eqn{\mathrm{Dm}}; \code{\link{ccc_pairwise_u_stat}} for a U-statistic
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
#' ccc_lmm_reml(dat, ry = "y", rind = "id", rmet = "method")
#'
#' # 95% CI container (limits currently NA by design)
#' ccc_lmm_reml(dat, ry = "y", rind = "id", rmet = "method", ci = TRUE)
#'
#' #--------------------------------------------------------------------
#' ## Two methods × time (balanced 2x2), with and without interaction
#' #--------------------------------------------------------------------
#' dat$time <- factor(rep(rep(c("t1","t2"), each = n_subj/2), times = 2))
#' ccc_lmm_reml(dat, "y", "id", rmet = "method", rtime = "time",
#'              interaction = FALSE)
#' ccc_lmm_reml(dat, "y", "id", rmet = "method", rtime = "time",
#'              interaction = TRUE)
#'
#' #--------------------------------------------------------------------
#' ## Three methods — pairwise CCCs
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
#' @author Thiago de Paula Oliveira
#' @importFrom stats as.formula model.matrix setNames
#' @export
ccc_lmm_reml <- function(data, ry, rind,
                         rmet = NULL, rtime = NULL, interaction = TRUE,
                         max_iter = 100, tol = 1e-6,
                         Dmat = NULL, ci = FALSE, conf_level = 0.95) {

  df <- as.data.frame(data)
  df[[ry]]   <- as.numeric(df[[ry]])
  df[[rind]] <- factor(df[[rind]])
  if (!is.null(rmet))  df[[rmet]]  <- factor(df[[rmet]])
  if (!is.null(rtime)) df[[rtime]] <- factor(df[[rtime]])

  #-----------------------------------------------------------------------------
  ## ---- fixed-effects formula ----
  #-----------------------------------------------------------------------------
  rhs <- "1"
  if (!is.null(rmet))  rhs <- paste(rhs, "+", rmet)
  if (!is.null(rtime)) rhs <- paste(rhs, "+", rtime)
  if (!is.null(rmet) && !is.null(rtime) && isTRUE(interaction))
    rhs <- paste(rhs, "+", paste0(rmet, ":", rtime))
  fml <- as.formula(paste("~", rhs))

  #-----------------------------------------------------------------------------
  ## Case 1: no (or single) method -> overall CCC (1x1)
  #-----------------------------------------------------------------------------
  if (is.null(rmet) || nlevels(df[[rmet]]) < 2L) {
    X <- model.matrix(fml, data = df)
    Laux <- build_L_Dm(df_sub = df, fml_sub = fml,
                       rmet = rmet, rtime = rtime,
                       # nm will be 0
                       Dmat_global = Dmat, has_interaction = interaction)

    # integer codes for method/time (0-length if absent or single level)
    method_int <- if (!is.null(rmet)  && nlevels(df[[rmet]])  >= 2L) {
      as.integer(df[[rmet]])
    } else integer(0)
    time_int   <- if (!is.null(rtime) && nlevels(df[[rtime]]) >= 2L) {
      as.integer(df[[rtime]])
    } else integer(0)

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
        Dmr = if (is.null(Laux$Dm)) NULL else unname(Laux$Dm)
      ),
      error = function(e) {
        stop("ccc_vc_cpp failed on this dataset (near-singular tiny data): ",
             conditionMessage(e), call. = FALSE)
      }
    )

    lab <- "Overall"
    est_mat <- matrix(unname(ans$ccc), 1, 1, dimnames = list(lab, lab))

    if (isTRUE(ci)) {
      out <- list(
        est    = est_mat,
        lwr.ci = matrix(if (!is.null(ans$lwr)) {
          unname(ans$lwr)
        } else NA_real_, 1, 1, dimnames = list(lab, lab)),
        upr.ci = matrix(if (!is.null(ans$upr)) {
          unname(ans$upr)
        } else NA_real_, 1, 1, dimnames = list(lab, lab))
      )
      attr(out, "method")     <-
        "Variance Components REML (LMM, cccrm-style random effects)"
      attr(out, "description")<- "Lin's CCC from random-effects LMM"
      attr(out, "package")    <- "matrixCorr"
      attr(out, "conf.level") <- conf_level
      class(out) <- c("ccc", "ccc_ci")
      return(out)
    } else {
      est <- est_mat
      attr(est, "method")     <-
        "Variance Components REML (LMM, cccrm-style random effects)"
      attr(est, "description")<- "Lin's CCC from random-effects LMM"
      attr(est, "package")    <- "matrixCorr"
      class(est) <- c("ccc", "matrix")
      return(est)
    }
  }

  #-----------------------------------------------------------------------------
  ## Case 2: pairwise by method (L >= 2)
  #-----------------------------------------------------------------------------
  df[[rmet]] <- droplevels(df[[rmet]])
  method_levels <- levels(df[[rmet]])
  Lm <- length(method_levels)

  est_mat <- matrix(1,  Lm, Lm, dimnames = list(method_levels, method_levels))
  if (isTRUE(ci)) {
    lwr_mat <- matrix(NA_real_, Lm, Lm,
                      dimnames = list(method_levels, method_levels))
    upr_mat <- matrix(NA_real_, Lm, Lm,
                      dimnames = list(method_levels, method_levels))
  }

  for (i in 1:(Lm - 1L)) {
    for (j in (i + 1L):Lm) {
      m1 <- method_levels[i]; m2 <- method_levels[j]
      sub <- df[df[[rmet]] %in% c(m1, m2), , drop = FALSE]
      sub[[rmet]] <- droplevels(sub[[rmet]])
      if (!is.null(rtime)) sub[[rtime]] <- droplevels(sub[[rtime]])

      Xp   <- model.matrix(fml, data = sub)

      Laux <- build_L_Dm(df_sub = sub, fml_sub = fml,
                           rmet = rmet, rtime = rtime,
                           Dmat_global = Dmat, has_interaction = interaction)

      method_int <- if (nlevels(sub[[rmet]])  >= 2L) {
        as.integer(sub[[rmet]])
      } else integer(0)
      time_int   <-
        if (!is.null(rtime) && nlevels(sub[[rtime]]) >= 2L) {
          as.integer(sub[[rtime]])
        } else integer(0)

      ans <- tryCatch(
        ccc_vc_cpp(
          Xr = unname(Xp),
          yr = sub[[ry]],
          subject = as.integer(sub[[rind]]),
          method  = method_int,
          time    = time_int,
          nm = Laux$nm, nt = Laux$nt,
          max_iter = max_iter, tol = tol,
          conf_level = conf_level,
          Lr   = if (is.null(Laux$L))    NULL else unname(Laux$L),
          Dmr = if (is.null(Laux$Dm)) NULL else unname(Laux$Dm)
        ),
        error = function(e) {
          warning(sprintf("ccc_vc_cpp failed for pair (%s, %s): %s",
                          m1, m2, conditionMessage(e)))
          NULL
        }
      )

      val <- if (is.null(ans)) NA_real_ else unname(ans$ccc)
      est_mat[i, j] <- est_mat[j, i] <- val

      if (isTRUE(ci)) {
        lwr_mat[i, j] <-
          lwr_mat[j, i] <-
          if (!is.null(ans) && !is.null(ans$lwr)) unname(ans$lwr) else NA_real_
        upr_mat[i, j] <-
          upr_mat[j, i] <-
          if (!is.null(ans) && !is.null(ans$upr)) unname(ans$upr) else NA_real_
      }
    }
  }

  diag(est_mat) <- 1
  if (isTRUE(ci)) {
    diag(lwr_mat) <- NA_real_
    diag(upr_mat) <- NA_real_
    out <- list(est = est_mat, lwr.ci = lwr_mat, upr.ci = upr_mat)
    attr(out, "method")     <-
      "Variance Components REML (LMM, cccrm-style random effects) — pairwise"
    attr(out, "description")<-
      "Lin's CCC per method pair from random-effects LMM"
    attr(out, "package")    <- "matrixCorr"
    attr(out, "conf.level") <- conf_level
    class(out) <- c("ccc", "ccc_ci")
    return(out)
  } else {
    est <- est_mat
    attr(est, "method")     <-
      "Variance Components REML (LMM, cccrm-style random effects) — pairwise"
    attr(est, "description")<-
      "Lin's CCC per method pair from random-effects LMM"
    attr(est, "package")    <- "matrixCorr"
    class(est) <- c("ccc", "matrix")
    return(est)
  }
}


#' Build per-time method contrast matrix L and time-weighting D
#'
#' @description
#' Helper that constructs the contrast matrix \code{L} (aligned to the columns
#' of the *actual* design matrix produced by
#' \code{model.matrix(fml_sub, df_sub)}) and the time-weighting matrix
#' \code{Dm} used to compute the fixed-effect dispersion term \eqn{S_B} in the
#' CCC. For the common case of exactly two methods (\eqn{nm = 2}), \code{L} is
#' built directly from column names, guaranteeing perfect alignment and
#' avoiding grid/contrast mismatches. For \eqn{nm \ge 3},
#' a grid-based construction is used but re-ordered to match the real design.
#'
#' @param df_sub A \code{data.frame} subset corresponding to the data being
#' analyzed (e.g., overall or a pair of methods). \code{df_sub[[rmet]]} and
#' \code{df_sub[[rtime]]} must already be factors with the intended levels
#' (use \code{droplevels()} upstream as needed).
#' @param fml_sub A right-hand-side formula (e.g. \code{~ 1 + rmet + rtime}
#' or \code{~ 1 + rmet + rtime + rmet:rtime}) that matches the fixed-effects
#' part of the model passed to 'C++'.
#' @param rmet Character scalar. Column name of the method factor in
#' \code{df_sub} (or \code{NULL} if no method factor).
#' @param rtime Character scalar. Column name of the time factor in
#' \code{df_sub} (or \code{NULL} if no time factor).
#' @param Dmat_global Optional numeric matrix of size \eqn{n_t \times n_t}
#' (where \eqn{n_t} is the number of time levels in the *original* data)
#' providing time weights for \eqn{S_B}. If \code{NULL}, the identity is used.
#' When working on a subset, this function subsets \code{Dmat_global} to the
#' time levels present in \code{df_sub}.
#' @param has_interaction Logical. Whether the fixed-effects design includes the
#' \code{rmet:rtime} interaction term. This determines the form of the method
#' contrasts when \eqn{nm = 2}: with interaction, \eqn{\Delta(t_j) =
#' \beta_{\text{met2}} + \beta_{\text{met2:time}_j}}; without interaction,
#' \eqn{\Delta(t_j) = \beta_{\text{met2}}} for all \eqn{j}.
#'
#' @details
#' Let \eqn{X = \texttt{model.matrix}(fml\_sub, df\_sub)} with \eqn{p} columns.
#' This function returns:
#'
#' \itemize{
#'   \item \strong{Two-method fast path} (\eqn{nm = 2}):
#'     \itemize{
#'       \item If there is no time factor (\eqn{nt = 0}), \code{L} is a
#'         \eqn{p \times 1} matrix that picks the \code{met2} column in \eqn{X}
#'         (i.e., the treatment-coded difference
#'         \eqn{\beta_{\text{met2}} = \mu_2 - \mu_1}).
#'       \item If \eqn{nt \ge 1}, \code{L} is \eqn{p \times nt} with columns
#'         representing the method difference at each time level. With
#'         interaction, column \eqn{j} encodes \eqn{\beta_{\text{met2}} +
#'         \beta_{\text{met2:time}_j}} for \eqn{j > 1} and
#'         \eqn{\beta_{\text{met2}}} for the baseline;
#'         without interaction, all columns equal \eqn{\beta_{\text{met2}}}.
#'       \item \code{Dm} is the time-weighting matrix of size \eqn{nt \times nt}
#'         (\code{diag(nt)} if \code{Dmat_global} is \code{NULL}).
#'     }
#'   \item \strong{General path} (\eqn{nm \ge 3}) represents
#'     a method–time grid, where pairwise method differences are encoded
#'     within each time level, yielding \code{L} of size
#'     \eqn{p \times (nd \cdot \max(nt,1))}, with
#'     \eqn{nd = nm(nm - 1)/2}. The grid design matrix is
#'     column-reordered to match the actual \eqn{X} columns (critical to keep
#'     \eqn{S_B} correct). \code{Dm} is \eqn{\mathrm{kronecker}(D, I_{nd})}
#'     with \code{D} the time-weighting matrix.
#' }
#'
#' All constructions are based on the \emph{actual} column names of \eqn{X}
#' to ensure \code{L} aligns perfectly with the design passed to the 'C++' core.
#' This makes the \eqn{S_B} term numerically stable, especially in small or
#' unbalanced designs.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{L}: numeric matrix of contrasts, dimension
#'         \eqn{p \times q}, where \eqn{q = nt} if \eqn{nm = 2} and
#'         \eqn{q = nd \cdot \max(nt,1)} otherwise.
#'   \item \code{Dm}: numeric time-weighting matrix of size \eqn{q \times q}.
#'   \item \code{nm}: effective number of method levels used (\eqn{\ge 0}).
#'   \item \code{nt}: effective number of time levels used (\eqn{\ge 0}).
#' }
#'
#' @section Errors and checks:
#' \itemize{
#'   \item If the required method or interaction columns cannot be found in
#'         \code{colnames(model.matrix(fml_sub, df_sub))}, an error is thrown.
#'   \item If \code{Dmat_global} is supplied but its dimension does not match
#'         the total number of time levels in the original data, an error is
#'         thrown.
#' }
#'
#' @seealso \code{\link{ccc_lmm_reml}} for the high-level wrapper that calls
#' this helper, and the 'C++' core that consumes \code{L}/\code{Dm}.
#' @author Thiago de Paula Oliveira
#' @keywords internal
build_L_Dm <- function(df_sub, fml_sub, rmet, rtime, Dmat_global,
                       has_interaction) {
  nm_levels <- if (is.null(rmet)) 0L else nlevels(df_sub[[rmet]])
  nt_levels <- if (is.null(rtime)) 0L else nlevels(df_sub[[rtime]])
  nm_eff <- if (nm_levels >= 2L) nm_levels else 0L
  nt_eff <- if (nt_levels >= 2L) nt_levels else 0L

  #-----------------------------------------------------------------------------
  # If no method contrasts are needed
  #-----------------------------------------------------------------------------
  if (nm_eff == 0L) return(list(L = NULL, Dm = NULL, nm = nm_eff, nt = nt_eff))

  #-----------------------------------------------------------------------------
  # Use the actual design you'll pass to 'C++' (column names must match)
  #-----------------------------------------------------------------------------
  Xsub <- model.matrix(fml_sub, data = df_sub)
  cn   <- colnames(Xsub)
  p    <- ncol(Xsub)

  #-----------------------------------------------------------------------------
  # ---- Fast, exact path for two methods ----
  #-----------------------------------------------------------------------------
  if (nm_eff == 2L) {
    lev_met  <- levels(df_sub[[rmet]])
    met2_name <- paste0(rmet, lev_met[2L])
    met2_idx  <- match(met2_name, cn)
    if (is.na(met2_idx)) stop("Could not find column for '", met2_name,
                              "' in model matrix.")

    # Time weighting matrix
    if (nt_eff == 0L) {
      Dsub <- matrix(1, 1, 1)
    } else if (is.null(Dmat_global)) {
      Dsub <- diag(nt_eff)
    } else {
      lev_all <- levels(df_sub[[rtime]])
      if (!all(dim(Dmat_global) == c(length(lev_all), length(lev_all))))
        stop("Dmat dimension mismatch: expected ", length(lev_all), " x ",
             length(lev_all))
      idx_t <- match(levels(df_sub[[rtime]]), lev_all)
      Dsub  <- as.matrix(Dmat_global[idx_t, idx_t, drop = FALSE])
    }

    # Build L directly against Xsub
    if (nt_eff == 0L) {
      # Methods-only: one contrast (met2 - met1)
      L <- matrix(0, nrow = p, ncol = 1); colnames(L) <- "method"
      L[met2_idx, 1L] <- 1
    } else {
      # Methods x time uses per-time diffs
      lev_time <- levels(df_sub[[rtime]])
      L <- matrix(0, nrow = p, ncol = nt_eff); colnames(L) <- lev_time

      # baseline time (first level): delta = beta_met2
      L[met2_idx, 1L] <- 1

      if (isTRUE(has_interaction)) {
        # other times uses delta(t_j) = beta_met2 + beta_met2:time_j
        for (j in 2:nt_eff) {
          L[met2_idx, j] <- 1
          inter_name1 <- paste0(met2_name, ":", rtime, lev_time[j])
          inter_name2 <- paste0(rtime, lev_time[j], ":", met2_name) # fallback
          inter_idx   <- match(inter_name1, cn, nomatch = NA_integer_)
          if (is.na(inter_idx)) inter_idx <- match(inter_name2, cn,
                                                   nomatch = NA_integer_)
          if (is.na(inter_idx))
            stop("Could not find interaction column for '", met2_name, ":",
                 rtime, lev_time[j], "' in model matrix.")
          L[inter_idx, j] <- 1
        }
      } else {
        # No interaction: method difference is the same at all times
        for (j in 2:nt_eff) L[met2_idx, j] <- 1
      }
    }

    return(list(L = L, Dm = Dsub, nm = nm_eff, nt = nt_eff))
  }

  #-----------------------------------------------------------------------------
  # ---- General path (>=3 methods) ----
  #-----------------------------------------------------------------------------
  if (nt_eff == 0L) {
    grid <- data.frame(setNames(list(levels(df_sub[[rmet]])), rmet),
                       stringsAsFactors = FALSE)
  } else {
    grid <- do.call(expand.grid, c(
      setNames(list(levels(df_sub[[rtime]])), rtime),  # time outer
      setNames(list(levels(df_sub[[rmet]])),  rmet),   # method inner
      list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    ))
  }
  grid[[rmet]]  <- factor(grid[[rmet]],  levels = levels(df_sub[[rmet]]))
  if (nt_eff > 0L)
    grid[[rtime]] <- factor(grid[[rtime]], levels = levels(df_sub[[rtime]]))

  gridX <- model.matrix(fml_sub, data = grid)
  #-----------------------------------------------------------------------------
  # Reorder to match the actual Xsub columns (CRITICAL)
  #-----------------------------------------------------------------------------
  if (!identical(colnames(gridX), cn)) {
    gridX <- gridX[, cn, drop = FALSE]
  }

  nd <- nm_eff * (nm_eff - 1L) / 2L
  Diffs <- matrix(0, nrow = nd * max(nt_eff, 1L),
                  ncol = nm_eff * max(nt_eff, 1L))
  row <- 0L
  idx <- function(m, t) (t - 1L) * nm_eff + m
  for (t in seq_len(max(nt_eff, 1L))) {
    for (i in 1:(nm_eff - 1L)) for (j in (i + 1L):nm_eff) {
      row <- row + 1L
      Diffs[row, idx(i, t)] <-  1
      Diffs[row, idx(j, t)] <- -1
    }
  }
  L <- t(gridX) %*% t(Diffs)

  if (nt_eff == 0L) {
    Dsub <- matrix(1, 1, 1)
  } else if (is.null(Dmat_global)) {
    Dsub <- diag(nt_levels)
  } else {
    lev_all <- levels(df_sub[[rtime]])
    idx_t <- match(levels(df_sub[[rtime]]), lev_all)
    Dsub  <- as.matrix(Dmat_global[idx_t, idx_t, drop = FALSE])
  }
  Dm <- kronecker(Dsub, diag(nd))
  list(L = L, Dm = Dm, nm = nm_eff, nt = nt_eff)
}
