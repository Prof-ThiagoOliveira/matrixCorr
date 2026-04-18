#' @title Bland-Altman for repeated measurements
#' @description
#' Repeated-measures Bland-Altman (BA) analysis for method comparison based on a
#' mixed-effects model fitted to **subject-time matched paired differences**.
#' The fitted model includes a subject-specific random intercept and, optionally,
#' an AR(1) residual correlation structure within subject.
#'
#' The function accepts either exactly two methods or \eqn{\ge 3} methods.
#' With exactly two methods it returns a single fitted BA object. With
#' \eqn{\ge 3} methods it fits the same model to every unordered method pair and
#' returns pairwise matrices of results.
#'
#' **Required variables**
#' \itemize{
#'   \item \code{response}: numeric measurements.
#'   \item \code{subject}: subject identifier.
#'   \item \code{method}: method label with at least two distinct levels.
#'   \item \code{time}: replicate/time key used to form within-subject pairs.
#' }
#'
#' For any analysed pair of methods, only records where both methods are present
#' for the same \code{subject} and integer-coerced \code{time} contribute to the
#' fit. Rows with missing values in any required field are excluded for that
#' analysed pair.
#'
#' @param data Optional data frame-like object. If supplied, `response`,
#'   `subject`, `method`, and `time` may be column names. Objects not already
#'   inheriting from `data.frame` are first coerced with `as.data.frame()`.
#' @param response Numeric response vector, or a single character string naming
#'   the response column in `data`.
#' @param subject Subject identifier vector (integer, numeric, or factor), or a
#'   single character string naming the subject column in `data`.
#' @param method Method label vector (character, factor, integer, or numeric), or
#'   a single character string naming the method column in `data`. At least two
#'   distinct method levels are required.
#' @param time Replicate/time index vector (integer or numeric), or a single
#'   character string naming the time column in `data`. Values are coerced to
#'   integer before pairing and before AR(1) contiguity checks.
#' @param loa_multiplier Positive scalar giving the SD multiplier used to form
#'   the limits of agreement. Default is `1.96`. In the exported `ba_rm()`
#'   interface this default is fixed and does not depend on `conf_level`.
#' @param conf_level Confidence level for Wald confidence intervals for the
#'   reported bias and both LoA endpoints. Must lie in `(0, 1)`. Default `0.95`.
#' @param n_threads Integer \eqn{\geq 1}. Number of OpenMP threads. Defaults to
#'   \code{getOption("matrixCorr.threads", 1L)}.
#' @param include_slope Logical. If `TRUE`, the model includes the paired mean as
#'   a fixed effect and estimates a proportional-bias slope. The reported BA
#'   centre remains the fitted mean difference at the centred reference paired
#'   mean used internally by the backend; the returned LoA remain horizontal
#'   bands and are not regression-adjusted curves.
#' @param use_ar1 Logical. If `TRUE`, request an AR(1) residual structure within
#'   subject over contiguous integer time blocks.
#' @param ar1_rho Optional AR(1) parameter. Must satisfy `abs(ar1_rho) < 0.999`
#'   when supplied. If `NA` and `use_ar1 = TRUE`, the backend estimates `rho`
#'   separately for each analysed pair.
#' @param max_iter Maximum number of EM/GLS iterations used by the backend.
#' @param tol Convergence tolerance for the backend EM/GLS iterations.
#' @param verbose Logical. If `TRUE`, print progress for the pairwise
#'   \eqn{\ge 3}-method path. It has no material effect in the exactly-two-method
#'   path.
#'
#' @return
#' Either a \code{"ba_repeated"} object (exactly two methods) or a
#' \code{"ba_repeated_matrix"} object (\eqn{\ge 3} methods).
#'
#' \strong{If exactly two methods are supplied}, the returned
#' \code{"ba_repeated"} object is a list with components:
#' \itemize{
#'   \item \code{means}: numeric vector of paired means
#'   \eqn{(y_1 + y_2)/2} used for plotting helpers.
#'   \item \code{diffs}: numeric vector of paired differences
#'   \eqn{y_2 - y_1} used for plotting helpers.
#'   \item \code{n_obs}: integer number of complete subject-time pairs used.
#'   \item \code{based.on}: compatibility alias for \code{n_obs}.
#'   \item \code{mean.diffs}: scalar model-based BA centre. When
#'   `include_slope = FALSE`, this is the fitted intercept of the paired-
#'   difference model. When `include_slope = TRUE`, this is the fitted mean
#'   difference at the centred reference paired mean used internally by the
#'   backend.
#'   \item \code{lower.limit}, \code{upper.limit}: scalar limits of agreement,
#'   computed as
#'   \eqn{\mu_0 \pm \texttt{loa\_multiplier}\sqrt{\sigma_u^2 + \sigma_e^2}}.
#'   \item \code{lines}: named numeric vector with entries `lower`, `mean`, and
#'   `upper`.
#'   \item \code{CI.lines}: named numeric vector containing Wald confidence
#'   interval bounds for the bias and both LoA endpoints:
#'   `mean.diff.ci.lower`, `mean.diff.ci.upper`,
#'   `lower.limit.ci.lower`, `lower.limit.ci.upper`,
#'   `upper.limit.ci.lower`, `upper.limit.ci.upper`.
#'   \item \code{loa_multiplier}: scalar LoA multiplier actually used.
#'   \item \code{critical.diff}: scalar LoA half-width
#'   \eqn{\texttt{loa\_multiplier} \times \texttt{sd\_loa}}.
#'   \item \code{include_slope}: logical, copied from the call.
#'   \item \code{beta_slope}: proportional-bias slope on the original paired-mean
#'   scale when `include_slope = TRUE`; otherwise `NA`.
#'   \item \code{sigma2_subject}: estimated variance of the subject-level random
#'   intercept on paired differences.
#'   \item \code{sigma2_resid}: estimated residual variance on paired differences.
#'   \item \code{use_ar1}: logical, copied from the call.
#'   \item \code{residual_model}: either `"ar1"` or `"iid"`, indicating the final
#'   residual structure actually used.
#'   \item \code{ar1_rho}: AR(1) correlation actually used in the final fit when
#'   `residual_model == "ar1"`; otherwise `NA`.
#'   \item \code{ar1_estimated}: logical indicating whether `ar1_rho` was
#'   estimated internally (`TRUE`) or supplied by the user (`FALSE`) when the
#'   final residual model is AR(1); otherwise `NA`.
#' }
#'
#' The confidence level is stored as `attr(x, "conf.level")`.
#'
#' \strong{If \eqn{\ge 3} methods are supplied}, the returned
#' \code{"ba_repeated_matrix"} object is a list with components:
#' \itemize{
#'   \item \code{bias}: numeric \eqn{m \times m} matrix of model-based BA centres.
#'   For indices \eqn{(j, k)} with \eqn{j < k}, \code{bias[j, k]} estimates
#'   \eqn{\mathrm{method}_k - \mathrm{method}_j}. Thus the matrix orientation is
#'   **column minus row**, not row minus column. The diagonal is `NA`.
#'   \item \code{sd_loa}: numeric \eqn{m \times m} matrix of LoA SDs,
#'   \eqn{\sqrt{\sigma_u^2 + \sigma_e^2}}. This matrix is symmetric.
#'   \item \code{loa_lower}, \code{loa_upper}: numeric \eqn{m \times m} matrices
#'   of LoA endpoints corresponding to `bias`. These satisfy
#'   \eqn{\texttt{loa\_lower}[j,k] = -\texttt{loa\_upper}[k,j]} and
#'   \eqn{\texttt{loa\_upper}[j,k] = -\texttt{loa\_lower}[k,j]}.
#'   \item \code{width}: numeric \eqn{m \times m} matrix of LoA widths,
#'   \code{loa_upper - loa_lower}. This matrix is symmetric.
#'   \item \code{n}: integer \eqn{m \times m} matrix giving the number of complete
#'   subject-time pairs used for each analysed contrast. Pairs with fewer than
#'   two complete matches are left as `NA` in the estimate matrices.
#'   \item \code{mean_ci_low}, \code{mean_ci_high}: numeric \eqn{m \times m}
#'   matrices of Wald confidence interval bounds for `bias`.
#'   \item \code{loa_lower_ci_low}, \code{loa_lower_ci_high}: numeric
#'   \eqn{m \times m} matrices of Wald confidence interval bounds for the lower
#'   LoA.
#'   \item \code{loa_upper_ci_low}, \code{loa_upper_ci_high}: numeric
#'   \eqn{m \times m} matrices of Wald confidence interval bounds for the upper
#'   LoA.
#'   \item \code{slope}: optional numeric \eqn{m \times m} matrix of
#'   proportional-bias slopes on the original paired-mean scale when
#'   `include_slope = TRUE`; otherwise `NULL`. This matrix is antisymmetric in
#'   sign because each fitted contrast is reversed across the transpose.
#'   \item \code{methods}: character vector of method levels defining matrix row
#'   and column order.
#'   \item \code{loa_multiplier}: scalar LoA multiplier actually used.
#'   \item \code{conf_level}: scalar confidence level used for the reported Wald
#'   intervals.
#'   \item \code{use_ar1}: logical, copied from the call.
#'   \item \code{ar1_rho}: scalar equal to the user-supplied common `ar1_rho`
#'   when `use_ar1 = TRUE` and a value was supplied; otherwise `NA`. This field
#'   does \emph{not} store the per-pair estimated AR(1) parameters.
#'   \item \code{residual_model}: character \eqn{m \times m} matrix whose entries
#'   are `"ar1"`, `"iid"`, or `NA`, indicating the final residual structure used
#'   for each pair.
#'   \item \code{sigma2_subject}: numeric \eqn{m \times m} matrix of estimated
#'   subject-level random-intercept variances.
#'   \item \code{sigma2_resid}: numeric \eqn{m \times m} matrix of estimated
#'   residual variances.
#'   \item \code{ar1_rho_pair}: optional numeric \eqn{m \times m} matrix giving
#'   the AR(1) correlation actually used for each pair when the final residual
#'   model is `"ar1"`; otherwise `NA` for that entry. Present only when
#'   `use_ar1 = TRUE`.
#'   \item \code{ar1_estimated}: optional logical \eqn{m \times m} matrix
#'   indicating whether the pair-specific `ar1_rho_pair` was estimated internally
#'   (`TRUE`) or supplied by the user (`FALSE`) for entries whose final residual
#'   model is `"ar1"`; otherwise `NA`. Present only when `use_ar1 = TRUE`.
#'   \item \code{data_long}: canonical long-form data frame with columns
#'   \code{.response}, \code{.subject}, \code{.method}, and \code{.time}
#'   retained for pairwise plotting helpers.
#'   \item \code{mapping}: named list mapping \code{response}, \code{subject},
#'   \code{method}, and \code{time} to the canonical stored columns in
#'   \code{data_long}.
#' }
#'
#' @details
#' For a selected pair of methods \eqn{(a, b)}, the backend first forms complete
#' within-subject pairs at matched \code{subject} and integer-coerced
#' \code{time}. Let
#' \deqn{
#' d_{it} = y_{itb} - y_{ita},
#' \qquad
#' m_{it} = \frac{y_{ita} + y_{itb}}{2},
#' }
#' where \eqn{d_{it}} is the paired difference and \eqn{m_{it}} is the paired
#' mean for subject \eqn{i} at time/replicate \eqn{t}. Only complete
#' subject-time matches contribute to that pairwise fit.
#'
#' If multiple rows are present for the same `subject`-`time`-`method`
#' combination within an analysed pair, the backend keeps the last encountered
#' value for that combination when forming the pair. The function therefore
#' implicitly assumes at most one observation per `subject`-`time`-`method`
#' cell for each analysed contrast.
#'
#' The fitted model for each analysed pair is
#' \deqn{
#' d_{it} = \beta_0 + \beta_1 x_{it} + u_i + \varepsilon_{it},
#' }
#' where \eqn{x_{it} = m_{it}} if `include_slope = TRUE` and the term is omitted
#' otherwise; \eqn{u_i \sim \mathcal{N}(0, \sigma_u^2)} is a subject-specific
#' random intercept; and the within-subject residual vector satisfies
#' \eqn{\mathrm{Cov}(\varepsilon_i) = \sigma_e^2 R_i}.
#'
#' When `use_ar1 = FALSE`, \eqn{R_i = I}. When `use_ar1 = TRUE`, the backend
#' works with the residual \emph{precision} matrix \eqn{C_i = R_i^{-1}} over
#' contiguous time blocks within subject and uses \eqn{\sigma_e^2 C_i^{-1}} as
#' the residual covariance.
#'
#' \subsection{AR(1) residual structure}{
#' Within each subject, paired observations are ordered by integer-coerced
#' `time`. AR(1) correlation is applied only over strictly contiguous runs
#' satisfying \eqn{t_{k+1} = t_k + 1}. Gaps break the run. Negative times, and
#' any isolated positions not belonging to a contiguous run, are treated as
#' independent singletons.
#'
#' For a contiguous run of length \eqn{L} and correlation parameter \eqn{\rho},
#' the block precision matrix is
#' \deqn{
#' C
#' =
#' \frac{1}{1-\rho^2}
#' \begin{bmatrix}
#' 1      & -\rho &        &        &  \\
#' -\rho  & 1+\rho^2 & -\rho &      &  \\
#'         & \ddots & \ddots & \ddots & \\
#'         &        & -\rho & 1+\rho^2 & -\rho \\
#'         &        &       & -\rho & 1
#' \end{bmatrix},
#' }
#' with a very small ridge added to the diagonal for numerical stability.
#'
#' If `use_ar1 = TRUE` and `ar1_rho` is supplied, that value is used after
#' validation and clipping to the admissible numerical range handled by the
#' backend.
#'
#' If `use_ar1 = TRUE` and `ar1_rho = NA`, the backend estimates `rho`
#' separately for each analysed pair by:
#' \enumerate{
#'   \item fitting the corresponding iid model;
#'   \item computing a moments-based lag-1 estimate from detrended residuals
#'   within contiguous blocks, used only as a seed; and
#'   \item refining that seed by a short profile search over `rho` using the
#'   profiled REML log-likelihood.
#' }
#'
#' In the exported `ba_rm()` wrapper, if an AR(1) fit for a given analysed pair
#' fails specifically because the backend EM/GLS routine did not converge to
#' admissible finite variance-component estimates, the wrapper retries that pair
#' with iid residuals. If the iid refit succeeds, the final reported residual
#' model for that pair is `"iid"` and a warning is issued. Other AR(1) failures
#' are not simplified and are propagated as errors.
#' }
#'
#' \subsection{Internal centring and scaling for the proportional-bias slope}{
#' When `include_slope = TRUE`, the paired mean regressor is centred and scaled
#' internally before fitting. Let \eqn{\bar m} be the mean of the observed paired
#' means. The backend chooses a scaling denominator from:
#' \itemize{
#'   \item the sample SD;
#'   \item the IQR-based scale \eqn{\mathrm{IQR}(m)/1.349};
#'   \item the MAD-based scale \eqn{1.4826\,\mathrm{MAD}(m)}.
#' }
#' It uses the first of these that is not judged near-zero relative to the
#' largest finite positive candidate scale, under a threshold proportional to
#' \eqn{\sqrt{\epsilon_{\mathrm{mach}}}}. If all candidate scales are treated as
#' near-zero, the fit stops with an error because the proportional-bias slope is
#' not estimable on the observed paired-mean scale.
#'
#' The returned `beta_slope` is back-transformed to the original paired-mean
#' scale. The returned BA centre is the fitted mean difference at the centred
#' reference paired mean \eqn{\bar m}, not the original-scale intercept
#' coefficient.
#' }
#'
#' \subsection{Estimation}{
#' The backend uses a stabilised EM/GLS scheme.
#'
#' Conditional on current variance components, the fixed effects are updated by
#' GLS using the marginal precision of the paired differences after integrating
#' out the random subject intercept. The resulting fixed-effect covariance used
#' in the confidence-interval calculations is the GLS covariance
#' \deqn{
#' \mathrm{Var}(\hat\beta \mid \hat\theta)
#' =
#' \left( \sum_i X_i^\top V_i^{-1} X_i \right)^{-1}.
#' }
#'
#' Given updated fixed effects, the variance components are refreshed by EM using
#' the conditional moments of the subject random intercept and the residual
#' quadratic forms. Variance updates are ratio-damped and clipped to admissible
#' ranges for numerical stability.
#' }
#'
#' \subsection{Reported BA centre and limits of agreement}{
#' The reported BA centre is always model-based.
#'
#' When `include_slope = FALSE`, it is the fitted intercept of the paired-
#' difference mixed model.
#'
#' When `include_slope = TRUE`, it is the fitted mean difference at the centred
#' reference paired mean used internally by the backend.
#'
#' The reported limits of agreement are
#' \deqn{
#' \mu_0 \pm \texttt{loa\_multiplier}\sqrt{\sigma_u^2 + \sigma_e^2},
#' }
#' where \eqn{\mu_0} is the reported model-based BA centre. These LoA are for a
#' single new paired difference from a random subject under the fitted model.
#'
#' Under the implemented parameterisation, AR(1) correlation affects the
#' off-diagonal within-subject covariance structure and therefore the estimation
#' of the model parameters and their uncertainty, but not the marginal variance
#' of a single paired difference. Consequently `rho` does not appear explicitly
#' in the LoA point-estimate formula.
#' }
#'
#' \subsection{Confidence intervals}{
#' The backend returns Wald confidence intervals for the reported BA centre and
#' for both LoA endpoints.
#'
#' These intervals combine:
#' \itemize{
#'   \item the conditional GLS uncertainty in the fixed effects at the fitted
#'   covariance parameters; and
#'   \item a delta-method propagation of covariance-parameter uncertainty from
#'   the observed information matrix of the profiled REML log-likelihood.
#' }
#'
#' The covariance-parameter vector is profiled on transformed scales:
#' log-variances for \eqn{\sigma_u^2} and \eqn{\sigma_e^2}, and, when `rho` is
#' estimated internally under AR(1), a transformed correlation parameter
#' mapped back by \eqn{\rho = 0.95\tanh(z)}.
#'
#' Numerical central finite differences are used to approximate both the
#' observed Hessian of the profiled REML log-likelihood and the gradients of the
#' reported derived quantities. The resulting variances are combined and the
#' final intervals are formed with the normal quantile corresponding to
#' `conf_level`.
#' }
#'
#' \subsection{Exactly two methods versus \eqn{\ge 3} methods}{
#' With exactly two methods, at least two complete subject-time pairs are
#' required; otherwise the function errors.
#'
#' With \eqn{\ge 3} methods, the function analyses every unordered pair of method
#' levels. For a given pair with fewer than two complete subject-time matches,
#' that contrast is skipped and the corresponding matrix entries remain `NA`.
#'
#' For a fitted contrast between methods in matrix positions \eqn{(j, k)} with
#' \eqn{j < k}, the stored orientation is:
#' \deqn{
#' \texttt{bias}[j,k] \approx \mathrm{method}_k - \mathrm{method}_j.
#' }
#' Hence the transposed entry changes sign, while `sd_loa` and `width` are
#' symmetric.
#' }
#'
#' \subsection{Identifiability and safeguards}{
#' Separate estimation of the residual and subject-level variance components
#' requires sufficient complete within-subject replication after pairing. If the
#' paired data are not adequate to separate these components, the fit stops with
#' an identifiability error.
#'
#' If the model is conceptually estimable but no finite positive pooled
#' within-subject variance can be formed during initialisation, the backend uses
#' \eqn{0.5 \times v_{\mathrm{ref}}} only as a temporary positive starting value
#' for the EM routine and records a warning string in the backend output. The
#' exported wrapper does not otherwise modify the final estimates.
#'
#' If the EM/GLS routine fails to reach admissible finite variance-component
#' estimates, the backend throws an explicit convergence error rather than
#' returning fallback estimates.
#' }
#'
#' @examples
#' # -------- Simulate repeated-measures data --------
#' set.seed(1)
#'
#' # design (no AR)
#' # subjects
#' S   <- 30L
#' # replicates per subject
#' Tm  <- 15L
#' subj <- rep(seq_len(S), each = Tm)
#' time <- rep(seq_len(Tm), times = S)
#'
#' # subject signal centered at 0 so BA "bias" won't be driven by the mean level
#' mu_s  <- rnorm(S, mean = 0, sd = 8)
#' # constant within subject across replicates
#' true  <- mu_s[subj]
#'
#' # common noise (no AR, i.i.d.)
#' sd_e <- 2
#' e0   <- rnorm(length(true), 0, sd_e)
#'
#' # --- Methods ---
#' # M1: signal + noise
#' y1 <- true + e0
#'
#' # M2: same precision as M1; here identical so M3 can be
#' #     almost perfectly the inverse of both M1 and M2
#' y2 <- y1 + rnorm(length(true), 0, 0.01)
#'
#' # M3: perfect inverse of M1 and M2
#' y3 <- -y1   # = -(true + e0)
#'
#' # M4: unrelated to all others (pure noise, different scale)
#' y4 <- rnorm(length(true), 3, 6)
#'
#' data <- rbind(
#'   data.frame(y = y1, subject = subj, method = "M1", time = time),
#'   data.frame(y = y2, subject = subj, method = "M2", time = time),
#'   data.frame(y = y3, subject = subj, method = "M3", time = time),
#'   data.frame(y = y4, subject = subj, method = "M4", time = time)
#' )
#' data$method <- factor(data$method, levels = c("M1","M2","M3","M4"))
#'
#' # quick sanity checks
#' with(data, {
#'   Y <- split(y, method)
#'   round(cor(cbind(M1 = Y$M1, M2 = Y$M2, M3 = Y$M3, M4 = Y$M4)), 3)
#' })
#'
#' # Run BA (no AR)
#' ba4 <- ba_rm(
#'   data = data,
#'   response = "y", subject = "subject", method = "method", time = "time",
#'   loa_multiplier = 1.96, conf_level = 0.95,
#'   include_slope = FALSE, use_ar1 = FALSE
#' )
#' summary(ba4)
#' plot(ba4)
#'
#' # -------- Simulate repeated-measures with AR(1) data --------
#' set.seed(123)
#' S <- 40L                      # subjects
#' Tm <- 50L                     # replicates per subject
#' methods <- c("A","B","C")     # N = 3 methods
#' rho <- 0.4                    # AR(1) within-subject across time
#'
#' ar1_sim <- function(n, rho, sd = 1) {
#'   z <- rnorm(n)
#'   e <- numeric(n)
#'   e[1] <- z[1] * sd
#'   if (n > 1) for (t in 2:n) e[t] <- rho * e[t-1] + sqrt(1 - rho^2) * z[t] * sd
#'   e
#' }
#'
#' # Subject baseline + time trend (latent "true" signal)
#' subj <- rep(seq_len(S), each = Tm)
#' time <- rep(seq_len(Tm), times = S)
#' # subject effects
#' mu_s  <- rnorm(S, 50, 7)
#' trend <- rep(seq_len(Tm) - mean(seq_len(Tm)), times = S) * 0.8
#' true  <- mu_s[subj] + trend
#'
#' # Method-specific biases (B has +1.5 constant; C has slight proportional bias)
#' bias  <- c(A = 0, B = 1.5, C = -0.5)
#' # proportional component on "true"
#' prop  <- c(A = 0.00, B = 0.00, C = 0.10)
#'
#' # Build long data: for each method, add AR(1) noise within subject over time
#' make_method <- function(meth, sd = 3) {
#'   e <- unlist(lapply(split(seq_along(time), subj),
#'                      function(ix) ar1_sim(length(ix), rho, sd)))
#'   y <- true * (1 + prop[meth]) + bias[meth] + e
#'   data.frame(y = y, subject = subj, method = meth, time = time,
#'              check.names = FALSE)
#' }
#'
#' data <- do.call(rbind, lapply(methods, make_method))
#' data$method <- factor(data$method, levels = methods)
#'
#' # -------- Repeated BA (pairwise matrix) ---------------------
#' baN <- ba_rm(
#'   response = data$y, subject = data$subject, method = data$method, time = data$time,
#'   loa_multiplier = 1.96, conf_level = 0.95,
#'   include_slope = FALSE,         # estimate proportional bias per pair
#'   use_ar1 = TRUE, ar1_rho = rho
#' )
#'
#' # Matrices (row - column orientation)
#' print(baN)
#' summary(baN)
#'
#' # Faceted BA scatter by pair
#' plot(baN, smoother = "lm", facet_scales = "free_y")
#'
#' # -------- Two-method AR(1) path (A vs B only) ------------------------------
#' data_AB <- subset(data, method %in% c("A","B"))
#' baAB <- ba_rm(
#'   response = data_AB$y, subject = data_AB$subject,
#'   method = droplevels(data_AB$method), time = data_AB$time,
#'   include_slope = FALSE, use_ar1 = TRUE, ar1_rho = 0.4
#' )
#' print(baAB)
#' plot(baAB)
#'
#' @author Thiago de Paula Oliveira
#' @export
ba_rm <- function(data = NULL, response, subject, method, time,
                                  conf_level = 0.95,
                                  n_threads = getOption("matrixCorr.threads", 1L),
                                  verbose = FALSE,
                                  loa_multiplier = 1.96,
                                  include_slope = FALSE,
                                  use_ar1 = FALSE, ar1_rho = NA_real_,
                                  max_iter = 200L, tol = 1e-6) {

  # legacy positional signature support: ba_rm(response, subject, method, time, ...)
  if (!missing(data) && !inherits(data, "data.frame") && missing(time)) {
    time <- method
    method <- subject
    subject <- response
    response <- data
    data <- NULL
  }

  if (!is.null(data) && !inherits(data, "data.frame")) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  # --- resolve columns if 'data' provided and names given ---
  if (!is.null(data)) {
    if (!inherits(data, "data.frame")) {
      abort_bad_arg("data",
        message = "must be a data.frame or data.table."
      )
    }
    pull <- function(x) {
      if (is.character(x) && length(x) == 1L) {
        if (!x %in% names(data)) {
          abort_bad_arg("data",
            message = "Column '{x}' not found in `data`.",
            x = x
          )
        }
        data[[x]]
      } else x
    }
    response <- pull(response)
    subject  <- pull(subject)
    method   <- pull(method)
    time     <- pull(time)
    mapping <- list(
      response =
        if (is.character(substitute(response)) && length(substitute(response)) == 1L &&
                     is.character(match.call()$response)) as.character(match.call()$response) else ".response",
      subject  = if (is.character(match.call()$subject)) as.character(match.call()$subject) else ".subject",
      method   = if (is.character(match.call()$method))  as.character(match.call()$method)  else ".method",
      time     = if (is.character(match.call()$time))    as.character(match.call()$time)    else ".time"
    )

    # always use canonical internal names in the stored table
    data_long <- setNames(
      data.frame(response, subject, method, time, check.names = FALSE),
      c(".response", ".subject", ".method", ".time")
    )
    mapping <- list(
      response = ".response",
      subject  = ".subject",
      method   = ".method",
      time     = ".time"
    )
  } else {
    # vectors supplied; build canonical internal columns
    mapping <- list(
      response = ".response",
      subject  = ".subject",
      method   = ".method",
      time     = ".time"
    )
    data_long <- data.frame(
      ".response" = response,
      ".subject"  = subject,
      ".method"   = method,
      ".time"     = time,
      check.names = FALSE
    )
  }

  # ---- validate / normalise types (work with local copies) ----
  y <- data_long[[mapping$response]]
  s <- data_long[[mapping$subject ]]
  m <- data_long[[mapping$method  ]]
  t <- data_long[[mapping$time    ]]

  if (!is.numeric(y) || length(y) < 2L) {
    abort_bad_arg("response",
      message = "must be numeric with length > 1."
    )
  }
  if (!(is.integer(s) || is.factor(s) || is.numeric(s))) {
    abort_bad_arg("subject",
      message = "must be integer, factor, or numeric."
    )
  }
  s <- as.integer(as.factor(s))
  if (!is.factor(m)) m <- as.factor(m)
  m <- droplevels(m)
  mlev <- levels(m)
  if (length(mlev) < 2L) {
    abort_bad_arg("method",
      message = "Need at least 2 distinct methods in `method`."
    )
  }
  if (!(is.integer(t) || is.numeric(t))) {
    abort_bad_arg("time",
      message = "must be integer or numeric."
    )
  }
  t <- as.integer(t)

  if (!is.numeric(loa_multiplier) || length(loa_multiplier) != 1L || !is.finite(loa_multiplier) || loa_multiplier <= 0) {
    abort_bad_arg("loa_multiplier",
      message = "`loa_multiplier` must be a positive scalar."
    )
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    abort_bad_arg("conf_level",
      message = "`conf_level` must be in (0,1)."
    )
  }
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")
  check_bool(include_slope, arg = "include_slope")
  check_bool(use_ar1, arg = "use_ar1")
  check_bool(verbose, arg = "verbose")
  prev_threads <- get_omp_threads()
  on.exit(set_omp_threads(as.integer(prev_threads)), add = TRUE)
  if (isTRUE(use_ar1)) {
    if (!is.na(ar1_rho)) {
      ar1_rho <- check_ar1_rho(ar1_rho, arg = "ar1_rho", bound = 0.999)
    }
  } else {
    ar1_rho <- NA_real_
  }

  if (isTRUE(verbose)) {
    cat("Using", n_threads, "OpenMP threads\n")
  }

  # ---- two-method path -------------------------------------------------------
  if (length(mlev) == 2L) {
    idx <- m %in% mlev & is.finite(y) & !is.na(s) & !is.na(t)
    res1 <- .ba_rep_two_methods(
      response = y[idx],
      subject  = s[idx],
      method12 = as.integer(m[idx] == mlev[2L]) + 1L,
      time     = t[idx],
      loa_multiplier = loa_multiplier, conf_level = conf_level, n_threads = n_threads,
      include_slope = include_slope,
      use_ar1 = use_ar1, ar1_rho = ar1_rho, max_iter = max_iter, tol = tol
    )
    attr(res1, "conf.level") <- conf_level
    return(res1)
  }

  # ---- N-method pairwise path ------------------------------------------------
  if (isTRUE(verbose)) {
    cat("Repeated BA (pairwise):", length(mlev), "methods ->",
        choose(length(mlev), 2L), "pairs\n")
  }

  methods <- mlev
  mm <- length(methods)

  bias         <- .make_named_matrix_ba(methods)
  sd_loa       <- .make_named_matrix_ba(methods)
  loa_lower    <- .make_named_matrix_ba(methods)
  loa_upper    <- .make_named_matrix_ba(methods)
  width        <- .make_named_matrix_ba(methods)
  n_mat        <- .make_named_matrix_ba(methods, NA_integer_, "integer")
  mean_ci_low  <- .make_named_matrix_ba(methods)
  mean_ci_high <- .make_named_matrix_ba(methods)
  lo_ci_low    <- .make_named_matrix_ba(methods)
  lo_ci_high   <- .make_named_matrix_ba(methods)
  hi_ci_low    <- .make_named_matrix_ba(methods)
  hi_ci_high   <- .make_named_matrix_ba(methods)
  slope_mat    <- if (isTRUE(include_slope)) .make_named_matrix_ba(methods) else NULL
  vc_subject   <- .make_named_matrix_ba(methods)
  vc_resid     <- .make_named_matrix_ba(methods)
  ar1_rho_mat   <- if (isTRUE(use_ar1)) .make_named_matrix_ba(methods) else NULL
  residual_model <- matrix(NA_character_, mm, mm, dimnames = list(methods, methods))
  ar1_estimated <- if (isTRUE(use_ar1)) {
    z <- matrix(NA, mm, mm, dimnames = list(methods, methods))
    storage.mode(z) <- "logical"
    z
  } else NULL
  ar1_simplified_pairs <- character()

  .recompose_pair <- function(fit) {
    list(
      md    = as.numeric(fit$bias_mu0),
      sd    = as.numeric(fit$sd_loa),
      lo    = as.numeric(fit$loa_lower),
      hi    = as.numeric(fit$loa_upper),
      bias_l= as.numeric(fit$bias_lwr),
      bias_u= as.numeric(fit$bias_upr),
      lo_l  = as.numeric(fit$loa_lower_lwr),
      lo_u  = as.numeric(fit$loa_lower_upr),
      hi_l  = as.numeric(fit$loa_upper_lwr),
      hi_u  = as.numeric(fit$loa_upper_upr)
    )
  }

  for (j in 1:(mm - 1L)) for (k in (j + 1L):mm) {
    lev_j <- methods[j]; lev_k <- methods[k]
    sel <- m %in% c(lev_j, lev_k) & is.finite(y) & !is.na(s) & !is.na(t)
    if (!any(sel)) next
    m12 <- ifelse(m[sel] == lev_j, 1L, 2L)
    n_pair <- .count_ba_rm_complete_pairs(y[sel], s[sel], m12, t[sel])
    n_mat[j,k] <- n_mat[k,j] <- n_pair

    if (n_pair < 2L) {
      bias[j,k]      <- bias[k,j]      <- NA_real_
      sd_loa[j,k]    <- sd_loa[k,j]    <- NA_real_
      loa_lower[j,k] <- loa_lower[k,j] <- NA_real_
      loa_upper[j,k] <- loa_upper[k,j] <- NA_real_
      width[j,k]     <- width[k,j]     <- NA_real_

      mean_ci_low[j,k]  <- mean_ci_low[k,j]   <- NA_real_
      mean_ci_high[j,k] <- mean_ci_high[k,j]  <- NA_real_
      lo_ci_low[j,k]    <- lo_ci_low[k,j]     <- NA_real_
      lo_ci_high[j,k]   <- lo_ci_high[k,j]    <- NA_real_
      hi_ci_low[j,k]    <- hi_ci_low[k,j]     <- NA_real_
      hi_ci_high[j,k]   <- hi_ci_high[k,j]    <- NA_real_

      vc_subject[j,k] <- vc_subject[k,j] <- NA_real_
      vc_resid[j,k]   <- vc_resid[k,j]   <- NA_real_
      residual_model[j,k] <- residual_model[k,j] <- NA_character_

      if (!is.null(slope_mat)) {
        slope_mat[j,k] <- slope_mat[k,j] <- NA_real_
      }
      if (isTRUE(use_ar1)) {
        ar1_rho_mat[j,k] <- ar1_rho_mat[k,j] <- NA_real_
        ar1_estimated[j,k] <- ar1_estimated[k,j] <- NA
      }

      inform_if_verbose(
        "Skipping Bland-Altman pair {.val {lev_j}} vs {.val {lev_k}}: need at least 2 complete subject-time pairs (found {n_pair}).",
        .verbose = verbose
      )
      next
    }

    pair_label <- paste(lev_j, lev_k, sep = "-")
    fit_info <- .fit_ba_rm_pair_model(
      response = y[sel], subject = s[sel], method12 = m12, time = t[sel],
      n_threads = n_threads,
      include_slope = include_slope, use_ar1 = use_ar1, ar1_rho = ar1_rho,
      max_iter = max_iter, tol = tol, conf_level = conf_level,
      loa_multiplier = loa_multiplier, pair_label = pair_label
    )
    fit <- fit_info$fit
    comp <- .recompose_pair(fit)

    bias[j,k]      <- comp$md;  bias[k,j]      <- -comp$md
    sd_loa[j,k]    <- comp$sd;  sd_loa[k,j]    <-  comp$sd
    loa_lower[j,k] <- comp$lo;  loa_lower[k,j] <- -comp$hi
    loa_upper[j,k] <- comp$hi;  loa_upper[k,j] <- -comp$lo
    width[j,k]     <- comp$hi - comp$lo; width[k,j] <- width[j,k]

    mean_ci_low[j,k]  <- comp$bias_l; mean_ci_low[k,j]  <- -comp$bias_u
    mean_ci_high[j,k] <- comp$bias_u; mean_ci_high[k,j] <- -comp$bias_l
    lo_ci_low[j,k]    <- comp$lo_l;   lo_ci_low[k,j]    <- -comp$hi_u
    lo_ci_high[j,k]   <- comp$lo_u;   lo_ci_high[k,j]   <- -comp$hi_l
    hi_ci_low[j,k]    <- comp$hi_l;   hi_ci_low[k,j]    <- -comp$lo_u
    hi_ci_high[j,k]   <- comp$hi_u;   hi_ci_high[k,j]   <- -comp$lo_l

    vc_subject[j,k] <- vc_subject[k,j] <- as.numeric(fit$sigma2_subject)
    vc_resid[j,k]   <- vc_resid[k,j]   <- as.numeric(fit$sigma2_resid)
    residual_model[j,k] <- residual_model[k,j] <- fit_info$residual_model

    if (isTRUE(use_ar1)) {
      ar1_rho_used <- if (identical(fit_info$residual_model, "ar1")) as.numeric(fit$ar1_rho) else NA_real_
      ar1_rho_mat[j,k] <- ar1_rho_mat[k,j] <- ar1_rho_used
      ar1_estimated[j,k] <- ar1_estimated[k,j] <- if (identical(fit_info$residual_model, "ar1")) isTRUE(fit$ar1_estimated) else NA
      if (isTRUE(fit_info$ar1_simplified)) {
        ar1_simplified_pairs <- c(ar1_simplified_pairs, pair_label)
      }
    }

    if (!is.null(slope_mat)) {
      slope <- as.numeric(fit$beta_slope)
      slope_mat[j,k] <-  slope
      slope_mat[k,j] <- -slope
    }

    if (isTRUE(verbose)) {
      cat(sprintf(" pair %s - %s: n=%d, bias=%.4f, sd=%.4f\n",
                  lev_j, lev_k, as.integer(fit$n_pairs), comp$md, comp$sd))
    }
  }

  ba_repeated <- list(
    bias = bias,
    sd_loa = sd_loa,
    loa_lower = loa_lower,
    loa_upper = loa_upper,
    width = width,
    n = n_mat,
    mean_ci_low = mean_ci_low,
    mean_ci_high = mean_ci_high,
    loa_lower_ci_low = lo_ci_low,
    loa_lower_ci_high = lo_ci_high,
    loa_upper_ci_low = hi_ci_low,
    loa_upper_ci_high = hi_ci_high,
    slope = slope_mat,
    methods = methods,
    loa_multiplier = loa_multiplier,
    conf_level = conf_level,
    use_ar1 = use_ar1,
    ar1_rho = if (use_ar1) ar1_rho else NA_real_,
    residual_model = residual_model,
    sigma2_subject = vc_subject,
    sigma2_resid   = vc_resid,
    ar1_rho_pair   = if (use_ar1) ar1_rho_mat else NULL,
    ar1_estimated  = if (use_ar1) ar1_estimated else NULL,
    data_long = data_long,
    mapping = mapping
  )
  ba_repeated <- structure(ba_repeated, class = c("ba_repeated_matrix","list"))
  attr(ba_repeated, "conf.level") <- conf_level
  if (length(ar1_simplified_pairs)) {
    .warn_ba_rm_ar1_fallback(ar1_simplified_pairs)
  }
  ba_repeated
}




#' @keywords internal
.count_ba_rm_complete_pairs <- function(response, subject, method12, time) {
  ba_rm_complete_pairs_cpp(response, subject, method12, time)
}


.ba_rm_is_ar1_convergence_failure <- function(cnd) {
  grepl(
    "failed to converge to admissible finite variance-component estimates",
    conditionMessage(cnd),
    fixed = TRUE
  )
}

.warn_ba_rm_ar1_fallback <- function(pair_labels = NULL) {
  where <- if (length(pair_labels)) {
    sprintf(" for pair(s): %s", paste(unique(pair_labels), collapse = ", "))
  } else {
    " for this fit"
  }

  warning(
    sprintf(
      "Requested AR(1) residual structure could not be fit%s; using iid residuals instead.",
      where
    ),
    call. = FALSE
  )
}


.fit_ba_rm_pair_model <- function(response, subject, method12, time,
                                  n_threads,
                                  include_slope, use_ar1, ar1_rho,
                                  max_iter, tol, conf_level, loa_multiplier,
                                  pair_label = NULL) {
  run_backend <- function(use_ar1_fit, ar1_rho_fit = NA_real_) {
    bland_altman_repeated_em_ext_cpp(
      y = response, subject = subject, method = method12, time = time,
      include_slope = include_slope,
      use_ar1 = use_ar1_fit, ar1_rho = ar1_rho_fit,
      max_iter = max_iter, tol = tol, conf_level = conf_level,
      loa_multiplier_arg = loa_multiplier,
      n_threads = n_threads
    )
  }

  if (!isTRUE(use_ar1)) {
    return(list(
      fit = run_backend(FALSE),
      residual_model = "iid",
      ar1_simplified = FALSE,
      ar1_simplification_reason = NA_character_
    ))
  }

  fit_ar1 <- tryCatch(
    run_backend(TRUE, ar1_rho),
    error = identity
  )
  if (!inherits(fit_ar1, "error")) {
    return(list(
      fit = fit_ar1,
      residual_model = "ar1",
      ar1_simplified = FALSE,
      ar1_simplification_reason = NA_character_
    ))
  }

  if (!.ba_rm_is_ar1_convergence_failure(fit_ar1)) {
    stop(fit_ar1)
  }

  fit_iid <- tryCatch(
    run_backend(FALSE),
    error = identity
  )
  if (inherits(fit_iid, "error")) {
    stop(fit_iid)
  }

  list(
    fit = fit_iid,
    residual_model = "iid",
    ar1_simplified = TRUE,
    ar1_simplification_reason = "ar1_fit_failed_then_iid_succeeded",
    pair_label = pair_label
  )
}


#' two-method helper
#' @keywords internal
.ba_rep_two_methods <- function(response, subject, method12, time,
                                loa_multiplier, conf_level, n_threads, include_slope,
                                use_ar1, ar1_rho, max_iter, tol) {
  n_pairs <- .count_ba_rm_complete_pairs(response, subject, method12, time)
  if (n_pairs < 2L) {
    abort_bad_arg("response",
      message = "must provide at least two subject-time matched pairs after removing missing observations; only {n_pairs} available.",
      n_pairs = n_pairs
    )
  }

  fit_info <- .fit_ba_rm_pair_model(
    response = response, subject = subject, method12 = method12, time = time,
    n_threads = n_threads,
    include_slope = include_slope, use_ar1 = use_ar1, ar1_rho = ar1_rho,
    max_iter = max_iter, tol = tol, conf_level = conf_level,
    loa_multiplier = loa_multiplier
  )
  fit <- fit_info$fit
  n_pairs <- as.integer(fit$n_pairs)

  md  <- as.numeric(fit$bias_mu0)
  sdL <- as.numeric(fit$sd_loa)

  loa_lower <- as.numeric(fit$loa_lower)
  loa_upper <- as.numeric(fit$loa_upper)

  CI.lines <- c(
    "mean.diff.ci.lower"   = as.numeric(fit$bias_lwr),
    "mean.diff.ci.upper"   = as.numeric(fit$bias_upr),
    "lower.limit.ci.lower" = as.numeric(fit$loa_lower_lwr),
    "lower.limit.ci.upper" = as.numeric(fit$loa_lower_upr),
    "upper.limit.ci.lower" = as.numeric(fit$loa_upper_lwr),
    "upper.limit.ci.upper" = as.numeric(fit$loa_upper_upr)
  )

  means <- as.numeric(fit$pairs_mean)
  diffs <- as.numeric(fit$pairs_diff)

  ba_repeated <- list(
    means         = means,
    diffs         = diffs,
    n_obs         = n_pairs,
    based.on      = n_pairs,
    lower.limit   = loa_lower,
    mean.diffs    = md,
    upper.limit   = loa_upper,
    lines         = c(lower = loa_lower, mean = md, upper = loa_upper),
    CI.lines      = CI.lines,
    loa_multiplier = loa_multiplier,
    critical.diff = loa_multiplier * sdL,
    include_slope = include_slope,
    beta_slope    = if (include_slope) as.numeric(fit$beta_slope) else NA_real_,
    sigma2_subject= as.numeric(fit$sigma2_subject),
    sigma2_resid  = as.numeric(fit$sigma2_resid),
    use_ar1       = use_ar1,
    residual_model = fit_info$residual_model,
    ar1_rho       = if (identical(fit_info$residual_model, "ar1")) as.numeric(fit$ar1_rho) else NA_real_,
    ar1_estimated = if (identical(fit_info$residual_model, "ar1")) isTRUE(fit$ar1_estimated) else NA
  )
  ba_repeated <- structure(ba_repeated, class = c("ba_repeated","list"))
  attr(ba_repeated, "conf.level") <- conf_level
  if (isTRUE(fit_info$ar1_simplified)) {
    .warn_ba_rm_ar1_fallback()
  }
  ba_repeated
}




# ---------- printers ----------
#' @rdname ba_rm
#' @method print ba_repeated
#' @param x A \code{"ba_repeated"} object.
#' @param digits Number of digits for estimates (default 3).
#' @param ci_digits Number of digits for CI bounds (default 3).
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to `getOption("width")`.
#' @param show_ci One of `"yes"` or `"no"`.
#' @param ... Unused.
#' @export
print.ba_repeated <- function(x, digits = 3, ci_digits = 3,
                              n = NULL, topn = NULL,
                              max_vars = NULL, width = NULL,
                              show_ci = NULL, ...) {
  check_inherits(x, "ba_repeated")
  show_ci <- .mc_resolve_show_ci(show_ci, context = "print")
  display_width <- width
  n   <- as.integer(x$based.on)
  loa_multiplier <- as.numeric(x$loa_multiplier)
  cl  <- suppressWarnings(as.numeric(attr(x, "conf.level")))
  if (!is.finite(cl)) cl <- NA_real_

  md    <- as.numeric(x$mean.diffs)
  loaL  <- as.numeric(x$lower.limit)
  loaU  <- as.numeric(x$upper.limit)
  sd_d  <- as.numeric(x$critical.diff) / loa_multiplier
  loa_width <- loaU - loaL

  cil <- function(nm) as.numeric(x$CI.lines[[nm]])
  bias_l <- cil("mean.diff.ci.lower"); bias_u <- cil("mean.diff.ci.upper")
  lo_l   <- cil("lower.limit.ci.lower"); lo_u <- cil("lower.limit.ci.upper")
  hi_l   <- cil("upper.limit.ci.lower"); hi_u <- cil("upper.limit.ci.upper")

  .mc_print_named_digest(
    c(
      pairs = n,
      loa_rule = sprintf("mean +/- %.3g * SD", loa_multiplier),
      if (identical(show_ci, "yes") && is.finite(cl)) c(ci = sprintf("%g%%", 100 * cl)),
      sd_single_pair = formatC(sd_d, format = "f", digits = digits),
      width = formatC(loa_width, format = "f", digits = digits)
    ),
    header = "Repeated-measures Bland-Altman preview:"
  )
  cat("\n")

  df <- data.frame(
    quantity = c("Mean difference", "Lower LoA", "Upper LoA"),
    estimate = c(md, loaL, loaU),
    lwr      = c(bias_l, lo_l, hi_l),
    upr      = c(bias_u, lo_u, hi_u),
    check.names = FALSE
  )
  df$estimate <- formatC(df$estimate, format = "f", digits = digits)
  if (identical(show_ci, "yes")) {
    df$lwr <- formatC(df$lwr, format = "f", digits = ci_digits)
    df$upr <- formatC(df$upr, format = "f", digits = ci_digits)
  } else {
    df$lwr <- NULL
    df$upr <- NULL
  }
  .mc_print_preview_table(
    df,
    n = .mc_coalesce(n, .mc_display_option("print_max_rows", 20L)),
    topn = .mc_coalesce(topn, .mc_display_option("print_topn", 5L)),
    max_vars = .mc_coalesce(max_vars, .mc_display_option("print_max_vars", NULL)),
    width = .mc_coalesce(display_width, getOption("width", 80L)),
    context = "print",
    full_hint = TRUE,
    summary_hint = TRUE,
    ...
  )

  if (isTRUE(x$include_slope) && is.finite(x$beta_slope)) {
    cat(sprintf("\nProportional bias slope (vs pair mean): %s\n",
                formatC(x$beta_slope, format = "f", digits = digits)))
  } else {
    cat("\n")
  }
  invisible(x)
}

#' @rdname ba_rm
#' @method print ba_repeated_matrix
#' @param x A \code{"ba_repeated_matrix"} object.
#' @param digits Number of digits for estimates (default 3).
#' @param ci_digits Number of digits for CI bounds (default 3).
#' @param style Show as pairs or matrix format?
#' @param ... Unused.
#' @export
print.ba_repeated_matrix <- function(x,
                                     digits = 3,
                                     ci_digits = 3,
                                     n = NULL,
                                     topn = NULL,
                                     max_vars = NULL,
                                     width = NULL,
                                     show_ci = NULL,
                                     style = c("pairs","matrices"),
                                     ...) {
  check_inherits(x, "ba_repeated_matrix")
  style <- match.arg(style)
  show_ci <- .mc_resolve_show_ci(show_ci, context = "print")
  cl <- suppressWarnings(as.numeric(attr(x, "conf.level")))
  has_ci <- all(c("mean_ci_low","mean_ci_high",
                  "loa_lower_ci_low","loa_lower_ci_high",
                  "loa_upper_ci_low","loa_upper_ci_high") %in% names(x))

  if (style == "matrices") {
    .mc_print_named_digest(
      c(
        methods = nrow(x$bias),
        if (identical(show_ci, "yes") && is.finite(cl)) c(ci = sprintf("%g%%", 100 * cl)),
        components = "bias, sd_loa, loa_low, loa_up, width, n"
      ),
      header = "Repeated-measures Bland-Altman matrix preview:"
    )
    cat("\nPrimary component: bias\n\n")
    .mc_print_preview_matrix(
      x$bias,
      digits = digits,
      n = .mc_coalesce(n, .mc_display_option("print_max_rows", 20L)),
      topn = .mc_coalesce(topn, .mc_display_option("print_topn", 5L)),
      max_vars = .mc_coalesce(max_vars, .mc_display_option("print_max_vars", NULL)),
      width = .mc_coalesce(width, getOption("width", 80L)),
      context = "print",
      full_hint = TRUE,
      summary_hint = TRUE,
      ...
    )
    return(invisible(x))
  }

  sm <- summary(x, digits = digits, ci_digits = ci_digits)
  cols <- c("method1", "method2", "bias", "sd_loa", "loa_low", "loa_up", "width", "n")
  if ("slope" %in% names(sm)) cols <- c(cols, "slope")
  if (identical(show_ci, "yes") && has_ci) {
    cols <- c(cols, "bias_lwr", "bias_upr", "lo_lwr", "lo_upr", "up_lwr", "up_upr")
  }
  df <- as.data.frame(sm)[, cols[cols %in% names(sm)], drop = FALSE]

  header <- .mc_header_with_ci("Bland-Altman (row \u2212 column)", cl, if (has_ci) show_ci else "no")
  .mc_print_summary_table(
    df,
    header = header,
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

#' @method summary ba_repeated
#' @export
summary.ba_repeated <- function(object,
                                digits = 3,
                                ci_digits = 3,
                                ...) {
  check_inherits(object, "ba_repeated")
  cl <- suppressWarnings(as.numeric(attr(object, "conf.level")))
  n  <- as.integer(object$based.on)

  ba_repeated <- data.frame(
    method1   = "method 1",
    method2   = "method 2",
    bias      = round(num_or_na_ba(object$mean.diffs), digits),
    sd_loa    = round(num_or_na_ba(object$critical.diff) / num_or_na_ba(object$loa_multiplier), digits),
    loa_low   = round(num_or_na_ba(object$lower.limit), digits),
    loa_up    = round(num_or_na_ba(object$upper.limit), digits),
    width     = round(num_or_na_ba(object$upper.limit - object$lower.limit), digits),
    n_obs     = n,
    stringsAsFactors = FALSE, check.names = FALSE
  )

  if (isTRUE(object$include_slope) && is.finite(num_or_na_ba(object$beta_slope))) {
    ba_repeated$slope <- round(num_or_na_ba(object$beta_slope), digits)
  }

  cil <- function(nm) num_or_na_ba(object$CI.lines[[nm]])
  if (!any(is.na(c(cil("mean.diff.ci.lower"), cil("mean.diff.ci.upper"),
                   cil("lower.limit.ci.lower"), cil("lower.limit.ci.upper"),
                   cil("upper.limit.ci.lower"), cil("upper.limit.ci.upper"))))) {
    ba_repeated$bias_lwr <- round(cil("mean.diff.ci.lower"), ci_digits)
    ba_repeated$bias_upr <- round(cil("mean.diff.ci.upper"), ci_digits)
    ba_repeated$lo_lwr   <- round(cil("lower.limit.ci.lower"), ci_digits)
    ba_repeated$lo_upr   <- round(cil("lower.limit.ci.upper"), ci_digits)
    ba_repeated$up_lwr   <- round(cil("upper.limit.ci.lower"), ci_digits)
    ba_repeated$up_upr   <- round(cil("upper.limit.ci.upper"), ci_digits)
  }

  ba_repeated$sigma2_subject <- round(num_or_na_ba(object$sigma2_subject), digits)
  ba_repeated$sigma2_resid   <- round(num_or_na_ba(object$sigma2_resid),   digits)
  ba_repeated$use_ar1        <- isTRUE(object$use_ar1)
  ba_repeated$residual_model <- if (!is.null(object$residual_model)) object$residual_model else if (isTRUE(object$use_ar1)) "ar1" else "iid"
  if (identical(ba_repeated$residual_model, "ar1")) {
    ba_repeated$ar1_rho       <- round(num_or_na_ba(object$ar1_rho), digits)
    ba_repeated$ar1_estimated <- object$ar1_estimated
  }

  ba_repeated <- .mc_finalize_summary_df(
    ba_repeated,
    class_name = "summary.ba_repeated",
    repeated = TRUE
  )
  attr(ba_repeated, "conf.level") <- cl
  ba_repeated
}

#' @method summary ba_repeated_matrix
#' @export
summary.ba_repeated_matrix <- function(object,
                                       digits = 3,
                                       ci_digits = 3,
                                       ...) {
  check_inherits(object, "ba_repeated_matrix")
  cl <- suppressWarnings(as.numeric(attr(object, "conf.level")))
  methods <- object$methods
  m <- length(methods)

  has_ci <- all(c("mean_ci_low","mean_ci_high",
                  "loa_lower_ci_low","loa_lower_ci_high",
                  "loa_upper_ci_low","loa_upper_ci_high") %in% names(object))

  rows <- vector("list", m * (m - 1L) / 2L)
  k <- 0L
  for (i in 1:(m - 1L)) for (j in (i + 1L):m) {
    k <- k + 1L
    row <- list(
      method1 = methods[i],
      method2 = methods[j],
      bias    = round(object$bias[i, j], digits),
      sd_loa  = round(object$sd_loa[i, j], digits),
      loa_low = round(object$loa_lower[i, j], digits),
      loa_up  = round(object$loa_upper[i, j], digits),
      width   = round(object$width[i, j], digits),
      n_obs   = suppressWarnings(as.integer(object$n[i, j])),
      sigma2_subject = round(object$sigma2_subject[i, j], digits),
      sigma2_resid   = round(object$sigma2_resid[i, j],   digits)
    )

    if (!is.null(object$slope)) {
      row$slope <- round(object$slope[i, j], digits)
    }

    if (has_ci && is.finite(cl)) {
      row$bias_lwr <- round(object$mean_ci_low[i, j],       ci_digits)
      row$bias_upr <- round(object$mean_ci_high[i, j],      ci_digits)
      row$lo_lwr   <- round(object$loa_lower_ci_low[i, j],  ci_digits)
      row$lo_upr   <- round(object$loa_lower_ci_high[i, j], ci_digits)
      row$up_lwr   <- round(object$loa_upper_ci_low[i, j],  ci_digits)
      row$up_upr   <- round(object$loa_upper_ci_high[i, j], ci_digits)
    }
    if (!is.null(object$residual_model)) {
      row$residual_model <- object$residual_model[i, j]
    }
    if (isTRUE(object$use_ar1)) {
      rho_ij <- if (!is.null(object$ar1_rho_pair)) object$ar1_rho_pair[i, j] else NA_real_
      est_ij <- if (!is.null(object$ar1_estimated)) object$ar1_estimated[i, j] else NA
      row$ar1_rho       <- round(rho_ij, digits)
      row$ar1_estimated <- est_ij
    }
    rows[[k]] <- row
  }
  df <- .mc_finalize_summary_df(
    do.call(rbind.data.frame, rows),
    class_name = "summary.ba_repeated_matrix",
    repeated = TRUE
  )
  if ("ar1_rho" %in% names(df) && all(is.na(df$ar1_rho))) {
    df$ar1_rho <- NULL
  }
  if ("ar1_estimated" %in% names(df) && all(is.na(df$ar1_estimated))) {
    df$ar1_estimated <- NULL
  }
  attr(df, "conf.level") <- cl
  df
}

.print_ba_summary_blocks <- function(x, ...) {
  .mc_print_sectioned_table(
    x,
    sections = list(
      list(
        title = "Agreement estimates",
        cols = c("item1", "item2", "n_obs", "bias", "sd_loa",
                 "loa_low", "loa_up", "width", "slope")
      ),
      list(
        title = "Confidence intervals",
        cols = c("bias_lwr", "bias_upr", "lo_lwr", "lo_upr", "up_lwr", "up_upr")
      ),
      list(
        title = "Model details",
        cols = c("sigma2_subject", "sigma2_resid",
                 "residual_model", "ar1_rho", "ar1_estimated")
      )
    ),
    header = NULL,
    ...
  )
}


#' @method print summary.ba_repeated
#' @export
print.summary.ba_repeated <- function(x, digits = NULL, n = NULL,
                                      topn = NULL, max_vars = NULL,
                                      width = NULL, show_ci = NULL, ...) {
  cl <- suppressWarnings(as.numeric(attr(x, "conf.level")))
  show_ci <- .mc_resolve_show_ci(show_ci, context = "summary")
  header <- .mc_header_with_ci("Bland-Altman (two methods)", cl, show_ci)
  .mc_print_sectioned_table(
    x,
    sections = list(
      list(
        title = "Agreement estimates",
        cols = c("item1", "item2", "n_obs", "bias", "sd_loa",
                 "loa_low", "loa_up", "width", "slope")
      ),
      list(
        title = "Confidence intervals",
        cols = c("bias_lwr", "bias_upr", "lo_lwr", "lo_upr", "up_lwr", "up_upr")
      ),
      list(
        title = "Model details",
        cols = c("sigma2_subject", "sigma2_resid",
                 "residual_model", "ar1_rho", "ar1_estimated")
      )
    ),
    header = header,
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

#' @method print summary.ba_repeated_matrix
#' @export
print.summary.ba_repeated_matrix <- function(x, digits = NULL, n = NULL,
                                             topn = NULL, max_vars = NULL,
                                             width = NULL, show_ci = NULL, ...) {
  cl <- suppressWarnings(as.numeric(attr(x, "conf.level")))
  show_ci <- .mc_resolve_show_ci(show_ci, context = "summary")
  header <- .mc_header_with_ci("Bland-Altman (pairwise)", cl, show_ci)
  .mc_print_sectioned_table(
    x,
    sections = list(
      list(
        title = "Agreement estimates",
        cols = c("item1", "item2", "n_obs", "bias", "sd_loa",
                 "loa_low", "loa_up", "width", "slope")
      ),
      list(
        title = "Confidence intervals",
        cols = c("bias_lwr", "bias_upr", "lo_lwr", "lo_upr", "up_lwr", "up_upr")
      ),
      list(
        title = "Model details",
        cols = c("sigma2_subject", "sigma2_resid",
                 "residual_model", "ar1_rho", "ar1_estimated")
      )
    ),
    header = header,
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

#-------- Plots ------------

#' @rdname ba_rm
#' @method plot ba_repeated
#' @param title Plot title (character scalar). Defaults to
#'   `"Bland-Altman (repeated measurements)"` for two methods and
#'   `"Bland-Altman (repeated, pairwise)"` for the faceted matrix plot.
#' @param subtitle Optional subtitle (character scalar). If `NULL`, a compact
#'   summary is shown using the fitted object.
#' @param point_alpha Numeric in `[0, 1]`. Transparency for scatter points
#'   drawn at (pair mean, pair difference) when point data are available.
#'   Passed to `ggplot2::geom_point(alpha = ...)`. Default `0.7`.
#'
#' @param point_size Positive numeric. Size of scatter points; passed to
#'   `ggplot2::geom_point(size = ...)`. Default `2.2`.
#'
#' @param line_size Positive numeric. Line width for horizontal bands
#'   (bias and both LoA) and, when requested, the proportional-bias line.
#'   Passed to `ggplot2::geom_hline(linewidth = ...)` (and `geom_abline`).
#'   Default `0.8`.
#'
#' @param shade_ci Logical. If `TRUE` and confidence intervals are available in
#'   the object (`CI.lines` for two methods; `*_ci_*` matrices for the pairwise
#'   case), semi-transparent rectangles are drawn to indicate CI bands for the
#'   bias and each LoA. If `FALSE`, dashed horizontal CI lines are drawn instead.
#'   Has no effect if CIs are not present. Default `TRUE`.
#'
#' @param shade_alpha Numeric in `[0, 1]`. Opacity of the CI shading
#'   rectangles when `shade_ci = TRUE`. Passed to `ggplot2::annotate("rect", alpha = ...)`.
#'   Default `0.08`.
#'
#' @param smoother One of `"none"`, `"loess"`, or `"lm"`. Adds an overlaid trend
#'   for differences vs means when points are drawn, to visualise proportional
#'   bias. `"lm"` fits a straight line with no SE ribbon; `"loess"` draws a
#'   locally-smoothed curve (span 0.9) with no SE ribbon; `"none"` draws no
#'   smoother. Ignored if `show_points = FALSE` or if no point data are available.
#'
#' @param symmetrize_y Logical (two-method plot only). If `TRUE`, the y-axis is
#'   centred at the estimated bias and expanded symmetrically to cover all
#'   elements used to compute the range (bands, CIs, and points if shown).
#'   Default `TRUE`.
#'
#' @param show_points Logical. If `TRUE`, per-pair points are drawn for the
#'   exactly-two-method path from the stored `means` and `diffs` vectors.
#'   Pairwise matrix plots draw points reconstructed from the canonical
#'   long-form columns stored by `ba_rm()`. If `FALSE` or if point data are
#'   unavailable, only the bands (and optional CI indicators) are drawn.
#'   Default `TRUE`.
#' @param ... Additional theme adjustments passed to `ggplot2::theme(...)`
#'   (e.g., `plot.title.position = "plot"`, `axis.title.x = element_text(size=11)`).
#' @param show_value Logical; included for a consistent plotting interface.
#'   Repeated-measures Bland-Altman plots do not overlay numeric cell values,
#'   so this argument currently has no effect.
#' @importFrom graphics abline lines par rect plot
#' @export
plot.ba_repeated <- function(x,
                             title = "Bland-Altman (repeated measurements)",
                             subtitle = NULL,
                             point_alpha = 0.7,
                             point_size  = 2.2,
                             line_size   = 0.8,
                             shade_ci    = TRUE,
                             shade_alpha = 0.08,
                             smoother    = c("none", "loess", "lm"),
                             symmetrize_y = TRUE,
                             show_points = TRUE,
                             show_value = TRUE,
                             ...) {
  check_inherits(x, "ba_repeated")
  check_bool(show_value, arg = "show_value")
  smoother <- match.arg(smoother)

  # Prefer pre-computed points if present
  has_pts <- !is.null(x$means) && !is.null(x$diffs)
  means <- if (has_pts) as.numeric(x$means) else numeric()
  diffs <- if (has_pts) as.numeric(x$diffs) else numeric()
  if (!has_pts) show_points <- FALSE

  # Bands
  md    <- as.numeric(x$mean.diffs)
  loaL  <- as.numeric(x$lower.limit)
  loaU  <- as.numeric(x$upper.limit)
  loa_multiplier <- as.numeric(x$loa_multiplier)
  n     <- as.integer(x$based.on)
  cl    <- suppressWarnings(as.numeric(attr(x, "conf.level")))
  ci_val <- function(nm) if (!is.null(x$CI.lines)) suppressWarnings(as.numeric(x$CI.lines[[nm]])) else NA_real_
  has_ci <- !is.null(x$CI.lines) &&
    all(is.finite(c(ci_val("mean.diff.ci.lower"), ci_val("mean.diff.ci.upper"),
                    ci_val("lower.limit.ci.lower"), ci_val("lower.limit.ci.upper"),
                    ci_val("upper.limit.ci.lower"), ci_val("upper.limit.ci.upper"))))

  if (is.null(subtitle)) {
    subtitle <- if (is.finite(cl)) {
      sprintf("pairs = %d |  mean diff = %.3g |  LoA = [%.3g, %.3g] |  %g%% CI",
              n, md, loaL, loaU, 100*cl)
    } else {
      sprintf("pairs = %d |  mean diff = %.3g |  LoA = [%.3g, %.3g]",
              n, md, loaL, loaU)
    }
  }

  # y-range
  y_bits <- c(loaL, loaU, md)
  if (has_ci) {
    y_bits <- c(y_bits,
                ci_val("mean.diff.ci.lower"), ci_val("mean.diff.ci.upper"),
                ci_val("lower.limit.ci.lower"), ci_val("lower.limit.ci.upper"),
                ci_val("upper.limit.ci.lower"), ci_val("upper.limit.ci.upper"))
  }
  if (show_points && length(diffs)) y_bits <- c(y_bits, diffs)
  y_rng <- range(y_bits, na.rm = TRUE)
  if (isTRUE(symmetrize_y) && all(is.finite(c(y_rng, md)))) {
    half <- max(abs(c(y_rng[1] - md, y_rng[2] - md)))
    y_rng <- c(md - half, md + half)
  }

  p <- ggplot2::ggplot()

  if (isTRUE(has_ci) && shade_ci) {
    p <- p +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                        ymin = ci_val("mean.diff.ci.lower"), ymax = ci_val("mean.diff.ci.upper"),
                        alpha = shade_alpha) +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                        ymin = ci_val("lower.limit.ci.lower"), ymax = ci_val("lower.limit.ci.upper"),
                        alpha = shade_alpha) +
      ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                        ymin = ci_val("upper.limit.ci.lower"), ymax = ci_val("upper.limit.ci.upper"),
                        alpha = shade_alpha)
  } else if (isTRUE(has_ci)) {
    p <- p + ggplot2::geom_hline(yintercept = c(
      ci_val("mean.diff.ci.lower"),  ci_val("mean.diff.ci.upper"),
      ci_val("lower.limit.ci.lower"),ci_val("lower.limit.ci.upper"),
      ci_val("upper.limit.ci.lower"),ci_val("upper.limit.ci.upper")),
      linetype = "dashed"
    )
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dotted", colour = "grey40") +
    ggplot2::geom_hline(yintercept = md,   linewidth = line_size) +
    ggplot2::geom_hline(yintercept = loaL, linewidth = line_size) +
    ggplot2::geom_hline(yintercept = loaU, linewidth = line_size)

  if (show_points && length(means) && length(diffs)) {
    df <- data.frame(means = means, diffs = diffs)
    p <- p + ggplot2::geom_point(data = df, ggplot2::aes(x = means, y = diffs),
                                 alpha = point_alpha, size = point_size)
    if (smoother == "lm") {
      p <- p + ggplot2::geom_smooth(data = df, ggplot2::aes(x = means, y = diffs),
                                    method = "lm", se = FALSE, linewidth = 0.7)
    } else if (smoother == "loess") {
      p <- p + ggplot2::geom_smooth(data = df, ggplot2::aes(x = means, y = diffs),
                                    method = "loess", se = FALSE, linewidth = 0.7, span = 0.9)
    }
    if (isTRUE(x$include_slope) && is.finite(x$beta_slope)) {
      p <- p + ggplot2::geom_abline(intercept = x$mean.diffs, slope = x$beta_slope, linewidth = line_size)
    }
  }

  p +
    ggplot2::coord_cartesian(ylim = y_rng, expand = TRUE) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), ...) +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = if (show_points) "Pair mean" else NULL,
                  y = "Difference (method 2 \u2212 method 1)")
}

#' @rdname ba_rm
#' @method plot ba_repeated_matrix
#' @param pairs (Faceted pairwise plot only.) Optional character vector of
#'   labels specifying which method contrasts to display. Labels must match the
#'   "row - column" convention used by `print()`/`summary()` (e.g., `"B - A"`).
#'   Defaults to all upper-triangle pairs.
#'
#' @param against (Faceted pairwise plot only.) Optional single method name.
#'   If supplied, facets are restricted to contrasts of the chosen method
#'   against all others. Ignored when `pairs` is provided.
#'
#' @param facet_scales (Faceted pairwise plot only.) Either `"free_y"` (default)
#'   to allow each facet its own y-axis limits, or `"fixed"` for a common scale
#'   across facets. Passed to `ggplot2::facet_wrap(scales = ...)`.
#' @importFrom ggplot2 ggplot annotate geom_hline geom_point geom_smooth
#' @importFrom ggplot2 coord_cartesian theme_minimal theme labs facet_wrap
#' @importFrom ggplot2 element_blank
#' @export
plot.ba_repeated_matrix <- function(
    x,
    pairs = NULL, against = NULL,
    facet_scales = c("free_y","fixed"),
    title = "Bland-Altman (repeated, pairwise)",
    point_alpha = 0.6, point_size = 1.8,
    line_size = 0.7, shade_ci = TRUE, shade_alpha = 0.08,
    smoother = c("none","loess","lm"),
    show_points = TRUE,
    show_value = TRUE,
    ...
) {
  check_inherits(x, "ba_repeated_matrix")
  check_bool(show_value, arg = "show_value")
  facet_scales <- match.arg(facet_scales)
  smoother <- match.arg(smoother)
  methods <- x$methods
  m <- length(methods)
  lab_pair <- function(i,j) paste(methods[i], "\u2212", methods[j])

  idx_upper <- which(upper.tri(matrix(NA_real_, m, m)), arr.ind = TRUE)
  all_pairs <- data.frame(
    j   = idx_upper[,1],
    k   = idx_upper[,2],
    lab = lab_pair(idx_upper[,1], idx_upper[,2]),
    stringsAsFactors = FALSE
  )

  if (!is.null(against)) {
    if (!against %in% methods) {
      abort_bad_arg("against",
        message = "must be one of: {methods*?{, }{ and }}.",
        methods = methods
      )
    }
    js <- match(against, methods)
    all_pairs <- subset(all_pairs, j == js | k == js)
  } else if (!is.null(pairs)) {
    all_pairs <- subset(all_pairs, lab %in% pairs)
    if (!nrow(all_pairs)) {
      abort_bad_arg("pairs",
        message = "None of the requested pairs matched available contrasts.",
        .hint = "Use the \"row - column\" labels shown in print/summary output."
      )
    }
  }
  pairs_order <- all_pairs$lab

  bands <- data.frame(
    pair = lab_pair(idx_upper[,1], idx_upper[,2]),
    md   = as.vector(x$bias[idx_upper]),
    loaL = as.vector(x$loa_lower[idx_upper]),
    loaU = as.vector(x$loa_upper[idx_upper]),
    md_l = as.vector(x$mean_ci_low[idx_upper]),
    md_u = as.vector(x$mean_ci_high[idx_upper]),
    lo_l = as.vector(x$loa_lower_ci_low[idx_upper]),
    lo_u = as.vector(x$loa_lower_ci_high[idx_upper]),
    hi_l = as.vector(x$loa_upper_ci_low[idx_upper]),
    hi_u = as.vector(x$loa_upper_ci_high[idx_upper]),
    stringsAsFactors = FALSE
  )
  bands <- bands[bands$pair %in% pairs_order, , drop = FALSE]
  bands$pair <- factor(bands$pair, levels = pairs_order)

  has_ci <- with(bands, all(is.finite(c(md_l, md_u, lo_l, lo_u, hi_l, hi_u))))

  # Build point data from stored mapping (no hard-coded names)
  pts <- NULL
  if (isTRUE(show_points) && !is.null(x$data_long) && !is.null(x$mapping)) {
    df  <- as.data.frame(x$data_long, check.names = FALSE)
    map <- x$mapping
    mY <- .resolve_map_col_ba(df, map, "response")
    mS <- .resolve_map_col_ba(df, map, "subject")
    mM <- .resolve_map_col_ba(df, map, "method")
    mT <- .resolve_map_col_ba(df, map, "time")

    build_panel <- function(lev_j, lev_k, lab) {
      a <- df[df[[mM]] == lev_j, c(mS, mT, mY)]
      b <- df[df[[mM]] == lev_k, c(mS, mT, mY)]
      names(a) <- c("subject","time","yA"); names(b) <- c("subject","time","yB")
      ab <- merge(a, b, by = c("subject","time"), all = FALSE)
      if (!nrow(ab)) return(NULL)
      data.frame(pair = lab,
                 means = (ab$yA + ab$yB)/2,
                 diffs =  ab$yB - ab$yA,
                 stringsAsFactors = FALSE)
    }

    pts <- do.call(rbind, lapply(seq_len(nrow(all_pairs)), function(i) {
      build_panel(methods[all_pairs$j[i]], methods[all_pairs$k[i]], all_pairs$lab[i])
    }))
    if (!is.null(pts) && nrow(pts)) {
      pts$pair <- factor(pts$pair, levels = pairs_order)
    } else {
      pts <- NULL
    }
  }

  p <- ggplot2::ggplot() + ggplot2::facet_wrap(~ pair, scales = facet_scales)

  # Ensure y-scale trained even without points
  p <- p +
    ggplot2::geom_blank(data = bands, ggplot2::aes(x = 0, y = md))   +
    ggplot2::geom_blank(data = bands, ggplot2::aes(x = 0, y = loaL)) +
    ggplot2::geom_blank(data = bands, ggplot2::aes(x = 0, y = loaU))

  if (has_ci) {
    if (shade_ci) {
      p <- p +
        ggplot2::geom_rect(data = bands, ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = md_l, ymax = md_u),
                           inherit.aes = FALSE, alpha = shade_alpha) +
        ggplot2::geom_rect(data = bands, ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = lo_l, ymax = lo_u),
                           inherit.aes = FALSE, alpha = shade_alpha) +
        ggplot2::geom_rect(data = bands, ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = hi_l, ymax = hi_u),
                           inherit.aes = FALSE, alpha = shade_alpha)
    } else {
      p <- p +
        ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = md_l)) +
        ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = md_u)) +
        ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = lo_l)) +
        ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = lo_u)) +
        ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = hi_l)) +
        ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = hi_u))
    }
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dotted", colour = "grey40") +
    ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = md),   linewidth = line_size) +
    ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = loaL), linewidth = line_size) +
    ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = loaU), linewidth = line_size)

  if (!is.null(pts)) {
    p <- p + ggplot2::geom_point(data = pts, ggplot2::aes(x = means, y = diffs),
                                 alpha = point_alpha, size = point_size)
    if (smoother == "lm") {
      p <- p + ggplot2::geom_smooth(data = pts, ggplot2::aes(x = means, y = diffs),
                                    method = "lm", se = FALSE, linewidth = 0.7)
    } else if (smoother == "loess") {
      p <- p + ggplot2::geom_smooth(data = pts, ggplot2::aes(x = means, y = diffs),
                                    method = "loess", se = FALSE, linewidth = 0.7, span = 0.9)
    }
  }

  p +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), ...) +
    ggplot2::labs(title = title,
                  x = if (!is.null(pts)) "Pair mean" else NULL,
                  y = "Difference (row \u2212 column)")
}

#' @keywords internal
num_or_na_ba <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y[!is.finite(y)] <- NA_real_
  y
}

#' @keywords internal
.resolve_map_col_ba <- function(df, mapping, key) {
  candidate <- mapping[[key]] %||% switch(key,
                                          response = ".response",
                                          subject  = ".subject",
                                          method   = ".method",
                                          time     = ".time"
  )
  if (candidate %in% names(df)) return(candidate)

  canonical <- switch(key,
                      response = ".response",
                      subject  = ".subject",
                      method   = ".method",
                      time     = ".time"
  )
  if (canonical %in% names(df)) return(canonical)

  abort_bad_arg("mapping",
    message = "Column for '{key}' not found in stored data_long.",
    key = key
  )
}

#' @keywords internal
.make_named_matrix_ba <- function(methods, fill = NA_real_, storage = "double") {
  m <- length(methods)
  out <- matrix(fill, m, m, dimnames = list(methods, methods))
  storage.mode(out) <- storage
  out
}
