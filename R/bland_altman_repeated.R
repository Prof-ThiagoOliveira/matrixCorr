#' @title Bland-Altman for repeated measurements
#' @description
#' Repeated-measures Bland–Altman (BA) for method comparison based on a
#' mixed-effects model fitted to **subject–time matched paired differences**.
#' A subject-specific random intercept accounts for clustering, and (optionally)
#' an AR(1) process captures serial correlation across replicates within subject.
#' The function returns bias (mean difference), limits of agreement (LoA),
#' confidence intervals, and variance components, for either two methods or
#' all pairwise contrasts when \eqn{\ge}3 methods are supplied.
#'
#' **Required columns / vectors**
#' \itemize{
#'   \item \code{response}: numeric measurements.
#'   \item \code{subject}: subject identifier (integer/factor/numeric).
#'   \item \code{method}: method label with \eqn{\ge}2 levels (factor/character/integer).
#'   \item \code{time}: replicate index used to form pairs; only records where
#'         both methods are present for the same \code{subject} and \code{time}
#'         contribute to a given pairwise BA.
#' }
#'
#' @param data Optional `data.frame`/`data.table` containing required columns.
#' @param response Numeric vector (stacked outcomes) **or** a single character
#'   string giving the column name in `data`.
#' @param subject Subject ID (integer/factor/numeric) **or** a single character
#'   string giving the column name in `data`.
#' @param method Method label (factor/character/integer; N >= 2 levels) **or**
#'   a single character string giving the column name in `data`.
#' @param time Integer/numeric replicate/time index (pairs within subject) **or**
#'   a single character string giving the column name in `data`.
#' @param two Positive scalar; LoA multiple of SD (default 1.96).
#' @param conf_level Confidence level for CIs (default 0.95).
#' @param include_slope Logical; if TRUE, estimates proportional bias per pair.
#' @param use_ar1 Logical; AR(1) within-subject residual correlation.
#' @param ar1_rho AR(1) parameter (|rho|<1).
#' @param max_iter,tol EM control for the backend (defaults 200, 1e-6).
#' @param verbose Logical; print brief progress.
#'
#' @return
#' Either a \code{"ba_repeated"} object (exactly two methods) or a
#' \code{"ba_repeated_matrix"} object (pairwise results when \eqn{\ge}3 methods).
#'
#' \strong{If \code{"ba_repeated_matrix"} (N\eqn{\ge}3 methods)}, outputs are:
#' \itemize{
#'   \item \code{bias} \eqn{(m \times m)}; estimated mean difference (row − column).
#'          Diagonal is \code{NA}.
#'   \item \code{sd_loa} \eqn{(m \times m)}; estimated SD of a \emph{single new} paired
#'         difference for the (row, column) methods, accounting for the random
#'         subject intercept and (if enabled) AR(1) correlation.
#'   \item \code{loa_lower}, \code{loa_upper} \eqn{(m \times m)}; limits of agreement
#'         for a single measurement pair, computed as
#'         \eqn{\mathrm{bias} \pm \mathrm{two}\times \mathrm{sd\_loa}}.
#'         Signs follow the row − column convention
#'         (e.g., \code{loa_lower[j,i] = -loa_upper[i,j]}).
#'   \item \code{width} \eqn{(m \times m)}; LoA width,
#'         \code{loa_upper − loa_lower} (= \code{2 * two * sd_loa}).
#'   \item \code{n} \eqn{(m \times m)}; number of subject–time pairs used in each
#'         pairwise BA (complete cases where both methods are present).
#'
#'   \item \strong{CI matrices at nominal \code{conf_level}} (delta method):
#'         \itemize{
#'           \item \code{mean_ci_low}, \code{mean_ci_high}; CI for the bias.
#'           \item \code{loa_lower_ci_low}, \code{loa_lower_ci_high}; CI for the lower LoA.
#'           \item \code{loa_upper_ci_low}, \code{loa_upper_ci_high}; CI for the upper LoA.
#'         }
#'
#'   \item \code{slope} (\eqn{m \times m}; optional); proportional-bias slope (difference
#'         vs pair mean) when \code{include_slope = TRUE}.
#'
#'   \item \code{sigma2_subject} \eqn{(m \times m)}; estimated variance of the
#'         subject-level random intercept (on differences).
#'   \item \code{sigma2_resid} \eqn{(m \times m)}; estimated residual variance of a
#'         single difference after accounting for the random intercept (and AR(1), if used).
#'
#'   \item \code{use_ar1} \emph{(scalar logical)}; whether AR(1) modeling was requested.
#'   \item \code{ar1_rho} \emph{(scalar numeric or \code{NA})}; user-supplied \eqn{\rho}
#'         if a single common value was provided; \code{NA} otherwise.
#'   \item \code{ar1_rho_pair} (\eqn{m \times m}; optional); \eqn{\rho} actually used per pair
#'         (may be estimated from data or equal to the supplied value).
#'   \item \code{ar1_estimated} (\eqn{m \times m}; optional logical); for each pair, \code{TRUE}
#'         if \eqn{\rho} was estimated internally; \code{FALSE} if supplied.
#'
#'   \item \code{methods} \emph{(character)}; method level names; matrix rows/columns
#'         follow this order.
#'   \item \code{two} \emph{(scalar)}; LoA multiplier used (default \code{1.96}).
#'   \item \code{conf_level} \emph{(scalar)}; nominal confidence level used for CIs.
#'
#'   \item \code{data_long} \emph{(data.frame)}; the long data used for fitting
#'         (\code{response}, \code{subject}, \code{method}, \code{time}). Included to
#'         facilitate plotting/reproducibility; not required for summary methods.
#' }
#'
#' \strong{If \code{"ba_repeated"} (exactly two methods)}, outputs are:
#' \itemize{
#'   \item \code{mean.diffs} \emph{(scalar)}; estimated bias (method 2 − method 1).
#'   \item \code{lower.limit}, \code{upper.limit} \emph{(scalars)}; LoA
#'         \eqn{\mu \pm \mathrm{two}\times \mathrm{SD}} for a single new pair.
#'   \item \code{critical.diff} \emph{(scalar)}; \code{two * SD}; LoA half-width.
#'   \item \code{two}, \code{conf_level} \emph{(scalars)}; as above.
#'   \item \code{CI.lines} \emph{(named numeric)}; CI bounds for bias and both LoA
#'         (\code{*.ci.lower}, \code{*.ci.upper}) at \code{conf_level}.
#'   \item \code{means}, \code{diffs} \emph{(vectors)}; per-pair means and differences
#'         used by plotting helpers.
#'   \item \code{based.on} \emph{(integer)}; number of subject–time pairs used.
#'   \item \code{include_slope}, \code{beta_slope}; whether a proportional-bias slope
#'         was estimated and its value (if requested).
#'   \item \code{sigma2_subject}, \code{sigma2_resid}; variance components as above.
#'   \item \code{use_ar1}, \code{ar1_rho}, \code{ar1_estimated}; AR(1) settings/results
#'         as above (scalars for the two-method fit).
#' }
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
#' ba4 <- bland_altman_repeated(
#'   data = data,
#'   response = "y", subject = "subject", method = "method", time = "time",
#'   two = 1.96, conf_level = 0.95,
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
#' baN <- bland_altman_repeated(
#'   response = data$y, subject = data$subject, method = data$method, time = data$time,
#'   two = 1.96, conf_level = 0.95,
#'   include_slope = TRUE,         # estimate proportional bias per pair
#'   use_ar1 = TRUE # model AR(1) within-subject
#' )
#'
#' # Matrices (row - column orientation)
#' print(baN)
#'
#' # Faceted BA scatter by pair
#' plot(baN, smoother = "lm", facet_scales = "free_y")
#'
#' # -------- Two-method path (A vs B only) -----------------------------------
#' data_AB <- subset(data, method %in% c("A","B"))
#' baAB <- bland_altman_repeated(
#'   response = data_AB$y, subject = data_AB$subject,
#'   method = droplevels(data_AB$method), time = data_AB$time,
#'   include_slope = TRUE, use_ar1 = TRUE, ar1_rho = 0.4
#' )
#' print(baAB)
#' plot(baAB)
#'
#' @author Thiago de Paula Oliveira
#' @export
bland_altman_repeated <- function(data = NULL, response, subject, method, time,
                                  two = 1.96, conf_level = 0.95,
                                  include_slope = FALSE,
                                  use_ar1 = FALSE, ar1_rho = NA_real_,
                                  max_iter = 200L, tol = 1e-6,
                                  verbose = FALSE) {

  # ---- resolve columns if 'data' provided and names given ----
  if (!is.null(data)) {
    if (!inherits(data, "data.frame"))
      stop("`data` must be a data.frame or data.table.")
    # helper to pull a column by name
    get_col <- function(x) {
      if (is.character(x) && length(x) == 1L) {
        if (!x %in% names(data))
          stop("Column '", x, "' not found in `data`.")
        data[[x]]
      } else {
        x
      }
    }
    response <- get_col(response)
    subject  <- get_col(subject)
    method   <- get_col(method)
    time     <- get_col(time)
  }

  # ---- original code (unchanged) ----
  stopifnot(is.numeric(response), length(response) > 1L)
  if (!(is.integer(subject) || is.factor(subject) || is.numeric(subject)))
    stop("'subject' must be integer/factor/numeric.")
  subject <- as.integer(as.factor(subject))

  # normalize method to factor with stable levels
  if (!is.factor(method)) method <- as.factor(method)
  method <- droplevels(method)
  mlev <- levels(method)
  if (length(mlev) < 2L) stop("Need at least 2 distinct methods in 'method'.")

  if (!(is.integer(time) || is.numeric(time))) stop("'time' must be integer/numeric.")
  time <- as.integer(time)
  if (!is.numeric(two) || length(two) != 1L || two <= 0) stop("'two' must be a positive scalar.")
  if (!(is.numeric(conf_level) && length(conf_level) == 1L && conf_level > 0 && conf_level < 1))
    stop("'conf_level' must be in (0,1).")
  if (isTRUE(use_ar1)) {
    if (!is.na(ar1_rho)) {
      if (!is.finite(ar1_rho)) stop("'ar1_rho' must be finite or NA when use_ar1=TRUE.")
      if (abs(ar1_rho) >= 0.999) stop("'ar1_rho' must be in (-0.999, 0.999).")
    }
  } else {
    ar1_rho <- NA_real_
  }

  # Two-method path
  if (length(mlev) == 2L) {
    idx <- method %in% mlev
    res1 <- .ba_rep_two_methods(
      response[idx], subject[idx], as.integer(method[idx] == mlev[1L]) + 1L,
      time[idx], two, conf_level, include_slope, use_ar1, ar1_rho, max_iter, tol
    )
    return(res1)
  }

  # N-method pairwise path
  if (isTRUE(verbose)) cat("Repeated BA (pairwise):", length(mlev), "methods ->",
                           choose(length(mlev), 2L), "pairs\n")

  methods <- mlev
  m <- length(methods)
  mk <- function(val = NA_real_, storage = "double") {
    out <- matrix(val, m, m)
    dimnames(out) <- list(methods, methods)
    storage.mode(out) <- storage
    out
  }
  bias         <- mk()
  sd_loa       <- mk()
  loa_lower    <- mk()
  loa_upper    <- mk()
  width        <- mk()
  n_mat        <- mk(NA_integer_, storage = "integer")
  mean_ci_low  <- mk()
  mean_ci_high <- mk()
  lo_ci_low    <- mk()
  lo_ci_high   <- mk()
  hi_ci_low    <- mk()
  hi_ci_high   <- mk()
  slope_mat    <- if (isTRUE(include_slope)) mk() else NULL
  vc_subject   <- mk()
  vc_resid     <- mk()
  ar1_rho_mat   <- if (isTRUE(use_ar1)) mk() else NULL
  ar1_estimated <- if (isTRUE(use_ar1)) {
    mtx <- matrix(NA, m, m, dimnames = list(methods, methods))
    storage.mode(mtx) <- "logical"
    mtx
  } else NULL

  # helper: from backend object LoA (±two*sd) and CIs at 'two'
  .recompose_pair <- function(fit) {
    list(
      md = as.numeric(fit$bias_mu0),
      sd = as.numeric(fit$sd_loa),
      lo = as.numeric(fit$loa_lower),
      hi = as.numeric(fit$loa_upper),
      bias_l = as.numeric(fit$bias_lwr),
      bias_u = as.numeric(fit$bias_upr),
      lo_l   = as.numeric(fit$loa_lower_lwr),
      lo_u   = as.numeric(fit$loa_lower_upr),
      hi_l   = as.numeric(fit$loa_upper_lwr),
      hi_u   = as.numeric(fit$loa_upper_upr)
    )
  }

  # loop pairs
  for (j in 1:(m-1)) for (k in (j+1):m) {
    lev_j <- methods[j]; lev_k <- methods[k]
    sel <- method %in% c(lev_j, lev_k) & is.finite(response) & !is.na(subject) & !is.na(time)
    if (!any(sel)) next
    m12 <- ifelse(method[sel] == lev_j, 1L, 2L)

    fit <- bland_altman_repeated_em_ext_cpp(
      y = response[sel], subject = subject[sel], method = m12, time = time[sel],
      include_slope = include_slope,
      use_ar1 = use_ar1, ar1_rho = ar1_rho,
      max_iter = max_iter, tol = tol, conf_level = conf_level,
      two_arg = two
    )

    comp <- .recompose_pair(fit)

    bias[j,k]      <- comp$md;     bias[k,j]      <- -comp$md
    sd_loa[j,k]    <- comp$sd;     sd_loa[k,j]    <- comp$sd
    loa_lower[j,k] <- comp$lo;     loa_lower[k,j] <- -comp$hi
    loa_upper[j,k] <- comp$hi;     loa_upper[k,j] <- -comp$lo
    width[j,k]     <- comp$hi - comp$lo; width[k,j] <- width[j,k]
    n_mat[j,k]     <- as.integer(fit$n_pairs); n_mat[k,j] <- n_mat[j,k]

    mean_ci_low[j,k]  <- comp$bias_l; mean_ci_low[k,j]  <- -comp$bias_u
    mean_ci_high[j,k] <- comp$bias_u; mean_ci_high[k,j] <- -comp$bias_l
    lo_ci_low[j,k]    <- comp$lo_l;   lo_ci_low[k,j]    <- -comp$hi_u
    lo_ci_high[j,k]   <- comp$lo_u;   lo_ci_high[k,j]   <- -comp$hi_l
    hi_ci_low[j,k]    <- comp$hi_l;   hi_ci_low[k,j]    <- -comp$lo_u
    hi_ci_high[j,k]   <- comp$hi_u;   hi_ci_high[k,j]   <- -comp$lo_l

    vc_subject[j,k]   <- vc_subject[k,j] <- as.numeric(fit$sigma2_subject)
    vc_resid[j,k]     <- vc_resid[k,j]   <- as.numeric(fit$sigma2_resid)

    if (isTRUE(use_ar1)) {
      ar1_rho_used <- as.numeric(fit$ar1_rho)
      ar1_rho_mat[j,k] <- ar1_rho_mat[k,j] <- ar1_rho_used
      ar1_estimated[j,k] <- ar1_estimated[k,j] <- isTRUE(fit$ar1_estimated)
    }

    if (!is.null(slope_mat)) {
      slope <- as.numeric(fit$beta_slope)
      slope_mat[j,k] <- slope
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
    two = two,
    conf_level = conf_level,
    use_ar1 = use_ar1,
    ar1_rho = if (use_ar1) ar1_rho else NA_real_,
    sigma2_subject = vc_subject,
    sigma2_resid   = vc_resid,
    ar1_rho_pair   = if (use_ar1) ar1_rho_mat else NULL,
    ar1_estimated  = if (use_ar1) ar1_estimated else NULL,
    data_long      = data.frame(
      response = response, subject = subject, method = method, time = time,
      check.names = FALSE
    )
  )
  class(ba_repeated) <- c("ba_repeated_matrix","list")
  attr(ba_repeated, "conf.level") <- conf_level
  ba_repeated
}



#' two-method helper
#' @keywords internal
.ba_rep_two_methods <- function(response, subject, method12, time,
                                two, conf_level, include_slope,
                                use_ar1, ar1_rho, max_iter, tol) {
  fit <- bland_altman_repeated_em_ext_cpp(
    y = response, subject = subject, method = method12, time = time,
    include_slope = include_slope,
    use_ar1 = use_ar1, ar1_rho = ar1_rho,
    max_iter = max_iter, tol = tol, conf_level = conf_level,
    two_arg = two
  )

  md  <- as.numeric(fit$bias_mu0)
  sdL <- as.numeric(fit$sd_loa)

  loa_lower <- as.numeric(fit$loa_lower)
  loa_upper <- as.numeric(fit$loa_upper)

  # CI lines straight from the backend
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
    groups        = data.frame(subject = subject[!is.na(response)], time = time[!is.na(response)])[seq_along(means), , drop = FALSE],
    based_on      = as.integer(fit$n_pairs),  # (kept as-is if you already used 'based.on')
    based.on      = as.integer(fit$n_pairs),
    lower.limit   = loa_lower,
    mean.diffs    = md,
    upper.limit   = loa_upper,
    lines         = c(lower = loa_lower, mean = md, upper = loa_upper),
    CI.lines      = CI.lines,
    two           = two,
    critical.diff = two * sdL,
    conf_level    = conf_level,
    include_slope = include_slope,
    beta_slope    = if (include_slope) as.numeric(fit$beta_slope) else NA_real_,
    sigma2_subject= as.numeric(fit$sigma2_subject),
    sigma2_resid  = as.numeric(fit$sigma2_resid),
    use_ar1       = use_ar1,
    ar1_rho       = if (use_ar1) as.numeric(fit$ar1_rho) else NA_real_,
    ar1_estimated = if (use_ar1) isTRUE(fit$ar1_estimated) else NA
  )
  class(ba_repeated) <- c("ba_repeated","list")
  attr(ba_repeated, "conf.level") <- conf_level
  ba_repeated
}


# ---------- printers ----------
#' @rdname bland_altman_repeated
#' @method print ba_repeated
#' @param x A \code{"ba_repeated"} object.
#' @param digits Number of digits for estimates (default 3).
#' @param ci_digits Number of digits for CI bounds (default 3).
#' @param ... Unused.
#' @export
print.ba_repeated <- function(x, digits = 3, ci_digits = 3, ...) {
  stopifnot(inherits(x, "ba_repeated"))
  n   <- as.integer(x$based.on)
  two <- as.numeric(x$two)
  cl  <- suppressWarnings(as.numeric(attr(x, "conf.level")))
  if (!is.finite(cl)) cl <- NA_real_

  md    <- as.numeric(x$mean.diffs)
  loaL  <- as.numeric(x$lower.limit)
  loaU  <- as.numeric(x$upper.limit)
  sd_d  <- as.numeric(x$critical.diff) / two
  width <- loaU - loaL

  cil <- function(nm) as.numeric(x$CI.lines[[nm]])
  bias_l <- cil("mean.diff.ci.lower"); bias_u <- cil("mean.diff.ci.upper")
  lo_l   <- cil("lower.limit.ci.lower"); lo_u <- cil("lower.limit.ci.upper")
  hi_l   <- cil("upper.limit.ci.lower"); hi_u <- cil("upper.limit.ci.upper")

  head <- sprintf("Repeated-measures Bland-Altman (pairs = %d) - LoA = mean +/- %.3g * SD%s\n\n",
                  n, two, if (is.finite(cl)) sprintf(", %g%% CI", 100*cl) else "")
  cat(head)

  df <- data.frame(
    quantity = c("Mean difference", "Lower LoA", "Upper LoA"),
    estimate = c(md, loaL, loaU),
    lwr      = c(bias_l, lo_l, hi_l),
    upr      = c(bias_u, lo_u, hi_u),
    check.names = FALSE,
    row.names   = NULL
  )
  df$estimate <- formatC(df$estimate, format = "f", digits = digits)
  df$lwr      <- formatC(df$lwr,      format = "f", digits = ci_digits)
  df$upr      <- formatC(df$upr,      format = "f", digits = ci_digits)
  print(df, row.names = FALSE, right = FALSE)

  cat(sprintf("\nSD(single-pair differences): %s   LoA width: %s",
              formatC(sd_d, format = "f", digits = digits),
              formatC(width, format = "f", digits = digits)))
  if (isTRUE(x$include_slope) && is.finite(x$beta_slope)) {
    cat(sprintf("\nProportional bias slope (vs pair mean): %s\n",
                formatC(x$beta_slope, format = "f", digits = digits)))
  } else {
    cat("\n")
  }
  invisible(x)
}

#' @export
print.ba_repeated_matrix <- function(x,
                                     digits = 3,
                                     ci_digits = 3,
                                     style = c("pairs","matrices"),
                                     ...) {
  stopifnot(inherits(x, "ba_repeated_matrix"))
  style <- match.arg(style)

  # Back-compat: allow old 'show' argument to force matrix view
  dots <- list(...)
  if (!is.null(dots$show)) style <- "matrices"

  cl <- suppressWarnings(as.numeric(attr(x, "conf.level")))
  has_ci <- all(c("mean_ci_low","mean_ci_high",
                  "loa_lower_ci_low","loa_lower_ci_high",
                  "loa_upper_ci_low","loa_upper_ci_high") %in% names(x))

  if (style == "matrices") {
    # old behavior (unchanged)
    show <- if (is.null(dots$show)) c("bias","width","sd_loa") else dots$show
    NextMethod("print", x, digits = digits, show = show, ...)
    return(invisible(x))
  }

  methods <- x$methods
  m <- length(methods)
  rows <- vector("list", m * (m-1) / 2L)
  k <- 0L
  for (i in 1:(m-1)) for (j in (i+1):m) {
    k <- k + 1L
    bias <- x$bias[i,j]
    sd   <- x$sd_loa[i,j]
    lo   <- x$loa_lower[i,j]
    hi   <- x$loa_upper[i,j]
    w    <- x$width[i,j]

    row <- list(
      method1 = methods[i],
      method2 = methods[j],
      bias    = round(bias, digits),
      sd_loa  = round(sd,   digits),
      loa_low = round(lo,   digits),
      loa_up  = round(hi,   digits),
      width   = round(w,    digits),
      n       = suppressWarnings(as.integer(x$n[i,j]))
    )

    if (has_ci && is.finite(cl)) {
      row$bias_lwr <- round(x$mean_ci_low[i,j],        ci_digits)
      row$bias_upr <- round(x$mean_ci_high[i,j],       ci_digits)
      row$lo_lwr   <- round(x$loa_lower_ci_low[i,j],   ci_digits)
      row$lo_upr   <- round(x$loa_lower_ci_high[i,j],  ci_digits)
      row$up_lwr   <- round(x$loa_upper_ci_low[i,j],   ci_digits)
      row$up_upr   <- round(x$loa_upper_ci_high[i,j],  ci_digits)
    }
    rows[[k]] <- row
  }

  df <- do.call(rbind.data.frame, rows)
  if (is.finite(cl)) cat(sprintf("Bland–Altman (row − column), %g%% CI\n\n", 100*cl))
  else               cat("Bland–Altman (row − column)\n\n")
  print(df, row.names = FALSE, right = FALSE)
  invisible(x)
}


#' @method summary ba_repeated
#' @export
summary.ba_repeated <- function(object,
                                digits = 3,
                                ci_digits = 3,
                                ...) {
  stopifnot(inherits(object, "ba_repeated"))
  cl <- suppressWarnings(as.numeric(attr(object, "conf.level")))
  n  <- as.integer(object$based.on)

  out <- data.frame(
    method1   = "method 1",
    method2   = "method 2",
    bias      = round(num_or_na(object$mean.diffs), digits),
    sd_loa    = round(num_or_na(object$critical.diff) / num_or_na(object$two), digits),
    loa_low   = round(num_or_na(object$lower.limit), digits),
    loa_up    = round(num_or_na(object$upper.limit), digits),
    width     = round(num_or_na(object$upper.limit - object$lower.limit), digits),
    n         = n,
    stringsAsFactors = FALSE, check.names = FALSE
  )

  cil <- function(nm) num_or_na(object$CI.lines[[nm]])
  if (!any(is.na(c(cil("mean.diff.ci.lower"), cil("mean.diff.ci.upper"),
                   cil("lower.limit.ci.lower"), cil("lower.limit.ci.upper"),
                   cil("upper.limit.ci.lower"), cil("upper.limit.ci.upper"))))) {
    out$bias_lwr <- round(cil("mean.diff.ci.lower"), ci_digits)
    out$bias_upr <- round(cil("mean.diff.ci.upper"), ci_digits)
    out$lo_lwr   <- round(cil("lower.limit.ci.lower"), ci_digits)
    out$lo_upr   <- round(cil("lower.limit.ci.upper"), ci_digits)
    out$up_lwr   <- round(cil("upper.limit.ci.lower"), ci_digits)
    out$up_upr   <- round(cil("upper.limit.ci.upper"), ci_digits)
  }

  out$sigma2_subject <- round(num_or_na(object$sigma2_subject), digits)
  out$sigma2_resid   <- round(num_or_na(object$sigma2_resid),   digits)
  out$use_ar1        <- isTRUE(object$use_ar1)
  out$ar1_rho        <- if (isTRUE(object$use_ar1)) round(num_or_na(object$ar1_rho), digits) else NA_real_
  out$ar1_estimated  <- if (isTRUE(object$use_ar1)) isTRUE(object$ar1_estimated) else NA

  attr(out, "conf.level") <- cl
  class(out) <- c("summary.ba_repeated","data.frame")
  out
}

#' @method summary ba_repeated_matrix
#' @export
summary.ba_repeated_matrix <- function(object,
                                       digits = 3,
                                       ci_digits = 3,
                                       ...) {
  stopifnot(inherits(object, "ba_repeated_matrix"))
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
      n       = suppressWarnings(as.integer(object$n[i, j]))
    )
    if (has_ci && is.finite(cl)) {
      row$bias_lwr <- round(object$mean_ci_low[i, j],       ci_digits)
      row$bias_upr <- round(object$mean_ci_high[i, j],      ci_digits)
      row$lo_lwr   <- round(object$loa_lower_ci_low[i, j],  ci_digits)
      row$lo_upr   <- round(object$loa_lower_ci_high[i, j], ci_digits)
      row$up_lwr   <- round(object$loa_upper_ci_low[i, j],  ci_digits)
      row$up_upr   <- round(object$loa_upper_ci_high[i, j], ci_digits)
    }

    row$sigma2_subject <- round(object$sigma2_subject[i, j], digits)
    row$sigma2_resid   <- round(object$sigma2_resid[i, j],   digits)
    if (isTRUE(object$use_ar1)) {
      rho_ij <- if (!is.null(object$ar1_rho_pair)) object$ar1_rho_pair[i, j] else NA_real_
      est_ij <- if (!is.null(object$ar1_estimated)) isTRUE(object$ar1_estimated[i, j]) else NA
      row$ar1_rho       <- round(rho_ij, digits)
      row$ar1_estimated <- est_ij
    }

    rows[[k]] <- row
  }
  df <- do.call(rbind.data.frame, rows)
  attr(df, "conf.level") <- cl
  class(df) <- c("summary.ba_repeated_matrix", "data.frame")
  df
}


#' @method print summary.ba_repeated
#' @export
print.summary.ba_repeated <- function(x, ...) {
  cl <- suppressWarnings(as.numeric(attr(x, "conf.level")))
  if (is.finite(cl)) cat(sprintf("Bland–Altman (two methods), %g%% CI\n\n", 100*cl))
  else               cat("Bland–Altman (two methods)\n\n")
  print.data.frame(x, row.names = FALSE, right = FALSE, ...)
  invisible(x)
}

#' @method print summary.ba_repeated_matrix
#' @export
print.summary.ba_repeated_matrix <- function(x, ...) {
  cl <- suppressWarnings(as.numeric(attr(x, "conf.level")))
  if (is.finite(cl)) cat(sprintf("Bland–Altman (pairwise), %g%% CI\n\n", 100*cl))
  else               cat("Bland–Altman (pairwise)\n\n")
  print.data.frame(x, row.names = FALSE, right = FALSE, ...)
  invisible(x)
}

#-------- Plots ------------

#' @rdname bland_altman_repeated
#' @method plot ba_repeated
#' @param title Plot title.
#' @param subtitle Optional subtitle; if NULL, shows n and LoA summary.
#' @param point_alpha,point_size,line_size,shade_ci,shade_alpha Same semantics as \code{plot.ba()}.
#' @param smoother One of "none","loess","lm" to visualize proportional bias (only if points shown).
#' @param symmetrize_y Logical; if TRUE, y-axis centered at bias with symmetric limits.
#' @param show_points Logical; if FALSE or if no point data stored, draw bands-only.
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
                             ...) {
  stopifnot(inherits(x, "ba_repeated"))
  smoother <- match.arg(smoother)

  # Stored point data (may be absent in some objects)
  has_pts <- !is.null(x$means) && !is.null(x$diffs)
  if (!has_pts) show_points <- FALSE

  means <- if (has_pts) as.numeric(x$means) else numeric()
  diffs <- if (has_pts) as.numeric(x$diffs) else numeric()

  # Bands
  md    <- as.numeric(x$mean.diffs)
  loaL  <- as.numeric(x$lower.limit)
  loaU  <- as.numeric(x$upper.limit)
  two   <- as.numeric(x$two)
  n     <- as.integer(x$based.on)
  sd_d  <- as.numeric(x$critical.diff) / two
  cl    <- suppressWarnings(as.numeric(attr(x, "conf.level")))
  ci_val <- function(nm) if (!is.null(x$CI.lines)) suppressWarnings(as.numeric(x$CI.lines[[nm]])) else NA_real_
  has_ci <- !is.null(x$CI.lines) &&
    all(is.finite(c(ci_val("mean.diff.ci.lower"), ci_val("mean.diff.ci.upper"),
                    ci_val("lower.limit.ci.lower"), ci_val("lower.limit.ci.upper"),
                    ci_val("upper.limit.ci.lower"), ci_val("upper.limit.ci.upper"))))

  if (is.null(subtitle)) {
    subtitle <- if (is.finite(cl)) {
      sprintf("pairs = %d  •  mean diff = %.2f  •  LoA = [%.2f, %.2f]  •  %g%% CI",
              n, md, loaL, loaU, 100*cl)
    } else {
      sprintf("pairs = %d  •  mean diff = %.2f  •  LoA = [%.2f, %.2f]",
              n, md, loaL, loaU)
    }
  }

  # Choose y-range: prefer points if present; else from bands/CI
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
  # ggplot
  p <- ggplot()

  # CI shading or dashed bands
  if (isTRUE(has_ci) && shade_ci) {
    p <- p +
      annotate("rect", xmin = -Inf, xmax = Inf,
               ymin = ci_val("mean.diff.ci.lower"), ymax = ci_val("mean.diff.ci.upper"),
               alpha = shade_alpha) +
      annotate("rect", xmin = -Inf, xmax = Inf,
               ymin = ci_val("lower.limit.ci.lower"), ymax = ci_val("lower.limit.ci.upper"),
               alpha = shade_alpha) +
      annotate("rect", xmin = -Inf, xmax = Inf,
               ymin = ci_val("upper.limit.ci.lower"), ymax = ci_val("upper.limit.ci.upper"),
               alpha = shade_alpha)
  } else if (isTRUE(has_ci)) {
    p <- p + geom_hline(yintercept = c(
      ci_val("mean.diff.ci.lower"),  ci_val("mean.diff.ci.upper"),
      ci_val("lower.limit.ci.lower"),ci_val("lower.limit.ci.upper"),
      ci_val("upper.limit.ci.lower"),ci_val("upper.limit.ci.upper")),
      linetype = "dashed"
    )
  }

  # Zero and BA bands
  p <- p +
    geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dotted", color = "grey40") +
    geom_hline(yintercept = md,   linewidth = line_size) +
    geom_hline(yintercept = loaL, linewidth = line_size) +
    geom_hline(yintercept = loaU, linewidth = line_size)

  # Points + optional smoother (only if we have them and user wants)
  if (show_points && length(means) && length(diffs)) {
    df <- data.frame(means = means, diffs = diffs)
    p <- p + geom_point(data = df, aes(x = means, y = diffs),
                        alpha = point_alpha, size = point_size)
    if (smoother == "lm") {
      p <- p + geom_smooth(data = df, aes(x = means, y = diffs),
                           method = "lm", se = FALSE, linewidth = 0.7)
    } else if (smoother == "loess") {
      p <- p + geom_smooth(data = df, aes(x = means, y = diffs),
                           method = "loess", se = FALSE, linewidth = 0.7, span = 0.9)
    }
    # Optional proportional-bias line (requires x support)
    if (isTRUE(x$include_slope) && is.finite(x$beta_slope)) {
      p <- p + geom_abline(intercept = x$mean.diffs, slope = x$beta_slope, linewidth = line_size)
    }
  }

  p +
    coord_cartesian(ylim = y_rng, expand = TRUE) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), ...) +
    labs(title = title, subtitle = subtitle,
         x = if (show_points) "Pair mean" else NULL,
         y = "Difference (method 2 - method 1)")
}


#' Faceted Bland-Altman plot for repeated (pairwise, NO external data)
#' @param x A "ba_repeated_matrix" from bland_altman_repeated().
#' @param pairs Optional character vector of "row − column" labels to plot; default all upper-tri pairs.
#' @param against Optional single method to plot against all others.
#' @param facet_scales "free_y" (default) or "fixed".
#' @export
#' Faceted Bland-Altman plot for repeated (pairwise, NO external data)
#' @param x A "ba_repeated_matrix" from bland_altman_repeated().
#' @param pairs Optional character vector of "row − column" labels to plot; default all upper-tri pairs.
#' @param against Optional single method to plot against all others.
#' @param facet_scales "free_y" (default) or "fixed".
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
    ...
) {
  stopifnot(inherits(x, "ba_repeated_matrix"))
  facet_scales <- match.arg(facet_scales)
  smoother <- match.arg(smoother)
  methods <- x$methods
  m <- length(methods)
  lab_pair <- function(i,j) paste(methods[i], "\u2212", methods[j])

  ## ---------- build pair list (upper triangle) ----------
  idx_upper <- which(upper.tri(matrix(NA_real_, m, m)), arr.ind = TRUE)
  all_pairs <- data.frame(
    j   = idx_upper[,1],
    k   = idx_upper[,2],
    lab = lab_pair(idx_upper[,1], idx_upper[,2]),
    stringsAsFactors = FALSE
  )

  if (!is.null(against)) {
    if (!against %in% methods)
      stop("`against` must be one of: ", paste(methods, collapse=", "))
    js <- match(against, methods)
    all_pairs <- subset(all_pairs, j == js | k == js)
  } else if (!is.null(pairs)) {
    all_pairs <- subset(all_pairs, lab %in% pairs)
    if (!nrow(all_pairs)) stop("None of requested `pairs` matched.")
  }

  pairs_order <- all_pairs$lab

  ## ---------- bands per panel from x ----------
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

  ## ---------- point data (optional) ----------
  pts <- NULL
  if (isTRUE(show_points) && !is.null(x$data_long)) {
    df <- as.data.frame(x$data_long, check.names = FALSE)

    # --- robust column detection (y/response) ---
    subj_col <- "subject"
    meth_col <- "method"
    time_col <- "time"
    if (!all(c(subj_col, meth_col, time_col) %in% names(df))) {
      stop("`data_long` must contain columns 'subject', 'method', and 'time'.")
    }
    y_col <- if ("y" %in% names(df)) {
      "y"
    } else if ("response" %in% names(df)) {
      "response"
    } else {
      # last resort: pick the first numeric column not among the keys
      cand <- setdiff(names(df), c(subj_col, meth_col, time_col))
      cand_num <- cand[vapply(df[cand], is.numeric, logical(1))]
      if (length(cand_num) == 0L)
        stop("Could not find outcome column in `data_long` (expected 'y' or 'response').")
      cand_num[1]
    }

    build_panel <- function(lev_j, lev_k, lab) {
      a <- df[df[[meth_col]] == lev_j, c(subj_col, time_col, y_col)]
      b <- df[df[[meth_col]] == lev_k, c(subj_col, time_col, y_col)]
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

  ## ---------- ggplot ----------
  p <- ggplot2::ggplot() + ggplot2::facet_wrap(~ pair, scales = facet_scales)

  # Ensure y-scale is trained even without points
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

  # Zero and BA bands
  p <- p +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dotted", colour = "grey40") +
    ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = md),   linewidth = line_size) +
    ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = loaL), linewidth = line_size) +
    ggplot2::geom_hline(data = bands, ggplot2::aes(yintercept = loaU), linewidth = line_size)

  # Optional points + smoother
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
    ggplot2::labs(title = title, x = if (!is.null(pts)) "Pair mean" else NULL,
                  y = "Difference (row - column)")
}
