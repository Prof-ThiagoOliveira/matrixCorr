# Internal helpers for the Poisson GLMM implementation.
.mc_ccc_glmm_ghq_rule <- function(n = 40L) {
  n <- as.integer(n)
  i <- seq_len(n - 1L)
  J <- matrix(0, n, n)
  off <- sqrt(i / 2)
  J[cbind(i, i + 1L)] <- off
  J[cbind(i + 1L, i)] <- off
  eig <- eigen(J, symmetric = TRUE)
  ord <- order(eig$values)
  nodes <- eig$values[ord]
  weights <- sqrt(pi) * eig$vectors[1L, ord]^2
  list(nodes = nodes, weights = weights)
}

.mc_ccc_glmm_replicates_per_cell <- function(subject, method, replicate) {
  if (!length(replicate)) return(1L)
  tab <- table(subject, method)
  cells <- as.numeric(tab[tab > 0])
  if (!length(cells)) return(1L)
  max(1L, as.integer(round(mean(cells))))
}

.mc_ccc_glmm_vcov_from_hessian <- function(hessian) {
  if (is.null(hessian)) return(NULL)
  hessian <- 0.5 * (hessian + t(hessian))
  if (any(!is.finite(hessian))) return(NULL)

  eig <- tryCatch(
    eigen(hessian, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NULL
  )
  if (is.null(eig) ||
      any(!is.finite(eig)) ||
      min(eig) <= sqrt(.Machine$double.eps)) {
    return(NULL)
  }

  tryCatch(
    chol2inv(chol(hessian)),
    error = function(e) NULL
  )
}

.mc_num_grad <- function(f, x, rel_step = 1e-5) {
  g <- numeric(length(x))
  for (i in seq_along(x)) {
    h <- rel_step * max(abs(x[[i]]), 1)
    xp <- x
    xm <- x
    xp[[i]] <- xp[[i]] + h
    xm[[i]] <- xm[[i]] - h
    fp <- f(xp)
    fm <- f(xm)
    if (!is.finite(fp) || !is.finite(fm)) {
      g[[i]] <- NA_real_
    } else {
      g[[i]] <- (fp - fm) / (2 * h)
    }
  }
  g
}

.mc_delta_ci <- function(estimate,
                         f,
                         theta,
                         vcov_theta,
                         conf_level = 0.95,
                         transform = c("fisher", "wald")) {
  transform <- match.arg(transform)
  out <- list(
    estimate = estimate,
    se = NA_real_,
    lwr.ci = NA_real_,
    upr.ci = NA_real_,
    conf_level = conf_level,
    ci_method = paste0("delta_", transform)
  )

  if (is.null(vcov_theta) ||
      any(!is.finite(vcov_theta)) ||
      !is.finite(estimate)) {
    return(out)
  }

  grad <- .mc_num_grad(f, theta)
  if (any(!is.finite(grad))) return(out)

  var_est <- as.numeric(t(grad) %*% vcov_theta %*% grad)
  if (!is.finite(var_est) || var_est < 0) return(out)

  se <- sqrt(var_est)
  alpha <- 1 - conf_level
  zcrit <- stats::qnorm(1 - alpha / 2)

  if (identical(transform, "wald")) {
    lwr <- estimate - zcrit * se
    upr <- estimate + zcrit * se
  } else {
    estimate_safe <- max(min(estimate, 1 - 1e-12), -1 + 1e-12)
    z_est <- atanh(estimate_safe)
    se_z <- se / (1 - estimate_safe^2)
    lwr <- tanh(z_est - zcrit * se_z)
    upr <- tanh(z_est + zcrit * se_z)
  }

  out$se <- se
  out$lwr.ci <- lwr
  out$upr.ci <- upr
  out
}

.mc_ccc_glmm_delta_functions <- function(include_subject_method,
                                         phi,
                                         m_reps) {
  clamp01 <- function(x) pmax(0, pmin(1, x))

  rho_total <- function(theta) {
    beta0 <- theta[[1L]]
    beta_method <- theta[[2L]]
    sigma2_subject <- exp(2 * theta[[3L]])
    sigma2_subject_method <- if (isTRUE(include_subject_method)) exp(2 * theta[[4L]]) else 0
    sigma2_method <- 0.5 * beta_method * beta_method
    sigma2_total <- sigma2_subject + sigma2_method + sigma2_subject_method
    mu <- exp(beta0 + 0.5 * beta_method + 0.5 * sigma2_total)
    clamp01(
      mu * (exp(sigma2_subject) - 1) /
        (mu * (exp(sigma2_total) - 1) + phi)
    )
  }

  rho_inter <- function(theta) {
    beta0 <- theta[[1L]]
    beta_method <- theta[[2L]]
    sigma2_subject <- exp(2 * theta[[3L]])
    sigma2_subject_method <- if (isTRUE(include_subject_method)) exp(2 * theta[[4L]]) else 0
    sigma2_method <- 0.5 * beta_method * beta_method
    sigma2_total <- sigma2_subject + sigma2_method + sigma2_subject_method
    mu <- exp(beta0 + 0.5 * beta_method + 0.5 * sigma2_total)
    clamp01(
      mu * (exp(sigma2_subject) - 1) /
        (mu * (exp(sigma2_total) - 1) + phi / m_reps)
    )
  }

  rho_intra_method1 <- function(theta) {
    beta0 <- theta[[1L]]
    sigma2_subject <- exp(2 * theta[[3L]])
    sigma2_subject_method <- if (isTRUE(include_subject_method)) exp(2 * theta[[4L]]) else 0
    sigma2_intra <- sigma2_subject + sigma2_subject_method
    mu1 <- exp(beta0 + 0.5 * sigma2_intra)
    clamp01(
      mu1 * (exp(sigma2_intra) - 1) /
        (mu1 * (exp(sigma2_intra) - 1) + phi)
    )
  }

  rho_intra_method2 <- function(theta) {
    beta0 <- theta[[1L]]
    beta_method <- theta[[2L]]
    sigma2_subject <- exp(2 * theta[[3L]])
    sigma2_subject_method <- if (isTRUE(include_subject_method)) exp(2 * theta[[4L]]) else 0
    sigma2_intra <- sigma2_subject + sigma2_subject_method
    mu2 <- exp(beta0 + beta_method + 0.5 * sigma2_intra)
    clamp01(
      mu2 * (exp(sigma2_intra) - 1) /
        (mu2 * (exp(sigma2_intra) - 1) + phi)
    )
  }

  list(
    rho_ccc = rho_total,
    rho_ccc_inter = rho_inter,
    rho_ccc_intra_method1 = rho_intra_method1,
    rho_ccc_intra_method2 = rho_intra_method2
  )
}

.mc_ccc_glmm_fit <- function(y,
                             subject,
                             method,
                             n_subjects,
                             include_subject_method = FALSE,
                             max_iter = 1000L,
                             tol = 1e-8,
                             compute_vcov = FALSE) {
  gh <- .mc_ccc_glmm_ghq_rule(if (isTRUE(include_subject_method)) 11L else 40L)
  mean_a <- mean(y[method == 1L])
  mean_b <- mean(y[method == 2L])
  beta0_start <- log(mean_a + 0.1)
  beta_method_start <- log(mean_b + 0.1) - log(mean_a + 0.1)
  subject_means <- tapply(y, subject, mean)
  raw_subject_var <- stats::var(log(subject_means + 0.1))
  sigma_start <- sqrt(max(raw_subject_var, 0.05, na.rm = TRUE))
  sm_start <- if (isTRUE(include_subject_method)) {
    subject_method <- interaction(subject, method, drop = TRUE)
    sm_means <- tapply(y, subject_method, mean)
    sqrt(max(stats::var(log(sm_means + 0.1)) - raw_subject_var, 0.02, na.rm = TRUE))
  } else {
    NULL
  }

  starts <- if (isTRUE(include_subject_method)) {
    list(
      c(beta0_start, beta_method_start, log(sigma_start), log(sm_start)),
      c(beta0_start, beta_method_start, log(0.25), log(0.20)),
      c(beta0_start, beta_method_start, log(0.50), log(0.25)),
      c(beta0_start, beta_method_start, log(1.00), log(0.30)),
      c(beta0_start, beta_method_start, log(0.50), log(0.50)),
      c(log(mean(y) + 0.1), 0, log(0.50), log(0.25))
    )
  } else {
    list(
      c(beta0_start, beta_method_start, log(sigma_start)),
      c(beta0_start, beta_method_start, log(0.25)),
      c(beta0_start, beta_method_start, log(0.50)),
      c(beta0_start, beta_method_start, log(1.00)),
      c(beta0_start, beta_method_start, log(1.50)),
      c(beta0_start, beta_method_start, log(2.00)),
      c(log(mean(y) + 0.1), 0, log(0.50))
    )
  }

  objective <- function(par) {
    ccc_glmm_poisson_ghq_nll_cpp(
      par = par,
      y = y,
      subject = subject,
      method_code = method,
      n_subjects = n_subjects,
      include_subject_method = include_subject_method,
      gh_nodes = gh$nodes,
      gh_weights = gh$weights
    )
  }

  fits <- lapply(starts, function(start) {
    tryCatch(
      stats::optim(
        par = start,
        fn = objective,
        method = "BFGS",
        control = list(maxit = max_iter, reltol = tol)
      ),
      error = function(e) NULL
    )
  })

  fits <- Filter(Negate(is.null), fits)
  fits <- Filter(function(z) is.finite(z$value), fits)
  if (!length(fits)) return(NULL)

  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  sigma <- exp(best$par[[3L]])
  sigma_sm <- if (isTRUE(include_subject_method)) exp(best$par[[4L]]) else 0
  if (!identical(best$convergence, 0L) ||
      !is.finite(sigma) ||
      sigma <= 1e-6 ||
      (isTRUE(include_subject_method) && (!is.finite(sigma_sm) || sigma_sm <= 1e-6))) {
    return(NULL)
  }

  hessian <- NULL
  vcov_par <- NULL
  if (isTRUE(compute_vcov)) {
    hessian <- tryCatch(
      stats::optimHess(best$par, objective),
      error = function(e) NULL
    )
    vcov_par <- .mc_ccc_glmm_vcov_from_hessian(hessian)
  }

  list(
    par = best$par,
    value = best$value,
    converged = TRUE,
    convergence = best$convergence,
    message = best$message %||% "",
    counts = best$counts,
    gh = gh,
    include_subject_method = include_subject_method,
    hessian = hessian,
    vcov_par = vcov_par
  )
}

.mc_ccc_glmm_random_effect_hat <- function(y, subject, method, par, gh, include_subject_method) {
  beta0 <- par[[1L]]
  beta_method <- par[[2L]]
  sigma <- exp(par[[3L]])
  sigma_sm <- if (isTRUE(include_subject_method)) exp(par[[4L]]) else 0
  n_subjects <- max(subject)
  alpha_hat <- numeric(n_subjects)
  gamma1_hat <- numeric(n_subjects)
  gamma2_hat <- numeric(n_subjects)
  log_sqrt_pi <- 0.5 * log(pi)

  for (s in seq_len(n_subjects)) {
    idx <- which(subject == s)
    if (!isTRUE(include_subject_method)) {
      alpha <- sqrt(2) * sigma * gh$nodes
      log_terms <- log(gh$weights) - log_sqrt_pi
      for (q in seq_along(alpha)) {
        eta <- beta0 + ifelse(method[idx] == 2L, beta_method, 0) + alpha[[q]]
        log_terms[[q]] <- log_terms[[q]] + sum(y[idx] * eta - exp(eta) - lgamma(y[idx] + 1))
      }
      mx <- max(log_terms)
      w <- exp(log_terms - mx)
      alpha_hat[[s]] <- sum(w * alpha) / sum(w)
    } else {
      vals <- expand.grid(qa = seq_along(gh$nodes), q1 = seq_along(gh$nodes), q2 = seq_along(gh$nodes))
      alpha <- sqrt(2) * sigma * gh$nodes[vals$qa]
      gamma1 <- sqrt(2) * sigma_sm * gh$nodes[vals$q1]
      gamma2 <- sqrt(2) * sigma_sm * gh$nodes[vals$q2]
      log_terms <- log(gh$weights[vals$qa]) + log(gh$weights[vals$q1]) +
        log(gh$weights[vals$q2]) - 3 * log_sqrt_pi
      for (q in seq_along(alpha)) {
        eta <- beta0 + ifelse(method[idx] == 2L, beta_method + gamma2[[q]], gamma1[[q]]) + alpha[[q]]
        log_terms[[q]] <- log_terms[[q]] + sum(y[idx] * eta - exp(eta) - lgamma(y[idx] + 1))
      }
      mx <- max(log_terms)
      w <- exp(log_terms - mx)
      denom <- sum(w)
      alpha_hat[[s]] <- sum(w * alpha) / denom
      gamma1_hat[[s]] <- sum(w * gamma1) / denom
      gamma2_hat[[s]] <- sum(w * gamma2) / denom
    }
  }
  list(alpha = alpha_hat, gamma1 = gamma1_hat, gamma2 = gamma2_hat)
}

.mc_ccc_glmm_metrics_from_fit <- function(y,
                                          subject,
                                          method,
                                          replicate,
                                          fit,
                                          overdispersion) {
  beta0 <- fit$par[[1L]]
  beta_method <- fit$par[[2L]]
  sigma2_subject <- exp(2 * fit$par[[3L]])
  sigma2_subject_method <- if (isTRUE(fit$include_subject_method)) exp(2 * fit$par[[4L]]) else 0
  sigma2_method <- 0.5 * beta_method * beta_method
  sigma2_total <- sigma2_subject + sigma2_method + sigma2_subject_method
  eta1 <- beta0
  eta2 <- beta0 + beta_method
  eta_bar <- 0.5 * (eta1 + eta2)
  mu <- exp(eta_bar + 0.5 * sigma2_total)
  m_reps <- .mc_ccc_glmm_replicates_per_cell(subject, method, replicate)

  phi <- 1
  if (identical(overdispersion, "pearson")) {
    re_hat <- .mc_ccc_glmm_random_effect_hat(
      y = y,
      subject = subject,
      method = method,
      par = fit$par,
      gh = fit$gh,
      include_subject_method = fit$include_subject_method
    )
    eta_hat <- beta0 + ifelse(method == 2L, beta_method + re_hat$gamma2[subject], re_hat$gamma1[subject]) +
      re_hat$alpha[subject]
    mu_hat <- pmax(exp(eta_hat), .Machine$double.eps)
    pearson <- sum((y - mu_hat)^2 / mu_hat)
    n_re <- max(subject) + if (isTRUE(fit$include_subject_method)) 2L * max(subject) else 0L
    df <- max(length(y) - 2L - n_re, 1L)
    phi <- max(pearson / df, .Machine$double.eps)
  }

  cov_subject <- mu^2 * (exp(sigma2_subject) - 1)
  var_marginal_mean <- mu^2 * (exp(sigma2_total) - 1)
  residual_var <- phi * mu
  rho_ccc <- cov_subject / (var_marginal_mean + residual_var)
  rho_inter <- cov_subject / (var_marginal_mean + residual_var / m_reps)

  sigma2_intra <- sigma2_subject + sigma2_subject_method
  mu1 <- exp(eta1 + 0.5 * sigma2_intra)
  mu2 <- exp(eta2 + 0.5 * sigma2_intra)
  cov_intra1 <- mu1^2 * (exp(sigma2_intra) - 1)
  cov_intra2 <- mu2^2 * (exp(sigma2_intra) - 1)
  rho_intra1 <- cov_intra1 / (cov_intra1 + phi * mu1)
  rho_intra2 <- cov_intra2 / (cov_intra2 + phi * mu2)
  precision <- mu * (exp(sigma2_subject) - 1) /
    (mu * (exp(sigma2_subject + sigma2_subject_method) - 1) + phi)
  accuracy <- if (is.finite(precision) && precision > 0) rho_ccc / precision else NA_real_

  clamp01 <- function(x) pmax(0, pmin(1, x))
  list(
    rho_ccc = clamp01(rho_ccc),
    rho_ccc_inter = clamp01(rho_inter),
    rho_ccc_intra_method1 = clamp01(rho_intra1),
    rho_ccc_intra_method2 = clamp01(rho_intra2),
    sigma2_subject = sigma2_subject,
    sigma2_method = sigma2_method,
    sigma2_subject_method = sigma2_subject_method,
    phi = phi,
    mu = mu,
    precision = clamp01(precision),
    accuracy = clamp01(accuracy),
    n_obs = length(y),
    n_subjects = max(subject),
    m_reps = m_reps,
    beta0 = beta0,
    beta_method = beta_method,
    nll = fit$value,
    logLik = -fit$value,
    fit_engine = "matrixCorr_ghq",
    approximation = if (isTRUE(fit$include_subject_method)) "tensor_ghq" else "ghq",
    convergence_code = fit$convergence,
    optimiser_message = fit$message,
    converged = TRUE,
    iterations = as.integer(fit$counts[["function"]] %||% NA_integer_),
    fit_status = "converged"
  )
}

#' Poisson GLMM concordance correlation for count agreement
#'
#' @description
#' Computes pairwise generalized concordance correlation coefficients (CCC) for
#' count outcomes measured repeatedly by two or more methods, observers, or
#' devices. The current implementation fits Poisson-log generalized linear
#' mixed models and reports total, inter-method, and intra-method agreement
#' summaries from the fitted mean and variance components.
#'
#' @param data A data frame containing the measurements.
#' @param response Character. Name of the non-negative integer-like count column.
#' @param subject Character. Name of the subject identifier column.
#' @param method Character. Name of the method/observer column.
#' @param replicate Optional character. Name of the replicate column. If
#'   \code{NULL}, one replicate per subject-method cell is assumed.
#' @param family Character. Currently only \code{"poisson"} is implemented.
#' @param link Character. Currently only \code{"log"} is implemented.
#' @param overdispersion Character. \code{"none"} fixes \eqn{\phi = 1};
#'   \code{"pearson"} estimates \eqn{\phi} from post-fit Pearson residuals.
#' @param include_subject_method Logical. If \code{TRUE}, include an additional
#'   subject-by-method random intercept variance component.
#' @param ci Logical. If \code{TRUE}, compute delta-method confidence intervals
#'   for the total, inter-method, and intra-method CCC estimates.
#' @param conf_level Confidence level for delta-method confidence intervals.
#' @param max_iter Positive integer. Maximum optimiser iterations.
#' @param tol Positive numeric convergence tolerance.
#' @param n_threads Integer \eqn{\geq 1}. Number of OpenMP threads to make
#'   available to compiled code.
#' @param verbose Logical. If \code{TRUE}, emit progress messages.
#'
#' @details
#' The fitted model for a pair of methods is
#' \deqn{
#' Y_{ijl} \mid \alpha_i, \gamma_{ij}
#'   \sim \mathrm{Poisson}(\mu_{ijl}),
#' \qquad
#' \log(\mu_{ijl}) = \eta_j + \alpha_i + \gamma_{ij},
#' }
#' where \eqn{i} indexes subjects, \eqn{j} indexes methods, and \eqn{l}
#' indexes replicate readings. The subject effect \eqn{\alpha_i} captures
#' between-subject heterogeneity. When \code{include_subject_method = TRUE},
#' \eqn{\gamma_{ij}} captures subject-specific method departures; otherwise it
#' is fixed at zero. The fixed method difference contributes to disagreement
#' through \code{sigma2_method}.
#'
#' The model is fitted by marginal maximum likelihood using Gauss-Hermite
#' quadrature. The random-intercept model uses 40 quadrature points. The
#' subject-by-method model uses tensor-product quadrature with fewer points per
#' dimension because it integrates over a three-dimensional subject block. The
#' reported CCC quantities are then computed from the fitted Poisson-log mean
#' and variance components.
#'
#' The main matrix value, \code{rho_ccc}, is the total agreement coefficient.
#' Use it as the primary overall count-agreement summary when each individual
#' reading is the unit of inference. It penalizes lack of between-subject
#' signal, systematic method bias, subject-by-method disagreement, and Poisson
#' residual variation.
#'
#' \code{rho_ccc_inter} is the inter-method agreement coefficient for averages
#' over replicated readings. It is useful when decisions are based on the mean
#' of \eqn{m} replicate readings per subject-method cell. Replication reduces
#' only the residual count variation term; it does not dilute systematic method
#' disagreement or between-subject variation.
#'
#' \code{rho_ccc_intra_method1} and \code{rho_ccc_intra_method2} are
#' method-specific repeatability coefficients. Use them to diagnose whether one
#' method is internally more repeatable than the other on the count scale. These
#' are not direct method-comparison coefficients; they describe within-method
#' reproducibility after accounting for the fitted mean and random effects.
#'
#' \code{precision} isolates the share of non-systematic variation attributable
#' to subject ranking/heterogeneity, while \code{accuracy} is the ratio
#' \code{rho_ccc / precision}. Low accuracy with reasonable precision indicates
#' that disagreement is driven mainly by method bias or extra method-specific
#' variation rather than poor subject discrimination.
#'
#' \code{overdispersion = "none"} fixes \eqn{\phi = 1}, the Poisson model.
#' \code{overdispersion = "pearson"} replaces the residual term by a Pearson
#' dispersion estimate. Use the Pearson adjustment as a sensitivity analysis
#' when counts appear more variable than the Poisson model allows; it should
#' generally reduce CCC when extra-Poisson variation is present.
#'
#' ## Confidence intervals
#' When \code{ci = TRUE}, large-sample confidence intervals are computed for
#' \code{rho_ccc}, \code{rho_ccc_inter}, \code{rho_ccc_intra_method1}, and
#' \code{rho_ccc_intra_method2}. The implementation uses a delta-method
#' standard error based on the fitted GLMM parameter vector and the inverse
#' Hessian of the marginal negative log-likelihood:
#' \deqn{
#' \mathrm{Var}\{g(\hat\theta)\}
#'   \approx
#'   \nabla g(\hat\theta)^\top
#'   \mathrm{Var}(\hat\theta)
#'   \nabla g(\hat\theta),
#' }
#' where \eqn{g(\theta)} is the relevant CCC function. Gradients are evaluated
#' numerically by central finite differences. The reported point estimates are
#' not replaced by \eqn{g(\hat\theta)}; they remain the values from the standard
#' point-estimate path.
#'
#' Intervals are formed on a Fisher-Z transformed CCC scale,
#' \eqn{z = 0.5\log\{(1+\rho)/(1-\rho)\}}, and then back-transformed to the
#' CCC scale. This is the same broad strategy used elsewhere for CCC-style
#' intervals because it is usually more stable than a raw Wald interval near
#' the boundaries. CI limits are not forcibly truncated to \eqn{[0, 1]}.
#'
#' For \code{overdispersion = "pearson"}, \eqn{\phi} is a post-fit Pearson
#' dispersion estimate rather than a likelihood parameter. Delta-method
#' confidence intervals therefore treat \eqn{\phi} as fixed and a warning is
#' issued.
#'
#' This function currently implements the Poisson-log count-data case. The
#' \code{family} and \code{link} arguments are present for API stability and
#' future extensions, but only \code{family = "poisson"} and
#' \code{link = "log"} are currently supported. The returned estimates are
#' variance-component CCCs constrained to \eqn{[0, 1]}; they are not Lin's raw
#' moment CCC and should not be expected to produce negative values.
#'
#' @return
#' A symmetric numeric matrix of class \code{c("ccc_glmm", "ccc")} containing
#' \code{rho_ccc}. Additional pairwise matrices are stored as attributes:
#' \itemize{
#'   \item \code{rho_ccc_inter}: agreement for replicated method averages.
#'   \item \code{rho_ccc_intra_method1}, \code{rho_ccc_intra_method2}:
#'     method-specific repeatability coefficients.
#'   \item \code{sigma2_subject}, \code{sigma2_method},
#'     \code{sigma2_subject_method}: fitted variance/disagreement components.
#'   \item \code{phi}, \code{mu}, \code{precision}, \code{accuracy}: fitted
#'     count-scale diagnostics used in the CCC decomposition.
#'   \item \code{beta0}, \code{beta_method}, \code{nll}, \code{logLik},
#'     \code{convergence_code}: fitting diagnostics.
#'   \item \code{n_obs}, \code{n_subjects}, \code{m_reps}: design diagnostics.
#'   \item when \code{ci = TRUE}, standard error and confidence-limit matrices
#'     for \code{rho_ccc}, \code{rho_ccc_inter},
#'     \code{rho_ccc_intra_method1}, and \code{rho_ccc_intra_method2}.
#' }
#'
#' @references
#' Carrasco JL (2010). A generalized concordance correlation coefficient based
#' on the variance components generalized linear mixed models for overdispersed
#' count ddata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAbElEQVR4Xs2RQQrAMAgEfZgf7W9LAguybljJpR3wEse5JOL3ZObDb4x1loDhHbBOFU6i2Ddnw2KNiXcdAXygJlwE8OFVBHDgKrLgSInN4WMe9iXiqIVsTMjH7z/GhNTEibOxQswcYIWYOR/zAjBJfiXh3jZ6AAAAAElFTkSuQmCCata. \emph{Biometrics}.
#'
#' @examples
#' # Example 1: overall agreement for two Poisson count methods.
#' # Here the methods have similar means, so rho_ccc is mainly driven by
#' # between-subject count heterogeneity versus Poisson residual noise.
#' set.seed(1)
#' df1 <- expand.grid(
#'   subject = factor(seq_len(12)),
#'   method = factor(c("A", "B")),
#'   replicate = factor(seq_len(2))
#' )
#' subject_eff <- rnorm(12, 0, 0.5)
#' df1$eta <- 1.1 + subject_eff[as.integer(df1$subject)]
#' df1$y <- rpois(nrow(df1), exp(df1$eta))
#' fit1 <- ccc_glmm(df1, "y", "subject", "method", replicate = "replicate",
#' ci = TRUE)
#' fit1
#' summary(fit1)
#'
#' # Example 2: method bias lowers total agreement.
#' # This tests whether a systematic method shift is reflected in rho_ccc
#' # and the accuracy component.
#' set.seed(2)
#' df2 <- expand.grid(
#'   subject = factor(seq_len(12)),
#'   method = factor(c("A", "B")),
#'   replicate = factor(seq_len(2))
#' )
#' subject_eff <- rnorm(12, 0, 0.6)
#' df2$eta <- 1.0 + subject_eff[as.integer(df2$subject)] +
#'   ifelse(df2$method == "B", 0.6, 0)
#' df2$y <- rpois(nrow(df2), exp(df2$eta))
#' fit2 <- ccc_glmm(df2, "y", "subject", "method", replicate = "replicate")
#' summary(fit2)
#'
#' # Example 3: subject-by-method variation.
#' # This tests whether individual subjects respond differently by method.
#' # The subject-method component is useful when disagreement is not explained
#' # by a single fixed method bias.
#' set.seed(3)
#' df3 <- expand.grid(
#'   subject = factor(seq_len(10)),
#'   method = factor(c("A", "B")),
#'   replicate = factor(seq_len(2))
#' )
#' subject_eff <- rnorm(10, 0, 0.4)
#' subject_method_eff <- matrix(rnorm(20, 0, 0.25), nrow = 10)
#' method_id <- as.integer(df3$method)
#' df3$eta <- 1.1 + subject_eff[as.integer(df3$subject)] +
#'   subject_method_eff[cbind(as.integer(df3$subject), method_id)]
#' df3$y <- rpois(nrow(df3), exp(df3$eta))
#' fit3 <- ccc_glmm(
#'   df3, "y", "subject", "method",
#'   replicate = "replicate",
#'   include_subject_method = TRUE
#' )
#' attr(fit3, "sigma2_subject_method")
#'
#' # Example 4: four methods.
#' # This tests pairwise agreement across several count methods and helps
#' # identify which method pairs have the strongest total CCC.
#' set.seed(4)
#' df4 <- expand.grid(
#'   subject = factor(seq_len(10)),
#'   method = factor(c("A", "B", "C", "D")),
#'   replicate = factor(seq_len(2))
#' )
#' subject_eff <- rnorm(10, 0, 0.5)
#' method_bias <- c(A = 0, B = 0.1, C = 0.4, D = -0.2)
#' df4$eta <- 1.0 + subject_eff[as.integer(df4$subject)] +
#'   method_bias[as.character(df4$method)]
#' df4$y <- rpois(nrow(df4), exp(df4$eta))
#' fit4 <- ccc_glmm(df4, "y", "subject", "method", replicate = "replicate")
#' fit4
#' summary(fit4, n = 3)
#'
#' @export
ccc_glmm <- function(data,
                     response,
                     subject,
                     method,
                     replicate = NULL,
                     family = "poisson",
                     link = "log",
                     overdispersion = c("none", "pearson"),
                     include_subject_method = FALSE,
                     ci = FALSE,
                     conf_level = 0.95,
                     max_iter = 1000,
                     tol = 1e-8,
                     n_threads = getOption("matrixCorr.threads", 1L),
                     verbose = FALSE) {
  overdispersion <- match_arg(
    overdispersion,
    values = c("none", "pearson"),
    arg_name = "overdispersion"
  )

  if (!identical(family, "poisson")) {
    abort_bad_arg("family", message = "only {.val poisson} is implemented.")
  }
  if (!identical(link, "log")) {
    abort_bad_arg("link", message = "only {.val log} is implemented for Poisson GLMM CCC.")
  }

  check_bool(ci, arg = "ci")
  if (isTRUE(ci)) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
    if (identical(overdispersion, "pearson")) {
      cli::cli_warn(
        "Delta-method confidence intervals for {.arg overdispersion = \"pearson\"} treat the Pearson dispersion estimate as fixed."
      )
    }
  }
  if (!isTRUE(ci)) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  }
  max_iter <- check_scalar_int_pos(max_iter, arg = "max_iter")
  check_scalar_nonneg(tol, arg = "tol", strict = TRUE)
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")
  check_bool(verbose, arg = "verbose")
  check_bool(include_subject_method, arg = "include_subject_method")

  df <- as.data.frame(data)
  req_cols <- c(response, subject, method, replicate)
  req_cols <- req_cols[!vapply(req_cols, is.null, logical(1))]
  check_required_cols(df, req_cols, df_arg = "data")

  complete_mask <- stats::complete.cases(df[, req_cols, drop = FALSE])
  df <- df[complete_mask, , drop = FALSE]
  if (!nrow(df)) {
    abort_bad_arg("data", message = "must contain at least one complete row.")
  }

  if (!is.numeric(df[[response]]) || any(!is.finite(df[[response]]))) {
    abort_bad_arg("response", message = "must reference a finite numeric count column.")
  }
  if (any(df[[response]] < 0)) {
    abort_bad_arg("response", message = "must contain non-negative counts.")
  }
  if (any(abs(df[[response]] - round(df[[response]])) > sqrt(.Machine$double.eps))) {
    abort_bad_arg("response", message = "must contain integer-like counts.")
  }

  df[[subject]] <- droplevels(factor(df[[subject]]))
  df[[method]] <- droplevels(factor(df[[method]]))
  if (!is.null(replicate)) {
    df[[replicate]] <- droplevels(factor(df[[replicate]]))
  }

  method_levels <- levels(df[[method]])
  L <- length(method_levels)
  if (L < 2L) {
    abort_bad_arg("method", message = "must have at least two levels.")
  }

  prev_threads <- .mc_prepare_omp_threads(
    n_threads,
    n_threads_missing = missing(n_threads)
  )
  if (!is.null(prev_threads)) {
    on.exit(.mc_exit_omp_threads(prev_threads), add = TRUE)
  }
  inform_if_verbose(
    "Using {get_omp_threads()} OpenMP thread{?s}.",
    .verbose = verbose
  )

  dn <- list(method_levels, method_levels)
  est_mat <- matrix(NA_real_, L, L, dimnames = dn)
  rho_inter_mat <- matrix(NA_real_, L, L, dimnames = dn)
  rho_intra_1_mat <- matrix(NA_real_, L, L, dimnames = dn)
  rho_intra_2_mat <- matrix(NA_real_, L, L, dimnames = dn)
  sigma2_subject_mat <- matrix(NA_real_, L, L, dimnames = dn)
  sigma2_method_mat <- matrix(NA_real_, L, L, dimnames = dn)
  sigma2_subject_method_mat <- matrix(NA_real_, L, L, dimnames = dn)
  phi_mat <- matrix(NA_real_, L, L, dimnames = dn)
  mu_mat <- matrix(NA_real_, L, L, dimnames = dn)
  precision_mat <- matrix(NA_real_, L, L, dimnames = dn)
  accuracy_mat <- matrix(NA_real_, L, L, dimnames = dn)
  n_obs_mat <- matrix(NA_integer_, L, L, dimnames = dn)
  n_subjects_mat <- matrix(NA_integer_, L, L, dimnames = dn)
  m_reps_mat <- matrix(NA_integer_, L, L, dimnames = dn)
  converged_mat <- matrix(NA, L, L, dimnames = dn)
  iterations_mat <- matrix(NA_integer_, L, L, dimnames = dn)
  fit_status_mat <- matrix(NA_character_, L, L, dimnames = dn)
  beta0_mat <- matrix(NA_real_, L, L, dimnames = dn)
  beta_method_mat <- matrix(NA_real_, L, L, dimnames = dn)
  nll_mat <- matrix(NA_real_, L, L, dimnames = dn)
  loglik_mat <- matrix(NA_real_, L, L, dimnames = dn)
  convergence_code_mat <- matrix(NA_integer_, L, L, dimnames = dn)
  optimiser_message_mat <- matrix(NA_character_, L, L, dimnames = dn)
  rho_ccc_se_mat <- rho_ccc_lwr_mat <- rho_ccc_upr_mat <- NULL
  rho_inter_se_mat <- rho_inter_lwr_mat <- rho_inter_upr_mat <- NULL
  rho_intra_1_se_mat <- rho_intra_1_lwr_mat <- rho_intra_1_upr_mat <- NULL
  rho_intra_2_se_mat <- rho_intra_2_lwr_mat <- rho_intra_2_upr_mat <- NULL
  if (isTRUE(ci)) {
    rho_ccc_se_mat <- matrix(NA_real_, L, L, dimnames = dn)
    rho_ccc_lwr_mat <- matrix(NA_real_, L, L, dimnames = dn)
    rho_ccc_upr_mat <- matrix(NA_real_, L, L, dimnames = dn)
    rho_inter_se_mat <- matrix(NA_real_, L, L, dimnames = dn)
    rho_inter_lwr_mat <- matrix(NA_real_, L, L, dimnames = dn)
    rho_inter_upr_mat <- matrix(NA_real_, L, L, dimnames = dn)
    rho_intra_1_se_mat <- matrix(NA_real_, L, L, dimnames = dn)
    rho_intra_1_lwr_mat <- matrix(NA_real_, L, L, dimnames = dn)
    rho_intra_1_upr_mat <- matrix(NA_real_, L, L, dimnames = dn)
    rho_intra_2_se_mat <- matrix(NA_real_, L, L, dimnames = dn)
    rho_intra_2_lwr_mat <- matrix(NA_real_, L, L, dimnames = dn)
    rho_intra_2_upr_mat <- matrix(NA_real_, L, L, dimnames = dn)
  }

  for (i in seq_len(L - 1L)) {
    for (j in (i + 1L):L) {
      pair_levels <- method_levels[c(i, j)]
      df_sub <- df[df[[method]] %in% pair_levels, , drop = FALSE]
      df_sub[[method]] <- droplevels(factor(df_sub[[method]], levels = pair_levels))
      df_sub[[subject]] <- droplevels(factor(df_sub[[subject]]))
      if (!is.null(replicate)) {
        df_sub[[replicate]] <- droplevels(factor(df_sub[[replicate]]))
      }

      if (nlevels(df_sub[[subject]]) < 2L) {
        cli::cli_warn(
          "ccc_glmm Poisson fit skipped for pair ({.val {pair_levels[1]}}, {.val {pair_levels[2]}}): fewer than two subjects."
        )
        next
      }

      y <- as.numeric(df_sub[[response]])
      subject_int <- as.integer(df_sub[[subject]])
      method_int <- as.integer(df_sub[[method]])
      replicate_int <- if (is.null(replicate)) integer(0) else as.integer(df_sub[[replicate]])

      inform_if_verbose(
        "Fitting ccc_glmm pair {.val {pair_levels[1]}} vs {.val {pair_levels[2]}}.",
        .verbose = verbose
      )

      ans <- tryCatch(
        {
          fit <- .mc_ccc_glmm_fit(
            y = y,
            subject = subject_int,
            method = method_int,
            n_subjects = nlevels(df_sub[[subject]]),
            include_subject_method = include_subject_method,
            max_iter = max_iter,
            tol = tol,
            compute_vcov = ci
          )
          if (is.null(fit)) {
            NULL
          } else {
            metrics <- .mc_ccc_glmm_metrics_from_fit(
              y = y,
              subject = subject_int,
              method = method_int,
              replicate = replicate_int,
              fit = fit,
              overdispersion = overdispersion
            )
            list(fit = fit, metrics = metrics)
          }
        },
        error = function(e) {
          cli::cli_warn(
            "ccc_glmm Poisson fit failed for pair ({.val {pair_levels[1]}}, {.val {pair_levels[2]}}): {conditionMessage(e)}"
          )
          NULL
        }
      )

      if (!is.null(ans)) {
        fit_pair <- ans$fit
        ans <- ans$metrics
        boundary_subject <- is.finite(num_or_na(ans$sigma2_subject)) &&
          num_or_na(ans$sigma2_subject) <= 1e-7
        if (!isTRUE(ans$converged) || isTRUE(boundary_subject)) {
          reason <- if (isTRUE(boundary_subject)) "reached the subject-variance boundary" else "did not converge"
          cli::cli_warn(
            "ccc_glmm Poisson fit {reason} for pair ({.val {pair_levels[1]}}, {.val {pair_levels[2]}}); returning NA for this pair."
          )
          converged_mat[i, j] <- converged_mat[j, i] <- FALSE
          iterations_mat[i, j] <- iterations_mat[j, i] <- as.integer(ans$iterations %||% NA_integer_)
          convergence_code_mat[i, j] <- convergence_code_mat[j, i] <- as.integer(ans$convergence_code %||% NA_integer_)
          optimiser_message_mat[i, j] <- optimiser_message_mat[j, i] <- ans$optimiser_message %||% reason
          fit_status_mat[i, j] <- fit_status_mat[j, i] <- if (isTRUE(boundary_subject)) {
            "boundary_subject_variance"
          } else {
            "nonconverged"
          }
          next
        }
        est_mat[i, j] <- est_mat[j, i] <- num_or_na(ans$rho_ccc)
        rho_inter_mat[i, j] <- rho_inter_mat[j, i] <- num_or_na(ans$rho_ccc_inter)
        rho_intra_1_mat[i, j] <- num_or_na(ans$rho_ccc_intra_method1)
        rho_intra_2_mat[i, j] <- num_or_na(ans$rho_ccc_intra_method2)
        rho_intra_1_mat[j, i] <- num_or_na(ans$rho_ccc_intra_method2)
        rho_intra_2_mat[j, i] <- num_or_na(ans$rho_ccc_intra_method1)

        sigma2_subject_mat[i, j] <- sigma2_subject_mat[j, i] <- num_or_na(ans$sigma2_subject)
        sigma2_method_mat[i, j] <- sigma2_method_mat[j, i] <- num_or_na(ans$sigma2_method)
        sigma2_subject_method_mat[i, j] <- sigma2_subject_method_mat[j, i] <- num_or_na(ans$sigma2_subject_method)
        phi_mat[i, j] <- phi_mat[j, i] <- num_or_na(ans$phi)
        mu_mat[i, j] <- mu_mat[j, i] <- num_or_na(ans$mu)
        precision_mat[i, j] <- precision_mat[j, i] <- num_or_na(ans$precision)
        accuracy_mat[i, j] <- accuracy_mat[j, i] <- num_or_na(ans$accuracy)
        n_obs_mat[i, j] <- n_obs_mat[j, i] <- as.integer(ans$n_obs)
        n_subjects_mat[i, j] <- n_subjects_mat[j, i] <- as.integer(ans$n_subjects)
        m_reps_mat[i, j] <- m_reps_mat[j, i] <- as.integer(ans$m_reps)
        converged_mat[i, j] <- converged_mat[j, i] <- isTRUE(ans$converged)
        iterations_mat[i, j] <- iterations_mat[j, i] <- as.integer(ans$iterations)
        beta0_mat[i, j] <- beta0_mat[j, i] <- num_or_na(ans$beta0)
        beta_method_mat[i, j] <- num_or_na(ans$beta_method)
        beta_method_mat[j, i] <- -num_or_na(ans$beta_method)
        nll_mat[i, j] <- nll_mat[j, i] <- num_or_na(ans$nll)
        loglik_mat[i, j] <- loglik_mat[j, i] <- num_or_na(ans$logLik)
        convergence_code_mat[i, j] <- convergence_code_mat[j, i] <- as.integer(ans$convergence_code)
        optimiser_message_mat[i, j] <- optimiser_message_mat[j, i] <- ans$optimiser_message %||% ""
        fit_status_mat[i, j] <- fit_status_mat[j, i] <- ans$fit_status %||% "converged"

        if (isTRUE(ci)) {
          funcs <- .mc_ccc_glmm_delta_functions(
            include_subject_method = include_subject_method,
            phi = ans$phi,
            m_reps = ans$m_reps
          )
          ci_total <- .mc_delta_ci(
            estimate = ans$rho_ccc,
            f = funcs$rho_ccc,
            theta = fit_pair$par,
            vcov_theta = fit_pair$vcov_par,
            conf_level = conf_level,
            transform = "fisher"
          )
          ci_inter <- .mc_delta_ci(
            estimate = ans$rho_ccc_inter,
            f = funcs$rho_ccc_inter,
            theta = fit_pair$par,
            vcov_theta = fit_pair$vcov_par,
            conf_level = conf_level,
            transform = "fisher"
          )
          ci_intra1 <- .mc_delta_ci(
            estimate = ans$rho_ccc_intra_method1,
            f = funcs$rho_ccc_intra_method1,
            theta = fit_pair$par,
            vcov_theta = fit_pair$vcov_par,
            conf_level = conf_level,
            transform = "fisher"
          )
          ci_intra2 <- .mc_delta_ci(
            estimate = ans$rho_ccc_intra_method2,
            f = funcs$rho_ccc_intra_method2,
            theta = fit_pair$par,
            vcov_theta = fit_pair$vcov_par,
            conf_level = conf_level,
            transform = "fisher"
          )

          rho_ccc_se_mat[i, j] <- rho_ccc_se_mat[j, i] <- num_or_na(ci_total$se)
          rho_ccc_lwr_mat[i, j] <- rho_ccc_lwr_mat[j, i] <- num_or_na(ci_total$lwr.ci)
          rho_ccc_upr_mat[i, j] <- rho_ccc_upr_mat[j, i] <- num_or_na(ci_total$upr.ci)
          rho_inter_se_mat[i, j] <- rho_inter_se_mat[j, i] <- num_or_na(ci_inter$se)
          rho_inter_lwr_mat[i, j] <- rho_inter_lwr_mat[j, i] <- num_or_na(ci_inter$lwr.ci)
          rho_inter_upr_mat[i, j] <- rho_inter_upr_mat[j, i] <- num_or_na(ci_inter$upr.ci)

          rho_intra_1_se_mat[i, j] <- num_or_na(ci_intra1$se)
          rho_intra_1_lwr_mat[i, j] <- num_or_na(ci_intra1$lwr.ci)
          rho_intra_1_upr_mat[i, j] <- num_or_na(ci_intra1$upr.ci)
          rho_intra_2_se_mat[i, j] <- num_or_na(ci_intra2$se)
          rho_intra_2_lwr_mat[i, j] <- num_or_na(ci_intra2$lwr.ci)
          rho_intra_2_upr_mat[i, j] <- num_or_na(ci_intra2$upr.ci)

          rho_intra_1_se_mat[j, i] <- num_or_na(ci_intra2$se)
          rho_intra_1_lwr_mat[j, i] <- num_or_na(ci_intra2$lwr.ci)
          rho_intra_1_upr_mat[j, i] <- num_or_na(ci_intra2$upr.ci)
          rho_intra_2_se_mat[j, i] <- num_or_na(ci_intra1$se)
          rho_intra_2_lwr_mat[j, i] <- num_or_na(ci_intra1$lwr.ci)
          rho_intra_2_upr_mat[j, i] <- num_or_na(ci_intra1$upr.ci)
        }
      }
    }
  }

  diag(est_mat) <- 1
  diag(rho_inter_mat) <- 1
  diag(rho_intra_1_mat) <- 1
  diag(rho_intra_2_mat) <- 1
  diag(sigma2_subject_mat) <- 0
  diag(sigma2_method_mat) <- 0
  diag(sigma2_subject_method_mat) <- 0
  diag(precision_mat) <- 1
  diag(accuracy_mat) <- 1
  diag(phi_mat) <- if (identical(overdispersion, "none")) 1 else NA_real_

  out <- structure(
    est_mat,
    class = c("ccc_glmm", "ccc")
  )
  attr(out, "method") <- "GLMM variance-components CCC for count data"
  attr(out, "description") <- "Poisson GLMM concordance correlation for count data"
  attr(out, "package") <- "matrixCorr"
  attr(out, "family") <- family
  attr(out, "link") <- link
  attr(out, "overdispersion") <- overdispersion
  attr(out, "include_subject_method") <- include_subject_method
  attr(out, "rho_ccc") <- est_mat
  attr(out, "rho_ccc_inter") <- rho_inter_mat
  attr(out, "rho_ccc_intra_method1") <- rho_intra_1_mat
  attr(out, "rho_ccc_intra_method2") <- rho_intra_2_mat
  attr(out, "sigma2_subject") <- sigma2_subject_mat
  attr(out, "sigma2_method") <- sigma2_method_mat
  attr(out, "sigma2_subject_method") <- sigma2_subject_method_mat
  attr(out, "phi") <- phi_mat
  attr(out, "mu") <- mu_mat
  attr(out, "precision") <- precision_mat
  attr(out, "accuracy") <- accuracy_mat
  attr(out, "n_obs") <- n_obs_mat
  attr(out, "n_subjects") <- n_subjects_mat
  attr(out, "m_reps") <- m_reps_mat
  attr(out, "converged") <- converged_mat
  attr(out, "iterations") <- iterations_mat
  attr(out, "fit_status") <- fit_status_mat
  attr(out, "beta0") <- beta0_mat
  attr(out, "beta_method") <- beta_method_mat
  attr(out, "nll") <- nll_mat
  attr(out, "logLik") <- loglik_mat
  attr(out, "fit_engine") <- "matrixCorr_ghq"
  attr(out, "approximation") <- if (isTRUE(include_subject_method)) "tensor_ghq" else "ghq"
  attr(out, "convergence_code") <- convergence_code_mat
  attr(out, "optimiser_message") <- optimiser_message_mat
  if (isTRUE(ci)) {
    attr(out, "rho_ccc_se") <- rho_ccc_se_mat
    attr(out, "rho_ccc_lwr.ci") <- rho_ccc_lwr_mat
    attr(out, "rho_ccc_upr.ci") <- rho_ccc_upr_mat
    attr(out, "rho_ccc_inter_se") <- rho_inter_se_mat
    attr(out, "rho_ccc_inter_lwr.ci") <- rho_inter_lwr_mat
    attr(out, "rho_ccc_inter_upr.ci") <- rho_inter_upr_mat
    attr(out, "rho_ccc_intra_method1_se") <- rho_intra_1_se_mat
    attr(out, "rho_ccc_intra_method1_lwr.ci") <- rho_intra_1_lwr_mat
    attr(out, "rho_ccc_intra_method1_upr.ci") <- rho_intra_1_upr_mat
    attr(out, "rho_ccc_intra_method2_se") <- rho_intra_2_se_mat
    attr(out, "rho_ccc_intra_method2_lwr.ci") <- rho_intra_2_lwr_mat
    attr(out, "rho_ccc_intra_method2_upr.ci") <- rho_intra_2_upr_mat
    attr(out, "ci") <- TRUE
    attr(out, "conf_level") <- conf_level
    attr(out, "conf.level") <- conf_level
    attr(out, "ci_method") <- "delta_fisher"
  } else {
    attr(out, "ci") <- FALSE
  }
  out
}

#' @export
#' @method print ccc_glmm
print.ccc_glmm <- function(x,
                           digits = 4,
                           ci_digits = 4,
                           n = NULL,
                           topn = NULL,
                           max_vars = NULL,
                           width = NULL,
                           show_ci = NULL,
                           ...) {
  est <- if (is.list(x) && !is.null(x$est)) as.matrix(x$est) else as.matrix(x)
  .mc_print_corr_matrix(
    x,
    header = "GLMM concordance correlation matrix",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    mat = est,
    ...
  )
  invisible(x)
}

#' @export
#' @method summary ccc_glmm
summary.ccc_glmm <- function(object,
                             digits = 4,
                             ci_digits = 2,
                             n = NULL,
                             topn = NULL,
                             max_vars = NULL,
                             width = NULL,
                             show_ci = NULL,
                             ...) {
  base <- summary.ccc(
    object,
    digits = digits,
    ci_digits = ci_digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
  if (all(c("item1", "item2") %in% names(base))) {
    names(base)[match(c("item1", "item2"), names(base))] <- c("method1", "method2")
  }
  est <- as.matrix(object)
  extract_pair <- function(mat) {
    out <- numeric(nrow(est) * (ncol(est) - 1L) / 2L)
    k <- 0L
    for (i in seq_len(nrow(est) - 1L)) {
      for (j in (i + 1L):ncol(est)) {
        k <- k + 1L
        out[k] <- suppressWarnings(as.numeric(mat[i, j]))
      }
    }
    out
  }

  attrs <- c(
    "rho_ccc_inter",
    "rho_ccc_intra_method1",
    "rho_ccc_intra_method2",
    "rho_ccc_inter_se",
    "rho_ccc_inter_lwr.ci",
    "rho_ccc_inter_upr.ci",
    "rho_ccc_intra_method1_se",
    "rho_ccc_intra_method1_lwr.ci",
    "rho_ccc_intra_method1_upr.ci",
    "rho_ccc_intra_method2_se",
    "rho_ccc_intra_method2_lwr.ci",
    "rho_ccc_intra_method2_upr.ci",
    "sigma2_subject",
    "sigma2_method",
    "sigma2_subject_method",
    "phi",
    "mu",
    "precision",
    "accuracy",
    "beta0",
    "beta_method",
    "nll",
    "logLik",
    "convergence_code",
    "n_subjects",
    "n_obs",
    "m_reps"
  )
  for (nm in attrs) {
    val <- attr(object, nm, exact = TRUE)
    if (is.matrix(val)) {
      base[[nm]] <- extract_pair(val)
    }
  }
  if (isTRUE(attr(object, "ci", exact = TRUE))) {
    ci_map <- c(
      se = "rho_ccc_se",
      lwr.ci = "rho_ccc_lwr.ci",
      upr.ci = "rho_ccc_upr.ci"
    )
    for (nm in names(ci_map)) {
      val <- attr(object, ci_map[[nm]], exact = TRUE)
      if (is.matrix(val)) {
        base[[nm]] <- extract_pair(val)
      }
    }
  }

  num_cols <- names(base)[vapply(base, is.numeric, logical(1))]
  for (nm in num_cols) {
    if (!nm %in% c("n_subjects", "n_obs", "m_reps")) {
      base[[nm]] <- round(base[[nm]], digits)
    }
  }
  if ("n_subjects" %in% names(base)) base$n_subjects <- as.integer(base$n_subjects)
  if ("n_obs" %in% names(base)) base$n_obs <- as.integer(base$n_obs)
  if ("m_reps" %in% names(base)) base$m_reps <- as.integer(base$m_reps)

  class(base) <- c("summary.ccc_glmm", class(base))
  attr(base, "method") <- attr(object, "method", exact = TRUE)
  attr(base, "family") <- attr(object, "family", exact = TRUE)
  attr(base, "link") <- attr(object, "link", exact = TRUE)
  attr(base, "overdispersion") <- attr(object, "overdispersion", exact = TRUE)
  attr(base, "has_ci") <- isTRUE(attr(object, "ci", exact = TRUE))
  attr(base, "conf.level") <- attr(object, "conf.level", exact = TRUE)
  attr(base, "ci_method") <- attr(object, "ci_method", exact = TRUE)
  attr(base, "digits") <- digits
  attr(base, "ci_digits") <- ci_digits
  attr(base, "print_n") <- n
  attr(base, "print_topn") <- topn
  base
}

#' @export
#' @method print summary.ccc_glmm
print.summary.ccc_glmm <- function(x,
                                   digits = NULL,
                                   n = NULL,
                                   topn = NULL,
                                   max_vars = NULL,
                                   width = NULL,
                                   show_ci = NULL,
                                   style = c("auto", "compact", "pairwise"),
                                   ...) {
  style <- match.arg(style)
  digits <- .mc_coalesce(digits, attr(x, "digits", exact = TRUE))
  if (is.null(digits) || !is.finite(digits)) digits <- 4L
  n <- .mc_coalesce(n, attr(x, "print_n", exact = TRUE))
  topn <- .mc_coalesce(topn, attr(x, "print_topn", exact = TRUE))

  if (identical(style, "auto")) {
    style <- if (nrow(x) <= 1L) "pairwise" else "compact"
  }

  if (identical(style, "compact")) {
    rows <- seq_len(nrow(x))
    if (!is.null(n)) {
      n <- check_scalar_int_pos(n, arg = "n")
      rows <- rows[seq_len(min(length(rows), n))]
    }
    omitted <- nrow(x) - length(rows)

    print_section <- function(title, cols) {
      cols <- cols[cols %in% names(x)]
      if (!length(cols)) return(invisible(NULL))
      cli::cat_line("")
      cli::cat_line(title)
      cli::cat_line("")
      print.data.frame(x[rows, cols, drop = FALSE], row.names = FALSE, right = FALSE, ...)
      invisible(NULL)
    }

    cli::cat_line("GLMM concordance correlation summary")
    if (omitted > 0L) {
      cli::cat_line("Showing ", length(rows), " of ", nrow(x), " method pairs.")
    }

    print_section(
      "Concordance estimates",
      c("method1", "method2", "estimate", "se", "lwr.ci", "upr.ci", "precision", "accuracy")
    )
    print_section(
      "Inter-method agreement for replicate means",
      c(
        "method1", "method2",
        "rho_ccc_inter",
        "rho_ccc_inter_se",
        "rho_ccc_inter_lwr.ci",
        "rho_ccc_inter_upr.ci"
      )
    )
    print_section(
      "Intra-method repeatability",
      c(
        "method1", "method2",
        "rho_ccc_intra_method1",
        "rho_ccc_intra_method1_se",
        "rho_ccc_intra_method1_lwr.ci",
        "rho_ccc_intra_method1_upr.ci",
        "rho_ccc_intra_method2",
        "rho_ccc_intra_method2_se",
        "rho_ccc_intra_method2_lwr.ci",
        "rho_ccc_intra_method2_upr.ci"
      )
    )
    print_section(
      "Variance components and mean",
      c("method1", "method2", "sigma2_subject", "sigma2_method", "sigma2_subject_method", "mu", "phi")
    )
    print_section(
      "Fit diagnostics",
      c("method1", "method2", "beta0", "beta_method", "nll", "convergence_code")
    )
    print_section(
      "Design information",
      c("method1", "method2", "n_subjects", "n_obs", "m_reps")
    )
    if (omitted > 0L) {
      cli::cat_line(
        "",
        "... ", omitted, " more method pair", if (omitted == 1L) "" else "s",
        " not shown. Use a larger `n` or `as.data.frame()` to inspect all rows."
      )
    }
    return(invisible(x))
  }

  rows <- seq_len(nrow(x))
  if (!is.null(n)) {
    n <- check_scalar_int_pos(n, arg = "n")
    rows <- rows[seq_len(min(length(rows), n))]
  }

  value <- function(row, name) {
    if (!name %in% names(row)) return(NA)
    out <- row[[name]]
    if (is.numeric(out) && name != "convergence_code") {
      out <- round(out, digits)
    }
    out
  }

  make_table <- function(row, labels, names) {
    data.frame(
      quantity = labels,
      value = vapply(names, function(nm) as.character(value(row, nm)), character(1)),
      check.names = FALSE
    )
  }

  make_concordance_table <- function(row, method1, method2) {
    labels <- c(
      "Total CCC (single count)",
      "Inter-CCC (replicate mean)",
      paste0("Intra CCC (", method1, ")"),
      paste0("Intra CCC (", method2, ")"),
      "Precision",
      "Accuracy"
    )
    estimate_names <- c(
      "estimate",
      "rho_ccc_inter",
      "rho_ccc_intra_method1",
      "rho_ccc_intra_method2",
      "precision",
      "accuracy"
    )
    out <- data.frame(
      quantity = labels,
      estimate = vapply(estimate_names, function(nm) as.character(value(row, nm)), character(1)),
      check.names = FALSE
    )

    if (isTRUE(attr(x, "has_ci", exact = TRUE))) {
      se_names <- c(
        "se",
        "rho_ccc_inter_se",
        "rho_ccc_intra_method1_se",
        "rho_ccc_intra_method2_se",
        NA_character_,
        NA_character_
      )
      lwr_names <- c(
        "lwr.ci",
        "rho_ccc_inter_lwr.ci",
        "rho_ccc_intra_method1_lwr.ci",
        "rho_ccc_intra_method2_lwr.ci",
        NA_character_,
        NA_character_
      )
      upr_names <- c(
        "upr.ci",
        "rho_ccc_inter_upr.ci",
        "rho_ccc_intra_method1_upr.ci",
        "rho_ccc_intra_method2_upr.ci",
        NA_character_,
        NA_character_
      )
      pull <- function(nms) {
        vapply(nms, function(nm) {
          if (is.na(nm)) "" else as.character(value(row, nm))
        }, character(1))
      }
      out$se <- pull(se_names)
      out$lwr.ci <- pull(lwr_names)
      out$upr.ci <- pull(upr_names)
    }

    out
  }

  cli::cat_line("GLMM concordance correlation summary")
  if (length(rows) < nrow(x)) {
    cli::cat_line("Showing ", length(rows), " of ", nrow(x), " method pairs.")
  }

  for (idx in rows) {
    row <- x[idx, , drop = FALSE]
    method1 <- as.character(row$method1)
    method2 <- as.character(row$method2)

    cli::cat_line("")
    cli::cat_line("Pair: ", method1, " vs ", method2)
    cli::cat_line("")

    concordance <- make_concordance_table(row, method1, method2)
    cli::cat_line("Concordance")
    print.data.frame(concordance, row.names = FALSE, right = FALSE, ...)

    variance <- make_table(
      row,
      labels = c(
        "Subject variance",
        "Method disagreement variance",
        "Subject-method variance",
        "Marginal mean",
        "Dispersion phi"
      ),
      names = c(
        "sigma2_subject",
        "sigma2_method",
        "sigma2_subject_method",
        "mu",
        "phi"
      )
    )
    cli::cat_line("")
    cli::cat_line("Variance components and mean")
    print.data.frame(variance, row.names = FALSE, right = FALSE, ...)

    fit <- make_table(
      row,
      labels = c(
        paste0("Intercept (", method1, ")"),
        paste0("Method effect (", method2, " - ", method1, ")"),
        "Negative log-likelihood",
        "Log-likelihood",
        "Convergence code",
        "Subjects",
        "Observations",
        "Replicates per subject-method"
      ),
      names = c(
        "beta0",
        "beta_method",
        "nll",
        "logLik",
        "convergence_code",
        "n_subjects",
        "n_obs",
        "m_reps"
      )
    )
    cli::cat_line("")
    cli::cat_line("Fit and design")
    print.data.frame(fit, row.names = FALSE, right = FALSE, ...)
  }

  invisible(x)
}
