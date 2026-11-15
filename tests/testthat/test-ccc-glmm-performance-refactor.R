make_ccc_glmm_perf_data <- function(n_subjects = 20,
                                    reps = 2,
                                    method_bias = 0.2,
                                    subject_sd = 0.7,
                                    seed = 1) {
  set.seed(seed)
  df <- expand.grid(
    subject = factor(seq_len(n_subjects)),
    method = factor(c("A", "B")),
    replicate = factor(seq_len(reps))
  )
  subj_eff <- stats::rnorm(n_subjects, 0, subject_sd)
  df$eta <- 1.5 + subj_eff[as.integer(df$subject)] +
    ifelse(df$method == "B", method_bias, 0)
  df$y <- stats::rpois(nrow(df), exp(df$eta))
  df
}

reference_subject_blocks <- function(y, subject, method) {
  ids <- sort(unique(subject))
  out <- lapply(ids, function(id) {
    keep <- subject == id
    y_id <- y[keep]
    m_id <- method[keep]
    list(
      y1 = sum(y_id[m_id == 1L]),
      y2 = sum(y_id[m_id == 2L]),
      n1 = sum(m_id == 1L),
      n2 = sum(m_id == 2L),
      log_factorial = sum(lgamma(y_id + 1))
    )
  })
  list(
    y1 = vapply(out, `[[`, numeric(1), "y1"),
    y2 = vapply(out, `[[`, numeric(1), "y2"),
    n1 = vapply(out, `[[`, numeric(1), "n1"),
    n2 = vapply(out, `[[`, numeric(1), "n2"),
    log_factorial = vapply(out, `[[`, numeric(1), "log_factorial")
  )
}

reference_poisson_ghq_nll <- function(theta, y, subject, method, include_subject_method, nodes, weights) {
  beta0 <- theta[[1L]]
  beta_method <- theta[[2L]]
  sigma_subject <- exp(theta[[3L]])
  sigma_sm <- if (isTRUE(include_subject_method)) exp(theta[[4L]]) else 0
  log_w <- log(weights)
  log_sqrt_pi <- 0.5 * log(pi)
  total <- 0
  for (id in sort(unique(subject))) {
    keep <- subject == id
    y_id <- y[keep]
    m_id <- method[keep]
    if (!include_subject_method) {
      log_terms <- vapply(seq_along(nodes), function(q) {
        alpha <- sqrt(2) * sigma_subject * nodes[[q]]
        eta <- ifelse(m_id == 1L, beta0 + alpha, beta0 + beta_method + alpha)
        log_w[[q]] - log_sqrt_pi + sum(stats::dpois(y_id, exp(eta), log = TRUE))
      }, numeric(1))
    } else {
      log_terms <- numeric(length(nodes)^3)
      pos <- 0L
      for (qa in seq_along(nodes)) {
        alpha <- sqrt(2) * sigma_subject * nodes[[qa]]
        for (q1 in seq_along(nodes)) {
          gamma1 <- sqrt(2) * sigma_sm * nodes[[q1]]
          for (q2 in seq_along(nodes)) {
            gamma2 <- sqrt(2) * sigma_sm * nodes[[q2]]
            eta <- ifelse(
              m_id == 1L,
              beta0 + alpha + gamma1,
              beta0 + beta_method + alpha + gamma2
            )
            pos <- pos + 1L
            log_terms[[pos]] <- log_w[[qa]] + log_w[[q1]] + log_w[[q2]] - 3 * log_sqrt_pi +
              sum(stats::dpois(y_id, exp(eta), log = TRUE))
          }
        }
      }
    }
    m <- max(log_terms)
    total <- total + m + log(sum(exp(log_terms - m)))
  }
  -total
}

test_that("block Poisson GLMM likelihood matches an independent Gauss-Hermite reference", {
  df <- make_ccc_glmm_perf_data(n_subjects = 8, reps = 2, method_bias = 0.4, seed = 101)
  y <- as.numeric(df$y)
  subject <- as.integer(droplevels(factor(df$subject)))
  method <- as.integer(droplevels(factor(df$method)))
  n_subjects <- max(subject)

  blocks <- matrixCorr:::ccc_glmm_poisson_prepare_blocks_cpp(y, subject, method, n_subjects)
  ref_blocks <- reference_subject_blocks(y, subject, method)
  expect_equal(blocks$y1, ref_blocks$y1, tolerance = 1e-12)
  expect_equal(blocks$y2, ref_blocks$y2, tolerance = 1e-12)
  expect_equal(blocks$n1, ref_blocks$n1, tolerance = 1e-12)
  expect_equal(blocks$n2, ref_blocks$n2, tolerance = 1e-12)
  expect_equal(blocks$log_factorial, ref_blocks$log_factorial, tolerance = 1e-12)

  gh_simple <- matrixCorr:::.mc_ccc_glmm_ghq_rule(11L)
  gh_sm <- matrixCorr:::.mc_ccc_glmm_ghq_rule(7L)
  theta_simple <- list(
    c(1.2, 0.15, log(0.35)),
    c(1.4, -0.2, log(0.65)),
    c(0.8, 0.6, log(1.1))
  )
  theta_sm <- list(
    c(1.2, 0.15, log(0.35), log(0.25)),
    c(1.4, -0.2, log(0.65), log(0.40)),
    c(0.8, 0.6, log(1.1), log(0.55))
  )
  check_case <- function(theta, include_sm, gh) {
    prepared <- matrixCorr:::ccc_glmm_poisson_ghq_nll_blocks_cpp(
      theta,
      blocks$y1, blocks$y2, blocks$n1, blocks$n2, blocks$log_factorial,
      include_sm, gh$nodes, gh$weights
    )
    raw <- reference_poisson_ghq_nll(theta, y, subject, method, include_sm, gh$nodes, gh$weights)
    expect_equal(prepared, raw, tolerance = 1e-12)
  }
  for (theta in theta_simple) check_case(theta, FALSE, gh_simple)
  for (theta in theta_sm) check_case(theta, TRUE, gh_sm)
})

analytic_ccc_glmm_grad <- function(theta, include_subject_method, phi, m_reps) {
  nms <- c("t1", "t2", "t3", if (include_subject_method) "t4")
  a <- "exp(2 * t3)"
  b <- if (include_subject_method) "exp(2 * t4)" else "0"
  cc <- "(0.5 * t2 * t2)"
  total <- sprintf("(%s + %s + %s)", a, b, cc)
  intra <- sprintf("(%s + %s)", a, b)
  mu <- sprintf("exp(t1 + 0.5 * t2 + 0.5 * %s)", total)
  numer <- sprintf("(%s * (exp(%s) - 1))", mu, a)
  mu1 <- sprintf("exp(t1 + 0.5 * %s)", intra)
  mu2 <- sprintf("exp(t1 + t2 + 0.5 * %s)", intra)
  ni <- function(mu_i) sprintf("(%s * (exp(%s) - 1))", mu_i, intra)
  specs <- c(
    rho_ccc = sprintf("%s / (%s * (exp(%s) - 1) + phi)", numer, mu, total),
    rho_ccc_inter = sprintf("%s / (%s * (exp(%s) - 1) + phi / m_reps)", numer, mu, total),
    rho_ccc_intra_method1 = sprintf("%s / (%s + phi)", ni(mu1), ni(mu1)),
    rho_ccc_intra_method2 = sprintf("%s / (%s + phi)", ni(mu2), ni(mu2))
  )
  env <- list2env(list(phi = phi, m_reps = m_reps), parent = baseenv())
  for (nm in nms) env[[nm]] <- theta[[match(nm, nms)]]
  grad <- matrix(NA_real_, 4L, length(theta), dimnames = list(names(specs), NULL))
  metrics <- setNames(rep(NA_real_, length(specs)), names(specs))
  for (nm in names(specs)) {
    val <- eval(stats::deriv(parse(text = specs[[nm]])[[1L]], nms), env)
    metrics[[nm]] <- as.numeric(val)
    grad[nm, ] <- attr(val, "gradient")
  }
  list(metrics = metrics, gradient = grad)
}

fisher_delta_from_grad <- function(estimate, grad_row, vcov_theta, conf_level) {
  var_est <- as.numeric(t(grad_row) %*% vcov_theta %*% grad_row)
  se <- sqrt(var_est)
  zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)
  estimate_safe <- max(min(estimate, 1 - 1e-12), -1 + 1e-12)
  z_est <- atanh(estimate_safe)
  se_z <- se / (1 - estimate_safe^2)
  c(se = se, lwr.ci = tanh(z_est - zcrit * se_z), upr.ci = tanh(z_est + zcrit * se_z))
}

test_that("vectorized ccc_glmm delta CIs match an analytic gradient reference", {
  cases <- list(
    list(
      include_subject_method = FALSE,
      theta = c(1.1, 0.25, log(0.55)),
      vcov = matrix(c(
        0.02, 0.001, 0.0005,
        0.001, 0.03, 0.002,
        0.0005, 0.002, 0.015
      ), 3L, 3L)
    ),
    list(
      include_subject_method = TRUE,
      theta = c(1.1, 0.25, log(0.55), log(0.30)),
      vcov = matrix(c(
        0.02, 0.001, 0.0004, 0.0002,
        0.001, 0.03, 0.001, 0.0003,
        0.0004, 0.001, 0.015, 0.0005,
        0.0002, 0.0003, 0.0005, 0.012
      ), 4L, 4L)
    )
  )
  phi <- 1.2
  m_reps <- 2L
  for (case in cases) {
    theta <- case$theta
    include_sm <- case$include_subject_method
    estimates <- matrixCorr:::.mc_ccc_glmm_metric_vector(theta, include_sm, phi, m_reps)
    expect_true(all(estimates > 0.05 & estimates < 0.95))
    live_grad <- matrixCorr:::.mc_num_grad_vec(
      function(x) matrixCorr:::.mc_ccc_glmm_metric_vector(x, include_sm, phi, m_reps),
      theta
    )
    analytic <- analytic_ccc_glmm_grad(theta, include_sm, phi, m_reps)
    expect_equal(unname(estimates), unname(analytic$metrics), tolerance = 1e-8)
    # Central steps are 1e-5, so the leading error is O(step^2) = 1e-10 times
    # a moderate third derivative. This 1e-8 check compares the numerical
    # gradient with an analytical reference. It is not a tighter form of the
    # old scalar-versus-vector identity check.
    expect_equal(live_grad, analytic$gradient, tolerance = 1e-8)
    vec <- matrixCorr:::.mc_delta_ci_vec(
      estimates,
      function(x) matrixCorr:::.mc_ccc_glmm_metric_vector(x, include_sm, phi, m_reps),
      theta,
      case$vcov,
      conf_level = 0.9,
      transform = "fisher"
    )
    for (nm in names(estimates)) {
      ref <- fisher_delta_from_grad(analytic$metrics[[nm]], analytic$gradient[nm, ], case$vcov, 0.9)
      expect_equal(vec[nm, "se"], unname(ref[["se"]]), tolerance = 1e-8, info = nm)
      expect_equal(vec[nm, "lwr.ci"], unname(ref[["lwr.ci"]]), tolerance = 1e-8, info = nm)
      expect_equal(vec[nm, "upr.ci"], unname(ref[["upr.ci"]]), tolerance = 1e-8, info = nm)
    }
  }
})

test_that("ccc_glmm optimized path preserves CI object structure", {
  df <- make_ccc_glmm_perf_data(n_subjects = 16, reps = 2, method_bias = 0.5, seed = 102)
  expect_warning(
    fit <- ccc_glmm(
      df,
      "y", "subject", "method",
      replicate = "replicate",
      overdispersion = "pearson",
      ci = TRUE,
      max_iter = 500
    ),
    "treat the Pearson dispersion estimate as fixed"
  )

  expect_s3_class(fit, "ccc_glmm")
  expect_identical(class(fit), c("ccc_glmm", "ccc"))
  expect_identical(attr(fit, "ci", exact = TRUE), TRUE)
  expect_identical(attr(fit, "ci_method", exact = TRUE), "delta_fisher")
  for (nm in c(
    "rho_ccc_se", "rho_ccc_lwr.ci", "rho_ccc_upr.ci",
    "rho_ccc_inter_se", "rho_ccc_intra_method1_se",
    "rho_ccc_intra_method2_se"
  )) {
    val <- attr(fit, nm, exact = TRUE)
    expect_identical(dim(val), dim(fit), info = nm)
    expect_identical(dimnames(val), dimnames(fit), info = nm)
  }

  expect_glmm_symbolic_assembly <- function(include_subject_method) {
    assembled <- ccc_glmm(
      df,
      "y", "subject", "method",
      replicate = "replicate",
      overdispersion = "pearson",
      include_subject_method = include_subject_method,
      ci = TRUE,
      conf_level = 0.95,
      max_iter = 500
    )
    y <- as.numeric(df$y)
    subject <- as.integer(droplevels(factor(df$subject)))
    method <- as.integer(droplevels(factor(df$method)))
    fitted <- matrixCorr:::.mc_ccc_glmm_fit(
      y = y,
      subject = subject,
      method = method,
      n_subjects = max(subject),
      include_subject_method = include_subject_method,
      max_iter = 500L,
      tol = 1e-8,
      compute_vcov = TRUE
    )
    phi <- attr(assembled, "phi", exact = TRUE)["A", "B"]
    m_reps <- attr(assembled, "m_reps", exact = TRUE)["A", "B"]
    analytic <- analytic_ccc_glmm_grad(
      fitted$par,
      include_subject_method,
      phi,
      m_reps
    )
    stored <- c(
      rho_ccc = unname(assembled["A", "B"]),
      rho_ccc_inter = unname(attr(assembled, "rho_ccc_inter", exact = TRUE)["A", "B"]),
      rho_ccc_intra_method1 = unname(attr(assembled, "rho_ccc_intra_method1", exact = TRUE)["A", "B"]),
      rho_ccc_intra_method2 = unname(attr(assembled, "rho_ccc_intra_method2", exact = TRUE)["A", "B"])
    )
    expect_equal(stored, analytic$metrics, tolerance = 1e-8, info = include_subject_method)
    limits <- list(
      rho_ccc = c("rho_ccc_se", "rho_ccc_lwr.ci", "rho_ccc_upr.ci"),
      rho_ccc_inter = c("rho_ccc_inter_se", "rho_ccc_inter_lwr.ci", "rho_ccc_inter_upr.ci"),
      rho_ccc_intra_method1 = c("rho_ccc_intra_method1_se", "rho_ccc_intra_method1_lwr.ci", "rho_ccc_intra_method1_upr.ci"),
      rho_ccc_intra_method2 = c("rho_ccc_intra_method2_se", "rho_ccc_intra_method2_lwr.ci", "rho_ccc_intra_method2_upr.ci")
    )
    for (nm in names(limits)) {
      ref <- fisher_delta_from_grad(
        analytic$metrics[[nm]],
        analytic$gradient[nm, ],
        fitted$vcov_par,
        0.95
      )
      got <- unname(vapply(limits[[nm]], function(attr_nm) {
        unname(attr(assembled, attr_nm, exact = TRUE)["A", "B"])
      }, numeric(1)))
      expect_equal(got[[1L]], unname(ref[["se"]]), tolerance = 1e-8, info = nm)
      expect_equal(got[[2L]], unname(ref[["lwr.ci"]]), tolerance = 1e-8, info = nm)
      expect_equal(got[[3L]], unname(ref[["upr.ci"]]), tolerance = 1e-8, info = nm)
    }
  }
  expect_warning(
    expect_glmm_symbolic_assembly(FALSE),
    "treat the Pearson dispersion estimate as fixed"
  )
  expect_warning(
    expect_glmm_symbolic_assembly(TRUE),
    "treat the Pearson dispersion estimate as fixed"
  )
})
