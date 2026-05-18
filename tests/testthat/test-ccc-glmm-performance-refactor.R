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

test_that("block-based Poisson GLMM likelihood matches raw likelihood", {
  df <- make_ccc_glmm_perf_data(n_subjects = 8, reps = 2, method_bias = 0.4, seed = 101)
  y <- as.numeric(df$y)
  subject <- as.integer(droplevels(factor(df$subject)))
  method <- as.integer(droplevels(factor(df$method)))
  n_subjects <- max(subject)

  blocks <- matrixCorr:::ccc_glmm_poisson_prepare_blocks_cpp(y, subject, method, n_subjects)
  gh_simple <- matrixCorr:::.mc_ccc_glmm_ghq_rule(11L)
  gh_sm <- matrixCorr:::.mc_ccc_glmm_ghq_rule(7L)

  theta_simple <- list(
    c(1.2, 0.15, log(0.35)),
    c(1.4, -0.2, log(0.65)),
    c(0.8, 0.6, log(1.1))
  )
  for (theta in theta_simple) {
    raw <- matrixCorr:::ccc_glmm_poisson_ghq_nll_cpp(
      theta, y, subject, method, n_subjects, FALSE,
      gh_simple$nodes, gh_simple$weights
    )
    prepared <- matrixCorr:::ccc_glmm_poisson_ghq_nll_blocks_cpp(
      theta,
      blocks$y1, blocks$y2, blocks$n1, blocks$n2, blocks$log_factorial,
      FALSE, gh_simple$nodes, gh_simple$weights
    )
    expect_equal(prepared, raw, tolerance = 1e-12)
  }

  theta_sm <- list(
    c(1.2, 0.15, log(0.35), log(0.25)),
    c(1.4, -0.2, log(0.65), log(0.40)),
    c(0.8, 0.6, log(1.1), log(0.55))
  )
  for (theta in theta_sm) {
    raw <- matrixCorr:::ccc_glmm_poisson_ghq_nll_cpp(
      theta, y, subject, method, n_subjects, TRUE,
      gh_sm$nodes, gh_sm$weights
    )
    prepared <- matrixCorr:::ccc_glmm_poisson_ghq_nll_blocks_cpp(
      theta,
      blocks$y1, blocks$y2, blocks$n1, blocks$n2, blocks$log_factorial,
      TRUE, gh_sm$nodes, gh_sm$weights
    )
    expect_equal(prepared, raw, tolerance = 1e-12)
  }
})

test_that("vectorized ccc_glmm delta CIs match scalar delta helper", {
  theta <- c(1.1, 0.25, log(0.55))
  vcov_theta <- matrix(
    c(0.02, 0.001, 0.0005,
      0.001, 0.03, 0.002,
      0.0005, 0.002, 0.015),
    nrow = 3L
  )
  include_subject_method <- FALSE
  phi <- 1.2
  m_reps <- 2L

  funcs <- matrixCorr:::.mc_ccc_glmm_delta_functions(include_subject_method, phi, m_reps)
  estimates <- matrixCorr:::.mc_ccc_glmm_metric_vector(theta, include_subject_method, phi, m_reps)
  vec <- matrixCorr:::.mc_delta_ci_vec(
    estimates,
    function(x) matrixCorr:::.mc_ccc_glmm_metric_vector(x, include_subject_method, phi, m_reps),
    theta,
    vcov_theta,
    conf_level = 0.9,
    transform = "fisher"
  )

  scalar <- lapply(names(estimates), function(nm) {
    matrixCorr:::.mc_delta_ci(
      estimate = estimates[[nm]],
      f = funcs[[nm]],
      theta = theta,
      vcov_theta = vcov_theta,
      conf_level = 0.9,
      transform = "fisher"
    )
  })
  names(scalar) <- names(estimates)

  for (nm in names(estimates)) {
    expect_equal(vec[nm, "se"], scalar[[nm]]$se, tolerance = 1e-12, info = nm)
    expect_equal(vec[nm, "lwr.ci"], scalar[[nm]]$lwr.ci, tolerance = 1e-12, info = nm)
    expect_equal(vec[nm, "upr.ci"], scalar[[nm]]$upr.ci, tolerance = 1e-12, info = nm)
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
})
