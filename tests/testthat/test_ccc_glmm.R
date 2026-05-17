make_ccc_glmm_data <- function(n_subjects = 20,
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
  subj_eff <- rnorm(n_subjects, 0, subject_sd)
  df$eta <- 1.5 + subj_eff[as.integer(df$subject)] +
    ifelse(df$method == "B", method_bias, 0)
  df$y <- rpois(nrow(df), exp(df$eta))
  df
}

test_that("ccc_glmm returns matrixCorr CCC object", {
  df <- make_ccc_glmm_data()

  fit <- ccc_glmm(
    df,
    response = "y",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    overdispersion = "none"
  )

  expect_s3_class(fit, "ccc_glmm")
  expect_s3_class(fit, "ccc")
  expect_identical(class(fit), c("ccc_glmm", "ccc"))
  expect_true(is.matrix(fit))
  expect_equal(diag(fit), c(A = 1, B = 1), ignore_attr = TRUE)
})

test_that("ccc_glmm exposes symmetric bounded CCC matrices and components", {
  fit <- ccc_glmm(
    make_ccc_glmm_data(),
    response = "y",
    subject = "subject",
    method = "method",
    replicate = "replicate"
  )

  expect_true(isSymmetric(unclass(fit), check.attributes = FALSE))
  expect_equal(diag(fit), c(A = 1, B = 1), ignore_attr = TRUE)
  expect_true(is.finite(fit[1, 2]))
  expect_gte(fit[1, 2], 0)
  expect_lte(fit[1, 2], 1)

  for (nm in c("rho_ccc_inter", "rho_ccc_intra_method1", "rho_ccc_intra_method2")) {
    val <- attr(fit, nm, exact = TRUE)
    expect_true(is.matrix(val), info = nm)
    expect_true(is.finite(val[1, 2]), info = nm)
    expect_gte(val[1, 2], 0)
    expect_lte(val[1, 2], 1)
  }
})

test_that("ccc_glmm handles overdispersion phi modes", {
  df <- make_ccc_glmm_data()
  fit_none <- ccc_glmm(
    df,
    response = "y",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    overdispersion = "none"
  )
  fit_pearson <- ccc_glmm(
    df,
    response = "y",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    overdispersion = "pearson"
  )

  expect_equal(attr(fit_none, "phi", exact = TRUE)[1, 2], 1)
  phi <- attr(fit_pearson, "phi", exact = TRUE)[1, 2]
  expect_true(is.finite(phi))
  expect_gt(phi, 0)
})

test_that("larger fixed method bias lowers GLMM CCC", {
  low_bias <- ccc_glmm(
    make_ccc_glmm_data(method_bias = 0, seed = 2),
    "y", "subject", "method",
    replicate = "replicate",
    include_subject_method = FALSE
  )
  high_bias <- ccc_glmm(
    make_ccc_glmm_data(method_bias = 1.5, seed = 2),
    "y", "subject", "method",
    replicate = "replicate",
    include_subject_method = FALSE
  )

  expect_lt(high_bias[1, 2], low_bias[1, 2])
})

test_that("ccc_glmm validates count data and implemented model choices", {
  df <- make_ccc_glmm_data()
  bad_non_integer <- df
  bad_non_integer$y[1] <- 1.25
  bad_negative <- df
  bad_negative$y[1] <- -1

  expect_error(
    ccc_glmm(bad_non_integer, "y", "subject", "method"),
    "integer-like"
  )
  expect_error(
    ccc_glmm(bad_negative, "y", "subject", "method"),
    "non-negative"
  )
  expect_error(
    ccc_glmm(df, "y", "subject", "method", family = "binomial"),
    "poisson"
  )
  expect_error(
    ccc_glmm(df, "y", "subject", "method", link = "identity"),
    "log"
  )
})

test_that("summary.ccc_glmm includes GLMM quantities", {
  fit <- ccc_glmm(
    make_ccc_glmm_data(),
    "y", "subject", "method",
    replicate = "replicate",
    overdispersion = "pearson"
  )
  sm <- summary(fit)

  expect_s3_class(sm, "summary.ccc_glmm")
  expect_named(
    sm,
    c(
      "method1", "method2", "estimate",
      "rho_ccc_inter",
      "rho_ccc_intra_method1",
      "rho_ccc_intra_method2",
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
    ),
    ignore.order = TRUE
  )
  expect_true(is.finite(sm$rho_ccc_inter[1]))
  expect_true(is.finite(sm$phi[1]))
})

test_that("ccc_glmm fits subject-method variance component when requested", {
  df <- make_ccc_glmm_data(n_subjects = 16, reps = 2, seed = 11)
  fit <- ccc_glmm(
    df,
    "y", "subject", "method",
    replicate = "replicate",
    include_subject_method = TRUE,
    max_iter = 500
  )

  expect_s3_class(fit, "ccc_glmm")
  expect_true(is.finite(fit[1, 2]))
  expect_true(is.finite(attr(fit, "sigma2_subject_method", exact = TRUE)[1, 2]))
  expect_gte(attr(fit, "sigma2_subject_method", exact = TRUE)[1, 2], 0)
  expect_identical(attr(fit, "include_subject_method", exact = TRUE), TRUE)
})

test_that("ccc_glmm summary includes fitted beta and variance components", {
  fit <- ccc_glmm(
    make_ccc_glmm_data(),
    "y", "subject", "method",
    replicate = "replicate"
  )
  sm <- summary(fit)

  expect_true(is.finite(sm$beta0[1]))
  expect_true(is.finite(sm$beta_method[1]))
  expect_true(is.finite(sm$sigma2_subject[1]))
  expect_true(is.finite(sm$mu[1]))
  expect_true(is.finite(sm$nll[1]))
  expect_true(is.finite(sm$logLik[1]))
})


test_that("plot.ccc works for ccc_glmm through inherited ccc class", {
  fit <- ccc_glmm(
    make_ccc_glmm_data(n_subjects = 12),
    "y", "subject", "method",
    replicate = "replicate"
  )
  p <- plot(fit)
  expect_s3_class(p, "ggplot")
})

test_that("ccc_glmm delta CIs do not alter point estimates", {
  df <- make_ccc_glmm_data(n_subjects = 18, reps = 2, seed = 21)
  fit_no_ci <- ccc_glmm(
    df,
    "y", "subject", "method",
    replicate = "replicate",
    ci = FALSE
  )
  fit_ci <- ccc_glmm(
    df,
    "y", "subject", "method",
    replicate = "replicate",
    ci = TRUE
  )

  expect_equal(
    unclass(as.matrix(fit_ci)),
    unclass(as.matrix(fit_no_ci)),
    tolerance = 0,
    ignore_attr = TRUE
  )
  point_attrs <- c(
    "rho_ccc", "rho_ccc_inter",
    "rho_ccc_intra_method1", "rho_ccc_intra_method2",
    "precision", "accuracy",
    "sigma2_subject", "sigma2_method", "sigma2_subject_method",
    "phi", "mu"
  )
  for (nm in point_attrs) {
    expect_equal(
      attr(fit_ci, nm, exact = TRUE),
      attr(fit_no_ci, nm, exact = TRUE),
      tolerance = 0,
      info = nm
    )
  }
})

test_that("ccc_glmm delta CI attributes exist and are shaped like estimates", {
  df <- make_ccc_glmm_data(n_subjects = 18, reps = 3, method_bias = 0.7, seed = 22)
  fit <- ccc_glmm(
    df,
    "y", "subject", "method",
    replicate = "replicate",
    ci = TRUE
  )

  ci_attrs <- c(
    "rho_ccc_se",
    "rho_ccc_lwr.ci",
    "rho_ccc_upr.ci",
    "rho_ccc_inter_se",
    "rho_ccc_intra_method1_se",
    "rho_ccc_intra_method2_se"
  )
  for (nm in ci_attrs) {
    val <- attr(fit, nm, exact = TRUE)
    expect_true(is.matrix(val), info = nm)
    expect_identical(dim(val), dim(fit), info = nm)
    expect_identical(dimnames(val), dimnames(fit), info = nm)
    expect_true(is.finite(val[1, 2]), info = nm)
  }
  expect_true(isTRUE(attr(fit, "ci", exact = TRUE)))
  expect_identical(attr(fit, "ci_method", exact = TRUE), "delta_fisher")
  sm <- summary(fit)
  expect_true(all(c("se", "lwr.ci", "upr.ci") %in% names(sm)))
  expect_true(is.finite(sm$se[1]))
  expect_true(is.finite(sm$lwr.ci[1]))
  expect_true(is.finite(sm$upr.ci[1]))
})

test_that("ccc_glmm ci = FALSE remains lean", {
  fit <- ccc_glmm(
    make_ccc_glmm_data(seed = 23),
    "y", "subject", "method",
    replicate = "replicate",
    ci = FALSE
  )

  expect_null(attr(fit, "rho_ccc_se", exact = TRUE))
  expect_null(attr(fit, "rho_ccc_lwr.ci", exact = TRUE))
  expect_false(isTRUE(attr(fit, "ci", exact = TRUE)))
})

test_that("ccc_glmm delta CI handles unavailable covariance safely", {
  ci <- .mc_delta_ci(
    estimate = 0.5,
    f = function(theta) theta[[1L]],
    theta = 0,
    vcov_theta = NULL
  )
  expect_true(is.na(ci$se))
  expect_true(is.na(ci$lwr.ci))
  expect_true(is.na(ci$upr.ci))
})

test_that("ccc_glmm inter and intra CIs are metric-specific", {
  df <- make_ccc_glmm_data(n_subjects = 20, reps = 4, method_bias = 1.0, seed = 24)
  fit <- ccc_glmm(
    df,
    "y", "subject", "method",
    replicate = "replicate",
    ci = TRUE
  )

  expect_false(isTRUE(all.equal(
    attr(fit, "rho_ccc_se", exact = TRUE)[1, 2],
    attr(fit, "rho_ccc_inter_se", exact = TRUE)[1, 2]
  )))
  expect_false(isTRUE(all.equal(
    attr(fit, "rho_ccc_intra_method1_se", exact = TRUE)[1, 2],
    attr(fit, "rho_ccc_intra_method2_se", exact = TRUE)[1, 2]
  )))
})

test_that("ccc_glmm Pearson overdispersion CIs warn and compute", {
  df <- make_ccc_glmm_data(n_subjects = 18, reps = 2, seed = 25)
  expect_warning(
    fit <- ccc_glmm(
      df,
      "y", "subject", "method",
      replicate = "replicate",
      overdispersion = "pearson",
      ci = TRUE
    ),
    "treat the Pearson dispersion estimate as fixed"
  )
  expect_true(is.finite(attr(fit, "rho_ccc_se", exact = TRUE)[1, 2]))
})
