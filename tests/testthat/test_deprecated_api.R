test_that("bland_altman wrapper warns and preserves `two` mapping", {
  set.seed(1)
  x <- rnorm(40)
  y <- x + rnorm(40, sd = 0.5)

  expect_warning(
    bland_altman(x, y, two = 2),
    "removed in 2.0.0"
  )

  old_fit <- suppressWarnings(bland_altman(x, y, two = 2))
  new_fit <- ba(x, y, loa_multiplier = 2)

  expect_equal(old_fit$mean.diffs, new_fit$mean.diffs, tolerance = 1e-12)
  expect_equal(old_fit$lower.limit, new_fit$lower.limit, tolerance = 1e-12)
  expect_equal(old_fit$upper.limit, new_fit$upper.limit, tolerance = 1e-12)
  expect_equal(old_fit$loa_multiplier, new_fit$loa_multiplier)
})

test_that("repeated-measures wrappers warn and preserve outputs", {
  set.seed(2)
  n_subject <- 8L
  n_time <- 3L
  subject <- rep(seq_len(n_subject), each = n_time)
  time <- rep(seq_len(n_time), times = n_subject)
  truth <- rnorm(n_subject, 10, 1)[subject]
  mA <- truth + rnorm(length(truth), sd = 0.3)
  mB <- truth + 0.2 + rnorm(length(truth), sd = 0.4)

  dat <- rbind(
    data.frame(y = mA, subject = subject, method = "A", time = time),
    data.frame(y = mB, subject = subject, method = "B", time = time)
  )

  expect_warning(
    bland_altman_repeated(
      data = dat,
      response = "y",
      subject = "subject",
      method = "method",
      time = "time"
    ),
    "removed in 2.0.0"
  )

  old_ba <- suppressWarnings(
    bland_altman_repeated(
      data = dat,
      response = "y",
      subject = "subject",
      method = "method",
      time = "time",
      two = 2
    )
  )
  new_ba <- ba_rm(
    data = dat,
    response = "y",
    subject = "subject",
    method = "method",
    time = "time",
    loa_multiplier = 2
  )

  expect_equal(old_ba$mean.diffs, new_ba$mean.diffs, tolerance = 1e-12)
  expect_equal(old_ba$lower.limit, new_ba$lower.limit, tolerance = 1e-12)
  expect_equal(old_ba$upper.limit, new_ba$upper.limit, tolerance = 1e-12)

  df <- expand.grid(subject = 1:10, time = 1:2, method = c("A", "B", "C"))
  df$y <- rnorm(nrow(df), mean = match(df$method, c("A", "B", "C")), sd = 1)

  expect_warning(
    ccc_pairwise_u_stat(
      df,
      response = "y",
      method = "method",
      subject = "subject",
      time = "time"
    ),
    "removed in 2.0.0"
  )

  old_u <- suppressWarnings(
    ccc_pairwise_u_stat(
      df,
      response = "y",
      method = "method",
      subject = "subject",
      time = "time"
    )
  )
  new_u <- ccc_rm_ustat(
    df,
    response = "y",
    method = "method",
    subject = "subject",
    time = "time"
  )

  expect_equal(unclass(old_u), unclass(new_u), tolerance = 1e-12)

  expect_warning(
    ccc_lmm_reml(
      df,
      response = "y",
      rind = "subject",
      method = "method",
      time = "time",
      ci = TRUE
    ),
    "removed in 2.0.0"
  )

  old_reml <- suppressWarnings(
    ccc_lmm_reml(
      df,
      response = "y",
      rind = "subject",
      method = "method",
      time = "time",
      ci = TRUE
    )
  )
  new_reml <- ccc_rm_reml(
    df,
    response = "y",
    rind = "subject",
    method = "method",
    time = "time",
    ci = TRUE
  )

  expect_equal(old_reml$est, new_reml$est, tolerance = 1e-12)
})

test_that("correlation wrappers warn and preserve outputs", {
  set.seed(123)
  X <- matrix(rnorm(120), ncol = 4)
  X[1, 1] <- 8

  expect_warning(biweight_mid_corr(X), "removed in 2.0.0")
  old_bicor <- suppressWarnings(biweight_mid_corr(X))
  new_bicor <- bicor(X)
  expect_equal(unclass(old_bicor), unclass(new_bicor), tolerance = 1e-12)

  expect_warning(distance_corr(X), "removed in 2.0.0")
  old_dcor <- suppressWarnings(distance_corr(X))
  new_dcor <- dcor(X)
  expect_equal(unclass(old_dcor), unclass(new_dcor), tolerance = 1e-12)
})

test_that("partial_correlation wrapper preserves the old default method", {
  set.seed(42)
  X <- matrix(rnorm(100 * 5), ncol = 5)

  expect_warning(partial_correlation(X), "pre-1.0 default `method = \"oas\"`")

  old_fit <- suppressWarnings(partial_correlation(X))
  new_fit <- pcorr(X, method = "oas")

  expect_equal(old_fit$pcor, new_fit$pcor, tolerance = 1e-12)
  expect_equal(old_fit$method, "oas")
})
