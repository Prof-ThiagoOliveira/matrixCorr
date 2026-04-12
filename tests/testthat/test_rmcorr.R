manual_rmcorr <- function(x, y, subject, conf_level = 0.95) {
  dat <- data.frame(x = x, y = y, subject = subject, check.names = FALSE)
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  if (!nrow(dat)) {
    return(list(valid = FALSE))
  }

  counts <- table(dat$subject)
  keep <- names(counts[counts >= 2L])
  dat <- dat[dat$subject %in% keep, , drop = FALSE]
  dat$subject <- droplevels(factor(dat$subject))

  n_complete <- nrow(dat)
  n_subjects <- nlevels(dat$subject)
  if (n_subjects < 2L) {
    return(list(valid = FALSE, n_complete = n_complete, n_subjects = n_subjects))
  }

  xc <- dat$x - ave(dat$x, dat$subject)
  yc <- dat$y - ave(dat$y, dat$subject)
  sxx <- sum(xc * xc)
  syy <- sum(yc * yc)
  sxy <- sum(xc * yc)
  r <- sxy / sqrt(sxx * syy)
  slope <- sxy / sxx
  df <- n_complete - n_subjects - 1L
  sse <- syy - (sxy * sxy) / sxx
  if (sse < 0 && sse > -1e-12) {
    sse <- 0
  }
  mse <- sse / df
  se_slope <- sqrt(mse / sxx)
  t_value <- if (se_slope == 0) {
    if (slope == 0) 0 else Inf
  } else {
    slope / se_slope
  }
  p_value <- 2 * stats::pt(-abs(t_value), df = df)

  ci <- c(lower = NA_real_, upper = NA_real_)
  eff_n <- df + 2
  if (abs(r) >= 1) {
    ci[] <- r
  } else if (eff_n > 3) {
    zr <- atanh(r)
    zcrit <- stats::qnorm(0.5 * (1 + conf_level))
    se_z <- 1 / sqrt(eff_n - 3)
    ci <- c(
      lower = tanh(zr - zcrit * se_z),
      upper = tanh(zr + zcrit * se_z)
    )
  }

  list(
    valid = TRUE,
    r = r,
    slope = slope,
    p_value = p_value,
    conf_int = ci,
    df = df,
    n_complete = n_complete,
    n_subjects = n_subjects
  )
}

sim_rmcorr_long <- function(seed = 2026L) {
  set.seed(seed)
  n_subjects <- 18L
  n_rep <- 4L
  subject <- factor(rep(seq_len(n_subjects), each = n_rep))
  subj_x <- rnorm(n_subjects, sd = 1.2)
  subj_y <- rnorm(n_subjects, sd = 1.5)
  signal <- rnorm(n_subjects * n_rep)

  data.frame(
    subject = subject,
    x = subj_x[subject] + signal + rnorm(n_subjects * n_rep, sd = 0.15),
    y = subj_y[subject] + 0.7 * signal + rnorm(n_subjects * n_rep, sd = 0.2),
    z = subj_y[subject] - 0.4 * signal + rnorm(n_subjects * n_rep, sd = 0.25),
    check.names = FALSE
  )
}

test_that("rmcorr pair interface matches manual within-subject formulas", {
  dat <- sim_rmcorr_long()
  ref <- manual_rmcorr(dat$x, dat$y, dat$subject)

  fit <- rmcorr(dat, response = c("x", "y"), subject = "subject")
  fit_positional <- rmcorr(dat[c("x", "y")], dat$subject)

  expect_s3_class(fit, "rmcorr")
  expect_equal(fit$r, ref$r, tolerance = 1e-10)
  expect_equal(fit$estimate, ref$r, tolerance = 1e-10)
  expect_equal(fit$slope, ref$slope, tolerance = 1e-10)
  expect_equal(fit$p_value, ref$p_value, tolerance = 1e-10)
  expect_equal(fit$conf_int, ref$conf_int, tolerance = 1e-10)
  expect_equal(c(fit$lwr, fit$upr), unname(ref$conf_int), tolerance = 1e-10)
  expect_equal(fit$df, ref$df)
  expect_equal(fit$based.on, ref$n_complete)
  expect_equal(fit$n_obs, ref$n_complete)
  expect_equal(fit$n_subjects, ref$n_subjects)
  expect_equal(fit$r, fit_positional$r, tolerance = 1e-10)
  expect_equal(fit$slope, fit_positional$slope, tolerance = 1e-10)

  sm <- summary(fit)
  expect_s3_class(sm, "summary.rmcorr")
  expect_s3_class(sm, "summary.matrixCorr")
  expect_true(any(grepl("Repeated-measures correlation", capture.output(print(fit)))))
  sum_out <- capture.output(matrixCorr:::print.summary.rmcorr(sm))
  expect_true(any(grepl("^Repeated-measures correlation summary:", sum_out)))
  expect_true(any(grepl("obs_per_subject", sum_out, fixed = TRUE)))
  expect_true(any(grepl("ci_method", sum_out, fixed = TRUE)))
  expect_named(sm, c(
    "class", "method", "description", "responses", "n_obs", "n_subjects", "df",
    "estimate", "slope", "p_value", "lwr", "upr", "conf_level",
    "obs_per_subject_min", "obs_per_subject_max", "n_intercepts",
    "ci_method", "ci_width", "valid"
  ))
  expect_equal(fit$r, fit$estimate, tolerance = 1e-12)
  expect_equal(fit$conf_int, c(lower = fit$lwr, upper = fit$upr), tolerance = 1e-12)
  expect_equal(fit$based.on, fit$n_obs)
  expect_false(any(c("r", "conf_int", "based.on", "data_long", "intercepts", "fitted") %in% names(unclass(fit))))
  expect_null(attr(fit, "source_data", exact = TRUE))
  expect_error(plot(fit), "keep_data = TRUE", fixed = TRUE)

  fit_keep <- rmcorr(dat, response = c("x", "y"), subject = "subject", keep_data = TRUE)
  expect_s3_class(plot(fit_keep), "ggplot")
  expect_true(is.list(attr(fit_keep, "source_data", exact = TRUE)))
})

test_that("rmcorr matrix output is symmetric and agrees with pair fits", {
  dat <- sim_rmcorr_long(seed = 44L)

  fit_xy <- rmcorr(dat, response = c("x", "y"), subject = "subject")
  fit_xz <- rmcorr(dat, response = c("x", "z"), subject = "subject")
  fit_mat <- rmcorr(dat, response = c("x", "y", "z"), subject = "subject", n_threads = 1L)

  expect_s3_class(fit_mat, "rmcorr_matrix")
  expect_equal(dim(fit_mat), c(3L, 3L))
  expect_equal(unname(diag(fit_mat)), c(1, 1, 1), tolerance = 1e-12)
  expect_equal(unclass(fit_mat), t(unclass(fit_mat)), tolerance = 1e-10)
  expect_equal(unname(fit_mat["x", "y"]), fit_xy$r, tolerance = 1e-10)
  expect_equal(unname(fit_mat["x", "z"]), fit_xz$r, tolerance = 1e-10)

  diag_attr <- attr(fit_mat, "diagnostics")
  expect_true(is.list(diag_attr))
  expect_equal(unname(diag_attr$slope["x", "y"]), fit_xy$slope, tolerance = 1e-10)
  expect_equal(unname(diag_attr$n_complete["x", "y"]), fit_xy$based.on)
  expect_equal(unname(diag_attr$n_subjects["x", "z"]), fit_xz$n_subjects)

  sm <- summary(fit_mat)
  expect_s3_class(sm, "summary.rmcorr_matrix")
  expect_s3_class(sm, "summary.matrixCorr")
  expect_true(any(grepl("Repeated-measures correlation matrix", capture.output(print(fit_mat)))))
  expect_s3_class(plot(fit_mat), "ggplot")
  expect_null(attr(fit_mat, "source_data", exact = TRUE))
})

test_that("rmcorr matrix stores compact source data only when requested", {
  dat <- sim_rmcorr_long(seed = 77L)
  fit_mat <- rmcorr(
    dat,
    response = c("x", "y", "z"),
    subject = "subject",
    keep_data = TRUE,
    n_threads = 1L
  )

  source_data <- attr(fit_mat, "source_data", exact = TRUE)
  expect_true(is.list(source_data))
  expect_equal(dim(source_data$response), c(nrow(dat), 3L))
  expect_null(colnames(source_data$response))
  expect_equal(source_data$response_names, c("x", "y", "z"))
  expect_equal(source_data$subject_code, as.integer(dat$subject))
  expect_equal(source_data$subject_levels, levels(dat$subject))
  expect_equal(source_data$subject_name, "subject")
  expect_equal(source_data$conf_level, 0.95)
})

test_that("rmcorr pair keeps compact source data only when requested", {
  dat <- sim_rmcorr_long(seed = 79L)
  fit_base <- rmcorr(dat, response = c("x", "y"), subject = "subject")
  fit_keep <- rmcorr(dat, response = c("x", "y"), subject = "subject", keep_data = TRUE)

  expect_null(attr(fit_base, "source_data", exact = TRUE))
  expect_true(is.list(attr(fit_keep, "source_data", exact = TRUE)))
  expect_lt(as.numeric(object.size(fit_base)), as.numeric(object.size(fit_keep)))
})

test_that("rmcorr check_na = FALSE uses pairwise complete cases", {
  dat <- data.frame(
    subject = c(1, 1, 2, 2, 2, 3, 3),
    x = c(1, 2, 3, 4, 5, 1, 2),
    y = c(2, 4, NA, 8, 10, 2, 4),
    z = c(1, 3, 3, 5, NA, 2, 5)
  )

  expect_error(
    rmcorr(dat, response = c("x", "y"), subject = "subject", na_method = "error"),
    "missing, NaN, or infinite values"
  )

  ref_xy <- manual_rmcorr(dat$x, dat$y, dat$subject)
  ref_xz <- manual_rmcorr(dat$x, dat$z, dat$subject)

  expect_warning(
    fit_xy <- rmcorr(dat, response = c("x", "y"), subject = "subject", check_na = FALSE),
    "deprecated"
  )
  expect_warning(
    fit_mat <- rmcorr(dat, response = c("x", "y", "z"), subject = "subject", check_na = FALSE),
    "deprecated"
  )

  expect_equal(fit_xy$based.on, ref_xy$n_complete)
  expect_equal(fit_xy$n_obs, ref_xy$n_complete)
  expect_equal(fit_xy$n_subjects, ref_xy$n_subjects)
  expect_equal(fit_xy$r, ref_xy$r, tolerance = 1e-10)

  diag_attr <- attr(fit_mat, "diagnostics")
  expect_equal(unname(diag_attr$n_complete["x", "y"]), ref_xy$n_complete)
  expect_equal(unname(diag_attr$n_complete["x", "z"]), ref_xz$n_complete)
  expect_equal(unname(fit_mat["x", "z"]), ref_xz$r, tolerance = 1e-10)
})

test_that("rmcorr canonical na_method matches legacy check_na behavior", {
  dat <- data.frame(
    subject = c(1, 1, 2, 2, 2, 3, 3),
    x = c(1, 2, 3, 4, 5, 1, 2),
    y = c(2, 4, NA, 8, 10, 2, 4)
  )

  fit_new <- rmcorr(dat, response = c("x", "y"), subject = "subject", na_method = "pairwise")
  expect_warning(
    fit_old <- rmcorr(dat, response = c("x", "y"), subject = "subject", check_na = FALSE),
    "deprecated"
  )

  expect_equal(fit_new$estimate, fit_old$estimate, tolerance = 1e-12)
  expect_equal(fit_new$n_obs, fit_old$n_obs)
  expect_equal(fit_new$n_subjects, fit_old$n_subjects)
})

test_that("rmcorr validates subject structure", {
  dat <- data.frame(
    subject = rep(1, 4),
    x = c(1, 2, 3, 4),
    y = c(2, 3, 4, 5)
  )

  expect_error(
    rmcorr(dat, response = c("x", "y"), subject = "subject"),
    "at least two distinct subjects"
  )
})
