.test_icc_reference <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  k <- ncol(X)
  stopifnot(k == 2L, n >= 2L)

  row_means <- rowMeans(X)
  grand_mean <- mean(X)
  col_means <- colMeans(X)

  ss_subject <- k * sum((row_means - grand_mean)^2)
  ss_rater <- n * sum((col_means - grand_mean)^2)
  ss_total <- sum((X - grand_mean)^2)
  ss_error <- ss_total - ss_subject - ss_rater
  ss_within <- ss_total - ss_subject

  ms_subject <- ss_subject / (n - 1)
  ms_rater <- ss_rater / (k - 1)
  ms_error <- ss_error / ((n - 1) * (k - 1))
  ms_within <- ss_within / (n * (k - 1))

  icc1 <- (ms_subject - ms_within) / (ms_subject + (k - 1) * ms_within)
  icc1k <- (ms_subject - ms_within) / ms_subject
  icc2 <- (ms_subject - ms_error) /
    (ms_subject + (k - 1) * ms_error + k * (ms_rater - ms_error) / n)
  icc2k <- (ms_subject - ms_error) / (ms_subject + (ms_rater - ms_error) / n)
  icc3 <- (ms_subject - ms_error) / (ms_subject + (k - 1) * ms_error)
  icc3k <- (ms_subject - ms_error) / ms_subject

  f11 <- ms_subject / ms_within
  df11n <- n - 1
  df11d <- n * (k - 1)
  f1l <- f11 / qf(0.975, df11n, df11d)
  f1u <- f11 * qf(0.975, df11d, df11n)
  l1 <- (f1l - 1) / (f1l + (k - 1))
  u1 <- (f1u - 1) / (f1u + (k - 1))

  f21 <- ms_subject / ms_error
  df21n <- n - 1
  df21d <- (n - 1) * (k - 1)
  f3l <- f21 / qf(0.975, df21n, df21d)
  f3u <- f21 * qf(0.975, df21d, df21n)
  l3_single <- (f3l - 1) / (f3l + (k - 1))
  u3_single <- (f3u - 1) / (f3u + (k - 1))

  fj <- ms_rater / ms_error
  vn <- (k - 1) * (n - 1) * ((k * icc2 * fj + n * (1 + (k - 1) * icc2) - k * icc2)^2)
  vd <- (n - 1) * k^2 * icc2^2 * fj^2 + (n * (1 + (k - 1) * icc2) - k * icc2)^2
  v <- vn / vd
  f2u <- qf(0.975, n - 1, v)
  f2l <- qf(0.975, v, n - 1)
  l2_single <- n * (ms_subject - f2u * ms_error) /
    (f2u * (k * ms_rater + (k * n - k - n) * ms_error) + n * ms_subject)
  u2_single <- n * (f2l * ms_subject - ms_error) /
    (k * ms_rater + (k * n - k - n) * ms_error + n * f2l * ms_subject)
  l2_avg <- l2_single * k / (1 + l2_single * (k - 1))
  u2_avg <- u2_single * k / (1 + u2_single * (k - 1))

  list(
    ICC1 = icc1,
    ICC1k = icc1k,
    ICC2 = icc2,
    ICC2k = icc2k,
    ICC3 = icc3,
    ICC3k = icc3k,
    ICC1_lwr = l1,
    ICC1_upr = u1,
    ICC2_lwr = l2_single,
    ICC2_upr = u2_single,
    ICC2k_lwr = l2_avg,
    ICC2k_upr = u2_avg,
    ICC3_lwr = l3_single,
    ICC3_upr = u3_single
  )
}

.test_icc_overall_reference <- function(X, conf_level = 0.95) {
  X <- as.matrix(X)
  n <- nrow(X)
  k <- ncol(X)
  stopifnot(n >= 2L, k >= 2L)

  row_means <- rowMeans(X)
  grand_mean <- mean(X)
  col_means <- colMeans(X)

  ss_subject <- k * sum((row_means - grand_mean)^2)
  ss_rater <- n * sum((col_means - grand_mean)^2)
  ss_total <- sum((X - grand_mean)^2)
  ss_error <- ss_total - ss_subject - ss_rater
  ss_within <- ss_total - ss_subject

  ms_subject <- ss_subject / (n - 1)
  ms_rater <- ss_rater / (k - 1)
  ms_error <- ss_error / ((n - 1) * (k - 1))
  ms_within <- ss_within / (n * (k - 1))

  icc1 <- (ms_subject - ms_within) / (ms_subject + (k - 1) * ms_within)
  icc2 <- (ms_subject - ms_error) /
    (ms_subject + (k - 1) * ms_error + k * (ms_rater - ms_error) / n)
  icc3 <- (ms_subject - ms_error) / (ms_subject + (k - 1) * ms_error)
  icc1k <- (ms_subject - ms_within) / ms_subject
  icc2k <- (ms_subject - ms_error) / (ms_subject + (ms_rater - ms_error) / n)
  icc3k <- (ms_subject - ms_error) / ms_subject

  f11 <- ms_subject / ms_within
  df11n <- n - 1
  df11d <- n * (k - 1)
  p11 <- pf(f11, df11n, df11d, lower.tail = FALSE)

  f21 <- ms_subject / ms_error
  df21n <- n - 1
  df21d <- (n - 1) * (k - 1)
  p21 <- pf(f21, df21n, df21d, lower.tail = FALSE)

  alpha <- 1 - conf_level
  f1l <- f11 / qf(1 - alpha / 2, df11n, df11d)
  f1u <- f11 * qf(1 - alpha / 2, df11d, df11n)
  l1 <- (f1l - 1) / (f1l + (k - 1))
  u1 <- (f1u - 1) / (f1u + (k - 1))

  f3l <- f21 / qf(1 - alpha / 2, df21n, df21d)
  f3u <- f21 * qf(1 - alpha / 2, df21d, df21n)
  l3 <- (f3l - 1) / (f3l + (k - 1))
  u3 <- (f3u - 1) / (f3u + (k - 1))

  fj <- ms_rater / ms_error
  vn <- (k - 1) * (n - 1) * ((k * icc2 * fj + n * (1 + (k - 1) * icc2) - k * icc2)^2)
  vd <- (n - 1) * k^2 * icc2^2 * fj^2 + (n * (1 + (k - 1) * icc2) - k * icc2)^2
  v <- vn / vd
  f2u <- qf(1 - alpha / 2, n - 1, v)
  f2l <- qf(1 - alpha / 2, v, n - 1)
  l2 <- n * (ms_subject - f2u * ms_error) /
    (f2u * (k * ms_rater + (k * n - k - n) * ms_error) + n * ms_subject)
  u2 <- n * (f2l * ms_subject - ms_error) /
    (k * ms_rater + (k * n - k - n) * ms_error + n * f2l * ms_subject)

  data.frame(
    coefficient = c("ICC1", "ICC2", "ICC3", "ICC1k", "ICC2k", "ICC3k"),
    estimate = c(icc1, icc2, icc3, icc1k, icc2k, icc3k),
    statistic = c(f11, f21, f21, f11, f21, f21),
    df1 = c(df11n, df21n, df21n, df11n, df21n, df21n),
    df2 = c(df11d, df21d, df21d, df11d, df21d, df21d),
    p_value = c(p11, p21, p21, p11, p21, p21),
    lwr = c(l1, l2, l3, 1 - 1 / f1l, l2 * k / (1 + l2 * (k - 1)), 1 - 1 / f3l),
    upr = c(u1, u2, u3, 1 - 1 / f1u, u2 * k / (1 + u2 * (k - 1)), 1 - 1 / f3u),
    stringsAsFactors = FALSE
  )
}

test_that("classical pairwise ICC matches the standard two-rater ANOVA formulas", {

  set.seed(101)
  n <- 80L
  subj <- rnorm(n, sd = 1.1)
  A <- subj + rnorm(n, sd = 0.3)
  B <- 0.4 + subj + rnorm(n, sd = 0.5)
  C <- -0.2 + 0.9 * subj + rnorm(n, sd = 0.4)
  X <- cbind(A = A, B = B, C = C)

  pairs <- list(
    c("A", "B"),
    c("A", "C"),
    c("B", "C")
  )

  forms <- list(
    list(model = "oneway", type = "consistency", unit = "single", ref = "ICC1"),
    list(model = "oneway", type = "consistency", unit = "average", ref = "ICC1k"),
    list(model = "twoway_random", type = "agreement", unit = "single", ref = "ICC2"),
    list(model = "twoway_random", type = "agreement", unit = "average", ref = "ICC2k"),
    list(model = "twoway_random", type = "consistency", unit = "single", ref = "ICC3"),
    list(model = "twoway_random", type = "consistency", unit = "average", ref = "ICC3k"),
    list(model = "twoway_mixed", type = "agreement", unit = "single", ref = "ICC2"),
    list(model = "twoway_mixed", type = "consistency", unit = "single", ref = "ICC3")
  )

  for (spec in forms) {
    fit <- icc(X,
      model = spec$model,
      type = spec$type,
      unit = spec$unit
    )
    for (pair in pairs) {
      ref <- .test_icc_reference(X[, pair, drop = FALSE])
      expect_equal(
        unname(fit[pair[1], pair[2]]),
        unname(ref[[spec$ref]]),
        tolerance = 1e-8
      )
    }
  }
})

test_that("classical ICC confidence intervals align with the standard reference output", {
  set.seed(202)
  n <- 60L
  subj <- rnorm(n, sd = 1)
  X <- data.frame(
    A = subj + rnorm(n, sd = 0.4),
    B = 0.3 + subj + rnorm(n, sd = 0.6)
  )

  fit <- icc(X,
    model = "twoway_random",
    type = "agreement",
    unit = "single",
    ci = TRUE
  )
  ref <- .test_icc_reference(X)

  expect_s3_class(fit, "icc_ci")
  expect_equal(fit$est["A", "B"], ref$ICC2, tolerance = 1e-8)
  expect_equal(fit$lwr.ci["A", "B"], ref$ICC2_lwr, tolerance = 1e-6)
  expect_equal(fit$upr.ci["A", "B"], ref$ICC2_upr, tolerance = 1e-6)
})

test_that("pairwise ICC stores CI payload only when requested", {
  set.seed(205)
  X <- data.frame(
    A = rnorm(30),
    B = rnorm(30),
    C = rnorm(30)
  )

  fit <- icc(X,
    model = "twoway_random",
    type = "agreement",
    unit = "single",
    ci = FALSE
  )
  fit_ci <- icc(X,
    model = "twoway_random",
    type = "agreement",
    unit = "single",
    ci = TRUE
  )

  expect_false(is.list(fit))
  expect_s3_class(fit_ci, "icc_ci")
  expect_identical(names(fit_ci), c("est", "lwr.ci", "upr.ci"))
  expect_false("n_complete" %in% names(fit_ci))
  expect_true(is.matrix(attr(fit, "diagnostics")$n_complete))
  expect_true(is.matrix(attr(fit_ci, "diagnostics")$n_complete))
})

test_that("classical ICC argument combinations and matrix behaviour are consistent", {
  set.seed(303)
  n <- 40L
  subj <- rnorm(n, sd = 1.2)
  X <- data.frame(
    A = subj + rnorm(n, sd = 0.2),
    B = 0.8 + subj + rnorm(n, sd = 0.2),
    C = rep(4, n)
  )

  expect_error(
    icc(X[, 1:2], model = "oneway", type = "agreement"),
    "cannot be.*agreement",
    ignore.case = TRUE
  )

  cons_single <- icc(X, model = "twoway_random", type = "consistency", unit = "single")
  agree_single <- icc(X[, 1:2], model = "twoway_random", type = "agreement", unit = "single")
  cons_average <- icc(X[, 1:2], model = "twoway_random", type = "consistency", unit = "average")

  expect_gt(cons_single["A", "B"], agree_single["A", "B"])
  expect_gte(cons_average["A", "B"], cons_single["A", "B"])
  expect_true(is.na(cons_single["A", "C"]))
  expect_equal(as.numeric(diag(cons_single)), c(1, 1, 1))
  expect_true(isSymmetric(unclass(cons_single)))

  fit_ci <- icc(X[, 1:2], model = "twoway_random", type = "agreement", unit = "single", ci = TRUE)
  sm <- summary(fit_ci)
  expect_s3_class(sm, "summary.icc")
  expect_true(all(c("item1", "item2", "estimate", "lwr", "upr") %in% names(sm)))

  txt_fit <- capture.output(print(fit_ci))
  txt_sm <- capture.output(print(sm))
  expect_true(any(grepl("^Intraclass correlation matrix$", txt_fit)))
  expect_true(any(grepl("^Intraclass correlation summary$", txt_sm)))
})

test_that("classical ICC respects pairwise-complete handling", {
  X <- data.frame(
    A = c(1, 2, 3, 4, 5, NA),
    B = c(1.2, 2.2, 3.1, 4.3, NA, 6.0),
    C = c(0.9, 2.0, 2.8, 4.1, 5.2, 6.1)
  )

  expect_error(icc(X, na_method = "error"), "Missing")

  fit <- icc(X, na_method = "pairwise")
  diag_n <- attr(fit, "diagnostics")$n_complete
  expect_equal(diag_n["A", "B"], 4L)
  expect_true(is.finite(fit["A", "B"]))
  expect_equal(as.numeric(diag(fit)), c(1, 1, 1))
})

test_that("overall ICC reproduces the standard all-rater ANOVA table", {
  sf <- matrix(c(
    9, 2, 5, 8,
    6, 1, 3, 2,
    8, 4, 6, 8,
    7, 1, 2, 6,
    10, 5, 6, 9,
    6, 2, 4, 7
  ), ncol = 4, byrow = TRUE)
  colnames(sf) <- paste0("J", 1:4)
  rownames(sf) <- paste0("S", 1:6)

  fit <- icc(
    sf,
    model = "twoway_random",
    type = "agreement",
    unit = "single",
    scope = "overall",
    ci = TRUE
  )
  ref <- .test_icc_overall_reference(sf)

  expect_s3_class(fit, "icc_overall")
  expect_true(all(c("coefficients", "anova", "mean_squares") %in% names(fit)))
  expect_identical(attr(fit, "selected_coefficient"), "ICC2")

  got <- fit$coefficients
  expect_equal(got$estimate, ref$estimate, tolerance = 1e-8)
  expect_equal(got$statistic, ref$statistic, tolerance = 1e-8)
  expect_equal(got$df1, ref$df1, tolerance = 1e-8)
  expect_equal(got$df2, ref$df2, tolerance = 1e-8)
  expect_equal(got$p_value, ref$p_value, tolerance = 1e-8)
  expect_equal(got$lwr, ref$lwr, tolerance = 1e-8)
  expect_equal(got$upr, ref$upr, tolerance = 1e-8)
  expect_true(got$selected[got$coefficient == "ICC2"])
})

test_that("overall ICC summary and display follow the native object layout", {
  set.seed(404)
  n <- 24
  subj <- rnorm(n, sd = 1)
  X <- data.frame(
    A = subj + rnorm(n, sd = 0.2),
    B = 0.2 + subj + rnorm(n, sd = 0.4),
    C = -0.1 + subj + rnorm(n, sd = 0.5)
  )

  fit <- icc(
    X,
    model = "twoway_mixed",
    type = "consistency",
    unit = "average",
    scope = "overall",
    ci = TRUE
  )
  sm <- summary(fit)

  expect_s3_class(sm, "summary.icc_overall")
  expect_true(all(c(
    "coefficient", "label", "estimate", "statistic",
    "df1", "df2", "p_value", "lwr", "upr", "selected"
  ) %in% names(sm)))
  expect_identical(attr(sm, "selected_coefficient"), "ICC3k")
  expect_true(is.data.frame(attr(sm, "anova")))

  txt_fit <- capture.output(print(fit))
  txt_sm <- capture.output(print(sm))
  expect_true(any(grepl("^Overall intraclass correlation$", txt_fit)))
  expect_true(any(grepl("^Overall intraclass correlation summary$", txt_sm)))
  expect_true(any(grepl("^ICC coefficients$", txt_sm)))
  expect_true(any(grepl("^ANOVA decomposition$", txt_sm)))
})

test_that("overall ICC drops incomplete rows when requested", {
  X <- data.frame(
    A = c(1, 2, 3, 4, 5, NA),
    B = c(1.2, 2.2, 3.1, 4.3, NA, 6.0),
    C = c(0.9, 2.0, 2.8, 4.1, 5.2, 6.1)
  )

  expect_error(
    icc(X, scope = "overall", na_method = "error"),
    "Missing"
  )

  fit <- icc(X, scope = "overall", na_method = "pairwise")
  expect_identical(attr(fit, "diagnostics")$n_complete, 4L)
  expect_identical(attr(fit, "diagnostics")$dropped_rows, 2L)
})

sim_icc_rm_long <- function(seed, sig_subject = 1.2, sig_method = 0.2,
                            sig_error = 0.4, bias_b = 0.0,
                            n_subj = 120L, n_time = 4L) {
  set.seed(seed)
  id <- factor(rep(seq_len(n_subj), each = 2L * n_time))
  method <- factor(rep(rep(c("A", "B"), each = n_time), times = n_subj))
  time <- factor(rep(rep(seq_len(n_time), times = 2L), times = n_subj))
  subj <- rnorm(n_subj, 0, sqrt(sig_subject))
  subj_method <- matrix(rnorm(n_subj * 2L, 0, sqrt(sig_method)), n_subj, 2L)
  y <- numeric(length(id))

  for (s in seq_len(n_subj)) {
    for (m in 1:2) {
      idx <- which(as.integer(id) == s & as.integer(method) == m)
      y[idx] <- subj[s] + subj_method[s, m] + (m == 2L) * bias_b + rnorm(n_time, 0, sqrt(sig_error))
    }
  }

  data.frame(y = y, id = id, method = method, time = time)
}

test_that("repeated-measures ICC behaves sensibly as variance components change", {
  dat_high <- sim_icc_rm_long(seed = 1, sig_subject = 1.4, sig_method = 0.1, sig_error = 0.15)
  dat_low <- sim_icc_rm_long(seed = 2, sig_subject = 0.6, sig_method = 0.3, sig_error = 1.2)

  fit_high <- icc_rm_reml(
    dat_high, "y", "id",
    method = "method", time = "time",
    ci = TRUE, ci_mode = "raw"
  )
  fit_low <- icc_rm_reml(
    dat_low, "y", "id",
    method = "method", time = "time",
    ci = TRUE, ci_mode = "raw"
  )

  expect_s3_class(fit_high, "icc_rm_reml")
  expect_gt(fit_high$est["A", "B"], 0.7)
  expect_lt(fit_low$est["A", "B"], fit_high$est["A", "B"])
  expect_true(fit_high$lwr.ci["A", "B"] <= fit_high$est["A", "B"])
  expect_true(fit_high$upr.ci["A", "B"] >= fit_high$est["A", "B"])
})

test_that("repeated-measures ICC uses the subject argument name", {
  dat <- sim_icc_rm_long(seed = 11, sig_subject = 1.0, sig_method = 0.2, sig_error = 0.3)

  fit_named <- icc_rm_reml(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    ci = FALSE
  )
  fit_positional <- icc_rm_reml(
    dat,
    "y",
    "id",
    method = "method",
    time = "time",
    ci = FALSE
  )

  expect_equal(unname(fit_named["A", "B"]), unname(fit_positional["A", "B"]))
})

test_that("repeated-measures ICC honors n_threads without changing estimates", {
  dat <- sim_icc_rm_long(seed = 12, sig_subject = 1.0, sig_method = 0.2, sig_error = 0.3)

  fit1 <- icc_rm_reml(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    n_threads = 1L,
    ci = FALSE
  )
  fit2 <- icc_rm_reml(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    n_threads = 2L,
    ci = FALSE
  )

  expect_equal(unname(fit1["A", "B"]), unname(fit2["A", "B"]), tolerance = 1e-12)
})

test_that("repeated-measures ICC matches the fitted variance-component denominator", {
  dat <- sim_icc_rm_long(seed = 3, sig_subject = 1.0, sig_method = 0.35, sig_error = 0.4, n_time = 1L)

  fit_cons <- icc_rm_reml(
    dat,
    "y",
    "id",
    method = "method",
    time = NULL,
    type = "consistency",
    vc_select = "none",
    include_subj_method = TRUE,
    ci = FALSE
  )
  fit_cons_nosm <- icc_rm_reml(
    dat,
    "y",
    "id",
    method = "method",
    time = NULL,
    type = "consistency",
    vc_select = "none",
    include_subj_method = FALSE,
    ci = FALSE
  )
  fit_agree <- icc_rm_reml(
    dat,
    "y",
    "id",
    method = "method",
    time = NULL,
    type = "agreement",
    vc_select = "none",
    include_subj_method = TRUE,
    ci = FALSE
  )

  sa <- attr(fit_cons, "sigma2_subject")["A", "B"]
  sab <- attr(fit_cons, "sigma2_subject_method")["A", "B"]
  se <- attr(fit_cons, "sigma2_error")["A", "B"]
  sa_nosm <- attr(fit_cons_nosm, "sigma2_subject")["A", "B"]
  se_nosm <- attr(fit_cons_nosm, "sigma2_error")["A", "B"]
  sb <- attr(fit_agree, "SB")["A", "B"]

  expect_equal(fit_cons["A", "B"], sa / (sa + sab + se), tolerance = 1e-6)
  expect_equal(fit_cons_nosm["A", "B"], sa_nosm / (sa_nosm + se_nosm), tolerance = 1e-6)
  expect_equal(fit_agree["A", "B"], sa / (sa + sab + se + sb), tolerance = 1e-6)
  expect_gte(fit_cons_nosm["A", "B"], fit_cons["A", "B"])
  expect_gte(fit_cons["A", "B"], fit_agree["A", "B"])
})

test_that("repeated-measures ICC summaries expose the native REML payload", {
  dat <- sim_icc_rm_long(seed = 4, sig_subject = 0.9, sig_method = 0.2, sig_error = 0.3, bias_b = 0.4)

  fit <- icc_rm_reml(
    dat,
    "y",
    "id",
    method = "method",
    time = "time",
    type = "agreement",
    ci = TRUE,
    ar = "none"
  )
  sm <- summary(fit)

  expect_s3_class(sm, "summary.icc_rm_reml")
  expect_true(all(c(
    "item1", "item2", "estimate", "lwr", "upr", "n_subjects", "n_obs", "se_icc",
    "sigma2_subject", "sigma2_subject_method", "sigma2_subject_time",
    "sigma2_error", "SB", "residual_model"
  ) %in% names(sm)))
  expect_identical(sm$residual_model[1], "iid")

  txt_fit <- capture.output(print(fit))
  txt_sm <- capture.output(print(sm))
  expect_true(any(grepl("^Repeated-measures intraclass correlation matrix$", txt_fit)))
  expect_true(any(grepl("^ICC estimates$", txt_sm)))
  expect_true(any(grepl("^Variance components$", txt_sm)))
})

test_that("repeated-measures ICC carries AR(1) diagnostics consistently", {
  dat <- sim_icc_rm_long(seed = 5, sig_subject = 1.1, sig_method = 0.1, sig_error = 0.35, n_time = 5L)
  fit <- icc_rm_reml(
    dat,
    "y",
    "id",
    method = "method",
    time = "time",
    ar = "ar1",
    ar_rho = 0.4,
    ci = FALSE
  )

  expect_true(is.matrix(attr(fit, "ar_rho")))
  expect_equal(attr(fit, "ar_rho")["A", "B"], 0.4, tolerance = 1e-12)
  sm <- summary(fit)
  expect_true("ar1_rho" %in% names(sm))
})
