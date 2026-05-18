make_cia_rm_perf_data <- function(n_subjects = 10L, n_methods = 3L, n_times = 4L) {
  subjects <- sprintf("s%02d", seq_len(n_subjects))
  methods <- LETTERS[seq_len(n_methods)]
  times <- sprintf("t%02d", seq_len(n_times))
  dat <- expand.grid(
    subject = subjects,
    method = methods,
    time = times,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  s <- match(dat$subject, subjects)
  m <- match(dat$method, methods)
  tt <- match(dat$time, times)
  dat$value <- 5 + 0.3 * s + 0.2 * m + 0.15 * tt + sin(s + m * tt) / 20
  dat
}

test_that("cia_rm no-CI output remains structurally stable", {
  dat2 <- make_cia_rm_perf_data(n_methods = 2L)
  dat3 <- make_cia_rm_perf_data(n_methods = 3L)

  for (estimator in c("vc_constrained", "mom_unconstrained")) {
    for (homogeneous in c(FALSE, TRUE)) {
      fit2 <- cia_rm(dat2, "value", "subject", "method", "time",
                     estimator = estimator, homogeneous = homogeneous,
                     use_message = FALSE)
      fit3 <- cia_rm(dat3, "value", "subject", "method", "time",
                     estimator = estimator, homogeneous = homogeneous,
                     use_message = FALSE)

      expect_s3_class(fit2, "cia_rm")
      expect_identical(class(fit2), c("cia_rm", "cia"))
      expect_named(fit2, c("est", "common", "condition"))
      expect_identical(attr(fit2, "estimator", exact = TRUE), estimator)
      expect_identical(attr(fit2, "homogeneous", exact = TRUE), homogeneous)

      expect_s3_class(fit3, "cia_rm")
      expect_true(all(c("overall", "overall.common") %in% names(fit3)))
      expect_identical(dimnames(fit3$est)[[1L]], attr(fit3, "method_levels", exact = TRUE))
      expect_identical(names(fit3$overall), attr(fit3, "time_levels", exact = TRUE))
    }
  }
})

test_that("cia_rm delta CI remains structurally stable", {
  dat2 <- make_cia_rm_perf_data(n_methods = 2L)
  dat3 <- make_cia_rm_perf_data(n_methods = 3L)

  for (estimator in c("vc_constrained", "mom_unconstrained")) {
    for (homogeneous in c(FALSE, TRUE)) {
      fit2 <- cia_rm(dat2, "value", "subject", "method", "time",
                     ci = TRUE, inference = "delta", estimator = estimator,
                     homogeneous = homogeneous, use_message = FALSE)
      fit3 <- cia_rm(dat3, "value", "subject", "method", "time",
                     ci = TRUE, inference = "delta", estimator = estimator,
                     homogeneous = homogeneous, use_message = FALSE)

      expect_identical(class(fit2), c("cia_rm", "cia_ci", "cia"))
      expect_identical(attr(fit2, "ci.method", exact = TRUE), "delta_normal")
      expect_true(all(c("se", "lwr.ci", "upr.ci", "common.se") %in% names(fit2)))
      expect_true(all(c("overall.se", "overall.lwr.ci", "overall.upr.ci") %in% names(fit3)))
      expect_identical(dim(fit2$se), dim(fit2$est))
      expect_identical(dimnames(fit2$lwr.ci), dimnames(fit2$est))
    }
  }
})

test_that("cia_rm bootstrap CI is reproducible with fixed seed", {
  dat <- make_cia_rm_perf_data(n_methods = 3L)

  for (estimator in c("vc_constrained", "mom_unconstrained")) {
    for (homogeneous in c(FALSE, TRUE)) {
      fit1 <- cia_rm(dat, "value", "subject", "method", "time",
                     ci = TRUE, inference = "bootstrap", estimator = estimator,
                     homogeneous = homogeneous, B = 49L, seed = 77L,
                     use_message = FALSE)
      fit2 <- cia_rm(dat, "value", "subject", "method", "time",
                     ci = TRUE, inference = "bootstrap", estimator = estimator,
                     homogeneous = homogeneous, B = 49L, seed = 77L,
                     use_message = FALSE)

      expect_equal(fit1, fit2, tolerance = 1e-12)
      expect_identical(attr(fit1, "ci.method", exact = TRUE), "subject_bootstrap_percentile")
      expect_identical(attr(fit1, "bootstrap", exact = TRUE), list(n_successful = 49L, B = 49L))
      expect_true(all(c("overall.lwr.ci", "overall.upr.ci", "overall.common.lwr.ci") %in% names(fit1)))
    }
  }
})

test_that("cia_rm time metadata remains stable for categorical numeric and Date time", {
  categorical <- make_cia_rm_perf_data()
  numeric <- categorical
  numeric$time <- match(numeric$time, sort(unique(numeric$time))) * 5
  date <- categorical
  date$time <- as.Date("2024-01-01") + match(date$time, sort(unique(date$time))) - 1L

  fit_cat <- cia_rm(categorical, "value", "subject", "method", "time", use_message = FALSE)
  fit_num <- cia_rm(numeric, "value", "subject", "method", "time", use_message = FALSE)
  fit_date <- cia_rm(date, "value", "subject", "method", "time", use_message = FALSE)

  expect_identical(attr(fit_cat, "time_scale", exact = TRUE), "categorical")
  expect_identical(attr(fit_num, "time_scale", exact = TRUE), "continuous")
  expect_type(attr(fit_num, "time_plot_values", exact = TRUE), "double")
  expect_identical(attr(fit_date, "time_scale", exact = TRUE), "continuous")
  expect_s3_class(attr(fit_date, "time_plot_values", exact = TRUE), "Date")
})

test_that("cia_rm error behaviour remains stable", {
  dat <- make_cia_rm_perf_data()
  missing_response <- dat
  missing_response$value[[1L]] <- NA_real_
  missing_key <- dat
  missing_key$subject[[1L]] <- NA_character_
  one_method <- subset(dat, method == "A")
  one_time <- subset(dat, time == "t01")
  one_subject <- subset(dat, subject == "s01")
  unbalanced <- dat[-1L, , drop = FALSE]
  duplicated <- rbind(dat, dat[1L, , drop = FALSE])

  expect_error(cia_rm(missing_response, "value", "subject", "method", "time", use_message = FALSE),
               "finite numeric")
  expect_error(cia_rm(missing_key, "value", "subject", "method", "time", use_message = FALSE),
               "must not contain missing")
  expect_error(cia_rm(one_method, "value", "subject", "method", "time", use_message = FALSE),
               "at least two method")
  expect_error(cia_rm(one_time, "value", "subject", "method", "time", use_message = FALSE),
               "at least two time")
  expect_error(cia_rm(one_subject, "value", "subject", "method", "time", use_message = FALSE),
               "at least two subjects")
  expect_error(cia_rm(unbalanced, "value", "subject", "method", "time", use_message = FALSE),
               "each subject-method-time cell must contain exactly one observation")
  expect_error(cia_rm(duplicated, "value", "subject", "method", "time", use_message = FALSE),
               "each subject-method-time cell must contain exactly one observation")
  expect_error(cia_rm(dat, "value", "subject", "method", "time", ci = TRUE,
                      inference = "none", use_message = FALSE),
               "must not be")
  expect_error(cia_rm(dat, "value", "subject", "method", "time", ci = TRUE,
                      inference = "bootstrap", B = 1L, use_message = FALSE),
               "must be >= 2")
})

test_that("cia_rm print summary and plot compatibility is retained", {
  dat <- make_cia_rm_perf_data(n_methods = 3L)
  fit <- cia_rm(dat, "value", "subject", "method", "time",
                ci = TRUE, inference = "delta", homogeneous = TRUE,
                use_message = FALSE)

  expect_output(print(fit), "Repeated-measures")
  sm <- summary(fit)
  expect_s3_class(sm, "summary.cia_rm")
  expect_output(print(sm), "Repeated-measures")

  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  }, add = TRUE)
  expect_silent(plot(fit))
})
