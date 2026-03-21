test_that("rmcorr viewer helpers preserve source data and rebuild pair fits", {
  set.seed(2026)
  n_subjects <- 12
  n_rep <- 4
  subject <- rep(seq_len(n_subjects), each = n_rep)
  subj_eff_x <- rnorm(n_subjects, sd = 1.2)
  subj_eff_y <- rnorm(n_subjects, sd = 1.5)
  within_signal <- rnorm(n_subjects * n_rep)

  dat <- data.frame(
    subject = subject,
    x = subj_eff_x[subject] + within_signal + rnorm(n_subjects * n_rep, sd = 0.2),
    y = subj_eff_y[subject] + 0.7 * within_signal + rnorm(n_subjects * n_rep, sd = 0.3),
    z = subj_eff_y[subject] - 0.4 * within_signal + rnorm(n_subjects * n_rep, sd = 0.4)
  )

  fit_mat <- rmcorr(
    dat,
    response = c("x", "y", "z"),
    subject = "subject",
    keep_data = TRUE
  )
  prepared <- matrixCorr:::`.mc_prepare_rmcorr_inputs`(fit_mat)

  expect_length(prepared, 1L)
  expect_true(inherits(prepared$default$matrix, "matrix"))
  expect_null(colnames(prepared$default$source_data$response))
  expect_equal(prepared$default$source_data$response_names, c("x", "y", "z"))
  expect_equal(prepared$default$source_data$subject_code, as.integer(dat$subject))
  expect_equal(prepared$default$source_data$subject_levels, as.character(seq_len(n_subjects)))
  expect_equal(prepared$default$n_subjects, n_subjects)

  fit_xy_direct <- rmcorr(dat, response = c("x", "y"), subject = "subject")
  fit_xy_view <- matrixCorr:::`.mc_rmcorr_view_pair_fit`(
    prepared$default,
    x_var = "x",
    y_var = "y"
  )

  expect_s3_class(fit_xy_view, "rmcorr")
  expect_equal(fit_xy_view$responses, c("x", "y"))
  expect_equal(fit_xy_view$r, fit_xy_direct$r, tolerance = 1e-12)
  expect_equal(fit_xy_view$slope, fit_xy_direct$slope, tolerance = 1e-12)
  expect_equal(fit_xy_view$p_value, fit_xy_direct$p_value, tolerance = 1e-12)
  expect_equal(fit_xy_view$based.on, fit_xy_direct$based.on)
  expect_equal(fit_xy_view$n_subjects, fit_xy_direct$n_subjects)
})

test_that("rmcorr viewer parser rejects matrix objects without retained data", {
  dat <- data.frame(
    subject = rep(1:3, each = 3),
    x = rnorm(9),
    y = rnorm(9),
    z = rnorm(9)
  )
  bad <- rmcorr(dat, response = c("x", "y", "z"), subject = "subject")

  expect_error(
    matrixCorr:::`.mc_parse_rmcorr_object`(bad, label = "bad"),
    "keep_data = TRUE"
  )
})
