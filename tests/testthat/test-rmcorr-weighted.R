reference_rmcorr_weighted <- function(data, x, y, subject) {
  dat <- data.frame(
    x = data[[x]],
    y = data[[y]],
    subject = data[[subject]]
  )
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  dat$subject <- factor(dat$subject)

  counts <- table(dat$subject)
  keep <- names(counts[counts >= 2L])
  dat <- dat[dat$subject %in% keep, , drop = FALSE]
  dat$subject <- droplevels(dat$subject)

  if (nlevels(dat$subject) < 2L) {
    return(list(
      estimate = NA_real_,
      slope = NA_real_,
      weighted_ss_x = NA_real_,
      weighted_ss_e = NA_real_,
      n_complete = nrow(dat),
      n_subjects = nlevels(dat$subject),
      valid = FALSE
    ))
  }

  df <- nrow(dat) - nlevels(dat$subject) - 1L
  if (df <= 0L) {
    return(list(
      estimate = NA_real_,
      slope = NA_real_,
      weighted_ss_x = NA_real_,
      weighted_ss_e = NA_real_,
      n_complete = nrow(dat),
      n_subjects = nlevels(dat$subject),
      valid = FALSE
    ))
  }

  fit_null <- stats::lm(y ~ subject, data = dat)
  fit_full <- stats::lm(y ~ x + subject, data = dat)

  beta <- unname(stats::coef(fit_full)[["x"]])

  e_null <- stats::residuals(fit_null)
  e_full <- stats::residuals(fit_full)

  dat$qx <- e_null^2 - e_full^2
  dat$qe <- e_full^2

  by_subject <- stats::aggregate(
    cbind(qx, qe) ~ subject,
    data = dat,
    FUN = mean
  )

  weighted_ss_x <- sum(by_subject$qx)
  weighted_ss_e <- sum(by_subject$qe)
  denom <- weighted_ss_x + weighted_ss_e

  if (!is.finite(beta) ||
      !is.finite(weighted_ss_x) ||
      !is.finite(weighted_ss_e) ||
      denom <= 0) {
    est <- NA_real_
    valid <- FALSE
  } else {
    if (weighted_ss_x < 0 && weighted_ss_x > -sqrt(.Machine$double.eps)) {
      weighted_ss_x <- 0
    }
    if (weighted_ss_x < 0) {
      est <- NA_real_
      valid <- FALSE
    } else {
      est <- sign(beta) * sqrt(weighted_ss_x / denom)
      est <- max(-1, min(1, est))
      valid <- is.finite(est)
    }
  }

  list(
    estimate = est,
    slope = beta,
    weighted_ss_x = weighted_ss_x,
    weighted_ss_e = weighted_ss_e,
    n_complete = nrow(dat),
    n_subjects = nlevels(dat$subject),
    df = df,
    valid = valid
  )
}

make_balanced_complete_data <- function() {
  data.frame(
    subject = rep(1:4, each = 4),
    x = rep(c(-1, 0, 1, 2), times = 4) + rep(c(0.0, 0.2, -0.1, 0.15), each = 4),
    y = rep(c(10, 11, 12, 13), each = 4) +
      1.25 * rep(c(-1, 0, 1, 2), times = 4) +
      rep(c(0.10, -0.05, 0.03, -0.02), times = 4),
    z = rep(c(4, 6, 5, 7), each = 4) -
      0.75 * rep(c(-1, 0, 1, 2), times = 4)
  )
}

make_unbalanced_weight_data <- function() {
  data.frame(
    subject = c(
      rep("A", 8),
      rep("B", 2),
      rep("C", 2)
    ),
    x = c(1:8, 1, 2, 1, 2),
    y = c(1:8, 2, 1, 2, 1)
  )
}

test_that("rmcorr default ANCOVA output matches explicit ancova configuration", {
  dat <- make_balanced_complete_data()

  fit_default <- rmcorr(dat, response = c("x", "y"), subject = "subject")
  fit_ancova <- rmcorr(
    dat,
    response = c("x", "y"),
    subject = "subject",
    estimator = "ancova",
    ci_method = "auto"
  )

  expect_identical(unclass(fit_default), unclass(fit_ancova))
  expect_identical(attributes(fit_default), attributes(fit_ancova))

  mat_default <- rmcorr(dat, response = c("x", "y", "z"), subject = "subject", n_threads = 1L)
  mat_ancova <- rmcorr(
    dat,
    response = c("x", "y", "z"),
    subject = "subject",
    estimator = "ancova",
    ci_method = "auto",
    n_threads = 1L
  )
  expect_equal(unclass(mat_default), unclass(mat_ancova), tolerance = 1e-12)
  expect_equal(attr(mat_default, "diagnostics"), attr(mat_ancova, "diagnostics"))
})

test_that("weighted rmcorr matches lm reference implementation exactly", {
  dat <- data.frame(
    subject = c(
      rep("S1", 5),
      rep("S2", 3),
      rep("S3", 2),
      rep("S4", 4)
    ),
    x = c(1, 2, 3, 4, 5, 2, 3, 4, 3, 5, 1, 2, 4, 6),
    y = c(2, 4, 6, 7, 10, 4, 5, 8, 5, 9, 2, 4, 7, 11)
  )

  ref <- reference_rmcorr_weighted(dat, "x", "y", "subject")
  fit <- rmcorr(
    dat,
    response = c("x", "y"),
    subject = "subject",
    estimator = "weighted",
    ci_method = "none"
  )

  expect_equal(fit$estimate, ref$estimate, tolerance = 1e-10)
  expect_equal(fit$slope, ref$slope, tolerance = 1e-10)
  expect_equal(fit$weighted_ss_x, ref$weighted_ss_x, tolerance = 1e-10)
  expect_equal(fit$weighted_ss_e, ref$weighted_ss_e, tolerance = 1e-10)
  expect_equal(fit$n_complete, ref$n_complete)
  expect_equal(fit$n_subjects, ref$n_subjects)
  expect_equal(fit$df, ref$df)
  expect_identical(fit$valid, ref$valid)
})

test_that("weighted equals ancova on balanced complete repeated data", {
  dat <- make_balanced_complete_data()
  fit_ancova <- rmcorr(
    dat,
    response = c("x", "y"),
    subject = "subject",
    estimator = "ancova",
    ci_method = "none"
  )
  fit_weighted <- rmcorr(
    dat,
    response = c("x", "y"),
    subject = "subject",
    estimator = "weighted",
    ci_method = "none"
  )

  expect_equal(fit_weighted$estimate, fit_ancova$estimate, tolerance = 1e-12)
})

test_that("weighted differs from ancova on unbalanced subject counts", {
  dat <- make_unbalanced_weight_data()

  fit_weighted <- rmcorr(
    dat,
    response = c("x", "y"),
    subject = "subject",
    estimator = "weighted",
    ci_method = "none"
  )
  fit_ancova <- rmcorr(
    dat,
    response = c("x", "y"),
    subject = "subject",
    estimator = "ancova",
    ci_method = "none"
  )
  ref <- reference_rmcorr_weighted(dat, "x", "y", "subject")

  expect_equal(fit_weighted$estimate, ref$estimate, tolerance = 1e-10)
  expect_equal(fit_weighted$weighted_ss_x, ref$weighted_ss_x, tolerance = 1e-10)
  expect_equal(fit_weighted$weighted_ss_e, ref$weighted_ss_e, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(fit_weighted$estimate, fit_ancova$estimate, tolerance = 1e-8)))
})

test_that("weighted matrix mode applies pairwise filtering per response pair", {
  dat <- data.frame(
    subject = rep(1:4, each = 3),
    x = c(1, 2, 3, 1, 3, 5, 2, 4, 6, 1, 2, 3),
    y = c(2, 4, NA, 1, 2, 3, 4, 6, 8, 2, NA, 4),
    z = c(NA, 3, 4, 2, 4, 6, 3, 5, NA, 1, 2, 4)
  )

  fit <- rmcorr(
    dat,
    response = c("x", "y", "z"),
    subject = "subject",
    estimator = "weighted",
    na_method = "pairwise",
    ci_method = "none",
    n_threads = 1L
  )
  diag_attr <- attr(fit, "diagnostics")

  ref_xy <- reference_rmcorr_weighted(dat, "x", "y", "subject")
  ref_xz <- reference_rmcorr_weighted(dat, "x", "z", "subject")
  ref_yz <- reference_rmcorr_weighted(dat, "y", "z", "subject")

  expect_equal(unname(fit["x", "y"]), ref_xy$estimate, tolerance = 1e-10)
  expect_equal(unname(fit["x", "z"]), ref_xz$estimate, tolerance = 1e-10)
  expect_equal(unname(fit["y", "z"]), ref_yz$estimate, tolerance = 1e-10)

  expect_equal(unname(diag_attr$n_complete["x", "y"]), ref_xy$n_complete)
  expect_equal(unname(diag_attr$n_complete["x", "z"]), ref_xz$n_complete)
  expect_equal(unname(diag_attr$n_complete["y", "z"]), ref_yz$n_complete)

  expect_equal(unname(diag_attr$n_subjects["x", "y"]), ref_xy$n_subjects)
  expect_equal(unname(diag_attr$n_subjects["x", "z"]), ref_xz$n_subjects)
  expect_equal(unname(diag_attr$n_subjects["y", "z"]), ref_yz$n_subjects)

  pair_n <- c(
    unname(diag_attr$n_complete["x", "y"]),
    unname(diag_attr$n_complete["x", "z"]),
    unname(diag_attr$n_complete["y", "z"])
  )
  expect_gt(length(unique(pair_n)), 1L)
})

test_that("weighted estimator respects na_method error and pairwise", {
  dat <- data.frame(
    subject = c(1, 1, 2, 2, 3, 3),
    x = c(1, 2, 1, NA, 2, 3),
    y = c(2, 4, 2, 3, 3, 5)
  )

  expect_error(
    rmcorr(dat, response = c("x", "y"), subject = "subject", estimator = "weighted", na_method = "error"),
    "missing, NaN, or infinite values"
  )
  expect_no_error(
    rmcorr(dat, response = c("x", "y"), subject = "subject", estimator = "weighted", na_method = "pairwise")
  )
})

test_that("subjects with fewer than two complete pairs are dropped", {
  dat <- data.frame(
    subject = c(1, 1, 1, 2, 2, 3, 3),
    x = c(1, 2, 3, 1, 2, 1, 2),
    y = c(2, 4, 6, 3, 6, 5, NA)
  )
  fit <- rmcorr(
    dat,
    response = c("x", "y"),
    subject = "subject",
    estimator = "weighted",
    na_method = "pairwise",
    ci_method = "none"
  )

  expect_equal(fit$n_subjects, 2L)
  expect_equal(fit$n_complete, 5L)
  expect_equal(fit$obs_per_subject_min, 2L)
  expect_equal(fit$obs_per_subject_max, 3L)
})

test_that("weighted estimator preserves negative association sign", {
  dat <- data.frame(
    subject = rep(1:3, each = 4),
    x = rep(c(1, 2, 3, 4), times = 3),
    y = rep(c(10, 12, 14), each = 4) - 1.5 * rep(c(1, 2, 3, 4), times = 3)
  )
  ref <- reference_rmcorr_weighted(dat, "x", "y", "subject")
  fit <- rmcorr(
    dat,
    response = c("x", "y"),
    subject = "subject",
    estimator = "weighted",
    ci_method = "none"
  )

  expect_lt(fit$slope, 0)
  expect_lt(fit$estimate, 0)
  expect_equal(abs(fit$estimate), abs(ref$estimate), tolerance = 1e-10)
})

test_that("weighted estimator handles zero within-subject x variance", {
  dat <- data.frame(
    subject = c(1, 1, 2, 2, 3, 3),
    x = c(1, 1, 2, 2, 3, 3),
    y = c(1, 2, 3, 4, 5, 6)
  )
  fit <- rmcorr(
    dat,
    response = c("x", "y"),
    subject = "subject",
    estimator = "weighted",
    ci_method = "none"
  )

  expect_true(is.na(fit$estimate))
  expect_true(is.na(fit$slope))
  expect_false(fit$valid)
})

test_that("weighted estimator handles zero within-subject y variance", {
  dat <- data.frame(
    subject = c(1, 1, 2, 2, 3, 3),
    x = c(1, 2, 1, 2, 1, 2),
    y = c(5, 5, 7, 7, 9, 9)
  )
  fit <- rmcorr(
    dat,
    response = c("x", "y"),
    subject = "subject",
    estimator = "weighted",
    ci_method = "none"
  )

  expect_true(is.na(fit$estimate))
  expect_false(fit$valid)
})

test_that("weighted bootstrap CI is reproducible with fixed seed", {
  dat <- make_balanced_complete_data()

  fit1 <- rmcorr(
    dat,
    response = c("x", "y"),
    subject = "subject",
    estimator = "weighted",
    ci_method = "bootstrap",
    n_boot = 199L,
    seed = 2026L
  )
  fit2 <- rmcorr(
    dat,
    response = c("x", "y"),
    subject = "subject",
    estimator = "weighted",
    ci_method = "bootstrap",
    n_boot = 199L,
    seed = 2026L
  )

  expect_equal(c(fit1$lwr, fit1$upr), c(fit2$lwr, fit2$upr))
  expect_true(all(is.finite(c(fit1$lwr, fit1$upr))))
  expect_identical(fit1$ci_method, "bootstrap")
})

test_that("weighted bootstrap resamples subjects as blocks with reassigned IDs", {
  dat <- data.frame(
    .response1 = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    .response2 = c(2, 3, 2, 5, 7, 6, 8, 9, 11),
    .subject = factor(c("A", "A", "B", "B", "B", "C", "C", "C", "C"))
  )

  sampled <- matrixCorr:::`.mc_rmcorr_resample_subject_blocks`(
    dat,
    draw = c(2L, 2L, 1L),
    keep_source = TRUE
  )$data

  expect_equal(as.integer(table(sampled$.subject)), c(3L, 3L, 2L))
  expect_identical(as.character(unique(sampled$.subject_source[sampled$.subject == 1L])), "B")
  expect_identical(as.character(unique(sampled$.subject_source[sampled$.subject == 2L])), "B")
  expect_identical(as.character(unique(sampled$.subject_source[sampled$.subject == 3L])), "A")

  orig_A <- dat[dat$.subject == "A", c(".response1", ".response2")]
  orig_B <- dat[dat$.subject == "B", c(".response1", ".response2")]
  new_1 <- sampled[sampled$.subject == 1L, c(".response1", ".response2")]
  new_2 <- sampled[sampled$.subject == 2L, c(".response1", ".response2")]
  new_3 <- sampled[sampled$.subject == 3L, c(".response1", ".response2")]

  expect_equal(new_1, orig_B)
  expect_equal(new_2, orig_B)
  expect_equal(new_3, orig_A)
})

test_that("weighted matrix output is symmetric with symmetric diagnostics", {
  dat <- data.frame(
    subject = rep(1:5, each = 3),
    x = c(1, 2, 3, 2, 4, 6, 1, 3, 5, 2, 3, 4, 4, 5, 6),
    y = c(2, 3, 5, 3, 5, 7, 1, 2, 4, 3, 4, 5, 5, 6, 8),
    z = c(5, 4, 3, 4, 3, 2, 6, 5, 4, 3, 2, 1, 7, 6, 5)
  )
  fit <- rmcorr(
    dat,
    response = c("x", "y", "z"),
    subject = "subject",
    estimator = "weighted",
    ci_method = "none",
    n_threads = 1L
  )
  diag_attr <- attr(fit, "diagnostics")

  expect_equal(unclass(fit), t(unclass(fit)), tolerance = 1e-12)
  expect_equal(rownames(fit), c("x", "y", "z"))
  expect_equal(colnames(fit), c("x", "y", "z"))

  matrix_fields <- c(
    "slope", "p_value", "df", "n_complete", "n_subjects",
    "conf_low", "conf_high", "weighted_ss_x", "weighted_ss_e",
    "obs_per_subject_min", "obs_per_subject_max", "valid"
  )
  for (nm in matrix_fields) {
    expect_equal(diag_attr[[nm]], t(diag_attr[[nm]]), tolerance = 1e-12)
  }

  valid_diag <- diag(diag_attr$valid)
  est_diag <- diag(as.matrix(fit))
  idx <- which(!is.na(valid_diag) & valid_diag)
  if (length(idx)) {
    expect_equal(est_diag[idx], rep(1, length(idx)), tolerance = 1e-12)
  }
})

test_that("weighted outputs remain compatible with rmcorr/rmcorr_matrix methods", {
  dat <- make_balanced_complete_data()

  fit_pair <- rmcorr(
    dat,
    response = c("x", "y"),
    subject = "subject",
    estimator = "weighted",
    ci_method = "none",
    keep_data = TRUE
  )
  expect_s3_class(fit_pair, "rmcorr")
  expect_equal(fit_pair$r, fit_pair$estimate, tolerance = 1e-12)
  expect_equal(fit_pair$conf_int, c(lower = fit_pair$lwr, upper = fit_pair$upr), tolerance = 1e-12)
  expect_no_error(print(fit_pair))
  expect_no_error(summary(fit_pair))
  expect_s3_class(plot(fit_pair), "ggplot")

  fit_matrix <- rmcorr(
    dat,
    response = c("x", "y", "z"),
    subject = "subject",
    estimator = "weighted",
    ci_method = "none",
    n_threads = 1L
  )
  expect_s3_class(fit_matrix, "rmcorr_matrix")
  expect_no_error(summary(fit_matrix))
  expect_s3_class(plot(fit_matrix), "ggplot")
})

test_that("weighted matrix results are thread-consistent", {
  dat <- data.frame(
    subject = rep(1:10, each = 4),
    x = rep(1:4, times = 10) + rep(seq_len(10) / 10, each = 4),
    y = rep(c(2, 3, 5, 6), times = 10) + rep(seq_len(10) / 20, each = 4),
    z = rep(c(6, 5, 4, 3), times = 10) - rep(seq_len(10) / 25, each = 4),
    w = rep(c(1, 4, 2, 5), times = 10)
  )
  dat$y[c(3, 17, 29)] <- NA_real_
  dat$z[c(6, 14, 31)] <- NA_real_

  fit1 <- rmcorr(
    dat,
    response = c("x", "y", "z", "w"),
    subject = "subject",
    estimator = "weighted",
    na_method = "pairwise",
    ci_method = "none",
    n_threads = 1L
  )
  fit2 <- rmcorr(
    dat,
    response = c("x", "y", "z", "w"),
    subject = "subject",
    estimator = "weighted",
    na_method = "pairwise",
    ci_method = "none",
    n_threads = 2L
  )

  expect_equal(unclass(fit1), unclass(fit2), tolerance = 1e-12)

  d1 <- attr(fit1, "diagnostics")
  d2 <- attr(fit2, "diagnostics")
  expect_equal(names(d1), names(d2))
  for (nm in names(d1)) {
    if (is.matrix(d1[[nm]]) || is.array(d1[[nm]])) {
      expect_equal(d1[[nm]], d2[[nm]], tolerance = 1e-12)
    } else {
      expect_equal(d1[[nm]], d2[[nm]])
    }
  }
})
