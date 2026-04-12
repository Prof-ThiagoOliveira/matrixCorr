biserial_manual_test <- function(x, y) {
  y <- as.integer(as.logical(y))
  keep <- stats::complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]

  p <- mean(y == 1L)
  q <- 1 - p
  z_p <- stats::qnorm(p)
  phi_z <- stats::dnorm(z_p)
  x1 <- mean(x[y == 1L])
  x0 <- mean(x[y == 0L])
  sx <- stats::sd(x)

  ((x1 - x0) / sx) * (p * q / phi_z)
}

# Brent optimization in the C++ latent estimators currently targets 1e-8,
# so allow for small platform/compiler drift around these regression values.
latent_reg_tol <- 1e-7

expect_equal_abs <- function(object, expected, tolerance = latent_reg_tol) {
  expect_lte(
    max(abs(as.numeric(object) - as.numeric(expected))),
    tolerance
  )
}

polyserial_reference_old <- function(x, y, check_na = TRUE) {
  if (!check_na) {
    keep <- stats::complete.cases(x, y)
    x <- x[keep]
    y <- y[keep]
  }

  if (length(x) < 2L || length(unique(y[!is.na(y)])) < 2L) {
    return(NA_real_)
  }

  x <- as.numeric(x)
  y <- as.integer(as.factor(y))
  z <- as.numeric(scale(x))
  tab <- table(y)
  n <- sum(tab)
  s <- length(tab)
  if (s < 2L || sum(tab != 0L) < 2L) {
    return(NA_real_)
  }

  cuts <- stats::qnorm(cumsum(tab) / n)[-s]
  rho <- sqrt((n - 1) / n) * stats::sd(y) * stats::cor(x, y) /
    sum(stats::dnorm(cuts))
  maxcor <- 0.9999

  if (abs(rho) > maxcor) {
    rho <- sign(rho) * maxcor
  }

  out <- suppressWarnings(
    tryCatch(
      stats::optim(
        c(rho, cuts),
        function(pars) matrixCorr:::matrixCorr_polyserial_negloglik_cpp(z, y, pars, maxcor = maxcor)
      )$par[1L],
      error = function(e) NA_real_
    )
  )

  if (is.na(out)) {
    return(NA_real_)
  }
  if (out > 1) {
    out <- maxcor
  } else if (out < -1) {
    out <- -maxcor
  }
  as.numeric(out)
}

latent_estimate_template <- function(fit) {
  if (is.null(dim(fit))) {
    return(as.numeric(fit))
  }

  est_mat <- as.matrix(fit)
  attributes(est_mat) <- attributes(est_mat)[c("dim", "dimnames")]
  est_mat
}

expect_same_shape <- function(object, template) {
  if (is.null(dim(template))) {
    expect_length(object, length(template))
    return(invisible(NULL))
  }

  expect_equal(dim(object), dim(template))
  expect_identical(dimnames(object), dimnames(template))
  invisible(NULL)
}

expect_latent_wald_payloads <- function(fit,
                                        ci_method,
                                        inference_method) {
  est <- latent_estimate_template(fit)
  inf <- attr(fit, "inference", exact = TRUE)
  ci <- attr(fit, "ci", exact = TRUE)

  if (is.list(inf)) {
    expect_identical(inf$method, inference_method)
    expect_identical(inf$alternative, "two.sided")
    expect_same_shape(inf$estimate, est)
    expect_same_shape(inf$statistic, est)
    expect_same_shape(inf$parameter, est)
    expect_same_shape(inf$p_value, est)
    expect_same_shape(inf$n_obs, est)
    expect_equal(as.numeric(inf$estimate), as.numeric(est), tolerance = 1e-12)
    expect_true(all(is.na(as.numeric(inf$parameter))))
    expect_true(all(is.na(as.numeric(inf$statistic)) | is.finite(as.numeric(inf$statistic))))
    expect_true(all(is.na(as.numeric(inf$p_value)) |
      (as.numeric(inf$p_value) >= 0 & as.numeric(inf$p_value) <= 1)))
    expect_true(all(is.na(as.numeric(inf$n_obs)) | as.numeric(inf$n_obs) >= 0))
  }

  if (is.list(ci)) {
    expect_identical(ci$ci.method, ci_method)
    expect_same_shape(ci$est, est)
    expect_same_shape(ci$lwr.ci, est)
    expect_same_shape(ci$upr.ci, est)
    expect_equal(as.numeric(ci$est), as.numeric(est), tolerance = 1e-12)
    expect_equal(attr(fit, "conf.level", exact = TRUE), ci$conf.level, tolerance = 1e-12)
    expect_true(all(is.na(as.numeric(ci$lwr.ci)) |
      (as.numeric(ci$lwr.ci) >= -1 & as.numeric(ci$lwr.ci) <= 1)))
    expect_true(all(is.na(as.numeric(ci$upr.ci)) |
      (as.numeric(ci$upr.ci) >= -1 & as.numeric(ci$upr.ci) <= 1)))
    expect_true(all(
      is.na(as.numeric(ci$lwr.ci)) |
      is.na(as.numeric(ci$upr.ci)) |
      as.numeric(ci$lwr.ci) <= as.numeric(ci$upr.ci)
    ))
    keep <- is.finite(as.numeric(est)) &
      is.finite(as.numeric(ci$lwr.ci)) &
      is.finite(as.numeric(ci$upr.ci))
    expect_true(all(
      as.numeric(est)[keep] >= as.numeric(ci$lwr.ci)[keep] &
        as.numeric(est)[keep] <= as.numeric(ci$upr.ci)[keep]
    ))
  }
}

test_that("tetrachoric pair and matrix modes match fixed regression values", {
  set.seed(1001)
  x <- sample(c(FALSE, TRUE), 300, replace = TRUE)
  y <- sample(c(FALSE, TRUE), 300, replace = TRUE)

  est_pair <- tetrachoric(x, y)
  expect_equal_abs(as.numeric(est_pair), -0.272117497540896)
  expect_true(is.list(attr(est_pair, "diagnostics")))
  expect_true(is.list(attr(est_pair, "thresholds")))

  dat <- data.frame(
    a = x,
    b = y,
    c = sample(c(FALSE, TRUE), 300, replace = TRUE)
  )
  est_mat <- tetrachoric(dat)
  ref_mat <- matrix(
    c(
      1, -0.272117497540896, 0.0624107955993122,
      -0.272117497540896, 1, 0.0348797202693802,
      0.0624107955993122, 0.0348797202693802, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(names(dat), names(dat))
  )

  expect_s3_class(est_mat, "tetrachoric_corr")
  expect_equal_abs(as.numeric(est_mat), as.numeric(ref_mat))
  expect_equal_abs(as.numeric(est_mat["a", "b"]), as.numeric(est_pair))
})

test_that("polychoric pair and matrix modes match fixed regression values", {
  set.seed(1002)
  x <- ordered(sample(1:4, 400, replace = TRUE))
  y <- ordered(sample(1:5, 400, replace = TRUE))

  est_pair <- polychoric(x, y)
  expect_equal_abs(as.numeric(est_pair), 0.0034270943899109)

  dat <- data.frame(
    x = x,
    y = y,
    z = ordered(sample(1:3, 400, replace = TRUE))
  )
  est_mat <- polychoric(dat)
  ref_mat <- matrix(
    c(
      1, 0.0034270943899109, 0.00263097305519399,
      0.0034270943899109, 1, 0.0231228372629091,
      0.00263097305519399, 0.0231228372629091, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(names(dat), names(dat))
  )

  expect_s3_class(est_mat, "polychoric_corr")
  expect_equal_abs(as.numeric(est_mat), as.numeric(ref_mat))
  expect_equal_abs(as.numeric(est_mat["x", "y"]), as.numeric(est_pair))
})

test_that("polyserial pair and matrix modes match fixed regression values", {
  set.seed(1003)
  x <- rnorm(350)
  y <- ordered(sample(1:4, 350, replace = TRUE))

  est_pair <- polyserial(x, y)
  expect_equal(as.numeric(est_pair), 0.0418913522648746, tolerance = 1e-12)

  X <- data.frame(x1 = rnorm(350), x2 = rnorm(350))
  Y <- data.frame(
    y1 = y,
    y2 = ordered(sample(1:5, 350, replace = TRUE))
  )
  est_mat <- polyserial(X, Y)
  expect_s3_class(est_mat, "polyserial_corr")
  expect_equal(
    as.numeric(est_mat),
    c(
      -0.0127484729795922, -0.0159091216414217,
      0.0706339051579041, -0.0149611262803363
    ),
    tolerance = 1e-12
  )
})

test_that("polyserial matches the previous R optim reference path", {
  set.seed(2468)
  x <- rnorm(240)
  y <- ordered(sample(c(1, 2, 4), 240, replace = TRUE))

  expect_equal(
    polyserial(x, y),
    polyserial_reference_old(x, y),
    tolerance = 1e-12
  )

  x[c(5, 19)] <- NA_real_
  y[c(11, 19)] <- NA
  expect_equal(
    polyserial(x, y, na_method = "pairwise"),
    polyserial_reference_old(x, y, check_na = FALSE),
    tolerance = 1e-12
  )

  X <- data.frame(x1 = rnorm(240), x2 = rnorm(240))
  Y <- data.frame(
    y1 = ordered(sample(1:4, 240, replace = TRUE)),
    y2 = ordered(sample(c(1, 3, 4), 240, replace = TRUE))
  )
  ref_mat <- outer(
    seq_len(ncol(X)),
    seq_len(ncol(Y)),
    Vectorize(function(j, k) polyserial_reference_old(X[[j]], Y[[k]]))
  )
  dimnames(ref_mat) <- list(names(X), names(Y))

  expect_equal(
    as.numeric(polyserial(X, Y)),
    as.numeric(ref_mat),
    tolerance = 1e-12
  )
})

test_that("tetrachoric Wald inference preserves the established estimates", {
  set.seed(1101)
  x <- sample(c(FALSE, TRUE), 320, replace = TRUE)
  y <- sample(c(FALSE, TRUE), 320, replace = TRUE)
  dat <- data.frame(
    a = x,
    b = y,
    c = sample(c(FALSE, TRUE), 320, replace = TRUE)
  )

  fit_pair <- tetrachoric(x, y, ci = TRUE, p_value = TRUE)
  fit_mat <- tetrachoric(dat, ci = TRUE, p_value = TRUE)

  expect_equal_abs(as.numeric(fit_pair), as.numeric(tetrachoric(x, y)))
  expect_equal_abs(as.numeric(fit_mat), as.numeric(tetrachoric(dat)))
  expect_latent_wald_payloads(
    fit_pair,
    ci_method = "wald_information_tetrachoric",
    inference_method = "wald_z_tetrachoric"
  )
  expect_latent_wald_payloads(
    fit_mat,
    ci_method = "wald_information_tetrachoric",
    inference_method = "wald_z_tetrachoric"
  )
})

test_that("polychoric Wald inference preserves the established estimates", {
  set.seed(1102)
  x <- ordered(sample(1:4, 360, replace = TRUE))
  y <- ordered(sample(1:5, 360, replace = TRUE))
  dat <- data.frame(
    x = x,
    y = y,
    z = ordered(sample(1:3, 360, replace = TRUE))
  )

  fit_pair <- polychoric(x, y, ci = TRUE, p_value = TRUE)
  fit_mat <- polychoric(dat, ci = TRUE, p_value = TRUE)

  expect_equal_abs(as.numeric(fit_pair), as.numeric(polychoric(x, y)))
  expect_equal_abs(as.numeric(fit_mat), as.numeric(polychoric(dat)))
  expect_latent_wald_payloads(
    fit_pair,
    ci_method = "wald_information_polychoric",
    inference_method = "wald_z_polychoric"
  )
  expect_latent_wald_payloads(
    fit_mat,
    ci_method = "wald_information_polychoric",
    inference_method = "wald_z_polychoric"
  )
})

test_that("polyserial Wald inference preserves the established estimates", {
  set.seed(1103)
  X <- data.frame(x1 = rnorm(280), x2 = rnorm(280))
  Y <- data.frame(
    y1 = ordered(sample(1:4, 280, replace = TRUE)),
    y2 = ordered(sample(c(1, 3, 4), 280, replace = TRUE))
  )

  fit_pair <- polyserial(X$x1, Y$y1, ci = TRUE, p_value = TRUE)
  fit_mat <- polyserial(X, Y, ci = TRUE, p_value = TRUE)

  expect_equal_abs(as.numeric(fit_pair), as.numeric(polyserial(X$x1, Y$y1)))
  expect_equal_abs(as.numeric(fit_mat), as.numeric(polyserial(X, Y)))
  expect_latent_wald_payloads(
    fit_pair,
    ci_method = "wald_information_polyserial",
    inference_method = "wald_z_polyserial"
  )
  expect_latent_wald_payloads(
    fit_mat,
    ci_method = "wald_information_polyserial",
    inference_method = "wald_z_polyserial"
  )
})

test_that("latent inference payloads are attached only when requested", {
  set.seed(1104)
  x_bin <- sample(c(FALSE, TRUE), 180, replace = TRUE)
  y_bin <- sample(c(FALSE, TRUE), 180, replace = TRUE)
  x_ord <- ordered(sample(1:4, 180, replace = TRUE))
  y_ord <- ordered(sample(1:5, 180, replace = TRUE))
  X <- data.frame(x1 = rnorm(180), x2 = rnorm(180))
  Y <- data.frame(y1 = x_ord, y2 = y_ord)

  tetra_0 <- tetrachoric(data.frame(a = x_bin, b = y_bin), ci = FALSE, p_value = FALSE)
  tetra_ci <- tetrachoric(data.frame(a = x_bin, b = y_bin), ci = TRUE, p_value = FALSE)
  tetra_p <- tetrachoric(data.frame(a = x_bin, b = y_bin), ci = FALSE, p_value = TRUE)
  tetra_both <- tetrachoric(data.frame(a = x_bin, b = y_bin), ci = TRUE, p_value = TRUE)

  poly_0 <- polychoric(data.frame(x = x_ord, y = y_ord), ci = FALSE, p_value = FALSE)
  poly_ci <- polychoric(data.frame(x = x_ord, y = y_ord), ci = TRUE, p_value = FALSE)
  poly_p <- polychoric(data.frame(x = x_ord, y = y_ord), ci = FALSE, p_value = TRUE)
  poly_both <- polychoric(data.frame(x = x_ord, y = y_ord), ci = TRUE, p_value = TRUE)

  ps_0 <- polyserial(X, Y, ci = FALSE, p_value = FALSE)
  ps_ci <- polyserial(X, Y, ci = TRUE, p_value = FALSE)
  ps_p <- polyserial(X, Y, ci = FALSE, p_value = TRUE)
  ps_both <- polyserial(X, Y, ci = TRUE, p_value = TRUE)

  for (fit in list(tetra_0, poly_0, ps_0)) {
    expect_null(attr(fit, "ci", exact = TRUE))
    expect_null(attr(fit, "inference", exact = TRUE))
  }

  for (fit in list(tetra_ci, poly_ci, ps_ci)) {
    expect_true(is.list(attr(fit, "ci", exact = TRUE)))
    expect_null(attr(fit, "inference", exact = TRUE))
  }

  for (fit in list(tetra_p, poly_p, ps_p)) {
    expect_true(is.list(attr(fit, "inference", exact = TRUE)))
    expect_null(attr(fit, "ci", exact = TRUE))
  }

  for (fit in list(tetra_both, poly_both, ps_both)) {
    expect_true(is.list(attr(fit, "ci", exact = TRUE)))
    expect_true(is.list(attr(fit, "inference", exact = TRUE)))
  }
})

test_that("latent core estimates are unchanged when CI and p-values are disabled", {
  set.seed(1105)
  dat_bin <- data.frame(
    a = sample(c(FALSE, TRUE), 160, replace = TRUE),
    b = sample(c(FALSE, TRUE), 160, replace = TRUE),
    c = sample(c(FALSE, TRUE), 160, replace = TRUE)
  )
  dat_ord <- data.frame(
    x = ordered(sample(1:4, 160, replace = TRUE)),
    y = ordered(sample(1:5, 160, replace = TRUE)),
    z = ordered(sample(1:3, 160, replace = TRUE))
  )
  X <- data.frame(x1 = rnorm(160), x2 = rnorm(160))
  Y <- data.frame(
    y1 = ordered(sample(1:4, 160, replace = TRUE)),
    y2 = ordered(sample(c(1, 3, 4), 160, replace = TRUE))
  )

  fit_t_base <- tetrachoric(dat_bin)
  fit_t_false <- tetrachoric(dat_bin, ci = FALSE, p_value = FALSE)
  fit_p_base <- polychoric(dat_ord)
  fit_p_false <- polychoric(dat_ord, ci = FALSE, p_value = FALSE)
  fit_ps_base <- polyserial(X, Y)
  fit_ps_false <- polyserial(X, Y, ci = FALSE, p_value = FALSE)

  expect_equal(unclass(fit_t_base), unclass(fit_t_false), tolerance = 1e-12)
  expect_equal(unclass(fit_p_base), unclass(fit_p_false), tolerance = 1e-12)
  expect_equal(unclass(fit_ps_base), unclass(fit_ps_false), tolerance = 1e-12)
  expect_identical(attributes(fit_t_base), attributes(fit_t_false))
  expect_identical(attributes(fit_p_base), attributes(fit_p_false))
  expect_identical(attributes(fit_ps_base), attributes(fit_ps_false))
})

test_that("latent inference summaries switch to the pairwise view only when requested", {
  set.seed(1106)
  dat_bin <- data.frame(
    a = sample(c(FALSE, TRUE), 150, replace = TRUE),
    b = sample(c(FALSE, TRUE), 150, replace = TRUE),
    c = sample(c(FALSE, TRUE), 150, replace = TRUE)
  )
  dat_ord <- data.frame(
    x = ordered(sample(1:4, 150, replace = TRUE)),
    y = ordered(sample(1:5, 150, replace = TRUE)),
    z = ordered(sample(1:3, 150, replace = TRUE))
  )
  X <- data.frame(x1 = rnorm(150), x2 = rnorm(150))
  Y <- data.frame(
    y1 = ordered(sample(1:4, 150, replace = TRUE)),
    y2 = ordered(sample(c(1, 3, 4), 150, replace = TRUE))
  )

  fits0 <- list(
    tetrachoric(dat_bin),
    polychoric(dat_ord),
    polyserial(X, Y)
  )
  fits1 <- list(
    tetrachoric(dat_bin, ci = TRUE, p_value = TRUE),
    polychoric(dat_ord, ci = TRUE, p_value = TRUE),
    polyserial(X, Y, ci = TRUE, p_value = TRUE)
  )
  summary_classes <- c("summary.tetrachoric_corr", "summary.polychoric_corr", "summary.polyserial_corr")
  summary_titles <- c("Tetrachoric correlation summary", "Polychoric correlation summary", "Polyserial correlation summary")

  for (fit in fits0) {
    expect_s3_class(summary(fit), "summary.latent_corr")
    expect_s3_class(summary(fit), "summary.matrixCorr")
  }

  for (i in seq_along(fits1)) {
    sm <- summary(fits1[[i]])
    expect_s3_class(sm, summary_classes[[i]])
    expect_true(all(c("item1", "item2", "estimate", "lwr", "upr", "statistic", "df", "p_value", "n_complete") %in% names(sm)))
    out <- capture.output(print(sm, digits = 3))
    expect_true(any(grepl(summary_titles[[i]], out, fixed = TRUE)))
  }
})

test_that("biserial matches the closed-form estimator in pair and matrix mode", {
  set.seed(1004)
  x <- rnorm(320)
  y <- sample(c(FALSE, TRUE), 320, replace = TRUE)

  est_pair <- biserial(x, y)
  expect_equal(est_pair, biserial_manual_test(x, y), tolerance = 1e-12)

  X <- data.frame(x1 = rnorm(320), x2 = rnorm(320))
  Y <- data.frame(
    b1 = y,
    b2 = sample(c(FALSE, TRUE), 320, replace = TRUE)
  )
  est_mat <- biserial(X, Y)
  ref_mat <- outer(
    seq_len(ncol(X)),
    seq_len(ncol(Y)),
    Vectorize(function(j, k) biserial_manual_test(X[[j]], Y[[k]]))
  )
  dimnames(ref_mat) <- list(names(X), names(Y))

  expect_s3_class(est_mat, "biserial_corr")
  expect_equal(as.numeric(est_mat), as.numeric(ref_mat), tolerance = 1e-12)
})

test_that("biserial large-sample inference follows the established reference behaviour", {
  set.seed(160)
  n <- 180
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  X <- data.frame(
    bx1 = 0.7 * z1 + rnorm(n, sd = 0.6),
    bx2 = -0.5 * z2 + 0.2 * z1 + rnorm(n, sd = 0.7)
  )
  Y <- data.frame(
    g1 = z1 > stats::qnorm(0.58),
    g2 = (0.4 * z1 + 0.8 * z2 + rnorm(n, sd = 0.5)) > stats::qnorm(0.52)
  )

  fit <- biserial(X, Y, ci = TRUE, p_value = TRUE)
  inf <- attr(fit, "inference", exact = TRUE)
  ci <- attr(fit, "ci", exact = TRUE)
  est_mat <- fit
  attributes(est_mat) <- attributes(est_mat)[c("dim", "dimnames")]

  ref_estimate <- structure(
    c(0.792193280766738, 0.103575171678435, 0.260351318125141, -0.330748623730774),
    dim = c(2L, 2L),
    dimnames = list(names(X), names(Y))
  )
  ref_statistic <- structure(
    c(17.3188455939721, 1.38933751928313, 3.59758636169522, -4.67590196928005),
    dim = c(2L, 2L),
    dimnames = list(names(X), names(Y))
  )
  ref_p_value <- structure(
    c(4.99830673160238e-40, 0.166465648114298, 0.000416008709530862, 5.76492513948137e-06),
    dim = c(2L, 2L),
    dimnames = list(names(X), names(Y))
  )
  ref_lwr <- structure(
    c(0.730581588076149, -0.0433447862124336, 0.118604440134517, -0.455000651378134),
    dim = c(2L, 2L),
    dimnames = list(names(X), names(Y))
  )
  ref_upr <- structure(
    c(0.841009788517119, 0.246110105580895, 0.391698804407742, -0.193863734545158),
    dim = c(2L, 2L),
    dimnames = list(names(X), names(Y))
  )
  ref_n_obs <- structure(rep(180L, 4L), dim = c(2L, 2L), dimnames = list(names(X), names(Y)))

  expect_equal(est_mat, ref_estimate, tolerance = 1e-12)
  expect_equal(inf$statistic, ref_statistic, tolerance = 1e-12)
  expect_equal(inf$p_value, ref_p_value, tolerance = 1e-12)
  expect_equal(ci$lwr.ci, ref_lwr, tolerance = 1e-12)
  expect_equal(ci$upr.ci, ref_upr, tolerance = 1e-12)
  expect_equal(inf$n_obs, ref_n_obs, tolerance = 0)
  expect_equal(attr(fit, "diagnostics", exact = TRUE)$n_complete, ref_n_obs, tolerance = 0)
})

test_that("biserial returns CI and p-value payloads only when requested", {
  set.seed(1007)
  x <- rnorm(150)
  y <- x + rnorm(150, sd = 0.7)
  g <- y > stats::qnorm(0.55)
  X <- data.frame(x1 = x, x2 = rnorm(150))
  Y <- data.frame(g1 = g, g2 = sample(c(FALSE, TRUE), 150, replace = TRUE))

  fit0 <- biserial(X, Y, ci = FALSE, p_value = FALSE)
  fit_ci <- biserial(X, Y, ci = TRUE, p_value = FALSE)
  fit_p <- biserial(X, Y, ci = FALSE, p_value = TRUE)
  fit_both <- biserial(X, Y, ci = TRUE, p_value = TRUE)

  expect_null(attr(fit0, "ci", exact = TRUE))
  expect_null(attr(fit0, "inference", exact = TRUE))

  expect_true(is.list(attr(fit_ci, "ci", exact = TRUE)))
  expect_null(attr(fit_ci, "inference", exact = TRUE))

  expect_true(is.list(attr(fit_p, "inference", exact = TRUE)))
  expect_null(attr(fit_p, "ci", exact = TRUE))

  expect_true(is.list(attr(fit_both, "ci", exact = TRUE)))
  expect_true(is.list(attr(fit_both, "inference", exact = TRUE)))
})

test_that("biserial confidence intervals are ordered and bounded", {
  set.seed(1008)
  z <- rnorm(170)
  X <- data.frame(
    x1 = 0.8 * z + rnorm(170, sd = 0.5),
    x2 = rnorm(170)
  )
  Y <- data.frame(
    g1 = z > stats::qnorm(0.6),
    g2 = (0.5 * z + rnorm(170, sd = 0.8)) > stats::qnorm(0.45)
  )

  fit <- biserial(X, Y, ci = TRUE)
  ci <- attr(fit, "ci", exact = TRUE)

  expect_true(all(ci$lwr.ci <= ci$upr.ci, na.rm = TRUE))
  expect_true(all(ci$lwr.ci >= -1, na.rm = TRUE))
  expect_true(all(ci$upr.ci <= 1, na.rm = TRUE))
  expect_true(all(unclass(fit) >= ci$lwr.ci & unclass(fit) <= ci$upr.ci, na.rm = TRUE))
})

test_that("biserial core estimate is unchanged when CI and p-values are disabled", {
  set.seed(1009)
  X <- data.frame(
    x1 = rnorm(140),
    x2 = rnorm(140)
  )
  Y <- data.frame(
    g1 = sample(c(FALSE, TRUE), 140, replace = TRUE),
    g2 = sample(c(FALSE, TRUE), 140, replace = TRUE)
  )

  fit_base <- biserial(X, Y)
  fit_flags_false <- biserial(X, Y, ci = FALSE, p_value = FALSE)

  expect_equal(unclass(fit_base), unclass(fit_flags_false), tolerance = 1e-12)
  expect_equal(attributes(fit_base)[c("dim", "dimnames", "class", "method", "description", "package", "diagnostics", "thresholds", "correct")],
               attributes(fit_flags_false)[c("dim", "dimnames", "class", "method", "description", "package", "diagnostics", "thresholds", "correct")])
})

test_that("biserial summary switches to pairwise inference view when requested", {
  set.seed(1010)
  z <- rnorm(160)
  X <- data.frame(
    x1 = 0.6 * z + rnorm(160, sd = 0.6),
    x2 = rnorm(160)
  )
  Y <- data.frame(
    g1 = z > stats::qnorm(0.57),
    g2 = sample(c(FALSE, TRUE), 160, replace = TRUE)
  )

  fit0 <- biserial(X, Y)
  fit1 <- biserial(X, Y, ci = TRUE, p_value = TRUE)

  sm0 <- summary(fit0)
  sm1 <- summary(fit1)

  expect_s3_class(sm0, "summary.latent_corr")
  expect_s3_class(sm0, "summary.matrixCorr")
  expect_s3_class(sm1, "summary.biserial_corr")
  expect_true(all(c("item1", "item2", "estimate", "lwr", "upr", "statistic", "df", "p_value", "n_complete") %in% names(sm1)))

  out <- capture.output(print(sm1, digits = 3))
  expect_true(any(grepl("Biserial correlation summary", out, fixed = TRUE)))
})

test_that("latent correlation methods support pairwise complete cases", {
  x_bin <- c(FALSE, TRUE, NA, FALSE, TRUE, FALSE)
  y_bin <- c(TRUE, FALSE, TRUE, NA, FALSE, TRUE)
  expect_false(is.na(tetrachoric(x_bin, y_bin, na_method = "pairwise")))

  x_cont <- c(0.1, -0.4, NA, 0.7, 1.2, -0.5)
  y_ord <- ordered(c(1, 2, 3, NA, 2, 1))
  expect_false(is.na(polyserial(x_cont, y_ord, na_method = "pairwise")))
  expect_false(is.na(biserial(x_cont, x_bin, na_method = "pairwise")))
})

test_that("latent correlation print, plot, and summary methods work", {
  skip_if_not_installed("ggplot2")

  set.seed(1005)
  dat <- data.frame(
    a = sample(c(FALSE, TRUE), 120, replace = TRUE),
    b = sample(c(FALSE, TRUE), 120, replace = TRUE),
    c = sample(c(FALSE, TRUE), 120, replace = TRUE)
  )
  tc <- tetrachoric(dat)

  out <- capture.output(print(tc, digits = 3, max_rows = 2, max_cols = 2))
  expect_true(any(grepl("omitted", out)))

  p <- plot(tc, title = "Tetra plot", low_color = "blue", high_color = "red",
            mid_color = "white", value_text_size = 3)
  expect_s3_class(p, "ggplot")
  scale <- p$scales$get_scales("fill")
  expect_equal(scale$limits, c(-1, 1))

  sm <- summary(tc)
  expect_s3_class(sm, "summary.latent_corr")
  expect_s3_class(sm, "summary.matrixCorr")
  sum_out <- capture.output(print(sm, digits = 3))
  expect_true(any(grepl("Latent correlation summary", sum_out)))
  expect_true(any(grepl("thresholds", sum_out)))
  expect_false(any(grepl("zero_cells", sum_out)))
})

test_that("latent summaries keep all pairs for square but non-symmetric outputs", {
  set.seed(1006)
  X <- data.frame(x1 = rnorm(200), x2 = rnorm(200))
  Y <- data.frame(
    g1 = sample(c(FALSE, TRUE), 200, replace = TRUE),
    g2 = sample(c(FALSE, TRUE), 200, replace = TRUE)
  )

  bs <- biserial(X, Y)
  expect_false(isSymmetric(unclass(as.matrix(bs))))

  sm <- summary(bs)

  expect_s3_class(sm, "summary.latent_corr")
  expect_s3_class(sm, "summary.matrixCorr")
  expect_identical(sm$n_pairs, 4L)
  expect_s3_class(sm$top_results, "data.frame")
  expect_identical(nrow(sm$top_results), 4L)

  out <- capture.output(print(sm, digits = 3))
  expect_true(any(grepl("x1", out, fixed = TRUE)))
  expect_true(any(grepl("x2", out, fixed = TRUE)))
  expect_true(any(grepl("g1", out, fixed = TRUE)))
  expect_true(any(grepl("g2", out, fixed = TRUE)))
})

test_that("regular table-input latent fits match known regression values", {
  tab_2x2 <- as.table(matrix(c(61661, 1610, 85, 20), 2, 2))
  tab_3x3 <- as.table(matrix(c(13, 69, 41, 6, 113, 132, 0, 22, 104), 3, 3))

  expect_equal_abs(as.numeric(tetrachoric(tab_2x2, correct = 0.5)), 0.348073131513184)
  expect_equal_abs(as.numeric(polychoric(tab_2x2, correct = 0.5)), 0.348073131513616)

  pc_3x3 <- polychoric(tab_3x3, correct = 0)
  expect_equal_abs(as.numeric(pc_3x3), 0.491371943011718)
  thr <- attr(pc_3x3, "thresholds")
  expect_equal(unname(thr$row), c(-1.77438191034496, -0.135773931302112), tolerance = 1e-12)
  expect_equal(unname(thr$col), c(-0.687131286795469, 0.668209299725723), tolerance = 1e-12)
})

test_that("table-input latent fits attach optional inference payloads", {
  tab_2x2 <- as.table(matrix(c(61661, 1610, 85, 20), 2, 2))
  tab_3x3 <- as.table(matrix(c(13, 69, 41, 6, 113, 132, 0, 22, 104), 3, 3))

  tc <- tetrachoric(tab_2x2, ci = TRUE, p_value = TRUE)
  pc <- polychoric(tab_3x3, correct = 0, ci = TRUE, p_value = TRUE)

  expect_true(is.list(attr(tc, "ci", exact = TRUE)))
  expect_true(is.list(attr(tc, "inference", exact = TRUE)))
  expect_true(is.list(attr(pc, "ci", exact = TRUE)))
  expect_true(is.list(attr(pc, "inference", exact = TRUE)))
  expect_equal(as.numeric(tc), as.numeric(tetrachoric(tab_2x2, correct = 0.5)), tolerance = 1e-12)
  expect_equal(as.numeric(pc), as.numeric(polychoric(tab_3x3, correct = 0)), tolerance = 1e-12)
  expect_latent_wald_payloads(
    tc,
    ci_method = "wald_information_tetrachoric",
    inference_method = "wald_z_tetrachoric"
  )
  expect_latent_wald_payloads(
    pc,
    ci_method = "wald_information_polychoric",
    inference_method = "wald_z_polychoric"
  )
})

test_that("sparse zero-cell latent fits expose diagnostics and thresholds", {
  sparse_2x2 <- as.table(matrix(c(44268, 193, 14, 0), 2, 2))

  tc0 <- tetrachoric(sparse_2x2, correct = 0)
  pc0 <- polychoric(sparse_2x2, correct = 0)
  tc5 <- tetrachoric(sparse_2x2, correct = 0.5)

  expect_true(is.list(attr(tc0, "diagnostics")))
  expect_true(is.list(attr(pc0, "diagnostics")))
  expect_true(is.list(attr(tc0, "thresholds")))
  expect_true(is.list(attr(pc0, "thresholds")))

  expect_true(isTRUE(attr(tc0, "diagnostics")$boundary))
  expect_true(isTRUE(attr(pc0, "diagnostics")$boundary))
  expect_true(isTRUE(attr(tc5, "diagnostics")$corrected))
  expect_true(attr(tc0, "diagnostics")$zero_cells > 0L)
  expect_true(attr(pc0, "diagnostics")$zero_cells > 0L)
  expect_equal(attr(tc0, "diagnostics")$n_complete, sum(sparse_2x2))
  expect_equal(attr(pc0, "diagnostics")$n_complete, sum(sparse_2x2))

  expand_table_pairs <- function(tab) {
    tab <- as.matrix(tab)
    x <- integer(sum(tab))
    y <- integer(sum(tab))
    pos <- 1L
    for (i in seq_len(nrow(tab))) {
      for (j in seq_len(ncol(tab))) {
        nij <- tab[i, j]
        if (nij > 0L) {
          idx <- pos:(pos + nij - 1L)
          x[idx] <- i
          y[idx] <- j
          pos <- pos + nij
        }
      }
    }
    data.frame(
      a = ordered(x),
      b = ordered(y)
    )
  }
  sparse_df <- expand_table_pairs(sparse_2x2)
  pc_mat <- polychoric(sparse_df, na_method = "error")
  pc_mat0 <- polychoric(sparse_df, correct = 0, na_method = "error")
  diag_pc <- attr(pc_mat, "diagnostics")
  diag_pc0 <- attr(pc_mat0, "diagnostics")
  expect_true(is.list(diag_pc))
  expect_true(is.matrix(diag_pc$boundary))
  expect_true(is.matrix(diag_pc$zero_cells))
  expect_true(isTRUE(diag_pc$corrected[1, 2]))
  expect_true(isTRUE(diag_pc0$boundary[1, 2]))
  expect_true(is.list(attr(pc_mat, "thresholds")))
  expect_true(is.numeric(attr(pc_mat, "correct")))
})
