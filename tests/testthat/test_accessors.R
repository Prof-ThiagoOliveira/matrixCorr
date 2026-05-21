test_that("estimate returns a plain numeric matrix for dense correlation results", {
  set.seed(101)
  X <- cbind(
    a = rnorm(30),
    b = rnorm(30),
    c = rnorm(30)
  )

  fit <- pearson_corr(X, ci = TRUE)
  est <- estimate(fit)

  expect_true(is.matrix(est))
  expect_type(est, "double")
  expect_false(inherits(est, "corr_matrix"))
  expect_equal(est, unclass(as.matrix(fit)), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(dim(est), dim(fit))
  expect_equal(dimnames(est), dimnames(fit))
  expect_equal(coef(fit), est, tolerance = 1e-12)
})

test_that("tidy returns unique upper-triangle pairs by default", {
  set.seed(102)
  X <- cbind(
    a = rnorm(30),
    b = rnorm(30),
    c = rnorm(30)
  )

  fit <- pearson_corr(X, ci = TRUE)
  td <- tidy(fit)

  expect_s3_class(td, "data.frame")
  expect_true(all(c("item1", "item2", "estimate") %in% names(td)))
  expect_false(any(td$item1 == td$item2))
  expect_equal(nrow(td), choose(ncol(X), 2L))
  expect_false(any(duplicated(paste(pmin(td$item1, td$item2), pmax(td$item1, td$item2)))))
})

test_that("tidy triangle and diagonal controls are predictable", {
  set.seed(103)
  X <- cbind(
    a = rnorm(20),
    b = rnorm(20),
    c = rnorm(20)
  )

  fit <- pearson_corr(X, ci = TRUE)
  p <- ncol(X)

  expect_equal(nrow(tidy(fit)), choose(p, 2L))
  expect_equal(nrow(tidy(fit, diag = TRUE)), choose(p, 2L) + p)
  expect_equal(nrow(tidy(fit, triangle = "full")), p * (p - 1L))
  expect_equal(nrow(tidy(fit, triangle = "full", diag = TRUE)), p * p)
  expect_true(any(tidy(fit, diag = TRUE)$item1 == tidy(fit, diag = TRUE)$item2))
})

test_that("confint and ci expose confidence intervals cleanly", {
  set.seed(104)
  X <- cbind(
    a = rnorm(30),
    b = rnorm(30),
    c = rnorm(30)
  )

  fit <- pearson_corr(X, ci = TRUE, conf_level = 0.95)
  ci_df <- confint(fit)

  expect_identical(ci(fit), attr(fit, "ci", exact = TRUE))
  expect_s3_class(ci_df, "data.frame")
  expect_true(all(c("item1", "item2", "lwr", "upr") %in% names(ci_df)))
  expect_false("conf.level" %in% names(ci_df))
  expect_false("conf_level" %in% names(ci_df))
  expect_false("ci.method" %in% names(ci_df))
  expect_no_error(confint(fit, level = attr(fit, "conf.level", exact = TRUE)))
  expect_error(
    confint(fit, level = 0.90),
    "This object contains pre-computed confidence intervals at conf.level"
  )
})

test_that("ci and confint fail clearly when intervals are absent", {
  set.seed(105)
  X <- cbind(
    a = rnorm(20),
    b = rnorm(20),
    c = rnorm(20)
  )

  fit <- pearson_corr(X, ci = FALSE)

  expect_error(ci(fit), "No confidence intervals are available for this object.")
  expect_error(confint(fit), "does not contain confidence limits")
})

test_that("estimate handles edge-list outputs predictably", {
  set.seed(106)
  X <- cbind(
    a = rnorm(20),
    b = rnorm(20),
    c = rnorm(20)
  )

  fit <- pearson_corr(X, output = "edge_list")
  est <- estimate(fit)

  expect_s3_class(est, "data.frame")
  expect_equal(names(est), c("item1", "item2", "estimate"))
  expect_false(inherits(est, "corr_edge_list"))
})

test_that("generic tidy preserves directed matrix orientation", {
  x <- seq(-1, 1, length.out = 80)
  X <- cbind(x = x, x2 = x^2)

  fit <- xi_corr(X, tie_method = "first", ci = FALSE)
  td <- tidy(fit)

  expect_equal(estimate(fit), unclass(as.matrix(fit)), tolerance = 1e-12, ignore_attr = TRUE)
  expect_true(all(c("x", "x2") %in% td$item1))
  expect_true(all(c("x", "x2") %in% td$item2))
  expect_equal(
    td$estimate[td$item1 == "x" & td$item2 == "x2"],
    unname(fit["x", "x2"]),
    tolerance = 1e-12
  )
  expect_equal(
    td$estimate[td$item1 == "x2" & td$item2 == "x"],
    unname(fit["x2", "x"]),
    tolerance = 1e-12
  )
})

test_that("generic accessors work for scalar CI results", {
  set.seed(107)
  x <- runif(120)
  y <- sin(2 * pi * x) + rnorm(120, sd = 0.2)

  fit <- xi_corr(
    x,
    y,
    ci = TRUE,
    ci_method = "dette_kroll",
    bootstrap_reps = 9,
    seed = 107,
    tie_method = "first"
  )

  td <- tidy(fit)
  ci_df <- confint(fit)

  expect_equal(estimate(fit), unname(as.numeric(fit)), tolerance = 1e-12)
  expect_true(all(c("estimate", "lwr", "upr", "m", "bootstrap_reps") %in% names(td)))
  expect_false("ci.method" %in% names(td))
  expect_false("n_complete" %in% names(td))
  expect_true(all(c("lwr", "upr") %in% names(ci_df)))
  expect_false("conf.level" %in% names(ci_df))
  expect_false("ci.method" %in% names(ci_df))
  expect_false("conf_level" %in% names(ci_df))
  expect_false("ci_method" %in% names(ci_df))
  expect_equal(ci_df$lwr, td$lwr, tolerance = 1e-12)
  expect_equal(ci_df$upr, td$upr, tolerance = 1e-12)
})

test_that("generic accessors work for scalar agreement results", {
  ratings <- data.frame(
    r1 = c(1, 1, 2, 2, 3, 3, 1, 2),
    r2 = c(1, 2, 2, 2, 3, 1, 1, 2)
  )

  fit <- cohen_kappa(ratings$r1, ratings$r2, ci = TRUE)
  td <- tidy(fit)
  ci_df <- confint(fit)

  expect_equal(estimate(fit), unname(as.numeric(fit)), tolerance = 1e-12)
  expect_true(all(c("estimate", "lwr", "upr") %in% names(td)))
  expect_true(all(c("lwr", "upr") %in% names(ci_df)))
})

test_that("generic accessors work for overall CIA results", {
  set.seed(110)
  subjects <- sprintf("s%02d", 1:10)
  methods <- c("A", "B", "C")
  replicates <- sprintf("r%02d", 1:4)
  dat <- expand.grid(
    subject = subjects,
    method = methods,
    replicate = replicates,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  subject_effect <- stats::rnorm(length(subjects), sd = 3)
  method_shift <- c(A = 0, B = 0.15, C = -0.10)
  dat$value <- subject_effect[match(dat$subject, subjects)] +
    method_shift[dat$method] +
    stats::rnorm(nrow(dat), sd = 0.4)

  fit <- cia(
    dat,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    scope = "overall"
  )
  td <- tidy(fit)

  expect_equal(estimate(fit), fit$cia, tolerance = 1e-12)
  expect_equal(estimate(summary(fit)), fit$cia, tolerance = 1e-12)
  expect_true(all(c("estimate", "cia_raw", "tau2") %in% names(td)))
  expect_equal(td$estimate, fit$cia, tolerance = 1e-12)
})

test_that("generic tidy handles list-style CI outputs", {
  X <- cbind(
    a = c(1, 2, 3, 4, 5, 6),
    b = c(1.1, 2.1, 2.9, 4.2, 4.9, 6.1),
    c = c(6, 5, 4, 3, 2, 1)
  )

  fit <- icc(X, ci = TRUE)
  est <- estimate(fit)
  td <- tidy(fit)
  ci_df <- confint(fit)

  expect_true(is.matrix(est))
  expect_false(inherits(est, "icc"))
  expect_equal(dim(est), dim(fit$est))
  expect_equal(dimnames(est), dimnames(fit$est))
  expect_true(all(c("item1", "item2", "estimate") %in% names(td)))
  expect_true(all(c("item1", "item2", "lwr", "upr") %in% names(ci_df)))
  expect_false("n_complete" %in% names(td))
  expect_false(any(td$item1 == td$item2))
  expect_equal(nrow(td), choose(ncol(X), 2L))
})

test_that("estimate preserves matrix shape for matrixCorr matrix-like classes", {
  X <- cbind(
    m1 = c(1, 2, 3, 4, 5, 6),
    m2 = c(1.1, 2.1, 2.9, 4.2, 4.9, 6.1),
    m3 = c(6, 5, 4, 3, 2, 1)
  )

  fit <- icc(
    X,
    model = "twoway_random",
    type = "agreement",
    unit = "single",
    scope = "pairwise"
  )
  est <- estimate(fit)

  expect_true(is.matrix(est))
  expect_false(inherits(est, "icc"))
  expect_equal(dim(est), dim(fit))
  expect_equal(dimnames(est), dimnames(fit))
})

test_that("ccc_glmm accessors return matrices and CI payloads", {
  set.seed(109)
  n_subjects <- 10L
  reps <- 2L
  dat <- data.frame(
    subject = factor(rep(seq_len(n_subjects), each = 2L * reps)),
    method = factor(rep(rep(c("A", "B"), each = reps), times = n_subjects)),
    replicate = factor(rep(seq_len(reps), times = 2L * n_subjects))
  )
  subject_eff <- rnorm(n_subjects, 0, 0.5)
  dat$eta <- 1.1 + subject_eff[as.integer(dat$subject)] +
    ifelse(dat$method == "B", 0.1, 0)
  dat$y <- rpois(nrow(dat), exp(dat$eta))

  fit <- ccc_glmm(
    dat,
    "y",
    "subject",
    "method",
    replicate = "replicate",
    ci = TRUE
  )
  est <- estimate(fit)
  payload <- ci(fit)

  expect_true(is.matrix(est))
  expect_false(inherits(est, "ccc_glmm"))
  expect_equal(dim(est), dim(fit))
  expect_equal(dimnames(est), dimnames(fit))
  expect_equal(coef(fit), est, tolerance = 1e-12)
  expect_true(is.list(payload))
  expect_false(is.logical(payload))
  expect_true(all(c("lwr", "upr", "se", "total", "inter", "intra_method1", "intra_method2") %in% names(payload)))
  expect_true(is.matrix(payload$lwr))
  expect_true(is.matrix(payload$upr))
  expect_equal(dim(payload$lwr), dim(fit))
  expect_equal(dim(payload$upr), dim(fit))
})

test_that("ci works for agreement list and Bland-Altman matrix outputs", {
  set.seed(108)
  Sigma <- matrix(c(1, 0.5, 0.3, 0.5, 1, 0.4, 0.3, 0.4, 1), nrow = 3)
  X <- MASS::mvrnorm(n = 60, mu = rep(0, 3), Sigma = Sigma)

  ccc_fit <- ccc(X, ci = TRUE)
  ccc_ci <- ci(ccc_fit)
  ccc_td <- tidy(ccc_fit)
  ccc_limits <- confint(ccc_fit)

  expect_true(all(c("lwr", "upr") %in% names(ccc_ci)))
  expect_false(any(ccc_td$item1 == ccc_td$item2))
  expect_equal(nrow(ccc_td), choose(ncol(X), 2L))
  expect_equal(nrow(ccc_limits), choose(ncol(X), 2L))

  wide <- data.frame(
    ref = rnorm(50, 100, 8),
    m2 = rnorm(50, 101, 8),
    m3 = rnorm(50, 99, 9)
  )
  ba_fit <- ba(wide)
  ba_ci <- ci(ba_fit)
  ba_limits <- confint(ba_fit)

  expect_true(all(c("bias", "loa_low", "loa_up", "conf.level") %in% names(ba_ci)))
  expect_true(all(c("bias_lwr", "bias_upr", "lo_lwr", "lo_upr", "up_lwr", "up_upr") %in% names(ba_limits)))
  expect_false("conf.level" %in% names(ba_limits))
})
