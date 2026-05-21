manual_xi_no_y_ties <- function(x, y) {
  o <- order(x)
  r <- rank(y, ties.method = "first")
  n <- length(x)
  1 - 3 * sum(abs(diff(r[o]))) / (n^2 - 1)
}

test_that("chatterjee_xi matches deterministic monotone values without y ties", {
  n <- 25
  expected <- (n - 2) / (n + 1)

  expect_equal(
    xi_corr(1:n, 1:n, tie_method = "first"),
    expected,
    tolerance = 1e-12
  )
  expect_equal(
    xi_corr(1:n, n:1, tie_method = "first"),
    expected,
    tolerance = 1e-12
  )
})

test_that("chatterjee_xi captures directed non-monotone deterministic dependence", {
  n <- 401
  x <- seq(-1, 1, length.out = n)
  y <- x^2
  xy <- xi_corr(x, y, tie_method = "first")
  yx <- xi_corr(y, x, tie_method = "first")

  expect_gt(xy, 0.75)
  expect_lt(yx, xy - 0.20)
})

test_that("chatterjee_xi is near zero under independence", {
  set.seed(12)
  x <- rnorm(2000)
  y <- rnorm(2000)
  expect_lt(abs(xi_corr(x, y, tie_method = "first")), 0.08)
})

test_that("chatterjee_xi handles ties and optional bias correction", {
  expect_true(is.na(xi_corr(1:5, rep(1, 5), tie_method = "first")))

  x <- c(1, 1, 1, 2, 2, 3, 3, 3)
  y <- c(2, 4, 1, 3, 5, 8, 7, 6)
  first_1 <- xi_corr(x, y, tie_method = "first")
  first_2 <- xi_corr(x, y, tie_method = "first")
  expect_identical(first_1, first_2)

  random_1 <- xi_corr(x, y, tie_method = "random", seed = 99)
  random_2 <- xi_corr(x, y, tie_method = "random", seed = 99)
  expect_identical(random_1, random_2)

  n <- 20
  expected <- (n - 2) / (n + 1)
  expect_equal(
    xi_corr(1:n, 1:n, tie_method = "first", bias_correction = "none"),
    expected,
    tolerance = 1e-12
  )
  expect_equal(
    xi_corr(1:n, 1:n, tie_method = "first", bias_correction = "upper_bound"),
    1,
    tolerance = 1e-12
  )
})

test_that("chatterjee_xi matrix orientation is row predictor and column response", {
  x <- seq(-1, 1, length.out = 101)
  X <- cbind(x = x, square = x^2, reverse = rev(x))
  fit <- xi_corr(X, tie_method = "first")

  expect_equal(
    unname(fit["x", "square"]),
    xi_corr(X[, "x"], X[, "square"], tie_method = "first"),
    tolerance = 1e-12
  )
  expect_equal(
    unname(fit["square", "x"]),
    xi_corr(X[, "square"], X[, "x"], tie_method = "first"),
    tolerance = 1e-12
  )
  expect_false(isSymmetric(unclass(fit)))
  expect_false(isTRUE(attr(fit, "corr_symmetric", exact = TRUE)))
})

test_that("chatterjee_xi missing-data modes follow package contracts", {
  X <- cbind(
    a = c(1, 2, 3, 4, NA, 6),
    b = c(1, 4, 9, 16, 25, 36),
    c = c(6, 5, Inf, 3, 2, 1)
  )

  expect_error(xi_corr(X, na_method = "error"), "Missing values|NA|Inf")

  complete <- xi_corr(X, na_method = "complete", tie_method = "first")
  d_complete <- attr(complete, "diagnostics", exact = TRUE)
  expect_identical(d_complete$na_method, "complete")
  expect_identical(d_complete$n_complete, 4L)

  pairwise <- xi_corr(X, na_method = "pairwise", tie_method = "first")
  d_pairwise <- attr(pairwise, "diagnostics", exact = TRUE)
  expect_true(is.matrix(d_pairwise$n_complete))
  expect_identical(unname(d_pairwise$n_complete["a", "b"]), 5L)
  expect_identical(unname(d_pairwise$n_complete["a", "c"]), 4L)
})

test_that("chatterjee_xi confidence intervals are finite and auto-select methods", {
  set.seed(42)
  X <- cbind(a = rnorm(80), b = rnorm(80), c = rnorm(80))
  fit <- xi_corr(
    X,
    ci = TRUE,
    ci_method = "dette_kroll",
    bootstrap_reps = 19,
    seed = 1,
    tie_method = "first"
  )
  ci <- attr(fit, "ci", exact = TRUE)
  expect_true(is.finite(ci$lwr.ci["a", "b"]))
  expect_true(is.finite(ci$upr.ci["a", "b"]))
  expect_identical(as.character(ci$ci.method["a", "b"]), "dette_kroll")

  big <- cbind(a = rnorm(1001), b = rnorm(1001))
  fit_big <- xi_corr(
    big,
    ci = TRUE,
    ci_method = "n_choose_m",
    bootstrap_reps = 9,
    seed = 2,
    tie_method = "first"
  )
  ci_big <- attr(fit_big, "ci", exact = TRUE)
  expect_true(is.finite(ci_big$lwr.ci["a", "b"]))
  expect_true(is.finite(ci_big$upr.ci["a", "b"]))
  expect_lte(ci_big$lwr.ci["a", "b"], fit_big["a", "b"])
  expect_gte(ci_big$upr.ci["a", "b"], fit_big["a", "b"])
  expect_identical(as.character(ci_big$ci.method["a", "b"]), "n_choose_m")

  set.seed(1299)
  x <- runif(1200, -1, 1)
  y <- sin(2 * pi * x) + rnorm(1200, sd = 0.2)
  scalar_big <- xi_corr(
    x,
    y,
    ci = TRUE,
    ci_method = "n_choose_m",
    bootstrap_reps = 19,
    seed = 123,
    tie_method = "first"
  )
  scalar_ci <- attr(scalar_big, "ci", exact = TRUE)
  scalar_est <- unname(as.numeric(scalar_big))
  expect_lte(scalar_ci$lwr, scalar_est)
  expect_gte(scalar_ci$upr, scalar_est)
  expect_s3_class(scalar_big, "chatterjee_xi_scalar")
  scalar_print <- capture.output(print(scalar_big))
  expect_true(any(grepl("Chatterjee xi directed-dependence estimate", scalar_print, fixed = TRUE)))
  expect_true(any(grepl("interval", scalar_print, fixed = TRUE)))
  scalar_summary <- summary(scalar_big)
  expect_s3_class(scalar_summary, "summary.chatterjee_xi_scalar")
  expect_equal(scalar_summary$estimate, scalar_est, tolerance = 1e-12)
  expect_equal(scalar_summary$lwr, scalar_ci$lwr, tolerance = 1e-12)
  summary_print <- capture.output(print(scalar_summary))
  expect_true(any(grepl("Chatterjee xi directed-dependence estimate", summary_print, fixed = TRUE)))
  expect_true(any(grepl("ci_method", summary_print, fixed = TRUE)))

  small_auto <- xi_corr(
    cbind(a = rnorm(1000), b = rnorm(1000)),
    ci = TRUE,
    ci_method = "auto",
    bootstrap_reps = 9,
    seed = 3,
    tie_method = "first"
  )
  large_auto <- xi_corr(
    cbind(a = rnorm(1001), b = rnorm(1001)),
    ci = TRUE,
    ci_method = "auto",
    bootstrap_reps = 9,
    seed = 3,
    tie_method = "first"
  )
  expect_identical(as.character(attr(small_auto, "ci")$ci.method["a", "b"]), "dette_kroll")
  expect_identical(as.character(attr(large_auto, "ci")$ci.method["a", "b"]), "n_choose_m")

  rep1 <- xi_corr(X, ci = TRUE, bootstrap_reps = 9, seed = 10, tie_method = "first")
  rep2 <- xi_corr(X, ci = TRUE, bootstrap_reps = 9, seed = 10, tie_method = "first")
  expect_equal(attr(rep1, "ci", exact = TRUE), attr(rep2, "ci", exact = TRUE))

  sm <- summary(fit)
  expect_s3_class(sm, "summary.chatterjee_xi")
  expect_true(all(c("lwr", "upr", "ci_method", "m", "se", "bootstrap_reps") %in% names(sm)))
  expect_identical(attr(sm, "ci_method", exact = TRUE), "dette_kroll")
})

test_that("chatterjee_xi output classes and directed output modes are correct", {
  x <- seq(-1, 1, length.out = 51)
  X <- cbind(a = x, b = x^2, c = rev(x))
  dense <- xi_corr(X, tie_method = "first")
  mat <- as.matrix(dense)

  expect_s3_class(dense, "corr_matrix")
  expect_s3_class(dense, "chatterjee_xi")
  expect_identical(attr(dense, "method", exact = TRUE), "chatterjee_xi")
  expect_false(isTRUE(attr(dense, "corr_symmetric", exact = TRUE)))

  edge <- xi_corr(X, tie_method = "first", output = "edge_list", threshold = 0, diag = FALSE)
  expect_s3_class(edge, "corr_edge_list")
  edge_df <- as.data.frame(edge, stringsAsFactors = FALSE)
  expect_equal(nrow(edge_df), ncol(X) * (ncol(X) - 1L))
  expect_equal(
    edge_df$value[edge_df$row == "b" & edge_df$col == "a"],
    unname(mat["b", "a"])
  )

  sparse <- xi_corr(X, tie_method = "first", output = "sparse", threshold = 0, diag = FALSE)
  expect_s4_class(sparse, "sparseMatrix")
  expect_false(isTRUE(attr(sparse, "corr_symmetric", exact = TRUE)))
  expect_equal(as.matrix(sparse)["a", "b"], mat["a", "b"], tolerance = 1e-12)
  expect_equal(as.matrix(sparse)["b", "a"], mat["b", "a"], tolerance = 1e-12)
})
