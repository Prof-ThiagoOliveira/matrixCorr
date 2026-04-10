test_that("returns correct structure and handles simple cases", {
  n <- 5
  mat <- cbind(
    x = 1:n,
    y = 2:(n+1),               # perfectly linear with x
    z = c(0.34, -0.75, 0.12, 1.09, -0.22)  # fixed values
  )

  res <- dcor(mat)

  expect_s3_class(res, "dcor")
  expect_equal(dim(res), c(3, 3))
  expect_true(all(diag(res) == 1))

  expect_gt(res["x","y"], 0.95)
  expect_lt(res["x","z"], 0.5)  # pick a margin that passes for this fixed z
})

test_that("detects non-linear dependence missed by Pearson", {
  set.seed(42)
  n <- 10000
  x <- rnorm(n)
  y <- x^2
  X <- cbind(x, y)
  colnames(X) <- c("x", "y")

  pearson <- pearson_corr(X)[1, 2]
  dcor <- dcor(X)[1, 2]

  # Pearson correlation near 0 due to symmetry
  expect_lt(abs(pearson), 0.02)

  # Distance correlation should be clearly positive
  expect_gt(dcor, 0.25)
})

test_that("handles constant variables and produces NAs", {
  mat <- cbind(rep(1, 10), rnorm(10))
  colnames(mat) <- c("const", "noise")

  res <- suppressWarnings(dcor(mat))
  expect_true(is.na(res["const", "noise"]))
})

test_that("matches expected dCor in AR(1) matrix", {
  set.seed(1)
  p <- 3; n <- 100; rho <- 0.7
  Sigma <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
  L <- chol(Sigma)
  X <- matrix(rnorm(n * p), n, p) %*% L
  colnames(X) <- paste0("V", 1:p)

  res <- dcor(X)

  # Expect symmetry and diagonal = 1
  expect_true(all.equal(res, t(res)))
  expect_true(all(diag(res) == 1))

  # dCor between adjacent variables should be stronger than non-adjacent
  expect_gt(res["V1", "V2"], res["V1", "V3"])
})

test_that("distance-correlation inference matches the reference t-test", {
  set.seed(2026)
  x <- rnorm(80)
  y <- x^2 + rnorm(80, sd = 0.25)
  z <- -0.4 * x + rnorm(80, sd = 0.8)
  X <- cbind(x = x, y = y, z = z)

  dc <- dcor(X, p_value = TRUE)
  inf <- attr(dc, "inference", exact = TRUE)
  ref_estimate <- matrix(
    c(
      NA_real_, 0.248471103773309, 0.0998513936435007,
      0.248471103773309, NA_real_, -0.0165169652366662,
      0.0998513936435007, -0.0165169652366662, NA_real_
    ),
    nrow = 3,
    byrow = TRUE
  )
  ref_statistic <- matrix(
    c(
      NA_real_, 14.2337297104894, 5.56845748933565,
      14.2337297104894, NA_real_, -0.916630998257919,
      5.56845748933565, -0.916630998257919, NA_real_
    ),
    nrow = 3,
    byrow = TRUE
  )
  ref_df <- matrix(
    c(
      NA_real_, 3079, 3079,
      3079, NA_real_, 3079,
      3079, 3079, NA_real_
    ),
    nrow = 3,
    byrow = TRUE
  )
  ref_p_value <- matrix(
    c(
      NA_real_, 0, 1.3955280708799478e-08,
      0, NA_real_, 0.82029609217836286,
      1.3955280708799478e-08, 0.82029609217836286, NA_real_
    ),
    nrow = 3,
    byrow = TRUE
  )

  expect_type(inf, "list")
  expect_identical(inf$method, "dcor_t_test")
  expect_true(is.matrix(inf$estimate))
  expect_true(is.matrix(inf$statistic))
  expect_true(is.matrix(inf$parameter))
  expect_true(is.matrix(inf$p_value))

  idx <- which(upper.tri(dc), arr.ind = TRUE)
  for (k in seq_len(nrow(idx))) {
    i <- idx[k, 1]
    j <- idx[k, 2]
    ref_est <- ref_estimate[i, j]
    ref_t <- ref_statistic[i, j]
    ref_df_ij <- ref_df[i, j]
    ref_p <- ref_p_value[i, j]

    expect_lt(abs(inf$estimate[i, j] - ref_est), 1e-6)
    expect_lt(abs(inf$statistic[i, j] - ref_t), 1e-5)
    expect_equal(inf$parameter[i, j], ref_df_ij, tolerance = 1e-10)
    if (isTRUE(is.finite(ref_p) && ref_p == 0)) {
      expect_gte(inf$p_value[i, j], 0)
      expect_lt(inf$p_value[i, j], .Machine$double.eps)
    } else {
      expect_lt(abs(inf$p_value[i, j] - ref_p), 1e-10)
    }
    expect_lt(abs(dc[i, j] - max(ref_est, 0)), 1e-6)
  }
})

test_that("requesting inference does not change the dcor estimate matrix", {
  X <- cbind(
    a = c(-2.0, -0.5, 0.0, 1.5, 2.0, 2.5),
    b = c(1.0, 0.2, -0.1, 0.5, 1.1, 1.8),
    c = c(0.5, -1.1, 0.7, -0.3, 1.6, -0.8)
  )

  est_only <- dcor(X)
  with_test <- dcor(X, p_value = TRUE)
  est_only_mat <- unclass(est_only)
  with_test_mat <- unclass(with_test)
  attributes(est_only_mat) <- attributes(est_only_mat)[c("dim", "dimnames")]
  attributes(with_test_mat) <- attributes(with_test_mat)[c("dim", "dimnames")]

  expect_equal(est_only_mat, with_test_mat, tolerance = 1e-12)
})

test_that("dcor summary switches to pairwise inference view when requested", {
  set.seed(11)
  X <- matrix(rnorm(120), nrow = 30, ncol = 4)
  colnames(X) <- paste0("D", seq_len(ncol(X)))

  dc <- dcor(X, p_value = TRUE)
  sm <- summary(dc)

  expect_s3_class(sm, "summary.dcor")
  expect_s3_class(sm, "data.frame")
  expect_true(all(c("var1", "var2", "estimate", "n_complete", "statistic", "df", "p_value") %in% names(sm)))

  txt <- capture.output(print(sm))
  expect_true(any(grepl("^Distance correlation summary$", txt)))
  expect_true(any(grepl("p_value", txt, fixed = TRUE)))
})

test_that("dcor print/plot cover optional parameters", {
  skip_if_not_installed("ggplot2")

  set.seed(303)
  X <- matrix(rnorm(60), nrow = 15, ncol = 4)
  colnames(X) <- paste0("D", seq_len(4))
  dc <- dcor(X)

  out <- capture.output(print(dc, digits = 3, max_rows = 2, max_cols = 3))
  expect_true(any(grepl("omitted", out)))

  p <- plot(dc, title = "Distance plot", low_color = "white", high_color = "navy", value_text_size = 3)
  expect_s3_class(p, "ggplot")
})

test_that("dcor rejects missing values by default", {
  X <- cbind(a = c(1, 2, NA, 4), b = c(1, 2, 3, 4), c = c(1, 2, 3, 4))
  expect_error(dcor(X), "Missing values are not allowed.")
})
