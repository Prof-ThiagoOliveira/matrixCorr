test_that("returns correct structure and handles simple cases", {
  n <- 5
  mat <- cbind(
    x = 1:n,
    y = 2:(n+1),               # perfectly linear with x
    z = c(0.34, -0.75, 0.12, 1.09, -0.22)  # fixed values
  )

  res <- distance_corr(mat)

  expect_s3_class(res, "distance_corr")
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
  dcor <- distance_corr(X)[1, 2]

  # Pearson correlation near 0 due to symmetry
  expect_lt(abs(pearson), 0.02)

  # Distance correlation should be clearly positive
  expect_gt(dcor, 0.25)
})

test_that("handles constant variables and produces NAs", {
  mat <- cbind(rep(1, 10), rnorm(10))
  colnames(mat) <- c("const", "noise")

  res <- suppressWarnings(distance_corr(mat))
  expect_true(is.na(res["const", "noise"]))
})

test_that("matches expected dCor in AR(1) matrix", {
  set.seed(1)
  p <- 3; n <- 100; rho <- 0.7
  Sigma <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
  L <- chol(Sigma)
  X <- matrix(rnorm(n * p), n, p) %*% L
  colnames(X) <- paste0("V", 1:p)

  res <- distance_corr(X)

  # Expect symmetry and diagonal = 1
  expect_true(all.equal(res, t(res)))
  expect_true(all(diag(res) == 1))

  # dCor between adjacent variables should be stronger than non-adjacent
  expect_gt(res["V1", "V2"], res["V1", "V3"])
})

test_that("distance_corr print/plot cover optional parameters", {
  skip_if_not_installed("ggplot2")

  set.seed(303)
  X <- matrix(rnorm(60), nrow = 15, ncol = 4)
  colnames(X) <- paste0("D", seq_len(4))
  dc <- distance_corr(X)

  out <- capture.output(print(dc, digits = 3, max_rows = 2, max_cols = 3))
  expect_true(any(grepl("omitted", out)))

  p <- plot(dc, title = "Distance plot", low_color = "white", high_color = "navy", value_text_size = 3)
  expect_s3_class(p, "ggplot")
})

test_that("distance_corr rejects missing values by default", {
  X <- cbind(a = c(1, 2, NA, 4), b = c(1, 2, 3, 4), c = c(1, 2, 3, 4))
  expect_error(distance_corr(X), "Missing values are not allowed.")
})
