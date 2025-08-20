test_that("returns correct structure and handles simple cases", {
  mat <- matrix(c(1:5, 2:6, rnorm(5)), ncol = 3)
  colnames(mat) <- c("x", "y", "z")  # x and y are perfectly linearly related

  res <- distance_corr(mat)
  expect_s3_class(res, "distance_corr")
  expect_equal(dim(res), c(3, 3))
  expect_equal(colnames(res), c("x", "y", "z"))
  expect_equal(rownames(res), c("x", "y", "z"))

  # Diagonal values must be 1
  expect_true(all(diag(res) == 1))

  # Distance correlation between perfectly linear x and y should be near 1
  expect_gt(res["x", "y"], 0.95)

  # Distance correlation with noise should be significantly lower
  expect_lt(res["x", "z"], 0.9)
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
