test_that("robust_dcor returns basic shape and class", {
  set.seed(1)
  X <- matrix(rnorm(100 * 5), ncol = 5)
  colnames(X) <- paste0("V", 1:5)

  R <- robust_dcor(X)

  expect_s3_class(R, "robust_dcor")
  expect_s3_class(R, "corr_matrix")
  expect_equal(dim(R), c(5L, 5L))
  expect_true(isSymmetric(unclass(R)))
  expect_equal(unname(diag(R)), rep(1, 5), tolerance = 1e-12)
})

test_that("robust_dcor routes output modes", {
  set.seed(1)
  X <- matrix(rnorm(100 * 5), ncol = 5)
  colnames(X) <- paste0("V", 1:5)

  dense <- robust_dcor(X, output = "matrix")
  sparse <- robust_dcor(X, output = "sparse", threshold = 0.2)
  edges <- robust_dcor(X, output = "edge_list", threshold = 0.2)

  expect_s3_class(dense, "corr_matrix")
  expect_s4_class(sparse, "sparseMatrix")
  expect_identical(attr(sparse, "corr_output_class", exact = TRUE), "corr_sparse")
  expect_s3_class(edges, "corr_edge_list")
})

test_that("robust_dcor rejects missing values with na_method error", {
  set.seed(1)
  X <- matrix(rnorm(100 * 5), ncol = 5)
  colnames(X) <- paste0("V", 1:5)
  X[1, 1] <- NA_real_

  expect_error(robust_dcor(X, na_method = "error"))
})

test_that("robust_dcor supports pairwise missing data", {
  set.seed(1)
  X <- matrix(rnorm(100 * 5), ncol = 5)
  colnames(X) <- paste0("V", 1:5)
  X[1, 1] <- NA_real_

  Rna <- robust_dcor(X, na_method = "pairwise")

  expect_s3_class(Rna, "robust_dcor")
  expect_true(!is.null(attr(Rna, "diagnostics")$n_complete))
})

test_that("robust_dcor detects nonlinear dependence", {
  set.seed(2)
  x <- rnorm(200)
  y <- x^2 + rnorm(200, sd = 0.2)

  R <- robust_dcor(cbind(x = x, y = y))

  expect_true(R["x", "y"] > 0.1)
})

test_that("robust_dcor is less outlier-sensitive than dcor", {
  set.seed(11)
  x <- rnorm(100)
  y <- rnorm(100)
  x[1] <- 50
  y[1] <- 50

  Dc <- dcor(cbind(x = x, y = y))
  Dr <- robust_dcor(cbind(x = x, y = y))

  expect_true(Dr["x", "y"] < Dc["x", "y"])
})

test_that("robust_dcor uses shared S3 methods", {
  skip_if_not_installed("ggplot2")

  set.seed(2)
  x <- rnorm(200)
  y <- x^2 + rnorm(200, sd = 0.2)
  R <- robust_dcor(cbind(x = x, y = y))

  expect_s3_class(summary(R), "summary.corr_result")
  expect_s3_class(plot(R), "ggplot")
  expect_true(length(capture.output(print(R))) > 0L)
})

test_that("robust_dcor supports permutation p-values", {
  set.seed(4)
  x <- rnorm(30)
  y <- x^2 + rnorm(30, sd = 0.2)

  R <- robust_dcor(cbind(x = x, y = y), p_value = TRUE, inference = "permutation", n_perm = 19L, seed = 1L)
  inf <- attr(R, "inference", exact = TRUE)

  expect_identical(inf$method, "permutation")
  expect_true(is.matrix(inf$p_value))
  expect_true(is.finite(inf$p_value["x", "y"]))
})

test_that("robust_dcor infers p-values from permutation inference", {
  set.seed(5)
  x <- rnorm(24)
  y <- x^2 + rnorm(24, sd = 0.3)

  R <- robust_dcor(cbind(x = x, y = y), inference = "permutation", n_perm = 9L, seed = 2L)
  inf <- attr(R, "inference", exact = TRUE)

  expect_identical(inf$method, "permutation")
  expect_true(is.matrix(inf$p_value))
  expect_true(is.finite(inf$p_value["x", "y"]))
})

test_that("robust_dcor warns for large permutation workloads", {
  expect_warning(
    .mc_warn_large_robust_dcor_permutation_request(n_cols = 6L, n_perm = 10L, warn_at = 100L),
    class = "matrixCorr_permutation_warning"
  )
})
