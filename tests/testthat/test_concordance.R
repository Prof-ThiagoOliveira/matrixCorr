test_that("ccc_cpp and lins_ccc_detail_cpp produce expected values", {
  # Test 1: Perfect agreement
  x1 <- c(1, 2, 3, 4, 5)
  mat1 <- cbind(x1, x1)
  ccc_mat <- ccc_cpp(mat1)
  expect_equal(ccc_mat[1, 2], 1.0, tolerance = 1e-12)

  # Test 2: Perfect negative agreement
  x2 <- c(1, 2, 3, 4, 5)
  y2 <- c(5, 4, 3, 2, 1)
  mat2 <- cbind(x2, y2)
  ccc_val2 <- ccc_cpp(mat2)[1, 2]
  expected_r2 <- -1
  mean_diff2 <- mean(x2) - mean(y2)
  var_x2 <- var(x2) * 4 / 5  # biased var
  var_y2 <- var(y2) * 4 / 5
  expected_sxy2 <- expected_r2 * sqrt(var_x2 * var_y2)
  expected_ccc2 <- 2 * expected_sxy2 / (var_x2 + var_y2 + mean_diff2^2)
  expect_equal(ccc_val2, expected_ccc2, tolerance = 1e-12)

  # Test 3: Manually calculated CCC
  x3 <- c(1, 2, 3, 4, 5)
  y3 <- c(1.1, 2.1, 2.9, 3.8, 5.2)
  mean_x3 <- mean(x3)
  mean_y3 <- mean(y3)
  var_x3 <- var(x3) * 4 / 5
  var_y3 <- var(y3) * 4 / 5
  cov_xy3 <- sum((x3 - mean_x3) * (y3 - mean_y3)) / 5
  r3 <- cov_xy3 / sqrt(var_x3 * var_y3)
  sxy3 <- r3 * sqrt(var_x3 * var_y3)
  expected_ccc3 <- 2 * sxy3 / (var_x3 + var_y3 + (mean_x3 - mean_y3)^2)

  mat3 <- cbind(x3, y3)
  fast_ccc3 <- ccc_cpp(mat3)[1, 2]
  expect_equal(fast_ccc3, expected_ccc3, tolerance = 1e-10)
})

test_that("CCC matrix is symmetric", {
  mat <- matrix(rnorm(100 * 4), ncol = 4)
  result <- ccc_cpp(mat)
  expect_equal(result, t(result), tolerance = 1e-12)
})

test_that("CCC detects lack of agreement despite perfect correlation", {
  x <- rnorm(100)
  y <- x * 2 + 5
  mat <- cbind(x, y)
  result <- ccc_cpp(mat)
  expect_true(result[1, 2] < 1)
})

test_that("CCC with CI returns correctly structured result", {
  mat <- matrix(rnorm(100 * 3), ncol = 3)
  result <- ccc_with_ci_cpp(mat)
  expect_named(result, c("est", "lwr.ci", "upr.ci"))
  expect_equal(dim(result$est), dim(result$lwr.ci))
  expect_equal(dim(result$est), dim(result$upr.ci))
})

test_that("ccc public API honors n_threads without changing estimates", {
  set.seed(44)
  X <- matrix(rnorm(160), ncol = 4)
  colnames(X) <- paste0("V", 1:4)

  fit1 <- ccc(X, n_threads = 1L)
  fit2 <- ccc(X, n_threads = 2L)
  fit1_ci <- ccc(X, ci = TRUE, n_threads = 1L)
  fit2_ci <- ccc(X, ci = TRUE, n_threads = 2L)

  expect_equal(unname(fit1), unname(fit2), tolerance = 1e-12)
  expect_equal(unname(fit1_ci$est), unname(fit2_ci$est), tolerance = 1e-12)
})

test_that("ccc keeps only canonical CI components when intervals are requested", {
  set.seed(45)
  X <- matrix(rnorm(120), ncol = 3)
  colnames(X) <- c("A", "B", "C")

  fit <- ccc(X)
  fit_ci <- ccc(X, ci = TRUE)

  expect_false(is.list(fit))
  expect_s3_class(fit, "ccc")
  expect_s3_class(fit_ci, "ccc_ci")
  expect_identical(names(fit_ci), c("est", "lwr.ci", "upr.ci"))
  expect_false(any(c("estimate", "lwr", "upr") %in% names(fit_ci)))
  expect_true(is.matrix(attr(fit, "diagnostics")$n_complete))
  expect_true(is.matrix(attr(fit_ci, "diagnostics")$n_complete))
})

test_that("ccc missing-data modes match Pearson/Spearman semantics", {
  X <- cbind(
    a = c(1, 2, 3, 4, 5, NA),
    b = c(1, 2, NA, 4, 5, 6),
    c = c(6, 5, 4, 3, 2, 1)
  )

  expect_error(ccc(X, na_method = "error"), "Missing")

  keep_all <- apply(is.finite(X), 1L, all)
  fit_complete <- ccc(X, na_method = "complete")
  ref_complete <- ccc(X[keep_all, , drop = FALSE])
  got_complete <- unclass(as.matrix(fit_complete))
  ref_complete_mat <- unclass(as.matrix(ref_complete))
  attributes(got_complete) <- attributes(got_complete)[c("dim", "dimnames")]
  attributes(ref_complete_mat) <- attributes(ref_complete_mat)[c("dim", "dimnames")]
  expect_equal(got_complete, ref_complete_mat, tolerance = 1e-12)
  expect_equal(attr(fit_complete, "diagnostics")$n_complete, matrix(sum(keep_all), 3, 3, dimnames = dimnames(fit_complete)))
  expect_identical(attr(fit_complete, "diagnostics")$na_method, "complete")

  fit_pairwise <- ccc(X, na_method = "pairwise")
  diag_n <- attr(fit_pairwise, "diagnostics")$n_complete
  expect_equal(diag_n["a", "b"], sum(is.finite(X[, "a"]) & is.finite(X[, "b"])))
  expect_equal(diag_n["a", "c"], sum(is.finite(X[, "a"]) & is.finite(X[, "c"])))
  expect_equal(diag_n["b", "c"], sum(is.finite(X[, "b"]) & is.finite(X[, "c"])))

  ref_ab <- ccc(X[is.finite(X[, "a"]) & is.finite(X[, "b"]), c("a", "b"), drop = FALSE])
  expect_equal(fit_pairwise["a", "b"], ref_ab["a", "b"], tolerance = 1e-12)
})

test_that("ccc pairwise CIs and negative estimates are supported across output modes", {
  X <- cbind(
    x = c(1, 2, 3, 4, 5, NA),
    y = c(5.1, 3.8, 3.2, 2.1, 0.9, 0),
    z = c(1, 2, 1, 2, NA, 3)
  )

  fit_ci <- ccc(X, ci = TRUE, na_method = "pairwise")
  expect_s3_class(fit_ci, "ccc_ci")
  expect_identical(attr(fit_ci, "ci.method", exact = TRUE), "lin_delta_fisher_z")
  expect_true(is.finite(fit_ci$lwr.ci["x", "y"]))
  expect_lt(fit_ci$est["x", "y"], 0)

  edge <- ccc(X, na_method = "pairwise", output = "edge_list", threshold = 0.5)
  expect_s3_class(edge, "corr_edge_list")
  expect_true(any(edge$row == "x" & edge$col == "y" & edge$value < 0))

  edge_ci <- ccc(X, ci = TRUE, na_method = "pairwise", output = "edge_list", threshold = 0.5)
  expect_s3_class(edge_ci, "corr_edge_list")
  ci_attr <- attr(edge_ci, "ci", exact = TRUE)
  expect_true(is.matrix(ci_attr$lwr.ci))
  expect_identical(ci_attr$ci.method, "lin_delta_fisher_z")
  expect_true(any(edge_ci$row == "x" & edge_ci$col == "y" & edge_ci$value < 0))
})

