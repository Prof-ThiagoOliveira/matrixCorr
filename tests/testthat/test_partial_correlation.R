test_that("partial_correlation returns expected components for each method", {
  set.seed(111)
  X <- matrix(rnorm(200), nrow = 40, ncol = 5)
  colnames(X) <- paste0("V", seq_len(5))

  samp <- partial_correlation(X, method = "sample", return_cov_precision = TRUE)
  expect_s3_class(samp, "partial_corr")
  expect_equal(dim(samp$pcor), c(5, 5))
  expect_true(all(diag(samp$pcor) == 1))
  expect_true(is.matrix(samp$cov))
  expect_true(is.matrix(samp$precision))
  expect_equal(samp$method, "sample")
  expect_true(is.na(samp$lambda))
  expect_true(is.na(samp$rho))

  ridge <- partial_correlation(X, method = "ridge", lambda = 1e-2, return_cov_precision = TRUE)
  expect_equal(ridge$method, "ridge")
  expect_equal(ridge$lambda, 1e-2)
  expect_true(is.na(ridge$rho))

  oas <- partial_correlation(X, method = "oas", return_cov_precision = TRUE)
  expect_equal(oas$method, "oas")
  expect_true(is.na(oas$lambda))
  expect_true(oas$rho >= 0 && oas$rho <= 1)
})

test_that("partial_correlation print and plot methods cover options", {
  skip_if_not_installed("ggplot2")

  set.seed(222)
  X <- matrix(rnorm(150), nrow = 30, ncol = 5)
  colnames(X) <- paste0("G", seq_len(5))
  pc <- partial_correlation(X, method = "ridge", lambda = 5e-3, return_cov_precision = TRUE)

  out1 <- capture.output(print(pc, digits = 4, max_rows = 3, max_cols = 4))
  expect_true(any(grepl("omitted", out1)))

  out2 <- capture.output(print(pc, show_method = FALSE, max_rows = 2, max_cols = 2))
  expect_true(any(grepl("Partial correlation matrix", out2)))

  p <- plot(pc, mask_diag = FALSE, reorder = TRUE, value_text_size = 3,
            low_color = "blue", high_color = "red", mid_color = "white")
  expect_s3_class(p, "ggplot")
})

test_that("partial_correlation validates lambda", {
  set.seed(1)
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)
  expect_error(partial_correlation(X, method = "ridge", lambda = -1), "must be >=")
})

test_that("partial_correlation exposes shrinkage metadata without cov/precision", {
  set.seed(333)
  X <- matrix(rnorm(120), nrow = 30, ncol = 4)

  oas <- partial_correlation(X, method = "oas")
  expect_true(is.numeric(oas$rho))
  expect_true(isTRUE(oas$rho >= 0 && oas$rho <= 1))
  expect_true(is.numeric(oas$jitter))
  expect_false(is.na(oas$jitter))
  expect_null(oas$cov)
  expect_null(oas$precision)

  samp <- partial_correlation(X, method = "sample")
  expect_true(is.na(samp$rho))
  expect_true(is.na(samp$lambda))
  expect_true(is.numeric(samp$jitter))
})
