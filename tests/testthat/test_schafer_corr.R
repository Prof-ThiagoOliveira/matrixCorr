test_that("schafer_corr returns shrinkage matrix with metadata", {
  set.seed(123)
  X <- matrix(rnorm(80), nrow = 20, ncol = 4)
  colnames(X) <- paste0("V", seq_len(4))

  R <- schafer_corr(X)

  expect_s3_class(R, "schafer_corr")
  expect_equal(dim(R), c(4, 4))
  expect_true(all(diag(R) == 1))
  expect_equal(attr(R, "method"), "schafer_shrinkage")
  expect_equal(attr(R, "package"), "matrixCorr")

  # shrinkage pulls towards zero relative to sample correlation
  R_raw <- stats::cor(X)
  off <- upper.tri(R_raw, diag = FALSE)
  expect_true(mean(abs(R[off])) <= mean(abs(R_raw[off])) + 1e-10)
})

test_that("schafer_corr flags zero-variance columns as NA", {
  X <- cbind(const = rep(2, 6), var = rnorm(6))

  R <- schafer_corr(X)

  expect_true(all(is.na(R["const", ])))
  expect_true(all(is.na(R[, "const"])))
})

test_that("print and plot methods cover optional arguments", {
  skip_if_not_installed("ggplot2")

  set.seed(456)
  X <- matrix(rnorm(60), nrow = 15, ncol = 4)
  colnames(X) <- paste0("G", seq_len(4))
  R <- schafer_corr(X)

  # print truncation path
  out <- capture.output(print(R, digits = 3, max_rows = 2, max_cols = 3))
  expect_true(any(grepl("omitted", out)))

  p1 <- plot(R, cluster = TRUE, triangle = "lower", show_values = TRUE, value_text_limit = 10)
  expect_s3_class(p1, "ggplot")

  if (requireNamespace("viridisLite", quietly = TRUE)) {
    p2 <- plot(R, cluster = FALSE, triangle = "upper", palette = "viridis")
    expect_s3_class(p2, "ggplot")
  }
})
