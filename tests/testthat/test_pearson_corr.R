test_that("pearson_corr matches base::cor(..., method = 'pearson')", {
  for (i in 1:10) {
    n <- sample(c(10, 50, 100, 500), 1)
    p <- sample(2:6, 1)
    mat <- replicate(p, rnorm(n))
    colnames(mat) <- paste0("P", seq_len(p))

    base_cor <- cor(mat, method = "pearson")
    fast_cor <- suppressWarnings(pearson_corr(mat))

    attributes(fast_cor) <- NULL
    attributes(base_cor) <- NULL

    expect_equal(
      fast_cor,
      base_cor,
      tolerance = 1e-10,
      info = paste("Mismatch on Pearson test dataset", i, "n =", n, "p =", p)
    )
  }
})

test_that("pearson_corr matches stats::cor on transformed data", {
  set.seed(123)
  X <- matrix(rnorm(200), ncol = 2)

  X_affine <- X
  X_affine[,1] <-  3 * X_affine[,1] + 5
  X_affine[,2] <- -2 * X_affine[,2] + 7

  base <- unname(cor(X_affine, method = "pearson"))
  fast <- unname(as.matrix(unclass(pearson_corr(X_affine))))

  expect_equal(as.numeric(fast),
               as.numeric(unname(base)),
               tolerance = 1e-12)
})

test_that("pearson_corr returns NA when a column is constant", {
  X <- cbind(a = rnorm(20), b = rep(1, 20))
  pc <- pearson_corr(X)

  expect_true(all(is.na(pc["b", ])))
  expect_true(all(is.na(pc[, "b"])))
})

test_that("pearson_corr matches base::cor on a known toy dataset", {
  X <- matrix(c(1, 2, 3, 4,
                2, 4, 6, 8), ncol = 2)
  colnames(X) <- c("x", "y")

  base_cor <- cor(X, method = "pearson")
  fast_cor <- pearson_corr(X)
  expect_equal(as.numeric(fast_cor),
               as.numeric(unname(base_cor)),
               tolerance = 1e-12)
})

test_that("pearson_corr recovers true correlation in large MVN sample", {
  set.seed(42)
  p <- 2; n <- 5000; rho <- 0.7
  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  Z <- MASS::mvrnorm(n, mu = c(0,0), Sigma = Sigma)

  est <- pearson_corr(Z)[1,2]

  expect_equal(est, rho, tolerance = 0.02,
               info = paste("Estimated =", round(est,3),
                            "Expected =", rho))
})

test_that("pearson_corr print/plot methods cover optional arguments", {
  skip_if_not_installed("ggplot2")

  set.seed(202)
  X <- matrix(rnorm(60), nrow = 15, ncol = 4)
  colnames(X) <- paste0("P", seq_len(4))
  pc <- pearson_corr(X)

  out <- capture.output(print(pc, digits = 3, max_rows = 2, max_cols = 3))
  expect_true(any(grepl("omitted", out)))

  p <- plot(pc, title = "Pearson plot", low_color = "blue", high_color = "red", mid_color = "white", value_text_size = 3)
  expect_s3_class(p, "ggplot")
  scale <- p$scales$get_scales("fill")
  expect_equal(scale$limits, c(-1, 1))
})

test_that("pearson_corr rejects missing values by default", {
  X <- cbind(a = c(1, 2, NA, 4), b = c(1, 2, 3, 4))
  expect_error(pearson_corr(X), "Missing values are not allowed.")
})

test_that("pearson_corr pairwise mode uses complete-case overlap per pair", {
  X <- cbind(
    x = c(1, 2, 3, 4, NA, 6),
    y = c(2, 4, 6, 8, 10, NA),
    z = c(1, NA, 2, 3, 4, 5)
  )

  R <- pearson_corr(X, check_na = FALSE)
  diag_attr <- attr(R, "diagnostics", exact = TRUE)

  keep_xy <- stats::complete.cases(X[, "x"], X[, "y"])
  keep_xz <- stats::complete.cases(X[, "x"], X[, "z"])

  expect_equal(unname(R["x", "y"]), stats::cor(X[keep_xy, "x"], X[keep_xy, "y"]))
  expect_equal(unname(R["x", "z"]), stats::cor(X[keep_xz, "x"], X[keep_xz, "z"]))
  expect_identical(diag_attr$n_complete["x", "y"], as.integer(sum(keep_xy)))
  expect_identical(diag_attr$n_complete["x", "z"], as.integer(sum(keep_xz)))
})

test_that("pearson_corr Fisher-z confidence intervals match manual calculation", {
  X <- cbind(
    x = c(1, 2, 3, 4, 5, 6, NA, 8),
    y = c(1.4, 2.1, 3.2, 4.0, 5.4, 5.9, 7.1, NA),
    z = c(8, 7, 6, 5, 4, 3, 2, 1)
  )

  R <- pearson_corr(X, check_na = FALSE, ci = TRUE)
  ci <- attr(R, "ci", exact = TRUE)
  diag_attr <- attr(R, "diagnostics", exact = TRUE)

  keep_xy <- stats::complete.cases(X[, "x"], X[, "y"])
  x_xy <- X[keep_xy, "x"]
  y_xy <- X[keep_xy, "y"]
  r_xy <- stats::cor(x_xy, y_xy)
  n_xy <- sum(keep_xy)
  z_xy <- atanh(r_xy)
  crit <- stats::qnorm(0.975)
  se_xy <- 1 / sqrt(n_xy - 3)
  lwr_xy <- tanh(z_xy - crit * se_xy)
  upr_xy <- tanh(z_xy + crit * se_xy)

  expect_equal(unname(R["x", "y"]), r_xy, tolerance = 1e-12)
  expect_identical(diag_attr$n_complete["x", "y"], as.integer(n_xy))
  expect_equal(unname(ci$lwr.ci["x", "y"]), lwr_xy, tolerance = 1e-12)
  expect_equal(unname(ci$upr.ci["x", "y"]), upr_xy, tolerance = 1e-12)
  expect_equal(ci$conf.level, 0.95)
})

test_that("pearson_corr CI-aware print, summary, and plot follow the CI style", {
  skip_if_not_installed("ggplot2")

  X <- cbind(
    A = c(1, 2, 3, 4, 5, 6),
    B = c(1.2, 2.1, 3.2, 3.9, 5.1, 5.8),
    C = c(6, 5, 4, 3, 2, 1)
  )

  R <- pearson_corr(X, ci = TRUE)
  sm <- summary(R)
  txt_print <- capture.output(print(R, show_ci = "yes"))
  txt_summary <- capture.output(print(sm))
  p <- plot(R, ci_text_size = 2)

  expect_s3_class(sm, "summary.pearson_corr")
  expect_equal(nrow(sm), choose(ncol(X), 2))
  expect_true(all(c("var1", "var2", "estimate", "n_complete", "lwr", "upr") %in% names(sm)))
  expect_match(paste(txt_print, collapse = "\n"), "Pearson correlation matrix")
  expect_false(any(grepl("Pearson correlation summary", txt_print, fixed = TRUE)))
  expect_match(paste(txt_summary, collapse = "\n"), "Pearson correlation summary")
  expect_match(paste(txt_summary, collapse = "\n"), "ci_width")
  expect_match(paste(txt_summary, collapse = "\n"), "n_complete")
  expect_s3_class(p, "ggplot")
})

