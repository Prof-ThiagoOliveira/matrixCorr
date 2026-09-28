local_v_dcor <- function(x, y, squared = FALSE) {
  ax <- abs(outer(x, x, "-"))
  ay <- abs(outer(y, y, "-"))
  Ax <- sweep(sweep(ax, 1L, rowMeans(ax), "-"), 2L, colMeans(ax), "-") + mean(ax)
  Ay <- sweep(sweep(ay, 1L, rowMeans(ay), "-"), 2L, colMeans(ay), "-") + mean(ay)
  xy <- mean(Ax * Ay)
  x2 <- mean(Ax * Ax)
  y2 <- mean(Ay * Ay)
  r2 <- xy / sqrt(x2 * y2)
  if (!is.finite(r2)) {
    return(NA_real_)
  }
  r2 <- min(1, max(0, r2))
  if (isTRUE(squared)) r2 else sqrt(r2)
}

local_u_center <- function(x) {
  a <- abs(outer(x, x, "-"))
  n <- nrow(a)
  out <- matrix(0, n, n)
  off <- row(a) != col(a)
  row_sum <- rowSums(a)
  col_sum <- colSums(a)
  total <- sum(a)
  out[off] <- a[off] -
    row_sum[row(a)[off]] / (n - 2) -
    col_sum[col(a)[off]] / (n - 2) +
    total / ((n - 1) * (n - 2))
  out
}

local_bcdcor <- function(x, y) {
  Ax <- local_u_center(x)
  Ay <- local_u_center(y)
  xy <- sum(Ax * Ay)
  x2 <- sum(Ax * Ax)
  y2 <- sum(Ay * Ay)
  xy / sqrt(x2 * y2)
}

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
  expect_lt(res["x","z"], 0.7)  # finite-sample conventional dCor for this fixed z
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

test_that("bcdcor inference matches the signed reference t-test", {
  set.seed(2026)
  x <- rnorm(80)
  y <- x^2 + rnorm(80, sd = 0.25)
  z <- -0.4 * x + rnorm(80, sd = 0.8)
  X <- cbind(x = x, y = y, z = z)

  bc <- bcdcor(X, p_value = TRUE)
  inf <- attr(bc, "inference", exact = TRUE)
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
  expect_identical(inf$method, "bcdcor_t_test")
  expect_true(is.matrix(inf$bcdcor))
  expect_true(is.matrix(inf$statistic))
  expect_true(is.matrix(inf$parameter))
  expect_true(is.matrix(inf$p_value))

  idx <- which(upper.tri(bc), arr.ind = TRUE)
  for (k in seq_len(nrow(idx))) {
    i <- idx[k, 1]
    j <- idx[k, 2]
    ref_est <- ref_estimate[i, j]
    ref_t <- ref_statistic[i, j]
    ref_df_ij <- ref_df[i, j]
    ref_p <- ref_p_value[i, j]

    expect_lt(abs(bc[i, j] - ref_est), 1e-6)
    expect_lt(abs(inf$bcdcor[i, j] - ref_est), 1e-6)
    expect_lt(abs(inf$statistic[i, j] - ref_t), 1e-5)
    expect_equal(inf$parameter[i, j], ref_df_ij, tolerance = 1e-10)
    if (isTRUE(is.finite(ref_p) && ref_p == 0)) {
      expect_gte(inf$p_value[i, j], 0)
      expect_lt(inf$p_value[i, j], .Machine$double.eps)
    } else {
      expect_lt(abs(inf$p_value[i, j] - ref_p), 1e-10)
    }
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

test_that("dcor stores inference payload only when requested", {
  X <- cbind(
    a = c(1, 2, 3, 4, 5, 6),
    b = c(2, 1, 4, 3, 6, 5),
    c = c(1, 4, 2, 5, 3, 6)
  )

  fit <- dcor(X)
  fit_p <- dcor(X, p_value = TRUE)

  expect_null(attr(fit, "inference", exact = TRUE))
  expect_null(attr(fit, "diagnostics", exact = TRUE))

  inf <- attr(fit_p, "inference", exact = TRUE)
  expect_type(inf, "list")
  expect_identical(names(inf), c("method", "bcdcor", "statistic", "parameter", "p_value", "alternative"))
  expect_true(is.matrix(inf$bcdcor))
  expect_true(is.matrix(inf$statistic))
  expect_true(is.matrix(inf$parameter))
  expect_true(is.matrix(inf$p_value))
})

test_that("dcor honors n_threads without changing estimates", {
  set.seed(246)
  X <- matrix(rnorm(240), nrow = 40, ncol = 6)
  colnames(X) <- paste0("D", seq_len(ncol(X)))

  fit1 <- dcor(X, n_threads = 1L)
  fit2 <- dcor(X, n_threads = 2L)
  fit1_p <- dcor(X, p_value = TRUE, n_threads = 1L)
  fit2_p <- dcor(X, p_value = TRUE, n_threads = 2L)

  expect_equal(unclass(fit1), unclass(fit2), tolerance = 1e-12)
  expect_equal(attr(fit1_p, "inference", exact = TRUE), attr(fit2_p, "inference", exact = TRUE), tolerance = 1e-12)
})

test_that("dcor summary switches to pairwise inference view when requested", {
  set.seed(11)
  X <- matrix(rnorm(120), nrow = 30, ncol = 4)
  colnames(X) <- paste0("D", seq_len(ncol(X)))

  dc_no_p <- dcor(X)
  sm_no_p <- summary(dc_no_p)
  expect_false(isTRUE(attr(sm_no_p, "has_p", exact = TRUE)))
  expect_false("p_value" %in% names(sm_no_p))

  dc <- dcor(X, p_value = TRUE)
  sm <- summary(dc)

  expect_s3_class(sm, "summary.dcor")
  expect_s3_class(sm, "data.frame")
  expect_true(all(c("item1", "item2", "estimate", "n_complete", "bcdcor", "statistic", "df", "p_value") %in% names(sm)))

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
  expect_error(bcdcor(X), "Missing values are not allowed.")
})

test_that("dcor matches a conventional V-statistic reference and squared identity", {
  set.seed(100)
  n <- 400
  rho <- 0.7
  x <- rnorm(n)
  y <- rho * x + sqrt(1 - rho^2) * rnorm(n)
  X <- cbind(x = x, y = y)

  mc_R <- dcor(X)["x", "y"]
  mc_R2 <- dcor(X, squared = TRUE)["x", "y"]
  mc_bc <- bcdcor(X)["x", "y"]
  ref_R <- local_v_dcor(x, y)
  ref_R2 <- local_v_dcor(x, y, squared = TRUE)
  ref_bc <- local_bcdcor(x, y)

  expect_equal(mc_R, ref_R, tolerance = 1e-10)
  expect_equal(mc_R2, ref_R2, tolerance = 1e-10)
  expect_equal(unname(mc_bc), unname(ref_bc), tolerance = 1e-10)
  expect_equal(mc_R2, mc_R^2, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(mc_bc, mc_R^2, tolerance = 1e-4)))
})

test_that("dcor and bcdcor match self-contained references across Gaussian correlations", {
  n <- 350
  for (rho in c(0, 0.2, 0.4, 0.6, 0.8)) {
    set.seed(900 + as.integer(100 * rho))
    x <- rnorm(n)
    y <- rho * x + sqrt(1 - rho^2) * rnorm(n)
    X <- cbind(x = x, y = y)

    mc_R <- dcor(X)["x", "y"]
    mc_R2 <- dcor(X, squared = TRUE)["x", "y"]
    mc_bc <- bcdcor(X)["x", "y"]
    ref_R <- local_v_dcor(x, y)
    ref_R2 <- local_v_dcor(x, y, squared = TRUE)
    ref_bc <- local_bcdcor(x, y)

    expect_equal(mc_R, ref_R, tolerance = 1e-10, info = paste("rho", rho))
    expect_equal(mc_R2, ref_R2, tolerance = 1e-10, info = paste("rho", rho))
    expect_equal(mc_R2, mc_R^2, tolerance = 1e-12, info = paste("rho", rho))
    expect_equal(unname(mc_bc), unname(ref_bc), tolerance = 1e-10, info = paste("rho", rho))
  }
})

test_that("independent data can produce negative signed bcdcor", {
  set.seed(3)
  x <- rnorm(80)
  y <- rnorm(80)
  X <- cbind(x = x, y = y)

  mc_R <- dcor(X)["x", "y"]
  mc_bc <- bcdcor(X)["x", "y"]
  ref_bc <- local_bcdcor(x, y)

  expect_gte(mc_R, -1e-12)
  expect_lte(mc_R, 1 + 1e-12)
  expect_lt(mc_bc, 0)
  expect_equal(unname(mc_bc), unname(ref_bc), tolerance = 1e-10)
})

test_that("perfect and negative linear dependence remain nonnegative for dcor", {
  x <- seq(-2, 2, length.out = 80)

  pos <- cbind(x = x, y = 2 * x + 1)
  neg <- cbind(x = x, y = -2 * x)

  expect_equal(dcor(pos)["x", "y"], 1, tolerance = 1e-12)
  expect_equal(dcor(pos, squared = TRUE)["x", "y"], 1, tolerance = 1e-12)
  expect_equal(bcdcor(pos)["x", "y"], 1, tolerance = 1e-12)
  expect_equal(dcor(neg)["x", "y"], 1, tolerance = 1e-12)
  expect_equal(dcor(neg, squared = TRUE)["x", "y"], 1, tolerance = 1e-12)
})

test_that("dcor is bounded and squared entries equal squared dcor", {
  set.seed(44)
  X <- matrix(rnorm(300), nrow = 100, ncol = 3)
  colnames(X) <- letters[1:3]

  R <- dcor(X)
  R2 <- dcor(X, squared = TRUE)
  vals <- R[upper.tri(R)]

  expect_true(all(vals >= -1e-12, na.rm = TRUE))
  expect_true(all(vals <= 1 + 1e-12, na.rm = TRUE))
  expect_equal(R2[upper.tri(R2)], R[upper.tri(R)]^2, tolerance = 1e-12)
})

test_that("dcor and bcdcor support missing-data modes", {
  X <- cbind(
    a = c(1, 2, NA, 4, 5, 6, 7),
    b = c(2, 4, 6, 8, NA, 12, 14),
    c = c(7, 6, 5, 4, 3, 2, 1)
  )
  Xcc <- X[stats::complete.cases(X) & apply(is.finite(X), 1L, all), , drop = FALSE]

  expect_error(dcor(X, na_method = "error"), "Missing values are not allowed.")
  expect_error(bcdcor(X, na_method = "error"), "Missing values are not allowed.")

  dcor_complete <- unclass(dcor(X, na_method = "complete"))
  dcor_manual <- unclass(dcor(Xcc, na_method = "error"))
  attributes(dcor_complete) <- attributes(dcor_complete)[c("dim", "dimnames")]
  attributes(dcor_manual) <- attributes(dcor_manual)[c("dim", "dimnames")]
  expect_equal(
    dcor_complete,
    dcor_manual,
    tolerance = 1e-12
  )
  bcdcor_complete <- unclass(bcdcor(X, na_method = "complete"))
  bcdcor_manual <- unclass(bcdcor(Xcc, na_method = "error"))
  attributes(bcdcor_complete) <- attributes(bcdcor_complete)[c("dim", "dimnames")]
  attributes(bcdcor_manual) <- attributes(bcdcor_manual)[c("dim", "dimnames")]
  expect_equal(
    bcdcor_complete,
    bcdcor_manual,
    tolerance = 1e-12
  )

  expect_s3_class(dcor(X, na_method = "pairwise"), "dcor")
  expect_s3_class(bcdcor(X, na_method = "pairwise"), "bcdcor")
})

test_that("bcdcor output modes preserve signed threshold semantics", {
  set.seed(3)
  X <- cbind(x = rnorm(80), y = rnorm(80), z = rnorm(80))
  dense <- bcdcor(X)
  sparse <- bcdcor(X, output = "sparse", threshold = 0.001, diag = FALSE)
  edges <- bcdcor(X, output = "edge_list", threshold = 0.001, diag = FALSE)

  expect_s3_class(dense, "bcdcor")
  expect_true(isTRUE(attr(sparse, "corr_result")))
  expect_s3_class(edges, "corr_edge_list")
  expect_true(any(edges$value < 0))
  expect_true(all(abs(edges$value) >= 0.001))
})
