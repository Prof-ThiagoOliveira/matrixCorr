manual_kendall_fieller_ci_R <- function(x, y, conf_level = 0.95) {
  keep <- stats::complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  tau <- suppressWarnings(stats::cor(x, y, method = "kendall"))

  out <- list(estimate = tau, lwr = NA_real_, upr = NA_real_, n = n)
  if (!is.finite(tau) || n <= 4L) return(out)

  z <- atanh(max(min(tau, 1 - 1e-15), -1 + 1e-15))
  crit <- stats::qnorm(0.5 * (1 + conf_level))
  se <- sqrt(0.437 / (n - 4))
  out$lwr <- max(-1, tanh(z - crit * se))
  out$upr <- min(1, tanh(z + crit * se))
  out
}

manual_kendall_brown_benedetti_ci_R <- function(x, y, conf_level = 0.95) {
  keep <- stats::complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)

  out <- list(estimate = NA_real_, lwr = NA_real_, upr = NA_real_, n = n)
  if (n < 2L) return(out)

  tab <- table(x, y)
  if (nrow(tab) < 2L || ncol(tab) < 2L) return(out)

  pi_c <- matrix(0, nrow = nrow(tab), ncol = ncol(tab))
  pi_d <- matrix(0, nrow = nrow(tab), ncol = ncol(tab))

  for (i in seq_len(nrow(tab))) {
    i_lt <- if (i > 1L) seq_len(i - 1L) else integer(0)
    i_gt <- if (i < nrow(tab)) seq.int(i + 1L, nrow(tab)) else integer(0)
    for (j in seq_len(ncol(tab))) {
      j_lt <- if (j > 1L) seq_len(j - 1L) else integer(0)
      j_gt <- if (j < ncol(tab)) seq.int(j + 1L, ncol(tab)) else integer(0)
      pi_c[i, j] <- sum(tab[i_lt, j_lt, drop = FALSE]) +
        sum(tab[i_gt, j_gt, drop = FALSE])
      pi_d[i, j] <- sum(tab[i_lt, j_gt, drop = FALSE]) +
        sum(tab[i_gt, j_lt, drop = FALSE])
    }
  }

  C <- sum(pi_c * tab) / 2
  D <- sum(pi_d * tab) / 2
  n0 <- n * (n - 1) / 2
  ti <- rowSums(tab)
  uj <- colSums(tab)
  n1 <- sum(ti * (ti - 1) / 2)
  n2 <- sum(uj * (uj - 1) / 2)
  den <- sqrt((n0 - n1) * (n0 - n2))
  if (!(den > 0)) return(out)

  tau_b <- (C - D) / den
  pi <- tab / sum(tab)
  pdiff <- (pi_c - pi_d) / sum(tab)
  Pdiff <- 2 * (C - D) / sum(tab)^2
  rowsum <- rowSums(pi)
  colsum <- colSums(pi)
  rowmat <- matrix(rep(rowsum, ncol(tab)), ncol = ncol(tab))
  colmat <- matrix(rep(colsum, nrow(tab)), nrow = nrow(tab), byrow = TRUE)
  delta1 <- sqrt(1 - sum(rowsum^2))
  delta2 <- sqrt(1 - sum(colsum^2))
  if (!(delta1 > 0) || !(delta2 > 0)) {
    out$estimate <- tau_b
    return(out)
  }

  tauphi <- (2 * pdiff + Pdiff * colmat) * delta2 * delta1 +
    (Pdiff * rowmat * delta2) / delta1
  sigma2 <- ((sum(pi * tauphi^2) - sum(pi * tauphi)^2) /
    (delta1 * delta2)^4) / n
  if (sigma2 < .Machine$double.eps * 10) sigma2 <- 0
  if (!is.finite(sigma2) || sigma2 < 0) {
    out$estimate <- tau_b
    return(out)
  }

  crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  ci <- tau_b + crit * sqrt(sigma2) * c(-1, 1)
  out$estimate <- tau_b
  out$lwr <- max(ci[1], -1)
  out$upr <- min(ci[2], 1)
  out
}

manual_kendall_ifel_ci_R <- function(x, y, conf_level = 0.95) {
  keep <- stats::complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  tau <- suppressWarnings(stats::cor(x, y, method = "kendall"))

  out <- list(estimate = tau, lwr = NA_real_, upr = NA_real_, n = n)
  if (!is.finite(tau) || n < 3L) return(out)

  s_sum <- numeric(n)
  a_sum <- numeric(n)
  b_sum <- numeric(n)
  s_total <- 0
  a_total <- 0
  b_total <- 0

  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      sx <- sign(x[i] - x[j])
      sy <- sign(y[i] - y[j])
      s_ij <- sx * sy
      a_ij <- as.numeric(sx != 0)
      b_ij <- as.numeric(sy != 0)

      s_sum[i] <- s_sum[i] + s_ij
      s_sum[j] <- s_sum[j] + s_ij
      a_sum[i] <- a_sum[i] + a_ij
      a_sum[j] <- a_sum[j] + a_ij
      b_sum[i] <- b_sum[i] + b_ij
      b_sum[j] <- b_sum[j] + b_ij

      s_total <- s_total + s_ij
      a_total <- a_total + a_ij
      b_total <- b_total + b_ij
    }
  }

  if (!(a_total > 0) || !(b_total > 0)) return(out)

  n0 <- n * (n - 1) / 2
  tau_hat <- s_total / sqrt(a_total * b_total)
  s_bar <- s_total / n0
  a_bar <- a_total / n0
  b_bar <- b_total / n0
  if (!(a_bar > 0) || !(b_bar > 0)) return(out)

  den <- sqrt(a_bar * b_bar)
  s_i <- s_sum / (n - 1)
  a_i <- a_sum / (n - 1)
  b_i <- b_sum / (n - 1)
  psi <- 2 * ((s_i - s_bar) / den -
    0.5 * tau_hat * ((a_i - a_bar) / a_bar + (b_i - b_bar) / b_bar))
  w <- tau_hat + psi

  if (!all(is.finite(w))) return(out)
  spread <- mean((w - tau_hat)^2)
  if (!is.finite(spread)) return(out)
  if (spread <= 1e-14) {
    out$estimate <- tau_hat
    out$lwr <- tau_hat
    out$upr <- tau_hat
    return(out)
  }

  elr <- function(theta) {
    v <- w - theta
    if (max(abs(v)) <= 1e-12) return(0)
    if (!(min(v) < 0 && max(v) > 0)) return(Inf)

    lower <- max(-1 / v[v > 0]) + 1e-12
    upper <- min(-1 / v[v < 0]) - 1e-12
    g <- function(lambda) sum(v / (1 + lambda * v))
    lambda <- if (abs(sum(v)) <= 1e-12) {
      0
    } else {
      stats::uniroot(g, interval = c(lower, upper), tol = 1e-11)$root
    }
    den_l <- 1 + lambda * v
    if (any(den_l <= 0)) return(Inf)
    2 * sum(log(den_l))
  }

  crit <- stats::qchisq(conf_level, df = 1)
  f <- function(theta) elr(theta) - crit

  left_edge <- max(-1, min(w))
  right_edge <- min(1, max(w))

  if (tau_hat <= left_edge + 1e-10) {
    out$lwr <- left_edge
  } else {
    left_start <- min(tau_hat, left_edge + 1e-10)
    if (!is.finite(f(left_start)) || f(left_start) > 0) {
      out$lwr <- stats::uniroot(f, interval = c(left_start, tau_hat), tol = 1e-11)$root
    } else {
      out$lwr <- left_edge
    }
  }

  if (tau_hat >= right_edge - 1e-10) {
    out$upr <- right_edge
  } else {
    right_start <- max(tau_hat, right_edge - 1e-10)
    if (!is.finite(f(right_start)) || f(right_start) > 0) {
      out$upr <- stats::uniroot(f, interval = c(tau_hat, right_start), tol = 1e-11)$root
    } else {
      out$upr <- right_edge
    }
  }

  out$estimate <- tau_hat
  out
}

test_that("kendall_tau matches base::cor(..., method = 'kendall')", {
  for (i in 1:10) {
    n <- sample(c(10, 50, 100, 500), 1)
    p <- sample(2:6, 1)
    mat <- replicate(p, rnorm(n))
    colnames(mat) <- paste0("V", seq_len(p))

    base_cor <- cor(mat, method = "kendall")
    fast_cor <- suppressWarnings(kendall_tau(mat))

    attributes(fast_cor) <- NULL
    attributes(base_cor) <- NULL

    expect_equal(
      fast_cor,
      base_cor,
      tolerance = 1e-8,
      info = paste("Mismatch on test dataset", i, "n =", n, "p =", p)
    )
  }
})

test_that("kendall_tau handles ties correctly and matches base::cor", {
  for (i in 1:5) {
    n <- sample(c(50, 100, 200), 1)
    p <- sample(2:5, 1)

    mat <- replicate(p, sample(rep(1:5, length.out = n)))
    colnames(mat) <- paste0("T", seq_len(p))

    base_cor <- cor(mat, method = "kendall")
    fast_cor <- suppressWarnings(kendall_tau(mat))

    attributes(fast_cor) <- NULL
    attributes(base_cor) <- NULL

    expect_equal(
      fast_cor,
      base_cor,
      tolerance = 1e-8,
      info = paste("Mismatch on tied test dataset", i, "n =", n, "p =", p)
    )
  }
})

test_that("kendall_tau numerics match stats::cor(Kendall)", {
  set.seed(123)

  X <- matrix(rnorm(200), ncol = 2)
  X_mono <- cbind(exp(X[, 1]), log1p(exp(X[, 2])))

  base_cor <- unname(cor(X, method = "kendall"))

  fast_raw <- kendall_tau(X_mono)
  fast_cor <- unname(as.matrix(fast_raw))
  class(fast_cor) <- NULL
  dimnames(fast_cor) <- NULL

  expect_equal(dim(fast_cor), dim(base_cor))
  expect_true(isSymmetric(matrix(fast_cor, nrow = nrow(base_cor))))
  expect_equal(diag(fast_cor), rep(1, ncol(X)))

  expect_equal(as.numeric(fast_cor),
    as.numeric(unname(cor(X, method = "kendall"))),
    tolerance = 1e-12
  )

  x <- X[, 1]
  y <- X[, 2]
  base12 <- cor(x, y, method = "kendall")
  fast12 <- as.numeric(as.matrix(kendall_tau(cbind(x, y)))[1, 2])
  expect_equal(fast12, base12, tolerance = 1e-12)

  y_rev <- -y
  expect_equal(
    as.numeric(as.matrix(kendall_tau(cbind(x, y_rev)))[1, 2]),
    -fast12,
    tolerance = 1e-12
  )
})

test_that("kendall_tau honors n_threads without changing estimates", {
  set.seed(987)
  X <- matrix(rnorm(240), nrow = 40, ncol = 6)
  colnames(X) <- paste0("K", seq_len(ncol(X)))

  fit1 <- kendall_tau(X, n_threads = 1L)
  fit2 <- kendall_tau(X, n_threads = 2L)
  fit1_ci <- kendall_tau(X, ci = TRUE, n_threads = 1L)
  fit2_ci <- kendall_tau(X, ci = TRUE, n_threads = 2L)

  expect_equal(unclass(fit1), unclass(fit2), tolerance = 1e-12)
  expect_equal(attr(fit1_ci, "ci", exact = TRUE), attr(fit2_ci, "ci", exact = TRUE), tolerance = 1e-12)
})

test_that("kendall_tau two-vector mode returns a scalar matching base::cor", {
  set.seed(456)
  x <- rnorm(300)
  y <- -0.4 * x + rnorm(300)
  kt <- kendall_tau(x, y)
  expect_type(kt, "double")
  expect_length(kt, 1L)
  expect_equal(
    kt,
    cor(x, y, method = "kendall"),
    tolerance = 1e-12
  )
})

test_that("kendall_tau returns NA when a column is constant", {
  X <- cbind(a = rnorm(20), b = rep(1, 20))
  kt <- kendall_tau(X)

  expect_true(all(is.na(kt["b", 1])))
  expect_true(all(is.na(kt[1, "b"])))
})

test_that("kendall_tau matches base::cor on a known toy dataset", {
  X <- matrix(c(
    1, 2, 3, 4,
    4, 3, 2, 1
  ), ncol = 2)
  colnames(X) <- c("x", "y")

  base_cor <- cor(X, method = "kendall")
  fast_cor <- kendall_tau(X)

  expect_equal(as.numeric(fast_cor),
    as.numeric(unname(base_cor)),
    tolerance = 1e-12
  )
})

test_that("kendall_tau rejects NA input in two-column mode", {
  X <- cbind(a = c(1, 2, NA, 4), b = c(0, 1, 2, 3))
  expect_error(kendall_tau(X), "Missing values are not allowed.")
})

test_that("kendall_tau estimates agree with theoretical BVN relationship", {
  set.seed(321)
  rhos <- c(-0.8, -0.4, 0, 0.4, 0.8)
  n <- 2000

  for (r in rhos) {
    Sigma <- matrix(c(1, r, r, 1), 2, 2)
    Z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma)

    est <- kendall_tau(Z)[1, 2]
    theory <- (2 / pi) * asin(r)

    expect_equal(est, theory, tolerance = 0.08,
      info = paste("Mismatch for true rho =", r)
    )
  }
})

test_that("kendall_tau default CI behaviour remains estimate-only", {
  set.seed(2026)
  X <- matrix(rnorm(80), nrow = 20, ncol = 4)
  colnames(X) <- paste0("K", 1:4)

  kt <- kendall_tau(X)
  sm <- summary(kt)

  expect_s3_class(sm, "summary.matrixCorr")
  expect_s3_class(sm, "summary.corr_matrix")
  expect_null(attr(kt, "ci", exact = TRUE))
  expect_null(attr(kt, "conf.level", exact = TRUE))
  expect_null(attr(kt, "ci.method", exact = TRUE))
})

test_that("kendall_tau rejects CI requests in two-vector mode", {
  expect_error(
    kendall_tau(1:6, 6:1, ci = TRUE),
    "Confidence intervals are not available"
  )
})

test_that("kendall_tau Fieller CI matches reference calculation for positive association", {
  x <- c(1, 2, 3, 4, 5, 6, 7, 8)
  y <- c(1.1, 2.2, 2.8, 4.0, 5.1, 5.9, 7.0, 8.3)
  ref <- manual_kendall_fieller_ci_R(x, y, conf_level = 0.95)
  kt <- kendall_tau(cbind(x = x, y = y), ci = TRUE)
  ci <- attr(kt, "ci", exact = TRUE)

  expect_equal(unname(kt["x", "y"]), ref$estimate, tolerance = 1e-12)
  expect_equal(unname(ci$lwr.ci["x", "y"]), ref$lwr, tolerance = 1e-10)
  expect_equal(unname(ci$upr.ci["x", "y"]), ref$upr, tolerance = 1e-10)
  expect_identical(ci$ci.method, "fieller")
})

test_that("kendall_tau Fieller CI matches reference calculation for negative association", {
  x <- c(1, 2, 3, 4, 5, 6, 7, 8)
  y <- c(8, 7, 6, 5, 4, 2, 3, 1)
  ref <- manual_kendall_fieller_ci_R(x, y, conf_level = 0.95)
  kt <- kendall_tau(cbind(x = x, y = y), ci = TRUE)
  ci <- attr(kt, "ci", exact = TRUE)

  expect_equal(unname(kt["x", "y"]), ref$estimate, tolerance = 1e-12)
  expect_equal(unname(ci$lwr.ci["x", "y"]), ref$lwr, tolerance = 1e-10)
  expect_equal(unname(ci$upr.ci["x", "y"]), ref$upr, tolerance = 1e-10)
})

test_that("kendall_tau Fieller CI matches reference calculation near zero", {
  x <- 1:10
  y <- c(2, 7, 1, 10, 4, 9, 3, 8, 6, 5)
  ref <- manual_kendall_fieller_ci_R(x, y, conf_level = 0.95)
  kt <- kendall_tau(cbind(x = x, y = y), ci = TRUE)
  ci <- attr(kt, "ci", exact = TRUE)

  expect_equal(unname(kt["x", "y"]), ref$estimate, tolerance = 1e-12)
  expect_equal(unname(ci$lwr.ci["x", "y"]), ref$lwr, tolerance = 1e-10)
  expect_equal(unname(ci$upr.ci["x", "y"]), ref$upr, tolerance = 1e-10)
})

test_that("kendall_tau custom conf_level changes Fieller interval width", {
  x <- 1:20
  y <- c(1, 3, 2, 4, 5, 7, 6, 8, 10, 9, 11, 13, 12, 14, 15, 17, 16, 18, 20, 19)

  kt90 <- kendall_tau(cbind(x = x, y = y), ci = TRUE, conf_level = 0.90)
  kt99 <- kendall_tau(cbind(x = x, y = y), ci = TRUE, conf_level = 0.99)

  ci90 <- attr(kt90, "ci", exact = TRUE)
  ci99 <- attr(kt99, "ci", exact = TRUE)
  width90 <- ci90$upr.ci["x", "y"] - ci90$lwr.ci["x", "y"]
  width99 <- ci99$upr.ci["x", "y"] - ci99$lwr.ci["x", "y"]

  expect_equal(attr(kt90, "conf.level", exact = TRUE), 0.90)
  expect_equal(attr(kt99, "conf.level", exact = TRUE), 0.99)
  expect_true(width99 > width90)
})

test_that("kendall_tau Fieller CI handles ties and pairwise-complete sample sizes", {
  X <- cbind(
    x = c(1, 1, 2, 2, 3, 3, NA, 4),
    y = c(1, 2, 2, 3, 3, 4, 4, NA),
    z = c(4, 4, 3, 3, 2, NA, 1, 1)
  )

  kt <- kendall_tau(X, na_method = "pairwise", ci = TRUE)
  diag_attr <- attr(kt, "diagnostics", exact = TRUE)
  ci <- attr(kt, "ci", exact = TRUE)
  ref_xy <- manual_kendall_fieller_ci_R(X[, "x"], X[, "y"], conf_level = 0.95)

  expect_identical(diag_attr$n_complete["x", "y"], as.integer(ref_xy$n))
  expect_equal(unname(kt["x", "y"]), ref_xy$estimate, tolerance = 1e-12)
  expect_equal(unname(ci$lwr.ci["x", "y"]), ref_xy$lwr, tolerance = 1e-10)
  expect_equal(unname(ci$upr.ci["x", "y"]), ref_xy$upr, tolerance = 1e-10)
})

test_that("kendall_tau Brown-Benedetti CI matches reference calculation for tied positive association", {
  x <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
  y <- c(1, 2, 2, 3, 3, 4, 4, 5, 5, 5)
  ref <- manual_kendall_brown_benedetti_ci_R(x, y, conf_level = 0.95)
  kt <- kendall_tau(cbind(x = x, y = y), ci = TRUE, ci_method = "brown_benedetti")
  ci <- attr(kt, "ci", exact = TRUE)

  expect_equal(unname(kt["x", "y"]), ref$estimate, tolerance = 1e-10)
  expect_equal(unname(ci$lwr.ci["x", "y"]), ref$lwr, tolerance = 1e-10)
  expect_equal(unname(ci$upr.ci["x", "y"]), ref$upr, tolerance = 1e-10)
  expect_identical(ci$ci.method, "brown_benedetti")
  expect_identical(attr(kt, "ci.method", exact = TRUE), "brown_benedetti")
})

test_that("kendall_tau Brown-Benedetti CI matches reference calculation for tied negative association", {
  x <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
  y <- c(5, 5, 4, 4, 3, 3, 2, 2, 1, 1)
  ref <- manual_kendall_brown_benedetti_ci_R(x, y, conf_level = 0.95)
  kt <- kendall_tau(cbind(x = x, y = y), ci = TRUE, ci_method = "brown_benedetti")
  ci <- attr(kt, "ci", exact = TRUE)

  expect_equal(unname(kt["x", "y"]), ref$estimate, tolerance = 1e-10)
  expect_equal(unname(ci$lwr.ci["x", "y"]), ref$lwr, tolerance = 1e-10)
  expect_equal(unname(ci$upr.ci["x", "y"]), ref$upr, tolerance = 1e-10)
})

test_that("kendall_tau Brown-Benedetti custom conf_level changes interval width", {
  x <- c(rep(1:6, each = 2), 7, 7, 8, 8)
  y <- c(1, 1, 2, 3, 3, 4, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8)

  kt90 <- kendall_tau(cbind(x = x, y = y), ci = TRUE, conf_level = 0.90, ci_method = "brown_benedetti")
  kt99 <- kendall_tau(cbind(x = x, y = y), ci = TRUE, conf_level = 0.99, ci_method = "brown_benedetti")

  ci90 <- attr(kt90, "ci", exact = TRUE)
  ci99 <- attr(kt99, "ci", exact = TRUE)
  width90 <- ci90$upr.ci["x", "y"] - ci90$lwr.ci["x", "y"]
  width99 <- ci99$upr.ci["x", "y"] - ci99$lwr.ci["x", "y"]

  expect_equal(attr(kt90, "conf.level", exact = TRUE), 0.90)
  expect_equal(attr(kt99, "conf.level", exact = TRUE), 0.99)
  expect_true(width99 > width90)
})

test_that("kendall_tau IF-EL CI matches reference implementation for positive association", {
  x <- 1:8
  y <- c(1.0, 2.2, 3.1, 4.0, 5.0, 5.7, 7.2, 8.0)
  ref <- manual_kendall_ifel_ci_R(x, y, conf_level = 0.95)
  kt <- kendall_tau(cbind(x = x, y = y), ci = TRUE, ci_method = "if_el")
  ci <- attr(kt, "ci", exact = TRUE)

  expect_equal(unname(kt["x", "y"]), ref$estimate, tolerance = 1e-12)
  expect_equal(unname(ci$lwr.ci["x", "y"]), ref$lwr, tolerance = 1e-8)
  expect_equal(unname(ci$upr.ci["x", "y"]), ref$upr, tolerance = 1e-8)
  expect_identical(ci$ci.method, "if_el")
  expect_identical(attr(kt, "ci.method", exact = TRUE), "if_el")
})

test_that("kendall_tau IF-EL CI matches reference implementation for negative association", {
  x <- 1:9
  y <- c(9, 8, 7, 5, 6, 4, 3, 2, 1)
  ref <- manual_kendall_ifel_ci_R(x, y, conf_level = 0.95)
  kt <- kendall_tau(cbind(x = x, y = y), ci = TRUE, ci_method = "if_el")
  ci <- attr(kt, "ci", exact = TRUE)

  expect_equal(unname(kt["x", "y"]), ref$estimate, tolerance = 1e-12)
  expect_equal(unname(ci$lwr.ci["x", "y"]), ref$lwr, tolerance = 1e-8)
  expect_equal(unname(ci$upr.ci["x", "y"]), ref$upr, tolerance = 1e-8)
})

test_that("kendall_tau IF-EL CI matches reference implementation near zero", {
  x <- 1:10
  y <- c(3, 9, 1, 10, 5, 7, 4, 8, 2, 6)
  ref <- manual_kendall_ifel_ci_R(x, y, conf_level = 0.95)
  kt <- kendall_tau(cbind(x = x, y = y), ci = TRUE, ci_method = "if_el")
  ci <- attr(kt, "ci", exact = TRUE)

  expect_equal(unname(kt["x", "y"]), ref$estimate, tolerance = 1e-12)
  expect_equal(unname(ci$lwr.ci["x", "y"]), ref$lwr, tolerance = 1e-8)
  expect_equal(unname(ci$upr.ci["x", "y"]), ref$upr, tolerance = 1e-8)
})

test_that("kendall_tau IF-EL CI handles ties and small samples", {
  x <- c(1, 1, 2, 2, 3, 3)
  y <- c(1, 2, 2, 3, 3, 4)
  ref <- manual_kendall_ifel_ci_R(x, y, conf_level = 0.95)
  kt <- kendall_tau(cbind(x = x, y = y), ci = TRUE, ci_method = "if_el")
  ci <- attr(kt, "ci", exact = TRUE)

  expect_equal(unname(kt["x", "y"]), ref$estimate, tolerance = 1e-12)
  expect_equal(unname(ci$lwr.ci["x", "y"]), ref$lwr, tolerance = 1e-8)
  expect_equal(unname(ci$upr.ci["x", "y"]), ref$upr, tolerance = 1e-8)
})

test_that("kendall_tau IF-EL returns NA bounds when the interval cannot be formed", {
  X <- cbind(
    x = c(1, 2),
    y = c(2, 1)
  )

  kt <- kendall_tau(X, ci = TRUE, ci_method = "if_el")
  ci <- attr(kt, "ci", exact = TRUE)

  expect_true(is.na(ci$lwr.ci["x", "y"]))
  expect_true(is.na(ci$upr.ci["x", "y"]))
})

test_that("kendall_tau CI-aware print, summary, and plot follow the CI contract", {
  skip_if_not_installed("ggplot2")

  X <- cbind(
    A = c(1, 2, 3, 4, 5, 6, 7),
    B = c(1, 3, 2, 4, 6, 5, 7),
    C = c(7, 6, 5, 4, 3, 2, 1)
  )

  kt <- kendall_tau(X, ci = TRUE, ci_method = "fieller")
  sm <- summary(kt)
  txt_print <- capture.output(print(kt, show_ci = "yes"))
  txt_summary <- capture.output(print(sm))
  p <- plot(kt, ci_text_size = 2)

  expect_s3_class(sm, "summary.kendall_matrix")
  expect_equal(nrow(sm), choose(ncol(X), 2))
  expect_true(all(c("item1", "item2", "estimate", "n_complete", "lwr", "upr") %in% names(sm)))
  expect_match(paste(txt_print, collapse = "\n"), "Kendall correlation matrix")
  expect_false(any(grepl("Kendall correlation summary", txt_print, fixed = TRUE)))
  expect_match(paste(txt_summary, collapse = "\n"), "Kendall correlation summary")
  expect_match(paste(txt_summary, collapse = "\n"), "fieller")
  expect_match(paste(txt_summary, collapse = "\n"), "n_complete")
  expect_s3_class(p, "ggplot")
})

test_that("kendall_tau CI-enabled objects keep the same outer contract as other methods", {
  X <- cbind(
    A = c(1, 2, 3, 4, 5, 6),
    B = c(1, 2, 4, 3, 5, 6),
    C = c(6, 5, 4, 3, 2, 1)
  )

  k_fieller <- kendall_tau(X, ci = TRUE, ci_method = "fieller")
  k_bb <- kendall_tau(X, ci = TRUE, ci_method = "brown_benedetti")
  k_ifel <- kendall_tau(X, ci = TRUE, ci_method = "if_el")
  p <- pearson_corr(X, ci = TRUE)
  s <- spearman_rho(X, ci = TRUE)

  expect_true(all(c("est", "lwr.ci", "upr.ci", "conf.level") %in% names(attr(k_fieller, "ci", exact = TRUE))))
  expect_true(all(c("est", "lwr.ci", "upr.ci", "conf.level") %in% names(attr(k_bb, "ci", exact = TRUE))))
  expect_true(all(c("est", "lwr.ci", "upr.ci", "conf.level") %in% names(attr(p, "ci", exact = TRUE))))
  expect_true(all(c("est", "lwr.ci", "upr.ci", "conf.level") %in% names(attr(s, "ci", exact = TRUE))))
  expect_equal(names(summary(k_fieller)), names(summary(k_bb)))
  expect_equal(names(summary(k_fieller)), names(summary(k_ifel)))
})
