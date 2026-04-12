manual_spearman_jel_ci_R <- function(x, y, conf_level = 0.95) {
  keep <- stats::complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  est <- suppressWarnings(stats::cor(x, y, method = "spearman"))

  out <- list(estimate = est, lwr = NA_real_, upr = NA_real_, n = n)
  if (!is.finite(est) || n < 3L) return(out)

  z <- vapply(seq_len(n), function(i) {
    suppressWarnings(n * est - (n - 1) * stats::cor(x[-i], y[-i], method = "spearman"))
  }, numeric(1))
  if (any(!is.finite(z))) return(out)

  s_u <- mean((z - est)^2)
  if (!is.finite(s_u)) return(out)
  if (s_u <= 1e-14) {
    out$lwr <- est
    out$upr <- est
    return(out)
  }

  crit <- stats::qchisq(conf_level, df = 1)
  g <- function(theta) {
    s_theta <- mean((z - theta)^2)
    if (!(s_theta > 0)) {
      return(if ((est - theta)^2 <= 1e-18) -crit else Inf)
    }
    n * (est - theta)^2 / s_theta - crit
  }

  if (est <= -1) {
    out$lwr <- est
  } else if (!is.finite(g(-1)) || g(-1) <= 0) {
    out$lwr <- -1
  } else {
    out$lwr <- uniroot(g, interval = c(-1, est), tol = 1e-11)$root
  }

  if (est >= 1) {
    out$upr <- est
  } else if (!is.finite(g(1)) || g(1) <= 0) {
    out$upr <- 1
  } else {
    out$upr <- uniroot(g, interval = c(est, 1), tol = 1e-11)$root
  }

  out
}

test_that("spearman_rho matches base::cor(..., method = 'spearman')", {
  for (i in 1:10) {
    n <- sample(c(10, 50, 100, 500), 1)
    p <- sample(2:6, 1)
    mat <- replicate(p, rnorm(n))
    colnames(mat) <- paste0("S", seq_len(p))

    base_cor <- cor(mat, method = "spearman")
    fast_cor <- suppressWarnings(spearman_rho(mat))

    attributes(fast_cor) <- NULL
    attributes(base_cor) <- NULL

    expect_equal(
      fast_cor,
      base_cor,
      tolerance = 1e-8,
      info = paste("Mismatch on Spearman test dataset", i, "n =", n, "p =", p)
    )
  }
})

test_that("spearman_rho handles ties correctly and matches base::cor", {
  for (i in 1:5) {
    n <- sample(c(50, 100, 200), 1)
    p <- sample(2:5, 1)

    mat <- replicate(p, sample(rep(1:5, length.out = n)))
    colnames(mat) <- paste0("T", seq_len(p))

    base_cor <- cor(mat, method = "spearman")
    fast_cor <- suppressWarnings(spearman_rho(mat))

    attributes(fast_cor) <- NULL
    attributes(base_cor) <- NULL

    expect_equal(
      fast_cor,
      base_cor,
      tolerance = 1e-8,
      info = paste("Mismatch on tied Spearman dataset", i, "n =", n, "p =", p)
    )
  }
})

test_that("spearman_rho is invariant to strictly monotone transformations", {
  set.seed(123)
  X <- matrix(rnorm(200), ncol = 2)
  X_mono <- X
  X_mono[,1] <- exp(X_mono[,1])     # monotone transform
  X_mono[,2] <- log1p(exp(X_mono[,2]))  # monotone transform

  base_cor <- cor(X, method = "spearman")
  fast_cor <- spearman_rho(X_mono)

  expect_equal(as.numeric(fast_cor),
               as.numeric(unname(base_cor)),
               tolerance = 1e-12)
})

test_that("spearman_rho returns NA when a column is constant", {
  X <- cbind(a = rnorm(20), b = rep(1, 20))
  sp <- spearman_rho(X)

  expect_true(all(is.na(sp["b", ])))
  expect_true(all(is.na(sp[, "b"])))
})

test_that("spearman_rho agrees with base::cor on known values (small matrix)", {
  X <- matrix(c(1, 2, 3, 4,
                1, 4, 9, 16), ncol = 2)
  colnames(X) <- c("x", "y")

  base_cor <- cor(X, method = "spearman")
  fast_cor <- spearman_rho(X)

  expect_equal(as.numeric(fast_cor),
               as.numeric(unname(base_cor)),
               tolerance = 1e-12)
})

test_that("theoretical BVN formula for Spearman holds approximately", {
  set.seed(321)
  rhos <- c(-0.8, -0.4, 0, 0.4, 0.8)
  n <- 2000

  for (r in rhos) {
    Sigma <- matrix(c(1, r, r, 1), 2, 2)
    Z <- MASS::mvrnorm(n, mu = c(0,0), Sigma = Sigma)

    est <- spearman_rho(Z)[1,2]
    theory <- (6/pi) * asin(r/2)

    expect_equal(est, theory, tolerance = 0.08,
                 info = paste("Mismatch for true rho =", r))
  }
})

test_that("spearman_rho print/plot methods cover options", {
  skip_if_not_installed("ggplot2")

  set.seed(101)
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)
  colnames(X) <- paste0("S", seq_len(4))
  sp <- spearman_rho(X)

  out <- capture.output(print(sp, digits = 3, max_rows = 2, max_cols = 3))
  expect_true(any(grepl("omitted", out)))

  p <- plot(sp, title = "Spearman plot", low_color = "blue", high_color = "red", mid_color = "white", value_text_size = 3)
  expect_s3_class(p, "ggplot")
  scale <- p$scales$get_scales("fill")
  expect_equal(scale$limits, c(-1, 1))
})

test_that("spearman_rho rejects missing values when check_na = TRUE", {
  X <- cbind(a = c(1, 2, NA, 4), b = c(1, 2, 3, 4))
  expect_error(spearman_rho(X), "Missing values are not allowed.")
})

test_that("spearman_rho default CI behaviour remains estimate-only", {
  set.seed(2026)
  X <- matrix(rnorm(80), nrow = 20, ncol = 4)
  colnames(X) <- paste0("S", 1:4)

  sp <- spearman_rho(X)
  sm <- summary(sp)

  expect_s3_class(sm, "summary.matrixCorr")
  expect_s3_class(sm, "summary.corr_matrix")
  expect_null(attr(sp, "ci", exact = TRUE))
  expect_null(attr(sp, "conf.level", exact = TRUE))
})

test_that("spearman_rho honors n_threads without changing estimates", {
  set.seed(654)
  X <- matrix(rnorm(240), nrow = 40, ncol = 6)
  colnames(X) <- paste0("S", seq_len(ncol(X)))

  fit1 <- spearman_rho(X, n_threads = 1L)
  fit2 <- spearman_rho(X, n_threads = 2L)
  fit1_ci <- spearman_rho(X, ci = TRUE, n_threads = 1L)
  fit2_ci <- spearman_rho(X, ci = TRUE, n_threads = 2L)

  expect_equal(unclass(fit1), unclass(fit2), tolerance = 1e-12)
  expect_equal(attr(fit1_ci, "ci", exact = TRUE), attr(fit2_ci, "ci", exact = TRUE), tolerance = 1e-12)
})

test_that("spearman_rho CI matches reference JEL implementation for positive association", {
  x <- c(1, 2, 3, 4, 5, 6, 7, 8)
  y <- c(1.0, 1.8, 2.9, 4.1, 4.9, 6.2, 7.1, 8.0)
  ref <- manual_spearman_jel_ci_R(x, y, conf_level = 0.95)
  sp <- spearman_rho(cbind(x = x, y = y), ci = TRUE)
  ci <- attr(sp, "ci", exact = TRUE)

  expect_equal(unname(sp["x", "y"]), ref$estimate, tolerance = 1e-12)
  expect_equal(unname(ci$lwr.ci["x", "y"]), ref$lwr, tolerance = 1e-8)
  expect_equal(unname(ci$upr.ci["x", "y"]), ref$upr, tolerance = 1e-8)
})

test_that("spearman_rho CI matches reference JEL implementation for negative association", {
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
  y <- c(9, 7, 8, 6, 5, 4, 3, 2, 1)
  ref <- manual_spearman_jel_ci_R(x, y, conf_level = 0.95)
  sp <- spearman_rho(cbind(x = x, y = y), ci = TRUE)
  ci <- attr(sp, "ci", exact = TRUE)

  expect_equal(unname(sp["x", "y"]), ref$estimate, tolerance = 1e-12)
  expect_equal(unname(ci$lwr.ci["x", "y"]), ref$lwr, tolerance = 1e-8)
  expect_equal(unname(ci$upr.ci["x", "y"]), ref$upr, tolerance = 1e-8)
})

test_that("spearman_rho CI matches reference JEL implementation near zero", {
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  y <- c(3, 8, 1, 9, 4, 10, 2, 7, 5, 6)
  ref <- manual_spearman_jel_ci_R(x, y, conf_level = 0.95)
  sp <- spearman_rho(cbind(x = x, y = y), ci = TRUE)
  ci <- attr(sp, "ci", exact = TRUE)

  expect_equal(unname(sp["x", "y"]), ref$estimate, tolerance = 1e-12)
  expect_equal(unname(ci$lwr.ci["x", "y"]), ref$lwr, tolerance = 1e-8)
  expect_equal(unname(ci$upr.ci["x", "y"]), ref$upr, tolerance = 1e-8)
})

test_that("spearman_rho custom conf_level changes JEL interval width", {
  x <- 1:20
  y <- c(1, 3, 2, 4, 5, 7, 6, 8, 10, 9, 11, 13, 12, 14, 15, 17, 16, 18, 20, 19)

  sp90 <- spearman_rho(cbind(x = x, y = y), ci = TRUE, conf_level = 0.90)
  sp99 <- spearman_rho(cbind(x = x, y = y), ci = TRUE, conf_level = 0.99)

  ci90 <- attr(sp90, "ci", exact = TRUE)
  ci99 <- attr(sp99, "ci", exact = TRUE)
  width90 <- ci90$upr.ci["x", "y"] - ci90$lwr.ci["x", "y"]
  width99 <- ci99$upr.ci["x", "y"] - ci99$lwr.ci["x", "y"]

  expect_equal(attr(sp90, "conf.level", exact = TRUE), 0.90)
  expect_equal(attr(sp99, "conf.level", exact = TRUE), 0.99)
  expect_true(width99 > width90)
})

test_that("spearman_rho CI works for pairwise-complete input and reports n_complete", {
  X <- cbind(
    x = c(1, 2, 3, 4, 5, 6, NA, 8),
    y = c(2, 1, 4, 3, 6, 5, 7, NA),
    z = c(8, 7, 6, 5, 4, NA, 2, 1)
  )

  sp <- spearman_rho(X, na_method = "pairwise", ci = TRUE)
  diag_attr <- attr(sp, "diagnostics", exact = TRUE)
  ci <- attr(sp, "ci", exact = TRUE)
  ref_xy <- manual_spearman_jel_ci_R(X[, "x"], X[, "y"], conf_level = 0.95)

  expect_identical(diag_attr$n_complete["x", "y"], as.integer(ref_xy$n))
  expect_equal(unname(sp["x", "y"]), ref_xy$estimate, tolerance = 1e-12)
  expect_equal(unname(ci$lwr.ci["x", "y"]), ref_xy$lwr, tolerance = 1e-8)
  expect_equal(unname(ci$upr.ci["x", "y"]), ref_xy$upr, tolerance = 1e-8)
})

test_that("spearman_rho CI matches reference implementation with ties and small sample", {
  x <- c(1, 1, 2, 2, 3, 3)
  y <- c(1, 2, 2, 3, 3, 4)
  ref <- manual_spearman_jel_ci_R(x, y, conf_level = 0.95)
  sp <- spearman_rho(cbind(x = x, y = y), ci = TRUE)
  ci <- attr(sp, "ci", exact = TRUE)

  expect_equal(unname(sp["x", "y"]), ref$estimate, tolerance = 1e-12)
  expect_equal(unname(ci$lwr.ci["x", "y"]), ref$lwr, tolerance = 1e-8)
  expect_equal(unname(ci$upr.ci["x", "y"]), ref$upr, tolerance = 1e-8)
})

test_that("spearman_rho CI-aware print, summary, and plot follow the CI contract", {
  skip_if_not_installed("ggplot2")

  X <- cbind(
    A = c(1, 2, 3, 4, 5, 6, 7),
    B = c(1, 3, 2, 4, 6, 5, 7),
    C = c(7, 6, 5, 4, 3, 2, 1)
  )

  sp <- spearman_rho(X, ci = TRUE)
  sm <- summary(sp)
  txt_print <- capture.output(print(sp, show_ci = "yes"))
  txt_summary <- capture.output(print(sm))
  p <- plot(sp, ci_text_size = 2)

  expect_s3_class(sm, "summary.spearman_rho")
  expect_equal(nrow(sm), choose(ncol(X), 2))
  expect_true(all(c("item1", "item2", "estimate", "n_complete", "lwr", "upr") %in% names(sm)))
  expect_match(paste(txt_print, collapse = "\n"), "Spearman correlation matrix")
  expect_false(any(grepl("Spearman correlation summary", txt_print, fixed = TRUE)))
  expect_match(paste(txt_summary, collapse = "\n"), "Spearman correlation summary")
  expect_match(paste(txt_summary, collapse = "\n"), "jackknife_euclidean_likelihood")
  expect_match(paste(txt_summary, collapse = "\n"), "n_complete")
  expect_s3_class(p, "ggplot")
})
