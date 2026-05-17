hsic_gaussian_raw_R <- function(x, y) {
  med_bw <- function(z) {
    D <- abs(outer(z, z, "-"))
    d <- sort(as.numeric(D[upper.tri(D, diag = FALSE)]))
    d <- d[is.finite(d)]
    if (!length(d)) return(0.001)
    sigma <- d[floor(length(d) / 2) + 1L] / sqrt(2)
    if (!is.finite(sigma) || sigma <= 0) return(0.001)
    sigma
  }
  gram <- function(z) {
    bw <- med_bw(z)
    D <- outer(z, z, "-")
    exp(-(D^2) / (2 * bw^2))
  }
  n <- length(x)
  H <- diag(n) - matrix(1 / n, n, n)
  Kc <- H %*% gram(x) %*% H
  Lc <- H %*% gram(y) %*% H
  sum(Kc * Lc) / n^2
}

test_that("hsic detects independence and non-linear dependence", {
  set.seed(1)
  x0 <- rnorm(200)
  y0 <- rnorm(200)
  H0 <- hsic(cbind(x = x0, y = y0))
  expect_s3_class(H0, "hsic")
  expect_lt(H0[1, 2], 0.10)

  set.seed(2)
  x1 <- rnorm(200)
  y1 <- x1^2 + rnorm(200, sd = 0.1)
  H1 <- hsic(cbind(x = x1, y = y1))
  expect_gt(H1[1, 2], 0.20)
})

test_that("hsic output is symmetric with expected diagonal", {
  set.seed(3)
  X <- matrix(rnorm(120), ncol = 3)
  colnames(X) <- c("a", "b", "c")
  H <- hsic(X)

  expect_equal(unname(diag(H)), rep(1, ncol(X)))
  expect_equal(H, t(H), tolerance = 1e-12)
  expect_equal(dimnames(H), list(colnames(X), colnames(X)))
  expect_true(is.matrix(attr(H, "hsic_raw")))
})

test_that("biased HSIC matches a direct R implementation", {
  x <- c(-1.2, -0.4, 0.1, 0.9, 1.4)
  y <- c(0.2, 0.5, -0.1, 1.0, 1.8)
  H <- hsic(cbind(x = x, y = y), normalise = FALSE, estimator = "biased")

  expect_equal(
    as.numeric(H[1, 2]),
    hsic_gaussian_raw_R(x, y),
    tolerance = 1e-12
  )
})

test_that("unbiased HSIC retains signed raw estimates", {
  set.seed(4)
  X <- cbind(a = rnorm(30), b = rnorm(30), c = rnorm(30))
  H <- hsic(X, estimator = "unbiased")
  raw <- attr(H, "hsic_raw")

  expect_true(is.matrix(raw))
  expect_equal(dim(raw), dim(H))
  expect_true(all(is.na(H[upper.tri(H)]) | (H[upper.tri(H)] >= 0 & H[upper.tri(H)] <= 1)))
})

test_that("hsic missing-data modes follow dcor-style handling", {
  X <- cbind(
    a = c(1, 2, 3, 4, 5, 6),
    b = c(1, NA, 2, 4, 5, 6),
    c = c(6, 5, 4, 3, 2, 1)
  )

  expect_error(hsic(X, na_method = "error"), "missing|finite|NA|NaN|infinite")

  H_complete <- hsic(X, na_method = "complete")
  expect_s3_class(H_complete, "hsic")
  expect_equal(attr(H_complete, "diagnostics")$n_complete, 5L)

  H_pairwise <- hsic(X, na_method = "pairwise")
  n_complete <- attr(H_pairwise, "diagnostics")$n_complete
  expect_equal(n_complete["a", "b"], 5)
  expect_equal(n_complete["a", "c"], 6)
})

test_that("hsic permutation inference is reproducible and stored as metadata", {
  set.seed(5)
  x <- rnorm(60)
  y <- x^2 + rnorm(60, sd = 0.2)
  H1 <- hsic(cbind(x = x, y = y), p_value = TRUE, B = 19L, seed = 11L)
  H2 <- hsic(cbind(x = x, y = y), p_value = TRUE, B = 19L, seed = 11L)

  inf1 <- attr(H1, "inference")
  inf2 <- attr(H2, "inference")
  expect_equal(inf1$p_value, inf2$p_value)
  expect_equal(inf1$B, 19L)
  expect_equal(inf1$alternative, "greater")
  expect_true(is.finite(inf1$p_value[1, 2]))
})

test_that("hsic supports sparse and edge-list output conversion", {
  set.seed(6)
  X <- matrix(rnorm(80), ncol = 4)
  colnames(X) <- paste0("v", 1:4)

  sp <- hsic(X, output = "sparse", threshold = 0, diag = FALSE)
  ed <- hsic(X, output = "edge_list", threshold = 0, diag = FALSE)

  expect_s4_class(sp, "sparseMatrix")
  expect_equal(attr(sp, "corr_output_class"), "corr_sparse")
  expect_s3_class(ed, "corr_edge_list")
  expect_true(all(c("row", "col", "value") %in% names(ed)))
})
