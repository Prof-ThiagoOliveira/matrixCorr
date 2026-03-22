pcor_from_precision <- function(Theta) {
  d <- diag(Theta)
  S <- -Theta / sqrt(outer(d, d))
  diag(S) <- 1
  S
}

oas_shrink_R <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  mu <- colMeans(X)
  Xc <- sweep(X, 2, mu, check.margin = FALSE)
  S_mle <- crossprod(Xc) / n
  trS <- sum(diag(S_mle))
  trS2 <- sum(S_mle * S_mle)
  mu_sc <- trS / p
  num <- (1 - 2 / p) * trS2 + trS^2
  den <- (n + 1 - 2 / p) * (trS2 - trS^2 / p)
  rho <- if (den > 0) num / den else 1
  rho <- max(0, min(1, rho))
  Sigma <- (1 - rho) * S_mle
  diag(Sigma) <- diag(Sigma) + rho * mu_sc
  list(Sigma = Sigma, rho = rho)
}

test_that("pcorr returns expected components for each method", {
  set.seed(111)
  X <- matrix(rnorm(200), nrow = 40, ncol = 5)
  colnames(X) <- paste0("V", seq_len(5))

  samp <- pcorr(X, method = "sample", return_cov_precision = TRUE)
  expect_s3_class(samp, "partial_corr")
  expect_equal(dim(samp$pcor), c(5, 5))
  expect_true(all(diag(samp$pcor) == 1))
  expect_true(is.matrix(samp$cov))
  expect_true(is.matrix(samp$precision))
  expect_equal(samp$method, "sample")
  expect_true(is.na(samp$lambda))
  expect_true(is.na(samp$rho))

  ridge <- pcorr(X, method = "ridge", lambda = 1e-2, return_cov_precision = TRUE)
  expect_equal(ridge$method, "ridge")
  expect_equal(ridge$lambda, 1e-2)
  expect_true(is.na(ridge$rho))

  oas <- pcorr(X, method = "oas", return_cov_precision = TRUE)
  expect_equal(oas$method, "oas")
  expect_true(is.na(oas$lambda))
  expect_true(oas$rho >= 0 && oas$rho <= 1)
})

test_that("pcorr print and plot methods cover options", {
  skip_if_not_installed("ggplot2")

  set.seed(222)
  X <- matrix(rnorm(150), nrow = 30, ncol = 5)
  colnames(X) <- paste0("G", seq_len(5))
  pc <- pcorr(X, method = "ridge", lambda = 5e-3, return_cov_precision = TRUE)

  out1 <- capture.output(print(pc, digits = 4, max_rows = 3, max_cols = 4))
  expect_true(any(grepl("omitted", out1)))

  out2 <- capture.output(print(pc, show_method = FALSE, max_rows = 2, max_cols = 2))
  expect_true(any(grepl("Partial correlation matrix", out2)))

  p <- plot(pc, mask_diag = FALSE, reorder = TRUE, value_text_size = 3,
            low_color = "blue", high_color = "red", mid_color = "white")
  expect_s3_class(p, "ggplot")
})

test_that("pcorr validates lambda", {
  set.seed(1)
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)
  expect_error(pcorr(X, method = "ridge", lambda = -1), "must be >=")
})

test_that("pcorr exposes shrinkage metadata without cov/precision", {
  set.seed(333)
  X <- matrix(rnorm(120), nrow = 30, ncol = 4)

  oas <- pcorr(X, method = "oas")
  expect_true(is.numeric(oas$rho))
  expect_true(isTRUE(oas$rho >= 0 && oas$rho <= 1))
  expect_true(is.numeric(oas$jitter))
  expect_false(is.na(oas$jitter))
  expect_null(oas$cov)
  expect_null(oas$precision)

  samp <- pcorr(X, method = "sample")
  expect_true(is.na(samp$rho))
  expect_true(is.na(samp$lambda))
  expect_true(is.numeric(samp$jitter))
})

test_that("sample partial correlation is close to truth in a structured MVN model", {
  skip_if_not_installed("MASS")
  set.seed(20240819)

  p <- 8
  alpha <- 0.3
  Theta <- diag(p)
  for (j in 1:(p - 1)) Theta[j, j + 1] <- Theta[j + 1, j] <- -alpha
  Sigma <- solve(Theta)

  n <- 1500
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(X) <- paste0("V", seq_len(p))

  ours <- pcorr(X, method = "sample", return_cov_precision = FALSE)$pcor
  truth <- pcor_from_precision(Theta)

  ut <- upper.tri(ours, diag = FALSE)
  z_hat <- atanh(ours[ut])
  z_true <- atanh(truth[ut])
  se_z <- 1 / sqrt(n - p - 1)

  expect_lt(max(abs(z_hat - z_true)), 3 * se_z)

  non_edge <- matrix(TRUE, p, p)
  diag(non_edge) <- FALSE
  non_edge[abs(row(non_edge) - col(non_edge)) == 1] <- FALSE
  expect_lt(max(abs(ours[non_edge])), 0.07)
  expect_true(isSymmetric(ours, tol = 1e-12))
  expect_equal(as.numeric(diag(ours)), rep(1, p), tolerance = 1e-12)
})

test_that("sample method equals the base-R precision construction", {
  set.seed(123)
  n <- 200
  p <- 12
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", seq_len(p))

  ours <- pcorr(X, method = "sample", return_cov_precision = TRUE)

  S_unb <- stats::cov(X)
  Theta <- solve(S_unb)
  ref <- pcor_from_precision(Theta)

  expect_equal(ours$pcor, ref, tolerance = 1e-10)
  expect_true(isSymmetric(ours$pcor))
  expect_equal(as.numeric(diag(ours$pcor)), rep(1, p))
})

test_that("ridge method equals the base-R ridge construction", {
  set.seed(99)
  n <- 180
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", seq_len(p))
  lambda <- 5e-3

  ours <- pcorr(X, method = "ridge", lambda = lambda,
                              return_cov_precision = TRUE)

  S_unb <- stats::cov(X)
  diag(S_unb) <- diag(S_unb) + lambda
  Theta <- solve(S_unb)
  ref <- pcor_from_precision(Theta)

  expect_equal(ours$pcor, ref, tolerance = 1e-10)
  I <- diag(p)
  dimnames(I) <- dimnames(S_unb)
  expect_equal(S_unb %*% Theta, I, tolerance = 1e-8)
})

test_that("OAS method matches an R implementation of the same formula", {
  set.seed(7)
  n <- 120
  p <- 25
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", seq_len(p))

  ours <- pcorr(X, method = "oas", return_cov_precision = TRUE)

  oas <- oas_shrink_R(X)
  Sigma <- oas$Sigma
  Theta <- solve(Sigma)
  ref <- pcor_from_precision(Theta)

  expect_equal(ours$pcor, ref, tolerance = 1e-10)
  expect_true(isSymmetric(ours$pcor))
  expect_equal(as.numeric(diag(ours$pcor)), rep(1, p))
})

test_that("p >> n: OAS returns a finite, well-formed matrix", {
  set.seed(4242)
  n <- 40
  p <- 80
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", seq_len(p))

  oas <- pcorr(X, method = "oas", return_cov_precision = TRUE)
  M <- oas$pcor
  expect_true(is.matrix(M))
  expect_true(isSymmetric(M, tol = 1e-12))
  expect_equal(as.numeric(diag(M)), rep(1, p), tolerance = 1e-12)
  expect_true(all(is.finite(M)))
  expect_true(all(abs(M[upper.tri(M, diag = FALSE)]) <= 1 + 1e-10))

  I <- diag(p)
  dimnames(I) <- dimnames(oas$cov)
  expect_equal(oas$cov %*% oas$precision, I, tolerance = 1e-6)
})

test_that("non-numeric columns are ignored and dimnames propagate", {
  set.seed(1)
  X <- data.frame(
    a = rnorm(50),
    b = rnorm(50),
    f = factor(sample(letters[1:3], 50, TRUE)),
    c = rnorm(50),
    s = sample(c("x", "y", "z"), 50, TRUE),
    l = sample(c(TRUE, FALSE), 50, TRUE)
  )
  cols <- c("a", "b", "c")

  pc <- pcorr(X, method = "oas", return_cov_precision = FALSE)$pcor
  expect_equal(dim(pc), c(3L, 3L))
  expect_equal(dimnames(pc), list(cols, cols))
})
