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

    # Create tied data
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

test_that("kendall_tau is invariant to strictly monotone transformations", {
  set.seed(123)
  X <- matrix(rnorm(200), ncol = 2)
  X_mono <- X
  X_mono[,1] <- exp(X_mono[,1])           # monotone transform
  X_mono[,2] <- log1p(exp(X_mono[,2]))    # monotone transform

  base_cor <- cor(X, method = "kendall")
  fast_cor <- kendall_tau(X_mono)

  attributes(fast_cor) <- NULL
  expect_equal(fast_cor, base_cor, tolerance = 1e-8)
})

test_that("kendall_tau returns NA when a column is constant", {
  X <- cbind(a = rnorm(20), b = rep(1, 20))
  kt <- kendall_tau(X)

  expect_true(all(is.na(kt["b", ])))
  expect_true(all(is.na(kt[, "b"])))
})

test_that("kendall_tau matches base::cor on a known toy dataset", {
  X <- matrix(c(1, 2, 3, 4,
                4, 3, 2, 1), ncol = 2)
  colnames(X) <- c("x", "y")

  base_cor <- cor(X, method = "kendall")
  fast_cor <- kendall_tau(X)

  attributes(fast_cor) <- NULL
  expect_equal(fast_cor, base_cor, tolerance = 1e-12)
})

test_that("kendall_tau estimates agree with theoretical BVN relationship", {
  # For BVN, tau = (2/pi) * arcsin(r)
  set.seed(321)
  rhos <- c(-0.8, -0.4, 0, 0.4, 0.8)
  n <- 2000

  for (r in rhos) {
    Sigma <- matrix(c(1, r, r, 1), 2, 2)
    Z <- MASS::mvrnorm(n, mu = c(0,0), Sigma = Sigma)

    est <- kendall_tau(Z)[1,2]
    theory <- (2 / pi) * asin(r)

    expect_equal(est, theory, tolerance = 0.05,
                 info = paste("Mismatch for true rho =", r))
  }
})
