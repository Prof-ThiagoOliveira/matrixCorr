schafer_shrink_base <- function(X) {
  stopifnot(is.matrix(X))
  n <- nrow(X)
  R <- stats::cor(X)

  zero <- is.na(diag(R))
  if (any(zero)) {
    R[zero, ] <- NA_real_
    R[, zero] <- NA_real_
  }

  off <- upper.tri(R, diag = FALSE)
  r <- R[off]
  r <- r[is.finite(r)]
  sum_sq <- sum(r^2)

  if (sum_sq > 0) {
    v <- ((1 - r^2)^2) / (n - 1)
    lambda <- sum(v) / sum_sq
    lambda <- min(max(lambda, 0), 1)
  } else {
    lambda <- 1
  }

  Rsh <- R
  Rsh[off] <- (1 - lambda) * R[off]
  Rsh[lower.tri(Rsh)] <- t(Rsh)[lower.tri(Rsh)]
  diag(Rsh) <- ifelse(is.na(diag(R)), NA_real_, 1)
  attr(Rsh, "lambda") <- lambda
  Rsh
}

drop_schafer_attrs <- function(x) {
  m <- as.matrix(x)
  attributes(m) <- attributes(m)[c("dim", "dimnames")]
  m
}

test_that("shrinkage_corr returns shrinkage matrix with metadata", {
  set.seed(123)
  X <- matrix(rnorm(80), nrow = 20, ncol = 4)
  colnames(X) <- paste0("V", seq_len(4))

  R <- shrinkage_corr(X)

  expect_s3_class(R, "shrinkage_corr")
  expect_s3_class(R, "schafer_corr")
  expect_equal(dim(R), c(4, 4))
  expect_true(all(diag(R) == 1))
  expect_equal(attr(R, "method"), "schafer_shrinkage")
  expect_equal(attr(R, "package"), "matrixCorr")

  R_raw <- stats::cor(X)
  off <- upper.tri(R_raw, diag = FALSE)
  expect_true(mean(abs(R[off])) <= mean(abs(R_raw[off])) + 1e-10)
})

test_that("base-R reference agrees within tolerance and attributes are sane", {
  set.seed(123)
  n <- 80
  p <- 30
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", seq_len(p))

  R_pkg <- shrinkage_corr(X)
  R_base <- schafer_shrink_base(X)

  expect_s3_class(R_pkg, "shrinkage_corr")
  expect_s3_class(R_pkg, "schafer_corr")
  expect_true(is.matrix(R_pkg))
  expect_true(isSymmetric(R_pkg))
  expect_identical(attr(R_pkg, "method"), "schafer_shrinkage")
  expect_true(grepl("Scha", attr(R_pkg, "description"), fixed = TRUE))
  expect_equal(drop_schafer_attrs(R_pkg), drop_schafer_attrs(R_base), tolerance = 2e-4)
  expect_true(all(diag(R_pkg) == 1))
})

test_that("deterministic small example matches the base reference tightly", {
  X <- matrix(
    c(
      1.1, 0.3, -0.7,
      0.9, 0.4, -0.5,
      1.2, 0.2, -0.6,
      0.8, 0.5, -0.4,
      1.0, 0.1, -0.8
    ),
    nrow = 5, ncol = 3, byrow = TRUE
  )
  colnames(X) <- c("A", "B", "C")

  R_pkg <- shrinkage_corr(X)
  R_base <- schafer_shrink_base(X)

  expect_equal(unname(drop_schafer_attrs(R_pkg)),
               unname(drop_schafer_attrs(R_base)),
               tolerance = 1e-10)
})

test_that("shrinkage_corr flags zero-variance columns as NA", {
  X <- cbind(const = rep(2, 6), var = rnorm(6))

  R <- shrinkage_corr(X)

  expect_true(all(is.na(R["const", ])))
  expect_true(all(is.na(R[, "const"])))
})

test_that("data.frame path ignores non-numeric columns and preserves names", {
  set.seed(7)
  n <- 40
  df <- data.frame(
    a = rnorm(n),
    b = rnorm(n),
    grp = rep(letters[1:4], length.out = n),
    flag = rep(c(TRUE, FALSE), length.out = n)
  )
  out <- shrinkage_corr(df)

  expect_equal(dim(out), c(2, 2))
  expect_setequal(colnames(out), c("a", "b"))
  expect_true(isSymmetric(out))
  expect_true(all(diag(out) == 1))
})

test_that("independent columns are shrunk closer to identity than raw correlations", {
  set.seed(99)
  n <- 60
  p <- 40
  X <- matrix(rnorm(n * p), n, p)

  R_raw <- stats::cor(X)
  R_shr <- shrinkage_corr(X)
  I <- diag(p)

  err_raw <- sqrt(sum((R_raw - I)^2))
  err_shr <- sqrt(sum((R_shr - I)^2))
  expect_lt(err_shr, err_raw)

  off <- upper.tri(R_raw, diag = FALSE)
  r <- R_raw[off]
  lambda_hat <- {
    r2 <- r^2
    sum_sq <- sum(r2[is.finite(r2)])
    if (sum_sq > 0) {
      v <- ((1 - r^2)^2) / (n - 1)
      min(max(sum(v[is.finite(v)]) / sum_sq, 0), 1)
    } else {
      1
    }
  }
  expect_equal(unname(as.vector(R_shr[off])),
               unname(as.vector((1 - lambda_hat) * r)),
               tolerance = 5e-5)
})

test_that("p >> n returns a valid symmetric PSD correlation matrix", {
  set.seed(2024)
  n <- 30
  p <- 120
  X <- matrix(rnorm(n * p), n, p)

  R_shr <- shrinkage_corr(X)

  expect_equal(dim(R_shr), c(p, p))
  expect_true(isSymmetric(R_shr))
  ev <- eigen(R_shr, symmetric = TRUE, only.values = TRUE)$values
  expect_gte(min(ev), -1e-8)
})

test_that("print and plot methods cover optional arguments", {
  skip_if_not_installed("ggplot2")

  set.seed(456)
  X <- matrix(rnorm(60), nrow = 15, ncol = 4)
  colnames(X) <- paste0("G", seq_len(4))
  R <- shrinkage_corr(X)

  out <- capture.output(print(R, digits = 3, max_rows = 2, max_cols = 3))
  expect_true(any(grepl("omitted", out)))

  p1 <- plot(R, cluster = TRUE, triangle = "lower", show_values = TRUE, value_text_limit = 10)
  expect_s3_class(p1, "ggplot")

  if (requireNamespace("viridisLite", quietly = TRUE)) {
    p2 <- plot(R, cluster = FALSE, triangle = "upper", palette = "viridis")
    expect_s3_class(p2, "ggplot")
  }
})

test_that("schafer_corr remains a compatibility alias", {
  set.seed(55)
  X <- matrix(rnorm(50), nrow = 10, ncol = 5)

  new_fit <- shrinkage_corr(X)
  old_fit <- schafer_corr(X)

  expect_identical(unclass(new_fit), unclass(old_fit))
  expect_s3_class(old_fit, "shrinkage_corr")
  expect_s3_class(old_fit, "schafer_corr")
})

test_that("shrinkage_corr honors n_threads without changing estimates", {
  set.seed(808)
  X <- matrix(rnorm(240), nrow = 40, ncol = 6)
  colnames(X) <- paste0("G", seq_len(ncol(X)))

  fit1 <- shrinkage_corr(X, n_threads = 1L)
  fit2 <- shrinkage_corr(X, n_threads = 2L)
  alias_fit <- schafer_corr(X, n_threads = 2L)

  expect_equal(unclass(fit1), unclass(fit2), tolerance = 1e-12)
  expect_equal(unclass(fit2), unclass(alias_fit), tolerance = 1e-12)
})
