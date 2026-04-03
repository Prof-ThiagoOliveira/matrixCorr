pcor_from_precision <- function(Theta) {
  d <- diag(Theta)
  S <- -Theta / sqrt(outer(d, d))
  diag(S) <- 1
  S
}

strip_matrix_metadata <- function(x) {
  attrs <- attributes(x)
  keep <- intersect(names(attrs), c("dim", "dimnames"))
  attributes(x) <- attrs[keep]
  x
}

manual_partial_fisher_ci_R <- function(X, conf_level = 0.95) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  c_vars <- p - 2L
  S <- stats::cov(X)
  Theta <- solve(S)
  pcor <- pcor_from_precision(Theta)

  lwr <- matrix(NA_real_, p, p, dimnames = dimnames(pcor))
  upr <- matrix(NA_real_, p, p, dimnames = dimnames(pcor))
  diag(lwr) <- 1
  diag(upr) <- 1

  se <- 1 / sqrt(n - 3 - c_vars)
  crit <- stats::qnorm(0.5 * (1 + conf_level))
  eps <- sqrt(.Machine$double.eps)

  for (j in seq_len(p - 1L)) {
    for (i in (j + 1L):p) {
      r <- pcor[j, i]
      if (!is.finite(r) || abs(r) >= 1) next
      z <- atanh(max(min(r, 1 - eps), -1 + eps))
      lo <- tanh(z - crit * se)
      hi <- tanh(z + crit * se)
      lwr[j, i] <- lwr[i, j] <- lo
      upr[j, i] <- upr[i, j] <- hi
    }
  }

  list(
    pcor = pcor,
    lwr = lwr,
    upr = upr,
    n_complete = n,
    n_conditioning = c_vars
  )
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

glasso_reference_case <- function() {
  vars <- paste0("V", seq_len(6))

  Sigma <- structure(
    c(0.892577511863956, -3.47460427750079e-05, -0.048912146675666,
      0.00145719602008398, 0.00184184223658341, 7.32384941200259e-08,
      -3.47460427750079e-05, 1.18178348379229, 0.000607657901650682,
      -1.81034106257749e-05, -0.0211063204776342, -0.00249098993197849,
      -0.048912146675666, 0.000607657901650682, 0.855402516097313,
      -0.0254842452589996, -0.0322111497962843, -1.28083505639143e-06,
      0.00145719602008398, -1.81034106257749e-05, -0.0254842452589996,
      1.08365876399129, 0.000959638095557692, 3.81587780019788e-08,
      0.00184184223658341, -0.0211063204776342, -0.0322111497962843,
      0.000959638095557692, 1.11881841527388, 4.44883792437081e-05,
      7.32384941200259e-08, -0.00249098993197849, -1.28083505639143e-06,
      3.81587780019788e-08, 4.44883792437081e-05, 0.894638069650798),
    dim = c(6L, 6L),
    dimnames = list(vars, vars)
  )

  Theta <- structure(
    c(1.12387243015938, 0, 0.0642633288004597, 0, 0, 0,
      0, 0.846468847136127, 0, 0, 0.0159683981780265, 0.0023560756458958,
      0.0642633288004597, 0, 1.1748032714163, 0.0275114249784953,
      0.033693582191095, 0, 0, 0, 0.0275114249784953, 0.923446698493756, 0, 0,
      0, 0.0159683981780265, 0.033693582191095, 0, 0.895071380200203, 0,
      0, 0.0023560756458958, 0, 0, 0, 1.11777701272096),
    dim = c(6L, 6L),
    dimnames = list(vars, vars)
  )

  list(lambda = 0.08, Sigma = Sigma, Theta = Theta)
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

  glasso <- pcorr(X, method = "glasso", lambda = 0.05, return_cov_precision = TRUE)
  expect_equal(glasso$method, "glasso")
  expect_equal(glasso$lambda, 0.05)
  expect_true(is.na(glasso$rho))
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

test_that("pcorr rejects non-finite matrix inputs without an R-side copy", {
  X <- matrix(rnorm(16), nrow = 4)
  X[3] <- Inf
  expect_error(pcorr(X, method = "sample"), "finite values")
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
  expect_null(samp$p_value)
})

test_that("pcorr default CI behaviour remains estimate-only", {
  set.seed(334)
  X <- matrix(rnorm(120), nrow = 30, ncol = 4)

  fit <- pcorr(X, method = "sample")

  expect_null(fit$ci)
  expect_null(attr(fit, "ci", exact = TRUE))
  expect_null(attr(fit, "conf.level", exact = TRUE))
  expect_true(is.list(fit$diagnostics))
})

test_that("pcorr Fisher-z CI matches manual partial-correlation calculation", {
  set.seed(335)
  X <- matrix(rnorm(240), nrow = 60, ncol = 4)
  colnames(X) <- paste0("V", 1:4)

  fit <- pcorr(X, method = "sample", ci = TRUE, return_p_value = TRUE)
  ref <- manual_partial_fisher_ci_R(X, conf_level = 0.95)
  ci <- attr(fit, "ci", exact = TRUE)

  expect_true(is.list(ci))
  expect_identical(ci$ci.method, "fisher_z_partial")
  expect_equal(ci$conf.level, 0.95)
  expect_equal(
    strip_matrix_metadata(fit$pcor),
    strip_matrix_metadata(ref$pcor),
    tolerance = 1e-10
  )
  expect_equal(ci$lwr.ci, ref$lwr, tolerance = 1e-10)
  expect_equal(ci$upr.ci, ref$upr, tolerance = 1e-10)
  expect_identical(unique(as.integer(fit$diagnostics$n_complete)), as.integer(ref$n_complete))
  expect_identical(unique(as.integer(fit$diagnostics$n_conditioning[upper.tri(fit$diagnostics$n_conditioning)])), as.integer(ref$n_conditioning))
  expect_true(is.matrix(fit$p_value))
})

test_that("pcorr custom conf_level changes CI width", {
  set.seed(336)
  X <- matrix(rnorm(320), nrow = 80, ncol = 4)
  colnames(X) <- paste0("V", 1:4)

  fit90 <- pcorr(X, method = "sample", ci = TRUE, conf_level = 0.90)
  fit99 <- pcorr(X, method = "sample", ci = TRUE, conf_level = 0.99)

  ci90 <- attr(fit90, "ci", exact = TRUE)
  ci99 <- attr(fit99, "ci", exact = TRUE)
  width90 <- ci90$upr.ci["V1", "V2"] - ci90$lwr.ci["V1", "V2"]
  width99 <- ci99$upr.ci["V1", "V2"] - ci99$lwr.ci["V1", "V2"]

  expect_equal(attr(fit90, "conf.level", exact = TRUE), 0.90)
  expect_equal(attr(fit99, "conf.level", exact = TRUE), 0.99)
  expect_true(width99 > width90)
})

test_that("pcorr CI integrates with the existing summary contract", {
  set.seed(337)
  X <- matrix(rnorm(240), nrow = 60, ncol = 4)
  colnames(X) <- paste0("V", 1:4)

  fit <- pcorr(X, method = "sample", ci = TRUE)
  sm <- summary(fit)
  txt <- capture.output(print(sm, show_ci = "yes"))

  expect_s3_class(sm, "summary_partial_corr")
  expect_s3_class(sm, "data.frame")
  expect_true(isTRUE(attr(sm, "has_ci")))
  expect_true(all(c("estimate", "lwr", "upr") %in% names(sm)))
  expect_match(paste(txt, collapse = "\n"), "Partial correlation summary")
  expect_match(paste(txt, collapse = "\n"), "ci_method")
  expect_match(paste(txt, collapse = "\n"), "ci_width")
  expect_match(paste(txt, collapse = "\n"), "Strongest pairs by \\|estimate\\|")
  expect_true(any(grepl("\\blwr\\b", txt)))
  expect_true(any(grepl("\\bupr\\b", txt)))
})

test_that("sample partial-correlation p-values match the reference example", {
  y.data <- data.frame(
    hl = c(7, 15, 19, 15, 21, 22, 57, 15, 20, 18),
    disp = c(0.000, 0.964, 0.000, 0.000, 0.921, 0.000, 0.000, 1.006, 0.000, 1.011),
    deg = c(9, 2, 3, 4, 1, 3, 1, 3, 6, 1),
    BC = c(1.78e-02, 1.05e-06, 1.37e-05, 7.18e-03, 0.00e+00,
           0.00e+00, 0.00e+00, 4.48e-03, 2.10e-06, 0.00e+00)
  )

  ours <- pcorr(y.data, method = "sample", return_p_value = TRUE)

  vars <- c("hl", "disp", "deg", "BC")
  expected_pcor <- structure(
    c(
      1.0000000, -0.6720863, -0.6161163,  0.1148459,
     -0.6720863,  1.0000000, -0.7215522,  0.2855420,
     -0.6161163, -0.7215522,  1.0000000,  0.6940953,
      0.1148459,  0.2855420,  0.6940953,  1.0000000
    ),
    dim = c(4L, 4L),
    dimnames = list(vars, vars)
  )

  expected_p_value <- structure(
    c(
      0.00000000, 0.06789202, 0.10383620, 0.78654997,
      0.06789202, 0.00000000, 0.04332869, 0.49299871,
      0.10383620, 0.04332869, 0.00000000, 0.05615021,
      0.78654997, 0.49299871, 0.05615021, 0.00000000
    ),
    dim = c(4L, 4L),
    dimnames = list(vars, vars)
  )

  expect_equal(
    strip_matrix_metadata(ours$pcor),
    strip_matrix_metadata(expected_pcor),
    tolerance = 1e-7
  )
  expect_equal(ours$p_value, expected_p_value, tolerance = 1e-8)
  expect_equal(unname(diag(ours$p_value)), rep(0, ncol(y.data)))
})

test_that("pcorr p-values are restricted to the classical sample setting", {
  X <- matrix(rnorm(60), nrow = 20, ncol = 3)

  expect_error(
    pcorr(X, method = "ridge", return_p_value = TRUE),
    "available only for"
  )

  X_wide <- matrix(rnorm(36), nrow = 6, ncol = 6)
  expect_error(
    pcorr(X_wide, method = "sample", return_p_value = TRUE),
    "requires .*n > p"
  )
})

test_that("pcorr CI is restricted to the classical sample setting", {
  X <- matrix(rnorm(80), nrow = 20, ncol = 4)

  expect_error(
    pcorr(X, method = "ridge", ci = TRUE),
    "available only for the classical .*method = \"sample\""
  )
  expect_error(
    pcorr(X, method = "oas", ci = TRUE),
    "available only for the classical .*method = \"sample\""
  )
  expect_error(
    pcorr(X, method = "glasso", ci = TRUE),
    "available only for the classical .*method = \"sample\""
  )
})

test_that("pcorr CI requires positive Fisher-z residual degrees of freedom", {
  X <- matrix(rnorm(20), nrow = 5, ncol = 4)

  expect_error(
    pcorr(X, method = "sample", ci = TRUE),
    "requires .*n > p \\+ 1"
  )
})

test_that("pcorr CI warns and returns NA bounds when the sample covariance needs jitter", {
  x1 <- 1:20
  x2 <- x1
  x3 <- rep(c(0, 1), each = 10)
  X <- cbind(x1 = x1, x2 = x2, x3 = x3)

  expect_warning(
    fit <- pcorr(X, method = "sample", ci = TRUE),
    class = "matrixCorr_ci_warning"
  )

  ci <- attr(fit, "ci", exact = TRUE)
  expect_true(all(is.na(ci$lwr.ci[upper.tri(ci$lwr.ci)])))
  expect_true(all(is.na(ci$upr.ci[upper.tri(ci$upr.ci)])))
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

  expect_equal(
    strip_matrix_metadata(ours$pcor),
    strip_matrix_metadata(ref),
    tolerance = 1e-10
  )
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

  expect_equal(
    strip_matrix_metadata(ours$pcor),
    strip_matrix_metadata(ref),
    tolerance = 1e-10
  )
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

  expect_equal(
    strip_matrix_metadata(ours$pcor),
    strip_matrix_metadata(ref),
    tolerance = 1e-10
  )
  expect_true(isSymmetric(ours$pcor))
  expect_equal(as.numeric(diag(ours$pcor)), rep(1, p))
})

test_that("glasso method matches fixed graphical-lasso reference values", {
  set.seed(1234)
  n <- 120
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("V", seq_len(p))
  ref <- glasso_reference_case()
  lambda <- ref$lambda

  ours <- pcorr(X, method = "glasso", lambda = lambda, return_cov_precision = TRUE)
  ref_pcor <- pcor_from_precision(ref$Theta)
  dimnames(ref_pcor) <- dimnames(ours$pcor)

  expect_equal(ours$precision, ref$Theta, tolerance = 2e-4)
  expect_equal(ours$cov, ref$Sigma, tolerance = 2e-4)
  expect_equal(
    strip_matrix_metadata(ours$pcor),
    strip_matrix_metadata(ref_pcor),
    tolerance = 2e-4
  )
  expect_true(isSymmetric(ours$pcor, tol = 1e-12))
  expect_equal(as.numeric(diag(ours$pcor)), rep(1, p), tolerance = 1e-12)
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
