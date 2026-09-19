robust_ccc_reference <- function(x, y, alpha = 0.75, seed = NULL,
                                 nsamp = 500L) {
  evaluate <- function() {
    fit <- tryCatch(
      suppressWarnings(
        robustbase::covMcd(
          cbind(x, y),
          alpha = alpha,
          raw.only = FALSE,
          nsamp = nsamp
        )
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(NA_real_)
    }
    mu <- fit$center
    scatter <- fit$cov
    denominator <- scatter[1, 1] + scatter[2, 2] + (mu[1] - mu[2])^2
    if (!is.finite(denominator) || denominator <= 0 ||
        !is.finite(scatter[1, 1]) || !is.finite(scatter[2, 2]) ||
        scatter[1, 1] <= 0 || scatter[2, 2] <= 0) {
      return(NA_real_)
    }
    value <- 2 * scatter[1, 2] / denominator
    if (is.finite(value)) value else NA_real_
  }

  if (is.null(seed)) {
    return(evaluate())
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  evaluate()
}

robust_ccc_percentile_reference <- function(values, conf_level) {
  values <- sort(values[is.finite(values)])
  if (!length(values)) {
    return(c(NA_real_, NA_real_))
  }
  alpha <- 1 - conf_level
  low <- floor((alpha / 2) * length(values) + 0.5)
  high <- floor((1 - alpha / 2) * length(values) + 0.5)
  low <- min(max(low, 1L), length(values))
  high <- min(max(high, 1L), length(values))
  c(values[low], values[high])
}

test_that("robust_ccc matches the joint reweighted MCD formula", {
  x <- c(-2.2, -1.9, -1.5, -1.1, -0.8, -0.4, -0.1, 0.3, 0.7, 1.0,
         1.4, 1.8, 2.1, 2.5, 6.5, -7.0)
  y <- c(-2.0, -1.7, -1.4, -1.0, -0.7, -0.5, 0.0, 0.4, 0.8, 1.2,
         1.3, 1.9, 2.0, 2.6, -6.0, 7.5)

  expected <- robust_ccc_reference(x, y, alpha = 0.75, seed = 541L)
  fit <- robust_ccc(cbind(device_a = x, device_b = y), seed = 541L)

  expect_equal(unname(fit["device_a", "device_b"]), unname(expected), tolerance = 1e-12)
  expect_identical(
    class(fit),
    c("corr_matrix", "robust_ccc", "corr_result", "matrix")
  )
  expect_identical(dimnames(fit), list(c("device_a", "device_b"), c("device_a", "device_b")))
  expect_equal(unname(fit), t(unname(fit)), tolerance = 0)
  expect_true(all(is.na(fit) | (fit >= -1 & fit <= 1)))
  expect_equal(unname(diag(fit)), c(1, 1), tolerance = 0)
  expect_identical(attr(fit, "method"), "mcd_robust_concordance")
  expect_identical(attr(fit, "alpha"), 0.75)
  expect_identical(attr(fit, "mcd_estimate"), "reweighted")
})

test_that("robust_ccc uses reweighted rather than raw MCD estimates", {
  x <- seq(-2, 2, length.out = 40)
  y <- x + rep(c(-0.08, 0.04, 0.02, -0.03), 10)
  x <- c(x, 8, -9, 10)
  y <- c(y, -8, 9, -11)

  set.seed(901)
  mcd <- suppressWarnings(
    robustbase::covMcd(cbind(x, y), alpha = 0.75, raw.only = FALSE, nsamp = 500)
  )
  coefficient <- function(mu, scatter) {
    2 * scatter[1, 2] /
      (scatter[1, 1] + scatter[2, 2] + (mu[1] - mu[2])^2)
  }
  expected_reweighted <- coefficient(mcd$center, mcd$cov)
  expected_raw <- coefficient(mcd$raw.center, mcd$raw.cov)

  fit <- robust_ccc(cbind(x = x, y = y), alpha = 0.75, seed = 901L)

  expect_gt(abs(expected_reweighted - expected_raw), 1e-5)
  expect_equal(unname(fit["x", "y"]), unname(expected_reweighted), tolerance = 1e-12)
  expect_gt(abs(unname(fit["x", "y"]) - expected_raw), 1e-5)
})

test_that("robust_ccc is stable under discordant leverage contamination", {
  x <- seq(-2, 2, length.out = 45)
  y <- x + 0.04 * sin(seq_along(x))
  clean <- cbind(x = x, y = y)
  contaminated <- rbind(clean, c(12, -11), c(-13, 12))

  classical_clean <- unname(ccc(clean)["x", "y"])
  classical_contaminated <- unname(ccc(contaminated)["x", "y"])
  robust_clean <- unname(robust_ccc(clean, seed = 101L)["x", "y"])
  robust_contaminated <- unname(robust_ccc(contaminated, seed = 101L)["x", "y"])

  expect_gt(abs(classical_contaminated - classical_clean), 0.5)
  expect_lt(abs(robust_contaminated - robust_clean), 0.05)
  expect_gt(robust_contaminated, classical_contaminated)
})

test_that("robust_ccc validates alpha and seeded point estimation", {
  set.seed(12)
  X <- cbind(a = rnorm(35), b = rnorm(35), c = rnorm(35))
  X[c(1, 3), 1] <- c(9, -8)
  X[c(1, 3), 2] <- c(-7, 10)

  fit_half <- robust_ccc(X, alpha = 0.5, seed = 19L)
  fit_one <- robust_ccc(X, alpha = 1, seed = 19L)
  expect_true(any(is.finite(fit_half[upper.tri(fit_half)])))
  expect_true(any(is.finite(fit_one[upper.tri(fit_one)])))
  expect_false(isTRUE(all.equal(
    unclass(fit_half),
    unclass(fit_one),
    tolerance = 1e-8,
    check.attributes = FALSE
  )))

  expect_error(robust_ccc(X, alpha = 0.49), class = "matrixCorr_arg_error")
  expect_error(robust_ccc(X, alpha = 1.01), class = "matrixCorr_arg_error")
  expect_error(robust_ccc(X, alpha = NA_real_), class = "matrixCorr_arg_error")
  expect_error(robust_ccc(X, seed = "bad"), class = "matrixCorr_arg_error")

  first <- robust_ccc(X, seed = 733L)
  second <- robust_ccc(X, seed = 733L)
  expect_identical(first, second)
})

test_that("robust_ccc diagonal and degenerate cases follow robust scatter", {
  X <- cbind(
    varying = seq_len(14),
    constant = rep(3, 14),
    matching = seq_len(14)
  )
  fit <- robust_ccc(X, seed = 32L)

  expect_identical(unname(fit["varying", "varying"]), 1)
  expect_identical(unname(fit["matching", "matching"]), 1)
  expect_true(is.na(fit["constant", "constant"]))
  expect_true(all(is.na(fit["constant", ])))
  expect_equal(unname(fit["varying", "matching"]), 1, tolerance = 1e-12)

  too_short <- robust_ccc(
    cbind(a = c(1, 2), b = c(1, 3)),
    na_method = "pairwise",
    seed = 9L
  )
  expect_true(all(is.na(too_short)))
})

test_that("robust_ccc implements error complete and pairwise NA semantics", {
  x <- c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, NA, 3.5)
  y <- c(-1.9, -1.4, -0.8, -0.6, 0.1, 0.6, 1.2, 1.4, 2.2, NA, 3.1, 3.4)
  z <- c(1, 2, 3, 4, 5, 6, 7, 8, NA, 10, 11, 12)
  X <- cbind(x = x, y = y, z = z)

  expect_error(robust_ccc(X, na_method = "error"), "Missing values are not allowed")
  X_inf <- X
  X_inf[1, 1] <- Inf
  expect_error(robust_ccc(X_inf, na_method = "error"), "non-finite values")

  pair_ok <- is.finite(x) & is.finite(y)
  pairwise <- robust_ccc(X, na_method = "pairwise", seed = 411L)
  expect_equal(
    unname(pairwise["x", "y"]),
    unname(robust_ccc_reference(x[pair_ok], y[pair_ok], seed = 411L)),
    tolerance = 1e-12
  )

  common_ok <- apply(is.finite(X), 1L, all)
  complete <- robust_ccc(X, na_method = "complete", seed = 411L)
  expect_equal(
    unname(complete["x", "y"]),
    unname(robust_ccc_reference(x[common_ok], y[common_ok], seed = 411L)),
    tolerance = 1e-12
  )
  expect_identical(attr(complete, "diagnostics")$n_complete, as.integer(sum(common_ok)))
  expect_false(isTRUE(all.equal(
    unname(pairwise["x", "y"]),
    unname(complete["x", "y"]),
    tolerance = 0
  )))
})

test_that("robust_ccc bootstrap is paired percentile inference and reproducible", {
  x <- c(-2.1, -1.8, -1.4, -1.0, -0.7, -0.3, 0.1, 0.4, 0.8, 1.1,
         1.5, 1.9, 2.2, 2.6)
  y <- c(-1.9, -1.7, -1.2, -1.1, -0.5, -0.4, 0.2, 0.5, 0.9, 1.0,
         1.6, 1.8, 2.4, 2.5)
  n_boot <- 16L
  conf_level <- 0.90
  seed <- 612L

  set.seed(seed)
  bootstrap_values <- vapply(
    seq_len(n_boot),
    function(b) {
      index <- sample.int(length(x), size = length(x), replace = TRUE)
      robust_ccc_reference(x[index], y[index], alpha = 0.75, seed = NULL)
    },
    numeric(1L)
  )
  expected <- robust_ccc_percentile_reference(bootstrap_values, conf_level)

  first <- robust_ccc(
    cbind(x = x, y = y),
    ci = TRUE,
    conf_level = conf_level,
    n_boot = n_boot,
    seed = seed
  )
  second <- robust_ccc(
    cbind(x = x, y = y),
    ci = TRUE,
    conf_level = conf_level,
    n_boot = n_boot,
    seed = seed
  )
  interval <- attr(first, "ci")
  diagnostics <- attr(first, "diagnostics")

  expect_identical(first, second)
  expect_identical(interval$ci.method, "percentile_bootstrap")
  expect_equal(unname(interval$lwr.ci["x", "y"]), expected[1], tolerance = 1e-12)
  expect_equal(unname(interval$upr.ci["x", "y"]), expected[2], tolerance = 1e-12)
  expect_identical(
    unname(diagnostics$n_boot_success["x", "y"]),
    as.integer(sum(is.finite(bootstrap_values)))
  )
  expect_identical(dim(interval$lwr.ci), dim(first))
  expect_identical(dimnames(interval$lwr.ci), dimnames(first))
  expect_identical(dim(interval$upr.ci), dim(first))
  expect_identical(dimnames(interval$upr.ci), dimnames(first))
  expect_equal(unname(diag(interval$lwr.ci)), c(1, 1), tolerance = 0)
  expect_equal(unname(diag(interval$upr.ci)), c(1, 1), tolerance = 0)
  expect_true(all(interval$lwr.ci >= -1, na.rm = TRUE))
  expect_true(all(interval$upr.ci <= 1, na.rm = TRUE))
})

test_that("robust_ccc bootstrap honours pairwise missingness and fit failures", {
  x <- c(-2, -1.4, -0.9, -0.3, 0.2, 0.8, 1.3, 1.9, NA, NA)
  y <- c(-1.8, -1.2, -1.0, -0.1, 0.3, 0.7, 1.5, 2.0, 2.5, 3.0)
  z <- c(NA, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3)
  fit <- robust_ccc(
    cbind(x = x, y = y, z = z),
    na_method = "pairwise",
    ci = TRUE,
    conf_level = 0.90,
    n_boot = 12L,
    seed = 202L
  )
  expect_identical(attr(fit, "diagnostics")$n_complete["x", "y"], 8L)
  expect_no_error(ci(fit))

  pair_ok <- is.finite(x) & is.finite(y)
  set.seed(202L)
  pair_boot <- vapply(
    seq_len(12L),
    function(b) {
      index <- sample.int(sum(pair_ok), size = sum(pair_ok), replace = TRUE)
      robust_ccc_reference(x[pair_ok][index], y[pair_ok][index], seed = NULL)
    },
    numeric(1L)
  )
  expected_pair_ci <- robust_ccc_percentile_reference(pair_boot, 0.90)
  expect_equal(unname(ci(fit)$lwr.ci["x", "y"]), expected_pair_ci[1], tolerance = 1e-12)
  expect_equal(unname(ci(fit)$upr.ci["x", "y"]), expected_pair_ci[2], tolerance = 1e-12)

  small <- robust_ccc(
    cbind(a = c(0, 1, 0, 2), b = c(0, 0, 2, 1)),
    ci = TRUE,
    n_boot = 30L,
    seed = 303L
  )
  successes <- attr(small, "diagnostics")$n_boot_success["a", "b"]
  expect_true(successes >= 0L && successes < 30L)
  expect_no_error(ci(small))
})

test_that("robust_ccc output representations filter the same seeded matrix", {
  set.seed(44)
  X <- cbind(
    a = rnorm(30),
    b = rnorm(30),
    c = rnorm(30),
    d = rnorm(30)
  )
  X[, 2] <- 0.8 * X[, 1] + 0.2 * X[, 2]
  threshold <- 0.25

  dense <- robust_ccc(X, seed = 77L)
  sparse <- robust_ccc(
    X,
    seed = 77L,
    output = "sparse",
    threshold = threshold,
    diag = FALSE
  )
  edge <- robust_ccc(
    X,
    seed = 77L,
    output = "edge_list",
    threshold = threshold,
    diag = FALSE
  )

  expected_sparse <- matrix(0, ncol(X), ncol(X), dimnames = dimnames(dense))
  index <- upper.tri(dense, diag = FALSE)
  keep <- index & is.finite(dense) & abs(dense) >= threshold
  expected_sparse[keep] <- dense[keep]
  expected_sparse <- expected_sparse + t(expected_sparse)

  expected_edge <- which(keep, arr.ind = TRUE)
  expected_edge <- data.frame(
    row = rownames(dense)[expected_edge[, 1]],
    col = colnames(dense)[expected_edge[, 2]],
    value = as.numeric(dense[keep]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  expect_equal(as.matrix(sparse), expected_sparse, tolerance = 0)
  expect_equal(
    as.data.frame(edge, stringsAsFactors = FALSE),
    expected_edge,
    tolerance = 0
  )
  expect_identical(attr(sparse, "alpha"), 0.75)
  expect_identical(attr(edge, "mcd_estimate"), "reweighted")
})

test_that("robust_ccc inherits standard accessors and representation dispatch", {
  skip_if_not_installed("ggplot2")
  set.seed(55)
  X <- cbind(a = rnorm(24), b = rnorm(24), c = rnorm(24))
  X[, 2] <- X[, 1] + rnorm(24, sd = 0.3)

  fit <- robust_ccc(X, ci = TRUE, n_boot = 8L, seed = 808L)
  estimate_matrix <- estimate(fit)

  expect_true(is.matrix(estimate_matrix))
  expect_false(inherits(estimate_matrix, "corr_result"))
  expect_equal(coef(fit), estimate_matrix, tolerance = 0)
  expect_identical(ci(fit), attr(fit, "ci"))
  expect_true(all(c("item1", "item2", "estimate", "lwr", "upr") %in% names(tidy(fit))))
  expect_true(all(c("item1", "item2", "lwr", "upr") %in% names(confint(fit))))

  summary_fit <- summary(fit)
  expect_s3_class(summary_fit, "summary.corr_result")
  expect_s3_class(summary_fit, "summary.corr_matrix")
  expect_s3_class(plot(fit, show_value = FALSE), "ggplot")

  sparse <- robust_ccc(X, seed = 808L, output = "sparse", threshold = 0.2, diag = FALSE)
  edge <- robust_ccc(X, seed = 808L, output = "edge_list", threshold = 0.2, diag = FALSE)
  expect_s3_class(summary(sparse), "summary.corr_sparse")
  expect_s3_class(summary(edge), "summary.corr_edge_list")
  expect_s3_class(plot(sparse, show_value = FALSE), "ggplot")
  expect_s3_class(plot(edge, show_value = FALSE), "ggplot")

  expect_null(getS3method("summary", "robust_ccc", optional = TRUE))
  expect_null(getS3method("plot", "robust_ccc", optional = TRUE))
  expect_null(getS3method("estimate", "robust_ccc", optional = TRUE))
  expect_null(getS3method("coef", "robust_ccc", optional = TRUE))
  expect_null(getS3method("ci", "robust_ccc", optional = TRUE))
  expect_null(getS3method("tidy", "robust_ccc", optional = TRUE))
})

test_that("robust_ccc validates the MCD search count without moving existing arguments", {
  expect_identical(tail(names(formals(robust_ccc)), 1L), "mcd_nsamp")
  expect_identical(formals(robust_ccc)$mcd_nsamp, 500L)
  set.seed(616)
  X <- matrix(rnorm(80), ncol = 2)
  expect_identical(robust_ccc(X, seed = 43L),
                   robust_ccc(X, seed = 43L, mcd_nsamp = 500L))

  invalid <- list(NULL, numeric(), NA_real_, NaN, Inf, -Inf, 0, -1, 1.5,
                  c(10, 20), TRUE, "500", "best", "exact", "deterministic",
                  .Machine$integer.max + 1)
  for (value in invalid) {
    expect_error(robust_ccc(X, mcd_nsamp = value),
                 class = "matrixCorr_arg_error")
  }
})

test_that("robust_ccc matches direct MCD with a user selected search count", {
  set.seed(782)
  X <- cbind(a = rnorm(70), b = rnorm(70))
  X[1:3, ] <- cbind(c(9, -10, 11), c(-8, 12, -9))

  for (count in c(1L, 41L, 1000L)) {
    ref <- robust_ccc_reference(X[, 1], X[, 2], seed = 917L, nsamp = count)
    fit <- robust_ccc(X, seed = 917L, mcd_nsamp = count)
    expect_equal(unname(fit[1, 2]), unname(ref), tolerance = 1e-12)
    expect_identical(attr(fit, "mcd_nsamp"), count)
    expect_identical(attr(fit, "mcd_estimate"), "reweighted")
    expect_identical(fit, robust_ccc(X, seed = 917L, mcd_nsamp = count))
  }
})

test_that("robust_ccc forwards the count and fixed reweighting to every MCD fit", {
  original <- robustbase::covMcd
  counts <- integer()
  raw_flags <- logical()
  testthat::local_mocked_bindings(
    covMcd = function(x, alpha, raw.only, nsamp, ...) {
      counts <<- c(counts, nsamp)
      raw_flags <<- c(raw_flags, raw.only)
      original(x, alpha = alpha, raw.only = raw.only, nsamp = nsamp, ...)
    },
    .package = "robustbase"
  )
  set.seed(825)
  X <- matrix(rnorm(120), ncol = 3)
  robust_ccc(X, ci = TRUE, n_boot = 6L, seed = 823L, mcd_nsamp = 73L)

  expected_fits <- choose(3, 2) * (1L + 6L) + 3L
  expect_identical(counts, rep(73L, expected_fits))
  expect_identical(raw_flags, rep(FALSE, expected_fits))
})

test_that("custom MCD counts preserve NA handling bootstrap and output equivalence", {
  set.seed(983)
  X <- matrix(rnorm(120), ncol = 3, dimnames = list(NULL, c("a", "b", "c")))
  count <- 41L
  seed <- 211L
  n_boot <- 10L
  for (mode in c("error", "complete", "pairwise")) {
    Z <- X
    if (mode != "error") {
      Z[1, 1] <- NA
      Z[2, 2] <- Inf
      Z[3, 3] <- NaN
      Z[4, 1] <- -Inf
    }
    rows <- if (mode == "complete") apply(is.finite(Z), 1L, all) else
      is.finite(Z[, 1]) & is.finite(Z[, 2])
    x <- Z[rows, 1]
    y <- Z[rows, 2]
    ref <- robust_ccc_reference(x, y, seed = seed, nsamp = count)

    set.seed(seed)
    boot <- vapply(seq_len(n_boot), function(b) {
      idx <- sample.int(length(x), length(x), replace = TRUE)
      robust_ccc_reference(x[idx], y[idx], nsamp = count)
    }, numeric(1L))
    limits <- robust_ccc_percentile_reference(boot, 0.95)

    set.seed(273L)
    initial_rng <- .Random.seed
    dense <- robust_ccc(Z, na_method = mode, ci = TRUE, n_boot = n_boot,
                        seed = seed, mcd_nsamp = count)
    expect_identical(.Random.seed, initial_rng)
    expect_equal(unname(dense[1, 2]), unname(ref), tolerance = 1e-12)
    expect_equal(unname(ci(dense)$lwr.ci[1, 2]), limits[1], tolerance = 1e-12)
    expect_equal(unname(ci(dense)$upr.ci[1, 2]), limits[2], tolerance = 1e-12)
    expect_identical(attr(dense, "diagnostics")$n_boot_success[1, 2],
                     as.integer(sum(is.finite(boot))))
    expect_identical(dense, robust_ccc(Z, na_method = mode, ci = TRUE,
      n_boot = n_boot, seed = seed, mcd_nsamp = count))

    for (out in c("sparse", "edge_list")) {
      converted <- robust_ccc(Z, na_method = mode, ci = TRUE, n_boot = n_boot,
        seed = seed, mcd_nsamp = count, output = out, threshold = 0.2, diag = FALSE)
      expect_identical(attr(converted, "mcd_nsamp"), count)
      expect_identical(attr(converted, "ci"), attr(dense, "ci"))
      expect_identical(class(converted), class(robust_ccc(
        X, seed = seed, output = out, threshold = 0.2, diag = FALSE)))
      keep <- upper.tri(dense) & is.finite(dense) & abs(dense) >= 0.2
      if (out == "sparse") {
        expected <- matrix(0, ncol(Z), ncol(Z), dimnames = dimnames(dense))
        expected[keep] <- dense[keep]
        expect_equal(as.matrix(converted), expected + t(expected), tolerance = 0)
      } else {
        expect_identical(converted$value, as.numeric(dense[keep]))
      }
    }
  }
})
