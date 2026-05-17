upper_median_bandwidth_ref <- function(x) {
  D <- abs(outer(x, x, "-"))
  d <- sort(as.numeric(D[upper.tri(D, diag = FALSE)]))

  if (!length(d)) {
    return(0.001)
  }

  sigma <- d[floor(length(d) / 2) + 1L] / sqrt(2)

  if (!is.finite(sigma) || sigma <= 0) {
    return(0.001)
  }

  sigma
}

silverman_bandwidth_ref <- function(x) {
  n <- length(x)
  s <- stats::sd(x)
  if (!is.finite(s) || s <= 1e-12) {
    return(1.0)
  }
  max(1.06 * s * n^(-1 / 5), 1e-12)
}

scott_bandwidth_ref <- function(x) {
  n <- length(x)
  s <- stats::sd(x)
  if (!is.finite(s) || s <= 1e-12) {
    return(1.0)
  }
  max(s * n^(-1 / 5), 1e-12)
}

bandwidth_ref <- function(x, bandwidth) {
  switch(
    bandwidth,
    median = upper_median_bandwidth_ref(x),
    silverman = silverman_bandwidth_ref(x),
    scott = scott_bandwidth_ref(x),
    stop("unsupported bandwidth")
  )
}

gram_ref <- function(x,
                     kernel = c("gaussian", "linear", "laplace", "polynomial"),
                     bandwidth = c("median", "silverman", "scott")) {
  kernel <- match.arg(kernel)
  bandwidth <- match.arg(bandwidth)

  if (identical(kernel, "linear")) {
    return(outer(x, x, "*"))
  }

  if (identical(kernel, "polynomial")) {
    return((outer(x, x, "*") + 1)^2)
  }

  bw <- bandwidth_ref(x, bandwidth)
  D <- abs(outer(x, x, "-"))

  if (identical(kernel, "gaussian")) {
    return(exp(-(D^2) / (2 * bw^2)))
  }

  if (identical(kernel, "laplace")) {
    return(exp(-D / bw))
  }

  stop("unsupported kernel")
}

center_ref <- function(K) {
  n <- nrow(K)
  H <- diag(n) - matrix(1 / n, n, n)
  H %*% K %*% H
}

hsic_biased_ref <- function(x, y, kernel, bandwidth) {
  K <- gram_ref(x, kernel = kernel, bandwidth = bandwidth)
  L <- gram_ref(y, kernel = kernel, bandwidth = bandwidth)
  Kc <- center_ref(K)
  Lc <- center_ref(L)
  sum(Kc * Lc) / length(x)^2
}

hsic_kcor_ref <- function(x, y, kernel, bandwidth) {
  K <- gram_ref(x, kernel = kernel, bandwidth = bandwidth)
  L <- gram_ref(y, kernel = kernel, bandwidth = bandwidth)

  Kc <- center_ref(K)
  Lc <- center_ref(L)

  xy <- sum(Kc * Lc)
  xx <- sum(Kc * Kc)
  yy <- sum(Lc * Lc)

  if (
    !is.finite(xy) || !is.finite(xx) || !is.finite(yy) ||
      xx <= 1e-12 || yy <= 1e-12
  ) {
    return(NA_real_)
  }

  out <- xy / sqrt(xx * yy)
  max(0, min(1, out))
}

hsic_unbiased_ref <- function(x, y, kernel, bandwidth) {
  n <- length(x)
  K <- gram_ref(x, kernel = kernel, bandwidth = bandwidth)
  L <- gram_ref(y, kernel = kernel, bandwidth = bandwidth)

  diag(K) <- 0
  diag(L) <- 0

  rK <- rowSums(K)
  rL <- rowSums(L)

  (
    sum(K * L) +
      sum(K) * sum(L) / ((n - 1) * (n - 2)) -
      2 * sum(rK * rL) / (n - 2)
  ) / (n * (n - 3))
}

hsic_math_scenario <- function(name, n) {
  set.seed(
    switch(
      name,
      independent_normal = 101L,
      linear_dependence = 102L,
      quadratic_dependence = 103L,
      duplicate_columns = 104L,
      duplicated_observations = 105L,
      constant_column = 106L
    ) + n
  )

  if (identical(name, "independent_normal")) {
    return(list(x = rnorm(n), y = rnorm(n)))
  }
  if (identical(name, "linear_dependence")) {
    x <- rnorm(n)
    return(list(x = x, y = 0.8 * x + rnorm(n, sd = 0.2)))
  }
  if (identical(name, "quadratic_dependence")) {
    x <- rnorm(n)
    return(list(x = x, y = x^2 + rnorm(n, sd = 0.1)))
  }
  if (identical(name, "duplicate_columns")) {
    x <- rnorm(n)
    return(list(x = x, y = x))
  }
  if (identical(name, "duplicated_observations")) {
    x <- as.numeric(sample(c(-1, 0, 1), n, replace = TRUE))
    return(list(x = x, y = x + rnorm(n, sd = 0.1)))
  }
  if (identical(name, "constant_column")) {
    return(list(x = rnorm(n), y = rep(1, n)))
  }

  stop("unknown scenario")
}

test_that("hsic matches direct biased and normalised mathematical references", {
  kernels <- c("gaussian", "linear", "laplace", "polynomial")
  bandwidths <- c("median", "silverman", "scott")
  scenarios <- c(
    "independent_normal",
    "linear_dependence",
    "quadratic_dependence",
    "duplicate_columns",
    "duplicated_observations",
    "constant_column"
  )

  for (scenario in scenarios) {
    for (n in c(5L, 100L)) {
      dat <- hsic_math_scenario(scenario, n)
      X <- cbind(x = dat$x, y = dat$y)

      for (kernel in kernels) {
        for (bandwidth in bandwidths) {
          raw <- hsic(
            X,
            kernel = kernel,
            bandwidth = bandwidth,
            estimator = "biased",
            normalise = FALSE
          )
          norm <- hsic(
            X,
            kernel = kernel,
            bandwidth = bandwidth,
            estimator = "biased",
            normalise = TRUE
          )

          expect_equal(
            as.numeric(raw[1, 2]),
            hsic_biased_ref(dat$x, dat$y, kernel, bandwidth),
            tolerance = 1e-10,
            info = paste(scenario, n, kernel, bandwidth, "raw")
          )
          expect_equal(
            as.numeric(norm[1, 2]),
            hsic_kcor_ref(dat$x, dat$y, kernel, bandwidth),
            tolerance = 1e-10,
            info = paste(scenario, n, kernel, bandwidth, "normalised")
          )
        }
      }
    }
  }
})

test_that("hsic unbiased estimator matches direct mathematical reference", {
  kernels <- c("gaussian", "linear", "laplace", "polynomial")
  bandwidths <- c("median", "silverman", "scott")
  dat <- hsic_math_scenario("quadratic_dependence", 20L)
  X <- cbind(x = dat$x, y = dat$y)

  for (kernel in kernels) {
    for (bandwidth in bandwidths) {
      raw <- hsic(
        X,
        kernel = kernel,
        bandwidth = bandwidth,
        estimator = "unbiased",
        normalise = FALSE
      )
      expect_equal(
        as.numeric(raw[1, 2]),
        hsic_unbiased_ref(dat$x, dat$y, kernel, bandwidth),
        tolerance = 1e-10,
        info = paste(kernel, bandwidth)
      )
    }
  }
})

test_that("linear and polynomial HSIC ignore bandwidth", {
  dat <- hsic_math_scenario("linear_dependence", 50L)
  X <- cbind(x = dat$x, y = dat$y)

  for (kernel in c("linear", "polynomial")) {
    vals <- vapply(
      c("median", "silverman", "scott"),
      function(bandwidth) {
        as.numeric(hsic(X, kernel = kernel, bandwidth = bandwidth, normalise = FALSE)[1, 2])
      },
      numeric(1)
    )
    expect_equal(unname(vals), rep(vals[[1]], length(vals)), tolerance = 1e-12)
  }
})

test_that("gaussian and laplace HSIC respond to bandwidth", {
  dat <- hsic_math_scenario("linear_dependence", 100L)
  X <- cbind(x = dat$x, y = dat$y)

  for (kernel in c("gaussian", "laplace")) {
    vals <- vapply(
      c("median", "silverman", "scott"),
      function(bandwidth) {
        as.numeric(hsic(X, kernel = kernel, bandwidth = bandwidth, normalise = FALSE)[1, 2])
      },
      numeric(1)
    )
    expect_true(length(unique(round(vals, 14))) > 1)
  }
})

test_that("hsic seed does not alter global R RNG state", {
  set.seed(42)
  before <- .Random.seed
  X <- cbind(x = rnorm(30), y = rnorm(30))
  set.seed(42)
  before <- .Random.seed
  invisible(hsic(X, p_value = TRUE, B = 19L, seed = 10L))
  after <- .Random.seed
  expect_identical(after, before)
})
