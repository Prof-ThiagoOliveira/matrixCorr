ratings_to_counts <- function(x, levels, min_raters = 2L) {
  x <- as.matrix(x)
  keep <- rowSums(!is.na(x)) >= min_raters
  x <- x[keep, , drop = FALSE]
  out <- matrix(0L, nrow = nrow(x), ncol = length(levels))
  colnames(out) <- as.character(levels)
  for (i in seq_len(nrow(x))) {
    xi <- x[i, ]
    xi <- xi[!is.na(xi)]
    out[i, ] <- tabulate(match(as.character(xi), as.character(levels)), nbins = ncol(out))
  }
  out
}

nominal_example <- function() {
  matrix(c(
    1, 2, 3, 3, 2, 1, 4, 1, 2, NA, NA, NA,
    1, 2, 3, 3, 2, 2, 4, 1, 2, 5, NA, 3,
    NA, 3, 3, 3, 2, 3, 4, 2, 2, 5, 1, NA,
    1, 2, 3, 3, 2, 4, 4, 1, 2, 5, 1, NA
  ), nrow = 12, ncol = 4)
}

interval_example <- function() {
  matrix(c(
    10, 10, 12, 11,
    20, 19, 21, 20,
    15, 15, NA, 16,
    30, 31, 29, 30,
    40, 39, 40, 41,
    50, 52, 51, 50,
    80, 78, 79, NA,
    100, 101, 99, 100
  ), byrow = TRUE, ncol = 4)
}

ratio_example <- function() {
  matrix(c(
    1, 1, 2, 1,
    2, 2, 2, 4,
    4, 4, 8, 4,
    8, 8, 8, 16,
    16, 16, 32, NA,
    32, 32, 32, 64,
    64, 64, 64, 64,
    128, 128, NA, 256
  ), byrow = TRUE, ncol = 4)
}

test_that("krippendorff_alpha customary nominal matches classic benchmark", {
  nominal <- nominal_example()
  fit <- krippendorff_alpha(
    nominal,
    level = "nominal",
    method = "customary",
    na_method = "available"
  )

  expect_s3_class(fit, "krippendorff_alpha")
  expect_equal(fit$alpha, 0.7434, tolerance = 5e-4)
  expect_equal(fit$n_items, 11L)
  expect_equal(fit$n_categories, 5L)
})

test_that("krippendorff_alpha smoke-test benchmarks match frozen estimate values", {
  nominal <- nominal_example()
  interval_dat <- interval_example()
  ratio_dat <- ratio_example()

  fit_nominal_customary <- krippendorff_alpha(
    nominal,
    level = "nominal",
    method = "customary",
    na_method = "available"
  )
  fit_nominal_analytical <- krippendorff_alpha(
    nominal,
    level = "nominal",
    method = "analytical",
    na_method = "available",
    ci = TRUE,
    p_value = TRUE
  )
  fit_interval_customary <- krippendorff_alpha(
    interval_dat,
    level = "interval",
    method = "customary",
    na_method = "available"
  )
  fit_interval_analytical <- krippendorff_alpha(
    interval_dat,
    level = "interval",
    method = "analytical",
    na_method = "available",
    ci = TRUE,
    p_value = TRUE
  )
  fit_ratio_customary <- krippendorff_alpha(
    ratio_dat,
    level = "ratio",
    method = "customary",
    na_method = "available"
  )
  fit_ratio_analytical <- krippendorff_alpha(
    ratio_dat,
    level = "ratio",
    method = "analytical",
    na_method = "available",
    ci = TRUE,
    p_value = TRUE
  )

  expect_equal(as.numeric(fit_nominal_customary$alpha[[1L]]), 0.74342105, tolerance = 5e-7)
  expect_equal(as.numeric(fit_nominal_analytical$alpha[[1L]]), 0.75711893, tolerance = 5e-7)
  expect_equal(as.numeric(fit_interval_customary$alpha[[1L]]), 0.99919521, tolerance = 5e-7)
  expect_equal(as.numeric(fit_interval_analytical$alpha[[1L]]), 0.99927285, tolerance = 5e-7)
  expect_equal(as.numeric(fit_ratio_customary$alpha[[1L]]), 0.88950587, tolerance = 5e-7)
  expect_equal(as.numeric(fit_ratio_analytical$alpha[[1L]]), 0.89909781, tolerance = 5e-7)

  expect_equal(as.numeric(fit_nominal_analytical$lwr.ci[[1L]]), 0.23035062, tolerance = 5e-7)
  expect_equal(as.numeric(fit_nominal_analytical$upr.ci[[1L]]), 0.95178774, tolerance = 5e-7)
  expect_equal(as.numeric(fit_interval_analytical$lwr.ci[[1L]]), 0.99625147, tolerance = 5e-7)
  expect_equal(as.numeric(fit_interval_analytical$upr.ci[[1L]]), 0.99985920, tolerance = 5e-7)
  expect_equal(as.numeric(fit_ratio_analytical$lwr.ci[[1L]]), 0.82735250, tolerance = 5e-7)
  expect_equal(as.numeric(fit_ratio_analytical$upr.ci[[1L]]), 0.94246775, tolerance = 5e-7)
})

test_that("ratings input and counts input return identical customary nominal alpha", {
  nominal <- nominal_example()
  lev <- 1:5
  counts <- ratings_to_counts(nominal, levels = lev, min_raters = 2L)

  fit_ratings <- krippendorff_alpha(
    nominal,
    levels = lev,
    input = "ratings",
    level = "nominal",
    method = "customary",
    na_method = "available"
  )
  fit_counts <- krippendorff_alpha(
    counts,
    levels = lev,
    input = "counts",
    level = "nominal",
    method = "customary"
  )

  expect_equal(fit_ratings$alpha, fit_counts$alpha, tolerance = 1e-12)
  expect_equal(fit_ratings$observed_disagreement, fit_counts$observed_disagreement, tolerance = 1e-12)
  expect_equal(fit_ratings$expected_disagreement, fit_counts$expected_disagreement, tolerance = 1e-12)
})

test_that("na_method = error rejects missing ratings", {
  nominal <- nominal_example()
  expect_error(
    krippendorff_alpha(
      nominal,
      level = "nominal",
      method = "customary",
      na_method = "error"
    ),
    "missing ratings"
  )
})

test_that("na_method = complete drops incomplete items", {
  nominal <- nominal_example()
  fit <- krippendorff_alpha(
    nominal,
    level = "nominal",
    method = "customary",
    na_method = "complete"
  )
  expect_equal(fit$n_items, sum(stats::complete.cases(nominal)))
})

test_that("na_method = available keeps rows with at least min_raters ratings", {
  nominal <- nominal_example()
  fit <- krippendorff_alpha(
    nominal,
    level = "nominal",
    method = "customary",
    na_method = "available",
    min_raters = 3L
  )
  expect_equal(fit$n_items, sum(rowSums(!is.na(nominal)) >= 3L))
})

test_that("explicit levels preserve unused categories", {
  x <- data.frame(
    r1 = factor(c("A", "A", "B", "C"), levels = c("A", "B", "C", "D")),
    r2 = factor(c("A", "B", "B", "C"), levels = c("A", "B", "C", "D")),
    r3 = factor(c("A", "A", "B", "C"), levels = c("A", "B", "C", "D"))
  )
  fit <- krippendorff_alpha(
    x,
    levels = c("A", "B", "C", "D"),
    level = "nominal",
    method = "customary",
    na_method = "available"
  )

  expect_equal(fit$n_categories, 4L)
  expect_identical(as.character(attr(fit, "categories", exact = TRUE)), c("A", "B", "C", "D"))
})

test_that("ordered factor + ordinal respects factor order", {
  lev <- c("low", "mid", "high")
  x <- data.frame(
    r1 = ordered(c("low", "mid", "high", "mid", "low"), levels = lev),
    r2 = ordered(c("low", "high", "high", "mid", "low"), levels = lev),
    r3 = ordered(c("mid", "mid", "high", "mid", "low"), levels = lev)
  )

  fit_factor <- krippendorff_alpha(
    x,
    level = "ordinal",
    method = "customary",
    na_method = "available"
  )
  fit_levels <- krippendorff_alpha(
    data.frame(lapply(x, as.character), stringsAsFactors = FALSE),
    levels = lev,
    level = "ordinal",
    method = "customary",
    na_method = "available"
  )

  expect_equal(fit_factor$alpha, fit_levels$alpha, tolerance = 1e-12)
})

test_that("character ordinal data without levels errors", {
  x <- data.frame(
    r1 = c("low", "mid", "high"),
    r2 = c("low", "high", "high"),
    r3 = c("mid", "mid", "high"),
    stringsAsFactors = FALSE
  )
  expect_error(
    krippendorff_alpha(
      x,
      level = "ordinal",
      method = "customary",
      na_method = "available"
    ),
    "must be supplied"
  )
})

test_that("interval alpha requires numeric categories", {
  x <- data.frame(
    r1 = c("low", "mid", "high"),
    r2 = c("low", "high", "high"),
    r3 = c("mid", "mid", "high"),
    stringsAsFactors = FALSE
  )
  expect_error(
    krippendorff_alpha(
      x,
      level = "interval",
      method = "customary",
      na_method = "available"
    ),
    "numeric"
  )
})

test_that("ratio alpha rejects negative categories", {
  x <- data.frame(
    r1 = c(-1, 0, 1, 2),
    r2 = c(-1, 1, 1, 2),
    r3 = c(0, 0, 1, 2)
  )
  expect_error(
    krippendorff_alpha(
      x,
      level = "ratio",
      method = "customary",
      na_method = "available"
    ),
    "negative"
  )
})

test_that("analytical estimator returns finite jackknife CI and p-value", {
  nominal <- nominal_example()
  fit <- krippendorff_alpha(
    nominal,
    level = "nominal",
    method = "analytical",
    na_method = "available",
    ci = TRUE,
    p_value = TRUE
  )

  expect_true(is.finite(fit$alpha))
  expect_true(is.finite(fit$lwr.ci))
  expect_true(is.finite(fit$upr.ci))
  expect_true(is.finite(fit$se))
  expect_true(is.finite(fit$statistic))
  expect_true(is.finite(fit$p_value))
  expect_identical(attr(fit, "se.method", exact = TRUE), "jackknife")
})

test_that("customary alpha rejects p_value = TRUE", {
  nominal <- nominal_example()
  expect_error(
    krippendorff_alpha(
      nominal,
      level = "nominal",
      method = "customary",
      na_method = "available",
      p_value = TRUE
    ),
    "P-values are available only"
  )
})

test_that("customary bootstrap CI is reproducible with fixed seed", {
  nominal <- nominal_example()
  fit1 <- krippendorff_alpha(
    nominal,
    level = "nominal",
    method = "customary",
    na_method = "available",
    ci = TRUE,
    n_boot = 100L,
    seed = 11L
  )
  fit2 <- krippendorff_alpha(
    nominal,
    level = "nominal",
    method = "customary",
    na_method = "available",
    ci = TRUE,
    n_boot = 100L,
    seed = 11L
  )

  expect_equal(fit1$lwr.ci, fit2$lwr.ci, tolerance = 0)
  expect_equal(fit1$upr.ci, fit2$upr.ci, tolerance = 0)
  expect_equal(fit1$se, fit2$se, tolerance = 0)
  expect_identical(attr(fit1, "se.method", exact = TRUE), "bootstrap")
})

test_that("return_matrices attaches counts, coincidence, expected coincidence and delta2", {
  nominal <- nominal_example()
  fit <- krippendorff_alpha(
    nominal,
    level = "nominal",
    method = "customary",
    na_method = "available",
    return_matrices = TRUE
  )
  mats <- attr(fit, "matrices", exact = TRUE)

  expect_true(is.list(mats))
  expect_true(all(c("counts", "coincidence", "expected_coincidence", "delta2") %in% names(mats)))
  expect_true(is.matrix(mats$counts))
  expect_true(is.matrix(mats$coincidence))
  expect_true(is.matrix(mats$expected_coincidence))
  expect_true(is.matrix(mats$delta2))
  expect_equal(ncol(mats$counts), fit$n_categories)
})

test_that("print, summary, and plot methods do not error", {
  skip_if_not_installed("ggplot2")

  nominal <- nominal_example()
  fit <- krippendorff_alpha(
    nominal,
    level = "nominal",
    method = "customary",
    na_method = "available",
    ci = TRUE,
    n_boot = 50L,
    seed = 1L,
    return_matrices = TRUE
  )

  expect_no_error(capture.output(print(fit)))
  sm <- summary(fit)
  expect_s3_class(sm, "summary.krippendorff_alpha")
  expect_no_error(capture.output(print(sm)))
  expect_s3_class(plot(fit, type = "estimate"), "ggplot")
  expect_s3_class(plot(fit, type = "item_disagreement"), "ggplot")
  expect_s3_class(plot(fit, type = "category_proportion"), "ggplot")
  expect_s3_class(plot(fit, type = "coincidence"), "ggplot")
})
