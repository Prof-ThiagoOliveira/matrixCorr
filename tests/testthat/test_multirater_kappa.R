manual_multirater_kappa_counts <- function(counts, method = c("fleiss", "randolph")) {
  method <- match.arg(method)
  counts <- as.matrix(counts)
  n_items <- nrow(counts)
  n_categories <- ncol(counts)
  if (n_items < 2L || n_categories < 2L) {
    return(list(kappa = NA_real_, observed = NA_real_, expected = NA_real_))
  }

  n_i <- rowSums(counts)
  item_agreement <- rowSums(counts * (counts - 1)) / (n_i * (n_i - 1))
  observed <- mean(item_agreement)
  expected <- if (identical(method, "fleiss")) {
    p_j <- colSums(counts) / sum(counts)
    sum(p_j^2)
  } else {
    1 / n_categories
  }
  if (!is.finite(1 - expected) || (1 - expected) <= 0) {
    kappa <- NA_real_
  } else {
    kappa <- (observed - expected) / (1 - expected)
  }

  list(
    kappa = kappa,
    observed = observed,
    expected = expected,
    item_agreement = item_agreement,
    n_i = n_i,
    category_counts = colSums(counts),
    category_proportions = colSums(counts) / sum(counts)
  )
}

manual_multirater_kappa_exact_ratings <- function(ratings) {
  ratings <- as.matrix(ratings)
  ratings <- ratings[stats::complete.cases(ratings), , drop = FALSE]
  n_items <- nrow(ratings)
  n_raters <- ncol(ratings)
  lev <- unique(as.character(ratings))
  counts <- vapply(seq_len(n_items), function(i) {
    tab <- table(factor(ratings[i, ], levels = lev))
    as.numeric(tab)
  }, numeric(length(lev)))
  counts <- t(counts)
  observed <- mean(rowSums(counts * (counts - 1)) / (n_raters * (n_raters - 1)))
  chance_base <- sum(colSums(counts)^2) / (n_items * n_raters)^2
  rtab <- vapply(seq_len(n_raters), function(j) {
    tab <- table(factor(ratings[, j], levels = lev))
    as.numeric(tab) / n_items
  }, numeric(length(lev)))
  chance_exact <- chance_base - sum(apply(rtab, 1, stats::var)) / n_raters
  if (!is.finite(1 - chance_exact) || (1 - chance_exact) <= 0) {
    return(NA_real_)
  }
  (observed - chance_exact) / (1 - chance_exact)
}

counts_to_ratings <- function(counts, categories = colnames(counts)) {
  counts <- as.matrix(counts)
  max_raters <- max(rowSums(counts))
  out <- matrix(NA_character_, nrow(counts), max_raters)
  for (i in seq_len(nrow(counts))) {
    labs <- rep(categories, times = counts[i, ])
    out[i, seq_along(labs)] <- labs
  }
  out
}

collapse_counts_for_category_r <- function(counts, j) {
  counts <- as.matrix(counts)
  cbind(counts[, j], rowSums(counts) - counts[, j])
}

test_that("multirater_kappa reproduces the standard Fleiss count example", {
  counts <- rbind(
    c(0, 0, 0, 0, 14),
    c(0, 2, 6, 4, 2),
    c(0, 0, 3, 5, 6),
    c(0, 3, 9, 2, 0),
    c(2, 2, 8, 1, 1),
    c(7, 7, 0, 0, 0),
    c(3, 2, 6, 3, 0),
    c(2, 5, 3, 2, 2),
    c(6, 5, 2, 1, 0),
    c(0, 2, 2, 3, 7)
  )
  colnames(counts) <- paste0("c", 1:5)

  fit <- multirater_kappa(counts, input = "counts", method = "fleiss")
  expect_s3_class(fit, c("multirater_kappa", "agreement_result", "data.frame"))
  expect_equal(fit$kappa[[1L]], 0.2099307, tolerance = 1e-6)
})

test_that("multirater_kappa ratings and counts inputs agree", {
  counts <- rbind(
    c(0, 0, 0, 0, 14),
    c(0, 2, 6, 4, 2),
    c(0, 0, 3, 5, 6),
    c(0, 3, 9, 2, 0),
    c(2, 2, 8, 1, 1),
    c(7, 7, 0, 0, 0),
    c(3, 2, 6, 3, 0),
    c(2, 5, 3, 2, 2),
    c(6, 5, 2, 1, 0),
    c(0, 2, 2, 3, 7)
  )
  colnames(counts) <- paste0("c", 1:5)
  ratings <- counts_to_ratings(counts, categories = colnames(counts))

  fit_counts <- multirater_kappa(counts, input = "counts")
  fit_ratings <- multirater_kappa(ratings, input = "ratings")
  expect_equal(fit_ratings$kappa[[1L]], fit_counts$kappa[[1L]], tolerance = 1e-12)
})

test_that("multirater_kappa reproduces Randolph's free-marginal kappa", {
  counts <- rbind(
    c(3, 0, 0),
    c(2, 1, 0),
    c(0, 3, 0),
    c(0, 1, 2)
  )
  manual <- manual_multirater_kappa_counts(counts, method = "randolph")
  fit <- multirater_kappa(counts, input = "counts", method = "randolph")

  expect_equal(fit$observed_agreement[[1L]], manual$observed, tolerance = 1e-12)
  expect_equal(fit$expected_agreement[[1L]], manual$expected, tolerance = 1e-12)
  expect_equal(fit$kappa[[1L]], manual$kappa, tolerance = 1e-12)
})

test_that("multirater_kappa exact Fleiss matches the irr-compatible exact formula", {
  ratings <- data.frame(
    r1 = c("A", "A", "B", "B", "C", "A", "B", "C"),
    r2 = c("A", "B", "B", "B", "C", "A", "B", "C"),
    r3 = c("A", "A", "B", "C", "C", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  fit_std <- multirater_kappa(ratings, method = "fleiss")
  fit_exact <- multirater_kappa(ratings, method = "fleiss", exact = TRUE)
  expected_exact <- manual_multirater_kappa_exact_ratings(ratings)

  expect_equal(fit_exact$kappa[[1L]], expected_exact, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(fit_std$kappa[[1L]], fit_exact$kappa[[1L]], tolerance = 1e-12)))
  expect_true(isTRUE(attr(fit_exact, "exact")))
})

test_that("multirater_kappa balanced 3-rater regression examples match known estimate and SE", {
  set.seed(20260514)
  n <- 250
  latent <- rnorm(n)
  mk <- function(s) factor(
    cut(latent + rnorm(n, sd = s), breaks = c(-Inf, -0.5, 0.5, Inf), labels = c("A", "B", "C")),
    levels = c("A", "B", "C")
  )
  dat <- data.frame(r1 = mk(0.8), r2 = mk(0.9), r3 = mk(0.85))

  fit_std <- multirater_kappa(dat, method = "fleiss", p_value = TRUE)
  fit_exact <- multirater_kappa(
    dat,
    method = "fleiss",
    exact = TRUE,
    p_value = TRUE,
    se_method = "jackknife"
  )

  expect_equal(fit_std$kappa[[1L]], 0.3159716, tolerance = 1e-7)
  expect_equal(fit_std$se[[1L]], 0.02586107388050949, tolerance = 1e-12)
  expect_equal(fit_std$statistic[[1L]], 12.2180397, tolerance = 1e-7)
  expect_identical(attr(fit_std, "se.method"), "asymptotic")

  expect_equal(fit_exact$kappa[[1L]], 0.3163428, tolerance = 1e-7)
  expect_equal(fit_exact$se[[1L]], 0.03317357033930296, tolerance = 1e-12)
  expect_equal(fit_exact$statistic[[1L]], 9.5359894, tolerance = 1e-7)
  expect_identical(attr(fit_exact, "se.method"), "jackknife")
})

test_that("multirater_kappa resolves categories for character, factor, and supplied levels", {
  char_ratings <- data.frame(
    r1 = c("low", "mid", "high", "low"),
    r2 = c("low", "high", "high", "mid"),
    r3 = c("low", "mid", "high", "mid"),
    stringsAsFactors = FALSE
  )
  fit_char <- multirater_kappa(char_ratings)
  expect_equal(attr(fit_char, "categories"), c("low", "mid", "high"))

  fac_ratings <- data.frame(
    r1 = factor(c("A", "B", "A", "C"), levels = c("A", "B", "C", "D")),
    r2 = factor(c("A", "B", "C", "C"), levels = c("A", "B", "C", "D")),
    r3 = factor(c("A", "B", "A", "D"), levels = c("A", "B", "D", "C"))
  )
  fit_fac <- multirater_kappa(fac_ratings)
  expect_equal(attr(fit_fac, "categories"), c("A", "B", "C", "D"))

  fit_unobs <- multirater_kappa(char_ratings, levels = c("low", "mid", "high", "very_high"))
  expect_equal(attr(fit_unobs, "categories"), c("low", "mid", "high", "very_high"))

  expect_error(
    multirater_kappa(char_ratings, levels = c("low", "mid")),
    "must contain every non-missing observed category"
  )
})

test_that("multirater_kappa requires explicit levels for mixed factor and non-factor ratings", {
  mixed <- data.frame(
    r1 = factor(c("mid", "low", "high"), levels = c("low", "mid", "high")),
    r2 = c("mid", "low", "high"),
    r3 = c("low", "low", "high"),
    stringsAsFactors = FALSE
  )

  expect_error(
    multirater_kappa(mixed),
    "factor and non-factor"
  )

  fit <- multirater_kappa(mixed, levels = c("low", "mid", "high"))
  expect_identical(attr(fit, "categories"), c("low", "mid", "high"))
})

test_that("multirater_kappa missing-data modes follow the documented rules", {
  raters <- data.frame(
    r1 = c("A", "A", "B", NA, "C", "A"),
    r2 = c("A", "B", "B", "B", NA, "A"),
    r3 = c("A", "A", NA, "B", "C", "A"),
    r4 = c("A", NA, "B", "B", "C", "A"),
    stringsAsFactors = FALSE
  )

  expect_error(multirater_kappa(raters, na_method = "error"), "contains missing ratings")

  fit_complete <- multirater_kappa(raters, na_method = "complete")
  keep_complete <- stats::complete.cases(raters)
  expect_equal(fit_complete$n_items[[1L]], sum(keep_complete))

  fit_available <- multirater_kappa(raters, na_method = "available", min_raters = 3L)
  keep_available <- rowSums(!is.na(raters)) >= 3L
  expect_equal(fit_available$n_items[[1L]], sum(keep_available))
  expect_true(!is.null(attr(fit_available, "diagnostics")$n_raters_by_item))

  expect_error(
    multirater_kappa(raters, na_method = "pairwise"),
    "panel-level multi-rater kappa"
  )
  expect_error(
    multirater_kappa(raters, na_method = "available", exact = TRUE),
    "must be FALSE when"
  )
  fit_available_auto <- multirater_kappa(raters, na_method = "available", p_value = TRUE)
  expect_identical(attr(fit_available_auto, "se.method"), "jackknife")
  expect_true(is.finite(fit_available_auto$se[[1L]]) || is.na(fit_available_auto$se[[1L]]))
  expect_error(
    multirater_kappa(raters, na_method = "available", p_value = TRUE, se_method = "asymptotic"),
    "asymptotic"
  )
  fit_available_jk <- multirater_kappa(raters, na_method = "available", p_value = TRUE, se_method = "jackknife")
  expect_identical(attr(fit_available_jk, "se.method"), "jackknife")
})

test_that("multirater_kappa counts input validates counts and category labels explicitly", {
  counts <- matrix(
    c(
      2, 1, 0,
      0, 3, 0,
      0, 0, 1,
      1, 1, 1
    ),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("A", "B", "C"))
  )

  fit <- multirater_kappa(counts, input = "counts", levels = c("low", "mid", "high"), min_raters = 2L)
  expect_identical(attr(fit, "categories"), c("low", "mid", "high"))
  expect_equal(fit$n_items[[1L]], 3L)
  expect_equal(attr(fit, "diagnostics")$row_sums, c(3, 3, 1, 3))

  expect_error(
    multirater_kappa(counts, input = "counts", levels = c("A", "B")),
    "length equal to ncol"
  )
  expect_error(
    multirater_kappa(replace(counts, 1, NA), input = "counts"),
    "missing counts"
  )
  expect_error(
    multirater_kappa(replace(counts, 1, -1), input = "counts"),
    "non-negative"
  )
  expect_error(
    multirater_kappa(replace(counts, 1, 1.5), input = "counts"),
    "integer-like"
  )
})

test_that("multirater_kappa validates inference-only arguments only when inference is requested", {
  raters <- data.frame(
    r1 = c("A", "A", "B", "B"),
    r2 = c("A", "B", "B", "B"),
    r3 = c("A", "A", "B", "A"),
    stringsAsFactors = FALSE
  )

  fit <- multirater_kappa(raters, exact = TRUE, se_method = "not-a-method", conf_level = NA_real_)
  expect_true(isTRUE(attr(fit, "exact")))
  expect_identical(attr(fit, "se.method"), "none")

  expect_error(
    multirater_kappa(raters, p_value = TRUE, se_method = "not-a-method"),
    class = "rlang_error"
  )
  expect_error(
    multirater_kappa(raters, p_value = TRUE, conf_level = NA_real_),
    class = "matrixCorr_arg_error"
  )
})

test_that("multirater_kappa records unbalanced available-rater data", {
  raters <- data.frame(
    r1 = c("A", "A", "B", "C", "A"),
    r2 = c("A", "B", "B", NA, "A"),
    r3 = c("A", NA, "B", NA, "C"),
    r4 = c("A", "B", NA, NA, "A"),
    stringsAsFactors = FALSE
  )

  fit <- multirater_kappa(raters, na_method = "available", min_raters = 2L)
  expect_false(fit$balanced[[1L]])
  expect_equal(fit$n_raters_min[[1L]], 3L)
  expect_equal(fit$n_raters_max[[1L]], 4L)
})

test_that("multirater_kappa handles degenerate validation and denominator cases", {
  expect_error(
    multirater_kappa(matrix(c("A", "A"), ncol = 2), input = "ratings"),
    "at least two categories"
  )

  expect_error(
    multirater_kappa(matrix(c(1, 1, 1, 1), ncol = 1), input = "counts"),
    "at least two category columns"
  )
  expect_error(
    multirater_kappa(matrix(c(1, 1, 1, 1), ncol = 2), input = "counts", exact = TRUE),
    "must be FALSE when"
  )
  expect_error(
    multirater_kappa(data.frame(r1 = c("A", "B"), r2 = c("A", "B")), method = "randolph", exact = TRUE),
    "must be FALSE unless"
  )

  raters <- data.frame(
    r1 = rep("A", 4),
    r2 = rep("A", 4),
    r3 = rep("A", 4),
    stringsAsFactors = FALSE
  )
  fit <- multirater_kappa(raters, levels = c("A", "B"))
  expect_true(is.na(fit$kappa[[1L]]))

  raters_miss <- data.frame(
    r1 = c("A", NA, "B"),
    r2 = c("A", NA, "B"),
    r3 = c(NA, NA, "B"),
    stringsAsFactors = FALSE
  )
  fit_available <- multirater_kappa(raters_miss, na_method = "available", min_raters = 2L)
  expect_equal(fit_available$n_items[[1L]], 2L)
})

test_that("multirater_kappa attaches by-category binary-collapsed results", {
  counts <- rbind(
    c(3, 0, 0),
    c(2, 1, 0),
    c(0, 3, 0),
    c(0, 1, 2)
  )
  colnames(counts) <- c("A", "B", "C")
  fit <- multirater_kappa(counts, input = "counts", by_category = TRUE)
  by_cat <- attr(fit, "by_category")

  expect_true(is.data.frame(by_cat))
  expect_equal(nrow(by_cat), ncol(counts))
  for (j in seq_len(ncol(counts))) {
    manual <- manual_multirater_kappa_counts(collapse_counts_for_category_r(counts, j), method = "fleiss")
    expect_equal(by_cat$kappa[[j]], manual$kappa, tolerance = 1e-12)
  }
})

test_that("multirater_kappa jackknife inference behaves as documented", {
  counts <- rbind(
    c(0, 0, 0, 0, 14),
    c(0, 2, 6, 4, 2),
    c(0, 0, 3, 5, 6),
    c(0, 3, 9, 2, 0),
    c(2, 2, 8, 1, 1),
    c(7, 7, 0, 0, 0),
    c(3, 2, 6, 3, 0),
    c(2, 5, 3, 2, 2),
    c(6, 5, 2, 1, 0),
    c(0, 2, 2, 3, 7)
  )
  fit <- multirater_kappa(counts, input = "counts", ci = TRUE, p_value = TRUE)
  expect_true(is.finite(fit$se[[1L]]) || is.na(fit$se[[1L]]))
  expect_true(is.finite(fit$lwr.ci[[1L]]))
  expect_true(is.finite(fit$upr.ci[[1L]]))
  expect_true(fit$lwr.ci[[1L]] >= -1 && fit$upr.ci[[1L]] <= 1)
  expect_identical(attr(fit, "se.method"), "asymptotic")

  fit_jk <- multirater_kappa(counts, input = "counts", ci = TRUE, p_value = TRUE, se_method = "jackknife")
  expect_identical(attr(fit_jk, "se.method"), "jackknife")
  expect_true(is.finite(fit_jk$se[[1L]]) || is.na(fit_jk$se[[1L]]))

  expect_error(
    multirater_kappa(counts, input = "counts", ci = TRUE, se_method = "none"),
    "must not be 'none'"
  )
  fit_exact_auto <- multirater_kappa(
    data.frame(r1 = c("A", "A", "B"), r2 = c("A", "B", "B"), r3 = c("A", "A", "B")),
    exact = TRUE,
    p_value = TRUE
  )
  expect_identical(attr(fit_exact_auto, "se.method"), "jackknife")
  expect_error(
    multirater_kappa(
      data.frame(r1 = c("A", "A", "B"), r2 = c("A", "B", "B"), r3 = c("A", "A", "B")),
      exact = TRUE,
      p_value = TRUE,
      se_method = "asymptotic"
    ),
    "asymptotic"
  )
  expect_error(
    multirater_kappa(counts, input = "counts", exact = TRUE, p_value = TRUE),
    "must be FALSE when"
  )
  expect_error(
    multirater_kappa(
      data.frame(r1 = c("A", "A", "B"), r2 = c("A", "B", "B"), r3 = c("A", "A", "B")),
      method = "randolph",
      exact = TRUE,
      p_value = TRUE
    ),
    "must be FALSE unless"
  )
})

test_that("multirater_kappa defaults to jackknife when asymptotic inference is unavailable", {
  raters <- data.frame(
    r1 = c("A", "A", "B", "C", "A", "B"),
    r2 = c("A", "B", "B", "C", "A", "B"),
    r3 = c("A", "A", "B", "B", "A", "C"),
    stringsAsFactors = FALSE
  )

  fit_randolph <- multirater_kappa(raters, method = "randolph", ci = TRUE, p_value = TRUE)
  expect_identical(attr(fit_randolph, "se.method"), "jackknife")
  expect_true(is.finite(fit_randolph$se[[1L]]) || is.na(fit_randolph$se[[1L]]))

  expect_error(
    multirater_kappa(raters, method = "randolph", ci = TRUE, p_value = TRUE, se_method = "asymptotic"),
    "asymptotic"
  )
})

test_that("multirater_kappa S3 methods return the expected object types", {
  skip_if_not_installed("ggplot2")

  raters <- data.frame(
    r1 = c("A", "A", "B", "C", "A", "B"),
    r2 = c("A", "B", "B", "C", "A", "B"),
    r3 = c("A", "A", "B", "B", "A", "C"),
    stringsAsFactors = FALSE
  )
  fit <- multirater_kappa(raters, by_category = TRUE)
  expect_s3_class(fit, c("multirater_kappa", "agreement_result", "data.frame"))
  expect_true(length(capture.output(print(fit))) > 0L)

  sm <- summary(fit)
  expect_s3_class(sm, "summary.multirater_kappa")

  expect_s3_class(plot(fit, type = "item_agreement"), "ggplot")
  expect_s3_class(plot(fit, type = "category_proportion"), "ggplot")
  expect_error(plot(multirater_kappa(raters), type = "by_category"), "by_category")
  expect_s3_class(plot(fit, type = "by_category"), "ggplot")
})

test_that("multirater_kappa C++ ratings and counts backends are consistent", {
  counts <- matrix(
    c(
      2, 1, 0,
      0, 3, 0,
      1, 1, 1,
      0, 1, 2
    ),
    ncol = 3,
    byrow = TRUE
  )
  ratings <- counts_to_ratings(counts, categories = c("A", "B", "C"))
  encoded <- matrix(match(ratings, c("A", "B", "C")), nrow = nrow(ratings))

  fit_counts_1 <- multirater_kappa_counts_cpp(counts, method_code = 1L, n_threads = 1L)
  fit_counts_2 <- multirater_kappa_counts_cpp(counts, method_code = 1L, n_threads = 2L)
  fit_ratings <- multirater_kappa_ratings_cpp(encoded, n_levels = 3L, method_code = 1L, na_code = 1L, n_threads = 1L)

  expect_equal(fit_counts_1$estimate, fit_ratings$estimate, tolerance = 1e-12)
  expect_equal(fit_counts_1$estimate, fit_counts_2$estimate, tolerance = 1e-12)

  fit_ratings_exact <- multirater_kappa_ratings_cpp(
    encoded,
    n_levels = 3L,
    method_code = 1L,
    na_code = 1L,
    n_threads = 1L,
    exact = TRUE
  )
  expect_equal(
    fit_ratings_exact$estimate,
    manual_multirater_kappa_exact_ratings(ratings),
    tolerance = 1e-12
  )
})
