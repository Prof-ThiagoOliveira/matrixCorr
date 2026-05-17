manual_gwet_pair <- function(x, y, levels = NULL, weights = NULL) {
  keep <- stats::complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  if (is.null(levels)) {
    levels <- unique(c(as.character(x), as.character(y)))
  }
  q <- length(levels)
  if (q < 2L || length(x) < 2L) {
    return(list(ac = NA_real_, observed = NA_real_, expected = NA_real_, n = length(x)))
  }

  if (is.null(weights)) {
    weights <- diag(q)
  }
  tab <- table(
    factor(as.character(x), levels = as.character(levels)),
    factor(as.character(y), levels = as.character(levels))
  )
  p <- unclass(tab) / sum(tab)
  pi_k <- (rowSums(p) + colSums(p)) / 2
  observed <- sum(weights * p)
  expected <- sum(weights) * sum(pi_k * (1 - pi_k)) / (q * (q - 1))
  ac <- if (!is.finite(1 - expected) || (1 - expected) <= 0) {
    NA_real_
  } else {
    (observed - expected) / (1 - expected)
  }
  list(ac = ac, observed = observed, expected = expected, n = sum(tab))
}

manual_gwet_counts <- function(counts, weights = NULL) {
  counts <- as.matrix(counts)
  n_items <- nrow(counts)
  q <- ncol(counts)
  if (n_items < 2L || q < 2L) {
    return(list(ac = NA_real_, observed = NA_real_, expected = NA_real_))
  }
  if (is.null(weights)) {
    weights <- diag(q)
  }

  n_i <- rowSums(counts)
  item_agreement <- vapply(seq_len(n_items), function(i) {
    ri <- n_i[[i]]
    if (ri < 2L) {
      return(NA_real_)
    }
    sum(counts[i, ] * (drop(weights %*% counts[i, ]) - 1)) / (ri * (ri - 1))
  }, numeric(1))
  observed <- mean(item_agreement, na.rm = TRUE)
  pi_k <- colMeans(counts / n_i)
  expected <- sum(weights) * sum(pi_k * (1 - pi_k)) / (q * (q - 1))
  ac <- if (!is.finite(1 - expected) || (1 - expected) <= 0) {
    NA_real_
  } else {
    (observed - expected) / (1 - expected)
  }

  list(
    ac = ac,
    observed = observed,
    expected = expected,
    item_agreement = item_agreement,
    category_proportions = colSums(counts) / sum(counts),
    pi_k = pi_k
  )
}

collapse_counts_for_category_ac1 <- function(counts, j) {
  counts <- as.matrix(counts)
  cbind(counts[, j], rowSums(counts) - counts[, j])
}

ratings_to_counts <- function(ratings, levels = NULL) {
  ratings <- as.matrix(ratings)
  if (is.null(levels)) {
    levels <- unique(as.character(ratings[!is.na(ratings)]))
  }
  out <- matrix(0L, nrow(ratings), length(levels), dimnames = list(NULL, levels))
  for (i in seq_len(nrow(ratings))) {
    tab <- table(factor(as.character(ratings[i, !is.na(ratings[i, ])]), levels = levels))
    out[i, ] <- as.integer(tab)
  }
  out
}

test_that("gwet_ac reproduces a manual two-rater binary AC1 example", {
  x <- c(rep("A", 25), rep("B", 25))
  y <- c(rep("A", 20), rep("B", 5), rep("A", 10), rep("B", 15))
  manual <- manual_gwet_pair(x, y)
  fit <- gwet_ac(x, y)

  expect_s3_class(fit, "gwet_ac")
  expect_equal(as.numeric(fit), manual$ac, tolerance = 1e-12)
  expect_equal(attr(fit, "diagnostics")$observed_agreement, manual$observed, tolerance = 1e-12)
  expect_equal(attr(fit, "diagnostics")$expected_agreement, manual$expected, tolerance = 1e-12)
})

test_that("gwet_ac reproduces a manual two-rater three-category AC1 example", {
  tab <- matrix(
    c(
      12, 3, 1,
      4, 10, 2,
      2, 5, 11
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(x = c("A", "B", "C"), y = c("A", "B", "C"))
  )
  x <- rep(rownames(tab), times = rowSums(tab))
  y <- unlist(lapply(seq_len(nrow(tab)), function(i) {
    rep(colnames(tab), times = tab[i, ])
  }), use.names = FALSE)

  manual <- manual_gwet_pair(x, y, levels = c("A", "B", "C"))
  fit <- gwet_ac(
    factor(x, levels = c("A", "B", "C")),
    factor(y, levels = c("A", "B", "C"))
  )

  expect_equal(attr(fit, "diagnostics")$observed_agreement, manual$observed, tolerance = 1e-12)
  expect_equal(attr(fit, "diagnostics")$expected_agreement, manual$expected, tolerance = 1e-12)
  expect_equal(as.numeric(fit), manual$ac, tolerance = 1e-12)
})

test_that("gwet_ac reproduces a manual weighted AC2 example", {
  x <- factor(c("low", "low", "mid", "mid", "high", "high"), levels = c("low", "mid", "high"))
  y <- factor(c("low", "mid", "mid", "high", "high", "mid"), levels = c("low", "mid", "high"))
  W <- matrix(
    c(
      1.0, 0.5, 0.0,
      0.5, 1.0, 0.5,
      0.0, 0.5, 1.0
    ),
    nrow = 3,
    byrow = TRUE
  )
  manual <- manual_gwet_pair(x, y, levels = levels(x), weights = W)
  fit <- gwet_ac(x, y, weights = W)

  expect_equal(as.numeric(fit), manual$ac, tolerance = 1e-12)
  expect_identical(attr(fit, "coefficient"), "AC2")
  expect_equal(attr(fit, "diagnostics")$observed_agreement, manual$observed, tolerance = 1e-12)
  expect_equal(attr(fit, "diagnostics")$expected_agreement, manual$expected, tolerance = 1e-12)
})

test_that("gwet_ac rejects ambiguous multi-value weight schemes", {
  x <- c("A", "A", "B", "B")
  y <- c("A", "B", "B", "A")

  expect_error(
    gwet_ac(x, y, weights = c("linear", "quadratic")),
    class = "matrixCorr_arg_error"
  )
  expect_error(
    .mc_resolve_gwet_weights(c("linear", "quadratic"), c("A", "B")),
    class = "matrixCorr_arg_error"
  )
  expect_s3_class(gwet_ac(x, y), "gwet_ac")
})

test_that("gwet_ac ordinal weights respect explicit and implicit category ordering", {
  x <- factor(c("low", "high", "mid", "low"), levels = c("high", "mid", "low"))
  y <- factor(c("high", "high", "mid", "low"), levels = c("high", "mid", "low"))
  fit_factor <- gwet_ac(x, y, weights = "linear")

  expect_identical(attr(fit_factor, "levels"), c("high", "mid", "low"))
  expect_identical(dimnames(attr(fit_factor, "weights"))[[1L]], c("high", "mid", "low"))
  expect_equal(attr(fit_factor, "weights")["high", "low"], 0)
  expect_equal(attr(fit_factor, "weights")["high", "mid"], 0.5)

  fit_explicit <- gwet_ac(
    as.character(x),
    as.character(y),
    weights = "linear",
    levels = c("low", "mid", "high")
  )
  expect_identical(attr(fit_explicit, "levels"), c("low", "mid", "high"))
  expect_identical(dimnames(attr(fit_explicit, "weights"))[[1L]], c("low", "mid", "high"))
  expect_equal(attr(fit_explicit, "weights")["low", "high"], 0)

  counts <- matrix(
    c(2, 0, 1, 0, 1, 2),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("high", "mid", "low"))
  )
  fit_counts <- gwet_ac(counts, input = "counts", weights = "linear")
  expect_identical(attr(fit_counts, "categories"), c("high", "mid", "low"))
  expect_identical(dimnames(attr(fit_counts, "weights"))[[1L]], c("high", "mid", "low"))
})

test_that("gwet_ac differs from Cohen's kappa in a high-prevalence example", {
  x <- c(rep("A", 90), rep("B", 10))
  y <- c(rep("A", 85), rep("B", 5), rep("A", 5), rep("B", 5))

  fit_gwet <- gwet_ac(x, y)
  fit_cohen <- cohen_kappa(x, y)
  manual <- manual_gwet_pair(x, y)

  expect_equal(as.numeric(fit_gwet), manual$ac, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(as.numeric(fit_gwet), as.numeric(fit_cohen), tolerance = 1e-6)))
})

test_that("gwet_ac matrix mode is symmetric and matches scalar pairwise values", {
  raters <- data.frame(
    r1 = c("A", "A", "B", "B", "C", "C", "A", "B"),
    r2 = c("A", "B", "B", "B", "C", "A", "A", "B"),
    r3 = c(1, 1, 2, 2, 3, 3, 1, 2),
    r4 = c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  fit <- gwet_ac(raters)
  mat <- as.matrix(fit)
  expect_true(isSymmetric(mat))
  expect_equal(unname(diag(mat)), rep(1, ncol(raters)))

  for (j in seq_len(ncol(raters) - 1L)) {
    for (k in (j + 1L):ncol(raters)) {
      scalar <- gwet_ac(raters[[j]], raters[[k]])
      expect_equal(mat[j, k], as.numeric(scalar), tolerance = 1e-12)
    }
  }
})

test_that("gwet_ac output modes integrate with corr_result helpers", {
  raters <- data.frame(
    a = c("A", "A", "B", "B", "C", "A", "B", "C"),
    b = c("A", "B", "B", "B", "C", "A", "B", "C"),
    c = c("A", "A", "B", "C", "C", "A", "B", "B"),
    d = c("A", "B", "A", "B", "C", "C", "B", "C"),
    stringsAsFactors = FALSE
  )

  fit_matrix <- gwet_ac(raters, output = "matrix")
  expect_s3_class(fit_matrix, c("corr_matrix", "gwet_ac", "corr_result", "matrix"))

  edge <- gwet_ac(raters, output = "edge_list", threshold = 0.2, diag = FALSE)
  edge_df <- as.data.frame(edge, stringsAsFactors = FALSE)
  expect_true(all(c("row", "col", "value") %in% names(edge_df)))

  sparse <- gwet_ac(raters, output = "sparse", threshold = 0.2, diag = FALSE)
  expect_s4_class(sparse, "sparseMatrix")

  expect_error(
    gwet_ac(raters, input = "ratings", output = "edge_list", threshold = 0.1),
    "must be .*matrix"
  )
})

test_that("gwet_ac missing-data modes behave as documented", {
  raters <- data.frame(
    a = c("A", "A", "B", NA, "C", "A"),
    b = c("A", "B", "B", "B", NA, "A"),
    c = c("A", "A", NA, "B", "C", "A"),
    stringsAsFactors = FALSE
  )

  expect_error(gwet_ac(raters, na_method = "error"), "contains missing values")

  fit_complete <- gwet_ac(raters, na_method = "complete")
  keep <- stats::complete.cases(raters)
  expected_complete <- matrix(sum(keep), ncol(raters), ncol(raters))
  diag(expected_complete) <- rep(sum(keep), ncol(raters))
  dimnames(expected_complete) <- dimnames(as.matrix(fit_complete))
  expect_equal(attr(fit_complete, "diagnostics")$n_complete, expected_complete)

  fit_pairwise <- gwet_ac(raters, na_method = "pairwise")
  n_complete <- attr(fit_pairwise, "diagnostics")$n_complete
  expect_equal(n_complete["a", "b"], sum(stats::complete.cases(raters[c("a", "b")])))
  expect_equal(n_complete["a", "c"], sum(stats::complete.cases(raters[c("a", "c")])))
  expect_equal(n_complete["b", "c"], sum(stats::complete.cases(raters[c("b", "c")])))
})

test_that("gwet_ac attaches asymptotic CI and inference metadata", {
  x <- factor(c("A", "A", "B", "B", "C", "A", "B", "C", "A", "C"))
  y <- factor(c("A", "B", "B", "B", "C", "A", "B", "C", "A", "B"))

  fit <- gwet_ac(x, y, ci = TRUE, p_value = TRUE)
  ci_attr <- attr(fit, "ci")
  inf_attr <- attr(fit, "inference")

  expect_true(is.list(ci_attr))
  expect_identical(ci_attr$ci.method, "asymptotic")
  expect_true(all(is.finite(c(ci_attr$lwr.ci, ci_attr$upr.ci))))
  expect_true(is.list(inf_attr))
  expect_true(all(is.finite(c(inf_attr$se, inf_attr$statistic, inf_attr$p_value))))
})

test_that("gwet_ac returns NA for a degenerate custom-weight denominator", {
  x <- factor(c("A", "A", "B", "B"), levels = c("A", "B"))
  y <- factor(c("A", "B", "A", "B"), levels = c("A", "B"))
  W <- matrix(1, 2, 2)
  fit <- gwet_ac(x, y, weights = W, ci = TRUE, p_value = TRUE)

  expect_true(is.na(as.numeric(fit)))
  expect_true(is.na(attr(fit, "inference")$se))
  expect_true(is.na(attr(fit, "ci")$lwr.ci))
  expect_true(is.na(attr(fit, "ci")$upr.ci))
})

test_that("gwet_ac multi-rater ratings matches a manual counts implementation", {
  ratings <- data.frame(
    r1 = c("A", "A", "B", "C", "A", "B"),
    r2 = c("A", "B", "B", "C", "A", "B"),
    r3 = c("A", "A", "B", "B", "A", "C"),
    r4 = c("A", "A", "B", "C", "B", "B"),
    stringsAsFactors = FALSE
  )

  counts <- ratings_to_counts(ratings, levels = c("A", "B", "C"))
  manual <- manual_gwet_counts(counts)
  fit <- gwet_ac(ratings, input = "ratings")

  expect_equal(fit$observed_agreement[[1L]], manual$observed, tolerance = 1e-12)
  expect_equal(fit$expected_agreement[[1L]], manual$expected, tolerance = 1e-12)
  expect_equal(fit$ac[[1L]], manual$ac, tolerance = 1e-12)
})

test_that("gwet_ac weighted multi-rater AC2 matches a manual counts implementation", {
  ratings <- data.frame(
    r1 = c("L", "L", "M", "H", "M", "H"),
    r2 = c("L", "M", "M", "H", "H", "H"),
    r3 = c("L", "L", "H", "H", "M", "M"),
    stringsAsFactors = FALSE
  )
  counts <- ratings_to_counts(ratings, levels = c("L", "M", "H"))
  W <- matrix(
    c(
      1.0, 0.5, 0.0,
      0.5, 1.0, 0.5,
      0.0, 0.5, 1.0
    ),
    nrow = 3,
    byrow = TRUE
  )
  manual <- manual_gwet_counts(counts, weights = W)
  fit <- gwet_ac(ratings, input = "ratings", weights = W)

  expect_identical(attr(fit, "coefficient"), "AC2")
  expect_equal(fit$observed_agreement[[1L]], manual$observed, tolerance = 1e-12)
  expect_equal(fit$expected_agreement[[1L]], manual$expected, tolerance = 1e-12)
  expect_equal(fit$ac[[1L]], manual$ac, tolerance = 1e-12)
})

test_that("gwet_ac ratings and counts panel inputs agree", {
  ratings <- data.frame(
    r1 = c("A", "A", "B", "C", "A", "B"),
    r2 = c("A", "B", "B", "C", "A", "B"),
    r3 = c("A", "A", "B", "B", "A", "C"),
    r4 = c("A", "A", "B", "C", "B", "B"),
    stringsAsFactors = FALSE
  )
  counts <- ratings_to_counts(ratings, levels = c("A", "B", "C"))
  colnames(counts) <- c("A", "B", "C")

  fit_ratings <- gwet_ac(ratings, input = "ratings")
  fit_counts <- gwet_ac(counts, input = "counts")
  expect_equal(fit_ratings$ac[[1L]], fit_counts$ac[[1L]], tolerance = 1e-12)
  expect_equal(fit_ratings$observed_agreement[[1L]], fit_counts$observed_agreement[[1L]], tolerance = 1e-12)
  expect_equal(fit_ratings$expected_agreement[[1L]], fit_counts$expected_agreement[[1L]], tolerance = 1e-12)
})

test_that("gwet_ac by-category panel results match manual binary collapses", {
  counts <- rbind(
    c(3, 0, 0),
    c(2, 1, 0),
    c(0, 3, 0),
    c(0, 1, 2)
  )
  colnames(counts) <- c("A", "B", "C")

  fit <- gwet_ac(counts, input = "counts", by_category = TRUE)
  by_cat <- attr(fit, "by_category")

  expect_true(is.data.frame(by_cat))
  expect_equal(nrow(by_cat), ncol(counts))
  for (j in seq_len(ncol(counts))) {
    manual <- manual_gwet_counts(collapse_counts_for_category_ac1(counts, j))
    expect_equal(by_cat$ac[[j]], manual$ac, tolerance = 1e-12)
  }
})

test_that("gwet_ac panel inference keeps asymptotic primary and jackknife optional", {
  counts <- rbind(
    c(3, 0, 0),
    c(2, 1, 0),
    c(0, 3, 0),
    c(0, 1, 2),
    c(1, 2, 0)
  )
  colnames(counts) <- c("A", "B", "C")

  fit_default <- gwet_ac(counts, input = "counts", ci = TRUE, p_value = TRUE)
  fit_jk <- gwet_ac(counts, input = "counts", ci = TRUE, p_value = TRUE, se_method = "jackknife")

  expect_identical(attr(fit_default, "se.method"), "asymptotic")
  expect_identical(attr(fit_jk, "se.method"), "jackknife")
  expect_equal(fit_default$ac[[1L]], fit_jk$ac[[1L]], tolerance = 1e-12)
  expect_true(is.finite(fit_default$se[[1L]]) || is.na(fit_default$se[[1L]]))
  expect_true(is.finite(fit_jk$se[[1L]]) || is.na(fit_jk$se[[1L]]))
})

test_that("gwet_ac rejects by-category AC2 output", {
  counts <- rbind(
    c(3, 0, 0),
    c(2, 1, 0),
    c(0, 3, 0)
  )
  colnames(counts) <- c("L", "M", "H")

  expect_error(
    gwet_ac(counts, input = "counts", weights = "linear", by_category = TRUE),
    "must be FALSE for weighted AC2 fits"
  )
})

test_that("gwet_ac thread count does not change estimates", {
  raters <- data.frame(
    a = c("A", "A", "B", "B", "C", "A", "B", "C", "A", "C"),
    b = c("A", "B", "B", "B", "C", "A", "B", "C", "A", "B"),
    c = c("A", "A", "B", "C", "C", "A", "B", "B", "A", "C"),
    d = c("A", "B", "A", "B", "C", "C", "B", "C", "A", "C"),
    stringsAsFactors = FALSE
  )

  fit1 <- gwet_ac(raters, n_threads = 1L)
  fit2 <- gwet_ac(raters, n_threads = 2L)
  expect_equal(as.matrix(fit1), as.matrix(fit2), tolerance = 1e-12)

  counts <- ratings_to_counts(raters, levels = c("A", "B", "C"))
  fit_panel_1 <- gwet_ac(counts, input = "counts", ci = TRUE, p_value = TRUE, n_threads = 1L)
  fit_panel_2 <- gwet_ac(counts, input = "counts", ci = TRUE, p_value = TRUE, n_threads = 2L)
  expect_equal(fit_panel_1$ac[[1L]], fit_panel_2$ac[[1L]], tolerance = 1e-12)
  expect_equal(fit_panel_1$se[[1L]], fit_panel_2$se[[1L]], tolerance = 1e-12)
})

test_that("gwet_ac print summary and plot methods work across result types", {
  skip_if_not_installed("ggplot2")

  x <- factor(c("A", "A", "B", "B", "A", "B"))
  y <- factor(c("A", "B", "B", "B", "A", "A"))
  fit_scalar <- gwet_ac(x, y, ci = TRUE, p_value = TRUE)

  expect_true(length(capture.output(print(fit_scalar))) > 0L)
  sm_scalar <- summary(fit_scalar)
  expect_s3_class(sm_scalar, "summary.gwet_ac")
  expect_s3_class(plot(fit_scalar, show_value = FALSE), "ggplot")

  raters <- data.frame(
    a = c("A", "A", "B", "B", "C", "A", "B", "C"),
    b = c("A", "B", "B", "B", "C", "A", "B", "C"),
    c = c("A", "A", "B", "C", "C", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  fit_matrix <- gwet_ac(raters)
  expect_true(length(capture.output(print(fit_matrix))) > 0L)
  expect_s3_class(summary(fit_matrix), "summary.corr_result")
  expect_s3_class(plot(fit_matrix, show_value = FALSE), "ggplot")
  expect_error(
    plot(fit_matrix, type = "estimate"),
    "not used for matrix-style correlation heatmaps"
  )

  counts <- ratings_to_counts(raters, levels = c("A", "B", "C"))
  fit_panel <- gwet_ac(counts, input = "counts", by_category = TRUE, ci = TRUE)
  expect_true(length(capture.output(print(fit_panel, show_by_category = TRUE))) > 0L)
  sm_panel <- summary(fit_panel)
  expect_s3_class(sm_panel, "summary.gwet_ac")
  expect_s3_class(plot(fit_panel, type = "estimate"), "ggplot")
  expect_s3_class(plot(fit_panel, type = "by_category"), "ggplot")
})

test_that("gwet_ac matches fixed benchmark values for estimates and CIs", {
  benchmark_tol <- 1e-6

  pair1_x <- c(rep("A", 12), rep("A", 3), rep("A", 1), rep("B", 4), rep("B", 10), rep("B", 2), rep("C", 2), rep("C", 5), rep("C", 11))
  pair1_y <- c(rep("A", 12), rep("B", 3), rep("C", 1), rep("A", 4), rep("B", 10), rep("C", 2), rep("A", 2), rep("B", 5), rep("C", 11))
  fit_pair1 <- gwet_ac(pair1_x, pair1_y, ci = TRUE)
  expect_equal(as.numeric(fit_pair1), 0.49010198, tolerance = benchmark_tol)
  expect_equal(attr(fit_pair1, "ci")$lwr.ci, 0.2882755, tolerance = benchmark_tol)
  expect_equal(attr(fit_pair1, "ci")$upr.ci, 0.6919285, tolerance = benchmark_tol)

  pair2_x <- c(rep("A", 90), rep("B", 10))
  pair2_y <- c(rep("A", 85), rep("B", 5), rep("A", 5), rep("B", 5))
  fit_pair2 <- gwet_ac(pair2_x, pair2_y, ci = TRUE)
  expect_equal(as.numeric(fit_pair2), 0.87804878, tolerance = benchmark_tol)
  expect_equal(attr(fit_pair2, "ci")$lwr.ci, 0.7984957, tolerance = benchmark_tol)
  expect_equal(attr(fit_pair2, "ci")$upr.ci, 0.9576018, tolerance = benchmark_tol)

  counts_u1 <- rbind(
    c(5, 0, 0),
    c(2, 1, 0),
    c(0, 3, 1),
    c(0, 1, 2),
    c(1, 1, 0),
    c(0, 0, 4)
  )
  colnames(counts_u1) <- c("A", "B", "C")
  fit_u1 <- gwet_ac(counts_u1, input = "counts", ci = TRUE)
  expect_equal(fit_u1$ac[[1L]], 0.29228101, tolerance = benchmark_tol)
  expect_equal(fit_u1$lwr.ci[[1L]], -0.3393303, tolerance = benchmark_tol)
  expect_equal(fit_u1$upr.ci[[1L]], 0.9238923, tolerance = benchmark_tol)

  counts_u2 <- rbind(
    c(4, 1, 0, 0),
    c(1, 2, 0, 0),
    c(0, 2, 2, 0),
    c(0, 1, 1, 1),
    c(0, 0, 2, 3),
    c(1, 0, 0, 1)
  )
  colnames(counts_u2) <- c("W", "X", "Y", "Z")
  fit_u2 <- gwet_ac(counts_u2, input = "counts", ci = TRUE)
  expect_equal(fit_u2$ac[[1L]], 0.03861956, tolerance = benchmark_tol)
  expect_equal(fit_u2$lwr.ci[[1L]], -0.2970932, tolerance = benchmark_tol)
  expect_equal(fit_u2$upr.ci[[1L]], 0.3743323, tolerance = benchmark_tol)
})

test_that("gwet_ac AC2 matches fixed benchmark values for estimates and CIs", {
  benchmark_tol <- 1e-6

  levels3 <- c("L", "M", "H")
  W <- matrix(
    c(
      1.0, 0.5, 0.0,
      0.5, 1.0, 0.5,
      0.0, 0.5, 1.0
    ),
    nrow = 3,
    byrow = TRUE
  )

  pair_x <- factor(c("L", "L", "M", "M", "H", "H"), levels = levels3)
  pair_y <- factor(c("L", "M", "M", "H", "H", "M"), levels = levels3)
  fit_pair <- gwet_ac(pair_x, pair_y, weights = W, ci = TRUE)
  expect_equal(as.numeric(fit_pair), 0.4517767, tolerance = benchmark_tol)
  expect_equal(attr(fit_pair, "ci")$lwr.ci, -0.1044028, tolerance = benchmark_tol)
  expect_equal(attr(fit_pair, "ci")$upr.ci, 1.0000000, tolerance = benchmark_tol)

  counts <- rbind(
    c(3, 0, 0),
    c(2, 1, 0),
    c(0, 2, 1),
    c(0, 0, 3),
    c(0, 2, 1),
    c(0, 1, 2)
  )
  colnames(counts) <- levels3
  fit_panel <- gwet_ac(counts, input = "counts", weights = W, ci = TRUE)
  expect_equal(fit_panel$ac[[1L]], 0.5057208, tolerance = benchmark_tol)
  expect_equal(fit_panel$lwr.ci[[1L]], 0.1044950, tolerance = benchmark_tol)
  expect_equal(fit_panel$upr.ci[[1L]], 0.9069466, tolerance = benchmark_tol)
})
