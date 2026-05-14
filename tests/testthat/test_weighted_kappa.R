manual_weighted_kappa_weights <- function(type, K) {
  idx <- outer(seq_len(K), seq_len(K), "-")
  switch(
    type,
    unweighted = diag(1, K),
    linear = if (K <= 1L) matrix(1, K, K) else 1 - abs(idx) / (K - 1),
    quadratic = if (K <= 1L) matrix(1, K, K) else 1 - (idx^2) / ((K - 1)^2)
  )
}

manual_weighted_kappa <- function(x, y, levels, weight_type = "quadratic") {
  keep <- stats::complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 2L) {
    return(NA_real_)
  }

  lev_chr <- as.character(levels)
  tab <- table(
    factor(as.character(x), levels = lev_chr),
    factor(as.character(y), levels = lev_chr)
  )
  P <- unclass(tab) / sum(tab)
  W <- manual_weighted_kappa_weights(weight_type, length(levels))
  Ao <- sum(W * P)
  Ae <- sum(W * outer(rowSums(P), colSums(P)))
  if (!is.finite(1 - Ae) || (1 - Ae) <= 0) {
    return(NA_real_)
  }
  (Ao - Ae) / (1 - Ae)
}

expected_weighted_edge_df <- function(mat, threshold = 0, diag = TRUE) {
  idx <- upper.tri(mat, diag = isTRUE(diag))
  vals <- mat[idx]
  keep <- is.finite(vals) & abs(vals) >= threshold
  data.frame(
    row = rownames(mat)[row(mat)[idx][keep]],
    col = colnames(mat)[col(mat)[idx][keep]],
    value = as.numeric(vals[keep]),
    stringsAsFactors = FALSE
  )
}

expected_weighted_sparse_dense <- function(mat, threshold = 0, diag = TRUE) {
  out <- matrix(0, nrow(mat), ncol(mat), dimnames = dimnames(mat))
  idx <- upper.tri(mat, diag = isTRUE(diag))
  vals <- mat[idx]
  keep <- is.finite(vals) & abs(vals) >= threshold
  out[idx] <- ifelse(keep, vals, 0)
  out[lower.tri(out)] <- t(out)[lower.tri(out)]
  out
}

test_that("weighted_kappa reproduces a manual two-vector example", {
  lev <- c("low", "mid", "high")
  x <- ordered(c("low", "low", "mid", "mid", "high", "high", "high"), levels = lev)
  y <- ordered(c("low", "mid", "mid", "high", "high", "mid", "high"), levels = lev)

  expected <- manual_weighted_kappa(x, y, levels = lev, weight_type = "linear")
  fit <- weighted_kappa(x, y, weights = "linear")

  expect_type(fit, "double")
  expect_equal(fit, expected, tolerance = 1e-12)
})

test_that("weighted_kappa two-rater regression example matches known estimate and SE", {
  lev <- c("low", "mid", "high")
  x <- ordered(
    c("low", "low", "mid", "mid", "high", "high", "high", "mid", "low", "high", "mid", "low"),
    levels = lev
  )
  y <- ordered(
    c("low", "mid", "mid", "high", "high", "mid", "high", "mid", "low", "high", "low", "mid"),
    levels = lev
  )

  fit <- weighted_kappa(x, y, weights = "quadratic")
  fit_inf <- weighted_kappa(data.frame(a = x, b = y), weights = "quadratic", p_value = TRUE)
  inf <- attr(fit_inf, "inference")

  expect_equal(fit, 0.6666666666666667, tolerance = 1e-12)
  expect_equal(inf$se["a", "b"], 0.1385799032138497, tolerance = 1e-12)
  expect_equal(inf$statistic["a", "b"], 4.810702354423639, tolerance = 1e-12)
})

test_that("weighted_kappa resolves linear and quadratic weights and aliases", {
  W_linear <- .mc_resolve_weighted_kappa_weights("linear", 4L)
  attr(W_linear, "weight_type") <- NULL
  W_quadratic <- .mc_resolve_weighted_kappa_weights("quadratic", 4L)
  attr(W_quadratic, "weight_type") <- NULL

  expect_equal(
    W_linear,
    manual_weighted_kappa_weights("linear", 4L),
    tolerance = 1e-12
  )
  expect_equal(
    W_quadratic,
    manual_weighted_kappa_weights("quadratic", 4L),
    tolerance = 1e-12
  )

  expect_equal(
    .mc_resolve_weighted_kappa_weights("equal", 4L),
    .mc_resolve_weighted_kappa_weights("linear", 4L),
    tolerance = 1e-12
  )
  expect_equal(
    .mc_resolve_weighted_kappa_weights("equal_spacing", 4L),
    .mc_resolve_weighted_kappa_weights("linear", 4L),
    tolerance = 1e-12
  )
  expect_equal(
    .mc_resolve_weighted_kappa_weights("equal-spacing", 4L),
    .mc_resolve_weighted_kappa_weights("linear", 4L),
    tolerance = 1e-12
  )

  expect_equal(
    .mc_resolve_weighted_kappa_weights("squared", 4L),
    .mc_resolve_weighted_kappa_weights("quadratic", 4L),
    tolerance = 1e-12
  )
  expect_equal(
    .mc_resolve_weighted_kappa_weights("fleiss_cohen", 4L),
    .mc_resolve_weighted_kappa_weights("quadratic", 4L),
    tolerance = 1e-12
  )
  expect_equal(
    .mc_resolve_weighted_kappa_weights("fleiss-cohen", 4L),
    .mc_resolve_weighted_kappa_weights("quadratic", 4L),
    tolerance = 1e-12
  )
})

test_that("unweighted weighted_kappa matches cohen_kappa", {
  x <- factor(c("A", "A", "B", "B", "C", "A", "B"))
  y <- factor(c("A", "B", "B", "B", "C", "A", "C"))

  expect_equal(
    weighted_kappa(x, y, weights = "unweighted", levels = c("A", "B", "C")),
    as.numeric(cohen_kappa(x, y)),
    tolerance = 1e-12
  )
})

test_that("weighted_kappa matrix mode is symmetric and matches manual pairwise values", {
  lev <- c("low", "mid", "high")
  raters <- data.frame(
    r1 = ordered(c("low", "low", "mid", "mid", "high", "high", "low", "high"), levels = lev),
    r2 = ordered(c("low", "mid", "mid", "high", "high", "mid", "low", "high"), levels = lev),
    r3 = ordered(c("low", "low", "mid", "high", "high", "high", "low", "mid"), levels = lev),
    r4 = ordered(c("mid", "mid", "mid", "high", "high", "high", "low", "mid"), levels = lev)
  )

  fit <- weighted_kappa(raters)
  mat <- as.matrix(fit)
  expect_true(isSymmetric(mat))
  expect_equal(unname(diag(mat)), rep(1, ncol(raters)))

  expected <- matrix(NA_real_, ncol(raters), ncol(raters), dimnames = dimnames(mat))
  diag(expected) <- 1
  for (j in seq_len(ncol(raters) - 1L)) {
    for (k in (j + 1L):ncol(raters)) {
      expected[j, k] <- manual_weighted_kappa(raters[[j]], raters[[k]], levels = lev)
      expected[k, j] <- expected[j, k]
    }
  }

  expect_equal(as.vector(mat), as.vector(expected), tolerance = 1e-12)
})

test_that("weighted_kappa category order handling follows the documented rules", {
  lev <- c("low", "mid", "high")
  ordered_df <- data.frame(
    a = ordered(c("low", "mid", "high", "mid"), levels = lev),
    b = ordered(c("low", "mid", "mid", "high"), levels = lev)
  )
  ordered_fit <- weighted_kappa(ordered_df)
  expect_equal(attr(ordered_fit, "diagnostics")$levels, lev)

  x_num <- c(2, 1, 3, 2, 1, 3)
  y_num <- c(2, 2, 3, 1, 1, 3)
  numeric_fit <- weighted_kappa(x_num, y_num, weights = "linear")
  expect_equal(
    numeric_fit,
    manual_weighted_kappa(x_num, y_num, levels = sort(unique(c(x_num, y_num))), weight_type = "linear"),
    tolerance = 1e-12
  )

  bad_nominal <- data.frame(
    a = c("low", "mid", "high"),
    b = c("low", "high", "mid"),
    stringsAsFactors = FALSE
  )
  expect_error(weighted_kappa(bad_nominal), "must be supplied")

  fit_levels <- weighted_kappa(bad_nominal, levels = lev)
  diag_attr <- attr(fit_levels, "diagnostics")
  expect_equal(diag_attr$levels, lev)
  expect_equal(dim(diag_attr$weights), c(3L, 3L))

  fit_unobs <- weighted_kappa(bad_nominal, levels = c("low", "mid", "high", "very_high"))
  expect_equal(dim(attr(fit_unobs, "diagnostics")$weights), c(4L, 4L))
})

test_that("weighted_kappa missing-data modes behave as documented", {
  lev <- c("low", "mid", "high")
  raters <- data.frame(
    a = ordered(c("low", "low", "mid", NA, "high", "low"), levels = lev),
    b = ordered(c("low", "mid", "mid", "high", NA, "low"), levels = lev),
    c = ordered(c("low", "low", NA, "high", "high", "low"), levels = lev)
  )

  expect_error(
    weighted_kappa(raters, na_method = "error"),
    "contains missing values"
  )

  fit_complete <- weighted_kappa(raters, na_method = "complete")
  keep <- stats::complete.cases(raters)
  expected_complete <- matrix(sum(keep), ncol(raters), ncol(raters))
  diag(expected_complete) <- rep(sum(keep), ncol(raters))
  dimnames(expected_complete) <- dimnames(as.matrix(fit_complete))
  expect_equal(attr(fit_complete, "diagnostics")$n_complete, expected_complete)

  fit_pairwise <- weighted_kappa(raters, na_method = "pairwise")
  n_complete <- attr(fit_pairwise, "diagnostics")$n_complete
  expect_equal(n_complete["a", "b"], sum(stats::complete.cases(raters[c("a", "b")])))
  expect_equal(n_complete["a", "c"], sum(stats::complete.cases(raters[c("a", "c")])))
  expect_equal(n_complete["b", "c"], sum(stats::complete.cases(raters[c("b", "c")])))
})

test_that("weighted_kappa output modes route through corr_result helpers", {
  lev <- c("low", "mid", "high")
  raters <- data.frame(
    a = ordered(c("low", "low", "mid", "mid", "high", "high", "low", "high"), levels = lev),
    b = ordered(c("low", "mid", "mid", "high", "high", "mid", "low", "high"), levels = lev),
    c = ordered(c("low", "low", "mid", "high", "high", "high", "low", "mid"), levels = lev),
    d = ordered(c("mid", "mid", "mid", "high", "high", "high", "low", "mid"), levels = lev)
  )

  fit_matrix <- weighted_kappa(raters, output = "matrix")
  expect_s3_class(fit_matrix, c("corr_matrix", "weighted_kappa", "corr_result", "matrix"))

  edge <- weighted_kappa(raters, output = "edge_list", threshold = 0.2, diag = FALSE)
  edge_df <- as.data.frame(edge, stringsAsFactors = FALSE)
  expect_true(all(c("row", "col", "value") %in% names(edge_df)))

  sparse <- weighted_kappa(raters, output = "sparse", threshold = 0.2, diag = FALSE)
  expect_s4_class(sparse, "sparseMatrix")

  expect_error(
    weighted_kappa(raters, output = "matrix", threshold = 0.2),
    "must be 0 when"
  )
})

test_that("weighted_kappa direct threshold outputs match dense filtering", {
  lev <- c("low", "mid", "high")
  raters <- data.frame(
    a = ordered(c("low", "low", "mid", "mid", "high", "high", "low", "high"), levels = lev),
    b = ordered(c("low", "mid", "mid", "high", "high", "mid", "low", "high"), levels = lev),
    c = ordered(c("low", "low", "mid", "high", "high", "high", "low", "mid"), levels = lev),
    d = ordered(c("mid", "mid", "mid", "high", "high", "high", "low", "mid"), levels = lev)
  )

  base <- weighted_kappa(raters, na_method = "error", ci = FALSE, p_value = FALSE)
  mat <- as.matrix(base)

  edge <- weighted_kappa(
    raters,
    na_method = "error",
    ci = FALSE,
    p_value = FALSE,
    output = "edge_list",
    threshold = 0.2,
    diag = FALSE
  )
  edge_df <- .mc_corr_as_edge_df(edge)
  edge_df <- edge_df[order(edge_df$col, edge_df$row), , drop = FALSE]
  expected_edge <- expected_weighted_edge_df(mat, threshold = 0.2, diag = FALSE)
  expected_edge <- expected_edge[order(expected_edge$col, expected_edge$row), , drop = FALSE]
  rownames(edge_df) <- NULL
  rownames(expected_edge) <- NULL
  expect_equal(edge_df, expected_edge, tolerance = 1e-12)

  sparse <- weighted_kappa(
    raters,
    na_method = "error",
    ci = FALSE,
    p_value = FALSE,
    output = "sparse",
    threshold = 0.2,
    diag = FALSE
  )
  expect_equal(
    as.matrix(sparse),
    expected_weighted_sparse_dense(mat, threshold = 0.2, diag = FALSE),
    tolerance = 1e-12
  )
})

test_that("weighted_kappa attaches CI and inference metadata", {
  lev <- c("low", "mid", "high")
  raters <- data.frame(
    a = ordered(c("low", "low", "mid", "mid", "high", "high", "low", "high", "mid", "high"), levels = lev),
    b = ordered(c("low", "mid", "mid", "high", "high", "mid", "low", "high", "mid", "high"), levels = lev),
    c = ordered(c("low", "low", "mid", "high", "high", "high", "low", "mid", "mid", "high"), levels = lev)
  )

  fit_ci <- weighted_kappa(raters, ci = TRUE)
  ci_attr <- attr(fit_ci, "ci")
  expect_true(is.list(ci_attr))
  expect_true(all(c("est", "lwr.ci", "upr.ci", "conf.level", "ci.method") %in% names(ci_attr)))
  expect_identical(ci_attr$ci.method, "delta")

  fit_p <- weighted_kappa(raters, p_value = TRUE)
  inf_attr <- attr(fit_p, "inference")
  expect_true(is.list(inf_attr))
  expect_true(all(c("estimate", "se", "statistic", "p_value", "n_obs") %in% names(inf_attr)))

  degenerate <- data.frame(
    a = ordered(rep("low", 6), levels = lev),
    b = ordered(rep("low", 6), levels = lev),
    c = ordered(c("low", "mid", "low", "mid", "low", "mid"), levels = lev)
  )
  fit_deg <- weighted_kappa(degenerate, ci = TRUE, p_value = TRUE)
  ci_deg <- attr(fit_deg, "ci")
  inf_deg <- attr(fit_deg, "inference")
  expect_true(is.na(ci_deg$lwr.ci["a", "b"]))
  expect_true(is.na(inf_deg$se["a", "b"]))
})

test_that("weighted_kappa handles degenerate and low-information cases", {
  lev <- c("low", "mid", "high")
  raters <- data.frame(
    a = ordered(rep("low", 6), levels = lev),
    b = ordered(rep("low", 6), levels = lev),
    c = ordered(c("low", "mid", "low", "mid", "low", "mid"), levels = lev)
  )
  fit <- weighted_kappa(raters)
  mat <- as.matrix(fit)
  expect_true(is.na(mat["a", "b"]))
  expect_equal(unname(diag(mat)), c(1, 1, 1))

  low_info <- data.frame(
    a = ordered(c("low", NA), levels = lev),
    b = ordered(c("low", "mid"), levels = lev),
    c = ordered(c("low", "low"), levels = lev)
  )
  fit_low <- weighted_kappa(low_info, na_method = "pairwise")
  mat_low <- as.matrix(fit_low)
  expect_true(is.na(mat_low["a", "b"]))
  expect_true(is.na(mat_low["a", "a"]))
  expect_equal(mat_low["b", "b"], 1)
  expect_equal(mat_low["c", "c"], 1)
})

test_that("weighted_kappa uses corr_result S3 summary and plotting", {
  skip_if_not_installed("ggplot2")

  lev <- c("low", "mid", "high")
  raters <- data.frame(
    a = ordered(c("low", "low", "mid", "mid", "high", "high", "low", "high"), levels = lev),
    b = ordered(c("low", "mid", "mid", "high", "high", "mid", "low", "high"), levels = lev),
    c = ordered(c("low", "low", "mid", "high", "high", "high", "low", "mid"), levels = lev)
  )

  fit <- weighted_kappa(raters)
  expect_true(length(capture.output(print(fit))) > 0L)
  expect_s3_class(summary(fit), "summary.corr_result")
  expect_s3_class(plot(fit, show_value = FALSE), "ggplot")
})
