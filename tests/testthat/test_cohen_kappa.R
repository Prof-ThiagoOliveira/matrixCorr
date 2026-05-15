manual_cohen_kappa <- function(x, y) {
  keep <- stats::complete.cases(x, y)
  x <- as.character(x[keep])
  y <- as.character(y[keep])
  if (length(x) < 2L) {
    return(NA_real_)
  }
  lev <- unique(c(x, y))
  tab <- table(
    factor(x, levels = lev),
    factor(y, levels = lev)
  )
  p <- unclass(tab) / sum(tab)
  po <- sum(diag(p))
  pe <- sum(rowSums(p) * colSums(p))
  if (!is.finite(1 - pe) || (1 - pe) <= 0) {
    return(NA_real_)
  }
  (po - pe) / (1 - pe)
}

test_that("cohen_kappa reproduces the manual 2x2 example", {
  x <- c(rep("c1", 25), rep("c2", 25))
  y <- c(
    rep("c1", 20), rep("c2", 5),
    rep("c1", 10), rep("c2", 15)
  )

  fit <- cohen_kappa(x, y)
  expect_s3_class(fit, "cohen_kappa")
  expect_equal(as.numeric(fit), 0.4, tolerance = 1e-12)
  expect_equal(attr(fit, "diagnostics")$observed_agreement, 0.7, tolerance = 1e-12)
  expect_equal(attr(fit, "diagnostics")$expected_agreement, 0.5, tolerance = 1e-12)
})

test_that("cohen_kappa returns the verified 2-category estimate and delta CI", {
  x <- factor(c(rep("c1", 25), rep("c2", 25)))
  y <- factor(c(
    rep("c1", 20), rep("c2", 5),
    rep("c1", 10), rep("c2", 15)
  ))

  fit <- cohen_kappa(x, y, ci = TRUE, p_value = TRUE)
  ci_attr <- attr(fit, "ci")
  inf_attr <- attr(fit, "inference")

  expect_equal(as.numeric(fit), 0.4, tolerance = 1e-12)
  expect_equal(ci_attr$lwr.ci, 0.1510923, tolerance = 1e-6)
  expect_equal(ci_attr$upr.ci, 0.6489077, tolerance = 1e-6)
  expect_equal(inf_attr$se, 0.1269961, tolerance = 1e-6)
})

test_that("cohen_kappa returns the verified 3-category estimate and delta CI", {
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

  fit <- cohen_kappa(
    factor(x, levels = rownames(tab)),
    factor(y, levels = colnames(tab)),
    ci = TRUE,
    p_value = TRUE
  )
  ci_attr <- attr(fit, "ci")
  inf_attr <- attr(fit, "inference")

  expect_equal(as.numeric(fit), 0.4916268, tolerance = 1e-6)
  expect_equal(ci_attr$lwr.ci, 0.2964395, tolerance = 1e-6)
  expect_equal(ci_attr$upr.ci, 0.6868141, tolerance = 1e-6)
  expect_equal(inf_attr$se, 0.0995872, tolerance = 1e-6)
})

test_that("cohen_kappa matrix mode is symmetric and matches manual pairwise values", {
  raters <- data.frame(
    r1 = factor(c("A", "A", "B", "B", "C", "C", "A", "B")),
    r2 = c("A", "B", "B", "B", "C", "A", "A", "B"),
    r3 = c(1, 1, 2, 2, 3, 3, 1, 2),
    r4 = c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  fit <- cohen_kappa(raters)
  mat <- as.matrix(fit)
  expect_true(isSymmetric(mat))
  expect_equal(unname(diag(mat)), rep(1, ncol(raters)))

  expected <- matrix(NA_real_, ncol(raters), ncol(raters), dimnames = dimnames(mat))
  diag(expected) <- 1
  for (j in seq_len(ncol(raters) - 1L)) {
    for (k in (j + 1L):ncol(raters)) {
      expected[j, k] <- manual_cohen_kappa(raters[[j]], raters[[k]])
      expected[k, j] <- expected[j, k]
    }
  }

  expect_identical(dimnames(mat), dimnames(expected))
  expect_equal(as.vector(mat), as.vector(expected), tolerance = 1e-12)
})

test_that("cohen_kappa missing-data modes behave as documented", {
  raters <- data.frame(
    a = c("A", "A", "B", NA, "C", "A"),
    b = c("A", "B", "B", "B", NA, "A"),
    c = c("A", "A", NA, "B", "C", "A"),
    stringsAsFactors = FALSE
  )

  expect_error(
    cohen_kappa(raters, na_method = "error"),
    "contains missing values"
  )

  fit_complete <- cohen_kappa(raters, na_method = "complete")
  keep <- stats::complete.cases(raters)
  expected_complete <- matrix(sum(keep), ncol(raters), ncol(raters))
  diag(expected_complete) <- rep(sum(keep), ncol(raters))
  dimnames(expected_complete) <- dimnames(as.matrix(fit_complete))
  expect_equal(attr(fit_complete, "diagnostics")$n_complete, expected_complete)
  expect_equal(
    as.vector(as.matrix(fit_complete)),
    as.vector(as.matrix(cohen_kappa(raters[keep, , drop = FALSE], na_method = "error"))),
    tolerance = 1e-12
  )

  fit_pairwise <- cohen_kappa(raters, na_method = "pairwise")
  n_complete <- attr(fit_pairwise, "diagnostics")$n_complete
  expect_equal(n_complete["a", "b"], sum(stats::complete.cases(raters[c("a", "b")])))
  expect_equal(n_complete["a", "c"], sum(stats::complete.cases(raters[c("a", "c")])))
  expect_equal(n_complete["b", "c"], sum(stats::complete.cases(raters[c("b", "c")])))
})

test_that("cohen_kappa scalar mode has dedicated print, summary, and plot methods", {
  skip_if_not_installed("ggplot2")

  x <- factor(c("A", "A", "B", "B", "A", "B"))
  y <- factor(c("A", "B", "B", "B", "A", "A"))
  fit <- cohen_kappa(x, y, ci = TRUE, p_value = TRUE)

  txt <- capture.output(print(fit))
  expect_true(any(grepl("Cohen's kappa agreement summary", txt, fixed = TRUE)))

  sm <- summary(fit)
  expect_s3_class(sm, "summary.cohen_kappa")
  expect_true(all(c("estimate", "n_complete", "observed_agreement", "expected_agreement") %in% names(sm)))
  expect_true(all(c("lwr", "upr", "se", "statistic", "p_value") %in% names(sm)))

  p <- plot(fit, show_value = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("cohen_kappa attaches CI and inference metadata", {
  raters <- data.frame(
    a = c("A", "A", "B", "B", "C", "A", "B", "C", "A", "C"),
    b = c("A", "B", "B", "B", "C", "A", "B", "C", "A", "B"),
    c = c("A", "A", "B", "C", "C", "A", "B", "B", "A", "C"),
    stringsAsFactors = FALSE
  )

  fit_ci <- cohen_kappa(raters, ci = TRUE)
  ci_attr <- attr(fit_ci, "ci")
  expect_true(is.list(ci_attr))
  expect_true(all(c("est", "lwr.ci", "upr.ci", "conf.level", "ci.method") %in% names(ci_attr)))
  expect_identical(ci_attr$ci.method, "delta")

  fit_p <- cohen_kappa(raters, p_value = TRUE)
  inf_attr <- attr(fit_p, "inference")
  expect_true(is.list(inf_attr))
  expect_true(all(c("statistic", "p_value") %in% names(inf_attr)))
})

test_that("cohen_kappa handles degenerate and low-information cases", {
  raters <- data.frame(
    a = rep("A", 6),
    b = rep("A", 6),
    c = c("A", "B", "A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )
  fit <- cohen_kappa(raters)
  mat <- as.matrix(fit)
  expect_true(is.na(mat["a", "b"]))
  expect_equal(unname(diag(mat)), c(1, 1, 1))

  low_info <- data.frame(
    a = c("A", NA),
    b = c("A", "B"),
    c = c("A", "A"),
    stringsAsFactors = FALSE
  )
  fit_low <- cohen_kappa(low_info, na_method = "pairwise")
  mat_low <- as.matrix(fit_low)
  expect_true(is.na(mat_low["a", "b"]))
  expect_true(is.na(mat_low["a", "a"]))
  expect_equal(mat_low["b", "b"], 1)
  expect_equal(mat_low["c", "c"], 1)
})
