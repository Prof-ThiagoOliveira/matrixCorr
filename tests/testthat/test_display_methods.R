test_that("matrix previews truncate rows and columns without mutating objects", {
  set.seed(1)
  X <- matrix(rnorm(20 * 12), nrow = 20, ncol = 12)
  colnames(X) <- paste0("V", seq_len(ncol(X)))

  obj <- pearson_corr(X)
  original <- unclass(obj)

  txt <- capture.output(print(obj, n = 8, topn = 2, max_vars = 4, show_ci = "no"))

  expect_true(any(grepl("Pearson correlation matrix", txt, fixed = TRUE)))
  expect_true(any(grepl("... 8 more rows not shown", txt, fixed = TRUE)))
  expect_true(any(grepl("... 8 more variables not shown", txt, fixed = TRUE)))
  expect_true(any(grepl("Use summary() for a richer digest.", txt, fixed = TRUE)))
  expect_equal(unclass(obj), original)
})

test_that("summary output remains bounded and preserves full underlying data", {
  set.seed(2)
  X <- matrix(rnorm(30 * 12), nrow = 30, ncol = 12)
  colnames(X) <- paste0("S", seq_len(ncol(X)))

  obj <- spearman_rho(X)
  sm <- summary(obj, topn = 6)
  txt <- capture.output(print(sm, n = 10, topn = 3, max_vars = 6))

  expect_s3_class(sm, "summary_corr_matrix")
  expect_equal(nrow(sm$top_results), 6L)
  expect_true(any(grepl("^Correlation summary$", txt)))
  expect_true(any(grepl("^Strongest pairs by \\|estimate\\|$", txt)))
  expect_true(any(grepl("Use as.data.frame()/tidy()/as.matrix() to inspect the full result.", txt, fixed = TRUE)))
})

test_that("width-aware truncation responds to width and max_vars overrides", {
  set.seed(3)
  X <- matrix(rnorm(25 * 10), nrow = 25, ncol = 10)
  colnames(X) <- paste0("LongVariable", seq_len(ncol(X)))

  obj <- kendall_tau(X)

  txt_narrow <- capture.output(print(obj, width = 40, max_vars = NULL, show_ci = "no"))
  txt_wide <- capture.output(print(obj, width = 220, max_vars = NULL, show_ci = "no"))
  txt_override <- capture.output(print(obj, width = 220, max_vars = 4, show_ci = "no"))

  expect_true(any(grepl("more variables not shown", txt_narrow, fixed = TRUE)))
  expect_false(any(grepl("more variables not shown", txt_wide, fixed = TRUE)))
  expect_true(any(grepl("... 6 more variables not shown", txt_override, fixed = TRUE)))
})

test_that("binary show_ci controls bounded CI digest display for CI-enabled summaries", {
  set.seed(4)
  X <- matrix(rnorm(40 * 5), nrow = 40, ncol = 5)
  colnames(X) <- paste0("P", seq_len(ncol(X)))

  obj <- pearson_corr(X, ci = TRUE)
  sm <- summary(obj, show_ci = "yes")

  txt_yes <- capture.output(print(sm, show_ci = "yes"))
  txt_no <- capture.output(print(sm, show_ci = "no"))

  expect_true(any(grepl("^\\s*ci\\s*:\\s*95%", txt_yes)))
  expect_true(any(grepl("^\\s*ci_width\\s*:", txt_yes)))
  expect_false(any(grepl("^\\s*ci\\s*:\\s*95%", txt_no)))
  expect_false(any(grepl("^\\s*ci_width\\s*:", txt_no)))
  expect_false(any(grepl("\\blwr\\b", txt_no)))
  expect_false(any(grepl("\\bupr\\b", txt_no)))
})

test_that("display helpers reject auto and invalid display arguments", {
  set.seed(5)
  X <- matrix(rnorm(30 * 6), nrow = 30, ncol = 6)
  colnames(X) <- paste0("T", seq_len(ncol(X)))

  obj <- spearman_rho(X, ci = TRUE)
  sm <- summary(obj, show_ci = "yes")

  expect_error(print(obj, show_ci = "auto"), class = "matrixCorr_arg_error")
  expect_error(summary(obj, show_ci = "auto"), class = "matrixCorr_arg_error")
  expect_error(print(sm, show_ci = "auto"), class = "matrixCorr_arg_error")

  expect_error(print(obj, n = 1), class = "matrixCorr_arg_error")
  expect_error(print(obj, topn = 0), class = "matrixCorr_arg_error")
  expect_error(print(obj, max_vars = 0), class = "matrixCorr_arg_error")
  expect_error(print(obj, width = 0), class = "matrixCorr_arg_error")
  expect_warning(print(obj, n = 6, topn = 4), class = "matrixCorr_display_warning")
  expect_no_warning(print(obj, n = 6))
})

test_that("package display options control defaults", {
  set.seed(6)
  X <- matrix(rnorm(30 * 12), nrow = 30, ncol = 12)
  colnames(X) <- paste0("O", seq_len(ncol(X)))

  obj <- pearson_corr(X, ci = TRUE)
  old <- options(
    matrixCorr.print_max_rows = 8L,
    matrixCorr.print_topn = 2L,
    matrixCorr.print_max_vars = 4L,
    matrixCorr.print_show_ci = "no"
  )
  on.exit(options(old), add = TRUE)

  txt <- capture.output(print(obj))

  expect_true(any(grepl("... 8 more rows not shown", txt, fixed = TRUE)))
  expect_true(any(grepl("... 8 more variables not shown", txt, fixed = TRUE)))
  expect_false(any(grepl("Pearson correlation summary", txt, fixed = TRUE)))
})

test_that("non-correlation S3 classes use the same bounded preview philosophy", {
  set.seed(7)
  x <- rnorm(40)
  y <- x + rnorm(40, sd = 0.5)

  fit <- ba(x, y)

  txt_print <- capture.output(print(fit, show_ci = "no"))
  txt_summary <- capture.output(print(summary(fit), show_ci = "no"))

  expect_true(any(grepl("^Bland-Altman preview:", txt_print)))
  expect_true(any(grepl("^Bland-Altman \\(two methods\\)$", txt_summary)))
  expect_false(any(grepl("bias_lwr", txt_summary, fixed = TRUE)))
  expect_false(any(grepl("bias_upr", txt_summary, fixed = TRUE)))
})

test_that("CI-enabled correlation print methods keep a single primary body", {
  set.seed(8)
  X <- matrix(rnorm(60), nrow = 20, ncol = 3)
  colnames(X) <- c("A", "B", "C")

  p_txt <- capture.output(print(pearson_corr(X, ci = TRUE), show_ci = "yes"))
  s_txt <- capture.output(print(spearman_rho(X, ci = TRUE), show_ci = "yes"))
  k_txt <- capture.output(print(kendall_tau(X, ci = TRUE), show_ci = "yes"))

  expect_false(any(grepl("Pearson correlation summary", p_txt, fixed = TRUE)))
  expect_false(any(grepl("Spearman correlation summary", s_txt, fixed = TRUE)))
  expect_false(any(grepl("Kendall correlation summary", k_txt, fixed = TRUE)))
})

test_that("plot methods honour show_value for numeric overlays", {
  skip_if_not_installed("ggplot2")

  geom_names <- function(p) vapply(p$layers, function(layer) class(layer$geom)[1], character(1))

  set.seed(9)
  X <- matrix(rnorm(80), nrow = 20, ncol = 4)
  colnames(X) <- paste0("V", seq_len(ncol(X)))

  p1 <- plot(pearson_corr(X, ci = TRUE), show_value = FALSE)
  p2 <- plot(pbcor(X), show_value = FALSE)
  p3 <- plot(schafer_corr(X), show_value = FALSE)
  p4 <- plot(schafer_corr(X), show_value = TRUE)

  expect_false(any(geom_names(p1) == "GeomText"))
  expect_false(any(geom_names(p2) == "GeomText"))
  expect_false(any(geom_names(p3) == "GeomText"))
  expect_true(any(geom_names(p4) == "GeomText"))
})
