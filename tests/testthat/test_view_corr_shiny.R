test_that("viewer helpers coerce package outputs", {
  skip_if_not_installed("stats")

  pearson <- pearson_corr(mtcars)
  spearman <- spearman_rho(mtcars)
  distance <- dcor(mtcars)
  partial <- pcorr(mtcars, method = "ridge", lambda = 1e-2)
  kendall <- kendall_tau(mtcars)
  biweight <- bicor(mtcars)
  schafer <- schafer_corr(mtcars)

  res <- matrixCorr:::`.mc_prepare_corr_inputs`(list(
    Pearson = pearson,
    Spearman = spearman,
    Distance = distance,
    Partial = partial,
    Kendall = kendall,
    Biweight = biweight,
    Schafer = schafer
  ))

  expect_true(all(vapply(res, function(x) is.matrix(x$matrix), logical(1))))
  expect_equal(res$Distance$signed, FALSE)
  expect_true(res$Pearson$signed)
  expect_equal(colnames(res$Partial$matrix), colnames(partial$pcor))
  expect_equal(res$Kendall$class, paste(class(kendall), collapse = ", "))
  expect_equal(res$Biweight$class, paste(class(biweight), collapse = ", "))
  expect_equal(res$Schafer$class, paste(class(schafer), collapse = ", "))
})

test_that("heatmap helper returns ggplot when plotly missing", {
  mat <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
  p <- matrixCorr:::`.mc_build_heatmap`(
    mat = mat,
    signed = TRUE,
    show_values = TRUE,
    use_abs = FALSE,
    use_plotly = FALSE
  )
  expect_s3_class(p, "ggplot")
})

test_that("distribution helpers use unique off-diagonal pairs", {
  mat <- matrix(
    c(1, 0.4, -0.2,
      0.4, 1, 0.1,
      -0.2, 0.1, 1),
    nrow = 3,
    byrow = TRUE
  )

  vals_signed <- matrixCorr:::`.mc_corr_distribution_values`(mat, use_abs = FALSE)
  vals_abs <- matrixCorr:::`.mc_corr_distribution_values`(mat, use_abs = TRUE)

  expect_equal(vals_signed, c(0.4, -0.2, 0.1))
  expect_equal(vals_abs, c(0.4, 0.2, 0.1))
})

test_that("distribution summary table reports visible pair statistics", {
  tbl <- matrixCorr:::`.mc_corr_summary_table`(c(0.4, -0.2, 0.1))

  expect_equal(tbl$Metric[[1]], "Pairs shown")
  expect_equal(tbl$Value[[1]], "3")
  expect_true("Median" %in% tbl$Metric)
  expect_true("Mean" %in% tbl$Metric)
})

test_that("distribution axis limits adapt to visible correlation range", {
  limits_signed <- matrixCorr:::`.mc_corr_axis_limits`(c(0.1, 0.3, 0.6), use_abs = FALSE)
  limits_abs <- matrixCorr:::`.mc_corr_axis_limits`(c(0.1, 0.3, 0.6), use_abs = TRUE)

  expect_true(limits_signed[[1]] > -1)
  expect_true(limits_signed[[2]] < 1)
  expect_true(limits_signed[[1]] <= 0.1)
  expect_true(limits_signed[[2]] >= 0.6)
  expect_true(limits_abs[[1]] >= 0)
  expect_true(limits_abs[[2]] <= 1)
  expect_true(limits_abs[[1]] <= 0.1)
  expect_true(limits_abs[[2]] >= 0.6)
})

test_that("matrix reorder helper handles absolute and signed modes", {
  mat <- matrix(c(1, 0.2, 0.8,
                  0.2, 1, 0.3,
                  0.8, 0.3, 1), nrow = 3, byrow = TRUE)
  res_abs <- matrixCorr:::`.mc_reorder_matrix`(mat, mode = "abs", method = "complete")
  expect_null(res_abs$message)
  expect_equal(dim(res_abs$matrix), dim(mat))
  expect_true(isSymmetric(res_abs$matrix, tol = 1e-12))
  expect_equal(diag(res_abs$matrix), diag(mat)[res_abs$order])
  expect_equal(sort(res_abs$order), seq_len(nrow(mat)))

  res_signed <- matrixCorr:::`.mc_reorder_matrix`(mat, mode = "signed", method = "average")
  expect_null(res_signed$message)
  expect_equal(dim(res_signed$matrix), dim(mat))
  expect_true(isSymmetric(res_signed$matrix, tol = 1e-12))
  expect_equal(diag(res_signed$matrix), diag(mat)[res_signed$order])
  expect_equal(sort(res_signed$order), seq_len(nrow(mat)))
})

test_that("correlation range filter subsets variables in signed mode", {
  mat <- matrix(
    c(1, 0.8, 0.1, -0.4,
      0.8, 1, 0.3, -0.2,
      0.1, 0.3, 1, 0.05,
      -0.4, -0.2, 0.05, 1),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))
  )

  res <- matrixCorr:::`.mc_filter_matrix_by_range`(
    mat = mat,
    corr_range = c(0.25, 0.85),
    use_abs = FALSE
  )

  expect_equal(colnames(res$matrix), c("A", "B", "C"))
  expect_match(res$message, "kept 3 of 4 variables")
})

test_that("correlation range filter subsets variables in absolute mode", {
  mat <- matrix(
    c(1, 0.8, 0.1, -0.4,
      0.8, 1, 0.3, -0.2,
      0.1, 0.3, 1, 0.05,
      -0.4, -0.2, 0.05, 1),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C", "D"), c("A", "B", "C", "D"))
  )

  res <- matrixCorr:::`.mc_filter_matrix_by_range`(
    mat = mat,
    corr_range = c(0.35, 0.5),
    use_abs = TRUE
  )

  expect_equal(colnames(res$matrix), c("A", "D"))
  expect_match(res$message, "kept 2 of 4 variables")
})

test_that("correlation range helper clamps out-of-bounds ranges", {
  expect_equal(matrixCorr:::`.mc_clamp_corr_range`(c(-2, 0.4), use_abs = TRUE), c(0, 0.4))
  expect_equal(matrixCorr:::`.mc_clamp_corr_range`(c(0.8, -2), use_abs = FALSE), c(-1, 0.8))
})

test_that("viewer helpers resolve logo asset and sanitize filenames", {
  expect_true(file.exists(matrixCorr:::`.mc_logo_path`()))
  expect_equal(matrixCorr:::`.mc_sanitize_filename`("Repeated measures / matrix"), "Repeated-measures-matrix")
})

test_that("signed plotly palette keeps blue at positive and red at negative", {
  scale <- matrixCorr:::`.mc_plotly_colorscale`(TRUE)

  expect_equal(scale[[1]], list(0, "#FF6A6A"))
  expect_equal(scale[[2]], list(0.5, "white"))
  expect_equal(scale[[3]], list(1, "#63B8FF"))
})

test_that("heatmap click helpers recover selected variable pairs", {
  mat <- matrix(
    c(1, 0.4,
      0.4, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("A", "B"), c("A", "B"))
  )
  df <- matrixCorr:::`.mc_heatmap_df`(mat)

  expect_equal(
    matrixCorr:::`.mc_plotly_click_pair`(data.frame(x = "B", y = "A"), vars = c("A", "B")),
    c("B", "A")
  )

  expect_equal(
    matrixCorr:::`.mc_static_click_pair`(df, list(x = "B", y = "A")),
    c("B", "A")
  )
})
