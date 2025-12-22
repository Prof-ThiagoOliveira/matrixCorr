test_that("viewer helpers coerce package outputs", {
  skip_if_not_installed("stats")

  pearson <- pearson_corr(mtcars)
  spearman <- spearman_rho(mtcars)
  distance <- distance_corr(mtcars)
  partial <- partial_correlation(mtcars, method = "ridge", lambda = 1e-2)
  kendall <- kendall_tau(mtcars)
  biweight <- biweight_mid_corr(mtcars)
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

test_that("matrix reorder helper is stable", {
  mat <- matrix(c(1, 0.2, 0.8,
                  0.2, 1, 0.3,
                  0.8, 0.3, 1), nrow = 3, byrow = TRUE)
  ord <- matrixCorr:::`.mc_reorder_matrix`(mat, signed = TRUE)
  expect_equal(dim(ord), dim(mat))
  expect_true(isSymmetric(ord, tol = 1e-12))
})
