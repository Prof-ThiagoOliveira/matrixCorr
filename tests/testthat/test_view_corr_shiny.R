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
