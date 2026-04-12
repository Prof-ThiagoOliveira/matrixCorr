test_that("matrix-style correlation summaries use the standard compact format", {
  set.seed(123)
  X <- matrix(rnorm(120), nrow = 30, ncol = 4)
  colnames(X) <- paste0("V", seq_len(ncol(X)))

  objs <- list(
    pearson_corr(X),
    spearman_rho(X),
    kendall_tau(X),
    dcor(X),
    shrinkage_corr(X),
    bicor(X),
    pbcor(X),
    wincor(X),
    skipped_corr(X)
  )

  for (obj in objs) {
    sm <- summary(obj)
    expect_s3_class(sm, "summary.matrixCorr")
    expect_s3_class(sm, "summary.corr_matrix")
    expect_identical(sm$n_rows, nrow(obj))
    expect_identical(sm$n_cols, ncol(obj))
    expect_true(isTRUE(sm$symmetric))

    txt <- capture.output(matrixCorr:::print.summary.matrixCorr(sm))
    expect_true(any(grepl("^Correlation summary$", txt)))
    expect_true(any(grepl("pairs", txt, fixed = TRUE)))
    expect_true(any(grepl("estimate", txt, fixed = TRUE)))
    expect_true(any(grepl("Strongest pairs by \\|estimate\\|", txt)))
  }
})

test_that("partial correlation summary follows the same matrix-style contract", {
  set.seed(456)
  X <- matrix(rnorm(160), nrow = 40, ncol = 4)
  colnames(X) <- paste0("P", seq_len(ncol(X)))

  pc <- pcorr(X, method = "ridge", lambda = 1e-2)
  sm <- summary(pc)

  expect_s3_class(sm, "summary.partial_corr")
  expect_s3_class(sm, "summary.matrixCorr")
  expect_s3_class(sm, "summary.corr_matrix")
  expect_identical(sm$class, "partial_corr")
  expect_identical(sm$method, "ridge")
  expect_equal(sm$lambda, 1e-2)
  expect_identical(sm$n_rows, 4L)
  expect_identical(sm$n_cols, 4L)

  txt <- capture.output(matrixCorr:::print.summary.partial_corr(sm))
  expect_true(any(grepl("^Correlation summary$", txt)))
  expect_true(any(grepl("lambda", txt, fixed = TRUE)))
})

test_that("latent summaries retain the latent header", {
  set.seed(789)
  X <- data.frame(
    x1 = rnorm(50),
    x2 = rnorm(50)
  )
  Y <- data.frame(
    y1 = ordered(sample(1:3, 50, TRUE)),
    y2 = ordered(sample(1:4, 50, TRUE))
  )

  ps <- polyserial(X, Y)
  sm <- summary(ps)

  expect_s3_class(sm, "summary.latent_corr")
  expect_s3_class(sm, "summary.matrixCorr")
  expect_s3_class(sm, "summary.corr_matrix")
  txt <- capture.output(matrixCorr:::print.summary.latent_corr(sm))
  expect_true(any(grepl("^Latent correlation summary$", txt)))
})
