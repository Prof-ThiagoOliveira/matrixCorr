test_that("na_method complete matches manual listwise deletion", {
  set.seed(1)
  X <- matrix(rnorm(25 * 5), 25, 5)
  X[sample(length(X), 8)] <- NA_real_
  colnames(X) <- paste0("V", 1:5)
  Xcc <- X[stats::complete.cases(X) & apply(is.finite(X), 1L, all), , drop = FALSE]

  funs <- list(
    pearson_corr = function(z, na_method) pearson_corr(z, na_method = na_method, ci = FALSE),
    spearman_rho = function(z, na_method) spearman_rho(z, na_method = na_method, ci = FALSE),
    kendall_tau = function(z, na_method) kendall_tau(z, na_method = na_method, ci = FALSE),
    dcor = function(z, na_method) dcor(z, na_method = na_method, p_value = FALSE),
    robust_dcor = function(z, na_method) robust_dcor(z, na_method = na_method, p_value = FALSE),
    wincor = function(z, na_method) wincor(z, na_method = na_method, ci = FALSE, p_value = FALSE),
    pbcor = function(z, na_method) pbcor(z, na_method = na_method, ci = FALSE, p_value = FALSE),
    bicor = function(z, na_method) bicor(z, na_method = na_method, ci = FALSE),
    skipped_corr = function(z, na_method) skipped_corr(z, na_method = na_method, ci = FALSE, p_value = FALSE)
  )

  for (nm in names(funs)) {
    R_complete <- funs[[nm]](X, "complete")
    R_manual <- funs[[nm]](Xcc, "error")
    mat_complete <- unclass(as.matrix(R_complete))
    mat_manual <- unclass(as.matrix(R_manual))
    attributes(mat_complete) <- attributes(mat_complete)[c("dim", "dimnames")]
    attributes(mat_manual) <- attributes(mat_manual)[c("dim", "dimnames")]
    expect_equal(
      mat_complete,
      mat_manual,
      tolerance = 1e-8,
      info = nm
    )
  }
})

test_that("complete diagnostics record the common sample", {
  set.seed(1)
  X <- matrix(rnorm(25 * 5), 25, 5)
  X[sample(length(X), 8)] <- NA_real_
  colnames(X) <- paste0("V", 1:5)
  Xcc <- X[stats::complete.cases(X) & apply(is.finite(X), 1L, all), , drop = FALSE]

  R <- pearson_corr(X, na_method = "complete")
  diag <- attr(R, "diagnostics")

  expect_equal(diag$na_method, "complete")
  expect_equal(diag$n_original, nrow(X))
  expect_equal(diag$n_complete, nrow(Xcc))
  expect_equal(diag$n_removed, nrow(X) - nrow(Xcc))
  expect_true(diag$common_sample)
})

test_that("error and pairwise behaviours remain available", {
  set.seed(2)
  X <- matrix(rnorm(50), 10, 5)
  colnames(X) <- paste0("V", 1:5)
  X_na <- X
  X_na[1, 1] <- NA_real_

  expect_equal(
    pearson_corr(X, na_method = "error"),
    pearson_corr(X, na_method = "error")
  )
  R_pairwise <- pearson_corr(X_na, na_method = "pairwise")
  expect_true(!is.null(attr(R_pairwise, "diagnostics")$n_complete))
})

test_that("complete errors when too few rows remain", {
  X_bad <- matrix(NA_real_, 5, 3)
  X_bad[1, ] <- c(1, 2, 3)

  expect_error(
    pearson_corr(X_bad, na_method = "complete"),
    "complete rows"
  )
})

test_that("pcorr complete matches manual listwise deletion", {
  set.seed(3)
  X <- matrix(rnorm(25 * 5), 25, 5)
  X[sample(length(X), 8)] <- NA_real_
  colnames(X) <- paste0("V", 1:5)
  Xcc <- X[stats::complete.cases(X) & apply(is.finite(X), 1L, all), , drop = FALSE]

  R_complete <- pcorr(X, na_method = "complete")
  R_manual <- pcorr(Xcc, na_method = "error")

  expect_equal(
    unclass(R_complete$pcor),
    unclass(R_manual$pcor),
    tolerance = 1e-8
  )
})
