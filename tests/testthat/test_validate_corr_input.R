test_that("validate_corr_input rejects non-finite real matrix values", {
  X <- cbind(a = c(1, 2, Inf, 4), b = c(1, 2, 3, 4))

  expect_error(
    validate_corr_input(X, check_na = TRUE),
    "Missing values are not allowed.*non-finite values are not allowed"
  )

  expect_error(
    validate_corr_input(replace(X, 3, -Inf), check_na = TRUE),
    "Missing values are not allowed.*non-finite values are not allowed"
  )
})

test_that("validate_corr_input rejects non-finite real data-frame values", {
  df <- data.frame(a = c(1, 2, NaN, 4), b = c(1, 2, Inf, 4))

  expect_error(
    validate_corr_input(df, check_na = TRUE),
    "Missing values are not allowed.*non-finite values are not allowed"
  )
})

test_that("validate_corr_input permits non-finite values only when requested", {
  X <- cbind(a = c(1, 2, Inf, 4), b = c(1, 2, 3, 4))

  out <- validate_corr_input(X, check_na = FALSE)

  expect_true(is.matrix(out))
  expect_true(is.infinite(out[3, "a"]))
})
