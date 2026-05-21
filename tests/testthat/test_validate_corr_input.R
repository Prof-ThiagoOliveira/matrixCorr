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

test_that("integer NA is preserved when check_na = FALSE", {
  X <- matrix(c(1L, NA_integer_, 3L, 4L), ncol = 2)

  out <- validate_corr_input(X, check_na = FALSE)

  expect_true(is.na(out[2, 1]))
  expect_false(is.finite(out[2, 1]))
  expect_false(identical(out[2, 1], -2147483648))
})

test_that("integer NA is preserved in data frames when check_na = FALSE", {
  df <- data.frame(
    a = c(1L, NA_integer_, 3L),
    b = c(4L, 5L, 6L),
    c = letters[1:3]
  )

  out <- validate_corr_input(df, check_na = FALSE)

  expect_true(is.na(out[2, "a"]))
  expect_false(is.finite(out[2, "a"]))
  expect_false(identical(out[2, "a"], -2147483648))
})
