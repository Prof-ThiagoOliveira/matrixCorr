test_that("check_bool validates logical flags", {
  expect_invisible(check_bool(TRUE))
  expect_invisible(check_bool(FALSE))

  expect_error(check_bool(c(TRUE, FALSE)), class = "matrixCorr_arg_error")
  expect_error(check_bool(1), class = "matrixCorr_arg_error")
})

test_that("check_scalar_numeric enforces bounds and types", {
  expect_invisible(check_scalar_numeric(1))
  expect_invisible(check_scalar_numeric(0, lower = 0, upper = 1))

  expect_error(check_scalar_numeric("a"), class = "matrixCorr_arg_error")
  expect_error(check_scalar_numeric(-1, lower = 0), class = "matrixCorr_arg_error")
  expect_error(check_scalar_numeric(2, upper = 1), class = "matrixCorr_arg_error")
})

test_that("check_required_cols detects missing columns", {
  df <- data.frame(a = 1, b = 2)
  expect_invisible(check_required_cols(df, c("a")))
  expect_error(check_required_cols(df, c("a", "c")), class = "matrixCorr_arg_error")
})

test_that("check_same_length enforces matching lengths", {
  expect_invisible(check_same_length(1:3, 4:6))
  expect_error(check_same_length(1:3, 4:7), class = "matrixCorr_arg_error")
})

test_that("check_prob_scalar covers open/closed intervals", {
  expect_invisible(check_prob_scalar(0))
  expect_invisible(check_prob_scalar(1, open_ends = FALSE))
  expect_error(check_prob_scalar(-0.1), class = "matrixCorr_arg_error")
  expect_error(check_prob_scalar(1, open_ends = TRUE), class = "matrixCorr_arg_error")
})

test_that("check_weights validates length and non-negativity", {
  expect_invisible(check_weights(NULL, n = 2))
  expect_invisible(check_weights(c(0.2, 0.8), n = 2))
  expect_error(check_weights(c(1, -1), n = 2), class = "matrixCorr_arg_error")
  expect_error(check_weights(1:3, n = 2), class = "matrixCorr_arg_error")
})

test_that("match_arg wraps rlang::arg_match", {
  expect_equal(match_arg("a", c("a", "b")), "a")
  expect_error(match_arg("c", c("a", "b")), class = "rlang_error")
})

test_that("inform_if_verbose obeys verbosity flag", {
  expect_silent(inform_if_verbose("quiet", .verbose = FALSE))
  expect_message(inform_if_verbose("loud", .verbose = TRUE), "loud")
})
