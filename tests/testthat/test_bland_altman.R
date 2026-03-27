test_that("structure, classes and lengths are correct", {
  set.seed(1)
  x <- rnorm(20, 100, 10)
  y <- x + rnorm(20, 0, 8)

  ba <- ba(x, y)

  expect_s3_class(ba, "ba")
  expect_true(all(c("means","diffs","groups","based.on",
                    "lower.limit","mean.diffs","upper.limit",
                    "lines","CI.lines","loa_multiplier","critical.diff") %in% names(ba)))

  # basic size checks
  expect_equal(length(ba$means), ba$based.on)
  expect_equal(length(ba$diffs), ba$based.on)
  expect_equal(nrow(ba$groups),  ba$based.on)

  # named vectors present
  expect_true(all(c("lower.limit","mean.diffs","upper.limit") %in% names(ba$lines)))
  expect_true(all(c("lower.limit.ci.lower","lower.limit.ci.upper",
                    "mean.diff.ci.lower","mean.diff.ci.upper",
                    "upper.limit.ci.lower","upper.limit.ci.upper")
                  %in% names(ba$CI.lines)))
})

test_that("matches known output on the 12-row example (exact numbers)", {
  example.data <- data.frame(
    sex = gl(2, 6, labels = c("f", "m")),
    m1  = c(16,10,14,18,16,15,18,19,14,11,11,17),
    m2  = c(18, 9,15,19,19,13,19,20,14,11,13,17)
  )
  ba <- ba(example.data$m1, example.data$m2)

  # reference values supplied in prompt (tolerance 1e-6)
  expect_equal(as.numeric(ba$based.on), 12L)
  expect_equal(as.numeric(ba$mean.diffs),  -2/3,      tolerance = 1e-6)
  expect_equal(as.numeric(ba$lower.limit), -3.353217, tolerance = 1e-6)
  expect_equal(as.numeric(ba$upper.limit),  2.019883, tolerance = 1e-6)
  expect_equal(as.numeric(ba$critical.diff), 2.68655, tolerance = 1e-5)

  # CI lines
  exp_CI <- c(
    "lower.limit.ci.lower" = -4.8616497,
    "lower.limit.ci.upper" = -1.8447839,
    "mean.diff.ci.lower"   = -1.5375608,
    "mean.diff.ci.upper"   =  0.2042275,
    "upper.limit.ci.lower" =  0.5114506,
    "upper.limit.ci.upper" =  3.5283163
  )
  # compare by name, same tolerance
  for (nm in names(exp_CI)) {
    expect_equal(as.numeric(ba$CI.lines[[nm]]), exp_CI[[nm]], tolerance = 1e-6)
  }
})

test_that("summary.ba returns a compact summary object and prints grouped sections", {
  set.seed(11)
  x <- rnorm(40, 100, 10)
  y <- x + rnorm(40, 0, 8)

  fit <- ba(x, y)
  sm <- summary(fit)

  expect_s3_class(sm, "summary.ba")
  expect_true(all(c("n", "bias", "sd_loa", "loa_low", "loa_up", "width",
                    "loa_multiplier", "bias_lwr", "bias_upr",
                    "lo_lwr", "lo_upr", "up_lwr", "up_upr") %in% names(sm)))

  out <- capture.output(print(sm))
  expect_true(any(grepl("^Agreement estimates$", out)))
  expect_true(any(grepl("^Confidence intervals$", out)))
})

test_that("'loa_multiplier' scales the LoA correctly", {
  set.seed(2)
  x <- rnorm(50, 100, 10)
  y <- x + rnorm(50, 0, 8)

  a <- ba(x, y, loa_multiplier = 1.96)
  b <- ba(x, y, loa_multiplier = 2.00)

  # the LoA width is 2 * loa_multiplier * sd_diffs -> so proportional to 'loa_multiplier'
  width_a <- as.numeric(a$upper.limit - a$lower.limit)
  width_b <- as.numeric(b$upper.limit - b$lower.limit)
  expect_equal(width_b / width_a, 2.00 / 1.96, tolerance = 1e-12)

  # 'critical.diff' equals loa_multiplier * sd_diffs
  sd_a <- as.numeric(a$critical.diff) / 1.96
  sd_b <- as.numeric(b$critical.diff) / 2.00
  expect_equal(sd_a, sd_b, tolerance = 1e-12)
})

test_that("mode=2 flips the sign as expected", {
  set.seed(3)
  x <- rnorm(40, 100, 10)
  y <- x + rnorm(40, 0, 8)

  m1 <- ba(x, y, mode = 1)  # diffs = x - y
  m2 <- ba(x, y, mode = 2)  # diffs = y - x

  # means are the same, diffs flip sign
  expect_equal(as.numeric(m1$means), as.numeric(m2$means))
  expect_equal(as.numeric(m1$diffs), -as.numeric(m2$diffs))

  # scalar summaries flip sign accordingly
  expect_equal(as.numeric(m1$mean.diffs), -as.numeric(m2$mean.diffs))
  expect_equal(as.numeric(m1$lower.limit), -as.numeric(m2$upper.limit), tolerance = 1e-12)
  expect_equal(as.numeric(m1$upper.limit), -as.numeric(m2$lower.limit), tolerance = 1e-12)

  # CI lines for mean difference flip sign; LoA CIs swap and flip
  expect_equal(as.numeric(m1$CI.lines["mean.diff.ci.lower"]),
               -as.numeric(m2$CI.lines["mean.diff.ci.upper"]), tolerance = 1e-12)
  expect_equal(as.numeric(m1$CI.lines["mean.diff.ci.upper"]),
               -as.numeric(m2$CI.lines["mean.diff.ci.lower"]), tolerance = 1e-12)
})

test_that("pairwise NA removal works and 'based.on' reflects that", {
  x <- c(1, 2, NA, 4, 5, 6)
  y <- c(1, NA, 3, 4, NA, 7)
  # valid pairs: (1,1), (4,4), (6,7) => n = 3
  ba <- ba(x, y)
  expect_equal(as.integer(ba$based.on), 3L)
  expect_equal(length(ba$means), 3L)
  expect_equal(nrow(ba$groups), 3L)
})

test_that("constant difference yields zero SD and zero-width LoA/CI", {
  set.seed(4)
  g1 <- rnorm(30, 50, 5)
  g2 <- g1 + 2                       # constant +2
  ba <- ba(g1, g2)         # diffs = g1 - g2 = -2

  expect_equal(as.numeric(ba$mean.diffs), -2, tolerance = 1e-12)
  expect_equal(as.numeric(ba$critical.diff), 0, tolerance = 1e-12)
  expect_equal(as.numeric(ba$lower.limit), -2, tolerance = 1e-12)
  expect_equal(as.numeric(ba$upper.limit), -2, tolerance = 1e-12)

  # all CI lines collapse to the corresponding estimates
  expect_equal(as.numeric(ba$CI.lines["mean.diff.ci.lower"]), -2, tolerance = 1e-12)
  expect_equal(as.numeric(ba$CI.lines["mean.diff.ci.upper"]), -2, tolerance = 1e-12)
  expect_equal(as.numeric(ba$CI.lines["lower.limit.ci.lower"]), -2, tolerance = 1e-12)
  expect_equal(as.numeric(ba$CI.lines["lower.limit.ci.upper"]), -2, tolerance = 1e-12)
  expect_equal(as.numeric(ba$CI.lines["upper.limit.ci.lower"]), -2, tolerance = 1e-12)
  expect_equal(as.numeric(ba$CI.lines["upper.limit.ci.upper"]), -2, tolerance = 1e-12)
})

test_that("errors on invalid inputs", {
  expect_error(ba(1:3, 1:2), "same length")
  expect_error(ba(letters[1:3], 1:3), "numeric", ignore.case = TRUE)
  expect_error(ba(1:3, 1:3, loa_multiplier = 0), "positive")
  expect_error(ba(1:3, 1:3, mode = 0), "1 or 2")

  expect_error(ba(1:3, 1:3, conf_level = 1.5), "in \\(0, 1\\)")
})

test_that("requires at least two complete pairs", {
  x <- c(1, NA)
  y <- c(1, 2)
  expect_error(
    ba(x, y),
    "at least two complete pairs",
    fixed = FALSE
  )
})

test_that("plot.ba returns a ggplot without line-size deprecation warnings", {
  skip_if_not_installed("ggplot2")

  set.seed(5)
  x <- rnorm(40, 100, 10)
  y <- x + rnorm(40, 0, 8)
  fit <- ba(x, y)

  expect_no_warning(
    expect_s3_class(plot(fit), "ggplot")
  )
})
