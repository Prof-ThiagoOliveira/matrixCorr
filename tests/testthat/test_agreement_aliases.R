test_that("ba aliases bland_altman", {
  set.seed(1)
  x <- rnorm(40)
  y <- x + rnorm(40, sd = 0.5)

  ref <- bland_altman(x, y)
  ali <- ba(x, y)

  expect_equal(ali$mean.diffs, ref$mean.diffs, tolerance = 1e-12)
  expect_equal(ali$lower.limit, ref$lower.limit, tolerance = 1e-12)
  expect_equal(ali$upper.limit, ref$upper.limit, tolerance = 1e-12)
  expect_equal(ali$based.on, ref$based.on)
})

test_that("ba_rm aliases bland_altman_repeated", {
  set.seed(2)
  S <- 8
  Tm <- 3
  subj <- rep(seq_len(S), each = Tm)
  time <- rep(seq_len(Tm), times = S)
  truth <- rnorm(S, 10, 1)[subj]
  mA <- truth + rnorm(length(truth), sd = 0.3)
  mB <- truth + 0.2 + rnorm(length(truth), sd = 0.4)

  dat <- rbind(
    data.frame(y = mA, subject = subj, method = "A", time = time),
    data.frame(y = mB, subject = subj, method = "B", time = time)
  )

  ref <- bland_altman_repeated(
    data = dat,
    response = "y",
    subject = "subject",
    method = "method",
    time = "time"
  )
  ali <- ba_rm(
    data = dat,
    response = "y",
    subject = "subject",
    method = "method",
    time = "time"
  )

  expect_equal(ali$mean.diffs, ref$mean.diffs, tolerance = 1e-12)
  expect_equal(ali$lower.limit, ref$lower.limit, tolerance = 1e-12)
  expect_equal(ali$upper.limit, ref$upper.limit, tolerance = 1e-12)
  expect_equal(ali$based.on, ref$based.on)
})

test_that("ccc repeated aliases preserve outputs", {
  set.seed(123)
  df <- expand.grid(subject = 1:10, time = 1:2, method = c("A", "B", "C"))
  df$y <- rnorm(nrow(df), mean = match(df$method, c("A", "B", "C")), sd = 1)

  ref_u <- ccc_pairwise_u_stat(
    df,
    response = "y",
    method = "method",
    subject = "subject",
    time = "time"
  )
  ali_u <- ccc_rm_ustat(
    df,
    response = "y",
    method = "method",
    subject = "subject",
    time = "time"
  )

  expect_equal(unclass(ali_u), unclass(ref_u), tolerance = 1e-12)

  ref_reml <- ccc_lmm_reml(
    df,
    response = "y",
    rind = "subject",
    method = "method",
    time = "time"
  )
  ali_reml <- ccc_rm_reml(
    df,
    response = "y",
    rind = "subject",
    method = "method",
    time = "time"
  )

  expect_equal(ali_reml$est, ref_reml$est, tolerance = 1e-12)
})
