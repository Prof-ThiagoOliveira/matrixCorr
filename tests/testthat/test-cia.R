test_that("cia matches a simple two-method hand calculation", {
  dat <- data.frame(
    subject = rep(c("s1", "s2"), each = 4),
    method = rep(c("A", "A", "B", "B"), times = 2),
    replicate = rep(c("r1", "r2"), times = 4),
    value = c(0, 2, 0, 0, 10, 12, 10, 10)
  )

  overall <- cia(
    dat,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    scope = "overall"
  )
  pairwise <- suppressWarnings(cia(
    dat,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    scope = "pairwise",
    estimator = "vc_constrained"
  ))

  expect_s3_class(overall, "cia_overall")
  expect_equal(overall$cia[[1]], 1)
  expect_false(any(c("lwr.ci", "upr.ci") %in% names(overall)))
  expect_s3_class(pairwise, "cia")
  expect_equal(unname(pairwise["A", "B"]), 1)
  expect_equal(unname(pairwise["B", "A"]), 1)
  expect_equal(unname(diag(pairwise)), c(1, 1))
  expect_false("fisher_z" %in% names(summary(pairwise)))
})

test_that("no-reference cia retains raw ratios and reports constrained defaults", {
  dat <- data.frame(
    subject = rep(c("s1", "s2"), each = 6),
    method = rep(c("A", "A", "B", "B", "C", "C"), times = 2),
    replicate = rep(c("r1", "r2"), times = 6),
    value = c(
      0, 2, 0, 0, 1, 1,
      10, 12, 10, 10, 11, 11
    )
  )

  overall <- cia(
    dat,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    scope = "overall"
  )
  pairwise <- suppressWarnings(cia(
    dat,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    scope = "pairwise"
  ))
  pairwise_constrained <- suppressWarnings(cia(
    dat,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    scope = "pairwise",
    estimator = "vc_constrained"
  ))

  expect_equal(overall$cia[[1]], 1, tolerance = 1e-12)
  expect_equal(overall$cia_raw[[1]], 1, tolerance = 1e-12)
  expect_equal(unname(pairwise["A", "B"]), 1, tolerance = 1e-12)
  expect_equal(unname(pairwise["A", "C"]), 2, tolerance = 1e-12)
  expect_equal(unname(pairwise["B", "C"]), 0, tolerance = 1e-12)
  expect_equal(unname(attr(pairwise, "cia_raw")["A", "C"]), 2, tolerance = 1e-12)
  expect_equal(unname(attr(pairwise, "tau2_raw")["A", "C"]), -1, tolerance = 1e-12)
  expect_true(isTRUE(attr(pairwise, "boundary")["A", "C"]))
  expect_equal(unname(pairwise_constrained["A", "C"]), 1, tolerance = 1e-12)
})

test_that("cia supports raw and constrained reference-scaled estimation", {
  dat <- data.frame(
    subject = rep(c("s1", "s2"), each = 6),
    method = rep(c("A", "A", "B", "B", "C", "C"), times = 2),
    replicate = rep(c("r1", "r2"), times = 6),
    value = c(
      0, 2, 1, 1, 2, 2,
      10, 12, 11, 11, 14, 14
    )
  )

  overall <- cia(
    dat,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    reference = "A",
    scope = "overall"
  )
  overall_raw <- cia(
    dat,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    reference = "A",
    scope = "overall",
    estimator = "mom_unconstrained"
  )
  pairwise <- suppressWarnings(cia(
    dat,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    reference = "A",
    scope = "pairwise"
  ))
  pairwise_raw <- suppressWarnings(cia(
    dat,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    reference = "A",
    scope = "pairwise",
    estimator = "mom_unconstrained"
  ))

  expect_equal(overall$cia[[1]], 8 / 7, tolerance = 1e-12)
  expect_equal(overall$cia_raw[[1]], 8 / 7, tolerance = 1e-12)
  expect_equal(overall$tau2_raw[[1]], -0.25, tolerance = 1e-12)
  expect_true(isTRUE(overall$boundary[[1]]))
  expect_equal(overall$numerator_term[[1]], 4, tolerance = 1e-12)
  expect_equal(overall$denominator_term[[1]], 3.5, tolerance = 1e-12)
  expect_equal(overall_raw$cia[[1]], 8 / 7, tolerance = 1e-12)
  expect_equal(unname(pairwise["A", "B"]), 4, tolerance = 1e-12)
  expect_equal(unname(pairwise["A", "C"]), 4 / 6, tolerance = 1e-12)
  expect_equal(unname(attr(pairwise, "cia_raw")["A", "B"]), 4, tolerance = 1e-12)
  expect_true(isTRUE(attr(pairwise, "boundary")["A", "B"]))
  expect_equal(unname(pairwise_raw["A", "B"]), 4, tolerance = 1e-12)
  expect_equal(unname(pairwise["B", "C"]), NA_real_)
})

test_that("cia summaries and printing expose raw boundary diagnostics", {
  dat <- data.frame(
    subject = rep(c("s1", "s2"), each = 6),
    method = rep(c("A", "A", "B", "B", "C", "C"), times = 2),
    replicate = rep(c("r1", "r2"), times = 6),
    value = c(
      0, 2, 1, 1, 2, 2,
      10, 12, 11, 11, 14, 14
    )
  )

  fit <- suppressWarnings(cia(
    dat,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    reference = "A",
    scope = "pairwise",
    estimator = "vc_constrained"
  ))
  smry <- summary(fit)

  expect_true(all(c("cia_raw", "tau2_raw", "tau2", "boundary") %in% names(smry)))
  expect_true(any(smry$boundary, na.rm = TRUE))
  expect_match(
    paste(capture.output(print(fit)), collapse = "\n"),
    "negative inter-method variability"
  )
  expect_match(
    paste(capture.output(print(smry)), collapse = "\n"),
    "negative inter-method variability"
  )
})

test_that("cia requires within-method replication", {
  dat <- data.frame(
    subject = rep(c("s1", "s2"), each = 2),
    method = rep(c("A", "B"), times = 2),
    replicate = "r1",
    value = c(1, 1, 2, 2)
  )

  expect_error(
    cia(
      dat,
      response = "value",
      subject = "subject",
      method = "method",
      replicate = "replicate",
      scope = "overall"
    ),
    "balanced-replication|replicated readings"
  )
})
