pa_sim_data <- function(n = 120L) {
  set.seed(11)
  dat <- data.frame(
    age = c(runif(n, 0, 50), runif(n, 0, 50)),
    generation = rep(c("Gen1", "Gen2"), each = n)
  )
  eta <- ifelse(dat$generation == "Gen1", 3.8 - 0.05 * dat$age, 3.5 - 0.045 * dat$age)
  dat$pass <- rbinom(2L * n, 1L, plogis(eta))
  dat
}

test_that("prob_agree fits logit reliability curves from data", {
  dat <- pa_sim_data()
  fit <- prob_agree(
    dat,
    response = "pass",
    predictor = "age",
    group = "generation",
    delta = 0.05,
    link = "logit",
    newdata = data.frame(age = c(0, 25, 50)),
    ci = TRUE
  )

  expect_s3_class(fit, "prob_agree")
  expect_equal(names(fit), c("group1", "group2", "age", "prob_agree", "lwr.ci", "upr.ci"))
  expect_equal(fit$age, c(0, 25, 50))
  expect_equal(unique(fit$group1), "Gen1")
  expect_equal(unique(fit$group2), "Gen2")
  expect_true(all(fit$prob_agree >= 0 & fit$prob_agree <= 1))
  expect_true(all(fit$lwr.ci >= 0 & fit$lwr.ci <= 1))
  expect_true(all(fit$upr.ci >= 0 & fit$upr.ci <= 1))
  expect_equal(attr(fit, "group.levels"), c("Gen1", "Gen2"))
  expect_equal(attr(fit, "pairs")$pair, "Gen1 vs Gen2")
  expect_equal(names(attr(fit, "pair_models")[[1]]$coefficients), c("Gen1", "Gen2"))
})

test_that("prob_agree supports probit link and auto grid", {
  dat <- pa_sim_data()
  fit <- prob_agree(
    dat,
    response = "pass",
    predictor = "age",
    group = "generation",
    delta = 0.10,
    link = "probit",
    grid_size = 7,
    ci = FALSE
  )

  expect_equal(nrow(fit), 7)
  expect_false("lwr.ci" %in% names(fit))
  expect_equal(attr(fit, "link"), "probit")
})

test_that("prob_agree supports asymmetric limits", {
  dat <- pa_sim_data()
  fit <- prob_agree(
    dat,
    response = "pass",
    predictor = "age",
    group = "generation",
    limits = c(-0.05, 0.10),
    newdata = c(10, 20),
    ci = FALSE
  )

  expect_equal(attr(fit, "tolerance")$limits, c(-0.05, 0.10))
})

test_that("prob_agree accepts asymmetric delta", {
  dat <- pa_sim_data()
  fit <- prob_agree(
    dat,
    response = "pass",
    predictor = "age",
    group = "generation",
    delta = c(-0.03, 0.05),
    newdata = c(10, 20),
    ci = FALSE
  )

  expect_equal(names(fit), c("group1", "group2", "age", "prob_agree"))
  expect_equal(attr(fit, "tolerance")$limits, c(-0.03, 0.05))
})

test_that("prob_agree computes all pairwise comparisons for more than two groups", {
  set.seed(12)
  n <- 70L
  dat <- data.frame(
    age = rep(runif(n, 0, 45), times = 4),
    generation = rep(paste0("Gen", 1:4), each = n)
  )
  gen_id <- match(dat$generation, paste0("Gen", 1:4))
  dat$pass <- rbinom(
    nrow(dat),
    1L,
    plogis(c(3.8, 3.6, 3.4, 3.2)[gen_id] - c(0.04, 0.042, 0.044, 0.046)[gen_id] * dat$age)
  )

  fit <- prob_agree(
    dat,
    response = "pass",
    predictor = "age",
    group = "generation",
    delta = 0.05,
    newdata = c(0, 20, 40),
    ci = FALSE
  )

  expect_equal(nrow(fit), 18)
  expect_equal(nrow(attr(fit, "pairs")), 6)
  expect_equal(unique(paste(fit$group1, fit$group2, sep = " vs ")), attr(fit, "pairs")$pair)
  expect_true(all(fit$prob_agree >= 0 & fit$prob_agree <= 1))
  skip_if_not_installed("ggplot2")
  expect_warning(p <- plot(fit), NA)
  expect_s3_class(p, "ggplot")
  expect_s3_class(plot(fit, style = "curve"), "ggplot")
  expect_s3_class(plot(fit, style = "facet"), "ggplot")
  expect_s3_class(plot(fit, style = "heatmap"), "ggplot")
})

test_that("prob_agree validates data and tolerances", {
  dat <- pa_sim_data()
  expect_error(prob_agree(dat, "pass", "age", "generation"), "tolerance")
  expect_error(
    prob_agree(dat, "pass", "age", "generation", delta = 0.1, limits = c(-0.1, 0.1)),
    "mutually exclusive"
  )
  expect_error(
    prob_agree(dat, "pass", "age", "generation", limits = c(0.1, -0.1)),
    class = "matrixCorr_arg_error"
  )
  expect_error(
    prob_agree(dat, "pass", "age", "generation", delta = -0.1),
    class = "matrixCorr_arg_error"
  )
  expect_error(
    prob_agree(transform(dat, generation = "one"), "pass", "age", "generation", delta = 0.1),
    class = "matrixCorr_arg_error"
  )
  expect_error(
    prob_agree(transform(dat, pass = 2), "pass", "age", "generation", delta = 0.1),
    class = "matrixCorr_arg_error"
  )
})

test_that("prob_agree S3 methods work", {
  skip_if_not_installed("ggplot2")
  dat <- pa_sim_data()
  fit <- prob_agree(
    dat,
    response = "pass",
    predictor = "age",
    group = "generation",
    delta = 0.05,
    newdata = c(0, 20, 40),
    ci = TRUE
  )

  expect_true(length(capture.output(print(fit))) > 0L)
  sm <- summary(fit)
  expect_s3_class(sm, "summary.prob_agree")
  expect_true(all(c(
    "group1", "group2", "age", "prob_agree", "lwr.ci", "upr.ci"
  ) %in% names(sm)))
  expect_s3_class(plot(fit), "ggplot")
  expect_s3_class(plot(fit, style = "curve"), "ggplot")
  expect_s3_class(plot(fit, style = "facet"), "ggplot")
})
