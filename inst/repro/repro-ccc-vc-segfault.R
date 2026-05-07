#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(matrixCorr))

run_case <- function(label, expr) {
  cat(sprintf("\n== %s ==\n", label))
  gc()
  flush.console()
  force(expr)
  cat(sprintf("OK: %s\n", label))
}

run_case("Case 1: pairwise, no time", {
  set.seed(123)
  n_subj <- 500L
  id <- factor(rep(seq_len(n_subj), each = 2L))
  method <- factor(rep(c("A", "B"), times = n_subj))

  sigA <- 1
  sigE <- 0.5
  biasB <- 0.2

  u <- rnorm(n_subj, 0, sqrt(sigA))[as.integer(id)]
  e <- rnorm(n_subj * 2L, 0, sqrt(sigE))
  y <- (method == "B") * biasB + u + e

  df <- data.frame(y, id, method)

  cfit <- ccc_rm_reml(
    df,
    response = "y",
    subject = "id",
    method = "method",
    ci = TRUE
  )

  print(cfit)
})

run_case("Case 2: method + time + vc_select = auto", {
  set.seed(102)
  n_subj <- 60L
  n_time <- 8L

  id <- factor(rep(seq_len(n_subj), each = 2L * n_time))
  time <- factor(rep(rep(seq_len(n_time), times = 2L), times = n_subj))
  method <- factor(rep(rep(c("A", "B"), each = n_time), times = n_subj))

  sigA <- 0.6
  sigAM <- 0.3
  sigAT <- 0.5
  sigE <- 0.4
  biasB <- 0.2

  u_i <- rnorm(n_subj, 0, sqrt(sigA))
  u <- u_i[as.integer(id)]

  sm <- interaction(id, method, drop = TRUE)
  w_im_lv <- rnorm(nlevels(sm), 0, sqrt(sigAM))
  w_im <- w_im_lv[as.integer(sm)]

  st <- interaction(id, time, drop = TRUE)
  g_it_lv <- rnorm(nlevels(st), 0, sqrt(sigAT))
  g_it <- g_it_lv[as.integer(st)]

  e <- rnorm(length(id), 0, sqrt(sigE))
  y <- (method == "B") * biasB + u + w_im + g_it + e

  dat_both <- data.frame(y, id, method, time)

  fit_both <- ccc_rm_reml(
    dat_both,
    "y",
    "id",
    method = "method",
    time = "time",
    vc_select = "auto",
    verbose = TRUE
  )

  print(fit_both)
})

run_case("Vignette ICC consistency case", {
  set.seed(50)
  n_id <- 14L
  n_time <- 4L

  dat <- expand.grid(
    id = factor(seq_len(n_id)),
    time = factor(seq_len(n_time)),
    method = factor(c("A", "B"))
  )

  dat$time_index <- as.integer(dat$time)

  subj <- rnorm(n_id, sd = 1.0)[dat$id]
  subject_method <- rnorm(n_id * 2L, sd = 0.25)
  sm <- subject_method[(as.integer(dat$id) - 1L) * 2L + as.integer(dat$method)]
  subject_time <- rnorm(n_id * n_time, sd = 0.75)
  st <- subject_time[(as.integer(dat$id) - 1L) * n_time + as.integer(dat$time)]

  dat$y <- subj + sm + st + 0.35 * (dat$method == "B") +
    rnorm(nrow(dat), sd = 0.35)

  fit_icc_cons <- icc_rm_reml(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time_index",
    type = "consistency",
    ci = TRUE
  )

  print(fit_icc_cons)
})

run_case("Vignette CCC no-CI time_index case", {
  set.seed(50)
  n_id <- 14L
  n_time <- 4L

  dat <- expand.grid(
    id = factor(seq_len(n_id)),
    time = factor(seq_len(n_time)),
    method = factor(c("A", "B"))
  )

  dat$time_index <- as.integer(dat$time)

  subj <- rnorm(n_id, sd = 1.0)[dat$id]
  subject_method <- rnorm(n_id * 2L, sd = 0.25)
  sm <- subject_method[(as.integer(dat$id) - 1L) * 2L + as.integer(dat$method)]
  subject_time <- rnorm(n_id * n_time, sd = 0.75)
  st <- subject_time[(as.integer(dat$id) - 1L) * n_time + as.integer(dat$time)]

  dat$y <- subj + sm + st + 0.35 * (dat$method == "B") +
    rnorm(nrow(dat), sd = 0.35)

  fit_ccc_reml <- ccc_rm_reml(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time_index",
    ci = FALSE
  )

  print(fit_ccc_reml)
})
