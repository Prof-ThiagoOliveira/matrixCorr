# AR(1) simulator for within-subject residuals
ar1_sim <- function(n, rho, sd = 1) {
  z <- rnorm(n)
  e <- numeric(n)
  e[1] <- z[1] * sd
  if (n > 1) for (t in 2:n) e[t] <- rho * e[t - 1] + sqrt(1 - rho^2) * z[t] * sd
  e
}

# Two-method data with known difference model:
# d_{s,t} = (y2 - y1) = mu + b_s + e_{s,t},
# with b_s ~ N(0, sig_s^2), e ~ N(0, sig_e^2) or AR(1)
# ==> truth: bias = mu; sd_loa = sqrt(sig_s^2 + sig_e^2)
sim_two_method_known <- function(S = 60L, Tm = 25L,
                                 mu = 1.25, sig_s = 1.7, sig_e = 2.3,
                                 rho = NULL, seed = 1L) {
  set.seed(seed)
  subj <- rep(seq_len(S), each = Tm)
  time <- rep(seq_len(Tm), times = S)
  b_s  <- rnorm(S, 0, sig_s)
  b    <- b_s[subj]
  e    <- if (is.null(rho)) rnorm(S * Tm, 0, sig_e) else
    unlist(lapply(split(seq_along(time), subj), function(ix) ar1_sim(length(ix), rho, sig_e)))

  # Baseline cancels in differences; include to mimic realistic responses
  baseline <- rnorm(S * Tm, 10, 3)
  y1 <- baseline
  y2 <- baseline + mu + b + e

  data.frame(
    y = c(y1, y2),
    subject = rep(subj, 2),
    method  = factor(rep(c("M1", "M2"), each = S * Tm), levels = c("M1", "M2")),
    time    = rep(time, 2),
    check.names = FALSE
  )
}

# Multi-method data for matrix invariants
# M2 ~= M1 + small noise; M3 = -M1; M4 = unrelated
sim_multi_method <- function(S = 30L, Tm = 15L, seed = 11L) {
  set.seed(seed)
  subj <- rep(seq_len(S), each = Tm)
  time <- rep(seq_len(Tm), times = S)

  mu_s <- rnorm(S, mean = 0, sd = 8)
  true <- mu_s[subj]
  e0   <- rnorm(length(true), 0, 2)

  y1 <- true + e0
  y2 <- y1 + rnorm(length(true), 0, 0.01)
  y3 <- -y1
  y4 <- rnorm(length(true), 3, 6)

  out <- rbind(
    data.frame(y = y1, subject = subj, method = "M1", time = time),
    data.frame(y = y2, subject = subj, method = "M2", time = time),
    data.frame(y = y3, subject = subj, method = "M3", time = time),
    data.frame(y = y4, subject = subj, method = "M4", time = time),
    make.row.names = FALSE
  )
  out$method <- factor(out$method, levels = c("M1", "M2", "M3", "M4"))
  out
}

sim_two_method_from_pairs <- function(means, diffs, S, Tm) {
  stopifnot(length(means) == length(diffs), length(means) == S * Tm)
  subj <- rep(seq_len(S), each = Tm)
  time <- rep(seq_len(Tm), times = S)
  y1 <- means - diffs / 2
  y2 <- means + diffs / 2

  data.frame(
    y = c(y1, y2),
    subject = rep(subj, 2),
    method = factor(rep(c("A", "B"), each = length(y1)), levels = c("A", "B")),
    time = rep(time, 2),
    check.names = FALSE
  )
}

test_that("two-method contrast matches factor order (method2 - method1)", {
  dat <- sim_two_method_known(S = 50, Tm = 20, mu = 1.0, sig_s = 1.2, sig_e = 1.5, seed = 1)
  dat$method <- factor(dat$method, levels = c("M1","M2"))
  fit <- ba_rm(data = dat, response="y", subject="subject",
                               method="method", time="time", use_ar1=FALSE)
  expect_gt(as.numeric(fit$mean.diffs), 0)  # should be ~ +1.0
})

test_that("Two-method BA recovers bias and sd_loa under i.i.d. residuals", {
  dat <- sim_two_method_known(S = 80, Tm = 50, mu = 1.25, sig_s = 1.7,
                              sig_e = 2.3, rho = NULL, seed = 100)
  truth_sd <- sqrt(1.7^2 + 2.3^2)

  fit <- ba_rm(
    data = dat,
    response = "y", subject = "subject", method = "method", time = "time",
    include_slope = FALSE, use_ar1 = FALSE, loa_multiplier = 1.96, conf_level = 0.95
  )

  expect_s3_class(fit, "ba_repeated")

  # Bias and sd_loa close to truth
  expect_equal(as.numeric(fit$mean.diffs), 1.25, tolerance = 0.08)
  expect_equal(as.numeric(fit$critical.diff) / as.numeric(fit$loa_multiplier), truth_sd, tolerance = 0.08)

  # LoA consistent with sd and 'loa_multiplier'
  sd_hat <- as.numeric(fit$critical.diff) / as.numeric(fit$loa_multiplier)
  expect_equal(as.numeric(fit$upper.limit) - as.numeric(fit$mean.diffs), 1.96 * sd_hat, tolerance = 1e-6)
  expect_equal(as.numeric(fit$mean.diffs) - as.numeric(fit$lower.limit), 1.96 * sd_hat, tolerance = 1e-6)

  # Summary object sanity
  sm <- summary(fit)
  expect_s3_class(sm, "summary.ba_repeated")
  expect_true(all(c("bias","sd_loa","loa_low","loa_up","width","n") %in% names(sm)))
  expect_false(any(c("ar1_rho", "ar1_estimated") %in% names(sm)))
})

test_that("Two-method BA with AR(1) honours supplied rho and recovers truth", {
  skip_on_cran()

  mu    <- 0.8
  sig_s <- 1.2
  sig_e <- 1.5
  rho   <- 0.45
  dat <- sim_two_method_known(S = 70, Tm = 40, mu = mu, sig_s = sig_s, sig_e = sig_e, rho = rho, seed = 200)
  truth_sd <- sqrt(sig_s^2 + sig_e^2)  # definition of sd_loa in docs

  fit <- ba_rm(
    data = dat,
    response = "y", subject = "subject", method = "method", time = "time",
    include_slope = FALSE, use_ar1 = TRUE, ar1_rho = rho, loa_multiplier = 2.0, conf_level = 0.9
  )

  expect_true(isTRUE(fit$use_ar1))
  expect_false(is.na(fit$ar1_rho))
  expect_equal(as.numeric(fit$ar1_rho), rho, tolerance = 1e-12)
  expect_false(isTRUE(fit$ar1_estimated))  # supplied, so not estimated

  expect_equal(as.numeric(fit$mean.diffs), mu,   tolerance = 0.08)
  expect_equal(as.numeric(fit$critical.diff) / 2.0, truth_sd, tolerance = 0.09)
})

test_that("Two-method BA estimates positive AR(1) dependence on short panels", {
  skip_on_cran()

  mu <- 0.35
  sig_s <- sqrt(0.28^2 + 0.38^2)
  sig_e <- sqrt(0.45^2 + 0.60^2)
  rho <- 0.60
  dat <- sim_two_method_known(
    S = 120,
    Tm = 10,
    mu = mu,
    sig_s = sig_s,
    sig_e = sig_e,
    rho = rho,
    seed = 20260401L
  )

  fit <- ba_rm(
    data = dat,
    response = "y", subject = "subject", method = "method", time = "time",
    include_slope = FALSE, use_ar1 = TRUE, ar1_rho = NA_real_,
    loa_multiplier = 1.96, conf_level = 0.95
  )

  expect_true(isTRUE(fit$use_ar1))
  expect_true(isTRUE(fit$ar1_estimated))
  expect_gt(as.numeric(fit$ar1_rho), 0.30)
  expect_equal(as.numeric(fit$ar1_rho), rho, tolerance = 0.20)
  expect_true(abs(as.numeric(fit$mean.diffs) - mu) < 0.08)
  expect_equal(as.numeric(fit$critical.diff) / 1.96, sqrt(sig_s^2 + sig_e^2), tolerance = 0.08)
})

test_that("AR(1) requests simplify to iid with a warning when needed for some pairs", {
  set.seed(123)
  S <- 40L
  Tm <- 50L
  methods <- c("A", "B", "C")
  rho <- 0.4

  ar1_sim <- function(n, rho, sd = 1) {
    z <- rnorm(n)
    e <- numeric(n)
    e[1] <- z[1] * sd
    if (n > 1) for (t in 2:n) e[t] <- rho * e[t - 1] + sqrt(1 - rho^2) * z[t] * sd
    e
  }

  subj <- rep(seq_len(S), each = Tm)
  time <- rep(seq_len(Tm), times = S)
  mu_s <- rnorm(S, 50, 7)
  trend <- rep(seq_len(Tm) - mean(seq_len(Tm)), times = S) * 0.8
  true <- mu_s[subj] + trend
  bias <- c(A = 0, B = 1.5, C = -0.5)
  prop <- c(A = 0.00, B = 0.00, C = 0.10)

  make_method <- function(meth, sd = 3) {
    e <- unlist(lapply(split(seq_along(time), subj),
                       function(ix) ar1_sim(length(ix), rho, sd)))
    y <- true * (1 + prop[meth]) + bias[meth] + e
    data.frame(y = y, subject = subj, method = meth, time = time,
               check.names = FALSE)
  }

  dat <- do.call(rbind, lapply(methods, make_method))
  dat$method <- factor(dat$method, levels = methods)

 fit <- ba_rm(
     response = dat$y, subject = dat$subject, method = dat$method, time = dat$time,
    include_slope = FALSE, use_ar1 = TRUE, ar1_rho = rho)


  expect_s3_class(fit, "ba_repeated_matrix")
  simplified <- fit$residual_model == "iid"
  expect_true(all(is.na(fit$ar1_rho_pair[simplified])))

  sm <- summary(fit)
  expect_true(all(c("residual_model", "ar1_rho") %in% names(sm)))
})

test_that("AR(1) fallback warning text is explicit", {
  expect_warning(
    matrixCorr:::.warn_ba_rm_ar1_fallback(),
    "Requested AR\\(1\\) residual structure could not be fit for this fit; using iid residuals instead\\."
  )

  expect_warning(
    matrixCorr:::.warn_ba_rm_ar1_fallback(c("A-B", "B-C", "A-B")),
    "Requested AR\\(1\\) residual structure could not be fit for pair\\(s\\): A-B, B-C; using iid residuals instead\\."
  )
})

test_that("Edge case: constant difference gives zero sd_loa and degenerate LoA", {
  # Small jitter to avoid numerical singularities in optimisers; still near-zero SD
  set.seed(1)
  S <- 25; Tm <- 8; mu <- 5
  subj <- rep(seq_len(S), each = Tm)
  time <- rep(seq_len(Tm), times = S)
  base <- rnorm(S * Tm, 12, 2)
  y1 <- base
  y2 <- base + mu + rnorm(S * Tm, 0, 1e-6)

  dat <- data.frame(y = c(y1,y2), subject = rep(subj, 2), method = rep(c("A","B"), each = length(y1)),
                    time = rep(time, 2))
  dat$method <- factor(dat$method, levels = c("A","B"))

  fit <- ba_rm(
    data = dat, response = "y", subject = "subject", method = "method", time = "time",
    include_slope = FALSE, use_ar1 = FALSE
  )

  expect_equal(as.numeric(fit$mean.diffs), mu, tolerance = 1e-3)
  expect_lt(as.numeric(fit$critical.diff), 1e-2)  # ~ zero
  expect_equal(as.numeric(fit$lower.limit), as.numeric(fit$upper.limit), tolerance = 1e-2)
})

test_that("Temporary positive residual-variance initialization remains only an EM starting heuristic", {
  S <- 10L
  Tm <- 2L
  subject_means <- seq(12, 30, length.out = S)
  means <- rep(subject_means, each = Tm)
  slope_true <- 0.25
  diffs <- rep(0.8 + slope_true * (subject_means - mean(subject_means)), each = Tm)
  dat <- sim_two_method_from_pairs(means, diffs, S = S, Tm = Tm)

  fit_cpp <- matrixCorr:::bland_altman_repeated_em_ext_cpp(
    y = dat$y,
    subject = dat$subject,
    method = as.integer(dat$method == "B") + 1L,
    time = dat$time,
    include_slope = TRUE,
    use_ar1 = FALSE
  )
  expect_match(
    fit_cpp$warn,
    "Constrained boundary fit used with sigma2_subject fixed at 0.",
    perl = TRUE
  )
  expect_true(isTRUE(fit_cpp$converged))
  expect_true(is.finite(as.numeric(fit_cpp$sigma2_subject)))
  expect_true(is.finite(as.numeric(fit_cpp$sigma2_resid)))

  fit <- ba_rm(
    data = dat,
    response = "y", subject = "subject", method = "method", time = "time",
    include_slope = TRUE, use_ar1 = FALSE
  )
  expect_true(is.finite(as.numeric(fit$beta_slope)))
  expect_equal(as.numeric(fit$beta_slope), slope_true, tolerance = 0.05)
})

test_that("Slope scale selection uses SD when paired means have ordinary spread", {
  normal_m <- seq(1, 8)
  info_normal <- matrixCorr:::ba_rm_slope_scale_cpp(normal_m)
  expect_identical(info_normal$source, "sd")
  expect_equal(info_normal$s_sd, stats::sd(normal_m), tolerance = 1e-15)
  expect_equal(info_normal$s_iqr, stats::IQR(normal_m) / 1.349, tolerance = 1e-15)
  expect_equal(
    info_normal$s_mad,
    stats::mad(normal_m, center = stats::median(normal_m), constant = 1.4826, na.rm = TRUE),
    tolerance = 1e-15
  )
  expect_equal(info_normal$s_m_star, info_normal$s_sd, tolerance = 1e-15)
})

test_that("Slope scale selection falls back to IQR when SD is near zero", {
  near_constant_m <- c(rep(0, 16), rep(5e-324, 16))
  info_near <- matrixCorr:::ba_rm_slope_scale_cpp(near_constant_m)

  expect_identical(info_near$source, "iqr")
  expect_equal(info_near$s_sd, 0, tolerance = 0)
  expect_true(is.finite(info_near$s_iqr) && info_near$s_iqr > 0)
  expect_equal(info_near$s_m_star, info_near$s_iqr, tolerance = 0)
  expect_true(info_near$s_sd <= info_near$tau * info_near$s_ref)
})

test_that("Slope fit is recovered when pair means have ordinary spread", {
  set.seed(321)
  S <- 18L
  Tm <- 4L
  n <- S * Tm
  means <- seq(10, 40, length.out = n)
  subj_eff <- rep(rnorm(S, 0, 0.15), each = Tm)
  slope_true <- 0.35
  diffs <- 0.75 + slope_true * (means - mean(means)) + subj_eff + rnorm(n, 0, 0.04)
  dat <- sim_two_method_from_pairs(means, diffs, S = S, Tm = Tm)

  fit <- ba_rm(
    data = dat,
    response = "y", subject = "subject", method = "method", time = "time",
    include_slope = TRUE, use_ar1 = FALSE
  )

  expect_true(is.finite(as.numeric(fit$beta_slope)))
  expect_equal(as.numeric(fit$beta_slope), slope_true, tolerance = 0.05)
})

test_that("Slope scaling errors when paired means are fully degenerate", {
  expect_error(
    matrixCorr:::ba_rm_slope_scale_cpp(rep(25, 32)),
    "The proportional-bias slope is not estimable because the paired means are near-degenerate on the observed data scale.",
    fixed = TRUE
  )

  set.seed(323)
  S <- 14L
  Tm <- 4L
  n <- S * Tm
  means <- rep(25, n)
  diffs <- 1.1 + rep(rnorm(S, 0, 0.12), each = Tm) + rnorm(n, 0, 0.03)
  dat <- sim_two_method_from_pairs(means, diffs, S = S, Tm = Tm)

  expect_error(
    ba_rm(
      data = dat,
      response = "y", subject = "subject", method = "method", time = "time",
      include_slope = TRUE, use_ar1 = FALSE
    ),
    "The proportional-bias slope is not estimable because the paired means are near-degenerate on the observed data scale.",
    fixed = TRUE
  )
})

test_that("Ordinary repeated-measures cases still converge successfully", {
  dat <- sim_two_method_known(S = 40L, Tm = 12L, mu = 1.1, sig_s = 1.0, sig_e = 0.8, seed = 913)

  fit_cpp <- matrixCorr:::bland_altman_repeated_em_ext_cpp(
    y = dat$y,
    subject = dat$subject,
    method = as.integer(dat$method == "M2") + 1L,
    time = dat$time,
    include_slope = FALSE,
    use_ar1 = FALSE
  )

  expect_true(isTRUE(fit_cpp$converged))
  expect_true(is.finite(as.numeric(fit_cpp$sigma2_subject)))
  expect_true(is.finite(as.numeric(fit_cpp$sigma2_resid)))

  fit <- ba_rm(
    data = dat,
    response = "y", subject = "subject", method = "method", time = "time",
    include_slope = FALSE, use_ar1 = FALSE
  )

  expect_s3_class(fit, "ba_repeated")
  expect_true(is.finite(as.numeric(fit$mean.diffs)))
  expect_true(is.finite(as.numeric(fit$critical.diff)))
  expect_gt(as.numeric(fit$critical.diff), 0)
  expect_true(is.finite(as.numeric(fit$sigma2_subject)))
  expect_true(is.finite(as.numeric(fit$sigma2_resid)))
})

test_that("Ordinary slope cases still converge successfully", {
  set.seed(914)
  S <- 18L
  Tm <- 4L
  n <- S * Tm
  means <- seq(10, 40, length.out = n)
  subj_eff <- rep(rnorm(S, 0, 0.15), each = Tm)
  slope_true <- 0.35
  diffs <- 0.75 + slope_true * (means - mean(means)) + subj_eff + rnorm(n, 0, 0.04)
  dat <- sim_two_method_from_pairs(means, diffs, S = S, Tm = Tm)

  fit_cpp <- matrixCorr:::bland_altman_repeated_em_ext_cpp(
    y = dat$y,
    subject = dat$subject,
    method = as.integer(dat$method == "B") + 1L,
    time = dat$time,
    include_slope = TRUE,
    use_ar1 = FALSE
  )

  expect_true(isTRUE(fit_cpp$converged))
  expect_true(is.finite(as.numeric(fit_cpp$beta_slope)))

  fit <- ba_rm(
    data = dat,
    response = "y", subject = "subject", method = "method", time = "time",
    include_slope = TRUE, use_ar1 = FALSE
  )

  expect_true(is.finite(as.numeric(fit$beta_slope)))
  expect_equal(as.numeric(fit$beta_slope), slope_true, tolerance = 0.05)
})

test_that("Pairwise matrix algebraic invariants hold", {
  dat <- sim_multi_method(S = 24, Tm = 12, seed = 99)

  fit <- ba_rm(
    data = dat,
    response = "y", subject = "subject", method = "method", time = "time",
    include_slope = FALSE, use_ar1 = FALSE, loa_multiplier = 1.96, conf_level = 0.95
  )

  expect_s3_class(fit, "ba_repeated_matrix")
  m <- length(fit$methods)

  # Symmetry / antisymmetry conditions
  expect_true(all(is.na(diag(fit$bias))))
  expect_true(all(abs(fit$bias + t(fit$bias)) < 1e-8, na.rm = TRUE))

  expect_true(all(abs(fit$sd_loa - t(fit$sd_loa)) < 1e-12, na.rm = TRUE))
  expect_true(all(abs(fit$width  - t(fit$width )) < 1e-12, na.rm = TRUE))

  # LoA sign conventions (row - column)
  for (i in 1:(m-1)) for (j in (i+1):m) {
    expect_equal(fit$loa_lower[j,i], -fit$loa_upper[i,j], tolerance = 1e-10)
    expect_equal(fit$loa_upper[j,i], -fit$loa_lower[i,j], tolerance = 1e-10)
  }

  # Width equals 2 * loa_multiplier * sd_loa
  expect_equal(fit$width, 2 * fit$loa_multiplier * fit$sd_loa, tolerance = 1e-10)
})

test_that("Pairwise result equals two-method fit for the same pair", {
  dat <- sim_multi_method(S = 30, Tm = 15, seed = 7)
  # Compare M1 vs M2
  dat12 <- subset(dat, method %in% c("M1","M2"))

  fit2 <- ba_rm(
    data = dat12, response = "y", subject = "subject", method = "method", time = "time",
    include_slope = FALSE, use_ar1 = FALSE, loa_multiplier = 1.96, conf_level = 0.95
  )

  fitN <- ba_rm(
    data = dat, response = "y", subject = "subject", method = "method", time = "time",
    include_slope = FALSE, use_ar1 = FALSE, loa_multiplier = 1.96, conf_level = 0.95
  )

  i <- match("M1", fitN$methods); j <- match("M2", fitN$methods)
  expect_equal(as.numeric(fit2$mean.diffs), fitN$bias[i,j], tolerance = 1e-8)
  expect_equal(as.numeric(fit2$critical.diff) / as.numeric(fit2$loa_multiplier), fitN$sd_loa[i,j], tolerance = 1e-8)
  expect_equal(as.numeric(fit2$lower.limit), fitN$loa_lower[i,j], tolerance = 1e-8)
  expect_equal(as.numeric(fit2$upper.limit), fitN$loa_upper[i,j], tolerance = 1e-8)
})

test_that("n matrix equals number of subject-time complete pairs per contrast", {
  # Create imbalance: drop some method-time rows
  dat <- sim_multi_method(S = 12, Tm = 10, seed = 123)
  # Remove random 15% of M3 rows to induce missingness
  set.seed(123)
  drop_idx <- which(dat$method == "M3")
  rm <- sample(drop_idx, size = floor(0.15 * length(drop_idx)))
  dat2 <- dat[-rm, ]

  fit <- ba_rm(
    data = dat2, response = "y", subject = "subject", method = "method", time = "time",
    include_slope = FALSE, use_ar1 = FALSE
  )

  # Manual complete-case counts for M1 vs M3
  cc <- merge(
    subset(dat2, method == "M1")[, c("subject","time")],
    subset(dat2, method == "M3")[, c("subject","time")],
    by = c("subject","time"),
    all = FALSE
  )
  i <- match("M1", fit$methods); j <- match("M3", fit$methods)
  expect_equal(as.integer(fit$n[i,j]), nrow(cc))
})

test_that("native complete-pair counter matches the expected matched-pair total", {
  dat <- data.frame(
    y = c(1, 2, 3, 4, 5, NA, 7, 8),
    subject = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
    method = c(1L, 2L, 1L, 1L, 1L, 2L, 2L, 1L),
    time = c(1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L)
  )

  expect_equal(
    matrixCorr:::ba_rm_complete_pairs_cpp(dat$y, dat$subject, dat$method, dat$time),
    2L
  )
})

test_that("Column mapping works both with data= and with vectors", {
  dat <- sim_multi_method(S = 10, Tm = 6, seed = 42)
  # With data=
  fit1 <- ba_rm(
    data = dat, response = "y", subject = "subject", method = "method", time = "time"
  )
  expect_true(inherits(fit1, "ba_repeated_matrix"))

  # With vectors (no data=)
  fit2 <- ba_rm(
    response = dat$y, subject = dat$subject, method = dat$method, time = dat$time
  )
  expect_true(inherits(fit2, "ba_repeated_matrix"))

  # Stored data_long/mapping present for plotting
  expect_true(!is.null(fit1$data_long) && !is.null(fit1$mapping))
  expect_true(!is.null(fit2$data_long) && !is.null(fit2$mapping))
})

test_that("summary/print produce expected classes and do not error", {
  dat <- sim_multi_method(S = 10, Tm = 6, seed = 2)

  fit_mat <- ba_rm(
    data = dat, response = "y", subject = "subject", method = "method", time = "time"
  )
  out_fit_mat <- capture.output(print(fit_mat))
  expect_false(any(grepl("bias_lwr|sigma2_subject|ar1_rho|Confidence intervals|Model details",
                         out_fit_mat)))

  sm_mat <- summary(fit_mat)
  expect_s3_class(sm_mat, "summary.ba_repeated_matrix")
  expect_true(all(c("method1","method2","bias","sd_loa","loa_low","loa_up","width","n") %in% names(sm_mat)))
  expect_false(any(c("ar1_rho", "ar1_estimated") %in% names(sm_mat)))
  out_mat <- capture.output(print(sm_mat))
  expect_true(any(grepl("^Agreement estimates$", out_mat)))
  expect_true(any(grepl("^Confidence intervals$", out_mat)))
  expect_true(any(grepl("^Model details$", out_mat)))
  expect_false(any(grepl("ar1_rho|ar1_estimated", out_mat)))

  # Two-method
  dat12 <- subset(dat, method %in% c("M1","M2"))
  fit2  <- ba_rm(
    data = dat12, response = "y", subject = "subject", method = droplevels(dat12$method), time = "time"
  )
  sm2 <- summary(fit2)
  expect_s3_class(sm2, "summary.ba_repeated")
  expect_true(all(c("bias","sd_loa","loa_low","loa_up","width","n") %in% names(sm2)))
  expect_false(any(c("ar1_rho", "ar1_estimated") %in% names(sm2)))
  out2 <- capture.output(print(sm2))
  expect_true(any(grepl("^Agreement estimates$", out2)))
  expect_true(any(grepl("^Confidence intervals$", out2)))
  expect_true(any(grepl("^Model details$", out2)))
  expect_false(any(grepl("ar1_rho|ar1_estimated", out2)))
})

test_that("two-method summary omits AR columns when final residual model is iid", {
  fit <- structure(
    list(
      based.on = 12L,
      mean.diffs = 0.8,
      critical.diff = 1.96,
      loa_multiplier = 1.96,
      lower.limit = -1.16,
      upper.limit = 2.76,
      CI.lines = c(
        "mean.diff.ci.lower" = 0.5,
        "mean.diff.ci.upper" = 1.1,
        "lower.limit.ci.lower" = -1.6,
        "lower.limit.ci.upper" = -0.7,
        "upper.limit.ci.lower" = 2.3,
        "upper.limit.ci.upper" = 3.2
      ),
      include_slope = FALSE,
      sigma2_subject = 0.4,
      sigma2_resid = 1.1,
      use_ar1 = TRUE,
      residual_model = "iid",
      ar1_rho = NA_real_,
      ar1_estimated = FALSE
    ),
    class = "ba_repeated"
  )
  attr(fit, "conf.level") <- 0.95

  sm <- summary(fit)

  expect_false(any(c("ar1_rho", "ar1_estimated") %in% names(sm)))
})

test_that("plot methods return a ggplot object and do not error", {
  skip_if_not_installed("ggplot2")
  dat <- sim_multi_method(S = 12, Tm = 8, seed = 5)

  # Matrix plot
  fit_mat <- ba_rm(
    data = dat, response = "y", subject = "subject", method = "method", time = "time",
    include_slope = FALSE
  )
  p_mat <- plot(fit_mat, smoother = "lm", facet_scales = "free_y")
  expect_s3_class(p_mat, "ggplot")

  # Two-method plot
  dat12 <- subset(dat, method %in% c("M1","M2"))
  fit2  <- ba_rm(
    data = dat12, response = "y", subject = "subject", method = droplevels(dat12$method), time = "time",
    include_slope = TRUE
  )
  p2 <- plot(fit2, smoother = "loess", symmetrize_y = TRUE)
  expect_s3_class(p2, "ggplot")
})

test_that("Input argument validation works", {
  dat <- sim_multi_method(S = 6, Tm = 4, seed = 1)

  # bad 'loa_multiplier'
  expect_error(
    ba_rm(data = dat, response = "y", subject = "subject", method = "method", time = "time",
                          loa_multiplier = -1),
    "`loa_multiplier` must be a positive scalar.", fixed = TRUE
  )

  # bad conf_level
  expect_error(
    ba_rm(data = dat, response = "y", subject = "subject", method = "method", time = "time",
                          conf_level = 1.1),
    "`conf_level` must be in (0,1).", fixed = TRUE
  )

  # AR1 rho bounds when use_ar1 = TRUE
  expect_error(
    ba_rm(data = dat, response = "y", subject = "subject", method = "method", time = "time",
                          use_ar1 = TRUE, ar1_rho = 1.2),
    "`ar1_rho` must be in (-0.999, 0.999).", fixed = TRUE
  )

  # Need at least 2 methods
  dat1 <- subset(dat, method == "M1")
  expect_error(
    ba_rm(data = dat1, response = "y", subject = "subject", method = "method", time = "time"),
    "Need at least 2 distinct methods in `method`.", fixed = TRUE
  )

  # Wrong column name
  expect_error(
    ba_rm(data = dat, response = "y", subject = "subject", method = "NOPE", time = "time"),
    "Column 'NOPE' not found in `data`.", fixed = TRUE
  )
})

test_that("two-method path requires at least two matched pairs", {
  dat <- data.frame(
    y = c(10, 11, 12),
    subject = c(1L, 1L, 1L),
    method = factor(c("A", "B", "A"), levels = c("A", "B")),
    time = c(1L, 1L, 2L),
    check.names = FALSE
  )

  expect_error(
    ba_rm(dat$y, dat$subject, dat$method, dat$time,
                          include_slope = FALSE, use_ar1 = FALSE),
    "at least two subject-time matched pairs",
    fixed = FALSE
  )
})

test_that("pairwise matrix drops contrasts with <2 matched pairs", {
  dat <- data.frame(
    y = c(1, 2, 3, 4, 5, 6, 7),
    subject = c(1L, 1L, 1L, 1L, 1L, 2L, 2L),
    method = factor(c("A", "B", "C", "A", "C", "A", "C"), levels = c("A", "B", "C")),
    time = c(1L, 1L, 1L, 2L, 2L, 1L, 1L),
    check.names = FALSE
  )

  fit <- ba_rm(dat$y, dat$subject, dat$method, dat$time,
                               include_slope = FALSE, use_ar1 = FALSE)

  expect_true(is.na(fit$bias["A", "B"]))
  expect_equal(fit$n["A", "B"], 1L)
  expect_true(is.na(fit$sd_loa["A", "B"]))
  expect_true(is.na(fit$loa_lower["A", "B"]))
  expect_true(is.na(fit$loa_upper["A", "B"]))
  expect_equal(fit$n["A", "C"], 3L)
  expect_true(is.finite(fit$bias["A", "C"]))
})

test_that("backend refactor preserves key model-based BA estimands", {
  extract_fit <- function(fit) {
    with(fit, list(
      bias_mu0 = bias_mu0,
      bias_lwr = bias_lwr,
      bias_upr = bias_upr,
      sigma2_subject = sigma2_subject,
      sigma2_resid = sigma2_resid,
      sd_loa = sd_loa,
      loa_lower = loa_lower,
      loa_upper = loa_upper,
      loa_lower_lwr = loa_lower_lwr,
      loa_lower_upr = loa_lower_upr,
      loa_upper_lwr = loa_upper_lwr,
      loa_upper_upr = loa_upper_upr,
      ar1_rho = ar1_rho,
      ar1_estimated = ar1_estimated
    ))
  }

  expect_close_to_reference <- function(observed, expected, tol = 5e-6) {
    for (nm in names(expected)) {
      if (is.logical(expected[[nm]])) {
        expect_identical(observed[[nm]], expected[[nm]], info = nm)
      } else if (is.na(expected[[nm]])) {
        expect_true(is.na(observed[[nm]]), info = nm)
      } else {
        expect_equal(observed[[nm]], expected[[nm]], tolerance = tol, info = nm)
      }
    }
  }

  expected <- list(
    iid_boundary = list(
      bias_mu0 = 1.10880351648566,
      bias_lwr = 0.881755926331974,
      bias_upr = 1.33585110663935,
      sigma2_subject = 0,
      sigma2_resid = 2.65706881134442,
      sd_loa = 1.63005178179849,
      loa_lower = -2.08609797583937,
      loa_upper = 4.30370500881069,
      loa_lower_lwr = -2.57353079579423,
      loa_lower_upr = -1.5986651558845,
      loa_upper_lwr = 3.81627218885742,
      loa_upper_upr = 4.79113782876397,
      ar1_rho = NA_real_,
      ar1_estimated = FALSE
    ),
    ar1_fixed_slope = list(
      bias_mu0 = 0.101168854430398,
      bias_lwr = -0.232691880183468,
      bias_upr = 0.435029589044263,
      sigma2_subject = 0.603133877867929,
      sigma2_resid = 0.604988223925724,
      sd_loa = 1.09914607845984,
      loa_lower = -2.09712330248929,
      loa_upper = 2.29946101135009,
      loa_lower_lwr = -2.61556143639825,
      loa_lower_upr = -1.57868516858033,
      loa_upper_lwr = 1.78104737272341,
      loa_upper_upr = 2.81787464997676,
      ar1_rho = 0.35,
      ar1_estimated = FALSE
    ),
    ar1_profiled = list(
      bias_mu0 = 0.35051291798968,
      bias_lwr = 0.024251325329714,
      bias_upr = 0.676774510649645,
      sigma2_subject = 0.247592888510583,
      sigma2_resid = 1.13252478725102,
      sd_loa = 1.17478409750967,
      loa_lower = -1.95206391312926,
      loa_upper = 2.65308974910862,
      loa_lower_lwr = -2.45695199981367,
      loa_lower_upr = -1.44717582644486,
      loa_upper_lwr = 2.14433296879168,
      loa_upper_upr = 3.16184652942557,
      ar1_rho = 0.555555555555556,
      ar1_estimated = TRUE
    ),
    exact_boundary = list(
      bias_mu0 = 0.799999999999998,
      bias_lwr = 0.799999561738728,
      bias_upr = 0.800000438261268,
      sigma2_subject = 0,
      sigma2_resid = 1e-12,
      sd_loa = 1e-06,
      loa_lower = 0.799998039999998,
      loa_upper = 0.800001959999998,
      loa_lower_lwr = 0.799997601738728,
      loa_lower_upr = 0.799998478261268,
      loa_upper_lwr = 0.800001521738728,
      loa_upper_upr = 0.800002398261268,
      ar1_rho = NA_real_,
      ar1_estimated = FALSE
    )
  )

  dat_iid <- sim_two_method_known(S = 22L, Tm = 9L, mu = 1.1, sig_s = 0.9, sig_e = 1.3, seed = 42L)
  fit_iid <- matrixCorr:::bland_altman_repeated_em_ext_cpp(
    y = dat_iid$y,
    subject = dat_iid$subject,
    method = as.integer(dat_iid$method == "M2") + 1L,
    time = dat_iid$time,
    include_slope = FALSE,
    use_ar1 = FALSE,
    conf_level = 0.95,
    loa_multiplier_arg = 1.96
  )
  expect_close_to_reference(extract_fit(fit_iid), expected$iid_boundary, tol = 5e-6)

  dat_ar1_fixed <- sim_two_method_known(S = 18L, Tm = 8L, mu = 0.6, sig_s = 1.0, sig_e = 0.8, rho = 0.35, seed = 314L)
  fit_ar1_fixed <- matrixCorr:::bland_altman_repeated_em_ext_cpp(
    y = dat_ar1_fixed$y,
    subject = dat_ar1_fixed$subject,
    method = as.integer(dat_ar1_fixed$method == "M2") + 1L,
    time = dat_ar1_fixed$time,
    include_slope = TRUE,
    use_ar1 = TRUE,
    ar1_rho = 0.35,
    conf_level = 0.90,
    loa_multiplier_arg = 2.0
  )
  expect_close_to_reference(extract_fit(fit_ar1_fixed), expected$ar1_fixed_slope, tol = 5e-6)

  dat_ar1_profiled <- sim_two_method_known(S = 24L, Tm = 7L, mu = 0.25, sig_s = 0.7, sig_e = 1.1, rho = 0.55, seed = 20260401L)
  fit_ar1_profiled <- matrixCorr:::bland_altman_repeated_em_ext_cpp(
    y = dat_ar1_profiled$y,
    subject = dat_ar1_profiled$subject,
    method = as.integer(dat_ar1_profiled$method == "M2") + 1L,
    time = dat_ar1_profiled$time,
    include_slope = FALSE,
    use_ar1 = TRUE,
    ar1_rho = NA_real_,
    conf_level = 0.95,
    loa_multiplier_arg = 1.96
  )
  expect_close_to_reference(extract_fit(fit_ar1_profiled), expected$ar1_profiled, tol = 5e-6)

  S <- 10L
  Tm <- 2L
  subject_means <- seq(12, 30, length.out = S)
  means <- rep(subject_means, each = Tm)
  diffs <- rep(0.8 + 0.25 * (subject_means - mean(subject_means)), each = Tm)
  dat_boundary <- sim_two_method_from_pairs(means, diffs, S = S, Tm = Tm)
  fit_boundary <- matrixCorr:::bland_altman_repeated_em_ext_cpp(
    y = dat_boundary$y,
    subject = dat_boundary$subject,
    method = as.integer(dat_boundary$method == "B") + 1L,
    time = dat_boundary$time,
    include_slope = TRUE,
    use_ar1 = FALSE,
    conf_level = 0.95,
    loa_multiplier_arg = 1.96
  )
  expect_close_to_reference(extract_fit(fit_boundary), expected$exact_boundary, tol = 1e-8)
  expect_equal(as.numeric(fit_boundary$sigma2_subject), 0, tolerance = 1e-12)
})
