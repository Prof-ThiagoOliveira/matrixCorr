make_cia_rm_dat2 <- function() {
  dat <- expand.grid(
    id = factor(c("s1", "s2", "s3")),
    method = factor(c("A", "B")),
    time = factor(c("t1", "t2", "t3")),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
    )
  dat$y <- c(
    1.0, 2.0, 1.5,
    2.5, 2.8, 3.2,
    3.0, 4.5, 4.2,
    1.7, 2.8, 2.0,
    2.9, 3.4, 4.1,
    4.0, 5.8, 4.9
  )
  dat
}

make_cia_rm_dat3 <- function() {
  dat <- expand.grid(
    id = factor(c("s1", "s2", "s3")),
    method = factor(c("A", "B", "C")),
    time = factor(c("t1", "t2", "t3")),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
    )
  dat$y <- c(
    1.0, 2.0, 1.5,
    2.5, 2.8, 3.2,
    3.0, 4.5, 4.2,
    1.7, 2.8, 2.0,
    2.9, 3.4, 4.1,
    4.0, 5.8, 4.9,
    0.9, 1.6, 1.2,
    2.2, 2.7, 3.0,
    2.6, 4.0, 4.1
  )
  dat
}

make_cia_rm_dat_near1 <- function() {
  set.seed(42)
  dat <- expand.grid(
    id = factor(sprintf("s%02d", 1:8)),
    method = factor(c("Sensor_A", "Sensor_B")),
    time = 1:5,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  subj_eff <- rnorm(8, sd = 0.8)[dat$id]
  time_eff <- c(0, 0.1, 0.2, 0.3, 0.4)[dat$time]
  method_eff <- c(Sensor_A = 0, Sensor_B = 0.03)[dat$method]
  interaction <- ifelse(dat$method == "Sensor_B" & dat$time >= 4, 0.02, 0)
  dat$y <- subj_eff + time_eff + method_eff + interaction + rnorm(nrow(dat), sd = 0.18)
  dat
}

make_cia_rm_dat3_near1 <- function() {
  set.seed(122)
  dat <- expand.grid(
    id = factor(sprintf("s%02d", 1:8)),
    method = factor(c("Sensor_A", "Sensor_B", "Sensor_C")),
    time = 1:5,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  subj_eff <- rnorm(8, sd = 0.75)[dat$id]
  time_eff <- c(0, 0.08, 0.17, 0.26, 0.35)[dat$time]
  method_eff <- c(Sensor_A = 0, Sensor_B = 0.03, Sensor_C = -0.02)[dat$method]
  interaction <- ifelse(dat$method == "Sensor_B" & dat$time >= 4, 0.02, 0) +
    ifelse(dat$method == "Sensor_C" & dat$time <= 2, -0.01, 0)
  dat$y <- subj_eff + time_eff + method_eff + interaction + rnorm(nrow(dat), sd = 0.25)
  dat
}

manual_cia_rm_pair <- function(dat, method_a, method_b, constrain_vc = FALSE, homogeneous = FALSE) {
  dat <- dat[dat$method %in% c(method_a, method_b), , drop = FALSE]
  dat$method <- factor(dat$method, levels = c(method_a, method_b))
  dat$time <- factor(dat$time)
  dat$id <- factor(dat$id)

  ids <- levels(dat$id)
  times <- levels(dat$time)
  n <- length(ids)
  K <- length(times)

  Z <- array(NA_real_, dim = c(n, 2L, K), dimnames = list(ids, c(method_a, method_b), times))
  for (i in seq_along(ids)) {
    for (r in 1:2) {
      for (k in seq_along(times)) {
        Z[i, r, k] <- dat$y[dat$id == ids[[i]] & dat$method == levels(dat$method)[[r]] & dat$time == times[[k]]]
      }
    }
  }

  grand <- mean(Z)
  Zi <- apply(Z, 1, mean)
  Zr <- apply(Z, 2, mean)
  Zk <- apply(Z, 3, mean)
  Zir <- apply(Z, c(1, 2), mean)
  Zik <- apply(Z, c(1, 3), mean)
  Zrk <- apply(Z, c(2, 3), mean)

  ss_subject_method <- K * sum((Zir - Zi - rep(Zr, each = n) + grand)^2)
  ss_method_time <- n * sum((Zrk - matrix(Zr, nrow = 2L, ncol = K) - matrix(Zk, nrow = 2L, ncol = K, byrow = TRUE) + grand)^2)

  ss_error <- 0
  for (i in seq_len(n)) {
    for (r in 1:2) {
      for (k in seq_len(K)) {
        resid <- Z[i, r, k] - Zir[i, r] - Zik[i, k] - Zrk[r, k] + Zi[i] + Zr[r] + Zk[k] - grand
        ss_error <- ss_error + resid^2
      }
    }
  }

  df_subject_method <- n - 1
  df_method_time <- K - 1
  df_error <- (n - 1) * (K - 1)
  ms_subject_method <- ss_subject_method / df_subject_method
  ms_method_time <- ss_method_time / df_method_time
  ms_error <- ss_error / df_error

  sigma2_error <- ms_error
  sigma2_subject_method <- (ms_subject_method - ms_error) / K
  if (isTRUE(constrain_vc)) {
    sigma2_error <- max(sigma2_error, 0)
    sigma2_subject_method <- max(sigma2_subject_method, 0)
  }

  d_k <- vapply(seq_len(K), function(k) mean(Z[, 1, k] - Z[, 2, k]), numeric(1))
  denom <- d_k^2 + 2 * sigma2_subject_method + 2 * sigma2_error
  psi <- ifelse(is.finite(denom) & denom > 0, 2 * sigma2_error / denom, NA_real_)

  F_stat <- ms_method_time / ms_error
  p_val <- stats::pf(F_stat, df1 = df_method_time, df2 = df_error, lower.tail = FALSE)

  out <- list(
    ss_subject_method = unname(ss_subject_method),
    ss_method_time = unname(ss_method_time),
    ss_error = unname(ss_error),
    ms_subject_method = unname(ms_subject_method),
    ms_method_time = unname(ms_method_time),
    ms_error = unname(ms_error),
    sigma2_error = unname(sigma2_error),
    sigma2_subject_method = unname(sigma2_subject_method),
    repeatability = unname(1.96 * sqrt(2 * sigma2_error)),
    d_k = unname(d_k),
    psi = unname(psi),
    F = unname(F_stat),
    p = unname(p_val),
    df_method_time = unname(df_method_time),
    df_error = unname(df_error)
  )

  if (isTRUE(homogeneous)) {
    ms_error_reduced <- (ss_method_time + ss_error) / (df_method_time + df_error)
    sigma2_error_common <- ms_error_reduced
    sigma2_subject_method_common <- (ms_subject_method - ms_error_reduced) / K
    if (isTRUE(constrain_vc)) {
      sigma2_error_common <- max(sigma2_error_common, 0)
      sigma2_subject_method_common <- max(sigma2_subject_method_common, 0)
    }
    d_common <- mean(vapply(seq_len(K), function(k) Z[, 1, k] - Z[, 2, k], numeric(n)))
    denom_common <- d_common^2 + 2 * sigma2_subject_method_common + 2 * sigma2_error_common
    out$common <- unname(if (is.finite(denom_common) && denom_common > 0) {
      2 * sigma2_error_common / denom_common
    } else {
      NA_real_
    })
    out$common_sigma2_error <- unname(sigma2_error_common)
    out$common_sigma2_subject_method <- unname(sigma2_subject_method_common)
  }

  out
}

manual_cia_rm_overall <- function(dat, constrain_vc = FALSE, homogeneous = FALSE) {
  dat$method <- factor(dat$method)
  dat$time <- factor(dat$time)
  dat$id <- factor(dat$id)

  ids <- levels(dat$id)
  methods <- levels(dat$method)
  times <- levels(dat$time)
  n <- length(ids)
  J <- length(methods)
  K <- length(times)

  Z <- array(NA_real_, dim = c(n, J, K), dimnames = list(ids, methods, times))
  for (i in seq_along(ids)) {
    for (j in seq_along(methods)) {
      for (k in seq_along(times)) {
        Z[i, j, k] <- dat$y[
          dat$id == ids[[i]] &
            dat$method == methods[[j]] &
            dat$time == times[[k]]
        ]
      }
    }
  }

  grand <- mean(Z)
  Zi <- apply(Z, 1, mean)
  Zj <- apply(Z, 2, mean)
  Zk <- apply(Z, 3, mean)
  Zij <- apply(Z, c(1, 2), mean)
  Zik <- apply(Z, c(1, 3), mean)
  Zjk <- apply(Z, c(2, 3), mean)

  ss_subject_method <- K * sum(vapply(
    seq_len(n),
    function(i) {
      sum((Zij[i, ] - Zi[i] - Zj + grand)^2)
    },
    numeric(1)
  ))
  ss_method_time <- n * sum(vapply(
    seq_len(J),
    function(j) {
      sum((Zjk[j, ] - Zj[j] - Zk + grand)^2)
    },
    numeric(1)
  ))

  ss_error <- 0
  for (i in seq_len(n)) {
    for (j in seq_len(J)) {
      for (k in seq_len(K)) {
        resid <- Z[i, j, k] - Zij[i, j] - Zik[i, k] - Zjk[j, k] + Zi[i] + Zj[j] + Zk[k] - grand
        ss_error <- ss_error + resid^2
      }
    }
  }

  df_subject_method <- (n - 1) * (J - 1)
  df_method_time <- (J - 1) * (K - 1)
  df_error <- (n - 1) * (J - 1) * (K - 1)

  ms_subject_method <- ss_subject_method / df_subject_method
  ms_method_time <- ss_method_time / df_method_time
  ms_error <- ss_error / df_error

  sigma2_error <- ms_error
  sigma2_subject_method <- (ms_subject_method - ms_error) / K
  if (isTRUE(constrain_vc)) {
    sigma2_error <- max(sigma2_error, 0)
    sigma2_subject_method <- max(sigma2_subject_method, 0)
  }

  pair_idx <- utils::combn(seq_len(J), 2)
  mean_pair_d2 <- vapply(
    seq_len(K),
    function(k) {
      mean(vapply(
        seq_len(ncol(pair_idx)),
        function(p) {
          idx <- pair_idx[, p]
          (Zjk[idx[1], k] - Zjk[idx[2], k])^2
        },
        numeric(1)
      ))
    },
    numeric(1)
  )
  denom <- mean_pair_d2 + 2 * sigma2_subject_method + 2 * sigma2_error
  overall <- ifelse(is.finite(denom) & denom > 0, 2 * sigma2_error / denom, NA_real_)

  out <- list(
    overall = unname(overall),
    sigma2_error = unname(sigma2_error),
    sigma2_subject_method = unname(sigma2_subject_method),
    repeatability = unname(1.96 * sqrt(2 * sigma2_error)),
    ms_subject_method = unname(ms_subject_method),
    ms_method_time = unname(ms_method_time),
    ms_error = unname(ms_error),
    ss_subject_method = unname(ss_subject_method),
    ss_method_time = unname(ss_method_time),
    ss_error = unname(ss_error),
    df_subject_method = unname(df_subject_method),
    df_method_time = unname(df_method_time),
    df_error = unname(df_error),
    F = unname(ms_method_time / ms_error),
    p = unname(stats::pf(ms_method_time / ms_error, df1 = df_method_time, df2 = df_error, lower.tail = FALSE))
  )

  if (isTRUE(homogeneous)) {
    ms_error_reduced <- (ss_method_time + ss_error) / (df_method_time + df_error)
    common_sigma2_error <- ms_error_reduced
    common_sigma2_subject_method <- (ms_subject_method - ms_error_reduced) / K
    if (isTRUE(constrain_vc)) {
      common_sigma2_error <- max(common_sigma2_error, 0)
      common_sigma2_subject_method <- max(common_sigma2_subject_method, 0)
    }
    mean_pair_d2_common <- mean(vapply(
      seq_len(ncol(pair_idx)),
      function(p) {
        idx <- pair_idx[, p]
        (Zj[idx[1]] - Zj[idx[2]])^2
      },
      numeric(1)
    ))
    denom_common <- mean_pair_d2_common + 2 * common_sigma2_subject_method + 2 * common_sigma2_error
    out$overall_common <- unname(if (is.finite(denom_common) && denom_common > 0) {
      2 * common_sigma2_error / denom_common
    } else {
      NA_real_
    })
    out$common_sigma2_error <- unname(common_sigma2_error)
    out$common_sigma2_subject_method <- unname(common_sigma2_subject_method)
    out$common_repeatability <- unname(1.96 * sqrt(2 * common_sigma2_error))
  }

  out
}

test_that("cia_rm wrapper validates required repeated-measures structure", {
  dat <- make_cia_rm_dat2()

  expect_error(cia_rm(dat, "y", "id", method = NULL, time = "time"), "method")
  expect_error(cia_rm(dat, "y", "id", method = "method", time = NULL), "time")

  dat_one_method <- subset(dat, method == "A")
  expect_error(cia_rm(dat_one_method, "y", "id", method = "method", time = "time"), "two method levels")

  dat_one_time <- subset(dat, time == "t1")
  expect_error(cia_rm(dat_one_time, "y", "id", method = "method", time = "time"), "two time/condition levels")

  dat_one_subject <- subset(dat, id == "s1")
  expect_error(cia_rm(dat_one_subject, "y", "id", method = "method", time = "time"), "two subjects")

  dat_dup <- rbind(dat, dat[1, , drop = FALSE])
  expect_error(
    cia_rm(dat_dup, "y", "id", method = "method", time = "time"),
    "each subject-method-time cell must contain exactly one observation"
  )

  dat_miss <- dat[-1, , drop = FALSE]
  expect_error(
    cia_rm(dat_miss, "y", "id", method = "method", time = "time"),
    "each subject-method-time cell must contain exactly one observation"
  )

  dat_bad_y <- dat
  dat_bad_y$y[1] <- Inf
  expect_error(cia_rm(dat_bad_y, "y", "id", method = "method", time = "time"), "finite numeric values")
})

test_that("cia_rm_anova_cpp returns the hand-computed ANOVA quantities", {
  dat <- make_cia_rm_dat2()
  manual <- manual_cia_rm_pair(dat, "A", "B")
  cpp_fun <- get("cia_rm_anova_cpp", envir = asNamespace("matrixCorr"))

  raw <- cpp_fun(
    y = dat$y,
    subject = as.integer(factor(dat$id)),
    method = as.integer(factor(dat$method, levels = c("A", "B"))),
    time = as.integer(factor(dat$time, levels = c("t1", "t2", "t3"))),
    n_subjects = 3L,
    n_methods = 2L,
    n_times = 3L,
    homogeneous = TRUE,
    constrain_vc = FALSE,
    n_threads = 1L
  )

  expect_equal(raw$ss_subject_method[1, 2], manual$ss_subject_method, tolerance = 1e-12)
  expect_equal(raw$ss_method_time[1, 2], manual$ss_method_time, tolerance = 1e-12)
  expect_equal(raw$ss_error[1, 2], manual$ss_error, tolerance = 1e-12)
  expect_equal(raw$ms_subject_method[1, 2], manual$ms_subject_method, tolerance = 1e-12)
  expect_equal(raw$ms_method_time[1, 2], manual$ms_method_time, tolerance = 1e-12)
  expect_equal(raw$ms_error[1, 2], manual$ms_error, tolerance = 1e-12)
  expect_equal(raw$sigma2_error[1, 2], manual$sigma2_error, tolerance = 1e-12)
  expect_equal(raw$sigma2_subject_method[1, 2], manual$sigma2_subject_method, tolerance = 1e-12)
})

test_that("cia_rm condition-specific CIA matches the repeated-measures formula", {
  dat <- make_cia_rm_dat2()
  manual <- manual_cia_rm_pair(dat, "A", "B", constrain_vc = TRUE)

  fit <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time"
  )

  expect_s3_class(fit, "cia_rm")
  expect_equal(unname(fit$est["A", "B", ]), manual$psi, tolerance = 1e-12)
  expect_equal(unname(attr(fit, "method_time_diff")["A", "B", ]), manual$d_k, tolerance = 1e-12)
})

test_that("cia_rm homogeneity test matches the method-time interaction F test", {
  dat <- make_cia_rm_dat2()
  manual <- manual_cia_rm_pair(dat, "A", "B")

  fit <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time"
  )

  expect_equal(unname(attr(fit, "homogeneity_F")["A", "B"]), manual$F, tolerance = 1e-12)
  expect_equal(unname(attr(fit, "homogeneity_p")["A", "B"]), manual$p, tolerance = 1e-12)
  expect_equal(unname(attr(fit, "df_method_time")["A", "B"]), manual$df_method_time, tolerance = 1e-12)
  expect_equal(unname(attr(fit, "df_error")["A", "B"]), manual$df_error, tolerance = 1e-12)
})

test_that("cia_rm common estimate matches the reduced homogeneous model", {
  dat <- make_cia_rm_dat2()
  manual <- manual_cia_rm_pair(dat, "A", "B", constrain_vc = TRUE, homogeneous = TRUE)

  fit <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    homogeneous = TRUE
  )

  expect_equal(unname(fit$common["A", "B"]), manual$common, tolerance = 1e-12)
  expect_equal(unname(attr(fit, "common_sigma2_error")["A", "B"]), manual$common_sigma2_error, tolerance = 1e-12)
  expect_equal(unname(attr(fit, "common_sigma2_subject_method")["A", "B"]), manual$common_sigma2_subject_method, tolerance = 1e-12)
})

test_that("cia_rm computes each method pair independently when there are three methods", {
  dat <- make_cia_rm_dat3()
  fit <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    homogeneous = TRUE
  )

  for (pair in list(c("A", "B"), c("A", "C"), c("B", "C"))) {
    manual <- manual_cia_rm_pair(dat, pair[[1]], pair[[2]], constrain_vc = TRUE, homogeneous = TRUE)
    expect_equal(unname(fit$est[pair[[1]], pair[[2]], ]), manual$psi, tolerance = 1e-12)
    expect_equal(unname(fit$common[pair[[1]], pair[[2]]]), manual$common, tolerance = 1e-12)
  }
})

test_that("global overall J=2 estimate matches the existing pairwise estimate", {
  dat <- make_cia_rm_dat2()
  cpp_fun <- get("cia_rm_anova_cpp", envir = asNamespace("matrixCorr"))

  raw <- cpp_fun(
    y = dat$y,
    subject = as.integer(factor(dat$id)),
    method = as.integer(factor(dat$method, levels = c("A", "B"))),
    time = as.integer(factor(dat$time, levels = c("t1", "t2", "t3"))),
    n_subjects = 3L,
    n_methods = 2L,
    n_times = 3L,
    homogeneous = TRUE,
    constrain_vc = FALSE,
    n_threads = 1L
  )

  expect_equal(as.numeric(raw$overall), c(raw$est[1, 2, ]), tolerance = 1e-12)

  fit <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    homogeneous = TRUE
  )
  expect_null(fit$overall)
  expect_null(fit$overall.common)
})

test_that("cia_rm overall multi-method estimate is computed from the global ANOVA model", {
  dat <- make_cia_rm_dat3()
  manual <- manual_cia_rm_overall(dat, constrain_vc = TRUE, homogeneous = TRUE)
  fit <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    homogeneous = TRUE
  )

  pair_mean <- vapply(
    seq_len(length(levels(dat$time))),
    function(k) mean(c(
      fit$est["A", "B", k],
      fit$est["A", "C", k],
      fit$est["B", "C", k]
    )),
    numeric(1)
  )

  expect_equal(unname(fit$overall), manual$overall, tolerance = 1e-12)
  expect_gt(max(abs(unname(fit$overall) - pair_mean)), 1e-6)
  expect_equal(unname(attr(fit, "overall_sigma2_error")), manual$sigma2_error, tolerance = 1e-12)
  expect_equal(unname(attr(fit, "overall_sigma2_subject_method")), manual$sigma2_subject_method, tolerance = 1e-12)
  expect_equal(unname(attr(fit, "overall_repeatability")), manual$repeatability, tolerance = 1e-12)
  expect_equal(unname(attr(fit, "overall_homogeneity_F")), manual$F, tolerance = 1e-12)
  expect_equal(unname(attr(fit, "overall_homogeneity_p")), manual$p, tolerance = 1e-12)
})

test_that("cia_rm overall common estimate uses the reduced homogeneous model", {
  dat <- make_cia_rm_dat3()
  manual <- manual_cia_rm_overall(dat, constrain_vc = TRUE, homogeneous = TRUE)
  fit <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    homogeneous = TRUE
  )

  expect_equal(unname(fit$overall.common), manual$overall_common, tolerance = 1e-12)
  expect_equal(unname(attr(fit, "overall_common_sigma2_error")), manual$common_sigma2_error, tolerance = 1e-12)
  expect_equal(unname(attr(fit, "overall_common_sigma2_subject_method")), manual$common_sigma2_subject_method, tolerance = 1e-12)
  expect_equal(unname(attr(fit, "overall_common_repeatability")), manual$common_repeatability, tolerance = 1e-12)
})

test_that("cia_rm bootstrap intervals resample subjects and preserve array shapes", {
  dat <- make_cia_rm_dat3()

  fit <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    ci = TRUE,
    inference = "bootstrap",
    homogeneous = TRUE,
    B = 50L,
    seed = 1L
  )

  expect_s3_class(fit, "cia_rm")
  expect_identical(attr(fit, "ci.method", exact = TRUE), "subject_bootstrap_percentile")
  expect_identical(dim(fit$est), dim(fit$lwr.ci))
  expect_identical(dim(fit$est), dim(fit$upr.ci))
  expect_identical(dim(fit$common), dim(fit$common.lwr.ci))
  expect_identical(dim(fit$common), dim(fit$common.upr.ci))
})

test_that("cia_rm uses delta-method confidence intervals by default", {
  dat <- make_cia_rm_dat3()

  fit <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    ci = TRUE,
    homogeneous = TRUE
  )

  expect_identical(attr(fit, "ci.method", exact = TRUE), "delta_normal")
  expect_identical(attr(fit, "estimator", exact = TRUE), "vc_constrained")
  expect_true(!is.null(fit$se))
  expect_true(!is.null(fit$lwr.ci))
  expect_true(!is.null(fit$upr.ci))
  expect_identical(dim(fit$se), dim(fit$est))
  expect_identical(dim(fit$lwr.ci), dim(fit$est))
  expect_identical(dim(fit$upr.ci), dim(fit$est))
  expect_identical(dim(fit$common.se), dim(fit$common))
  expect_identical(dim(fit$common.lwr.ci), dim(fit$common))
  expect_identical(dim(fit$common.upr.ci), dim(fit$common))
})

test_that("cia_rm unconstrained and constrained estimators differ only by the boundary rule", {
  dat <- make_cia_rm_dat_near1()

  fit_raw <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    ci = TRUE,
    inference = "delta",
    estimator = "mom_unconstrained"
  )
  fit_bnd <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    ci = TRUE,
    inference = "delta",
    estimator = "vc_constrained"
  )

  expect_true(any(fit_raw$est > 1, na.rm = TRUE) || any(fit_raw$est < 0, na.rm = TRUE) ||
                any(fit_raw$upr.ci > 1, na.rm = TRUE) || any(fit_raw$lwr.ci < 0, na.rm = TRUE))
  expect_true(all(fit_bnd$est[is.finite(fit_bnd$est)] <= 1))
  expect_true(all(fit_bnd$est[is.finite(fit_bnd$est)] >= 0))
  expect_true(any(fit_raw$upr.ci > 1, na.rm = TRUE) || any(fit_raw$lwr.ci < 0, na.rm = TRUE))
  expect_true(all(fit_bnd$upr.ci[is.finite(fit_bnd$upr.ci)] <= 1))
  expect_true(all(fit_bnd$lwr.ci[is.finite(fit_bnd$lwr.ci)] >= 0))
})

test_that("cia_rm overall constrained estimator truncates negative subject-method variance and bounds CIs", {
  dat <- make_cia_rm_dat3_near1()
  manual_raw <- manual_cia_rm_overall(dat, constrain_vc = FALSE, homogeneous = TRUE)
  manual_bnd <- manual_cia_rm_overall(dat, constrain_vc = TRUE, homogeneous = TRUE)

  fit_raw <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    ci = TRUE,
    inference = "delta",
    estimator = "mom_unconstrained",
    homogeneous = TRUE
  )
  fit_bnd <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    ci = TRUE,
    inference = "delta",
    estimator = "vc_constrained",
    homogeneous = TRUE
  )

  expect_lt(manual_raw$sigma2_subject_method, 0)
  expect_equal(unname(attr(fit_bnd, "overall_sigma2_subject_method")), 0, tolerance = 1e-12)
  expect_equal(unname(attr(fit_bnd, "overall_sigma2_subject_method")), manual_bnd$sigma2_subject_method, tolerance = 1e-12)
  expect_equal(unname(fit_bnd$overall), manual_bnd$overall, tolerance = 1e-10)
  expect_true(any(fit_raw$overall.upr.ci > 1, na.rm = TRUE) || any(fit_raw$overall.lwr.ci < 0, na.rm = TRUE))
  expect_true(all(fit_bnd$overall.upr.ci[is.finite(fit_bnd$overall.upr.ci)] <= 1))
  expect_true(all(fit_bnd$overall.lwr.ci[is.finite(fit_bnd$overall.lwr.ci)] >= 0))
  expect_true(is.finite(fit_bnd$overall.common))
  expect_true(is.finite(fit_bnd$overall.common.lwr.ci))
  expect_true(is.finite(fit_bnd$overall.common.upr.ci))
  expect_true(fit_bnd$overall.common.lwr.ci >= 0)
  expect_true(fit_bnd$overall.common.upr.ci <= 1)
})

test_that("cia_rm delta gradient matches an independent finite-difference check", {
  dat <- make_cia_rm_dat2()
  prep <- matrixCorr:::.mc_cia_rm_prepare(dat, "y", "id", "method", "time")
  cube <- matrixCorr:::.mc_cia_rm_cube(prep)
  D <- matrixCorr:::.mc_cia_rm_pair_diff_matrix(cube, 1L, 2L)
  stats <- matrixCorr:::.mc_cia_rm_pair_moment_stats(D)
  m <- c(stats$d, as.vector(stats$M))

  grad_pkg <- matrixCorr:::.mc_cia_rm_num_grad(
    function(mm) matrixCorr:::.mc_cia_rm_pair_est_from_moments(mm, K = stats$K, n = stats$n)[[1L]],
    m
  )

  grad_ref <- numeric(length(m))
  for (idx in seq_along(m)) {
    h <- 1e-7 * max(1, abs(m[[idx]]))
    plus <- m
    minus <- m
    plus[[idx]] <- plus[[idx]] + h
    minus[[idx]] <- minus[[idx]] - h
    grad_ref[[idx]] <-
      (matrixCorr:::.mc_cia_rm_pair_est_from_moments(plus, K = stats$K, n = stats$n)[[1L]] -
         matrixCorr:::.mc_cia_rm_pair_est_from_moments(minus, K = stats$K, n = stats$n)[[1L]]) / (2 * h)
  }

  expect_equal(grad_pkg, grad_ref, tolerance = 1e-6)
})

test_that("cia_rm returns common delta intervals when homogeneous is requested", {
  dat <- make_cia_rm_dat2()

  fit <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    ci = TRUE,
    homogeneous = TRUE
  )

  expect_true(is.matrix(fit$common.se))
  expect_true(is.matrix(fit$common.lwr.ci))
  expect_true(is.matrix(fit$common.upr.ci))
  expect_true(is.finite(fit$common.se["A", "B"]))
})

test_that("cia_rm returns overall multi-method confidence intervals when there are at least three methods", {
  dat <- make_cia_rm_dat3()

  fit_delta <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    ci = TRUE,
    homogeneous = TRUE
  )
  fit_boot <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    ci = TRUE,
    inference = "bootstrap",
    homogeneous = TRUE,
    B = 40L,
    seed = 1L
  )

  expect_length(fit_delta$overall, 3L)
  expect_length(fit_delta$overall.se, 3L)
  expect_length(fit_delta$overall.lwr.ci, 3L)
  expect_length(fit_delta$overall.upr.ci, 3L)
  expect_true(all(names(fit_delta$overall) == levels(factor(dat$time))))
  expect_true(is.finite(fit_delta$overall.common.se))
  expect_length(fit_boot$overall.lwr.ci, 3L)
  expect_length(fit_boot$overall.upr.ci, 3L)
  expect_true(is.finite(fit_boot$overall.common.lwr.ci))
  expect_true(is.finite(fit_boot$overall.common.upr.ci))
})

test_that("cia_rm objects carry repeated-measures diagnostics and summary structure", {
  dat <- make_cia_rm_dat3()
  fit <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    homogeneous = TRUE
  )

  expect_identical(attr(fit, "method"), "Repeated-measures coefficient of individual agreement")
  expect_identical(attr(fit, "description"), "Pairwise repeated-measures CIA from categorical-condition ANOVA")
  expect_identical(attr(fit, "package"), "matrixCorr")
  expect_true(is.matrix(attr(fit, "sigma2_error")))
  expect_true(is.matrix(attr(fit, "sigma2_subject_method")))
  expect_true(is.matrix(attr(fit, "repeatability")))
  expect_true(is.matrix(attr(fit, "homogeneity_F")))
  expect_true(is.matrix(attr(fit, "homogeneity_p")))
  expect_length(fit$overall, 3L)
  expect_true(is.numeric(attr(fit, "overall_sigma2_error")))
  expect_true(is.numeric(attr(fit, "overall_sigma2_subject_method")))
  expect_true(is.numeric(attr(fit, "overall_repeatability")))
  expect_true(is.numeric(attr(fit, "overall_homogeneity_F")))
  expect_true(is.numeric(attr(fit, "overall_homogeneity_p")))

  sm <- summary(fit)
  expect_s3_class(sm, "summary.cia_rm")
  expect_true(all(c(
    "method_1", "method_2", "time", "cia",
    "sigma2_error", "sigma2_subject_method",
    "repeatability", "homogeneity_F", "homogeneity_p"
  ) %in% names(sm)))
  expect_true(any(sm$time == "common"))
  expect_true(is.data.frame(attr(sm, "overall_section", exact = TRUE)))
  expect_true(is.data.frame(attr(sm, "overall_common_section", exact = TRUE)))
})

test_that("cia_rm print is compact while summary print includes diagnostics", {
  dat <- make_cia_rm_dat3()
  fit <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    homogeneous = TRUE
  )

  out_print <- paste(capture.output(print(fit)), collapse = "\n")
  out_summary <- paste(capture.output(print(summary(fit))), collapse = "\n")

  expect_match(out_print, "Condition-specific CIA estimates")
  expect_match(out_print, "Overall multi-method CIA")
  expect_match(out_print, "Common overall CIA")
  expect_no_match(out_print, "Pair diagnostics")
  expect_match(out_summary, "CIA estimates")
  expect_match(out_summary, "Overall multi-method CIA")
  expect_match(out_summary, "Common overall CIA")
  expect_match(out_summary, "Pair diagnostics")
})

test_that("cia_rm print preview respects digits", {
  dat <- make_cia_rm_dat2()
  fit <- cia_rm(
    dat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    ci = TRUE,
    inference = "bootstrap",
    B = 20L,
    seed = 1L
  )

  out_print <- paste(capture.output(print(fit, digits = 3)), collapse = "\n")
  expect_no_match(out_print, "0\\.[0-9]{10,}")
})

test_that("cia_rm plot detects categorical versus continuous condition scales", {
  dat_cat <- make_cia_rm_dat2()
  fit_cat <- cia_rm(
    dat_cat,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    homogeneous = TRUE
  )

  expect_identical(attr(fit_cat, "time_scale", exact = TRUE), "categorical")
  p_cat <- plot(fit_cat)
  expect_s3_class(p_cat, "ggplot")

  dat_num <- make_cia_rm_dat2()
  dat_num$time <- c(1, 2, 3)[match(dat_num$time, levels(dat_num$time))]
  fit_num <- cia_rm(
    dat_num,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    homogeneous = TRUE
  )

  expect_identical(attr(fit_num, "time_scale", exact = TRUE), "continuous")
  p_num <- plot(fit_num, facet_by_pair = TRUE)
  expect_s3_class(p_num, "ggplot")
  expect_true(inherits(p_num$facet, "FacetWrap"))
})

test_that("cia_rm plot handles confidence intervals for continuous condition scales", {
  dat_num <- make_cia_rm_dat3()
  dat_num$time <- c(1, 2, 3)[match(dat_num$time, levels(dat_num$time))]

  fit_num <- cia_rm(
    dat_num,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    ci = TRUE,
    inference = "bootstrap",
    homogeneous = TRUE,
    B = 20L,
    seed = 1L
  )

  p_num <- plot(fit_num, facet_by_pair = TRUE)
  expect_s3_class(p_num, "ggplot")
})
