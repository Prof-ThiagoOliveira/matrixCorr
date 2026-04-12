test_that("ccc_rm_ustat: basic structure, symmetry, CI container", {
  set.seed(123)
  n_subj <- 600L
  id     <- factor(rep(seq_len(n_subj), each = 2L))
  method <- factor(rep(c("A","B"), times = n_subj))

  # model: y_A = u + e ; y_B = b + u + e
  sigA <- 1.0   # subject variance
  sigE <- 0.5   # error variance
  biasB <- 0.2  # fixed method bias for B
  u <- rnorm(n_subj, 0, sqrt(sigA))[as.integer(id)]
  e <- rnorm(n_subj * 2L, 0, sqrt(sigE))
  y <- (method == "B") * biasB + u + e
  df <- data.frame(y, id, method)

  # Theoretical CCC for 2 methods, T = 1, no AxM/AxT:
  # N = sA; D = sA + S_B + sE
  ccc_theory <- sigA / (sigA + biasB^2 + sigE)

  # estimates only
  c1 <- ccc_rm_ustat(df, response = "y", subject = "id", method = "method")
  expect_s3_class(c1, "ccc")
  expect_true(is.matrix(c1) && all(rownames(c1) == c("A","B")))
  expect_equal(as.numeric(diag(c1)), c(1,1))
  expect_equal(c1["A","B"], c1["B","A"])
  expect_true(c1["A","B"] > 0 && c1["A","B"] < 1)

  # within Monte Carlo tolerance
  expect_lt(abs(c1["A","B"] - ccc_theory), 0.05)

  # with CI container
  c2 <- ccc_rm_ustat(df, response = "y", subject = "id", method = "method", ci = TRUE, conf_level = 0.95)
  expect_s3_class(c2, "ccc_ci")
  expect_named(c2, c("est","lwr.ci","upr.ci"))
  expect_equal(dim(c2$est), c(2,2))
  expect_true(is.na(diag(c2$lwr.ci))[1] && is.na(diag(c2$upr.ci))[1])
  est <- c2$est["A","B"]; lwr <- c2$lwr.ci["A","B"]; upr <- c2$upr.ci["A","B"]
  expect_true(lwr <= est && est <= upr)
})

test_that("ccc_rm_ustat uses the aligned repeated-measures argument order", {
  set.seed(124)
  n_subj <- 80L
  n_time <- 3L
  id <- factor(rep(seq_len(n_subj), each = 2L * n_time))
  method <- factor(rep(rep(c("A", "B"), each = n_time), times = n_subj))
  time <- factor(rep(rep(seq_len(n_time), times = 2L), times = n_subj))
  u <- rnorm(n_subj, 0, 1)[as.integer(id)]
  y <- u + 0.2 * (method == "B") + rnorm(length(id), 0, 0.4)
  df <- data.frame(y, id, method, time)

  fit_named <- ccc_rm_ustat(
    df,
    response = "y",
    subject = "id",
    method = "method",
    time = "time"
  )
  fit_positional <- ccc_rm_ustat(
    df,
    "y",
    "id",
    "method",
    "time"
  )

  expect_equal(unname(fit_named["A", "B"]), unname(fit_positional["A", "B"]))
})

test_that("ccc_rm_reml (pairwise, no time): matches simple theory and returns VCs", {
  set.seed(123)
  n_subj <- 500L
  id     <- factor(rep(seq_len(n_subj), each = 2L))
  method <- factor(rep(c("A","B"), times = n_subj))

  sigA <- 1.0
  sigE <- 0.5
  biasB <- 0.2
  u <- rnorm(n_subj, 0, sqrt(sigA))[as.integer(id)]
  e <- rnorm(n_subj * 2L, 0, sqrt(sigE))
  y <- (method == "B") * biasB + u + e
  df <- data.frame(y, id, method)

  cfit <- ccc_rm_reml(df, response = "y", subject = "id",
                       method = "method", ci = TRUE)
  expect_s3_class(cfit, "ccc_rm_reml")
  expect_named(cfit, c("est","lwr.ci","upr.ci"))
  expect_equal(rownames(cfit$est), c("A","B"))
  expect_equal(colnames(cfit$est), c("A","B"))
  expect_equal(as.numeric(diag(cfit$est)), c(1,1))

  # theory for T=1 (no AxM/AxT): sA / (sA + b^2 + sE)
  ccc_theory <- sigA / (sigA + biasB^2 + sigE)
  est <- cfit$est["A","B"]
  expect_lt(abs(est - ccc_theory), 0.05)

  # variance-component attributes exist and are numeric matrices
  for (nm in c("sigma2_subject","sigma2_subject_method","sigma2_subject_time",
               "sigma2_error","SB","se_ccc")) {
    v <- attr(cfit, nm)
    expect_true(is.matrix(v))
    expect_equal(dim(v), c(2,2))
  }

  out_print <- capture.output(print(cfit))
  expect_false(any(grepl("sigma2_subject|sigma2_subject_method|sigma2_subject_time|sigma2_error|SB|se_ccc|Variance components|AR\\(1\\) diagnostics",
                         out_print)))

  # summary data frame columns present
  sm <- summary(cfit, show_ci = "yes", digits = 4)
  expect_s3_class(sm, "summary.ccc_rm_reml")
  expect_true(all(c("item1","item2","estimate","lwr","upr","n_subjects","n_obs",
                    "sigma2_subject","sigma2_subject_method","sigma2_subject_time",
                    "sigma2_error","SB","se_ccc") %in% names(sm)))
  out <- capture.output(print(sm))
  expect_true(any(grepl("^Concordance estimates$", out)))
  expect_true(any(grepl("^Variance components$", out)))
})

# helper to center time within subject
center_by_id <- function(x, id) ave(x, id, FUN = function(v) v - mean(v))

sim_ccc_rm_dat <- function(seed, rho = 0, n_subj = 120L, n_time = 6L,
                           sig_subject = 1, sig_error = 0.7, bias_b = 0.2) {
  set.seed(seed)
  id <- factor(rep(seq_len(n_subj), each = 2L * n_time))
  time <- factor(rep(rep(seq_len(n_time), times = 2L), times = n_subj))
  method <- factor(rep(rep(c("A", "B"), each = n_time), times = n_subj))
  subj_eff <- rnorm(n_subj, 0, sig_subject)
  y <- numeric(length(id))

  for (s in seq_len(n_subj)) {
    for (m in c("A", "B")) {
      idx <- which(id == levels(id)[s] & method == m)
      eps <- if (abs(rho) < 1e-12) {
        rnorm(n_time, 0, sig_error)
      } else {
        as.numeric(stats::arima.sim(list(ar = rho), n = n_time, sd = sig_error))
      }
      y[idx] <- subj_eff[s] + (m == "B") * bias_b + eps
    }
  }

  data.frame(y = y, id = id, method = method, time = time)
}

sim_ccc_rm_irregular <- function(seed, rho = 0, n_subj = 30L, n_time = 6L, n_method = 3L,
                                 miss = 0.35, sig_subject = 0.5, sig_error = 0.8,
                                 bias_scale = 0.5, time_effect = 1.0,
                                 sig_subject_time = 0.7) {
  set.seed(seed)
  ids <- seq_len(n_subj)
  methods <- LETTERS[seq_len(n_method)]
  rows <- vector("list", n_subj * n_method)
  kk <- 0L
  subj_eff <- rnorm(n_subj, 0, sig_subject)
  subj_time_eff <- matrix(rnorm(n_subj * n_time, 0, sig_subject_time), n_subj, n_time)

  for (s in ids) {
    for (m in seq_along(methods)) {
      keep_t <- sort(sample(seq_len(n_time), size = max(2L, ceiling((1 - miss) * n_time))))
      eps <- if (abs(rho) < 1e-12) {
        rnorm(length(keep_t), 0, sig_error)
      } else {
        as.numeric(stats::arima.sim(list(ar = rho), n = length(keep_t), sd = sig_error))
      }
      y <- subj_eff[s] + subj_time_eff[s, keep_t] +
        time_effect * (keep_t - mean(seq_len(n_time))) / n_time +
        (m - 1) * bias_scale + eps
      kk <- kk + 1L
      rows[[kk]] <- data.frame(
        y = y,
        id = factor(s),
        method = factor(methods[m], levels = methods),
        time = factor(keep_t, levels = seq_len(n_time))
      )
    }
  }

  do.call(rbind, rows[seq_len(kk)])
}

ccc_rm_raw_ar1_recommend <- function(row) {
  scalar_or_na <- function(x) {
    if (is.null(x) || !length(x)) return(NA_real_)
    suppressWarnings(as.numeric(x[[1L]]))
  }

  sag_hat <- scalar_or_na(row[["sigma2_subject_time"]])
  se_hat  <- scalar_or_na(row[["sigma2_error"]])
  rho_hat <- scalar_or_na(row[["ar1_rho_lag1"]])
  pval    <- scalar_or_na(row[["ar1_pval"]])

  sag_share <- if (is.finite(sag_hat) && is.finite(se_hat)) {
    sag_hat / max(1e-12, sag_hat + se_hat)
  } else {
    0
  }
  thr   <- if (sag_share > 0.25) 0.20 else 0.10
  p_thr <- if (sag_share > 0.25) 0.01 else 0.05

  is.finite(rho_hat) && is.finite(pval) && rho_hat >= thr && pval < p_thr
}

find_ccc_rm_false_positive_candidate <- function() {
  param_grid <- list(
    list(n_subj = 30L, n_time = 6L, n_method = 3L, miss = 0.35, sig_subject = 0.5, sig_error = 0.8, bias_scale = 0.5, time_effect = 1.0, sig_subject_time = 0.7),
    list(n_subj = 24L, n_time = 6L, n_method = 3L, miss = 0.25, sig_subject = 0.4, sig_error = 0.9, bias_scale = 0.6, time_effect = 0.9, sig_subject_time = 0.8),
    list(n_subj = 36L, n_time = 5L, n_method = 3L, miss = 0.30, sig_subject = 0.6, sig_error = 0.7, bias_scale = 0.5, time_effect = 1.1, sig_subject_time = 0.6)
  )

  for (params in param_grid) {
    for (seed in 1:80) {
      dat <- do.call(sim_ccc_rm_irregular, c(list(seed = seed, rho = 0), params))
      fit <- suppressMessages(
        ccc_rm_reml(dat, "y", "id", method = "method", time = "time", ar = "none")
      )
      sm <- as.data.frame(summary(fit))
      raw_reco <- vapply(seq_len(nrow(sm)), function(i) ccc_rm_raw_ar1_recommend(sm[i, , drop = FALSE]), logical(1))
      final_reco <- sm$use_ar1 %in% TRUE
      hit <- which(raw_reco & !final_reco)
      if (length(hit)) {
        return(list(
          seed = seed,
          params = params,
          summary = sm[hit[1], , drop = FALSE]
        ))
      }
    }
  }

  NULL
}

test_that("Dmat_type affects CCC as expected when biases flip over time", {
  set.seed(42)
  n_subj <- 200L
  n_time <- 6L
  id     <- factor(rep(seq_len(n_subj), each = 2L * n_time))
  time   <- factor(rep(rep(seq_len(n_time), times = 2L), times = n_subj))
  method <- factor(rep(rep(c("A","B"), each = n_time), times = n_subj))

  sigA  <- 0.8
  sigAT <- 0.3
  sigE  <- 0.4

  # time-varying bias for B: +b, -b, +b, -b, ...
  b0 <- 0.35
  tnum <- as.integer(time)
  bias_t <- ifelse(tnum %% 2L == 1L,  b0, -b0)
  bias  <- ifelse(method == "B", bias_t, 0)

  u  <- rnorm(n_subj, 0, sqrt(sigA))[as.integer(id)]
  gIT <- rnorm(n_subj * n_time, 0, sqrt(sigAT))
  g  <- gIT[(as.integer(id) - 1L) * n_time + as.integer(time)]
  y  <- bias + u + g + rnorm(length(id), 0, sqrt(sigE))
  df <- data.frame(y, id, method, time)

  fit_avg <- ccc_rm_reml(df, "y", "id", method = "method", time = "time",
                          Dmat_type = "time-avg")
  fit_typ <- ccc_rm_reml(df, "y", "id", method = "method", time = "time",
                          Dmat_type = "typical-visit")

  # With alternating biases, squared-average is ~0, average of squares > 0,
  # so CCC(time-avg) should be >= CCC(typical-visit)
  c_avg <- fit_avg["A","B"]; c_typ <- fit_typ["A","B"]
  expect_gte(c_avg, c_typ)
})

test_that("Supplying a wrong-sized Dmat errors clearly", {
  set.seed(1)
  n_subj <- 10L; n_time <- 3L
  id     <- factor(rep(seq_len(n_subj), each = 2L * n_time))
  time   <- factor(rep(rep(seq_len(n_time), times = 2L), times = n_subj))
  method <- factor(rep(rep(c("A","B"), each = n_time), times = n_subj))
  y <- rnorm(length(id))
  df <- data.frame(y, id, method, time)

  badD <- diag(5L)   # wrong dimension
  expect_error(
    ccc_rm_reml(df, "y", "id", method = "method", time = "time", Dmat = badD),
    "Dmat has incompatible dimension", fixed = FALSE
  )
})

test_that("summary adds sigma2_extra* columns when slope is enabled", {
  set.seed(2)
  n_subj <- 120L; n_time <- 4L
  id     <- factor(rep(seq_len(n_subj), each = 2L * n_time))
  time   <- factor(rep(rep(seq_len(n_time), times = 2L), times = n_subj))
  method <- factor(rep(rep(c("A","B"), each = n_time), times = n_subj))

  # Build a genuine slope effect so sigma2_extra > 0
  tnum <- as.integer(time)
  t_c  <- center_by_id(tnum, id)
  slope_subj <- rnorm(n_subj, 0, 0.15)[as.integer(id)]

  y <- 0.3 * (method == "B") +
    slope_subj * t_c +
    rnorm(length(id), 0, 0.6)

  df <- data.frame(y, id, method, time, t_c = t_c)

  fit_slope <- ccc_rm_reml(df, "y", "id", method = "method", time = "time",
                            slope = "subject", slope_var = "t_c", ci = TRUE)

  # Summary should expose sigma2_extra columns
  sm <- summary(fit_slope, show_ci = "yes")
  extra_cols <- grep("^sigma2_extra", names(sm), value = TRUE)
  expect_true(length(extra_cols) >= 1)
  expect_true(any(is.finite(sm[[extra_cols[1]]])))

  # In the no-slope case, those columns should not be present
  fit_noslope <- ccc_rm_reml(df, "y", "id", method = "method", time = "time", ci = TRUE)
  sm0 <- summary(fit_noslope, show_ci = "yes")
  expect_length(grep("^sigma2_extra", names(sm0), value = TRUE), 0)
})

test_that("AR(1) path: fixed rho is carried in attributes", {
  set.seed(10)
  n_subj <- 50L; n_time <- 6L
  id     <- factor(rep(seq_len(n_subj), each = 2L * n_time))
  time   <- factor(rep(rep(seq_len(n_time), times = 2L), times = n_subj))
  method <- factor(rep(rep(c("A","B"), each = n_time), times = n_subj))

  rho_true <- 0.6
  # simulate AR(1) per (subject,method) run
  y <- numeric(length(id))
  for (s in seq_len(n_subj)) {
    for (m in c("A","B")) {
      idx <- which(id == levels(id)[s] & method == m)
      y[idx] <- 0.2 * (m == "B") +
        as.numeric(stats::arima.sim(list(ar = rho_true), n = n_time, sd = 0.7))
    }
  }
  df <- data.frame(y, id, method, time)

  fit_ar <- ccc_rm_reml(df, "y", "id", method = "method", time = "time",
                         ar = "ar1", ar_rho = rho_true)
  rr <- attr(fit_ar, "ar_rho")
  expect_true(is.matrix(rr))
  expect_equal(rr["A","B"], rho_true, tolerance = 1e-12)

  # summary returns AR1 diagnostics columns (may have NAs but should exist)
  sm <- summary(fit_ar)
  expect_true("ar1_rho" %in% names(sm) || "ar1_rho_lag1" %in% names(sm))
  out <- capture.output(print(sm))
  expect_true(any(grepl("^AR\\(1\\) diagnostics$", out)))
})

test_that("AR(1) recommendation distinguishes IID from positive serial correlation", {
  dat_iid <- sim_ccc_rm_dat(seed = 7, rho = 0)
  fit_iid <- expect_no_message(
    ccc_rm_reml(dat_iid, "y", "id", method = "method", time = "time", ar = "none")
  )
  sm_iid <- as.data.frame(summary(fit_iid))
  expect_false(isTRUE(sm_iid$use_ar1[1]))
  expect_lt(sm_iid$ar1_rho_lag1[1], 0)

  dat_ar1 <- sim_ccc_rm_dat(seed = 1, rho = 0.6)
  expect_message(
    fit_ar1_diag <- ccc_rm_reml(dat_ar1, "y", "id", method = "method", time = "time", ar = "none"),
    "Positive lag-1 residual correlation detected"
  )
  sm_ar1 <- as.data.frame(summary(fit_ar1_diag))
  expect_true(isTRUE(sm_ar1$use_ar1[1]))
  expect_gt(sm_ar1$ar1_rho_lag1[1], 0.1)
})

test_that("AR(1) fallback message aligns with ba_rm wording", {
  dat <- sim_ccc_rm_dat(seed = 11, rho = 0)
  subj_time <- rep(c(1L, 2L), length.out = nlevels(dat$id))
  dat$time <- factor(subj_time[as.integer(dat$id)], levels = 1:2)

  expect_message(
    ccc_rm_reml(dat, "y", "id", method = "method", time = "time", ar = "ar1", verbose = TRUE),
    "Requested AR\\(1\\) residual structure could not be fit for pair\\(s\\): A vs B; using iid residuals instead\\."
  )
})

test_that("AR(1) recommendation is not triggered by an IID irregular-panel false positive", {
  cand <- find_ccc_rm_false_positive_candidate()
  skip_if(is.null(cand), "No IID irregular-panel raw false-positive candidate found on this platform/configuration.")

  row <- cand$summary[1, , drop = FALSE]
  expect_true(ccc_rm_raw_ar1_recommend(row))
  expect_false(isTRUE(row$use_ar1))
})


ccc_lin <- function(x, y, na.rm = TRUE) {
  stopifnot(length(x) == length(y))
  if (na.rm) {
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]; y <- y[keep]
  }
  if (length(x) < 2L) return(NA_real_)
  mx <- mean(x); my <- mean(y)
  vx <- stats::var(x); vy <- stats::var(y)
  sxy <- stats::cov(x, y)
  den <- vx + vy + (mx - my)^2
  if (!is.finite(den) || den <= 0) return(NA_real_)
  2 * sxy / den
}

test_that("ccc_rm_ustat reduces to Lin's CCC when T = 1", {
  set.seed(123)
  n <- 1500L
  u  <- rnorm(n, 0, 1.0)
  eA <- rnorm(n, 0, 0.7)
  eB <- rnorm(n, 0, 0.7)
  bias <- 0.3
  xA <- u + eA
  xB <- bias + u + eB

  # Lin on paired data
  c_lin <- ccc_lin(xA, xB)

  df <- data.frame(
    id     = factor(rep(seq_len(n), each = 2L)),
    method = factor(rep(c("A","B"), times = n)),
    time   = factor(rep(1L, 2L * n))
  )

  # INTERLEAVE A,B per subject (A1,B1,A2,B2,...)
  df$y <- c(rbind(xA, xB))

  c_us <- ccc_rm_ustat(df, response = "y", subject = "id", method = "method")
  c_us_AB <- unname(c_us["A","B"])

  # should match Lin very closely
  expect_equal(c_us_AB, c_lin, tolerance = 5e-3)
})

test_that("ccc_rm_reml (no time) approximates Lin's CCC when T = 1", {
  set.seed(124)
  n <- 1500L
  u  <- rnorm(n, 0, 1.0)
  eA <- rnorm(n, 0, 0.7)
  eB <- rnorm(n, 0, 0.7)
  bias <- 0.25
  xA <- u + eA
  xB <- bias + u + eB
  c_lin <- ccc_lin(xA, xB)

  df <- data.frame(
    id     = factor(rep(seq_len(n), each = 2L)),
    method = factor(rep(c("A","B"), times = n)),
    y      = c(rbind(xA, xB))
  )

  fit <- ccc_rm_reml(df, response = "y", subject = "id", method = "method")
  c_reml <- unname(fit["A","B"])

  expect_equal(c_reml, c_lin, tolerance = 1e-2)
})

test_that("repeated-measures CCC uses the subject argument name", {
  set.seed(432)
  n_subj <- 40L
  n_time <- 3L
  id <- factor(rep(seq_len(n_subj), each = 2L * n_time))
  method <- factor(rep(rep(c("A", "B"), each = n_time), times = n_subj))
  time <- factor(rep(rep(seq_len(n_time), times = 2L), times = n_subj))
  u <- rnorm(n_subj, 0, 1)[as.integer(id)]
  y <- u + 0.15 * (method == "B") + rnorm(length(id), 0, 0.3)
  df <- data.frame(y, id, method, time)

  fit_named <- ccc_rm_reml(
    df,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    ci = FALSE
  )
  fit_positional <- ccc_rm_reml(
    df,
    "y",
    "id",
    method = "method",
    time = "time",
    ci = FALSE
  )

  expect_equal(unname(fit_named["A", "B"]), unname(fit_positional["A", "B"]))
})

test_that("repeated-measures CCC honors n_threads without changing estimates", {
  set.seed(433)
  df <- sim_ccc_rm_dat(seed = 9, rho = 0.2, n_subj = 50L, n_time = 4L)

  fit1 <- ccc_rm_reml(
    df,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    n_threads = 1L,
    ci = FALSE
  )
  fit2 <- ccc_rm_reml(
    df,
    response = "y",
    subject = "id",
    method = "method",
    time = "time",
    n_threads = 2L,
    ci = FALSE
  )

  expect_equal(unname(fit1["A", "B"]), unname(fit2["A", "B"]), tolerance = 1e-12)
})
