dat_ref2 <- data.frame(
  subject = c(
    "s1", "s1", "s1", "s1",
    "s2", "s2", "s2", "s2"
  ),
  method = c(
    "R", "R", "A", "A",
    "R", "R", "A", "A"
  ),
  replicate = c(
    "r1", "r2", "r1", "r2",
    "r1", "r2", "r1", "r2"
  ),
  value = c(
    0, 2, 0, 0,
    10, 12, 10, 10
  )
)

dat_ref3 <- data.frame(
  subject = c(
    "s1", "s1", "s1", "s1", "s1", "s1",
    "s2", "s2", "s2", "s2", "s2", "s2"
  ),
  method = c(
    "R", "R", "A", "A", "B", "B",
    "R", "R", "A", "A", "B", "B"
  ),
  replicate = c(
    "r1", "r2", "r1", "r2", "r1", "r2",
    "r1", "r2", "r1", "r2", "r1", "r2"
  ),
  value = c(
    0, 2, 0, 0, 4, 4,
    10, 12, 10, 10, 14, 14
  )
)

dat_overall_bal <- data.frame(
  subject = c(
    rep("s1", 6),
    rep("s2", 6),
    rep("s3", 6)
  ),
  method = rep(c("R", "R", "A", "A", "B", "B"), times = 3),
  replicate = rep(c("r1", "r2"), times = 9),
  value = c(
    0, 2, 0, 0, 4, 4,
    10, 12, 10, 10, 14, 14,
    20, 24, 20, 20, 25, 25
  )
)

manual_overall_cia <- function(dat, reference = NULL, conf_level = 0.95) {
  subjects <- unique(dat$subject)
  methods <- unique(dat$method)
  n <- length(subjects)
  J <- length(methods)
  counts <- xtabs(~ subject + method, data = dat)
  stopifnot(length(unique(as.integer(counts))) == 1L)
  K <- unique(as.integer(counts))[[1L]]

  Ybar <- matrix(NA_real_, n, J, dimnames = list(subjects, methods))
  A <- Ybar
  for (i in seq_along(subjects)) {
    for (j in seq_along(methods)) {
      vals <- dat$value[dat$subject == subjects[[i]] & dat$method == methods[[j]]]
      mu <- mean(vals)
      Ybar[i, j] <- mu
      A[i, j] <- sum((vals - mu)^2) / (K - 1)
    }
  }

  if (is.null(reference)) {
    A_i <- rowMeans(A)
    Y_i_dot_dot <- rowMeans(Ybar)
    B_i <- rowSums((Ybar - Y_i_dot_dot)^2) / (J - 1) + (1 - 1 / K) * A_i
    A_bar <- mean(A_i)
    B_bar <- mean(B_i)
    psi <- A_bar / B_bar
    var_A <- stats::var(A_i) / n
    var_B <- stats::var(B_i) / n
    cov_A_B <- stats::cov(A_i, B_i) / n
    var_psi <-
      (1 / (B_bar^2)) * var_A +
      (A_bar^2 / B_bar^4) * var_B -
      (2 * A_bar / B_bar^3) * cov_A_B
    numerator_term <- A_bar
    denominator_term <- B_bar
  } else {
    ref_idx <- match(reference, methods)
    nonref_idx <- setdiff(seq_len(J), ref_idx)
    A_i <- A[, ref_idx]
    B_i <- rowSums((Ybar[, nonref_idx, drop = FALSE] - Ybar[, ref_idx])^2) / (J - 1) +
      (1 - 1 / K) * rowSums(A[, nonref_idx, drop = FALSE]) / (J - 1) +
      (1 - 1 / K) * A[, ref_idx]
    A_bar <- mean(A_i)
    B_bar <- mean(B_i)
    psi <- 2 * A_bar / B_bar
    var_A <- stats::var(A_i) / n
    var_B <- stats::var(B_i) / n
    cov_A_B <- stats::cov(A_i, B_i) / n
    var_ratio <-
      (1 / (B_bar^2)) * var_A +
      (A_bar^2 / B_bar^4) * var_B -
      (2 * A_bar / B_bar^3) * cov_A_B
    var_psi <- 4 * var_ratio
    numerator_term <- 2 * A_bar
    denominator_term <- B_bar
  }

  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  se <- sqrt(max(var_psi, 0))
  list(
    psi = psi,
    lwr = psi - z * se,
    upr = psi + z * se,
    se = se,
    A_i = A_i,
    B_i = B_i,
    A_bar = mean(A_i),
    B_bar = mean(B_i),
    numerator_term = numerator_term,
    denominator_term = denominator_term
  )
}

test_that("pairwise reference CIA uses the unhalved reference replicate disagreement", {
  fit <- suppressWarnings(cia(
    dat_ref2,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    reference = "R",
    scope = "pairwise",
    estimator = "mom_unconstrained"
  ))

  mat <- as.matrix(fit)
  diag_attr <- attr(fit, "diagnostics", exact = TRUE)

  expect_equal(mat["A", "R"], 2, tolerance = 1e-12)
  expect_equal(unname(diag_attr$within_msd["R"]), 4, tolerance = 1e-12)
  expect_equal(unname(diag_attr$sigma_within["R"]), 2, tolerance = 1e-12)
  expect_equal(unname(diag_attr$between_msd["A", "R"]), 2, tolerance = 1e-12)
})

test_that("overall reference CIA uses 2 * A_bar / B_bar", {
  manual <- manual_overall_cia(dat_overall_bal, reference = "R")
  fit <- cia(
    dat_overall_bal,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    reference = "R",
    scope = "overall",
    estimator = "mom_unconstrained"
  )

  expect_equal(fit$cia[[1]], manual$psi, tolerance = 1e-12)
  expect_equal(fit$numerator_term[[1]], manual$numerator_term, tolerance = 1e-12)
  expect_equal(fit$denominator_term[[1]], manual$denominator_term, tolerance = 1e-12)
  expect_equal(fit$within_msd[[1]], manual$numerator_term, tolerance = 1e-12)
  expect_equal(fit$between_msd[[1]], manual$denominator_term, tolerance = 1e-12)
})

test_that("overall no-reference CIA uses mean(A_i_dot) / mean(B_N_i)", {
  manual <- manual_overall_cia(dat_overall_bal, reference = NULL)
  fit <- cia(
    dat_overall_bal,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    scope = "overall",
    estimator = "mom_unconstrained"
  )

  expect_equal(fit$cia[[1]], manual$psi, tolerance = 1e-12)
  expect_equal(fit$numerator_term[[1]], manual$numerator_term, tolerance = 1e-12)
  expect_equal(fit$denominator_term[[1]], manual$denominator_term, tolerance = 1e-12)
})

test_that("pairwise reference CIA is correctly scaled for all non-reference methods", {
  fit <- suppressWarnings(cia(
    dat_ref3,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    reference = "R",
    scope = "pairwise",
    estimator = "mom_unconstrained"
  ))

  mat <- as.matrix(fit)

  expect_equal(mat["A", "R"], 2, tolerance = 1e-12)
  expect_equal(mat["B", "R"], 0.4, tolerance = 1e-12)
})

test_that("pairwise no-reference CIA uses the two-method subject-level estimator", {
  fit <- suppressWarnings(cia(
    dat_ref2,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    reference = NULL,
    scope = "pairwise",
    estimator = "mom_unconstrained"
  ))

  mat <- as.matrix(fit)

  expect_equal(mat["A", "R"], 1, tolerance = 1e-12)
})

test_that("C++ pairwise reference estimate uses the unhalved reference numerator", {
  cpp_fun <- get("cia_moments_cpp", envir = asNamespace("matrixCorr"))

  y <- c(
    0, 2, 0, 0, 4, 4,
    10, 12, 10, 10, 14, 14
  )
  subject_code <- c(
    1, 1, 1, 1, 1, 1,
    2, 2, 2, 2, 2, 2
  )
  method_code <- c(
    1, 1, 2, 2, 3, 3,
    1, 1, 2, 2, 3, 3
  )
  replicate_code <- c(
    1, 2, 1, 2, 1, 2,
    1, 2, 1, 2, 1, 2
  )

  raw <- cpp_fun(
    y = y,
    subject = subject_code,
    method = method_code,
    replicate = replicate_code,
    n_methods = 3L,
    reference_method = 1L,
    has_reference = TRUE,
    pairwise = TRUE,
    n_threads = 1L
  )

  expect_equal(raw$within_msd[[1]], 4, tolerance = 1e-12)
  expect_equal(raw$sigma_within[[1]], 2, tolerance = 1e-12)
  expect_equal(raw$between_msd[2, 1], 2, tolerance = 1e-12)
  expect_equal(raw$between_msd[3, 1], 10, tolerance = 1e-12)
  expect_equal(raw$estimate[2, 1], 2, tolerance = 1e-12)
  expect_equal(raw$estimate[3, 1], 0.4, tolerance = 1e-12)
})

test_that("unbalanced replication uses equal subject weighting in pairwise CIA", {
  dat_unbalanced <- data.frame(
    subject = c(
      rep("s1", 5),
      rep("s2", 4)
    ),
    method = c(
      "A", "A", "A", "B", "B",
      "A", "A", "B", "B"
    ),
    replicate = c(
      "r1", "r2", "r3", "r1", "r2",
      "r1", "r2", "r1", "r2"
    ),
    value = c(
      0, 2, 4, 0, 0,
      10, 12, 10, 10
    )
  )

  fit <- suppressWarnings(cia(
    dat_unbalanced,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    scope = "pairwise",
    estimator = "mom_unconstrained"
  ))

  expect_equal(as.matrix(fit)["A", "B"], 9 / 13, tolerance = 1e-12)
})

test_that("pairwise CIA uses delta normal CI by default", {
  dat_big <- do.call(
    rbind,
    lapply(seq_len(12), function(i) {
      data.frame(
        subject = rep(sprintf("s%02d", i), 4),
        method = c("A", "A", "B", "B"),
        replicate = c("r1", "r2", "r1", "r2"),
        value = c(i, i + 1, i + 0.1, i + 1.1)
      )
    })
  )

  fit <- cia(
    dat_big,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    scope = "pairwise",
    ci = TRUE
  )

  expect_identical(attr(fit, "ci.method", exact = TRUE), "delta_normal")
  expect_true(all(c("est", "lwr.ci", "upr.ci") %in% names(fit)))
})

test_that("pairwise bootstrap CIA remains available", {
  dat_big <- do.call(
    rbind,
    lapply(seq_len(12), function(i) {
      data.frame(
        subject = rep(sprintf("s%02d", i), 4),
        method = c("A", "A", "B", "B"),
        replicate = c("r1", "r2", "r1", "r2"),
        value = c(i, i + 1, i + 0.1, i + 1.1)
      )
    })
  )

  fit <- cia(
    dat_big,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    scope = "pairwise",
    ci = TRUE,
    inference = "bootstrap",
    B = 20L,
    seed = 1L
  )

  expect_identical(attr(fit, "ci.method", exact = TRUE), "subject_bootstrap_percentile")
  expect_true(all(c("est", "lwr.ci", "upr.ci") %in% names(fit)))
})

test_that("overall reference delta CI follows the paper's ratio delta method", {
  manual <- manual_overall_cia(dat_overall_bal, reference = "R")
  fit <- cia(
    dat_overall_bal,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    reference = "R",
    scope = "overall",
    estimator = "mom_unconstrained",
    ci = TRUE,
    inference = "delta"
  )

  expect_identical(attr(fit, "ci.method", exact = TRUE), "delta_normal")
  expect_equal(fit$cia[[1]], manual$psi, tolerance = 1e-12)
  expect_equal(fit$se[[1]], manual$se, tolerance = 1e-12)
  expect_equal(fit$lwr.ci[[1]], manual$lwr, tolerance = 1e-12)
  expect_equal(fit$upr.ci[[1]], manual$upr, tolerance = 1e-12)
})

test_that("overall no-reference delta CI follows the paper's ratio delta method", {
  manual <- manual_overall_cia(dat_overall_bal, reference = NULL)
  fit <- cia(
    dat_overall_bal,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    scope = "overall",
    estimator = "mom_unconstrained",
    ci = TRUE,
    inference = "delta"
  )

  expect_identical(attr(fit, "ci.method", exact = TRUE), "delta_normal")
  expect_equal(fit$cia[[1]], manual$psi, tolerance = 1e-12)
  expect_equal(fit$se[[1]], manual$se, tolerance = 1e-12)
  expect_equal(fit$lwr.ci[[1]], manual$lwr, tolerance = 1e-12)
  expect_equal(fit$upr.ci[[1]], manual$upr, tolerance = 1e-12)
})

test_that("overall CIA rejects unbalanced replication", {
  dat_bad <- dat_overall_bal[-1, , drop = FALSE]

  expect_error(
    cia(
      dat_bad,
      response = "value",
      subject = "subject",
      method = "method",
      replicate = "replicate",
      scope = "overall"
    ),
    "balanced-replication"
  )
})

test_that("overall bootstrap recomputes the corrected moment estimator", {
  fit <- cia(
    dat_overall_bal,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    reference = "R",
    scope = "overall",
    estimator = "mom_unconstrained",
    ci = TRUE,
    inference = "bootstrap",
    B = 20L,
    seed = 1L
  )
  fit_plain <- cia(
    dat_overall_bal,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    reference = "R",
    scope = "overall",
    estimator = "mom_unconstrained"
  )

  expect_identical(attr(fit, "ci.method", exact = TRUE), "subject_bootstrap_percentile")
  expect_equal(fit$cia[[1]], fit_plain$cia[[1]], tolerance = 1e-12)
  expect_true(all(c("lwr.ci", "upr.ci") %in% names(fit)))
})

test_that("constrained overall estimator rejects delta CI", {
  expect_error(
    cia(
      dat_overall_bal,
      response = "value",
      subject = "subject",
      method = "method",
      replicate = "replicate",
      scope = "overall",
      estimator = "vc_constrained",
      ci = TRUE,
      inference = "delta"
    ),
    "unconstrained moment estimator"
  )
})

test_that("reference pairwise CIA can exceed 1 under the raw MOM estimator", {
  dat_gt1 <- data.frame(
    subject = c(
      "s1", "s1", "s1", "s1",
      "s2", "s2", "s2", "s2"
    ),
    method = c(
      "R", "R", "A", "A",
      "R", "R", "A", "A"
    ),
    replicate = c(
      "r1", "r2", "r1", "r2",
      "r1", "r2", "r1", "r2"
    ),
    value = c(
      0, 2, 1, 1,
      10, 12, 11, 11
    )
  )

  fit <- suppressWarnings(cia(
    dat_gt1,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    reference = "R",
    scope = "pairwise",
    estimator = "mom_unconstrained"
  ))

  expect_gt(as.matrix(fit)["A", "R"], 1)
})

test_that("pairwise CIA warns when fewer than 10 eligible subjects are used", {
  expect_warning(
    cia(
      dat_ref2,
      response = "value",
      subject = "subject",
      method = "method",
      replicate = "replicate",
      scope = "pairwise",
      estimator = "mom_unconstrained"
    ),
    "fewer than 10 eligible subjects"
  )
})
