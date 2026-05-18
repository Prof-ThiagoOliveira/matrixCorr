make_cia_perf_data <- function(n_subjects = 12L, n_methods = 3L, n_replicates = 3L) {
  subjects <- sprintf("s%02d", seq_len(n_subjects))
  methods <- LETTERS[seq_len(n_methods)]
  reps <- sprintf("r%02d", seq_len(n_replicates))
  out <- expand.grid(
    subject = subjects,
    method = methods,
    replicate = reps,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  s_idx <- match(out$subject, subjects)
  m_idx <- match(out$method, methods)
  r_idx <- match(out$replicate, reps)
  out$value <- 4 + 0.7 * s_idx + 0.2 * m_idx + 0.05 * r_idx +
    sin(s_idx + m_idx + r_idx) / 10
  out
}

cia_perf_prep <- function(dat, reference = NULL) {
  matrixCorr:::.mc_cia_prepare_long(
    data = dat,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    reference = reference
  )
}

cia_perf_diagnostics <- function(prep) {
  raw <- matrixCorr:::cia_moments_cpp(
    y = prep$y,
    subject = prep$subject_code,
    method = prep$method_code,
    replicate = prep$replicate_code,
    n_methods = prep$n_methods,
    reference_method = prep$reference_index,
    has_reference = prep$has_reference,
    pairwise = TRUE,
    n_threads = 1L
  )
  matrixCorr:::.mc_cia_build_diagnostics(prep, raw)
}

test_that("pairwise payload reuse preserves matrix output", {
  dat <- make_cia_perf_data()
  prep <- cia_perf_prep(dat)
  diagnostics <- cia_perf_diagnostics(prep)
  payload <- matrixCorr:::.mc_cia_pairwise_payload(
    prep,
    diagnostics = diagnostics,
    estimator = "mom_unconstrained",
    n_threads = 1L
  )

  old_path <- matrixCorr:::.mc_cia_build_pairwise_matrix(
    prep,
    diagnostics,
    estimator = "mom_unconstrained",
    n_threads = 1L
  )
  new_path <- matrixCorr:::.mc_cia_build_pairwise_matrix(
    prep,
    diagnostics,
    estimator = "mom_unconstrained",
    n_threads = 1L,
    payload = payload
  )

  expect_identical(class(new_path), class(old_path))
  expect_identical(dimnames(new_path), dimnames(old_path))
  expect_identical(attributes(new_path), attributes(old_path))
  expect_equal(unclass(new_path), unclass(old_path), tolerance = 1e-12)
})

test_that("pairwise payload reuse preserves delta CI output with reference", {
  dat <- make_cia_perf_data()
  prep <- cia_perf_prep(dat, reference = "A")
  diagnostics <- cia_perf_diagnostics(prep)
  payload <- matrixCorr:::.mc_cia_pairwise_payload(
    prep,
    diagnostics = diagnostics,
    estimator = "vc_constrained",
    conf_level = 0.9,
    return_delta = TRUE,
    n_threads = 1L
  )

  old_path <- matrixCorr:::.mc_cia_build_pairwise_ci(
    prep,
    diagnostics,
    estimator = "vc_constrained",
    inference = "delta",
    conf_level = 0.9,
    n_threads = 1L
  )
  new_path <- matrixCorr:::.mc_cia_build_pairwise_ci(
    prep,
    diagnostics,
    estimator = "vc_constrained",
    inference = "delta",
    conf_level = 0.9,
    n_threads = 1L,
    payload = payload
  )

  expect_identical(class(new_path), class(old_path))
  expect_identical(names(new_path), names(old_path))
  expect_identical(attributes(new_path), attributes(old_path))
  expect_equal(new_path, old_path, tolerance = 1e-12)
})

test_that("pairwise bootstrap remains reproducible and reports retained resamples", {
  dat <- make_cia_perf_data()
  fit1 <- cia(
    dat,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    scope = "pairwise",
    ci = TRUE,
    inference = "bootstrap",
    B = 49L,
    seed = 42L
  )
  fit2 <- cia(
    dat,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    scope = "pairwise",
    ci = TRUE,
    inference = "bootstrap",
    B = 49L,
    seed = 42L
  )

  expect_equal(fit1, fit2, tolerance = 1e-12)
  expect_identical(attr(fit1, "bootstrap", exact = TRUE)$B, 49L)
  expect_identical(attr(fit1, "bootstrap", exact = TRUE)$n_successful, 49L)
})

test_that("reference pairwise bootstrap remains reproducible", {
  dat <- make_cia_perf_data()
  fit1 <- cia(
    dat,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    reference = "A",
    scope = "pairwise",
    estimator = "vc_constrained",
    ci = TRUE,
    inference = "bootstrap",
    B = 49L,
    seed = 42L
  )
  fit2 <- cia(
    dat,
    response = "value",
    subject = "subject",
    method = "method",
    replicate = "replicate",
    reference = "A",
    scope = "pairwise",
    estimator = "vc_constrained",
    ci = TRUE,
    inference = "bootstrap",
    B = 49L,
    seed = 42L
  )

  expect_equal(fit1, fit2, tolerance = 1e-12)
  expect_identical(attr(fit1, "bootstrap", exact = TRUE), list(n_successful = 49L, B = 49L))
})

test_that("overall CIA public paths remain stable", {
  dat <- make_cia_perf_data()
  plain <- cia(dat, "value", "subject", "method", "replicate", scope = "overall")
  ref <- cia(dat, "value", "subject", "method", "replicate", reference = "A", scope = "overall")
  delta <- cia(dat, "value", "subject", "method", "replicate", scope = "overall", ci = TRUE, inference = "delta")
  boot <- cia(dat, "value", "subject", "method", "replicate", scope = "overall", ci = TRUE,
              inference = "bootstrap", B = 49L, seed = 42L)

  expect_s3_class(plain, "cia_overall")
  expect_s3_class(ref, "cia_overall")
  expect_identical(attr(delta, "ci.method", exact = TRUE), "delta_normal")
  expect_identical(attr(boot, "bootstrap", exact = TRUE), list(n_successful = 49L, B = 49L))
  expect_named(
    plain,
    c("coefficient", "estimator", "cia", "cia_raw", "tau2", "tau2_raw",
      "boundary", "reference", "n_methods", "n_subjects", "n_obs", "K",
      "numerator_term", "denominator_term", "within_msd", "between_msd")
  )
})

test_that("CIA error and warning behaviour remains unchanged", {
  dat <- make_cia_perf_data()
  dup <- rbind(dat, dat[1L, , drop = FALSE])
  one_rep <- subset(dat, replicate == "r01")
  unbalanced <- dat[-1L, , drop = FALSE]
  small <- make_cia_perf_data(n_subjects = 8L)

  expect_error(
    cia(dup, "value", "subject", "method", "replicate"),
    "duplicated subject-method-replicate"
  )
  expect_error(
    cia(one_rep, "value", "subject", "method", "replicate"),
    "at least two replicated readings"
  )
  expect_error(
    cia(unbalanced, "value", "subject", "method", "replicate", scope = "overall"),
    "balanced-replication"
  )
  expect_error(
    cia(dat, "value", "subject", "method", "replicate", reference = "missing"),
    "reference"
  )
  expect_error(
    cia(dat, "value", "subject", "method", "replicate", ci = TRUE, inference = "none"),
    "must not be"
  )
  expect_error(
    cia(dat, "value", "subject", "method", "replicate", scope = "overall",
        estimator = "vc_constrained", ci = TRUE, inference = "delta"),
    "unconstrained moment estimator only"
  )
  expect_warning(
    cia(small, "value", "subject", "method", "replicate", scope = "pairwise"),
    "fewer than 10 eligible subjects"
  )
})
