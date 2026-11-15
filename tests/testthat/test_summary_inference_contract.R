pearson_fisher_stat <- function(estimate, n_complete) {
  z <- suppressWarnings(atanh(estimate))
  z[!is.finite(z)] <- NA_real_
  stat <- rep(NA_real_, length(estimate))
  p <- rep(NA_real_, length(estimate))
  eligible <- is.finite(n_complete) & n_complete > 3 & is.finite(z)
  se <- 1 / sqrt(n_complete[eligible] - 3)
  stat[eligible] <- z[eligible] / se
  p[eligible] <- 2 * stats::pnorm(abs(stat[eligible]), lower.tail = FALSE)
  list(fisher_z = z, statistic = stat, p_value = p)
}

pair_field <- function(sm, item_a, item_b, field) {
  hit <- (sm$item1 == item_a & sm$item2 == item_b) |
    (sm$item1 == item_b & sm$item2 == item_a)
  sm[[field]][which(hit)[[1L]]]
}

expect_no_invented_inference <- function(sm) {
  expect_false("statistic" %in% names(sm))
  expect_false("p_value" %in% names(sm))
  expect_false("fisher_z" %in% names(sm))
  expect_false("z_statistic" %in% names(sm))
  expect_false(isTRUE(attr(sm, "has_p", exact = TRUE)))
}

test_that("Pearson summary uses the per-pair Fisher normal approximation only for n > 3", {
  set.seed(20261114)
  wide <- matrix(rnorm(40), nrow = 10, ncol = 4)
  colnames(wide) <- paste0("V", seq_len(ncol(wide)))
  narrow <- wide[1:3, 1:2, drop = FALSE]
  mixed <- wide
  mixed[4:10, 3] <- NA_real_

  for (output in c("matrix", "sparse", "edge_list")) {
    fit <- pearson_corr(
      wide,
      ci = FALSE,
      na_method = "pairwise",
      output = output,
      threshold = 0,
      diag = TRUE
    )
    sm <- summary(fit)
    expect_s3_class(sm, "summary.corr_result")
    expect_true("fisher_z" %in% names(sm))
    expect_true(isTRUE(attr(sm, "has_p", exact = TRUE)))
    ref <- pearson_fisher_stat(sm$estimate, sm$n_complete)
    expect_equal(sm$fisher_z, ref$fisher_z)
    expect_equal(sm$statistic, ref$statistic)
    expect_equal(sm$p_value, ref$p_value)
    expect_true(all(sm$n_complete > 3))
    expect_true(all(is.finite(sm$p_value)))
  }

  fit_n3 <- pearson_corr(narrow, ci = FALSE, na_method = "pairwise")
  sm_n3 <- summary(fit_n3)
  expect_identical(unique(sm_n3$n_complete), 3L)
  expect_true("fisher_z" %in% names(sm_n3))
  expect_false("statistic" %in% names(sm_n3))
  expect_false("p_value" %in% names(sm_n3))
  expect_false(isTRUE(attr(sm_n3, "has_p", exact = TRUE)))

  fit_mixed <- pearson_corr(mixed, ci = FALSE, na_method = "pairwise")
  sm_mixed <- summary(fit_mixed)
  ref_mixed <- pearson_fisher_stat(sm_mixed$estimate, sm_mixed$n_complete)
  expect_equal(sm_mixed$statistic, ref_mixed$statistic)
  expect_equal(sm_mixed$p_value, ref_mixed$p_value)
  expect_true(any(is.na(sm_mixed$p_value)))
  expect_true(any(is.finite(sm_mixed$p_value)))
  expect_true(isTRUE(attr(sm_mixed, "has_p", exact = TRUE)))
  expect_true(is.na(pair_field(sm_mixed, "V1", "V3", "p_value")))
  expect_true(is.finite(pair_field(sm_mixed, "V1", "V2", "p_value")))
})

test_that("stored Pearson inference is retrieved and not replaced by the Fisher fallback", {
  set.seed(11)
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  colnames(X) <- c("A", "B", "C")
  fit <- pearson_corr(X, ci = TRUE)
  stat <- matrix(c(0, 4.5, 1.25, 4.5, 0, -2, 1.25, -2, 0), 3, 3)
  p <- matrix(c(1, 0.02, 0.4, 0.02, 1, 0.11, 0.4, 0.11, 1), 3, 3)
  df <- matrix(7, 3, 3)
  dimnames(stat) <- dimnames(p) <- dimnames(df) <- dimnames(fit)
  attr(fit, "inference") <- list(statistic = stat, p_value = p, df = df)
  sm <- summary(fit)
  expect_identical(pair_field(sm, "A", "B", "statistic"), 4.5)
  expect_identical(pair_field(sm, "A", "B", "p_value"), 0.02)
  expect_identical(pair_field(sm, "A", "B", "df"), 7)
  expect_identical(pair_field(sm, "A", "B", "fisher_z"), atanh(pair_field(sm, "A", "B", "estimate")))
  expect_false(isTRUE(all.equal(
    pair_field(sm, "A", "B", "statistic"),
    atanh(pair_field(sm, "A", "B", "estimate")) / (1 / sqrt(10 - 3))
  )))
})

test_that("summaries do not manufacture inference when the estimator stored none", {
  set.seed(12)
  X <- matrix(rnorm(80), nrow = 20, ncol = 4)
  colnames(X) <- paste0("V", 1:4)

  expect_no_invented_inference(summary(spearman_rho(X, ci = FALSE)))
  expect_no_invented_inference(summary(kendall_tau(X, ci = FALSE)))
  expect_no_invented_inference(summary(pbcor(X, ci = FALSE, p_value = FALSE)))
  expect_no_invented_inference(summary(wincor(X, ci = FALSE, p_value = FALSE)))
  expect_no_invented_inference(summary(skipped_corr(
    X, method = "pearson", ci = FALSE, p_value = FALSE, n_boot = 20, seed = 1,
    return_masks = FALSE
  )))

  for (method in c("sample", "ridge", "oas", "glasso")) {
    fit <- pcorr(X, method = method, ci = FALSE, return_p_value = FALSE)
    expect_no_invented_inference(summary(fit))
  }

  for (fit in list(
    dcor(X, p_value = FALSE),
    bcdcor(X, p_value = FALSE),
    hsic(X),
    xi_corr(X, ci = FALSE),
    ccc(X, ci = FALSE)
  )) {
    sm <- summary(fit)
    expect_false("fisher_z" %in% names(sm))
  }
})

test_that("bicor summary keeps stored T, df, p and exposes Z as z_statistic", {
  set.seed(13)
  X <- matrix(rnorm(60), nrow = 20, ncol = 3)
  colnames(X) <- c("A", "B", "C")
  fit <- bicor(X, ci = TRUE)
  inf <- attr(fit, "inference", exact = TRUE)
  sm <- summary(fit)
  expect_false("fisher_z" %in% names(sm))
  expect_identical(pair_field(sm, "A", "B", "statistic"), inf$statistic[1, 2])
  expect_identical(pair_field(sm, "A", "B", "p_value"), inf$p_value[1, 2])
  expect_identical(pair_field(sm, "A", "B", "df"), inf$parameter[1, 2])
  expect_identical(pair_field(sm, "A", "B", "z_statistic"), inf$Z[1, 2])
  expect_false(isTRUE(all.equal(
    pair_field(sm, "A", "B", "statistic"),
    pair_field(sm, "A", "B", "z_statistic")
  )))
  ci <- attr(fit, "ci", exact = TRUE)
  expect_identical(pair_field(sm, "A", "B", "lwr"), ci$lwr.ci[1, 2])
  expect_identical(pair_field(sm, "A", "B", "upr"), ci$upr.ci[1, 2])
})

test_that("an explicit stored p-value column stays present when every value is NA", {
  set.seed(14)
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)
  colnames(X) <- paste0("V", seq_len(ncol(X)))
  fit <- spearman_rho(X, ci = FALSE)
  p <- matrix(NA_real_, 4, 4, dimnames = dimnames(fit))
  stat <- matrix(NA_real_, 4, 4, dimnames = dimnames(fit))
  attr(fit, "inference") <- list(
    method = "stored_unavailable_numeric",
    statistic = stat,
    p_value = p
  )
  sm <- summary(fit)
  expect_true("p_value" %in% names(sm))
  expect_true(all(is.na(sm$p_value)))
  expect_true(isTRUE(attr(sm, "has_p", exact = TRUE)))
  expect_identical(sm$inference_method, "stored_unavailable_numeric")
  expect_identical(attr(sm, "overview", exact = TRUE)$has_p, TRUE)

  empty <- spearman_rho(
    X, ci = FALSE, output = "edge_list", threshold = 1, diag = FALSE
  )
  attr(empty, "inference") <- attr(fit, "inference", exact = TRUE)
  sm_empty <- summary(empty)
  expect_identical(nrow(sm_empty), 0L)
  expect_true("p_value" %in% names(sm_empty))
  expect_true(isTRUE(attr(sm_empty, "has_p", exact = TRUE)))
})

test_that("CI-only non-Pearson summaries keep stored intervals and do not invent tests", {
  set.seed(15)
  X <- matrix(rnorm(36), nrow = 12, ncol = 3)
  colnames(X) <- c("A", "B", "C")
  fit_for <- function(output) {
    list(
      spearman_rho(X, ci = TRUE, output = output),
      kendall_tau(X, ci = TRUE, output = output),
      pbcor(X, ci = TRUE, p_value = FALSE, n_boot = 30L, seed = 15, output = output),
      wincor(X, ci = TRUE, p_value = FALSE, n_boot = 30L, seed = 15, output = output),
      skipped_corr(
        X, method = "pearson", ci = TRUE, p_value = FALSE,
        n_boot = 20L, seed = 15, return_masks = FALSE, output = output
      )
    )
  }
  for (output in c("matrix", "sparse", "edge_list")) {
    for (fit in fit_for(output)) {
      sm <- summary(fit)
      expect_no_invented_inference(sm)
      expect_true(all(c("lwr", "upr") %in% names(sm)))
      expect_true(isTRUE(attr(sm, "has_ci", exact = TRUE)))
      ci <- attr(fit, "ci", exact = TRUE)
      dn <- dimnames(ci$lwr.ci)[[1L]]
      for (k in seq_len(nrow(sm))) {
        i <- match(sm$item1[[k]], dn)
        j <- match(sm$item2[[k]], dn)
        expect_identical(sm$lwr[[k]], ci$lwr.ci[i, j])
        expect_identical(sm$upr[[k]], ci$upr.ci[i, j])
      }
    }
  }
})

test_that("a stored robust p-value payload is copied into the summary without recalculation", {
  set.seed(16)
  X <- matrix(rnorm(36), nrow = 12, ncol = 3)
  colnames(X) <- c("A", "B", "C")
  fit <- pbcor(X, ci = FALSE, p_value = TRUE, n_boot = 40L, seed = 16)
  inf <- attr(fit, "inference", exact = TRUE)
  sm <- summary(fit)
  expect_true(isTRUE(attr(sm, "has_p", exact = TRUE)))
  expect_identical(sm$inference_method, inf$method)
  dn <- dimnames(inf$p_value)[[1L]]
  for (k in seq_len(nrow(sm))) {
    i <- match(sm$item1[[k]], dn)
    j <- match(sm$item2[[k]], dn)
    expect_identical(sm$p_value[[k]], inf$p_value[i, j])
    expect_identical(sm$statistic[[k]], inf$statistic[i, j])
    expect_identical(sm$df[[k]], inf$parameter[i, j])
  }
})

test_that("empty non-Pearson summaries keep baseline attribute order", {
  expected_attr_order <- c("names", "row.names", "class")

  X <- cbind(A = 1:4, B = 4:1)
  for (output in c("edge_list", "sparse")) {
    sm <- summary(shrinkage_corr(
      X, output = output, threshold = 1.01, diag = FALSE
    ))
    expect_identical(nrow(sm), 0L)
    expect_identical(head(names(attributes(sm)), 3L), expected_attr_order)
    expect_identical(names(sm), c("item1", "item2", "estimate"))
  }
})

test_that("summary reuses the diagonal-inclusive pair table", {
  x <- seq(-1, 1, length.out = 16)
  directed <- xi_corr(
    cbind(x = x, x2 = x^2),
    tie_method = "first",
    ci = FALSE,
    output = "edge_list",
    threshold = 0,
    diag = TRUE
  )
  directed_summary <- summary(directed)
  directed_overview <- attr(directed_summary, "overview", exact = TRUE)
  expect_false(any(directed_summary$item1 == directed_summary$item2))
  expect_gt(directed_overview$n_pairs, nrow(directed_summary))

  dn <- list(c("a", "b"), c("a", "b"))
  mat <- matrix(c(0.5, 0.2, 0.2, NaN), 2, 2, dimnames = dn)
  cross <- matrixCorr:::.mc_new_corr_matrix(
    mat,
    estimator_class = "biserial_corr",
    method = "biserial",
    description = "cross-axis diagonal",
    symmetric = TRUE
  )
  same <- matrixCorr:::.mc_new_corr_matrix(
    mat,
    estimator_class = "pearson_corr",
    method = "pearson",
    description = "same-axis diagonal",
    symmetric = TRUE
  )
  cross_summary <- summary(cross)
  same_summary <- summary(same)
  expect_true(any(cross_summary$item1 == cross_summary$item2 & is.finite(cross_summary$estimate)))
  expect_false(any(same_summary$item1 == same_summary$item2))
  expect_false(any(is.nan(cross_summary$estimate)))

  sparse <- Matrix::sparseMatrix(
    i = 2L, j = 1L, x = 0.2, dims = c(2L, 2L), dimnames = dn, symmetric = TRUE
  )
  sparse <- matrixCorr:::.mc_new_corr_sparse(
    sparse,
    estimator_class = "pearson_corr",
    method = "pearson",
    description = "stored lower only",
    symmetric = TRUE
  )
  sparse_tidy <- tidy(sparse)
  expect_identical(nrow(sparse_tidy), 1L)
  expect_identical(sparse_tidy$estimate, 0.2)
  expect_false(any(sparse_tidy$estimate == 0))
})

test_that("symmetric pair counts keep rows with missing positions", {
  missing_pos <- data.frame(
    row = c("1", NA_character_, "2"),
    col = c("2", "2", "2"),
    value = c(0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )
  missing_pos <- matrixCorr:::.mc_new_corr_edge_list(
    missing_pos,
    estimator_class = "pearson_corr",
    method = "pearson",
    description = "missing position",
    threshold = 0,
    diag = TRUE,
    source_dim = c(3L, 3L),
    source_dimnames = list(NULL, NULL),
    symmetric = TRUE
  )
  overview <- matrixCorr:::.mc_corr_overview(missing_pos)
  expect_identical(overview$n_pairs, 2L)
  expect_identical(overview$n_pairs_retained, 2L)
})
