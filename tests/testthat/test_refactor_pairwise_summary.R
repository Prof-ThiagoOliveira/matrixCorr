expect_pairwise_summary_equal <- function(actual, expected) {
  expect_identical(names(actual), names(expected))
  expect_identical(class(actual), class(expected))
  expect_identical(setdiff(names(attributes(actual)), "row.names"),
                   setdiff(names(attributes(expected)), "row.names"))
  expect_identical(rownames(actual), rownames(expected))

  for (nm in names(expected)) {
    if (is.numeric(expected[[nm]])) {
      expect_equal(actual[[nm]], expected[[nm]], tolerance = 1e-12, info = nm)
    } else {
      expect_identical(actual[[nm]], expected[[nm]], info = nm)
    }
  }

  actual_attrs <- attributes(actual)
  expected_attrs <- attributes(expected)
  actual_attrs$row.names <- NULL
  expected_attrs$row.names <- NULL
  expect_equal(actual_attrs, expected_attrs, tolerance = 1e-12)
}

reference_pairwise_summary <- function(object,
                                       class_name,
                                       digits = 4,
                                       ci_digits = 3,
                                       p_digits = 4,
                                       show_ci = "yes",
                                       ci = attr(object, "ci", exact = TRUE),
                                       inf = attr(object, "inference", exact = TRUE),
                                       diagnostics = attr(object, "diagnostics", exact = TRUE),
                                       include_ci = identical(show_ci, "yes") && !is.null(ci),
                                       include_p = FALSE,
                                       extra = NULL,
                                       overview = matrixCorr:::.mc_summary_corr_matrix(object),
                                       attrs = list()) {
  est <- as.matrix(object)
  rn <- rownames(est)
  cn <- colnames(est)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(est)))
  if (is.null(cn)) cn <- as.character(seq_len(ncol(est)))

  rows <- vector("list", nrow(est) * (ncol(est) - 1L) / 2L)
  k <- 0L
  for (i in seq_len(nrow(est) - 1L)) {
    for (j in (i + 1L):ncol(est)) {
      k <- k + 1L
      rec <- list(
        var1 = rn[[i]],
        var2 = cn[[j]],
        estimate = round(est[i, j], digits)
      )
      if (is.list(diagnostics) && is.matrix(diagnostics$n_complete)) {
        rec$n_complete <- as.integer(diagnostics$n_complete[i, j])
      }
      if (include_ci) {
        rec$lwr <- if (!is.null(ci$lwr.ci) && is.finite(ci$lwr.ci[i, j])) {
          round(ci$lwr.ci[i, j], ci_digits)
        } else {
          NA_real_
        }
        rec$upr <- if (!is.null(ci$upr.ci) && is.finite(ci$upr.ci[i, j])) {
          round(ci$upr.ci[i, j], ci_digits)
        } else {
          NA_real_
        }
      }
      if (is.function(extra)) {
        rec <- c(rec, extra(i, j, inf, diagnostics))
      }
      rows[[k]] <- rec
    }
  }

  df <- do.call(rbind.data.frame, rows)
  rownames(df) <- NULL
  for (nm in intersect(c("estimate", "lwr", "upr", "statistic", "df",
                        "fisher_z", "p_value", "p_value_adjusted",
                        "skipped_prop"), names(df))) {
    df[[nm]] <- as.numeric(df[[nm]])
  }
  for (nm in intersect(c("n_complete", "skipped_n"), names(df))) {
    df[[nm]] <- as.integer(df[[nm]])
  }

  out <- matrixCorr:::.mc_finalize_summary_df(df, class_name = class_name)
  attr(out, "overview") <- overview
  attr(out, "has_ci") <- include_ci
  if (!is.null(include_p)) attr(out, "has_p") <- include_p
  attr(out, "conf.level") <- if (is.null(ci)) NA_real_ else ci$conf.level
  attr(out, "digits") <- digits
  attr(out, "ci_digits") <- ci_digits
  if (!is.null(p_digits)) attr(out, "p_digits") <- p_digits
  for (nm in names(attrs)) attr(out, nm) <- attrs[[nm]]
  out
}

test_that("shared pairwise summaries preserve estimator-specific data frames", {
  set.seed(20260522)
  X <- matrix(rnorm(72), nrow = 18, ncol = 4)
  colnames(X) <- paste0("S", seq_len(ncol(X)))

  pearson <- pearson_corr(X, ci = TRUE)
  expect_pairwise_summary_equal(
    matrixCorr:::.mc_pearson_pairwise_summary(pearson),
    reference_pairwise_summary(pearson, "summary.pearson_corr", p_digits = NULL, include_p = NULL)
  )

  spearman <- spearman_rho(X, ci = TRUE)
  expect_pairwise_summary_equal(
    matrixCorr:::.mc_spearman_pairwise_summary(spearman),
    reference_pairwise_summary(spearman, "summary.spearman_rho", p_digits = NULL, include_p = NULL)
  )

  kendall <- kendall_tau(X, ci = TRUE)
  expect_pairwise_summary_equal(
    matrixCorr:::.mc_kendall_pairwise_summary(kendall),
    reference_pairwise_summary(
      kendall,
      "summary.kendall_matrix",
      p_digits = NULL,
      include_p = NULL,
      attrs = list(ci.method = attr(kendall, "ci", exact = TRUE)$ci.method)
    )
  )

  bic <- bicor(X, ci = TRUE)
  bic_extra <- function(i, j, inf, diagnostics) {
    list(
      statistic = if (is.finite(inf$statistic[i, j])) round(inf$statistic[i, j], 4) else NA_real_,
      fisher_z = if (is.finite(inf$Z[i, j])) round(inf$Z[i, j], 4) else NA_real_,
      p_value = if (is.finite(inf$p_value[i, j])) round(inf$p_value[i, j], 4) else NA_real_
    )
  }
  expect_pairwise_summary_equal(
    matrixCorr:::.mc_bicor_pairwise_summary(bic),
    reference_pairwise_summary(
      bic,
      "summary.bicor",
      include_p = TRUE,
      extra = bic_extra,
      attrs = list(inference_method = attr(bic, "inference", exact = TRUE)$method)
    )
  )

  dc <- dcor(X, p_value = TRUE)
  dc_extra <- function(i, j, inf, diagnostics) {
    list(
      statistic = if (is.matrix(inf$statistic) && is.finite(inf$statistic[i, j])) round(inf$statistic[i, j], 4) else NA_real_,
      df = if (is.matrix(inf$parameter) && is.finite(inf$parameter[i, j])) round(inf$parameter[i, j], 4) else NA_real_,
      p_value = if (is.matrix(inf$p_value) && is.finite(inf$p_value[i, j])) round(inf$p_value[i, j], 4) else NA_real_
    )
  }
  expect_pairwise_summary_equal(
    matrixCorr:::.mc_dcor_pairwise_summary(dc),
    reference_pairwise_summary(
      dc,
      "summary.dcor",
      ci_digits = NULL,
      include_ci = FALSE,
      include_p = TRUE,
      extra = dc_extra,
      attrs = list(inference_method = attr(dc, "inference", exact = TRUE)$method)
    )
  )
})

test_that("shared robust and partial summaries preserve extra columns", {
  set.seed(20260523)
  X <- matrix(rnorm(72), nrow = 18, ncol = 4)
  colnames(X) <- paste0("R", seq_len(ncol(X)))

  robust_extra <- function(i, j, inf, diagnostics) {
    list(
      statistic = if (is.finite(inf$statistic[i, j])) round(inf$statistic[i, j], 4) else NA_real_,
      p_value = if (is.finite(inf$p_value[i, j])) round(inf$p_value[i, j], 4) else NA_real_
    )
  }

  pb <- pbcor(X, ci = TRUE, p_value = TRUE, n_boot = 20L, seed = 1L)
  expect_pairwise_summary_equal(
    matrixCorr:::.mc_pbcor_pairwise_summary(pb),
    reference_pairwise_summary(
      pb,
      "summary.pbcor",
      include_p = TRUE,
      extra = robust_extra,
      attrs = list(
        inference_method = attr(pb, "inference", exact = TRUE)$method,
        n_boot = attr(pb, "n_boot", exact = TRUE)
      )
    )
  )

  win <- wincor(X, ci = TRUE, p_value = TRUE, n_boot = 20L, seed = 2L)
  expect_pairwise_summary_equal(
    matrixCorr:::.mc_wincor_pairwise_summary(win),
    reference_pairwise_summary(
      win,
      "summary.wincor",
      include_p = TRUE,
      extra = robust_extra,
      attrs = list(
        inference_method = attr(win, "inference", exact = TRUE)$method,
        n_boot = attr(win, "n_boot", exact = TRUE)
      )
    )
  )

  pc <- pcorr(X, method = "sample", ci = TRUE, return_p_value = TRUE, return_details = TRUE)
  pc_extra <- function(i, j, inf, diagnostics) {
    list(p_value = pc$p_value[i, j])
  }
  expect_pairwise_summary_equal(
    matrixCorr:::.mc_partial_corr_pairwise_summary(pc),
    reference_pairwise_summary(
      pc$pcor,
      "summary.partial_corr",
      p_digits = NULL,
      ci = pc$ci,
      diagnostics = pc$diagnostics,
      include_p = NULL,
      extra = pc_extra,
      overview = matrixCorr:::.mc_summary_corr_matrix(pc$pcor)
    )
  )

  skipped <- skipped_corr(
    X,
    ci = TRUE,
    p_value = TRUE,
    n_boot = 20L,
    n_mc = 20L,
    seed = 3L
  )
  skip_extra <- function(i, j, inf, diagnostics) {
    out <- list(
      p_value = if (!is.null(inf$p_value) && is.finite(inf$p_value[i, j])) round(inf$p_value[i, j], 4) else NA_real_
    )
    if (!is.null(inf$p_value_adjusted)) {
      out$p_value_adjusted <- if (is.finite(inf$p_value_adjusted[i, j])) round(inf$p_value_adjusted[i, j], 4) else NA_real_
    }
    if (!is.null(inf$reject)) out$reject <- isTRUE(inf$reject[i, j])
    if (is.matrix(diagnostics$skipped_n)) out$skipped_n <- as.integer(diagnostics$skipped_n[i, j])
    if (is.matrix(diagnostics$skipped_prop)) out$skipped_prop <- round(diagnostics$skipped_prop[i, j], 4)
    out
  }
  expect_pairwise_summary_equal(
    matrixCorr:::.mc_skipcor_pairwise_summary(skipped),
    reference_pairwise_summary(
      skipped,
      "summary.skipped_corr",
      ci_digits = 2,
      include_p = TRUE,
      extra = skip_extra,
      attrs = list(
        inference_method = attr(skipped, "inference", exact = TRUE)$method,
        p_adjust = attr(skipped, "inference", exact = TRUE)$p_adjust,
        critical_p_value = attr(skipped, "inference", exact = TRUE)$critical_p_value %||% NA_real_
      )
    )
  )
})
