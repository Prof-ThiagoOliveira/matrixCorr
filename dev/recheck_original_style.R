#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pkgload)
  library(bench)
})

pkgload::load_all(".", quiet = TRUE, export_all = FALSE)

set.seed(20260413)

make_rmcorr_data <- function(n_subject = 200L, n_rep = 5L, p = 8L) {
  subject <- rep(seq_len(n_subject), each = n_rep)
  n <- length(subject)
  subj_eff <- matrix(rnorm(n_subject * p, sd = 0.8), n_subject, p)
  W <- matrix(rnorm(n * p), n, p)
  X <- matrix(0, n, p)
  for (j in seq_len(p)) {
    X[, j] <- subj_eff[subject, j] + 0.7 * W[, 1] + 0.3 * W[, j] + rnorm(n, sd = 0.2)
  }
  colnames(X) <- paste0("V", seq_len(p))
  out <- as.data.frame(X)
  out$subject <- subject
  out
}

make_ccc_data <- function(n_subject = 160L, n_time = 4L, methods = c("A", "B", "C")) {
  subj <- rep(seq_len(n_subject), each = n_time)
  tim <- rep(seq_len(n_time), times = n_subject)
  latent <- rnorm(n_subject, mean = 10, sd = 2)[subj] + 0.3 * tim
  make_method <- function(m) {
    bias <- switch(m, A = 0.0, B = 0.15, C = -0.10, 0.0)
    noise <- switch(m, A = 0.35, B = 0.45, C = 0.40, 0.4)
    data.frame(
      y = latent + bias + rnorm(length(latent), sd = noise),
      subject = subj,
      method = m,
      time = tim,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, lapply(methods, make_method))
}

make_latent_data <- function(n = 1600L) {
  z1 <- rnorm(n)
  z2 <- 0.55 * z1 + sqrt(1 - 0.55^2) * rnorm(n)
  z3 <- -0.35 * z1 + sqrt(1 - 0.35^2) * rnorm(n)
  z4 <- 0.45 * z2 + sqrt(1 - 0.45^2) * rnorm(n)
  data.frame(
    O1 = cut(z1, breaks = quantile(z1, probs = c(0, 0.3, 0.7, 1)), include.lowest = TRUE, ordered_result = TRUE),
    O2 = cut(z2, breaks = quantile(z2, probs = c(0, 0.25, 0.6, 1)), include.lowest = TRUE, ordered_result = TRUE),
    O3 = cut(z3, breaks = quantile(z3, probs = c(0, 0.4, 0.8, 1)), include.lowest = TRUE, ordered_result = TRUE),
    O4 = cut(z4, breaks = quantile(z4, probs = c(0, 0.2, 0.5, 0.8, 1)), include.lowest = TRUE, ordered_result = TRUE)
  )
}

n <- 1800L
p <- 32L
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("X", seq_len(p))
X <- X %*% chol(0.35^abs(outer(seq_len(p), seq_len(p), "-")))
X_na <- X
na_idx <- sample.int(length(X_na), size = floor(0.06 * length(X_na)))
X_na[na_idx] <- NA_real_
ccc <- make_ccc_data()
lat_ord <- make_latent_data()
rmc <- make_rmcorr_data()

# Match the previous harness call order that happened before benchmarks.
invisible(icc(X[, 1:8], scope = "pairwise", ci = TRUE, na_method = "error"))
invisible(icc(X[, 1:8], scope = "overall", ci = TRUE, na_method = "error"))
invisible(rmcorr(rmc, response = paste0("V", 1:6), subject = "subject", na_method = "error"))
invisible(ccc_rm_reml(ccc, response = "y", subject = "subject", method = "method", time = "time", ci = TRUE, ar = "none"))
invisible(ccc_rm_ustat(ccc, response = "y", subject = "subject", method = "method", time = "time", ci = TRUE))
invisible(dcor(X[, 1:20]))
invisible(dcor(X_na[, 1:20], na_method = "pairwise", p_value = TRUE))
invisible(pearson_corr(X_na[, 1:20], na_method = "pairwise", ci = TRUE))
invisible(polychoric(lat_ord, ci = TRUE, p_value = TRUE))
invisible(summary(dcor(X[, 1:20]), topn = 25L))

run_case <- function(label, expr, iterations) {
  expr_sub <- substitute(expr)
  expr_env <- parent.frame()
  gc(full = TRUE)
  bm <- bench::mark(
    result = eval(expr_sub, envir = expr_env),
    iterations = iterations,
    check = FALSE,
    memory = TRUE
  )
  data.frame(
    case = label,
    median_now = as.numeric(bm$median[[1]]),
    mem_now = as.numeric(bm$mem_alloc[[1]]),
    stringsAsFactors = FALSE
  )
}

now_tbl <- do.call(rbind, list(
  run_case("polychoric_ci_p", polychoric(lat_ord, ci = TRUE, p_value = TRUE), 4L),
  run_case("pearson_pairwise_ci", pearson_corr(X_na[, 1:20], na_method = "pairwise", ci = TRUE), 4L),
  run_case("ccc_rm_reml_ci", ccc_rm_reml(ccc, response = "y", subject = "subject", method = "method", time = "time", ci = TRUE, ar = "none"), 4L),
  run_case("ccc_rm_ustat_ci", ccc_rm_ustat(ccc, response = "y", subject = "subject", method = "method", time = "time", ci = TRUE), 4L),
  run_case("dcor_pairwise_inf", dcor(X_na[, 1:20], na_method = "pairwise", p_value = TRUE), 4L),
  run_case("summary_dcor_top", summary(dcor(X[, 1:20]), topn = 25L), 5L),
  run_case("dcor_matrix", dcor(X[, 1:20]), 5L)
))

baseline <- data.frame(
  case = c(
    "polychoric_ci_p",
    "pearson_pairwise_ci",
    "ccc_rm_reml_ci",
    "ccc_rm_ustat_ci",
    "dcor_pairwise_inf",
    "summary_dcor_top",
    "dcor_matrix"
  ),
  median_base = c(
    0.02940575,
    0.00352075,
    0.58700355,
    0.09991880,
    0.00518725,
    0.11176740,
    0.07924120
  ),
  mem_base = c(
    2801952,
    564560,
    5944736,
    2870456,
    314784,
    1472696,
    298544
  ),
  stringsAsFactors = FALSE
)

out <- merge(baseline, now_tbl, by = "case", all.x = TRUE, sort = FALSE)
out$speedup_vs_original <- out$median_base / out$median_now
out$mem_ratio_vs_original <- out$mem_now / out$mem_base

print(
  out[, c("case", "median_base", "median_now", "speedup_vs_original", "mem_base", "mem_now", "mem_ratio_vs_original")],
  row.names = FALSE
)

cat("\nOpenMP threads after run:", matrixCorr:::get_omp_threads(), "\n")
