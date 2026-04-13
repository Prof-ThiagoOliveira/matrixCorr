#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pkgload)
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
  ord <- data.frame(
    O1 = cut(z1, breaks = quantile(z1, probs = c(0, 0.3, 0.7, 1)), include.lowest = TRUE, ordered_result = TRUE),
    O2 = cut(z2, breaks = quantile(z2, probs = c(0, 0.25, 0.6, 1)), include.lowest = TRUE, ordered_result = TRUE),
    O3 = cut(z3, breaks = quantile(z3, probs = c(0, 0.4, 0.8, 1)), include.lowest = TRUE, ordered_result = TRUE),
    O4 = cut(z4, breaks = quantile(z4, probs = c(0, 0.2, 0.5, 0.8, 1)), include.lowest = TRUE, ordered_result = TRUE)
  )
  cont <- data.frame(
    X1 = scale(z1 + rnorm(n, sd = 0.3))[, 1],
    X2 = scale(z2 + rnorm(n, sd = 0.3))[, 1],
    X3 = scale(z3 + rnorm(n, sd = 0.3))[, 1],
    X4 = scale(z4 + rnorm(n, sd = 0.3))[, 1]
  )
  bin <- data.frame(
    B1 = as.integer(z1 > median(z1)),
    B2 = as.integer(z2 > median(z2))
  )
  list(ord = ord, cont = cont, bin = bin)
}

n <- 1800L
p <- 32L
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("X", seq_len(p))
X <- X %*% chol(0.35^abs(outer(seq_len(p), seq_len(p), "-")))
X_na <- X
na_idx <- sample.int(length(X_na), size = floor(0.06 * length(X_na)))
X_na[na_idx] <- NA_real_

mix_df <- as.data.frame(X)
mix_df$grp <- sample(letters[1:4], nrow(mix_df), replace = TRUE)
mix_df$flag <- sample(c(TRUE, FALSE), nrow(mix_df), replace = TRUE)
mix_df$date <- as.Date("2024-01-01") + seq_len(nrow(mix_df))

rmc <- make_rmcorr_data()
ccc <- make_ccc_data()
lat <- make_latent_data()

call_and_trace <- function(label, expr) {
  cat(sprintf("%-20s before=%d\n", label, matrixCorr:::get_omp_threads()))
  force(expr)
  cat(sprintf("%-20s after =%d\n", label, matrixCorr:::get_omp_threads()))
}

cat("initial threads:", matrixCorr:::get_omp_threads(), "\n")

call_and_trace("validate", matrixCorr:::validate_corr_input(mix_df, check_na = FALSE))
call_and_trace("icc_pair", icc(X[, 1:8], scope = "pairwise", ci = TRUE, na_method = "error"))
call_and_trace("icc_overall", icc(X[, 1:8], scope = "overall", ci = TRUE, na_method = "error"))
call_and_trace("rmcorr", rmcorr(rmc, response = paste0("V", 1:6), subject = "subject", na_method = "error"))
call_and_trace("ccc_reml", ccc_rm_reml(ccc, response = "y", subject = "subject", method = "method", time = "time", ci = TRUE, ar = "none"))
call_and_trace("ccc_ustat", ccc_rm_ustat(ccc, response = "y", subject = "subject", method = "method", time = "time", ci = TRUE))
call_and_trace("dcor", dcor(X[, 1:20]))
call_and_trace("dcor_inf", dcor(X_na[, 1:20], na_method = "pairwise", p_value = TRUE))
call_and_trace("pearson", pearson_corr(X_na[, 1:20], na_method = "pairwise", ci = TRUE))
call_and_trace("spearman", spearman_rho(X_na[, 1:18], na_method = "pairwise", ci = TRUE))
call_and_trace("polychoric", polychoric(lat$ord, ci = TRUE, p_value = TRUE))
call_and_trace("polyserial", polyserial(lat$cont, lat$ord[, 1:2], ci = TRUE, p_value = TRUE))
call_and_trace("biserial", biserial(lat$cont, lat$bin[, 1:2], ci = TRUE, p_value = TRUE))
call_and_trace("summary_dcor", summary(dcor(X[, 1:20]), topn = 25L))

cat("final threads:", matrixCorr:::get_omp_threads(), "\n")
