#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pkgload)
})

pkgload::load_all(".", quiet = TRUE, export_all = FALSE)

set.seed(20260413)

n <- 300L
p <- 8L
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("V", seq_len(p))
X_na <- X
X_na[sample.int(length(X_na), size = floor(0.04 * length(X_na)))] <- NA_real_

subj <- rep(seq_len(80L), each = 4L)
tim <- rep(seq_len(4L), times = 80L)
latent <- rnorm(80L, 10, 2)[subj] + 0.2 * tim
mk <- function(m) {
  bias <- switch(m, A = 0.0, B = 0.1, C = -0.08, 0.0)
  noise <- switch(m, A = 0.35, B = 0.45, C = 0.4, 0.4)
  data.frame(
    y = latent + bias + rnorm(length(latent), sd = noise),
    subject = subj,
    method = m,
    time = tim,
    stringsAsFactors = FALSE
  )
}
ccc_df <- do.call(rbind, lapply(c("A", "B", "C"), mk))
ba2_df <- subset(ccc_df, method %in% c("A", "B"))

probe <- function(name, expr) {
  matrixCorr:::set_omp_threads(16L)
  before <- matrixCorr:::get_omp_threads()
  ok <- TRUE
  err <- ""
  tryCatch(force(expr), error = function(e) {
    ok <<- FALSE
    err <<- conditionMessage(e)
  })
  after <- matrixCorr:::get_omp_threads()
  data.frame(
    method = name,
    before = as.integer(before),
    after = as.integer(after),
    restored = identical(as.integer(before), as.integer(after)),
    ok = ok,
    error = err,
    stringsAsFactors = FALSE
  )
}

res <- do.call(rbind, list(
  probe("pearson_corr", pearson_corr(X_na[, 1:6], na_method = "pairwise", n_threads = 1L)),
  probe("spearman_rho", spearman_rho(X_na[, 1:6], na_method = "pairwise", n_threads = 1L)),
  probe("kendall_tau", kendall_tau(X_na[, 1:6], na_method = "pairwise", n_threads = 1L)),
  probe("dcor", dcor(X_na[, 1:6], na_method = "pairwise", n_threads = 1L)),
  probe("bicor", bicor(X[, 1:6], n_threads = 1L)),
  probe("wincor", wincor(X[, 1:6], n_threads = 1L)),
  probe("pbcor", pbcor(X[, 1:6], n_threads = 1L)),
  probe("skipped_corr", skipped_corr(X[, 1:6], n_threads = 1L, n_boot = 40L)),
  probe("schafer_corr", schafer_corr(X[, 1:6], n_threads = 1L)),
  probe("ccc", ccc(X[, 1:6], n_threads = 1L)),
  probe("icc_pairwise", icc(X[, 1:6], scope = "pairwise", n_threads = 1L)),
  probe("rmcorr_matrix", {
    rmc <- data.frame(X[, 1:4], subject = rep(seq_len(75L), each = 4L))
    rmcorr(rmc, response = c("V1", "V2", "V3", "V4"), subject = "subject", n_threads = 1L)
  }),
  probe("ccc_rm_ustat", ccc_rm_ustat(ccc_df, response = "y", subject = "subject", method = "method", time = "time", n_threads = 1L)),
  probe("ccc_rm_reml", ccc_rm_reml(ccc_df, response = "y", subject = "subject", method = "method", time = "time", n_threads = 1L)),
  probe("ba", ba(X[, 1], X[, 2], n_threads = 1L)),
  probe("ba_rm", ba_rm(ba2_df, response = "y", subject = "subject", method = "method", time = "time", n_threads = 1L))
))

print(res, row.names = FALSE)
if (any(!res$restored | !res$ok)) {
  quit(status = 1L)
}
