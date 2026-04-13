#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pkgload)
  library(bench)
})

repo_root <- normalizePath(".", winslash = "\\", mustWork = TRUE)
pkgload::load_all(repo_root, quiet = TRUE, export_all = FALSE)

set.seed(20260413)

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

n <- 1800L
p <- 32L
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("X", seq_len(p))
X <- X %*% chol(0.35^abs(outer(seq_len(p), seq_len(p), "-")))
X_na <- X
na_idx <- sample.int(length(X_na), size = floor(0.06 * length(X_na)))
X_na[na_idx] <- NA_real_
ccc <- make_ccc_data()

probe <- function(label, expr, iters = 3L) {
  expr_sub <- substitute(expr)
  expr_env <- parent.frame()
  gc(full = TRUE)
  bm <- bench::mark(
    result = eval(expr_sub, envir = expr_env),
    iterations = iters,
    check = FALSE,
    memory = TRUE
  )
  data.frame(
    case = label,
    median = as.numeric(bm$median[[1]]),
    mem_alloc = as.numeric(bm$mem_alloc[[1]]),
    stringsAsFactors = FALSE
  )
}

run_set <- function(n_threads) {
  rbind(
    probe(
      sprintf("dcor_matrix_t%d", n_threads),
      dcor(X[, 1:20], n_threads = n_threads),
      iters = 3L
    ),
    probe(
      sprintf("dcor_pairwise_inf_t%d", n_threads),
      dcor(X_na[, 1:20], na_method = "pairwise", p_value = TRUE, n_threads = n_threads),
      iters = 2L
    ),
    probe(
      sprintf("pearson_pairwise_ci_t%d", n_threads),
      pearson_corr(X_na[, 1:20], na_method = "pairwise", ci = TRUE, n_threads = n_threads),
      iters = 2L
    ),
    probe(
      sprintf("ccc_rm_reml_ci_t%d", n_threads),
      ccc_rm_reml(
        ccc,
        response = "y",
        subject = "subject",
        method = "method",
        time = "time",
        ci = TRUE,
        ar = "none",
        n_threads = n_threads
      ),
      iters = 2L
    ),
    probe(
      sprintf("ccc_rm_ustat_ci_t%d", n_threads),
      ccc_rm_ustat(
        ccc,
        response = "y",
        subject = "subject",
        method = "method",
        time = "time",
        ci = TRUE,
        n_threads = n_threads
      ),
      iters = 2L
    ),
    probe(
      sprintf("summary_dcor_top_t%d", n_threads),
      summary(dcor(X[, 1:20], n_threads = n_threads), topn = 25L),
      iters = 2L
    )
  )
}

res <- rbind(run_set(1L), run_set(16L))
print(res, row.names = FALSE)
