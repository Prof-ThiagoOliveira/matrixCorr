args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) >= 1L) normalizePath(args[[1L]]) else .libPaths()[1L]
iterations <- if (length(args) >= 2L) as.integer(args[[2L]]) else 5L

.libPaths(c(lib, .libPaths()))

library(bench)
library(matrixCorr)

set.seed(1)

X20 <- matrix(rnorm(250 * 20), 250, 20)
colnames(X20) <- paste0("V", seq_len(ncol(X20)))

X30 <- matrix(rnorm(300 * 30), 300, 30)
colnames(X30) <- paste0("V", seq_len(ncol(X30)))

X50 <- matrix(rnorm(400 * 50), 400, 50)
colnames(X50) <- paste0("V", seq_len(ncol(X50)))

Xsh <- matrix(rnorm(80 * 120), 80, 120)
colnames(Xsh) <- paste0("V", seq_len(ncol(Xsh)))

Xbin <- data.frame(
  b1 = sample(0:1, 120, TRUE),
  b2 = sample(0:1, 120, TRUE),
  b3 = sample(0:1, 120, TRUE)
)

Xord <- data.frame(
  o1 = sample(1:4, 120, TRUE),
  o2 = sample(1:5, 120, TRUE),
  o3 = sample(1:4, 120, TRUE)
)

Xcont <- matrix(rnorm(120 * 3), 120, 3)
colnames(Xcont) <- c("x1", "x2", "x3")
y_bin <- sample(0:1, 120, TRUE)
y_ord <- sample(1:4, 120, TRUE)

g1 <- rnorm(400)
g2 <- g1 + rnorm(400, sd = 0.5)

n_subjects <- 30L
n_rep <- 4L
subject <- rep(seq_len(n_subjects), each = n_rep)
within_signal <- rep(rep(c(-1, 1), each = n_rep / 2L), n_subjects)
subj_eff_x <- rnorm(n_subjects)
subj_eff_y <- 0.6 * subj_eff_x + rnorm(n_subjects, sd = 0.4)

dat_rm <- data.frame(
  subject = subject,
  x = subj_eff_x[subject] + within_signal + rnorm(n_subjects * n_rep, sd = 0.2),
  y = subj_eff_y[subject] + 0.8 * within_signal + rnorm(n_subjects * n_rep, sd = 0.3),
  z = subj_eff_y[subject] - 0.4 * within_signal + rnorm(n_subjects * n_rep, sd = 0.4)
)

dat_agree <- expand.grid(
  id = seq_len(24),
  time = seq_len(2),
  method = factor(c("A", "B"))
)
subj_eff <- rnorm(24)[dat_agree$id]
time_eff <- c(0, 0.15)[dat_agree$time]
method_eff <- c(0, 0.2)[as.integer(dat_agree$method)]
dat_agree$y <- subj_eff + time_eff + method_eff + rnorm(nrow(dat_agree), sd = 0.35)

bench_cases <- list(
  pearson_corr = quote(pearson_corr(X50)),
  spearman_rho = quote(spearman_rho(X50)),
  kendall_tau = quote(kendall_tau(X20)),
  dcor = quote(dcor(X20)),
  pcorr = quote(pcorr(X20[, 1:10])),
  bicor = quote(bicor(X30)),
  pbcor = quote(pbcor(X30)),
  wincor = quote(wincor(X30)),
  skipped_corr = quote(skipped_corr(X20[, 1:8])),
  schafer_corr = quote(schafer_corr(Xsh)),
  tetrachoric = quote(tetrachoric(Xbin)),
  polychoric = quote(polychoric(Xord)),
  biserial = quote(biserial(Xcont, y_bin)),
  polyserial = quote(polyserial(Xcont, y_ord)),
  rmcorr = quote(rmcorr(dat_rm, response = c("x", "y", "z"), subject = "subject")),
  ccc = quote(ccc(X30)),
  ba = quote(ba(g1, g2)),
  ba_rm = quote(ba_rm(dat_agree, response = "y", subject = "id", method = "method", time = "time")),
  ccc_rm_ustat = quote(ccc_rm_ustat(dat_agree, response = "y", subject = "id", method = "method", time = "time")),
  ccc_rm_reml = quote(ccc_rm_reml(dat_agree, "y", "id", method = "method", time = "time"))
)

results <- lapply(names(bench_cases), function(name) {
  b <- bench::mark(eval(bench_cases[[name]], envir = environment()), check = FALSE, iterations = iterations)
  s <- summary(b)[1L, c("median", "mem_alloc")]
  data.frame(
    expression = name,
    median = as.character(s$median),
    mem_alloc = as.character(s$mem_alloc),
    stringsAsFactors = FALSE
  )
})

out <- do.call(rbind, results)
rownames(out) <- NULL
print(out)
