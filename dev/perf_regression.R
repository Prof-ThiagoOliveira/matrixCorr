#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pkgload)
  library(bench)
})

mode <- if (length(commandArgs(trailingOnly = TRUE)) >= 1L) {
  commandArgs(trailingOnly = TRUE)[1L]
} else {
  "baseline"
}

if (!mode %in% c("baseline", "compare")) {
  stop("Mode must be 'baseline' or 'compare'.")
}

repo_root <- normalizePath(".", winslash = "\\", mustWork = TRUE)
baseline_path <- file.path(repo_root, "dev", "perf_baseline.rds")
compare_path <- file.path(repo_root, "dev", "perf_compare.rds")

pkgload::load_all(repo_root, quiet = TRUE, export_all = FALSE)

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

make_inputs <- function() {
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

  list(
    X = X,
    X_na = X_na,
    mix_df = mix_df,
    rmc = rmc,
    ccc = ccc,
    lat = lat
  )
}

inputs <- make_inputs()

compute_outputs <- function(inp) {
  dcor_inf <- dcor(inp$X_na[, 1:20], na_method = "pairwise", p_value = TRUE)
  list(
    validate = matrixCorr:::validate_corr_input(inp$mix_df, check_na = FALSE),
    icc_pair = icc(inp$X[, 1:8], scope = "pairwise", ci = TRUE, na_method = "error"),
    icc_overall = icc(inp$X[, 1:8], scope = "overall", ci = TRUE, na_method = "error"),
    rmcorr_mat = rmcorr(inp$rmc, response = paste0("V", 1:6), subject = "subject", na_method = "error"),
    ccc_reml = ccc_rm_reml(
      inp$ccc,
      response = "y",
      subject = "subject",
      method = "method",
      time = "time",
      ci = TRUE,
      ar = "none"
    ),
    ccc_ustat = ccc_rm_ustat(
      inp$ccc,
      response = "y",
      subject = "subject",
      method = "method",
      time = "time",
      ci = TRUE
    ),
    dcor = dcor(inp$X[, 1:20]),
    dcor_inf = dcor_inf,
    pearson_pw = pearson_corr(inp$X_na[, 1:20], na_method = "pairwise", ci = TRUE),
    spearman_pw = spearman_rho(inp$X_na[, 1:18], na_method = "pairwise", ci = TRUE),
    polychoric = polychoric(inp$lat$ord, ci = TRUE, p_value = TRUE),
    polyserial = polyserial(inp$lat$cont, inp$lat$ord[, 1:2], ci = TRUE, p_value = TRUE),
    biserial = biserial(inp$lat$cont, inp$lat$bin[, 1:2], ci = TRUE, p_value = TRUE),
    summary_dcor = summary(dcor(inp$X[, 1:20]), topn = 25L)
  )
}

run_case <- function(label, expr, iterations = 3L) {
  expr_sub <- substitute(expr)
  expr_env <- parent.frame()
  runner <- function() eval(expr_sub, envir = expr_env)
  gc(full = TRUE)
  bm <- bench::mark(
    result = runner(),
    iterations = iterations,
    check = FALSE,
    memory = TRUE
  )
  data.frame(
    case = label,
    iterations = iterations,
    median = bm$median[[1]],
    itr_sec = bm$`itr/sec`[[1]],
    mem_alloc = bm$mem_alloc[[1]],
    stringsAsFactors = FALSE
  )
}

compute_benchmarks <- function(inp) {
  bench_rows <- list(
    run_case("validate_df_no_na", matrixCorr:::validate_corr_input(inp$mix_df, check_na = FALSE), iterations = 5L),
    run_case("icc_pairwise_ci", icc(inp$X[, 1:8], scope = "pairwise", ci = TRUE, na_method = "error"), iterations = 3L),
    run_case("icc_overall_ci", icc(inp$X[, 1:8], scope = "overall", ci = TRUE, na_method = "error"), iterations = 3L),
    run_case("rmcorr_matrix", rmcorr(inp$rmc, response = paste0("V", 1:6), subject = "subject", na_method = "error"), iterations = 3L),
    run_case("ccc_rm_reml_ci", ccc_rm_reml(inp$ccc, response = "y", subject = "subject", method = "method", time = "time", ci = TRUE, ar = "none"), iterations = 2L),
    run_case("ccc_rm_ustat_ci", ccc_rm_ustat(inp$ccc, response = "y", subject = "subject", method = "method", time = "time", ci = TRUE), iterations = 2L),
    run_case("dcor_matrix", dcor(inp$X[, 1:20]), iterations = 3L),
    run_case("dcor_pairwise_inf", dcor(inp$X_na[, 1:20], na_method = "pairwise", p_value = TRUE), iterations = 2L),
    run_case("pearson_pairwise_ci", pearson_corr(inp$X_na[, 1:20], na_method = "pairwise", ci = TRUE), iterations = 2L),
    run_case("spearman_pairwise_ci", spearman_rho(inp$X_na[, 1:18], na_method = "pairwise", ci = TRUE), iterations = 2L),
    run_case("polychoric_ci_p", polychoric(inp$lat$ord, ci = TRUE, p_value = TRUE), iterations = 2L),
    run_case("polyserial_ci_p", polyserial(inp$lat$cont, inp$lat$ord[, 1:2], ci = TRUE, p_value = TRUE), iterations = 2L),
    run_case("summary_dcor_top", summary(dcor(inp$X[, 1:20]), topn = 25L), iterations = 3L)
  )
  do.call(rbind, bench_rows)
}

compare_outputs <- function(before, after) {
  nms <- intersect(names(before), names(after))
  rows <- lapply(nms, function(nm) {
    ae <- all.equal(before[[nm]], after[[nm]], tolerance = 1e-12, check.attributes = TRUE)
    data.frame(
      object = nm,
      identical = isTRUE(ae),
      detail = if (isTRUE(ae)) "OK" else paste(ae, collapse = " | "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

compare_bench <- function(before, after) {
  merged <- merge(before, after, by = "case", suffixes = c("_before", "_after"), all = FALSE)
  transform(
    merged,
    speedup = as.numeric(median_before) / as.numeric(median_after),
    mem_ratio = as.numeric(mem_alloc_after) / as.numeric(mem_alloc_before)
  )
}

cat("Running mode:", mode, "\n")
outputs <- compute_outputs(inputs)
bench_tbl <- compute_benchmarks(inputs)

if (mode == "baseline") {
  saveRDS(
    list(
      outputs = outputs,
      benchmarks = bench_tbl,
      session = sessionInfo()
    ),
    baseline_path
  )
  cat("Baseline saved to:", baseline_path, "\n")
} else {
  if (!file.exists(baseline_path)) {
    stop("Baseline not found at: ", baseline_path, ". Run baseline mode first.")
  }
  baseline <- readRDS(baseline_path)
  out_cmp <- compare_outputs(baseline$outputs, outputs)
  bench_cmp <- compare_bench(baseline$benchmarks, bench_tbl)

  saveRDS(
    list(
      outputs_compare = out_cmp,
      benchmarks_before = baseline$benchmarks,
      benchmarks_after = bench_tbl,
      benchmarks_compare = bench_cmp,
      session = sessionInfo()
    ),
    compare_path
  )

  cat("\nOutput Equality\n")
  print(out_cmp, row.names = FALSE)

  cat("\nBenchmark Compare (speedup > 1 faster, mem_ratio < 1 lower allocation)\n")
  print(
    bench_cmp[, c("case", "speedup", "mem_ratio")],
    row.names = FALSE
  )
  cat("\nCompare report saved to:", compare_path, "\n")
}
