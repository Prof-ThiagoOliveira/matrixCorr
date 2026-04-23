#!/usr/bin/env Rscript

baseline_sha <- "9c47c4fcff7ff22328213077ce7eab3684c95566"

required_pkgs <- c("Rcpp", "bench", "profmem", "pkgload")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Install required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  stop("Run this script from the package root (directory containing DESCRIPTION).", call. = FALSE)
}

git_show_lines <- function(sha, path, repo) {
  old <- setwd(repo)
  on.exit(setwd(old), add = TRUE)
  out <- system2("git", c("show", sprintf("%s:%s", sha, path)),
                 stdout = TRUE, stderr = TRUE)
  st <- attr(out, "status")
  if (!is.null(st) && st != 0L) {
    stop("Failed to read ", path, " from commit ", sha, ".\n",
         paste(out, collapse = "\n"),
         call. = FALSE)
  }
  out
}

current_file_lines <- function(path, repo) {
  readLines(file.path(repo, path), warn = FALSE)
}

build_impl_env <- function(read_path_lines, ns) {
  tmp <- tempfile("matrixCorr_baseline_")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)

  cpp_files <- c(
    "src/matrixCorr_detail.h",
    "src/matrixCorr_omp.h",
    "src/threshold_triplets.h",
    "src/spearman.cpp",
    "src/kendall_corr.cpp"
  )
  for (f in cpp_files) {
    txt <- read_path_lines(f)
    writeLines(txt, con = file.path(tmp, basename(f)), useBytes = TRUE)
  }

  cpp_env <- new.env(parent = baseenv())
  Rcpp::sourceCpp(file.path(tmp, "spearman.cpp"), env = cpp_env, rebuild = TRUE, showOutput = FALSE)
  Rcpp::sourceCpp(file.path(tmp, "kendall_corr.cpp"), env = cpp_env, rebuild = TRUE, showOutput = FALSE)

  wrappers_env <- new.env(parent = ns)
  eval(parse(text = paste(read_path_lines("R/spearman_rho.R"), collapse = "\n")),
       envir = wrappers_env)
  eval(parse(text = paste(read_path_lines("R/kendall_corr.R"), collapse = "\n")),
       envir = wrappers_env)

  wrappers_env$spearman_matrix_cpp <- cpp_env$spearman_matrix_cpp
  wrappers_env$spearman_matrix_pairwise_cpp <- cpp_env$spearman_matrix_pairwise_cpp
  wrappers_env$spearman_threshold_triplets_cpp <- cpp_env$spearman_threshold_triplets_cpp
  wrappers_env$kendall_matrix_cpp <- cpp_env$kendall_matrix_cpp
  wrappers_env$kendall_tau2_cpp <- cpp_env$kendall_tau2_cpp
  wrappers_env$kendall_tau2_from_mat_cpp <- cpp_env$kendall_tau2_from_mat_cpp
  wrappers_env$kendall_matrix_pairwise_cpp <- cpp_env$kendall_matrix_pairwise_cpp

  wrappers_env
}

validate_contract <- function(before, after) {
  bm <- as.matrix(before)
  am <- as.matrix(after)

  checks <- list(
    same_dim = identical(dim(bm), dim(am)),
    same_dimnames = identical(dimnames(bm), dimnames(am)),
    same_storage_mode = identical(storage.mode(bm), storage.mode(am)),
    same_class = identical(class(before), class(after)),
    same_diagonal = identical(diag(bm), diag(am)),
    same_symmetry = identical(isSymmetric(bm), isSymmetric(am))
  )

  id <- identical(before, after)
  ae <- isTRUE(all.equal(before, after, tolerance = 1e-12, check.attributes = TRUE))
  list(
    checks = checks,
    identical = id,
    all_equal_1e12 = ae
  )
}

bench_one <- function(label, before_expr, after_expr, envir, iterations = 5L) {
  b <- bench::mark(
    before = eval(before_expr, envir = envir),
    after = eval(after_expr, envir = envir),
    iterations = iterations,
    check = FALSE,
    memory = TRUE,
    time_unit = "ms",
    gc = TRUE
  )

  expr_chr <- as.character(b$expression)
  before_row <- b[expr_chr == "before", , drop = FALSE]
  after_row <- b[expr_chr == "after", , drop = FALSE]
  if (nrow(before_row) != 1L || nrow(after_row) != 1L) {
    stop("Unexpected bench::mark result rows for label: ", label, call. = FALSE)
  }

  before_ms <- as.numeric(before_row$median)
  after_ms <- as.numeric(after_row$median)
  pct <- 100 * (before_ms - after_ms) / before_ms

  data.frame(
    label = label,
    before_ms = before_ms,
    after_ms = after_ms,
    pct_improvement = pct,
    before_mem_alloc = as.numeric(before_row$mem_alloc),
    after_mem_alloc = as.numeric(after_row$mem_alloc),
    stringsAsFactors = FALSE
  )
}

alloc_bytes <- function(expr, envir) {
  p <- profmem::profmem(eval(expr, envir = envir))
  sum(p$bytes[p$bytes > 0], na.rm = TRUE)
}

cat("Loading current package from:", repo_root, "\n")
pkgload::load_all(
  repo_root,
  export_all = FALSE,
  helpers = FALSE,
  quiet = TRUE,
  reset = TRUE,
  compile = FALSE
)
ns <- asNamespace("matrixCorr")

cat("Building baseline kernels/wrappers from commit:", baseline_sha, "\n")
baseline_env <- build_impl_env(
  read_path_lines = function(path) git_show_lines(baseline_sha, path, repo_root),
  ns = ns
)

cat("Building current kernels/wrappers from working tree...\n")
current_env <- build_impl_env(
  read_path_lines = function(path) current_file_lines(path, repo_root),
  ns = ns
)

set.seed(20260423)

spearman_no_na <- matrix(rnorm(7000 * 80), nrow = 7000, ncol = 80)
colnames(spearman_no_na) <- paste0("S", seq_len(ncol(spearman_no_na)))
spearman_pairwise <- spearman_no_na
spearman_pairwise[matrix(runif(length(spearman_pairwise)) < 0.08,
                         nrow = nrow(spearman_pairwise))] <- NA_real_

kendall_no_na <- matrix(rnorm(4500 * 56), nrow = 4500, ncol = 56)
colnames(kendall_no_na) <- paste0("K", seq_len(ncol(kendall_no_na)))
kendall_pairwise <- kendall_no_na
kendall_pairwise[matrix(runif(length(kendall_pairwise)) < 0.08,
                        nrow = nrow(kendall_pairwise))] <- NA_real_

cases <- list(
  list(
    label = "spearman_no_na",
    before = quote(baseline_env$spearman_rho(spearman_no_na)),
    after = quote(current_env$spearman_rho(spearman_no_na))
  ),
  list(
    label = "spearman_pairwise_na",
    before = quote(baseline_env$spearman_rho(spearman_pairwise, na_method = "pairwise")),
    after = quote(current_env$spearman_rho(spearman_pairwise, na_method = "pairwise"))
  ),
  list(
    label = "kendall_no_na",
    before = quote(baseline_env$kendall_tau(kendall_no_na)),
    after = quote(current_env$kendall_tau(kendall_no_na))
  ),
  list(
    label = "kendall_pairwise_na",
    before = quote(baseline_env$kendall_tau(kendall_pairwise, na_method = "pairwise")),
    after = quote(current_env$kendall_tau(kendall_pairwise, na_method = "pairwise"))
  )
)

results <- vector("list", length(cases))
contracts <- vector("list", length(cases))
memory_rows <- vector("list", length(cases))

for (i in seq_along(cases)) {
  cs <- cases[[i]]
  cat("\nRunning:", cs$label, "\n")

  before_out <- eval(cs$before, envir = environment())
  after_out <- eval(cs$after, envir = environment())
  contracts[[i]] <- c(list(label = cs$label), validate_contract(before_out, after_out))

  results[[i]] <- bench_one(
    label = cs$label,
    before_expr = cs$before,
    after_expr = cs$after,
    envir = environment(),
    iterations = 5L
  )

  before_alloc <- alloc_bytes(cs$before, envir = environment())
  after_alloc <- alloc_bytes(cs$after, envir = environment())
  memory_rows[[i]] <- data.frame(
    case = cs$label,
    before_alloc_bytes = before_alloc,
    after_alloc_bytes = after_alloc,
    alloc_pct_improvement = 100 * (before_alloc - after_alloc) / before_alloc,
    stringsAsFactors = FALSE
  )
}

timing_df <- do.call(rbind, results)
memory_df <- do.call(rbind, memory_rows)

cat("\n=== Runtime / bench::mark (median) ===\n")
print(timing_df, row.names = FALSE)

cat("\n=== Allocation / profmem ===\n")
print(memory_df, row.names = FALSE)

cat("\n=== Output Equivalence Checks ===\n")
for (x in contracts) {
  cat("\nCase:", x$label, "\n")
  cat("  identical():", x$identical, "\n")
  cat("  all.equal(tol=1e-12):", x$all_equal_1e12, "\n")
  chk <- x$checks
  cat("  same_dim:", chk$same_dim, "\n")
  cat("  same_dimnames:", chk$same_dimnames, "\n")
  cat("  same_storage_mode:", chk$same_storage_mode, "\n")
  cat("  same_class:", chk$same_class, "\n")
  cat("  same_diagonal:", chk$same_diagonal, "\n")
  cat("  same_symmetry:", chk$same_symmetry, "\n")
}
