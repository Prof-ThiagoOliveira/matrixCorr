#!/usr/bin/env Rscript

baseline_sha <- "9c47c4fcff7ff22328213077ce7eab3684c95566"

required_pkgs <- c("Rcpp", "pkgload")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Install required packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

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

build_baseline_env <- function(sha, repo) {
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
    writeLines(git_show_lines(sha, f, repo), file.path(tmp, basename(f)), useBytes = TRUE)
  }

  cpp_env <- new.env(parent = baseenv())
  Rcpp::sourceCpp(file.path(tmp, "spearman.cpp"), env = cpp_env, rebuild = TRUE, showOutput = FALSE)
  Rcpp::sourceCpp(file.path(tmp, "kendall_corr.cpp"), env = cpp_env, rebuild = TRUE, showOutput = FALSE)

  ns <- asNamespace("matrixCorr")
  wrappers_env <- new.env(parent = ns)
  eval(parse(text = paste(git_show_lines(sha, "R/spearman_rho.R", repo), collapse = "\n")),
       envir = wrappers_env)
  eval(parse(text = paste(git_show_lines(sha, "R/kendall_corr.R", repo), collapse = "\n")),
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

check_result_pair <- function(name, before, after) {
  if (is.atomic(before) && length(before) == 1L && is.null(dim(before))) {
    same_scalar <- identical(before, after) ||
      isTRUE(all.equal(before, after, tolerance = 1e-12, check.attributes = TRUE))
    if (!same_scalar) {
      stop("[", name, "] scalar mismatch: before=", before, " after=", after, call. = FALSE)
    }
    cat("[PASS]", name, "(scalar)\n")
    return(invisible(NULL))
  }

  bm <- as.matrix(before)
  am <- as.matrix(after)

  hard_checks <- list(
    same_dim = identical(dim(bm), dim(am)),
    same_dimnames = identical(dimnames(bm), dimnames(am)),
    same_storage_mode = identical(storage.mode(bm), storage.mode(am)),
    same_class = identical(class(before), class(after)),
    same_diag = identical(diag(bm), diag(am)),
    same_symmetry = identical(isSymmetric(bm), isSymmetric(am))
  )
  failed <- names(hard_checks)[!unlist(hard_checks, use.names = FALSE)]
  if (length(failed)) {
    stop("[", name, "] structural mismatch: ", paste(failed, collapse = ", "), call. = FALSE)
  }

  if (identical(before, after)) {
    cat("[PASS]", name, "(identical)\n")
    return(invisible(NULL))
  }

  if (!isTRUE(all.equal(before, after, tolerance = 1e-12, check.attributes = TRUE))) {
    stop("[", name, "] value mismatch (not identical, not all.equal tol=1e-12).", call. = FALSE)
  }

  cat("[PASS]", name, "(all.equal tol=1e-12)\n")
}

cat("Loading current package...\n")
pkgload::load_all(repo_root, export_all = FALSE, helpers = FALSE, quiet = TRUE, reset = TRUE)
cat("Loading baseline implementation from commit", baseline_sha, "...\n")
baseline_env <- build_baseline_env(baseline_sha, repo_root)

set.seed(20260423)

# Spearman data
sp_no_na <- matrix(rnorm(500 * 14), nrow = 500, ncol = 14)
colnames(sp_no_na) <- paste0("SP", seq_len(ncol(sp_no_na)))
sp_ties <- matrix(sample(1:25, 420 * 12, replace = TRUE), nrow = 420, ncol = 12)
colnames(sp_ties) <- paste0("SPT", seq_len(ncol(sp_ties)))
sp_na <- sp_no_na
sp_na[matrix(runif(length(sp_na)) < 0.12, nrow = nrow(sp_na))] <- NA_real_

# Kendall data
kd_no_na <- matrix(rnorm(420 * 12), nrow = 420, ncol = 12)
colnames(kd_no_na) <- paste0("KD", seq_len(ncol(kd_no_na)))
kd_ties <- matrix(sample(1:18, 420 * 10, replace = TRUE), nrow = 420, ncol = 10)
colnames(kd_ties) <- paste0("KDT", seq_len(ncol(kd_ties)))
kd_na <- kd_no_na
kd_na[matrix(runif(length(kd_na)) < 0.12, nrow = nrow(kd_na))] <- NA_real_

checks <- list(
  list(
    name = "spearman_no_na_matrix",
    before = quote(baseline_env$spearman_rho(sp_no_na)),
    after = quote(spearman_rho(sp_no_na))
  ),
  list(
    name = "spearman_ties_matrix",
    before = quote(baseline_env$spearman_rho(sp_ties)),
    after = quote(spearman_rho(sp_ties))
  ),
  list(
    name = "spearman_pairwise_na_matrix",
    before = quote(baseline_env$spearman_rho(sp_na, na_method = "pairwise")),
    after = quote(spearman_rho(sp_na, na_method = "pairwise"))
  ),
  list(
    name = "spearman_ci_pairwise",
    before = quote(baseline_env$spearman_rho(sp_na[, 1:6], na_method = "pairwise", ci = TRUE)),
    after = quote(spearman_rho(sp_na[, 1:6], na_method = "pairwise", ci = TRUE))
  ),
  list(
    name = "spearman_sparse_threshold",
    before = quote(baseline_env$spearman_rho(sp_no_na, output = "sparse", threshold = 0.2, diag = FALSE)),
    after = quote(spearman_rho(sp_no_na, output = "sparse", threshold = 0.2, diag = FALSE))
  ),
  list(
    name = "kendall_no_na_matrix",
    before = quote(baseline_env$kendall_tau(kd_no_na)),
    after = quote(kendall_tau(kd_no_na))
  ),
  list(
    name = "kendall_ties_matrix",
    before = quote(baseline_env$kendall_tau(kd_ties)),
    after = quote(kendall_tau(kd_ties))
  ),
  list(
    name = "kendall_pairwise_na_matrix",
    before = quote(baseline_env$kendall_tau(kd_na, na_method = "pairwise")),
    after = quote(kendall_tau(kd_na, na_method = "pairwise"))
  ),
  list(
    name = "kendall_ci_fieller_pairwise",
    before = quote(baseline_env$kendall_tau(kd_na[, 1:6], na_method = "pairwise", ci = TRUE, ci_method = "fieller")),
    after = quote(kendall_tau(kd_na[, 1:6], na_method = "pairwise", ci = TRUE, ci_method = "fieller"))
  ),
  list(
    name = "kendall_ci_brown_benedetti_pairwise",
    before = quote(baseline_env$kendall_tau(kd_na[, 1:5], na_method = "pairwise", ci = TRUE, ci_method = "brown_benedetti")),
    after = quote(kendall_tau(kd_na[, 1:5], na_method = "pairwise", ci = TRUE, ci_method = "brown_benedetti"))
  ),
  list(
    name = "kendall_ci_ifel_pairwise",
    before = quote(baseline_env$kendall_tau(kd_na[, 1:5], na_method = "pairwise", ci = TRUE, ci_method = "if_el")),
    after = quote(kendall_tau(kd_na[, 1:5], na_method = "pairwise", ci = TRUE, ci_method = "if_el"))
  ),
  list(
    name = "kendall_two_vector_mode",
    before = quote(baseline_env$kendall_tau(kd_no_na[, 1], kd_no_na[, 2])),
    after = quote(kendall_tau(kd_no_na[, 1], kd_no_na[, 2]))
  )
)

for (cs in checks) {
  before_out <- eval(cs$before, envir = environment())
  after_out <- eval(cs$after, envir = environment())
  check_result_pair(cs$name, before_out, after_out)
}

cat("\nAll equivalence checks passed.\n")
