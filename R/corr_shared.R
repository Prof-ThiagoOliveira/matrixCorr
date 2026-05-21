# Shared correlation-result helpers.
#
# These helpers are used by several correlation and agreement estimators. Keep
# them outside estimator-specific files so result construction does not depend
# on latent-correlation implementation details.

.mc_square_dimnames <- function(names) {
  list(names, names)
}

.mc_set_matrix_dimnames <- function(x, row_names = NULL, col_names = row_names) {
  if (is.null(row_names) && is.null(col_names)) {
    return(x)
  }
  dimnames(x) <- list(row_names, col_names)
  x
}

.mc_resolve_corr_na <- function(na_method,
                                dots = list(),
                                na_method_missing = FALSE,
                                allowed = c("error", "pairwise", "complete"),
                                warn = TRUE) {
  if (!length(dots) && isTRUE(na_method_missing)) {
    return(list(na_method = "error", check_na = TRUE))
  }

  legacy_args <- .mc_extract_legacy_aliases(dots, allowed = "check_na")
  resolve_na_args(
    na_method = na_method,
    check_na = legacy_args$check_na %||% NULL,
    na_method_missing = na_method_missing,
    allowed = allowed,
    warn = warn
  )
}

.mc_prepare_corr_input <- function(data,
                                   na_cfg,
                                   min_n = 2L,
                                   arg = "data") {
  numeric_data <- validate_corr_input(data, check_na = na_cfg$check_na)
  diagnostics <- NULL

  if (identical(na_cfg$na_method, "complete")) {
    cc <- .mc_complete_case_matrix(numeric_data, min_n = min_n, arg = arg)
    numeric_data <- cc$data
    diagnostics <- cc$diagnostics
  }

  col_names <- colnames(numeric_data)
  list(
    data = numeric_data,
    colnames = col_names,
    dimnames = if (is.null(col_names)) NULL else .mc_square_dimnames(col_names),
    diagnostics = diagnostics
  )
}

.mc_with_omp_threads <- function(n_threads,
                                 n_threads_missing = FALSE,
                                 code) {
  prev_threads <- .mc_prepare_omp_threads(
    n_threads,
    n_threads_missing = n_threads_missing
  )
  if (!is.null(prev_threads)) {
    on.exit(.mc_exit_omp_threads(prev_threads), add = TRUE)
  }
  force(code)
}

.mc_prepare_corr_output <- function(output = c("matrix", "sparse", "edge_list"),
                                    threshold = 0,
                                    diag = TRUE,
                                    thresholded = FALSE) {
  validate <- if (isTRUE(thresholded)) {
    .mc_validate_thresholded_output_request
  } else {
    .mc_validate_output_args
  }
  validate(output = output, threshold = threshold, diag = diag)
}

.mc_structure_corr_matrix <- function(mat, class_name, method, description,
                                      diagnostics = NULL, thresholds = NULL,
                                      correct = NULL, dimnames = NULL,
                                      symmetric = NULL,
                                      extra_attrs = NULL,
                                      classes = c(class_name, "matrix")) {
  if (!is.null(dimnames)) {
    dimnames(mat) <- dimnames
  }
  symmetric_flag <- symmetric
  if (is.null(symmetric_flag)) {
    symmetric_flag <- isTRUE(nrow(mat) == ncol(mat)) &&
      isTRUE(isSymmetric(mat, check.attributes = FALSE))
  }
  keep_classes <- setdiff(
    classes,
    c(class_name, "matrix", "corr_matrix", "corr_result")
  )
  .mc_new_corr_matrix(
    mat = mat,
    estimator_class = class_name,
    method = method,
    description = description,
    output = "matrix",
    threshold = 0,
    diag = TRUE,
    diagnostics = diagnostics,
    ci = attr(mat, "ci", exact = TRUE),
    conf.level = attr(mat, "conf.level", exact = TRUE),
    symmetric = symmetric_flag,
    extra_attrs = c(
      list(
        thresholds = thresholds,
        correct = correct
      ),
      extra_attrs %||% list()
    ),
    extra_classes = keep_classes
  )
}

.mc_finalize_corr_result <- function(mat,
                                     class_name,
                                     method,
                                     description,
                                     output_cfg,
                                     diagnostics = NULL,
                                     dimnames = NULL,
                                     symmetric = NULL,
                                     extra_attrs = NULL,
                                     classes = c(class_name, "matrix")) {
  out <- .mc_structure_corr_matrix(
    mat = mat,
    class_name = class_name,
    method = method,
    description = description,
    diagnostics = diagnostics,
    dimnames = dimnames,
    symmetric = symmetric,
    extra_attrs = extra_attrs,
    classes = classes
  )

  .mc_finalize_corr_output_fast(
    out,
    output = output_cfg$output,
    threshold = output_cfg$threshold,
    diag = output_cfg$diag
  )
}
