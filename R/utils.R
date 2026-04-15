#' Internal validation utilities
#'
#' Centralised argument and consistency checks using cli/rlang.
#' These helpers are not exported and are intended for internal use.
#'
#' @keywords internal
#' @name matrixCorr-utils
#' @noRd
#' @importFrom cli cli_abort cli_warn cli_inform
#' @importFrom rlang is_bool is_scalar_integerish is_scalar_character arg_match
NULL

#' Abort for internal errors (should not happen)
#' @keywords internal
abort_internal <- function(message, ...,
                           .class = character()) {
  cli::cli_abort(
    message,
    ...,
    class = c(.class, "matrixCorr_error", "matrixCorr_internal_error")
  )
}

#' Abort with a standardised argument error
#' @keywords internal
abort_bad_arg <- function(arg,
                          message,
                          ...,
                          .hint = NULL,
                          .class = character()) {
  env <- rlang::env_clone(rlang::caller_env())
  rlang::env_bind(env, arg = arg)
  bullets <- c(
    "Problem with {.arg {arg}}.",
    "x" = message
  )
  if (!is.null(.hint)) {
    bullets <- c(bullets, "i" = .hint)
  }

  cli::cli_abort(
    bullets,
    arg = arg,
    ...,
    .envir = env,
    class = c(.class, "matrixCorr_error", "matrixCorr_arg_error")
  )
}

#' Check a single logical flag
#' @keywords internal
check_bool <- function(x,
                       arg = as.character(substitute(x))) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    abort_bad_arg(arg,
      message = "must be a single TRUE or FALSE."
    )
  }
  invisible(x)
}

#' Check scalar numeric (with optional bounds)
#' @keywords internal
check_scalar_numeric <- function(x,
                                 arg = as.character(substitute(x)),
                                 finite = TRUE,
                                 lower = -Inf,
                                 upper = Inf,
                                 allow_na = FALSE,
                                 closed_lower = TRUE,
                                 closed_upper = TRUE) {
  ok_len <- length(x) == 1L
  ok_num <- is.numeric(x)
  if (!ok_len || !ok_num) {
    abort_bad_arg(arg,
      message = "must be a single numeric value."
    )
  }

  if (!allow_na && (is.na(x) || is.nan(x))) {
    abort_bad_arg(arg,
      message = "must not be NA or NaN."
    )
  }

  if (finite && !(allow_na && is.na(x)) && any(!is.finite(x))) {
    abort_bad_arg(arg,
      message = "must be finite."
    )
  }

  if (!is.na(lower)) {
    if (closed_lower && x < lower) {
      abort_bad_arg(arg,
        message = "must be >= {lower}.",
        lower = lower
      )
    } else if (!closed_lower && x <= lower) {
      abort_bad_arg(arg,
        message = "must be > {lower}.",
        lower = lower
      )
    }
  }

  if (!is.na(upper)) {
    if (closed_upper && x > upper) {
      abort_bad_arg(arg,
        message = "must be <= {upper}.",
        upper = upper
      )
    } else if (!closed_upper && x >= upper) {
      abort_bad_arg(arg,
        message = "must be < {upper}.",
        upper = upper
      )
    }
  }

  invisible(x)
}

#' Check scalar character (non-empty)
#' @keywords internal
check_scalar_character <- function(x,
                                   arg = as.character(substitute(x)),
                                   allow_empty = FALSE) {
  if (!rlang::is_scalar_character(x)) {
    abort_bad_arg(arg,
      message = "must be a single character string."
    )
  }
  if (!allow_empty && !nzchar(x)) {
    abort_bad_arg(arg,
      message = "must be a non-empty character string."
    )
  }
  invisible(x)
}

#' Check that a data frame contains required columns
#' @keywords internal
check_required_cols <- function(df,
                                cols,
                                df_arg = "data") {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0L) {
    abort_bad_arg(df_arg,
      message = "is missing required column(s): {paste(missing, collapse = ', ')}.",
      .hint   = "Check spelling and ensure all required variables are present.",
      missing = missing
    )
  }
  invisible(df)
}

#' Check that two vectors have the same length
#' @keywords internal
check_same_length <- function(x, y,
                              arg_x = as.character(substitute(x)),
                              arg_y = as.character(substitute(y))) {
  if (length(x) != length(y)) {
    cli::cli_abort(
      c(
        "Length mismatch: {.arg {arg_x}} and {.arg {arg_y}} must have the same length.",
        "x" = "{.arg {arg_x}} has length {length(x)}.",
        "x" = "{.arg {arg_y}} has length {length(y)}."
      ),
      arg_x = arg_x,
      arg_y = arg_y,
      class = c("matrixCorr_error", "matrixCorr_arg_error")
    )
  }
  invisible(NULL)
}

#' Check matrix dimensions
#' @keywords internal
check_matrix_dims <- function(M,
                              nrow = NULL,
                              ncol = NULL,
                              arg = as.character(substitute(M))) {
  if (!is.matrix(M)) {
    abort_bad_arg(arg,
      message = "must be a matrix."
    )
  }
  dm <- dim(M)
  if (!is.null(nrow) && dm[1L] != nrow) {
    abort_bad_arg(arg,
      message = "must have {nrow} row{?s}, not {dm[1L]}.",
      nrow = nrow, dm = dm
    )
  }
  if (!is.null(ncol) && dm[2L] != ncol) {
    abort_bad_arg(arg,
      message = "must have {ncol} column{?s}, not {dm[2L]}.",
      ncol = ncol, dm = dm
    )
  }
  invisible(M)
}

#' Check matrix symmetry
#' @keywords internal
check_symmetric_matrix <- function(M,
                                   tol = .Machine$double.eps^0.5,
                                   arg = as.character(substitute(M))) {
  check_matrix_dims(M, arg = arg)
  if (!isSymmetric(M, tol = tol)) {
    abort_bad_arg(arg,
      message = "must be symmetric within tolerance {tol}.",
      tol = tol
    )
  }
  invisible(M)
}

#' Check a non-negative scalar (>= 0)
#' @keywords internal
check_scalar_nonneg <- function(x,
                                arg = as.character(substitute(x)),
                                strict = FALSE) {
  if (strict) {
    # Strictly positive: custom message to surface "positive"
    if (!is.numeric(x) || length(x) != 1L) {
      abort_bad_arg(arg, message = "must be a single numeric value.")
    }
    if (is.na(x) || is.nan(x)) {
      abort_bad_arg(arg, message = "must not be NA or NaN.")
    }
    if (!is.finite(x)) {
      abort_bad_arg(arg, message = "must be finite.")
    }
    if (x <= 0) {
      abort_bad_arg(arg, message = "must be a positive value (> 0).")
    }
  } else {
    check_scalar_numeric(x, arg = arg, lower = 0, closed_lower = TRUE)
  }
  invisible(x)
}

#' Check strictly positive scalar integer
#' @keywords internal
check_scalar_int_pos <- function(x,
                                 arg = as.character(substitute(x))) {
  ok <- length(x) == 1L &&
    is.numeric(x) &&
    !is.na(x) &&
    is.finite(x) &&
    x > 0 &&
    abs(x - round(x)) <= sqrt(.Machine$double.eps)
  if (!ok) {
    abort_bad_arg(arg,
      message = "must be a positive integer."
    )
  }
  as.integer(x)
}

#' Check probability in unit interval
#' @keywords internal
check_prob_scalar <- function(x,
                              arg = as.character(substitute(x)),
                              open_ends = FALSE) {
  if (open_ends) {
    check_scalar_numeric(x, arg = arg,
                         lower = 0, closed_lower = FALSE,
                         upper = 1, closed_upper = FALSE)
  } else {
    check_scalar_numeric(x, arg = arg,
                         lower = 0, closed_lower = TRUE,
                         upper = 1, closed_upper = TRUE)
  }
  invisible(x)
}

#' Check AR(1) correlation parameter
#' @keywords internal
check_ar1_rho <- function(x,
                          arg = as.character(substitute(x)),
                          bound = 0.999) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x <= -bound || x >= bound) {
    abort_bad_arg(arg,
      message = "`{arg}` must be in (-{bound}, {bound}).",
      bound = bound
    )
  }
  as.numeric(x)
}

#' Check weights vector (non-negative, finite, correct length)
#' @keywords internal
check_weights <- function(w,
                          n,
                          arg = as.character(substitute(w))) {
  if (is.null(w)) return(invisible(NULL))

  w <- as.numeric(w)
  if (length(w) != n) {
    abort_bad_arg(arg,
      message = "must have length {n}, not {length(w)}.",
      n = n, w = w
    )
  }
  if (anyNA(w) || any(!is.finite(w)) || any(w < 0)) {
    abort_bad_arg(arg,
      message = "must be non-negative, finite, and contain no missing values."
    )
  }
  invisible(w)
}

#' Check object class
#' @keywords internal
check_inherits <- function(x,
                           class,
                           arg = as.character(substitute(x))) {
  if (!inherits(x, class)) {
    abort_bad_arg(arg,
      message = "must inherit from class{?es}: {cli::qty(class)}{class*?{, }{ and }}.",
      class = class
    )
  }
  invisible(x)
}

#' Match argument to allowed values
#' @keywords internal
match_arg <- function(arg,
                      values,
                      arg_name = as.character(substitute(arg))) {
  if (is.character(arg) &&
      length(arg) == 1L &&
      !is.na(arg) &&
      identical(match(arg, values, nomatch = 0L) > 0L, TRUE)) {
    return(arg)
  }
  rlang::arg_match(arg = arg, values = values, multiple = FALSE, error_arg = arg_name)
}

#' Validate the package-wide missing-data method selector
#' @keywords internal
#' @noRd
validate_na_method <- function(na_method,
                               arg = as.character(substitute(na_method))) {
  if (is.logical(na_method) && length(na_method) == 1L && !is.na(na_method)) {
    return(if (isTRUE(na_method)) "error" else "pairwise")
  }
  if (is.character(na_method) &&
      length(na_method) == 1L &&
      !is.na(na_method) &&
      (identical(na_method, "error") || identical(na_method, "pairwise"))) {
    return(na_method)
  }
  match_arg(
    na_method,
    values = c("error", "pairwise"),
    arg_name = arg
  )
}

#' Resolve canonical/deprecated missing-data arguments
#' @keywords internal
#' @noRd
resolve_na_args <- function(na_method = "error",
                            check_na = NULL,
                            na_method_missing = FALSE,
                            arg_na_method = "na_method",
                            arg_check_na = "check_na",
                            warn = TRUE) {
  if (is.null(check_na) && isTRUE(na_method_missing)) {
    return(list(
      na_method = "error",
      check_na = TRUE
    ))
  }
  if (!is.null(check_na)) {
    check_bool(check_na, arg = arg_check_na)
    if (!isTRUE(na_method_missing)) {
      abort_bad_arg(
        arg_na_method,
        message = "and {.arg check_na} cannot both be supplied.",
        .hint = "Use only {.arg na_method}; {.arg check_na} is deprecated."
      )
    }
    if (isTRUE(warn)) {
      replacement <- if (isTRUE(check_na)) "error" else "pairwise"
      cli::cli_warn(
        c(
          "{.arg check_na} is deprecated.",
          "i" = "Use {.arg na_method} = {.val {replacement}} instead."
        ),
        replacement = replacement,
        class = c("matrixCorr_warning", "matrixCorr_deprecated_warning")
      )
    }
    resolved <- if (isTRUE(check_na)) "error" else "pairwise"
    return(list(
      na_method = resolved,
      check_na = identical(resolved, "error")
    ))
  }

  if (!isTRUE(na_method_missing) && rlang::is_bool(na_method)) {
    if (isTRUE(warn)) {
      replacement <- if (isTRUE(na_method)) "error" else "pairwise"
      cli::cli_warn(
        c(
          "Logical {.arg na_method} is deprecated compatibility input.",
          "i" = "Use {.arg na_method} = {.val {replacement}} instead."
        ),
        replacement = replacement,
        class = c("matrixCorr_warning", "matrixCorr_deprecated_warning")
      )
    }
  }

  resolved <- validate_na_method(na_method, arg = arg_na_method)
  list(
    na_method = resolved,
    check_na = identical(resolved, "error")
  )
}

#' Set OpenMP threads only when a change is required
#' @keywords internal
#' @noRd
.mc_enter_omp_threads <- function(n_threads) {
  current <- as.integer(get_omp_threads())
  target <- as.integer(n_threads)
  if (identical(current, target)) {
    return(NULL)
  }
  set_omp_threads(target)
  current
}

#' Restore OpenMP thread count after `.mc_enter_omp_threads()`
#' @keywords internal
#' @noRd
.mc_exit_omp_threads <- function(prev_threads) {
  if (!is.null(prev_threads)) {
    set_omp_threads(as.integer(prev_threads))
  }
  invisible(NULL)
}

#' Prepare OpenMP thread state for a wrapper call
#' @keywords internal
#' @noRd
.mc_prepare_omp_threads <- function(n_threads,
                                    n_threads_missing = FALSE,
                                    arg = "n_threads") {
  if (isTRUE(n_threads_missing) &&
      identical(getOption("matrixCorr.threads", 1L), 1L)) {
    return(NULL)
  }
  .mc_enter_omp_threads(check_scalar_int_pos(n_threads, arg = arg))
}

#' Extract deprecated compatibility aliases from dots
#' @keywords internal
#' @noRd
.mc_extract_legacy_aliases <- function(dots,
                                       allowed = character(),
                                       dots_arg = "...") {
  dots <- if (is.null(dots)) list() else dots
  nms <- names(dots)
  if (length(dots) && (is.null(nms) || anyNA(nms) || any(!nzchar(nms)))) {
    abort_bad_arg(
      dots_arg,
      message = "does not allow unnamed extra arguments."
    )
  }
  unknown <- setdiff(nms, allowed)
  if (length(unknown)) {
    abort_bad_arg(
      dots_arg,
      message = "contains unsupported argument(s): {paste(unknown, collapse = ', ')}."
    )
  }
  dots
}

#' Validate common correlation output arguments
#' @keywords internal
#' @noRd
.mc_validate_output_args <- function(output = c("matrix", "sparse", "edge_list"),
                                     threshold = 0,
                                     diag = TRUE,
                                     arg_output = "output",
                                     arg_threshold = "threshold",
                                     arg_diag = "diag") {
  output <- match_arg(
    output,
    values = c("matrix", "sparse", "edge_list"),
    arg_name = arg_output
  )
  check_scalar_numeric(
    threshold,
    arg = arg_threshold,
    lower = 0,
    closed_lower = TRUE
  )
  check_bool(diag, arg = arg_diag)

  threshold <- as.numeric(threshold)
  if (identical(output, "matrix") && threshold > 0) {
    abort_bad_arg(
      arg_threshold,
      message = "must be 0 when {.arg output} is {.val matrix}.",
      .hint = "Use {.arg output} = {.val sparse} or {.val edge_list} for thresholded output."
    )
  }

  list(
    output = output,
    threshold = threshold,
    diag = isTRUE(diag)
  )
}

#' Validate thresholded-output routing request
#' @keywords internal
#' @noRd
.mc_validate_thresholded_output_request <- function(output = c("matrix", "sparse", "edge_list"),
                                                    threshold = 0,
                                                    diag = TRUE,
                                                    arg_output = "output",
                                                    arg_threshold = "threshold",
                                                    arg_diag = "diag") {
  cfg <- .mc_validate_output_args(
    output = output,
    threshold = threshold,
    diag = diag,
    arg_output = arg_output,
    arg_threshold = arg_threshold,
    arg_diag = arg_diag
  )
  cfg$thresholded <- identical(cfg$output, "sparse") ||
    identical(cfg$output, "edge_list")
  cfg$thresholded <- isTRUE(cfg$thresholded) && isTRUE(cfg$threshold > 0)
  cfg
}

#' Report whether direct thresholded path is eligible
#' @keywords internal
#' @noRd
.mc_supports_direct_threshold_path <- function(method,
                                               na_method = "error",
                                               ci = FALSE,
                                               output = "matrix",
                                               threshold = 0,
                                               symmetric = TRUE,
                                               exact = TRUE,
                                               pairwise = FALSE,
                                               has_ci = FALSE,
                                               has_inference = FALSE,
                                               weighted = FALSE,
                                               ...) {
  method <- tolower(as.character(method %||% ""))
  if (!method %in% c("pearson", "spearman", "ccc", "bicor", "pbcor", "wincor")) {
    return(FALSE)
  }
  if (!identical(na_method, "error")) {
    return(FALSE)
  }
  if (!identical(output, "sparse") && !identical(output, "edge_list")) {
    return(FALSE)
  }
  if (!(is.numeric(threshold) && length(threshold) == 1L && is.finite(threshold) && threshold > 0)) {
    return(FALSE)
  }
  if (!isTRUE(symmetric) || !isTRUE(exact)) {
    return(FALSE)
  }
  if (isTRUE(ci) || isTRUE(has_ci) || isTRUE(has_inference)) {
    return(FALSE)
  }
  if (isTRUE(pairwise) || isTRUE(weighted)) {
    return(FALSE)
  }
  TRUE
}

#' First estimator class from a correlation result
#' @keywords internal
#' @noRd
.mc_corr_estimator_class <- function(x, default = "corr_estimator") {
  cls <- class(x)
  if (is.null(cls) || !length(cls)) {
    return(default)
  }
  skip <- c(
    "corr_matrix", "corr_sparse", "corr_packed_upper", "corr_edge_list",
    "corr_result", "matrix", "array", "data.frame", "list",
    "sparseMatrix", "Matrix", "mMatrix", "dMatrix", "symmetricMatrix",
    "CsparseMatrix", "dsparseMatrix", "dCsparseMatrix", "generalMatrix"
  )
  keep <- setdiff(cls, skip)
  if (!length(keep)) {
    return(default)
  }
  keep[[1L]]
}

#' Collect non-structural attributes for propagation across outputs
#' @keywords internal
#' @noRd
.mc_corr_extra_attrs <- function(x) {
  attrs <- attributes(x)
  if (is.null(attrs) || !length(attrs)) {
    return(list())
  }
  drop <- c(
    "dim", "dimnames", "class", "names", "row.names",
    "method", "description", "package", "output", "threshold", "diag",
    "diagnostics", "ci", "conf.level", "corr_result", "corr_output_class",
    "corr_estimator_class", "corr_dim", "corr_dimnames", "corr_symmetric"
  )
  attrs[setdiff(names(attrs), drop)]
}

#' Attach arbitrary attributes from a named list
#' @keywords internal
#' @noRd
.mc_set_attrs <- function(x, attrs) {
  if (is.null(attrs) || !length(attrs)) {
    return(x)
  }
  nm <- names(attrs)
  if (is.null(nm) || any(!nzchar(nm))) {
    return(x)
  }
  for (i in seq_along(attrs)) {
    attr(x, nm[[i]]) <- attrs[[i]]
  }
  x
}

#' New dense correlation result
#' @keywords internal
#' @noRd
.mc_new_corr_matrix <- function(mat,
                                estimator_class,
                                method,
                                description,
                                output = "matrix",
                                threshold = 0,
                                diag = TRUE,
                                diagnostics = NULL,
                                ci = NULL,
                                conf.level = NULL,
                                symmetric = NULL,
                                package_name = "matrixCorr",
                                extra_attrs = list(),
                                extra_classes = character()) {
  check_matrix_dims(mat, arg = "mat")
  cls <- unique(c(
    "corr_matrix",
    estimator_class,
    "corr_result",
    extra_classes,
    "matrix"
  ))
  out <- structure(
    mat,
    class = cls,
    method = method,
    description = description,
    package = package_name,
    output = output,
    threshold = as.numeric(threshold),
    diag = isTRUE(diag),
    diagnostics = diagnostics,
    ci = ci,
    conf.level = conf.level,
    corr_result = TRUE,
    corr_output_class = "corr_matrix",
    corr_estimator_class = estimator_class,
    corr_dim = as.integer(dim(mat)),
    corr_dimnames = dimnames(mat),
    corr_symmetric = if (is.null(symmetric)) {
      isTRUE(nrow(mat) == ncol(mat)) && isTRUE(isSymmetric(mat, check.attributes = FALSE))
    } else {
      isTRUE(symmetric)
    }
  )
  .mc_set_attrs(out, extra_attrs)
}

#' New sparse correlation result
#' @keywords internal
#' @noRd
.mc_new_corr_sparse <- function(x,
                                estimator_class,
                                method,
                                description,
                                output = "sparse",
                                threshold = 0,
                                diag = TRUE,
                                diagnostics = NULL,
                                ci = NULL,
                                conf.level = NULL,
                                symmetric = NULL,
                                package_name = "matrixCorr",
                                extra_attrs = list()) {
  out <- x
  attr(out, "method") <- method
  attr(out, "description") <- description
  attr(out, "package") <- package_name
  attr(out, "output") <- output
  attr(out, "threshold") <- as.numeric(threshold)
  attr(out, "diag") <- isTRUE(diag)
  attr(out, "diagnostics") <- diagnostics
  attr(out, "ci") <- ci
  attr(out, "conf.level") <- conf.level
  attr(out, "corr_result") <- TRUE
  attr(out, "corr_output_class") <- "corr_sparse"
  attr(out, "corr_estimator_class") <- estimator_class
  attr(out, "corr_dim") <- as.integer(dim(out))
  attr(out, "corr_dimnames") <- dimnames(out)
  attr(out, "corr_symmetric") <- if (is.null(symmetric)) {
    isTRUE(nrow(out) == ncol(out))
  } else {
    isTRUE(symmetric)
  }
  .mc_set_attrs(out, extra_attrs)
}

#' New packed-upper correlation result
#' @keywords internal
#' @noRd
.mc_new_corr_packed_upper <- function(df,
                                      estimator_class,
                                      method,
                                      description,
                                      threshold = 0,
                                      diag = TRUE,
                                      diagnostics = NULL,
                                      ci = NULL,
                                      conf.level = NULL,
                                      source_dim = NULL,
                                      source_dimnames = NULL,
                                      symmetric = TRUE,
                                      package_name = "matrixCorr",
                                      extra_attrs = list()) {
  out <- data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  class(out) <- c("corr_packed_upper", estimator_class, "corr_result", "data.frame")
  attr(out, "method") <- method
  attr(out, "description") <- description
  attr(out, "package") <- package_name
  attr(out, "output") <- "packed_upper"
  attr(out, "threshold") <- as.numeric(threshold)
  attr(out, "diag") <- isTRUE(diag)
  attr(out, "diagnostics") <- diagnostics
  attr(out, "ci") <- ci
  attr(out, "conf.level") <- conf.level
  attr(out, "corr_result") <- TRUE
  attr(out, "corr_output_class") <- "corr_packed_upper"
  attr(out, "corr_estimator_class") <- estimator_class
  attr(out, "corr_dim") <- as.integer(source_dim %||% c(NA_integer_, NA_integer_))
  attr(out, "corr_dimnames") <- source_dimnames
  attr(out, "corr_symmetric") <- isTRUE(symmetric)
  .mc_set_attrs(out, extra_attrs)
}

#' New edge-list correlation result
#' @keywords internal
#' @noRd
.mc_new_corr_edge_list <- function(df,
                                   estimator_class,
                                   method,
                                   description,
                                   threshold = 0,
                                   diag = TRUE,
                                   diagnostics = NULL,
                                   ci = NULL,
                                   conf.level = NULL,
                                   source_dim = NULL,
                                   source_dimnames = NULL,
                                   symmetric = NULL,
                                   package_name = "matrixCorr",
                                   extra_attrs = list()) {
  out <- data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  class(out) <- c("corr_edge_list", estimator_class, "corr_result", "data.frame")
  attr(out, "method") <- method
  attr(out, "description") <- description
  attr(out, "package") <- package_name
  attr(out, "output") <- "edge_list"
  attr(out, "threshold") <- as.numeric(threshold)
  attr(out, "diag") <- isTRUE(diag)
  attr(out, "diagnostics") <- diagnostics
  attr(out, "ci") <- ci
  attr(out, "conf.level") <- conf.level
  attr(out, "corr_result") <- TRUE
  attr(out, "corr_output_class") <- "corr_edge_list"
  attr(out, "corr_estimator_class") <- estimator_class
  attr(out, "corr_dim") <- as.integer(source_dim %||% c(NA_integer_, NA_integer_))
  attr(out, "corr_dimnames") <- source_dimnames
  attr(out, "corr_symmetric") <- if (is.null(symmetric)) NA else isTRUE(symmetric)
  .mc_set_attrs(out, extra_attrs)
}

#' Convert C++ triplets to sparse matrix
#' @keywords internal
#' @noRd
.mc_triplets_to_sparse <- function(triplets,
                                   dim,
                                   dimnames = NULL,
                                   symmetric = TRUE) {
  check_same_length(triplets$i, triplets$j, arg_x = "triplets$i", arg_y = "triplets$j")
  check_same_length(triplets$i, triplets$x, arg_x = "triplets$i", arg_y = "triplets$x")
  dm <- as.integer(dim)
  if (length(dm) != 2L) {
    abort_bad_arg("dim", message = "must be length 2.")
  }

  i <- as.integer(triplets$i)
  j <- as.integer(triplets$j)
  x <- as.numeric(triplets$x)

  keep <- is.finite(i) & is.finite(j) & is.finite(x) &
    i >= 1L & i <= dm[[1L]] &
    j >= 1L & j <= dm[[2L]]
  i <- i[keep]
  j <- j[keep]
  x <- x[keep]

  out <- Matrix::sparseMatrix(
    i = i,
    j = j,
    x = x,
    dims = dm,
    dimnames = dimnames
  )
  if (isTRUE(symmetric)) {
    return(Matrix::forceSymmetric(out, uplo = "U"))
  }
  out
}

#' Convert C++ triplets to internal edge-list payload
#' @keywords internal
#' @noRd
.mc_triplets_to_edge_list <- function(triplets) {
  check_same_length(triplets$i, triplets$j, arg_x = "triplets$i", arg_y = "triplets$j")
  check_same_length(triplets$i, triplets$x, arg_x = "triplets$i", arg_y = "triplets$x")
  data.frame(
    i = as.integer(triplets$i),
    j = as.integer(triplets$j),
    x = as.numeric(triplets$x),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Build sparse/edge output directly from retained triplets
#' @keywords internal
#' @noRd
.mc_finalize_triplets_output <- function(triplets,
                                         output = c("sparse", "edge_list"),
                                         estimator_class,
                                         method,
                                         description,
                                         threshold,
                                         diag,
                                         source_dim,
                                         source_dimnames = NULL,
                                         diagnostics = NULL,
                                         ci = NULL,
                                         conf.level = NULL,
                                         symmetric = TRUE,
                                         package_name = "matrixCorr",
                                         extra_attrs = list()) {
  output <- match.arg(output)
  source_dim <- as.integer(source_dim)

  meta <- list(
    source_class = estimator_class,
    method = method,
    description = description,
    package = package_name,
    diagnostics = diagnostics
  )

  if (identical(output, "sparse")) {
    out <- .mc_new_corr_sparse(
      x = .mc_triplets_to_sparse(
        triplets,
        dim = source_dim,
        dimnames = source_dimnames,
        symmetric = symmetric
      ),
      estimator_class = estimator_class,
      method = method,
      description = description,
      output = "sparse",
      threshold = threshold,
      diag = diag,
      diagnostics = diagnostics,
      ci = ci,
      conf.level = conf.level,
      symmetric = symmetric,
      package_name = package_name,
      extra_attrs = extra_attrs
    )
    attr(out, "matrixCorr_meta") <- meta
    return(out)
  }

  out <- .mc_new_corr_edge_list(
    df = .mc_triplets_to_edge_list(triplets),
    estimator_class = estimator_class,
    method = method,
    description = description,
    threshold = threshold,
    diag = diag,
    diagnostics = diagnostics,
    ci = ci,
    conf.level = conf.level,
    source_dim = source_dim,
    source_dimnames = source_dimnames,
    symmetric = symmetric,
    package_name = package_name,
    extra_attrs = extra_attrs
  )
  attr(out, "matrixCorr_meta") <- meta
  out
}

#' Convert a symmetric matrix to an upper-triangle edge frame
#' @keywords internal
#' @noRd
.mc_matrix_to_edge_list <- function(mat,
                                    threshold = 0,
                                    diag = TRUE) {
  check_matrix_dims(mat, arg = "mat")
  check_scalar_numeric(threshold, arg = "threshold", lower = 0, closed_lower = TRUE)
  check_bool(diag, arg = "diag")

  if (nrow(mat) == 0L || ncol(mat) == 0L) {
    return(data.frame(
      row = character(),
      col = character(),
      value = numeric(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  same_dimnames <- isTRUE(identical(rownames(mat), colnames(mat))) ||
    (is.null(rownames(mat)) && is.null(colnames(mat)))
  use_upper <- isTRUE(nrow(mat) == ncol(mat)) &&
    same_dimnames &&
    isTRUE(isSymmetric(mat, check.attributes = FALSE))
  if (use_upper) {
    idx <- upper.tri(mat, diag = isTRUE(diag))
  } else {
    idx <- matrix(TRUE, nrow = nrow(mat), ncol = ncol(mat))
    if (!isTRUE(diag)) {
      d <- seq_len(min(nrow(mat), ncol(mat)))
      idx[cbind(d, d)] <- FALSE
    }
  }
  i <- row(mat)[idx]
  j <- col(mat)[idx]
  vals <- mat[idx]
  keep <- is.finite(vals) & (abs(vals) >= as.numeric(threshold))

  rn <- rownames(mat)
  cn <- colnames(mat)
  row_out <- if (is.null(rn)) as.character(i[keep]) else rn[i[keep]]
  col_out <- if (is.null(cn)) as.character(j[keep]) else cn[j[keep]]

  data.frame(
    row = row_out,
    col = col_out,
    value = as.numeric(vals[keep]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Convert a symmetric matrix to packed upper-triangle records
#' @keywords internal
#' @noRd
.mc_matrix_to_packed_upper <- function(mat,
                                       threshold = 0,
                                       diag = TRUE) {
  .mc_matrix_to_edge_list(mat, threshold = threshold, diag = diag)
}

#' Convert a symmetric matrix to thresholded sparse matrix
#' @keywords internal
#' @noRd
.mc_matrix_to_sparse_thresholded <- function(mat,
                                             threshold = 0,
                                             diag = TRUE) {
  check_matrix_dims(mat, arg = "mat")
  check_scalar_numeric(threshold, arg = "threshold", lower = 0, closed_lower = TRUE)
  check_bool(diag, arg = "diag")

  dn <- dimnames(mat)
  same_dimnames <- isTRUE(identical(rownames(mat), colnames(mat))) ||
    (is.null(rownames(mat)) && is.null(colnames(mat)))
  use_upper <- isTRUE(nrow(mat) == ncol(mat)) &&
    same_dimnames &&
    isTRUE(isSymmetric(mat, check.attributes = FALSE))
  if (use_upper) {
    idx <- upper.tri(mat, diag = isTRUE(diag))
  } else {
    idx <- matrix(TRUE, nrow = nrow(mat), ncol = ncol(mat))
    if (!isTRUE(diag)) {
      d <- seq_len(min(nrow(mat), ncol(mat)))
      idx[cbind(d, d)] <- FALSE
    }
  }
  i <- row(mat)[idx]
  j <- col(mat)[idx]
  vals <- mat[idx]
  keep <- is.finite(vals) & (abs(vals) >= as.numeric(threshold))

  out <- Matrix::sparseMatrix(
    i = as.integer(i[keep]),
    j = as.integer(j[keep]),
    x = as.numeric(vals[keep]),
    dims = dim(mat),
    dimnames = dn
  )
  if (use_upper) {
    return(Matrix::forceSymmetric(out, uplo = "U"))
  }
  out
}

#' Route a computed correlation matrix object to requested output shape
#' @keywords internal
#' @noRd
.mc_finalize_corr_output <- function(x,
                                     output = c("matrix", "sparse", "edge_list"),
                                     threshold = 0,
                                     diag = TRUE) {
  cfg <- .mc_validate_output_args(output = output, threshold = threshold, diag = diag)
  mat <- as.matrix(x)
  estimator_class <- .mc_corr_estimator_class(x)
  method <- attr(x, "method", exact = TRUE) %||% estimator_class
  description <- attr(x, "description", exact = TRUE) %||% "Correlation result"
  diagnostics <- attr(x, "diagnostics", exact = TRUE)
  ci <- attr(x, "ci", exact = TRUE)
  conf.level <- attr(x, "conf.level", exact = TRUE)
  package_name <- attr(x, "package", exact = TRUE) %||% "matrixCorr"
  extra_attrs <- .mc_corr_extra_attrs(x)
  meta <- list(
    source_class = class(x),
    method = method,
    description = description,
    package = package_name,
    diagnostics = diagnostics
  )
  source_dim <- dim(mat)
  source_dimnames <- dimnames(mat)
  source_symmetric <- isTRUE(nrow(mat) == ncol(mat)) &&
    (isTRUE(identical(rownames(mat), colnames(mat))) ||
      (is.null(rownames(mat)) && is.null(colnames(mat)))) &&
    isTRUE(isSymmetric(mat, check.attributes = FALSE))

  if (identical(cfg$output, "matrix")) {
    if (inherits(x, "corr_matrix")) {
      attr(x, "output") <- "matrix"
      attr(x, "threshold") <- 0
      attr(x, "diag") <- TRUE
      return(x)
    }
    return(.mc_new_corr_matrix(
      mat = mat,
      estimator_class = estimator_class,
      method = method,
      description = description,
      output = "matrix",
      threshold = 0,
      diag = TRUE,
      diagnostics = diagnostics,
      ci = ci,
      conf.level = conf.level,
      symmetric = source_symmetric,
      package_name = package_name,
      extra_attrs = extra_attrs,
      extra_classes = setdiff(class(x), c(estimator_class, "matrix"))
    ))
  }

  if (identical(cfg$output, "packed_upper")) {
    out <- .mc_new_corr_packed_upper(
      df = .mc_matrix_to_packed_upper(mat, threshold = cfg$threshold, diag = cfg$diag),
      estimator_class = estimator_class,
      method = method,
      description = description,
      threshold = cfg$threshold,
      diag = cfg$diag,
      diagnostics = diagnostics,
      ci = ci,
      conf.level = conf.level,
      source_dim = source_dim,
      source_dimnames = source_dimnames,
      symmetric = source_symmetric,
      package_name = package_name,
      extra_attrs = extra_attrs
    )
    attr(out, "matrixCorr_meta") <- meta
    return(out)
  }

  if (identical(cfg$output, "edge_list")) {
    out <- .mc_new_corr_edge_list(
      df = .mc_matrix_to_edge_list(mat, threshold = cfg$threshold, diag = cfg$diag),
      estimator_class = estimator_class,
      method = method,
      description = description,
      threshold = cfg$threshold,
      diag = cfg$diag,
      diagnostics = diagnostics,
      ci = ci,
      conf.level = conf.level,
      source_dim = source_dim,
      source_dimnames = source_dimnames,
      symmetric = source_symmetric,
      package_name = package_name,
      extra_attrs = extra_attrs
    )
    attr(out, "matrixCorr_meta") <- meta
    return(out)
  }

  out <- .mc_new_corr_sparse(
    x = .mc_matrix_to_sparse_thresholded(mat, threshold = cfg$threshold, diag = cfg$diag),
    estimator_class = estimator_class,
    method = method,
    description = description,
    output = "sparse",
    threshold = cfg$threshold,
    diag = cfg$diag,
    diagnostics = diagnostics,
    ci = ci,
    conf.level = conf.level,
    symmetric = source_symmetric,
    package_name = package_name,
    extra_attrs = extra_attrs
  )
  attr(out, "matrixCorr_meta") <- meta
  out
}

#' Correlation result dimensions
#' @keywords internal
#' @noRd
.mc_corr_dim <- function(x) {
  dm <- attr(x, "corr_dim", exact = TRUE)
  if (!is.null(dm) && length(dm) == 2L) {
    return(as.integer(dm))
  }
  as.integer(dim(x))
}

#' Correlation result dimnames
#' @keywords internal
#' @noRd
.mc_corr_dimnames <- function(x) {
  dn <- attr(x, "corr_dimnames", exact = TRUE)
  if (is.list(dn) && length(dn) == 2L) {
    return(dn)
  }
  dimnames(x)
}

#' Coerce correlation object to edge data frame
#' @keywords internal
#' @noRd
.mc_corr_as_edge_df <- function(x) {
  output <- attr(x, "output", exact = TRUE) %||% "matrix"

  if (inherits(x, "corr_packed_upper") || inherits(x, "corr_edge_list")) {
    out <- x
    class(out) <- "data.frame"
    out <- as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
    if (all(c("row", "col", "value") %in% names(out))) {
      return(out[, c("row", "col", "value"), drop = FALSE])
    }
    if (all(c("i", "j", "x") %in% names(out))) {
      dn <- .mc_corr_dimnames(x)
      ii <- as.integer(out$i)
      jj <- as.integer(out$j)
      vv <- as.numeric(out$x)
      row_lab <- if (is.null(dn[[1L]])) as.character(ii) else dn[[1L]][ii]
      col_lab <- if (is.null(dn[[2L]])) as.character(jj) else dn[[2L]][jj]
      return(data.frame(
        row = row_lab,
        col = col_lab,
        value = vv,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ))
    }
    abort_bad_arg("x", "must contain columns `row`, `col`, and `value` (or internal `i`, `j`, `x`).")
  }

  if (identical(output, "sparse") || inherits(x, "sparseMatrix")) {
    sx <- Matrix::summary(x)
    ii <- sx$i
    jj <- sx$j
    vv <- sx$x
    dn <- .mc_corr_dimnames(x)
    row_lab <- if (is.null(dn[[1L]])) as.character(ii) else dn[[1L]][ii]
    col_lab <- if (is.null(dn[[2L]])) as.character(jj) else dn[[2L]][jj]
    return(data.frame(
      row = row_lab,
      col = col_lab,
      value = as.numeric(vv),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }

  .mc_matrix_to_edge_list(
    as.matrix(x),
    threshold = if (identical(output, "matrix")) 0 else (attr(x, "threshold", exact = TRUE) %||% 0),
    diag = if (identical(output, "matrix")) TRUE else (attr(x, "diag", exact = TRUE) %||% TRUE)
  )
}

#' Coerce correlation object to dense matrix
#' @keywords internal
#' @noRd
.mc_corr_as_dense_matrix <- function(x) {
  if (is.matrix(x) || inherits(x, "corr_matrix")) {
    return(as.matrix(x))
  }
  if (identical(attr(x, "output", exact = TRUE), "sparse") || inherits(x, "sparseMatrix")) {
    return(as.matrix(x))
  }

  dm <- .mc_corr_dim(x)
  if (length(dm) != 2L || any(!is.finite(dm))) {
    abort_bad_arg("x", "does not expose matrix dimensions.")
  }
  dn <- .mc_corr_dimnames(x)
  mat <- matrix(NA_real_, nrow = dm[[1L]], ncol = dm[[2L]], dimnames = dn)
  edf <- .mc_corr_as_edge_df(x)
  if (!nrow(edf)) {
    return(mat)
  }

  rmap <- dn[[1L]]
  cmap <- dn[[2L]]
  ii <- if (is.null(rmap)) as.integer(edf$row) else match(edf$row, rmap)
  jj <- if (is.null(cmap)) as.integer(edf$col) else match(edf$col, cmap)
  keep <- is.finite(ii) & is.finite(jj)
  ii <- ii[keep]
  jj <- jj[keep]
  vv <- as.numeric(edf$value[keep])
  mat[cbind(ii, jj)] <- vv

  symmetric <- isTRUE(attr(x, "corr_symmetric", exact = TRUE))
  if (symmetric) {
    mat[cbind(jj, ii)] <- vv
  }
  mat
}

#' Correlation pairs abstraction
#' @keywords internal
#' @noRd
.mc_corr_pairs <- function(x, diag = FALSE) {
  check_bool(diag, arg = "diag")
  out <- .mc_corr_as_edge_df(x)
  if (!isTRUE(diag) && nrow(out)) {
    out <- out[out$row != out$col, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

#' Correlation overview abstraction
#' @keywords internal
#' @noRd
.mc_corr_overview <- function(x) {
  dm <- .mc_corr_dim(x)
  edges <- .mc_corr_pairs(x, diag = TRUE)
  method <- attr(x, "method", exact = TRUE) %||% .mc_corr_estimator_class(x)
  output <- attr(x, "output", exact = TRUE) %||% "matrix"
  thr <- attr(x, "threshold", exact = TRUE) %||% 0
  use_diag <- attr(x, "diag", exact = TRUE) %||% TRUE
  symmetric <- isTRUE(attr(x, "corr_symmetric", exact = TRUE))
  has_ci <- !is.null(attr(x, "ci", exact = TRUE))
  inf_attr <- attr(x, "inference", exact = TRUE)
  has_p <- is.list(inf_attr) && is.matrix(inf_attr$p_value)
  conf.level <- attr(x, "conf.level", exact = TRUE) %||%
    (attr(x, "ci", exact = TRUE)$conf.level %||% NA_real_)
  ci_method <- attr(x, "ci", exact = TRUE)$ci.method %||%
    attr(x, "ci.method", exact = TRUE) %||%
    NA_character_
  inference_method <- if (is.list(inf_attr) && is.character(inf_attr$method) &&
                          length(inf_attr$method) == 1L) {
    inf_attr$method
  } else {
    NA_character_
  }

  vals <- as.numeric(edges$value)
  vals <- vals[is.finite(vals)]
  m <- .mc_corr_as_dense_matrix(x)
  selector <- if (symmetric && dm[[1L]] == dm[[2L]]) {
    upper.tri(m, diag = FALSE)
  } else {
    matrix(TRUE, nrow = dm[[1L]], ncol = dm[[2L]])
  }
  pick_range <- function(mat) {
    if (!is.matrix(mat) || !identical(dim(mat), dm)) {
      return(c(NA_real_, NA_real_))
    }
    vv <- as.numeric(mat[selector])
    vv <- vv[is.finite(vv)]
    if (!length(vv)) {
      return(c(NA_real_, NA_real_))
    }
    c(min(vv), max(vv))
  }
  diag_attr <- attr(x, "diagnostics", exact = TRUE)
  skipped_n <- if (is.list(diag_attr)) pick_range(diag_attr$skipped_n) else c(NA_real_, NA_real_)
  skipped_prop <- if (is.list(diag_attr)) pick_range(diag_attr$skipped_prop) else c(NA_real_, NA_real_)
  n_missing <- if (is.matrix(m)) sum(is.na(m)) else NA_integer_
  threshold_sets <- {
    thr <- attr(x, "thresholds", exact = TRUE)
    if (is.list(thr)) length(thr) else NA_integer_
  }
  thr_num <- if (is.list(thr)) numeric() else suppressWarnings(as.numeric(thr))
  if (!length(thr_num) || !is.finite(thr_num[[1L]])) {
    thr_num <- 0
  } else {
    thr_num <- thr_num[[1L]]
  }
  conf_num <- suppressWarnings(as.numeric(conf.level))
  if (!length(conf_num) || !is.finite(conf_num[[1L]])) {
    conf_num <- NA_real_
  } else {
    conf_num <- conf_num[[1L]]
  }
  if (length(ci_method) != 1L || is.list(ci_method)) {
    ci_method <- NA_character_
  }
  if (length(inference_method) != 1L || is.list(inference_method)) {
    inference_method <- NA_character_
  }

  list(
    class = .mc_corr_estimator_class(x),
    method = method,
    output = output,
    threshold = thr_num,
    diag = isTRUE(use_diag),
    n_variables = if (symmetric) as.integer(dm[[1L]]) else NA_integer_,
    n_rows = as.integer(dm[[1L]]),
    n_cols = as.integer(dm[[2L]]),
    n_pairs = nrow(edges),
    n_pairs_retained = nrow(edges),
    symmetric = symmetric,
    n_missing = n_missing,
    threshold_sets = threshold_sets,
    estimate_min = if (length(vals)) min(vals) else NA_real_,
    estimate_max = if (length(vals)) max(vals) else NA_real_,
    has_ci = has_ci,
    has_p = has_p,
    conf.level = conf_num,
    ci_method = ci_method,
    inference_method = inference_method,
    skipped_n_min = skipped_n[[1L]],
    skipped_n_max = skipped_n[[2L]],
    skipped_prop_min = skipped_prop[[1L]],
    skipped_prop_max = skipped_prop[[2L]]
  )
}

#' Inform only when verbose
#' @keywords internal
inform_if_verbose <- function(...,
                              .verbose = getOption("matrixCorr.verbose", TRUE)) {
  if (isTRUE(.verbose)) {
    cli::cli_inform(...)
  }
  invisible(NULL)
}


