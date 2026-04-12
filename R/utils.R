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

#' Inform only when verbose
#' @keywords internal
inform_if_verbose <- function(...,
                              .verbose = getOption("matrixCorr.verbose", TRUE)) {
  if (isTRUE(.verbose)) {
    cli::cli_inform(...)
  }
  invisible(NULL)
}
