# Correlation output validation, constructors, conversion, and pair views.

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
  if (!method %in% c("pearson", "spearman", "ccc", "bicor", "pbcor", "wincor", "cohen_kappa", "weighted_kappa")) {
    return(FALSE)
  }
  if (identical(method, "cohen_kappa")) {
    if (!na_method %in% c("error", "complete")) {
      return(FALSE)
    }
  } else if (identical(method, "weighted_kappa")) {
    if (!identical(na_method, "error")) {
      return(FALSE)
    }
  } else if (!identical(na_method, "error")) {
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
    "corr_matrix", "corr_sparse", "corr_edge_list",
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
  pos <- attr(df, "corr_pair_pos", exact = TRUE)
  if (!is.null(pos) && !is.null(rownames(pos)) && !is.null(rownames(df))) {
    idx <- match(rownames(df), rownames(pos))
    if (anyNA(idx) || length(idx) != nrow(df)) {
      abort_internal("positional pair mapping does not align with edge rows.")
    }
    pos <- pos[idx, , drop = FALSE]
  }
  out <- data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  native_i <- NULL
  native_j <- NULL
  native_rn <- NULL
  native_cn <- NULL
  if (!all(c("row", "col", "value") %in% names(out)) &&
      all(c("i", "j", "x") %in% names(out))) {
    dn <- source_dimnames
    loc <- .mc_validated_native_positions(out$i, out$j, source_dim)
    ii <- loc$i
    jj <- loc$j
    rn <- if (is.list(dn) && length(dn) == 2L) dn[[1L]] else NULL
    cn <- if (is.list(dn) && length(dn) == 2L) dn[[2L]] else NULL
    native_i <- ii
    native_j <- jj
    native_rn <- rn
    native_cn <- cn
    pos <- NULL
    out <- data.frame(
      row = if (is.null(rn)) as.character(ii) else rn[ii],
      col = if (is.null(cn)) as.character(jj) else cn[jj],
      value = as.numeric(out$x),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  if (!all(c("row", "col", "value") %in% names(out))) {
    abort_bad_arg("df", message = "must contain columns `row`, `col`, and `value`.")
  }
  out <- out[, c("row", "col", "value"), drop = FALSE]
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
  out <- .mc_set_attrs(out, extra_attrs)
  if (!is.null(native_i)) {
    out <- .mc_stamp_pair_pos(out, native_i, native_j, native_rn, native_cn)
  } else if (!is.null(pos)) {
    if (!is.matrix(pos) || nrow(pos) != nrow(out) || !all(c("i", "j") %in% colnames(pos))) {
      abort_internal("positional pair mapping does not match the edge list.")
    }
    rownames(pos) <- rownames(out)
    attr(out, "corr_pair_pos") <- pos
  }
  out
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

#' Duplicate names on either axis cannot be recovered by label matching.
#' A single missing or blank label remains matchable and is not ambiguous.
#' @keywords internal
#' @noRd
.mc_axis_names_ambiguous <- function(row_names, col_names) {
  (is.character(row_names) && anyDuplicated(row_names) > 0L) ||
    (is.character(col_names) && anyDuplicated(col_names) > 0L)
}

#' Positional plot keys when labels are duplicated, missing, or blank
#' @keywords internal
#' @noRd
.mc_axis_names_positional <- function(row_names, col_names) {
  needs <- function(nm) {
    is.character(nm) && (
      anyDuplicated(nm) > 0L || anyNA(nm) || any(!nzchar(nm))
    )
  }
  needs(row_names) || needs(col_names)
}

#' Continuous-by-discrete estimators pair different variables on each axis
#' @keywords internal
#' @noRd
.mc_cross_axis_estimator <- function(x) {
  cls <- attr(x, "corr_estimator_class", exact = TRUE) %||% .mc_corr_estimator_class(x)
  cls[[1L]] %in% c("biserial_corr", "polyserial_corr")
}

#' Square result whose axes are the same variables, including directed squares.
#' Named-alike biserial and polyserial squares stay cross-axis results.
#' @keywords internal
#' @noRd
.mc_same_axis_square <- function(x) {
  if (.mc_cross_axis_estimator(x)) {
    return(FALSE)
  }
  dm <- .mc_corr_dim(x)
  if (length(dm) != 2L || !isTRUE(dm[[1L]] == dm[[2L]])) {
    return(FALSE)
  }
  if (isTRUE(attr(x, "corr_symmetric", exact = TRUE))) {
    return(TRUE)
  }
  dn <- .mc_corr_dimnames(x)
  identical(dn[[1L]], dn[[2L]]) || (is.null(dn[[1L]]) && is.null(dn[[2L]]))
}

#' Attach positions only when duplicate names make labels ambiguous
#' @keywords internal
#' @noRd
.mc_stamp_pair_pos <- function(df, i, j, row_names, col_names) {
  if (!nrow(df) || !.mc_axis_names_ambiguous(row_names, col_names)) {
    return(df)
  }
  if (length(i) != nrow(df) || length(j) != nrow(df)) {
    abort_internal("positional pair mapping does not match the edge list.")
  }
  pos <- cbind(i = as.integer(i), j = as.integer(j))
  rownames(pos) <- as.character(seq_len(nrow(df)))
  rownames(df) <- rownames(pos)
  attr(df, "corr_pair_pos") <- pos
  df
}

.mc_abort_ambiguous_pairs <- function() {
  cli::cli_abort(
    c(
      "Duplicate variable names have no positional pair mapping.",
      "i" = "Pair identity cannot be recovered by matching labels.",
      "i" = "Rebuild the edge list from the correlation result."
    ),
    class = c("matrixCorr_ambiguous_pairs", "matrixCorr_error")
  )
}

.mc_abort_stale_pair_pos <- function() {
  cli::cli_abort(
    "Stored positional pair mapping does not cover these rows or no longer matches their labels and dimensions.",
    class = c("matrixCorr_stale_pair_pos", "matrixCorr_error")
  )
}

#' Whole positions inside source dimensions, before any truncation
#' @keywords internal
#' @noRd
.mc_validated_native_positions <- function(i, j, dm) {
  whole <- function(v) {
    is.numeric(v) && all(is.finite(v)) && all(v == trunc(v))
  }
  dm <- as.integer(dm)
  if (!whole(i) || !whole(j) || length(dm) != 2L || any(!is.finite(dm))) {
    .mc_abort_stale_pair_pos()
  }
  ii <- as.integer(i)
  jj <- as.integer(j)
  if (any(ii < 1L | ii > dm[[1L]] | jj < 1L | jj > dm[[2L]])) {
    .mc_abort_stale_pair_pos()
  }
  list(i = ii, j = jj)
}

#' Validate a row-name keyed position matrix against the current edge rows
#' @keywords internal
#' @noRd
.mc_validated_pair_pos <- function(pos, df, dm, dn) {
  if (!is.matrix(pos) || !all(c("i", "j") %in% colnames(pos)) || is.null(rownames(pos))) {
    .mc_abort_stale_pair_pos()
  }
  pos_keys <- rownames(pos)
  if (anyNA(pos_keys) || anyDuplicated(pos_keys) > 0L) {
    .mc_abort_stale_pair_pos()
  }
  keys <- rownames(df)
  if (is.null(keys) || anyNA(keys) || anyDuplicated(keys) > 0L || anyNA(match(keys, pos_keys))) {
    .mc_abort_stale_pair_pos()
  }
  ii_raw <- pos[keys, "i"]
  jj_raw <- pos[keys, "j"]
  whole <- function(v) {
    is.numeric(v) && all(is.finite(v)) && all(v == trunc(v))
  }
  if (!whole(ii_raw) || !whole(jj_raw)) {
    .mc_abort_stale_pair_pos()
  }
  ii <- as.integer(ii_raw)
  jj <- as.integer(jj_raw)
  if (length(dm) != 2L || any(!is.finite(dm)) ||
      anyNA(ii) || anyNA(jj) ||
      any(ii < 1L | ii > dm[[1L]] | jj < 1L | jj > dm[[2L]])) {
    .mc_abort_stale_pair_pos()
  }
  rn <- dn[[1L]]
  cn <- dn[[2L]]
  if (!is.null(rn) && !identical(as.character(rn[ii]), as.character(df$row))) {
    .mc_abort_stale_pair_pos()
  }
  if (!is.null(cn) && !identical(as.character(cn[jj]), as.character(df$col))) {
    .mc_abort_stale_pair_pos()
  }
  list(i = ii, j = jj)
}

#' Convert C++ triplets to public edge-list payload
#' @keywords internal
#' @noRd
.mc_triplets_to_edge_list <- function(triplets,
                                      dimnames = NULL) {
  check_same_length(triplets$i, triplets$j, arg_x = "triplets$i", arg_y = "triplets$j")
  check_same_length(triplets$i, triplets$x, arg_x = "triplets$i", arg_y = "triplets$x")
  ii <- as.integer(triplets$i)
  jj <- as.integer(triplets$j)

  rn <- NULL
  cn <- NULL
  if (is.list(dimnames) && length(dimnames) == 2L) {
    rn <- dimnames[[1L]]
    cn <- dimnames[[2L]]
  }

  df <- data.frame(
    row = if (is.null(rn)) as.character(ii) else rn[ii],
    col = if (is.null(cn)) as.character(jj) else cn[jj],
    value = as.numeric(triplets$x),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  .mc_stamp_pair_pos(df, ii, jj, rn, cn)
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
    df = .mc_triplets_to_edge_list(triplets, dimnames = source_dimnames),
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
.mc_matrix_selected_entries <- function(mat, threshold = 0, diag = TRUE, symmetric = NULL) {
  check_matrix_dims(mat, arg = "mat")
  check_scalar_numeric(threshold, arg = "threshold", lower = 0, closed_lower = TRUE)
  check_bool(diag, arg = "diag")

  empty <- list(
    i = integer(),
    j = integer(),
    value = numeric(),
    row = character(),
    col = character(),
    row_names = rownames(mat),
    col_names = colnames(mat)
  )
  if (nrow(mat) == 0L || ncol(mat) == 0L) {
    return(empty)
  }

  same_dimnames <- isTRUE(identical(rownames(mat), colnames(mat))) ||
    (is.null(rownames(mat)) && is.null(colnames(mat)))
  use_upper <- if (is.null(symmetric)) {
    isTRUE(nrow(mat) == ncol(mat)) &&
      same_dimnames &&
      isTRUE(isSymmetric(mat, check.attributes = FALSE))
  } else {
    isTRUE(symmetric) &&
      isTRUE(nrow(mat) == ncol(mat)) &&
      same_dimnames
  }
  if (use_upper) {
    idx <- upper.tri(mat, diag = isTRUE(diag))
  } else {
    idx <- matrix(TRUE, nrow = nrow(mat), ncol = ncol(mat))
    if (!isTRUE(diag)) {
      d <- seq_len(min(nrow(mat), ncol(mat)))
      idx[cbind(d, d)] <- FALSE
    }
  }
  ij <- which(idx, arr.ind = TRUE)
  vals <- mat[idx]
  keep <- is.finite(vals) & (abs(vals) >= as.numeric(threshold))
  if (!any(keep)) {
    return(empty)
  }
  i <- ij[keep, 1L]
  j <- ij[keep, 2L]
  rn <- rownames(mat)
  cn <- colnames(mat)
  list(
    i = i,
    j = j,
    value = as.numeric(vals[keep]),
    row = if (is.null(rn)) as.character(i) else rn[i],
    col = if (is.null(cn)) as.character(j) else cn[j],
    row_names = rn,
    col_names = cn
  )
}

.mc_matrix_to_edge_list <- function(mat,
                                    threshold = 0,
                                    diag = TRUE,
                                    symmetric = NULL) {
  entries <- .mc_matrix_selected_entries(
    mat,
    threshold = threshold,
    diag = diag,
    symmetric = symmetric
  )
  df <- data.frame(
    row = entries$row,
    col = entries$col,
    value = entries$value,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  .mc_stamp_pair_pos(df, entries$i, entries$j, entries$row_names, entries$col_names)
}

#' Convert a symmetric matrix to thresholded sparse matrix
#' @keywords internal
#' @noRd
.mc_matrix_to_sparse_thresholded <- function(mat,
                                             threshold = 0,
                                             diag = TRUE,
                                             symmetric = NULL) {
  check_matrix_dims(mat, arg = "mat")
  check_scalar_numeric(threshold, arg = "threshold", lower = 0, closed_lower = TRUE)
  check_bool(diag, arg = "diag")

  dn <- dimnames(mat)
  same_dimnames <- isTRUE(identical(rownames(mat), colnames(mat))) ||
    (is.null(rownames(mat)) && is.null(colnames(mat)))
  use_upper <- if (is.null(symmetric)) {
    isTRUE(nrow(mat) == ncol(mat)) &&
      same_dimnames &&
      isTRUE(isSymmetric(mat, check.attributes = FALSE))
  } else {
    isTRUE(symmetric) &&
      isTRUE(nrow(mat) == ncol(mat)) &&
      same_dimnames
  }
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

#' Fast route for already-validated output routing
#' @keywords internal
#' @noRd
.mc_finalize_corr_output_fast <- function(x,
                                          output = "matrix",
                                          threshold = 0,
                                          diag = TRUE) {
  if (identical(output, "matrix") &&
    isTRUE(is.finite(threshold)) && isTRUE(as.numeric(threshold) == 0) &&
    isTRUE(diag) &&
    inherits(x, "corr_matrix")) {
    attr(x, "output") <- "matrix"
    attr(x, "threshold") <- 0
    attr(x, "diag") <- TRUE
    return(x)
  }
  .mc_finalize_corr_output(
    x = x,
    output = output,
    threshold = threshold,
    diag = diag
  )
}

#' Route a computed correlation matrix object to requested output shape
#' @keywords internal
#' @noRd
.mc_finalize_corr_output <- function(x,
                                     output = c("matrix", "sparse", "edge_list"),
                                     threshold = 0,
                                     diag = TRUE) {
  if (inherits(x, "corr_matrix") &&
    is.character(output) &&
    length(output) == 1L &&
    identical(output, "matrix") &&
    is.numeric(threshold) &&
    length(threshold) == 1L &&
    !is.na(threshold) &&
    as.numeric(threshold) == 0 &&
    isTRUE(diag)) {
    attr(x, "output") <- "matrix"
    attr(x, "threshold") <- 0
    attr(x, "diag") <- TRUE
    return(x)
  }
  cfg <- .mc_validate_output_args(output = output, threshold = threshold, diag = diag)
  if (identical(cfg$output, "matrix") && inherits(x, "corr_matrix")) {
    attr(x, "output") <- "matrix"
    attr(x, "threshold") <- 0
    attr(x, "diag") <- TRUE
    return(x)
  }
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
  source_symmetric_attr <- attr(x, "corr_symmetric", exact = TRUE)
  source_symmetric <- if (is.logical(source_symmetric_attr) &&
    length(source_symmetric_attr) == 1L &&
    !is.na(source_symmetric_attr)) {
    isTRUE(source_symmetric_attr)
  } else {
    isTRUE(nrow(mat) == ncol(mat)) &&
      (isTRUE(identical(rownames(mat), colnames(mat))) ||
        (is.null(rownames(mat)) && is.null(colnames(mat)))) &&
      isTRUE(isSymmetric(mat, check.attributes = FALSE))
  }

  if (identical(cfg$output, "matrix")) {
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

  if (identical(cfg$output, "edge_list")) {
    out <- .mc_new_corr_edge_list(
      df = .mc_matrix_to_edge_list(
        mat,
        threshold = cfg$threshold,
        diag = cfg$diag,
        symmetric = source_symmetric
      ),
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
    x = .mc_matrix_to_sparse_thresholded(
      mat,
      threshold = cfg$threshold,
      diag = cfg$diag,
      symmetric = source_symmetric
    ),
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
  if (inherits(x, "corr_edge_list")) {
    return(list(NULL, NULL))
  }
  dimnames(x)
}

#' Coerce correlation object to edge data frame
#' @keywords internal
#' @noRd
.mc_corr_as_edge_df <- function(x) {
  output <- attr(x, "output", exact = TRUE) %||% "matrix"

  if (inherits(x, "corr_edge_list")) {
    out <- x
    class(out) <- "data.frame"
    out <- as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
    if (all(c("row", "col", "value") %in% names(out))) {
      return(out[, c("row", "col", "value"), drop = FALSE])
    }
    if (all(c("i", "j", "x") %in% names(out))) {
      dn <- .mc_corr_dimnames(x)
      loc <- .mc_validated_native_positions(out$i, out$j, .mc_corr_dim(x))
      ii <- loc$i
      jj <- loc$j
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
    diag = if (identical(output, "matrix")) TRUE else (attr(x, "diag", exact = TRUE) %||% TRUE),
    symmetric = attr(x, "corr_symmetric", exact = TRUE)
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
  edf <- .mc_edge_indexed_pairs(x)
  if (!nrow(edf)) {
    return(mat)
  }
  keep <- is.finite(edf$i) & is.finite(edf$j)
  ii <- as.integer(edf$i[keep])
  jj <- as.integer(edf$j[keep])
  vv <- as.numeric(edf$value[keep])
  if (!length(ii)) {
    return(mat)
  }
  mat[cbind(ii, jj)] <- vv

  symmetric <- isTRUE(attr(x, "corr_symmetric", exact = TRUE))
  if (symmetric) {
    mat[cbind(jj, ii)] <- vv
  }
  mat
}

#' Edge rows with validated source positions
#' @keywords internal
#' @noRd
.mc_edge_indexed_pairs <- function(x) {
  if (!all(c("row", "col", "value") %in% names(x)) &&
      all(c("i", "j", "x") %in% names(x))) {
    dm <- .mc_corr_dim(x)
    dn <- .mc_corr_dimnames(x)
    loc <- .mc_validated_native_positions(x$i, x$j, dm)
    ii <- loc$i
    jj <- loc$j
    df <- data.frame(
      row = if (is.null(dn[[1L]])) as.character(ii) else as.character(dn[[1L]][ii]),
      col = if (is.null(dn[[2L]])) as.character(jj) else as.character(dn[[2L]][jj]),
      value = as.numeric(x$x),
      i = ii,
      j = jj,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    return(df)
  }
  if (!all(c("row", "col", "value") %in% names(x))) {
    df <- .mc_corr_as_edge_df(x)
  } else {
    df <- data.frame(
      row = as.character(x$row),
      col = as.character(x$col),
      value = as.numeric(x$value),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  dm <- .mc_corr_dim(x)
  dn <- .mc_corr_dimnames(x)
  if (!nrow(df)) {
    df$i <- integer()
    df$j <- integer()
    return(df)
  }
  if (!.mc_axis_names_ambiguous(dn[[1L]], dn[[2L]])) {
    ii <- if (is.null(dn[[1L]])) suppressWarnings(as.integer(df$row)) else match(df$row, dn[[1L]])
    jj <- if (is.null(dn[[2L]])) suppressWarnings(as.integer(df$col)) else match(df$col, dn[[2L]])
    df$i <- as.integer(ii)
    df$j <- as.integer(jj)
    return(df)
  }
  pos <- attr(x, "corr_pair_pos", exact = TRUE)
  if (is.null(pos)) {
    .mc_abort_ambiguous_pairs()
  }
  rownames(df) <- rownames(x)
  loc <- .mc_validated_pair_pos(pos, df, dm, dn)
  df$i <- loc$i
  df$j <- loc$j
  df
}

#' Indexed pairs for matrix, sparse, and edge results
#' @keywords internal
#' @noRd
.mc_resolve_indexed_pairs <- function(x) {
  if (inherits(x, "corr_edge_list")) {
    return(.mc_edge_indexed_pairs(x))
  }
  if (identical(attr(x, "output", exact = TRUE), "sparse") || inherits(x, "sparseMatrix")) {
    sx <- Matrix::summary(x)
    ii <- as.integer(sx$i)
    jj <- as.integer(sx$j)
    dn <- .mc_corr_dimnames(x)
    return(data.frame(
      row = if (is.null(dn[[1L]])) as.character(ii) else as.character(dn[[1L]][ii]),
      col = if (is.null(dn[[2L]])) as.character(jj) else as.character(dn[[2L]][jj]),
      value = as.numeric(sx$x),
      i = ii,
      j = jj,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  output <- attr(x, "output", exact = TRUE) %||% "matrix"
  .mc_matrix_edge_frame(
    as.matrix(x),
    threshold = if (identical(output, "matrix")) 0 else (attr(x, "threshold", exact = TRUE) %||% 0),
    diag = if (identical(output, "matrix")) TRUE else (attr(x, "diag", exact = TRUE) %||% TRUE),
    symmetric = attr(x, "corr_symmetric", exact = TRUE)
  )
}

#' Matrix entries in edge order, retaining positions
#' @keywords internal
#' @noRd
.mc_matrix_edge_frame <- function(mat, threshold = 0, diag = TRUE, symmetric = NULL) {
  entries <- .mc_matrix_selected_entries(
    mat,
    threshold = threshold,
    diag = diag,
    symmetric = symmetric
  )
  data.frame(
    row = entries$row,
    col = entries$col,
    value = entries$value,
    i = as.integer(entries$i),
    j = as.integer(entries$j),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Correlation pairs abstraction
#' @keywords internal
#' @noRd
.mc_corr_pairs <- function(x, diag = FALSE, indexed_pairs = NULL) {
  check_bool(diag, arg = "diag")
  out <- if (is.null(indexed_pairs)) .mc_resolve_indexed_pairs(x) else indexed_pairs
  if (!isTRUE(diag) && nrow(out) && .mc_same_axis_square(x)) {
    keep <- out$i != out$j
    if (!isTRUE(all(keep))) {
      out <- out[keep, , drop = FALSE]
    }
  }
  rownames(out) <- NULL
  out
}

