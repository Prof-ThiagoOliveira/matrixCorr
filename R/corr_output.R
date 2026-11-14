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
  out <- data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("row", "col", "value") %in% names(out)) &&
      all(c("i", "j", "x") %in% names(out))) {
    dn <- source_dimnames
    ii <- as.integer(out$i)
    jj <- as.integer(out$j)
    rn <- if (is.list(dn) && length(dn) == 2L) dn[[1L]] else NULL
    cn <- if (is.list(dn) && length(dn) == 2L) dn[[2L]] else NULL
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

  data.frame(
    row = if (is.null(rn)) as.character(ii) else rn[ii],
    col = if (is.null(cn)) as.character(jj) else cn[jj],
    value = as.numeric(triplets$x),
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
.mc_matrix_to_edge_list <- function(mat,
                                    threshold = 0,
                                    diag = TRUE,
                                    symmetric = NULL) {
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

