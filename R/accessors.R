#' Access MatrixCorr Estimates, Confidence Intervals, and Tidy Results
#'
#' @description
#' Lightweight accessors for matrixCorr result objects. These functions do not
#' change the stored object structure; they just provide a stable way to extract the
#' estimate matrix, confidence intervals, or a pairwise data frame without
#' reading attributes directly.
#'
#' @param x,object A matrixCorr result, summary object, or scalar result.
#' @param ... Additional arguments passed to methods.
#' @param diag Logical; include diagonal entries for matrix-style results.
#'   Default is \code{FALSE}.
#' @param triangle For matrix-style correlation results, which triangle to
#'   return from \code{tidy()}: \code{"upper"} (default), \code{"lower"}, or
#'   \code{"full"}. The default avoids duplicate pairs for symmetric matrices;
#'   directed matrices use all off-diagonal entries unless \code{triangle} is
#'   supplied explicitly.
#' @param parm Ignored; included for compatibility with \code{stats::confint()}.
#' @param level Confidence level requested by \code{stats::confint()}. MatrixCorr
#'   objects store pre-computed intervals; if \code{level} differs from the
#'   stored level, recompute the original object with the desired level.
#'
#' @return
#' \code{estimate()} returns the primary estimate in its natural shape: a
#' matrix-like result returns a matrix, an edge-list result returns a data frame,
#' and a scalar result returns a numeric value.
#'
#' \code{tidy()} returns a data frame with columns such as \code{item1},
#' \code{item2}, \code{estimate}, \code{lwr}, \code{upr}, diagnostics, and
#' inferential quantities when available.
#'
#' \code{confint()} returns a data frame containing confidence limits when they
#' are available.
#'
#' \code{ci()} returns the stored confidence-interval payload. It is mainly a
#' structured alternative to reading \code{attr(x, "ci")} directly.
#'
#' @examples
#' X <- cbind(a = 1:10, b = 1:10 + rnorm(10), c = rnorm(10))
#' fit <- pearson_corr(X, ci = TRUE)
#'
#' estimate(fit)
#' tidy(fit)
#' ci(fit)
#' confint(fit)
#'
#' @author Thiago de Paula Oliveira
#' @export
estimate <- function(x, ...) {
  UseMethod("estimate")
}

#' @importFrom generics tidy
#' @rdname estimate
#' @export
generics::tidy

#' @rdname estimate
#' @export
ci <- function(x, ...) {
  UseMethod("ci")
}

#' @rdname estimate
#' @export
estimate.corr_result <- function(x, ...) {
  if (inherits(x, "corr_edge_list")) {
    out <- .mc_corr_as_edge_df(x)
    names(out)[names(out) == "row"] <- "item1"
    names(out)[names(out) == "col"] <- "item2"
    names(out)[names(out) == "value"] <- "estimate"
    out <- out[, intersect(c("item1", "item2", "estimate"), names(out)), drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  .mc_plain_numeric_matrix(.mc_corr_as_dense_matrix(x))
}

#' @rdname estimate
#' @export
estimate.summary.corr_result <- function(x, ...) {
  df <- .mc_tidy_summary_df(x)
  if ("estimate" %in% names(df)) {
    return(df$estimate)
  }
  abort_bad_arg("x", message = "does not contain an `estimate` column.")
}

#' @rdname estimate
#' @export
estimate.summary.matrixCorr <- estimate.summary.corr_result

#' @rdname estimate
#' @export
coef.corr_result <- function(object, ...) {
  estimate(object, ...)
}

#' @rdname estimate
#' @export
ci.corr_result <- function(x, ...) {
  out <- attr(x, "ci", exact = TRUE)
  if (is.null(out)) {
    cli::cli_abort("No confidence intervals are available for this object.")
  }
  out
}

#' @rdname estimate
#' @export
ci.default <- function(x, ...) {
  .mc_ci_payload(x)
}

#' @rdname estimate
#' @export
ci.ba <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.ba_matrix <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.ba_repeated <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.ba_repeated_matrix <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.ccc <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.ccc_ci <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.ccc_glmm <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.cia <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.cia_ci <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.cia_rm <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.cohen_kappa <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.gwet_ac <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.icc <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.icc_overall <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.icc_rm_reml <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.krippendorff_alpha <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.multirater_kappa <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.partial_corr <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.prob_agree <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.rmcorr <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.rmcorr_matrix <- function(x, ...) .mc_ci_payload(x)

#' @rdname estimate
#' @export
ci.weighted_kappa <- function(x, ...) .mc_ci_payload(x)

.mc_ci_payload <- function(x) {
  if (inherits(x, "ccc_glmm")) {
    return(.mc_ccc_glmm_ci_payload(x))
  }
  out <- attr(x, "ci", exact = TRUE)
  if (!is.null(out) && !is.logical(out)) {
    return(out)
  }
  if (inherits(x, c("ba", "ba_repeated"))) {
    return(.mc_ba_lines_ci_payload(x))
  }
  if (inherits(x, c("ba_matrix", "ba_repeated_matrix"))) {
    return(.mc_ba_matrix_ci_payload(x))
  }
  if (.mc_is_matrixcorr_like(x) && is.list(x)) {
    out <- x[["ci", exact = TRUE]]
    if (!is.null(out)) {
      return(out)
    }
    lwr <- x[["lwr", exact = TRUE]] %||% x[["lwr.ci", exact = TRUE]]
    upr <- x[["upr", exact = TRUE]] %||% x[["upr.ci", exact = TRUE]]
    if (!is.null(lwr) && !is.null(upr)) {
      return(list(
        lwr = lwr,
        upr = upr,
        conf.level = x[["conf.level", exact = TRUE]] %||%
          x[["conf_level", exact = TRUE]] %||%
          attr(x, "conf.level", exact = TRUE),
        ci.method = x[["ci.method", exact = TRUE]] %||%
          x[["ci_method", exact = TRUE]] %||%
          attr(x, "ci.method", exact = TRUE)
      ))
    }
  }
  sm <- tryCatch(summary(x), error = function(e) NULL)
  if (is.data.frame(sm) && all(c("lwr", "upr") %in% names(sm))) {
    return(list(
      lwr = sm$lwr,
      upr = sm$upr,
      conf.level = attr(sm, "conf.level", exact = TRUE),
      ci.method = attr(sm, "ci.method", exact = TRUE)
    ))
  }
  cli::cli_abort("No confidence intervals are available for this object.")
}

#' @rdname estimate
#' @export
estimate.partial_corr <- function(x, ...) {
  if (!is.null(x$pcor)) {
    return(as.matrix(x$pcor))
  }
  estimate.default(x, ...)
}

#' @rdname estimate
#' @export
estimate.ba <- function(x, ...) .mc_ba_estimate_columns(tidy(x, ...))

#' @rdname estimate
#' @export
estimate.ba_matrix <- estimate.ba

#' @rdname estimate
#' @export
estimate.ba_repeated <- estimate.ba

#' @rdname estimate
#' @export
estimate.ba_repeated_matrix <- estimate.ba

#' @rdname estimate
#' @export
estimate.ccc_glmm <- function(x, ...) .mc_plain_numeric_matrix(x)

#' @rdname estimate
#' @export
coef.ccc_glmm <- function(object, ...) estimate(object, ...)

#' @rdname estimate
#' @export
estimate.default <- function(x, ...) {
  if (is.data.frame(x) && .mc_is_matrixcorr_like(x)) {
    for (nm in c("estimate", "kappa", "alpha", "ac", "rho", "ccc", "cia", "prob_agree")) {
      if (nm %in% names(x)) {
        return(x[[nm]])
      }
    }
  }
  est <- if (is.list(x) && .mc_is_matrixcorr_like(x)) x[["est", exact = TRUE]] else NULL
  if (!is.null(est)) {
    if (is.matrix(est) || is.data.frame(est)) {
      return(.mc_plain_numeric_matrix(est))
    }
    return(as.numeric(est)[1L])
  }
  estimate_value <- if (is.list(x) && .mc_is_matrixcorr_like(x)) x[["estimate", exact = TRUE]] else NULL
  if (!is.null(estimate_value)) {
    return(estimate_value)
  }
  ac <- if (is.list(x) && .mc_is_matrixcorr_like(x)) x[["ac", exact = TRUE]] else NULL
  if (!is.null(ac)) {
    return(ac)
  }
  est_attr <- attr(x, "estimate", exact = TRUE)
  if (!is.null(est_attr)) {
    return(est_attr)
  }
  if (is.matrix(x) && .mc_is_matrixcorr_like(x)) {
    return(.mc_plain_numeric_matrix(x))
  }
  if (is.numeric(x) && length(x) >= 1L && .mc_is_matrixcorr_like(x)) {
    return(unname(as.numeric(x)))
  }
  sm <- tryCatch(summary(x), error = function(e) NULL)
  if (is.data.frame(sm) && .mc_is_matrixcorr_like(sm) && "estimate" %in% names(sm)) {
    return(sm$estimate)
  }
  abort_bad_arg("x", message = "does not expose a matrixCorr estimate.")
}

#' @rdname estimate
#' @export
tidy.corr_result <- function(x,
                             diag = FALSE,
                             triangle = c("upper", "lower", "full"),
                             ...) {
  check_bool(diag, arg = "diag")
  triangle_missing <- missing(triangle)
  triangle <- match.arg(triangle)
  out <- .mc_tidy_corr_result(
    x,
    diag = diag,
    triangle = triangle,
    triangle_missing = triangle_missing
  )
  rownames(out) <- NULL
  out
}

#' @rdname estimate
#' @export
tidy.summary.corr_result <- function(x, ...) {
  .mc_tidy_summary_df(x)
}

#' @rdname estimate
#' @export
tidy.summary.matrixCorr <- tidy.summary.corr_result

#' @rdname estimate
#' @export
tidy.partial_corr <- function(x,
                              diag = FALSE,
  triangle = c("upper", "lower", "full"),
                              ...) {
  if (!is.null(x$pcor)) {
    triangle <- match.arg(triangle)
    pcor_obj <- .mc_structure_corr_matrix(
      x$pcor,
      class_name = "partial_corr_matrix",
      method = x$method %||% attr(x, "method", exact = TRUE) %||% "partial_corr",
      description = "Partial correlation matrix",
      classes = c("partial_corr_matrix", "matrix"),
      extra_attrs = list(
        ci = x$ci %||% attr(x, "ci", exact = TRUE),
        diagnostics = x$diagnostics %||% attr(x, "diagnostics", exact = TRUE),
        inference = x$inference %||% attr(x, "inference", exact = TRUE),
        conf.level = x$conf.level %||% attr(x, "conf.level", exact = TRUE),
        ci.method = x$ci.method %||% attr(x, "ci.method", exact = TRUE)
      )
    )
    return(tidy.corr_result(pcor_obj, diag = diag, triangle = triangle, ...))
  }
  tidy.default(x, ...)
}

#' @rdname estimate
#' @export
tidy.ba <- function(x, ...) .mc_tidy_summary_with_conf_level(summary(x, ...))

#' @rdname estimate
#' @export
tidy.ba_matrix <- tidy.ba

#' @rdname estimate
#' @export
tidy.ba_repeated <- tidy.ba

#' @rdname estimate
#' @export
tidy.ba_repeated_matrix <- tidy.ba

#' @rdname estimate
#' @export
tidy.default <- function(x,
                         diag = FALSE,
                         triangle = c("upper", "lower", "full"),
                         ...) {
  check_bool(diag, arg = "diag")
  triangle <- match.arg(triangle)
  if (is.list(x) && .mc_is_matrixcorr_like(x) && !is.null(x[["est", exact = TRUE]])) {
    return(.mc_tidy_estimate_list(x, diag = diag, triangle = triangle))
  }
  if (is.list(x) && .mc_is_matrixcorr_like(x) && !is.null(x[["estimate", exact = TRUE]])) {
    return(.mc_tidy_named_estimate_list(x))
  }
  sm <- tryCatch(summary(x), error = function(e) NULL)
  if (is.data.frame(sm) && .mc_is_matrixcorr_like(sm)) {
    out <- .mc_tidy_summary_df(sm)
    if ("estimate" %in% names(out)) {
      return(out)
    }
  }
  if (is.numeric(x) && length(x) >= 1L && .mc_is_matrixcorr_like(x)) {
    return(.mc_tidy_scalar_result(x))
  }
  if (is.list(sm) && .mc_is_matrixcorr_like(sm) && !is.null(sm[["estimate", exact = TRUE]])) {
    return(.mc_tidy_named_estimate_list(sm))
  }
  abort_bad_arg("x", message = "does not expose a tidy matrixCorr representation.")
}

#' @rdname estimate
#' @export
confint.corr_result <- function(object, parm, level = NULL, ...) {
  .mc_confint_from_object(object, level = level, ...)
}

#' @rdname estimate
#' @export
confint.summary.corr_result <- function(object, parm, level = NULL, ...) {
  .mc_check_confint_level_value(attr(object, "conf.level", exact = TRUE), level = level)
  .mc_confint_from_tidy(tidy(object, ...), level = NULL)
}

#' @rdname estimate
#' @export
confint.summary.matrixCorr <- confint.summary.corr_result

#' @rdname estimate
#' @export
confint.ba <- function(object, parm, level = NULL, ...) .mc_confint_ba_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.ba_repeated <- function(object, parm, level = NULL, ...) .mc_confint_ba_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.ba_matrix <- function(object, parm, level = NULL, ...) .mc_confint_ba_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.ba_repeated_matrix <- function(object, parm, level = NULL, ...) .mc_confint_ba_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.ccc <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.ccc_ci <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.ccc_glmm <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.chatterjee_xi_scalar <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.cia <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.cia_ci <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.cia_rm <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.cohen_kappa <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.gwet_ac <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.icc <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.icc_overall <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.icc_rm_reml <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.krippendorff_alpha <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.multirater_kappa <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.prob_agree <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.partial_corr <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.rmcorr <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.rmcorr_matrix <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

#' @rdname estimate
#' @export
confint.weighted_kappa <- function(object, parm, level = NULL, ...) .mc_confint_from_object(object, level = level, ...)

.mc_is_matrixcorr_like <- function(x) {
  cls <- class(x)
  isTRUE(attr(x, "corr_result", exact = TRUE)) ||
    identical(attr(x, "package", exact = TRUE), "matrixCorr") ||
    any(cls %in% c(
      "corr_result", "agreement_result", "matrixCorr_ccc_ci",
      "ba", "ba_matrix", "ba_repeated", "ba_repeated_matrix",
      "ccc", "ccc_ci", "ccc_glmm", "cia", "cia_ci", "cia_rm",
      "cohen_kappa", "weighted_kappa", "gwet_ac", "icc", "icc_overall",
      "icc_rm_reml", "krippendorff_alpha", "multirater_kappa",
      "partial_corr", "prob_agree", "rmcorr", "rmcorr_matrix",
      "chatterjee_xi_scalar", "summary.matrixCorr"
    ))
}

.mc_plain_numeric_matrix <- function(x) {
  x <- as.matrix(x)
  matrix(
    as.numeric(x),
    nrow = nrow(x),
    ncol = ncol(x),
    dimnames = dimnames(x)
  )
}

.mc_tidy_index_matrix <- function(mat, diag = FALSE, triangle = c("upper", "lower", "full")) {
  triangle <- match.arg(triangle)
  idx <- switch(
    triangle,
    upper = upper.tri(mat, diag = isTRUE(diag)),
    lower = lower.tri(mat, diag = isTRUE(diag)),
    full = matrix(TRUE, nrow = nrow(mat), ncol = ncol(mat))
  )
  if (!isTRUE(diag) && identical(triangle, "full")) {
    d <- seq_len(min(nrow(mat), ncol(mat)))
    idx[cbind(d, d)] <- FALSE
  }
  idx
}

.mc_complete_symmetric_matrix <- function(mat, symmetric = FALSE) {
  if (!is.matrix(mat) || !isTRUE(symmetric) || nrow(mat) != ncol(mat)) {
    return(mat)
  }
  out <- mat
  miss <- is.na(out)
  transpose <- t(out)
  out[miss & !is.na(transpose)] <- transpose[miss & !is.na(transpose)]
  out
}

.mc_standard_conf_level <- function(x) {
  out <- suppressWarnings(as.numeric(x)[1L])
  if (length(out) && is.finite(out)) out else NA_real_
}

.mc_check_confint_level <- function(df, level = NULL) {
  stored <- .mc_standard_conf_level(df[["conf.level", exact = TRUE]] %||% attr(df, "conf.level", exact = TRUE))
  .mc_check_confint_level_value(stored, level = level)
}

.mc_object_conf_level <- function(object) {
  if (is.list(object)) {
    object[["conf.level", exact = TRUE]] %||%
      object[["conf_level", exact = TRUE]] %||%
      attr(object, "conf.level", exact = TRUE)
  } else {
    attr(object, "conf.level", exact = TRUE)
  }
}

.mc_check_confint_level_value <- function(stored, level = NULL) {
  stored <- .mc_standard_conf_level(stored)
  if (!is.null(level)) {
    requested <- .mc_standard_conf_level(level)
    if (!is.finite(requested)) {
      abort_bad_arg("level", message = "must be a finite confidence level.")
    }
    if (is.finite(stored) && !isTRUE(all.equal(requested, stored, tolerance = sqrt(.Machine$double.eps)))) {
      cli::cli_abort(sprintf(
        "This object contains pre-computed confidence intervals at conf.level = %s. Recompute the correlation object to use a different level.",
        format(stored)
      ))
    }
  }
  invisible(stored)
}

.mc_confint_from_object <- function(object, level = NULL, ...) {
  .mc_check_confint_level_value(.mc_object_conf_level(object), level = level)
  .mc_confint_from_tidy(tidy(object, ...), level = NULL)
}

.mc_confint_ba_from_object <- function(object, level = NULL, ...) {
  .mc_check_confint_level_value(.mc_object_conf_level(object), level = level)
  .mc_confint_ba(tidy(object, ...), level = NULL)
}

.mc_clean_extractor_metadata <- function(df) {
  if (!is.data.frame(df)) {
    return(df)
  }

  conf_level <- df[["conf.level", exact = TRUE]] %||%
    df[["conf_level", exact = TRUE]] %||%
    attr(df, "conf.level", exact = TRUE)
  conf_level <- .mc_standard_conf_level(conf_level)
  df$conf.level <- NULL
  df$conf_level <- NULL
  attr(df, "conf.level") <- NULL

  df$ci.method <- NULL
  df$ci_method <- NULL
  attr(df, "ci.method") <- NULL
  attr(df, "ci_method") <- NULL

  if ("n_complete" %in% names(df)) {
    n_complete <- unique(stats::na.omit(df$n_complete))
    if (length(n_complete) <= 1L) {
      df$n_complete <- NULL
    }
  }

  attributes(df) <- attributes(df)[intersect(names(attributes(df)), c("names", "row.names", "class"))]
  df
}

.mc_tidy_corr_result <- function(x,
                                 diag = FALSE,
                                 triangle = c("upper", "lower", "full"),
                                 triangle_missing = FALSE) {
  triangle <- match.arg(triangle)
  mat <- .mc_plain_numeric_matrix(.mc_corr_as_dense_matrix(x))
  dm <- dim(mat)
  dn <- dimnames(mat)
  symmetric <- isTRUE(attr(x, "corr_symmetric", exact = TRUE)) &&
    length(dm) == 2L &&
    dm[[1L]] == dm[[2L]] &&
    identical(dn[[1L]], dn[[2L]])
  use_triangle <- if (isTRUE(symmetric)) triangle else if (isTRUE(triangle_missing)) "full" else triangle
  cross_matrix <- !identical(dn[[1L]], dn[[2L]]) && identical(use_triangle, "full")
  idx_mat <- .mc_tidy_index_matrix(
    mat,
    diag = isTRUE(diag) || isTRUE(cross_matrix),
    triangle = use_triangle
  )
  idx <- which(idx_mat, arr.ind = TRUE)
  rn <- dn[[1L]] %||% as.character(seq_len(dm[[1L]]))
  cn <- dn[[2L]] %||% as.character(seq_len(dm[[2L]]))
  out <- data.frame(
    item1 = rn[idx[, 1L]],
    item2 = cn[idx[, 2L]],
    estimate = as.numeric(mat[idx]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out <- out[is.finite(out$estimate), , drop = FALSE]
  idx <- list(i = idx[, 1L], j = idx[, 2L], ok = rep(TRUE, nrow(idx)))
  if (nrow(out) != length(idx$i)) {
    keep <- is.finite(as.numeric(mat[cbind(idx$i, idx$j)]))
    idx <- list(i = idx$i[keep], j = idx$j[keep], ok = rep(TRUE, sum(keep)))
  }

  ci <- attr(x, "ci", exact = TRUE)
  if (is.list(ci)) {
    lwr_mat <- .mc_complete_symmetric_matrix(ci$lwr.ci %||% ci$lwr, symmetric = symmetric)
    upr_mat <- .mc_complete_symmetric_matrix(ci$upr.ci %||% ci$upr, symmetric = symmetric)
    out$lwr <- .mc_extract_pair_matrix(lwr_mat, idx, dm)
    out$upr <- .mc_extract_pair_matrix(upr_mat, idx, dm)
    conf_level <- .mc_standard_conf_level(ci$conf.level %||% attr(x, "conf.level", exact = TRUE))
    if (is.finite(conf_level)) {
      attr(out, "conf.level") <- conf_level
    }
    ci_method <- ci$ci.method %||% ci$ci_method %||% attr(x, "ci.method", exact = TRUE) %||%
      attr(x, "ci_method", exact = TRUE)
    if (is.matrix(ci_method)) {
      out$ci.method <- .mc_extract_pair_matrix(
        .mc_complete_symmetric_matrix(ci_method, symmetric = symmetric),
        idx,
        dm,
        character = TRUE
      )
    } else if (!is.null(ci_method) && length(ci_method) >= 1L) {
      attr(out, "ci.method") <- as.character(ci_method[[1L]])
    }
  }

  diag_attr <- attr(x, "diagnostics", exact = TRUE)
  if (is.list(diag_attr)) {
    out$n_complete <- .mc_extract_pair_matrix(
      .mc_complete_symmetric_matrix(diag_attr$n_complete, symmetric = symmetric),
      idx,
      dm,
      integer = TRUE
    )
    out$m <- .mc_extract_pair_matrix(
      .mc_complete_symmetric_matrix(diag_attr$m, symmetric = symmetric),
      idx,
      dm,
      integer = TRUE
    )
    out$se <- .mc_extract_pair_matrix(
      .mc_complete_symmetric_matrix(diag_attr$se, symmetric = symmetric),
      idx,
      dm
    )
    if (!is.null(diag_attr$bootstrap_reps)) {
      out$bootstrap_reps <- as.integer(diag_attr$bootstrap_reps[[1L]])
    }
  }

  inf <- attr(x, "inference", exact = TRUE)
  if (is.list(inf)) {
    out$statistic <- .mc_extract_pair_matrix(
      .mc_complete_symmetric_matrix(inf$statistic, symmetric = symmetric),
      idx,
      dm
    )
    out$df <- .mc_extract_pair_matrix(
      .mc_complete_symmetric_matrix(inf$parameter %||% inf$df, symmetric = symmetric),
      idx,
      dm
    )
    out$p_value <- .mc_extract_pair_matrix(
      .mc_complete_symmetric_matrix(inf$p_value, symmetric = symmetric),
      idx,
      dm
    )
  }

  out <- .mc_clean_extractor_metadata(.mc_drop_all_na_columns(out))
  rownames(out) <- NULL
  out
}

.mc_pair_indices <- function(df, dm, dn) {
  if (!all(c("item1", "item2") %in% names(df)) || length(dm) != 2L) {
    return(list(i = integer(), j = integer(), ok = logical()))
  }
  rn <- dn[[1L]]
  cn <- dn[[2L]]
  ii <- if (is.null(rn)) suppressWarnings(as.integer(df$item1)) else match(df$item1, rn)
  jj <- if (is.null(cn)) suppressWarnings(as.integer(df$item2)) else match(df$item2, cn)
  ok <- is.finite(ii) & is.finite(jj) &
    ii >= 1L & ii <= dm[[1L]] &
    jj >= 1L & jj <= dm[[2L]]
  list(i = as.integer(ii), j = as.integer(jj), ok = ok)
}

.mc_extract_pair_matrix <- function(mat,
                                    idx,
                                    dm,
                                    integer = FALSE,
                                    character = FALSE) {
  n <- length(idx$ok)
  if (isTRUE(character)) {
    out <- rep(NA_character_, n)
  } else if (isTRUE(integer)) {
    out <- rep(NA_integer_, n)
  } else {
    out <- rep(NA_real_, n)
  }
  if (!is.matrix(mat) || !identical(dim(mat), dm) || !n) {
    return(out)
  }
  out[idx$ok] <- mat[cbind(idx$i[idx$ok], idx$j[idx$ok])]
  if (isTRUE(integer)) {
    return(as.integer(round(out)))
  }
  if (isTRUE(character)) {
    return(as.character(out))
  }
  as.numeric(out)
}

.mc_tidy_summary_df <- function(x) {
  conf_level <- attr(x, "conf.level", exact = TRUE)
  out <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  if ("row" %in% names(out) && !"item1" %in% names(out)) names(out)[names(out) == "row"] <- "item1"
  if ("col" %in% names(out) && !"item2" %in% names(out)) names(out)[names(out) == "col"] <- "item2"
  if ("method1" %in% names(out) && !"item1" %in% names(out)) names(out)[names(out) == "method1"] <- "item1"
  if ("method2" %in% names(out) && !"item2" %in% names(out)) names(out)[names(out) == "method2"] <- "item2"
  if ("value" %in% names(out) && !"estimate" %in% names(out)) names(out)[names(out) == "value"] <- "estimate"
  if ("ac" %in% names(out) && !"estimate" %in% names(out)) names(out)[names(out) == "ac"] <- "estimate"
  if ("kappa" %in% names(out) && !"estimate" %in% names(out)) names(out)[names(out) == "kappa"] <- "estimate"
  if ("alpha" %in% names(out) && !"estimate" %in% names(out)) names(out)[names(out) == "alpha"] <- "estimate"
  if ("rho" %in% names(out) && !"estimate" %in% names(out)) names(out)[names(out) == "rho"] <- "estimate"
  if ("ccc" %in% names(out) && !"estimate" %in% names(out)) names(out)[names(out) == "ccc"] <- "estimate"
  if ("cia" %in% names(out) && !"estimate" %in% names(out)) names(out)[names(out) == "cia"] <- "estimate"
  if ("prob_agree" %in% names(out) && !"estimate" %in% names(out)) names(out)[names(out) == "prob_agree"] <- "estimate"
  if ("lwr.ci" %in% names(out) && !"lwr" %in% names(out)) names(out)[names(out) == "lwr.ci"] <- "lwr"
  if ("upr.ci" %in% names(out) && !"upr" %in% names(out)) names(out)[names(out) == "upr.ci"] <- "upr"
  if ("conf_level" %in% names(out) && !"conf.level" %in% names(out)) names(out)[names(out) == "conf_level"] <- "conf.level"
  if ("ci_method" %in% names(out) && !"ci.method" %in% names(out)) names(out)[names(out) == "ci_method"] <- "ci.method"
  if ("conf_level" %in% names(out) && "conf.level" %in% names(out)) out$conf_level <- NULL
  if ("ci_method" %in% names(out) && "ci.method" %in% names(out)) out$ci_method <- NULL
  rownames(out) <- NULL
  attr(out, "conf.level") <- conf_level
  .mc_clean_extractor_metadata(out)
}

.mc_tidy_summary_with_conf_level <- function(x) {
  out <- .mc_tidy_summary_df(x)
  conf_level <- suppressWarnings(as.numeric(attr(x, "conf.level", exact = TRUE)))
  if (length(conf_level) && is.finite(conf_level[1L])) {
    attr(out, "conf.level") <- conf_level[1L]
  }
  .mc_clean_extractor_metadata(out)
}

.mc_ba_estimate_columns <- function(df) {
  keep <- intersect(
    c(
      "item1", "item2", "n_obs", "bias", "sd_loa", "loa_low", "loa_up",
      "width", "loa_multiplier", "sigma2_subject", "sigma2_resid",
      "residual_model", "ar1_rho"
    ),
    names(df)
  )
  if (!length(keep)) {
    abort_bad_arg("x", message = "does not contain Bland-Altman estimates.")
  }
  out <- df[, keep, drop = FALSE]
  rownames(out) <- NULL
  .mc_drop_all_na_columns(out)
}

.mc_confint_ba <- function(df, level = NULL) {
  df <- .mc_tidy_summary_df(df)
  keep <- intersect(
    c("item1", "item2", "bias_lwr", "bias_upr", "lo_lwr", "lo_upr", "up_lwr", "up_upr"),
    names(df)
  )
  if (!all(c("bias_lwr", "bias_upr", "lo_lwr", "lo_upr", "up_lwr", "up_upr") %in% names(df))) {
    abort_bad_arg("object", message = "does not contain Bland-Altman confidence limits.")
  }
  .mc_check_confint_level(df, level = level)
  out <- df[, keep, drop = FALSE]
  rownames(out) <- NULL
  attr(out, "conf.level") <- attr(df, "conf.level", exact = TRUE)
  out <- .mc_clean_extractor_metadata(out)
  attr(out, "conf.level") <- NULL
  out
}

.mc_ba_lines_ci_payload <- function(x) {
  ci_lines <- x[["CI.lines", exact = TRUE]]
  if (is.null(ci_lines)) {
    ci_lines <- x[["ci", exact = TRUE]]
  }
  if (is.null(ci_lines)) {
    cli::cli_abort("No confidence intervals are available for this object.")
  }
  get <- function(nm) suppressWarnings(as.numeric(ci_lines[[nm]]))[1L]
  out <- list(
    bias = c(lwr = get("mean.diff.ci.lower"), upr = get("mean.diff.ci.upper")),
    loa_low = c(lwr = get("lower.limit.ci.lower"), upr = get("lower.limit.ci.upper")),
    loa_up = c(lwr = get("upper.limit.ci.lower"), upr = get("upper.limit.ci.upper")),
    conf.level = attr(x, "conf.level", exact = TRUE),
    ci.method = attr(x, "ci.method", exact = TRUE)
  )
  if (all(!is.finite(unlist(out[1:3], use.names = FALSE)))) {
    cli::cli_abort("No confidence intervals are available for this object.")
  }
  out
}

.mc_ba_matrix_ci_payload <- function(x) {
  fields <- list(
    bias = c("mean_ci_low", "mean_ci_high"),
    loa_low = c("loa_lower_ci_low", "loa_lower_ci_high"),
    loa_up = c("loa_upper_ci_low", "loa_upper_ci_high")
  )
  if (!all(unlist(fields, use.names = FALSE) %in% names(x))) {
    cli::cli_abort("No confidence intervals are available for this object.")
  }
  out <- lapply(fields, function(nm) {
    list(lwr = x[[nm[[1L]], exact = TRUE]], upr = x[[nm[[2L]], exact = TRUE]])
  })
  out$conf.level <- attr(x, "conf.level", exact = TRUE)
  out$ci.method <- attr(x, "ci.method", exact = TRUE)
  out
}

.mc_ccc_glmm_ci_payload <- function(x) {
  lwr <- attr(x, "rho_ccc_lwr.ci", exact = TRUE)
  upr <- attr(x, "rho_ccc_upr.ci", exact = TRUE)
  se <- attr(x, "rho_ccc_se", exact = TRUE)
  if (is.null(lwr) || is.null(upr)) {
    cli::cli_abort("No confidence intervals are available for this object.")
  }

  out <- list(
    lwr = lwr,
    upr = upr,
    se = se,
    total = list(lwr = lwr, upr = upr, se = se),
    inter = list(
      lwr = attr(x, "rho_ccc_inter_lwr.ci", exact = TRUE),
      upr = attr(x, "rho_ccc_inter_upr.ci", exact = TRUE),
      se = attr(x, "rho_ccc_inter_se", exact = TRUE)
    ),
    intra_method1 = list(
      lwr = attr(x, "rho_ccc_intra_method1_lwr.ci", exact = TRUE),
      upr = attr(x, "rho_ccc_intra_method1_upr.ci", exact = TRUE),
      se = attr(x, "rho_ccc_intra_method1_se", exact = TRUE)
    ),
    intra_method2 = list(
      lwr = attr(x, "rho_ccc_intra_method2_lwr.ci", exact = TRUE),
      upr = attr(x, "rho_ccc_intra_method2_upr.ci", exact = TRUE),
      se = attr(x, "rho_ccc_intra_method2_se", exact = TRUE)
    )
  )
  conf_level <- attr(x, "conf.level", exact = TRUE) %||% attr(x, "conf_level", exact = TRUE)
  if (!is.null(conf_level)) {
    out$conf.level <- conf_level
  }
  ci_method <- attr(x, "ci.method", exact = TRUE) %||% attr(x, "ci_method", exact = TRUE)
  if (!is.null(ci_method)) {
    out$ci.method <- ci_method
  }
  out
}

.mc_tidy_estimate_list <- function(x,
                                   diag = FALSE,
                                   triangle = c("upper", "lower", "full")) {
  triangle <- match.arg(triangle)
  est <- x[["est", exact = TRUE]]
  if (is.matrix(est)) {
    out <- .mc_matrix_estimate_table(est, diag = diag, triangle = triangle)
    if (is.matrix(x$lwr.ci)) out$lwr <- x$lwr.ci[cbind(out$.i, out$.j)]
    if (is.matrix(x$upr.ci)) out$upr <- x$upr.ci[cbind(out$.i, out$.j)]
    out$.i <- out$.j <- NULL
    attr(out, "conf.level") <- x[["conf.level", exact = TRUE]] %||% attr(x, "conf.level", exact = TRUE)
    return(.mc_clean_extractor_metadata(.mc_drop_all_na_columns(out)))
  }
  out <- data.frame(
    estimate = as.numeric(est)[1L],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!is.null(x$lwr.ci)) out$lwr <- as.numeric(x$lwr.ci)[1L]
  if (!is.null(x$upr.ci)) out$upr <- as.numeric(x$upr.ci)[1L]
  attr(out, "conf.level") <- x[["conf.level", exact = TRUE]] %||% attr(x, "conf.level", exact = TRUE)
  .mc_clean_extractor_metadata(.mc_drop_all_na_columns(out))
}

.mc_tidy_named_estimate_list <- function(x) {
  out <- data.frame(
    estimate = as.numeric(x[["estimate", exact = TRUE]])[1L],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (nm in c(
    "lwr", "upr", "conf.level", "conf_level", "ci.method", "ci_method",
    "p_value", "statistic", "df", "slope", "n_obs", "n_subjects",
    "n_complete", "n_boot", "bootstrap_reps"
  )) {
    value <- x[[nm, exact = TRUE]]
    if (!is.null(value)) {
      out[[nm]] <- value[[1L]]
    }
  }
  .mc_clean_extractor_metadata(.mc_tidy_summary_df(.mc_drop_all_na_columns(out)))
}

.mc_matrix_estimate_table <- function(mat,
                                      diag = FALSE,
                                      triangle = c("upper", "lower", "full")) {
  triangle <- match.arg(triangle)
  mat <- as.matrix(mat)
  rn <- rownames(mat) %||% as.character(seq_len(nrow(mat)))
  cn <- colnames(mat) %||% as.character(seq_len(ncol(mat)))
  same_items <- identical(rn, cn) && nrow(mat) == ncol(mat)
  use_triangle <- if (same_items) triangle else "full"
  idx_mat <- .mc_tidy_index_matrix(
    mat,
    diag = isTRUE(diag) || !same_items,
    triangle = use_triangle
  )
  idx <- which(idx_mat, arr.ind = TRUE)
  data.frame(
    item1 = rn[idx[, 1L]],
    item2 = cn[idx[, 2L]],
    estimate = as.numeric(mat[idx]),
    .i = idx[, 1L],
    .j = idx[, 2L],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.mc_tidy_scalar_result <- function(x) {
  out <- data.frame(
    estimate = unname(as.numeric(x)[1L]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  ci <- attr(x, "ci", exact = TRUE)
  if (is.list(ci)) {
    out$lwr <- as.numeric(ci$lwr %||% ci$lwr.ci %||% NA_real_)[1L]
    out$upr <- as.numeric(ci$upr %||% ci$upr.ci %||% NA_real_)[1L]
    out$conf.level <- as.numeric(ci$conf.level %||% attr(x, "conf.level", exact = TRUE) %||% NA_real_)[1L]
    out$se <- as.numeric(ci$se %||% NA_real_)[1L]
    out$m <- as.integer(ci$m %||% NA_integer_)[1L]
    out$bootstrap_reps <- as.integer(ci$bootstrap_reps %||% NA_integer_)[1L]
    out$n_complete <- as.integer(ci$n_complete %||% NA_integer_)[1L]
  }
  .mc_clean_extractor_metadata(.mc_drop_all_na_columns(out))
}

.mc_confint_from_tidy <- function(df, level = NULL) {
  df <- .mc_tidy_summary_df(df)
  if (!all(c("lwr", "upr") %in% names(df))) {
    abort_bad_arg("object", message = "does not contain confidence limits.")
  }
  .mc_check_confint_level(df, level = level)
  keep <- intersect(c("item1", "item2", "lwr", "upr"), names(df))
  out <- df[, keep, drop = FALSE]
  rownames(out) <- NULL
  attr(out, "conf.level") <- attr(df, "conf.level", exact = TRUE)
  out <- .mc_clean_extractor_metadata(out)
  attr(out, "conf.level") <- NULL
  out
}

.mc_drop_all_na_columns <- function(df) {
  if (!is.data.frame(df) || !ncol(df)) {
    return(df)
  }
  keep <- vapply(df, function(z) {
    if (is.list(z)) {
      return(TRUE)
    }
    !all(is.na(z))
  }, logical(1L))
  df[, keep, drop = FALSE]
}
