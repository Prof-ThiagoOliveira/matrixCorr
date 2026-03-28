# Latent categorical / ordinal correlation wrappers

.mc_has_missing_or_bad <- function(x) {
  if (is.factor(x) || is.logical(x)) {
    return(anyNA(x))
  }
  anyNA(x) || any(!is.finite(x))
}

.mc_is_integerish_numeric <- function(x, tol = 1e-8) {
  is.numeric(x) && all(is.na(x) | abs(x - round(x)) <= tol)
}

.mc_as_data_frame <- function(data, arg = as.character(substitute(data))) {
  if (is.data.frame(data)) {
    return(data)
  }
  if (is.matrix(data)) {
    return(as.data.frame(data, stringsAsFactors = FALSE))
  }
  abort_bad_arg(
    arg,
    message = "must be a matrix or data frame."
  )
}

.mc_encode_ordinal_vector <- function(x, binary = FALSE) {
  if (is.logical(x)) {
    levels_x <- c("FALSE", "TRUE")
    code <- ifelse(is.na(x), NA_integer_, ifelse(x, 2L, 1L))
  } else if (is.factor(x)) {
    levels_x <- levels(x)
    code <- as.integer(x)
  } else if (.mc_is_integerish_numeric(x)) {
    vals <- sort(unique(as.numeric(x[!is.na(x)])))
    if (length(vals) == 0L) {
      return(NULL)
    }
    levels_x <- as.character(vals)
    code <- match(as.numeric(x), vals)
  } else {
    return(NULL)
  }

  if (length(levels_x) < 2L) {
    return(NULL)
  }
  if (binary && length(levels_x) != 2L) {
    return(NULL)
  }

  list(code = as.integer(code), levels = levels_x)
}

.mc_extract_discrete_columns <- function(data, kind = c("ordinal", "binary"),
                                         arg = as.character(substitute(data)),
                                         min_cols = 1L) {
  kind <- match.arg(kind)
  want_binary <- identical(kind, "binary")

  if (is.null(dim(data)) && is.atomic(data)) {
    enc <- .mc_encode_ordinal_vector(data, binary = want_binary)
    if (is.null(enc)) {
      msg <- if (want_binary) {
        "must be a binary (two-level) variable."
      } else {
        "must be an ordinal variable."
      }
      abort_bad_arg(arg, message = msg)
    }
    out <- list(enc)
    names(out) <- arg
    return(out)
  }

  df <- .mc_as_data_frame(data, arg = arg)

  keep <- vapply(df, function(col) {
    !is.null(.mc_encode_ordinal_vector(col, binary = want_binary))
  }, logical(1))

  df <- df[keep]
  if (ncol(df) < min_cols) {
    msg <- if (want_binary) {
      "must contain at least {min_cols} binary column{?s}."
    } else {
      paste0(
        "must contain at least {min_cols} ordinal column{?s} ",
        "(factor, ordered factor, logical, or integer-like numeric)."
      )
    }
    abort_bad_arg(arg, message = msg, min_cols = min_cols)
  }

  encoded <- lapply(df, .mc_encode_ordinal_vector, binary = want_binary)
  names(encoded) <- names(df)
  encoded
}

.mc_extract_continuous_matrix <- function(x,
                                          arg = as.character(substitute(x)),
                                          min_cols = 1L) {
  if (is.null(dim(x)) && is.atomic(x) && !is.factor(x)) {
    if (!is.numeric(x)) {
      abort_bad_arg(arg, message = "must be numeric.")
    }
    mat <- matrix(as.numeric(x), ncol = 1L)
    colnames(mat) <- arg
    return(mat)
  }

  df <- .mc_as_data_frame(x, arg = arg)
  keep <- vapply(df, is.numeric, logical(1))
  df <- df[keep]
  if (ncol(df) < min_cols) {
    abort_bad_arg(
      arg,
      message = "must contain at least {min_cols} numeric column{?s}.",
      min_cols = min_cols
    )
  }
  data.matrix(df)
}

.mc_check_latent_missing <- function(x, check_na, arg) {
  check_bool(check_na)
  if (!check_na) {
    return(invisible(NULL))
  }
  bad <- vapply(x, .mc_has_missing_or_bad, logical(1))
  if (any(bad)) {
    abort_bad_arg(
      arg,
      message = "contains missing, NaN, or infinite values.",
      .hint = "Set `check_na = FALSE` to use pairwise complete cases."
    )
  }
  invisible(NULL)
}

.mc_scalar_or_matrix <- function(mat, scalar) {
  if (scalar) {
    return(as.numeric(mat[1L, 1L]))
  }
  mat
}

.mc_attach_scalar_latent <- function(x, method, description,
                                     diagnostics = NULL, thresholds = NULL) {
  structure(
    as.numeric(x),
    method = method,
    description = description,
    package = "matrixCorr",
    diagnostics = diagnostics,
    thresholds = thresholds
  )
}

.mc_print_corr_matrix <- function(x, header, digits = 4, max_rows = NULL,
                                  max_cols = NULL, ...) {
  cat(header, "\n")
  m <- as.matrix(x)
  attributes(m) <- attributes(m)[c("dim", "dimnames")]

  if (!is.null(max_rows) || !is.null(max_cols)) {
    nr <- nrow(m)
    nc <- ncol(m)
    r <- if (is.null(max_rows)) nr else min(nr, max_rows)
    c <- if (is.null(max_cols)) nc else min(nc, max_cols)
    mm <- round(m[seq_len(r), seq_len(c), drop = FALSE], digits)
    print(mm, ...)
    if (nr > r || nc > c) {
      cat(sprintf("... omitted: %d rows, %d cols\n", nr - r, nc - c))
    }
  } else {
    print(round(m, digits), ...)
  }

  invisible(x)
}

.mc_plot_corr_matrix <- function(x, class_name, fill_name, title,
                                 low_color = "indianred1",
                                 high_color = "steelblue1",
                                 mid_color = "white",
                                 value_text_size = 4, ...) {
  check_inherits(x, class_name)

  mat <- as.matrix(x)
  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("Var1", "Var2", fill_name)
  df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))

  ggplot2::ggplot(df, ggplot2::aes(Var2, Var1, fill = .data[[fill_name]])) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", .data[[fill_name]])),
      size = value_text_size,
      color = "black"
    ) +
    ggplot2::scale_fill_gradient2(
      low = low_color,
      high = high_color,
      mid = mid_color,
      midpoint = 0,
      limits = c(-1, 1),
      name = fill_name
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      ...
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = NULL, y = NULL)
}

.mc_summary_corr_matrix <- function(object, header = "Correlation summary") {
  m <- as.matrix(object)
  symmetric <- isTRUE(nrow(m) == ncol(m)) && isTRUE(isSymmetric(unclass(m)))
  vals <- if (symmetric && nrow(m) > 1L) {
    m[upper.tri(m, diag = FALSE)]
  } else {
    as.vector(m)
  }
  vals <- vals[is.finite(vals)]
  n_pairs <- if (symmetric) {
    stats::setNames(nrow(m) * (nrow(m) - 1L) / 2L, NULL)
  } else {
    stats::setNames(nrow(m) * ncol(m), NULL)
  }

  out <- list(
    class = class(object)[1L],
    method = attr(object, "method"),
    description = attr(object, "description"),
    correct = attr(object, "correct"),
    n_variables = if (symmetric) nrow(m) else NA_integer_,
    n_pairs = as.integer(n_pairs),
    n_rows = nrow(m),
    n_cols = ncol(m),
    symmetric = symmetric,
    n_missing = sum(is.na(m)),
    estimate_min = if (length(vals)) min(vals) else NA_real_,
    estimate_max = if (length(vals)) max(vals) else NA_real_,
    threshold_sets = {
      thr <- attr(object, "thresholds")
      if (is.list(thr)) length(thr) else NA_integer_
    }
  )

  diag_attr <- attr(object, "diagnostics")
  if (is.list(diag_attr)) {
    selector <- if (symmetric && nrow(m) > 1L) upper.tri(m, diag = FALSE) else matrix(TRUE, nrow(m), ncol(m))
    pair_pick <- function(x, transform = identity, default = NA_integer_) {
      if (!is.matrix(x) || !identical(dim(x), dim(m))) return(default)
      vals <- transform(x[selector])
      vals <- vals[!is.na(vals)]
      if (!length(vals)) return(default)
      vals
    }

    zero_vals <- pair_pick(diag_attr$zero_cells)
    ncomp_vals <- pair_pick(diag_attr$n_complete)
    boundary_vals <- pair_pick(diag_attr$boundary, as.logical, NA)
    near_vals <- pair_pick(diag_attr$near_boundary, as.logical, NA)
    corrected_vals <- pair_pick(diag_attr$corrected, as.logical, NA)
    conv_vals <- pair_pick(diag_attr$converged, as.logical, NA)

    out$zero_cell_pairs <- if (is.atomic(zero_vals) && length(zero_vals)) sum(zero_vals > 0, na.rm = TRUE) else NA_integer_
    out$boundary_pairs <- if (is.atomic(boundary_vals) && length(boundary_vals)) sum(boundary_vals, na.rm = TRUE) else NA_integer_
    out$near_boundary_pairs <- if (is.atomic(near_vals) && length(near_vals)) sum(near_vals, na.rm = TRUE) else NA_integer_
    out$corrected_pairs <- if (is.atomic(corrected_vals) && length(corrected_vals)) sum(corrected_vals, na.rm = TRUE) else NA_integer_
    out$converged_pairs <- if (is.atomic(conv_vals) && length(conv_vals)) sum(conv_vals, na.rm = TRUE) else NA_integer_
    out$n_complete_min <- if (is.atomic(ncomp_vals) && length(ncomp_vals)) min(ncomp_vals, na.rm = TRUE) else NA_integer_
    out$n_complete_max <- if (is.atomic(ncomp_vals) && length(ncomp_vals)) max(ncomp_vals, na.rm = TRUE) else NA_integer_
    out$optimizer_tol <- if (is.null(diag_attr$optimizer_tol)) NA_real_ else diag_attr$optimizer_tol
  } else {
    out$zero_cell_pairs <- NA_integer_
    out$boundary_pairs <- NA_integer_
    out$near_boundary_pairs <- NA_integer_
    out$corrected_pairs <- NA_integer_
    out$converged_pairs <- NA_integer_
    out$n_complete_min <- NA_integer_
    out$n_complete_max <- NA_integer_
    out$optimizer_tol <- NA_real_
  }

  out$header <- header
  class(out) <- if (identical(header, "Latent correlation summary")) {
    c("summary_latent_corr", "summary_corr_matrix")
  } else {
    "summary_corr_matrix"
  }
  out
}

.mc_diag_discrete <- function(code) {
  keep <- !is.na(code)
  if (sum(keep) < 2L || length(unique(code[keep])) < 2L) {
    return(NA_real_)
  }
  1
}

.mc_pair_complete <- function(x, y, check_na) {
  if (check_na) {
    return(list(x = x, y = y))
  }
  keep <- stats::complete.cases(x, y)
  list(x = x[keep], y = y[keep])
}

.mc_fast_binary_table01 <- function(x, y) {
  matrix(
    tabulate(as.integer(x) + 2L * as.integer(y) + 1L, nbins = 4L),
    nrow = 2L, ncol = 2L
  )
}

.mc_fast_ordinal_table <- function(x, y, n_x, n_y) {
  matrix(
    tabulate(as.integer(x) + (as.integer(y) - 1L) * as.integer(n_x),
             nbins = as.integer(n_x) * as.integer(n_y)),
    nrow = as.integer(n_x),
    ncol = as.integer(n_y)
  )
}

.mc_binary_tau <- function(x01) {
  x01 <- x01[!is.na(x01)]
  -stats::qnorm(mean(x01))
}

.mc_global_cutpoints <- function(code, n_levels) {
  if (n_levels < 2L) {
    return(numeric())
  }
  code <- as.integer(code)
  code <- code[!is.na(code)]
  freq <- tabulate(code, nbins = as.integer(n_levels))
  probs <- freq / sum(freq)
  stats::qnorm(cumsum(probs)[-n_levels])
}

.mc_phi2_exact <- function(a, b, rho) {
  stats::integrate(
    function(z) stats::dnorm(z) * stats::pnorm((b - rho * z) / sqrt(1 - rho^2)),
    lower = -Inf,
    upper = a,
    rel.tol = 1e-12,
    abs.tol = 0
  )$value
}

.mc_binary_prob_exact <- function(rc, cc, rho) {
  p11 <- .mc_phi2_exact(rc, cc, rho)
  p21 <- stats::pnorm(rc) - p11
  p12 <- stats::pnorm(cc) - p11
  p22 <- 1 - stats::pnorm(rc) - p12
  matrix(c(p11, p21, p12, p22), 2L, 2L)
}

.mc_binary_prob_exact_poly <- function(rc, cc, rho) {
  p11 <- .mc_phi2_exact(rc, cc, rho)
  p12 <- stats::pnorm(rc) - p11
  p21 <- stats::pnorm(cc) - p11
  p22 <- 1 - stats::pnorm(rc) - p21
  matrix(c(p11, p21, p12, p22), 2L, 2L)
}

.mc_tetra_opt_exact <- function(tab, rc, cc) {
  fn <- function(rho) {
    P <- .mc_binary_prob_exact(rc, cc, rho)
    if (any(tab > 0 & (!is.finite(P) | P <= 0))) {
      return(Inf)
    }
    -sum(tab * log(P))
  }
  suppressWarnings(stats::optimize(fn, c(-1, 1))$minimum)
}

.mc_poly_opt_exact_binary <- function(tab, rc, cc) {
  fn <- function(rho) {
    P <- .mc_binary_prob_exact_poly(rc, cc, rho)
    P[P <= 0] <- NA_real_
    lP <- log(P)
    lP[!is.finite(lP)] <- NA_real_
    -sum(tab * lP, na.rm = TRUE)
  }
  suppressWarnings(stats::optimize(fn, c(-1, 1))$minimum)
}

.mc_tetra_table_thresholds <- function(tab, correct) {
  tab <- as.matrix(tab)
  if (correct > 0) {
    tab[tab == 0] <- correct
  }
  tab <- tab / sum(tab)
  list(
    row = if (nrow(tab) == 2L) stats::qnorm(rowSums(tab))[1L] else stats::qnorm(cumsum(rowSums(tab))[-nrow(tab)]),
    col = if (ncol(tab) == 2L) stats::qnorm(colSums(tab))[1L] else stats::qnorm(cumsum(colSums(tab))[-ncol(tab)])
  )
}

.mc_poly_table_thresholds <- function(tab, correct) {
  tab <- as.matrix(tab)
  tot <- sum(tab)
  tab <- tab / tot
  if (correct > 0) {
    tab[tab == 0] <- correct / tot
  }
  list(
    row = stats::qnorm(cumsum(rowSums(tab))[-nrow(tab)]),
    col = stats::qnorm(cumsum(colSums(tab))[-ncol(tab)])
  )
}

.mc_tetrachoric_table_psych <- function(tab, correct) {
  tab <- as.matrix(tab)
  if (all(dim(tab) == c(2L, 2L))) {
    if (correct > 0) {
      tab[tab == 0] <- correct
    }
    tab <- tab / sum(tab)
    rc <- stats::qnorm(colSums(tab))[1L]
    cc <- stats::qnorm(rowSums(tab))[1L]
    return(.mc_tetra_opt_exact(tab, rc = rc, cc = cc))
  }
  matrixCorr_tetrachoric_mle_cpp(unclass(tab), correct = correct)
}

.mc_polychoric_table_psych <- function(tab, correct) {
  tab <- as.matrix(tab)
  if (all(dim(tab) == c(2L, 2L))) {
    tot <- sum(tab)
    tab <- tab / tot
    if (correct > 0) {
      tab[tab == 0] <- correct / tot
    }
    cc <- stats::qnorm(cumsum(colSums(tab))[-ncol(tab)])
    rc <- stats::qnorm(cumsum(rowSums(tab))[-nrow(tab)])
    return(.mc_poly_opt_exact_binary(tab, rc = rc, cc = cc))
  }
  matrixCorr_polychoric_mle_cpp(unclass(tab), correct = correct)
}

.mc_tetrachoric_estimate_psych <- function(x, y, correct,
                                           tau_x = NULL, tau_y = NULL,
                                           global = FALSE) {
  x01 <- as.integer(x) - min(as.integer(x), na.rm = TRUE)
  y01 <- as.integer(y) - min(as.integer(y), na.rm = TRUE)
  tab <- .mc_fast_binary_table01(x01, y01)

  if ((sum(tab) > 1L) && (min(tab) == 0L) && (correct > 0)) {
    tab[tab == 0L] <- correct
  }

  if (any(tab == correct) || any(tab == 0L)) {
    if (global) {
      cc <- if (is.null(tau_x)) .mc_binary_tau(x01) else tau_x
      rc <- if (is.null(tau_y)) .mc_binary_tau(y01) else tau_y
      return(.mc_tetra_opt_exact(tab, rc = rc, cc = cc))
    }
    tab_local <- tab / sum(tab)
    rc <- stats::qnorm(colSums(tab_local))[1L]
    cc <- stats::qnorm(rowSums(tab_local))[1L]
    return(.mc_tetra_opt_exact(tab_local, rc = rc, cc = cc))
  }

  if (global) {
    cc <- if (is.null(tau_x)) .mc_binary_tau(x01) else tau_x
    rc <- if (is.null(tau_y)) .mc_binary_tau(y01) else tau_y
    return(matrixCorr_tetrachoric_fixed_cpp(unclass(tab), rc = rc, cc = cc, correct = correct))
  } else {
    return(matrixCorr_tetrachoric_mle_cpp(unclass(tab), correct = correct))
  }
}

.mc_polychoric_estimate_psych <- function(x, y, n_x, n_y, correct,
                                          tau_x = NULL, tau_y = NULL,
                                          global = FALSE) {
  tab <- .mc_fast_ordinal_table(x, y, n_x = n_x, n_y = n_y)
  tot <- sum(tab)
  if (tot <= 0) {
    return(NA_real_)
  }
  tab <- tab / tot

  if (n_x == 2L && n_y == 2L && any(tab == 0)) {
    if (global) {
      return(.mc_poly_opt_exact_binary(tab, rc = tau_x, cc = tau_y))
    }
    tab_local <- tab
    if (correct > 0 && any(tab_local == 0)) {
      tab_local[tab_local == 0] <- correct / tot
    }
    cc <- stats::qnorm(cumsum(colSums(tab_local))[-ncol(tab_local)])
    rc <- stats::qnorm(cumsum(rowSums(tab_local))[-nrow(tab_local)])
    return(.mc_poly_opt_exact_binary(tab_local, rc = rc, cc = cc))
  }

  if (global) {
    return(matrixCorr_polychoric_fixed_cpp(unclass(tab), rc_in = tau_x, cc_in = tau_y))
  } else {
    return(matrixCorr_polychoric_mle_cpp(unclass(tab), correct = correct))
  }
}

.mc_pair_tetrachoric <- function(x, y, correct, check_na,
                                 tau_x = NULL, tau_y = NULL,
                                 global = TRUE) {
  pair <- .mc_pair_complete(x, y, check_na)
  if (length(pair$x) < 2L) {
    return(NA_real_)
  }
  ux <- unique(pair$x)
  uy <- unique(pair$y)
  ux <- ux[!is.na(ux)]
  uy <- uy[!is.na(uy)]
  if (length(ux) < 2L || length(uy) < 2L) {
    return(NA_real_)
  }
  .mc_tetrachoric_estimate_psych(
    pair$x, pair$y, correct = correct,
    tau_x = tau_x, tau_y = tau_y, global = global
  )
}

.mc_pair_polychoric <- function(x, y, n_x, n_y, correct, check_na,
                                tau_x = NULL, tau_y = NULL,
                                global = identical(n_x, n_y)) {
  pair <- .mc_pair_complete(x, y, check_na)
  if (length(pair$x) < 2L) {
    return(NA_real_)
  }
  ux <- unique(pair$x)
  uy <- unique(pair$y)
  ux <- ux[!is.na(ux)]
  uy <- uy[!is.na(uy)]
  if (length(ux) < 2L || length(uy) < 2L) {
    return(NA_real_)
  }
  .mc_polychoric_estimate_psych(
    pair$x, pair$y,
    n_x = n_x, n_y = n_y,
    correct = correct,
    tau_x = tau_x, tau_y = tau_y,
    global = global
  )
}

.mc_pair_polyserial <- function(x, y, check_na) {
  pair <- .mc_pair_complete(x, y, check_na)
  if (length(pair$x) < 2L) {
    return(NA_real_)
  }
  if (length(unique(pair$y[!is.na(pair$y)])) < 2L) {
    return(NA_real_)
  }
  matrixCorr_polyserial_mle_cpp(
    x = as.numeric(pair$x),
    y = as.integer(as.factor(pair$y))
  )
}

.mc_pair_biserial <- function(x, y, check_na) {
  pair <- .mc_pair_complete(x, y, check_na)
  if (length(pair$x) < 2L) {
    return(NA_real_)
  }
  yy <- pair$y[!is.na(pair$y)]
  if (length(unique(yy)) < 2L) {
    return(NA_real_)
  }
  matrixCorr_biserial_latent_cpp(as.numeric(pair$x), as.logical(pair$y == 2L))
}

.mc_encode_list_to_matrix <- function(enc) {
  out <- matrix(NA_integer_, nrow = length(enc[[1L]]$code), ncol = length(enc))
  for (j in seq_along(enc)) {
    out[, j] <- as.integer(enc[[j]]$code)
  }
  colnames(out) <- names(enc)
  out
}

.mc_tau_matrix <- function(enc) {
  max_cuts <- max(vapply(enc, function(z) max(length(z$levels) - 1L, 0L), integer(1)))
  out <- matrix(NA_real_, nrow = max_cuts, ncol = length(enc))
  for (j in seq_along(enc)) {
    tau_j <- .mc_global_cutpoints(enc[[j]]$code, length(enc[[j]]$levels))
    if (length(tau_j) > 0L) {
      out[seq_along(tau_j), j] <- tau_j
    }
  }
  out
}

.mc_latent_pair_metadata <- function(est, tab, n_complete, correct, thresholds,
                                     optimizer_tol = 1e-8) {
  zero_cells <- sum(tab == 0)
  boundary <- isTRUE(correct == 0) && zero_cells > 0L
  near_boundary <- isTRUE(boundary) || (is.finite(est) && abs(est) >= 0.99)
  list(
    zero_cells = as.integer(zero_cells),
    correct = correct,
    corrected = isTRUE(correct > 0) && zero_cells > 0L,
    boundary = boundary,
    near_boundary = near_boundary,
    n_complete = as.integer(n_complete),
    converged = is.finite(est),
    optimizer_tol = optimizer_tol,
    thresholds = thresholds
  )
}

.mc_matrix_diagnostics_template <- function(p, dimnames) {
  mat_int <- matrix(NA_integer_, nrow = p, ncol = p, dimnames = dimnames)
  mat_lgl <- matrix(NA, nrow = p, ncol = p, dimnames = dimnames)
  list(
    zero_cells = mat_int,
    n_complete = mat_int,
    corrected = mat_lgl,
    boundary = mat_lgl,
    near_boundary = mat_lgl,
    converged = mat_lgl,
    optimizer_tol = 1e-8
  )
}

.mc_fill_pair_diag <- function(diag, j, k, est, tab, n_complete, correct) {
  zero_cells <- as.integer(sum(tab == 0))
  boundary <- isTRUE(correct == 0) && zero_cells > 0L
  near_boundary <- isTRUE(boundary) || (is.finite(est) && abs(est) >= 0.99)
  corrected <- isTRUE(correct > 0) && zero_cells > 0L
  converged <- is.finite(est)

  diag$zero_cells[j, k] <- zero_cells
  diag$zero_cells[k, j] <- zero_cells
  diag$n_complete[j, k] <- as.integer(n_complete)
  diag$n_complete[k, j] <- as.integer(n_complete)
  diag$corrected[j, k] <- corrected
  diag$corrected[k, j] <- corrected
  diag$boundary[j, k] <- boundary
  diag$boundary[k, j] <- boundary
  diag$near_boundary[j, k] <- near_boundary
  diag$near_boundary[k, j] <- near_boundary
  diag$converged[j, k] <- converged
  diag$converged[k, j] <- converged
  diag
}

.mc_structure_corr_matrix <- function(mat, class_name, method, description,
                                      diagnostics = NULL, thresholds = NULL,
                                      correct = NULL) {
  mat <- structure(mat, class = c(class_name, "matrix"))
  attr(mat, "method") <- method
  attr(mat, "description") <- description
  attr(mat, "package") <- "matrixCorr"
  attr(mat, "diagnostics") <- diagnostics
  attr(mat, "thresholds") <- thresholds
  attr(mat, "correct") <- correct
  mat
}

#' @title Pairwise Tetrachoric Correlation
#'
#' @description
#' Computes the tetrachoric correlation for either a pair of binary variables
#' or all pairwise combinations of binary columns in a matrix/data frame.
#'
#' @param data A binary vector, matrix, or data frame. In matrix/data-frame mode,
#' only binary columns are retained.
#' @param y Optional second binary vector. When supplied, the function returns a
#' single tetrachoric correlation estimate.
#' @param correct Non-negative continuity correction added to zero-count cells.
#' Default is \code{0.5}.
#' @param check_na Logical (default \code{TRUE}). If \code{TRUE}, missing values
#' are rejected. If \code{FALSE}, pairwise complete cases are used.
#'
#' @return
#' If \code{y} is supplied, a numeric scalar with attributes
#' \code{diagnostics} and \code{thresholds}. Otherwise a symmetric matrix of
#' class \code{tetrachoric_corr} with attributes \code{method},
#' \code{description}, \code{package = "matrixCorr"}, \code{diagnostics},
#' \code{thresholds}, and \code{correct}.
#'
#' @details
#' The tetrachoric correlation assumes that the observed binary variables arise
#' by dichotomising latent standard-normal variables. Let
#' \eqn{Z_1, Z_2 \sim N(0, 1)} with latent correlation \eqn{\rho}, and define
#' observed binary variables by thresholds \eqn{\tau_1, \tau_2}:
#' \deqn{
#' X = \mathbf{1}\{Z_1 > \tau_1\},
#' \qquad
#' Y = \mathbf{1}\{Z_2 > \tau_2\}.
#' }
#' If the observed \eqn{2 \times 2} table has counts
#' \eqn{n_{ij}} for \eqn{i,j \in \{0,1\}}, the marginal proportions determine
#' the thresholds:
#' \deqn{
#' \tau_1 = \Phi^{-1}\!\big(P(X = 0)\big),
#' \qquad
#' \tau_2 = \Phi^{-1}\!\big(P(Y = 0)\big).
#' }
#' The estimator returned here is the maximum-likelihood estimate of the latent
#' correlation \eqn{\rho}, obtained by maximizing the multinomial log-likelihood
#' built from the rectangle probabilities of the bivariate normal distribution:
#' \deqn{
#' \ell(\rho) = \sum_{i=0}^1 \sum_{j=0}^1 n_{ij}\log \pi_{ij}(\rho;\tau_1,\tau_2),
#' }
#' where \eqn{\pi_{ij}} are the four bivariate-normal cell probabilities implied
#' by \eqn{\rho} and the fixed thresholds. The implementation evaluates the
#' likelihood over \eqn{\rho \in (-1,1)} by a coarse search followed by Brent
#' refinement in C++.
#'
#' The argument \code{correct} adds a continuity correction only to zero-count
#' cells before threshold estimation and likelihood evaluation. This stabilises
#' the estimator for sparse tables and mirrors the conventional
#' \code{correct = 0.5} behaviour used in several psychometric implementations.
#' When \code{correct = 0} and the observed contingency table contains zero
#' cells, the fit is non-regular and may be boundary-driven. In those cases the
#' returned object stores sparse-fit diagnostics, including whether the fit was
#' classified as \code{boundary} or \code{near_boundary}.
#'
#' In matrix/data-frame mode, all pairwise tetrachoric correlations are computed
#' between binary columns. Diagonal entries are \code{1} for non-degenerate
#' columns and \code{NA} for columns with fewer than two observed levels.
#' Variable-specific latent thresholds are stored in the \code{thresholds}
#' attribute, and pairwise sparse-fit diagnostics are stored in
#' \code{diagnostics}.
#'
#' \strong{Computational complexity.} For \eqn{p} binary variables, the matrix
#' path evaluates \eqn{p(p-1)/2} pairwise likelihoods. Each pair uses a
#' one-dimensional optimisation with negligible memory overhead beyond the
#' output matrix.
#'
#' @references
#' Pearson, K. (1900). Mathematical contributions to the theory of evolution.
#' VII. On the correlation of characters not quantitatively measurable.
#' \emph{Philosophical Transactions of the Royal Society A}, 195, 1-47.
#'
#' Olsson, U. (1979). Maximum likelihood estimation of the polychoric
#' correlation coefficient. \emph{Psychometrika}, 44(4), 443-460.
#'
#' @examplesIf requireNamespace("mnormt", quietly = TRUE)
#' set.seed(123)
#' n <- 1000
#' Sigma <- matrix(c(
#'   1.00, 0.55, 0.35,
#'   0.55, 1.00, 0.45,
#'   0.35, 0.45, 1.00
#' ), 3, 3, byrow = TRUE)
#'
#' Z <- mnormt::rmnorm(n = n, mean = rep(0, 3), varcov = Sigma)
#' X <- data.frame(
#'   item1 = Z[, 1] > stats::qnorm(0.70),
#'   item2 = Z[, 2] > stats::qnorm(0.60),
#'   item3 = Z[, 3] > stats::qnorm(0.50)
#' )
#'
#' tc <- tetrachoric(X)
#' print(tc, digits = 3)
#' summary(tc)
#' plot(tc)
#'
#' # Interactive viewing (requires shiny)
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   view_corr_shiny(tc)
#' }
#'
#' # latent Pearson correlations used to generate the binary items
#' round(stats::cor(Z), 2)
#'
#' @author Thiago de Paula Oliveira
#' @export
tetrachoric <- function(data, y = NULL, correct = 0.5, check_na = TRUE) {
  check_bool(check_na)
  check_scalar_nonneg(correct, arg = "correct")

  if (is.null(y) && inherits(data, "table")) {
    tab <- as.matrix(data)
    est <- .mc_tetrachoric_table_psych(tab, correct = correct)
    thresholds <- .mc_tetra_table_thresholds(tab, correct = correct)
    diagnostics <- .mc_latent_pair_metadata(
      est = est,
      tab = tab,
      n_complete = sum(tab),
      correct = correct,
      thresholds = thresholds
    )
    return(.mc_attach_scalar_latent(
      est,
      method = "tetrachoric",
      description = "Pairwise tetrachoric correlation",
      diagnostics = diagnostics,
      thresholds = thresholds
    ))
  }

  if (!is.null(y)) {
    check_same_length(data, y, arg_x = "data", arg_y = "y")
    .mc_check_latent_missing(list(data = data, y = y), check_na = check_na, arg = "data")

    enc_x <- .mc_encode_ordinal_vector(data, binary = TRUE)
    enc_y <- .mc_encode_ordinal_vector(y, binary = TRUE)
    if (is.null(enc_x) || is.null(enc_y)) {
      abort_bad_arg(
        "data",
        message = "and {.arg y} must both be binary (two-level) variables."
      )
    }
    pair <- .mc_pair_complete(enc_x$code, enc_y$code, check_na)
    est <- .mc_pair_tetrachoric(pair$x, pair$y, correct, TRUE)
    x01 <- as.integer(pair$x) - min(as.integer(pair$x), na.rm = TRUE)
    y01 <- as.integer(pair$y) - min(as.integer(pair$y), na.rm = TRUE)
    tab <- .mc_fast_binary_table01(x01, y01)
    thresholds <- list(
      data = .mc_binary_tau(x01),
      y = .mc_binary_tau(y01)
    )
    diagnostics <- .mc_latent_pair_metadata(
      est = est,
      tab = tab,
      n_complete = length(pair$x),
      correct = correct,
      thresholds = thresholds
    )
    return(.mc_attach_scalar_latent(
      est,
      method = "tetrachoric",
      description = "Pairwise tetrachoric correlation",
      diagnostics = diagnostics,
      thresholds = thresholds
    ))
  }

  enc <- .mc_extract_discrete_columns(data, kind = "binary", arg = "data", min_cols = 2L)
  .mc_check_latent_missing(
    lapply(enc, `[[`, "code"),
    check_na = check_na,
    arg = "data"
  )

  p <- length(enc)
  x_mat <- .mc_encode_list_to_matrix(enc)
  tau <- vapply(enc, function(z) {
    code <- as.integer(z$code)
    .mc_binary_tau(code - min(code, na.rm = TRUE))
  }, numeric(1))
  out <- matrixCorr_tetrachoric_matrix_cpp(
    x = x_mat,
    tau = tau,
    correct = correct,
    pairwise_complete = !check_na
  )
  dimnames(out) <- list(names(enc), names(enc))
  diag_info <- .mc_matrix_diagnostics_template(p, dimnames(out))
  diag(diag_info$zero_cells) <- 0L
  diag(diag_info$n_complete) <- colSums(!is.na(x_mat))
  diag(diag_info$corrected) <- FALSE
  diag(diag_info$boundary) <- FALSE
  diag(diag_info$near_boundary) <- FALSE
  diag(diag_info$converged) <- is.finite(diag(out))
  for (j in seq_len(p - 1L)) {
    for (k in (j + 1L):p) {
      pair <- .mc_pair_complete(enc[[j]]$code, enc[[k]]$code, check_na)
      x01 <- as.integer(pair$x) - min(as.integer(pair$x), na.rm = TRUE)
      y01 <- as.integer(pair$y) - min(as.integer(pair$y), na.rm = TRUE)
      tab <- .mc_fast_binary_table01(x01, y01)
      diag_info <- .mc_fill_pair_diag(diag_info, j, k, out[j, k], tab, length(pair$x), correct)
    }
  }
  thresholds <- as.list(tau)
  names(thresholds) <- names(enc)

  .mc_structure_corr_matrix(
    out,
    class_name = "tetrachoric_corr",
    method = "tetrachoric",
    description = "Pairwise tetrachoric correlation matrix",
    diagnostics = diag_info,
    thresholds = thresholds,
    correct = correct
  )
}

#' @rdname tetrachoric
#' @method print tetrachoric_corr
#' @param x An object of class \code{tetrachoric_corr}.
#' @param digits Integer; number of decimal places to print.
#' @param max_rows Optional integer; maximum number of rows to display.
#' @param max_cols Optional integer; maximum number of columns to display.
#' @param ... Additional arguments passed to \code{print()}.
#' @export
print.tetrachoric_corr <- function(x, digits = 4, max_rows = NULL,
                                   max_cols = NULL, ...) {
  .mc_print_corr_matrix(
    x, header = "Tetrachoric correlation matrix:",
    digits = digits, max_rows = max_rows, max_cols = max_cols, ...
  )
}

#' @rdname tetrachoric
#' @method plot tetrachoric_corr
#' @param title Plot title. Default is \code{"Tetrachoric correlation heatmap"}.
#' @param low_color Color for the minimum correlation.
#' @param high_color Color for the maximum correlation.
#' @param mid_color Color for zero correlation.
#' @param value_text_size Font size used in tile labels.
#' @export
plot.tetrachoric_corr <- function(x, title = "Tetrachoric correlation heatmap",
                                  low_color = "indianred1",
                                  high_color = "steelblue1",
                                  mid_color = "white",
                                  value_text_size = 4, ...) {
  .mc_plot_corr_matrix(
    x, class_name = "tetrachoric_corr", fill_name = "Tetrachoric",
    title = title, low_color = low_color, high_color = high_color,
    mid_color = mid_color, value_text_size = value_text_size, ...
  )
}

#' @rdname tetrachoric
#' @method summary tetrachoric_corr
#' @param object An object of class \code{tetrachoric_corr}.
#' @export
summary.tetrachoric_corr <- function(object, ...) {
  .mc_summary_corr_matrix(object, header = "Latent correlation summary")
}

#' @title Pairwise Polychoric Correlation
#'
#' @description
#' Computes the polychoric correlation for either a pair of ordinal variables
#' or all pairwise combinations of ordinal columns in a matrix/data frame.
#'
#' @param data An ordinal vector, matrix, or data frame. Supported columns are
#' factors, ordered factors, logical values, or integer-like numerics.
#' In matrix/data-frame mode, only supported ordinal columns are retained.
#' @param y Optional second ordinal vector. When supplied, the function returns
#' a single polychoric correlation estimate.
#' @param correct Non-negative continuity correction added to zero-count cells.
#' Default is \code{0.5}.
#' @param check_na Logical (default \code{TRUE}). If \code{TRUE}, missing values
#' are rejected. If \code{FALSE}, pairwise complete cases are used.
#'
#' @return
#' If \code{y} is supplied, a numeric scalar with attributes
#' \code{diagnostics} and \code{thresholds}. Otherwise a symmetric matrix of
#' class \code{polychoric_corr} with attributes \code{method},
#' \code{description}, \code{package = "matrixCorr"}, \code{diagnostics},
#' \code{thresholds}, and \code{correct}.
#'
#' @details
#' The polychoric correlation generalises the tetrachoric model to ordered
#' categorical variables with more than two levels. It assumes latent
#' standard-normal variables \eqn{Z_1, Z_2} with correlation \eqn{\rho}, and
#' cut-points
#' \eqn{-\infty = \alpha_0 < \alpha_1 < \cdots < \alpha_R = \infty} and
#' \eqn{-\infty = \beta_0 < \beta_1 < \cdots < \beta_C = \infty} such that
#' \deqn{
#' X = r \iff \alpha_{r-1} < Z_1 \le \alpha_r,
#' \qquad
#' Y = c \iff \beta_{c-1} < Z_2 \le \beta_c.
#' }
#' For an observed \eqn{R \times C} contingency table with counts \eqn{n_{rc}},
#' the thresholds are estimated from the marginal cumulative proportions:
#' \deqn{
#' \alpha_r = \Phi^{-1}\!\Big(\sum_{k \le r} P(X = k)\Big),
#' \qquad
#' \beta_c = \Phi^{-1}\!\Big(\sum_{k \le c} P(Y = k)\Big).
#' }
#' Holding those thresholds fixed, the log-likelihood for the latent
#' correlation is
#' \deqn{
#' \ell(\rho) = \sum_{r=1}^{R}\sum_{c=1}^{C}
#' n_{rc} \log \Pr\!\big(
#' \alpha_{r-1} < Z_1 \le \alpha_r,\;
#' \beta_{c-1} < Z_2 \le \beta_c
#' \mid \rho \big),
#' }
#' and the estimator returned is the maximiser over \eqn{\rho \in (-1,1)}.
#' The C++ implementation performs a dense one-dimensional search followed by
#' Brent refinement.
#'
#' The argument \code{correct} adds a non-negative continuity correction to
#' empty cells before marginal threshold estimation and likelihood evaluation.
#' This avoids numerical failures for sparse tables with structurally zero cells.
#' When \code{correct = 0} and zero cells are present, the corresponding fit can
#' be boundary-driven rather than a regular interior maximum-likelihood problem.
#' The returned object stores sparse-fit diagnostics and the thresholds used for
#' estimation so those cases can be inspected explicitly.
#'
#' In matrix/data-frame mode, all pairwise polychoric correlations are computed
#' between supported ordinal columns. Diagonal entries are \code{1} for
#' non-degenerate columns and \code{NA} when a column has fewer than two
#' observed levels.
#'
#' \strong{Computational complexity.} For \eqn{p} ordinal variables, the matrix
#' path evaluates \eqn{p(p-1)/2} bivariate likelihoods. Each pair optimises a
#' single scalar parameter \eqn{\rho}, so the main cost is repeated evaluation
#' of bivariate normal rectangle probabilities.
#'
#' @references
#' Olsson, U. (1979). Maximum likelihood estimation of the polychoric
#' correlation coefficient. \emph{Psychometrika}, 44(4), 443-460.
#'
#' @examplesIf requireNamespace("mnormt", quietly = TRUE)
#' set.seed(124)
#' n <- 1200
#' Sigma <- matrix(c(
#'   1.00, 0.60, 0.40,
#'   0.60, 1.00, 0.50,
#'   0.40, 0.50, 1.00
#' ), 3, 3, byrow = TRUE)
#'
#' Z <- mnormt::rmnorm(n = n, mean = rep(0, 3), varcov = Sigma)
#' Y <- data.frame(
#'   y1 = ordered(cut(
#'     Z[, 1],
#'     breaks = c(-Inf, -0.7, 0.4, Inf),
#'     labels = c("low", "mid", "high")
#'   )),
#'   y2 = ordered(cut(
#'     Z[, 2],
#'     breaks = c(-Inf, -1.0, -0.1, 0.8, Inf),
#'     labels = c("1", "2", "3", "4")
#'   )),
#'   y3 = ordered(cut(
#'     Z[, 3],
#'     breaks = c(-Inf, -0.4, 0.2, 1.1, Inf),
#'     labels = c("A", "B", "C", "D")
#'   ))
#' )
#'
#' pc <- polychoric(Y)
#' print(pc, digits = 3)
#' summary(pc)
#' plot(pc)
#'
#' # Interactive viewing (requires shiny)
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   view_corr_shiny(pc)
#' }
#'
#' # latent Pearson correlations used to generate the ordinal variables
#' round(stats::cor(Z), 2)
#' @author Thiago de Paula Oliveira
#' @export
polychoric <- function(data, y = NULL, correct = 0.5, check_na = TRUE) {
  check_bool(check_na)
  check_scalar_nonneg(correct, arg = "correct")

  if (is.null(y) && inherits(data, "table")) {
    tab <- as.matrix(data)
    est <- .mc_polychoric_table_psych(tab, correct = correct)
    thresholds <- .mc_poly_table_thresholds(tab, correct = correct)
    diagnostics <- .mc_latent_pair_metadata(
      est = est,
      tab = tab,
      n_complete = sum(tab),
      correct = correct,
      thresholds = thresholds
    )
    return(.mc_attach_scalar_latent(
      est,
      method = "polychoric",
      description = "Pairwise polychoric correlation",
      diagnostics = diagnostics,
      thresholds = thresholds
    ))
  }

  if (!is.null(y)) {
    check_same_length(data, y, arg_x = "data", arg_y = "y")
    .mc_check_latent_missing(list(data = data, y = y), check_na = check_na, arg = "data")

    enc_x <- .mc_encode_ordinal_vector(data, binary = FALSE)
    enc_y <- .mc_encode_ordinal_vector(y, binary = FALSE)
    if (is.null(enc_x) || is.null(enc_y)) {
      abort_bad_arg(
        "data",
        message = "and {.arg y} must both be ordinal variables."
      )
    }
    n_x <- length(enc_x$levels)
    n_y <- length(enc_y$levels)
    global <- identical(n_x, n_y)
    pair <- .mc_pair_complete(enc_x$code, enc_y$code, check_na)
    est <- .mc_pair_polychoric(
      pair$x, pair$y,
      n_x = n_x, n_y = n_y,
      correct = correct, check_na = TRUE,
      tau_x = if (global) .mc_global_cutpoints(pair$x, n_x) else NULL,
      tau_y = if (global) .mc_global_cutpoints(pair$y, n_y) else NULL,
      global = global
    )
    tab <- .mc_fast_ordinal_table(pair$x, pair$y, n_x = n_x, n_y = n_y)
    thresholds <- list(
      data = .mc_global_cutpoints(pair$x, n_x),
      y = .mc_global_cutpoints(pair$y, n_y)
    )
    diagnostics <- .mc_latent_pair_metadata(
      est = est,
      tab = tab,
      n_complete = length(pair$x),
      correct = correct,
      thresholds = thresholds
    )
    return(.mc_attach_scalar_latent(
      est,
      method = "polychoric",
      description = "Pairwise polychoric correlation",
      diagnostics = diagnostics,
      thresholds = thresholds
    ))
  }

  enc <- .mc_extract_discrete_columns(data, kind = "ordinal", arg = "data", min_cols = 2L)
  .mc_check_latent_missing(
    lapply(enc, `[[`, "code"),
    check_na = check_na,
    arg = "data"
  )

  p <- length(enc)
  n_levels <- vapply(enc, function(z) length(z$levels), integer(1))
  global_all <- length(unique(n_levels)) == 1L
  tau_mat <- if (global_all) {
    .mc_tau_matrix(enc)
  } else {
    matrix(NA_real_, nrow = 0L, ncol = 0L)
  }
  x_mat <- .mc_encode_list_to_matrix(enc)
  out <- matrixCorr_polychoric_matrix_cpp(
    x = x_mat,
    n_levels = n_levels,
    tau_mat = tau_mat,
    global_all = global_all,
    correct = correct,
    pairwise_complete = !check_na
  )
  dimnames(out) <- list(names(enc), names(enc))
  diag_info <- .mc_matrix_diagnostics_template(p, dimnames(out))
  diag(diag_info$zero_cells) <- 0L
  diag(diag_info$n_complete) <- vapply(enc, function(z) sum(!is.na(z$code)), integer(1))
  diag(diag_info$corrected) <- FALSE
  diag(diag_info$boundary) <- FALSE
  diag(diag_info$near_boundary) <- FALSE
  diag(diag_info$converged) <- is.finite(diag(out))
  for (j in seq_len(p - 1L)) {
    for (k in (j + 1L):p) {
      pair <- .mc_pair_complete(enc[[j]]$code, enc[[k]]$code, check_na)
      tab <- .mc_fast_ordinal_table(
        pair$x, pair$y,
        n_x = length(enc[[j]]$levels),
        n_y = length(enc[[k]]$levels)
      )
      diag_info <- .mc_fill_pair_diag(diag_info, j, k, out[j, k], tab, length(pair$x), correct)
    }
  }
  thresholds <- lapply(enc, function(z) .mc_global_cutpoints(z$code, length(z$levels)))
  names(thresholds) <- names(enc)

  .mc_structure_corr_matrix(
    out,
    class_name = "polychoric_corr",
    method = "polychoric",
    description = "Pairwise polychoric correlation matrix",
    diagnostics = diag_info,
    thresholds = thresholds,
    correct = correct
  )
}

#' @rdname polychoric
#' @method print polychoric_corr
#' @param x An object of class \code{polychoric_corr}.
#' @param digits Integer; number of decimal places to print.
#' @param max_rows Optional integer; maximum number of rows to display.
#' @param max_cols Optional integer; maximum number of columns to display.
#' @param ... Additional arguments passed to \code{print()}.
#' @export
print.polychoric_corr <- function(x, digits = 4, max_rows = NULL,
                                  max_cols = NULL, ...) {
  .mc_print_corr_matrix(
    x, header = "Polychoric correlation matrix:",
    digits = digits, max_rows = max_rows, max_cols = max_cols, ...
  )
}

#' @rdname polychoric
#' @method plot polychoric_corr
#' @param title Plot title. Default is \code{"Polychoric correlation heatmap"}.
#' @param low_color Color for the minimum correlation.
#' @param high_color Color for the maximum correlation.
#' @param mid_color Color for zero correlation.
#' @param value_text_size Font size used in tile labels.
#' @export
plot.polychoric_corr <- function(x, title = "Polychoric correlation heatmap",
                                 low_color = "indianred1",
                                 high_color = "steelblue1",
                                 mid_color = "white",
                                 value_text_size = 4, ...) {
  .mc_plot_corr_matrix(
    x, class_name = "polychoric_corr", fill_name = "Polychoric",
    title = title, low_color = low_color, high_color = high_color,
    mid_color = mid_color, value_text_size = value_text_size, ...
  )
}

#' @rdname polychoric
#' @method summary polychoric_corr
#' @param object An object of class \code{polychoric_corr}.
#' @export
summary.polychoric_corr <- function(object, ...) {
  .mc_summary_corr_matrix(object, header = "Latent correlation summary")
}

#' @title Polyserial Correlation Between Continuous and Ordinal Variables
#'
#' @description
#' Computes polyserial correlations between continuous variables in \code{data}
#' and ordinal variables in \code{y}. Both pairwise vector mode and rectangular
#' matrix/data-frame mode are supported.
#'
#' @usage polyserial(data, y, check_na = TRUE)
#'
#' @param data A numeric vector, matrix, or data frame containing continuous
#' variables.
#' @param y An ordinal vector, matrix, or data frame containing ordinal
#' variables. Supported columns are factors, ordered factors, logical values,
#' or integer-like numerics.
#' @param check_na Logical (default \code{TRUE}). If \code{TRUE}, missing values
#' are rejected. If \code{FALSE}, pairwise complete cases are used.
#'
#' @return
#' If both \code{data} and \code{y} are vectors, a numeric scalar. Otherwise a
#' numeric matrix of class \code{polyserial_corr} with rows corresponding to
#' the continuous variables in \code{data} and columns to the ordinal variables
#' in \code{y}. Matrix outputs carry attributes \code{method},
#' \code{description}, and \code{package = "matrixCorr"}.
#'
#' @details
#' The polyserial correlation assumes a latent bivariate normal model between a
#' continuous variable and an unobserved continuous propensity underlying an
#' ordinal variable. Let
#' \eqn{(X, Z)^\top \sim N_2(0, \Sigma)} with
#' \eqn{\mathrm{corr}(X,Z)=\rho}, and suppose the observed ordinal response
#' \eqn{Y} is formed by cut-points
#' \eqn{-\infty = \beta_0 < \beta_1 < \cdots < \beta_K = \infty}:
#' \deqn{
#' Y = k \iff \beta_{k-1} < Z \le \beta_k.
#' }
#' After standardising the observed continuous variable \eqn{X}, the thresholds
#' are estimated from the marginal proportions of \eqn{Y}. Conditional on an
#' observed \eqn{x_i}, the category probability is
#' \deqn{
#' \Pr(Y_i = k \mid X_i = x_i, \rho)
#' =
#' \Phi\!\left(\frac{\beta_k - \rho x_i}{\sqrt{1-\rho^2}}\right)
#' -
#' \Phi\!\left(\frac{\beta_{k-1} - \rho x_i}{\sqrt{1-\rho^2}}\right).
#' }
#' The returned estimate maximises the log-likelihood
#' \deqn{
#' \ell(\rho) = \sum_{i=1}^{n}\log \Pr(Y_i = y_i \mid X_i = x_i, \rho)
#' }
#' over \eqn{\rho \in (-1,1)} via a one-dimensional Brent search in C++.
#'
#' In vector mode a single estimate is returned. In matrix/data-frame mode,
#' every numeric column of \code{data} is paired with every ordinal column of
#' \code{y}, producing a rectangular matrix of continuous-by-ordinal
#' polyserial correlations.
#'
#' \strong{Computational complexity.} If \code{data} has \eqn{p_x} continuous
#' columns and \code{y} has \eqn{p_y} ordinal columns, the matrix path computes
#' \eqn{p_x p_y} separate one-parameter likelihood optimisations.
#'
#' @references
#' Olsson, U., Drasgow, F., & Dorans, N. J. (1982). The polyserial
#' correlation coefficient. \emph{Psychometrika}, 47(3), 337-347.
#'
#' @examplesIf requireNamespace("mnormt", quietly = TRUE)
#' set.seed(125)
#' n <- 1000
#' Sigma <- matrix(c(
#'   1.00, 0.30, 0.55, 0.20,
#'   0.30, 1.00, 0.25, 0.50,
#'   0.55, 0.25, 1.00, 0.40,
#'   0.20, 0.50, 0.40, 1.00
#' ), 4, 4, byrow = TRUE)
#'
#' Z <- mnormt::rmnorm(n = n, mean = rep(0, 4), varcov = Sigma)
#' X <- data.frame(x1 = Z[, 1], x2 = Z[, 2])
#' Y <- data.frame(
#'   y1 = ordered(cut(
#'     Z[, 3],
#'     breaks = c(-Inf, -0.5, 0.7, Inf),
#'     labels = c("low", "mid", "high")
#'   )),
#'   y2 = ordered(cut(
#'     Z[, 4],
#'     breaks = c(-Inf, -1.0, 0.0, 1.0, Inf),
#'     labels = c("1", "2", "3", "4")
#'   ))
#' )
#'
#' ps <- polyserial(X, Y)
#' print(ps, digits = 3)
#' summary(ps)
#' plot(ps)
#' @author Thiago de Paula Oliveira
#' @export
polyserial <- function(data, y, check_na = TRUE) {
  check_bool(check_na)

  scalar <- is.null(dim(data)) && is.atomic(data) && !is.factor(data) &&
    is.null(dim(y)) && is.atomic(y)

  x_mat <- .mc_extract_continuous_matrix(data, arg = "data", min_cols = 1L)
  y_enc <- .mc_extract_discrete_columns(y, kind = "ordinal", arg = "y", min_cols = 1L)

  if (nrow(x_mat) != length(y_enc[[1L]]$code)) {
    abort_bad_arg(
      "data",
      message = "and {.arg y} must have the same number of observations."
    )
  }

  .mc_check_latent_missing(
    c(as.list(as.data.frame(x_mat)), lapply(y_enc, `[[`, "code")),
    check_na = check_na,
    arg = "data"
  )

  if (scalar) {
    return(.mc_pair_polyserial(x_mat[, 1L], y_enc[[1L]]$code, check_na))
  }

  out <- matrix(NA_real_, nrow = ncol(x_mat), ncol = length(y_enc),
                dimnames = list(colnames(x_mat), names(y_enc)))

  for (j in seq_len(ncol(x_mat))) {
    xj <- x_mat[, j]
    for (k in seq_along(y_enc)) {
      out[j, k] <- .mc_pair_polyserial(xj, y_enc[[k]]$code, check_na)
    }
  }

  .mc_structure_corr_matrix(
    out,
    class_name = "polyserial_corr",
    method = "polyserial",
    description = "Polyserial correlation matrix (continuous x ordinal)"
  )
}

#' @rdname polyserial
#' @method print polyserial_corr
#' @param x An object of class \code{polyserial_corr}.
#' @param digits Integer; number of decimal places to print.
#' @param max_rows Optional integer; maximum number of rows to display.
#' @param max_cols Optional integer; maximum number of columns to display.
#' @param ... Additional arguments passed to \code{print()}.
#' @export
print.polyserial_corr <- function(x, digits = 4, max_rows = NULL,
                                  max_cols = NULL, ...) {
  .mc_print_corr_matrix(
    x, header = "Polyserial correlation matrix:",
    digits = digits, max_rows = max_rows, max_cols = max_cols, ...
  )
}

#' @rdname polyserial
#' @method plot polyserial_corr
#' @param title Plot title. Default is \code{"Polyserial correlation heatmap"}.
#' @param low_color Color for the minimum correlation.
#' @param high_color Color for the maximum correlation.
#' @param mid_color Color for zero correlation.
#' @param value_text_size Font size used in tile labels.
#' @export
plot.polyserial_corr <- function(x, title = "Polyserial correlation heatmap",
                                 low_color = "indianred1",
                                 high_color = "steelblue1",
                                 mid_color = "white",
                                 value_text_size = 4, ...) {
  .mc_plot_corr_matrix(
    x, class_name = "polyserial_corr", fill_name = "Polyserial",
    title = title, low_color = low_color, high_color = high_color,
    mid_color = mid_color, value_text_size = value_text_size, ...
  )
}

#' @rdname polyserial
#' @method summary polyserial_corr
#' @param object An object of class \code{polyserial_corr}.
#' @export
summary.polyserial_corr <- function(object, ...) {
  .mc_summary_corr_matrix(object, header = "Latent correlation summary")
}

#' @title Biserial Correlation Between Continuous and Binary Variables
#'
#' @description
#' Computes biserial correlations between continuous variables in \code{data}
#' and binary variables in \code{y}. Both pairwise vector mode and rectangular
#' matrix/data-frame mode are supported.
#'
#' @usage biserial(data, y, check_na = TRUE)
#'
#' @param data A numeric vector, matrix, or data frame containing continuous
#' variables.
#' @param y A binary vector, matrix, or data frame. In data-frame mode, only
#' two-level columns are retained.
#' @param check_na Logical (default \code{TRUE}). If \code{TRUE}, missing values
#' are rejected. If \code{FALSE}, pairwise complete cases are used.
#'
#' @return
#' If both \code{data} and \code{y} are vectors, a numeric scalar. Otherwise a
#' numeric matrix of class \code{biserial_corr} with rows corresponding to
#' the continuous variables in \code{data} and columns to the binary variables
#' in \code{y}. Matrix outputs carry attributes \code{method},
#' \code{description}, and \code{package = "matrixCorr"}.
#'
#' @details
#' The biserial correlation is the special two-category case of the polyserial
#' model. It assumes that a binary variable \eqn{Y} arises by thresholding an
#' unobserved standard-normal variable \eqn{Z} that is jointly normal with a
#' continuous variable \eqn{X}. Writing \eqn{p = P(Y = 1)} and
#' \eqn{q = 1-p}, let \eqn{z_p = \Phi^{-1}(p)} and \eqn{\phi(z_p)} be the
#' standard-normal density evaluated at \eqn{z_p}. If \eqn{\bar x_1} and
#' \eqn{\bar x_0} denote the sample means of \eqn{X} in the two observed groups
#' and \eqn{s_x} is the sample standard deviation of \eqn{X}, the usual
#' biserial estimator is
#' \deqn{
#' r_b =
#' \frac{\bar x_1 - \bar x_0}{s_x}
#' \frac{pq}{\phi(z_p)}.
#' }
#' This is exactly the estimator implemented in the underlying C++ kernel.
#'
#' In vector mode a single biserial correlation is returned. In
#' matrix/data-frame mode, every numeric column of \code{data} is paired with every
#' binary column of \code{y}, producing a rectangular matrix of
#' continuous-by-binary biserial correlations.
#'
#' Unlike the point-biserial correlation, which is just Pearson correlation on a
#' 0/1 coding of the binary variable, the biserial coefficient explicitly
#' assumes an underlying latent normal threshold model and rescales the mean
#' difference accordingly.
#'
#' \strong{Computational complexity.} If \code{data} has \eqn{p_x} continuous
#' columns and \code{y} has \eqn{p_y} binary columns, the matrix path computes
#' \eqn{p_x p_y} closed-form estimates with negligible extra memory beyond the
#' output matrix.
#'
#' @references
#' Olsson, U., Drasgow, F., & Dorans, N. J. (1982). The polyserial
#' correlation coefficient. \emph{Psychometrika}, 47(3), 337-347.
#'
#' @examplesIf requireNamespace("mnormt", quietly = TRUE)
#' set.seed(126)
#' n <- 1000
#' Sigma <- matrix(c(
#'   1.00, 0.35, 0.50, 0.25,
#'   0.35, 1.00, 0.30, 0.55,
#'   0.50, 0.30, 1.00, 0.40,
#'   0.25, 0.55, 0.40, 1.00
#' ), 4, 4, byrow = TRUE)
#'
#' Z <- mnormt::rmnorm(n = n, mean = rep(0, 4), varcov = Sigma)
#' X <- data.frame(x1 = Z[, 1], x2 = Z[, 2])
#' Y <- data.frame(
#'   g1 = Z[, 3] > stats::qnorm(0.65),
#'   g2 = Z[, 4] > stats::qnorm(0.55)
#' )
#'
#' bs <- biserial(X, Y)
#' print(bs, digits = 3)
#' summary(bs)
#' plot(bs)
#' @author Thiago de Paula Oliveira
#' @export
biserial <- function(data, y, check_na = TRUE) {
  check_bool(check_na)

  scalar <- is.null(dim(data)) && is.atomic(data) && !is.factor(data) &&
    is.null(dim(y)) && is.atomic(y)

  x_mat <- .mc_extract_continuous_matrix(data, arg = "data", min_cols = 1L)
  y_enc <- .mc_extract_discrete_columns(y, kind = "binary", arg = "y", min_cols = 1L)

  if (nrow(x_mat) != length(y_enc[[1L]]$code)) {
    abort_bad_arg(
      "data",
      message = "and {.arg y} must have the same number of observations."
    )
  }

  .mc_check_latent_missing(
    c(as.list(as.data.frame(x_mat)), lapply(y_enc, `[[`, "code")),
    check_na = check_na,
    arg = "data"
  )

  out <- matrix(NA_real_, nrow = ncol(x_mat), ncol = length(y_enc),
                dimnames = list(colnames(x_mat), names(y_enc)))

  for (j in seq_len(ncol(x_mat))) {
    xj <- x_mat[, j]
    for (k in seq_along(y_enc)) {
      out[j, k] <- .mc_pair_biserial(xj, y_enc[[k]]$code, check_na)
    }
  }

  out <- .mc_scalar_or_matrix(out, scalar = scalar)
  if (scalar) {
    return(out)
  }

  .mc_structure_corr_matrix(
    out,
    class_name = "biserial_corr",
    method = "biserial",
    description = "Biserial correlation matrix (continuous x binary)"
  )
}

#' @rdname biserial
#' @method print biserial_corr
#' @param x An object of class \code{biserial_corr}.
#' @param digits Integer; number of decimal places to print.
#' @param max_rows Optional integer; maximum number of rows to display.
#' @param max_cols Optional integer; maximum number of columns to display.
#' @param ... Additional arguments passed to \code{print()}.
#' @export
print.biserial_corr <- function(x, digits = 4, max_rows = NULL,
                                max_cols = NULL, ...) {
  .mc_print_corr_matrix(
    x, header = "Biserial correlation matrix:",
    digits = digits, max_rows = max_rows, max_cols = max_cols, ...
  )
}

#' @rdname biserial
#' @method plot biserial_corr
#' @param title Plot title. Default is \code{"Biserial correlation heatmap"}.
#' @param low_color Color for the minimum correlation.
#' @param high_color Color for the maximum correlation.
#' @param mid_color Color for zero correlation.
#' @param value_text_size Font size used in tile labels.
#' @export
plot.biserial_corr <- function(x, title = "Biserial correlation heatmap",
                               low_color = "indianred1",
                               high_color = "steelblue1",
                               mid_color = "white",
                               value_text_size = 4, ...) {
  .mc_plot_corr_matrix(
    x, class_name = "biserial_corr", fill_name = "Biserial",
    title = title, low_color = low_color, high_color = high_color,
    mid_color = mid_color, value_text_size = value_text_size, ...
  )
}

#' @rdname biserial
#' @method summary biserial_corr
#' @param object An object of class \code{biserial_corr}.
#' @export
summary.biserial_corr <- function(object, ...) {
  .mc_summary_corr_matrix(object, header = "Latent correlation summary")
}

#' @title Summary Method for Correlation Matrices
#'
#' @description
#' Prints compact summary statistics returned by matrix-style
#' \code{summary()} methods in \pkg{matrixCorr}.
#'
#' @param x An object of class \code{summary_corr_matrix}.
#' @param digits Integer; number of decimal places to print.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.summary_corr_matrix <- function(x, digits = 4, ...) {
  cat(paste0(x$header %||% "Correlation summary", ":\n"))
  cat(sprintf("  class      : %s\n", x$class))
  cat(sprintf("  method     : %s\n", x$method))
  if (!is.null(x$lambda) && is.finite(x$lambda)) {
    cat(sprintf("  lambda     : %s\n", format(signif(x$lambda, digits = digits))))
  }
  if (!is.null(x$rho) && is.finite(x$rho)) {
    cat(sprintf("  rho        : %s\n", format(signif(x$rho, digits = digits))))
  }
  if (!is.null(x$jitter) && is.finite(x$jitter)) {
    cat(sprintf("  jitter     : %s\n", format(signif(x$jitter, digits = digits))))
  }
  cat(sprintf("  dimensions : %d x %d\n", x$n_rows, x$n_cols))
  if (!is.na(x$correct) && !is.null(x$correct)) {
    cat(sprintf("  correct    : %s\n", format(round(x$correct, digits), nsmall = digits)))
  }
  if (!is.na(x$n_variables)) {
    cat(sprintf("  n_variables: %d\n", x$n_variables))
  }
  cat(sprintf("  n_pairs    : %d\n", x$n_pairs))
  cat(sprintf("  symmetric  : %s\n", if (isTRUE(x$symmetric)) "yes" else "no"))
  cat(sprintf("  missing    : %d\n", x$n_missing))
  if (!is.na(x$n_complete_min) && !is.na(x$n_complete_max)) {
    cat(sprintf("  n_complete : %d to %d\n", x$n_complete_min, x$n_complete_max))
  }
  if (!is.na(x$zero_cell_pairs)) {
    cat(sprintf("  zero_cells : %d pair(s)\n", x$zero_cell_pairs))
  }
  if (!is.na(x$corrected_pairs)) {
    cat(sprintf("  corrected  : %d pair(s)\n", x$corrected_pairs))
  }
  if (!is.na(x$boundary_pairs)) {
    cat(sprintf("  boundary   : %d pair(s)\n", x$boundary_pairs))
  }
  if (!is.na(x$near_boundary_pairs)) {
    cat(sprintf("  near_bound : %d pair(s)\n", x$near_boundary_pairs))
  }
  if (!is.na(x$converged_pairs)) {
    cat(sprintf("  converged  : %d pair(s)\n", x$converged_pairs))
  }
  if (!is.na(x$optimizer_tol)) {
    cat(sprintf("  opt_tol    : %s\n", format(signif(x$optimizer_tol, digits = digits))))
  }
  if (!is.na(x$threshold_sets)) {
    cat(sprintf("  thresholds : %d set(s)\n", x$threshold_sets))
  }
  cat(sprintf("  min        : %s\n", format(round(x$estimate_min, digits), nsmall = digits)))
  cat(sprintf("  max        : %s\n", format(round(x$estimate_max, digits), nsmall = digits)))
  invisible(x)
}

#' @title Summary Method for Latent Correlation Matrices
#'
#' @description
#' Prints compact summary statistics returned by
#' \code{summary.tetrachoric_corr()}, \code{summary.polychoric_corr()},
#' \code{summary.polyserial_corr()}, and \code{summary.biserial_corr()}.
#'
#' @param x An object of class \code{summary_latent_corr}.
#' @param digits Integer; number of decimal places to print.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.summary_latent_corr <- function(x, digits = 4, ...) {
  print.summary_corr_matrix(x, digits = digits, ...)
}
