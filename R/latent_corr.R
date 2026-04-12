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
      .hint = "Use `na_method = \"pairwise\"` to use pairwise complete cases."
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

.mc_attach_scalar_latent <- function(x, method, description,
                                     diagnostics = NULL, thresholds = NULL) {
  extra_attrs <- list(estimate = as.numeric(x))
  if (is.list(diagnostics) && length(diagnostics$n_complete) == 1L) {
    extra_attrs$n_obs <- as.integer(diagnostics$n_complete)
  }
  do.call(
    structure,
    c(
      list(
        .Data = as.numeric(x),
        method = method,
        description = description,
        package = "matrixCorr",
        diagnostics = diagnostics,
        thresholds = thresholds
      ),
      extra_attrs
    )
  )
}

.mc_attach_latent_inference <- function(x, payload, conf_level = NULL) {
  attrs <- attributes(x)
  if (is.list(payload$inference)) {
    attrs$inference <- payload$inference
    if (length(payload$inference$estimate) == 1L) {
      attrs$estimate <- as.numeric(payload$inference$estimate)
    }
    if (length(payload$inference$n_obs) == 1L) {
      attrs$n_obs <- as.integer(payload$inference$n_obs)
    }
  }
  if (is.list(payload$ci)) {
    attrs$ci <- payload$ci
    attrs$estimate <- as.numeric(payload$ci$est)
    attrs$lwr <- as.numeric(payload$ci$lwr.ci)
    attrs$upr <- as.numeric(payload$ci$upr.ci)
    attrs$conf.level <- conf_level %||% payload$ci$conf.level
    attrs$conf_level <- conf_level %||% payload$ci$conf.level
  }
  attributes(x) <- attrs
  x
}

.mc_latent_ci_attr <- function(x) {
  attr(x, "ci", exact = TRUE)
}

.mc_latent_inference_attr <- function(x) {
  attr(x, "inference", exact = TRUE)
}

.mc_latent_scalar_payload <- function(fit,
                                      return_ci = FALSE,
                                      return_p = FALSE) {
  list(
    inference = if (isTRUE(return_p)) {
      list(
        method = as.character(fit$inference.method),
        estimate = as.numeric(fit$estimate),
        statistic = as.numeric(fit$statistic),
        parameter = as.numeric(fit$parameter),
        p_value = as.numeric(fit$p_value),
        n_obs = as.integer(fit$n_obs),
        alternative = "two.sided"
      )
    } else {
      NULL
    },
    ci = if (isTRUE(return_ci)) {
      list(
        est = as.numeric(fit$estimate),
        lwr.ci = as.numeric(fit$lwr),
        upr.ci = as.numeric(fit$upr),
        conf.level = as.numeric(fit$conf.level),
        ci.method = as.character(fit$ci.method)
      )
    } else {
      NULL
    }
  )
}

.mc_new_metric_matrix <- function(template, mode = c("double", "integer")) {
  mode <- match.arg(mode)
  matrix(
    if (identical(mode, "integer")) NA_integer_ else NA_real_,
    nrow = nrow(template),
    ncol = ncol(template),
    dimnames = dimnames(template)
  )
}

.mc_latent_matrix_payload <- function(estimate,
                                      statistic = NULL,
                                      parameter = NULL,
                                      p_value = NULL,
                                      n_obs = NULL,
                                      lwr = NULL,
                                      upr = NULL,
                                      conf_level = NULL,
                                      ci_method = NULL,
                                      inference_method = NULL,
                                      return_ci = FALSE,
                                      return_p = FALSE) {
  list(
    inference = if (isTRUE(return_p)) {
      list(
        method = inference_method,
        estimate = estimate,
        statistic = statistic,
        parameter = parameter,
        p_value = p_value,
        n_obs = n_obs,
        alternative = "two.sided"
      )
    } else {
      NULL
    },
    ci = if (isTRUE(return_ci)) {
      list(
        est = estimate,
        lwr.ci = lwr,
        upr.ci = upr,
        conf.level = conf_level,
        ci.method = ci_method
      )
    } else {
      NULL
    }
  )
}

.mc_recenter_wald_fit <- function(fit, estimate, conf_level) {
  fit$estimate <- as.numeric(estimate)
  se <- as.numeric(fit$se)
  est <- as.numeric(fit$estimate)

  if (is.finite(est) && is.finite(se) && se > 0) {
    z_crit <- stats::qnorm(0.5 * (1 + conf_level))
    fit$statistic <- est / se
    fit$p_value <- 2 * stats::pnorm(abs(fit$statistic), lower.tail = FALSE)
    fit$lwr <- max(-1, est - z_crit * se)
    fit$upr <- min(1, est + z_crit * se)
  } else {
    fit$statistic <- NA_real_
    fit$p_value <- NA_real_
    fit$lwr <- NA_real_
    fit$upr <- NA_real_
  }

  fit
}

.mc_tetrachoric_inference_one <- function(tab, correct, conf_level) {
  fit <- matrixCorr_tetrachoric_inference_cpp(
    tab = unclass(tab),
    correct = correct,
    conf_level = conf_level
  )
  .mc_recenter_wald_fit(
    fit,
    estimate = .mc_tetrachoric_table_estimate(tab, correct = correct),
    conf_level = conf_level
  )
}

.mc_polychoric_inference_one <- function(tab, correct, conf_level) {
  fit <- matrixCorr_polychoric_inference_cpp(
    tab = unclass(tab),
    correct = correct,
    conf_level = conf_level
  )
  .mc_recenter_wald_fit(
    fit,
    estimate = .mc_polychoric_table_estimate(tab, correct = correct),
    conf_level = conf_level
  )
}

.mc_polyserial_inference_one <- function(x, y, check_na, conf_level) {
  pair <- .mc_pair_complete(x, y, check_na)
  out <- list(
    estimate = NA_real_,
    statistic = NA_real_,
    parameter = NA_real_,
    p_value = NA_real_,
    n_obs = as.integer(length(pair$x)),
    lwr = NA_real_,
    upr = NA_real_,
    ci.method = "wald_information_polyserial",
    inference.method = "wald_z_polyserial"
  )

  if (length(pair$x) < 2L) {
    return(out)
  }
  yy <- pair$y[!is.na(pair$y)]
  if (length(unique(yy)) < 2L) {
    return(out)
  }

  fit <- matrixCorr_polyserial_inference_cpp(
    x = as.numeric(pair$x),
    y = as.integer(as.factor(pair$y)),
    conf_level = conf_level
  )
  fit$n_obs <- as.integer(fit$n_obs)
  fit
}

.mc_print_corr_matrix <- function(x,
                                  header,
                                  digits = 4,
                                  n = NULL,
                                  topn = NULL,
                                  max_vars = NULL,
                                  width = NULL,
                                  show_ci = NULL,
                                  mat = NULL,
                                  ...) {
  legacy <- .mc_extract_legacy_display_args(
    list(...),
    n = n,
    topn = topn,
    max_vars = max_vars
  )
  cfg <- .mc_resolve_display_args(
    context = "print",
    n = legacy$n,
    topn = legacy$topn,
    max_vars = legacy$max_vars,
    width = width,
    show_ci = show_ci
  )

  m <- as.matrix(.mc_coalesce(mat, x))
  attributes(m) <- attributes(m)[c("dim", "dimnames")]
  ci_present <- .mc_has_ci_payload(x)
  p_present <- .mc_has_p_payload(x)

  cat(header, "\n", sep = "")
  digest <- c(
      method = .mc_coalesce(attr(x, "method", exact = TRUE), NA_character_),
      dimensions = sprintf("%d x %d", nrow(m), ncol(m)),
      if (ci_present) c(ci = "yes"),
      if (p_present) c(p_values = "yes")
    )
  digest <- digest[!is.na(digest) & nzchar(digest)]
  .mc_print_named_digest(digest)
  cat("\n")
  do.call(
    .mc_print_preview_matrix,
    c(
      list(
        m = m,
        digits = digits,
        n = cfg$n,
        topn = cfg$topn,
        max_vars = cfg$max_vars,
        width = cfg$width,
        context = "print",
        full_hint = TRUE,
        summary_hint = TRUE
      ),
      legacy$dots
    )
  )

  invisible(x)
}

.mc_plot_corr_matrix <- function(x, class_name, fill_name, title,
                                 low_color = "indianred1",
                                 high_color = "steelblue1",
                                 mid_color = "white",
                                 value_text_size = 4,
                                 show_value = TRUE, ...) {
  check_inherits(x, class_name)
  check_bool(show_value, arg = "show_value")

  mat <- as.matrix(x)
  df <- as.data.frame(as.table(mat))
  colnames(df) <- c("Var1", "Var2", fill_name)
  df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))

  p <- ggplot2::ggplot(df, ggplot2::aes(Var2, Var1, fill = .data[[fill_name]])) +
    ggplot2::geom_tile(color = "white") +
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

  if (isTRUE(show_value) && !is.null(value_text_size) && is.finite(value_text_size)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", .data[[fill_name]])),
      size = value_text_size,
      color = "black"
    )
  }

  p
}

.mc_summary_corr_matrix <- function(object,
                                    header = "Correlation summary",
                                    topn = NULL) {
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
    },
    has_ci = .mc_has_ci_payload(object),
    has_p = .mc_has_p_payload(object)
  )
  extremes <- .mc_pair_extremes(object, digits = 4)
  out$most_negative <- extremes$most_negative
  out$most_positive <- extremes$most_positive

  diag_attr <- attr(object, "diagnostics")
  if (is.list(diag_attr)) {
    selector <- if (symmetric && nrow(m) > 1L) upper.tri(m, diag = FALSE) else matrix(TRUE, nrow(m), ncol(m))
    pair_pick <- function(x, transform = identity, default = NA_integer_) {
      if (!is.matrix(x) || !identical(dim(x), dim(m))) return(default)
      vals <- transform(x[selector])
      vals <- vals[!is.na(vals)]
      vals <- vals[is.finite(vals)]
      if (!length(vals)) return(default)
      vals
    }

    zero_vals <- pair_pick(diag_attr$zero_cells)
    ncomp_vals <- pair_pick(diag_attr$n_complete)
    skipped_n_vals <- pair_pick(diag_attr$skipped_n, default = integer())
    skipped_prop_vals <- pair_pick(diag_attr$skipped_prop, default = numeric())
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
    out$skipped_n_min <- if (is.atomic(skipped_n_vals) && length(skipped_n_vals)) min(skipped_n_vals, na.rm = TRUE) else NA_integer_
    out$skipped_n_max <- if (is.atomic(skipped_n_vals) && length(skipped_n_vals)) max(skipped_n_vals, na.rm = TRUE) else NA_integer_
    out$skipped_prop_min <- if (is.atomic(skipped_prop_vals) && length(skipped_prop_vals)) min(skipped_prop_vals, na.rm = TRUE) else NA_real_
    out$skipped_prop_max <- if (is.atomic(skipped_prop_vals) && length(skipped_prop_vals)) max(skipped_prop_vals, na.rm = TRUE) else NA_real_
      out$optimizer_tol <- if (is.null(diag_attr$optimizer_tol)) NA_real_ else diag_attr$optimizer_tol
    } else {
    out$zero_cell_pairs <- NA_integer_
    out$boundary_pairs <- NA_integer_
    out$near_boundary_pairs <- NA_integer_
    out$corrected_pairs <- NA_integer_
    out$converged_pairs <- NA_integer_
    out$n_complete_min <- NA_integer_
    out$n_complete_max <- NA_integer_
    out$skipped_n_min <- NA_integer_
    out$skipped_n_max <- NA_integer_
    out$skipped_prop_min <- NA_real_
    out$skipped_prop_max <- NA_real_
      out$optimizer_tol <- NA_real_
    }

  ci_attr <- attr(object, "ci", exact = TRUE)
  out$ci_conf_level <- NA_real_
  out$ci_method <- NA_character_
  out$ci_width_min <- NA_real_
  out$ci_width_max <- NA_real_
  out$ci_cross_zero <- NA_integer_
  if (is.list(ci_attr) &&
      is.matrix(ci_attr$lwr.ci) &&
      is.matrix(ci_attr$upr.ci) &&
      identical(dim(ci_attr$lwr.ci), dim(m)) &&
      identical(dim(ci_attr$upr.ci), dim(m))) {
    selector <- if (symmetric && nrow(m) > 1L) upper.tri(m, diag = FALSE) else matrix(TRUE, nrow(m), ncol(m))
    lwr_vals <- ci_attr$lwr.ci[selector]
    upr_vals <- ci_attr$upr.ci[selector]
    keep_ci <- is.finite(lwr_vals) & is.finite(upr_vals)
    out$ci_conf_level <- suppressWarnings(as.numeric(ci_attr$conf.level %||% attr(object, "conf.level", exact = TRUE)))
    out$ci_method <- as.character(ci_attr$ci.method %||% attr(object, "ci.method", exact = TRUE) %||% NA_character_)
    if (any(keep_ci)) {
      widths <- upr_vals[keep_ci] - lwr_vals[keep_ci]
      out$ci_width_min <- min(widths)
      out$ci_width_max <- max(widths)
      out$ci_cross_zero <- sum(lwr_vals[keep_ci] <= 0 & upr_vals[keep_ci] >= 0, na.rm = TRUE)
    }
  }
  
  out$header <- header
  out$top_results <- .mc_summary_top_pairs(
    object,
    digits = 4,
    topn = .mc_coalesce(topn, .mc_display_option("summary_topn", 5L))
  )
  out <- .mc_finalize_summary_list(
    out,
    class_name = if (identical(header, "Latent correlation summary")) {
      "summary.latent_corr"
    } else {
      "summary.matrixCorr"
    }
  )
  if (identical(class(out)[1L], "summary.matrixCorr")) {
    class(out) <- c("summary.matrixCorr", "summary.corr_matrix")
  } else {
    class(out) <- c("summary.latent_corr", "summary.matrixCorr", "summary.corr_matrix")
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

.mc_tetrachoric_table_estimate <- function(tab, correct) {
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

.mc_polychoric_table_estimate <- function(tab, correct) {
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

.mc_tetrachoric_estimate_internal <- function(x, y, correct,
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

.mc_polychoric_estimate_internal <- function(x, y, n_x, n_y, correct,
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
  .mc_tetrachoric_estimate_internal(
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
  .mc_polychoric_estimate_internal(
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

.mc_biserial_inference_one <- function(x, y, check_na,
                                       return_ci = FALSE,
                                       conf_level = 0.95) {
  pair <- .mc_pair_complete(x, y, check_na)
  n_complete <- length(pair$x)
  out <- list(
    estimate = NA_real_,
    statistic = NA_real_,
    parameter = NA_real_,
    p_value = NA_real_,
    n_obs = as.integer(n_complete),
    lwr = NA_real_,
    upr = NA_real_
  )

  if (n_complete < 2L) {
    return(out)
  }

  yy <- pair$y[!is.na(pair$y)]
  if (length(unique(yy)) < 2L) {
    return(out)
  }

  est <- matrixCorr_biserial_latent_cpp(
    as.numeric(pair$x),
    as.logical(pair$y == 2L)
  )
  out$estimate <- est

  df <- n_complete - 2L
  out$parameter <- if (df > 0L) as.numeric(df) else NA_real_
  if (is.finite(est) && df > 0L) {
    if (abs(est) >= 1) {
      statistic <- sign(est) * Inf
    } else {
      statistic <- est * sqrt(df / (1 - est^2))
    }
    out$statistic <- statistic
    out$p_value <- 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE)
  }

  if (isTRUE(return_ci) && is.finite(est) && n_complete > 3L) {
    if (abs(est) >= 1) {
      out$lwr <- est
      out$upr <- est
    } else {
      z_est <- atanh(est)
      z_crit <- stats::qnorm(0.5 * (1 + conf_level))
      z_half_width <- z_crit / sqrt(n_complete - 3)
      out$lwr <- tanh(z_est - z_half_width)
      out$upr <- tanh(z_est + z_half_width)
    }
  }

  out
}

.mc_biserial_ci_attr <- function(x) {
  attr(x, "ci", exact = TRUE)
}

.mc_biserial_inference_attr <- function(x) {
  attr(x, "inference", exact = TRUE)
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
                                      correct = NULL, dimnames = NULL,
                                      extra_attrs = NULL,
                                      classes = c(class_name, "matrix")) {
  if (!is.null(dimnames)) {
    dimnames(mat) <- dimnames
  }
  do.call(
    structure,
    c(
      list(
        .Data = mat,
        class = classes,
        method = method,
        description = description,
        package = "matrixCorr",
        diagnostics = diagnostics,
        thresholds = thresholds,
        correct = correct
      ),
      extra_attrs
    )
  )
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
#' @param na_method Character scalar controlling missing-data handling.
#'   \code{"error"} rejects missing values. \code{"pairwise"} uses pairwise
#'   complete cases.
#' @param ci Logical (default \code{FALSE}). If \code{TRUE}, attach
#' model-based large-sample Wald confidence intervals derived from the observed
#' information matrix of the latent-variable likelihood.
#' @param p_value Logical (default \code{FALSE}). If \code{TRUE}, attach
#' model-based large-sample Wald p-values and test statistics for each
#' estimated latent correlation.
#' @param conf_level Confidence level used when \code{ci = TRUE}. Default is
#' \code{0.95}.
#'
#' @return
#' If \code{y} is supplied, a numeric scalar with attributes
#' \code{diagnostics} and \code{thresholds}. Otherwise a symmetric matrix of
#' class \code{tetrachoric_corr} with attributes \code{method},
#' \code{description}, \code{package = "matrixCorr"}, \code{diagnostics},
#' \code{thresholds}, and \code{correct}. When \code{p_value = TRUE}, the
#' returned object also carries an \code{inference} attribute with elements
#' \code{estimate}, \code{statistic}, \code{parameter}, \code{p_value}, and
#' \code{n_obs}. When \code{ci = TRUE}, it also carries a \code{ci} attribute
#' with elements \code{est}, \code{lwr.ci}, \code{upr.ci}, \code{conf.level},
#' and \code{ci.method}, plus \code{attr(x, "conf.level")}. Scalar outputs keep
#' the same point estimate and gain the same metadata only when inference is
#' requested.
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
#' \code{correct = 0.5} continuity-correction behaviour used in several
#' latent-correlation implementations.
#' When \code{correct = 0} and the observed contingency table contains zero
#' cells, the fit is non-regular and may be boundary-driven. In those cases the
#' returned object stores sparse-fit diagnostics, including whether the fit was
#' classified as \code{boundary} or \code{near_boundary}.
#'
#' \strong{Assumptions.} The coefficient is appropriate when both observed
#' binary variables are viewed as thresholded versions of jointly normal latent
#' variables. The optional p-values and confidence intervals adopt this
#' latent-normal interpretation and use the same likelihood that defines the
#' tetrachoric estimate. These inferential quantities are therefore model-based
#' and should not be interpreted as distribution-free summaries.
#'
#' \strong{Inference.} When \code{ci = TRUE} or \code{p_value = TRUE}, the
#' function refits the pairwise tetrachoric model by maximum likelihood and
#' obtains the observed information matrix numerically in C++. The reported
#' confidence interval is a Wald interval
#' \eqn{\hat\rho \pm z_{1-\alpha/2}\operatorname{SE}(\hat\rho)}, and the
#' reported p-value is from the large-sample Wald \eqn{z}-test for
#' \eqn{H_0:\rho = 0}. These inferential quantities are only computed when
#' explicitly requested.
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
tetrachoric <- function(data,
                        y = NULL,
                        na_method = c("error", "pairwise"),
                        ci = FALSE,
                        p_value = FALSE,
                        conf_level = 0.95,
                        correct = 0.5,
                        ...) {
  if (missing(na_method) && isFALSE(ci) && isFALSE(p_value)) {
    check_scalar_nonneg(correct, arg = "correct")

    if (is.null(y) && inherits(data, "table")) {
      tab <- as.matrix(data)
      est <- .mc_tetrachoric_table_estimate(tab, correct = correct)
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
      .mc_check_latent_missing(list(data = data, y = y), check_na = TRUE, arg = "data")

      enc_x <- .mc_encode_ordinal_vector(data, binary = TRUE)
      enc_y <- .mc_encode_ordinal_vector(y, binary = TRUE)
      if (is.null(enc_x) || is.null(enc_y)) {
        abort_bad_arg(
          "data",
          message = "and {.arg y} must both be binary (two-level) variables."
        )
      }

      pair <- .mc_pair_complete(enc_x$code, enc_y$code, TRUE)
      est <- .mc_pair_tetrachoric(pair$x, pair$y, correct, TRUE)
      x01 <- as.integer(pair$x) - min(as.integer(pair$x), na.rm = TRUE)
      y01 <- as.integer(pair$y) - min(as.integer(pair$y), na.rm = TRUE)
      tab <- .mc_fast_binary_table01(x01, y01)
      thresholds <- list(data = .mc_binary_tau(x01), y = .mc_binary_tau(y01))
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
    .mc_check_latent_missing(lapply(enc, `[[`, "code"), check_na = TRUE, arg = "data")

    x_mat <- .mc_encode_list_to_matrix(enc)
    tau <- vapply(enc, function(z) {
      code <- as.integer(z$code)
      .mc_binary_tau(code - min(code, na.rm = TRUE))
    }, numeric(1))
    out <- matrixCorr_tetrachoric_matrix_cpp(
      x = x_mat,
      tau = tau,
      correct = correct,
      pairwise_complete = FALSE
    )
    dimnames(out) <- list(names(enc), names(enc))

    p <- length(enc)
    diag_info <- .mc_matrix_diagnostics_template(p, dimnames(out))
    diag(diag_info$zero_cells) <- 0L
    diag(diag_info$n_complete) <- colSums(!is.na(x_mat))
    diag(diag_info$corrected) <- FALSE
    diag(diag_info$boundary) <- FALSE
    diag(diag_info$near_boundary) <- FALSE
    diag(diag_info$converged) <- is.finite(diag(out))
    for (j in seq_len(p - 1L)) {
      for (k in (j + 1L):p) {
        pair <- .mc_pair_complete(enc[[j]]$code, enc[[k]]$code, TRUE)
        x01 <- as.integer(pair$x) - min(as.integer(pair$x), na.rm = TRUE)
        y01 <- as.integer(pair$y) - min(as.integer(pair$y), na.rm = TRUE)
        tab <- .mc_fast_binary_table01(x01, y01)
        diag_info <- .mc_fill_pair_diag(diag_info, j, k, out[j, k], tab, length(pair$x), correct)
      }
    }
    thresholds <- as.list(tau)
    names(thresholds) <- names(enc)
    return(.mc_structure_corr_matrix(
      out,
      class_name = "tetrachoric_corr",
      method = "tetrachoric",
      description = "Pairwise tetrachoric correlation matrix",
      diagnostics = diag_info,
      thresholds = thresholds,
      correct = correct
    ))
  }

  legacy_args <- .mc_extract_legacy_aliases(list(...), allowed = "check_na")
  na_cfg <- resolve_na_args(
    na_method = na_method,
    check_na = legacy_args$check_na %||% NULL,
    na_method_missing = missing(na_method)
  )
  check_na <- na_cfg$check_na
  check_bool(ci, arg = "ci")
  check_bool(p_value, arg = "p_value")
  check_scalar_nonneg(correct, arg = "correct")
  if (isTRUE(ci)) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  }

  if (is.null(y) && inherits(data, "table")) {
    tab <- as.matrix(data)
    fit <- if (isTRUE(ci) || isTRUE(p_value)) {
      .mc_tetrachoric_inference_one(tab, correct = correct, conf_level = conf_level)
    } else {
      NULL
    }
    est <- if (is.null(fit)) .mc_tetrachoric_table_estimate(tab, correct = correct) else as.numeric(fit$estimate)
    thresholds <- .mc_tetra_table_thresholds(tab, correct = correct)
    diagnostics <- .mc_latent_pair_metadata(
      est = est,
      tab = tab,
      n_complete = sum(tab),
      correct = correct,
      thresholds = thresholds
    )
    out <- .mc_attach_scalar_latent(
      est,
      method = "tetrachoric",
      description = "Pairwise tetrachoric correlation",
      diagnostics = diagnostics,
      thresholds = thresholds
    )
    if (isTRUE(ci) || isTRUE(p_value)) {
      out <- .mc_attach_latent_inference(
        out,
        .mc_latent_scalar_payload(fit, return_ci = ci, return_p = p_value),
        conf_level = conf_level
      )
    }
    return(out)
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
    x01 <- as.integer(pair$x) - min(as.integer(pair$x), na.rm = TRUE)
    y01 <- as.integer(pair$y) - min(as.integer(pair$y), na.rm = TRUE)
    tab <- .mc_fast_binary_table01(x01, y01)
    fit <- if (isTRUE(ci) || isTRUE(p_value)) {
      .mc_tetrachoric_inference_one(tab, correct = correct, conf_level = conf_level)
    } else {
      NULL
    }
    est <- .mc_pair_tetrachoric(pair$x, pair$y, correct, TRUE)
    if (!is.null(fit)) {
      fit <- .mc_recenter_wald_fit(fit, estimate = est, conf_level = conf_level)
    }
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
    out <- .mc_attach_scalar_latent(
      est,
      method = "tetrachoric",
      description = "Pairwise tetrachoric correlation",
      diagnostics = diagnostics,
      thresholds = thresholds
    )
    if (isTRUE(ci) || isTRUE(p_value)) {
      out <- .mc_attach_latent_inference(
        out,
        .mc_latent_scalar_payload(fit, return_ci = ci, return_p = p_value),
        conf_level = conf_level
      )
    }
    return(out)
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
  out <- if (!isTRUE(ci) && !isTRUE(p_value)) {
    matrixCorr_tetrachoric_matrix_cpp(
      x = x_mat,
      tau = tau,
      correct = correct,
      pairwise_complete = !check_na
    )
  } else {
    matrix(NA_real_, nrow = p, ncol = p)
  }
  dimnames(out) <- list(names(enc), names(enc))
  diag_info <- .mc_matrix_diagnostics_template(p, dimnames(out))
  diag(diag_info$zero_cells) <- 0L
  diag(diag_info$n_complete) <- colSums(!is.na(x_mat))
  diag(diag_info$corrected) <- FALSE
  diag(diag_info$boundary) <- FALSE
  diag(diag_info$near_boundary) <- FALSE
  diag(diag_info$converged) <- TRUE
  diag(out) <- ifelse(vapply(enc, function(z) length(unique(z$code[!is.na(z$code)])) >= 2L, logical(1)), 1, NA_real_)
  statistic <- parameter <- p_mat <- lwr <- upr <- NULL
  if (isTRUE(ci) || isTRUE(p_value)) {
    statistic <- if (isTRUE(p_value)) .mc_new_metric_matrix(out, "double") else NULL
    parameter <- if (isTRUE(p_value)) .mc_new_metric_matrix(out, "double") else NULL
    p_mat <- if (isTRUE(p_value)) .mc_new_metric_matrix(out, "double") else NULL
    lwr <- if (isTRUE(ci)) .mc_new_metric_matrix(out, "double") else NULL
    upr <- if (isTRUE(ci)) .mc_new_metric_matrix(out, "double") else NULL
  }
  for (j in seq_len(p - 1L)) {
    for (k in (j + 1L):p) {
      pair <- .mc_pair_complete(enc[[j]]$code, enc[[k]]$code, check_na)
      x01 <- as.integer(pair$x) - min(as.integer(pair$x), na.rm = TRUE)
      y01 <- as.integer(pair$y) - min(as.integer(pair$y), na.rm = TRUE)
      tab <- .mc_fast_binary_table01(x01, y01)
      est_jk <- .mc_pair_tetrachoric(
        enc[[j]]$code, enc[[k]]$code,
        correct = correct,
        check_na = check_na,
        tau_x = tau[j],
        tau_y = tau[k],
        global = TRUE
      )
      if (isTRUE(ci) || isTRUE(p_value)) {
        fit <- .mc_tetrachoric_inference_one(tab, correct = correct, conf_level = conf_level)
        fit <- .mc_recenter_wald_fit(fit, estimate = est_jk, conf_level = conf_level)
        out[j, k] <- est_jk
        out[k, j] <- out[j, k]
        if (isTRUE(p_value)) {
          statistic[j, k] <- as.numeric(fit$statistic)
          statistic[k, j] <- statistic[j, k]
          parameter[j, k] <- as.numeric(fit$parameter)
          parameter[k, j] <- parameter[j, k]
          p_mat[j, k] <- as.numeric(fit$p_value)
          p_mat[k, j] <- p_mat[j, k]
        }
        if (isTRUE(ci)) {
          lwr[j, k] <- as.numeric(fit$lwr)
          lwr[k, j] <- lwr[j, k]
          upr[j, k] <- as.numeric(fit$upr)
          upr[k, j] <- upr[j, k]
        }
      }
      diag_info <- .mc_fill_pair_diag(diag_info, j, k, out[j, k], tab, length(pair$x), correct)
    }
  }
  thresholds <- as.list(tau)
  names(thresholds) <- names(enc)

  out <- .mc_structure_corr_matrix(
    out,
    class_name = "tetrachoric_corr",
    method = "tetrachoric",
    description = "Pairwise tetrachoric correlation matrix",
    diagnostics = diag_info,
    thresholds = thresholds,
    correct = correct
  )
  if (isTRUE(ci) || isTRUE(p_value)) {
    out <- .mc_attach_latent_inference(
      out,
      .mc_latent_matrix_payload(
        estimate = unclass(out),
        statistic = statistic,
        parameter = parameter,
        p_value = p_mat,
        n_obs = diag_info$n_complete,
        lwr = lwr,
        upr = upr,
        conf_level = conf_level,
        ci_method = "wald_information_tetrachoric",
        inference_method = "wald_z_tetrachoric",
        return_ci = ci,
        return_p = p_value
      ),
      conf_level = conf_level
    )
  }
  out
}

#' @rdname tetrachoric
#' @method print tetrachoric_corr
#' @param x An object of class \code{tetrachoric_corr}.
#' @param digits Integer; number of decimal places to print.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param ... Additional arguments passed to \code{print()}.
#' @export
print.tetrachoric_corr <- function(x, digits = 4, n = NULL, topn = NULL,
                                   max_vars = NULL, width = NULL,
                                   show_ci = NULL, ...) {
  .mc_print_corr_matrix(
    x, header = "Tetrachoric correlation matrix",
    digits = digits, n = n, topn = topn,
    max_vars = max_vars, width = width, show_ci = show_ci, ...
  )
}

#' @rdname tetrachoric
#' @method plot tetrachoric_corr
#' @param title Plot title. Default is \code{"Tetrachoric correlation heatmap"}.
#' @param low_color Color for the minimum correlation.
#' @param high_color Color for the maximum correlation.
#' @param mid_color Color for zero correlation.
#' @param value_text_size Font size used in tile labels.
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
#' @export
plot.tetrachoric_corr <- function(x, title = "Tetrachoric correlation heatmap",
                                  low_color = "indianred1",
                                  high_color = "steelblue1",
                                  mid_color = "white",
                                  value_text_size = 4,
                                  show_value = TRUE, ...) {
  .mc_plot_corr_matrix(
    x, class_name = "tetrachoric_corr", fill_name = "Tetrachoric",
    title = title, low_color = low_color, high_color = high_color,
    mid_color = mid_color, value_text_size = value_text_size,
    show_value = show_value, ...
  )
}

#' @rdname tetrachoric
#' @method summary tetrachoric_corr
#' @param object An object of class \code{tetrachoric_corr}.
#' @param ci_digits Integer; digits for confidence limits in the pairwise
#'   summary.
#' @param p_digits Integer; digits for p-values in the pairwise summary.
#' @export
summary.tetrachoric_corr <- function(object, n = NULL, topn = NULL,
                                     max_vars = NULL, width = NULL,
                                     ci_digits = 3,
                                     p_digits = 4,
                                     show_ci = NULL, ...) {
  check_inherits(object, "tetrachoric_corr")
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )

  ci <- .mc_latent_ci_attr(object)
  inf <- .mc_latent_inference_attr(object)
  if (is.null(ci) && (is.null(inf) || is.null(inf$p_value))) {
    return(.mc_summary_corr_matrix(object, header = "Latent correlation summary", topn = topn))
  }

  .mc_latent_pairwise_summary(
    object,
    summary_class = "summary.tetrachoric_corr",
    ci_digits = ci_digits,
    p_digits = p_digits,
    show_ci = show_ci
  )
}

#' @rdname tetrachoric
#' @method print summary.tetrachoric_corr
#' @param x An object of class \code{summary.tetrachoric_corr}.
#' @export
print.summary.tetrachoric_corr <- function(x, digits = NULL, n = NULL,
                                           topn = NULL, max_vars = NULL,
                                           width = NULL, show_ci = NULL, ...) {
  .mc_print_latent_pairwise_summary(
    x,
    title = "Tetrachoric correlation summary",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
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
#' @param na_method Character scalar controlling missing-data handling.
#'   \code{"error"} rejects missing values. \code{"pairwise"} uses pairwise
#'   complete cases.
#' @param ci Logical (default \code{FALSE}). If \code{TRUE}, attach
#' model-based large-sample Wald confidence intervals derived from the observed
#' information matrix of the latent-variable likelihood.
#' @param p_value Logical (default \code{FALSE}). If \code{TRUE}, attach
#' model-based large-sample Wald p-values and test statistics for each
#' estimated latent correlation.
#' @param conf_level Confidence level used when \code{ci = TRUE}. Default is
#' \code{0.95}.
#'
#' @return
#' If \code{y} is supplied, a numeric scalar with attributes
#' \code{diagnostics} and \code{thresholds}. Otherwise a symmetric matrix of
#' class \code{polychoric_corr} with attributes \code{method},
#' \code{description}, \code{package = "matrixCorr"}, \code{diagnostics},
#' \code{thresholds}, and \code{correct}. When \code{p_value = TRUE}, the
#' returned object also carries an \code{inference} attribute with elements
#' \code{estimate}, \code{statistic}, \code{parameter}, \code{p_value}, and
#' \code{n_obs}. When \code{ci = TRUE}, it also carries a \code{ci} attribute
#' with elements \code{est}, \code{lwr.ci}, \code{upr.ci}, \code{conf.level},
#' and \code{ci.method}, plus \code{attr(x, "conf.level")}. Scalar outputs keep
#' the same point estimate and gain the same metadata only when inference is
#' requested.
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
#' \strong{Assumptions.} The coefficient is appropriate when both observed
#' ordinal variables are viewed as discretisations of jointly normal latent
#' variables. The optional p-values and confidence intervals adopt this
#' latent-normal interpretation and use the same likelihood that defines the
#' polychoric estimate. These inferential quantities are therefore model-based
#' and should not be interpreted as distribution-free summaries.
#'
#' \strong{Inference.} When \code{ci = TRUE} or \code{p_value = TRUE}, the
#' function refits the pairwise polychoric model by maximum likelihood and
#' obtains the observed information matrix numerically in C++. The reported
#' confidence interval is a Wald interval
#' \eqn{\hat\rho \pm z_{1-\alpha/2}\operatorname{SE}(\hat\rho)}, and the
#' reported p-value is from the large-sample Wald \eqn{z}-test for
#' \eqn{H_0:\rho = 0}. These inferential quantities are only computed when
#' explicitly requested.
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
polychoric <- function(data,
                       y = NULL,
                       na_method = c("error", "pairwise"),
                       ci = FALSE,
                       p_value = FALSE,
                       conf_level = 0.95,
                       correct = 0.5,
                       ...) {
  if (missing(na_method) && isFALSE(ci) && isFALSE(p_value)) {
    check_scalar_nonneg(correct, arg = "correct")

    if (is.null(y) && inherits(data, "table")) {
      tab <- as.matrix(data)
      est <- .mc_polychoric_table_estimate(tab, correct = correct)
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
      .mc_check_latent_missing(list(data = data, y = y), check_na = TRUE, arg = "data")

      enc_x <- .mc_encode_ordinal_vector(data, binary = FALSE)
      enc_y <- .mc_encode_ordinal_vector(y, binary = FALSE)
      if (is.null(enc_x) || is.null(enc_y)) {
        abort_bad_arg("data", message = "and {.arg y} must both be ordinal variables.")
      }
      n_x <- length(enc_x$levels)
      n_y <- length(enc_y$levels)
      global <- identical(n_x, n_y)
      pair <- .mc_pair_complete(enc_x$code, enc_y$code, TRUE)
      est <- .mc_pair_polychoric(
        pair$x, pair$y,
        n_x = n_x,
        n_y = n_y,
        correct = correct,
        check_na = TRUE,
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
    .mc_check_latent_missing(lapply(enc, `[[`, "code"), check_na = TRUE, arg = "data")

    p <- length(enc)
    n_levels <- vapply(enc, function(z) length(z$levels), integer(1))
    global_all <- length(unique(n_levels)) == 1L
    tau_mat <- if (global_all) .mc_tau_matrix(enc) else matrix(NA_real_, nrow = 0L, ncol = 0L)
    x_mat <- .mc_encode_list_to_matrix(enc)
    out <- matrixCorr_polychoric_matrix_cpp(
      x = x_mat,
      n_levels = n_levels,
      tau_mat = tau_mat,
      global_all = global_all,
      correct = correct,
      pairwise_complete = FALSE
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
        pair <- .mc_pair_complete(enc[[j]]$code, enc[[k]]$code, TRUE)
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
    return(.mc_structure_corr_matrix(
      out,
      class_name = "polychoric_corr",
      method = "polychoric",
      description = "Pairwise polychoric correlation matrix",
      diagnostics = diag_info,
      thresholds = thresholds,
      correct = correct
    ))
  }

  legacy_args <- .mc_extract_legacy_aliases(list(...), allowed = "check_na")
  na_cfg <- resolve_na_args(
    na_method = na_method,
    check_na = legacy_args$check_na %||% NULL,
    na_method_missing = missing(na_method)
  )
  check_na <- na_cfg$check_na
  check_bool(ci, arg = "ci")
  check_bool(p_value, arg = "p_value")
  check_scalar_nonneg(correct, arg = "correct")
  if (isTRUE(ci)) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  }

  if (is.null(y) && inherits(data, "table")) {
    tab <- as.matrix(data)
    fit <- if (isTRUE(ci) || isTRUE(p_value)) {
      .mc_polychoric_inference_one(tab, correct = correct, conf_level = conf_level)
    } else {
      NULL
    }
    est <- if (is.null(fit)) .mc_polychoric_table_estimate(tab, correct = correct) else as.numeric(fit$estimate)
    thresholds <- .mc_poly_table_thresholds(tab, correct = correct)
    diagnostics <- .mc_latent_pair_metadata(
      est = est,
      tab = tab,
      n_complete = sum(tab),
      correct = correct,
      thresholds = thresholds
    )
    out <- .mc_attach_scalar_latent(
      est,
      method = "polychoric",
      description = "Pairwise polychoric correlation",
      diagnostics = diagnostics,
      thresholds = thresholds
    )
    if (isTRUE(ci) || isTRUE(p_value)) {
      out <- .mc_attach_latent_inference(
        out,
        .mc_latent_scalar_payload(fit, return_ci = ci, return_p = p_value),
        conf_level = conf_level
      )
    }
    return(out)
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
    tab <- .mc_fast_ordinal_table(pair$x, pair$y, n_x = n_x, n_y = n_y)
    fit <- if (isTRUE(ci) || isTRUE(p_value)) {
      .mc_polychoric_inference_one(tab, correct = correct, conf_level = conf_level)
    } else {
      NULL
    }
    est <- .mc_pair_polychoric(
      pair$x, pair$y,
      n_x = n_x, n_y = n_y,
      correct = correct, check_na = TRUE,
      tau_x = if (global) .mc_global_cutpoints(pair$x, n_x) else NULL,
      tau_y = if (global) .mc_global_cutpoints(pair$y, n_y) else NULL,
      global = global
    )
    if (!is.null(fit)) {
      fit <- .mc_recenter_wald_fit(fit, estimate = est, conf_level = conf_level)
    }
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
    out <- .mc_attach_scalar_latent(
      est,
      method = "polychoric",
      description = "Pairwise polychoric correlation",
      diagnostics = diagnostics,
      thresholds = thresholds
    )
    if (isTRUE(ci) || isTRUE(p_value)) {
      out <- .mc_attach_latent_inference(
        out,
        .mc_latent_scalar_payload(fit, return_ci = ci, return_p = p_value),
        conf_level = conf_level
      )
    }
    return(out)
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
  global_thresholds <- lapply(enc, function(z) .mc_global_cutpoints(z$code, length(z$levels)))
  names(global_thresholds) <- names(enc)
  x_mat <- .mc_encode_list_to_matrix(enc)
  out <- if (!isTRUE(ci) && !isTRUE(p_value)) {
    matrixCorr_polychoric_matrix_cpp(
      x = x_mat,
      n_levels = n_levels,
      tau_mat = tau_mat,
      global_all = global_all,
      correct = correct,
      pairwise_complete = !check_na
    )
  } else {
    matrix(NA_real_, nrow = p, ncol = p)
  }
  dimnames(out) <- list(names(enc), names(enc))
  diag_info <- .mc_matrix_diagnostics_template(p, dimnames(out))
  diag(diag_info$zero_cells) <- 0L
  diag(diag_info$n_complete) <- vapply(enc, function(z) sum(!is.na(z$code)), integer(1))
  diag(diag_info$corrected) <- FALSE
  diag(diag_info$boundary) <- FALSE
  diag(diag_info$near_boundary) <- FALSE
  diag(diag_info$converged) <- TRUE
  diag(out) <- ifelse(vapply(enc, function(z) length(unique(z$code[!is.na(z$code)])) >= 2L, logical(1)), 1, NA_real_)
  statistic <- parameter <- p_mat <- lwr <- upr <- NULL
  if (isTRUE(ci) || isTRUE(p_value)) {
    statistic <- if (isTRUE(p_value)) .mc_new_metric_matrix(out, "double") else NULL
    parameter <- if (isTRUE(p_value)) .mc_new_metric_matrix(out, "double") else NULL
    p_mat <- if (isTRUE(p_value)) .mc_new_metric_matrix(out, "double") else NULL
    lwr <- if (isTRUE(ci)) .mc_new_metric_matrix(out, "double") else NULL
    upr <- if (isTRUE(ci)) .mc_new_metric_matrix(out, "double") else NULL
  }
  for (j in seq_len(p - 1L)) {
    for (k in (j + 1L):p) {
      pair <- .mc_pair_complete(enc[[j]]$code, enc[[k]]$code, check_na)
      tab <- .mc_fast_ordinal_table(
        pair$x, pair$y,
        n_x = length(enc[[j]]$levels),
        n_y = length(enc[[k]]$levels)
      )
      global <- if (global_all) TRUE else identical(n_levels[j], n_levels[k])
      est_jk <- .mc_pair_polychoric(
        enc[[j]]$code, enc[[k]]$code,
        n_x = n_levels[j],
        n_y = n_levels[k],
        correct = correct,
        check_na = check_na,
        tau_x = if (global) global_thresholds[[j]] else NULL,
        tau_y = if (global) global_thresholds[[k]] else NULL,
        global = global
      )
      if (isTRUE(ci) || isTRUE(p_value)) {
        fit <- .mc_polychoric_inference_one(tab, correct = correct, conf_level = conf_level)
        fit <- .mc_recenter_wald_fit(fit, estimate = est_jk, conf_level = conf_level)
        out[j, k] <- est_jk
        out[k, j] <- out[j, k]
        if (isTRUE(p_value)) {
          statistic[j, k] <- as.numeric(fit$statistic)
          statistic[k, j] <- statistic[j, k]
          parameter[j, k] <- as.numeric(fit$parameter)
          parameter[k, j] <- parameter[j, k]
          p_mat[j, k] <- as.numeric(fit$p_value)
          p_mat[k, j] <- p_mat[j, k]
        }
        if (isTRUE(ci)) {
          lwr[j, k] <- as.numeric(fit$lwr)
          lwr[k, j] <- lwr[j, k]
          upr[j, k] <- as.numeric(fit$upr)
          upr[k, j] <- upr[j, k]
        }
      }
      diag_info <- .mc_fill_pair_diag(diag_info, j, k, out[j, k], tab, length(pair$x), correct)
    }
  }
  thresholds <- global_thresholds

  out <- .mc_structure_corr_matrix(
    out,
    class_name = "polychoric_corr",
    method = "polychoric",
    description = "Pairwise polychoric correlation matrix",
    diagnostics = diag_info,
    thresholds = thresholds,
    correct = correct
  )
  if (isTRUE(ci) || isTRUE(p_value)) {
    out <- .mc_attach_latent_inference(
      out,
      .mc_latent_matrix_payload(
        estimate = unclass(out),
        statistic = statistic,
        parameter = parameter,
        p_value = p_mat,
        n_obs = diag_info$n_complete,
        lwr = lwr,
        upr = upr,
        conf_level = conf_level,
        ci_method = "wald_information_polychoric",
        inference_method = "wald_z_polychoric",
        return_ci = ci,
        return_p = p_value
      ),
      conf_level = conf_level
    )
  }
  out
}

#' @rdname polychoric
#' @method print polychoric_corr
#' @param x An object of class \code{polychoric_corr}.
#' @param digits Integer; number of decimal places to print.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param ... Additional arguments passed to \code{print()}.
#' @export
print.polychoric_corr <- function(x, digits = 4, n = NULL, topn = NULL,
                                  max_vars = NULL, width = NULL,
                                  show_ci = NULL, ...) {
  .mc_print_corr_matrix(
    x, header = "Polychoric correlation matrix",
    digits = digits, n = n, topn = topn,
    max_vars = max_vars, width = width, show_ci = show_ci, ...
  )
}

#' @rdname polychoric
#' @method plot polychoric_corr
#' @param title Plot title. Default is \code{"Polychoric correlation heatmap"}.
#' @param low_color Color for the minimum correlation.
#' @param high_color Color for the maximum correlation.
#' @param mid_color Color for zero correlation.
#' @param value_text_size Font size used in tile labels.
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
#' @export
plot.polychoric_corr <- function(x, title = "Polychoric correlation heatmap",
                                 low_color = "indianred1",
                                 high_color = "steelblue1",
                                 mid_color = "white",
                                 value_text_size = 4,
                                 show_value = TRUE, ...) {
  .mc_plot_corr_matrix(
    x, class_name = "polychoric_corr", fill_name = "Polychoric",
    title = title, low_color = low_color, high_color = high_color,
    mid_color = mid_color, value_text_size = value_text_size,
    show_value = show_value, ...
  )
}

#' @rdname polychoric
#' @method summary polychoric_corr
#' @param object An object of class \code{polychoric_corr}.
#' @param ci_digits Integer; digits for confidence limits in the pairwise
#'   summary.
#' @param p_digits Integer; digits for p-values in the pairwise summary.
#' @export
summary.polychoric_corr <- function(object, n = NULL, topn = NULL,
                                    max_vars = NULL, width = NULL,
                                    ci_digits = 3,
                                    p_digits = 4,
                                    show_ci = NULL, ...) {
  check_inherits(object, "polychoric_corr")
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )

  ci <- .mc_latent_ci_attr(object)
  inf <- .mc_latent_inference_attr(object)
  if (is.null(ci) && (is.null(inf) || is.null(inf$p_value))) {
    return(.mc_summary_corr_matrix(object, header = "Latent correlation summary", topn = topn))
  }

  .mc_latent_pairwise_summary(
    object,
    summary_class = "summary.polychoric_corr",
    ci_digits = ci_digits,
    p_digits = p_digits,
    show_ci = show_ci
  )
}

#' @rdname polychoric
#' @method print summary.polychoric_corr
#' @param x An object of class \code{summary.polychoric_corr}.
#' @export
print.summary.polychoric_corr <- function(x, digits = NULL, n = NULL,
                                          topn = NULL, max_vars = NULL,
                                          width = NULL, show_ci = NULL, ...) {
  .mc_print_latent_pairwise_summary(
    x,
    title = "Polychoric correlation summary",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
}

#' @title Polyserial Correlation Between Continuous and Ordinal Variables
#'
#' @description
#' Computes polyserial correlations between continuous variables in \code{data}
#' and ordinal variables in \code{y}. Both pairwise vector mode and rectangular
#' matrix/data-frame mode are supported.
#'
#' @usage polyserial(data, y, na_method = c("error", "pairwise"), ci = FALSE, p_value = FALSE,
#'   conf_level = 0.95, ...)
#'
#' @param data A numeric vector, matrix, or data frame containing continuous
#' variables.
#' @param y An ordinal vector, matrix, or data frame containing ordinal
#' variables. Supported columns are factors, ordered factors, logical values,
#' or integer-like numerics.
#' @param na_method Character scalar controlling missing-data handling.
#'   \code{"error"} rejects missing values. \code{"pairwise"} uses pairwise
#'   complete cases.
#' @param ci Logical (default \code{FALSE}). If \code{TRUE}, attach
#' model-based large-sample Wald confidence intervals derived from the observed
#' information matrix of the latent-variable likelihood.
#' @param p_value Logical (default \code{FALSE}). If \code{TRUE}, attach
#' model-based large-sample Wald p-values and test statistics for each
#' estimated latent correlation.
#' @param conf_level Confidence level used when \code{ci = TRUE}. Default is
#' \code{0.95}.
#'
#' @return
#' If both \code{data} and \code{y} are vectors, a numeric scalar. Otherwise a
#' numeric matrix of class \code{polyserial_corr} with rows corresponding to
#' the continuous variables in \code{data} and columns to the ordinal variables
#' in \code{y}. Matrix outputs carry attributes \code{method},
#' \code{description}, and \code{package = "matrixCorr"}. When
#' \code{p_value = TRUE}, the returned object also carries an \code{inference}
#' attribute with elements \code{estimate}, \code{statistic}, \code{parameter},
#' \code{p_value}, and \code{n_obs}. When \code{ci = TRUE}, it also carries a
#' \code{ci} attribute with elements \code{est}, \code{lwr.ci},
#' \code{upr.ci}, \code{conf.level}, and \code{ci.method}, plus
#' \code{attr(x, "conf.level")}. Scalar outputs keep the same point estimate
#' and gain the same metadata only when inference is requested.
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
#' \strong{Assumptions.} The coefficient is appropriate when the ordinal
#' variable is viewed as the discretised version of a latent normal variable
#' that is jointly normal with the observed continuous variable. The optional
#' p-values and confidence intervals adopt this latent-normal interpretation
#' and use the same likelihood that defines the polyserial estimate. These
#' inferential quantities are therefore model-based and should not be
#' interpreted as distribution-free summaries.
#'
#' \strong{Inference.} When \code{ci = TRUE} or \code{p_value = TRUE}, the
#' function refits the pairwise polyserial model by maximum likelihood and
#' obtains the observed information matrix numerically in C++. The reported
#' confidence interval is a Wald interval
#' \eqn{\hat\rho \pm z_{1-\alpha/2}\operatorname{SE}(\hat\rho)}, and the
#' reported p-value is from the large-sample Wald \eqn{z}-test for
#' \eqn{H_0:\rho = 0}. These inferential quantities are only computed when
#' explicitly requested.
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
#' @param ... Deprecated compatibility aliases. The legacy \code{check_na}
#'   argument is still accepted temporarily.
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
polyserial <- function(data,
                       y,
                       na_method = c("error", "pairwise"),
                       ci = FALSE,
                       p_value = FALSE,
                       conf_level = 0.95,
                       ...) {
  if (missing(na_method) && isFALSE(ci) && isFALSE(p_value)) {
    scalar <- is.null(dim(data)) && is.atomic(data) && !is.factor(data) &&
      is.null(dim(y)) && is.atomic(y)

    x_mat <- .mc_extract_continuous_matrix(data, arg = "data", min_cols = 1L)
    y_enc <- .mc_extract_discrete_columns(y, kind = "ordinal", arg = "y", min_cols = 1L)
    if (nrow(x_mat) != length(y_enc[[1L]]$code)) {
      abort_bad_arg("data", message = "and {.arg y} must have the same number of observations.")
    }
    .mc_check_latent_missing(
      c(as.list(as.data.frame(x_mat)), lapply(y_enc, `[[`, "code")),
      check_na = TRUE,
      arg = "data"
    )

    if (scalar) {
      return(.mc_pair_polyserial(x_mat[, 1L], y_enc[[1L]]$code, TRUE))
    }

    out <- matrix(
      NA_real_,
      nrow = ncol(x_mat),
      ncol = length(y_enc),
      dimnames = list(colnames(x_mat), names(y_enc))
    )
    for (j in seq_len(ncol(x_mat))) {
      xj <- x_mat[, j]
      for (k in seq_along(y_enc)) {
        out[j, k] <- .mc_pair_polyserial(xj, y_enc[[k]]$code, TRUE)
      }
    }
    return(.mc_structure_corr_matrix(
      out,
      class_name = "polyserial_corr",
      method = "polyserial",
      description = "Polyserial correlation matrix (continuous x ordinal)"
    ))
  }

  legacy_args <- .mc_extract_legacy_aliases(list(...), allowed = "check_na")
  na_cfg <- resolve_na_args(
    na_method = na_method,
    check_na = legacy_args$check_na %||% NULL,
    na_method_missing = missing(na_method)
  )
  check_na <- na_cfg$check_na
  check_bool(ci, arg = "ci")
  check_bool(p_value, arg = "p_value")
  if (isTRUE(ci)) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  }

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
    if (!isTRUE(ci) && !isTRUE(p_value)) {
      return(.mc_pair_polyserial(x_mat[, 1L], y_enc[[1L]]$code, check_na))
    }
    fit <- .mc_polyserial_inference_one(
      x_mat[, 1L],
      y_enc[[1L]]$code,
      check_na = check_na,
      conf_level = conf_level
    )
    out <- .mc_attach_scalar_latent(
      fit$estimate,
      method = "polyserial",
      description = "Polyserial correlation",
      diagnostics = list(n_complete = as.integer(fit$n_obs))
    )
    out <- .mc_attach_latent_inference(
      out,
      .mc_latent_scalar_payload(fit, return_ci = ci, return_p = p_value),
      conf_level = conf_level
    )
    return(out)
  }

  out <- matrix(NA_real_, nrow = ncol(x_mat), ncol = length(y_enc),
                dimnames = list(colnames(x_mat), names(y_enc)))
  n_complete <- NULL
  if (isTRUE(ci) || isTRUE(p_value)) {
    n_complete <- matrix(
      NA_integer_,
      nrow = nrow(out),
      ncol = ncol(out),
      dimnames = dimnames(out)
    )
  }
  statistic <- parameter <- p_mat <- lwr <- upr <- NULL
  if (isTRUE(ci) || isTRUE(p_value)) {
    statistic <- if (isTRUE(p_value)) .mc_new_metric_matrix(out, "double") else NULL
    parameter <- if (isTRUE(p_value)) .mc_new_metric_matrix(out, "double") else NULL
    p_mat <- if (isTRUE(p_value)) .mc_new_metric_matrix(out, "double") else NULL
    lwr <- if (isTRUE(ci)) .mc_new_metric_matrix(out, "double") else NULL
    upr <- if (isTRUE(ci)) .mc_new_metric_matrix(out, "double") else NULL
  }

  for (j in seq_len(ncol(x_mat))) {
    xj <- x_mat[, j]
    for (k in seq_along(y_enc)) {
      if (isTRUE(ci) || isTRUE(p_value)) {
        fit <- .mc_polyserial_inference_one(
          xj,
          y_enc[[k]]$code,
          check_na = check_na,
          conf_level = conf_level
        )
        out[j, k] <- as.numeric(fit$estimate)
        n_complete[j, k] <- as.integer(fit$n_obs)
        if (isTRUE(p_value)) {
          statistic[j, k] <- as.numeric(fit$statistic)
          parameter[j, k] <- as.numeric(fit$parameter)
          p_mat[j, k] <- as.numeric(fit$p_value)
        }
        if (isTRUE(ci)) {
          lwr[j, k] <- as.numeric(fit$lwr)
          upr[j, k] <- as.numeric(fit$upr)
        }
      } else {
        out[j, k] <- .mc_pair_polyserial(xj, y_enc[[k]]$code, check_na)
      }
    }
  }

  out <- .mc_structure_corr_matrix(
    out,
    class_name = "polyserial_corr",
    method = "polyserial",
    description = "Polyserial correlation matrix (continuous x ordinal)",
    diagnostics = if (isTRUE(ci) || isTRUE(p_value)) list(n_complete = n_complete) else NULL
  )
  if (isTRUE(ci) || isTRUE(p_value)) {
    out <- .mc_attach_latent_inference(
      out,
      .mc_latent_matrix_payload(
        estimate = unclass(out),
        statistic = statistic,
        parameter = parameter,
        p_value = p_mat,
        n_obs = n_complete,
        lwr = lwr,
        upr = upr,
        conf_level = conf_level,
        ci_method = "wald_information_polyserial",
        inference_method = "wald_z_polyserial",
        return_ci = ci,
        return_p = p_value
      ),
      conf_level = conf_level
    )
  }
  out
}

#' @rdname polyserial
#' @method print polyserial_corr
#' @param x An object of class \code{polyserial_corr}.
#' @param digits Integer; number of decimal places to print.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param ... Additional arguments passed to \code{print()}.
#' @export
print.polyserial_corr <- function(x, digits = 4, n = NULL, topn = NULL,
                                  max_vars = NULL, width = NULL,
                                  show_ci = NULL, ...) {
  .mc_print_corr_matrix(
    x, header = "Polyserial correlation matrix",
    digits = digits, n = n, topn = topn,
    max_vars = max_vars, width = width, show_ci = show_ci, ...
  )
}

#' @rdname polyserial
#' @method plot polyserial_corr
#' @param title Plot title. Default is \code{"Polyserial correlation heatmap"}.
#' @param low_color Color for the minimum correlation.
#' @param high_color Color for the maximum correlation.
#' @param mid_color Color for zero correlation.
#' @param value_text_size Font size used in tile labels.
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
#' @export
plot.polyserial_corr <- function(x, title = "Polyserial correlation heatmap",
                                 low_color = "indianred1",
                                 high_color = "steelblue1",
                                 mid_color = "white",
                                 value_text_size = 4,
                                 show_value = TRUE, ...) {
  .mc_plot_corr_matrix(
    x, class_name = "polyserial_corr", fill_name = "Polyserial",
    title = title, low_color = low_color, high_color = high_color,
    mid_color = mid_color, value_text_size = value_text_size,
    show_value = show_value, ...
  )
}

#' @rdname polyserial
#' @method summary polyserial_corr
#' @param object An object of class \code{polyserial_corr}.
#' @param ci_digits Integer; digits for confidence limits in the pairwise
#'   summary.
#' @param p_digits Integer; digits for p-values in the pairwise summary.
#' @export
summary.polyserial_corr <- function(object, n = NULL, topn = NULL,
                                    max_vars = NULL, width = NULL,
                                    ci_digits = 3,
                                    p_digits = 4,
                                    show_ci = NULL, ...) {
  check_inherits(object, "polyserial_corr")
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )

  ci <- .mc_latent_ci_attr(object)
  inf <- .mc_latent_inference_attr(object)
  if (is.null(ci) && (is.null(inf) || is.null(inf$p_value))) {
    return(.mc_summary_corr_matrix(object, header = "Latent correlation summary", topn = topn))
  }

  .mc_latent_pairwise_summary(
    object,
    summary_class = "summary.polyserial_corr",
    ci_digits = ci_digits,
    p_digits = p_digits,
    show_ci = show_ci
  )
}

#' @rdname polyserial
#' @method print summary.polyserial_corr
#' @param x An object of class \code{summary.polyserial_corr}.
#' @export
print.summary.polyserial_corr <- function(x, digits = NULL, n = NULL,
                                          topn = NULL, max_vars = NULL,
                                          width = NULL, show_ci = NULL, ...) {
  .mc_print_latent_pairwise_summary(
    x,
    title = "Polyserial correlation summary",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
}

#' @title Biserial Correlation Between Continuous and Binary Variables
#'
#' @description
#' Computes biserial correlations between continuous variables in \code{data}
#' and binary variables in \code{y}. Both pairwise vector mode and rectangular
#' matrix/data-frame mode are supported.
#'
#' @usage biserial(data, y, na_method = c("error", "pairwise"), ci = FALSE, p_value = FALSE,
#'   conf_level = 0.95, ...)
#'
#' @param data A numeric vector, matrix, or data frame containing continuous
#' variables.
#' @param y A binary vector, matrix, or data frame. In data-frame mode, only
#' two-level columns are retained.
#' @param na_method Character scalar controlling missing-data handling.
#'   \code{"error"} rejects missing values. \code{"pairwise"} uses pairwise
#'   complete cases.
#' @param ci Logical (default \code{FALSE}). If \code{TRUE}, attach
#' approximate large-sample confidence intervals derived from a Fisher
#' \eqn{z}-transformation of the biserial estimate.
#' @param p_value Logical (default \code{FALSE}). If \code{TRUE}, attach
#' model-based large-sample p-values, test statistics, and degrees of freedom
#' for each biserial estimate.
#' @param conf_level Confidence level used when \code{ci = TRUE}. Default is
#' \code{0.95}.
#'
#' @return
#' If both \code{data} and \code{y} are vectors, a numeric scalar. Otherwise a
#' numeric matrix of class \code{biserial_corr} with rows corresponding to
#' the continuous variables in \code{data} and columns to the binary variables
#' in \code{y}. Matrix outputs carry attributes \code{method},
#' \code{description}, and \code{package = "matrixCorr"}. When
#' \code{p_value = TRUE}, the object also carries an \code{inference}
#' attribute with matrices \code{estimate}, \code{statistic},
#' \code{parameter}, \code{p_value}, and \code{n_obs}. When \code{ci = TRUE},
#' it additionally carries a \code{ci} attribute with matrices
#' \code{lwr.ci} and \code{upr.ci}, plus \code{attr(x, "conf.level")}. Scalar
#' outputs keep the same point estimate and gain the same metadata only when
#' inference is requested.
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
#' \strong{Assumptions.} The biserial coefficient is appropriate when the
#' observed binary variable is viewed as a thresholded version of an
#' unobserved continuous latent variable that is jointly normal with the
#' observed continuous variable. The optional p-values and confidence
#' intervals adopt this latent-normal interpretation together with the usual
#' large-sample approximations used for correlation coefficients. These
#' inferential quantities are therefore model-based and should not be
#' interpreted as distribution-free summaries.
#'
#' \strong{Inference.} When \code{p_value = TRUE}, the package reports the
#' large-sample \eqn{t}-statistic
#' \deqn{
#' t = r_b \sqrt{\frac{n - 2}{1 - r_b^2}},
#' }
#' referenced to a Student \eqn{t}-distribution with \eqn{n - 2} degrees of
#' freedom. When \code{ci = TRUE}, the package forms an approximate Fisher
#' \eqn{z}-interval by transforming \eqn{r_b} with
#' \eqn{z = \operatorname{atanh}(r_b)}, using standard error
#' \eqn{1 / \sqrt{n - 3}}, and mapping the limits back with
#' \eqn{\tanh(\cdot)}. The CI is therefore an internal large-sample
#' extension and is only computed when explicitly requested.
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
#' Fisher, R. A. (1921). On the probable error of a coefficient of
#' correlation deduced from a small sample. \emph{Metron}, 1, 3-32.
#' @param ... Deprecated compatibility aliases. The legacy \code{check_na}
#'   argument is still accepted temporarily.
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
#' bs <- biserial(X, Y, ci = TRUE, p_value = TRUE)
#' print(bs, digits = 3)
#' summary(bs)
#' plot(bs)
#' @author Thiago de Paula Oliveira
#' @export
biserial <- function(data,
                     y,
                     na_method = c("error", "pairwise"),
                     ci = FALSE,
                     p_value = FALSE,
                     conf_level = 0.95,
                     ...) {
  if (missing(na_method) && isFALSE(ci) && isFALSE(p_value)) {
    scalar <- is.null(dim(data)) && is.atomic(data) && !is.factor(data) &&
      is.null(dim(y)) && is.atomic(y)

    x_mat <- .mc_extract_continuous_matrix(data, arg = "data", min_cols = 1L)
    y_enc <- .mc_extract_discrete_columns(y, kind = "binary", arg = "y", min_cols = 1L)
    if (nrow(x_mat) != length(y_enc[[1L]]$code)) {
      abort_bad_arg("data", message = "and {.arg y} must have the same number of observations.")
    }
    .mc_check_latent_missing(
      c(as.list(as.data.frame(x_mat)), lapply(y_enc, `[[`, "code")),
      check_na = TRUE,
      arg = "data"
    )

    out <- matrix(
      NA_real_,
      nrow = ncol(x_mat),
      ncol = length(y_enc),
      dimnames = list(colnames(x_mat), names(y_enc))
    )
    for (j in seq_len(ncol(x_mat))) {
      xj <- x_mat[, j]
      for (k in seq_along(y_enc)) {
        out[j, k] <- .mc_pair_biserial(xj, y_enc[[k]]$code, TRUE)
      }
    }
    out <- .mc_scalar_or_matrix(out, scalar = scalar)
    if (scalar) {
      return(out)
    }
    return(.mc_structure_corr_matrix(
      out,
      class_name = "biserial_corr",
      method = "biserial",
      description = "Biserial correlation matrix (continuous x binary)"
    ))
  }

  legacy_args <- .mc_extract_legacy_aliases(list(...), allowed = "check_na")
  na_cfg <- resolve_na_args(
    na_method = na_method,
    check_na = legacy_args$check_na %||% NULL,
    na_method_missing = missing(na_method)
  )
  check_na <- na_cfg$check_na
  check_bool(ci, arg = "ci")
  check_bool(p_value, arg = "p_value")
  if (isTRUE(ci)) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  }

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
  diagnostics <- NULL
  inference_attr <- NULL
  ci_attr <- NULL

  if (isTRUE(ci) || isTRUE(p_value)) {
    n_complete <- matrix(
      NA_integer_,
      nrow = nrow(out),
      ncol = ncol(out),
      dimnames = dimnames(out)
    )
    statistic <- if (isTRUE(p_value)) matrix(
      NA_real_,
      nrow = nrow(out),
      ncol = ncol(out),
      dimnames = dimnames(out)
    ) else NULL
    parameter <- if (isTRUE(p_value)) matrix(
      NA_real_,
      nrow = nrow(out),
      ncol = ncol(out),
      dimnames = dimnames(out)
    ) else NULL
    p_mat <- if (isTRUE(p_value)) matrix(
      NA_real_,
      nrow = nrow(out),
      ncol = ncol(out),
      dimnames = dimnames(out)
    ) else NULL
    lwr <- if (isTRUE(ci)) matrix(
      NA_real_,
      nrow = nrow(out),
      ncol = ncol(out),
      dimnames = dimnames(out)
    ) else NULL
    upr <- if (isTRUE(ci)) matrix(
      NA_real_,
      nrow = nrow(out),
      ncol = ncol(out),
      dimnames = dimnames(out)
    ) else NULL
  }

  for (j in seq_len(ncol(x_mat))) {
    xj <- x_mat[, j]
    for (k in seq_along(y_enc)) {
      if (isTRUE(ci) || isTRUE(p_value)) {
        fit <- .mc_biserial_inference_one(
          xj,
          y_enc[[k]]$code,
          check_na = check_na,
          return_ci = ci,
          conf_level = conf_level
        )
        out[j, k] <- fit$estimate
        n_complete[j, k] <- fit$n_obs
        if (isTRUE(p_value)) {
          statistic[j, k] <- fit$statistic
          parameter[j, k] <- fit$parameter
          p_mat[j, k] <- fit$p_value
        }
        if (isTRUE(ci)) {
          lwr[j, k] <- fit$lwr
          upr[j, k] <- fit$upr
        }
      } else {
        out[j, k] <- .mc_pair_biserial(xj, y_enc[[k]]$code, check_na)
      }
    }
  }

  estimate_mat <- out
  attributes(estimate_mat) <- attributes(estimate_mat)[c("dim", "dimnames")]
  out <- .mc_scalar_or_matrix(out, scalar = scalar)
  if (scalar) {
    if (!isTRUE(ci) && !isTRUE(p_value)) {
      return(out)
    }
    diagnostics <- list(n_complete = as.integer(n_complete[1L, 1L]))
    scalar_out <- .mc_attach_scalar_latent(
      out,
      method = "biserial",
      description = "Biserial correlation",
      diagnostics = diagnostics
    )
    if (isTRUE(p_value)) {
      attr(scalar_out, "inference") <- list(
        method = "biserial_t_test",
        estimate = as.numeric(out),
        statistic = as.numeric(statistic[1L, 1L]),
        parameter = as.numeric(parameter[1L, 1L]),
        p_value = as.numeric(p_mat[1L, 1L]),
        n_obs = as.integer(n_complete[1L, 1L]),
        alternative = "two.sided"
      )
    }
    if (isTRUE(ci)) {
      attr(scalar_out, "ci") <- list(
        est = as.numeric(out),
        lwr.ci = as.numeric(lwr[1L, 1L]),
        upr.ci = as.numeric(upr[1L, 1L]),
        conf.level = conf_level,
        ci.method = "fisher_z_biserial"
      )
      attr(scalar_out, "conf.level") <- conf_level
    }
    return(scalar_out)
  }

  diagnostics <- if (isTRUE(ci) || isTRUE(p_value)) {
    list(n_complete = n_complete)
  } else {
    NULL
  }

  out <- .mc_structure_corr_matrix(
    out,
    class_name = "biserial_corr",
    method = "biserial",
    description = "Biserial correlation matrix (continuous x binary)",
    diagnostics = diagnostics
  )
    if (isTRUE(p_value)) {
      inference_attr <- list(
        method = "biserial_t_test",
        estimate = estimate_mat,
        statistic = statistic,
        parameter = parameter,
        p_value = p_mat,
      n_obs = n_complete,
      alternative = "two.sided"
    )
    attr(out, "inference") <- inference_attr
  }
  if (isTRUE(ci)) {
    ci_attr <- list(
      est = estimate_mat,
      lwr.ci = lwr,
      upr.ci = upr,
      conf.level = conf_level,
      ci.method = "fisher_z_biserial"
    )
    attr(out, "ci") <- ci_attr
    attr(out, "conf.level") <- conf_level
  }
  out
}

#' @rdname biserial
#' @method print biserial_corr
#' @param x An object of class \code{biserial_corr}.
#' @param digits Integer; number of decimal places to print.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param ... Additional arguments passed to \code{print()}.
#' @export
print.biserial_corr <- function(x, digits = 4, n = NULL, topn = NULL,
                                max_vars = NULL, width = NULL,
                                show_ci = NULL, ...) {
  .mc_print_corr_matrix(
    x, header = "Biserial correlation matrix",
    digits = digits, n = n, topn = topn,
    max_vars = max_vars, width = width, show_ci = show_ci, ...
  )
}

#' @rdname biserial
#' @method plot biserial_corr
#' @param title Plot title. Default is \code{"Biserial correlation heatmap"}.
#' @param low_color Color for the minimum correlation.
#' @param high_color Color for the maximum correlation.
#' @param mid_color Color for zero correlation.
#' @param value_text_size Font size used in tile labels.
#' @param ci_text_size Text size for confidence intervals in the heatmap.
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
#' @export
plot.biserial_corr <- function(x, title = "Biserial correlation heatmap",
                               low_color = "indianred1",
                               high_color = "steelblue1",
                               mid_color = "white",
                               value_text_size = 4,
                               ci_text_size = 3,
                               show_value = TRUE, ...) {
  check_bool(show_value, arg = "show_value")
  ci <- .mc_biserial_ci_attr(x)
  if (is.null(ci) || is.null(ci$lwr.ci) || is.null(ci$upr.ci)) {
    return(.mc_plot_corr_matrix(
      x, class_name = "biserial_corr", fill_name = "Biserial",
      title = title, low_color = low_color, high_color = high_color,
      mid_color = mid_color, value_text_size = value_text_size,
      show_value = show_value, ...
    ))
  }

  est_mat <- as.matrix(x)
  df_est <- as.data.frame(as.table(est_mat))
  names(df_est) <- c("Var1", "Var2", "biserial")

  df_lwr <- as.data.frame(as.table(ci$lwr.ci))
  names(df_lwr)[3L] <- "lwr"
  df_upr <- as.data.frame(as.table(ci$upr.ci))
  names(df_upr)[3L] <- "upr"
  df <- Reduce(
    function(a, b) merge(a, b, by = c("Var1", "Var2"), all = TRUE),
    list(df_est, df_lwr, df_upr)
  )

  lev_row <- unique(df_est$Var1)
  lev_col <- unique(df_est$Var2)
  df$Var1 <- factor(df$Var1, levels = rev(lev_row))
  df$Var2 <- factor(df$Var2, levels = lev_col)
  df$ci_label <- ifelse(
    is.na(df$lwr) | is.na(df$upr),
    NA_character_,
    sprintf("[%.3f, %.3f]", df$lwr, df$upr)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Var2, y = Var1, fill = .data$biserial)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(
      low = low_color,
      high = high_color,
      mid = mid_color,
      midpoint = 0,
      limits = c(-1, 1),
      name = "Biserial"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      ...
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = NULL, y = NULL)

  if (isTRUE(show_value) && !is.null(value_text_size) && is.finite(value_text_size)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.na(biserial), "NA", sprintf("%.2f", biserial))),
      size = value_text_size,
      color = "black"
    )
  }
  if (isTRUE(show_value) && any(!is.na(df$ci_label))) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = ci_label, y = as.numeric(Var1) - 0.22),
      size = ci_text_size,
      color = "gray30",
      na.rm = TRUE
    )
  }

  p
}

.mc_latent_pairwise_summary <- function(object,
                                        summary_class,
                                        digits = 4,
                                        ci_digits = 3,
                                        p_digits = 4,
                                        show_ci = "yes") {
  est <- as.matrix(object)
  rn <- rownames(est)
  cn <- colnames(est)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(est)))
  if (is.null(cn)) cn <- as.character(seq_len(ncol(est)))

  ci <- .mc_latent_ci_attr(object)
  inf <- .mc_latent_inference_attr(object)
  diag_attr <- attr(object, "diagnostics", exact = TRUE)
  include_ci <- identical(show_ci, "yes") && is.list(ci)
  include_p <- is.list(inf) && !is.null(inf$p_value)
  symmetric <- isTRUE(nrow(est) == ncol(est)) && isTRUE(isSymmetric(unclass(est)))

  n_pairs <- if (symmetric) {
    nrow(est) * (ncol(est) - 1L) / 2L
  } else {
    nrow(est) * ncol(est)
  }
  rows <- vector("list", n_pairs)
  k <- 0L

  for (i in seq_len(nrow(est))) {
    j_start <- if (symmetric) i + 1L else 1L
    if (j_start > ncol(est)) {
      next
    }
    for (j in j_start:ncol(est)) {
      k <- k + 1L
      rec <- list(
        var1 = rn[i],
        var2 = cn[j],
        estimate = round(est[i, j], digits)
      )
      if (is.list(diag_attr) && is.matrix(diag_attr$n_complete)) {
        rec$n_complete <- as.integer(diag_attr$n_complete[i, j])
      }
      if (include_ci) {
        rec$lwr <- if (is.finite(ci$lwr.ci[i, j])) round(ci$lwr.ci[i, j], ci_digits) else NA_real_
        rec$upr <- if (is.finite(ci$upr.ci[i, j])) round(ci$upr.ci[i, j], ci_digits) else NA_real_
      }
      if (include_p) {
        rec$statistic <- if (is.finite(inf$statistic[i, j])) round(inf$statistic[i, j], digits) else NA_real_
        rec$df <- if (is.finite(inf$parameter[i, j])) round(inf$parameter[i, j], digits) else NA_real_
        rec$p_value <- if (is.finite(inf$p_value[i, j])) round(inf$p_value[i, j], p_digits) else NA_real_
      }
      rows[[k]] <- rec
    }
  }

  if (!k) {
    df <- data.frame(var1 = character(), var2 = character(), estimate = numeric())
  } else {
    rows <- rows[seq_len(k)]
    df <- do.call(rbind.data.frame, rows)
  }
  rownames(df) <- NULL
  num_cols <- intersect(c("estimate", "lwr", "upr", "statistic", "df", "p_value"), names(df))
  int_cols <- intersect(c("n_complete"), names(df))
  for (nm in num_cols) df[[nm]] <- as.numeric(df[[nm]])
  for (nm in int_cols) df[[nm]] <- as.integer(df[[nm]])

  out <- .mc_finalize_summary_df(df, class_name = summary_class)
  attr(out, "overview") <- .mc_summary_corr_matrix(object, header = "Latent correlation summary")
  attr(out, "has_ci") <- include_ci
  attr(out, "has_p") <- include_p
  attr(out, "conf.level") <- if (is.null(ci)) NA_real_ else ci$conf.level
  attr(out, "ci_method") <- if (is.null(ci)) NA_character_ else ci$ci.method
  attr(out, "digits") <- digits
  attr(out, "ci_digits") <- ci_digits
  attr(out, "p_digits") <- p_digits
  attr(out, "inference_method") <- if (is.null(inf)) NA_character_ else inf$method
  out
}

.mc_print_latent_pairwise_summary <- function(x,
                                              title,
                                              digits = NULL,
                                              n = NULL,
                                              topn = NULL,
                                              max_vars = NULL,
                                              width = NULL,
                                              show_ci = NULL,
                                              ...) {
  .mc_print_pairwise_summary_digest(
    x,
    title = title,
    digits = .mc_coalesce(digits, 4),
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ci_method = if (isTRUE(attr(x, "has_ci"))) attr(x, "ci_method", exact = TRUE) else NULL,
    extra_items = c(inference = attr(x, "inference_method", exact = TRUE)),
    ...
  )
  invisible(x)
}

.mc_biserial_pairwise_summary <- function(object,
                                          digits = 4,
                                          ci_digits = 3,
                                          p_digits = 4,
                                          show_ci = "yes") {
  check_inherits(object, "biserial_corr")

  est <- as.matrix(object)
  rn <- rownames(est)
  cn <- colnames(est)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(est)))
  if (is.null(cn)) cn <- as.character(seq_len(ncol(est)))

  ci <- .mc_biserial_ci_attr(object)
  inf <- .mc_biserial_inference_attr(object)
  diag_attr <- attr(object, "diagnostics", exact = TRUE)
  include_ci <- identical(show_ci, "yes") && is.list(ci)
  include_p <- is.list(inf) && !is.null(inf$p_value)

  rows <- vector("list", nrow(est) * ncol(est))
  k <- 0L
  for (i in seq_len(nrow(est))) {
    for (j in seq_len(ncol(est))) {
      k <- k + 1L
      rec <- list(
        var1 = rn[i],
        var2 = cn[j],
        estimate = round(est[i, j], digits)
      )
      if (is.list(diag_attr) && is.matrix(diag_attr$n_complete)) {
        rec$n_complete <- as.integer(diag_attr$n_complete[i, j])
      }
      if (include_ci) {
        rec$lwr <- if (is.finite(ci$lwr.ci[i, j])) round(ci$lwr.ci[i, j], ci_digits) else NA_real_
        rec$upr <- if (is.finite(ci$upr.ci[i, j])) round(ci$upr.ci[i, j], ci_digits) else NA_real_
      }
      if (include_p) {
        rec$statistic <- if (is.finite(inf$statistic[i, j])) round(inf$statistic[i, j], digits) else NA_real_
        rec$df <- if (is.finite(inf$parameter[i, j])) round(inf$parameter[i, j], digits) else NA_real_
        rec$p_value <- if (is.finite(inf$p_value[i, j])) round(inf$p_value[i, j], p_digits) else NA_real_
      }
      rows[[k]] <- rec
    }
  }

  df <- do.call(rbind.data.frame, rows)
  rownames(df) <- NULL
  num_cols <- intersect(c("estimate", "lwr", "upr", "statistic", "df", "p_value"), names(df))
  int_cols <- intersect(c("n_complete"), names(df))
  for (nm in num_cols) df[[nm]] <- as.numeric(df[[nm]])
  for (nm in int_cols) df[[nm]] <- as.integer(df[[nm]])

  out <- .mc_finalize_summary_df(df, class_name = "summary.biserial_corr")
  attr(out, "overview") <- .mc_summary_corr_matrix(object)
  attr(out, "has_ci") <- include_ci
  attr(out, "has_p") <- include_p
  attr(out, "conf.level") <- if (is.null(ci)) NA_real_ else ci$conf.level
  attr(out, "digits") <- digits
  attr(out, "ci_digits") <- ci_digits
  attr(out, "p_digits") <- p_digits
  attr(out, "inference_method") <- if (is.null(inf)) NA_character_ else inf$method
  out
}

#' @rdname biserial
#' @method summary biserial_corr
#' @param object An object of class \code{biserial_corr}.
#' @param ci_digits Integer; digits for biserial confidence limits in the
#'   pairwise summary.
#' @param p_digits Integer; digits for biserial p-values in the pairwise
#'   summary.
#' @export
summary.biserial_corr <- function(object, n = NULL, topn = NULL,
                                  max_vars = NULL, width = NULL,
                                  ci_digits = 3,
                                  p_digits = 4,
                                  show_ci = NULL, ...) {
  check_inherits(object, "biserial_corr")
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )

  ci <- .mc_biserial_ci_attr(object)
  inf <- .mc_biserial_inference_attr(object)
  if (is.null(ci) && (is.null(inf) || is.null(inf$p_value))) {
    return(.mc_summary_corr_matrix(object, header = "Latent correlation summary", topn = topn))
  }

  .mc_biserial_pairwise_summary(
    object,
    ci_digits = ci_digits,
    p_digits = p_digits,
    show_ci = show_ci
  )
}

#' @rdname biserial
#' @method print summary.biserial_corr
#' @param x An object of class \code{summary.biserial_corr}.
#' @export
print.summary.biserial_corr <- function(x, digits = NULL, n = NULL,
                                        topn = NULL, max_vars = NULL,
                                        width = NULL, show_ci = NULL, ...) {
  .mc_print_pairwise_summary_digest(
    x,
    title = "Biserial correlation summary",
    digits = .mc_coalesce(digits, 4),
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ci_method = if (isTRUE(attr(x, "has_ci"))) "fisher_z_biserial" else NULL,
    extra_items = c(inference = attr(x, "inference_method", exact = TRUE)),
    ...
  )
  invisible(x)
}

#' @title Summary Method for Correlation Matrices
#'
#' @description
#' Prints compact summary statistics returned by matrix-style
#' \code{summary()} methods in \pkg{matrixCorr}.
#'
#' @param x An object of class \code{summary.matrixCorr}.
#' @param digits Integer; number of decimal places to print.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#' @method print summary.matrixCorr
#' @export
print.summary.matrixCorr <- function(x,
                                     digits = 4,
                                     n = NULL,
                                     topn = NULL,
                                     max_vars = NULL,
                                     width = NULL,
                                     show_ci = NULL,
                                     ...) {
  cfg <- .mc_resolve_display_args(
    context = "summary",
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci
  )

  digest <- .mc_corr_summary_digest_items(x, digits = digits, show_ci = cfg$show_ci)
  if (isTRUE(x$has_ci) && identical(cfg$show_ci, "yes")) {
    if (is.finite(x$ci_conf_level)) {
      digest <- c(digest, ci_level = sprintf("%g%%", 100 * x$ci_conf_level))
    }
    if (!is.null(x$ci_method) && length(x$ci_method) == 1L &&
        !is.na(x$ci_method) && nzchar(x$ci_method)) {
      digest <- c(digest, ci_method = x$ci_method)
    }
    ci_width_txt <- .mc_format_scalar_or_range(
      x$ci_width_min,
      x$ci_width_max,
      digits = 3,
      integer = FALSE
    )
    if (!is.null(ci_width_txt)) {
      digest <- c(digest, ci_width = ci_width_txt)
    }
    if (is.finite(x$ci_cross_zero) && x$ci_cross_zero > 0L) {
      digest <- c(digest, cross_zero = sprintf("%s pair(s)", .mc_count_fmt(x$ci_cross_zero)))
    }
  }
  if (!is.null(x$lambda) && is.finite(x$lambda)) {
    digest <- c(digest, lambda = format(signif(x$lambda, digits = digits)))
  }
  if (!is.null(x$rho) && is.finite(x$rho)) {
    digest <- c(digest, rho = format(signif(x$rho, digits = digits)))
  }
  if (!is.null(x$jitter) && is.finite(x$jitter)) {
    digest <- c(digest, jitter = format(signif(x$jitter, digits = digits)))
  }

  .mc_print_named_digest(
    digest,
    header = .mc_coalesce(x$header, "Correlation summary")
  )

  if (is.data.frame(x$top_results) && nrow(x$top_results)) {
    .mc_print_ranked_pairs_preview(
      x$top_results,
      header = "Strongest pairs by |estimate|",
      topn = cfg$topn,
      max_vars = cfg$max_vars,
      width = cfg$width,
      show_ci = cfg$show_ci,
      ...
    )
  }
  invisible(x)
}

#' @method print summary.corr_matrix
#' @export
print.summary.corr_matrix <- function(x,
                                      digits = 4,
                                      n = NULL,
                                      topn = NULL,
                                      max_vars = NULL,
                                      width = NULL,
                                      show_ci = NULL,
                                      ...) {
  print.summary.matrixCorr(
    x,
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
}

print.summary_corr_matrix <- print.summary.corr_matrix

#' @title Summary Method for Latent Correlation Matrices
#'
#' @description
#' Prints compact summary statistics returned by
#' \code{summary.tetrachoric_corr()}, \code{summary.polychoric_corr()},
#' \code{summary.polyserial_corr()}, and \code{summary.biserial_corr()}.
#'
#' @param x An object of class \code{summary.latent_corr}.
#' @param digits Integer; number of decimal places to print.
#' @param ... Unused.
#'
#' @return Invisibly returns \code{x}.
#' @method print summary.latent_corr
#' @export
print.summary.latent_corr <- function(x, digits = 4, ...) {
  print.summary.matrixCorr(x, digits = digits, ...)
}

print.summary_latent_corr <- print.summary.latent_corr
