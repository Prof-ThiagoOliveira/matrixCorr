#' @title Pairwise Lin's concordance correlation coefficient
#'
#' @description
#' Computes all pairwise Lin's Concordance Correlation Coefficients (CCC)
#' from the numeric columns of a matrix or data frame. CCC measures both
#' precision (Pearson correlation) and accuracy (closeness to the 45-degree line).
#' This function is backed by a high-performance 'C++' implementation.
#'
#' Lin's CCC quantifies the concordance between a new test/measurement
#' and a gold-standard for the same variable. Like a correlation, CCC
#' ranges from -1 to 1 with perfect agreement at 1, and it cannot exceed the
#' absolute value of the Pearson correlation between variables. It can be
#' legitimately computed even with small samples (e.g., 10 observations),
#' and results are often similar to intraclass correlation coefficients.
#' CCC provides a single summary of agreement, but it may not capture
#' systematic bias; a Bland-Altman plot (differences vs. means) is recommended
#' to visualize bias, proportional trends, and heteroscedasticity (see
#' \code{\link{ba}}).
#'
#' @details
#' Lin's CCC is defined as
#' \deqn{
#' \rho_c \;=\; \frac{2\,\mathrm{cov}(X, Y)}
#'                  {\sigma_X^2 + \sigma_Y^2 + (\mu_X - \mu_Y)^2},
#' }
#' where \eqn{\mu_X,\mu_Y} are the means, \eqn{\sigma_X^2,\sigma_Y^2} the
#' variances, and \eqn{\mathrm{cov}(X,Y)} the covariance. Equivalently,
#' \deqn{
#' \rho_c \;=\; r \times C_b, \qquad
#' r \;=\; \frac{\mathrm{cov}(X,Y)}{\sigma_X \sigma_Y}, \quad
#' C_b \;=\; \frac{2 \sigma_X \sigma_Y}
#'                {\sigma_X^2 + \sigma_Y^2 + (\mu_X - \mu_Y)^2}.
#' }
#' Hence \eqn{|\rho_c| \le |r| \le 1}, \eqn{\rho_c = r} iff
#' \eqn{\mu_X=\mu_Y} and \eqn{\sigma_X=\sigma_Y}, and \eqn{\rho_c=1} iff, in
#' addition, \eqn{r=1}. CCC is symmetric in \eqn{(X,Y)} and penalises both
#' location and scale differences; unlike Pearson's \eqn{r}, it is not invariant
#' to affine transformations that change means or variances.
#'
#' When \code{ci = TRUE}, large-sample
#' confidence intervals for \eqn{\rho_c} are returned for each pair (delta-method
#' approximation). For speed, CIs are omitted when \code{ci = FALSE}.
#'
#'If either variable has zero variance, \eqn{\rho_c} is
#' undefined and \code{NA} is returned for that pair (including the diagonal).
#'
#' Missing values are not allowed; inputs must be numeric with at least two
#' distinct non-missing values per column.
#'
#' @param data A numeric matrix or data frame with at least two numeric columns.
#' Non-numeric columns will be ignored.
#' @param ci Logical; if TRUE, return lower and upper confidence bounds
#' @param conf_level Confidence level for CI, default = 0.95
#' @param n_threads Integer \eqn{\geq 1}. Number of OpenMP threads. Defaults to
#'   \code{getOption("matrixCorr.threads", 1L)}.
#' @param output Output representation for the computed estimates.
#'   \itemize{
#'   \item \code{"matrix"} (default): full dense matrix; best when you need
#'   matrix algebra, dense heatmaps, or full compatibility with existing code.
#'   \item \code{"sparse"}: sparse matrix from \pkg{Matrix} containing only
#'   retained entries; best when many values are dropped by thresholding.
#'   \item \code{"edge_list"}: long-form data frame with columns
#'   \code{row}, \code{col}, \code{value}; convenient for filtering, joins,
#'   and network-style workflows.
#'   }
#' @param threshold Non-negative absolute-value filter for non-matrix outputs:
#'   keep entries with \code{abs(value) >= threshold}. Use
#'   \code{threshold > 0} when you want only stronger associations (typically
#'   with \code{output = "sparse"} or \code{"edge_list"}). Keep
#'   \code{threshold = 0} to retain all values. Must be \code{0} when
#'   \code{output = "matrix"}.
#' @param diag Logical; whether to include diagonal entries in
#'   \code{"sparse"} and \code{"edge_list"} outputs.
#' @param verbose Logical; if TRUE, prints how many threads are used
#'
#' @return A symmetric numeric matrix with class \code{"ccc"} and attributes:
#' \itemize{
#'   \item \code{method}: The method used ("Lin's concordance")
#'   \item \code{description}: Description string
#' }
#'  If \code{ci = FALSE}, returns matrix of class \code{"ccc"}.
#'         If \code{ci = TRUE}, returns a list with elements: \code{est},
#'         \code{lwr.ci}, \code{upr.ci}.
#'
#' @seealso \code{\link{print.ccc}}, \code{\link{plot.ccc}},
#' \code{\link{ba}}
#'
#' @seealso For repeated measurements look at \code{\link{ccc_rm_reml}},
#' \code{\link{ccc_rm_ustat}} or \code{\link{ba_rm}}
#'
#' @examples
#' # Example with multivariate normal data
#' Sigma <- matrix(c(1, 0.5, 0.3,
#'                   0.5, 1, 0.4,
#'                   0.3, 0.4, 1), nrow = 3)
#' mu <- c(0, 0, 0)
#' set.seed(123)
#' mat_mvn <- MASS::mvrnorm(n = 100, mu = mu, Sigma = Sigma)
#' result_mvn <- ccc(mat_mvn)
#' print(result_mvn)
#' summary(result_mvn)
#' plot(result_mvn)
#'
#' # Interactive viewing (requires shiny)
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   view_corr_shiny(result_mvn)
#' }
#'
#' @importFrom stats var cov cor
#' @importFrom graphics plot
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2
#' @importFrom ggplot2 theme_minimal element_text coord_fixed labs theme
#' @author Thiago de Paula Oliveira
#' @references
#' Lin L (1989). A concordance correlation coefficient to evaluate
#' reproducibility. Biometrics 45: 255-268.
#' @references
#' Lin L (2000). A note on the concordance correlation coefficient.
#' Biometrics 56: 324-325.
#' @references
#' Bland J, Altman D (1986). Statistical methods for assessing agreement
#' between two methods of clinical measurement. The Lancet 327: 307-310.
#' @export
ccc <- function(data, ci = FALSE, conf_level = 0.95,
                n_threads = getOption("matrixCorr.threads", 1L),
                output = c("matrix", "sparse", "edge_list"),
                threshold = 0,
                diag = TRUE,
                verbose = FALSE) {
  output_cfg <- .mc_validate_thresholded_output_request(
    output = output,
    threshold = threshold,
    diag = diag
  )
  check_bool(ci, arg = "ci")
  if (isTRUE(ci)) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  }
  check_bool(verbose, arg = "verbose")

  numeric_data <- validate_corr_input(data)
  mat <- numeric_data
  colnames_data <- colnames(numeric_data)
  dn <- .mc_square_dimnames(colnames_data)
  diag_payload <- list(
    n_complete = matrix(
      as.integer(nrow(mat)),
      nrow = ncol(mat),
      ncol = ncol(mat),
      dimnames = dn
    )
  )

  prev_threads <- .mc_prepare_omp_threads(
    n_threads,
    n_threads_missing = missing(n_threads)
  )
  if (!is.null(prev_threads)) {
    on.exit(.mc_exit_omp_threads(prev_threads), add = TRUE)
  }

  if (verbose) cat("Using", openmp_threads(), "OpenMP threads\n")

  if (!isTRUE(ci) && .mc_supports_direct_threshold_path(
    method = "ccc",
    na_method = "error",
    ci = FALSE,
    output = output_cfg$output,
    threshold = output_cfg$threshold,
    pairwise = FALSE,
    has_ci = FALSE
  )) {
    trip <- ccc_threshold_triplets_cpp(
      mat,
      threshold = output_cfg$threshold,
      diag = output_cfg$diag
    )
    return(.mc_finalize_triplets_output(
      triplets = trip,
      output = output_cfg$output,
      estimator_class = "ccc",
      method = "Lin's concordance",
      description = "Pairwise Lin's concordance correlation matrix",
      threshold = output_cfg$threshold,
      diag = output_cfg$diag,
      source_dim = as.integer(c(ncol(mat), ncol(mat))),
      source_dimnames = dn,
      diagnostics = diag_payload,
      symmetric = TRUE
    ))
  }

  if (ci) {
    if (nrow(mat) <= 2L) {
      abort_bad_arg("data",
        message = "must provide at least three observations per variable when {.arg ci} = TRUE.",
        .hint   = "Increase the sample size or set ci = FALSE to obtain point estimates only."
      )
    }
    ccc_lin <- ccc_with_ci_cpp(mat, conf_level)
    diag(ccc_lin$est)    <- 1
    diag(ccc_lin$lwr.ci) <- 1
    diag(ccc_lin$upr.ci) <- 1
    ccc_lin$est <- .mc_set_matrix_dimnames(ccc_lin$est, colnames_data)
    ccc_lin$lwr.ci <- .mc_set_matrix_dimnames(ccc_lin$lwr.ci, colnames_data)
    ccc_lin$upr.ci <- .mc_set_matrix_dimnames(ccc_lin$upr.ci, colnames_data)

    ccc_lin <- structure(
      ccc_lin,
      class = c("ccc", "ccc_ci"),
      method = "Lin's concordance",
      description = "Pairwise Lin's concordance with confidence intervals",
      package = "matrixCorr",
      conf.level = conf_level,
      diagnostics = diag_payload
    )
    if (!identical(output_cfg$output, "matrix")) {
      ci_attr <- list(
        est = .mc_set_matrix_dimnames(unclass(ccc_lin$est), colnames_data),
        lwr.ci = .mc_set_matrix_dimnames(unclass(ccc_lin$lwr.ci), colnames_data),
        upr.ci = .mc_set_matrix_dimnames(unclass(ccc_lin$upr.ci), colnames_data),
        conf.level = conf_level
      )
      dense_obj <- .mc_structure_corr_matrix(
        ccc_lin$est,
        class_name = "ccc",
        method = "Lin's concordance",
        description = "Pairwise Lin's concordance with confidence intervals",
        diagnostics = diag_payload,
        dimnames = dn,
        classes = c("ccc", "matrix"),
        extra_attrs = list(ci = ci_attr, conf.level = conf_level)
      )
      return(.mc_finalize_corr_output(
        dense_obj,
        output = output_cfg$output,
        threshold = output_cfg$threshold,
        diag = output_cfg$diag
      ))
    }
  } else {
    est <- ccc_cpp(mat)
    ccc_lin <- .mc_structure_corr_matrix(
      est,
      class_name = "ccc",
      method = "Lin's concordance",
      description = "Pairwise Lin's concordance correlation matrix",
      diagnostics = diag_payload,
      dimnames = dn,
      classes = c("ccc", "matrix")
    )
    ccc_lin <- .mc_finalize_corr_output(
      ccc_lin,
      output = output_cfg$output,
      threshold = output_cfg$threshold,
      diag = output_cfg$diag
    )
  }

  ccc_lin
}


#' @rdname ccc
#' @method print ccc
#' @param digits Integer; decimals for CCC estimates (default 4).
#' @param ci_digits Integer; decimals for CI bounds (default 4).
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param ... Passed to \code{\link[base]{print.data.frame}}.
#' @export
print.ccc <- function(x,
                      digits = 4,
                      ci_digits = 4,
                      n = NULL,
                      topn = NULL,
                      max_vars = NULL,
                      width = NULL,
                      show_ci = NULL,
                      ...) {
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("print_show_ci", "yes")
  )

  # -- identify object type ----------------------------------------------------
  is_ci_obj <- inherits(x, "ccc_ci") ||
    (is.list(x) && all(c("est", "lwr.ci", "upr.ci") %in% names(x)))

  if (is_ci_obj) {
    est <- as.matrix(x$est)
    lwr <- as.matrix(x$lwr.ci)
    upr <- as.matrix(x$upr.ci)
  } else if (is.matrix(x)) {
    est <- as.matrix(x)
    lwr <- matrix(NA_real_, nrow(est), ncol(est), dimnames = dimnames(est))
    upr <- lwr
  } else {
    abort_bad_arg("x",
      message = "must be a matrix or a list with elements `est`, `lwr.ci`, and `upr.ci`."
    )
  }

  rn <- rownames(est); cn <- colnames(est)
  if (is.null(rn)) rn <- paste0("m", seq_len(nrow(est)))
  if (is.null(cn)) cn <- rn

  .mc_print_corr_matrix(
    x,
    header = "Lin's concordance correlation matrix",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    mat = est,
    ...
  )
  invisible(x)
}

#' @rdname ccc
#' @method summary ccc
#' @param object A \code{"ccc"} or \code{"ccc_ci"} object to summarize.
#' @param digits Integer; decimals for CCC estimates (default 4).
#' @param ci_digits Integer; decimals for CI bounds (default 2).
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param ... Ignored.
#' @return For \code{summary.ccc}, a data frame with columns
#'   \code{item1}, \code{item2}, \code{estimate}, and (optionally)
#'   \code{lwr}, \code{upr}, plus \code{n_complete} when available.
#' @export
summary.ccc <- function(object,
                        digits = 4,
                        ci_digits = 2,
                        n = NULL,
                        topn = NULL,
                        max_vars = NULL,
                        width = NULL,
                        show_ci = NULL,
                        ...) {
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )

  # detect CI container
  is_ci_obj <- inherits(object, "ccc_ci") ||
    (is.list(object) && all(c("est", "lwr.ci", "upr.ci") %in% names(object)))

  if (is_ci_obj) {
    est <- as.matrix(object$est)
    lwr <- as.matrix(object$lwr.ci)
    upr <- as.matrix(object$upr.ci)
    conf_level <- suppressWarnings(as.numeric(attr(object, "conf.level")))
  } else if (is.matrix(object)) {
    est <- as.matrix(object)
    lwr <- matrix(NA_real_, nrow(est), ncol(est), dimnames = dimnames(est))
    upr <- lwr
    conf_level <- NA_real_
  } else {
    abort_bad_arg("object",
      message = "must be a matrix or a list with elements `est`, `lwr.ci`, and `upr.ci`."
    )
  }

  # labels (fallback if missing)
  rn <- rownames(est); cn <- colnames(est)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(est)))
  if (is.null(cn)) cn <- as.character(seq_len(ncol(est)))

  # decide whether to include CI columns
  has_any_ci <- any(is.finite(lwr) | is.finite(upr))
  include_ci <- identical(show_ci, "yes") && has_any_ci
  diag_attr <- attr(object, "diagnostics", exact = TRUE)

  # 1x1 case
  if (nrow(est) == 1L && ncol(est) == 1L) {
    df <- data.frame(
      method1  = rn[1],
      method2  = cn[1],
      estimate = round(est[1, 1], digits),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (is.list(diag_attr) && is.matrix(diag_attr$n_complete)) {
      df$n_complete <- as.integer(diag_attr$n_complete[1, 1])
    }
    if (include_ci) {
      df$lwr <- if (is.na(lwr[1,1])) NA_real_ else round(lwr[1,1], ci_digits)
      df$upr <- if (is.na(upr[1,1])) NA_real_ else round(upr[1,1], ci_digits)
    }
  } else {
    # long table over i<j
    rows <- vector("list", nrow(est) * (ncol(est) - 1L) / 2L); k <- 0L
    for (i in seq_len(nrow(est) - 1L)) {
      for (j in (i + 1L):ncol(est)) {
        k <- k + 1L
        rec <- list(
          method1  = rn[i],
          method2  = cn[j],
          estimate = round(est[i, j], digits)
        )
        if (is.list(diag_attr) && is.matrix(diag_attr$n_complete)) {
          rec$n_complete <- as.integer(diag_attr$n_complete[i, j])
        }
        if (include_ci) {
          rec$lwr <- if (is.na(lwr[i, j])) NA_real_ else round(lwr[i, j], ci_digits)
          rec$upr <- if (is.na(upr[i, j])) NA_real_ else round(upr[i, j], ci_digits)
        }
        rows[[k]] <- rec
      }
    }
    df <- do.call(rbind.data.frame, rows)
    rownames(df) <- NULL
    # ensure proper column types
    num_cols <- c("estimate", if (include_ci) c("lwr","upr"))
    for (nm in num_cols) df[[nm]] <- as.numeric(df[[nm]])
  }

  # carry attrs for printing
  df <- .mc_finalize_summary_df(df, class_name = "summary.ccc")
  attr(df, "overview") <- .mc_summary_corr_matrix(est, topn = topn)
  attr(df, "conf.level") <- if (is.finite(conf_level)) conf_level else NA_real_
  attr(df, "has_ci")     <- isTRUE(include_ci)
  attr(df, "digits")     <- digits
  attr(df, "ci_digits")  <- ci_digits
  df
}

#' @rdname ccc
#' @method print summary.ccc
#' @param ... Passed to \code{\link[base]{print.data.frame}}.
#' @export
print.summary.ccc <- function(x, digits = NULL, n = NULL,
                              topn = NULL, max_vars = NULL,
                              width = NULL, show_ci = NULL, ...) {
  .mc_print_pairwise_summary_digest(
    x,
    title = "Lin's concordance summary",
    digits = .mc_coalesce(digits, 4),
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ci_method = "delta_method",
    ...
  )
  invisible(x)
}


#' @rdname ccc
#' @method plot ccc
#' @param x An object of class \code{"ccc"} (either a matrix or a list with CIs).
#' @param title Title for the plot.
#' @param low_color Color for low CCC values.
#' @param high_color Color for high CCC values.
#' @param mid_color Color for mid CCC values.
#' @param value_text_size Text size for CCC values in the heatmap.
#' @param ci_text_size Text size for confidence intervals.
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
#' @param ... Passed to \code{ggplot2::theme()}.
#' @export
plot.ccc <- function(x,
                     title = "Lin's Concordance Correlation Heatmap",
                     low_color = "indianred1",
                     high_color = "steelblue1",
                     mid_color = "white",
                     value_text_size = 4,
                     ci_text_size = 3,
                     show_value = TRUE,
                     ...) {

  check_inherits(x, "ccc")
  check_bool(show_value, arg = "show_value")

  # --- Build long data with proper alignment by (Var1, Var2) ---
  est_mat <- if (is.list(x) && !is.null(x$est)) x$est else unclass(x)
  df_est  <- as.data.frame(as.table(est_mat))
  names(df_est) <- c("Var1", "Var2", "CCC")

  if (is.list(x) && !is.null(x$lwr.ci) && !is.null(x$upr.ci)) {
    df_lwr <- as.data.frame(as.table(x$lwr.ci)); names(df_lwr)[3] <- "lwr"
    df_upr <- as.data.frame(as.table(x$upr.ci)); names(df_upr)[3] <- "upr"
    df <- Reduce(function(a, b) merge(a, b, by = c("Var1","Var2"), all = TRUE),
                 list(df_est, df_lwr, df_upr))

    # Blank CI on the diagonal (show 1.00 but no CI there)
    diag_idx <- df$Var1 == df$Var2
    df$lwr[diag_idx] <- NA_real_
    df$upr[diag_idx] <- NA_real_
    df$ci_label <- ifelse(is.na(df$lwr) | is.na(df$upr),
                          NA_character_,
                          sprintf("(%.2f, %.2f)", df$lwr, df$upr))
  } else {
    df <- df_est
    df$ci_label <- NA_character_
  }

  # Reverse Y axis for heatmap look, keep X in natural order
  lev_row <- unique(df_est$Var1)
  lev_col <- unique(df_est$Var2)
  df$Var1 <- factor(df$Var1, levels = rev(lev_row))
  df$Var2 <- factor(df$Var2, levels = lev_col)

  # Main estimate label
  df$label <- sprintf("%.2f", df$CCC)

  # --- Plot ---
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Var2, y = Var1, fill = CCC)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(
      low = low_color, high = high_color, mid = mid_color,
      midpoint = 0, limits = c(-1, 1), name = "CCC"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      ...
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL)

  if (isTRUE(show_value)) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = label), size = value_text_size)
  }

  if (isTRUE(show_value) && any(!is.na(df$ci_label))) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = ci_label, y = as.numeric(Var1) - 0.25),
      size = ci_text_size, color = "gray30", na.rm = TRUE
    )
  }

  p
}

