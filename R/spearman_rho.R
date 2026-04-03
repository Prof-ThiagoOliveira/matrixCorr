#' @title Pairwise Spearman's rank correlation
#'
#' @description
#' Computes pairwise Spearman's rank correlations for the numeric columns of a
#' matrix or data frame using a high-performance 'C++' backend. Optional
#' confidence intervals are available via a jackknife Euclidean-likelihood
#' method.
#'
#' @param data A numeric matrix or a data frame with at least two numeric
#' columns. All non-numeric columns will be excluded. Each column must have
#' at least two non-missing values.
#' @param check_na Logical (default \code{TRUE}). If \code{TRUE}, the input is
#' required to be free of \code{NA}/\code{NaN}/\code{Inf}. Set to
#' \code{FALSE} only when the caller already handled missingness.
#' @param ci Logical (default \code{FALSE}). If \code{TRUE}, attach
#' jackknife Euclidean-likelihood confidence intervals for the off-diagonal
#' Spearman correlations.
#' @param conf_level Confidence level used when \code{ci = TRUE}. Default is
#' \code{0.95}.
#'
#' @return A symmetric numeric matrix where the \code{(i, j)}-th element is
#' the Spearman correlation between the \code{i}-th and \code{j}-th
#' numeric columns of the input. When \code{ci = TRUE}, the object also
#' carries a \code{ci} attribute with elements \code{est}, \code{lwr.ci},
#' \code{upr.ci}, and \code{conf.level}. When pairwise-complete evaluation is
#' used, pairwise sample sizes are stored in \code{attr(x, "diagnostics")$n_complete}.
#'
#' @details
#' For each column \eqn{j=1,\ldots,p}, let
#' \eqn{R_{\cdot j} \in \{1,\ldots,n\}^n} denote the (mid-)ranks of
#' \eqn{X_{\cdot j}}, assigning average ranks to ties. The mean rank is
#' \eqn{\bar R_j = (n+1)/2} regardless of ties. Define the centred rank
#' vectors \eqn{\tilde R_{\cdot j} = R_{\cdot j} - \bar R_j \mathbf{1}},
#' where \eqn{\mathbf{1}\in\mathbb{R}^n} is the all-ones vector. The
#' Spearman correlation between columns \eqn{i} and \eqn{j} is the Pearson
#' correlation of their rank vectors:
#' \deqn{
#' \rho_S(i,j) \;=\;
#' \frac{\sum_{k=1}^n (R_{ki}-\bar R_i)(R_{kj}-\bar R_j)}
#'      {\sqrt{\sum_{k=1}^n (R_{ki}-\bar R_i)^2}\;
#'       \sqrt{\sum_{k=1}^n (R_{kj}-\bar R_j)^2}}.
#' }
#' In matrix form, with \eqn{R=[R_{\cdot 1},\ldots,R_{\cdot p}]},
#' \eqn{\mu=(n+1)\mathbf{1}_p/2} for \eqn{\mathbf{1}_p\in\mathbb{R}^p}, and
#' \eqn{S_R=\bigl(R-\mathbf{1}\mu^\top\bigr)^\top
#'            \bigl(R-\mathbf{1}\mu^\top\bigr)/(n-1)},
#' the Spearman correlation matrix is
#' \deqn{
#' \widehat{\rho}_S \;=\; D^{-1/2} S_R D^{-1/2}, \qquad
#' D \;=\; \mathrm{diag}(\mathrm{diag}(S_R)).
#' }
#' When there are no ties, the familiar rank-difference formula obtains
#' \deqn{
#' \rho_S(i,j) \;=\; 1 - \frac{6}{n(n^2-1)} \sum_{k=1}^n d_k^2,
#' \quad d_k \;=\; R_{ki}-R_{kj},
#' }
#' but this expression does \emph{not} hold under ties; computing Pearson on
#' mid-ranks (as above) is the standard tie-robust approach. Without ties,
#' \eqn{\mathrm{Var}(R_{\cdot j})=(n^2-1)/12}; with ties, the variance is
#' smaller.
#'
#' \eqn{\rho_S(i,j) \in [-1,1]} and \eqn{\widehat{\rho}_S} is symmetric
#' positive semi-definite by construction (up to floating-point error). The
#' implementation symmetrises the result to remove round-off asymmetry.
#' Spearman's correlation is invariant to strictly monotone transformations
#' applied separately to each variable.
#'
#' \strong{Computation.} Each column is ranked (mid-ranks) to form \eqn{R}.
#' The product \eqn{R^\top R} is computed via a 'BLAS' symmetric rank update
#' ('SYRK'), and centred using
#' \deqn{
#' (R-\mathbf{1}\mu^\top)^\top (R-\mathbf{1}\mu^\top)
#' \;=\; R^\top R \;-\; n\,\mu\mu^\top,
#' }
#' avoiding an explicit centred copy. Division by \eqn{n-1} yields the sample
#' covariance of ranks; standardising by \eqn{D^{-1/2}} gives \eqn{\widehat{\rho}_S}.
#' Columns with zero rank variance (all values equal) are returned as \code{NA}
#' along their row/column; the corresponding diagonal entry is also \code{NA}.
#'
#' When \code{check_na = FALSE}, each \eqn{(i,j)} estimate is recomputed on the
#' pairwise complete-case overlap of columns \eqn{i} and \eqn{j}. When
#' \code{ci = TRUE}, confidence intervals are computed in 'C++' using the
#' jackknife Euclidean-likelihood method of de Carvalho and Marques (2012).
#' For a pairwise estimate \eqn{U = \hat\rho_S}, delete-one jackknife
#' pseudo-values are formed as
#' \deqn{
#' Z_i = nU - (n-1)U_{(-i)}, \qquad i = 1,\ldots,n,
#' }
#' where \eqn{U_{(-i)}} is the Spearman correlation after removing observation
#' \eqn{i}. The confidence limits solve
#' \deqn{
#' \frac{n(U-\theta)^2}{n^{-1}\sum_{i=1}^n (Z_i - \theta)^2}
#' = \chi^2_{1,\;\texttt{conf\_level}}.
#' }
#'
#' Ranking costs
#' \eqn{O\!\bigl(p\,n\log n\bigr)}; forming and normalising
#' \eqn{R^\top R} costs \eqn{O\!\bigl(n p^2\bigr)} with \eqn{O(p^2)} additional
#' memory. The optional jackknife Euclidean-likelihood confidence intervals add
#' per-pair delete-one recomputation work and are intended for inference rather
#' than raw-matrix throughput.
#'
#' @note Missing values are not allowed when \code{check_na = TRUE}. Columns
#' with fewer than two observations are excluded.
#'
#' @references
#' Spearman, C. (1904). The proof and measurement of association between
#' two things. International Journal of Epidemiology, 39(5), 1137-1150.
#'
#' de Carvalho, M., & Marques, F. (2012). Jackknife Euclidean
#' likelihood-based inference for Spearman's rho. North American Actuarial
#' Journal, 16(4), 487-492.
#'
#' @examples
#' ## Monotone transformation invariance (Spearman is rank-based)
#' set.seed(123)
#' n <- 400; p <- 6; rho <- 0.6
#' Sigma <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
#' L <- chol(Sigma)
#' X <- matrix(rnorm(n * p), n, p) %*% L
#' colnames(X) <- paste0("V", seq_len(p))
#'
#' X_mono <- X
#' X_mono[, 1] <- exp(X_mono[, 1])
#' X_mono[, 2] <- log1p(exp(X_mono[, 2]))
#' X_mono[, 3] <- X_mono[, 3]^3
#'
#' sp_X <- spearman_rho(X)
#' sp_m <- spearman_rho(X_mono)
#' summary(sp_X)
#' round(max(abs(sp_X - sp_m)), 3)
#' plot(sp_X)
#'
#' ## Confidence intervals
#' sp_ci <- spearman_rho(X[, 1:3], ci = TRUE)
#' print(sp_ci, show_ci = "yes")
#' summary(sp_ci)
#'
#' ## Ties handled via mid-ranks
#' tied <- cbind(
#'   a = rep(1:5, each = 20),
#'   b = rep(5:1, each = 20) + rnorm(100, sd = 0.1),
#'   c = as.numeric(gl(10, 10))
#' )
#' sp_tied <- spearman_rho(tied, ci = TRUE)
#' print(sp_tied, digits = 2, show_ci = "yes")
#'
#' # Interactive viewing (requires shiny)
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   view_corr_shiny(sp_X)
#' }
#' @useDynLib matrixCorr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @seealso \code{\link{print.spearman_rho}}, \code{\link{plot.spearman_rho}}
#' @author Thiago de Paula Oliveira
#' @export
spearman_rho <- function(data, check_na = TRUE, ci = FALSE, conf_level = 0.95) {
  check_bool(ci, arg = "ci")
  if (isTRUE(ci)) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  }

  numeric_data <- validate_corr_input(data, check_na = check_na)
  colnames_data <- colnames(numeric_data)
  diagnostics <- NULL
  ci_attr <- NULL

  if (isTRUE(check_na) && !isTRUE(ci)) {
    result <- spearman_matrix_cpp(numeric_data)
  } else {
    pairwise <- spearman_matrix_pairwise_cpp(
      numeric_data,
      return_ci = ci,
      conf_level = conf_level
    )
    result <- pairwise$est
    diagnostics <- list(n_complete = pairwise$n_complete)
    dimnames(diagnostics$n_complete) <- list(colnames_data, colnames_data)
    if (isTRUE(ci)) {
      ci_attr <- list(
        est = unclass(result),
        lwr.ci = unclass(pairwise$lwr),
        upr.ci = unclass(pairwise$upr),
        conf.level = pairwise$conf_level
      )
      dimnames(ci_attr$est) <- list(colnames_data, colnames_data)
      dimnames(ci_attr$lwr.ci) <- list(colnames_data, colnames_data)
      dimnames(ci_attr$upr.ci) <- list(colnames_data, colnames_data)
    }
  }

  colnames(result) <- rownames(result) <- colnames_data
  out <- .mc_structure_corr_matrix(
    result,
    class_name = "spearman_rho",
    method = "spearman",
    description = "Pairwise Spearman's rank correlation matrix",
    diagnostics = diagnostics
  )
  if (!is.null(ci_attr)) {
    attr(out, "ci") <- ci_attr
    attr(out, "conf.level") <- conf_level
  }
  out
}

.mc_spearman_ci_attr <- function(x) {
  attr(x, "ci", exact = TRUE)
}

.mc_spearman_pairwise_summary <- function(object,
                                          digits = 4,
                                          ci_digits = 3,
                                          show_ci = NULL) {
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )
  check_inherits(object, "spearman_rho")

  est <- as.matrix(object)
  rn <- rownames(est); cn <- colnames(est)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(est)))
  if (is.null(cn)) cn <- as.character(seq_len(ncol(est)))

  ci <- .mc_spearman_ci_attr(object)
  diag_attr <- attr(object, "diagnostics", exact = TRUE)
  include_ci <- identical(show_ci, "yes") && !is.null(ci)

  rows <- vector("list", nrow(est) * (ncol(est) - 1L) / 2L)
  k <- 0L
  for (i in seq_len(nrow(est) - 1L)) {
    for (j in (i + 1L):ncol(est)) {
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
        rec$lwr <- if (!is.null(ci$lwr.ci) && is.finite(ci$lwr.ci[i, j])) round(ci$lwr.ci[i, j], ci_digits) else NA_real_
        rec$upr <- if (!is.null(ci$upr.ci) && is.finite(ci$upr.ci[i, j])) round(ci$upr.ci[i, j], ci_digits) else NA_real_
      }
      rows[[k]] <- rec
    }
  }

  df <- do.call(rbind.data.frame, rows)
  rownames(df) <- NULL
  if ("estimate" %in% names(df)) df$estimate <- as.numeric(df$estimate)
  if ("lwr" %in% names(df)) df$lwr <- as.numeric(df$lwr)
  if ("upr" %in% names(df)) df$upr <- as.numeric(df$upr)
  if ("n_complete" %in% names(df)) df$n_complete <- as.integer(df$n_complete)

  out <- structure(df, class = c("summary.spearman_rho", "data.frame"))
  attr(out, "overview") <- .mc_summary_corr_matrix(object)
  attr(out, "has_ci") <- include_ci
  attr(out, "conf.level") <- if (is.null(ci)) NA_real_ else ci$conf.level
  attr(out, "digits") <- digits
  attr(out, "ci_digits") <- ci_digits
  out
}

#' @rdname spearman_rho
#' @method print spearman_rho
#' @title Print Method for \code{spearman_rho} Objects
#'
#' @param x An object of class \code{spearman_rho}.
#' @param digits Integer; number of decimal places to print.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param ci_digits Integer; digits for Spearman confidence limits.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @return Invisibly returns the \code{spearman_rho} object.
#' @export
print.spearman_rho <- function(x,
                               digits = 4,
                               n = NULL,
                               topn = NULL,
                               max_vars = NULL,
                               width = NULL,
                               ci_digits = 3,
                               show_ci = NULL,
                               ...) {
  .mc_print_corr_matrix(
    x,
    header = "Spearman correlation matrix",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
  invisible(x)
}

#' @rdname spearman_rho
#' @method plot spearman_rho
#' @title Plot Method for \code{spearman_rho} Objects
#'
#' @param x An object of class \code{spearman_rho}.
#' @param title Plot title. Default is \code{"Spearman's rank correlation
#' heatmap"}.
#' @param low_color Color for the minimum rho value. Default is
#'  \code{"indianred1"}.
#' @param high_color Color for the maximum rho value. Default is
#' \code{"steelblue1"}.
#' @param mid_color Color for zero correlation. Default is \code{"white"}.
#' @param value_text_size Font size for displaying correlation values. Default
#' is \code{4}.
#' @param ci_text_size Text size for confidence intervals in the heatmap.
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
#' @param ... Additional arguments passed to \code{ggplot2::theme()} or other
#' \code{ggplot2} layers.
#'
#' @return A \code{ggplot} object representing the heatmap.
#' @import ggplot2
#' @export
plot.spearman_rho <-
  function(x, title = "Spearman's rank correlation heatmap",
           low_color = "indianred1", high_color = "steelblue1",
           mid_color = "white", value_text_size = 4,
           ci_text_size = 3, show_value = TRUE, ...) {
    check_bool(show_value, arg = "show_value")
    ci <- .mc_spearman_ci_attr(x)
    if (is.null(ci) || is.null(ci$lwr.ci) || is.null(ci$upr.ci)) {
      return(.mc_plot_corr_matrix(
        x, class_name = "spearman_rho", fill_name = "Rho",
        title = title, low_color = low_color, high_color = high_color,
        mid_color = mid_color, value_text_size = value_text_size,
        show_value = show_value, ...
      ))
    }

    est_mat <- as.matrix(x)
    df_est <- as.data.frame(as.table(est_mat))
    names(df_est) <- c("Var1", "Var2", "rho")

    df_lwr <- as.data.frame(as.table(ci$lwr.ci))
    names(df_lwr)[3] <- "lwr"
    df_upr <- as.data.frame(as.table(ci$upr.ci))
    names(df_upr)[3] <- "upr"
    df <- Reduce(
      function(a, b) merge(a, b, by = c("Var1", "Var2"), all = TRUE),
      list(df_est, df_lwr, df_upr)
    )

    diag_idx <- df$Var1 == df$Var2
    df$lwr[diag_idx] <- NA_real_
    df$upr[diag_idx] <- NA_real_
    df$ci_label <- ifelse(
      is.na(df$lwr) | is.na(df$upr),
      NA_character_,
      sprintf("[%.3f, %.3f]", df$lwr, df$upr)
    )

    lev_row <- unique(df_est$Var1)
    lev_col <- unique(df_est$Var2)
    df$Var1 <- factor(df$Var1, levels = rev(lev_row))
    df$Var2 <- factor(df$Var2, levels = lev_col)
    df$label <- sprintf("%.2f", df$rho)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = Var2, y = Var1, fill = .data$rho)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient2(
        low = low_color,
        high = high_color,
        mid = mid_color,
        midpoint = 0,
        limits = c(-1, 1),
        name = "Rho"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank(),
        ...
      ) +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = title, x = NULL, y = NULL)

    if (isTRUE(show_value)) {
      p <- p + ggplot2::geom_text(ggplot2::aes(label = label), size = value_text_size, color = "black")
    }

    if (isTRUE(show_value) && any(!is.na(df$ci_label))) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = ci_label, y = as.numeric(Var1) - 0.25),
        size = ci_text_size,
        color = "gray30",
        na.rm = TRUE
      )
    }

    p
  }

#' @rdname spearman_rho
#' @method summary spearman_rho
#' @param object An object of class \code{spearman_rho}.
#' @param ci_digits Integer; digits for Spearman confidence limits in the
#'   pairwise summary.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @export
summary.spearman_rho <- function(object,
                                 n = NULL,
                                 topn = NULL,
                                 max_vars = NULL,
                                 width = NULL,
                                 ci_digits = 3,
                                 show_ci = NULL,
                                 ...) {
  check_inherits(object, "spearman_rho")
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )
  if (is.null(.mc_spearman_ci_attr(object))) {
    return(.mc_summary_corr_matrix(object, topn = topn))
  }
  .mc_spearman_pairwise_summary(
    object,
    ci_digits = ci_digits,
    show_ci = show_ci
  )
}

#' @rdname spearman_rho
#' @method print summary.spearman_rho
#' @param x An object of class \code{summary.spearman_rho}.
#' @export
print.summary.spearman_rho <- function(x, digits = NULL, n = NULL,
                                       topn = NULL, max_vars = NULL,
                                       width = NULL, show_ci = NULL, ...) {
  .mc_print_pairwise_summary_digest(
    x,
    title = "Spearman correlation summary",
    digits = .mc_coalesce(digits, 4),
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ci_method = "jackknife_euclidean_likelihood",
    ...
  )
  invisible(x)
}
