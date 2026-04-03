#' @title Pairwise skipped correlation
#'
#' @description
#' Computes all pairwise skipped correlation coefficients for the numeric
#' columns of a matrix or data frame using a high-performance 'C++' backend.
#'
#' Skipped correlation detects bivariate outliers using a projection rule and
#' then computes Pearson or Spearman correlation on the retained observations.
#' It is designed for situations where marginally robust methods can still be
#' distorted by unusual points in the joint data cloud.
#'
#' @param data A numeric matrix or a data frame with at least two numeric
#' columns. All non-numeric columns will be excluded.
#' @param method Correlation computed after removing projected outliers. One of
#'   \code{"pearson"} (default) or \code{"spearman"}.
#' @param na_method One of \code{"error"} (default) or \code{"pairwise"}.
#'   With \code{"error"}, the function requires all retained numeric columns to
#'   be free of missing or non-finite values and aborts otherwise. This is the
#'   recommended setting when you want a single common sample size across all
#'   pairs, reproducible skipped-row diagnostics on the same rows, or bootstrap
#'   inference via \code{ci = TRUE} / \code{p_value = TRUE}. With
#'   \code{"pairwise"}, each variable pair is computed on its own overlap of
#'   finite rows. This is more permissive for incomplete data, but different
#'   pairs can be based on different effective samples and different skipped-row
#'   sets, so the resulting matrix is less directly comparable across entries.
#' @param stand Logical; if \code{TRUE} (default), each variable in the pair is
#'   centred by its median and divided by a robust scale estimate before the
#'   projection outlier search. The scale estimate is the MAD when positive,
#'   with fallback to \eqn{\mathrm{IQR}/1.34898} and then the usual sample
#'   standard deviation if needed. This standardisation affects only outlier
#'   detection, not the final correlation computed on the retained
#'   observations.
#' @param outlier_rule One of \code{"idealf"} (default) or \code{"mad"}.
#'   The default uses the ideal-fourths interquartile width of projected
#'   distances; \code{"mad"} uses the median absolute deviation of projected
#'   distances.
#' @param cutoff Positive numeric constant multiplying the projected spread in
#'   the outlier rule
#'   \eqn{\mathrm{med}(d_{i\cdot}) + cutoff \times s(d_{i\cdot})}. Larger
#'   values flag fewer observations as outliers; smaller values flag more.
#'   Default \code{sqrt(qchisq(0.975, df = 2))}.
#' @param n_threads Integer \eqn{\geq 1}. Number of OpenMP threads. Defaults to
#'   \code{getOption("matrixCorr.threads", 1L)}.
#' @param return_masks Logical; if \code{TRUE}, attach compact pairwise skipped-row
#'   indices as an attribute. Default \code{FALSE}.
#' @param ci Logical; if \code{TRUE}, attach percentile-bootstrap confidence
#'   intervals for each skipped correlation using the Wilcox (2015) B2
#'   resampling scheme. Default \code{FALSE}.
#' @param p_value Logical; if \code{TRUE}, attach bootstrap p-values for testing
#'   whether each skipped correlation is zero. Default \code{FALSE}.
#' @param conf_level Confidence level used when \code{ci = TRUE}. Default
#'   \code{0.95}.
#' @param n_boot Integer \eqn{\geq 2}. Number of bootstrap resamples used when
#'   \code{ci = TRUE} or \code{p_value = TRUE}. Default \code{2000}.
#' @param p_adjust One of \code{"none"} (default), \code{"hochberg"}, or
#'   \code{"ecp"}. Optional familywise-error procedure applied to the matrix of
#'   bootstrap p-values. \code{"hochberg"} corresponds to method H in Wilcox,
#'   Rousselet, and Pernet (2018); \code{"ecp"} corresponds to their simulated
#'   critical-p-value method ECP.
#' @param fwe_level Familywise-error level used when
#'   \code{p_adjust = "hochberg"} or \code{"ecp"}. Default \code{0.05}.
#' @param n_mc Integer \eqn{\geq 10}. Number of null Monte Carlo data sets used
#'   when \code{p_adjust = "ecp"} to estimate the critical p-value. Default
#'   \code{1000}.
#' @param seed Optional positive integer used to seed the bootstrap resampling
#'   when \code{ci = TRUE} or \code{p_value = TRUE}. If \code{NULL}, a fresh
#'   internal seed is generated.
#' @param x An object of class \code{skipped_corr}.
#' @param var1,var2 Optional column names or 1-based column indices used by
#'   [skipped_corr_masks()] to extract the skipped-row indices for one pair.
#' @param digits Integer; number of digits to print.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param ... Additional arguments passed to the underlying print or plot helper.
#' @param title Character; plot title.
#' @param low_color,high_color,mid_color Colors used in the heatmap.
#' @param value_text_size Numeric text size for overlaid cell values.
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
#'
#' @return A symmetric correlation matrix with class \code{skipped_corr} and
#'   attributes \code{method = "skipped_correlation"}, \code{description}, and
#'   \code{package = "matrixCorr"}. When \code{return_masks = TRUE}, the matrix
#'   also carries a \code{skipped_masks} attribute containing compact pairwise
#'   skipped-row indices. The \code{diagnostics} attribute stores per-pair
#'   complete-case counts and skipped-row counts/proportions. When
#'   \code{ci = TRUE} or \code{p_value = TRUE}, bootstrap inference matrices are
#'   attached via attributes.
#'
#' @details
#' Let \eqn{X \in \mathbb{R}^{n \times p}} be a numeric matrix with rows as
#' observations and columns as variables. For a given pair of columns
#' \eqn{(x, y)}, write the observed bivariate points as
#' \eqn{u_i = (x_i, y_i)^\top}, \eqn{i=1,\ldots,n}. If \code{stand = TRUE},
#' each margin is first centred by its median and divided by a robust scale
#' estimate before outlier detection; otherwise the original pair is used. The
#' robust scale is the MAD when positive, with fallback to
#' \eqn{\mathrm{IQR}/1.34898} and then the usual sample standard deviation if
#' needed. Let \eqn{\tilde u_i} denote the resulting points and let \eqn{c} be
#' the componentwise median center of the detection cloud.
#'
#' For each observation \eqn{i}, define the direction vector
#' \eqn{b_i = \tilde u_i - c}. When \eqn{\|b_i\| > 0}, all observations are
#' projected onto the line through \eqn{c} in direction \eqn{b_i}. The
#' projected distances are
#' \deqn{
#' d_{ij} \;=\; \frac{|(\tilde u_j - c)^\top b_i|}{\|b_i\|},
#' \qquad j=1,\ldots,n.
#' }
#' For each direction \eqn{i}, observation \eqn{j} is flagged as an outlier if
#' 
#' \deqn{
#' d_{ij} \;>\; \mathrm{med}(d_{i\cdot}) + g\, s(d_{i\cdot}),
#' \qquad g = \code{cutoff},
#' }
#' where \eqn{s(\cdot)} is either the ideal-fourths interquartile width
#' (\code{outlier_rule = "idealf"}) or the median absolute deviation
#' (\code{outlier_rule = "mad"}). An observation is removed if it is flagged
#' for at least one projection direction. The skipped correlation is then the
#' ordinary Pearson or Spearman correlation computed from the retained
#' observations:
#' \deqn{
#' r_{\mathrm{skip}}(x,y) \;=\;
#' \mathrm{cor}\!\left(x_{\mathcal{K}}, y_{\mathcal{K}}\right),
#' }
#' where \eqn{\mathcal{K}} is the index set of observations not flagged as
#' outliers.
#'
#' Unlike marginally robust methods such as \code{pbcor()}, \code{wincor()},
#' or \code{bicor()}, skipped correlation is explicitly pairwise because
#' outlier detection depends on the joint geometry of each variable pair. As a
#' result, the reported matrix need not be positive semidefinite, even with
#' complete data.
#'
#' \strong{Computational notes.} In the complete-data path, each column pair
#' requires a full bivariate projection search, so the dominant cost is higher
#' than for marginal robust methods. The implementation evaluates pairs in
#' 'C++'; where available, pairs are processed with 'OpenMP' parallelism. With
#' \code{na_method = "pairwise"}, each pair is recomputed on its overlap of
#' non-missing rows.
#'
#' \strong{Bootstrap inference.} When \code{ci = TRUE} or \code{p_value = TRUE},
#' the implementation uses the percentile-bootstrap strategy studied by Wilcox
#' (2015). Each bootstrap replicate resamples whole observation pairs with
#' replacement, reruns the skipped-correlation outlier detection on the
#' resampled data, and recomputes the skipped correlation on the retained
#' observations. This corresponds to Wilcox's B2 method and avoids the
#' statistically unsatisfactory shortcut of removing outliers only once before
#' bootstrapping. Bootstrap inference currently requires complete data
#' (\code{na_method = "error"}). When \code{p_adjust = "hochberg"}, the
#' bootstrap p-values are processed with Hochberg's step-up procedure (method H
#' in Wilcox, Rousselet, and Pernet, 2018). When \code{p_adjust = "ecp"}, the
#' package follows their ECP method and simulates \code{n_mc} null data sets
#' from a \eqn{p}-variate normal distribution with identity covariance,
#' recomputes the pairwise bootstrap p-values for each null data set, stores the
#' minimum p-value from each run, and estimates the \code{fwe_level} quantile of
#' that null distribution using the Harrell-Davis estimator. Hypotheses are then
#' rejected when their observed bootstrap p-values are less than or equal to the
#' estimated critical p-value. The calibrated H1 procedure from Wilcox,
#' Rousselet, and Pernet (2018) is not currently implemented.
#'
#' @references
#' Wilcox, R. R. (2004). Inferences based on a skipped correlation coefficient.
#' Journal of Applied Statistics, 31(2), 131-143.
#' \doi{10.1080/0266476032000148821}
#'
#' Wilcox, R. R. (2015). Inferences about the skipped correlation coefficient:
#' Dealing with heteroscedasticity and non-normality. Journal of Modern Applied
#' Statistical Methods, 14(1), 172-188.
#' \doi{10.22237/jmasm/1430453580}
#'
#' Wilcox, R. R., Rousselet, G. A., & Pernet, C. R. (2018). Improved methods
#' for making inferences about multiple skipped correlations. Journal of
#' Statistical Computation and Simulation, 88(16), 3116-3131.
#' \doi{10.1080/00949655.2018.1501051}
#'
#' @seealso [pbcor()], [wincor()], [bicor()]
#'
#' @examples
#' set.seed(12)
#' X <- matrix(rnorm(160 * 4), ncol = 4)
#' X[1, 1] <- 9
#' X[1, 2] <- -8
#'
#' R <- skipped_corr(X, method = "pearson")
#' print(R, digits = 2)
#' summary(R)
#' plot(R)
#'
#' Rm <- skipped_corr(X, method = "pearson", return_masks = TRUE)
#' skipped_corr_masks(Rm, 1, 2)
#'
#' # Example 1:
#' Xm <- as.matrix(datasets::mtcars[, c("mpg", "disp", "hp", "wt")])
#' Rm2 <- skipped_corr(Xm, method = "spearman")
#' print(Rm2, digits = 2)
#'
#' # Example 2:
#' Ri <- skipped_corr(Xm, method = "pearson", ci = TRUE, n_boot = 40, seed = 1)
#' Ri$ci
#'
#' # Interactive viewing (requires shiny)
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   view_corr_shiny(R)
#' }
#'
#' @author Thiago de Paula Oliveira
#' @export
skipped_corr <- function(data,
                    method = c("pearson", "spearman"),
                    na_method = c("error", "pairwise"),
                    stand = TRUE,
                    outlier_rule = c("idealf", "mad"),
                    cutoff = sqrt(stats::qchisq(0.975, df = 2)),
                    n_threads = getOption("matrixCorr.threads", 1L),
                    return_masks = FALSE,
                    ci = FALSE,
                    p_value = FALSE,
                    conf_level = 0.95,
                    n_boot = 2000L,
                    p_adjust = c("none", "hochberg", "ecp"),
                    fwe_level = 0.05,
                    n_mc = 1000L,
                    seed = NULL) {
  method <- match.arg(method)
  na_method <- match.arg(na_method)
  outlier_rule <- match.arg(outlier_rule)
  p_adjust <- match.arg(p_adjust)
  check_bool(stand, arg = "stand")
  check_bool(return_masks, arg = "return_masks")
  check_bool(ci, arg = "ci")
  check_bool(p_value, arg = "p_value")
  check_scalar_numeric(cutoff, arg = "cutoff", lower = 0, closed_lower = FALSE)
  check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  check_prob_scalar(fwe_level, arg = "fwe_level", open_ends = TRUE)
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")
  n_boot <- check_scalar_int_pos(n_boot, arg = "n_boot")
  n_mc <- check_scalar_int_pos(n_mc, arg = "n_mc")
  if (n_mc < 10L) {
    abort_bad_arg("n_mc", message = "must be >= 10.")
  }
  if (!is.null(seed)) {
    seed <- check_scalar_int_pos(seed, arg = "seed")
  }
  if ((isTRUE(ci) || isTRUE(p_value)) && na_method != "error") {
    abort_bad_arg(
      "na_method",
      message = "{.arg ci} and {.arg p_value} currently require {.code na_method = \"error\"}."
    )
  }
  if (!isTRUE(p_value) && !identical(p_adjust, "none")) {
    abort_bad_arg(
      "p_adjust",
      message = "{.arg p_adjust} can only be used when {.arg p_value} = TRUE."
    )
  }

  numeric_data <- if (na_method == "error") {
    validate_corr_input(data)
  } else {
    validate_corr_input(data, check_na = FALSE)
  }
  colnames_data <- colnames(numeric_data)

  method_int <- switch(method, pearson = 0L, spearman = 1L)
  use_mad <- identical(outlier_rule, "mad")
  res <- skipcor_matrix_cpp(
    numeric_data,
    method_int = method_int,
    stand = stand,
    use_mad = use_mad,
    gval = cutoff,
    min_n = 5L,
    n_threads = n_threads,
    return_masks = return_masks,
    return_inference = isTRUE(ci) || isTRUE(p_value),
    conf_level = conf_level,
    n_boot = n_boot,
    seed = if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else seed,
    multiple_method_int = switch(p_adjust, none = 0L, hochberg = 1L, ecp = 2L),
    fwe_level = fwe_level,
    n_mc = n_mc
  )

  mask_payload <- NULL
  if (isTRUE(return_masks)) {
    mask_payload <- structure(
      list(
        pair_i = as.integer(res$pair_i),
        pair_j = as.integer(res$pair_j),
        skipped_rows = unname(lapply(res$skipped_rows, as.integer)),
        n_rows = nrow(numeric_data),
        n_cols = ncol(numeric_data),
        colnames = colnames_data
      ),
      class = "skipped_corr_masks"
    )
  }

  diag_payload <- list(
    n_complete = unclass(res$n_complete),
    skipped_n = unclass(res$skipped_n),
    skipped_prop = unclass(res$skipped_prop)
  )
  res_cor <- res$cor
  colnames(res_cor) <- rownames(res_cor) <- colnames_data
  for (nm in names(diag_payload)) {
    dimnames(diag_payload[[nm]]) <- list(colnames_data, colnames_data)
  }

  infer_payload <- NULL
  if (isTRUE(ci) || isTRUE(p_value)) {
    infer_payload <- list(
      conf_level = conf_level,
      n_boot = n_boot,
      method = switch(
        p_adjust,
        none = "bootstrap_b2",
        hochberg = "method_h",
        ecp = "ecp"
      )
    )
    if (isTRUE(ci)) {
      lwr <- unclass(res$lwr)
      upr <- unclass(res$upr)
      dimnames(lwr) <- dimnames(upr) <- list(colnames_data, colnames_data)
      infer_payload$lwr <- lwr
      infer_payload$upr <- upr
    }
    if (isTRUE(p_value)) {
      p_mat <- unclass(res$p_value)
      dimnames(p_mat) <- list(colnames_data, colnames_data)
      infer_payload$p_value <- p_mat
      if (identical(p_adjust, "hochberg")) {
        padj <- unclass(res$p_value_adjusted)
        reject <- unclass(res$reject) > 0.5
        dimnames(padj) <- list(colnames_data, colnames_data)
        dimnames(reject) <- list(colnames_data, colnames_data)
        infer_payload$p_value_adjusted <- padj
        infer_payload$reject <- reject
        infer_payload$p_adjust <- "hochberg"
      } else if (identical(p_adjust, "ecp")) {
        reject <- unclass(res$reject) > 0.5
        dimnames(reject) <- list(colnames_data, colnames_data)
        infer_payload$reject <- reject
        infer_payload$critical_p_value <- unname(res$critical_p_value)
        infer_payload$p_adjust <- "ecp"
      } else {
        infer_payload$p_adjust <- "none"
      }
    }
  }

  out <- .mc_structure_corr_matrix(
    res_cor,
    class_name = "skipped_corr",
    method = "skipped_correlation",
    description = paste0(
      "Skipped correlation; base = ", method,
      "; rule = ", outlier_rule,
      "; standardise = ", stand,
      "; NA mode = ", na_method, "."
    ),
    diagnostics = diag_payload
  )
  if (isTRUE(return_masks)) attr(out, "skipped_masks") <- mask_payload
  if (!is.null(infer_payload)) {
    ci_attr <- NULL
    if (isTRUE(ci)) {
      ci_attr <- list(
        est = unclass(res_cor),
        lwr.ci = infer_payload$lwr,
        upr.ci = infer_payload$upr,
        conf.level = infer_payload$conf_level
      )
      attr(out, "ci") <- ci_attr
    }

    inference_attr <- list(
      method = infer_payload$method,
      n_boot = infer_payload$n_boot,
      conf.level = infer_payload$conf_level,
      fwe_level = fwe_level,
      p_adjust = infer_payload$p_adjust %||% "none"
    )
    if (isTRUE(ci)) {
      inference_attr$ci <- ci_attr
    }
    if (isTRUE(p_value)) {
      inference_attr$p_value <- infer_payload$p_value
      if (!is.null(infer_payload$p_value_adjusted)) {
        inference_attr$p_value_adjusted <- infer_payload$p_value_adjusted
      }
      if (!is.null(infer_payload$reject)) {
        inference_attr$reject <- infer_payload$reject
      }
      if (!is.null(infer_payload$critical_p_value)) {
        inference_attr$critical_p_value <- infer_payload$critical_p_value
      }
    }
    if (identical(p_adjust, "ecp")) inference_attr$n_mc <- n_mc
    attr(out, "inference") <- inference_attr
  }
  out
}

.mc_skipcor_ci_attr <- function(x) {
  ci <- attr(x, "ci", exact = TRUE)
  if (!is.null(ci)) return(ci)

  lwr <- attr(x, "lwr.ci", exact = TRUE)
  upr <- attr(x, "upr.ci", exact = TRUE)
  if (is.null(lwr) && is.null(upr)) return(NULL)

  list(
    est = unclass(as.matrix(x)),
    lwr.ci = lwr,
    upr.ci = upr,
    conf.level = attr(x, "conf.level", exact = TRUE)
  )
}

.mc_skipcor_inference_attr <- function(x) {
  inf <- attr(x, "inference", exact = TRUE)
  if (!is.null(inf)) return(inf)

  legacy <- list(
    method = attr(x, "inference_method", exact = TRUE),
    n_boot = attr(x, "n_boot", exact = TRUE),
    conf.level = attr(x, "conf.level", exact = TRUE),
    fwe_level = attr(x, "fwe_level", exact = TRUE),
    p_adjust = attr(x, "p_adjust", exact = TRUE),
    p_value = attr(x, "p_value", exact = TRUE),
    p_value_adjusted = attr(x, "p_value_adjusted", exact = TRUE),
    reject = attr(x, "reject", exact = TRUE),
    critical_p_value = attr(x, "critical_p_value", exact = TRUE)
  )
  if (all(vapply(legacy, is.null, logical(1)))) return(NULL)
  legacy$ci <- .mc_skipcor_ci_attr(x)
  legacy
}

.mc_skipcor_component_names <- function(x) {
  out <- c("estimate")
  if (!is.null(.mc_skipcor_ci_attr(x))) out <- c(out, "ci")
  inf <- .mc_skipcor_inference_attr(x)
  if (!is.null(inf)) {
    if (!is.null(inf$p_value)) out <- c(out, "p_value")
    if (!is.null(inf$p_value_adjusted)) out <- c(out, "p_value_adjusted")
    if (!is.null(inf$reject)) out <- c(out, "reject")
    if (!is.null(inf$critical_p_value)) out <- c(out, "critical_p_value")
    out <- c(out, "inference")
  }
  if (!is.null(attr(x, "diagnostics", exact = TRUE))) out <- c(out, "diagnostics")
  if (!is.null(attr(x, "skipped_masks", exact = TRUE))) out <- c(out, "skipped_masks")
  unique(out)
}

#' @export
names.skipped_corr <- function(x) {
  .mc_skipcor_component_names(x)
}

#' @export
`$.skipped_corr` <- function(x, name) {
  ci <- .mc_skipcor_ci_attr(x)
  inf <- .mc_skipcor_inference_attr(x)
  switch(name,
    estimate = unclass(as.matrix(x)),
    est = unclass(as.matrix(x)),
    cor = unclass(as.matrix(x)),
    ci = ci,
    p_value = if (is.null(inf)) NULL else inf$p_value,
    p_value_adjusted = if (is.null(inf)) NULL else inf$p_value_adjusted,
    reject = if (is.null(inf)) NULL else inf$reject,
    critical_p_value = if (is.null(inf)) NULL else inf$critical_p_value,
    inference = inf,
    diagnostics = attr(x, "diagnostics", exact = TRUE),
    skipped_masks = attr(x, "skipped_masks", exact = TRUE),
    NULL
  )
}

#' @export
`[[.skipped_corr` <- function(x, i, ...) {
  if (is.numeric(i)) {
    nms <- names(x)
    if (length(i) != 1L || is.na(i) || i < 1L || i > length(nms)) {
      abort_bad_arg("i", message = "must be a valid component index.")
    }
    i <- nms[[i]]
  }
  `$.skipped_corr`(x, i)
}

#' @rdname skipped_corr
#' @export
skipped_corr_masks <- function(x, var1 = NULL, var2 = NULL) {
  check_inherits(x, "skipped_corr")
  masks <- attr(x, "skipped_masks", exact = TRUE)
  if (is.null(masks)) return(NULL)
  if (is.null(var1) && is.null(var2)) return(masks)
  if (is.null(var1) || is.null(var2)) {
    abort_bad_arg("var1", message = "and {.arg var2} must both be supplied when extracting a specific pair mask.")
  }

  resolve_one <- function(var, arg) {
    if (is.character(var)) {
      if (length(var) != 1L || is.na(var)) {
        abort_bad_arg(arg, message = "must be a single non-missing column name or index.")
      }
      idx <- match(var, colnames(x))
      if (is.na(idx)) {
        abort_bad_arg(arg, message = "must match a column name in {.arg x}.")
      }
      idx
    } else {
      idx <- check_scalar_int_pos(var, arg = arg)
      if (idx > ncol(x)) {
        abort_bad_arg(arg, message = "must be <= ncol(x).")
      }
      idx
    }
  }

  i <- resolve_one(var1, "var1")
  j <- resolve_one(var2, "var2")
  if (i == j) return(integer())
  lo <- min(i, j)
  hi <- max(i, j)
  hit <- which(masks$pair_i == lo & masks$pair_j == hi)
  if (!length(hit)) integer() else masks$skipped_rows[[hit[[1L]]]]
}

#' @rdname skipped_corr
#' @method print skipped_corr
#' @param ci_digits Integer; digits for skipped-correlation confidence limits.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param show_p One of \code{"auto"}, \code{"yes"}, \code{"no"}. For
#'   \code{print()}, \code{"auto"} keeps the compact matrix-only display;
#'   use \code{"yes"} to also print pairwise p-values.
#' @export
print.skipped_corr <- function(x,
                               digits = 4,
                               n = NULL,
                               topn = NULL,
                               max_vars = NULL,
                               width = NULL,
                               ci_digits = 4,
                               show_ci = NULL,
                               show_p = c("auto", "yes", "no"),
                               ...) {
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("print_show_ci", "yes")
  )
  show_p <- match.arg(show_p)
  .mc_print_corr_matrix(
    x,
    header = "Skipped correlation matrix",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
  ci <- .mc_skipcor_ci_attr(x)
  inf <- .mc_skipcor_inference_attr(x)

  include_ci <- identical(show_ci, "yes") && !is.null(ci)
  include_p <- identical(show_p, "yes") && !is.null(inf) && !is.null(inf$p_value)
  if (include_p && !include_ci) {
    cli::cli_inform("Use {.code summary(x, show_ci = \"no\")} for bounded pairwise p-value output.")
  }
  invisible(x)
}

#' @rdname skipped_corr
#' @method plot skipped_corr
#' @param ci_text_size Text size for confidence intervals in the heatmap.
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
#' @export
plot.skipped_corr <- function(x,
                         title = "Skipped correlation heatmap",
                         low_color = "indianred1",
                         high_color = "steelblue1",
                         mid_color = "white",
                         value_text_size = 4,
                         ci_text_size = 3,
                         show_value = TRUE,
                         ...) {
  check_inherits(x, "skipped_corr")
  check_bool(show_value, arg = "show_value")

  est_mat <- as.matrix(x)
  df_est <- as.data.frame(as.table(est_mat))
  names(df_est) <- c("Var1", "Var2", "skipped_corr")

  ci <- .mc_skipcor_ci_attr(x)
  if (!is.null(ci) && !is.null(ci$lwr.ci) && !is.null(ci$upr.ci)) {
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
      sprintf("(%.2f, %.2f)", df$lwr, df$upr)
    )
  } else {
    df <- df_est
    df$ci_label <- NA_character_
  }

  lev_row <- unique(df_est$Var1)
  lev_col <- unique(df_est$Var2)
  df$Var1 <- factor(df$Var1, levels = rev(lev_row))
  df$Var2 <- factor(df$Var2, levels = lev_col)
  df$label <- sprintf("%.2f", df$skipped_corr)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Var2, y = Var1, fill = skipped_corr)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(
      low = low_color,
      high = high_color,
      mid = mid_color,
      midpoint = 0,
      limits = c(-1, 1),
      name = "skipped_corr"
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

.mc_skipcor_pairwise_summary <- function(object,
                                         digits = 4,
                                         ci_digits = 2,
                                         p_digits = 4,
                                         show_ci = NULL,
                                         show_p = c("auto", "yes", "no")) {
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )
  show_p <- match.arg(show_p)
  check_inherits(object, "skipped_corr")

  est <- as.matrix(object)
  rn <- rownames(est); cn <- colnames(est)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(est)))
  if (is.null(cn)) cn <- as.character(seq_len(ncol(est)))

  ci <- .mc_skipcor_ci_attr(object)
  inf <- .mc_skipcor_inference_attr(object)
  diag_attr <- attr(object, "diagnostics", exact = TRUE)

  include_ci <- identical(show_ci, "yes") && !is.null(ci)
  include_p <- switch(show_p, auto = !is.null(inf) && !is.null(inf$p_value), yes = TRUE, no = FALSE)

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
      if (include_ci) {
        rec$lwr <- if (!is.null(ci$lwr.ci) && is.finite(ci$lwr.ci[i, j])) round(ci$lwr.ci[i, j], ci_digits) else NA_real_
        rec$upr <- if (!is.null(ci$upr.ci) && is.finite(ci$upr.ci[i, j])) round(ci$upr.ci[i, j], ci_digits) else NA_real_
      }
      if (include_p) {
        rec$p_value <- if (!is.null(inf$p_value) && is.finite(inf$p_value[i, j])) round(inf$p_value[i, j], p_digits) else NA_real_
        if (!is.null(inf$p_value_adjusted)) {
          rec$p_value_adjusted <- if (is.finite(inf$p_value_adjusted[i, j])) round(inf$p_value_adjusted[i, j], p_digits) else NA_real_
        }
        if (!is.null(inf$reject)) rec$reject <- isTRUE(inf$reject[i, j])
      }
      if (is.list(diag_attr)) {
        if (is.matrix(diag_attr$skipped_n)) rec$skipped_n <- as.integer(diag_attr$skipped_n[i, j])
        if (is.matrix(diag_attr$skipped_prop)) rec$skipped_prop <- round(diag_attr$skipped_prop[i, j], p_digits)
        if (is.matrix(diag_attr$n_complete)) rec$n_complete <- as.integer(diag_attr$n_complete[i, j])
      }
      rows[[k]] <- rec
    }
  }

  df <- do.call(rbind.data.frame, rows)
  rownames(df) <- NULL
  num_cols <- intersect(c("estimate", "lwr", "upr", "p_value", "p_value_adjusted", "skipped_prop"), names(df))
  int_cols <- intersect(c("skipped_n", "n_complete"), names(df))
  for (nm in num_cols) df[[nm]] <- as.numeric(df[[nm]])
  for (nm in int_cols) df[[nm]] <- as.integer(df[[nm]])

  base <- .mc_summary_corr_matrix(object)
  out <- structure(df, class = c("summary.skipped_corr", "data.frame"))
  attr(out, "overview") <- base
  attr(out, "has_ci") <- include_ci
  attr(out, "has_p") <- include_p
  attr(out, "conf.level") <- if (is.null(ci)) NA_real_ else ci$conf.level
  attr(out, "digits") <- digits
  attr(out, "ci_digits") <- ci_digits
  attr(out, "p_digits") <- p_digits
  attr(out, "inference_method") <- if (is.null(inf)) NA_character_ else inf$method
  attr(out, "p_adjust") <- if (is.null(inf) || is.null(inf$p_adjust)) "none" else inf$p_adjust
  attr(out, "critical_p_value") <- if (is.null(inf)) NA_real_ else inf$critical_p_value %||% NA_real_
  out
}

#' @rdname skipped_corr
#' @method summary skipped_corr
#' @param object An object of class \code{skipped_corr}.
#' @export
summary.skipped_corr <- function(object, n = NULL, topn = NULL,
                                 max_vars = NULL, width = NULL,
                                 show_ci = NULL, ...) {
  check_inherits(object, "skipped_corr")
  ci <- .mc_skipcor_ci_attr(object)
  inf <- .mc_skipcor_inference_attr(object)

  if (is.null(ci) && (is.null(inf) || is.null(inf$p_value))) {
    return(.mc_summary_corr_matrix(object, topn = topn))
  }

  .mc_skipcor_pairwise_summary(
    object,
    show_ci = .mc_validate_yes_no(
      show_ci,
      arg = "show_ci",
      default = .mc_display_option("summary_show_ci", "yes")
    ),
    show_p = "auto"
  )
}

#' @rdname skipped_corr
#' @method print summary.skipped_corr
#' @param x An object of class \code{summary.skipped_corr}.
#' @export
print.summary.skipped_corr <- function(x, digits = NULL, n = NULL,
                                       topn = NULL, max_vars = NULL,
                                       width = NULL, show_ci = NULL, ...) {
  has_ci <- isTRUE(attr(x, "has_ci"))
  has_p <- isTRUE(attr(x, "has_p"))
  p_digits <- attr(x, "p_digits"); if (!is.numeric(p_digits)) p_digits <- 4
  extra_items <- c(
    if (has_p) c(inference = attr(x, "inference_method")),
    if (has_p) c(multiplicity = attr(x, "p_adjust"))
  )
  crit <- suppressWarnings(as.numeric(attr(x, "critical_p_value")))
  if (is.finite(crit)) {
    extra_items <- c(extra_items, critical_p = format(signif(crit, digits = p_digits)))
  }
  .mc_print_pairwise_summary_digest(
    x,
    title = "Skipped correlation summary",
    digits = .mc_coalesce(digits, 4),
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ci_method = if (has_ci) "bootstrap_b2" else NULL,
    extra_items = extra_items,
    ...
  )
  invisible(x)
}
