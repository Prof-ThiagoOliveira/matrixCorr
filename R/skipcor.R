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
#' @param x An object of class \code{skipped_corr}.
#' @param var1,var2 Optional column names or 1-based column indices used by
#'   [skipped_corr_masks()] to extract the skipped-row indices for one pair.
#' @param digits Integer; number of digits to print.
#' @param max_rows,max_cols Optional integers limiting the printed matrix size.
#' @param ... Additional arguments passed to the underlying print or plot helper.
#' @param title Character; plot title.
#' @param low_color,high_color,mid_color Colors used in the heatmap.
#' @param value_text_size Numeric text size for overlaid cell values.
#'
#' @return A symmetric correlation matrix with class \code{skipped_corr} and
#'   attributes \code{method = "skipped_correlation"}, \code{description}, and
#'   \code{package = "matrixCorr"}. When \code{return_masks = TRUE}, the matrix
#'   also carries a \code{skipped_masks} attribute containing compact pairwise
#'   skipped-row indices.
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
#' @references
#' Wilcox, R. R. (2004). Inferences based on a skipped correlation coefficient.
#' Journal of Applied Statistics, 31(2), 131-143.
#' \doi{10.1080/0266476032000148821}
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
                    return_masks = FALSE) {
  method <- match.arg(method)
  na_method <- match.arg(na_method)
  outlier_rule <- match.arg(outlier_rule)
  check_bool(stand, arg = "stand")
  check_bool(return_masks, arg = "return_masks")
  check_scalar_numeric(cutoff, arg = "cutoff", lower = 0, closed_lower = FALSE)
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")

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
    return_masks = return_masks
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

  res <- res$cor
  colnames(res) <- rownames(res) <- colnames_data
  out <- .mc_structure_corr_matrix(
    res,
    class_name = "skipped_corr",
    method = "skipped_correlation",
    description = paste0(
      "Skipped correlation; base = ", method,
      "; rule = ", outlier_rule,
      "; standardise = ", stand,
      "; NA mode = ", na_method, "."
    )
  )
  if (isTRUE(return_masks)) attr(out, "skipped_masks") <- mask_payload
  out
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
#' @export
print.skipped_corr <- function(x, digits = 4, max_rows = NULL, max_cols = NULL, ...) {
  .mc_print_corr_matrix(
    x,
    header = "Skipped correlation matrix:",
    digits = digits,
    max_rows = max_rows,
    max_cols = max_cols,
    ...
  )
}

#' @rdname skipped_corr
#' @method plot skipped_corr
#' @export
plot.skipped_corr <- function(x,
                         title = "Skipped correlation heatmap",
                         low_color = "indianred1",
                         high_color = "steelblue1",
                         mid_color = "white",
                         value_text_size = 4,
                       ...) {
  .mc_plot_corr_matrix(
    x,
    class_name = "skipped_corr",
    fill_name = "skipped_corr",
    title = title,
    low_color = low_color,
    high_color = high_color,
    mid_color = mid_color,
    value_text_size = value_text_size,
    ...
  )
}

#' @rdname skipped_corr
#' @method summary skipped_corr
#' @param object An object of class \code{skipped_corr}.
#' @export
summary.skipped_corr <- function(object, ...) {
  .mc_summary_corr_matrix(object)
}
