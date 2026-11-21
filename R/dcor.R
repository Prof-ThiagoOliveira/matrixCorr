#' @title Pairwise Distance Correlation (dCor)
#'
#' @description
#' Computes pairwise distance correlation for numeric columns using
#' a high-performance C++ backend. Population distance correlation is a nonnegative
#' measure of dependence that can capture nonlinear and nonmonotonic
#' relationships. By default, `dcor()` returns the conventional sample distance
#' correlation \eqn{R_n}; set `squared = TRUE` to return \eqn{R_n^2}.
#'
#' @param data A numeric matrix or a data frame with at least two numeric
#' columns. All non-numeric columns are dropped. Columns must be numeric.
#' @param squared Logical; if `FALSE` (default), return the correlation-scale estimate. If `TRUE`,
#'   return the squared-scale result described below.
#' @param estimator One of `v` (conventional, default), `u` (U-centred), or
#'   `improved` (bootstrap-weighted U/V combination).
#' @param B Number of smoothed bootstrap replicates per pair for `improved`
#'   (at least two). Use `set.seed()` for reproducibility.
#' @param bandwidth Positive smoothing bandwidth in observation units, either
#'   shared or one per retained numeric column. By default, `stats::bw.nrd0`
#'   is applied to each margin on the rows used for that pair.
#' @param na_method Character scalar controlling missing-data handling.
#'   \code{"error"} rejects missing, \code{NaN}, and infinite values.
#'   \code{"pairwise"} recomputes each association on its own pairwise
#'   complete-case overlap. \code{"complete"} performs listwise deletion once
#'   across the retained numeric columns.
#' @param p_value Logical (default \code{FALSE}). Compatibility option. If
#'   \code{TRUE}, attach the Székely-Rizzo t-test metadata computed from the
#'   signed bias-corrected statistic returned by `bcdcor()`, not from the
#'   estimates in the main matrix.
#' @param n_threads Integer \eqn{\geq 1}. Number of OpenMP threads. Defaults to
#'   \code{getOption("matrixCorr.threads", 1L)}.
#' @param output Output representation for the computed estimates:
#'   \code{"matrix"}, \code{"sparse"}, or \code{"edge_list"}.
#' @param threshold Non-negative absolute-value filter for non-matrix outputs:
#'   keep entries with \code{abs(value) >= threshold}. Must be \code{0} when
#'   \code{output = "matrix"}.
#' @param diag Logical; whether to include diagonal entries in
#'   \code{"sparse"} and \code{"edge_list"} outputs.
#' @param ... Compatibility arguments, including deprecated `check_na`.
#'
#' @return A symmetric correlation result where the \code{(i, j)} entry is
#' the selected estimate. The `estimator` attribute identifies the method.
#' For improved estimates, `attr(result, "diagnostics")$lambda`
#' gives the weight on U for each pair (NA where no weight was fitted). When `p_value = TRUE`, the `inference` attribute contains
#' `bcdcor`, `statistic`, `parameter`, and `p_value`; the t-test uses the
#' signed U-centred statistic in `inference$bcdcor`.
#'
#' @details
#' For distances \eqn{a_{ij}=|x_i-x_j|}, define the conventional doubly centred
#' distances
#' \deqn{A_{ij}=a_{ij}-\bar a_{i\cdot}-\bar a_{\cdot j}+\bar a_{\cdot\cdot}.}
#' Define \eqn{B_{ij}} analogously for \eqn{y}. The sample squared distance
#' covariance is
#' \deqn{V_n^2(x,y)=\frac{1}{n^2}\sum_{i,j} A_{ij}B_{ij},}
#' with \eqn{V_n^2(x,x)} and \eqn{V_n^2(y,y)} defined analogously. The
#' conventional squared distance correlation is
#' \deqn{R_n^2(x,y)=
#' \frac{V_n^2(x,y)}{\sqrt{V_n^2(x,x)V_n^2(y,y)}}.}
#' By default, `dcor()` returns \eqn{R_n=\sqrt{R_n^2}} unless `squared = TRUE`.
#'
#' For `estimator = "u"`, the squared-scale estimate is the signed ratio
#' returned by `bcdcor()`. On the correlation scale we use its signed square
#' root, following the software convention used in the authors' simulations.
#' Unbiasedness applies to squared distance covariance before normalization.
#'
#' For `estimator = "improved"`, Algorithm 1 of Monroy-Castillo et al. (2025)
#' resamples paired rows with replacement and adds independent Gaussian noise
#' to each margin. Bootstrap biases use the observed V estimate as reference.
#' We use sample variances and covariances (divisor B - 1) in equation (8)
#' and constrain the weight on U to the unit interval. When bootstrap
#' moments are undefined the estimate is NA. When the criterion is flat we
#' use V. The combination can be negative. With `squared = TRUE` we return
#' the ordinary square of the fitted combination. This argument changes the
#' reporting scale only for the improved estimator.
#'
#' Under the required finite first-moment conditions, population distance
#' correlation is zero if and only if the variables are independent. This
#' population property should not be read as a sample-level test by itself.
#'
#' @section Difference between `dcor()` and `bcdcor()`:
#' By default, `dcor()` returns conventional nonnegative distance correlation \eqn{R_n}.
#' `dcor(squared = TRUE)` returns conventional \eqn{R_n^2}. `bcdcor()` returns
#' a U-centred bias-corrected squared distance-correlation statistic that may be
#' negative in finite samples and is used by the bias-corrected dCor t-test. In
#' finite samples, `bcdcor(x)` is not generally equal to `dcor(x)^2`.
#'
#' @note Requires \eqn{n \ge 4}. Columns with zero distance variance yield
#' \code{NA} in their off-diagonal row/column entries.
#'
#' @references
#' Szekely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007).
#' Measuring and testing dependence by correlation of distances.
#' \emph{Annals of Statistics}, 35(6), 2769-2794.
#'
#' Szekely, G. J., & Rizzo, M. L. (2013).
#' The distance correlation t-test of independence.
#' \emph{Journal of Multivariate Analysis}, 117, 193-213.
#'
#' Monroy-Castillo, B. E., Jácome, M. A., & Cao, R. (2025). Improved distance
#' correlation estimation. Applied Intelligence, 55, 263. \doi{10.1007/s10489-024-05940-x}.
#'
#' @seealso [bcdcor()], [robust_dcor()]
#'
#' @examples
#' # Is eruption duration associated with waiting time to the next eruption?
#' X <- datasets::faithful
#' dcor(X)
#' dcor(X, squared = TRUE)
#' dcor(X, estimator = "u")
#'
#' \donttest{
#' set.seed(123)
#' n <- 300
#' x <- rnorm(n)
#' y <- x^2 + rnorm(n, sd = 0.25)
#' X <- cbind(x = x, y = y)
#'
#' cor(x, y)
#' dcor(X)
#'
#' # Pairwise distance-correlation matrix across columns
#' if (requireNamespace("MASS", quietly = TRUE)) {
#'   set.seed(7)
#'   p <- 5
#'   n <- 150
#'   rho <- 0.6
#'   Sigma <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
#'
#'   X <- MASS::mvrnorm(
#'     n,
#'     mu = rep(0, p),
#'     Sigma = Sigma
#'   )
#'   colnames(X) <- paste0("V", seq_len(p))
#'
#'   D <- dcor(X)
#'   D[1:3, 1:3]
#' }
#'
#' # Backward-compatible inference through dcor()
#' fit <- dcor(cbind(x = x, y = y), p_value = TRUE)
#' fit
#'
#' # The displayed matrix contains conventional dCor.
#' # The attached t-test is based on the signed bias-corrected statistic.
#' inf <- attr(fit, "inference")
#' inf$bcdcor
#' inf$p_value
#'
#' # Missing-data modes
#' X_miss <- cbind(
#'   a = c(1, 2, NA, 4, 5, 6),
#'   b = c(2, 3, 4, 5, NA, 7),
#'   c = c(6, 5, 4, 3, 2, 1)
#' )
#'
#' dcor(X_miss, na_method = "complete")
#' dcor(X_miss, na_method = "pairwise")
#'
#' # Edge-list output for larger matrices
#' set.seed(11)
#' X <- matrix(rnorm(500), nrow = 100, ncol = 5)
#' colnames(X) <- paste0("V", seq_len(ncol(X)))
#'
#' dcor(
#'   X,
#'   output = "edge_list",
#'   threshold = 0.2,
#'   diag = FALSE
#' )
#' }
#'
#' @author Thiago de Paula Oliveira
#'
#' @export
dcor <- function(data,
                 squared = FALSE,
                 na_method = c("error", "pairwise", "complete"),
                 p_value = FALSE,
                 n_threads = getOption("matrixCorr.threads", 1L),
                 output = c("matrix", "sparse", "edge_list"),
                 threshold = 0,
                 diag = TRUE,
                 ...,
                 estimator = c("v", "u", "improved"),
                 B = 1000L,
                 bandwidth = NULL) {
  estimator <- match.arg(estimator)
  check_bool(squared, arg = "squared")
  output_cfg <- .mc_prepare_corr_output(
    output = output,
    threshold = threshold,
    diag = diag
  )
  if (estimator == "v" && ...length() == 0L && missing(na_method) && isFALSE(p_value)) {
    input <- .mc_prepare_corr_input(
      data,
      na_cfg = list(na_method = "error", check_na = TRUE),
      min_n = 4L
    )
    return(.mc_with_omp_threads(
      n_threads,
      {
        .mc_finalize_corr_result(
          mat = dcor_matrix_cpp(input$data, squared = squared),
          class_name = "dcor",
          method = "distance_correlation",
          description = .mc_dcor_description(squared),
          output_cfg = output_cfg,
          dimnames = input$dimnames,
          symmetric = TRUE,
          extra_attrs = list(squared = squared, estimator = estimator)
        )
      }
    ))
  }

  na_cfg <- .mc_resolve_corr_na(
    na_method = na_method,
    dots = list(...),
    na_method_missing = missing(na_method),
    allowed = c("error", "pairwise", "complete")
  )
  if (!isFALSE(p_value)) {
    check_bool(p_value, arg = "p_value")
  } else if (!is.logical(p_value) || length(p_value) != 1L || is.na(p_value)) {
    check_bool(p_value, arg = "p_value")
  }
  input <- .mc_prepare_corr_input(data, na_cfg = na_cfg, min_n = 4L)
  diagnostics <- NULL
  inference_attr <- NULL

  .mc_with_omp_threads(
    n_threads,
    {
      pairwise <- identical(na_cfg$na_method, "pairwise")
      if (estimator != "v" || p_value) {
        if (pairwise || p_value) {
          bc_pairwise <- bcdcor_matrix_pairwise_cpp(input$data, return_inference = p_value)
          u2 <- if (p_value) bc_pairwise$estimate else bc_pairwise$est
        } else {
          u2 <- bcdcor_matrix_cpp(input$data)
        }
      }
      if (estimator != "u") {
        v_squared <- estimator == "v" && squared
        if (pairwise) {
          v_pair <- dcor_matrix_pairwise_cpp(input$data, squared = v_squared)
          dcor_matrix <- v_pair$est
        } else {
          dcor_matrix <- dcor_matrix_cpp(input$data, squared = v_squared)
        }
      }
      if (pairwise) {
        counts <- if (estimator == "v") v_pair$n_complete else bc_pairwise$n_complete
        diagnostics <- list(n_complete = .mc_set_matrix_dimnames(counts, input$colnames))
      }
      if (estimator == "u") {
        dcor_matrix <- if (squared) u2 else sign(u2) * sqrt(abs(u2))
      }
      if (estimator == "improved") {
        B <- check_scalar_int_pos(B, arg = "B")
        if (B < 2L) cli::cli_abort("{.arg B} must be at least two.")
        p <- ncol(input$data)
        if (!is.null(bandwidth) && (!is.numeric(bandwidth) ||
            !length(bandwidth) %in% c(1L, p) || any(!is.finite(bandwidth) | bandwidth <= 0))) {
          cli::cli_abort("{.arg bandwidth} must be positive and finite, shared or one per column.")
        }
        if (!is.null(bandwidth)) bandwidth <- rep_len(bandwidth, p)
        lambda <- matrix(NA_real_, p, p, dimnames = input$dimnames)
        for (j in seq.int(2L, p)) {
          for (i in seq_len(j - 1L)) {
            v <- dcor_matrix[i, j]
            u <- sign(u2[i, j]) * sqrt(abs(u2[i, j]))
            value <- NA_real_
            if (is.finite(v) && is.finite(u)) {
              keep <- is.finite(input$data[, i]) & is.finite(input$data[, j])
              weight <- .mc_dcor_weight(input$data[keep, i], input$data[keep, j],
                                        v, B, bandwidth[c(i, j)])
              lambda[i, j] <- lambda[j, i] <- weight
              value <- weight * u + (1 - weight) * v
            }
            dcor_matrix[i, j] <- dcor_matrix[j, i] <- value
          }
        }
        diagnostics$lambda <- lambda
        if (squared) dcor_matrix <- dcor_matrix^2
      }

      if (isTRUE(p_value)) {
        diagnostics <- .mc_merge_diagnostics(
          diagnostics,
          list(n_complete = .mc_set_matrix_dimnames(unclass(bc_pairwise$n_complete), input$colnames))
        )
        inference_attr <- list(
          method = "bcdcor_t_test",
          bcdcor = .mc_set_matrix_dimnames(unclass(bc_pairwise$estimate), input$colnames),
          statistic = .mc_set_matrix_dimnames(unclass(bc_pairwise$statistic), input$colnames),
          parameter = .mc_set_matrix_dimnames(unclass(bc_pairwise$parameter), input$colnames),
          p_value = .mc_set_matrix_dimnames(unclass(bc_pairwise$p_value), input$colnames),
          alternative = "greater"
        )
      }

      .mc_finalize_corr_result(
        mat = dcor_matrix,
        class_name = "dcor",
        method = "distance_correlation",
        description = .mc_dcor_description(squared, estimator),
        output_cfg = output_cfg,
        diagnostics = .mc_merge_diagnostics(diagnostics, input$diagnostics),
        dimnames = input$dimnames,
        symmetric = TRUE,
        extra_attrs = c(
          list(squared = squared, estimator = estimator),
          if (!is.null(inference_attr)) list(inference = inference_attr)
        )
      )
    }
  )
}

.mc_dcor_weight <- function(x, y, v, B, bandwidth) {
  if (is.null(bandwidth)) bandwidth <- c(stats::bw.nrd0(x), stats::bw.nrd0(y))
  n <- length(x)
  u_boot <- v_boot <- numeric(B)
  for (b in seq_len(B)) {
    rows <- sample.int(n, n, replace = TRUE)
    xb <- x[rows] + bandwidth[1L] * stats::rnorm(n)
    yb <- y[rows] + bandwidth[2L] * stats::rnorm(n)
    u2 <- bcdcor_pair_cpp(xb, yb)
    u_boot[b] <- sign(u2) * sqrt(abs(u2))
    v_boot[b] <- dcor_pair_cpp(xb, yb)
  }
  if (any(!is.finite(u_boot)) || any(!is.finite(v_boot))) return(NA_real_)
  # Equation (8), collecting differences before estimating the moments.
  delta <- v_boot - u_boot
  bias_diff <- mean(delta)
  denominator <- stats::var(delta) + bias_diff^2
  if (denominator == 0) return(0)
  numerator <- stats::cov(v_boot, delta) + (mean(v_boot) - v) * bias_diff
  max(0, min(1, numerator / denominator))
}

.mc_dcor_description <- function(squared = FALSE, estimator = "v") {
  if (estimator != "v") {
    scale <- if (squared) "squared distance correlation" else "distance correlation"
    return(sprintf("Pairwise %s matrix (%s estimator)", scale, estimator))
  }
  if (isTRUE(squared)) {
    "Pairwise squared distance correlation matrix"
  } else {
    "Pairwise distance correlation matrix"
  }
}

.mc_bcdcor_description <- function() {
  "Pairwise bias-corrected squared distance correlation matrix"
}

#' @title Pairwise Bias-Corrected Squared Distance Correlation
#'
#' @description
#' Computes the U-centred bias-corrected squared distance-correlation statistic
#' for numeric columns using a high-performance C++ backend. Unlike
#' conventional distance correlation, this finite-sample statistic may be
#' negative. It is the statistic used by the Székely-Rizzo bias-corrected
#' distance-correlation t-test.
#'
#' @inheritParams dcor
#' @param p_value Logical (default \code{FALSE}). If \code{TRUE}, attach
#'   pairwise p-values, test statistics, and degrees of freedom from the
#'   bias-corrected distance-correlation t-test.
#' @param low_color,mid_color,high_color Colours used in the `bcdcor()` heatmap.
#'
#' @return A symmetric correlation result where the \code{(i, j)} entry is the
#' signed U-centred bias-corrected squared distance-correlation statistic. The
#' object has class \code{bcdcor}. When `p_value = TRUE`, the `inference`
#' attribute contains `bcdcor`, `statistic`, `parameter`, and `p_value`.
#'
#' @details
#' Let \eqn{A^{(x)}} and \eqn{A^{(y)}} be U-centred distance matrices. The
#' U-centred construction provides an unbiased estimator of squared distance
#' covariance,
#' \deqn{V_{n,U}^2(x,y)=\frac{1}{n(n-3)}
#'       \sum_{i \ne j} A^{(x)}_{ij} A^{(y)}_{ij}.}
#' The bias-corrected squared distance-correlation statistic is
#' \deqn{R_n^\ast(x,y)=
#' \frac{V_{n,U}^2(x,y)}
#'      {\sqrt{V_{n,U}^2(x,x)V_{n,U}^2(y,y)}}.}
#' This ratio is signed and may be negative in finite samples.
#'
#' When `p_value = TRUE`, the t statistic is
#' \deqn{T = \sqrt{M - 1}\frac{R_n^\ast}
#'        {\sqrt{1 - (R_n^\ast)^2}}, \quad
#'        M=\frac{n(n-3)}{2},}
#' referenced to a Student \eqn{t} distribution with \eqn{M-1} degrees of
#' freedom. The signed statistic is not clipped before inference.
#'
#' @section Difference between `dcor()` and `bcdcor()`:
#' By default, `dcor()` returns conventional nonnegative distance correlation \eqn{R_n}.
#' `dcor(squared = TRUE)` returns conventional \eqn{R_n^2}. `bcdcor()` returns
#' the U-centred bias-corrected squared statistic used by the t-test. In finite
#' samples, `bcdcor(x)` is not generally equal to `dcor(x)^2`.
#'
#' @references
#' Szekely, G. J., & Rizzo, M. L. (2013).
#' The distance correlation t-test of independence.
#' \emph{Journal of Multivariate Analysis}, 117, 193-213.
#'
#' @seealso [dcor()], [robust_dcor()]
#'
#' @examples
#' # Is eruption duration associated with waiting time to the next eruption?
#' X <- datasets::faithful
#' dcor(X, squared = TRUE)
#' bcdcor(X)
#'
#' \donttest{
#' set.seed(3)
#' x <- rnorm(80)
#' y <- rnorm(80)
#' X <- cbind(x = x, y = y)
#'
#' bcdcor(X)
#'
#' # Bias-corrected dCor t-test
#' set.seed(2026)
#' n <- 80
#' x <- rnorm(n)
#' y <- x^2 + rnorm(n, sd = 0.25)
#' X <- cbind(x = x, y = y)
#'
#' fit <- bcdcor(X, p_value = TRUE)
#'
#' fit
#' summary(fit)
#'
#' inf <- attr(fit, "inference")
#' inf$bcdcor
#' inf$statistic
#' inf$parameter
#' inf$p_value
#'
#' # bcdcor() supports the same missing-data modes as dcor().
#' X_miss <- cbind(
#'   a = c(1, 2, NA, 4, 5, 6),
#'   b = c(2, 3, 4, 5, NA, 7),
#'   c = c(6, 5, 4, 3, 2, 1)
#' )
#'
#' bcdcor(X_miss, na_method = "complete")
#' bcdcor(X_miss, na_method = "pairwise")
#'
#' # Edge-list thresholding uses absolute magnitude and preserves signs.
#' set.seed(3)
#' X <- cbind(x = rnorm(80), y = rnorm(80), z = rnorm(80))
#'
#' bcdcor(
#'   X,
#'   output = "edge_list",
#'   threshold = 0.001,
#'   diag = FALSE
#' )
#' }
#'
#' @author Thiago de Paula Oliveira
#'
#' @export
bcdcor <- function(data,
                   na_method = c("error", "pairwise", "complete"),
                   p_value = FALSE,
                   n_threads = getOption("matrixCorr.threads", 1L),
                   output = c("matrix", "sparse", "edge_list"),
                   threshold = 0,
                   diag = TRUE,
                   ...) {
  output_cfg <- .mc_prepare_corr_output(
    output = output,
    threshold = threshold,
    diag = diag
  )
  if (...length() == 0L && missing(na_method) && isFALSE(p_value)) {
    input <- .mc_prepare_corr_input(
      data,
      na_cfg = list(na_method = "error", check_na = TRUE),
      min_n = 4L
    )
    return(.mc_with_omp_threads(
      n_threads,
      {
        .mc_finalize_corr_result(
          mat = bcdcor_matrix_cpp(input$data),
          class_name = "bcdcor",
          method = "bias_corrected_distance_correlation",
          description = .mc_bcdcor_description(),
          output_cfg = output_cfg,
          dimnames = input$dimnames,
          symmetric = TRUE
        )
      }
    ))
  }

  na_cfg <- .mc_resolve_corr_na(
    na_method = na_method,
    dots = list(...),
    na_method_missing = missing(na_method),
    allowed = c("error", "pairwise", "complete")
  )
  if (!isFALSE(p_value)) {
    check_bool(p_value, arg = "p_value")
  } else if (!is.logical(p_value) || length(p_value) != 1L || is.na(p_value)) {
    check_bool(p_value, arg = "p_value")
  }
  input <- .mc_prepare_corr_input(data, na_cfg = na_cfg, min_n = 4L)
  diagnostics <- NULL
  inference_attr <- NULL

  .mc_with_omp_threads(
    n_threads,
    {
      if (isTRUE(p_value) || identical(na_cfg$na_method, "pairwise")) {
        pairwise <- bcdcor_matrix_pairwise_cpp(
          input$data,
          return_inference = p_value
        )
        bc_matrix <- pairwise$est
        diagnostics <- list(
          n_complete = .mc_set_matrix_dimnames(unclass(pairwise$n_complete), input$colnames)
        )

        if (isTRUE(p_value)) {
          inference_attr <- list(
            method = "bcdcor_t_test",
            bcdcor = .mc_set_matrix_dimnames(unclass(pairwise$estimate), input$colnames),
            statistic = .mc_set_matrix_dimnames(unclass(pairwise$statistic), input$colnames),
            parameter = .mc_set_matrix_dimnames(unclass(pairwise$parameter), input$colnames),
            p_value = .mc_set_matrix_dimnames(unclass(pairwise$p_value), input$colnames),
            alternative = "greater"
          )
        }
      } else {
        bc_matrix <- bcdcor_matrix_cpp(input$data)
      }

      .mc_finalize_corr_result(
        mat = bc_matrix,
        class_name = "bcdcor",
        method = "bias_corrected_distance_correlation",
        description = .mc_bcdcor_description(),
        output_cfg = output_cfg,
        diagnostics = .mc_merge_diagnostics(diagnostics, input$diagnostics),
        dimnames = input$dimnames,
        symmetric = TRUE,
        extra_attrs = if (!is.null(inference_attr)) list(inference = inference_attr)
      )
    }
  )
}

.mc_dcor_pairwise_summary <- function(object,
                                      digits = 4,
                                      p_digits = 4) {
  check_inherits(object, c("dcor", "bcdcor"))

  inf <- .mc_inference_attr(object)
  has_p <- is.list(inf) && is.matrix(inf$p_value)
  estimator_class <- .mc_corr_estimator_class(object)
  extra_columns <- list()
  if (is.list(inf) && is.matrix(inf$bcdcor) && !inherits(object, "bcdcor")) {
    extra_columns$bcdcor <- list(matrix = inf$bcdcor, digits = digits)
  }
  if (has_p) {
    extra_columns$statistic <- list(matrix = inf$statistic, digits = digits)
    extra_columns$df <- list(matrix = inf$parameter, digits = digits)
    extra_columns$p_value <- list(matrix = inf$p_value, digits = p_digits)
  }
  .mc_pairwise_matrix_summary(
    object,
    class_name = paste0("summary.", estimator_class),
    digits = digits,
    p_digits = p_digits,
    include_ci = FALSE,
    include_p = has_p,
    extra_columns = extra_columns,
    extra_attrs = list(
      inference_method = if (isTRUE(has_p)) inf$method %||% NA_character_ else NA_character_,
      squared = attr(object, "squared", exact = TRUE) %||% FALSE
    )
  )
}

#' @rdname dcor
#' @method print dcor
#' @title Print Method for \code{dcor} Objects
#'
#' @param x An object of class \code{dcor}.
#' @param digits Integer; number of decimal places to print.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param ... Additional arguments passed to \code{print}.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.dcor <- function(x, digits = 4, n = NULL, topn = NULL,
                       max_vars = NULL, width = NULL,
                       show_ci = NULL, ...) {
  .mc_print_corr_matrix(
    x,
    header = if (isTRUE(attr(x, "squared", exact = TRUE))) {
      "Squared distance correlation (dCor^2) matrix"
    } else {
      "Distance correlation (dCor) matrix"
    },
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
}

#' @rdname dcor
#' @method plot dcor
#' @title Plot Method for \code{dcor} Objects
#'
#' @param x An object of class \code{dcor}.
#' @param title Plot title. Default is \code{"Distance correlation heatmap"}.
#' @param low_color Colour for zero correlation. Default is \code{"white"}.
#' @param high_color Colour for strong correlation. Default is \code{"steelblue1"}.
#' @param value_text_size Font size for displaying values. Default is \code{4}.
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
#' @param ... Additional arguments passed to \code{ggplot2::theme()} or other
#' \code{ggplot2} layers.
#'
#' @return A \code{ggplot} object representing the heatmap.
#' @import ggplot2
#' @export
plot.dcor <-
  function(x, title = NULL,
           low_color = "white", high_color = "steelblue1",
           value_text_size = 4, show_value = TRUE, ...) {

    check_inherits(x, "dcor")
    check_bool(show_value, arg = "show_value")
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg ggplot2} is required for plotting.")
    }

    mat <- as.matrix(x)
    df <- as.data.frame(as.table(mat))
    colnames(df) <- c("Var1", "Var2", "value")

    df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))
    fill_label <- if (isTRUE(attr(x, "squared", exact = TRUE))) "dCor^2" else "dCor"
    title <- title %||% if (isTRUE(attr(x, "squared", exact = TRUE))) {
      "Squared distance correlation heatmap"
    } else {
      "Distance correlation heatmap"
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(Var2, Var1, fill = .data$value)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient(
        low = low_color, high = high_color,
        limits = c(0, 1), name = fill_label
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
        ggplot2::aes(label = sprintf("%.2f", .data$value)),
        size = value_text_size,
        color = "black"
      )
    }

    p
  }

#' @rdname dcor
#' @method summary dcor
#' @param object An object of class \code{dcor}.
#' @export
summary.dcor <- function(object, n = NULL, topn = NULL,
                         max_vars = NULL, width = NULL,
                         show_ci = NULL, ...) {
  check_inherits(object, "dcor")
  inf <- .mc_inference_attr(object)
  if (!is.list(inf) || !is.matrix(inf$p_value)) {
    return(.mc_summary_corr_matrix(object, topn = topn))
  }
  .mc_dcor_pairwise_summary(object)
}

#' @rdname dcor
#' @method print summary.dcor
#' @param x An object of class \code{summary.dcor}.
#' @export
print.summary.dcor <- function(x, digits = NULL, n = NULL,
                               topn = NULL, max_vars = NULL,
                               width = NULL, show_ci = NULL, ...) {
  .mc_print_pairwise_summary_digest(
    x,
    title = if (isTRUE(attr(x, "squared", exact = TRUE))) {
      "Squared distance correlation summary"
    } else {
      "Distance correlation summary"
    },
    digits = .mc_coalesce(digits, 4),
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    extra_items = c(inference = attr(x, "inference_method", exact = TRUE)),
    ...
  )
  invisible(x)
}

#' @rdname bcdcor
#' @method print bcdcor
#' @export
print.bcdcor <- function(x, digits = 4, n = NULL, topn = NULL,
                         max_vars = NULL, width = NULL,
                         show_ci = NULL, ...) {
  .mc_print_corr_matrix(
    x,
    header = "Bias-corrected squared distance correlation matrix",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    ...
  )
}

#' @rdname bcdcor
#' @method plot bcdcor
#' @export
plot.bcdcor <-
  function(x, title = "Bias-corrected squared distance correlation heatmap",
           low_color = "firebrick3", mid_color = "white", high_color = "steelblue1",
           value_text_size = 4, show_value = TRUE, ...) {

    check_inherits(x, "bcdcor")
    check_bool(show_value, arg = "show_value")
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg ggplot2} is required for plotting.")
    }

    mat <- as.matrix(x)
    df <- as.data.frame(as.table(mat))
    colnames(df) <- c("Var1", "Var2", "value")

    df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))

    p <- ggplot2::ggplot(df, ggplot2::aes(Var2, Var1, fill = .data$value)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient2(
        low = low_color, mid = mid_color, high = high_color,
        midpoint = 0, name = "bc-dCor^2"
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
        ggplot2::aes(label = sprintf("%.2f", .data$value)),
        size = value_text_size,
        color = "black"
      )
    }

    p
  }

#' @rdname bcdcor
#' @method summary bcdcor
#' @export
summary.bcdcor <- function(object, n = NULL, topn = NULL,
                           max_vars = NULL, width = NULL,
                           show_ci = NULL, ...) {
  check_inherits(object, "bcdcor")
  inf <- .mc_inference_attr(object)
  if (!is.list(inf) || !is.matrix(inf$p_value)) {
    return(.mc_summary_corr_matrix(object, topn = topn))
  }
  .mc_dcor_pairwise_summary(object)
}

#' @rdname bcdcor
#' @method print summary.bcdcor
#' @export
print.summary.bcdcor <- function(x, digits = NULL, n = NULL,
                                 topn = NULL, max_vars = NULL,
                                 width = NULL, show_ci = NULL, ...) {
  .mc_print_pairwise_summary_digest(
    x,
    title = "Bias-corrected squared distance correlation summary",
    digits = .mc_coalesce(digits, 4),
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    extra_items = c(inference = attr(x, "inference_method", exact = TRUE)),
    ...
  )
  invisible(x)
}
