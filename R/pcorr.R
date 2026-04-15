#' @title Partial correlation matrix (sample / ridge / OAS / graphical lasso)
#'
#' @description
#' Computes Gaussian partial correlations for the numeric columns of a matrix
#' or data frame using a high-performance 'C++' backend. Covariance estimation
#' is available via the classical sample estimator, ridge regularisation, OAS
#' shrinkage, or graphical lasso. Optional p-values and Fisher-z confidence
#' intervals are available for the classical sample estimator in the ordinary
#' low-dimensional setting.
#'
#' @param data A numeric matrix or data frame with at least two numeric columns.
#'   Non-numeric columns are ignored.
#' @param method Character; one of \code{"sample"}, \code{"oas"},
#'   \code{"ridge"}, or \code{"glasso"}. Default \code{"sample"}.
#' @param lambda Numeric \eqn{\ge 0}; regularisation strength. For
#'   \code{method = "ridge"}, this is the penalty added to the covariance
#'   diagonal. For \code{method = "glasso"}, this is the off-diagonal
#'   precision-matrix \eqn{\ell_1} penalty. Ignored otherwise. Default
#'   \code{1e-3}.
#' @param return_cov_precision Logical; if \code{TRUE}, also return the
#'   covariance (\code{cov}) and precision (\code{precision}) matrices used to
#'   form the partial correlations. Default to \code{FALSE}
#' @param return_p_value Logical; if \code{TRUE}, also return the matrix of
#'   two-sided p-values for testing whether each sample partial correlation is
#'   zero. This option is available only for \code{method = "sample"} and
#'   requires \eqn{n > p}. Default to \code{FALSE}.
#' @param ci Logical (default \code{FALSE}). If \code{TRUE}, attach Fisher-z
#'   confidence intervals for the off-diagonal partial correlations. This
#'   option is available only for the classical \code{method = "sample"}
#'   estimator in the ordinary low-dimensional setting.
#' @param conf_level Confidence level used when \code{ci = TRUE}. Default is
#'   \code{0.95}.
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
#'
#' @return An object of class \code{"partial_corr"} (a list) with elements:
#'   \itemize{
#'     \item \code{pcor}: \eqn{p \times p} partial correlation matrix.
#'     \item \code{cov} (if requested): covariance matrix used.
#'     \item \code{precision} (if requested): precision matrix \eqn{\Theta}.
#'     \item \code{p_value} (if requested): matrix of two-sided p-values for
#'     the sample partial correlations.
#'     \item \code{ci} (if requested): a list with elements \code{est},
#'     \code{lwr.ci}, \code{upr.ci}, \code{conf.level}, and \code{ci.method}.
#'     \item \code{diagnostics}: metadata used for inference, including the
#'     effective complete-case sample size and number of conditioned variables.
#'     \item \code{method}: the estimator used (\code{"oas"}, \code{"ridge"},
#'     \code{"sample"}, or \code{"glasso"}).
#'     \item \code{lambda}: ridge or graphical-lasso penalty
#'     (or \code{NA_real_}).
#'     \item \code{rho}: OAS shrinkage weight in \eqn{[0,1]} (or \code{NA_real_}).
#'     \item \code{jitter}: diagonal jitter added (if any) to ensure positive
#'     definiteness.
#'   }
#'
#' @details
#' \strong{Statistical overview.} Given an \eqn{n \times p} data matrix \eqn{X}
#' (rows are observations, columns are variables), the routine estimates a
#' \emph{partial correlation} matrix via the precision (inverse covariance)
#' matrix. Let \eqn{\mu} be the vector of column means and
#' \deqn{S = (X - \mathbf{1}\mu)^\top (X - \mathbf{1}\mu)}
#' be the centred cross-product matrix computed without forming a centred copy
#' of \eqn{X}. Two conventional covariance scalings are formed:
#' \deqn{\hat\Sigma_{\mathrm{MLE}} = S/n, \qquad
#'       \hat\Sigma_{\mathrm{unb}} = S/(n-1).}
#'
#' \itemize{
#'   \item \emph{Sample:} \eqn{\Sigma = \hat\Sigma_{\mathrm{unb}}}.
#'   \item \emph{Ridge:} \eqn{\Sigma = \hat\Sigma_{\mathrm{unb}} + \lambda I_p}
#'         with user-supplied \eqn{\lambda \ge 0} (diagonal inflation).
#'   \item \emph{OAS (Oracle Approximating Shrinkage):}
#'         shrink \eqn{\hat\Sigma_{\mathrm{MLE}}} towards a scaled identity
#'         target \eqn{\mu_I I_p}, where \eqn{\mu_I = \mathrm{tr}(\hat\Sigma_{\mathrm{MLE}})/p}.
#'         The data-driven weight \eqn{\rho \in [0,1]} is
#'         \deqn{\rho = \min\!\left\{1,\max\!\left(0,\;
#'         \frac{(1-\tfrac{2}{p})\,\mathrm{tr}(\hat\Sigma_{\mathrm{MLE}}^2)
#'               + \mathrm{tr}(\hat\Sigma_{\mathrm{MLE}})^2}
#'              {(n + 1 - \tfrac{2}{p})
#'               \left[\mathrm{tr}(\hat\Sigma_{\mathrm{MLE}}^2)
#'               - \tfrac{\mathrm{tr}(\hat\Sigma_{\mathrm{MLE}})^2}{p}\right]}
#'         \right)\right\},}
#'         and
#'         \deqn{\Sigma = (1-\rho)\,\hat\Sigma_{\mathrm{MLE}} + \rho\,\mu_I I_p.}
#'   \item \emph{Graphical lasso:} estimate a sparse precision matrix
#'         \eqn{\Theta} by maximising
#'         \deqn{\log\det(\Theta) - \mathrm{tr}(\hat\Sigma_{\mathrm{MLE}}\Theta)
#'         - \lambda\sum_{i \ne j} |\theta_{ij}|,}
#'         with \eqn{\lambda \ge 0}. The returned covariance matrix is
#'         \eqn{\Sigma = \Theta^{-1}}.
#' }
#'
#' The method then ensures positive definiteness of \eqn{\Sigma} (adding a very
#' small diagonal \emph{jitter} only if necessary) and computes the precision
#' matrix \eqn{\Theta = \Sigma^{-1}}. Partial correlations are obtained by
#' standardising the off-diagonals of \eqn{\Theta}:
#' \deqn{\mathrm{pcor}_{ij} \;=\;
#'       -\,\frac{\theta_{ij}}{\sqrt{\theta_{ii}\,\theta_{jj}}}, \qquad
#'       \mathrm{pcor}_{ii}=1.}
#'
#' If \code{return_p_value = TRUE}, the function also reports the classical
#' two-sided test p-values for the sample partial correlations, using
#' \deqn{t_{ij} = \mathrm{pcor}_{ij}
#'   \sqrt{\frac{n - p}{1 - \mathrm{pcor}_{ij}^2}}}
#' with \eqn{n - p} degrees of freedom. These p-values are returned only for
#' \code{method = "sample"}, where they match the standard full-model partial
#' correlation test.
#'
#' When \code{ci = TRUE}, the function reports Fisher-\eqn{z} confidence
#' intervals for the sample partial correlations. For a partial correlation
#' \eqn{r_{xy \cdot Z}} conditioning on \eqn{c} variables, the transformed
#' statistic is \eqn{z = \operatorname{atanh}(r_{xy \cdot Z})} with standard
#' error
#' \deqn{\operatorname{SE}(z) = \frac{1}{\sqrt{n - 3 - c}},}
#' where \eqn{n} is the effective complete-case sample size used for the
#' estimate. The two-sided normal-theory interval is formed on the transformed
#' scale using \code{conf_level} and then mapped back with \code{tanh()}. In
#' the full matrix path implemented here, each off-diagonal entry conditions on
#' all remaining variables, so \eqn{c = p - 2} and the classical CI requires
#' \eqn{n > p + 1}. This inference is only supported for
#' \code{method = "sample"} without positive-definiteness repair; in
#' unsupported or numerically singular settings, CI bounds are returned as
#' \code{NA} with an informative \pkg{cli} warning or the request is rejected.
#'
#' \strong{Interpretation.} For Gaussian data, \eqn{\mathrm{pcor}_{ij}} equals
#' the correlation between residuals from regressing variable \eqn{i} and
#' variable \eqn{j} on all the remaining variables; equivalently, it encodes
#' conditional dependence in a Gaussian graphical model, where
#' \eqn{\mathrm{pcor}_{ij}=0} if variables \eqn{i} and \eqn{j} are
#' conditionally independent given the others. Partial correlations are
#' invariant to separate rescalings of each
#' variable; in particular, multiplying \eqn{\Sigma} by any positive scalar
#' leaves the partial correlations unchanged.
#'
#' \strong{Why shrinkage/regularisation?} When \eqn{p \ge n}, the sample
#' covariance is singular and inversion is ill-posed. Ridge and OAS both yield
#' well-conditioned \eqn{\Sigma}. Ridge adds a fixed \eqn{\lambda} on the
#' diagonal, whereas OAS shrinks adaptively towards \eqn{\mu_I I_p} with a
#' weight chosen to minimise (approximately) the Frobenius risk under a
#' Gaussian model, often improving mean-square accuracy in high dimension.
#'
#' \strong{Why glasso?} Glasso is useful when the goal is not just to
#' stabilise a covariance estimate, but to recover a manageable network of
#' direct relationships rather than a dense matrix of overall associations. In
#' Gaussian models, zeros in the precision matrix correspond to conditional
#' independences, so glasso can suppress indirect associations that are
#' explained by the other variables and return a smaller, more interpretable
#' conditional-dependence graph. This is especially practical in
#' high-dimensional settings, where the sample covariance may be unstable or
#' singular. Glasso yields a positive-definite precision estimate and supports
#' edge selection, graph recovery, and downstream network analysis.
#'
#' \strong{Computational notes.} The implementation forms \eqn{S} using 'BLAS'
#' \code{syrk} when available and constructs partial correlations by traversing
#' only the upper triangle with 'OpenMP' parallelism. Positive definiteness is
#' verified via a Cholesky factorisation; if it fails, a tiny diagonal jitter is
#' increased geometrically up to a small cap, at which point the routine
#' signals an error.
#'
#' @examples
#' ## Structured MVN with known partial correlations
#' set.seed(42)
#' p <- 12; n <- 1000
#'
#' ## Build a tri-diagonal precision (Omega) so the true partial correlations
#' ## are sparse
#' phi <- 0.35
#' Omega <- diag(p)
#' for (j in 1:(p - 1)) {
#'   Omega[j, j + 1] <- Omega[j + 1, j] <- -phi
#' }
#' ## Strict diagonal dominance
#' diag(Omega) <- 1 + 2 * abs(phi) + 0.05
#' Sigma <- solve(Omega)
#'
#' ## Upper Cholesky
#' L <- chol(Sigma)
#' Z <- matrix(rnorm(n * p), n, p)
#' X <- Z %*% L
#' colnames(X) <- sprintf("V%02d", seq_len(p))
#'
#' pc <- pcorr(X)
#' summary(pc)
#'
#' # Interactive viewing (requires shiny)
#' if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
#'   view_corr_shiny(pc)
#' }
#'
#' ## True partial correlation from Omega
#' pcor_true <- -Omega / sqrt(diag(Omega) %o% diag(Omega))
#' diag(pcor_true) <- 1
#'
#' ## Quick visual check (first 5x5 block)
#' round(pc$pcor[1:5, 1:5], 2)
#' round(pcor_true[1:5, 1:5], 2)
#'
#' ## Plot method
#' plot(pc)
#'
#' ## Graphical-lasso example
#' set.seed(100)
#' p <- 20; n <- 250
#' Theta_g <- diag(p)
#' Theta_g[cbind(1:5, 2:6)] <- -0.25
#' Theta_g[cbind(2:6, 1:5)] <- -0.25
#' Theta_g[cbind(8:11, 9:12)] <- -0.20
#' Theta_g[cbind(9:12, 8:11)] <- -0.20
#' diag(Theta_g) <- rowSums(abs(Theta_g)) + 0.2
#'
#' Sigma_g <- solve(Theta_g)
#' X_g <- matrix(rnorm(n * p), n, p) %*% chol(Sigma_g)
#' colnames(X_g) <- paste0("Node", seq_len(p))
#'
#' gfit_1 <- pcorr(X_g, method = "glasso", lambda = 0.02,
#'                 return_cov_precision = TRUE)
#' gfit_2 <- pcorr(X_g, method = "glasso", lambda = 0.08,
#'                 return_cov_precision = TRUE)
#'
#' ## Larger lambda gives a sparser conditional-dependence graph
#' edge_count <- function(M, tol = 1e-8) {
#'   sum(abs(M[upper.tri(M, diag = FALSE)]) > tol)
#' }
#'
#' c(edges_lambda_002 = edge_count(gfit_1$precision),
#'   edges_lambda_008 = edge_count(gfit_2$precision))
#'
#' ## Inspect strongest estimated conditional associations
#' pcor_g <- gfit_1$pcor
#' idx <- which(upper.tri(pcor_g), arr.ind = TRUE)
#' ord <- order(abs(pcor_g[idx]), decreasing = TRUE)
#' head(data.frame(
#'   i = rownames(pcor_g)[idx[ord, 1]],
#'   j = colnames(pcor_g)[idx[ord, 2]],
#'   pcor = round(pcor_g[idx][ord], 2)
#' ))
#'
#' ## High-dimensional case p >> n
#' set.seed(7)
#' n <- 60; p <- 120
#'
#' ar_block <- function(m, rho = 0.6) rho^abs(outer(seq_len(m), seq_len(m), "-"))
#'
#' ## Two AR(1) blocks on the diagonal
#' if (requireNamespace("Matrix", quietly = TRUE)) {
#'   Sigma_hd <- as.matrix(Matrix::bdiag(ar_block(60, 0.6), ar_block(60, 0.6)))
#' } else {
#'   Sigma_hd <- rbind(
#'     cbind(ar_block(60, 0.6), matrix(0, 60, 60)),
#'     cbind(matrix(0, 60, 60), ar_block(60, 0.6))
#'   )
#' }
#'
#' L <- chol(Sigma_hd)
#' X_hd <- matrix(rnorm(n * p), n, p) %*% L
#' colnames(X_hd) <- paste0("G", seq_len(p))
#'
#' pc_oas   <-
#'  pcorr(X_hd, method = "oas",   return_cov_precision = TRUE)
#' pc_ridge <-
#'  pcorr(X_hd, method = "ridge", lambda = 1e-2,
#'                      return_cov_precision = TRUE)
#' pc_samp  <-
#'  pcorr(X_hd, method = "sample", return_cov_precision = TRUE)
#' pc_glasso <-
#'  pcorr(X_hd, method = "glasso", lambda = 5e-3,
#'                      return_cov_precision = TRUE)
#'
#' ## Show how much diagonal regularisation was used
#' c(oas_jitter = pc_oas$jitter,
#'   ridge_lambda = pc_ridge$lambda,
#'   sample_jitter = pc_samp$jitter,
#'   glasso_lambda = pc_glasso$lambda)
#'
#' ## Compare conditioning of the estimated covariance matrices
#' c(kappa_oas = kappa(pc_oas$cov),
#'   kappa_ridge = kappa(pc_ridge$cov),
#'   kappa_sample = kappa(pc_samp$cov))
#'
#' ## Simple conditional-dependence graph from partial correlations
#' pcor <- pc_oas$pcor
#' vals <- abs(pcor[upper.tri(pcor, diag = FALSE)])
#' thresh <- quantile(vals, 0.98)  # top 2%
#' edges  <- which(abs(pcor) > thresh & upper.tri(pcor), arr.ind = TRUE)
#' head(data.frame(i = colnames(pcor)[edges[,1]],
#'                 j = colnames(pcor)[edges[,2]],
#'                 pcor = round(pcor[edges], 2)))
#'
#' @references
#' Chen, Y., Wiesel, A., & Hero, A. O. III (2011).
#' Robust Shrinkage Estimation of High-dimensional Covariance Matrices.
#' IEEE Transactions on Signal Processing.
#'
#' @references
#' Friedman, J., Hastie, T., & Tibshirani, R. (2007).
#' Sparse inverse covariance estimation with the graphical lasso.
#' Biostatistics.
#'
#' @references
#' Ledoit, O., & Wolf, M. (2004).
#' A well-conditioned estimator for large-dimensional covariance matrices.
#' Journal of Multivariate Analysis, 88(2), 365-411.
#'
#' @references
#' Schafer, J., & Strimmer, K. (2005).
#' A shrinkage approach to large-scale covariance matrix estimation and
#' implications for functional genomics.
#' Statistical Applications in Genetics and Molecular Biology, 4(1), Article 32.
#'
#' @export
pcorr <- function(data, method = c("sample","oas","ridge","glasso"),
                                ci = FALSE, conf_level = 0.95,
                                return_cov_precision = FALSE,
                                return_p_value = FALSE, lambda = 1e-3,
                                output = c("matrix", "sparse", "edge_list"),
                                threshold = 0,
                                diag = TRUE) {
  output_cfg <- .mc_validate_output_args(
    output = output,
    threshold = threshold,
    diag = diag
  )
  method <- match.arg(method)
  lambda <- check_scalar_numeric(lambda, arg = "lambda", lower = 0, closed_lower = TRUE)
  lambda <- as.numeric(lambda)
  check_bool(return_cov_precision, arg = "return_cov_precision")
  check_bool(return_p_value, arg = "return_p_value")
  check_bool(ci, arg = "ci")
  if (isTRUE(ci)) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  }
  if (!identical(output_cfg$output, "matrix") &&
      (isTRUE(return_cov_precision) || isTRUE(return_p_value) || isTRUE(ci))) {
    abort_bad_arg(
      "output",
      message = "non-matrix outputs currently support point estimates only.",
      .hint = "Set `.arg return_cov_precision = FALSE`, `.arg return_p_value = FALSE`, and `.arg ci = FALSE`."
    )
  }

  numeric_data <-
    if (is.matrix(data) && is.double(data)) {
      data
    } else {
      # drops non-numeric
      validate_corr_input(data)
    }

  if (isTRUE(return_p_value) && !identical(method, "sample")) {
    cli::cli_abort(
      "{.arg return_p_value} is available only for {.code method = \"sample\"}."
    )
  }
  if (isTRUE(return_p_value) && nrow(numeric_data) <= ncol(numeric_data)) {
    cli::cli_abort(
      "{.arg return_p_value} requires {.code n > p} so that the sample partial-correlation test has positive degrees of freedom."
    )
  }
  if (isTRUE(ci) && !identical(method, "sample")) {
    cli::cli_abort(
      "{.arg ci} is available only for the classical {.code method = \"sample\"} partial correlation."
    )
  }
  if (isTRUE(ci) && nrow(numeric_data) <= (ncol(numeric_data) + 1L)) {
    cli::cli_abort(
      "{.arg ci} requires {.code n > p + 1} so that the Fisher-z partial-correlation interval has positive degrees of freedom."
    )
  }

  res <- partial_correlation_cpp(
    numeric_data,
    method,
    lambda,
    return_cov_precision,
    return_p_value
  )

  # set dimnames (cheap; attributes only)
  dn <- list(colnames(numeric_data), colnames(numeric_data))
  if (!is.null(res$pcor)) {
    res$pcor <- structure(res$pcor, dimnames = dn)
    if (!is.null(res$p_value))   res$p_value <- structure(res$p_value, dimnames = dn)
    if (!is.null(res$cov))       res$cov <- structure(res$cov, dimnames = dn)
    if (!is.null(res$precision)) res$precision <- structure(res$precision, dimnames = dn)
  } else {
    res <- list(pcor = structure(res[[1L]], dimnames = dn))
  }

  diagnostics <- list(
    n_complete = matrix(
      as.integer(nrow(numeric_data)),
      nrow = ncol(numeric_data),
      ncol = ncol(numeric_data),
      dimnames = dn
    ),
    n_conditioning = matrix(
      as.integer(ncol(numeric_data) - 2L),
      nrow = ncol(numeric_data),
      ncol = ncol(numeric_data),
      dimnames = dn
    )
  )
  diag(diagnostics$n_conditioning) <- 0L
  ci_attr <- NULL
  if (isTRUE(ci)) {
    ci_needs_repair <- .mc_sample_covariance_needs_ci_repair(numeric_data)
    ci_attr <- .mc_partial_corr_fisher_ci(
      pcor = res$pcor,
      n_complete = nrow(numeric_data),
      n_conditioning = ncol(numeric_data) - 2L,
      conf_level = conf_level,
      jitter = res$jitter %||% NA_real_,
      needs_repair = ci_needs_repair
    )
  }

  res$method <- method
  res$lambda <- if (method %in% c("ridge", "glasso")) lambda else NA_real_
  res$rho    <- if (identical(method, "oas"))   res$rho %||% NA_real_ else NA_real_
  res$jitter <- res$jitter %||% NA_real_
  res$diagnostics <- diagnostics
  if (!is.null(ci_attr)) {
    res$ci <- ci_attr
  }
  out <- structure(
    res,
    class = c("partial_corr", "list"),
    method = method,
    ci = ci_attr,
    conf.level = if (!is.null(ci_attr)) conf_level else NULL,
    ci.method = if (!is.null(ci_attr)) ci_attr$ci.method else NULL
  )
  if (identical(output_cfg$output, "matrix")) {
    return(out)
  }

  pcor_obj <- .mc_structure_corr_matrix(
    out$pcor,
    class_name = "partial_corr_matrix",
    method = paste0("partial_correlation_", method),
    description = "Partial correlation matrix",
    diagnostics = diagnostics,
    dimnames = dn,
    classes = c("partial_corr_matrix", "matrix")
  )
  .mc_finalize_corr_output(
    pcor_obj,
    output = output_cfg$output,
    threshold = output_cfg$threshold,
    diag = output_cfg$diag
  )
}


# small helper for older R versions without %||%
`%||%` <- function(a, b) if (!is.null(a)) a else b

.mc_sample_covariance_needs_ci_repair <- function(x) {
  centered <- scale(as.matrix(x), center = TRUE, scale = FALSE)
  qr(centered)$rank < ncol(centered)
}

.mc_partial_corr_ci_attr <- function(x) {
  if (inherits(x, "partial_corr") && !is.null(x$ci)) {
    return(x$ci)
  }
  attr(x, "ci", exact = TRUE)
}

.mc_partial_corr_pairwise_summary <- function(object,
                                              digits = 4,
                                              ci_digits = 3,
                                              show_ci = NULL) {
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )
  check_inherits(object, "partial_corr")

  est <- as.matrix(object$pcor)
  rn <- rownames(est)
  cn <- colnames(est)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(est)))
  if (is.null(cn)) cn <- as.character(seq_len(ncol(est)))

  ci <- object$ci %||% attr(object, "ci", exact = TRUE)
  diag_attr <- object$diagnostics %||% attr(object, "diagnostics", exact = TRUE)
  p_value <- object$p_value %||% NULL
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
      if (!is.null(p_value) && is.matrix(p_value) && identical(dim(p_value), dim(est))) {
        rec$p_value <- p_value[i, j]
      }
      if (include_ci) {
        rec$lwr <- if (!is.null(ci$lwr.ci) && is.finite(ci$lwr.ci[i, j])) {
          round(ci$lwr.ci[i, j], ci_digits)
        } else {
          NA_real_
        }
        rec$upr <- if (!is.null(ci$upr.ci) && is.finite(ci$upr.ci[i, j])) {
          round(ci$upr.ci[i, j], ci_digits)
        } else {
          NA_real_
        }
      }
      rows[[k]] <- rec
    }
  }

  df <- do.call(rbind.data.frame, rows)
  rownames(df) <- NULL
  if ("estimate" %in% names(df)) df$estimate <- as.numeric(df$estimate)
  if ("lwr" %in% names(df)) df$lwr <- as.numeric(df$lwr)
  if ("upr" %in% names(df)) df$upr <- as.numeric(df$upr)
  if ("p_value" %in% names(df)) df$p_value <- as.numeric(df$p_value)
  if ("n_complete" %in% names(df)) df$n_complete <- as.integer(df$n_complete)

  out <- .mc_finalize_summary_df(df, class_name = "summary.partial_corr")
  attr(out, "overview") <- .mc_summary_corr_matrix(object$pcor)
  attr(out, "has_ci") <- include_ci
  attr(out, "conf.level") <- if (is.null(ci)) NA_real_ else ci$conf.level
  attr(out, "digits") <- digits
  attr(out, "ci_digits") <- ci_digits
  out
}

.mc_partial_corr_fisher_ci <- function(pcor,
                                       n_complete,
                                       n_conditioning,
                                       conf_level = 0.95,
                                       jitter = NA_real_,
                                       needs_repair = FALSE) {
  pcor <- as.matrix(pcor)
  p <- ncol(pcor)
  dn <- dimnames(pcor)
  lwr <- matrix(NA_real_, p, p, dimnames = dn)
  upr <- matrix(NA_real_, p, p, dimnames = dn)
  diag(lwr) <- 1
  diag(upr) <- 1

  out <- list(
    est = unclass(pcor),
    lwr.ci = lwr,
    upr.ci = upr,
    conf.level = conf_level,
    ci.method = "fisher_z_partial"
  )

  if (isTRUE(needs_repair) || (!is.na(jitter) && is.finite(jitter) && jitter > 0)) {
    cli::cli_warn(
      c(
        "Partial-correlation confidence intervals are unavailable for this fit.",
        "i" = "The sample covariance was rank-deficient or required positive-definiteness repair, so the classical Fisher-z partial-correlation interval is not identified.",
        "i" = "Returning {.code NA} confidence bounds."
      ),
      class = c("matrixCorr_warning", "matrixCorr_ci_warning")
    )
    return(out)
  }

  se_denom <- as.numeric(n_complete - 3L - n_conditioning)
  if (!is.finite(se_denom) || se_denom <= 0) {
    cli::cli_abort(
      "{.arg ci} requires positive Fisher-z residual degrees of freedom; got {.code n_complete - 3 - c = {se_denom}}."
    )
  }

  crit <- stats::qnorm(0.5 * (1 + conf_level))
  se <- 1 / sqrt(se_denom)
  eps <- sqrt(.Machine$double.eps)
  boundary_pairs <- 0L

  for (j in seq_len(p - 1L)) {
    for (i in (j + 1L):p) {
      r <- pcor[j, i]
      if (!is.finite(r)) next
      if (abs(r) >= 1) {
        boundary_pairs <- boundary_pairs + 1L
        next
      }
      r_safe <- max(min(r, 1 - eps), -1 + eps)
      z <- atanh(r_safe)
      lo <- tanh(z - crit * se)
      hi <- tanh(z + crit * se)
      lwr[j, i] <- lwr[i, j] <- max(-1, lo)
      upr[j, i] <- upr[i, j] <- min(1, hi)
    }
  }

  out$lwr.ci <- lwr
  out$upr.ci <- upr
  if (boundary_pairs > 0L) {
    cli::cli_warn(
      "{boundary_pairs} partial-correlation pair{?s} were at the boundary +/-1; returning {.code NA} confidence bounds for those pair{?s}.",
      boundary_pairs = boundary_pairs,
      class = c("matrixCorr_warning", "matrixCorr_ci_warning")
    )
  }
  out
}

#' @rdname pcorr
#' @title Print method for \code{partial_corr}
#'
#' @param x An object of class \code{partial_corr}.
#' @param digits Integer; number of decimal places for display (default 3).
#' @param show_method Logical; print a one-line header with \code{method}
#'   (and \code{lambda}/\code{rho} if available). Default \code{TRUE}.
#' @param n Optional row threshold for compact preview output.
#' @param topn Optional number of leading/trailing rows to show when truncated.
#' @param max_vars Optional maximum number of visible columns; `NULL` derives this
#'   from console width.
#' @param width Optional display width; defaults to \code{getOption("width")}.
#' @param show_ci One of \code{"yes"} or \code{"no"}.
#' @param ... Further arguments passed to \code{print.matrix()}.
#' @return Invisibly returns \code{x}.
#' @method print partial_corr
#' @importFrom utils capture.output
#' @export
print.partial_corr <- function(
    x,
    digits = 3,
    show_method = TRUE,
    n = NULL,
    topn = NULL,
    max_vars = NULL,
    width = NULL,
    show_ci = NULL,
    ...
) {
  check_inherits(x, "partial_corr")
  M <- x$pcor
  check_matrix_dims(M, arg = "x$pcor")
  M <- as.matrix(M)

  lines <- character()
  if (isTRUE(show_method)) {
    meth <- if (!is.null(x$method)) as.character(x$method) else NA_character_
    hdr <- switch(
      tolower(meth),
      "oas"   = {
        rho <- if (!is.null(x$rho) && is.finite(x$rho)) sprintf(", OAS rho=%.3f", x$rho) else ""
        paste0("Partial correlation (OAS", rho, ")")
      },
      "ridge" = {
        lam <- if (!is.null(x$lambda) && is.finite(x$lambda)) sprintf(", lambda=%.3g", x$lambda) else ""
        paste0("Partial correlation (ridge", lam, ")")
      },
      "glasso" = {
        lam <- if (!is.null(x$lambda) && is.finite(x$lambda)) sprintf(", lambda=%.3g", x$lambda) else ""
        paste0("Partial correlation (glasso", lam, ")")
      },
      "sample" = "Partial correlation (sample covariance)",
      "Partial correlation"
    )
    lines <- c(lines, hdr)
  } else {
    lines <- c(lines, "Partial correlation matrix:")
  }

  if (length(lines)) writeLines(lines)
  cat("\n")
  .mc_print_corr_matrix(
    x,
    header = "Partial correlation matrix",
    digits = digits,
    n = n,
    topn = topn,
    max_vars = max_vars,
    width = width,
    show_ci = show_ci,
    mat = M,
    ...
  )

  invisible(x)
}

#' @rdname pcorr
#' @method summary partial_corr
#' @title Summary Method for \code{partial_corr} Objects
#'
#' @param object An object of class \code{partial_corr}.
#' @param ... Unused.
#'
#' @return A compact summary object of class \code{summary.partial_corr}.
#' @export
summary.partial_corr <- function(object, n = NULL, topn = NULL,
                                 max_vars = NULL, width = NULL,
                                 show_ci = NULL, ...) {
  check_inherits(object, "partial_corr")
  show_ci <- .mc_validate_yes_no(
    show_ci,
    arg = "show_ci",
    default = .mc_display_option("summary_show_ci", "yes")
  )
  if (!is.null(.mc_partial_corr_ci_attr(object))) {
    return(.mc_partial_corr_pairwise_summary(
      object,
      show_ci = show_ci
    ))
  }
  out <- .mc_summary_corr_matrix(object$pcor, topn = topn)
  out$class <- "partial_corr"
  out$method <- object$method %||% attr(object, "method")
  out$lambda <- object$lambda %||% NA_real_
  out$rho <- object$rho %||% NA_real_
  out$jitter <- object$jitter %||% NA_real_
  out$header <- "Correlation summary"
  class(out) <- c("summary.partial_corr", "summary.matrixCorr", "summary.corr_matrix")
  out
}

#' @rdname pcorr
#' @method print summary.partial_corr
#' @export
print.summary.partial_corr <- function(x, digits = 4, n = NULL, topn = NULL,
                                       max_vars = NULL, width = NULL,
                                       show_ci = NULL, ...) {
  if (inherits(x, "data.frame")) {
    .mc_print_pairwise_summary_digest(
      x,
      title = "Partial correlation summary",
      digits = .mc_coalesce(digits, 4),
      n = n,
      topn = topn,
      max_vars = max_vars,
      width = width,
      show_ci = show_ci,
      ci_method = "fisher_z_partial",
      ...
    )
    return(invisible(x))
  }
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

print.summary_partial_corr <- print.summary.partial_corr


#' @rdname pcorr
#' @method plot partial_corr
#' @title Plot Method for \code{partial_corr} Objects
#'
#' @param x An object of class \code{partial_corr}.
#' @param title Plot title. By default, constructed from the estimator in
#' \code{x$method}.
#' @param low_color Colour for low (negative) values. Default
#' \code{"indianred1"}.
#' @param high_color Colour for high (positive) values. Default
#' \code{"steelblue1"}.
#' @param mid_color Colour for zero. Default \code{"white"}.
#' @param value_text_size Font size for cell labels. Default \code{4}.
#' @param show_value Logical; if \code{TRUE} (default), overlay numeric values
#'   on the heatmap tiles.
#' @param mask_diag Logical; if \code{TRUE}, the diagonal is masked
#' (set to \code{NA}) and not labelled. Default \code{TRUE}.
#' @param reorder Logical; if \code{TRUE}, variables are reordered by
#' hierarchical clustering of \eqn{1 - |pcor|}. Default \code{FALSE}.
#' @param ... Additional arguments passed to \code{ggplot2::theme()} or other
#'   \pkg{ggplot2} layers.
#'
#' @return A \code{ggplot} object.
#' @import ggplot2
#' @importFrom stats as.dist hclust
#' @export
plot.partial_corr <- function(
    x,
    title = NULL,
    low_color  = "indianred1",
    high_color = "steelblue1",
    mid_color  = "white",
    value_text_size = 4,
    show_value = TRUE,
    mask_diag = TRUE,
    reorder   = FALSE,
    ...
) {
  check_inherits(x, "partial_corr")
  check_bool(show_value, arg = "show_value")
  check_bool(mask_diag, arg = "mask_diag")
  check_bool(reorder,   arg = "reorder")
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ggplot2} is required for plotting.")
  }

  M <- x$pcor
  check_matrix_dims(M, arg = "x$pcor")

  # Ensure dimnames for labelling
  if (is.null(colnames(M))) colnames(M) <- paste0("V", seq_len(ncol(M)))
  if (is.null(rownames(M))) rownames(M) <- colnames(M)

  # Optional reordering by hierarchical clustering of 1 - |pcor|
  if (isTRUE(reorder) && nrow(M) >= 2L) {
    if (nrow(M) == ncol(M)) {
      dist_mat <- 1 - pmin(1, abs(M))
      dist_mat <- (dist_mat + t(dist_mat)) / 2
      diag(dist_mat) <- 0
      dist_mat[!is.finite(dist_mat)] <- 0
      dist_mat <- as.matrix(dist_mat)
      if (nrow(dist_mat) == ncol(dist_mat)) {
        D <- tryCatch(
          stats::as.dist(dist_mat),
          warning = function(w) {
            if (grepl("non-square matrix", conditionMessage(w), fixed = TRUE)) {
              return(NULL)
            }
            warning(w)
            NULL
          }
        )
        if (!is.null(D) && isTRUE(attr(D, "Size") >= 2)) {
          hc <- stats::hclust(D, method = "average")
          ord <- hc$order
          M <- M[ord, ord, drop = FALSE]
        }
      }
    }
  }

  # Default title constructed from x$method (+ tuning, if present)
  if (is.null(title)) {
    method <- tolower(as.character(x$method %||% ""))
    extra <- switch(
      method,
      "oas"   = if (is.finite(x$rho %||% NA_real_)) sprintf(" (OAS, rho=%.3f)", x$rho) else " (OAS)",
      "ridge" = if (is.finite(x$lambda %||% NA_real_)) sprintf(" (ridge, lambda=%.3g)", x$lambda) else " (ridge)",
      "glasso" = if (is.finite(x$lambda %||% NA_real_)) sprintf(" (glasso, lambda=%.3g)", x$lambda) else " (glasso)",
      "sample" = " (sample)",
      ""
    )
    title <- paste0("Partial correlation heatmap", extra)
  }

  # Build plotting frame
  df <- as.data.frame(as.table(M))
  names(df) <- c("Var1", "Var2", "PCor")

  # Mask diagonal if requested
  if (isTRUE(mask_diag)) {
    df$PCor[df$Var1 == df$Var2] <- NA_real_
  }

  # Reverse y-axis order for a tidy heatmap
  df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Var2, y = Var1, fill = PCor)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::scale_fill_gradient2(
      low = low_color, high = high_color, mid = mid_color,
      midpoint = 0, limits = c(-1, 1), name = "Partial r",
      na.value = "grey95"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid  = ggplot2::element_blank(),
      ...
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title, x = NULL, y = NULL)

  if (isTRUE(show_value) && !is.null(value_text_size) && is.finite(value_text_size)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.na(PCor), "", sprintf("%.2f", PCor))),
      size = value_text_size, colour = "black"
    )
  }

  p
}

