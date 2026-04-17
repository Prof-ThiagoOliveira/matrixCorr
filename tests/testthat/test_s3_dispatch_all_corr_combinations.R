case_output_args <- function(output, threshold_non_matrix = 0.2) {
  if (identical(output, "matrix")) {
    list(output = "matrix", threshold = 0, diag = TRUE)
  } else {
    list(output = output, threshold = threshold_non_matrix, diag = FALSE)
  }
}

assert_s3_dispatch <- function(obj,
                               case_label,
                               expected_summary_class,
                               expected_output = NULL,
                               expected_has_ci = NULL,
                               expected_has_p = NULL) {
  sm <- summary(obj)
  expect_s3_class(sm, expected_summary_class)

  if (inherits(sm, "summary.corr_result")) {
    expect_false(inherits(sm, "summaryDefault"), info = case_label)
    if (!is.null(expected_output)) {
      expect_identical(attr(sm, "output", exact = TRUE), expected_output, info = case_label)
    }
    if (!is.null(expected_has_ci)) {
      expect_identical(isTRUE(attr(sm, "has_ci", exact = TRUE)), expected_has_ci, info = case_label)
    }
    if (!is.null(expected_has_p)) {
      expect_identical(isTRUE(attr(sm, "has_p", exact = TRUE)), expected_has_p, info = case_label)
    }
  }

  gp <- plot(obj, show_value = FALSE)
  expect_s3_class(gp, "ggplot")
}

test_that("S3 dispatch works across all core correlation estimators and output combinations", {
  skip_if_not_installed("ggplot2")

  set.seed(20260417)
  X <- matrix(rnorm(180 * 5), nrow = 180, ncol = 5)
  colnames(X) <- paste0("V", seq_len(ncol(X)))

  Z <- matrix(rnorm(320 * 4), nrow = 320, ncol = 4)
  X_bin <- data.frame(
    b1 = Z[, 1] > -0.2,
    b2 = Z[, 2] > 0.1,
    b3 = Z[, 3] > 0.4,
    b4 = Z[, 4] > -0.5
  )
  X_ord <- data.frame(
    o1 = ordered(cut(Z[, 1], breaks = c(-Inf, -0.4, 0.3, Inf), labels = c("L", "M", "H"))),
    o2 = ordered(cut(Z[, 2], breaks = c(-Inf, -0.8, -0.1, 0.6, Inf), labels = c("1", "2", "3", "4"))),
    o3 = ordered(cut(Z[, 3], breaks = c(-Inf, -0.2, 0.5, Inf), labels = c("A", "B", "C"))),
    o4 = ordered(cut(Z[, 4], breaks = c(-Inf, -0.3, 0.2, Inf), labels = c("Q", "R", "S")))
  )

  outputs <- c("matrix", "sparse", "edge_list")

  ci_methods <- list(
    pearson_corr = function(output, ci) {
      do.call(
        pearson_corr,
        c(list(data = X, ci = ci), case_output_args(output))
      )
    },
    spearman_rho = function(output, ci) {
      do.call(
        spearman_rho,
        c(list(data = X, ci = ci), case_output_args(output))
      )
    },
    kendall_tau = function(output, ci) {
      do.call(
        kendall_tau,
        c(list(data = X, ci = ci), case_output_args(output))
      )
    },
    bicor = function(output, ci) {
      do.call(
        bicor,
        c(list(data = X, ci = ci), case_output_args(output))
      )
    },
    ccc = function(output, ci) {
      do.call(
        ccc,
        c(list(data = X, ci = ci), case_output_args(output))
      )
    },
    tetrachoric = function(output, ci) {
      do.call(
        tetrachoric,
        c(list(data = X_bin, ci = ci, p_value = ci), case_output_args(output, threshold_non_matrix = 0.1))
      )
    },
    polychoric = function(output, ci) {
      do.call(
        polychoric,
        c(list(data = X_ord, ci = ci, p_value = ci), case_output_args(output, threshold_non_matrix = 0.1))
      )
    }
  )

  for (nm in names(ci_methods)) {
    fn <- ci_methods[[nm]]
    for (output in outputs) {
      for (ci in c(FALSE, TRUE)) {
        case_label <- sprintf("%s output=%s ci=%s", nm, output, ci)
        obj <- fn(output, ci)
        expected_class <- if (identical(nm, "ccc") && identical(output, "matrix")) "summary.ccc" else "summary.corr_result"
        assert_s3_dispatch(
          obj = obj,
          case_label = case_label,
          expected_summary_class = expected_class,
          expected_output = if (identical(expected_class, "summary.corr_result")) output else NULL,
          expected_has_ci = if (identical(expected_class, "summary.corr_result")) ci else NULL
        )
      }
    }
  }

  for (output in outputs) {
    out_args <- case_output_args(output)
    for (p_value in c(FALSE, TRUE)) {
      case_label <- sprintf("dcor output=%s p_value=%s", output, p_value)
      obj <- do.call(
        dcor,
        c(list(data = X, p_value = p_value), out_args)
      )
      assert_s3_dispatch(
        obj = obj,
        case_label = case_label,
        expected_summary_class = "summary.corr_result",
        expected_output = output,
        expected_has_ci = FALSE,
        expected_has_p = p_value
      )
    }
  }

  robust_methods <- list(
    pbcor = function(output, infer) {
      do.call(
        pbcor,
        c(
          list(data = X, ci = infer, p_value = infer, n_boot = 20, seed = 42),
          case_output_args(output)
        )
      )
    },
    wincor = function(output, infer) {
      do.call(
        wincor,
        c(
          list(data = X, ci = infer, p_value = infer, n_boot = 20, seed = 42),
          case_output_args(output)
        )
      )
    },
    skipped_corr = function(output, infer) {
      do.call(
        skipped_corr,
        c(
          list(
            data = X,
            method = "pearson",
            ci = infer,
            p_value = infer,
            n_boot = 30,
            seed = 42,
            return_masks = FALSE
          ),
          case_output_args(output)
        )
      )
    }
  )

  for (nm in names(robust_methods)) {
    fn <- robust_methods[[nm]]
    for (output in outputs) {
      for (infer in c(FALSE, TRUE)) {
        case_label <- sprintf("%s output=%s infer=%s", nm, output, infer)
        obj <- fn(output, infer)
        assert_s3_dispatch(
          obj = obj,
          case_label = case_label,
          expected_summary_class = "summary.corr_result",
          expected_output = output,
          expected_has_ci = infer,
          expected_has_p = if (identical(nm, "skipped_corr")) NULL else infer
        )
      }
    }
  }

  aliases <- list(
    shrinkage_corr = function(output) {
      do.call(shrinkage_corr, c(list(data = X), case_output_args(output)))
    },
    schafer_corr = function(output) {
      do.call(schafer_corr, c(list(data = X), case_output_args(output)))
    }
  )

  for (nm in names(aliases)) {
    fn <- aliases[[nm]]
    for (output in outputs) {
      case_label <- sprintf("%s output=%s", nm, output)
      obj <- fn(output)
      assert_s3_dispatch(
        obj = obj,
        case_label = case_label,
        expected_summary_class = "summary.corr_result",
        expected_output = output,
        expected_has_ci = FALSE
      )
    }
  }
})

test_that("S3 dispatch works for partial correlation matrix and non-matrix outputs", {
  skip_if_not_installed("ggplot2")

  set.seed(20260418)
  X <- matrix(rnorm(220 * 5), nrow = 220, ncol = 5)
  colnames(X) <- paste0("P", seq_len(ncol(X)))

  matrix_fit <- pcorr(
    X,
    method = "sample",
    ci = TRUE,
    return_p_value = TRUE,
    output = "matrix"
  )
  sm_matrix <- summary(matrix_fit)
  expect_s3_class(sm_matrix, "summary.partial_corr")
  expect_s3_class(plot(matrix_fit, show_value = FALSE), "ggplot")

  sparse_fit <- pcorr(
    X,
    method = "sample",
    ci = FALSE,
    return_p_value = FALSE,
    output = "sparse",
    threshold = 0.2,
    diag = FALSE
  )
  assert_s3_dispatch(
    obj = sparse_fit,
    case_label = "pcorr sample sparse",
    expected_summary_class = "summary.corr_result",
    expected_output = "sparse",
    expected_has_ci = FALSE
  )

  edge_fit <- pcorr(
    X,
    method = "ridge",
    ci = FALSE,
    return_p_value = FALSE,
    output = "edge_list",
    threshold = 0.2,
    diag = FALSE
  )
  assert_s3_dispatch(
    obj = edge_fit,
    case_label = "pcorr ridge edge_list",
    expected_summary_class = "summary.corr_result",
    expected_output = "edge_list",
    expected_has_ci = FALSE
  )
})
