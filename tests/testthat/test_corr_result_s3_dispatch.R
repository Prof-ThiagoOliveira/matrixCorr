test_that("correlation outputs dispatch summary/plot by representation class", {
  skip_if_not_installed("ggplot2")

  set.seed(20260415)
  X <- matrix(rnorm(300 * 6), nrow = 300, ncol = 6)
  colnames(X) <- paste0("V", seq_len(ncol(X)))

  dense <- pearson_corr(X, ci = TRUE)
  expect_s3_class(dense, "corr_matrix")
  expect_s3_class(dense, "pearson_corr")
  expect_identical(attr(dense, "method", exact = TRUE), "pearson")
  expect_true(is.list(attr(dense, "ci", exact = TRUE)))

  sm_dense <- summary(dense)
  expect_s3_class(sm_dense, "summary.corr_result")
  expect_s3_class(sm_dense, "summary.corr_matrix")
  expect_identical(attr(sm_dense, "method", exact = TRUE), "pearson")
  expect_true(isTRUE(attr(sm_dense, "has_ci", exact = TRUE)))
  expect_false(inherits(sm_dense, "summaryDefault"))
  expect_s3_class(plot(dense, show_value = FALSE), "ggplot")

  sparse <- pearson_corr(X, output = "sparse", threshold = 0.30, diag = FALSE)
  expect_s4_class(sparse, "sparseMatrix")
  expect_identical(attr(sparse, "output", exact = TRUE), "sparse")
  expect_identical(attr(sparse, "method", exact = TRUE), "pearson")

  sm_sparse <- summary.corr_sparse(sparse)
  expect_s3_class(sm_sparse, "summary.corr_result")
  expect_s3_class(sm_sparse, "summary.corr_sparse")
  expect_identical(attr(sm_sparse, "output", exact = TRUE), "sparse")
  expect_false(inherits(sm_sparse, "summaryDefault"))
  expect_s3_class(plot.corr_sparse(sparse, show_value = FALSE), "ggplot")

  expect_error(pearson_corr(X, output = "packed_upper", threshold = 0.30, diag = FALSE))

  edge <- pearson_corr(X, output = "edge_list", threshold = 0.30, diag = FALSE)
  expect_s3_class(edge, "corr_edge_list")
  expect_identical(attr(edge, "output", exact = TRUE), "edge_list")
  expect_identical(attr(edge, "method", exact = TRUE), "pearson")

  sm_edge <- summary(edge)
  expect_s3_class(sm_edge, "summary.corr_result")
  expect_s3_class(sm_edge, "summary.corr_edge_list")
  expect_identical(attr(sm_edge, "output", exact = TRUE), "edge_list")
  expect_false(inherits(sm_edge, "summaryDefault"))
  expect_s3_class(plot(edge, show_value = FALSE), "ggplot")
  edge_all <- pearson_corr(X, output = "edge_list", threshold = 0, diag = TRUE)
  edge_print <- capture.output(print(edge_all, n = 8, topn = 2))
  expect_true(any(grepl("V\\d+\\s+V\\d+", edge_print)))
  expect_false(any(grepl("^X(\\.|\\s)", edge_print)))
})

test_that("correlation summary objects keep standardized overview metadata", {
  set.seed(42)
  X <- matrix(rnorm(160 * 4), nrow = 160, ncol = 4)
  colnames(X) <- paste0("M", seq_len(ncol(X)))

  out <- list(
    pearson_corr(X),
    pearson_corr(X, output = "sparse", threshold = 0.25, diag = FALSE),
    pearson_corr(X, output = "edge_list", threshold = 0.25, diag = FALSE)
  )

  for (obj in out) {
    sm <- if (inherits(obj, "sparseMatrix")) summary.corr_sparse(obj) else summary(obj)
    ov <- attr(sm, "overview", exact = TRUE)
    expect_true(is.list(ov))
    expect_identical(attr(sm, "method", exact = TRUE), "pearson")
    expect_true(attr(sm, "output", exact = TRUE) %in% c("matrix", "sparse", "edge_list"))
    expect_true(is.logical(attr(sm, "diag", exact = TRUE)))
    expect_true(is.numeric(attr(sm, "threshold", exact = TRUE)))
  }
})
