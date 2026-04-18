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

  sparse <- pearson_corr(X, output = "sparse", threshold = 0.30, diag = FALSE, ci = TRUE)
  expect_s4_class(sparse, "sparseMatrix")
  expect_identical(attr(sparse, "output", exact = TRUE), "sparse")
  expect_identical(attr(sparse, "method", exact = TRUE), "pearson")
  expect_true(is.list(attr(sparse, "ci", exact = TRUE)))

  sm_sparse <- summary(sparse)
  expect_s3_class(sm_sparse, "summary.corr_result")
  expect_s3_class(sm_sparse, "summary.corr_sparse")
  expect_identical(attr(sm_sparse, "output", exact = TRUE), "sparse")
  expect_true(isTRUE(attr(sm_sparse, "has_ci", exact = TRUE)))
  expect_false(inherits(sm_sparse, "summaryDefault"))
  expect_s3_class(plot(sparse, show_value = FALSE), "ggplot")

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
    sm <- summary(obj)
    ov <- attr(sm, "overview", exact = TRUE)
    expect_true(is.list(ov))
    expect_identical(attr(sm, "method", exact = TRUE), "pearson")
    expect_true(attr(sm, "output", exact = TRUE) %in% c("matrix", "sparse", "edge_list"))
    expect_true(is.logical(attr(sm, "diag", exact = TRUE)))
    expect_true(is.numeric(attr(sm, "threshold", exact = TRUE)))
  }
})

test_that("summary() dispatches CI-aware sparse outputs across correlation estimators", {
  set.seed(20260417)
  X <- matrix(rnorm(260 * 5), nrow = 260, ncol = 5)
  colnames(X) <- paste0("S", seq_len(ncol(X)))

  methods <- list(
    pearson_corr = function(...) pearson_corr(X, ci = TRUE, ...),
    spearman_rho = function(...) spearman_rho(X, ci = TRUE, ...),
    kendall_tau = function(...) kendall_tau(X, ci = TRUE, ...)
  )

  for (nm in names(methods)) {
    fn <- methods[[nm]]

    dense <- fn(output = "matrix")
    sm_dense <- summary(dense)
    expect_s3_class(sm_dense, "summary.corr_result")
    expect_s3_class(sm_dense, "summary.corr_matrix")
    expect_true(isTRUE(attr(sm_dense, "has_ci", exact = TRUE)), info = paste(nm, "dense"))

    sparse <- fn(output = "sparse", threshold = 0.20, diag = FALSE)
    sm_sparse <- summary(sparse)
    expect_s3_class(sm_sparse, "summary.corr_result")
    expect_s3_class(sm_sparse, "summary.corr_sparse")
    expect_true(isTRUE(attr(sm_sparse, "has_ci", exact = TRUE)), info = paste(nm, "sparse"))
    expect_false(inherits(sm_sparse, "summaryDefault"), info = paste(nm, "sparse"))

    edge <- fn(output = "edge_list", threshold = 0.20, diag = FALSE)
    sm_edge <- summary(edge)
    expect_s3_class(sm_edge, "summary.corr_result")
    expect_s3_class(sm_edge, "summary.corr_edge_list")
    expect_true(isTRUE(attr(sm_edge, "has_ci", exact = TRUE)), info = paste(nm, "edge_list"))
  }
})

test_that("corr_result heatmaps show CI labels when ci=TRUE for all output modes", {
  skip_if_not_installed("ggplot2")

  set.seed(20260418)
  x1 <- rnorm(260)
  x2 <- x1 + rnorm(260, sd = 0.15)
  x3 <- -x1 + rnorm(260, sd = 0.15)
  x4 <- rnorm(260)
  X <- cbind(x1 = x1, x2 = x2, x3 = x3, x4 = x4)

  count_text_layers <- function(p) {
    sum(vapply(p$layers, function(layer) inherits(layer$geom, "GeomText"), logical(1)))
  }

  matrix_fit <- pearson_corr(X, ci = TRUE)
  sparse_fit <- pearson_corr(X, output = "sparse", threshold = 0.6, diag = FALSE, ci = TRUE)
  edge_fit <- pearson_corr(X, output = "edge_list", threshold = 0.6, diag = FALSE, ci = TRUE)

  gp_matrix <- plot(matrix_fit, show_value = TRUE)
  gp_sparse <- plot(sparse_fit, show_value = TRUE)
  gp_edge <- plot(edge_fit, show_value = TRUE)

  expect_gte(count_text_layers(gp_matrix), 2L)
  expect_gte(count_text_layers(gp_sparse), 2L)
  expect_gte(count_text_layers(gp_edge), 2L)
})
