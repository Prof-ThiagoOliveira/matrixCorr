sort_unnamed_pairs <- function(td) {
  td[order(td$item1, td$item2), , drop = FALSE]
}

test_that("unnamed edge accessors use source coordinates", {
  set.seed(10)
  X <- matrix(rnorm(40), 10, 4)
  dense <- pearson_corr(X, output = "matrix", ci = TRUE, diag = FALSE)
  sparse <- pearson_corr(X, output = "sparse", ci = TRUE, diag = FALSE)
  edge <- pearson_corr(X, output = "edge_list", ci = TRUE, diag = FALSE)
  raw_edge <- edge
  expect_null(attr(edge, "corr_dimnames", exact = TRUE))
  expect_identical(names(edge), c("row", "col", "value"))
  expect_true(all(is.finite(edge$value)))
  expect_identical(matrixCorr:::.mc_corr_dimnames(edge), list(NULL, NULL))

  td_dense <- sort_unnamed_pairs(tidy(dense))
  td_sparse <- sort_unnamed_pairs(tidy(sparse))
  td_edge <- sort_unnamed_pairs(tidy(edge))
  expect_identical(td_edge$item1, td_dense$item1)
  expect_identical(td_edge$item2, td_dense$item2)
  expect_identical(td_edge$estimate, td_dense$estimate)
  expect_identical(td_edge$lwr, td_dense$lwr)
  expect_identical(td_edge$upr, td_dense$upr)
  expect_false(any(td_edge$item1 == td_edge$item2))
  expect_identical(td_sparse$estimate, td_dense$estimate)
  cols <- intersect(names(td_edge), names(td_dense))
  expect_identical(td_edge[, cols], td_dense[, cols])

  sm_edge <- summary(edge)
  sm_dense <- summary(dense)
  sm_cols <- intersect(c("item1", "item2", "estimate", "lwr", "upr", "n_complete"), names(sm_edge))
  expect_true(length(sm_cols) >= 3L)
  expect_identical(
    as.data.frame(sort_unnamed_pairs(sm_edge))[, sm_cols],
    as.data.frame(sort_unnamed_pairs(sm_dense))[, sm_cols]
  )
  if ("p_value" %in% names(sm_edge)) {
    expect_identical(sort_unnamed_pairs(sm_edge)$p_value, sort_unnamed_pairs(sm_dense)$p_value)
  }
  expect_identical(confint(edge), confint(dense))
  expect_identical(ci(edge)$conf.level, 0.95)
  expect_identical(ci(dense)$conf.level, ci(edge)$conf.level)
  expect_identical(edge, raw_edge)

  named <- X
  dimnames(named) <- list(NULL, as.character(seq_len(ncol(named))))
  named_edge <- pearson_corr(named, output = "edge_list", ci = TRUE, diag = FALSE)
  named_dense <- pearson_corr(named, output = "matrix", ci = TRUE, diag = FALSE)
  expect_identical(
    matrixCorr:::.mc_corr_dimnames(named_edge),
    matrixCorr:::.mc_corr_dimnames(named_dense)
  )
  plot_cols <- c("row", "col", "value", "lwr", "upr", "ci_label")
  edge_plot <- plot(edge)$data
  named_plot <- plot(named_edge)$data
  expect_true(all(plot_cols %in% names(edge_plot)))
  expect_identical(
    as.data.frame(edge_plot)[, plot_cols],
    as.data.frame(named_plot)[, plot_cols]
  )
})

test_that("empty and rectangular unnamed edges keep source positions", {
  set.seed(11)
  X <- matrix(rnorm(30), 10, 3)
  empty <- pearson_corr(X, output = "edge_list", threshold = 2, diag = FALSE)
  expect_null(attr(empty, "corr_dimnames", exact = TRUE))
  expect_identical(nrow(tidy(empty)), 0L)
  expect_identical(nrow(summary(empty)), 0L)

  legacy <- data.frame(i = c(1L, 2L), j = c(1L, 3L), x = c(0.2, -0.4))
  class(legacy) <- c("corr_edge_list", "pearson_corr", "corr_result", "data.frame")
  attr(legacy, "output") <- "edge_list"
  attr(legacy, "corr_dim") <- c(2L, 3L)
  attr(legacy, "corr_dimnames") <- NULL
  attr(legacy, "corr_symmetric") <- FALSE
  attr(legacy, "corr_estimator_class") <- "pearson_corr"
  attr(legacy, "diag") <- TRUE
  td <- tidy(legacy, triangle = "full", diag = TRUE)
  expect_identical(td$item1, c("1", "2"))
  expect_identical(td$item2, c("1", "3"))
  expect_identical(td$estimate, c(0.2, -0.4))
  expect_null(attr(legacy, "corr_dimnames", exact = TRUE))
})
