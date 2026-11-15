overview_threshold <- function(object) {
  sm <- summary(object, show_ci = "no")
  attr(sm, "overview", exact = TRUE)
}

test_that("public summaries keep filtering thresholds when latent threshold sets exist", {
  set.seed(5)
  dense_x <- matrix(rnorm(40), nrow = 10, ncol = 4)
  colnames(dense_x) <- paste0("y", seq_len(4))
  binary <- matrix(sample(0:1, 60, replace = TRUE), nrow = 20, ncol = 3)
  colnames(binary) <- paste0("b", seq_len(3))

  dense <- pearson_corr(dense_x, output = "matrix")
  attr(dense, "threshold") <- 0.35
  attr(dense, "thresholds") <- list(a = 0.1, b = -0.2)

  sparse <- tetrachoric(binary, output = "sparse", threshold = 0.25, diag = FALSE)
  edge <- tetrachoric(binary, output = "edge_list", threshold = 0.4, diag = FALSE)

  objects <- list(dense = dense, sparse = sparse, edge = edge)
  expected_threshold <- c(dense = 0.35, sparse = 0.25, edge = 0.4)

  for (nm in names(objects)) {
    overview <- overview_threshold(objects[[nm]])
    expect_identical(overview$threshold, expected_threshold[[nm]], info = nm)
    expect_true(is.finite(overview$threshold_sets), info = nm)
    expect_true(overview$threshold_sets > 0L, info = nm)
    printed <- capture.output(print(summary(objects[[nm]], show_ci = "no")))
    expect_true(
      any(grepl(formatC(expected_threshold[[nm]], format = "f", digits = 4), printed, fixed = TRUE)),
      info = nm
    )
    expect_true(any(grepl("set\\(s\\)", printed)), info = nm)
  }
})
