test_that("pearson_corr matches base::cor(..., method = 'pearson')", {
  for (i in 1:10) {
    n <- sample(c(10, 50, 100, 500), 1)
    p <- sample(2:6, 1)
    mat <- replicate(p, rnorm(n))
    colnames(mat) <- paste0("P", seq_len(p))

    base_cor <- cor(mat, method = "pearson")
    fast_cor <- suppressWarnings(pearson_corr(mat))

    attributes(fast_cor) <- NULL
    attributes(base_cor) <- NULL

    expect_equal(
      fast_cor,
      base_cor,
      tolerance = 1e-10,
      info = paste("Mismatch on Pearson test dataset", i, "n =", n, "p =", p)
    )
  }
})
