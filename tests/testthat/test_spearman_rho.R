test_that("spearman_rho matches base::cor(..., method = 'spearman')", {
  for (i in 1:10) {
    n <- sample(c(10, 50, 100, 500), 1)
    p <- sample(2:6, 1)
    mat <- replicate(p, rnorm(n))
    colnames(mat) <- paste0("S", seq_len(p))

    base_cor <- cor(mat, method = "spearman")
    fast_cor <- suppressWarnings(spearman_rho(mat))

    attributes(fast_cor) <- NULL
    attributes(base_cor) <- NULL

    expect_equal(
      fast_cor,
      base_cor,
      tolerance = 1e-8,
      info = paste("Mismatch on Spearman test dataset", i, "n =", n, "p =", p)
    )
  }
})

test_that("spearman_rho handles ties correctly and matches base::cor", {
  for (i in 1:5) {
    n <- sample(c(50, 100, 200), 1)
    p <- sample(2:5, 1)

    mat <- replicate(p, sample(rep(1:5, length.out = n)))
    colnames(mat) <- paste0("T", seq_len(p))

    base_cor <- cor(mat, method = "spearman")
    fast_cor <- suppressWarnings(spearman_rho(mat))

    attributes(fast_cor) <- NULL
    attributes(base_cor) <- NULL

    expect_equal(
      fast_cor,
      base_cor,
      tolerance = 1e-8,
      info = paste("Mismatch on tied Spearman dataset", i, "n =", n, "p =", p)
    )
  }
})
