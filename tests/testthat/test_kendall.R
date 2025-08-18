test_that("kendall_tau matches base::cor(..., method = 'kendall')", {
  set.seed(123)

  for (i in 1:10) {
    n <- sample(c(10, 50, 100, 500), 1)
    p <- sample(2:6, 1)
    mat <- replicate(p, rnorm(n))
    colnames(mat) <- paste0("V", seq_len(p))

    base_cor <- cor(mat, method = "kendall")
    fast_cor <- suppressWarnings(kendall_tau(mat))

    attributes(fast_cor) <- NULL
    attributes(base_cor) <- NULL

    expect_equal(
      fast_cor,
      base_cor,
      tolerance = 1e-8,
      info = paste("Mismatch on test dataset", i, "n =", n, "p =", p)
    )
  }
})

test_that("kendall_tau handles ties correctly and matches base::cor", {
  set.seed(456)

  for (i in 1:5) {
    n <- sample(c(50, 100, 200), 1)
    p <- sample(2:5, 1)

    # Create tied data
    mat <- replicate(p, sample(rep(1:5, length.out = n)))
    colnames(mat) <- paste0("T", seq_len(p))

    base_cor <- cor(mat, method = "kendall")
    fast_cor <- suppressWarnings(kendall_tau(mat))

    attributes(fast_cor) <- NULL
    attributes(base_cor) <- NULL

    expect_equal(
      fast_cor,
      base_cor,
      tolerance = 1e-8,
      info = paste("Mismatch on tied test dataset", i, "n =", n, "p =", p)
    )
  }
})
