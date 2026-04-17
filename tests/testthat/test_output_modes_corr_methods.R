expected_edge_df <- function(mat, threshold = 0, diag = TRUE) {
  idx <- upper.tri(mat, diag = diag)
  i <- row(mat)[idx]
  j <- col(mat)[idx]
  vals <- mat[idx]
  keep <- is.finite(vals) & abs(vals) >= threshold
  rn <- rownames(mat)
  cn <- colnames(mat)
  row_out <- if (is.null(rn)) as.character(i[keep]) else rn[i[keep]]
  col_out <- if (is.null(cn)) as.character(j[keep]) else cn[j[keep]]
  data.frame(
    row = row_out,
    col = col_out,
    value = as.numeric(vals[keep]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

expected_sparse_dense <- function(mat, threshold = 0, diag = TRUE) {
  out <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  idx <- upper.tri(mat, diag = diag)
  i <- row(mat)[idx]
  j <- col(mat)[idx]
  vals <- mat[idx]
  keep <- is.finite(vals) & abs(vals) >= threshold
  out[cbind(i[keep], j[keep])] <- vals[keep]
  out[cbind(j[keep], i[keep])] <- vals[keep]
  out
}

test_that("output modes are consistent across dense correlation methods", {
  set.seed(20260415)
  X <- matrix(rnorm(200 * 5), nrow = 200, ncol = 5)
  colnames(X) <- paste0("V", seq_len(ncol(X)))

  methods <- list(
    pearson = function(...) pearson_corr(X, na_method = "error", ci = FALSE, ...),
    spearman = function(...) spearman_rho(X, na_method = "error", ci = FALSE, ...),
    kendall = function(...) kendall_tau(X, na_method = "error", ci = FALSE, ...),
    bicor = function(...) bicor(X, na_method = "error", ci = FALSE, ...),
    dcor = function(...) dcor(X, na_method = "error", p_value = FALSE, ...),
    pbcor = function(...) pbcor(X, na_method = "error", ci = FALSE, p_value = FALSE, ...),
    wincor = function(...) wincor(X, na_method = "error", ci = FALSE, p_value = FALSE, ...),
    skipped = function(...) skipped_corr(
      X,
      method = "pearson",
      na_method = "error",
      ci = FALSE,
      p_value = FALSE,
      return_masks = FALSE,
      ...
    ),
    shrinkage = function(...) shrinkage_corr(X, ...),
    ccc = function(...) ccc(X, ci = FALSE, ...)
  )

  for (nm in names(methods)) {
    fn <- methods[[nm]]

    base <- fn()
    mat <- as.matrix(base)
    expect_equal(colnames(mat), colnames(X), info = nm)
    expect_equal(rownames(mat), colnames(X), info = nm)

    explicit_matrix <- fn(output = "matrix", threshold = 0, diag = TRUE)
    expect_equal(as.matrix(explicit_matrix), mat, tolerance = 1e-12, info = nm)

    expect_error(
      fn(output = "matrix", threshold = 0.25),
      "must be 0 when",
      info = nm
    )

    expect_error(fn(output = "packed_upper"), info = nm)

    thr <- 0.30
    edge <- fn(output = "edge_list", threshold = thr, diag = FALSE)
    expected_edge <- expected_edge_df(mat, threshold = thr, diag = FALSE)
    edge_df <- as.data.frame(edge, stringsAsFactors = FALSE)[, c("row", "col", "value")]
    expect_equal(edge_df, expected_edge, tolerance = 1e-12, info = nm)
    if (nrow(edge) > 0L) {
      expect_true(all(abs(edge$value) >= thr), info = nm)
    }

    sparse <- fn(output = "sparse", threshold = thr, diag = FALSE)
    expect_s4_class(sparse, "sparseMatrix")
    expect_equal(dimnames(sparse), dimnames(mat), info = nm)
    dense_sparse <- as.matrix(sparse)
    expect_true(isSymmetric(dense_sparse), info = nm)
    expect_equal(
      dense_sparse,
      expected_sparse_dense(mat, threshold = thr, diag = FALSE),
      tolerance = 1e-12,
      info = nm
    )
  }
})

test_that("eligible methods use direct triplet payloads and match dense filtering", {
  set.seed(20260416)
  X <- matrix(rnorm(320 * 6), nrow = 320, ncol = 6)
  colnames(X) <- paste0("T", seq_len(ncol(X)))

  methods <- list(
    pearson = function(...) pearson_corr(X, na_method = "error", ci = FALSE, ...),
    spearman = function(...) spearman_rho(X, na_method = "error", ci = FALSE, ...),
    ccc = function(...) ccc(X, ci = FALSE, ...),
    bicor = function(...) bicor(X, na_method = "error", ci = FALSE, ...),
    pbcor = function(...) pbcor(X, na_method = "error", ci = FALSE, p_value = FALSE, ...),
    wincor = function(...) wincor(X, na_method = "error", ci = FALSE, p_value = FALSE, ...)
  )

  thr <- 0.35
  for (nm in names(methods)) {
    fn <- methods[[nm]]
    dense <- fn(output = "matrix")
    mat <- as.matrix(dense)

    edge <- fn(output = "edge_list", threshold = thr, diag = FALSE)
    expect_true(all(c("i", "j", "x") %in% names(edge)), info = nm)
    edge_df <- .mc_corr_as_edge_df(edge)
    expect_equal(
      edge_df,
      expected_edge_df(mat, threshold = thr, diag = FALSE),
      tolerance = 1e-12,
      info = nm
    )

    sparse <- fn(output = "sparse", threshold = thr, diag = FALSE)
    expect_s4_class(sparse, "sparseMatrix")
    expect_equal(dimnames(sparse), dimnames(mat), info = nm)
    expect_true(isSymmetric(as.matrix(sparse)), info = nm)
    expect_equal(
      as.matrix(sparse),
      expected_sparse_dense(mat, threshold = thr, diag = FALSE),
      tolerance = 1e-12,
      info = nm
    )
  }
})

test_that("thresholded direct outputs are materially smaller than dense outputs", {
  set.seed(20260417)
  X <- matrix(rnorm(250 * 150), nrow = 250, ncol = 150)
  colnames(X) <- paste0("M", seq_len(ncol(X)))

  methods <- list(
    pearson = function(...) pearson_corr(X, na_method = "error", ci = FALSE, ...),
    spearman = function(...) spearman_rho(X, na_method = "error", ci = FALSE, ...),
    ccc = function(...) ccc(X, ci = FALSE, ...),
    bicor = function(...) bicor(X, na_method = "error", ci = FALSE, ...),
    pbcor = function(...) pbcor(X, na_method = "error", ci = FALSE, p_value = FALSE, ...),
    wincor = function(...) wincor(X, na_method = "error", ci = FALSE, p_value = FALSE, ...)
  )

  thr <- 0.90
  for (nm in names(methods)) {
    fn <- methods[[nm]]
    dense <- fn(output = "matrix")
    sparse <- fn(output = "sparse", threshold = thr, diag = FALSE)
    edge <- fn(output = "edge_list", threshold = thr, diag = FALSE)

    dense_size <- as.numeric(object.size(as.matrix(dense)))
    sparse_size <- as.numeric(object.size(sparse))
    edge_size <- as.numeric(object.size(edge))
    retained_ratio <- Matrix::nnzero(sparse) / (ncol(X) * ncol(X))

    expect_lt(retained_ratio, 0.20)

    expect_lt(sparse_size, dense_size * 0.90)
    expect_lt(edge_size, dense_size * 0.90)
  }
})

test_that("pearson edge_list threshold=0 keeps edge-list payload contract", {
  set.seed(20260417)
  X <- matrix(rnorm(240 * 8), nrow = 240, ncol = 8)
  colnames(X) <- paste0("E", seq_len(ncol(X)))

  mat <- as.matrix(pearson_corr(X, na_method = "error", ci = FALSE, output = "matrix"))
  edge <- pearson_corr(X, na_method = "error", ci = FALSE, output = "edge_list", threshold = 0, diag = FALSE)

  expect_true(all(c("row", "col", "value") %in% names(edge)))
  expect_false(any(c("i", "j", "x") %in% names(edge)))

  got <- as.data.frame(edge, stringsAsFactors = FALSE)[, c("row", "col", "value")]
  expected <- expected_edge_df(mat, threshold = 0, diag = FALSE)

  got <- got[order(got$col, got$row), , drop = FALSE]
  expected <- expected[order(expected$col, expected$row), , drop = FALSE]
  rownames(got) <- NULL
  rownames(expected) <- NULL

  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("pcorr supports output modes for point-estimate path", {
  set.seed(11)
  X <- matrix(rnorm(180 * 5), nrow = 180, ncol = 5)
  colnames(X) <- paste0("P", seq_len(ncol(X)))

  fit <- pcorr(X, method = "sample")
  mat <- as.matrix(fit$pcor)
  expect_equal(colnames(mat), colnames(X))
  expect_equal(rownames(mat), colnames(X))

  expect_error(
    pcorr(X, method = "sample", output = "matrix", threshold = 0.2),
    "must be 0 when"
  )
  expect_error(
    pcorr(X, method = "sample", output = "edge_list", return_cov_precision = TRUE),
    "point estimates only"
  )

  expect_error(pcorr(X, method = "sample", output = "packed_upper"))

  edge <- pcorr(X, method = "sample", output = "edge_list", threshold = 0.25, diag = FALSE)
  edge_df <- as.data.frame(edge, stringsAsFactors = FALSE)[, c("row", "col", "value")]
  expect_equal(edge_df, expected_edge_df(mat, threshold = 0.25, diag = FALSE), tolerance = 1e-12)

  sparse <- pcorr(X, method = "sample", output = "sparse", threshold = 0.25, diag = FALSE)
  expect_s4_class(sparse, "sparseMatrix")
  expect_equal(dimnames(sparse), dimnames(mat))
  expect_equal(
    as.matrix(sparse),
    expected_sparse_dense(mat, threshold = 0.25, diag = FALSE),
    tolerance = 1e-12
  )
})

test_that("latent symmetric methods support matrix/sparse/edge outputs", {
  set.seed(4242)
  Z <- matrix(rnorm(600 * 4), nrow = 600, ncol = 4)
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

  methods <- list(
    tetrachoric = function(...) tetrachoric(X_bin, ci = FALSE, p_value = FALSE, ...),
    polychoric = function(...) polychoric(X_ord, ci = FALSE, p_value = FALSE, ...)
  )

  for (nm in names(methods)) {
    fn <- methods[[nm]]
    base <- fn()
    mat <- as.matrix(base)

    expect_error(fn(output = "matrix", threshold = 0.2), "must be 0 when", info = nm)

    expect_error(fn(output = "packed_upper"), info = nm)

    edge <- fn(output = "edge_list", threshold = 0.25, diag = FALSE)
    edge_df <- as.data.frame(edge, stringsAsFactors = FALSE)[, c("row", "col", "value")]
    expect_equal(edge_df, expected_edge_df(mat, threshold = 0.25, diag = FALSE), tolerance = 1e-12, info = nm)

    sparse <- fn(output = "sparse", threshold = 0.25, diag = FALSE)
    expect_s4_class(sparse, "sparseMatrix")
    expect_equal(as.matrix(sparse), expected_sparse_dense(mat, threshold = 0.25, diag = FALSE), tolerance = 1e-12, info = nm)
  }
})

test_that("rectangular latent methods keep legacy interface without output arguments", {
  set.seed(5252)
  Z <- matrix(rnorm(400 * 4), nrow = 400, ncol = 4)
  X <- data.frame(x1 = Z[, 1], x2 = Z[, 2])
  Y_ord <- data.frame(
    y1 = ordered(cut(Z[, 3], breaks = c(-Inf, -0.3, 0.5, Inf), labels = c("L", "M", "H"))),
    y2 = ordered(cut(Z[, 4], breaks = c(-Inf, -0.8, 0.0, 0.8, Inf), labels = c("1", "2", "3", "4")))
  )
  Y_bin <- data.frame(
    g1 = Z[, 3] > 0,
    g2 = Z[, 4] > -0.2
  )

  expect_error(polyserial(X, Y_ord, output = "sparse"), "unsupported argument\\(s\\): output")
  expect_error(biserial(X, Y_bin, output = "sparse"), "unsupported argument\\(s\\): output")
})
