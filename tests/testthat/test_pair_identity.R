# Positional pair identity when dimnames are duplicated or overlap.

payload_from_positions <- function(est, ci, diagnostics, inference) {
  idx <- which(upper.tri(est, diag = FALSE), arr.ind = TRUE)
  nmat <- diagnostics$n_complete
  data.frame(
    estimate = est[idx],
    n_complete = if (is.matrix(nmat)) as.integer(round(nmat[idx])) else rep(NA_integer_, nrow(idx)),
    lwr = if (is.list(ci) && is.matrix(ci$lwr.ci)) as.numeric(ci$lwr.ci[idx]) else rep(NA_real_, nrow(idx)),
    upr = if (is.list(ci) && is.matrix(ci$upr.ci)) as.numeric(ci$upr.ci[idx]) else rep(NA_real_, nrow(idx)),
    p_value = if (is.list(inference) && is.matrix(inference$p_value)) {
      as.numeric(inference$p_value[idx])
    } else {
      rep(NA_real_, nrow(idx))
    },
    stringsAsFactors = FALSE
  )
}

expect_positional_summary <- function(fit, est, ci, diagnostics, inference = NULL) {
  sm <- summary(fit, topn = 1000L, show_ci = "yes")
  exp <- payload_from_positions(est, ci, diagnostics, inference)
  sort_cols <- intersect(c("estimate", "lwr", "upr", "n_complete", "p_value"), names(sm))
  sort_cols <- sort_cols[vapply(sm[sort_cols], is.atomic, logical(1))]
  ord <- do.call(order, sm[sort_cols])
  exp <- exp[do.call(order, exp[intersect(sort_cols, names(exp))]), , drop = FALSE]
  sm <- sm[ord, , drop = FALSE]
  expect_equal(sm$estimate, exp$estimate)
  expect_identical(as.integer(sm$n_complete), exp$n_complete)
  if (any(is.finite(exp$lwr))) {
    expect_equal(sm$lwr, exp$lwr)
    expect_equal(sm$upr, exp$upr)
  }
  if (any(is.finite(exp$p_value))) {
    expect_equal(sm$p_value, exp$p_value)
  }
  sm
}

test_that("duplicate names keep positional estimates, CI, n, and p", {
  set.seed(11)
  x <- matrix(rnorm(24 * 4), 24, 4)
  colnames(x) <- c("a", "a", "b", "b")
  x[1:6, 1] <- NA
  x[3:9, 3] <- NA

  dense <- pearson_corr(x, ci = TRUE, na_method = "pairwise")
  est <- as.matrix(dense)
  expect_positional_summary(
    dense, est, attr(dense, "ci"), attr(dense, "diagnostics")
  )
  sparse <- pearson_corr(
    x, ci = TRUE, na_method = "pairwise", output = "sparse", threshold = 0, diag = FALSE
  )
  edge <- pearson_corr(
    x, ci = TRUE, na_method = "pairwise", output = "edge_list", threshold = 0, diag = FALSE
  )
  expect_positional_summary(sparse, est, attr(sparse, "ci"), attr(sparse, "diagnostics"))
  expect_positional_summary(edge, est, attr(edge, "ci"), attr(edge, "diagnostics"))

  set.seed(12)
  xw <- matrix(rnorm(16 * 4), 16, 4)
  colnames(xw) <- c("a", "a", "b", "b")
  xw[1:4, 2] <- NA
  robust <- wincor(xw, ci = FALSE, p_value = TRUE, na_method = "pairwise")
  expect_positional_summary(
    robust,
    as.matrix(robust),
    attr(robust, "ci"),
    attr(robust, "diagnostics"),
    attr(robust, "inference")
  )
  robust_edge <- wincor(
    xw, ci = FALSE, p_value = TRUE, na_method = "pairwise",
    output = "edge_list", threshold = 0, diag = FALSE
  )
  expect_positional_summary(
    robust_edge,
    as.matrix(robust),
    attr(robust_edge, "ci"),
    attr(robust_edge, "diagnostics"),
    attr(robust_edge, "inference")
  )
})

test_that("edge positions survive reordering and round-trip to dense", {
  set.seed(13)
  x <- matrix(rnorm(20 * 3), 20, 3)
  colnames(x) <- c("a", "a", "b")
  dense <- pearson_corr(x, ci = TRUE)
  edge <- pearson_corr(x, ci = TRUE, output = "edge_list", diag = TRUE)
  expect_false(is.null(attr(edge, "corr_pair_pos", exact = TRUE)))
  rebuilt <- matrixCorr:::.mc_corr_as_dense_matrix(edge)
  plain <- function(m) {
    attributes(m) <- attributes(m)[c("dim", "dimnames")]
    m
  }
  expect_identical(plain(rebuilt), plain(as.matrix(dense)))

  reordered <- edge[rev(seq_len(nrow(edge))), , drop = FALSE]
  expect_identical(plain(matrixCorr:::.mc_corr_as_dense_matrix(reordered)), plain(as.matrix(dense)))

  keep <- c(1L, nrow(edge))
  partial <- edge[keep, , drop = FALSE]
  partial_mat <- matrixCorr:::.mc_corr_as_dense_matrix(partial)
  pos <- attr(edge, "corr_pair_pos", exact = TRUE)
  keys <- rownames(edge)[keep]
  expect_identical(
    partial_mat[cbind(pos[keys, "i"], pos[keys, "j"])],
    edge$value[keep]
  )
})

test_that("direct triplets and directed squares use positions, not label equality", {
  triplets <- list(i = c(1L, 2L), j = c(2L, 1L), x = c(0.25, -0.8))
  dn <- list(c("a", "a"), c("a", "a"))
  edge <- matrixCorr:::.mc_finalize_triplets_output(
    triplets = triplets,
    output = "edge_list",
    estimator_class = "directed_assoc",
    method = "directed",
    description = "directed triplet",
    threshold = 0,
    diag = FALSE,
    source_dim = c(2L, 2L),
    source_dimnames = dn,
    symmetric = FALSE
  )
  expect_null(attr(edge, "ci", exact = TRUE))
  mat <- matrixCorr:::.mc_corr_as_dense_matrix(edge)
  expect_identical(mat[1, 2], 0.25)
  expect_identical(mat[2, 1], -0.8)
  pairs <- matrixCorr:::.mc_corr_pairs(edge, diag = FALSE)
  expect_identical(nrow(pairs), 2L)
  expect_identical(sort(pairs$value), sort(c(0.25, -0.8)))
})

test_that("rectangular overlap is not treated as a same-axis diagonal", {
  mat <- matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), nrow = 2, byrow = TRUE)
  dimnames(mat) <- list(c("a", "b"), c("a", "a", "b"))
  dense <- matrixCorr:::.mc_new_corr_matrix(
    mat = mat,
    estimator_class = "cross_assoc",
    method = "cross",
    description = "rectangular",
    symmetric = FALSE
  )
  edge <- matrixCorr:::.mc_finalize_corr_output(dense, output = "edge_list", threshold = 0, diag = FALSE)
  pairs <- matrixCorr:::.mc_corr_pairs(edge, diag = FALSE)
  # (1, 2) shares the label "a" with the dropped positional diagonal (1, 1).
  expect_identical(nrow(pairs), 4L)
  expect_true(any(pairs$i == 1L & pairs$j == 2L))
  expect_false(any(pairs$i == pairs$j))
  rebuilt <- matrixCorr:::.mc_corr_as_dense_matrix(edge)
  expect_true(is.na(rebuilt[1, 1]))
  expect_true(is.na(rebuilt[2, 2]))
  expect_identical(rebuilt[1, 2], mat[1, 2])
  expect_identical(rebuilt[2, 3], mat[2, 3])
})

test_that("ambiguous edges without a mapping or with a stale mapping fail", {
  set.seed(14)
  x <- matrix(rnorm(15 * 3), 15, 3)
  colnames(x) <- c("a", "a", "b")
  edge <- pearson_corr(x, output = "edge_list")
  old <- edge
  attr(old, "corr_pair_pos") <- NULL
  expect_error(matrixCorr:::.mc_corr_as_dense_matrix(old), "positional pair mapping")
  expect_error(summary(old), "positional pair mapping")

  stale <- edge
  pos <- attr(stale, "corr_pair_pos", exact = TRUE)
  pos[1L, 1L] <- nrow(x) + 5L
  attr(stale, "corr_pair_pos") <- pos
  expect_error(summary(stale), "positional pair mapping")
})

test_that("unique-name edge views do not gain or lose the positional attribute", {
  x <- matrix(c(1, 2, 3, 4, 2, 1, 0, 5, 3, 0, 1, 6, 4, 5, 6, 1), 4, 4)
  colnames(x) <- c("w", "x", "y", "z")
  edge <- pearson_corr(x, ci = TRUE, output = "edge_list")
  expect_null(attr(edge, "corr_pair_pos", exact = TRUE))
  view <- as.data.frame(edge)
  expect_identical(names(view), c("row", "col", "value"))
  expect_null(attr(view, "corr_pair_pos", exact = TRUE))
  est <- estimate(edge)
  expect_identical(names(est), c("item1", "item2", "estimate"))
  expect_null(attr(est, "corr_pair_pos", exact = TRUE))
  df <- matrixCorr:::.mc_plot_ci_df(
    pearson_corr(x, ci = TRUE),
    attr(pearson_corr(x, ci = TRUE), "ci"),
    "pearson"
  )
  expect_null(attr(df, "mc_axis_labels", exact = TRUE))
  expect_identical(levels(df$Var2), colnames(x))
})

test_that("duplicate-name plots keep one cell per position and the original labels", {
  set.seed(15)
  x <- matrix(rnorm(18 * 3), 18, 3)
  colnames(x) <- c("a", "a", "b")
  fit <- pearson_corr(x, ci = TRUE)
  ci <- attr(fit, "ci", exact = TRUE)
  df <- matrixCorr:::.mc_plot_ci_df(fit, ci, "pearson")
  expect_identical(nrow(df), 9L)
  ii <- as.integer(as.character(df$Var1))
  jj <- as.integer(as.character(df$Var2))
  expect_identical(df$pearson, as.numeric(as.matrix(fit)))
  lwr <- as.numeric(ci$lwr.ci)
  lwr[ii == jj] <- NA_real_
  # df is column-major, matching as.numeric()
  expect_identical(df$lwr, lwr)

  p <- plot(fit)
  expect_identical(nrow(ggplot2::layer_data(p, 1)), 9L)
  built <- ggplot2::ggplot_build(p)
  params <- built$layout$panel_params[[1]]
  expect_identical(params$x$get_labels(), colnames(x))
  expect_identical(params$y$get_labels(), rev(colnames(x)))
})

distinct_xy <- function(p) {
  ld <- ggplot2::layer_data(p, 1)
  nrow(unique(ld[, c("x", "y")]))
}

test_that("public and legacy pearson plots keep one coordinate per cell", {
  set.seed(16)
  x <- matrix(rnorm(18 * 3), 18, 3)
  colnames(x) <- c("same", "same", "other")
  fit <- pearson_corr(x, ci = TRUE)
  fit0 <- pearson_corr(x, ci = FALSE)
  expect_identical(distinct_xy(plot(fit)), 9L)
  expect_identical(distinct_xy(plot(fit0)), 9L)
  expect_identical(distinct_xy(plot.pearson_corr(fit)), 9L)
  legacy <- plot.pearson_corr(fit0)
  expect_identical(distinct_xy(legacy), 9L)
  params <- ggplot2::ggplot_build(legacy)$layout$panel_params[[1]]
  expect_identical(params$x$get_labels(), colnames(x))
  expect_identical(params$y$get_labels(), rev(colnames(x)))
})

test_that("bicor plots keep positions, triangles, and clustering", {
  set.seed(17)
  z <- rnorm(30)
  unique_x <- cbind(z, z + rnorm(30, sd = 0.05), rnorm(30))
  colnames(unique_x) <- c("a", "b", "c")
  bu <- bicor(unique_x, ci = FALSE)
  full <- plot.bicor(bu, reorder = "none", triangle = "full", show_value = FALSE)
  expect_identical(full$data$bicor, as.numeric(as.matrix(bu)))
  lower <- plot.bicor(bu, reorder = "none", triangle = "lower", show_value = FALSE)
  expect_identical(nrow(lower$data), 6L)
  clustered <- plot.bicor(bu, reorder = "hclust", triangle = "full", show_value = FALSE)
  labels <- ggplot2::ggplot_build(clustered)$layout$panel_params[[1]]$x$get_labels()
  plain <- ggplot2::ggplot_build(full)$layout$panel_params[[1]]$x$get_labels()
  expect_identical(sort(labels), sort(plain))
  expect_false(identical(labels, plain))

  dup_x <- unique_x
  colnames(dup_x) <- c("same", "same", "other")
  for (use_ci in c(FALSE, TRUE)) {
    bd <- bicor(dup_x, ci = use_ci)
    expect_identical(distinct_xy(plot.bicor(bd, triangle = "full", reorder = "none")), 9L)
    expect_identical(distinct_xy(plot.bicor(bd, triangle = "lower", reorder = "none")), 6L)
    expect_identical(distinct_xy(plot.bicor(bd, triangle = "upper", reorder = "hclust")), 6L)
    params <- ggplot2::ggplot_build(
      plot.bicor(bd, triangle = "full", reorder = "none")
    )$layout$panel_params[[1]]
    expect_identical(params$x$get_labels(), colnames(dup_x))
    expect_identical(params$y$get_labels(), rev(colnames(dup_x)))
  }
})

test_that("fractional, duplicated, and stale positional maps are rejected", {
  set.seed(18)
  x <- matrix(rnorm(15 * 3), 15, 3)
  colnames(x) <- c("a", "a", "b")
  edge <- pearson_corr(x, output = "edge_list")

  fractional <- edge
  pos <- attr(fractional, "corr_pair_pos", exact = TRUE)
  pos[1L, "i"] <- 1.5
  attr(fractional, "corr_pair_pos") <- pos
  expect_error(summary(fractional), "positional pair mapping")

  duplicated <- edge
  pos <- attr(duplicated, "corr_pair_pos", exact = TRUE)
  rownames(pos)[2L] <- rownames(pos)[1L]
  attr(duplicated, "corr_pair_pos") <- pos
  expect_error(summary(duplicated), "positional pair mapping")

  stale_labels <- edge
  dn <- attr(stale_labels, "corr_dimnames", exact = TRUE)
  dn[[1L]][1L] <- "not-the-stored-label"
  attr(stale_labels, "corr_dimnames") <- dn
  expect_error(summary(stale_labels), "positional pair mapping")
})

dup_axis_payload <- function() {
  dn <- list(c("a", "a"), c("a", "a"))
  lwr <- matrix(c(0.11, 0.31, 0.21, 0.41), 2, 2)
  upr <- lwr + 0.05
  list(
    dn = dn,
    ci = list(lwr.ci = lwr, upr.ci = upr),
    diagnostics = list(n_complete = matrix(c(11, 13, 12, 14), 2, 2)),
    inference = list(p_value = matrix(c(0.9, 0.8, 0.01, 0.7), 2, 2))
  )
}

test_that("legacy i/j/x edges keep positions when names are duplicated", {
  payload <- dup_axis_payload()
  legacy <- structure(
    data.frame(i = 1L, j = 2L, x = 0.5),
    class = c("corr_edge_list", "pearson_corr", "corr_result", "data.frame"),
    corr_dim = c(2L, 2L),
    corr_dimnames = payload$dn,
    corr_symmetric = TRUE,
    output = "edge_list",
    corr_estimator_class = "pearson_corr",
    ci = payload$ci,
    diagnostics = payload$diagnostics,
    inference = payload$inference
  )
  pairs <- matrixCorr:::.mc_edge_indexed_pairs(legacy)
  expect_identical(pairs$i, 1L)
  expect_identical(pairs$j, 2L)
  expect_identical(pairs$value, 0.5)
  expect_identical(pairs$row, "a")
  expect_identical(pairs$col, "a")

  sm <- summary(legacy, topn = 10L, show_ci = "yes")
  expect_identical(sm$estimate, 0.5)
  expect_equal(sm$lwr, 0.21)
  expect_equal(sm$upr, 0.26)
  expect_identical(as.integer(sm$n_complete), 12L)
  expect_equal(sm$p_value, 0.01)

  rebuilt <- matrixCorr:::.mc_corr_as_dense_matrix(legacy)
  expect_identical(rebuilt[1, 2], 0.5)
  expect_identical(rebuilt[2, 1], 0.5)
  expect_true(is.na(rebuilt[1, 1]))
  expect_true(is.na(rebuilt[2, 2]))

  fractional <- legacy
  fractional$i <- 1.5
  expect_error(summary(fractional), "positional pair mapping")

  constructed <- matrixCorr:::.mc_new_corr_edge_list(
    data.frame(i = 1L, j = 2L, x = 0.5),
    estimator_class = "pearson_corr",
    method = "pearson",
    description = "legacy triplet",
    source_dim = c(2L, 2L),
    source_dimnames = payload$dn,
    symmetric = TRUE,
    ci = payload$ci,
    diagnostics = payload$diagnostics,
    extra_attrs = list(inference = payload$inference)
  )
  expect_identical(
    class(constructed),
    c("corr_edge_list", "pearson_corr", "corr_result", "data.frame")
  )
  expect_identical(names(constructed), c("row", "col", "value"))
  pos <- attr(constructed, "corr_pair_pos", exact = TRUE)
  expect_identical(pos[, "i"], 1L)
  expect_identical(pos[, "j"], 2L)
  sm2 <- summary(constructed, topn = 10L, show_ci = "yes")
  expect_identical(sm2$estimate, 0.5)
  expect_equal(sm2$lwr, 0.21)
  expect_identical(as.integer(sm2$n_complete), 12L)
  expect_equal(sm2$p_value, 0.01)
  rebuilt2 <- matrixCorr:::.mc_corr_as_dense_matrix(constructed)
  expect_identical(rebuilt2[1, 2], 0.5)
  expect_identical(rebuilt2[2, 1], 0.5)

  unique_edge <- matrixCorr:::.mc_new_corr_edge_list(
    data.frame(i = 1L, j = 2L, x = 0.5),
    estimator_class = "pearson_corr",
    method = "pearson",
    description = "unique legacy triplet",
    source_dim = c(2L, 2L),
    source_dimnames = list(c("a", "b"), c("a", "b")),
    symmetric = TRUE
  )
  expect_null(attr(unique_edge, "corr_pair_pos", exact = TRUE))
  expect_identical(unique_edge$row, "a")
  expect_identical(unique_edge$col, "b")
})

axis_labels_of <- function(p) {
  params <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]
  list(x = params$x$get_labels(), y = params$y$get_labels())
}

expect_positional_heatmap <- function(p, labels, fit, ci = NULL) {
  expect_identical(distinct_xy(p), 9L)
  labs <- axis_labels_of(p)
  expect_identical(labs$x, labels)
  expect_identical(labs$y, rev(labels))
  df <- p$data
  if (is.null(ci)) {
    return(invisible(df))
  }
  value_col <- intersect(c("value", "pearson", "bicor", "Pearson"), names(df))[[1L]]
  row_key <- if ("row" %in% names(df)) df$row else df$Var1
  col_key <- if ("col" %in% names(df)) df$col else df$Var2
  ii <- as.integer(as.character(row_key))
  jj <- as.integer(as.character(col_key))
  expect_false(anyNA(ii) || anyNA(jj))
  expect_equal(df[[value_col]], as.numeric(as.matrix(fit)[cbind(ii, jj)]))
  lwr_src <- as.numeric(ci$lwr.ci[cbind(ii, jj)])
  upr_src <- as.numeric(ci$upr.ci[cbind(ii, jj)])
  if ("ci_label" %in% names(df)) {
    off <- ii != jj
    shown <- !is.na(df$ci_label)
    expect_true(all(shown[off]))
    expect_false(any(shown[!off]))
    parsed <- sub("^\\[(.*), (.*)\\]$", "\\1", df$ci_label[off])
    expect_equal(as.numeric(parsed), lwr_src[off], tolerance = 1e-3)
  } else {
    expect_equal(df$lwr, lwr_src)
    expect_equal(df$upr, upr_src)
  }
  invisible(df)
}

test_that("one missing or blank label uses positional axes without renaming", {
  set.seed(19)
  base <- matrix(rnorm(24 * 3), 24, 3)
  cases <- list(
    na = c("a", NA, "c"),
    blank = c("", "b", "c")
  )
  for (labels in cases) {
    x <- base
    colnames(x) <- labels
    fit <- pearson_corr(x, ci = TRUE)
    fit0 <- pearson_corr(x, ci = FALSE)
    ci <- attr(fit, "ci", exact = TRUE)
    expect_positional_heatmap(plot(fit), labels, fit, ci)
    expect_positional_heatmap(plot(fit0), labels, fit0)
    expect_positional_heatmap(plot.pearson_corr(fit), labels, fit, ci)
    expect_positional_heatmap(plot.pearson_corr(fit0), labels, fit0)

    edge <- pearson_corr(x, output = "edge_list")
    expect_null(attr(edge, "corr_pair_pos", exact = TRUE))
    expect_s3_class(summary(edge, topn = 20L), "data.frame")

    for (use_ci in c(FALSE, TRUE)) {
      bd <- bicor(x, ci = use_ci)
      p <- plot.bicor(bd, reorder = "none", triangle = "full")
      if (isTRUE(use_ci)) {
        expect_positional_heatmap(p, labels, bd, attr(bd, "ci", exact = TRUE))
      } else {
        expect_positional_heatmap(p, labels, bd)
      }
      expect_identical(distinct_xy(plot.bicor(bd, triangle = "lower", reorder = "none")), 6L)
      expect_identical(distinct_xy(plot.bicor(bd, triangle = "upper", reorder = "hclust")), 6L)
    }
  }
})

expect_cross_square_payload <- function(fit) {
  est <- as.matrix(fit)
  idx <- which(matrix(TRUE, nrow(est), ncol(est)), arr.ind = TRUE)
  sm <- summary(fit, topn = 20L, show_ci = "yes")
  td <- tidy(fit)
  expect_identical(nrow(sm), 4L)
  expect_identical(nrow(td), 4L)
  sm_pos <- cbind(match(sm$item1, rownames(est)), match(sm$item2, colnames(est)))
  td_pos <- cbind(match(td$item1, rownames(est)), match(td$item2, colnames(est)))
  expect_false(anyNA(sm_pos) || anyNA(td_pos))
  expect_equal(sm$estimate, est[sm_pos])
  expect_equal(td$estimate, est[td_pos])
  expect_equal(td_pos, idx, ignore_attr = TRUE)
  ci <- attr(fit, "ci", exact = TRUE)
  inf <- attr(fit, "inference", exact = TRUE)
  nmat <- attr(fit, "diagnostics", exact = TRUE)$n_complete
  expect_equal(sm$lwr, as.numeric(ci$lwr.ci[sm_pos]))
  expect_equal(sm$upr, as.numeric(ci$upr.ci[sm_pos]))
  expect_equal(td$lwr, as.numeric(ci$lwr.ci[td_pos]))
  expect_equal(td$upr, as.numeric(ci$upr.ci[td_pos]))
  expect_equal(sm$p_value, as.numeric(inf$p_value[sm_pos]))
  expect_equal(td$p_value, as.numeric(inf$p_value[td_pos]))
  expect_identical(as.integer(sm$n_complete), as.integer(round(nmat[sm_pos])))
  expect_false("n_complete" %in% names(td))
  expect_identical(nrow(tidy(fit, triangle = "upper")), 1L)
  expect_identical(nrow(tidy(fit, triangle = "lower")), 1L)
  expect_identical(nrow(tidy(fit, triangle = "full")), 4L)
  p <- plot(fit)
  same <- as.character(p$data$row) == as.character(p$data$col)
  expect_true(any(same))
  expect_true(all(!is.na(p$data$ci_label)))
  parsed <- sub("^\\[(.*), (.*)\\]$", "\\1", p$data$ci_label[same])
  ii <- match(as.character(p$data$row[same]), rownames(est))
  jj <- match(as.character(p$data$col[same]), colnames(est))
  expect_false(anyNA(ii) || anyNA(jj))
  expect_equal(as.numeric(parsed), as.numeric(ci$lwr.ci[cbind(ii, jj)]), tolerance = 5e-4)
}

test_that("same-name biserial and polyserial squares keep every cross pair", {
  set.seed(44)
  x <- matrix(rnorm(80), 40, 2, dimnames = list(NULL, c("a", "b")))
  y_bin <- 1L * (matrix(rnorm(80), 40, 2) > 0)
  colnames(y_bin) <- c("a", "b")
  y_ord <- y_bin + 1L
  expect_cross_square_payload(biserial(x, y_bin, ci = TRUE, p_value = TRUE))
  expect_cross_square_payload(polyserial(x, y_ord, ci = TRUE, p_value = TRUE))

  pear <- pearson_corr(x, ci = TRUE)
  expect_identical(nrow(summary(pear, topn = 20L)), 1L)
  expect_identical(nrow(tidy(pear)), 1L)
  pear_plot <- plot(pear)
  pear_same <- as.character(pear_plot$data$row) == as.character(pear_plot$data$col)
  expect_true(any(pear_same))
  expect_true(all(is.na(pear_plot$data$ci_label[pear_same])))

  xi <- xi_corr(x)
  expect_identical(nrow(tidy(xi)), 2L)
  expect_identical(nrow(summary(xi, topn = 20L)), 2L)
  expect_false(any(tidy(xi)$item1 == tidy(xi)$item2))
})
