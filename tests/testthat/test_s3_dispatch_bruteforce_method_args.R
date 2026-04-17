grid_df <- function(...) {
  expand.grid(..., KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
}

dispatch_expectations <- function(obj, label, output) {
  sm <- summary(obj)

  if (inherits(obj, "partial_corr")) {
    expect_true(inherits(sm, "summary.partial_corr"), info = label)
  } else if (inherits(obj, "ccc") && identical(output, "matrix")) {
    expect_true(inherits(sm, "summary.ccc"), info = label)
  } else {
    expect_true(inherits(sm, "summary.corr_result"), info = label)
  }

  if (inherits(sm, "summary.corr_result")) {
    expect_false(inherits(sm, "summaryDefault"), info = label)
  }

  txt <- capture.output(print(sm))
  expect_true(length(txt) > 0L, info = label)
}

select_data_by_na <- function(x_error, x_pairwise, na_method) {
  if (identical(na_method, "pairwise")) x_pairwise else x_error
}

valid_threshold <- function(output, threshold) {
  !identical(output, "matrix") || identical(threshold, 0)
}

test_that("bruteforce dispatch: pearson/spearman/kendall/bicor/dcor", {
  set.seed(20260420)
  X <- matrix(rnorm(150 * 5), nrow = 150, ncol = 5)
  colnames(X) <- paste0("V", seq_len(ncol(X)))
  X_na <- X
  X_na[sample.int(length(X_na), size = 25)] <- NA_real_

  out_grid <- grid_df(
    output = c("matrix", "sparse", "edge_list"),
    threshold = c(0, 0.2),
    diag = c(TRUE, FALSE)
  )
  out_grid <- out_grid[vapply(
    seq_len(nrow(out_grid)),
    function(i) valid_threshold(out_grid$output[[i]], out_grid$threshold[[i]]),
    logical(1)
  ), , drop = FALSE]

  na_ci_grid <- grid_df(
    na_method = c("error", "pairwise"),
    ci = c(FALSE, TRUE),
    conf_level = c(0.90, 0.95),
    n_threads = c(1L, 2L)
  )

  for (i in seq_len(nrow(na_ci_grid))) {
    for (j in seq_len(nrow(out_grid))) {
      row <- na_ci_grid[i, , drop = FALSE]
      out <- out_grid[j, , drop = FALSE]
      x_use <- select_data_by_na(X, X_na, row$na_method[[1]])

      lbl <- sprintf(
        "pearson na=%s ci=%s cl=%.2f thr=%.2f diag=%s nt=%s out=%s",
        row$na_method[[1]], row$ci[[1]], row$conf_level[[1]],
        out$threshold[[1]], out$diag[[1]], row$n_threads[[1]], out$output[[1]]
      )
      obj <- pearson_corr(
        x_use,
        na_method = row$na_method[[1]],
        ci = row$ci[[1]],
        conf_level = row$conf_level[[1]],
        n_threads = row$n_threads[[1]],
        output = out$output[[1]],
        threshold = out$threshold[[1]],
        diag = out$diag[[1]]
      )
      dispatch_expectations(obj, lbl, out$output[[1]])

      lbl <- sprintf(
        "spearman na=%s ci=%s cl=%.2f thr=%.2f diag=%s nt=%s out=%s",
        row$na_method[[1]], row$ci[[1]], row$conf_level[[1]],
        out$threshold[[1]], out$diag[[1]], row$n_threads[[1]], out$output[[1]]
      )
      obj <- spearman_rho(
        x_use,
        na_method = row$na_method[[1]],
        ci = row$ci[[1]],
        conf_level = row$conf_level[[1]],
        n_threads = row$n_threads[[1]],
        output = out$output[[1]],
        threshold = out$threshold[[1]],
        diag = out$diag[[1]]
      )
      dispatch_expectations(obj, lbl, out$output[[1]])
    }
  }

  kend_grid <- grid_df(
    na_method = c("error", "pairwise"),
    ci = c(FALSE, TRUE),
    ci_method = c("fieller", "if_el", "brown_benedetti"),
    conf_level = c(0.90, 0.95)
  )

  for (i in seq_len(nrow(kend_grid))) {
    for (j in seq_len(nrow(out_grid))) {
      row <- kend_grid[i, , drop = FALSE]
      out <- out_grid[j, , drop = FALSE]
      x_use <- select_data_by_na(X, X_na, row$na_method[[1]])
      lbl <- sprintf(
        "kendall na=%s ci=%s cim=%s cl=%.2f thr=%.2f diag=%s out=%s",
        row$na_method[[1]], row$ci[[1]], row$ci_method[[1]], row$conf_level[[1]],
        out$threshold[[1]], out$diag[[1]], out$output[[1]]
      )
      obj <- kendall_tau(
        x_use,
        na_method = row$na_method[[1]],
        ci = row$ci[[1]],
        ci_method = row$ci_method[[1]],
        conf_level = row$conf_level[[1]],
        n_threads = 1L,
        output = out$output[[1]],
        threshold = out$threshold[[1]],
        diag = out$diag[[1]]
      )
      dispatch_expectations(obj, lbl, out$output[[1]])
    }
  }

  bicor_grid <- grid_df(
    na_method = c("error", "pairwise"),
    ci = c(FALSE, TRUE),
    c_const = c(6, 9),
    max_p_outliers = c(0.5, 1),
    pearson_fallback = c("hybrid", "none", "all"),
    mad_consistent = c(FALSE, TRUE)
  )

  for (i in seq_len(nrow(bicor_grid))) {
    for (j in seq_len(nrow(out_grid))) {
      row <- bicor_grid[i, , drop = FALSE]
      out <- out_grid[j, , drop = FALSE]
      x_use <- select_data_by_na(X, X_na, row$na_method[[1]])
      lbl <- sprintf(
        "bicor na=%s ci=%s c=%.1f mpo=%.2f pf=%s mad=%s thr=%.2f diag=%s out=%s",
        row$na_method[[1]], row$ci[[1]], row$c_const[[1]], row$max_p_outliers[[1]],
        row$pearson_fallback[[1]], row$mad_consistent[[1]],
        out$threshold[[1]], out$diag[[1]], out$output[[1]]
      )
      obj <- bicor(
        x_use,
        na_method = row$na_method[[1]],
        ci = row$ci[[1]],
        conf_level = 0.95,
        n_threads = 1L,
        output = out$output[[1]],
        threshold = out$threshold[[1]],
        diag = out$diag[[1]],
        c_const = row$c_const[[1]],
        max_p_outliers = row$max_p_outliers[[1]],
        pearson_fallback = row$pearson_fallback[[1]],
        mad_consistent = row$mad_consistent[[1]]
      )
      dispatch_expectations(obj, lbl, out$output[[1]])
    }
  }

  dcor_grid <- grid_df(
    na_method = c("error", "pairwise"),
    p_value = c(FALSE, TRUE),
    n_threads = c(1L, 2L)
  )

  for (i in seq_len(nrow(dcor_grid))) {
    for (j in seq_len(nrow(out_grid))) {
      row <- dcor_grid[i, , drop = FALSE]
      out <- out_grid[j, , drop = FALSE]
      x_use <- select_data_by_na(X, X_na, row$na_method[[1]])
      lbl <- sprintf(
        "dcor na=%s p=%s thr=%.2f diag=%s nt=%s out=%s",
        row$na_method[[1]], row$p_value[[1]], out$threshold[[1]],
        out$diag[[1]], row$n_threads[[1]], out$output[[1]]
      )
      obj <- dcor(
        x_use,
        na_method = row$na_method[[1]],
        p_value = row$p_value[[1]],
        n_threads = row$n_threads[[1]],
        output = out$output[[1]],
        threshold = out$threshold[[1]],
        diag = out$diag[[1]]
      )
      dispatch_expectations(obj, lbl, out$output[[1]])
    }
  }
})

test_that("bruteforce dispatch: pbcor/wincor/skipped with tuning arguments", {
  set.seed(20260421)
  X <- matrix(rnorm(150 * 5), nrow = 150, ncol = 5)
  colnames(X) <- paste0("R", seq_len(ncol(X)))
  X_na <- X
  X_na[sample.int(length(X_na), size = 25)] <- NA_real_

  out_grid <- grid_df(
    output = c("matrix", "sparse", "edge_list"),
    threshold = c(0, 0.2),
    diag = c(TRUE, FALSE)
  )
  out_grid <- out_grid[vapply(
    seq_len(nrow(out_grid)),
    function(i) valid_threshold(out_grid$output[[i]], out_grid$threshold[[i]]),
    logical(1)
  ), , drop = FALSE]

  infer_grid <- grid_df(
    na_method = c("error", "pairwise"),
    ci = c(FALSE, TRUE),
    p_value = c(FALSE, TRUE)
  )

  for (i in seq_len(nrow(infer_grid))) {
    row <- infer_grid[i, , drop = FALSE]
    x_use <- select_data_by_na(X, X_na, row$na_method[[1]])
    for (j in seq_len(nrow(out_grid))) {
      out <- out_grid[j, , drop = FALSE]

      lbl <- sprintf(
        "pbcor na=%s ci=%s p=%s thr=%.2f diag=%s out=%s",
        row$na_method[[1]], row$ci[[1]], row$p_value[[1]],
        out$threshold[[1]], out$diag[[1]], out$output[[1]]
      )
      obj <- pbcor(
        x_use,
        na_method = row$na_method[[1]],
        ci = row$ci[[1]],
        p_value = row$p_value[[1]],
        conf_level = 0.95,
        n_threads = 1L,
        beta = 0.2,
        n_boot = 20L,
        seed = 11L,
        output = out$output[[1]],
        threshold = out$threshold[[1]],
        diag = out$diag[[1]]
      )
      dispatch_expectations(obj, lbl, out$output[[1]])

      lbl <- sprintf(
        "wincor na=%s ci=%s p=%s thr=%.2f diag=%s out=%s",
        row$na_method[[1]], row$ci[[1]], row$p_value[[1]],
        out$threshold[[1]], out$diag[[1]], out$output[[1]]
      )
      obj <- wincor(
        x_use,
        na_method = row$na_method[[1]],
        ci = row$ci[[1]],
        p_value = row$p_value[[1]],
        conf_level = 0.95,
        n_threads = 1L,
        tr = 0.2,
        n_boot = 20L,
        seed = 11L,
        output = out$output[[1]],
        threshold = out$threshold[[1]],
        diag = out$diag[[1]]
      )
      dispatch_expectations(obj, lbl, out$output[[1]])
    }
  }

  skip_grid <- grid_df(
    method = c("pearson", "spearman"),
    na_method = c("error", "pairwise"),
    ci = c(FALSE, TRUE),
    p_value = c(FALSE, TRUE),
    stand = c(FALSE, TRUE),
    outlier_rule = c("idealf", "mad"),
    p_adjust = c("none", "hochberg")
  )

  for (i in seq_len(nrow(skip_grid))) {
    row <- skip_grid[i, , drop = FALSE]
    if ((row$ci[[1]] || row$p_value[[1]]) && !identical(row$na_method[[1]], "error")) {
      next
    }
    if (!row$p_value[[1]] && !identical(row$p_adjust[[1]], "none")) {
      next
    }
    x_use <- select_data_by_na(X, X_na, row$na_method[[1]])
    for (j in seq_len(nrow(out_grid))) {
      out <- out_grid[j, , drop = FALSE]
      lbl <- sprintf(
        "skip m=%s na=%s ci=%s p=%s st=%s or=%s pa=%s thr=%.2f diag=%s out=%s",
        row$method[[1]], row$na_method[[1]], row$ci[[1]], row$p_value[[1]],
        row$stand[[1]], row$outlier_rule[[1]], row$p_adjust[[1]],
        out$threshold[[1]], out$diag[[1]], out$output[[1]]
      )
      obj <- skipped_corr(
        x_use,
        method = row$method[[1]],
        na_method = row$na_method[[1]],
        ci = row$ci[[1]],
        p_value = row$p_value[[1]],
        conf_level = 0.95,
        n_threads = 1L,
        return_masks = FALSE,
        stand = row$stand[[1]],
        outlier_rule = row$outlier_rule[[1]],
        cutoff = sqrt(stats::qchisq(0.975, df = 2)),
        n_boot = 30L,
        p_adjust = row$p_adjust[[1]],
        fwe_level = 0.05,
        n_mc = 50L,
        seed = 11L,
        output = out$output[[1]],
        threshold = out$threshold[[1]],
        diag = out$diag[[1]]
      )
      dispatch_expectations(obj, lbl, out$output[[1]])
    }
  }
})

test_that("bruteforce dispatch: ccc/shrinkage/schafer and latent methods", {
  set.seed(20260422)
  X <- matrix(rnorm(150 * 5), nrow = 150, ncol = 5)
  colnames(X) <- paste0("C", seq_len(ncol(X)))

  out_grid <- grid_df(
    output = c("matrix", "sparse", "edge_list"),
    threshold = c(0, 0.2),
    diag = c(TRUE, FALSE)
  )
  out_grid <- out_grid[vapply(
    seq_len(nrow(out_grid)),
    function(i) valid_threshold(out_grid$output[[i]], out_grid$threshold[[i]]),
    logical(1)
  ), , drop = FALSE]

  for (i in seq_len(nrow(out_grid))) {
    out <- out_grid[i, , drop = FALSE]
    for (ci in c(FALSE, TRUE)) {
      lbl <- sprintf(
        "ccc ci=%s thr=%.2f diag=%s out=%s",
        ci, out$threshold[[1]], out$diag[[1]], out$output[[1]]
      )
      obj <- ccc(
        X,
        ci = ci,
        conf_level = 0.95,
        n_threads = 1L,
        output = out$output[[1]],
        threshold = out$threshold[[1]],
        diag = out$diag[[1]],
        verbose = FALSE
      )
      dispatch_expectations(obj, lbl, out$output[[1]])
    }

    lbl <- sprintf(
      "shrinkage thr=%.2f diag=%s out=%s",
      out$threshold[[1]], out$diag[[1]], out$output[[1]]
    )
    obj <- shrinkage_corr(
      X,
      n_threads = 1L,
      output = out$output[[1]],
      threshold = out$threshold[[1]],
      diag = out$diag[[1]]
    )
    dispatch_expectations(obj, lbl, out$output[[1]])

    lbl <- sprintf(
      "schafer thr=%.2f diag=%s out=%s",
      out$threshold[[1]], out$diag[[1]], out$output[[1]]
    )
    obj <- schafer_corr(
      X,
      n_threads = 1L,
      output = out$output[[1]],
      threshold = out$threshold[[1]],
      diag = out$diag[[1]]
    )
    dispatch_expectations(obj, lbl, out$output[[1]])
  }

  Z <- matrix(rnorm(260 * 4), nrow = 260, ncol = 4)
  X_bin <- data.frame(
    b1 = Z[, 1] > -0.2,
    b2 = Z[, 2] > 0.1,
    b3 = Z[, 3] > 0.4,
    b4 = Z[, 4] > -0.5
  )
  X_bin_na <- X_bin
  X_bin_na[sample.int(nrow(X_bin_na), 10), 1] <- NA
  X_bin_na[sample.int(nrow(X_bin_na), 10), 2] <- NA

  X_ord <- data.frame(
    o1 = ordered(cut(Z[, 1], breaks = c(-Inf, -0.4, 0.3, Inf), labels = c("L", "M", "H"))),
    o2 = ordered(cut(Z[, 2], breaks = c(-Inf, -0.8, -0.1, 0.6, Inf), labels = c("1", "2", "3", "4"))),
    o3 = ordered(cut(Z[, 3], breaks = c(-Inf, -0.2, 0.5, Inf), labels = c("A", "B", "C"))),
    o4 = ordered(cut(Z[, 4], breaks = c(-Inf, -0.3, 0.2, Inf), labels = c("Q", "R", "S")))
  )
  X_ord_na <- X_ord
  X_ord_na[sample.int(nrow(X_ord_na), 10), 1] <- NA
  X_ord_na[sample.int(nrow(X_ord_na), 10), 2] <- NA

  lat_grid <- grid_df(
    na_method = c("error", "pairwise"),
    ci = c(FALSE, TRUE),
    p_value = c(FALSE, TRUE)
  )

  for (i in seq_len(nrow(lat_grid))) {
    row <- lat_grid[i, , drop = FALSE]
    for (j in seq_len(nrow(out_grid))) {
      out <- out_grid[j, , drop = FALSE]

      x_bin_use <- select_data_by_na(X_bin, X_bin_na, row$na_method[[1]])
      lbl <- sprintf(
        "tetra na=%s ci=%s p=%s thr=%.2f diag=%s out=%s",
        row$na_method[[1]], row$ci[[1]], row$p_value[[1]],
        out$threshold[[1]], out$diag[[1]], out$output[[1]]
      )
      obj <- tetrachoric(
        x_bin_use,
        na_method = row$na_method[[1]],
        ci = row$ci[[1]],
        p_value = row$p_value[[1]],
        conf_level = 0.95,
        correct = 0.5,
        output = out$output[[1]],
        threshold = out$threshold[[1]],
        diag = out$diag[[1]]
      )
      dispatch_expectations(obj, lbl, out$output[[1]])

      x_ord_use <- select_data_by_na(X_ord, X_ord_na, row$na_method[[1]])
      lbl <- sprintf(
        "poly na=%s ci=%s p=%s thr=%.2f diag=%s out=%s",
        row$na_method[[1]], row$ci[[1]], row$p_value[[1]],
        out$threshold[[1]], out$diag[[1]], out$output[[1]]
      )
      obj <- polychoric(
        x_ord_use,
        na_method = row$na_method[[1]],
        ci = row$ci[[1]],
        p_value = row$p_value[[1]],
        conf_level = 0.95,
        correct = 0.5,
        output = out$output[[1]],
        threshold = out$threshold[[1]],
        diag = out$diag[[1]]
      )
      dispatch_expectations(obj, lbl, out$output[[1]])
    }
  }
})

test_that("bruteforce dispatch: pcorr method combinations and constraints", {
  set.seed(20260423)
  X <- matrix(rnorm(220 * 5), nrow = 220, ncol = 5)
  colnames(X) <- paste0("P", seq_len(ncol(X)))

  out_grid <- grid_df(
    output = c("matrix", "sparse", "edge_list"),
    threshold = c(0, 0.2),
    diag = c(TRUE, FALSE)
  )
  out_grid <- out_grid[vapply(
    seq_len(nrow(out_grid)),
    function(i) valid_threshold(out_grid$output[[i]], out_grid$threshold[[i]]),
    logical(1)
  ), , drop = FALSE]

  pc_grid <- grid_df(
    method = c("sample", "oas", "ridge", "glasso"),
    ci = c(FALSE, TRUE),
    return_cov_precision = c(FALSE, TRUE),
    return_p_value = c(FALSE, TRUE)
  )

  for (i in seq_len(nrow(pc_grid))) {
    row <- pc_grid[i, , drop = FALSE]
    for (j in seq_len(nrow(out_grid))) {
      out <- out_grid[j, , drop = FALSE]

      is_non_matrix <- !identical(out$output[[1]], "matrix")
      invalid <- FALSE
      if (is_non_matrix && (row$ci[[1]] || row$return_cov_precision[[1]] || row$return_p_value[[1]])) {
        invalid <- TRUE
      }
      if (row$return_p_value[[1]] && !identical(row$method[[1]], "sample")) {
        invalid <- TRUE
      }
      if (row$ci[[1]] && !identical(row$method[[1]], "sample")) {
        invalid <- TRUE
      }

      lbl <- sprintf(
        "pcorr m=%s ci=%s cov=%s p=%s thr=%.2f diag=%s out=%s",
        row$method[[1]], row$ci[[1]], row$return_cov_precision[[1]], row$return_p_value[[1]],
        out$threshold[[1]], out$diag[[1]], out$output[[1]]
      )

      call <- function() {
        pcorr(
          X,
          method = row$method[[1]],
          ci = row$ci[[1]],
          conf_level = 0.95,
          return_cov_precision = row$return_cov_precision[[1]],
          return_p_value = row$return_p_value[[1]],
          lambda = 1e-2,
          output = out$output[[1]],
          threshold = out$threshold[[1]],
          diag = out$diag[[1]]
        )
      }

      if (invalid) {
        expect_error(call(), info = lbl)
      } else {
        obj <- call()
        dispatch_expectations(obj, lbl, out$output[[1]])
      }
    }
  }
})

test_that("invalid output-mode combinations fail with informative errors", {
  set.seed(20260424)
  X <- matrix(rnorm(120 * 4), nrow = 120, ncol = 4)
  colnames(X) <- paste0("E", seq_len(ncol(X)))

  x <- rnorm(80)
  y <- x + rnorm(80, sd = 0.5)

  tab2 <- table(sample(c("0", "1"), 80, TRUE), sample(c("0", "1"), 80, TRUE))

  expect_error(kendall_tau(x, y, output = "sparse"), "non-matrix output")
  expect_error(kendall_tau(x, y, ci = TRUE), "not available when")

  expect_error(tetrachoric(tab2, output = "sparse"), "must be .*matrix")
  expect_error(polychoric(tab2, output = "edge_list"), "must be .*matrix")

  expect_error(pcorr(X, method = "sample", output = "sparse", ci = TRUE), "point estimates only")
  expect_error(pcorr(X, method = "sample", output = "edge_list", return_p_value = TRUE), "point estimates only")
})
