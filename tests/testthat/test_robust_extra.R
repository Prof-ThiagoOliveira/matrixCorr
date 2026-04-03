pbcor_manual_R <- function(x, y, beta = 0.2) {
  stopifnot(length(x) == length(y), !anyNA(x), !anyNA(y))

  pb_scores <- function(z) {
    med <- stats::median(z)
    omega <- sort(abs(z - med))[floor((1 - beta) * length(z))]
    if (!is.finite(omega) || omega <= 0) return(rep(NA_real_, length(z)))
    psi <- (z - med) / omega
    i1 <- sum(psi < -1)
    i2 <- sum(psi > 1)
    s <- z
    s[psi < -1 | psi > 1] <- 0
    denom <- length(z) - i1 - i2
    if (denom <= 0) return(rep(NA_real_, length(z)))
    pbos <- (sum(s) + omega * (i2 - i1)) / denom
    a <- (z - pbos) / omega
    a[a <= -1] <- -1
    a[a >= 1] <- 1
    a
  }

  ax <- pb_scores(x)
  ay <- pb_scores(y)
  if (anyNA(ax) || anyNA(ay)) return(NA_real_)
  sum(ax * ay) / sqrt(sum(ax^2) * sum(ay^2))
}

pbcor_manual_matrix_R <- function(X, beta = 0.2) {
  X <- as.matrix(X)
  p <- ncol(X)
  out <- diag(1, nrow = p, ncol = p)
  colnames(out) <- rownames(out) <- colnames(X)
  for (j in seq_len(p - 1L)) {
    for (k in seq.int(j + 1L, p)) {
      val <- pbcor_manual_R(X[, j], X[, k], beta = beta)
      out[j, k] <- val
      out[k, j] <- val
    }
  }
  out
}

plain_matrix_R <- function(X) {
  X <- as.matrix(X)
  matrix(
    data = as.vector(X),
    nrow = nrow(X),
    ncol = ncol(X),
    dimnames = dimnames(X)
  )
}

wincor_manual_R <- function(x, y, tr = 0.2) {
  stopifnot(length(x) == length(y), !anyNA(x), !anyNA(y))

  winval <- function(z) {
    zs <- sort(z)
    g <- floor(tr * length(z))
    ibot <- g + 1L
    itop <- length(z) - g
    zbot <- zs[ibot]
    ztop <- zs[itop]
    z[z <= zbot] <- zbot
    z[z >= ztop] <- ztop
    z
  }

  stats::cor(winval(x), winval(y))
}

idealf_iqr_R <- function(x) {
  n <- length(x)
  if (n < 2) return(0)
  j <- floor(n / 4 + 5 / 12)
  y <- sort(x)
  g <- (n / 4) - j + 5 / 12
  low <- (1 - g) * y[j] + g * y[j + 1]
  k <- n - j + 1
  up <- (1 - g) * y[k] + g * y[k - 1]
  up - low
}

skipcor_manual_R <- function(x, y,
                             method = c("pearson", "spearman"),
                             stand = TRUE,
                             outlier_rule = c("idealf", "mad"),
                             cutoff = sqrt(stats::qchisq(0.975, df = 2))) {
  method <- match.arg(method)
  outlier_rule <- match.arg(outlier_rule)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 5) return(NA_real_)

  X <- cbind(x, y)
  if (stand) {
    for (j in 1:2) {
      med <- stats::median(X[, j])
      sc <- stats::mad(X[, j], constant = 1)
      if (!is.finite(sc) || sc <= 0) {
        qs <- stats::quantile(X[, j], c(0.25, 0.75), names = FALSE, type = 7)
        sc <- diff(qs) / 1.3489795003921634
      }
      if (!is.finite(sc) || sc <= 0) sc <- 1
      X[, j] <- (X[, j] - med) / sc
    }
  }

  center <- apply(X, 2, stats::median)
  B <- sweep(X, 2, center, "-")
  bot <- rowSums(B^2)
  outlier <- rep(FALSE, nrow(B))

  for (i in seq_len(nrow(B))) {
    if (!(bot[i] > 0)) next
    dis <- abs(drop(B %*% B[i, ])) / sqrt(bot[i])
    med <- stats::median(dis)
    spread <- if (outlier_rule == "mad") {
      stats::mad(dis, center = med, constant = 1)
    } else {
      idealf_iqr_R(dis)
    }
    if (!is.finite(spread) || spread <= 0) next
    thresh <- med + cutoff * spread
    outlier <- outlier | (dis > thresh)
  }

  if (sum(!outlier) < 5) return(NA_real_)
  if (method == "spearman") {
    stats::cor(x[!outlier], y[!outlier], method = "spearman")
  } else {
    stats::cor(x[!outlier], y[!outlier], method = "pearson")
  }
}

skipcor_manual_mask_R <- function(x, y,
                                  stand = TRUE,
                                  outlier_rule = c("idealf", "mad"),
                                  cutoff = sqrt(stats::qchisq(0.975, df = 2))) {
  outlier_rule <- match.arg(outlier_rule)
  ok <- is.finite(x) & is.finite(y)
  x_ok <- x[ok]
  y_ok <- y[ok]
  idx_ok <- which(ok)
  if (length(x_ok) < 5) return(integer())

  X <- cbind(x_ok, y_ok)
  if (stand) {
    for (j in 1:2) {
      med <- stats::median(X[, j])
      sc <- stats::mad(X[, j], constant = 1)
      if (!is.finite(sc) || sc <= 0) {
        qs <- stats::quantile(X[, j], c(0.25, 0.75), names = FALSE, type = 7)
        sc <- diff(qs) / 1.3489795003921634
      }
      if (!is.finite(sc) || sc <= 0) sc <- 1
      X[, j] <- (X[, j] - med) / sc
    }
  }

  center <- apply(X, 2, stats::median)
  B <- sweep(X, 2, center, "-")
  bot <- rowSums(B^2)
  outlier <- rep(FALSE, nrow(B))

  for (i in seq_len(nrow(B))) {
    if (!(bot[i] > 0)) next
    dis <- abs(drop(B %*% B[i, ])) / sqrt(bot[i])
    med <- stats::median(dis)
    spread <- if (outlier_rule == "mad") {
      stats::mad(dis, center = med, constant = 1)
    } else {
      idealf_iqr_R(dis)
    }
    if (!is.finite(spread) || spread <= 0) next
    thresh <- med + cutoff * spread
    outlier <- outlier | (dis > thresh)
  }

  idx_ok[outlier]
}

skipcor_manual_matrix_R <- function(X,
                                    method = c("pearson", "spearman"),
                                    stand = TRUE,
                                    outlier_rule = c("idealf", "mad"),
                                    cutoff = sqrt(stats::qchisq(0.975, df = 2))) {
  method <- match.arg(method)
  outlier_rule <- match.arg(outlier_rule)
  X <- as.matrix(X)
  p <- ncol(X)
  out <- diag(1, nrow = p, ncol = p)
  colnames(out) <- rownames(out) <- colnames(X)
  for (j in seq_len(p - 1L)) {
    for (k in seq.int(j + 1L, p)) {
      val <- skipcor_manual_R(
        X[, j],
        X[, k],
        method = method,
        stand = stand,
        outlier_rule = outlier_rule,
        cutoff = cutoff
      )
      out[j, k] <- val
      out[k, j] <- val
    }
  }
  out
}

test_that("pbcor matches manual percentage bend calculation", {
  set.seed(1001)
  x <- rnorm(60)
  y <- 0.7 * x + rnorm(60, sd = 0.4)
  x[1] <- 9
  y[1] <- -8

  R <- pbcor(cbind(x = x, y = y))
  expect_s3_class(R, "pbcor")
  expect_equal(R["x", "y"], pbcor_manual_R(x, y), tolerance = 1e-12)
})

test_that("pbcor compiled backend matches a fixed known example", {
  x <- c(1, 2, 3, 4, 5, 20, 7, 8, 9, 10)
  y <- c(1.2, 2.1, 2.8, 4.2, 5.1, -18, 6.9, 8.3, 8.8, 10.4)
  X <- cbind(x = x, y = y)

  R_cpp <- pbcor_matrix_cpp(X, beta = 0.2, n_threads = 1L)
  expect_equal(R_cpp[1, 2], 0.5901576900261898, tolerance = 1e-12)
  expect_equal(R_cpp[1, 2], pbcor_manual_R(x, y, beta = 0.2), tolerance = 1e-12)
})

test_that("pbcor matches the manual percentage bend matrix on wider input", {
  set.seed(1005)
  X <- matrix(rnorm(120 * 6), nrow = 120, ncol = 6)
  X[, 2] <- 0.75 * X[, 1] + rnorm(120, sd = 0.35)
  X[, 4] <- -0.60 * X[, 3] + rnorm(120, sd = 0.40)
  X[c(2, 11, 77), 5] <- X[c(2, 11, 77), 5] + c(12, -11, 9)
  X[c(4, 39), 6] <- X[c(4, 39), 6] - c(13, 10)
  colnames(X) <- paste0("V", seq_len(ncol(X)))

  R <- pbcor(X, beta = 0.2)
  expect_equal(
    unname(plain_matrix_R(R)),
    unname(plain_matrix_R(pbcor_manual_matrix_R(X, beta = 0.2))),
    tolerance = 1e-12
  )
})

test_that("pbcor compiled backend matches the manual matrix on wider input", {
  set.seed(1005)
  X <- matrix(rnorm(120 * 6), nrow = 120, ncol = 6)
  X[, 2] <- 0.75 * X[, 1] + rnorm(120, sd = 0.35)
  X[, 4] <- -0.60 * X[, 3] + rnorm(120, sd = 0.40)
  X[c(2, 11, 77), 5] <- X[c(2, 11, 77), 5] + c(12, -11, 9)
  X[c(4, 39), 6] <- X[c(4, 39), 6] - c(13, 10)
  colnames(X) <- paste0("V", seq_len(ncol(X)))

  expect_equal(
    unname(pbcor_matrix_cpp(X, beta = 0.2, n_threads = 1L)),
    unname(pbcor_manual_matrix_R(X, beta = 0.2)),
    tolerance = 1e-12
  )
})

test_that("wincor matches a fixed known example", {
  x <- c(1, 2, 3, 4, 5, 20, 7, 8, 9, 10)
  y <- c(1.2, 2.1, 2.8, 4.2, 5.1, -18, 6.9, 8.3, 8.8, 10.4)
  X <- cbind(x = x, y = y)

  R <- wincor(X, tr = 0.2)
  expect_equal(R[1, 2], 0.6887508025909920, tolerance = 1e-12)
  expect_equal(R[1, 2], wincor_manual_R(x, y, tr = 0.2), tolerance = 1e-12)
})

test_that("wincor matches manual winsorized Pearson calculation", {
  set.seed(1002)
  x <- rnorm(80)
  y <- 0.6 * x + rnorm(80, sd = 0.5)
  x[c(2, 5)] <- c(-10, 11)
  y[c(2, 5)] <- c(-9, 10)

  R <- wincor(cbind(x = x, y = y), tr = 0.2)
  expect_s3_class(R, "wincor")
  expect_equal(R["x", "y"], wincor_manual_R(x, y, tr = 0.2), tolerance = 1e-12)
})

test_that("pbcor and wincor pairwise NA mode use overlap only", {
  x <- c(-3, -2, -1, 0, 1, 2, 3, NA)
  y <- c(-2.8, -2, -1.2, 0.1, 1.1, 2.2, 3.1, 4)
  z <- c(1, NA, 2, 3, 4, 5, 6, 7)
  M <- cbind(x = x, y = y, z = z)

  Rp <- pbcor(M, na_method = "pairwise")
  Rw <- wincor(M, na_method = "pairwise")

  idx_xy <- is.finite(x) & is.finite(y)
  idx_xz <- is.finite(x) & is.finite(z)
  expect_equal(Rp["x", "y"], pbcor_manual_R(x[idx_xy], y[idx_xy]), tolerance = 1e-12)
  expect_equal(Rw["x", "z"], wincor_manual_R(x[idx_xz], z[idx_xz]), tolerance = 1e-12)
})

test_that("skipped_corr matches fixed known examples", {
  x <- c(1, 2, 3, 4, 5, 20, 7, 8, 9, 10)
  y <- c(1.2, 2.1, 2.8, 4.2, 5.1, -18, 6.9, 8.3, 8.8, 10.4)
  X <- cbind(x = x, y = y)

  Rp <- skipped_corr(X, method = "pearson")
  Rs <- skipped_corr(X, method = "spearman")

  expect_equal(Rp[1, 2], 0.9978288433383867, tolerance = 1e-12)
  expect_equal(Rs[1, 2], 1.0000000000000000, tolerance = 1e-12)
  expect_equal(Rp[1, 2], skipcor_manual_R(x, y, method = "pearson"), tolerance = 1e-12)
  expect_equal(Rs[1, 2], skipcor_manual_R(x, y, method = "spearman"), tolerance = 1e-12)
})

test_that("skipped_corr matches manual skipped Pearson and Spearman calculations", {
  set.seed(1003)
  x <- rnorm(70)
  y <- 0.8 * x + rnorm(70, sd = 0.3)
  x[1] <- 7
  y[1] <- -7

  Rp <- skipped_corr(cbind(x = x, y = y), method = "pearson")
  Rs <- skipped_corr(cbind(x = x, y = y), method = "spearman")

  expect_s3_class(Rp, "skipped_corr")
  expect_equal(Rp["x", "y"], skipcor_manual_R(x, y, method = "pearson"), tolerance = 1e-12)
  expect_equal(Rs["x", "y"], skipcor_manual_R(x, y, method = "spearman"), tolerance = 1e-12)
  expect_gt(Rp["x", "y"], stats::cor(x, y))
})

test_that("skipped_corr matches the manual skipped-correlation matrix on wider input", {
  set.seed(1006)
  X <- matrix(rnorm(140 * 5), nrow = 140, ncol = 5)
  X[, 2] <- 0.65 * X[, 1] + rnorm(140, sd = 0.35)
  X[, 4] <- -0.70 * X[, 3] + rnorm(140, sd = 0.40)
  X[c(5, 17), 1] <- c(9, -8)
  X[c(5, 17), 2] <- c(-9, 8)
  X[c(33, 91), 4] <- X[c(33, 91), 4] + c(10, -11)
  colnames(X) <- paste0("V", seq_len(ncol(X)))

  R <- skipped_corr(X, method = "pearson")
  expect_equal(
    unname(plain_matrix_R(R)),
    unname(plain_matrix_R(skipcor_manual_matrix_R(X, method = "pearson"))),
    tolerance = 1e-12
  )
})

test_that("skipped_corr default output is unchanged and has no extra mask payload", {
  set.seed(1007)
  X <- matrix(rnorm(120), nrow = 30, ncol = 4)
  X[1, 1] <- 8
  X[1, 2] <- -8
  colnames(X) <- paste0("V", seq_len(ncol(X)))

  R_default <- skipped_corr(X, method = "pearson")
  R_masks <- skipped_corr(X, method = "pearson", return_masks = TRUE)

  expect_equal(unname(plain_matrix_R(R_default)), unname(plain_matrix_R(R_masks)), tolerance = 0)
  expect_null(attr(R_default, "skipped_masks", exact = TRUE))
  expect_true(inherits(attr(R_masks, "skipped_masks", exact = TRUE), "skipped_corr_masks"))
})

test_that("skipped_corr legacy example remains unchanged when inference is off", {
  set.seed(4242)
  X <- matrix(rnorm(80), ncol = 4)
  X[1, 1] <- 8
  X[1, 2] <- -8
  colnames(X) <- paste0("V", 1:4)

  R <- skipped_corr(X, method = "pearson")
  expected <- matrix(
    c(
      1.000000000000000,  0.442551271825158,  0.001069283910523, -0.056140658927032,
      0.442551271825158,  1.000000000000000, -0.193019542135606,  0.180008663965879,
      0.001069283910523, -0.193019542135606,  1.000000000000000,  0.072274267976520,
     -0.056140658927032,  0.180008663965879,  0.072274267976520,  1.000000000000000
    ),
    nrow = 4, byrow = TRUE,
    dimnames = list(colnames(X), colnames(X))
  )

  expect_equal(unname(plain_matrix_R(R)), unname(expected), tolerance = 1e-12)
  expect_null(attr(R, "lwr.ci", exact = TRUE))
  expect_null(attr(R, "p_value", exact = TRUE))
})

test_that("skipped_corr bootstrap inference attaches CI and p-value matrices", {
  set.seed(1101)
  X <- matrix(rnorm(60 * 3), nrow = 60, ncol = 3)
  X[, 2] <- 0.75 * X[, 1] + rnorm(60, sd = 0.35)
  X[1, 1] <- 9
  X[1, 2] <- -9
  colnames(X) <- c("A", "B", "C")

  R <- skipped_corr(
    X,
    method = "pearson",
    ci = TRUE,
    p_value = TRUE,
    n_boot = 200,
    seed = 123,
    p_adjust = "hochberg"
  )

  ci <- R$ci
  lwr <- ci$lwr.ci
  upr <- ci$upr.ci
  pval <- R$p_value
  padj <- R$p_value_adjusted
  rej <- R$reject

  expect_true(is.list(ci))
  expect_true(is.matrix(lwr))
  expect_true(is.matrix(upr))
  expect_true(is.matrix(pval))
  expect_true(is.matrix(padj))
  expect_true(is.matrix(rej))
  expect_equal(dim(lwr), dim(R))
  expect_equal(dim(upr), dim(R))
  expect_equal(dim(pval), dim(R))
  expect_equal(dim(padj), dim(R))
  expect_equal(unname(diag(lwr)), rep(1, ncol(X)), tolerance = 0)
  expect_equal(unname(diag(upr)), rep(1, ncol(X)), tolerance = 0)
  expect_equal(unname(diag(pval)), rep(0, ncol(X)), tolerance = 0)
  expect_equal(unname(diag(padj)), rep(0, ncol(X)), tolerance = 0)
  expect_true(isSymmetric(lwr))
  expect_true(isSymmetric(upr))
  expect_true(isSymmetric(pval))
  expect_true(isSymmetric(padj))
  expect_true(R["A", "B"] >= lwr["A", "B"] && R["A", "B"] <= upr["A", "B"])
  expect_identical(R$inference$method, "method_h")
  expect_identical(R$inference$p_adjust, "hochberg")

  idx <- upper.tri(pval, diag = FALSE) & is.finite(pval)
  expect_equal(
    padj[idx],
    stats::p.adjust(pval[idx], method = "hochberg"),
    tolerance = 1e-12
  )
  expect_identical(
    unname(rej[idx]),
    unname(padj[idx] <= R$inference$fwe_level)
  )
})

test_that("skipped_corr ECP inference returns a critical p-value and reject matrix", {
  set.seed(1103)
  X <- matrix(rnorm(50 * 3), nrow = 50, ncol = 3)
  X[, 2] <- 0.7 * X[, 1] + rnorm(50, sd = 0.3)
  X[1, 1] <- 8
  X[1, 2] <- -8
  colnames(X) <- c("A", "B", "C")

  R <- skipped_corr(
    X,
    method = "pearson",
    p_value = TRUE,
    n_boot = 80,
    p_adjust = "ecp",
    n_mc = 40,
    seed = 777
  )

  pval <- R$p_value
  reject <- R$reject
  critical <- R$critical_p_value

  expect_true(is.matrix(pval))
  expect_true(is.matrix(reject))
  expect_true(is.finite(critical))
  expect_identical(R$inference$method, "ecp")
  expect_identical(R$inference$p_adjust, "ecp")
  expect_identical(R$inference$n_mc, 40L)
  expect_equal(unname(diag(reject)), rep(FALSE, ncol(X)))
  idx <- upper.tri(pval, diag = FALSE) & is.finite(pval)
  expect_identical(unname(reject[idx]), unname(pval[idx] <= critical))
})

test_that("skipped_corr list-like accessors expose CI and inference components", {
  set.seed(1104)
  X <- matrix(rnorm(40 * 3), nrow = 40, ncol = 3)
  X[1, 1] <- 7
  X[1, 2] <- -7
  colnames(X) <- c("A", "B", "C")

  R <- skipped_corr(X, ci = TRUE, p_value = TRUE, n_boot = 60, seed = 99)

  expect_true("ci" %in% names(R))
  expect_true("p_value" %in% names(R))
  expect_true(is.list(R$ci))
  expect_true(is.matrix(R$ci$lwr.ci))
  expect_true(is.matrix(R[["p_value"]]))
  expect_true(is.list(R$inference))
})

test_that("skipped_corr plot supports CI overlays when intervals are available", {
  set.seed(1105)
  X <- matrix(rnorm(45 * 3), nrow = 45, ncol = 3)
  X[1, 1] <- 8
  X[1, 2] <- -8
  colnames(X) <- c("A", "B", "C")

  R <- skipped_corr(X, ci = TRUE, n_boot = 60, seed = 123)
  p <- plot(R, value_text_size = 2, ci_text_size = 2)

  expect_s3_class(p, "ggplot")
})

test_that("skipped_corr summary surfaces pairwise inference details when available", {
  set.seed(1106)
  X <- matrix(rnorm(40 * 4), nrow = 40, ncol = 4)
  X[1, 1] <- 8
  X[1, 2] <- -8
  colnames(X) <- paste0("V", 1:4)

  R_ci <- skipped_corr(X, ci = TRUE, n_boot = 60, seed = 123)
  sm_ci <- summary(R_ci)
  txt_ci <- capture.output(print(sm_ci, digits = 4))

  expect_s3_class(sm_ci, "summary.skipped_corr")
  expect_equal(nrow(sm_ci), choose(ncol(X), 2))
  expect_true(isTRUE(attr(sm_ci, "has_ci")))
  expect_false(isTRUE(attr(sm_ci, "has_p")))
  expect_match(paste(txt_ci, collapse = "\n"), "Skipped correlation summary")
  expect_match(paste(txt_ci, collapse = "\n"), "ci_width")

  R_p <- skipped_corr(X, ci = TRUE, p_value = TRUE, n_boot = 60, seed = 456)
  sm_p <- summary(R_p)
  txt_p <- capture.output(print(sm_p, digits = 4))

  expect_s3_class(sm_p, "summary.skipped_corr")
  expect_true(isTRUE(attr(sm_p, "has_ci")))
  expect_true(isTRUE(attr(sm_p, "has_p")))
  expect_match(paste(txt_p, collapse = "\n"), "inference")
  expect_match(paste(txt_p, collapse = "\n"), "bootstrap")
})

test_that("skipped_corr masks reconstruct the reported pairwise correlations", {
  set.seed(1008)
  X <- matrix(rnorm(90 * 3), nrow = 90, ncol = 3)
  X[, 2] <- 0.7 * X[, 1] + rnorm(90, sd = 0.25)
  X[c(5, 11), 1] <- c(9, -10)
  X[c(5, 11), 2] <- c(-9, 10)
  X[c(13, 24), 3] <- c(11, -12)
  X[7, 3] <- NA_real_
  colnames(X) <- c("A", "B", "C")

  R <- skipped_corr(X, method = "spearman", na_method = "pairwise", return_masks = TRUE)

  for (j in 1:(ncol(X) - 1L)) {
    for (k in (j + 1L):ncol(X)) {
      skipped <- skipped_corr_masks(R, j, k)
      keep <- is.finite(X[, j]) & is.finite(X[, k])
      if (length(skipped)) keep[skipped] <- FALSE
      expected <- if (sum(keep) < 5L) NA_real_ else stats::cor(X[keep, j], X[keep, k], method = "spearman")
      expect_equal(unname(R[j, k]), expected, tolerance = 1e-12)
    }
  }
})

test_that("skipped_corr mask accessor is symmetric by variable pair", {
  set.seed(1009)
  X <- matrix(rnorm(70 * 3), nrow = 70, ncol = 3)
  X[1, 1] <- 8
  X[1, 2] <- -8
  colnames(X) <- c("x", "y", "z")

  R <- skipped_corr(X, return_masks = TRUE)
  expect_identical(skipped_corr_masks(R, "x", "y"), skipped_corr_masks(R, "y", "x"))
})

test_that("skipped_corr two-variable masks expose expected skipped rows", {
  x <- c(1, 2, 3, 4, 5, 20, 7, 8, 9, 10)
  y <- c(1.2, 2.1, 2.8, 4.2, 5.1, -18, 6.9, 8.3, 8.8, 10.4)
  X <- cbind(x = x, y = y)

  R <- skipped_corr(X, method = "pearson", return_masks = TRUE)
  expect_identical(skipped_corr_masks(R, "x", "y"), 6L)
})

test_that("skipped_corr masks match manual skipped rows on pairwise NA input", {
  x <- c(-3, -2, -1, 0, 1, 2, 8, 3, NA)
  y <- c(-3.1, -2.2, -1.0, 0.2, 1.1, 2.0, -7.5, 2.9, 4)
  z <- c(1, 2, 3, 4, 5, 6, 7, NA, 9)
  M <- cbind(x = x, y = y, z = z)

  R <- skipped_corr(M, na_method = "pairwise", return_masks = TRUE)
  expect_identical(
    skipped_corr_masks(R, "x", "y"),
    skipcor_manual_mask_R(x, y)
  )
  expect_identical(
    skipped_corr_masks(R, "x", "z"),
    skipcor_manual_mask_R(x, z)
  )
})

test_that("new robust correlation classes support print and plot methods", {
  set.seed(1004)
  X <- matrix(rnorm(120), nrow = 30, ncol = 4)
  colnames(X) <- paste0("V", 1:4)

  objs <- list(
    pbcor(X),
    wincor(X),
    skipped_corr(X)
  )

  for (obj in objs) {
    txt <- capture.output(print(obj, digits = 2))
    expect_true(length(txt) > 0)
    p <- plot(obj, value_text_size = 2)
    expect_s3_class(p, "ggplot")
    sm <- summary(obj)
    expect_s3_class(sm, "summary_corr_matrix")
  }
})

test_that("skipped_corr summary reports skipped counts and proportions", {
  set.seed(1102)
  X <- matrix(rnorm(50 * 4), nrow = 50, ncol = 4)
  X[1, 1] <- 8
  X[1, 2] <- -8
  colnames(X) <- paste0("V", 1:4)

  sm <- summary(skipped_corr(X))
  txt <- capture.output(print(sm, digits = 4))

  expect_match(paste(txt, collapse = "\n"), "skipped_n")
  expect_match(paste(txt, collapse = "\n"), "skipped_prop")
  expect_true(is.finite(sm$skipped_n_min))
  expect_true(is.finite(sm$skipped_n_max))
  expect_true(is.finite(sm$skipped_prop_min))
  expect_true(is.finite(sm$skipped_prop_max))
})
