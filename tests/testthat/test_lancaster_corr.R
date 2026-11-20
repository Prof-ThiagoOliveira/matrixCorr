lancaster_components_ref <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  transform <- if (method == "rank") {
    function(z) qnorm(rank(z) / (length(z) + 1))
  } else {
    function(z) as.numeric(scale(z))
  }
  a <- transform(x[keep])
  b <- transform(y[keep])
  suppressWarnings(c(cor(a, b), cor(a^2, b^2)))
}

test_that("rank Lancaster uses the published plotting positions at the tails", {
  # The squared-score component dominates. Midpoint scores would give 0.6165093.
  fit <- lancaster_corr(cbind(x = 1:8, y = c(7, 2, 4, 3, 6, 5, 1, 8)))
  expect_equal(fit[1, 2], 0.6629908843567711, tolerance = 1e-12)
})

test_that("Lancaster uses the two published components on the pair's rows", {
  x <- rep(c(-3, -2, -1, 1, 2, 3), each = 2)
  y <- x * rep(c(-1, 1), 6)
  data <- cbind(x, y, z = c(NA, 3:10, 1, 1, Inf))
  for (method in c("moment", "rank")) {
    fit <- lancaster_corr(data, method, "pairwise")
    d <- attr(fit, "diagnostics")
    for (j in 1:2) for (k in (j + 1):3) {
      ref <- lancaster_components_ref(data[, j], data[, k], method)
      expect_equal(c(d$linear[j,k], d$quadratic[j,k]), ref)
      expect_equal(fit[j,k], max(abs(ref)))
    }
    expect_equal(c(d$linear[1,2], d$quadratic[1,2]), c(0, 1))
    cc <- data[apply(is.finite(data), 1, all), ]
    common <- lancaster_corr(data, method, "complete")
    expect_equal(estimate(common), estimate(lancaster_corr(cc, method)))
    expect_true(all(summary(common)$n_complete == nrow(cc)))
    expect_equal(length(summary(common)$n_complete), 3L)
    expect_error(lancaster_corr(data, method, "error"))
    bad <- lancaster_corr(cbind(one = rep(1, 6), two = rep(c(-1, 1), 3),
                                ok = seq_len(6)), method)
    expect_true(all(is.na(estimate(bad)[1:2, ])))
  }
  expect_equal(lancaster_corr(cbind(x * 1e-240, y * 1e260), "moment")[1,2], 1)
  expect_true(is.na(lancaster_corr(cbind(1:2, 2:1))[1,2]))
})

test_that("Lancaster inference matches exhaustive small-sample calculations", {
  x <- c(-3, -.5, 1, 4)
  y <- c(2, -3, 1, 5)
  idx <- as.matrix(expand.grid(rep(list(1:4), 4)))
  perms <- idx[apply(idx, 1, function(z) !anyDuplicated(z)), ]
  for (method in c("moment", "rank")) {
    ref <- lancaster_components_ref(x, y, method)
    size <- max(abs(ref))
    p_exact <- mean(apply(perms, 1, function(i)
      max(abs(lancaster_components_ref(x, y[i], method))) >= size - 1e-12))
    # All 4^4 ordered paired resamples have the same probability.
    boot <- apply(idx, 1, function(i) lancaster_components_ref(x[i], y[i], method))
    # A balanced two-point sample has constant squared magnitudes. Floating
    # point qnorm can distinguish their squares by one ulp, so exclude these
    # undefined samples explicitly as well as constant margins.
    valid <- apply(idx, 1, function(i) !any(table(i) == 4) &&
                       !(length(unique(i)) == 2 && all(table(i) == 2)))
    b <- boot[which.max(abs(ref)), valid]
    se <- sqrt(mean((b - mean(b))^2))
    fit <- lancaster_corr(cbind(x,y), method, ci = TRUE, p_value = TRUE,
                          conf_level = .6, ci_method = "regular", n_boot = 20000, n_perm = 20000, seed = 317)
    expected <- pmin(1, pmax(0, size + c(-1, 1) * qnorm(.8) * se))
    expect_lt(max(abs(unlist(confint(fit)[1, c("lwr", "upr")]) - expected)), .012)
    # Check equation (18) using the exact bootstrap covariance and an
    # independent conditional-normal integral for the maximum quantile.
    covar <- cov(t(boot[, valid]))
    sds <- sqrt(diag(covar))
    rho <- cov2cor(covar)[1,2] * if (prod(ref) < 0) -1 else 1
    prob <- function(q) integrate(function(z)
      dnorm(z) * pnorm((q/sds[2] - rho*z)/sqrt(1-rho^2)),
      -Inf, q/sds[1])$value
    critical <- uniroot(function(q) prob(q) - .8, c(0, 5 * max(sds)))$root
    con <- lancaster_corr(cbind(x,y), method, ci = TRUE, conf_level = .6,
                          n_boot = 20000, seed = 317)
    expect_lt(abs(ci(con)$lwr.ci[1,2] - max(0, size - critical)), .012)
    expect_equal(ci(con)$upr.ci, ci(fit)$upr.ci)
    expect_lt(abs(tidy(fit)$p_value - p_exact), .015)
    expect_lt(abs(attr(fit, "diagnostics")$n_boot_success[1,2]/20000 - mean(valid)), .015)
  }
})

test_that("Lancaster uses shared result methods across representations", {
  set.seed(18)
  X <- matrix(rnorm(200), 50, 4, dimnames = list(NULL, letters[1:4]))
  X[,2] <- X[,1] + rnorm(50, sd = .1)
  dense <- lancaster_corr(X, ci = TRUE, p_value = TRUE, n_perm = 39,
                          n_boot = 49, seed = 27)
  expect_identical(class(dense), c("corr_matrix","lancaster_corr","corr_result","matrix"))
  expect_identical(coef(dense), estimate(dense))
  expect_identical(abs(dense), abs(estimate(dense)))
  expect_identical(dense * 0, estimate(dense) * 0)
  expect_output(print(dense), "Lancaster correlation matrix")
  expect_s3_class(base::summary(dense), "summary.corr_result")
  expect_s3_class(plot(dense), "ggplot")
  for (out in c("sparse", "edge_list")) {
    fit <- lancaster_corr(X, output = out, threshold = .4, diag = FALSE,
                          ci = TRUE, p_value = TRUE, n_perm = 39,
                          n_boot = 49, n_threads = 2, seed = 27)
    expected <- subset(tidy(dense), estimate >= .4)
    tab <- tidy(fit)
    expect_equal(tab$estimate, expected$estimate)
    expect_equal(tab$p_value, expected$p_value)
    expect_equal(confint(fit)$lwr, expected$lwr)
    expect_equal(confint(fit)$upr, expected$upr)
    expect_false(is.null(ci(fit)))
    expect_s3_class(base::summary(fit), "summary.corr_result")
    expect_s3_class(plot(fit), "ggplot")
  }
})
