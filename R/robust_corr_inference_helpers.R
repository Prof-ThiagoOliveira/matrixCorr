.mc_corr_n_complete <- function(x) {
  x <- as.matrix(x)
  fin <- is.finite(x)
  n_complete <- crossprod(unclass(fin))
  storage.mode(n_complete) <- "integer"
  n_complete
}

.mc_eval_with_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(force(expr))
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  force(expr)
}

.mc_seed_offset <- function(seed, offset) {
  if (is.null(seed)) {
    return(NULL)
  }

  span <- .Machine$integer.max - 1
  # Reduce each term before adding. The sum of the two residues stays below
  # two spans, so that addition and the final remainder stay exact in double
  # precision whenever both truncated terms are exactly representable integers.
  seed_mod <- (trunc(as.numeric(seed)) - 1) %% span
  offset_mod <- trunc(as.numeric(offset)) %% span
  out <- ((seed_mod + offset_mod) %% span) + 1
  as.integer(out)
}

.mc_percentile_boot_ci <- function(values, conf_level = 0.95) {
  values <- sort(values[is.finite(values)])
  n_boot <- length(values)
  if (!n_boot) {
    return(c(NA_real_, NA_real_))
  }

  alpha <- 1 - conf_level
  ilow <- floor((alpha / 2) * n_boot + 0.5)
  ihi <- floor((1 - alpha / 2) * n_boot + 0.5)
  ilow <- min(max(ilow, 1L), n_boot)
  ihi <- min(max(ihi, 1L), n_boot)
  c(values[[ilow]], values[[ihi]])
}

.mc_robust_corr_pairwise_payload <- function(X,
                                             est,
                                             tune,
                                             ci = FALSE,
                                             p_value = FALSE,
                                             conf_level = 0.95,
                                             n_boot = 500L,
                                             seed = NULL,
                                             pair_inference,
                                             method) {
  est <- as.matrix(est)
  p <- ncol(est)
  dn <- dimnames(est)
  n_complete <- .mc_corr_n_complete(X)
  dimnames(n_complete) <- dn

  ci_lwr <- ci_upr <- NULL
  if (isTRUE(ci)) {
    ci_lwr <- matrix(NA_real_, p, p, dimnames = dn)
    ci_upr <- matrix(NA_real_, p, p, dimnames = dn)
    diag(ci_lwr) <- 1
    diag(ci_upr) <- 1
  }

  stat_mat <- p_mat <- df_mat <- NULL
  if (isTRUE(p_value)) {
    stat_mat <- matrix(NA_real_, p, p, dimnames = dn)
    p_mat <- matrix(NA_real_, p, p, dimnames = dn)
    df_mat <- matrix(NA_real_, p, p, dimnames = dn)
    diag(stat_mat) <- Inf
    diag(p_mat) <- 0
  }

  pair_id <- 0
  for (i in seq_len(p - 1L)) {
    for (j in seq.int(i + 1L, p)) {
      pair_id <- pair_id + 1
      ok <- is.finite(X[, i]) & is.finite(X[, j])
      n_ij <- sum(ok)
      n_complete[i, j] <- n_complete[j, i] <- n_ij
      if (n_ij < 5L) {
        next
      }

      info <- pair_inference(
        X[ok, i],
        X[ok, j],
        tune,
        ci = ci,
        p_value = p_value,
        conf_level = conf_level,
        n_boot = n_boot,
        seed = .mc_seed_offset(seed, pair_id - 1L)
      )

      if (isTRUE(ci)) {
        ci_lwr[i, j] <- ci_lwr[j, i] <- info$lwr
        ci_upr[i, j] <- ci_upr[j, i] <- info$upr
      }
      if (isTRUE(p_value)) {
        stat_mat[i, j] <- stat_mat[j, i] <- info$statistic
        p_mat[i, j] <- p_mat[j, i] <- info$p_value
        df_mat[i, j] <- df_mat[j, i] <- info$parameter
      }
    }
  }

  diagnostics <- list(n_complete = n_complete)
  ci_attr <- NULL
  if (isTRUE(ci)) {
    ci_attr <- list(
      est = unclass(est),
      lwr.ci = ci_lwr,
      upr.ci = ci_upr,
      conf.level = conf_level,
      ci.method = "percentile_bootstrap"
    )
  }

  inference_attr <- NULL
  if (isTRUE(p_value)) {
    inference_attr <- list(
      method = method,
      estimate = unclass(est),
      statistic = stat_mat,
      p_value = p_mat,
      parameter = df_mat,
      n_obs = n_complete,
      alternative = "two.sided"
    )
  }

  list(
    diagnostics = diagnostics,
    ci = ci_attr,
    inference = inference_attr
  )
}
