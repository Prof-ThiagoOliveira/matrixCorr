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

  span <- .Machine$integer.max - 1L
  out <- ((as.integer(seed) - 1L + as.integer(offset)) %% span) + 1L
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
