# Modular reference written independently of .mc_seed_offset:
# truncate toward zero, then a positive residue via floor division.
reference_seed_offset <- function(seed, offset) {
  span <- .Machine$integer.max - 1
  shifted <- trunc(as.numeric(seed)) - 1 + trunc(as.numeric(offset))
  residue <- shifted - span * floor(shifted / span)
  as.integer(residue + 1)
}

legacy_seed_offset <- function(seed, offset) {
  span <- .Machine$integer.max - 1L
  ((as.integer(seed) - 1L + as.integer(offset)) %% span) + 1L
}

legacy_seed_offset_defined <- function(seed, offset) {
  warned <- FALSE
  out <- withCallingHandlers(
    legacy_seed_offset(seed, offset),
    warning = function(w) {
      warned <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  !warned && length(out) == 1L && !is.na(out)
}

test_that("seed offsets match an independent modular reference", {
  span <- .Machine$integer.max - 1
  cases <- list(
    c(0, 0),
    c(0, 1),
    c(1, 0),
    c(1, 1),
    c(2.9, 0.2),
    c(-3.8, 4),
    c(11, -4.7),
    c(0, -1),
    c(span, 0),
    c(span, 1),
    c(1, span),
    c(.Machine$integer.max, 0),
    c(.Machine$integer.max, 1),
    c(1, .Machine$integer.max),
    c(.Machine$integer.max, .Machine$integer.max),
    c(0, .Machine$integer.max + 5),
    c(7, .Machine$integer.max + 5),
    c(3, 2 * .Machine$integer.max + 9),
    c(span, span + 4),
    c(-100, .Machine$integer.max * 3 + 2)
  )

  for (case in cases) {
    got <- matrixCorr:::.mc_seed_offset(case[[1L]], case[[2L]])
    expect_identical(got, reference_seed_offset(case[[1L]], case[[2L]]))
    expect_type(got, "integer")
    expect_length(got, 1L)
    expect_true(got >= 1L && got <= span)
  }
})

test_that("non-overflowing integer expression is preserved, including seed 0", {
  span <- .Machine$integer.max - 1L
  cases <- list(
    c(0, 0), c(0, 1), c(0, 8), c(1, 0), c(1, span - 1),
    c(4.9, 2.1), c(-2.2, 5), c(span, 0), c(10, -3),
    c(.Machine$integer.max, 0), c(100, 250)
  )
  for (case in cases) {
    expect_true(legacy_seed_offset_defined(case[[1L]], case[[2L]]))
    expect_identical(
      matrixCorr:::.mc_seed_offset(case[[1L]], case[[2L]]),
      legacy_seed_offset(case[[1L]], case[[2L]])
    )
  }
  expect_null(matrixCorr:::.mc_seed_offset(NULL, 4L))
  expect_null(matrixCorr:::.mc_seed_offset(NULL, .Machine$integer.max + 5))
})

test_that("fixed seeds reproduce percentage-bend and Winsorized bootstraps and restore RNG", {
  set.seed(21)
  X <- matrix(rnorm(48), nrow = 12, ncol = 4)
  colnames(X) <- paste0("v", seq_len(4))

  before <- get(".Random.seed", envir = .GlobalEnv)
  pb1 <- pbcor(X, ci = TRUE, n_boot = 20L, seed = 7L)
  after_pb <- get(".Random.seed", envir = .GlobalEnv)
  expect_identical(after_pb, before)
  pb2 <- pbcor(X, ci = TRUE, n_boot = 20L, seed = 7L)
  expect_identical(get(".Random.seed", envir = .GlobalEnv), before)
  expect_identical(as.matrix(pb1), as.matrix(pb2))
  expect_identical(attr(pb1, "ci"), attr(pb2, "ci"))

  win1 <- wincor(X, tr = 0.2, ci = TRUE, n_boot = 20L, seed = 11L)
  expect_identical(get(".Random.seed", envir = .GlobalEnv), before)
  win2 <- wincor(X, tr = 0.2, ci = TRUE, n_boot = 20L, seed = 11L)
  expect_identical(get(".Random.seed", envir = .GlobalEnv), before)
  expect_identical(as.matrix(win1), as.matrix(win2))
  expect_identical(attr(win1, "ci"), attr(win2, "ci"))

  set.seed(99)
  consumed_before <- get(".Random.seed", envir = .GlobalEnv)
  pbcor(X, ci = TRUE, n_boot = 15L, seed = NULL)
  consumed_after <- get(".Random.seed", envir = .GlobalEnv)
  expect_false(identical(consumed_before, consumed_after))
})

test_that("offsets at 2^53 - 1 match exact integer residues", {
  # Literals from Python arbitrary precision integers, span = 2^31 - 2.
  # .mc_seed_offset claims this only for exactly representable integer offsets.
  offset <- 2^53 - 1
  expected <- c(
    "0" = 8388607L,
    "1" = 8388608L,
    "2" = 8388609L,
    "-1" = 8388606L,
    "-100" = 8388507L
  )
  expected[[as.character(.Machine$integer.max - 1L)]] <- 8388607L
  expected[[as.character(.Machine$integer.max - 2L)]] <- 8388606L
  expected[[as.character(.Machine$integer.max)]] <- 8388608L

  for (nm in names(expected)) {
    got <- matrixCorr:::.mc_seed_offset(as.numeric(nm), offset)
    expect_identical(got, expected[[nm]], info = nm)
    expect_type(got, "integer")
  }
})

test_that("large seeds reproduce percentage-bend and Winsorized bootstraps and restore RNG", {
  set.seed(21)
  X <- matrix(rnorm(48), nrow = 12, ncol = 4)
  colnames(X) <- paste0("v", seq_len(4))
  large_seed <- .Machine$integer.max

  before <- get(".Random.seed", envir = .GlobalEnv)
  pb1 <- pbcor(X, ci = TRUE, n_boot = 20L, seed = large_seed)
  expect_identical(get(".Random.seed", envir = .GlobalEnv), before)
  pb2 <- pbcor(X, ci = TRUE, n_boot = 20L, seed = large_seed)
  expect_identical(get(".Random.seed", envir = .GlobalEnv), before)
  expect_identical(as.matrix(pb1), as.matrix(pb2))
  expect_identical(attr(pb1, "ci"), attr(pb2, "ci"))

  win1 <- wincor(X, tr = 0.2, ci = TRUE, n_boot = 20L, seed = large_seed - 1L)
  expect_identical(get(".Random.seed", envir = .GlobalEnv), before)
  win2 <- wincor(X, tr = 0.2, ci = TRUE, n_boot = 20L, seed = large_seed - 1L)
  expect_identical(get(".Random.seed", envir = .GlobalEnv), before)
  expect_identical(as.matrix(win1), as.matrix(win2))
  expect_identical(attr(win1, "ci"), attr(win2, "ci"))
})
