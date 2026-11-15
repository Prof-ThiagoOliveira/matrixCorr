dense_reference <- function(edge) {
  mat <- matrixCorr:::.mc_corr_as_dense_matrix(edge)
  matrixCorr:::.mc_new_corr_matrix(
    mat,
    estimator_class = attr(edge, "corr_estimator_class", exact = TRUE),
    method = attr(edge, "method", exact = TRUE),
    description = "dense reference",
    symmetric = attr(edge, "corr_symmetric", exact = TRUE),
    diagnostics = attr(edge, "diagnostics", exact = TRUE),
    ci = attr(edge, "ci", exact = TRUE),
    conf.level = attr(edge, "conf.level", exact = TRUE)
  )
}

make_edge <- function(i, j, value, dm, dn = list(NULL, NULL), symmetric = TRUE, diagnostics = NULL, ci = NULL) {
  row <- if (is.null(dn[[1L]])) as.character(i) else dn[[1L]][i]
  col <- if (is.null(dn[[2L]])) as.character(j) else dn[[2L]][j]
  df <- data.frame(row = row, col = col, value = value, stringsAsFactors = FALSE)
  edge <- matrixCorr:::.mc_new_corr_edge_list(
    df,
    estimator_class = "pearson_corr",
    method = "pearson",
    description = "edge fixture",
    threshold = 0,
    diag = TRUE,
    diagnostics = diagnostics,
    ci = ci,
    source_dim = dm,
    source_dimnames = dn,
    symmetric = symmetric
  )
  if (!is.null(dn[[1L]]) && (anyDuplicated(dn[[1L]]) > 0L || anyDuplicated(dn[[2L]]) > 0L)) {
    pos <- cbind(i = as.integer(i), j = as.integer(j))
    rownames(pos) <- rownames(df)
    attr(edge, "corr_pair_pos") <- pos
  }
  edge
}

expect_edge_matches_dense <- function(edge, ...) {
  expect_identical(tidy(edge, ...), tidy(dense_reference(edge), ...))
}

test_that("edge tidy matches the dense converter for triangles and direction", {
  set.seed(801)
  X <- cbind(a = rnorm(18), b = rnorm(18), c = rnorm(18))
  dense <- pearson_corr(X, ci = TRUE, conf_level = 0.9)
  edge <- pearson_corr(X, ci = TRUE, conf_level = 0.9, output = "edge_list", threshold = 0, diag = TRUE)
  expect_identical(tidy(edge), tidy(dense))
  for (triangle in c("upper", "lower", "full")) {
    for (include_diag in c(FALSE, TRUE)) {
      expect_identical(
        tidy(edge, triangle = triangle, diag = include_diag),
        tidy(dense, triangle = triangle, diag = include_diag),
        info = paste(triangle, include_diag)
      )
    }
  }

  x <- seq(-1, 1, length.out = 24)
  directed_x <- cbind(x = x, x2 = x^2)
  directed_dense <- xi_corr(directed_x, tie_method = "first", ci = FALSE)
  directed_edge <- xi_corr(
    directed_x,
    tie_method = "first",
    ci = FALSE,
    output = "edge_list",
    threshold = 0,
    diag = FALSE
  )
  expect_identical(tidy(directed_edge), tidy(directed_dense))
  expect_identical(
    tidy(directed_edge, triangle = "lower"),
    tidy(directed_dense, triangle = "lower")
  )
})

test_that("native edge coordinates keep last finite value and explicit zeros", {
  dn <- list(c("a", "b", "c"), c("a", "b", "c"))
  edge <- make_edge(
    i = c(1L, 2L, 1L, 1L, 3L, 2L),
    j = c(2L, 1L, 2L, 1L, 3L, 3L),
    value = c(0.2, 0.4, NA_real_, -0, Inf, 0),
    dm = c(3L, 3L),
    dn = dn,
    symmetric = TRUE
  )
  expect_edge_matches_dense(edge)
  expect_edge_matches_dense(edge, triangle = "lower", diag = TRUE)
  upper <- tidy(edge)
  expect_true(any(upper$estimate == 0))
  expect_false(any(is.na(upper$estimate)))

  labels <- list(c("a", "a"), c("a", "a"))
  duplicated <- make_edge(
    i = c(1L, 2L),
    j = c(2L, 1L),
    value = c(0.3, 0.8),
    dm = c(2L, 2L),
    dn = labels,
    symmetric = TRUE
  )
  expect_edge_matches_dense(duplicated, triangle = "full", diag = FALSE)

  unnamed <- make_edge(
    i = c(2L, 1L),
    j = c(1L, 3L),
    value = c(-0, 0.5),
    dm = c(3L, 4L),
    dn = list(NULL, NULL),
    symmetric = FALSE
  )
  td <- tidy(unnamed, triangle = "full", diag = TRUE)
  expect_identical(td$item1, c("2", "1"))
  expect_identical(td$item2, c("1", "3"))
  expect_true(identical(
    td$estimate,
    c(-0, 0.5),
    num.eq = FALSE,
    single.NA = FALSE,
    attrib.as.set = FALSE
  ))

  signed <- make_edge(
    i = c(1L, 1L),
    j = c(1L, 1L),
    value = c(0, -0),
    dm = c(2L, 2L),
    dn = list(c("a", "b"), c("a", "b")),
    symmetric = TRUE
  )
  signed_tidy <- tidy(signed, triangle = "upper", diag = TRUE)
  expect_true(identical(
    signed_tidy$estimate,
    -0,
    num.eq = FALSE,
    single.NA = FALSE,
    attrib.as.set = FALSE
  ))
  erased <- make_edge(
    i = c(1L, 1L, 1L),
    j = c(1L, 1L, 2L),
    value = c(0.4, NA_real_, 0.2),
    dm = c(2L, 2L),
    dn = list(c("a", "b"), c("a", "b")),
    symmetric = TRUE
  )
  expect_identical(tidy(erased)$estimate, 0.2)
  expect_false(any(tidy(erased)$item1 == tidy(erased)$item2))

  bad_names <- make_edge(1L, 2L, 0.3, c(2L, 2L), list(c("a", "b"), c("a", "b")), symmetric = TRUE)
  attr(bad_names, "corr_dimnames") <- list("only", c("a", "b"))
  expect_error(tidy(bad_names), "dimnames")
  expect_error(matrixCorr:::.mc_corr_as_dense_matrix(bad_names), "dimnames")
})

test_that("invalid edge coordinates use the dense subscript path", {
  bad <- make_edge(1L, 0L, 0.2, c(2L, 2L), symmetric = FALSE)
  expect_identical(tidy(bad), tidy(dense_reference(bad)))

  outside <- make_edge(3L, 1L, 0.2, c(2L, 2L), symmetric = FALSE)
  expect_error(tidy(outside), "subscript")
  expect_error(matrixCorr:::.mc_corr_as_dense_matrix(outside), "subscript")

  negative <- make_edge(-1L, 1L, 0.2, c(2L, 2L), symmetric = FALSE)
  expect_error(tidy(negative), "negative")

  rectangular <- make_edge(1L, 3L, 0.2, c(2L, 3L), symmetric = TRUE)
  expect_error(tidy(rectangular), "subscript")
  expect_error(matrixCorr:::.mc_corr_as_dense_matrix(rectangular), "subscript")

  missing_pos <- data.frame(row = c("a", "a"), col = c("a", "a"), value = c(0.1, 0.2))
  missing_pos <- matrixCorr:::.mc_new_corr_edge_list(
    missing_pos,
    estimator_class = "pearson_corr",
    method = "pearson",
    description = "ambiguous",
    threshold = 0,
    diag = TRUE,
    source_dim = c(2L, 2L),
    source_dimnames = list(c("a", "a"), c("a", "a")),
    symmetric = TRUE
  )
  expect_error(tidy(missing_pos), "positional pair mapping")

  unknown <- make_edge(1L, 2L, 0.4, c(2L, 2L), list(c("a", "b"), c("a", "b")), symmetric = FALSE)
  unknown$row <- "missing"
  expect_identical(ncol(tidy(unknown)), 0L)

  nonfinite <- make_edge(1L, 2L, 0.4, c(2L, 2L), symmetric = FALSE)
  nonfinite$i <- NA_real_
  attr(nonfinite, "corr_pair_pos") <- NULL
  nonfinite$row <- NA_character_
  expect_identical(ncol(tidy(nonfinite)), 0L)
})

test_that("edge confidence limits use the stored opposite cell", {
  dn <- list(c("a", "b"), c("a", "b"))
  lwr <- matrix(c(NA, 0.1, NA, NA), 2, 2, dimnames = dn)
  upr <- lwr + 0.2
  edge <- make_edge(
    i = 1L,
    j = 2L,
    value = 0.4,
    dm = c(2L, 2L),
    dn = dn,
    symmetric = TRUE,
    ci = list(lwr.ci = lwr, upr.ci = upr, conf.level = 0.95)
  )
  expect_edge_matches_dense(edge)
  expect_identical(tidy(edge)$lwr, 0.1)
})

test_that("empty edge and dense bootstrap diagnostics do not consume RNG", {
  x <- seq(-1, 1, length.out = 16)
  fit <- xi_corr(
    cbind(x = x, x2 = x^2),
    tie_method = "first",
    ci = TRUE,
    bootstrap_reps = 9L,
    seed = 3L,
    output = "edge_list",
    threshold = 2,
    diag = FALSE
  )
  seed <- .Random.seed
  td <- tidy(fit)
  sm <- summary(fit)
  expect_identical(.Random.seed, seed)
  expect_identical(nrow(td), 0L)
  expect_identical(ncol(td), 0L)
  expect_identical(nrow(sm), 0L)
  expect_error(confint(fit), "does not contain confidence limits")

  mat <- matrix(NA_real_, 2, 2, dimnames = list(c("a", "b"), c("a", "b")))
  dense <- matrixCorr:::.mc_new_corr_matrix(
    mat,
    estimator_class = "chatterjee_xi",
    method = "xi",
    description = "empty dense",
    symmetric = TRUE,
    diagnostics = list(bootstrap_reps = 9L)
  )
  expect_identical(nrow(tidy(dense)), 0L)
  expect_identical(nrow(summary(dense)), 0L)
})

test_that("skipped overview range ignores the retained edge filter", {
  dn <- list(c("a", "b", "c"), c("a", "b", "c"))
  skipped_n <- matrix(NA_real_, 3, 3, dimnames = dn)
  skipped_n[1, 2] <- 2
  skipped_n[1, 3] <- 99
  edge <- make_edge(
    i = 1L,
    j = 2L,
    value = 0.4,
    dm = c(3L, 3L),
    dn = dn,
    symmetric = TRUE,
    diagnostics = list(skipped_n = skipped_n)
  )
  overview <- matrixCorr:::.mc_corr_overview(edge)
  expect_identical(overview$n_pairs_retained, 1L)
  expect_identical(overview$skipped_n_min, 2)
  expect_identical(overview$skipped_n_max, 99)
  expect_identical(tidy(edge)$skipped_n, 2L)
})

test_that("large edge tidy does not allocate the estimate matrix", {
  edge <- make_edge(
    i = c(1L, 4L),
    j = c(2L, 8L),
    value = c(0, 0.25),
    dm = c(50000L, 50000L),
    dn = list(NULL, NULL),
    symmetric = FALSE
  )
  ns <- asNamespace("matrixCorr")
  trace(
    ".mc_corr_as_dense_matrix",
    tracer = quote(stop("edge densified", call. = FALSE)),
    print = FALSE,
    where = ns
  )
  on.exit(untrace(".mc_corr_as_dense_matrix", where = ns), add = TRUE)
  td <- tidy(edge, triangle = "full", diag = FALSE)
  expect_identical(td$item1, c("1", "4"))
  expect_identical(td$item2, c("2", "8"))
  expect_identical(td$estimate, c(0, 0.25))
})

test_that("thin edge projection keeps the last finite landing cell", {
  edge <- make_edge(
    i = c(1L, 2L, 1L, 1L, 3L, 2L),
    j = c(2L, 1L, 2L, 1L, 3L, 3L),
    value = c(0.2, 0.4, NA_real_, -0, Inf, 0),
    dm = c(100L, 100L),
    dn = list(NULL, NULL),
    symmetric = TRUE
  )
  ns <- asNamespace("matrixCorr")
  trace(
    ".mc_corr_as_dense_matrix",
    tracer = quote(stop("edge densified", call. = FALSE)),
    print = FALSE,
    where = ns
  )
  on.exit(untrace(".mc_corr_as_dense_matrix", where = ns), add = TRUE)
  upper <- tidy(edge)
  expect_identical(upper$item1, c("1", "2"))
  expect_identical(upper$item2, c("2", "3"))
  expect_true(identical(
    upper$estimate,
    c(0.4, 0),
    num.eq = FALSE,
    single.NA = FALSE,
    attrib.as.set = FALSE
  ))
  lower <- tidy(edge, triangle = "lower", diag = TRUE)
  expect_identical(lower$item1, c("1", "3"))
  expect_identical(lower$item2, c("1", "2"))
  expect_true(identical(
    lower$estimate,
    c(-0, 0),
    num.eq = FALSE,
    single.NA = FALSE,
    attrib.as.set = FALSE
  ))
  expect_false(any(lower$item1 == "2" & lower$item2 == "1"))
})

test_that("finite dense tidy keeps values and data-frame attribute order", {
  dn <- list(letters[1:3], letters[1:3])
  mat <- matrix(c(1, 0.2, NA, 0.2, 1, -0.3, NA, -0.3, 1), 3, 3, dimnames = dn)
  n_complete <- matrix(c(8, 8, 4, 8, 8, 8, 4, 4, 8), 3, 3, dimnames = dn)
  lwr <- mat - 0.1
  upr <- mat + 0.1
  mixed <- matrixCorr:::.mc_new_corr_matrix(
    mat,
    estimator_class = "pearson_corr",
    method = "pearson",
    description = "mixed finite",
    symmetric = TRUE,
    ci = list(lwr.ci = lwr, upr.ci = upr, conf.level = 0.95),
    diagnostics = list(n_complete = n_complete)
  )
  td <- tidy(mixed)
  expect_identical(td$item1, c("a", "b"))
  expect_identical(td$item2, c("b", "c"))
  expect_identical(td$estimate, c(0.2, -0.3))
  expect_identical(td$lwr, c(lwr[1, 2], lwr[2, 3]))
  expect_identical(td$upr, c(upr[1, 2], upr[2, 3]))
  expect_identical(td$n_complete, c(8L, 4L))
  expect_identical(names(attributes(td)), c("names", "row.names", "class"))

  finite <- mat
  finite[is.na(finite)] <- 0.1
  finite_obj <- matrixCorr:::.mc_new_corr_matrix(
    finite,
    estimator_class = "pearson_corr",
    method = "pearson",
    description = "all finite",
    symmetric = TRUE
  )
  finite_td <- tidy(finite_obj)
  expect_identical(names(attributes(finite_td)), c("names", "row.names", "class"))
  expect_identical(nrow(finite_td), 3L)

  signed <- finite
  signed[1, 2] <- -0
  signed[2, 1] <- -0
  signed[2, 3] <- NaN
  signed[3, 2] <- Inf
  signed_obj <- matrixCorr:::.mc_new_corr_matrix(
    signed,
    estimator_class = "pearson_corr",
    method = "pearson",
    description = "signed nonfinite",
    symmetric = TRUE
  )
  signed_td <- tidy(signed_obj)
  expect_true(identical(
    signed_td$estimate[signed_td$item1 == "a" & signed_td$item2 == "b"],
    -0,
    num.eq = FALSE,
    single.NA = FALSE,
    attrib.as.set = FALSE
  ))
  expect_false(any(is.nan(signed_td$estimate) | is.infinite(signed_td$estimate)))
})

test_that("a thin symmetric rectangle still reports a reverse bound error", {
  rectangular <- make_edge(1L, 50L, 0.2, c(20L, 50L), list(NULL, NULL), symmetric = TRUE)
  expect_error(tidy(rectangular), "subscript")
})
