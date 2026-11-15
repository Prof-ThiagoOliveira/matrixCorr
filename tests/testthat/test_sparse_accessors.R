pair_key <- function(df) {
  paste(df$item1, df$item2, sep = "\r")
}

aligned_columns <- function(reference, other, columns) {
  columns <- intersect(columns, intersect(names(reference), names(other)))
  other[match(pair_key(reference), pair_key(other)), columns, drop = FALSE]
}

expect_registered <- function(generic, class) {
  found <- utils::getS3method(generic, class, optional = TRUE)
  expect_false(
    is.null(found),
    info = sprintf("%s.%s is not registered", generic, class)
  )
  found
}

test_that("sparseMatrix accessors dispatch through metadata and leave unmarked Matrix alone", {
  set.seed(611)
  X <- cbind(a = rnorm(25), b = rnorm(25), c = rnorm(25))
  sp <- pearson_corr(X, ci = TRUE, conf_level = 0.95, output = "sparse", threshold = 0, diag = TRUE)
  raw <- as(Matrix::sparseMatrix(i = 1, j = 2, x = 1, dims = c(2, 2)), "generalMatrix")

  for (generic in c("estimate", "coef", "tidy", "confint")) {
    expect_registered(generic, "sparseMatrix")
    expect_null(utils::getS3method(generic, "dsCMatrix", optional = TRUE))
    expect_null(utils::getS3method(generic, "dgCMatrix", optional = TRUE))
    expect_null(utils::getS3method(generic, "Matrix", optional = TRUE))
  }

  expect_identical(estimate(sp), sp)
  expect_identical(coef(sp), sp)
  expect_s4_class(estimate(sp), "dsCMatrix")
  expect_identical(estimate(sp)@x, sp@x)
  expect_identical(estimate(sp)@i, sp@i)
  expect_identical(estimate(sp)@p, sp@p)
  expect_identical(estimate(sp)@Dim, sp@Dim)
  expect_identical(estimate(sp)@Dimnames, sp@Dimnames)
  expect_identical(attr(estimate(sp), "output", exact = TRUE), "sparse")
  expect_identical(attr(estimate(sp), "threshold", exact = TRUE), attr(sp, "threshold", exact = TRUE))
  expect_identical(attr(estimate(sp), "diag", exact = TRUE), attr(sp, "diag", exact = TRUE))

  expect_error(estimate(raw), "does not expose a matrixCorr estimate")
  expect_error(tidy(raw), "does not expose a tidy matrixCorr representation")
  expect_error(coef(raw), "\\$ operator not defined for this S4 class")
  expect_error(confint(raw), "\\$ operator not defined for this S4 class")
})

test_that("sparse tidy keeps explicit entries, including stored zeros, and does not densify", {
  set.seed(612)
  X <- cbind(a = rnorm(30), b = rnorm(30), c = rnorm(30), d = rnorm(30))
  sp <- pearson_corr(X, ci = FALSE, output = "sparse", threshold = 0.2, diag = FALSE)
  sp@x[1L] <- 0
  stored <- Matrix::summary(sp)
  ns <- asNamespace("matrixCorr")
  trace(
    ".mc_corr_as_dense_matrix",
    tracer = quote({
      if (inherits(x, "sparseMatrix") || identical(attr(x, "output", exact = TRUE), "sparse")) {
        stop("sparse densified", call. = FALSE)
      }
    }),
    print = FALSE,
    where = ns
  )
  on.exit(untrace(".mc_corr_as_dense_matrix", where = ns), add = TRUE)

  td <- tidy(sp)
  expect_s3_class(td, "data.frame")
  expect_true(any(td$estimate == 0))
  expect_false(any(td$item1 == td$item2))
  expect_identical(
    paste(td$item1, td$item2, td$estimate),
    paste(colnames(sp)[stored$i], colnames(sp)[stored$j], stored$x)
  )

  omitted <- expand.grid(i = seq_len(nrow(sp)), j = seq_len(ncol(sp)), KEEP.OUT.ATTRS = FALSE)
  omitted <- omitted[omitted$i < omitted$j, , drop = FALSE]
  stored_key <- paste(stored$i, stored$j, sep = "\r")
  omitted <- omitted[!paste(omitted$i, omitted$j, sep = "\r") %in% stored_key, , drop = FALSE]
  expect_true(nrow(omitted) > 0L)
  tidy_key <- paste(match(td$item1, colnames(sp)), match(td$item2, colnames(sp)), sep = "\r")
  expect_false(any(paste(omitted$i, omitted$j, sep = "\r") %in% tidy_key))

  diag_off <- pearson_corr(X, ci = FALSE, output = "sparse", threshold = 0, diag = FALSE)
  expect_false(any(tidy(diag_off, diag = TRUE)$item1 == tidy(diag_off, diag = TRUE)$item2))
})

test_that("symmetric sparse triangles mirror stored entries and directed sparse defaults to full", {
  set.seed(613)
  X <- cbind(a = rnorm(24), b = rnorm(24), c = rnorm(24))
  sp <- pearson_corr(X, ci = FALSE, output = "sparse", threshold = 0, diag = TRUE)
  stored <- Matrix::summary(sp)
  upper <- tidy(sp, triangle = "upper", diag = FALSE)
  lower <- tidy(sp, triangle = "lower", diag = FALSE)
  full <- tidy(sp, triangle = "full", diag = FALSE)
  stored_off <- stored[stored$i != stored$j, , drop = FALSE]

  expect_identical(
    pair_key(upper),
    paste(colnames(sp)[stored_off$i], colnames(sp)[stored_off$j], sep = "\r")
  )
  expect_identical(
    pair_key(lower),
    paste(colnames(sp)[stored_off$j], colnames(sp)[stored_off$i], sep = "\r")
  )
  expect_identical(pair_key(full), c(pair_key(upper), pair_key(lower)))
  expect_identical(upper$estimate, stored_off$x)
  expect_identical(lower$estimate, stored_off$x)
  expect_false(isTRUE(all(abs(diff(abs(full$estimate))) <= 0)))

  x <- seq(-1, 1, length.out = 40)
  directed <- xi_corr(
    cbind(x = x, x2 = x^2, x3 = x^3),
    tie_method = "first",
    ci = FALSE,
    output = "sparse",
    threshold = 0,
    diag = FALSE
  )
  directed_tidy <- tidy(directed)
  directed_stored <- Matrix::summary(directed)
  expect_s4_class(directed, "dgCMatrix")
  expect_identical(
    pair_key(directed_tidy),
    paste(colnames(directed)[directed_stored$i], colnames(directed)[directed_stored$j], sep = "\r")
  )
  expect_true(any(directed_tidy$item1 == "x" & directed_tidy$item2 == "x2"))
  expect_true(any(directed_tidy$item1 == "x2" & directed_tidy$item2 == "x"))
  expect_equal(nrow(tidy(directed, triangle = "upper")), sum(directed_stored$i < directed_stored$j))
})

test_that("duplicate labels and cross-axis sparse roles follow stored positions", {
  labels <- c("a", "a")
  sp <- Matrix::sparseMatrix(
    i = c(1L, 2L, 1L),
    j = c(1L, 2L, 2L),
    x = c(0.2, 0.4, 0),
    dims = c(2L, 2L),
    dimnames = list(labels, c("b", "b")),
    repr = "C"
  )
  attr(sp, "output") <- "sparse"
  attr(sp, "corr_output_class") <- "corr_sparse"
  attr(sp, "corr_symmetric") <- FALSE
  attr(sp, "corr_estimator_class") <- "polyserial_corr"
  attr(sp, "corr_dim") <- dim(sp)
  attr(sp, "corr_dimnames") <- dimnames(sp)
  td <- tidy(sp)
  stored <- Matrix::summary(sp)
  expect_identical(td$item1, rownames(sp)[stored$i])
  expect_identical(td$item2, colnames(sp)[stored$j])
  expect_identical(td$estimate, stored$x)
  expect_identical(td$item1, c("a", "a", "a"))
  expect_true(any(duplicated(td$item1)) && any(duplicated(td$item2)))

  same <- sp
  dimnames(same) <- list(labels, labels)
  attr(same, "corr_dimnames") <- dimnames(same)
  attr(same, "corr_estimator_class") <- "pearson_corr"
  plain <- tidy(same)
  expect_identical(plain$item1, "a")
  expect_identical(plain$item2, "a")
  expect_identical(plain$estimate, 0)

  cross_square <- same
  attr(cross_square, "corr_estimator_class") <- "biserial_corr"
  crossed <- tidy(cross_square)
  expect_identical(crossed$estimate, stored$x)
  expect_identical(nrow(crossed), nrow(stored))
})

test_that("sparse confint reads stored limits and reports unavailable intervals", {
  set.seed(614)
  X <- cbind(a = rnorm(22), b = rnorm(22), c = rnorm(22))
  sp <- pearson_corr(
    X,
    ci = TRUE,
    conf_level = 0.95,
    output = "sparse",
    threshold = 0,
    diag = FALSE
  )
  limits <- confint(sp)
  td <- tidy(sp)
  expect_identical(limits$item1, td$item1)
  expect_identical(limits$item2, td$item2)
  expect_identical(limits$lwr, td$lwr)
  expect_identical(limits$upr, td$upr)
  expect_no_error(confint(sp, level = 0.95))
  expect_error(confint(sp, level = 0.8), "pre-computed confidence intervals")

  bare <- pearson_corr(X, ci = FALSE, output = "sparse", threshold = 0, diag = FALSE)
  expect_error(confint(bare), "does not contain confidence limits")
  expect_identical(ci(sp), attr(sp, "ci", exact = TRUE))
})

test_that("dense, sparse, and edge accessors retain the same statistical values", {
  set.seed(615)
  X <- cbind(a = rnorm(28), b = rnorm(28), c = rnorm(28), d = rnorm(28))
  binary <- cbind(
    a = sample(0:1, 28, replace = TRUE),
    b = sample(0:1, 28, replace = TRUE),
    c = sample(0:1, 28, replace = TRUE),
    d = sample(0:1, 28, replace = TRUE)
  )
  methods <- list(
    pearson = function(...) pearson_corr(X, ci = TRUE, conf_level = 0.9, ...),
    spearman = function(...) spearman_rho(X, ci = TRUE, conf_level = 0.9, na_method = "error", ...),
    kendall = function(...) kendall_tau(X, ci = TRUE, conf_level = 0.9, na_method = "error", ...),
    bicor = function(...) bicor(X, ci = TRUE, conf_level = 0.9, na_method = "error", ...),
    dcor = function(...) dcor(X, na_method = "error", p_value = FALSE, ...),
    tetrachoric = function(...) tetrachoric(binary, ...),
    ccc = function(...) ccc(X, ci = TRUE, conf_level = 0.9, ...)
  )
  stat_columns <- c(
    "item1", "item2", "estimate", "lwr", "upr", "statistic", "df", "p_value",
    "n_complete", "se", "z_statistic"
  )

  for (nm in names(methods)) {
    dense <- methods[[nm]](output = "matrix")
    sparse <- methods[[nm]](output = "sparse", threshold = 0, diag = FALSE)
    edge <- methods[[nm]](output = "edge_list", threshold = 0, diag = FALSE)
    dense_td <- tidy(dense)
    sparse_td <- tidy(sparse)
    edge_td <- tidy(edge)
    expect_identical(estimate(sparse), sparse, info = nm)
    expect_identical(coef(sparse), sparse, info = nm)
    expect_identical(pair_key(sparse_td), pair_key(edge_td), info = nm)
    dense_cols <- dense_td[, intersect(stat_columns, names(dense_td)), drop = FALSE]
    rownames(dense_cols) <- NULL
    sparse_cols <- aligned_columns(dense_td, sparse_td, stat_columns)
    edge_cols <- aligned_columns(dense_td, edge_td, stat_columns)
    rownames(sparse_cols) <- NULL
    rownames(edge_cols) <- NULL
    expect_identical(pair_key(sparse_td), pair_key(dense_td), info = nm)
    expect_identical(pair_key(edge_td), pair_key(dense_td), info = nm)
    expect_identical(sparse_cols, dense_cols, info = nm)
    expect_identical(edge_cols, dense_cols, info = nm)
    if (all(c("lwr", "upr") %in% names(sparse_td))) {
      expect_identical(confint(sparse)[, c("lwr", "upr")], sparse_td[, c("lwr", "upr")], info = nm)
    }
  }
})

test_that("summary stays on the selected cell while tidy can read a missing opposite", {
  set.seed(616)
  X <- cbind(a = rnorm(20), b = rnorm(20), c = rnorm(20))
  fit <- pearson_corr(X, ci = FALSE)
  payload <- matrix(NA_real_, 3, 3, dimnames = dimnames(fit))
  payload[2, 1] <- 0.2
  payload[3, 1] <- 0.3
  payload[3, 2] <- 0.4
  attr(fit, "inference") <- list(p_value = payload)
  sm <- summary(fit)
  td <- tidy(fit)
  expect_true(all(is.na(sm$p_value)))
  expect_identical(td$p_value, c(0.2, 0.3, 0.4))
})

test_that("pair reader fills only a missing selected cell from the opposite entry", {
  mat <- matrix(
    c(NA, 2, Inf, NaN, NA, -1, 4, NA, 5),
    nrow = 3,
    byrow = TRUE
  )
  idx <- list(
    i = c(1L, 2L, 3L, 1L, 2L),
    j = c(1L, 1L, 2L, 3L, 3L),
    ok = c(TRUE, TRUE, TRUE, TRUE, FALSE)
  )
  expect_identical(
    matrixCorr:::.mc_extract_pair_matrix(mat, idx, c(3L, 3L), symmetric = TRUE),
    c(NA_real_, 2, -1, Inf, NA_real_)
  )
  expect_identical(
    matrixCorr:::.mc_extract_pair_matrix(mat, idx, c(3L, 3L), symmetric = FALSE),
    c(NA_real_, NaN, NA_real_, Inf, NA_real_)
  )

  chars <- matrix(c("keep", "down", NA_character_, "other"), nrow = 2)
  char_idx <- list(i = c(1L, 2L), j = c(2L, 2L), ok = c(TRUE, TRUE))
  expect_identical(
    matrixCorr:::.mc_extract_pair_matrix(
      chars, char_idx, c(2L, 2L), character = TRUE, symmetric = TRUE
    ),
    c("down", "other")
  )

  flags <- matrix(c(TRUE, FALSE, NA, TRUE), nrow = 2)
  expect_identical(
    matrixCorr:::.mc_extract_pair_matrix(
      flags, char_idx, c(2L, 2L), logical = TRUE, symmetric = TRUE
    ),
    c(FALSE, TRUE)
  )

  counts <- matrix(c(NA_real_, 2.5, NA_real_, NA_real_), nrow = 2)
  expect_identical(
    matrixCorr:::.mc_extract_pair_matrix(
      counts,
      list(i = c(1L, 2L), j = c(2L, 2L), ok = c(TRUE, TRUE)),
      c(2L, 2L),
      integer = TRUE,
      symmetric = TRUE
    ),
    c(2L, NA_integer_)
  )
  expect_identical(
    matrixCorr:::.mc_extract_pair_matrix(mat, idx, c(2L, 2L), symmetric = TRUE),
    rep(NA_real_, 5)
  )
  expect_identical(
    matrixCorr:::.mc_extract_pair_matrix(
      matrix(NA_real_, 3, 3),
      list(i = c(1L, 2L), j = c(2L, 1L), ok = c(TRUE, TRUE)),
      c(3L, 3L),
      symmetric = TRUE
    ),
    c(NA_real_, NA_real_)
  )
})

stamp_sparse <- function(x, symmetric = TRUE, estimator = "pearson_corr", ci = NULL) {
  attr(x, "output") <- "sparse"
  attr(x, "corr_output_class") <- "corr_sparse"
  attr(x, "corr_symmetric") <- symmetric
  attr(x, "corr_estimator_class") <- estimator
  attr(x, "corr_dim") <- dim(x)
  attr(x, "corr_dimnames") <- dimnames(x)
  if (!is.null(ci)) {
    attr(x, "ci") <- ci
    attr(x, "conf.level") <- 0.95
  }
  x
}

test_that("upper tidy recovers a lower-stored symmetric triangle", {
  labels <- c("a", "b", "c")
  sp <- Matrix::sparseMatrix(
    i = c(2L, 3L),
    j = c(1L, 1L),
    x = c(0, 0.7),
    dims = c(3L, 3L),
    dimnames = list(labels, labels),
    symmetric = TRUE
  )
  sp <- stamp_sparse(sp)
  expect_identical(sp@uplo, "L")
  expect_identical(estimate(sp), sp)
  expect_identical(coef(sp), sp)
  upper <- tidy(sp)
  expect_identical(upper$item1, c("a", "a"))
  expect_identical(upper$item2, c("b", "c"))
  expect_identical(upper$estimate, c(0, 0.7))

  general <- Matrix::sparseMatrix(
    i = c(2L, 3L),
    j = c(1L, 1L),
    x = c(0, 0.7),
    dims = c(3L, 3L),
    dimnames = list(labels, labels),
    repr = "C"
  )
  general <- stamp_sparse(general)
  expect_s4_class(general, "dgCMatrix")
  expect_identical(tidy(general)$estimate, c(0, 0.7))

  both <- Matrix::sparseMatrix(
    i = c(2L, 1L, 3L, 1L),
    j = c(1L, 2L, 1L, 3L),
    x = c(0, 0, 0.7, 0.7),
    dims = c(3L, 3L),
    dimnames = list(labels, labels),
    repr = "C"
  )
  both <- stamp_sparse(both)
  full <- tidy(both, triangle = "full")
  expect_identical(nrow(full), 4L)
  expect_identical(sum(full$estimate == 0), 2L)
  expect_identical(sum(full$estimate == 0.7), 2L)
})

test_that("empty sparse tidy keeps typed columns and stored intervals", {
  labels <- c("a", "b")
  bare <- Matrix::sparseMatrix(
    i = integer(),
    j = integer(),
    x = numeric(),
    dims = c(2L, 2L),
    dimnames = list(labels, labels),
    repr = "C"
  )
  bare <- stamp_sparse(bare)
  td <- tidy(bare)
  expect_identical(nrow(td), 0L)
  expect_identical(names(td), c("item1", "item2", "estimate"))
  expect_type(td$item1, "character")
  expect_type(td$item2, "character")
  expect_type(td$estimate, "double")
  expect_error(confint(bare), "does not contain confidence limits")

  limits <- matrix(c(NA, 0.1, 0.2, NA), 2, 2, dimnames = list(labels, labels))
  with_ci <- stamp_sparse(
    bare,
    ci = list(lwr.ci = limits, upr.ci = limits + 0.1, conf.level = 0.95)
  )
  td_ci <- tidy(with_ci)
  expect_identical(nrow(td_ci), 0L)
  expect_identical(names(td_ci), c("item1", "item2", "estimate", "lwr", "upr"))
  expect_type(td_ci$lwr, "double")
  expect_type(td_ci$upr, "double")
  expect_null(attr(td_ci, "conf.level", exact = TRUE))
  ci_df <- confint(with_ci)
  expect_identical(nrow(ci_df), 0L)
  expect_identical(names(ci_df), c("item1", "item2", "lwr", "upr"))

  lower_only <- stamp_sparse(bare, ci = list(lwr.ci = limits, conf.level = 0.95))
  expect_identical(names(tidy(lower_only)), c("item1", "item2", "estimate"))
  expect_error(confint(lower_only), "does not contain confidence limits")
  mismatch <- stamp_sparse(
    bare,
    ci = list(
      lwr.ci = limits,
      upr.ci = matrix(0.2, 3, 3),
      conf.level = 0.95
    )
  )
  expect_identical(names(tidy(mismatch)), c("item1", "item2", "estimate"))
  expect_error(confint(mismatch), "does not contain confidence limits")

  diag_only <- Matrix::sparseMatrix(
    i = c(1L, 2L),
    j = c(1L, 2L),
    x = c(1, 1),
    dims = c(2L, 2L),
    dimnames = list(labels, labels),
    repr = "C"
  )
  diag_only <- stamp_sparse(diag_only)
  expect_identical(nrow(tidy(diag_only, diag = FALSE)), 0L)
  expect_identical(tidy(diag_only, diag = TRUE)$estimate, c(1, 1))

  set.seed(618)
  X <- cbind(a = rnorm(12), b = rnorm(12))
  dropped <- pearson_corr(X, ci = FALSE, output = "sparse", threshold = 1.1, diag = FALSE)
  dropped_td <- tidy(dropped)
  expect_identical(nrow(dropped_td), 0L)
  expect_identical(names(dropped_td), c("item1", "item2", "estimate"))

  x <- seq(-1, 1, length.out = 16)
  empty_xi <- xi_corr(
    cbind(x = x, x2 = x^2),
    tie_method = "first",
    ci = TRUE,
    bootstrap_reps = 9L,
    seed = 3L,
    output = "sparse",
    threshold = 2,
    diag = FALSE
  )
  expect_true(is.list(attr(empty_xi, "diagnostics", exact = TRUE)))
  expect_true(!is.null(attr(empty_xi, "diagnostics", exact = TRUE)$bootstrap_reps))
  xi_td <- tidy(empty_xi)
  expect_identical(nrow(xi_td), 0L)
  expect_true(all(c("item1", "item2", "estimate", "lwr", "upr") %in% names(xi_td)))
  expect_identical(nrow(confint(empty_xi)), 0L)
})
