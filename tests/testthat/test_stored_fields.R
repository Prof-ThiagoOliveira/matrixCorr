stored_pair_value <- function(mat, item1, item2) {
  value <- unname(mat[item1, item2])
  if (is.integer(mat)) {
    as.integer(value)
  } else if (is.logical(mat)) {
    as.logical(value)
  } else {
    value
  }
}

expect_stored_column <- function(df, column, mat) {
  expect_true(column %in% names(df))
  got <- mapply(
    stored_pair_value,
    item1 = df$item1,
    item2 = df$item2,
    MoreArgs = list(mat = mat),
    USE.NAMES = FALSE
  )
  if (identical(column, "skipped_n")) {
    expect_type(df[[column]], "integer")
    got <- as.integer(round(got))
  }
  if (identical(column, "reject")) {
    expect_type(mat, "logical")
    expect_type(df[[column]], "logical")
  }
  expect_identical(df[[column]], got)
}

test_that("summary and tidy append stored skip fields without recomputing them", {
  set.seed(711)
  X <- cbind(a = rnorm(16), b = rnorm(16), c = rnorm(16))
  fits <- list(
    hochberg = skipped_corr(
      X,
      method = "pearson",
      na_method = "error",
      ci = FALSE,
      p_value = TRUE,
      p_adjust = "hochberg",
      n_boot = 20L,
      seed = 7L
    ),
    ecp = skipped_corr(
      X,
      method = "pearson",
      na_method = "error",
      ci = FALSE,
      p_value = TRUE,
      p_adjust = "ecp",
      n_boot = 20L,
      n_mc = 10L,
      seed = 8L
    ),
    none = skipped_corr(
      X,
      method = "pearson",
      na_method = "error",
      ci = FALSE,
      p_value = TRUE,
      p_adjust = "none",
      n_boot = 20L,
      seed = 9L
    )
  )

  for (nm in names(fits)) {
    base <- fits[[nm]]
    representations <- list(
      matrix = base,
      sparse = skipped_corr(
        X,
        method = "pearson",
        na_method = "error",
        ci = FALSE,
        p_value = TRUE,
        p_adjust = nm,
        n_boot = 20L,
        n_mc = 10L,
        seed = if (nm == "hochberg") 7L else if (nm == "ecp") 8L else 9L,
        output = "sparse",
        threshold = 0,
        diag = TRUE
      ),
      edge = skipped_corr(
        X,
        method = "pearson",
        na_method = "error",
        ci = FALSE,
        p_value = TRUE,
        p_adjust = nm,
        n_boot = 20L,
        n_mc = 10L,
        seed = if (nm == "hochberg") 7L else if (nm == "ecp") 8L else 9L,
        output = "edge_list",
        threshold = 0,
        diag = TRUE
      )
    )
    inf <- attr(base, "inference", exact = TRUE)
    diag_attr <- attr(base, "diagnostics", exact = TRUE)
    expect_true(is.matrix(diag_attr$skipped_n))
    expect_true(is.matrix(diag_attr$skipped_prop))
    expect_false(is.matrix(inf$p_value_adjusted) && nm != "hochberg")
    if (nm == "hochberg") {
      expect_true(is.matrix(inf$p_value_adjusted))
      expect_true(is.matrix(inf$reject))
    }
    if (nm == "ecp") {
      expect_null(inf$p_value_adjusted)
      expect_true(is.matrix(inf$reject))
    }
    if (nm == "none") {
      expect_null(inf$p_value_adjusted)
      expect_null(inf$reject)
    }

    compared <- list()
    for (rep_nm in names(representations)) {
      object <- representations[[rep_nm]]
      sm <- summary(object)
      td <- tidy(object, triangle = "upper", diag = FALSE)
      expect_stored_column(sm, "p_value", attr(object, "inference", exact = TRUE)$p_value)
      expect_stored_column(td, "p_value", attr(object, "inference", exact = TRUE)$p_value)
      raw_skipped_n <- attr(object, "diagnostics", exact = TRUE)$skipped_n
      expect_type(raw_skipped_n, "double")
      expect_type(sm$skipped_n, "integer")
      expect_type(td$skipped_n, "integer")
      expect_identical(
        sm$skipped_n,
        as.integer(round(mapply(
          function(a, b) raw_skipped_n[a, b],
          sm$item1, sm$item2,
          USE.NAMES = FALSE
        )))
      )
      expect_identical(
        td$skipped_n,
        as.integer(round(mapply(
          function(a, b) raw_skipped_n[a, b],
          td$item1, td$item2,
          USE.NAMES = FALSE
        )))
      )
      expect_stored_column(sm, "skipped_prop", attr(object, "diagnostics", exact = TRUE)$skipped_prop)
      expect_stored_column(td, "skipped_prop", attr(object, "diagnostics", exact = TRUE)$skipped_prop)
      if (nm == "hochberg") {
        expect_stored_column(sm, "p_value_adjusted", attr(object, "inference", exact = TRUE)$p_value_adjusted)
        expect_stored_column(td, "p_value_adjusted", attr(object, "inference", exact = TRUE)$p_value_adjusted)
        expect_stored_column(sm, "reject", attr(object, "inference", exact = TRUE)$reject)
        expect_stored_column(td, "reject", attr(object, "inference", exact = TRUE)$reject)
        expect_type(td$reject, "logical")
        expect_type(sm$reject, "logical")
      } else {
        expect_false("p_value_adjusted" %in% names(sm))
        expect_false("p_value_adjusted" %in% names(td))
      }
      if (nm == "ecp") {
        expect_stored_column(sm, "reject", attr(object, "inference", exact = TRUE)$reject)
        expect_stored_column(td, "reject", attr(object, "inference", exact = TRUE)$reject)
      }
      if (nm == "none") {
        expect_false("reject" %in% names(sm))
        expect_false("reject" %in% names(td))
      }
      compared[[rep_nm]] <- td[order(td$item1, td$item2), c("item1", "item2", "estimate", "p_value", "skipped_n", "skipped_prop")]
    }
    expect_identical(compared$sparse, compared$matrix)
    expect_identical(compared$edge, compared$matrix)
  }
})

test_that("tidy appends stored bicor Z and bcdcor without replacing the estimate", {
  set.seed(712)
  X <- cbind(a = rnorm(24), b = rnorm(24), c = rnorm(24))
  fit <- bicor(X, ci = TRUE, na_method = "error")
  td <- tidy(fit)
  inf <- attr(fit, "inference", exact = TRUE)
  expect_true(is.matrix(inf$Z))
  expect_identical(
    td$z_statistic,
    mapply(function(i, j) inf$Z[i, j], td$item1, td$item2, USE.NAMES = FALSE)
  )
  expect_false(identical(td$z_statistic, td$estimate))
  expect_false("bcdcor" %in% names(td))

  distance <- bcdcor(X, na_method = "error", p_value = TRUE)
  distance_td <- tidy(distance)
  distance_inf <- attr(distance, "inference", exact = TRUE)
  expect_true(is.matrix(distance_inf$bcdcor))
  expect_identical(
    distance_td$bcdcor,
    mapply(
      function(i, j) distance_inf$bcdcor[i, j],
      distance_td$item1,
      distance_td$item2,
      USE.NAMES = FALSE
    )
  )
  expect_identical(distance_td$estimate, mapply(
    function(i, j) unclass(as.matrix(distance))[i, j],
    distance_td$item1,
    distance_td$item2,
    USE.NAMES = FALSE
  ))
})
