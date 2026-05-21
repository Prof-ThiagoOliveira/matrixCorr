# Shared correlation-result helpers.
#
# These helpers are used by several correlation and agreement estimators. Keep
# them outside estimator-specific files so result construction does not depend
# on latent-correlation implementation details.

.mc_square_dimnames <- function(names) {
  list(names, names)
}

.mc_set_matrix_dimnames <- function(x, row_names = NULL, col_names = row_names) {
  if (is.null(row_names) && is.null(col_names)) {
    return(x)
  }
  dimnames(x) <- list(row_names, col_names)
  x
}

.mc_structure_corr_matrix <- function(mat, class_name, method, description,
                                      diagnostics = NULL, thresholds = NULL,
                                      correct = NULL, dimnames = NULL,
                                      symmetric = NULL,
                                      extra_attrs = NULL,
                                      classes = c(class_name, "matrix")) {
  if (!is.null(dimnames)) {
    dimnames(mat) <- dimnames
  }
  symmetric_flag <- symmetric
  if (is.null(symmetric_flag)) {
    symmetric_flag <- isTRUE(nrow(mat) == ncol(mat)) &&
      isTRUE(isSymmetric(mat, check.attributes = FALSE))
  }
  keep_classes <- setdiff(
    classes,
    c(class_name, "matrix", "corr_matrix", "corr_result")
  )
  .mc_new_corr_matrix(
    mat = mat,
    estimator_class = class_name,
    method = method,
    description = description,
    output = "matrix",
    threshold = 0,
    diag = TRUE,
    diagnostics = diagnostics,
    ci = attr(mat, "ci", exact = TRUE),
    conf.level = attr(mat, "conf.level", exact = TRUE),
    symmetric = symmetric_flag,
    extra_attrs = c(
      list(
        thresholds = thresholds,
        correct = correct
      ),
      extra_attrs %||% list()
    ),
    extra_classes = keep_classes
  )
}
