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

.mc_resolve_corr_na <- function(na_method,
                                dots = list(),
                                na_method_missing = FALSE,
                                allowed = c("error", "pairwise", "complete"),
                                warn = TRUE) {
  if (!length(dots) && isTRUE(na_method_missing)) {
    return(list(na_method = "error", check_na = TRUE))
  }

  legacy_args <- .mc_extract_legacy_aliases(dots, allowed = "check_na")
  resolve_na_args(
    na_method = na_method,
    check_na = legacy_args$check_na %||% NULL,
    na_method_missing = na_method_missing,
    allowed = allowed,
    warn = warn
  )
}

.mc_prepare_corr_input <- function(data,
                                   na_cfg,
                                   min_n = 2L,
                                   arg = "data") {
  numeric_data <- validate_corr_input(data, check_na = na_cfg$check_na)
  diagnostics <- NULL

  if (identical(na_cfg$na_method, "complete")) {
    cc <- .mc_complete_case_matrix(numeric_data, min_n = min_n, arg = arg)
    numeric_data <- cc$data
    diagnostics <- cc$diagnostics
  }

  col_names <- colnames(numeric_data)
  list(
    data = numeric_data,
    colnames = col_names,
    dimnames = if (is.null(col_names)) NULL else .mc_square_dimnames(col_names),
    diagnostics = diagnostics
  )
}

.mc_with_omp_threads <- function(n_threads,
                                 code) {
  prev_threads <- .mc_prepare_omp_threads(
    n_threads
  )
  if (!is.null(prev_threads)) {
    on.exit(.mc_exit_omp_threads(prev_threads), add = TRUE)
  }
  force(code)
}

.mc_prepare_corr_output <- function(output = c("matrix", "sparse", "edge_list"),
                                    threshold = 0,
                                    diag = TRUE,
                                    thresholded = FALSE) {
  validate <- if (isTRUE(thresholded)) {
    .mc_validate_thresholded_output_request
  } else {
    .mc_validate_output_args
  }
  validate(output = output, threshold = threshold, diag = diag)
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

.mc_finalize_corr_result <- function(mat,
                                     class_name,
                                     method,
                                     description,
                                     output_cfg,
                                     diagnostics = NULL,
                                     dimnames = NULL,
                                     symmetric = NULL,
                                     extra_attrs = NULL,
                                     classes = c(class_name, "matrix")) {
  out <- .mc_structure_corr_matrix(
    mat = mat,
    class_name = class_name,
    method = method,
    description = description,
    diagnostics = diagnostics,
    dimnames = dimnames,
    symmetric = symmetric,
    extra_attrs = extra_attrs,
    classes = classes
  )

  .mc_finalize_corr_output_fast(
    out,
    output = output_cfg$output,
    threshold = output_cfg$threshold,
    diag = output_cfg$diag
  )
}

.mc_finalize_edge_list_zero_threshold <- function(triplets,
                                                  estimator_class,
                                                  method,
                                                  description,
                                                  source_dim,
                                                  source_dimnames = NULL,
                                                  diag = TRUE,
                                                  diagnostics = NULL,
                                                  ci = NULL,
                                                  conf.level = NULL,
                                                  symmetric = TRUE,
                                                  package_name = "matrixCorr",
                                                  extra_attrs = list()) {
  ii <- as.integer(triplets$i)
  jj <- as.integer(triplets$j)
  vv <- as.numeric(triplets$x)

  rn <- NULL
  cn <- NULL
  if (is.list(source_dimnames) && length(source_dimnames) == 2L) {
    rn <- source_dimnames[[1L]]
    cn <- source_dimnames[[2L]]
  }

  row_out <- if (is.null(rn)) as.character(ii) else rn[ii]
  col_out <- if (is.null(cn)) as.character(jj) else cn[jj]

  out <- .mc_new_corr_edge_list(
    df = data.frame(
      row = row_out,
      col = col_out,
      value = vv,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    estimator_class = estimator_class,
    method = method,
    description = description,
    threshold = 0,
    diag = diag,
    diagnostics = diagnostics,
    ci = ci,
    conf.level = conf.level,
    source_dim = source_dim,
    source_dimnames = source_dimnames,
    symmetric = symmetric,
    package_name = package_name,
    extra_attrs = extra_attrs
  )

  attr(out, "matrixCorr_meta") <- list(
    source_class = estimator_class,
    method = method,
    description = description,
    package = package_name,
    diagnostics = diagnostics
  )
  out
}

.mc_try_direct_triplet_output <- function(input,
                                          output_cfg,
                                          kernel_threshold,
                                          estimator_class,
                                          method,
                                          description,
                                          na_method = "error",
                                          ci = FALSE,
                                          pairwise = FALSE,
                                          symmetric = TRUE,
                                          diagnostics = NULL,
                                          conf.level = NULL,
                                          extra_attrs = list()) {
  if (is.null(kernel_threshold)) {
    return(NULL)
  }

  source_dim <- as.integer(c(ncol(input$data), ncol(input$data)))
  can_use_error_no_ci <- identical(na_method, "error") &&
    isFALSE(ci) &&
    !isTRUE(pairwise)

  if (identical(output_cfg$output, "edge_list") &&
      isTRUE(output_cfg$threshold == 0) &&
      isTRUE(can_use_error_no_ci)) {
    trip <- kernel_threshold(
      input$data,
      threshold = 0,
      diag = output_cfg$diag
    )
    return(.mc_finalize_edge_list_zero_threshold(
      triplets = trip,
      estimator_class = estimator_class,
      method = method,
      description = description,
      source_dim = source_dim,
      source_dimnames = input$dimnames,
      diag = output_cfg$diag,
      diagnostics = diagnostics,
      conf.level = conf.level,
      symmetric = symmetric,
      extra_attrs = extra_attrs
    ))
  }

  if (.mc_supports_direct_threshold_path(
    method = method,
    na_method = na_method,
    ci = ci,
    output = output_cfg$output,
    threshold = output_cfg$threshold,
    pairwise = isTRUE(pairwise),
    has_ci = ci,
    symmetric = symmetric
  )) {
    trip <- kernel_threshold(
      input$data,
      threshold = output_cfg$threshold,
      diag = output_cfg$diag
    )
    return(.mc_finalize_triplets_output(
      triplets = trip,
      output = output_cfg$output,
      estimator_class = estimator_class,
      method = method,
      description = description,
      threshold = output_cfg$threshold,
      diag = output_cfg$diag,
      source_dim = source_dim,
      source_dimnames = input$dimnames,
      diagnostics = diagnostics,
      conf.level = conf.level,
      symmetric = symmetric,
      extra_attrs = extra_attrs
    ))
  }

  NULL
}

.mc_pairwise_ci_attr <- function(result,
                                 pairwise,
                                 colnames,
                                 conf_level_name = "conf_level") {
  list(
    est = .mc_set_matrix_dimnames(unclass(result), colnames),
    lwr.ci = .mc_set_matrix_dimnames(unclass(pairwise$lwr), colnames),
    upr.ci = .mc_set_matrix_dimnames(unclass(pairwise$upr), colnames),
    conf.level = pairwise[[conf_level_name]]
  )
}

.mc_corr_wrapper <- function(data,
                             dots,
                             na_method,
                             na_method_missing,
                             ci,
                             conf_level,
                             n_threads,
                             output,
                             threshold,
                             diag,
                             estimator_class,
                             method,
                             description,
                             kernel_matrix,
                             kernel_pairwise = NULL,
                             kernel_threshold = NULL,
                             min_n = 2L,
                             symmetric = TRUE,
                             thresholded_output = TRUE,
                             structure_matrix = NULL,
                             structure_matrix_fast = NULL) {
  output_cfg <- .mc_prepare_corr_output(
    output = output,
    threshold = threshold,
    diag = diag,
    thresholded = thresholded_output
  )

  na_cfg <- .mc_resolve_corr_na(
    na_method = na_method,
    dots = dots,
    na_method_missing = na_method_missing,
    allowed = c("error", "pairwise", "complete")
  )

  if (!isFALSE(ci)) {
    check_bool(ci, arg = "ci")
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
  } else if (!is.logical(ci) || length(ci) != 1L || is.na(ci)) {
    check_bool(ci, arg = "ci")
  }

  input <- .mc_prepare_corr_input(data, na_cfg = na_cfg, min_n = min_n)
  diagnostics <- NULL
  ci_attr <- NULL

  .mc_with_omp_threads(
    n_threads,
    {
      direct <- .mc_try_direct_triplet_output(
        input = input,
        output_cfg = output_cfg,
        kernel_threshold = kernel_threshold,
        estimator_class = estimator_class,
        method = method,
        description = description,
        na_method = na_cfg$na_method,
        ci = ci,
        pairwise = identical(na_cfg$na_method, "pairwise"),
        symmetric = symmetric
      )
      if (!is.null(direct)) {
        direct
      } else {
        if (!identical(na_cfg$na_method, "pairwise") && !isTRUE(ci)) {
          result <- kernel_matrix(input$data)
        } else {
          if (is.null(kernel_pairwise)) {
            abort_internal("Pairwise kernel is required for pairwise missingness or confidence intervals.")
          }
          pairwise <- kernel_pairwise(
            input$data,
            return_ci = ci,
            conf_level = conf_level
          )
          result <- pairwise$est
          diagnostics <- list(
            n_complete = .mc_set_matrix_dimnames(pairwise$n_complete, input$colnames)
          )
          if (isTRUE(ci)) {
            ci_attr <- .mc_pairwise_ci_attr(result, pairwise, input$colnames)
          }
        }

        diagnostics <- .mc_merge_diagnostics(diagnostics, input$diagnostics)

        if (identical(output_cfg$output, "matrix") &&
            is.null(diagnostics) &&
            is.null(ci_attr) &&
            !is.null(structure_matrix_fast)) {
          structure_matrix_fast(
            mat = result,
            dimnames = input$dimnames
          )
        } else {
          if (!is.null(structure_matrix)) {
            out <- structure_matrix(
              mat = result,
              dimnames = input$dimnames,
              diagnostics = diagnostics,
              ci_attr = ci_attr,
              conf_level = if (!is.null(ci_attr)) conf_level else NULL
            )
          } else {
            out <- .mc_structure_corr_matrix(
              result,
              class_name = estimator_class,
              method = method,
              description = description,
              symmetric = symmetric,
              diagnostics = diagnostics,
              dimnames = input$dimnames,
              extra_attrs = if (!is.null(ci_attr)) {
                list(
                  ci = ci_attr,
                  conf.level = conf_level
                )
              }
            )
          }

          .mc_finalize_corr_output_fast(
            out,
            output = output_cfg$output,
            threshold = output_cfg$threshold,
            diag = output_cfg$diag
          )
        }
      }
    }
  )
}
