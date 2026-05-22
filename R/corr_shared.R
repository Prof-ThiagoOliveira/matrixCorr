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
                                          has_inference = FALSE,
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
    isFALSE(has_inference) &&
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
    has_inference = has_inference,
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

.mc_ci_attr <- function(x) {
  attr(x, "ci", exact = TRUE)
}

.mc_inference_attr <- function(x) {
  attr(x, "inference", exact = TRUE)
}

.mc_round_finite <- function(x, digits) {
  x <- as.numeric(x)
  if (is.null(digits)) {
    return(x)
  }
  out <- rep(NA_real_, length(x))
  keep <- is.finite(x)
  out[keep] <- round(x[keep], digits)
  out
}

.mc_pairwise_matrix_summary <- function(object,
                                        class_name,
                                        digits,
                                        ci_digits = NULL,
                                        p_digits = NULL,
                                        show_ci = NULL,
                                        ci_attr = .mc_ci_attr(object),
                                        inference_attr = .mc_inference_attr(object),
                                        diagnostics_attr = attr(object, "diagnostics", exact = TRUE),
                                        include_n_complete = TRUE,
                                        include_ci = TRUE,
                                        include_p = TRUE,
                                        extra_columns = NULL,
                                        extra_attrs = list(),
                                        overview = NULL) {
  est <- as.matrix(object)
  nr <- nrow(est)
  nc <- ncol(est)
  rn <- rownames(est)
  cn <- colnames(est)
  if (is.null(rn)) rn <- as.character(seq_len(nr))
  if (is.null(cn)) cn <- as.character(seq_len(nc))

  n_pairs <- nr * (nc - 1L) / 2L
  var1 <- character(n_pairs)
  var2 <- character(n_pairs)
  estimate <- numeric(n_pairs)

  has_n_complete <- isTRUE(include_n_complete) &&
    is.list(diagnostics_attr) &&
    is.matrix(diagnostics_attr$n_complete)
  n_complete <- if (has_n_complete) integer(n_pairs) else NULL

  show_ci <- if (is.null(show_ci)) "yes" else show_ci
  has_ci <- isTRUE(include_ci) &&
    identical(show_ci, "yes") &&
    is.list(ci_attr)
  lwr <- if (has_ci) numeric(n_pairs) else NULL
  upr <- if (has_ci) numeric(n_pairs) else NULL

  extra_columns <- extra_columns %||% list()
  extra_values <- vector("list", length(extra_columns))
  names(extra_values) <- names(extra_columns)
  if (length(extra_columns)) {
    for (nm in names(extra_columns)) {
      spec <- extra_columns[[nm]]
      type <- spec$type %||% "numeric"
      extra_values[[nm]] <- switch(
        type,
        integer = integer(n_pairs),
        logical = logical(n_pairs),
        character = character(n_pairs),
        numeric(n_pairs)
      )
    }
  }

  k <- 0L
  for (i in seq_len(nr - 1L)) {
    for (j in (i + 1L):nc) {
      k <- k + 1L
      var1[[k]] <- rn[[i]]
      var2[[k]] <- cn[[j]]
      estimate[[k]] <- round(est[i, j], digits)
      if (!is.null(n_complete)) {
        n_complete[[k]] <- as.integer(diagnostics_attr$n_complete[i, j])
      }
      if (has_ci) {
        lwr[[k]] <- if (!is.null(ci_attr$lwr.ci) && is.finite(ci_attr$lwr.ci[i, j])) {
          round(ci_attr$lwr.ci[i, j], ci_digits)
        } else {
          NA_real_
        }
        upr[[k]] <- if (!is.null(ci_attr$upr.ci) && is.finite(ci_attr$upr.ci[i, j])) {
          round(ci_attr$upr.ci[i, j], ci_digits)
        } else {
          NA_real_
        }
      }
      if (length(extra_columns)) {
        for (nm in names(extra_columns)) {
          spec <- extra_columns[[nm]]
          mat <- spec$matrix
          value <- if (is.matrix(mat)) {
            mat[i, j]
          } else if (is.function(spec$value)) {
            spec$value(i, j)
          } else {
            NA
          }
          type <- spec$type %||% "numeric"
          extra_values[[nm]][[k]] <- switch(
            type,
            integer = as.integer(value),
            logical = isTRUE(value),
            character = as.character(value),
            .mc_round_finite(value, if ("digits" %in% names(spec)) spec$digits else digits)
          )
        }
      }
    }
  }

  df <- data.frame(
    var1 = var1,
    var2 = var2,
    estimate = as.numeric(estimate),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!is.null(n_complete)) df$n_complete <- as.integer(n_complete)
  if (!is.null(lwr)) df$lwr <- as.numeric(lwr)
  if (!is.null(upr)) df$upr <- as.numeric(upr)
  if (length(extra_values)) {
    for (nm in names(extra_values)) {
      df[[nm]] <- extra_values[[nm]]
    }
  }
  rownames(df) <- NULL

  out <- .mc_finalize_summary_df(df, class_name = class_name)
  attr(out, "overview") <- overview %||% .mc_summary_corr_matrix(object)
  attr(out, "has_ci") <- has_ci
  if (!is.null(include_p)) {
    attr(out, "has_p") <- isTRUE(include_p)
  }
  attr(out, "conf.level") <- if (is.null(ci_attr)) NA_real_ else ci_attr$conf.level
  attr(out, "digits") <- digits
  if (!is.null(ci_digits)) attr(out, "ci_digits") <- ci_digits
  if (!is.null(p_digits)) attr(out, "p_digits") <- p_digits
  if (length(extra_attrs)) {
    for (nm in names(extra_attrs)) {
      attr(out, nm) <- extra_attrs[[nm]]
    }
  }
  out
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
        has_inference = FALSE,
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

.mc_inference_corr_wrapper <- function(data,
                                       na_method,
                                       ci = FALSE,
                                       p_value = FALSE,
                                       conf_level = 0.95,
                                       n_threads = getOption("matrixCorr.threads", 1L),
                                       output = c("matrix", "sparse", "edge_list"),
                                       threshold = 0,
                                       diag = TRUE,
                                       estimator_class,
                                       method,
                                       description,
                                       kernel_matrix,
                                       kernel_pairwise = NULL,
                                       kernel_threshold = NULL,
                                       payload_builder = NULL,
                                       min_n = 2L,
                                       direct_method = method,
                                       symmetric = TRUE,
                                       extra_attrs = list(),
                                       n_boot = NULL,
                                       seed = NULL,
                                       arg_data = "data") {
  output_cfg <- .mc_validate_thresholded_output_request(
    output = output,
    threshold = threshold,
    diag = diag
  )
  na_method <- match.arg(na_method, c("error", "pairwise", "complete"))
  check_bool(ci, arg = "ci")
  check_bool(p_value, arg = "p_value")
  needs_ci <- isTRUE(ci)
  needs_p_value <- isTRUE(p_value)
  needs_bootstrap <- needs_ci
  if (needs_bootstrap) {
    check_prob_scalar(conf_level, arg = "conf_level", open_ends = TRUE)
    n_boot <- check_scalar_int_pos(n_boot %||% 500L, arg = "n_boot")
    if (!is.null(seed)) {
      seed <- check_scalar_int_pos(seed, arg = "seed")
    }
  } else {
    conf_level <- 0.95
    n_boot <- n_boot %||% 500L
    seed <- NULL
  }
  n_threads <- check_scalar_int_pos(n_threads, arg = "n_threads")

  numeric_data <- if (identical(na_method, "error")) {
    validate_corr_input(data)
  } else {
    validate_corr_input(data, check_na = FALSE)
  }
  diagnostics_extra <- NULL
  if (identical(na_method, "complete")) {
    cc <- .mc_complete_case_matrix(numeric_data, min_n = min_n, arg = arg_data)
    numeric_data <- cc$data
    diagnostics_extra <- cc$diagnostics
  }
  colnames_data <- colnames(numeric_data)
  dn <- if (is.null(colnames_data)) NULL else .mc_square_dimnames(colnames_data)

  .mc_with_omp_threads(
    n_threads,
    {
      direct <- .mc_try_direct_triplet_output(
        input = list(data = numeric_data, dimnames = dn),
        output_cfg = output_cfg,
        kernel_threshold = if (is.null(kernel_threshold)) {
          NULL
        } else {
          function(x, threshold, diag) {
            kernel_threshold(
              x,
              threshold = threshold,
              diag = diag,
              n_threads = n_threads
            )
          }
        },
        estimator_class = estimator_class,
        method = direct_method,
        description = description,
        na_method = na_method,
        ci = ci,
        has_inference = needs_p_value,
        pairwise = identical(na_method, "pairwise"),
        symmetric = symmetric
      )
      if (!is.null(direct)) {
        direct
      } else {
        use_pairwise_kernel <- identical(na_method, "pairwise")
        if (!isTRUE(use_pairwise_kernel)) {
          est <- kernel_matrix(numeric_data, n_threads = n_threads)
        } else {
          if (is.null(kernel_pairwise)) {
            abort_internal("Pairwise kernel is required when {.arg na_method} is {.val pairwise}.")
          }
          est <- kernel_pairwise(numeric_data, n_threads = n_threads)
        }

        est <- .mc_set_matrix_dimnames(est, colnames_data)
        payload <- NULL
        if ((isTRUE(ci) || isTRUE(p_value)) && !is.null(payload_builder)) {
          payload <- payload_builder(
            numeric_data,
            est = est,
            ci = ci,
            p_value = p_value,
            conf_level = conf_level,
            n_boot = n_boot,
            seed = seed
          )
        }

        diagnostics <- .mc_merge_diagnostics(
          if (is.null(payload)) NULL else payload$diagnostics,
          diagnostics_extra
        )

        out <- .mc_structure_corr_matrix(
          est,
          class_name = estimator_class,
          method = method,
          description = description,
          symmetric = symmetric,
          diagnostics = diagnostics,
          dimnames = dn,
          extra_attrs = c(
            extra_attrs,
            if (!is.null(payload) && !is.null(payload$ci)) {
              list(
                ci = payload$ci,
                conf.level = conf_level,
                n_boot = n_boot
              )
            },
            if (!is.null(payload) && !is.null(payload$inference)) {
              list(inference = payload$inference)
            }
          )
        )

        .mc_finalize_corr_output_fast(
          out,
          output = output_cfg$output,
          threshold = output_cfg$threshold,
          diag = output_cfg$diag
        )
      }
    }
  )
}
