#' Interactive Shiny viewer for matrixCorr objects
#'
#' Launches an interactive Shiny gadget that displays correlation heatmaps with
#' filtering, clustering, and hover inspection. The viewer accepts any
#' matrixCorr correlation result (for example the outputs from
#' [pearson_corr()], [spearman_rho()], [kendall_tau()], [bicor()],
#' [pbcor()], [wincor()], [skipped_corr()], [pcorr()], [dcor()], or
#' [shrinkage_corr()]), a plain
#' matrix, or a named list of such objects. When a list is supplied the gadget
#' offers a picker to switch between results.
#' @param x A correlation result, a numeric matrix, or a named list of those
#'   objects. Each element must be square with matching row/column names.
#' @param title Optional character title shown at the top of the gadget.
#' @param default_max_vars Integer; maximum number of variables pre-selected
#'   when the app opens. Defaults to 40 so very wide matrices start collapsed.
#' @return Invisibly returns `NULL`; the function is called for its side
#'   effect of launching a Shiny gadget.
#' @details This helper lives in `Suggests`; it requires the `shiny` and
#'   `shinyWidgets` packages at runtime and will optionally convert the plot to
#'   an interactive widget when `plotly` is installed. Variable selection uses
#'   a searchable picker, and clustering controls let you reorder variables via
#'   hierarchical clustering on either absolute or signed correlations with a
#'   choice of linkage methods.
#'
#' @importFrom rlang .data
#' @examples
#' if (interactive()) {
#'   data <- mtcars
#'   results <- list(
#'     Pearson = pearson_corr(data),
#'     Spearman = spearman_rho(data),
#'     Kendall = kendall_tau(data)
#'   )
#'   view_corr_shiny(results)
#' }
#'
#' @export
view_corr_shiny <- function(x, title = NULL, default_max_vars = 40L) {
  if (missing(x)) {
    abort_bad_arg("x", message = "must be supplied (correlation result or list).")
  }

  if (!requireNamespace("shiny", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg shiny} is required for {.fn view_corr_shiny}.")
  }
  if (!requireNamespace("shinyWidgets", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg shinyWidgets} is required for {.fn view_corr_shiny}.")
  }

  prepared <- .mc_prepare_corr_inputs(x)
  if (length(prepared) == 0L) {
    cli::cli_abort("{.fn view_corr_shiny} did not find a usable correlation matrix.")
  }

  use_plotly <- .mc_has_namespace("plotly")
  app_title <- title %||% "matrixCorr correlation viewer"
  logo_src <- .mc_register_logo_resource(.mc_logo_path())

  plot_widget <- if (use_plotly) {
    .mc_plotly_fn("plotlyOutput")("heatmap", height = "650px")
  } else {
    shiny::plotOutput("heatmap", height = "650px")
  }

  ui <- shiny::fluidPage(
    .mc_viewer_head(),
    shiny::div(
      class = "mc-app-shell",
      .mc_viewer_header(
        title = app_title,
        subtitle = "Explore matrix-style correlation outputs with searchable subsetting, responsive heatmaps, and focused visual inspection.",
        logo_src = logo_src,
        badges = list(
          .mc_header_badge("Correlation Viewer", tone = "neutral"),
          .mc_header_badge(
            if (use_plotly) "Interactive heatmap" else "Static heatmap",
            tone = if (use_plotly) "accent" else "neutral"
          )
        )
      ),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          width = 4,
          .mc_control_section(
            title = "Data",
            shinyWidgets::pickerInput(
              inputId = "matrix_choice",
              label = "Correlation object",
              choices = names(prepared),
              selected = names(prepared)[[1L]],
              multiple = FALSE,
              options = list(`live-search` = TRUE)
            ),
            shinyWidgets::pickerInput(
              inputId = "var_choice",
              label = "Variables (select to subset)",
              choices = colnames(prepared[[1L]]$matrix),
              selected = .mc_default_vars(prepared[[1L]]$matrix, default_max_vars),
              multiple = TRUE,
              options = list(`live-search` = TRUE, `actions-box` = TRUE)
            )
          ),
          .mc_control_section(
            title = "Display",
            shiny::sliderInput(
              inputId = "threshold",
              label = "Hide cells with |value| below",
              min = 0,
              max = 1,
              step = 0.05,
              value = 0
            ),
            shiny::sliderInput(
              inputId = "corr_range",
              label = "Keep variables with correlation in range",
              min = -1,
              max = 1,
              step = 0.05,
              value = c(-1, 1)
            ),
            shiny::checkboxInput("use_abs", "Colour by absolute value", value = FALSE),
            shiny::checkboxInput("mask_diag", "Hide diagonal", value = FALSE),
            shiny::selectInput(
              inputId = "cluster_mode",
              label = "Cluster variables",
              choices = c(
                "None" = "none",
                "Absolute correlation (|r|)" = "abs",
                "Signed correlation (r)" = "signed"
              ),
              selected = "none"
            ),
            shiny::conditionalPanel(
              condition = "input.cluster_mode !== 'none'",
              shiny::selectInput(
                inputId = "cluster_method",
                label = "Linkage method",
                choices = c(
                  "Complete" = "complete",
                  "Average" = "average",
                  "Single" = "single",
                  "Ward (D2)" = "ward.D2"
                ),
                selected = "complete"
              )
            ),
            shiny::checkboxInput("show_values", "Show cell labels", value = TRUE)
          ),
          .mc_control_section(
            title = "Actions",
            shiny::div(
              class = "mc-viewer-controls",
              shiny::actionButton("reset_controls", "Reset controls"),
              shiny::downloadButton("download_matrix", "Download displayed matrix")
            )
          )
        ),
        shiny::mainPanel(
          width = 8,
          shiny::div(
            class = "mc-tabset-shell",
            shiny::tabsetPanel(
              id = "corr_view",
              shiny::tabPanel(
                title = "Heatmap",
                shiny::div(
                  class = "mc-card mc-main-card",
                  .mc_card_header(
                    title = "Heatmap",
                    subtitle = "Filter, cluster, and inspect the currently displayed matrix."
                  ),
                  plot_widget,
                  shiny::uiOutput("status_alerts")
                ),
                shiny::div(
                  class = "mc-card mc-meta-card",
                  .mc_card_header(
                    title = "Summary",
                    subtitle = "Current object and display state."
                  ),
                  shiny::uiOutput("meta")
                )
              ),
              shiny::tabPanel(
                title = "Distribution",
                shiny::fluidRow(
                  shiny::column(
                    width = 7,
                    shiny::div(
                      class = "mc-card mc-main-card",
                      .mc_card_header(
                        title = "Correlation Distribution",
                        subtitle = "Boxplot of displayed off-diagonal values using unique variable pairs."
                      ),
                      shiny::plotOutput("dist_plot", height = "430px")
                    )
                  ),
                  shiny::column(
                    width = 5,
                    shiny::div(
                      class = "mc-card mc-meta-card",
                      .mc_card_header(
                        title = "Summary Table",
                        subtitle = "Computed from the currently displayed upper triangle."
                      ),
                      shiny::tableOutput("dist_summary")
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )

  server <- function(input, output, session) {
    shiny::observe({
      .mc_update_corr_range_slider(
        session = session,
        input_id = "corr_range",
        use_abs = isTRUE(input$use_abs),
        value = input$corr_range
      )
    })

    shiny::observeEvent(input$matrix_choice, {
      info <- prepared[[input$matrix_choice]]
      vars <- colnames(info$matrix)
      shinyWidgets::updatePickerInput(
        session,
        inputId = "var_choice",
        choices = vars,
        selected = .mc_default_vars(info$matrix, default_max_vars)
      )
      shiny::updateSliderInput(session, "threshold", min = 0, max = 1, value = min(input$threshold %||% 0, 1))
    }, ignoreNULL = FALSE)

    shiny::observeEvent(input$reset_controls, {
      info <- current_matrix()
      shinyWidgets::updatePickerInput(
        session,
        inputId = "var_choice",
        choices = colnames(info$matrix),
        selected = .mc_default_vars(info$matrix, default_max_vars)
      )
      shiny::updateSliderInput(session, "threshold", value = 0)
      shiny::updateCheckboxInput(session, "use_abs", value = FALSE)
      shiny::updateCheckboxInput(session, "mask_diag", value = FALSE)
      shiny::updateSelectInput(session, "cluster_mode", selected = "none")
      shiny::updateSelectInput(session, "cluster_method", selected = "complete")
      shiny::updateCheckboxInput(session, "show_values", value = TRUE)
      .mc_update_corr_range_slider(
        session = session,
        input_id = "corr_range",
        use_abs = FALSE,
        value = c(-1, 1)
      )
    })

    current_matrix <- shiny::reactive({
      prepared[[input$matrix_choice]]
    })

    filtered_matrix <- shiny::reactive({
      info <- current_matrix()
      vars <- input$var_choice
      if (is.null(vars) || !length(vars)) {
        vars <- colnames(info$matrix)
      }
      vars <- intersect(vars, colnames(info$matrix))
      shiny::validate(shiny::need(length(vars) >= 2L, "Select at least two variables."))
      M <- info$matrix[vars, vars, drop = FALSE]
      range_filter <- .mc_filter_matrix_by_range(
        mat = M,
        corr_range = input$corr_range,
        use_abs = isTRUE(input$use_abs)
      )
      M <- range_filter$matrix
      shiny::validate(shiny::need(
        ncol(M) >= 2L,
        "Correlation range filter left fewer than two variables. Widen the range or select more variables."
      ))
      diag_map <- diag(M)
      if (!is.null(colnames(M))) names(diag_map) <- colnames(M)
      cluster_mode <- input$cluster_mode %||% "none"
      cluster_method <- input$cluster_method %||% "complete"
      status_messages <- range_filter$message
      if (!identical(cluster_mode, "none")) {
        if (nrow(M) < 3L) {
          status_messages <- c(status_messages, "Clustering requires at least three variables.")
        } else {
          res <- .mc_reorder_matrix(
            M,
            mode = cluster_mode,
            method = cluster_method
          )
          M <- res$matrix
          if (length(diag_map)) {
            if (!is.null(res$order)) {
              diag_map <- diag_map[res$order]
            }
            if (!is.null(colnames(M))) {
              names(diag_map) <- colnames(M)
            }
          }
          status_messages <- c(status_messages, res$message)
        }
      }
      thr <- input$threshold %||% 0
      if (thr > 0) {
        suppress <- abs(M) < thr
        M[suppress] <- NA_real_
      }
      if (isTRUE(input$mask_diag)) {
        diag(M) <- NA_real_
      } else if (length(diag_map)) {
        if (!is.null(colnames(M))) {
          diag(M) <- diag_map[colnames(M)]
        } else {
          diag(M) <- diag_map
        }
      }
      if (isTRUE(input$show_values) && prod(dim(M)) > .mc_heatmap_label_limit(use_plotly = use_plotly)) {
        status_messages <- c(
          status_messages,
          "Cell labels are hidden automatically for larger heatmaps to keep the viewer responsive."
        )
      }
      list(
        matrix = M,
        meta = info,
        status_messages = status_messages[!is.na(status_messages) & nzchar(status_messages)]
      )
    })

    output$meta <- shiny::renderUI({
      res <- filtered_matrix()
      info <- res$meta
      dims <- dim(res$matrix)
      total_vars <- ncol(info$matrix)
      .mc_meta_grid(c(
        "Result" = info$label,
        "Object class" = info$class,
        "Displayed variables" = sprintf("%d of %d", dims[[1L]], total_vars),
        "Signed values" = if (info$signed) "Yes" else "No"
      ))
    })

    output$status_alerts <- shiny::renderUI({
      .mc_status_notices(filtered_matrix()$status_messages)
    })

    output$download_matrix <- shiny::downloadHandler(
      filename = function() {
        sprintf(
          "%s-displayed-matrix.csv",
          .mc_sanitize_filename(input$matrix_choice %||% "matrixCorr")
        )
      },
      content = function(file) {
        utils::write.csv(filtered_matrix()$matrix, file = file, row.names = TRUE)
      }
    )

    plot_data <- shiny::reactive({
      res <- filtered_matrix()
      mat <- res$matrix
      signed <- res$meta$signed
      list(
        plot = .mc_build_heatmap(
          mat = mat,
          signed = signed,
          show_values = isTRUE(input$show_values),
          use_abs = isTRUE(input$use_abs),
          use_plotly = use_plotly
        ),
        signed = signed
      )
    })

    distribution_values <- shiny::reactive({
      .mc_corr_distribution_values(
        mat = filtered_matrix()$matrix,
        use_abs = isTRUE(input$use_abs)
      )
    })

    if (use_plotly) {
      output$heatmap <- .mc_plotly_fn("renderPlotly")({
        plot_data()$plot
      })
    } else {
      output$heatmap <- shiny::renderPlot({
        plot_data()$plot
      })
    }

    output$dist_plot <- shiny::renderPlot({
      values <- distribution_values()
      shiny::validate(shiny::need(
        length(values) > 0L,
        "No off-diagonal correlations are visible after filtering."
      ))
      .mc_build_corr_boxplot(values = values, use_abs = isTRUE(input$use_abs))
    })

    output$dist_summary <- shiny::renderTable(
      {
        .mc_corr_summary_table(distribution_values())
      },
      striped = TRUE,
      bordered = FALSE,
      spacing = "s",
      width = "100%",
      rownames = FALSE
    )
  }

  shiny::runApp(shiny::shinyApp(ui = ui, server = server))
  invisible(NULL)
}

.mc_default_vars <- function(mat, max_vars) {
  vars <- colnames(mat)
  if (is.null(vars)) {
    vars <- paste0("V", seq_len(ncol(mat)))
  }
  vars[seq_len(min(length(vars), max_vars))]
}

.mc_prepare_corr_inputs <- function(x) {
  parsed_single <- try(.mc_parse_corr_object(x, label = "default"), silent = TRUE)
  if (!inherits(parsed_single, "try-error")) {
    return(list(default = parsed_single))
  }

  if (!(is.list(x) && !is.matrix(x) && !inherits(x, "Matrix"))) {
    return(list())
  }

  objects <- x
  if (is.null(names(objects))) {
    names(objects) <- paste0("Matrix ", seq_along(objects))
  }

  out <- list()
  for (nm in names(objects)) {
    obj <- objects[[nm]]
    parsed <- try(.mc_parse_corr_object(obj, label = nm), silent = TRUE)
    if (!inherits(parsed, "try-error")) {
      out[[parsed$name]] <- parsed
    }
  }
  out
}

.mc_is_corr_result_like <- function(x) {
  if (inherits(x, "corr_result")) {
    return(TRUE)
  }
  if (isTRUE(attr(x, "corr_result", exact = TRUE))) {
    return(TRUE)
  }
  output_class <- attr(x, "corr_output_class", exact = TRUE)
  if (is.character(output_class) && length(output_class) == 1L && nzchar(output_class)) {
    return(TRUE)
  }
  output <- attr(x, "output", exact = TRUE)
  is.character(output) && length(output) == 1L &&
    output %in% c("matrix", "sparse", "edge_list", "packed_upper")
}

.mc_parse_corr_object <- function(obj, label) {
  mat <- NULL
  signed <- TRUE
  desc <- NULL
  cls <- class(obj)
  if (inherits(obj, "partial_corr")) {
    mat <- obj$pcor
    desc <- attr(obj, "method", exact = TRUE) %||% "partial correlation"
  } else if (.mc_is_corr_result_like(obj)) {
    mat <- .mc_corr_as_dense_matrix(obj)
    desc <- attr(obj, "description", exact = TRUE) %||%
      attr(obj, "method", exact = TRUE) %||%
      label
  } else if (inherits(obj, "Matrix")) {
    mat <- as.matrix(obj)
  } else if (is.matrix(obj)) {
    mat <- obj
  } else if (is.list(obj) && !is.null(obj$pcor) && is.matrix(obj$pcor)) {
    mat <- obj$pcor
    desc <- "partial correlation"
  } else if (is.list(obj) && !is.null(obj$est) && is.matrix(obj$est)) {
    mat <- obj$est
    desc <- attr(obj, "description", exact = TRUE) %||%
      attr(obj, "method", exact = TRUE) %||%
      label
  } else {
    stop("Unsupported object class", call. = FALSE)
  }
  mat <- as.matrix(mat)
  if (is.null(colnames(mat))) {
    colnames(mat) <- rownames(mat) <- paste0("V", seq_len(ncol(mat)))
  }
  if (nrow(mat) != ncol(mat)) {
    stop("Correlation matrices must be square.", call. = FALSE)
  }
  if (is.null(desc)) {
    desc <- attr(obj, "description", exact = TRUE) %||%
      attr(obj, "method", exact = TRUE) %||% label
  }
  signed <- any(mat < 0, na.rm = TRUE)
  list(
    name = label,
    matrix = mat,
    signed = signed,
    label = desc,
    class = paste(cls, collapse = ", ")
  )
}

.mc_reorder_matrix <- function(mat, mode = c("abs", "signed"), method = "complete") {
  mode <- match.arg(mode, c("abs", "signed"))
  allowed_methods <- c(
    complete = "complete",
    average  = "average",
    single   = "single",
    ward.d   = "ward.D",
    ward.d2  = "ward.D2",
    mcquitty = "mcquitty",
    median   = "median",
    centroid = "centroid"
  )
  method_key <- match.arg(tolower(method), names(allowed_methods))
  method <- allowed_methods[[method_key]]

  base_mat <- pmax(pmin(mat, 1), -1)
  dist_mat <- switch(
    mode,
    abs = 1 - abs(base_mat),
    signed = 0.5 * (1 - base_mat)
  )
  dist_mat <- (dist_mat + t(dist_mat)) / 2
  dist_mat[dist_mat < 0] <- 0
  dist_mat[dist_mat > 1] <- 1
  diag(dist_mat) <- 0
  dist_mat[!is.finite(dist_mat)] <- 1

  dist_obj <- try(suppressWarnings(stats::as.dist(dist_mat)), silent = TRUE)
  if (inherits(dist_obj, "try-error") || is.null(dist_obj)) {
    return(list(
      matrix = mat,
      message = "Unable to compute a distance matrix for clustering.",
      order = seq_len(nrow(mat))
    ))
  }

  hc <- try(stats::hclust(dist_obj, method = method), silent = TRUE)
  if (inherits(hc, "try-error")) {
    cond <- attr(hc, "condition")
    msg <- if (inherits(cond, "condition")) {
      conditionMessage(cond)
    } else {
      as.character(hc)
    }
    return(list(
      matrix = mat,
      message = paste("Clustering failed:", msg),
      order = seq_len(nrow(mat))
    ))
  }

  ord <- hc$order
  list(matrix = mat[ord, ord, drop = FALSE], message = NULL, order = ord)
}

.mc_filter_matrix_by_range <- function(mat, corr_range = NULL, use_abs = FALSE) {
  corr_range <- .mc_clamp_corr_range(corr_range, use_abs = use_abs)
  values <- if (isTRUE(use_abs)) abs(mat) else mat
  diag(values) <- NA_real_

  keep <- apply(
    values,
    2,
    function(x) any(!is.na(x) & x >= corr_range[[1L]] & x <= corr_range[[2L]])
  )
  if (!any(keep)) {
    return(list(
      matrix = mat[FALSE, FALSE, drop = FALSE],
      kept = character(),
      message = sprintf(
        "No variables have %s within [%s, %s].",
        if (isTRUE(use_abs)) "|correlation|" else "correlation",
        format(signif(corr_range[[1L]], digits = 3)),
        format(signif(corr_range[[2L]], digits = 3))
      )
    ))
  }

  kept <- colnames(mat)[keep]
  message <- NULL
  if (sum(!keep) > 0L) {
    message <- sprintf(
      "Correlation range filter kept %d of %d variables using %s within [%s, %s].",
      sum(keep),
      length(keep),
      if (isTRUE(use_abs)) "|correlation|" else "correlation",
      format(signif(corr_range[[1L]], digits = 3)),
      format(signif(corr_range[[2L]], digits = 3))
    )
  }

  list(
    matrix = mat[keep, keep, drop = FALSE],
    kept = kept,
    message = message
  )
}

.mc_clamp_corr_range <- function(value, use_abs = FALSE) {
  bounds <- if (isTRUE(use_abs)) c(0, 1) else c(-1, 1)
  if (is.null(value) || length(value) != 2L || any(!is.finite(value))) {
    return(bounds)
  }
  value <- sort(value)
  c(
    max(bounds[[1L]], min(bounds[[2L]], value[[1L]])),
    max(bounds[[1L]], min(bounds[[2L]], value[[2L]]))
  )
}

.mc_update_corr_range_slider <- function(session, input_id, use_abs = FALSE, value = NULL) {
  bounds <- if (isTRUE(use_abs)) c(0, 1) else c(-1, 1)
  label <- if (isTRUE(use_abs)) {
    "Keep variables with |correlation| in range"
  } else {
    "Keep variables with correlation in range"
  }
  shiny::updateSliderInput(
    session = session,
    inputId = input_id,
    label = label,
    min = bounds[[1L]],
    max = bounds[[2L]],
    value = .mc_clamp_corr_range(value, use_abs = use_abs),
    step = 0.05
  )
}

.mc_build_heatmap <- function(mat, signed, show_values, use_abs, use_plotly,
                              plotly_source = NULL) {
  if (use_abs) {
    plot_mat <- abs(mat)
    signed <- FALSE
  } else {
    plot_mat <- mat
  }
  show_values <- isTRUE(show_values) && prod(dim(plot_mat)) <= .mc_heatmap_label_limit(use_plotly = use_plotly)
  if (use_plotly) {
    return(.mc_build_heatmap_plotly(
      mat = plot_mat,
      signed = signed,
      show_values = show_values,
      plotly_source = plotly_source
    ))
  }
  df <- .mc_heatmap_df(plot_mat)
  fill_scale <- if (signed) {
    ggplot2::scale_fill_gradient2(
      low = "#FF6A6A",
      mid = "white",
      high = "#63B8FF",
      limits = c(-1, 1),
      midpoint = 0,
      na.value = "grey90",
      name = "r"
    )
  } else {
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "steelblue1",
      limits = c(0, 1),
      na.value = "grey90",
      name = "|r|"
    )
  }
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$Var2, y = .data$Var1, fill = .data$Value)
  ) +
    ggplot2::geom_tile(color = "white") +
    fill_scale +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = NULL, y = NULL)
  if (show_values) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(
        label = ifelse(is.na(.data$Value), "", sprintf("%.2f", .data$Value))
      ),
      size = 3,
      colour = "black"
    )
  }
  p
}

.mc_heatmap_df <- function(mat) {
  row_levels <- rownames(mat)
  if (is.null(row_levels)) {
    row_levels <- sprintf("V%02d", seq_len(nrow(mat)))
    rownames(mat) <- row_levels
  }
  col_levels <- colnames(mat)
  if (is.null(col_levels)) {
    col_levels <- sprintf("V%02d", seq_len(ncol(mat)))
    colnames(mat) <- col_levels
  }
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(df) <- c("Var1", "Var2", "Value")
  df$Var1 <- factor(df$Var1, levels = rev(row_levels))
  df$Var2 <- factor(df$Var2, levels = col_levels)
  df
}

.mc_corr_distribution_values <- function(mat, use_abs = FALSE) {
  if (!is.matrix(mat) || nrow(mat) < 2L || ncol(mat) < 2L) {
    return(numeric())
  }
  values <- if (isTRUE(use_abs)) abs(mat) else mat
  values <- values[upper.tri(values, diag = FALSE)]
  values[is.finite(values)]
}

.mc_corr_summary_table <- function(values) {
  if (!length(values)) {
    return(data.frame(
      Metric = c("Pairs shown", "Status"),
      Value = c("0", "No off-diagonal correlations visible"),
      stringsAsFactors = FALSE
    ))
  }

  quants <- stats::quantile(values, probs = c(0.25, 0.5, 0.75), names = FALSE, na.rm = TRUE)
  data.frame(
    Metric = c("Pairs shown", "Min", "Q1", "Median", "Mean", "Q3", "Max", "SD"),
    Value = c(
      as.character(length(values)),
      .mc_format_metric(min(values)),
      .mc_format_metric(quants[[1L]]),
      .mc_format_metric(quants[[2L]]),
      .mc_format_metric(mean(values)),
      .mc_format_metric(quants[[3L]]),
      .mc_format_metric(max(values)),
      .mc_format_metric(stats::sd(values))
    ),
    stringsAsFactors = FALSE
  )
}

.mc_format_metric <- function(x, digits = 3L) {
  if (!length(x) || anyNA(x) || !is.finite(x)) {
    return("NA")
  }
  format(signif(x, digits = digits), trim = TRUE)
}

.mc_corr_axis_limits <- function(values, use_abs = FALSE) {
  bounds <- if (isTRUE(use_abs)) c(0, 1) else c(-1, 1)
  if (!length(values) || all(!is.finite(values))) {
    return(bounds)
  }

  rng <- range(values, na.rm = TRUE, finite = TRUE)
  span <- diff(rng)
  pad <- if (span <= 0) {
    max(0.04, abs(rng[[1L]]) * 0.12)
  } else {
    max(0.03, span * 0.12)
  }
  limits <- c(rng[[1L]] - pad, rng[[2L]] + pad)
  limits[[1L]] <- max(bounds[[1L]], limits[[1L]])
  limits[[2L]] <- min(bounds[[2L]], limits[[2L]])

  if (diff(limits) < 0.08) {
    centre <- mean(limits)
    limits <- c(centre - 0.04, centre + 0.04)
    limits[[1L]] <- max(bounds[[1L]], limits[[1L]])
    limits[[2L]] <- min(bounds[[2L]], limits[[2L]])
  }

  limits
}

.mc_build_corr_boxplot <- function(values, use_abs = FALSE) {
  df <- data.frame(Group = "Displayed pairs", Value = values)
  limits <- .mc_corr_axis_limits(values = values, use_abs = use_abs)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Group, y = .data$Value))
  if (!isTRUE(use_abs)) {
    p <- p + ggplot2::geom_hline(
      yintercept = 0,
      colour = "#C8D2DD",
      linewidth = 0.6,
      linetype = "dashed"
    )
  }
  p <- p +
    ggplot2::geom_boxplot(
      width = 0.22,
      outlier.shape = NA,
      fill = if (isTRUE(use_abs)) "#A7D8FF" else "#C6D8FF",
      colour = "#153A5B",
      linewidth = 0.8
    ) +
    ggplot2::geom_jitter(
      width = 0.08,
      height = 0,
      size = 1.7,
      alpha = 0.3,
      colour = if (isTRUE(use_abs)) "#1874CD" else "#24557A"
    ) +
    ggplot2::coord_cartesian(ylim = limits) +
    ggplot2::labs(
      x = NULL,
      y = if (isTRUE(use_abs)) "|r|" else "r"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(colour = "#274055", face = "bold"),
      axis.text.y = ggplot2::element_text(colour = "#425466"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
  p
}

.mc_plotly_click_pair <- function(event, vars = NULL) {
  if (is.null(event) || !NROW(event)) {
    return(NULL)
  }
  x <- event[["x"]][[1L]]
  y <- event[["y"]][[1L]]
  if (is.null(x) || is.null(y)) {
    return(NULL)
  }
  pair <- c(as.character(x), as.character(y))
  if (!is.null(vars) && !all(pair %in% vars)) {
    return(NULL)
  }
  pair
}

.mc_static_click_pair <- function(df, click, threshold = 15) {
  if (is.null(click) || is.null(df) || !nrow(df)) {
    return(NULL)
  }
  if (is.character(click[["x"]]) && is.character(click[["y"]])) {
    pair <- c(as.character(click[["x"]]), as.character(click[["y"]]))
    vars <- unique(c(as.character(df$Var1), as.character(df$Var2)))
    if (all(pair %in% vars)) {
      return(pair)
    }
  }
  hit <- shiny::nearPoints(
    df,
    click,
    xvar = "Var2",
    yvar = "Var1",
    threshold = threshold,
    maxpoints = 1,
    addDist = FALSE
  )
  if (!nrow(hit)) {
    return(NULL)
  }
  c(as.character(hit$Var2[[1L]]), as.character(hit$Var1[[1L]]))
}

.mc_heatmap_label_limit <- function(use_plotly = FALSE) {
  if (isTRUE(use_plotly)) 196L else 400L
}

.mc_build_heatmap_plotly <- function(mat, signed, show_values, plotly_source = NULL) {
  cols <- colnames(mat)
  if (is.null(cols)) {
    cols <- sprintf("V%02d", seq_len(ncol(mat)))
  }
  rows <- rownames(mat)
  if (is.null(rows)) {
    rows <- sprintf("V%02d", seq_len(nrow(mat)))
  }
  z <- mat[nrow(mat):1L, , drop = FALSE]
  y <- rev(rows)
  text_mat <- if (isTRUE(show_values)) {
    matrix(ifelse(is.na(z), "", sprintf("%.2f", z)), nrow = nrow(z), ncol = ncol(z))
  } else {
    NULL
  }
  plot_ly <- .mc_plotly_fn("plot_ly")
  layout_fn <- .mc_plotly_fn("layout")
  p <- plot_ly(
    x = cols,
    y = y,
    z = z,
    type = "heatmap",
    source = plotly_source,
    colorscale = .mc_plotly_colorscale(signed),
    reversescale = FALSE,
    zmin = if (signed) -1 else 0,
    zmax = 1,
    hoverongaps = FALSE,
    text = text_mat,
    texttemplate = if (isTRUE(show_values)) "%{text}" else NULL,
    textfont = if (isTRUE(show_values)) list(color = "black", size = 11) else NULL,
    hovertemplate = if (signed) {
      paste("x: %{x}<br>y: %{y}<br>r: %{z:.3f}<extra></extra>")
    } else {
      paste("x: %{x}<br>y: %{y}<br>|r|: %{z:.3f}<extra></extra>")
    },
    colorbar = list(title = if (signed) "r" else "|r|")
  )
  layout_fn(
    p,
    margin = list(l = 70, r = 30, b = 90, t = 20),
    xaxis = list(title = NULL, tickangle = 45),
    yaxis = list(title = NULL),
    paper_bgcolor = "white",
    plot_bgcolor = "white"
  )
}

.mc_plotly_colorscale <- function(signed) {
  if (isTRUE(signed)) {
    list(
      list(0.0, "#FF6A6A"),
      list(0.5, "white"),
      list(1.0, "#63B8FF")
    )
  } else {
    list(
      list(0.0, "white"),
      list(1.0, "#63B8FF")
    )
  }
}

.mc_has_namespace <- function(pkg) {
  !inherits(try(getNamespace(pkg), silent = TRUE), "try-error")
}

.mc_plotly_fn <- function(name) {
  if (!.mc_has_namespace("plotly")) {
    cli::cli_abort("Package {.pkg plotly} is required for {.arg use_plotly} = TRUE.")
  }
  getExportedValue("plotly", name)
}

.mc_viewer_header <- function(title, subtitle = NULL, badges = NULL, logo_src = NULL) {
  shiny::tags$div(
    class = "mc-hero",
    shiny::tags$div(
      class = "mc-hero-copy",
      shiny::tags$div(
        class = "mc-hero-title-row",
        if (!is.null(logo_src) && nzchar(logo_src)) shiny::tags$img(
          class = "mc-hero-logo",
          src = logo_src,
          alt = "matrixCorr logo"
        ),
        shiny::tags$div(class = "mc-hero-title", title)
      ),
      if (!is.null(subtitle) && nzchar(subtitle)) shiny::tags$p(class = "mc-hero-subtitle", subtitle)
    ),
    if (length(badges)) shiny::tags$div(class = "mc-hero-badges", badges)
  )
}

.mc_header_badge <- function(text, tone = c("neutral", "accent")) {
  tone <- match.arg(tone)
  shiny::tags$span(class = paste("mc-header-badge", paste0("mc-header-badge-", tone)), text)
}

.mc_control_section <- function(title, ...) {
  shiny::tags$div(
    class = "mc-control-section",
    shiny::tags$div(class = "mc-section-title", title),
    ...
  )
}

.mc_card_header <- function(title, subtitle = NULL) {
  shiny::tags$div(
    class = "mc-card-header",
    shiny::tags$div(class = "mc-card-title", title),
    if (!is.null(subtitle) && nzchar(subtitle)) shiny::tags$div(class = "mc-card-subtitle", subtitle)
  )
}

.mc_meta_grid <- function(items) {
  shiny::tags$div(
    class = "mc-meta-grid",
    lapply(seq_along(items), function(i) {
      shiny::tags$div(
        class = "mc-stat-card",
        shiny::tags$div(class = "mc-stat-label", names(items)[[i]]),
        shiny::tags$div(class = "mc-stat-value", as.character(items[[i]]))
      )
    })
  )
}

.mc_status_notices <- function(msgs) {
  if (is.null(msgs) || !length(msgs)) {
    return(NULL)
  }
  shiny::tagList(lapply(msgs, function(msg) {
    shiny::tags$div(
      class = "mc-status-note",
      shiny::tags$span(class = "mc-status-dot"),
      shiny::tags$span(msg)
    )
  }))
}

.mc_viewer_head <- function() {
  shiny::tags$head(
    shiny::tags$style(shiny::HTML(
      paste(
        "body { background: linear-gradient(180deg, #F3F7FB 0%, #ECF2F8 100%); color: #162535; }",
        ".container-fluid { max-width: 1580px; padding: 24px 24px 112px; }",
        ".mc-app-shell { position: relative; }",
        ".mc-hero { display: flex; justify-content: space-between; align-items: flex-start; gap: 24px; margin-bottom: 22px; padding: 28px 30px; border-radius: 26px; background: linear-gradient(135deg, #12263F 0%, #1C4E72 58%, #2B7A8E 100%); box-shadow: 0 18px 42px rgba(17, 37, 62, 0.18); color: #FFFFFF; }",
        ".mc-hero-copy { max-width: 860px; }",
        ".mc-hero-title-row { display: flex; align-items: center; gap: 16px; margin-bottom: 8px; }",
        ".mc-hero-logo { width: 74px; height: auto; display: block; flex: 0 0 auto; filter: drop-shadow(0 8px 18px rgba(5, 14, 24, 0.22)); }",
        ".mc-hero-title { font-size: 30px; line-height: 1.15; font-weight: 700; letter-spacing: -0.02em; margin: 0; }",
        ".mc-hero-subtitle { margin: 0; max-width: 760px; font-size: 14px; line-height: 1.6; color: rgba(255,255,255,0.82); }",
        ".mc-hero-badges { display: flex; flex-wrap: wrap; gap: 10px; justify-content: flex-end; }",
        ".mc-header-badge { display: inline-flex; align-items: center; padding: 8px 12px; border-radius: 999px; font-size: 12px; font-weight: 700; letter-spacing: 0.04em; text-transform: uppercase; }",
        ".mc-header-badge-neutral { background: rgba(255,255,255,0.16); color: #FFFFFF; border: 1px solid rgba(255,255,255,0.22); }",
        ".mc-header-badge-accent { background: #E7F5FA; color: #14405E; }",
        ".well { background: rgba(255,255,255,0.94); border: 1px solid #D9E4EE; border-radius: 24px; box-shadow: 0 14px 32px rgba(21, 37, 56, 0.08); padding: 22px; }",
        ".mc-control-section + .mc-control-section { margin-top: 20px; padding-top: 18px; border-top: 1px solid #E8EEF5; }",
        ".mc-section-title { margin: 0 0 12px; font-size: 11px; font-weight: 800; letter-spacing: 0.08em; text-transform: uppercase; color: #587086; }",
        ".mc-card { margin-bottom: 18px; padding: 22px 24px; border-radius: 24px; border: 1px solid #D9E4EE; background: rgba(255,255,255,0.96); box-shadow: 0 14px 34px rgba(21, 37, 56, 0.08); }",
        ".mc-card-header { margin-bottom: 16px; }",
        ".mc-card-title { font-size: 18px; font-weight: 700; color: #15283D; }",
        ".mc-card-subtitle { margin-top: 4px; font-size: 13px; line-height: 1.5; color: #60758A; }",
        ".mc-tabset-shell .nav-tabs { border-bottom: 0; margin-bottom: 18px; display: flex; gap: 10px; }",
        ".mc-tabset-shell .nav-tabs > li { float: none; }",
        ".mc-tabset-shell .nav-tabs > li > a { margin-right: 0; border: 0; border-radius: 999px; background: rgba(255,255,255,0.72); color: #33516A; font-weight: 700; padding: 10px 16px; }",
        ".mc-tabset-shell .nav-tabs > li.active > a, .mc-tabset-shell .nav-tabs > li.active > a:focus, .mc-tabset-shell .nav-tabs > li.active > a:hover { border: 0; background: #153A5B; color: #FFFFFF; box-shadow: 0 10px 20px rgba(21, 58, 91, 0.18); }",
        ".mc-meta-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 14px; }",
        ".mc-stat-card { padding: 16px 18px; border-radius: 18px; background: linear-gradient(180deg, #F9FBFD 0%, #F1F6FA 100%); border: 1px solid #E3EBF3; }",
        ".mc-stat-label { font-size: 11px; font-weight: 800; letter-spacing: 0.08em; text-transform: uppercase; color: #6E8193; }",
        ".mc-stat-value { margin-top: 8px; font-size: 16px; line-height: 1.4; font-weight: 700; color: #17324A; word-break: break-word; }",
        ".mc-status-note { display: flex; align-items: flex-start; gap: 10px; margin-top: 12px; padding: 12px 14px; border-radius: 16px; background: #FFF8E8; border: 1px solid #F2DEAE; color: #6F5515; }",
        ".mc-status-dot { width: 9px; height: 9px; border-radius: 50%; margin-top: 5px; flex: 0 0 auto; background: #D9A321; }",
        ".mc-viewer-controls { display: flex; flex-wrap: wrap; gap: 10px; }",
        ".mc-viewer-controls .btn, .mc-viewer-controls .btn-default { margin: 0; border-radius: 12px; border: 0; padding: 10px 14px; font-weight: 700; box-shadow: none; }",
        ".mc-viewer-controls .btn-default { background: #153A5B; color: #FFFFFF; }",
        ".mc-viewer-controls .btn-default:hover, .mc-viewer-controls .btn-default:focus { background: #1D4C77; color: #FFFFFF; }",
        ".bootstrap-select > .dropdown-toggle, .form-control { min-height: 44px; border-radius: 14px; border-color: #D8E2EC; box-shadow: none; }",
        ".bootstrap-select > .dropdown-toggle { background: #FFFFFF; }",
        ".control-label { font-size: 12px; font-weight: 700; letter-spacing: 0.02em; color: #274055; }",
        ".checkbox label { color: #274055; font-weight: 600; }",
        ".irs--shiny .irs-bar, .irs--shiny .irs-from, .irs--shiny .irs-to, .irs--shiny .irs-single { background: #1C5A7B; border-top-color: #1C5A7B; border-bottom-color: #1C5A7B; }",
        ".irs--shiny .irs-handle > i:first-child { background: #153A5B; }",
        ".table { margin-bottom: 0; }",
        ".table > thead > tr > th { border-bottom: 1px solid #D9E4EE; font-size: 11px; font-weight: 800; letter-spacing: 0.08em; text-transform: uppercase; color: #6E8193; }",
        ".table > tbody > tr > td { border-top: 1px solid #EDF2F7; color: #1C3249; font-weight: 600; }",
        "@media (max-width: 768px) {",
        "  .container-fluid { padding: 16px 14px 84px; }",
        "  .mc-hero { padding: 22px 18px; border-radius: 22px; }",
        "  .mc-hero-title-row { align-items: flex-start; }",
        "  .mc-hero-logo { width: 62px; }",
        "  .mc-hero-title { font-size: 24px; }",
        "  .mc-hero-badges { justify-content: flex-start; }",
        "  .mc-card { padding: 18px 16px; border-radius: 20px; }",
        "}",
        sep = "\n"
      )
    ))
  )
}

.mc_viewer_footer <- function(logo_src) {
  if (is.null(logo_src) || !nzchar(logo_src)) {
    return(NULL)
  }
  shiny::tags$div(
    class = "mc-viewer-footer",
    shiny::tags$img(src = logo_src, alt = "matrixCorr logo"),
    shiny::tags$div(class = "mc-viewer-footer-text", "matrixCorr")
  )
}

.mc_logo_path <- function() {
  candidates <- c(
    system.file("app-www", "logo.svg", package = "matrixCorr"),
    file.path("inst", "app-www", "logo.svg"),
    file.path("man", "figures", "logo.svg")
  )
  candidates <- candidates[nzchar(candidates)]
  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) {
    return(NULL)
  }
  normalizePath(existing[[1L]], winslash = "/", mustWork = FALSE)
}

.mc_register_logo_resource <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  prefix <- paste0("matrixcorr-assets-", as.integer(stats::runif(1L, min = 1, max = 1e9)))
  shiny::addResourcePath(prefix, dirname(path))
  sprintf("%s/%s", prefix, basename(path))
}

.mc_sanitize_filename <- function(x) {
  out <- gsub("[^A-Za-z0-9_-]+", "-", x)
  out <- gsub("(^-+|-+$)", "", out)
  if (!nzchar(out)) "matrixCorr" else out
}
