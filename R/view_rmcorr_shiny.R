#' Interactive Shiny Viewer for Repeated-Measures Correlation
#'
#' Launches a dedicated Shiny gadget for repeated-measures correlation matrix
#' objects of class \code{"rmcorr_matrix"}. The viewer combines the correlation
#' heatmap with a pairwise scatterplot panel that rebuilds the corresponding
#' two-variable \code{"rmcorr"} fit for user-selected variables.
#'
#' @param x An object of class \code{"rmcorr_matrix"} or a named list of such
#'   objects.
#' @param title Optional character title shown at the top of the gadget.
#' @param default_max_vars Integer; maximum number of variables pre-selected
#'   in the heatmap view when the app opens. Defaults to 40.
#'
#' @return Invisibly returns \code{NULL}; the function is called for its side
#'   effect of launching a Shiny gadget.
#'
#' @details This helper requires the \pkg{shiny} and \pkg{shinyWidgets}
#' packages at runtime and will optionally use \pkg{plotly} for the heatmap
#' when available. The pairwise panel reuses the package's regular
#' \code{plot.rmcorr()} method, so the Shiny scatterplot matches the standard
#' pairwise repeated-measures correlation plot. To rebuild pairwise fits from a
#' returned \code{"rmcorr_matrix"} object, the matrix must have been created
#' with \code{keep_data = TRUE}.
#'
#' @examples
#' if (interactive()) {
#'   set.seed(2026)
#'   n_subjects <- 20
#'   n_rep <- 4
#'   subject <- rep(seq_len(n_subjects), each = n_rep)
#'   subj_eff_x <- rnorm(n_subjects, sd = 1.5)
#'   subj_eff_y <- rnorm(n_subjects, sd = 2.0)
#'   within_signal <- rnorm(n_subjects * n_rep)
#'
#'   dat <- data.frame(
#'     subject = subject,
#'     x = subj_eff_x[subject] + within_signal + rnorm(n_subjects * n_rep, sd = 0.2),
#'     y = subj_eff_y[subject] + 0.8 * within_signal + rnorm(n_subjects * n_rep, sd = 0.3),
#'     z = subj_eff_y[subject] - 0.4 * within_signal + rnorm(n_subjects * n_rep, sd = 0.4)
#'   )
#'
#'   fit_mat <- rmcorr(
#'     dat,
#'     response = c("x", "y", "z"),
#'     subject = "subject",
#'     keep_data = TRUE
#'   )
#'   view_rmcorr_shiny(fit_mat)
#' }
#'
#' @export
view_rmcorr_shiny <- function(x, title = NULL, default_max_vars = 40L) {
  if (missing(x)) {
    abort_bad_arg("x", message = "must be supplied (an rmcorr_matrix result or list).")
  }

  if (!requireNamespace("shiny", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg shiny} is required for {.fn view_rmcorr_shiny}.")
  }
  if (!requireNamespace("shinyWidgets", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg shinyWidgets} is required for {.fn view_rmcorr_shiny}.")
  }

  prepared <- .mc_prepare_rmcorr_inputs(x)
  if (!length(prepared)) {
    cli::cli_abort("{.fn view_rmcorr_shiny} did not find a usable {.cls rmcorr_matrix} object.")
  }

  use_plotly <- .mc_has_namespace("plotly")
  app_title <- title %||% "matrixCorr repeated-measures correlation viewer"
  logo_src <- .mc_register_logo_resource(.mc_logo_path())

  heatmap_widget <- if (use_plotly) {
    .mc_plotly_fn("plotlyOutput")("heatmap", height = "650px")
  } else {
    shiny::plotOutput("heatmap", click = "heatmap_click", height = "650px")
  }

  ui <- shiny::fluidPage(
    .mc_viewer_head(),
    shiny::div(
      class = "mc-app-shell",
      .mc_viewer_header(
        title = app_title,
        subtitle = "Inspect repeated-measures correlation matrices with linked heatmap, pairwise fit, and distribution summaries.",
        logo_src = logo_src,
        badges = list(
          .mc_header_badge("Repeated Measures", tone = "neutral"),
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
              label = "Repeated-measures correlation object",
              choices = names(prepared),
              selected = names(prepared)[[1L]],
              multiple = FALSE,
              options = list(`live-search` = TRUE)
            ),
            shinyWidgets::pickerInput(
              inputId = "var_choice",
              label = "Heatmap variables",
              choices = colnames(prepared[[1L]]$matrix),
              selected = .mc_default_vars(prepared[[1L]]$matrix, default_max_vars),
              multiple = TRUE,
              options = list(`live-search` = TRUE, `actions-box` = TRUE)
            )
          ),
          .mc_control_section(
            title = "Heatmap",
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
            shiny::checkboxInput("show_values", "Show heatmap cell labels", value = TRUE)
          ),
          .mc_control_section(
            title = "Pair Plot",
            shiny::selectInput(
              inputId = "pair_x",
              label = "Pair plot X variable",
              choices = colnames(prepared[[1L]]$matrix),
              selected = colnames(prepared[[1L]]$matrix)[[1L]]
            ),
            shiny::selectInput(
              inputId = "pair_y",
              label = "Pair plot Y variable",
              choices = colnames(prepared[[1L]]$matrix),
              selected = colnames(prepared[[1L]]$matrix)[[2L]]
            ),
            shiny::checkboxInput("pair_show_legend", "Show subject legend", value = FALSE)
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
              id = "rmcorr_view",
              shiny::tabPanel(
                title = "Heatmap",
                shiny::div(
                  class = "mc-card mc-main-card",
                  .mc_card_header(
                    title = "Heatmap",
                    subtitle = "Click a cell to open the corresponding repeated-measures pair plot."
                  ),
                  heatmap_widget,
                  shiny::uiOutput("status_alerts")
                ),
                shiny::div(
                  class = "mc-card mc-meta-card",
                  .mc_card_header(
                    title = "Summary",
                    subtitle = "Current repeated-measures matrix and display state."
                  ),
                  shiny::uiOutput("heatmap_meta")
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
              ),
              shiny::tabPanel(
                title = "Pair Plot",
                shiny::div(
                  class = "mc-card mc-main-card",
                  .mc_card_header(
                    title = "Pair Plot",
                    subtitle = "Repeated-measures fit for the currently selected variable pair."
                  ),
                  shiny::plotOutput("pair_plot", height = "650px")
                ),
                shiny::div(
                  class = "mc-card mc-meta-card",
                  .mc_card_header(
                    title = "Pair Summary",
                    subtitle = "Model estimates for the selected repeated-measures pair."
                  ),
                  shiny::uiOutput("pair_meta")
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
      pair_default <- .mc_rmcorr_default_pair(vars)
      shinyWidgets::updatePickerInput(
        session,
        inputId = "var_choice",
        choices = vars,
        selected = .mc_default_vars(info$matrix, default_max_vars)
      )
      shiny::updateSelectInput(session, "pair_x", choices = vars, selected = pair_default[[1L]])
      shiny::updateSelectInput(session, "pair_y", choices = vars, selected = pair_default[[2L]])
      shiny::updateSliderInput(
        session,
        "threshold",
        min = 0,
        max = 1,
        value = min(input$threshold %||% 0, 1)
      )
    }, ignoreNULL = FALSE)

    shiny::observeEvent(input$reset_controls, {
      info <- current_info()
      vars <- colnames(info$matrix)
      pair_default <- .mc_rmcorr_default_pair(vars)
      shinyWidgets::updatePickerInput(
        session,
        inputId = "var_choice",
        choices = vars,
        selected = .mc_default_vars(info$matrix, default_max_vars)
      )
      shiny::updateSliderInput(session, "threshold", value = 0)
      shiny::updateCheckboxInput(session, "use_abs", value = FALSE)
      shiny::updateCheckboxInput(session, "mask_diag", value = FALSE)
      shiny::updateSelectInput(session, "cluster_mode", selected = "none")
      shiny::updateSelectInput(session, "cluster_method", selected = "complete")
      shiny::updateCheckboxInput(session, "show_values", value = TRUE)
      shiny::updateSelectInput(session, "pair_x", choices = vars, selected = pair_default[[1L]])
      shiny::updateSelectInput(session, "pair_y", choices = vars, selected = pair_default[[2L]])
      shiny::updateCheckboxInput(session, "pair_show_legend", value = FALSE)
      .mc_update_corr_range_slider(
        session = session,
        input_id = "corr_range",
        use_abs = FALSE,
        value = c(-1, 1)
      )
    })

    current_info <- shiny::reactive({
      prepared[[input$matrix_choice]]
    })

    filtered_heatmap <- shiny::reactive({
      info <- current_info()
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
          res <- .mc_reorder_matrix(M, mode = cluster_mode, method = cluster_method)
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

    pair_fit <- shiny::reactive({
      info <- current_info()
      vars <- colnames(info$matrix)
      x_var <- input$pair_x
      y_var <- input$pair_y
      shiny::validate(shiny::need(x_var %in% vars, "Choose a valid X variable."))
      shiny::validate(shiny::need(y_var %in% vars, "Choose a valid Y variable."))
      shiny::validate(shiny::need(!identical(x_var, y_var), "Choose two distinct variables."))
      .mc_rmcorr_view_pair_fit(info, x_var = x_var, y_var = y_var)
    })

    heatmap_click_pair <- shiny::reactive({
      if (use_plotly) {
        .mc_plotly_click_pair(
          .mc_plotly_fn("event_data")("plotly_click", source = "matrixcorr-rmcorr-heatmap"),
          vars = colnames(filtered_heatmap()$matrix)
        )
      } else {
        .mc_static_click_pair(
          .mc_heatmap_df(if (isTRUE(input$use_abs)) abs(filtered_heatmap()$matrix) else filtered_heatmap()$matrix),
          input$heatmap_click
        )
      }
    })

    shiny::observeEvent(heatmap_click_pair(), {
      pair <- heatmap_click_pair()
      if (is.null(pair) || length(pair) != 2L) {
        return()
      }
      shiny::updateSelectInput(session, "pair_x", selected = pair[[1L]])
      shiny::updateSelectInput(session, "pair_y", selected = pair[[2L]])
      shiny::updateTabsetPanel(session, "rmcorr_view", selected = "Pair Plot")
    }, ignoreInit = TRUE)

    output$heatmap_meta <- shiny::renderUI({
      res <- filtered_heatmap()
      info <- res$meta
      dims <- dim(res$matrix)
      total_vars <- ncol(info$matrix)
      diag_attr <- attr(info$matrix, "diagnostics", exact = TRUE)
      conf_level <- diag_attr$conf_level %||% info$conf_level
      .mc_meta_grid(c(
        "Result" = info$label,
        "Object class" = info$class,
        "Displayed variables" = sprintf("%d of %d", dims[[1L]], total_vars),
        "Subjects" = as.character(info$n_subjects),
        "Confidence level" = format(signif(conf_level, digits = 3))
      ))
    })

    output$status_alerts <- shiny::renderUI({
      .mc_status_notices(filtered_heatmap()$status_messages)
    })

    output$download_matrix <- shiny::downloadHandler(
      filename = function() {
        sprintf(
          "%s-displayed-matrix.csv",
          .mc_sanitize_filename(input$matrix_choice %||% "matrixCorr-rmcorr")
        )
      },
      content = function(file) {
        utils::write.csv(filtered_heatmap()$matrix, file = file, row.names = TRUE)
      }
    )

    heatmap_plot <- shiny::reactive({
      res <- filtered_heatmap()
      .mc_build_heatmap(
        mat = res$matrix,
        signed = TRUE,
        show_values = isTRUE(input$show_values),
        use_abs = isTRUE(input$use_abs),
        use_plotly = use_plotly,
        plotly_source = "matrixcorr-rmcorr-heatmap"
      )
    })

    distribution_values <- shiny::reactive({
      .mc_corr_distribution_values(
        mat = filtered_heatmap()$matrix,
        use_abs = isTRUE(input$use_abs)
      )
    })

    if (use_plotly) {
      output$heatmap <- .mc_plotly_fn("renderPlotly")({
        heatmap_plot()
      })
    } else {
      output$heatmap <- shiny::renderPlot({
        heatmap_plot()
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

    output$pair_plot <- shiny::renderPlot({
      print(plot(pair_fit(), show_legend = isTRUE(input$pair_show_legend)))
    })

    output$pair_meta <- shiny::renderUI({
      fit <- pair_fit()
      .mc_meta_grid(c(
        "Pair" = sprintf("%s vs %s", fit$responses[[1L]], fit$responses[[2L]]),
        "Correlation" = format(signif(fit$r, digits = 4)),
        "Slope" = format(signif(fit$slope, digits = 4)),
        "P-value" = format(signif(fit$p_value, digits = 4)),
        "Degrees of freedom" = format(signif(fit$df, digits = 4)),
        "Subjects" = as.character(fit$n_subjects),
        "Based on" = as.character(fit$based.on)
      ))
    })
  }

  shiny::runApp(shiny::shinyApp(ui = ui, server = server))
  invisible(NULL)
}

.mc_prepare_rmcorr_inputs <- function(x) {
  objects <-
    if (is.list(x) && !is.matrix(x) && !inherits(x, "Matrix") &&
        !inherits(x, "rmcorr_matrix")) {
      if (is.null(names(x))) {
        names(x) <- paste0("rmcorr ", seq_along(x))
      }
      x
    } else {
      list(default = x)
    }

  out <- list()
  for (nm in names(objects)) {
    obj <- objects[[nm]]
    parsed <- try(.mc_parse_rmcorr_object(obj, label = nm), silent = TRUE)
    if (!inherits(parsed, "try-error")) {
      out[[parsed$name]] <- parsed
    }
  }
  out
}

.mc_parse_rmcorr_object <- function(obj, label) {
  if (!inherits(obj, "rmcorr_matrix")) {
    stop("Unsupported object class", call. = FALSE)
  }

  mat <- as.matrix(obj)
  if (!is.matrix(mat) || nrow(mat) != ncol(mat)) {
    stop("Repeated-measures correlation matrices must be square.", call. = FALSE)
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- rownames(mat) <- paste0("V", seq_len(ncol(mat)))
  }

  source_data <- attr(obj, "source_data", exact = TRUE)
  if (!is.list(source_data)) {
    stop(
      "rmcorr_matrix objects need embedded source_data for pair plots. Refit with keep_data = TRUE.",
      call. = FALSE
    )
  }
  if (is.null(source_data$response) ||
      is.null(source_data$response_names) ||
      is.null(source_data$subject_code) ||
      is.null(source_data$subject_levels)) {
    stop("source_data is incomplete for rebuilding pair plots.", call. = FALSE)
  }

  response <- as.matrix(source_data$response)
  colnames(response) <- source_data$response_names
  if (!identical(colnames(response), colnames(mat))) {
    stop("source_data response columns do not match the matrix variables.", call. = FALSE)
  }
  if (nrow(response) != length(source_data$subject_code)) {
    stop("source_data subject length does not match the stored responses.", call. = FALSE)
  }

  list(
    name = label,
    matrix = mat,
    label = attr(obj, "description", exact = TRUE) %||% label,
    class = paste(class(obj), collapse = ", "),
    source_data = source_data,
    conf_level = attr(obj, "diagnostics", exact = TRUE)$conf_level %||% NA_real_,
    n_subjects = length(source_data$subject_levels)
  )
}

.mc_rmcorr_default_pair <- function(vars) {
  if (length(vars) < 2L) {
    stop("Repeated-measures correlation viewers require at least two variables.", call. = FALSE)
  }
  vars[1:2]
}

.mc_rmcorr_view_pair_fit <- function(info, x_var, y_var) {
  source_data <- info$source_data
  response <- as.matrix(source_data$response)
  colnames(response) <- source_data$response_names
  response_mat <- response[, c(x_var, y_var), drop = FALSE]
  subject <- factor(
    source_data$subject_levels[source_data$subject_code],
    levels = source_data$subject_levels
  )
  subj_info <- list(code = source_data$subject_code, levels = source_data$subject_levels)
  stats <- rmcorr_pair_cpp(
    x = response_mat[, 1L],
    y = response_mat[, 2L],
    subject = subj_info$code,
    conf_level = source_data$conf_level %||% info$conf_level %||% 0.95
  )
  .mc_rmcorr_build_pair_object(
    stats = stats,
    response_mat = response_mat,
    subject = subject,
    response_names = c(x_var, y_var),
    subject_name = source_data$subject_name %||% "subject",
    conf_level = source_data$conf_level %||% info$conf_level %||% 0.95
  )
}
