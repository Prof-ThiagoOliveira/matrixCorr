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

  use_plotly <- requireNamespace("plotly", quietly = TRUE)
  app_title <- title %||% "matrixCorr repeated-measures correlation viewer"

  heatmap_widget <- if (use_plotly) {
    .mc_plotly_fn("plotlyOutput")("heatmap", height = "650px")
  } else {
    shiny::plotOutput("heatmap", height = "650px")
  }

  ui <- shiny::fluidPage(
    shiny::titlePanel(app_title),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
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
        ),
        shiny::sliderInput(
          inputId = "threshold",
          label = "Hide cells with |value| below",
          min = 0,
          max = 1,
          step = 0.05,
          value = 0
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
        shiny::checkboxInput("show_values", "Show heatmap cell labels", value = TRUE),
        shiny::tags$hr(),
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
      shiny::mainPanel(
        shiny::tabsetPanel(
          id = "rmcorr_view",
          shiny::tabPanel(
            title = "Heatmap",
            heatmap_widget,
            shiny::uiOutput("cluster_alert"),
            shiny::htmlOutput("heatmap_meta")
          ),
          shiny::tabPanel(
            title = "Pair Plot",
            shiny::plotOutput("pair_plot", height = "650px"),
            shiny::htmlOutput("pair_meta")
          )
        )
      )
    )
  )

  server <- function(input, output, session) {
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
      diag_map <- diag(M)
      if (!is.null(colnames(M))) names(diag_map) <- colnames(M)
      cluster_mode <- input$cluster_mode %||% "none"
      cluster_method <- input$cluster_method %||% "complete"
      cluster_message <- NULL
      if (!identical(cluster_mode, "none")) {
        if (nrow(M) < 3L) {
          cluster_message <- "Clustering requires at least three variables."
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
          cluster_message <- res$message
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
      list(matrix = M, meta = info, cluster_message = cluster_message)
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

    output$heatmap_meta <- shiny::renderUI({
      info <- filtered_heatmap()$meta
      dims <- dim(info$matrix)
      diag_attr <- attr(info$matrix, "diagnostics", exact = TRUE)
      conf_level <- diag_attr$conf_level %||% info$conf_level
      shiny::HTML(sprintf(
        paste0(
          "<b>%s</b><br/>Class: %s<br/>Variables: %d<br/>",
          "Subjects: %d<br/>Confidence level: %s"
        ),
        info$label,
        info$class,
        dims[[1L]],
        info$n_subjects,
        format(signif(conf_level, digits = 3))
      ))
    })

    output$cluster_alert <- shiny::renderUI({
      msg <- filtered_heatmap()$cluster_message
      if (is.null(msg) || identical(msg, "")) {
        return(NULL)
      }
      shiny::div(class = "alert alert-warning", style = "margin-top: 10px;", msg)
    })

    heatmap_plot <- shiny::reactive({
      res <- filtered_heatmap()
      .mc_build_heatmap(
        mat = res$matrix,
        signed = TRUE,
        show_values = isTRUE(input$show_values),
        use_abs = isTRUE(input$use_abs),
        use_plotly = use_plotly
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

    output$pair_plot <- shiny::renderPlot({
      print(plot(pair_fit(), show_legend = isTRUE(input$pair_show_legend)))
    })

    output$pair_meta <- shiny::renderUI({
      fit <- pair_fit()
      shiny::HTML(sprintf(
        paste0(
          "<b>%s vs %s</b><br/>r = %s<br/>Slope = %s<br/>",
          "p-value = %s<br/>df = %s<br/>Subjects = %d<br/>Based on = %d"
        ),
        fit$responses[[1L]],
        fit$responses[[2L]],
        format(signif(fit$r, digits = 4)),
        format(signif(fit$slope, digits = 4)),
        format(signif(fit$p_value, digits = 4)),
        format(signif(fit$df, digits = 4)),
        fit$n_subjects,
        fit$based.on
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
