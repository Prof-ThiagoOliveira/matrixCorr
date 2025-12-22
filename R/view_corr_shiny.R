#' Interactive Shiny viewer for matrixCorr objects
#'
#' Launches an interactive Shiny gadget that displays correlation heatmaps with
#' filtering, clustering, and hover inspection. The viewer accepts any
#' matrixCorr correlation result (for example the outputs from
#' [pearson_corr()], [spearman_rho()], [kendall_tau()], [biweight_mid_corr()],
#' [partial_correlation()], [distance_corr()], or [schafer_corr()]), a plain
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
#'   a searchable picker, while clustering mirrors the package heatmap methods
#'   (complete-linkage on `1 - |r|` for signed matrices or `1 - r` for
#'   non-negative ones).
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

  use_plotly <- requireNamespace("plotly", quietly = TRUE)
  app_title <- title %||% "matrixCorr correlation viewer"

  plot_widget <- if (use_plotly) {
    plotly::plotlyOutput("heatmap", height = "650px")
  } else {
    shiny::plotOutput("heatmap", height = "650px")
  }

  ui <- shiny::fluidPage(
    shiny::titlePanel(app_title),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
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
        shiny::checkboxInput("reorder", "Reorder by hierarchical clustering", value = FALSE),
        shiny::checkboxInput("show_values", "Show cell labels", value = TRUE)
      ),
      shiny::mainPanel(
        plot_widget,
        shiny::htmlOutput("meta")
      )
    )
  )

  server <- function(input, output, session) {
    shiny::observeEvent(input$matrix_choice, {
      info <- prepared[[input$matrix_choice]]
      vars <- colnames(info$matrix)
      shinyWidgets::updatePickerInput(
        session,
        inputId = "var_choice",
        choices = vars,
        selected = .mc_default_vars(info$matrix, default_max_vars)
      )
      rng <- if (info$signed) c(0, 1) else c(0, 1)
      shiny::updateSliderInput(session, "threshold", min = rng[[1L]], max = rng[[2L]])
    }, ignoreNULL = FALSE)

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
      if (isTRUE(input$reorder) && nrow(M) > 2L) {
        M <- .mc_reorder_matrix(M, signed = info$signed)
      }
      thr <- input$threshold %||% 0
      if (thr > 0) {
        suppress <- abs(M) < thr
        M[suppress] <- NA_real_
      }
      if (isTRUE(input$mask_diag)) {
        diag(M) <- NA_real_
      }
      list(matrix = M, meta = info)
    })

    output$meta <- shiny::renderUI({
      info <- filtered_matrix()$meta
      dims <- dim(info$matrix)
      shiny::HTML(sprintf(
        "<b>%s</b><br/>Class: %s<br/>Variables: %d<br/>Signed: %s", 
        info$label, info$class, dims[[1L]], if (info$signed) "yes" else "no"
      ))
    })

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

    if (use_plotly) {
      output$heatmap <- plotly::renderPlotly({
        plot_data()$plot
      })
    } else {
      output$heatmap <- shiny::renderPlot({
        plot_data()$plot
      })
    }
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
  objects <-
    if (is.list(x) && !is.matrix(x) && !inherits(x, "Matrix")) {
      if (is.null(names(x))) {
        names(x) <- paste0("Matrix ", seq_along(x))
      }
      x
    } else {
      list(default = x)
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

.mc_parse_corr_object <- function(obj, label) {
  mat <- NULL
  signed <- TRUE
  desc <- NULL
  cls <- class(obj)
  if (inherits(obj, "partial_corr")) {
    mat <- obj$pcor
    desc <- attr(obj, "method", exact = TRUE) %||% "partial correlation"
  } else if (inherits(obj, "Matrix")) {
    mat <- as.matrix(obj)
  } else if (is.matrix(obj)) {
    mat <- obj
  } else if (is.list(obj) && !is.null(obj$pcor) && is.matrix(obj$pcor)) {
    mat <- obj$pcor
    desc <- "partial correlation"
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

.mc_reorder_matrix <- function(mat, signed) {
  dist_mat <- if (signed) 1 - pmin(1, abs(mat)) else 1 - pmin(1, mat)
  dist_mat <- (dist_mat + t(dist_mat)) / 2
  diag(dist_mat) <- 0
  dist_mat[!is.finite(dist_mat)] <- 0
  dist_obj <- try(suppressWarnings(stats::as.dist(dist_mat)), silent = TRUE)
  if (inherits(dist_obj, "try-error") || is.null(dist_obj)) {
    return(mat)
  }
  ord <- try(stats::hclust(dist_obj, method = "complete")$order, silent = TRUE)
  if (inherits(ord, "try-error")) {
    return(mat)
  }
  mat[ord, ord, drop = FALSE]
}

.mc_build_heatmap <- function(mat, signed, show_values, use_abs, use_plotly) {
  if (use_abs) {
    plot_mat <- abs(mat)
    signed <- FALSE
  } else {
    plot_mat <- mat
  }
  df <- as.data.frame(as.table(plot_mat), stringsAsFactors = FALSE)
  names(df) <- c("Var1", "Var2", "Value")
  df$Var1 <- factor(df$Var1, levels = rev(unique(df$Var1)))
  fill_scale <- if (signed) {
    ggplot2::scale_fill_gradient2(
      low = "indianred1",
      mid = "white",
      high = "steelblue1",
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
  if (use_plotly) {
    plotly::ggplotly(p, tooltip = c("x", "y", "fill"))
  } else {
    p
  }
}
