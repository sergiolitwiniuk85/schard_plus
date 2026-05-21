shiny_server <- function(input, output, session, data, qc_detected) {
  qc_metrics <- shiny::reactive({
    compute_qc_metrics(data()$obs, qc_detected)
  })

  output$qc_sliders_ui <- shiny::renderUI({
    shiny::req(qc_metrics())
    qc <- qc_metrics()

    sliders <- list()

    if (!all(is.na(qc$pct_mito))) {
      rng <- range(qc$pct_mito, na.rm = TRUE)
      sliders[[length(sliders) + 1]] <- shiny::sliderInput(
        "thresh_pct_mito", "Max Mito %",
        min = floor(rng[1]), max = ceiling(rng[2]),
        value = min(10, ceiling(rng[2])),
        step = 0.5
      )
    }

    if (!all(is.na(qc$n_genes))) {
      rng <- range(qc$n_genes, na.rm = TRUE)
      sliders[[length(sliders) + 1]] <- shiny::sliderInput(
        "thresh_n_genes", "Gene Count Range",
        min = floor(rng[1]), max = ceiling(rng[2]),
        value = c(max(500, floor(rng[1])), min(6000, ceiling(rng[2]))),
        step = 100
      )
    }

    if (!all(is.na(qc$n_UMI))) {
      rng <- range(qc$n_UMI, na.rm = TRUE)
      sliders[[length(sliders) + 1]] <- shiny::sliderInput(
        "thresh_n_UMI", "UMI Count Range",
        min = floor(rng[1]), max = ceiling(rng[2]),
        value = c(max(1000, floor(rng[1])), ceiling(rng[2])),
        step = 100
      )
    }

    if (!all(is.na(qc$doublet_score))) {
      rng <- range(qc$doublet_score, na.rm = TRUE)
      sliders[[length(sliders) + 1]] <- shiny::sliderInput(
        "thresh_doublet_score", "Max Doublet Score",
        min = 0, max = ceiling(rng[2] * 10) / 10,
        value = min(0.25, ceiling(rng[2] * 10) / 10),
        step = 0.05
      )
    }

    if (length(sliders) == 0) {
      return(shiny::p("No QC metrics detected. Use 'qc_cols' parameter to specify."))
    }

    sliders
  })

  thresholds <- shiny::reactive({
    list(
      pct_mito = input$thresh_pct_mito,
      n_genes = if (!is.null(input$thresh_n_genes)) input$thresh_n_genes else c(0, Inf),
      n_UMI = if (!is.null(input$thresh_n_UMI)) input$thresh_n_UMI else c(0, Inf),
      doublet_score = input$thresh_doublet_score
    )
  })

  filtered <- shiny::reactive({
    apply_filter_thresholds(qc_metrics(), thresholds())
  })

  output$qc_summary_text <- shiny::renderText({
    f <- filtered()
    n_total <- length(f$pass)
    n_pass <- sum(f$pass)
    sprintf("Passing: %d / %d (%.1f%%)", n_pass, n_total, n_pass / n_total * 100)
  })

  output$umap_plot <- plotly::renderPlotly({
    plotly::plot_ly() |>
      plotly::layout(title = "UMAP plot will render in Phase 6")
  })

  output$impact_cluster <- DT::renderDT({
    data.frame(Message = "Cluster impact table \u2014 Phase 6")
  })

  output$impact_replicate <- DT::renderDT({
    data.frame(Message = "Replicate impact table \u2014 Phase 7")
  })

  output$export_csv <- shiny::downloadHandler(
    filename = function() paste0("qc_filtered_cells_", Sys.Date(), ".csv"),
    content = function(file) {
      f <- filtered()
      write.csv(data.frame(
        cell_barcode = qc_metrics()$cell_barcode,
        pass = f$pass,
        fail_reason = f$fail_reason
      ), file, row.names = FALSE)
    }
  )

  output$export_html <- shiny::downloadHandler(
    filename = function() paste0("qc_report_", Sys.Date(), ".html"),
    content = function(file) {
      f <- filtered()
      results <- list(pass = f$pass, fail_reason = f$fail_reason, qc = qc_metrics())
      write_qc_html(results, file)
    }
  )
}
