shiny_server <- function(input, output, session, data, qc_detected, replicate_col = NULL, loss_threshold = 25) {
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

  shiny::observe({
    shiny::req(data())
    obs <- data()$obs
    char_cols <- names(obs)[vapply(obs, function(x) is.character(x) || is.factor(x), logical(1))]
    choices <- c("QC Status (Pass/Fail)" = "qc_status")
    for (col in char_cols) {
      choices[[col]] <- col
    }
    shiny::updateSelectInput(session, "color_by", choices = choices)
  })

  shiny::observe({
    shiny::req(data())
    obs <- data()$obs
    rep_candidates <- names(obs)[vapply(obs, function(x) {
      is.character(x) || is.factor(x)
    }, logical(1))]
    choices <- c("None" = "")
    for (col in rep_candidates) {
      choices[[col]] <- col
    }
    shiny::updateSelectInput(session, "replicate_select", choices = choices)
    if (!is.null(replicate_col()) && replicate_col() %in% rep_candidates) {
      shiny::updateSelectInput(session, "replicate_select", selected = replicate_col())
    }
  })

  output$umap_plot <- plotly::renderPlotly({
    shiny::req(data())
    obsm <- data()$obsm
    if (is.null(obsm) || length(obsm) == 0) {
      return(plotly::plot_ly() |>
        plotly::layout(title = "No reduced dimensions available"))
    }

    # Flexible reduction name detection: match by prefix / keyword
    obsm_names <- names(obsm)
    reduction_name <- NULL
    umap_idx <- grep("umap", obsm_names, ignore.case = TRUE)
    pca_idx  <- grep("pca",  obsm_names, ignore.case = TRUE)
    if (length(umap_idx) > 0) {
      reduction_name <- obsm_names[umap_idx[1]]
    } else if (length(pca_idx) > 0) {
      reduction_name <- obsm_names[pca_idx[1]]
    } else {
      reduction_name <- obsm_names[1]
    }

    coords <- obsm[[reduction_name]]
    shiny::req(NCOL(coords) >= 2)

    f <- filtered()
    qc <- qc_metrics()

    plot_df <- data.frame(
      UMAP1 = coords[, 1],
      UMAP2 = coords[, 2],
      cell_barcode = qc$cell_barcode,
      qc_status = ifelse(f$pass, "Pass", "Fail"),
      stringsAsFactors = FALSE
    )

    # Build the color vector and hover text
    if (input$color_by != "qc_status" &&
        input$color_by %in% colnames(data()$obs)) {
      color_vals <- as.character(data()$obs[[input$color_by]])
      color_title <- input$color_by
      plot_df$hover_text <- paste0(
        "Cell: ", plot_df$cell_barcode,
        "<br>", color_title, ": ", color_vals
      )
    } else {
      color_vals <- plot_df$qc_status
      color_title <- "QC Status"
      plot_df$hover_text <- paste0(
        "Cell: ", plot_df$cell_barcode,
        "<br>QC Status: ", color_vals
      )
    }

    plotly::plot_ly(
      data = plot_df,
      x = ~UMAP1,
      y = ~UMAP2,
      color = color_vals,
      colors = if (identical(color_title, "QC Status")) c("#A23B72", "#2E86AB") else NULL,
      type = "scatter",
      mode = "markers",
      marker = list(size = 3),
      text = ~hover_text,
      hoverinfo = "text",
      showlegend = TRUE
    ) |>
      plotly::toWebGL() |>
      plotly::layout(
        title = paste("UMAP \u2014", reduction_name),
        xaxis = list(title = "UMAP 1"),
        yaxis = list(title = "UMAP 2")
      )
  })

  shiny::observe({
    shiny::req(data(), input$replicate_select)
    if (input$replicate_select == "") {
      shiny::updateSelectInput(session, "replicate_filter", choices = c("All" = ""))
      return()
    }
    vals <- unique(data()$obs[[input$replicate_select]])
    choices <- c("All" = "")
    for (v in vals[!is.na(vals)]) {
      choices[[as.character(v)]] <- as.character(v)
    }
    shiny::updateSelectInput(session, "replicate_filter", choices = choices)
  })

  output$impact_cluster <- DT::renderDT({
    shiny::req(data(), filtered())
    obs <- data()$obs
    f <- filtered()

    possible <- c("cluster", "louvain", "leiden", "cell_type", "celltype",
                  "Cluster", "Louvain", "Leiden", "cell_type_main")
    cluster_col <- NULL
    for (col in possible) {
      if (col %in% colnames(obs)) {
        cluster_col <- col
        break
      }
    }

    if (is.null(cluster_col)) {
      return(DT::datatable(data.frame(Message = "No cluster column found in obs"),
                          options = list(dom = "t")))
    }

    if (!is.null(input$replicate_filter) && input$replicate_filter != "") {
      rep_filter <- input$replicate_filter
      rep_col <- input$replicate_select
      in_rep <- obs[[rep_col]] == rep_filter
      groups <- obs[[cluster_col]][in_rep]
      total <- as.numeric(table(groups))
      passing <- as.numeric(by(f$pass[in_rep], groups, sum))
      failing <- total - passing
      pct_pass <- round(passing / total * 100, 1)
    } else {
      groups <- obs[[cluster_col]]
      total <- as.numeric(table(groups))
      passing <- as.numeric(by(f$pass, groups, sum))
      failing <- total - passing
      pct_pass <- round(passing / total * 100, 1)
    }

    impact_df <- data.frame(
      Cluster = names(table(groups)),
      Total = total,
      Passing = passing,
      Failing = failing,
      `% Pass` = pct_pass,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    DT::datatable(impact_df, options = list(pageLength = 10, dom = "tip"),
                  rownames = FALSE) |>
      DT::formatStyle("Failing",
        background = DT::styleColorBar(c(0, max(impact_df$Failing)), "#A23B7244"))
  })

  output$impact_replicate <- DT::renderDT({
    shiny::req(data(), filtered())

    rep_col <- input$replicate_select
    if (is.null(rep_col) || rep_col == "") {
      return(DT::datatable(
        data.frame(Message = "Select a replicate column above to see per-replicate impact"),
        options = list(dom = "t")
      ))
    }

    obs <- data()$obs
    f <- filtered()

    impact <- impact_by_replicate(obs, rep_col, f$pass, loss_threshold = loss_threshold)

    impact$pct_pass <- paste0(impact$pct_pass, "%")
    impact$pct_lost <- paste0(impact$pct_lost, "%")

    impact$warning <- ifelse(impact$warning, "> HIGH LOSS", "OK")

    dt <- DT::datatable(impact,
      options = list(
        pageLength = 10,
        dom = "tip",
        rowCallback = DT::JS(
          "function(row, data) {",
          "  if (data[6] == '> HIGH LOSS') {",
          "    $(row).css('background-color', '#FFE4E1');",
          "  }",
          "}"
        )
      ),
      rownames = FALSE)

    dt
  })

  output$export_csv <- shiny::downloadHandler(
    filename = function() paste0("qc_filtered_cells_", Sys.Date(), ".csv"),
    content = function(file) {
      f <- filtered()
      qc <- qc_metrics()

      export_df <- data()$obs
      export_df$cell_barcode <- qc$cell_barcode
      export_df$pass <- f$pass
      export_df$fail_reason <- f$fail_reason

      col_order <- unique(c("cell_barcode", "pass", "fail_reason", colnames(export_df)))
      export_df <- export_df[, col_order[col_order %in% colnames(export_df)]]

      write.csv(export_df, file, row.names = FALSE)
    }
  )

  output$export_html <- shiny::downloadHandler(
    filename = function() paste0("qc_report_", Sys.Date(), ".html"),
    content = function(file) {
      f <- filtered()
      results <- list(pass = f$pass, fail_reason = f$fail_reason, qc = qc_metrics(), thresholds = thresholds())
      rep_data <- NULL
      if (!is.null(input$replicate_select) && input$replicate_select != "") {
        rep_data <- impact_by_replicate(data()$obs, input$replicate_select, f$pass, loss_threshold)
      }
      write_qc_html(results, file, replicate_data = rep_data)
    }
  )

  # Save filtered H5AD — writes to working directory (not downloadHandler)
  shiny::observeEvent(input$save_h5ad, {
    f <- filtered()
    n_pass <- sum(f$pass)
    if (n_pass == 0) {
      shiny::showNotification("No cells pass the current filters. Nothing to save.",
        type = "error", duration = 5)
      return()
    }

    shiny::withProgress(message = "Saving filtered H5AD", value = 0, {
      lst <- data()
      keep <- f$pass

      shiny::incProgress(0.2, detail = "Subsetting data...")

      filtered_list <- list(
        X = lst$X[, keep, drop = FALSE],
        obs = lst$obs[keep, , drop = FALSE],
        var = lst$var,
        obsm = lapply(lst$obsm, function(m) {
          if (is.matrix(m) || inherits(m, "Matrix")) {
            m[keep, , drop = FALSE]
          } else m
        })
      )

      shiny::incProgress(0.4, detail = "Writing H5AD file...")
      filename <- paste0("qc_filtered_", Sys.Date(), ".h5ad")

      tryCatch(
        {
          schard::write_h5ad(filtered_list, filename)
          shiny::incProgress(0.8, detail = "Finalizing...")
          shiny::showNotification(
            paste0("Saved: ", filename, " (", n_pass, " cells, ",
                   nrow(lst$var), " genes)"),
            type = "message", duration = 5)
        },
        error = function(e) {
          shiny::showNotification(paste0("Error saving H5AD: ", e$message),
            type = "error", duration = 10)
        }
      )

      shiny::incProgress(1.0, detail = "Done!")
    })
  })
}
