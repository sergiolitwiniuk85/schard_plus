shiny_ui <- function(data, qc_detected, replicate_col) {
  bslib::page_sidebar(
    title = "schardExplorer \u2014 Interactive QC Dashboard",
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),

    sidebar = bslib::sidebar(
      width = 350,
      bslib::card(
        bslib::card_header("QC Thresholds"),
        shiny::uiOutput("qc_sliders_ui")
      ),
      bslib::card(
        bslib::card_header("Display Options"),
        shiny::selectInput("color_by", "Color UMAP by:",
                    choices = c("QC Status (Pass/Fail)" = "qc_status")),
        shiny::selectInput("replicate_select", "Biological Replicate:",
                    choices = c("None" = ""),
                    selected = replicate_col %||% ""),
        shiny::selectInput("replicate_filter", "Filter cluster impact by replicate:",
                    choices = c("All" = "")),
        shiny::hr(),
        shiny::downloadButton("export_csv", "Export CSV"),
        shiny::downloadButton("export_html", "Export HTML Report"),
        shiny::hr(),
        shiny::actionButton("save_h5ad", "Save Filtered H5AD",
          class = "btn-primary",
          icon = shiny::icon("download"),
          onclick = "this.disabled=true; this.innerHTML='Saving...';"
        )
      )
    ),

    bslib::card(
      bslib::card_header("UMAP"),
      plotly::plotlyOutput("umap_plot", height = "500px")
    ),

    bslib::layout_columns(
      bslib::card(
        bslib::card_header("Impact by Cluster"),
        DT::DTOutput("impact_cluster")
      ),
      bslib::card(
        bslib::card_header("Impact by Replicate"),
        DT::DTOutput("impact_replicate")
      )
    ),

    bslib::card(
      bslib::card_header("QC Summary"),
      shiny::textOutput("qc_summary_text")
    )
  )
}
