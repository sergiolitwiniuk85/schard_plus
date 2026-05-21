#' Run the schardExplorer Shiny app
#'
#' @param data_list Named list with X, obs, var, obsm, uns (from normalize_to_list)
#' @param qc_cols Optional named vector for QC column detection
#' @param replicate_col Optional column name for biological replicates
#' @param loss_threshold Numeric, % loss warning threshold
#'
#' @return The Shiny app object
#' @keywords internal
run_dashboard_app <- function(data_list, qc_cols = NULL, replicate_col = NULL, loss_threshold = 25) {
  qc_detected <- detect_qc_cols(data_list$obs, qc_cols)

  data_reactive <- shiny::reactiveVal(data_list)
  replicate_col_reactive <- shiny::reactiveVal(replicate_col)

  shiny::shinyApp(
    ui = shiny_ui(data_list, qc_detected, replicate_col),
    server = function(input, output, session) {
      shiny_server(input, output, session,
             data = data_reactive,
             qc_detected = qc_detected,
             replicate_col = replicate_col_reactive,
             loss_threshold = loss_threshold)
    }
  )
}
