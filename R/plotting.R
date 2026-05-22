#' QC plots for batch reports
#'
#' @name qc-plots
#' @keywords internal
NULL

#' Plot UMAP with pass/fail coloring
#'
#' @param obsm Named list of reduced dimensions (from schard list)
#' @param pass Logical vector (TRUE = passes QC)
#' @param reduction Character, which reduction to plot (default "X_umap")
#' @return ggplot2 object
#' @keywords internal
plot_umap_qc <- function(obsm, pass, reduction = "X_umap") {
  coords <- obsm[[reduction]]
  if (is.null(coords) && length(obsm) > 0) {
    coords <- obsm[[1]]
    reduction <- names(obsm)[1]
  }
  if (is.null(coords) || nrow(coords) < 2) {
    return(ggplot2::ggplot() +
      ggplot2::labs(title = "No reduced dimensions available") +
      ggplot2::theme_void())
  }

  df <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    pass = ifelse(pass, "Pass", "Fail")
  )

  ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = pass)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.6) +
    ggplot2::scale_color_manual(values = c("Pass" = "#2E86AB", "Fail" = "#A23B72")) +
    ggplot2::labs(
      title = sprintf("UMAP (%s) \u2014 QC Filtering", reduction),
      x = "UMAP 1", y = "UMAP 2",
      color = "QC Status"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

#' Plot QC metric distributions
#'
#' Returns a list of individual ggplot2 histograms, one per metric.
#'
#' @param qc_metrics data.frame from compute_qc_metrics()
#' @param pass Logical vector (TRUE = passes)
#' @return Named list of ggplot2 objects
#' @keywords internal
plot_qc_distributions <- function(qc_metrics, pass) {
  metrics <- c("pct_mito", "n_genes", "n_UMI", "doublet_score")
  present <- metrics[metrics %in% colnames(qc_metrics)]
  present <- present[!vapply(present, function(m) all(is.na(qc_metrics[[m]])), logical(1))]

  if (length(present) == 0) return(list())

  plots <- list()
  for (m in present) {
    df <- data.frame(
      value = qc_metrics[[m]],
      status = ifelse(pass, "Pass", "Fail")
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = value, fill = status)) +
      ggplot2::geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
      ggplot2::scale_fill_manual(values = c("Pass" = "#2E86AB", "Fail" = "#A23B72")) +
      ggplot2::labs(title = m, x = m, y = "Count") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")
    plots[[m]] <- p
  }
  plots
}

#' Plot impact barplot
#'
#' @param impact data.frame from impact_by_replicate() or per-cluster impact
#' @param group_col Character, column name for grouping (e.g., "replicate")
#' @param title Character, plot title
#' @return ggplot2 object
#' @keywords internal
plot_impact_barplot <- function(impact, group_col = "replicate", title = "Impact by Group") {
  passing_df <- data.frame(
    group = impact[[group_col]],
    count = impact$passing,
    status = "passing",
    stringsAsFactors = FALSE
  )
  failing_df <- data.frame(
    group = impact[[group_col]],
    count = impact$failing,
    status = "failing",
    stringsAsFactors = FALSE
  )
  df_long <- rbind(passing_df, failing_df)

  ggplot2::ggplot(df_long, ggplot2::aes(x = group, y = count, fill = status)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = c("passing" = "#2E86AB", "failing" = "#A23B72")) +
    ggplot2::labs(title = title, x = group_col, y = "Cell count", fill = "Status") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
