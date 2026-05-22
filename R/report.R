#' Generate a batch QC report
#'
#' Runs QC filtering and generates static plots + CSV summary.
#' Does NOT require Shiny - runs in any R environment.
#'
#' @param object Single-cell object: list, SingleCellExperiment, Seurat, or .h5ad path.
#' @param output_dir Directory to write report files.
#' @param thresholds List of QC thresholds (see apply_filter_thresholds).
#' @param replicate_col Optional column name for biological replicates.
#' @param qc_cols Optional named vector for QC column detection.
#' @param loss_threshold Numeric, % loss warning threshold for replicates.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns the filtered object list.
#' @export
qc_report <- function(object, ...) {
  UseMethod("qc_report")
}

#' @rdname qc_report
#' @export
qc_report.default <- function(object, output_dir = ".",
                               thresholds = list(pct_mito = 10),
                               replicate_col = NULL,
                               qc_cols = NULL,
                               loss_threshold = 25,
                               ...) {
  data <- normalize_to_list(object)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop("Cannot create output directory: ", output_dir, call. = FALSE)
  }

  qc <- compute_qc_metrics(data$obs, qc_cols)
  filtered <- apply_filter_thresholds(data$obs, thresholds, qc_cols)
  n_total <- length(filtered$pass)
  n_pass <- sum(filtered$pass)

  cells_out <- data.frame(
    cell_barcode = qc$cell_barcode,
    pass = filtered$pass,
    fail_reason = filtered$fail_reason,
    qc[, !colnames(qc) %in% "cell_barcode", drop = FALSE],
    stringsAsFactors = FALSE
  )
  write.csv(cells_out, file.path(output_dir, "filtered_cells.csv"), row.names = FALSE)

  if (!is.null(data$obsm) && length(data$obsm) > 0) {
    p <- plot_umap_qc(data$obsm, filtered$pass)
    ggplot2::ggsave(file.path(output_dir, "umap_qc.png"), p, width = 8, height = 6, dpi = 150)
  }

  qc_plots <- plot_qc_distributions(qc, filtered$pass)
  if (length(qc_plots) > 0) {
    for (nm in names(qc_plots)) {
      ggplot2::ggsave(file.path(output_dir, paste0("qc_", nm, ".png")),
                      qc_plots[[nm]], width = 6, height = 4, dpi = 150)
    }
  }

  if (!is.null(replicate_col)) {
    impact <- impact_by_replicate(data$obs, replicate_col, filtered$pass, loss_threshold)
    write.csv(impact, file.path(output_dir, "impact_by_replicate.csv"), row.names = FALSE)

    p <- plot_impact_barplot(impact, "replicate", "Impact by Replicate")
    ggplot2::ggsave(file.path(output_dir, "impact_by_replicate.png"), p, width = 8, height = 5, dpi = 150)
  }

  results <- list(pass = filtered$pass, fail_reason = filtered$fail_reason, qc = qc, thresholds = thresholds)
  write_qc_html(results, file.path(output_dir, "qc_summary.html"))

  message(sprintf("QC Report: %d/%d cells pass (%.1f%%) \u2192 %s",
                  n_pass, n_total, n_pass / n_total * 100, output_dir))

  data$pass <- filtered$pass
  data$fail_reason <- filtered$fail_reason
  invisible(data)
}

#' @rdname qc_report
#' @export
qc_report.character <- function(object, output_dir = ".", thresholds = list(pct_mito = 10), replicate_col = NULL, qc_cols = NULL, loss_threshold = 25, ...) {
  NextMethod()
}

#' @rdname qc_report
#' @export
qc_report.list <- function(object, output_dir = ".", thresholds = list(pct_mito = 10), replicate_col = NULL, qc_cols = NULL, loss_threshold = 25, ...) {
  NextMethod()
}

#' @rdname qc_report
#' @export
qc_report.SingleCellExperiment <- function(object, output_dir = ".", thresholds = list(pct_mito = 10), replicate_col = NULL, qc_cols = NULL, loss_threshold = 25, ...) {
  NextMethod()
}

#' @rdname qc_report
#' @export
qc_report.Seurat <- function(object, output_dir = ".", thresholds = list(pct_mito = 10), replicate_col = NULL, qc_cols = NULL, loss_threshold = 25, ...) {
  NextMethod()
}

#' Write a self-contained HTML QC report
#'
#' @param results list with pass, fail_reason, qc, thresholds
#' @param output_path character, file path to write
#' @param replicate_data optional data.frame from impact_by_replicate
#' @return invisible output_path
#' @keywords internal
write_qc_html <- function(results, output_path, replicate_data = NULL) {
  n_total <- length(results$pass)
  n_pass <- sum(results$pass)
  n_fail <- n_total - n_pass
  pct_pass <- round(n_pass / n_total * 100, 1)
  pct_fail <- round(n_fail / n_total * 100, 1)

  qc_summary_rows <- ""
  qc_metrics <- c("pct_mito", "n_genes", "n_UMI", "doublet_score")
  for (m in qc_metrics) {
    if (!is.null(results$qc[[m]]) && !all(is.na(results$qc[[m]]))) {
      vals <- results$qc[[m]][!is.na(results$qc[[m]])]
      qc_summary_rows <- paste0(qc_summary_rows, sprintf(
        '<tr><td>%s</td><td>%.1f</td><td>%.1f</td><td>%.1f</td><td>%.1f</td></tr>\n',
        m, min(vals), mean(vals), stats::median(vals), max(vals)
      ))
    }
  }

  rep_rows <- ""
  if (!is.null(replicate_data) && nrow(replicate_data) > 0) {
    for (i in seq_len(nrow(replicate_data))) {
      r <- replicate_data[i, ]
      cl <- if (isTRUE(r$warning)) ' class="warning"' else ""
      rep_rows <- paste0(rep_rows, sprintf(
        '<tr%s><td>%s</td><td>%d</td><td>%d</td><td>%d</td><td>%.1f%%</td><td>%.1f%%</td></tr>\n',
        cl, r$replicate, r$total, r$passing, r$failing, r$pct_pass, r$pct_lost
      ))
    }
  }

  rep_section <- if (nzchar(rep_rows)) sprintf('
<h2>Impact by Replicate</h2>
<table>
<tr><th>Replicate</th><th>Total</th><th>Passing</th><th>Failing</th><th>%% Pass</th><th>%% Lost</th></tr>
%s
</table>', rep_rows) else ""

  thresholds <- results$thresholds
  params_rows <- ""
  if (!is.null(thresholds)) {
    for (nm in names(thresholds)) {
      v <- thresholds[[nm]]
      params_rows <- paste0(params_rows, sprintf(
        '<tr><td>%s</td><td>%s</td></tr>\n',
        nm, paste(v, collapse = " \u2013 ")
      ))
    }
  }
  params_section <- if (nzchar(params_rows)) sprintf('
<h2>Filtering Parameters</h2>
<table>
<tr><th>Parameter</th><th>Value</th></tr>
%s
</table>', params_rows) else ""

  html <- sprintf('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>QC Report \u2014 schardExplorer</title>
<style>
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; max-width: 960px; margin: 2em auto; padding: 0 1em; color: #333; }
h1 { color: #2E86AB; border-bottom: 2px solid #2E86AB; padding-bottom: 0.3em; }
h2 { color: #555; margin-top: 1.5em; }
table { border-collapse: collapse; width: 100%%; margin: 1em 0; }
th, td { border: 1px solid #ddd; padding: 8px 12px; text-align: left; }
th { background: #f5f5f5; font-weight: 600; }
tr:nth-child(even) { background: #fafafa; }
.pass { color: #2E86AB; font-weight: 600; }
.fail { color: #A23B72; font-weight: 600; }
.warning { background: #FFE4E1 !important; }
.summary-box { background: #f0f7fa; border-left: 4px solid #2E86AB; padding: 1em; margin: 1em 0; border-radius: 4px; }
.footer { margin-top: 2em; font-size: 0.85em; color: #999; border-top: 1px solid #eee; padding-top: 1em; }
</style>
</head>
<body>
<h1>schardExplorer QC Report</h1>
<p>Generated: %s</p>
<div class="summary-box">
<strong>Summary:</strong> <span class="pass">%d (%.1f%%) passing</span> | <span class="fail">%d (%.1f%%) failing</span> | %d total cells
</div>
<h2>QC Metric Summary</h2>
<table>
<tr><th>Metric</th><th>Min</th><th>Mean</th><th>Median</th><th>Max</th></tr>
%s
</table>
%s
%s
<div class="footer">
<p>Generated by schardExplorer v0.1.0</p>
</div>
</body>
</html>',
    Sys.time(), n_pass, pct_pass, n_fail, pct_fail, n_total,
    qc_summary_rows,
    rep_section,
    params_section)

  writeLines(html, output_path)
  invisible(output_path)
}
