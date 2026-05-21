#' Generate a batch QC report
#'
#' Runs QC filtering and generates static plots + CSV summary.
#' Does NOT require Shiny \u2014 runs in any R environment.
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

  results <- list(pass = filtered$pass, fail_reason = filtered$fail_reason, qc = qc)
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

#' Write a simple self-contained QC HTML report
#'
#' @param results List with elements: pass (logical), fail_reason (character), qc (data.frame)
#' @param output_path Character, path to write the HTML file
#' @return Invisibly returns output_path
#' @keywords internal
write_qc_html <- function(results, output_path, replicate_data = NULL) {
  n_total <- length(results$pass)
  n_pass <- sum(results$pass)
  n_fail <- n_total - n_pass
  pct_pass <- round(n_pass / n_total * 100, 1)
  pct_fail <- round(n_fail / n_total * 100, 1)

  rep_rows <- ""
  if (!is.null(replicate_data) && nrow(replicate_data) > 0) {
    for (i in seq_len(nrow(replicate_data))) {
      r <- replicate_data[i, ]
      cl <- if (isTRUE(r$warning)) ' class="warn"' else ""
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

  html <- sprintf('<!DOCTYPE html>
<html>
<head><title>QC Report</title>
<style>
body { font-family: -apple-system, sans-serif; max-width: 960px; margin: 2em auto; }
table { border-collapse: collapse; width: 100%%; }
th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
th { background: #f5f5f5; }
.pass { color: #2E86AB; }
.fail { color: #A23B72; }
.warn { background: #FFE4E1; }
</style></head>
<body>
<h1>QC Report</h1>
<p>Generated: %s</p>
<h2>Summary</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Total cells</td><td>%d</td></tr>
<tr><td>Passing</td><td class="pass">%d (%.1f%%%%)</td></tr>
<tr><td>Failing</td><td class="fail">%d (%.1f%%%%)</td></tr>
</table>
%s
</body></html>',
    Sys.time(), n_total, n_pass, pct_pass, n_fail, pct_fail, rep_section)

  writeLines(html, output_path)
  invisible(output_path)
}
