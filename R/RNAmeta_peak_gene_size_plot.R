.rnameta_peak_gene_size_table <- function(data, set_group_name = NULL) {
  shared <- attr(data, "peak_annotation", exact = TRUE)
  if (!is.null(shared$raw_peaks)) {
    result <- data.table::rbindlist(lapply(seq_along(shared$raw_peaks), function(i) {
      .rnameta_online_peak_gene_size_compute_lengths(shared$raw_peaks[[i]],
        shared$annotation_bundle$txdb_features$transcripts, shared$annotation_bundle$txlens,
        shared$samples[[i]])
    }), use.names = TRUE, fill = TRUE)
  } else {
    result <- .rnameta_preliminary_component(data, 5L, c("Sample", "tx_len"))
  }
  result <- .rnameta_clean_sample_names(result)
  result[, tx_len := suppressWarnings(as.numeric(tx_len))]
  result <- result[is.finite(tx_len) & tx_len > 0]
  result <- .rnameta_rename_samples(result, set_group_name)
  result
}

.rnameta_peak_gene_size_payload <- function(data, set_group_name = NULL) {
  table <- .rnameta_peak_gene_size_table(data, set_group_name)
  samples <- unique(as.character(table$Sample)); colors <- .rnameta_chart_palette()
  list(type = "boxplot", title = "Peak Gene Size", yLabel = "Transcript length (bp)", scaleTransform = "log2",
       groups = lapply(seq_along(samples), function(i) list(name = samples[[i]], originalName = samples[[i]], color = colors[((i - 1L) %% length(colors)) + 1L], values = unname(table[Sample == samples[[i]], tx_len]))))
}

#' Plot transcript lengths for peak-associated genes
#' @param preliminary_analysis_data Preliminary analysis output or prepared transcript table.
#' @param set_group_name Optional replacement sample labels.
#' @param width,height Optional widget dimensions.
#' @return An htmlwidget.
#' @export
RNAmeta_peak_gene_size_plot <- function(preliminary_analysis_data = NULL, set_group_name = NULL,
                                         width = NULL, height = NULL) {
  .rnameta_chart_widget(.rnameta_peak_gene_size_payload(preliminary_analysis_data, set_group_name), width, height)
}

#' Export a peak gene size chart
#' @param data Preliminary analysis output or prepared transcript table.
#' @param file PNG or PDF output path.
#' @inheritParams RNAmeta_peak_gene_size_plot
#' @param width Output width before PNG scaling.
#' @param height Optional exact output height.
#' @param aspect_ratio Width/height ratio used when height is NULL.
#' @param scale PNG resolution multiplier.
#' @return The normalized output path, invisibly.
#' @export
RNAmeta_peak_gene_size_export <- function(data, file, set_group_name = NULL, width = 1180,
                                          height = NULL, aspect_ratio = c(1180, 560), scale = 2) {
  .rnameta_chart_export(.rnameta_peak_gene_size_payload(data, set_group_name), file, width, height, aspect_ratio, scale)
}
