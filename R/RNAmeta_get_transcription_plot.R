#' Compute transcription boundary and heatmap data
#' @param preliminary_analysis_data Output from RNAmeta_get_preliminary_analysis_data().
#' @param seed Non-negative integer controlling transcript subsampling.
#' @return A list containing boundary hits and one heatmap per sample.
#' @export
get_transcription_data <- function(preliminary_analysis_data = NULL, seed = 1L) {
  seed <- .rnameta_validate_seed(seed)
  txdb_features <- .rnameta_preliminary_value(preliminary_analysis_data, "txdb_features", 7L); peak_gr <- .rnameta_preliminary_value(preliminary_analysis_data, "peak_gr_list", 2L)
  flank_size <- .rnameta_preliminary_value(preliminary_analysis_data, "flanking_size", 6L); gene_max_tx_lengths <- .rnameta_preliminary_value(preliminary_analysis_data, "gene_max_tx_lengths", 8L)
  annotation_key <- paste("transcription-v4", flank_size,
                          .rnameta_granges_fingerprint(txdb_features$cds),
                          .rnameta_object_fingerprint(gene_max_tx_lengths), seed, sep = ":")
  sample_results <- lapply(seq_along(peak_gr), function(i) {
    sample_names <- names(peak_gr); sample_name <- if (is.null(sample_names) || !nzchar(sample_names[[i]])) paste0("Group", i) else sample_names[[i]]
    sample_key <- paste(annotation_key, sample_name, .rnameta_granges_fingerprint(peak_gr[[i]]), sep = ":")
    .rnameta_site_profile_cached(sample_key, function() {
      query <- stats::setNames(list(peak_gr[[i]]), sample_name)
      sample_seed <- .rnameta_seed_offset(seed, i - 1L)
      starts <- .get_feature_boundary_coverage_new(1L, query, txdb_features$cds, flank_size, "Transcription_start_site", 10000, sample_seed)
      ends <- .get_feature_boundary_coverage_new(1L, query, txdb_features$cds, flank_size, "Transcription_end_site", 10000, sample_seed)
      list(coverage = data.table::rbindlist(list(starts, ends), use.names = TRUE, fill = TRUE),
           heatmap = .get_feature_boundary_heatmap(1L, query, txdb_features, gene_max_tx_lengths, flank_size))
    })
  })
  heatmaps <- stats::setNames(lapply(sample_results, `[[`, "heatmap"), names(peak_gr))
  list(data.table::rbindlist(lapply(sample_results, `[[`, "coverage"), use.names = TRUE, fill = TRUE), heatmaps)
}
#' Plot transcription site profiles
#' @param transcription_data Data returned by get_transcription_data().
#' @param set_group_name Optional sample labels.
#' @param heatmap_sample Sample name/index to display, or "all".
#' @param panel Panel name/index to display, or "all" to return separate plots for every panel.
#' @param width,height Optional widget dimensions.
#' @return One selected panel as an htmlwidget, or a named list of separate htmlwidgets when panel is "all".
#' @export
RNAmeta_get_transcription_plot <- function(transcription_data = NULL, set_group_name = NULL, heatmap_sample = "all", width = NULL, height = NULL, panel = "all") {
  .rnameta_site_profile_plot(transcription_data, set_group_name, "Transcription", width, height, heatmap_sample, panel)
}

#' Export transcription site profiles
#' @param data Data returned by get_transcription_data().
#' @param file PNG or PDF path.
#' @param set_group_name Optional sample labels.
#' @param heatmap_sample Sample name/index to export, or "all".
#' @param panel Panel name/index to export, or "all" to write one file per panel.
#' @param width Output width.
#' @param height Optional exact height.
#' @param aspect_ratio Width/height ratio when height is NULL.
#' @param scale PNG resolution multiplier.
#' @return One output path, or named output paths for all panels, invisibly.
#' @export
RNAmeta_transcription_export <- function(data, file, set_group_name = NULL, heatmap_sample = "all", width = 1180, height = NULL, aspect_ratio = c(1180, 560), scale = 2, panel = "all") {
  .rnameta_site_profile_export(data, file, set_group_name, "Transcription", width, height, aspect_ratio, scale, heatmap_sample, panel)
}
