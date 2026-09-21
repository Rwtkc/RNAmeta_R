#' Compute translation boundary and heatmap data
#' @param preliminary_analysis_data Output from RNAmeta_get_preliminary_analysis_data().
#' @return A list containing boundary hits and one heatmap per sample.
#' @export
get_translation_data <- function(preliminary_analysis_data = NULL) {
  txdb_features <- .rnameta_preliminary_value(preliminary_analysis_data, "txdb_features", 7L); peak_gr <- .rnameta_preliminary_value(preliminary_analysis_data, "peak_gr_list", 2L)
  flank_size <- .rnameta_preliminary_value(preliminary_analysis_data, "short_flanking_size", 9L); tx_lengths <- .rnameta_preliminary_value(preliminary_analysis_data, "tx_lengths", 10L)
  annotation_key <- paste("translation-v4", flank_size,
                          .rnameta_granges_fingerprint(txdb_features$exons),
                          .rnameta_object_fingerprint(tx_lengths), sep = ":")
  sample_results <- lapply(seq_along(peak_gr), function(i) {
    sample_names <- names(peak_gr); sample_name <- if (is.null(sample_names) || !nzchar(sample_names[[i]])) paste0("Group", i) else sample_names[[i]]
    sample_key <- paste(annotation_key, sample_name, .rnameta_granges_fingerprint(peak_gr[[i]]), sep = ":")
    .rnameta_site_profile_cached(sample_key, function() {
      query <- stats::setNames(list(peak_gr[[i]]), sample_name)
      list(coverage = .get_translation_boundary_coverage_new(1L, query, tx_lengths, txdb_features$exons, flank_size),
           heatmap = .get_translation_boundary_heatmap(1L, query, tx_lengths, txdb_features$exons, flank_size))
    })
  })
  heatmaps <- stats::setNames(lapply(sample_results, `[[`, "heatmap"), names(peak_gr))
  list(data.table::rbindlist(lapply(sample_results, `[[`, "coverage"), use.names = TRUE, fill = TRUE), heatmaps)
}
#' Plot translation site profiles
#' @param translation_data Data returned by get_translation_data().
#' @param set_group_name Optional sample labels.
#' @param heatmap_sample Sample name/index to display, or "all".
#' @param panel Panel name/index to display, or "all" to return separate plots for every panel.
#' @param width,height Optional widget dimensions.
#' @return One selected panel as an htmlwidget, or a named list of separate htmlwidgets when panel is "all".
#' @export
RNAmeta_get_translation_plot <- function(translation_data = NULL, set_group_name = NULL, heatmap_sample = "all", width = NULL, height = NULL, panel = "all") {
  .rnameta_site_profile_plot(translation_data, set_group_name, "Translation", width, height, heatmap_sample, panel)
}

#' Export translation site profiles
#' @param data Data returned by get_translation_data().
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
RNAmeta_translation_export <- function(data, file, set_group_name = NULL, heatmap_sample = "all", width = 1180, height = NULL, aspect_ratio = c(1180, 560), scale = 2, panel = "all") {
  .rnameta_site_profile_export(data, file, set_group_name, "Translation", width, height, aspect_ratio, scale, heatmap_sample, panel)
}
