.rnameta_splicesite_data <- function(preliminary_analysis_data) {
  peaks <- .rnameta_preliminary_value(preliminary_analysis_data, "peak_gr_list", 2L); exons <- .rnameta_preliminary_value(preliminary_analysis_data, "txdb_features", 7L)$exons; flank <- .rnameta_preliminary_value(preliminary_analysis_data, "short_flanking_size", 9L)
  annotation_key <- paste("splicesite-v4", flank, .rnameta_granges_fingerprint(exons), sep = ":")
  sample_results <- lapply(seq_along(peaks), function(i) {
    sample_names <- names(peaks); sample_name <- if (is.null(sample_names) || !nzchar(sample_names[[i]])) paste0("Group", i) else sample_names[[i]]
    sample_key <- paste(annotation_key, sample_name, .rnameta_granges_fingerprint(peaks[[i]]), sep = ":")
    .rnameta_site_profile_cached(sample_key, function() {
      query <- stats::setNames(list(peaks[[i]]), sample_name)
      .get_eejunct_new(1L, query, exons, flank)
    })
  })
  data.table::rbindlist(sample_results, use.names = TRUE, fill = TRUE)
}
#' Plot splice-site profiles
#' @param preliminary_analysis_data Preliminary analysis output or a prepared hit table.
#' @param set_group_name Optional sample labels.
#' @param panel Panel name/index to display, or "all" to return separate plots for every panel.
#' @param width,height Optional widget dimensions.
#' @return One selected panel as an htmlwidget, or a named list of separate htmlwidgets when panel is "all".
#' @export
RNAmeta_get_splicesite_plot <- function(preliminary_analysis_data = NULL, set_group_name = NULL, width = NULL, height = NULL, panel = "all") {
  hits <- if (is.data.frame(preliminary_analysis_data)) preliminary_analysis_data else .rnameta_splicesite_data(preliminary_analysis_data)
  .rnameta_site_profile_plot(list(hits), set_group_name, "Splicesite", width, height, panel = panel)
}

#' Export splice-site profiles
#' @param data Preliminary analysis output or prepared hit table.
#' @param file PNG or PDF path.
#' @param set_group_name Optional sample labels.
#' @param panel Panel name/index to export, or "all" to write one file per panel.
#' @param width Output width.
#' @param height Optional exact height.
#' @param aspect_ratio Width/height ratio when height is NULL.
#' @param scale PNG resolution multiplier.
#' @return One output path, or named output paths for all panels, invisibly.
#' @export
RNAmeta_splicesite_export <- function(data, file, set_group_name = NULL, width = 1180, height = NULL, aspect_ratio = c(1180, 560), scale = 2, panel = "all") {
  hits <- if (is.data.frame(data)) data else .rnameta_splicesite_data(data)
  .rnameta_site_profile_export(list(hits), file, set_group_name, "Splicesite", width, height, aspect_ratio, scale, panel = panel)
}
