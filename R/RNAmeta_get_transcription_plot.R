#' Compute transcription boundary and heatmap data
#' @param preliminary_analysis_data Output from RNAmeta_get_preliminary_analysis_data().
#' @param seed Non-negative integer controlling transcript subsampling.
#' @return A list containing boundary hits and one heatmap per sample.
#' @export
get_transcription_data <- function(preliminary_analysis_data = NULL, seed = 1L) {
  .rnameta_validate_seed(seed)
  .rnameta_online_profile_data(preliminary_analysis_data, "transcription")
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
