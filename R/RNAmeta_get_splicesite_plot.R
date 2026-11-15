.rnameta_splicesite_data <- function(preliminary_analysis_data) {
  .rnameta_online_profile_data(preliminary_analysis_data, "splicesite")
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
