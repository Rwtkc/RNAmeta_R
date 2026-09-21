.rnameta_peak_exon_size_table <- function(data, set_group_name = NULL) {
  result <- .rnameta_preliminary_component(data, 2L, c("Sample", "exonLoc", "length"))
  result <- .rnameta_clean_sample_names(result)
  result[, length := suppressWarnings(as.numeric(length))]
  result <- result[is.finite(length) & length > 0 & as.character(exonLoc) %in% c("First", "Middle", "Last")]
  result <- .rnameta_rename_samples(result, set_group_name)
  result[, exonLoc := factor(as.character(exonLoc), levels = c("First", "Middle", "Last"))]
  result
}

.rnameta_peak_exon_size_payload <- function(data, set_group_name = NULL) {
  table <- .rnameta_peak_exon_size_table(data, set_group_name)
  samples <- unique(as.character(table$Sample)); colors <- .rnameta_chart_palette()
  locations <- c("First", "Middle", "Last"); locations <- locations[locations %in% as.character(table$exonLoc)]
  list(type = "faceted_boxplot", title = "Peak Exon Size", yLabel = "Exon length (bp)", scaleTransform = "log2",
       facets = lapply(locations, function(location) list(name = location, groups = lapply(seq_along(samples), function(i) {
         values <- table[Sample == samples[[i]] & as.character(exonLoc) == location, length]
         list(name = samples[[i]], originalName = samples[[i]], color = colors[((i - 1L) %% length(colors)) + 1L], values = unname(values))
       }))))
}

#' Plot exon lengths by exon position
#' @param preliminary_analysis_data Preliminary analysis output or prepared exon table.
#' @param set_group_name Optional replacement sample labels.
#' @param width,height Optional widget dimensions.
#' @return An htmlwidget.
#' @export
RNAmeta_peak_exon_size_plot <- function(preliminary_analysis_data = NULL, set_group_name = NULL,
                                         width = NULL, height = NULL) {
  .rnameta_chart_widget(.rnameta_peak_exon_size_payload(preliminary_analysis_data, set_group_name), width, height)
}

#' Export a peak exon size chart
#' @param data Preliminary analysis output or prepared exon table.
#' @param file PNG or PDF output path.
#' @inheritParams RNAmeta_peak_exon_size_plot
#' @param width Output width before PNG scaling.
#' @param height Optional exact output height.
#' @param aspect_ratio Width/height ratio used when height is NULL.
#' @param scale PNG resolution multiplier.
#' @return The normalized output path, invisibly.
#' @export
RNAmeta_peak_exon_size_export <- function(data, file, set_group_name = NULL, width = 1180,
                                          height = NULL, aspect_ratio = c(1180, 560), scale = 2) {
  .rnameta_chart_export(.rnameta_peak_exon_size_payload(data, set_group_name), file, width, height, aspect_ratio, scale)
}
