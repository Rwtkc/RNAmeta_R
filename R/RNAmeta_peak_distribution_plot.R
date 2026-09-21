.rnameta_peak_distribution_table <- function(data, set_group_name = NULL,
                                              features_to_remove = NULL) {
  feature_order <- c("Promoter", "UTR5", "Start Codon", "CDS", "Stop Codon", "UTR3", "Intron", "Intergenic")
  result <- .rnameta_preliminary_component(data, 1L, c("Feature", "Frequency", "Sample"))
  result[, Feature := as.character(Feature)]
  result[, Frequency := suppressWarnings(as.numeric(Frequency))]
  result <- result[Feature %in% setdiff(feature_order, features_to_remove) & is.finite(Frequency)]
  result <- result[, .(Frequency = sum(Frequency)), by = .(Feature, Sample = as.character(Sample))]
  result <- .rnameta_rename_samples(result, set_group_name)
  result[, Feature := factor(Feature, levels = feature_order)]
  result[order(Feature, Sample)]
}

.rnameta_peak_distribution_payload <- function(data, set_group_name = NULL,
                                                features_to_remove = NULL,
                                                show_labels = TRUE,
                                                animate = TRUE) {
  table <- .rnameta_peak_distribution_table(data, set_group_name, features_to_remove)
  categories <- as.character(stats::na.omit(unique(table$Feature)))
  samples <- unique(as.character(table$Sample))
  colors <- .rnameta_chart_palette()
  values <- unlist(table$Frequency, use.names = FALSE)
  ymax <- max(c(1, values), na.rm = TRUE) * 1.16
  list(
    type = "grouped_bar", title = "Peak Distribution", yLabel = "Frequency",
    categories = categories, yDomain = c(0, ymax), yTicks = pretty(c(0, ymax), n = 5),
    showLabels = isTRUE(show_labels), animate = isTRUE(animate), hoverMinHeight = 12,
    series = lapply(seq_along(samples), function(index) {
      sample <- samples[[index]]
      sample_data <- table[Sample == sample]
      sample_values <- sample_data$Frequency[match(categories, as.character(sample_data$Feature))]
      sample_values[is.na(sample_values)] <- 0
      list(name = sample, originalName = sample, color = colors[((index - 1L) %% length(colors)) + 1L], values = unname(sample_values))
    })
  )
}

#' Plot peak distribution with the RNAmeta interactive chart
#' @param preliminary_analysis_data Preliminary analysis output or a prepared table with Feature, Frequency and Sample columns.
#' @param set_group_name Optional replacement sample labels.
#' @param features_to_remove Feature names to exclude.
#' @param show_labels Show frequency labels above bars.
#' @param animate Animate bars in interactive viewers.
#' @param width,height Optional widget dimensions.
#' @return An htmlwidget.
#' @export
RNAmeta_peak_distribution_plot <- function(preliminary_analysis_data = NULL,
                                            set_group_name = NULL,
                                            features_to_remove = NULL,
                                            show_labels = TRUE,
                                            animate = TRUE,
                                            width = NULL,
                                            height = NULL) {
  .rnameta_chart_widget(
    .rnameta_peak_distribution_payload(preliminary_analysis_data, set_group_name, features_to_remove, show_labels, animate),
    width = width, height = height
  )
}

#' Export a peak distribution chart
#' @inheritParams RNAmeta_peak_distribution_plot
#' @param data Preliminary analysis output or prepared table.
#' @param file PNG or PDF output path.
#' @param width Output width before PNG scaling.
#' @param height Optional exact output height.
#' @param aspect_ratio Width/height ratio used when height is NULL.
#' @param scale PNG resolution multiplier.
#' @return The normalized output path, invisibly.
#' @export
RNAmeta_peak_distribution_export <- function(data, file, set_group_name = NULL,
                                              features_to_remove = NULL,
                                              show_labels = TRUE,
                                              width = 1180, height = NULL,
                                              aspect_ratio = c(1180, 560), scale = 2) {
  config <- .rnameta_peak_distribution_payload(data, set_group_name, features_to_remove, show_labels, animate = FALSE)
  .rnameta_chart_export(config, file, width, height, aspect_ratio, scale)
}
