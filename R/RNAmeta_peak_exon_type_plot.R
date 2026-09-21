.rnameta_peak_exon_type_table <- function(data, set_group_name = NULL) {
  input <- .rnameta_preliminary_component(data, 2L, c("Sample", "exonLoc"))
  input <- .rnameta_clean_sample_names(input)
  if (all(c("Count", "Percentage") %in% names(input))) {
    result <- input[, .(Count = sum(as.numeric(Count)), Percentage = sum(as.numeric(Percentage))), by = .(Sample, exonLoc = as.character(exonLoc))]
  } else {
    result <- input[, .(Count = .N), by = .(Sample, exonLoc = as.character(exonLoc))]
    result[, Percentage := if (sum(Count) > 0) Count / sum(Count) * 100 else 0, by = Sample]
  }
  result <- .rnameta_rename_samples(result, set_group_name)
  result[, exonLoc := factor(exonLoc, levels = c("First", "Middle", "Last"))]
  result[order(Sample, exonLoc)]
}

.rnameta_peak_exon_type_payload <- function(data, set_group_name = NULL, animate = TRUE) {
  table <- .rnameta_peak_exon_type_table(data, set_group_name)
  samples <- unique(as.character(table$Sample))
  locations <- c("First", "Middle", "Last")
  locations <- locations[locations %in% as.character(table$exonLoc)]
  palette <- c(First = "#a30545", Middle = "#f26e42", Last = "#7ecba4")
  list(
    type = "stacked_bar", title = "Peak Exon Type", yLabel = "Percentage",
    categories = samples, yDomain = c(0, 100), yTicks = seq(0, 100, 20),
    showLabels = FALSE, animate = isTRUE(animate), tooltipValueLabel = "Percentage",
    series = lapply(locations, function(location) {
      current <- table[as.character(exonLoc) == location]
      values <- current$Percentage[match(samples, as.character(current$Sample))]
      values[is.na(values)] <- 0
      list(name = location, color = unname(palette[[location]]), values = unname(values))
    })
  )
}

#' Plot peak exon types with an interactive stacked bar chart
#' @param preliminary_analysis_data Preliminary analysis output or prepared exon data.
#' @param set_group_name Optional replacement sample labels.
#' @param animate Animate bars in interactive viewers.
#' @param width,height Optional widget dimensions.
#' @return An htmlwidget.
#' @export
RNAmeta_peak_exon_type_plot <- function(preliminary_analysis_data = NULL, set_group_name = NULL,
                                         animate = TRUE, width = NULL, height = NULL) {
  .rnameta_chart_widget(.rnameta_peak_exon_type_payload(preliminary_analysis_data, set_group_name, animate), width, height)
}

#' Export a peak exon type chart
#' @param data Preliminary analysis output or prepared exon data.
#' @param file PNG or PDF output path.
#' @inheritParams RNAmeta_peak_exon_type_plot
#' @param width Output width before PNG scaling.
#' @param height Optional exact output height.
#' @param aspect_ratio Width/height ratio used when height is NULL.
#' @param scale PNG resolution multiplier.
#' @return The normalized output path, invisibly.
#' @export
RNAmeta_peak_exon_type_export <- function(data, file, set_group_name = NULL, width = 1180,
                                          height = NULL, aspect_ratio = c(1180, 560), scale = 2) {
  .rnameta_chart_export(.rnameta_peak_exon_type_payload(data, set_group_name, FALSE), file, width, height, aspect_ratio, scale)
}
