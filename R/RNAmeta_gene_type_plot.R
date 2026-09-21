#' Gene Type plot
#'
#' Render Gene Type counts as the interactive RNAmeta grouped-bar widget.
#' The returned object works in RStudio Viewer, R Markdown, Quarto and Shiny.
#' On headless systems use `RNAmeta_gene_type_export()` instead of relying on a Viewer.
#' @param preliminary_analysis_data A legacy preliminary analysis list, or the
#'   counts table returned by `RNAmeta_gene_type_data()`.
#' @param set_group_name Optional replacement labels in current sample order.
#' @param show_labels Display count labels above bars.
#' @param animate Animate bars when rendered.
#' @param width,height Widget dimensions.
#' @param elementId Optional HTML element id.
#' @return An interactive htmlwidget.
#' @export
RNAmeta_gene_type_plot <- function(preliminary_analysis_data = NULL,
                                   set_group_name = NULL,
                                   show_labels = TRUE,
                                   animate = TRUE,
                                   width = NULL,
                                   height = 560,
                                   elementId = NULL) {
  .rnameta_gene_type_widget(
    preliminary_analysis_data,
    set_group_name = set_group_name,
    show_labels = show_labels,
    animate = animate,
    width = width,
    height = height,
    elementId = elementId
  )
}